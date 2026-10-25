"""
EDNiX Longitudinal Tools
========================
Dataset-agnostic functions for longitudinal morphometry analysis.

The only constraint on input data is a standard DataFrame format:

    subject  | session | age    | region  | hemisphere | atlas_level | <metric_col>
    ---------|---------|--------|---------|------------|-------------|-------------
    sub-01   | 1       | 25.3   | Isocortex | left     | 1           | 1234.5

Any dataset can be adapted to this format in a separate formatting script
before calling these functions.

Functions
---------
fit_gamm(df, metric_col, age_col)       → dict with spline-LMM fit results (historical name)
fit_lme(df, metric_col, age_col)        → dict with LME fit results
draw_longitudinal(ax, ...)               → trajectory panel into an existing axis
draw_pct_change(ax, ...)                 → % change panel into an existing axis
plot_longitudinal(...)                   → trajectory figure (one PNG)
plot_pct_change(...)                     → % change figure (one PNG)
prepare_longitudinal(...)                → filter + age merge + min-sessions
run_longitudinal(...)                    → main runner for one species/dataset
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── statsmodels ──────────────────────────────────────────────────────────────
try:
    import statsmodels.formula.api as smf
    from statsmodels.tools.sm_exceptions import ConvergenceWarning
    _HAS_SMF = True
except ImportError:
    _HAS_SMF = False
    warnings.warn("statsmodels not found — LME disabled. pip install statsmodels")

# ── pygam ────────────────────────────────────────────────────────────────────
try:
    from pygam import LinearGAM, s
    _HAS_GAM = True
except ImportError:
    _HAS_GAM = False
    # pygam is not used by fit_gamm (spline LMM, statsmodels/patsy only)

try:
    from Exemples.Study_EDNiX.analysis.ednix_bids_tools import PAPER_RC
except ImportError:          # the core package must not require the study scripts
    PAPER_RC = {}

opj = os.path.join

# ─────────────────────────────────────────────────────────────────────────────
# CONFIG (can be overridden at call site)
# ─────────────────────────────────────────────────────────────────────────────
MIN_SESSIONS  = 3    # minimum sessions per subject to include in model
GAM_N_SPLINES = None      # None = automatic (capped at 8 for short series)
GAM_LAM       = 5.0       # base penalty (higher = smoother)
import numpy as _np
GAM_LAM_GRID  = _np.logspace(0.5, 3.5, 12)  # search 3.2 … 3162 (smoother end)


# ══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def _safe_col(name: str) -> str:
    """Make column name safe for statsmodels formula."""
    import re
    return re.sub(r'[^a-zA-Z0-9_]', '_', f"v_{name}")


def _filter_df(df: pd.DataFrame,
               atlas_level: int = 1,
               hemisphere: str = 'left',
               regions: list = None) -> pd.DataFrame:
    """Apply atlas_level / hemisphere / region filters."""
    df = df.copy()
    if 'atlas_level' in df.columns:
        df = df[df['atlas_level'] == atlas_level]
    if hemisphere and hemisphere != 'all' and 'hemisphere' in df.columns:
        filt = df[df['hemisphere'] == hemisphere]
        df   = filt if not filt.empty else df
    if regions and 'region' in df.columns:
        df = df[df['region'].isin(regions)]
    return df


def _keep_longitudinal(df: pd.DataFrame,
                       age_col: str,
                       min_sessions: int = MIN_SESSIONS) -> pd.DataFrame:
    """Drop subjects with fewer than min_sessions valid age values."""
    n   = df.groupby('subject')[age_col].count()
    ok  = n[n >= min_sessions].index
    out = df[df['subject'].isin(ok)]
    dropped = df['subject'].nunique() - out['subject'].nunique()
    if dropped:
        print(f"    dropped {dropped} subject(s) with < {min_sessions} sessions")
    return out


# ══════════════════════════════════════════════════════════════════════════════
# MODELS
# ══════════════════════════════════════════════════════════════════════════════
def _ml_random_intercept(y, X, groups):
    """
    Exact ML fit of y = X b + u_subject + e (random intercept), by profiling the
    variance ratio lam = var(u)/var(e) on [0, 1e4] (boundary lam = 0 included,
    where statsmodels.MixedLM fails to invert the Hessian).
    Returns (llf, beta, cov_beta, description).
    """
    from scipy.optimize import minimize_scalar
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    _, inv = np.unique(groups, return_inverse=True)
    blocks = [np.flatnonzero(inv == g) for g in range(inv.max() + 1)]
    n, k = X.shape

    def parts(lam):
        XtViX = np.zeros((k, k)); XtViy = np.zeros(k); yViy = 0.0; logdet = 0.0
        for idx in blocks:
            Xg, yg, ng = X[idx], y[idx], len(idx)
            c = lam / (1.0 + lam * ng)
            sx, sy = Xg.sum(0), yg.sum()
            XtViX += Xg.T @ Xg - c * np.outer(sx, sx)
            XtViy += Xg.T @ yg - c * sx * sy
            yViy += yg @ yg - c * sy * sy
            logdet += np.log1p(lam * ng)
        beta = np.linalg.solve(XtViX, XtViy)
        sigma2 = max((yViy - beta @ XtViy) / n, 1e-300)
        llf = -0.5 * (n * np.log(2 * np.pi * sigma2) + logdet + n)
        return llf, beta, sigma2, XtViX

    res = minimize_scalar(lambda lam: -parts(lam)[0], bounds=(0.0, 1e4), method='bounded',
                          options={'xatol': 1e-8})
    lam = res.x if parts(res.x)[0] >= parts(0.0)[0] else 0.0
    llf, beta, sigma2, XtViX = parts(lam)
    cov = sigma2 * np.linalg.inv(XtViX)
    return float(llf), beta, cov, f'LMM random intercept (ML, var ratio {lam:.3g})'


def fit_gamm(df: pd.DataFrame,
             metric_col: str,
             age_col: str,
             min_sessions: int = MIN_SESSIONS,
             spline_df: int = 4) -> dict:
    """
    Mixed model with an unpenalised cubic B-spline age basis
    (metric ~ bs(age, df=spline_df) + (1|subject)), fitted by ML.
    NB: despite the historical name, this is not a penalised GAMM.

    Age effect: likelihood-ratio test against the intercept-only model
    (chi2, spline_df degrees of freedom). The previous version returned the
    minimum of the individual coefficient p-values (random-effect variance term
    included), which gave ~19% false positives at alpha = 0.05 under the null.

    Non-linearity: LRT spline vs linear-age model (same random intercept),
    spline_df - 1 degrees of freedom (p_nonlinear).

    Predictions use the spline basis of the fitted data (same knots) and the
    full fixed-effect covariance for the confidence band.
    """
    if not _HAS_SMF:
        warnings.warn("statsmodels not available - spline LMM disabled")
        return None

    from patsy import dmatrix, build_design_matrices
    from scipy.stats import chi2

    sub = df[[age_col, metric_col, 'subject']].dropna()
    if len(sub) < min_sessions:
        return None

    y = sub[metric_col].to_numpy(dtype=float)
    groups = sub['subject'].to_numpy()
    age_mu, age_sd = sub[age_col].mean(), sub[age_col].std() + 1e-9
    age_z = pd.DataFrame({'age_z': (sub[age_col].to_numpy(dtype=float) - age_mu) / age_sd})
    grid_z = np.linspace(age_z['age_z'].min(), age_z['age_z'].max(), 200)

    try:
        X = dmatrix(f"bs(age_z, df={spline_df}, degree=3)", age_z, return_type='dataframe')
        Xg = np.asarray(build_design_matrices([X.design_info], pd.DataFrame({'age_z': grid_z}))[0])
        ll1, beta, cov, kind = _ml_random_intercept(y, X.to_numpy(), groups)
        ll0, _, _, _ = _ml_random_intercept(y, np.ones((len(y), 1)), groups)
        ll_lin, _, _, _ = _ml_random_intercept(
            y, np.column_stack([np.ones(len(y)), age_z['age_z'].to_numpy()]), groups)
        dof = spline_df
        method = f'{kind.split(" (")[0]}, cubic B-spline age basis (df={spline_df}), LRT vs intercept-only'
    except Exception as e:
        warnings.warn(f"spline model failed — {metric_col}: {e}; linear model used")
        X = np.column_stack([np.ones(len(y)), age_z['age_z'].to_numpy()])
        Xg = np.column_stack([np.ones_like(grid_z), grid_z])
        ll1, beta, cov, kind = _ml_random_intercept(y, X, groups)
        ll0, _, _, _ = _ml_random_intercept(y, X[:, :1], groups)
        ll_lin = ll1
        dof = 1
        method = f'{kind.split(" (")[0]}, linear age (fallback), LRT vs intercept-only'

    lr = max(2.0 * (ll1 - ll0), 0.0)
    lr_nl = max(2.0 * (ll1 - ll_lin), 0.0)   # spline vs linear age (non-linearity), dof-1 df
    pred = Xg @ beta
    se = np.sqrt(np.clip(np.einsum('ij,jk,ik->i', Xg, cov, Xg), 0, None))
    return dict(
        age_grid=grid_z * age_sd + age_mu, pred=pred,
        ci_lower=pred - 1.96 * se, ci_upper=pred + 1.96 * se,
        p_age=float(chi2.sf(lr, dof)), lr_stat=lr,
        p_nonlinear=float(chi2.sf(lr_nl, dof - 1)) if dof > 1 else np.nan, lr_nonlinear=lr_nl,
        edf=float(dof), method=method, n=len(sub),
    )

def fit_lme(df: pd.DataFrame,
            metric_col: str,
            age_col: str,
            min_sessions: int = MIN_SESSIONS) -> dict:
    """
    Linear mixed model: metric ~ age_z + (1|subject).

    Returns
    -------
    dict with keys: slope (per SD of age), slope_per_unit (per age unit, e.g. per
    year), se_per_unit, p_slope, intercept, r2_marginal, n
    """
    if not _HAS_SMF:
        return {}
    sub = df[[age_col, metric_col, 'subject']].dropna()
    if len(sub) < min_sessions:
        return {}

    sc  = _safe_col(metric_col)
    tmp = sub.rename(columns={metric_col: sc})
    tmp['age_z'] = (tmp[age_col] - tmp[age_col].mean()) / (tmp[age_col].std() + 1e-9)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', ConvergenceWarning)
            md = smf.mixedlm(
                f"{sc} ~ age_z", tmp, groups=tmp['subject'],
            ).fit(reml=False, method='lbfgs', maxiter=300)

        slope = md.params.get('age_z', np.nan)
        pval  = md.pvalues.get('age_z', np.nan)
        age_sd = tmp[age_col].std() + 1e-9
        se_z  = md.bse.get('age_z', np.nan)
        inter = md.params.get('Intercept', np.nan)
        # Marginal R² (Nakagawa & Schielzeth 2013)
        # variance of the fixed-effect prediction X*beta (fittedvalues - resid = 2*fitted - y was used before)
        var_fe  = float(np.var(md.model.exog @ np.asarray(md.fe_params)))
        var_re  = float(md.cov_re.values[0][0]) if hasattr(md, 'cov_re') else 0.0
        var_res = float(md.scale)
        r2_m    = var_fe / (var_fe + var_re + var_res + 1e-9)
        return dict(slope=slope, slope_per_unit=slope / age_sd, se_per_unit=se_z / age_sd,
                    p_slope=pval, intercept=inter, r2_marginal=r2_m, n=len(sub))
    except Exception as e:
        warnings.warn(f"fit_lme error — {metric_col}: {e}")
        return {}


# ══════════════════════════════════════════════════════════════════════════════
# PLOTS
# ══════════════════════════════════════════════════════════════════════════════

def _fmt_p(p):
    return 'p<0.001' if p < 0.001 else f'p={p:.3f}'


def _unit_of(label):
    """'Volume (mm³)' -> 'mm³'; '' if no unit in parentheses."""
    return label.split('(')[-1].rstrip(')') if label and '(' in label else ''


def draw_longitudinal(ax,
                      df: pd.DataFrame,
                      metric_col: str,
                      age_col: str,
                      gamm_result: dict = None,
                      lme_result:  dict = None,
                      age_label: str = 'Age (years)',
                      metric_label: str = None,
                      stats_fontsize: float = 8,
                      inset: str = 'full'):
    """
    Draw individual trajectories + spline-LMM fit + 95% CI into `ax`.
    inset: 'full'  = spline-LMM LRT, n, linear-LME slope/p, R²m (per-panel PNGs)
           'slope' = paper figure: linear-LME slope per age unit and its p, or, when
                     the spline departs from linearity (p_nonlinear < 0.05), the
                     non-linearity LRT p instead (a linear slope would misdescribe it)
           None    = no text
    Subject colours: tab20 over the sorted subject list (same as draw_pct_change).
    """
    sub_df = df[[age_col, metric_col, 'subject']].dropna()
    if sub_df.empty:
        return False

    subjects = sorted(sub_df['subject'].unique())
    colors   = plt.cm.tab20(np.linspace(0, 1, max(len(subjects), 1)))
    sp_col   = dict(zip(subjects, colors))

    for subj in subjects:
        sdf = sub_df[sub_df['subject'] == subj].sort_values(age_col)
        kw  = dict(color=sp_col[subj], alpha=0.55, zorder=4)
        if len(sdf) < 2:
            ax.scatter(sdf[age_col], sdf[metric_col], s=30, **kw)
        else:
            ax.plot(sdf[age_col], sdf[metric_col],
                    'o-', markersize=4, linewidth=1, **kw)

    if gamm_result:
        ax.plot(gamm_result['age_grid'], gamm_result['pred'],
                color='black', lw=2.5, zorder=6, label='Spline LMM fit')
        ax.fill_between(gamm_result['age_grid'],
                        gamm_result['ci_lower'], gamm_result['ci_upper'],
                        color='black', alpha=0.15, zorder=5)

    lines = []
    if inset == 'full':
        if gamm_result:
            lines += [f"Spline LMM (LRT)  {_fmt_p(gamm_result['p_age'])}  df={gamm_result['edf']:.0f}",
                      f"n={gamm_result['n']} sessions, {sub_df['subject'].nunique()} subjects"]
        if lme_result:
            lines += [f"Linear LME slope={lme_result.get('slope', np.nan):.3f}/SD  "
                      f"{_fmt_p(lme_result.get('p_slope', np.nan))}",
                      f"R²_m={lme_result.get('r2_marginal', np.nan):.3f}"]
    elif (inset == 'slope' and gamm_result
          and gamm_result.get('p_nonlinear', np.nan) < 0.05):
        lines += ['non-linear', f"{_fmt_p(gamm_result['p_nonlinear'])} (LRT)"]
    elif inset == 'slope' and lme_result and 'slope_per_unit' in lme_result:
        unit = _unit_of(metric_label)
        t_unit = _unit_of(age_label).rstrip('s') or 'unit'   # years -> year
        v = lme_result['slope_per_unit']
        v_str = f"{v:.0f}" if abs(v) >= 100 else f"{v:.3g}"
        lines += [f"slope = {v_str} {unit}/{t_unit}",
                  _fmt_p(lme_result.get('p_slope', np.nan))]
    if lines:
        ax.text(0.03, 0.97, '\n'.join(lines), transform=ax.transAxes,
                va='top', fontsize=stats_fontsize, zorder=10,
                bbox=dict(facecolor='white', alpha=0.85, edgecolor='none'))

    ax.set_xlabel(age_label)
    ax.set_ylabel(metric_label or metric_col.replace('_', ' ').capitalize())
    return True


def plot_longitudinal(df: pd.DataFrame,
                      metric_col: str,
                      age_col: str,
                      region: str,
                      species: str,
                      output_path: str,
                      gamm_result: dict = None,
                      lme_result:  dict = None,
                      age_label: str = 'Age (years)',
                      metric_label: str = None):
    """
    Individual trajectories + spline-LMM fit + 95% CI + stats inset (one PNG).
    """
    if df[[age_col, metric_col, 'subject']].dropna().empty:
        return
    with plt.rc_context(PAPER_RC):
        fig, ax = plt.subplots(figsize=(8, 5))
        draw_longitudinal(ax, df, metric_col, age_col, gamm_result, lme_result,
                          age_label, metric_label)
        ax.xaxis.label.set_size(10)
        ax.yaxis.label.set_size(10)
        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, dpi=200, bbox_inches='tight')
        plt.close(fig)


def pct_change_table(df: pd.DataFrame,
                     metric_col: str,
                     age_col: str) -> pd.DataFrame:
    """
    % change from each subject's first session (subjects with >= 2 sessions).
    Adds columns pct_change and age_from_start.
    """
    sub_df = df[[age_col, metric_col, 'subject']].dropna().copy()
    pct_rows = []
    for subj, grp in sub_df.groupby('subject'):
        grp      = grp.sort_values(age_col)
        baseline = grp[metric_col].iloc[0]
        if abs(baseline) < 1e-9 or len(grp) < 2:
            continue
        grp = grp.copy()
        grp['pct_change']     = 100.0 * (grp[metric_col] - baseline) / abs(baseline)
        grp['age_from_start'] = grp[age_col] - grp[age_col].iloc[0]
        pct_rows.append(grp)
    return pd.concat(pct_rows, ignore_index=True) if pct_rows else pd.DataFrame()


def draw_pct_change(ax,
                    df: pd.DataFrame,
                    metric_col: str,
                    age_col: str,
                    age_label: str = 'Age (years)',
                    min_n_mean: int = 3):
    """
    Draw % change from each subject's first session into `ax`:
    individual curves + mean ± SEM on a common time-from-baseline grid.
    Mean/SEM are drawn only where >= min_n_mean subjects are observed
    (the tail of the grid otherwise rests on 1–2 subjects).
    """
    pct_df = pct_change_table(df, metric_col, age_col)
    if pct_df.empty:
        return False

    # same subject -> colour mapping as draw_longitudinal (all subjects of df)
    all_subjects = sorted(df[[age_col, metric_col, 'subject']].dropna()['subject'].unique())
    colors   = plt.cm.tab20(np.linspace(0, 1, max(len(all_subjects), 1)))
    sp_col   = dict(zip(all_subjects, colors))
    subjects = sorted(pct_df['subject'].unique())

    for subj in subjects:
        sdf = pct_df[pct_df['subject'] == subj].sort_values('age_from_start')
        ax.plot(sdf['age_from_start'], sdf['pct_change'],
                'o-', color=sp_col[subj], alpha=0.5, markersize=3.5, linewidth=1)

    # Mean ± SEM across subjects on a common grid
    grid = np.linspace(0, pct_df['age_from_start'].max(), 100)
    mat  = []
    for subj in subjects:
        sdf = pct_df[pct_df['subject'] == subj].sort_values('age_from_start')
        if len(sdf) < 2:
            continue
        mat.append(np.interp(grid, sdf['age_from_start'].values,
                             sdf['pct_change'].values, left=np.nan, right=np.nan))

    if len(mat) >= 2:
        mat    = np.array(mat)
        mean_y = np.nanmean(mat, axis=0)
        n_ok   = np.sum(~np.isnan(mat), axis=0)
        sem_y  = np.nanstd(mat, axis=0) / np.sqrt(np.maximum(n_ok, 1))
        low_n  = n_ok < min_n_mean
        mean_y[low_n] = np.nan
        sem_y[low_n]  = np.nan
        ax.plot(grid, mean_y, color='black', lw=2.5, zorder=6, label='Mean')
        ax.fill_between(grid, mean_y - sem_y, mean_y + sem_y,
                        color='black', alpha=0.2, zorder=5)

    ax.axhline(0, color='gray', lw=0.8, linestyle='--')
    unit = age_label.split('(')[-1].rstrip(')') if '(' in age_label else age_label
    ax.set_xlabel(f'Time from baseline ({unit})')
    ax.set_ylabel('% change from baseline')
    return True


def plot_pct_change(df: pd.DataFrame,
                    metric_col: str,
                    age_col: str,
                    region: str,
                    species: str,
                    output_path: str,
                    age_label: str = 'Age (years)'):
    """
    % change from each subject's first session (one PNG).
    Shows individual curves + mean ± SEM on a common time-from-baseline grid.
    """
    if pct_change_table(df, metric_col, age_col).empty:
        return
    with plt.rc_context(PAPER_RC):
        fig, ax = plt.subplots(figsize=(8, 5))
        draw_pct_change(ax, df, metric_col, age_col, age_label)
        ax.xaxis.label.set_size(10)
        ax.yaxis.label.set_size(10)
        ax.set_title(f'{species} — {region} — % change',
                     fontweight='bold', fontsize=11)
        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, dpi=200, bbox_inches='tight')
        plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# MAIN RUNNER
# ══════════════════════════════════════════════════════════════════════════════

def prepare_longitudinal(df_morph:     pd.DataFrame,
                         df_age:       pd.DataFrame,
                         age_col:      str,
                         atlas_level:  int  = 1,
                         hemisphere:   str  = 'left',
                         regions:      list = None,
                         min_sessions: int  = MIN_SESSIONS,
                         label:        str  = '') -> pd.DataFrame:
    """
    Filter (atlas level / hemisphere / regions), merge with the age table on
    subject (+ session when both tables have it), keep subjects with
    >= min_sessions sessions. Returns an empty DataFrame when nothing is left.
    """
    df = _filter_df(df_morph, atlas_level=atlas_level,
                    hemisphere=hemisphere, regions=regions)
    if df.empty:
        print(f"  [{label}] empty after filters — skipping")
        return pd.DataFrame()

    merge_on = ['subject']
    if 'session' in df_age.columns and 'session' in df.columns:
        merge_on.append('session')

    df['subject'] = df['subject'].astype(str)
    df_age        = df_age.copy()
    df_age['subject'] = df_age['subject'].astype(str)
    if 'session' in df_age.columns:
        df_age['session'] = df_age['session'].astype(str)
    if 'session' in df.columns:
        df['session'] = df['session'].astype(str)

    merged = pd.merge(df, df_age, on=merge_on, how='inner')
    if merged.empty:
        print(f"  [{label}] no rows after merge with age table")
        return pd.DataFrame()

    print(f"  [{label}] {len(merged)} rows  "
          f"{merged['subject'].nunique()} subjects  "
          f"{merged['region'].nunique() if 'region' in merged.columns else '?'} regions")

    merged = _keep_longitudinal(merged, age_col, min_sessions)
    if merged.empty:
        print(f"  [{label}] no subjects with ≥{min_sessions} sessions")
    return merged


def run_longitudinal(
    label:        str,
    df_morph:     pd.DataFrame,
    df_age:       pd.DataFrame,
    metric_col:   str,
    metric_label: str,
    age_col:      str,
    age_label:    str,
    out_dir:      str,
    atlas_level:  int   = 1,
    hemisphere:   str   = 'left',
    regions:      list  = None,
    min_sessions: int   = MIN_SESSIONS,
) -> pd.DataFrame:
    """
    Run longitudinal analysis for ONE species / modality.

    Parameters
    ----------
    label        : display name (e.g. 'Macaque_surface')
    df_morph     : tidy morphometry DataFrame from ednix_bids_tools
                   Required columns: subject, session, region, hemisphere,
                                     atlas_level, <metric_col>
    df_age       : standardised age table
                   Required columns: subject, <age_col>
                   Optional: session (for session-level age)
    metric_col   : column name of the morphometry value
    metric_label : y-axis label for plots
    age_col      : column name for age in df_age
    age_label    : x-axis label for plots
    out_dir      : output directory (created if needed)
    atlas_level  : which atlas hierarchy level to keep
    hemisphere   : 'left', 'right', 'bilateral', or 'all'
    regions      : list of region names to process (None = all)
    min_sessions : minimum sessions per subject

    Returns
    -------
    DataFrame with one row per region containing spline-LMM and linear-LME stats.
    """
    merged = prepare_longitudinal(df_morph, df_age, age_col, atlas_level=atlas_level,
                                  hemisphere=hemisphere, regions=regions,
                                  min_sessions=min_sessions, label=label)
    if merged.empty:
        return pd.DataFrame()

    os.makedirs(out_dir, exist_ok=True)
    all_regions = sorted(merged['region'].unique()) \
                  if 'region' in merged.columns else [label]
    stats_rows  = []

    for region in all_regions:
        rdf = merged[merged['region'] == region].copy() \
              if 'region' in merged.columns else merged.copy()

        if rdf[metric_col].dropna().empty:
            continue

        gamm = fit_gamm(rdf, metric_col, age_col, min_sessions)
        lme  = fit_lme(rdf,  metric_col, age_col, min_sessions)

        row = dict(label=label, region=region,
                   n_subjects=rdf['subject'].nunique(),
                   n_sessions=len(rdf))
        if gamm:
            row.update(gamm_p=gamm['p_age'], gamm_edf=gamm['edf'],
                       gamm_method=gamm['method'])
        if lme:
            row.update(lme_slope=lme.get('slope'),
                       lme_p=lme.get('p_slope'),
                       lme_r2=lme.get('r2_marginal'))
        stats_rows.append(row)

        safe = region.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')
        plot_longitudinal(
            rdf, metric_col, age_col, region, label,
            opj(out_dir, f'{safe}_trajectory.png'),
            gamm_result=gamm, lme_result=lme,
            age_label=age_label, metric_label=metric_label)

        plot_pct_change(
            rdf, metric_col, age_col, region, label,
            opj(out_dir, f'{safe}_pct_change.png'),
            age_label=age_label)

    df_stats = pd.DataFrame(stats_rows)
    if not df_stats.empty:
        df_stats.to_excel(opj(out_dir, 'stats.xlsx'), index=False)
        print(f"  [{label}] stats → {out_dir}/stats.xlsx")

    return df_stats