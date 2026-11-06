"""
longitudinal_ablation.py -- Supplementary Figure 15: what does the longitudinal mode buy?
=======================================================================================

Question (Results, "Tracking development and ageing within individuals", Fig. 3)
    The anatomical pipeline has a `coregistration_longitudinal` option
    (modalities/anat/_loop2.py, loop3/backtonative.py). With it, every session of a
    subject is registered to that subject's LATEST session (the within-subject
    reference), and only the reference is registered to the template; without it,
    every session is registered to the template on its own. Each dataset of Fig. 3
    was processed twice, once per mode, on the same images:

        Human    /scratch2/EDNiX/Human/BIDS_Park            longitudinal  (Fig. 3a-c)
                 /scratch2/EDNiX/Human/BIDS_Park_nolong     cross-sectional
        Macaque  /scratch2/EDNiX/Macaque/BIDS_Zhu_Garin     longitudinal  (Fig. 3d-i)
                 /scratch2/EDNiX/Macaque/*Zhu_Garin*non_long* cross-sectional

    Expected (Reuter et al. 2012 NeuroImage 61:1402): a shared path to group space
    removes registration noise between sessions -> smaller within-subject residual
    variance, smaller standard error of the age slope (= more power), and no change
    in the slope itself (no bias). EDNiX's reference is the LAST session, not an
    unbiased median template: asymmetric processing can bias change estimates
    (Reuter & Fischl 2011 NeuroImage 57:19), so the slope difference is tested too.

Unit of analysis and data (identical to Fig. 3, library/longitudinal.py)
    Isocortex, EDNiX atlas level 1, LEFT hemisphere, subjects with >= 2 sessions,
    ages from the ORIGINAL dataset (participants.tsv of BIDS_Park; Garin_macaque.xlsx).
    Only (subject, session) pairs present in BOTH pipelines are compared, so that the
    two modes are evaluated on exactly the same scans.

Statistics, per species x modality (surface area, thickness, volume)
    1. Integrity: the reference session is processed identically in both modes, so its
       values must agree (median |difference| ~ 0). A large difference there means
       the two folders were not processed with the same settings/atlas.
    2. Agreement of the two modes over all sessions: ICC(A,1) and Bland-Altman bias
       (% of the mean).
    3. Same linear mixed model as the Fig. 3 inset (metric ~ age + (1|subject), ML):
       slope per year, SE, P, within-subject residual SD (sigma_e) and between-subject
       SD (sigma_u), in each mode. Non-linearity LRT (spline vs linear, as Fig. 3).
    4. Paired cluster bootstrap over SUBJECTS (same resampled subjects in both modes,
       N_BOOT resamples): 95% CI and two-sided P of
         sigma_e ratio   cross-sectional / longitudinal   (> 1: longitudinal less noisy)
         SE ratio        cross-sectional / longitudinal   (> 1: longitudinal more precise)
         delta slope     longitudinal - cross-sectional   (0: no bias)
       Sample-size equivalence = SE ratio^2: number of subjects the cross-sectional mode
       would need to match the precision of the longitudinal mode, relative to n.
    5. Model-free check: per subject with >= 3 sessions, SD of the residuals around its
       own least-squares line; paired Wilcoxon signed-rank, longitudinal vs
       cross-sectional.
    6. Fig. 3g-i (macaque) / Supp. Fig. 7 (human): association vs primary paired t-test,
       recomputed in each mode (t, P, d_z).

Outputs (paper_stats.supp15 / paper_figures.supp15, run_paper.py supp15)
    <PAPER_DIR>/Supp_Fig_15/Supp_Fig_15.pdf|png
    <PAPER_DIR>/Supp_Fig_15/Supp_Fig_15_stats.xlsx   sheets summary, assoc_primary, paired_values
    <PAPER_DIR>/Supp_Fig_15/Supp_Fig_15_stats.txt
"""
import glob
import os
import warnings

import numpy as np
import pandas as pd
from scipy import stats

opj = os.path.join

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════
# One entry per species. "nolong" may be a path or a glob pattern; the pattern must
# match exactly one folder (find_folder raises otherwise, so a wrong folder is never
# used silently).
PAIRS = [
    dict(species='Human',
         long='/scratch2/EDNiX/Human/BIDS_Park',
         nolong='/scratch2/EDNiX/Human/BIDS_Park*no*long*'),
    dict(species='Macaque',
         long='/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin',
         nolong='/scratch2/EDNiX/Macaque/*Zhu_Garin*no*long*'),
]
REGION = 'Isocortex'     # Fig. 3 a-f
N_BOOT = 2000            # subject bootstrap resamples
SEED = 0
MODE_LONG, MODE_CROSS = 'longitudinal', 'cross-sectional'


def find_folder(pattern):
    """Path or glob -> the single matching folder (error if none or several)."""
    hits = sorted(p for p in glob.glob(pattern) if os.path.isdir(p))
    if len(hits) != 1:
        raise FileNotFoundError(
            f"{pattern!r} must match exactly one folder, found {hits or 'none'}; "
            "edit PAIRS in library/longitudinal_ablation.py")
    return hits[0]


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 -- load both modes with the Fig. 3 loaders
# ══════════════════════════════════════════════════════════════════════════════
def load_ablation_data(pairs=PAIRS):
    """[(species, {mode: morphometry dict}, age table, age column, {mode: folder})].
    Morphometry is read with the very function that feeds Fig. 3
    (longitudinal.load_morphometry); the age table comes from the longitudinal
    folder (the cross-sectional copy may lack participants.tsv), same subjects."""
    from Exemples.Study_EDNiX.analysis.library import longitudinal as lg
    rows = []
    for p in pairs:
        sp = p['species']
        folders = {MODE_LONG: find_folder(p['long']), MODE_CROSS: find_folder(p['nolong'])}
        atlas_paths, _ = lg._build_atlas_paths([sp])
        morph = {mode: lg.load_morphometry({sp: [d]}, atlas_paths, regions_of_interest=lg.REGIONS)
                 for mode, d in folders.items()}
        if sp == 'Macaque':
            age, age_col = lg.format_macaque_age(lg.MACAQUE_EXCEL), 'age_years'
        else:
            age, age_col = lg.format_human_age([folders[MODE_LONG]],
                                               age_prefix=lg.HUMAN_AGE_PREFIX), 'age'
        rows.append((sp, morph, age, age_col, folders))
    return rows


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 -- matched table: one row per (subject, session) present in both modes
# ══════════════════════════════════════════════════════════════════════════════
def matched_sessions(morph, age_df, age_col, mod, metric_col, region=REGION,
                     atlas_level=1, hemisphere='left', min_sessions=2):
    """Columns subject, session, <age_col>, long, cross, is_reference.
    Each mode is filtered exactly as Fig. 3 (prepare_longitudinal), then the two are
    inner-joined so that only scans present in both are kept. is_reference = latest
    session of the subject (the longitudinal reference, max_ses in _loop2.py)."""
    from Statistics.Group_fMRI.EDNiX_longitudinal import prepare_longitudinal
    per_mode = {}
    for mode in (MODE_LONG, MODE_CROSS):
        df = morph[mode].get(mod, pd.DataFrame())
        if df is None or df.empty:
            return pd.DataFrame()
        m = prepare_longitudinal(df, age_df, age_col, atlas_level=atlas_level,
                                 hemisphere=hemisphere, regions=[region],
                                 min_sessions=min_sessions, label=f'S15/{mode}/{mod}')
        if m.empty:
            return pd.DataFrame()
        m = m[m['region'] == region]
        # one value per scan (a duplicated row would weigh a session twice)
        per_mode[mode] = (m.groupby(['subject', 'session'], as_index=False)
                           .agg(**{age_col: (age_col, 'first'), mode: (metric_col, 'mean')}))
    out = per_mode[MODE_LONG].merge(per_mode[MODE_CROSS].drop(columns=age_col),
                                    on=['subject', 'session'], how='inner')
    out = out.dropna(subset=[MODE_LONG, MODE_CROSS, age_col])
    # keep subjects that still have >= min_sessions after matching
    n = out.groupby('subject')['session'].transform('size')
    out = out[n >= min_sessions].copy()
    # reference session = latest session number (sessions are numbered in time order)
    sess_num = pd.to_numeric(out['session'].astype(str).str.replace('ses-', ''), errors='coerce')
    out['_s'] = sess_num
    out['is_reference'] = out['_s'] == out.groupby('subject')['_s'].transform('max')
    return out.drop(columns='_s').rename(columns={MODE_LONG: 'long', MODE_CROSS: 'cross'})


# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 -- the linear mixed model of the Fig. 3 inset, fast and exact
# ══════════════════════════════════════════════════════════════════════════════
def lmm_random_intercept(y, age, groups):
    """ML fit of y = b0 + b1 * age + u_subject + e, u ~ N(0, s_u^2), e ~ N(0, s_e^2).
    Same model as Statistics/Group_fMRI/EDNiX_longitudinal.fit_lme (statsmodels MixedLM,
    reml=False); solved here by profiling lam = s_u^2 / s_e^2 (closed-form GLS for
    each lam), which is exact, handles the lam = 0 boundary and is ~100x faster, so it
    can be refitted thousands of times in the bootstrap. Age is in years, so b1 is the
    slope per year. Returns dict(slope, se, p, sigma_e, sigma_u, n_obs, n_subjects)."""
    from scipy.optimize import minimize_scalar
    y = np.asarray(y, float)
    X = np.column_stack([np.ones(len(y)), np.asarray(age, float)])
    _, inv = np.unique(np.asarray(groups), return_inverse=True)
    blocks = [np.flatnonzero(inv == g) for g in range(inv.max() + 1)]
    n = len(y)

    def parts(lam):
        # V = s_e^2 (I + lam J) per subject; Sherman-Morrison: V^-1 = (I - c J) / s_e^2
        XtViX, XtViy, yViy, logdet = np.zeros((2, 2)), np.zeros(2), 0.0, 0.0
        for idx in blocks:
            Xg, yg, ng = X[idx], y[idx], len(idx)
            c = lam / (1.0 + lam * ng)
            sx, sy = Xg.sum(0), yg.sum()
            XtViX += Xg.T @ Xg - c * np.outer(sx, sx)
            XtViy += Xg.T @ yg - c * sx * sy
            yViy += yg @ yg - c * sy * sy
            logdet += np.log1p(lam * ng)
        beta = np.linalg.solve(XtViX, XtViy)
        s2 = max((yViy - beta @ XtViy) / n, 1e-300)       # ML estimate of s_e^2
        llf = -0.5 * (n * np.log(2 * np.pi * s2) + logdet + n)
        return llf, beta, s2, XtViX

    res = minimize_scalar(lambda lam: -parts(lam)[0], bounds=(0.0, 1e4), method='bounded',
                          options={'xatol': 1e-8})
    lam = res.x if parts(res.x)[0] >= parts(0.0)[0] else 0.0
    _, beta, s2, XtViX = parts(lam)
    se = float(np.sqrt(s2 * np.linalg.inv(XtViX)[1, 1]))
    z = beta[1] / se if se > 0 else np.nan
    return dict(slope=float(beta[1]), intercept=float(beta[0]), se=se,
                p=float(2 * stats.norm.sf(abs(z))),            # Wald z, as statsmodels
                sigma_e=float(np.sqrt(s2)), sigma_u=float(np.sqrt(lam * s2)),
                n_obs=n, n_subjects=len(blocks))


# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 -- agreement and model-free statistics
# ══════════════════════════════════════════════════════════════════════════════
def icc_a1(x, y):
    """ICC(A,1), two-way, absolute agreement, single measure (McGraw & Wong 1996)."""
    Y = np.column_stack([x, y]).astype(float)
    n, k = Y.shape
    gm = Y.mean()
    msr = k * ((Y.mean(1) - gm) ** 2).sum() / (n - 1)
    msc = n * ((Y.mean(0) - gm) ** 2).sum() / (k - 1)
    mse = ((Y - Y.mean(1, keepdims=True) - Y.mean(0, keepdims=True) + gm) ** 2).sum() \
        / ((n - 1) * (k - 1))
    return float((msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n))


def within_subject_residual_sd(m, age_col, col):
    """{subject: SD of residuals around the subject's own OLS line} (>= 3 sessions)."""
    out = {}
    for s, g in m.groupby('subject'):
        if len(g) < 3:
            continue
        t, v = g[age_col].to_numpy(float), g[col].to_numpy(float)
        res = v - np.polyval(np.polyfit(t, v, 1), t)
        out[s] = float(np.sqrt((res ** 2).sum() / (len(v) - 2)))
    return out


def _boot_p(samples, null):
    """Two-sided bootstrap P that the statistic equals `null` (percentile method),
    floored at 1 / (B + 1)."""
    s = np.asarray(samples, float)
    s = s[np.isfinite(s)]
    if not len(s):
        return np.nan
    tail = min((s <= null).mean(), (s >= null).mean())
    return float(min(1.0, max(2 * tail, 1.0 / (len(s) + 1))))


# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 -- every number of one species x modality
# ══════════════════════════════════════════════════════════════════════════════
def compare_modes(m, age_col, n_boot=N_BOOT, seed=SEED):
    """m = matched_sessions(...). Returns one dict (one row of the summary sheet)."""
    out = dict(n_subjects=m['subject'].nunique(), n_sessions=len(m))
    # 1. integrity at the reference session, and differences elsewhere (% of the mean)
    d_pct = 100 * (m['long'] - m['cross']) / ((m['long'] + m['cross']) / 2)
    out['ref_abs_diff_pct_median'] = float(d_pct[m['is_reference']].abs().median())
    out['nonref_abs_diff_pct_median'] = float(d_pct[~m['is_reference']].abs().median())
    # 2. agreement over all sessions
    out['icc_a1'] = icc_a1(m['long'], m['cross'])
    out['bias_pct'] = float(d_pct.mean())
    out['loa_low'] = float(d_pct.mean() - 1.96 * d_pct.std(ddof=1))
    out['loa_high'] = float(d_pct.mean() + 1.96 * d_pct.std(ddof=1))
    # 3. same LMM in each mode, on the same scans
    fit = {}
    for key in ('long', 'cross'):
        fit[key] = lmm_random_intercept(m[key], m[age_col], m['subject'])
        for k in ('slope', 'se', 'p', 'sigma_e', 'sigma_u'):
            out[f'{key}_{k}'] = fit[key][k]
        # slope relative to the mean first-session value (% per year), as Results
        base = m.sort_values(age_col).groupby('subject')[key].first().mean()
        out[f'{key}_slope_pct_per_year'] = 100 * fit[key]['slope'] / base
        out[f'{key}_se_pct_per_year'] = 100 * fit[key]['se'] / abs(base)
    out['sigma_e_ratio'] = fit['cross']['sigma_e'] / fit['long']['sigma_e']
    out['se_ratio'] = fit['cross']['se'] / fit['long']['se']
    out['n_equivalent_ratio'] = out['se_ratio'] ** 2
    out['delta_slope'] = fit['long']['slope'] - fit['cross']['slope']
    # 4. paired cluster bootstrap over subjects: the SAME resampled subjects in both modes
    rng = np.random.default_rng(seed)
    subs = m['subject'].unique()
    by_sub = {s: g for s, g in m.groupby('subject')}
    bs = []
    for _ in range(n_boot):
        pick = rng.choice(subs, len(subs), replace=True)
        # a subject drawn twice becomes two clusters (new id), as required by the
        # cluster bootstrap (Field & Welsh 2007 JRSS-B 69:369)
        g = pd.concat([by_sub[s].assign(subject=f'{s}#{i}') for i, s in enumerate(pick)])
        try:
            fl = lmm_random_intercept(g['long'], g[age_col], g['subject'])
            fc = lmm_random_intercept(g['cross'], g[age_col], g['subject'])
        except np.linalg.LinAlgError:
            continue                          # degenerate resample (one age only)
        bs.append((fc['sigma_e'] / fl['sigma_e'], fc['se'] / fl['se'],
                   fl['slope'] - fc['slope']))
    bs = np.asarray(bs)
    for j, name, null in ((0, 'sigma_e_ratio', 1.0), (1, 'se_ratio', 1.0),
                          (2, 'delta_slope', 0.0)):
        out[f'{name}_ci_low'], out[f'{name}_ci_high'] = np.percentile(bs[:, j], [2.5, 97.5])
        out[f'{name}_p'] = _boot_p(bs[:, j], null)
    out['n_boot'] = len(bs)
    # 5. model-free within-subject residual SD (subjects with >= 3 sessions)
    rl = within_subject_residual_sd(m, age_col, 'long')
    rc = within_subject_residual_sd(m, age_col, 'cross')
    common = sorted(set(rl) & set(rc))
    out['n_subjects_3plus'] = len(common)
    if len(common) >= 5:
        a, b = np.array([rl[s] for s in common]), np.array([rc[s] for s in common])
        out['resid_sd_long_median'], out['resid_sd_cross_median'] = np.median(a), np.median(b)
        out['resid_sd_wilcoxon_p'] = float(stats.wilcoxon(a, b).pvalue)
        out['n_subjects_less_noisy_long'] = int((a < b).sum())
    return out


def nonlinearity(m, age_col):
    """Spline vs linear LRT P in each mode (fit_gamm of Fig. 3)."""
    from Statistics.Group_fMRI.EDNiX_longitudinal import fit_gamm
    out = {}
    for key in ('long', 'cross'):
        g = fit_gamm(m.rename(columns={key: 'v'}), 'v', age_col, 2) or {}
        out[f'{key}_p_nonlinear'] = g.get('p_nonlinear', np.nan)
        out[f'{key}_p_age_spline'] = g.get('p_age', np.nan)
    return out


def assoc_primary_both(morph, age_df, species):
    """Fig. 3g-i / Supp. Fig. 7 association - primary paired t-test in both modes."""
    from Exemples.Study_EDNiX.analysis.library import longitudinal as lg
    rows = []
    for mod, _, _, title in lg.FIG3_MODALITIES:
        for mode in (MODE_LONG, MODE_CROSS):
            res = lg.assoc_primary_trajectories(morph[mode], age_df, lg.PRIMARY_REGIONS,
                                                lg.ASSOCIATION_REGIONS, species_label=species,
                                                data_key=mod)
            if res is None:
                continue
            t = lg._assoc_primary_test(res[0], res[1])
            rows.append(dict(species=species, modality=mod, mode=mode, n=t['n'],
                             assoc_pct=t['assoc_pct'], primary_pct=t['primary_pct'],
                             t=t['t'], p=t['p'], dz=t['dz']))
    return rows


# ══════════════════════════════════════════════════════════════════════════════
# STEP 6 -- all species and modalities; text of Supp_Fig_15_stats.txt
# ══════════════════════════════════════════════════════════════════════════════
def ablation_stats(rows, n_boot=N_BOOT):
    """rows = load_ablation_data(). Returns dict(summary, assoc, paired, lines, matched)."""
    from Exemples.Study_EDNiX.analysis.library import longitudinal as lg
    summary, assoc, paired, matched = [], [], [], {}
    for sp, morph, age_df, age_col, folders in rows:
        for mod, col, label, title in lg.FIG3_MODALITIES:
            m = matched_sessions(morph, age_df, age_col, mod, col)
            if m.empty or m['subject'].nunique() < 3:
                warnings.warn(f'{sp} {mod}: fewer than 3 matched subjects, skipped')
                continue
            matched[(sp, mod)] = (m, age_col, label)
            r = dict(species=sp, modality=mod, folder_long=folders[MODE_LONG],
                     folder_cross=folders[MODE_CROSS])
            r.update(compare_modes(m, age_col, n_boot=n_boot))
            r.update(nonlinearity(m, age_col))
            summary.append(r)
            paired.append(m.assign(species=sp, modality=mod))
        assoc += assoc_primary_both(morph, age_df, sp)
    summary, assoc = pd.DataFrame(summary), pd.DataFrame(assoc)
    paired = pd.concat(paired, ignore_index=True) if paired else pd.DataFrame()
    return dict(summary=summary, assoc=assoc, paired=paired, matched=matched,
                lines=result_lines(summary, assoc))


def _p(p):
    return 'P < 0.001' if p < 0.001 else f'P = {p:.2g}'


def result_lines(summary, assoc):
    L = ['Supp. Fig. 15 -- longitudinal vs cross-sectional processing of the same scans',
         f'{REGION}, level 1, left hemisphere, subjects with >= 2 matched sessions; LMM '
         'metric ~ age + (1|subject), ML; ratios = cross-sectional / longitudinal; '
         'delta = longitudinal - cross-sectional; paired subject bootstrap 95% CI and P.', '']
    for r in summary.itertuples():
        L.append(
            f'{r.species} {r.modality} (n = {r.n_subjects} subjects, {r.n_sessions} sessions): '
            f'reference session |diff| median {r.ref_abs_diff_pct_median:.2f}%, other sessions '
            f'{r.nonref_abs_diff_pct_median:.2f}%; ICC(A,1) = {r.icc_a1:.3f}, bias '
            f'{r.bias_pct:+.2f}% (LoA {r.loa_low:+.2f} to {r.loa_high:+.2f}%); '
            f'slope/yr long {r.long_slope:.4g} +- {r.long_se:.2g} ({_p(r.long_p)}), cross '
            f'{r.cross_slope:.4g} +- {r.cross_se:.2g} ({_p(r.cross_p)}); '
            f'%/yr long {r.long_slope_pct_per_year:+.2f}, cross {r.cross_slope_pct_per_year:+.2f}; '
            f'sigma_e ratio {r.sigma_e_ratio:.2f} [{r.sigma_e_ratio_ci_low:.2f}, '
            f'{r.sigma_e_ratio_ci_high:.2f}] {_p(r.sigma_e_ratio_p)}; SE ratio {r.se_ratio:.2f} '
            f'[{r.se_ratio_ci_low:.2f}, {r.se_ratio_ci_high:.2f}] {_p(r.se_ratio_p)} '
            f'(n-equivalent x{r.n_equivalent_ratio:.2f}); delta slope {r.delta_slope:+.4g} '
            f'[{r.delta_slope_ci_low:+.4g}, {r.delta_slope_ci_high:+.4g}] {_p(r.delta_slope_p)}; '
            f'non-linear P long {r.long_p_nonlinear:.3g}, cross {r.cross_p_nonlinear:.3g}')
        if np.isfinite(getattr(r, 'resid_sd_wilcoxon_p', np.nan)):
            L.append(f'    model-free: {r.n_subjects_less_noisy_long}/{r.n_subjects_3plus} '
                     f'subjects less noisy in longitudinal mode, residual SD median '
                     f'{r.resid_sd_long_median:.4g} vs {r.resid_sd_cross_median:.4g}, '
                     f'Wilcoxon {_p(r.resid_sd_wilcoxon_p)}')
    if not assoc.empty:
        L += ['', 'Association - primary, paired t on last-session % change (Fig. 3g-i, Supp. Fig. 7):']
        for r in assoc.itertuples():
            L.append(f'  {r.species} {r.modality} {r.mode}: {r.assoc_pct:+.2f}% vs '
                     f'{r.primary_pct:+.2f}%, t({r.n - 1}) = {r.t:.2f}, {_p(r.p)}, dz = {r.dz:.2f}')
    return L


# ══════════════════════════════════════════════════════════════════════════════
# STEP 7 -- drawing
# ══════════════════════════════════════════════════════════════════════════════
def draw_supp15(res, out_base):
    """a-f: per-subject values in both modes (species rows, modality columns; solid =
    longitudinal, dashed grey = cross-sectional, same scans). g: within-subject residual
    SD ratio (cross / long) with bootstrap 95% CI. h: slope in % per year, both modes."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from Exemples.Study_EDNiX.analysis import paper_style as ps
    from Exemples.Study_EDNiX.analysis.library import longitudinal as lg
    summ = res['summary']
    species = list(dict.fromkeys(summ['species']))
    mods = lg.FIG3_MODALITIES
    with plt.rc_context(ps.PAPER_RC):
        fig = plt.figure(figsize=(ps.WIDTH_2COL, (48 * len(species) + 55) * ps.MM))
        gs = fig.add_gridspec(len(species) + 1, 6, hspace=0.75, wspace=1.4)
        letters = iter('abcdefghijklmnop')
        for i, sp in enumerate(species):
            col = ps.species_color(sp)
            for j, (mod, _, label, title) in enumerate(mods):
                ax = fig.add_subplot(gs[i, 2 * j:2 * j + 2])
                ps.panel_label(ax, next(letters), x=-0.3, y=1.02)
                if (sp, mod) not in res['matched']:
                    ax.set_axis_off()
                    continue
                m, age_col, _ = res['matched'][(sp, mod)]
                for _, g in m.sort_values(age_col).groupby('subject'):
                    ax.plot(g[age_col], g['cross'], color='0.55', lw=0.6, ls='--')
                    ax.plot(g[age_col], g['long'], color=col, lw=0.7)
                r = summ[(summ['species'] == sp) & (summ['modality'] == mod)].iloc[0]
                ax.text(0.03, 0.97, f'$\\sigma_e$ ratio {r.sigma_e_ratio:.2f}\n{_p(r.sigma_e_ratio_p)}',
                        transform=ax.transAxes, va='top', fontsize=5)
                ax.set_xlabel('Age (years)')
                ax.set_ylabel(label)
                ax.set_title(f'{sp}, {title.lower()}', fontsize=6)
        # g, h: forest plots, one row per species x modality
        lab = [f'{r.species} {r.modality}' for r in summ.itertuples()]
        y = np.arange(len(summ))[::-1]
        ax = fig.add_subplot(gs[-1, 0:3])
        ps.panel_label(ax, next(letters), x=-0.45, y=1.02)
        c = [ps.species_color(s) for s in summ['species']]
        ax.errorbar(summ['sigma_e_ratio'], y, fmt='none', ecolor=c, elinewidth=0.8,
                    xerr=[summ['sigma_e_ratio'] - summ['sigma_e_ratio_ci_low'],
                          summ['sigma_e_ratio_ci_high'] - summ['sigma_e_ratio']])
        ax.scatter(summ['sigma_e_ratio'], y, c=c, s=10, zorder=3)
        ax.axvline(1, color='0.6', lw=0.6, ls=':')
        ax.set_yticks(y)
        ax.set_yticklabels(lab, fontsize=5)
        ax.set_xlabel('Within-subject residual SD,\ncross-sectional / longitudinal')
        ax = fig.add_subplot(gs[-1, 3:6])
        ps.panel_label(ax, next(letters), x=-0.1, y=1.02)
        for key, dy, mk, fc in (('long', 0.15, 'o', None), ('cross', -0.15, 's', 'white')):
            ci = 1.96 * summ[f'{key}_se_pct_per_year']
            ax.errorbar(summ[f'{key}_slope_pct_per_year'], y + dy, xerr=ci, fmt='none',
                        ecolor=c, elinewidth=0.8)
            ax.scatter(summ[f'{key}_slope_pct_per_year'], y + dy, s=10, marker=mk, zorder=3,
                       facecolors=c if fc is None else fc, edgecolors=c,
                       label=MODE_LONG if key == 'long' else MODE_CROSS)
        ax.axvline(0, color='0.6', lw=0.6, ls=':')
        ax.set_yticks(y)
        ax.set_yticklabels([])
        ax.set_xlabel('Slope (% of baseline per year, 95% CI)')
        ax.legend(fontsize=5, frameon=False, loc='upper right')
        paths = ps.save(fig, out_base)
        plt.close(fig)
    return paths
