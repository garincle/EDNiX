"""
EDNiX — Primary Network FC Threshold Explorer & QC Classifier
==============================================================

Each run_*.csv file is treated as an independent observation.
One df_qc row = one BIDS session. score_all_runs() expands each row into
N run rows (one per run_* file found in the session directory).
Classification, fingerprint, and all figures operate on individual runs.

FIXED: All functions use the SAME run discovery logic to ensure consistency.
"""

import os, re, glob, argparse, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.cm as _mpl_cm
from matplotlib.colors import TwoSlopeNorm, to_rgba
from scipy import stats as _stats
from scipy.signal import find_peaks

# =============================================================================
# CONSTANTS
# =============================================================================

_CAT_COLORS = {"Specific": "#009E73", "Unspecific": "#E69F00",
               "No": "#D55E00", "Spurious": "#CC79A7"}
_PHYLO = ["Mouse", "Rat", "Mouselemur", "Bat", "Marmoset",
          "Macaque", "Human", "Dog", "Cat", "Pig"]
_PALETTE = ["#0072B2", "#E69F00", "#009E73", "#F0E442", "#56B4E9",
            "#D55E00", "#CC79A7", "#999999", "#332288", "#44AA99",
            "#117733", "#882255", "#DDCC77", "#88CCEE", "#AA4499"]
_GOLD = "#FFD700"
_GOLD_EDGE = "#B8860B"
CATS = ["Specific", "Unspecific", "No", "Spurious"]
_DEFAULT_BIDS = "/scratch2/EDNiX/{species}/{bids_dir}"

_PRIMARY_NETWORKS = {
    "motor": ["Motor_and_premotor"],
    "somatosensory": ["Somatosensory_cortex"],
    "visual_striate": ["Visual_striate_cortex"],
    "visual_extra": ["Visual_pre_and_extra_striate_cortex"],
    "auditory": ["Auditory_cortex"],
}

# =============================================================================
# HELPERS
# =============================================================================

def _sp_order(df):
    vals = df["species"].unique().tolist() if "species" in df.columns else []
    return [s for s in _PHYLO if s in vals] + [s for s in vals if s not in _PHYLO]

def _gcolors(groups):
    return {g: _PALETTE[i % len(_PALETTE)] for i, g in enumerate(groups)}

def _save_or_show(fig, path):
    if path:
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        fig.savefig(path, bbox_inches="tight", dpi=200)
        plt.close(fig)
        print(f"  [saved] {path}")
    else:
        plt.show()

def _load_matrix(path):
    try:
        from Plotting.ednix_bids_tools import load_corr_matrix
        return load_corr_matrix(path)
    except ImportError:
        df = pd.read_csv(path, index_col=0)
        return list(df.index), df.values.astype(float)

def _find_all_run_paths(bids_root, subject, session, fit_kind,
                        atlas_name="EDNIxCSC", atlas_level=2, use_lr=True):
    """
    Find ALL run matrix files - tries multiple directory structures.
    Returns all run_*_matrix.csv files found.
    """
    lr = "LR" if use_lr else ""
    excl = ("check_fit", "flattened", "pval", "tstat")
    fname_pat = f"{atlas_name}{lr}_{atlas_level}_run_*{fit_kind}_matrix.csv"

    sub_dir = os.path.join(bids_root, f"sub-{subject}")

    # Build list of candidate session directories to try
    ses_dirs = []
    if session:
        ses_dirs.append(os.path.join(sub_dir, f"ses-{session}"))
    ses_dirs.append(sub_dir)  # no session dir fallback

    all_files = []
    for ses_d in ses_dirs:
        # Try multiple intermediate path patterns
        for mid in [
            os.path.join("func", "acpc-func", "Stats", "Correl_matrix"),
            os.path.join("func", "Stats", "Correl_matrix"),
            os.path.join("func", "acpc-func", "Correl_matrix"),
            os.path.join("func", "Correl_matrix"),
        ]:
            pat = os.path.join(ses_d, mid, fname_pat)
            files = sorted([f for f in glob.glob(pat)
                           if not any(x in os.path.basename(f) for x in excl)])
            if files:
                all_files.extend(files)

    # Remove duplicates while preserving order
    seen = set()
    unique_files = []
    for f in all_files:
        if f not in seen:
            seen.add(f)
            unique_files.append(f)

    return unique_files

def _get_all_run_paths_for_row(fit_kind, row, bids_root_template,
                                atlas_name="EDNIxCSC", atlas_level=2, use_lr=True):
    """
    Return ALL run_* matrix paths for one df_qc row.
    Uses bids_root_template to find all runs.
    """
    sp = str(row.get("species", ""))
    bd = str(row.get("bids_dir", ""))
    sub = str(row.get("subject", ""))
    ses = str(row.get("session", ""))

    if not bids_root_template:
        return []

    root = bids_root_template.format(species=sp, bids_dir=bd)
    paths = _find_all_run_paths(root, sub, ses, fit_kind, atlas_name, atlas_level, use_lr)

    return [(p, sp, bd, sub, ses) for p in paths]

def _find_primary_rois(rois):
    found = {}
    for net, patterns in _PRIMARY_NETWORKS.items():
        l_idx = r_idx = None
        for p in patterns:
            for i, r in enumerate(rois):
                if p.lower() in r.lower():
                    if r.startswith(("L_", "l_")) and l_idx is None:
                        l_idx = i
                    elif r.startswith(("R_", "r_")) and r_idx is None:
                        r_idx = i
            if l_idx is not None and r_idx is not None:
                break
        if l_idx is not None and r_idx is not None:
            found[net] = {"L": l_idx, "R": r_idx,
                          "L_name": rois[l_idx], "R_name": rois[r_idx]}
    return found

def _is_homotopic(a, b):
    def _base(n):
        return re.sub(r"^[LlRr]_", "", str(n)).lower().strip()
    return _base(a) == _base(b)

def _kde_valley(values, n_pts=1000):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if len(v) < 5:
        return float(np.nanmedian(v))
    kde = _stats.gaussian_kde(v, bw_method="silverman")
    x = np.linspace(v.min(), v.max(), n_pts)
    y = kde(x)
    peaks, _ = find_peaks(y, prominence=y.max() * 0.05)
    if len(peaks) < 2:
        return float(np.nanpercentile(v, 30))
    top2 = peaks[np.argsort(y[peaks])[-2:]]
    lo, hi = sorted(top2)
    valley_idx = lo + int(np.argmin(y[lo:hi + 1]))
    return float(x[valley_idx])

# =============================================================================
# PRIMARY NETWORK SCORE
# =============================================================================
def compute_primary_score_from_matrix(mat, rois):
    primary = _find_primary_rois(rois)
    nets = list(primary.keys())
    if not nets:
        return dict(intra_means={}, mean_intra=np.nan,
                    mean_inter=np.nan, mean_inter_hetero=np.nan,
                    contrast_ratio=np.nan, networks_found=[])

    # Specific signal: homotopic same-network (L_Motor ↔ R_Motor)
    intra_means = {net: float(mat[v["L"], v["R"]])
                   if np.isfinite(mat[v["L"], v["R"]]) else np.nan
                   for net, v in primary.items()}
    mean_intra = float(np.nanmean(list(intra_means.values())))

    inter_vals_all    = []  # all cross-network pairs (kept for reference)
    inter_vals_hetero = []  # contralateral cross-network only: L_netA ↔ R_netB

    for i in range(len(nets)):
        for j in range(i + 1, len(nets)):
            na, nb = nets[i], nets[j]
            ia_L, ia_R = primary[na]["L"], primary[na]["R"]
            ib_L, ib_R = primary[nb]["L"], primary[nb]["R"]

            # All four cross-network pairs
            for ia, ib in [(ia_L, ib_L), (ia_L, ib_R), (ia_R, ib_L), (ia_R, ib_R)]:
                r = float(mat[ia, ib])
                if np.isfinite(r):
                    inter_vals_all.append(r)

            # Contralateral cross-network only (heterotopic):
            # L_netA ↔ R_netB  and  R_netA ↔ L_netB
            # excludes ipsilateral: L_netA ↔ L_netB  and  R_netA ↔ R_netB
            for ia, ib in [(ia_L, ib_R), (ia_R, ib_L)]:
                r = float(mat[ia, ib])
                if np.isfinite(r):
                    inter_vals_hetero.append(r)

    mean_inter        = float(np.nanmean(inter_vals_all))    if inter_vals_all    else np.nan
    mean_inter_hetero = float(np.nanmean(inter_vals_hetero)) if inter_vals_hetero else np.nan

    # contrast_ratio uses heterotopic inter as denominator
    denom = mean_inter_hetero
    if np.isnan(mean_intra) or np.isnan(denom):
        ratio = np.nan
    elif denom <= 0:
        ratio = 20.0
    else:
        ratio = mean_intra / denom

    return dict(intra_means=intra_means,
                mean_intra=mean_intra,
                mean_inter=mean_inter,
                mean_inter_hetero=mean_inter_hetero,
                contrast_ratio=ratio,
                networks_found=nets)

# =============================================================================
# STEP 1 — SCORE ALL RUNS
# =============================================================================

def score_all_runs(df_qc, fit_kind,
                   bids_root_template=_DEFAULT_BIDS,
                   atlas_name="EDNIxCSC", atlas_level=2,
                   use_lr=True, verbose=True):
    """
    Expand df_qc: each BIDS session row becomes N run rows (one per run file).
    Each run is scored independently — no averaging at any level.
    """
    all_runs = []
    n_scored = 0

    for idx, row in df_qc.iterrows():
        run_infos = _get_all_run_paths_for_row(
            fit_kind, row, bids_root_template, atlas_name, atlas_level, use_lr)

        for run_path, sp, bd, sub, ses in run_infos:
            try:
                rois, mat = _load_matrix(run_path)
                ps = compute_primary_score_from_matrix(mat, rois)

                run_row = {
                    "species": sp,
                    "bids_dir": bd,
                    "subject": sub,
                    "session": ses,
                    "run": os.path.basename(run_path),
                    "corr_matrix_path": run_path,
                    "primary_intra": ps["mean_intra"],
                    "primary_inter": ps["mean_inter_hetero"],
                    "contrast_ratio": ps["contrast_ratio"],
                }
                for net, val in ps["intra_means"].items():
                    run_row[f"per_net_intra_{net}"] = val

                # Carry over extra columns from original df_qc
                for col in df_qc.columns:
                    if col not in run_row:
                        run_row[col] = row.get(col, np.nan)

                all_runs.append(run_row)
                n_scored += 1

            except Exception as e:
                warnings.warn(f"[score] {run_path}: {e}")

    if not all_runs:
        raise ValueError("No runs could be loaded! Check bids_root_template.")

    df_scored = pd.DataFrame(all_runs)

    if verbose:
        print(f"  [score] {n_scored} runs scored from {len(df_qc)} session rows")
        pi = df_scored["primary_intra"].dropna()
        cr = df_scored["contrast_ratio"].dropna()
        if len(pi):
            print(f"  [score] primary_intra: {pi.min():.3f}--{pi.max():.3f} mean={pi.mean():.3f}")
        if len(cr):
            print(f"  [score] contrast_ratio: {cr.min():.2f}--{cr.max():.2f} mean={cr.mean():.2f}")

    return df_scored

# =============================================================================
# STEP 2 — THRESHOLDS FROM 2D CLASSIFICATION SPACE (GMM)
# =============================================================================
def suggest_thresholds(df_scored,
                       method="auto",
                       thresh_intra=None,
                       thresh_inter=None,
                       thresh_delta=None,
                       stringency=0.0,
                       verbose=True):
    """
    Thresholds for delta-based classifier:
        delta = primary_intra - primary_inter_hetero
        Specific   : delta > thresh_delta  AND  intra > thresh_intra
        No         : intra <= thresh_intra
        Unspecific : intra > thresh_intra  AND  delta <= thresh_delta
        Spurious   : inter_hetero > intra  (delta < 0, caught before the rest)

    thresh_delta defaults to the 25th-percentile of observed deltas.
    """
    inter_col = ("primary_inter_hetero"
                 if "primary_inter_hetero" in df_scored.columns
                 else "primary_inter")

    valid = (df_scored["primary_intra"].notna() &
             df_scored[inter_col].notna())
    intra_vals = df_scored.loc[valid, "primary_intra"].values
    inter_vals = df_scored.loc[valid, inter_col].values
    delta_vals  = intra_vals - inter_vals

    if len(intra_vals) < 10:
        t_intra = 0.2
        t_inter = 0.2   # kept for back-compat / display
        t_delta = 0.0
    else:
        s = np.clip(stringency, -1.0, 1.0)

        base_intra  = np.percentile(intra_vals, 25)
        strict_intra  = np.percentile(intra_vals, 50)
        lenient_intra = np.percentile(intra_vals, 10)

        base_inter  = np.percentile(inter_vals,  75)
        strict_inter  = np.percentile(inter_vals,  50)
        lenient_inter = np.percentile(inter_vals,  90)

        base_delta  = np.percentile(delta_vals,  25)
        strict_delta  = np.percentile(delta_vals,  50)
        lenient_delta = np.percentile(delta_vals,  10)

        if s >= 0:
            auto_intra = base_intra + s * (strict_intra  - base_intra)
            auto_inter = base_inter + s * (strict_inter  - base_inter)
            auto_delta = base_delta + s * (strict_delta  - base_delta)
        else:
            auto_intra = base_intra + abs(s) * (lenient_intra - base_intra)
            auto_inter = base_inter + abs(s) * (lenient_inter - base_inter)
            auto_delta = base_delta + abs(s) * (lenient_delta - base_delta)

        auto_intra = float(np.clip(auto_intra, 0.05, 0.8))
        auto_inter = float(np.clip(auto_inter, 0.05, 0.8))
        # delta can legitimately be negative → don't clip below 0
        auto_delta = float(auto_delta)

        t_intra = thresh_intra if thresh_intra is not None else auto_intra
        t_inter = thresh_inter if thresh_inter is not None else auto_inter
        t_delta = thresh_delta if thresh_delta is not None else auto_delta

    # Preview
    n_specific = (
        (intra_vals > t_intra) &
        ((intra_vals - inter_vals) > t_delta)
    ).sum()
    pct_specific = 100.0 * n_specific / len(intra_vals) if len(intra_vals) else 0

    if verbose:
        print(f"\n{'=' * 50}")
        print(f"  THRESHOLD SELECTION  (stringency={stringency:+.1f})")
        print(f"  intra  > {t_intra:.3f}")
        print(f"  delta  > {t_delta:.3f}  (intra − inter_hetero)")
        print(f"  [ref]  inter_hetero < {t_inter:.3f}  (informational)")
        print(f"\n  → {n_specific}/{len(intra_vals)} runs Specific ({pct_specific:.0f}%)")
        print(f"{'=' * 50}\n")

    return dict(thresh_intra=t_intra, thresh_inter=t_inter, thresh_delta=t_delta,
                auto_intra=auto_intra, auto_inter=auto_inter, auto_delta=auto_delta,
                stringency=stringency)

# =============================================================================
# STEP 3 — CLASSIFY
# =============================================================================

def classify_runs(df_scored, thresh_intra, thresh_inter=None,
                  thresh_delta=None, verbose=True):
    """
    Delta-based classification.

        delta = primary_intra − primary_inter_hetero

        inter_hetero > intra  →  Spurious   (negative delta, leakage / noise)
        intra <= thresh_intra →  No         (too weak bilateral signal)
        delta <= thresh_delta →  Unspecific (bilateral but no network segregation)
        else                  →  Specific
    """
    inter_col = ("primary_inter_hetero"
                 if "primary_inter_hetero" in df_scored.columns
                 else "primary_inter")

    if thresh_delta is None:
        thresh_delta = 0.0   # sensible fallback: delta must be positive

    cats = pd.Series("No", index=df_scored.index, dtype=object)
    intra_arr = df_scored["primary_intra"].values.astype(float)
    inter_arr = df_scored[inter_col].values.astype(float)

    for i, idx in enumerate(df_scored.index):
        intra = intra_arr[i]
        inter = inter_arr[i]
        if not np.isfinite(intra) or not np.isfinite(inter):
            continue
        delta = intra - inter
        if inter > intra:                  # delta < 0
            cats[idx] = "Spurious"
        elif intra <= thresh_intra:
            cats[idx] = "No"
        elif delta <= thresh_delta:
            cats[idx] = "Unspecific"
        elif delta > thresh_delta and intra > thresh_intra:
            cats[idx] = "Specific"
        else:
            print(
                f"Unhandled classification state: "
                f"intra={intra:.3f}, inter_hetero={inter:.3f}, "
                f"delta={delta:.3f}"
            )
            cats[idx] = 'Unclassified'

    if verbose:
        print(f"  Classification  thresh_intra={thresh_intra:.4f}  "
              f"thresh_delta={thresh_delta:.4f}  (col={inter_col})")
        if "species" in df_scored.columns:
            for sp in _sp_order(df_scored):
                mask = df_scored["species"] == sp
                if not mask.any():
                    continue
                c = cats[mask].value_counts()
                n = mask.sum()
                print(f"    {sp:14s} n={n:3d}   " +
                      "  ".join(f"{cat}={c.get(cat, 0)}" for cat in CATS))

    return cats

# =============================================================================
# FIGURE: THRESHOLD DISTRIBUTIONS
# =============================================================================

def plot_threshold_distributions(df_scored, thresh_intra=None, thresh_inter=None,
                                 output_path=None):
    if thresh_intra is None or thresh_inter is None:
        t = suggest_thresholds(df_scored, method="auto", verbose=False)
        thresh_intra = thresh_intra or t["thresh_intra"]
        thresh_inter = thresh_inter or t["thresh_inter"]

    sp_list = _sp_order(df_scored)
    gc = _gcolors(sp_list)
    markers = ["o", "s", "^", "D", "v", "P", "*", "h"]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(f"Threshold choice —\n"
                 f"thresh_intra={thresh_intra:.3f} thresh_inter={thresh_inter:.3f}",
                 fontsize=12, fontweight="bold")

    for ax, col, thresh, xlabel, want_above in [
        (axes[0], "primary_intra", thresh_intra, "Primary intra-network r (want > thresh)", True),
        (axes[1], "primary_inter", thresh_inter, "Primary inter-network r (want < thresh)", False),
    ]:
        for sp in sp_list:
            vals = df_scored.loc[df_scored["species"] == sp, col].dropna().values \
                   if "species" in df_scored.columns else df_scored[col].dropna().values
            vals = vals[np.isfinite(vals)]
            if len(vals) < 3:
                continue
            x = np.linspace(vals.min() - 0.02, vals.max() + 0.02, 400)
            try:
                k = _stats.gaussian_kde(vals, bw_method="silverman")
                ax.fill_between(x, k(x), alpha=0.12, color=gc[sp])
                ax.plot(x, k(x), color=gc[sp], lw=1.5, label=sp)
            except Exception:
                pass
        ls = "--" if want_above else ":"
        ax.axvline(thresh, color="#222", lw=2.0, ls=ls, label=f"thresh={thresh:.3f}")
        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel("Density", fontsize=9)
        title = "Bilateral connectivity" if want_above else "Network segregation"
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.legend(fontsize=7, frameon=False)
        ax.set_ylim(bottom=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    ax = axes[2]
    for i, sp in enumerate(sp_list):
        mask = df_scored["species"] == sp if "species" in df_scored.columns else pd.Series(True, index=df_scored.index)
        x = df_scored.loc[mask, "primary_intra"].values.astype(float)
        y = df_scored.loc[mask, "primary_inter"].values.astype(float)
        fin = np.isfinite(x) & np.isfinite(y)
        ax.scatter(x[fin], y[fin], color=gc[sp], marker=markers[i % len(markers)],
                   s=12, alpha=0.5, edgecolors="none", label=sp, zorder=3)

    xmin = float(np.nanmin(df_scored["primary_intra"].dropna())) - 0.01
    xmax = float(np.nanmax(df_scored["primary_intra"].dropna())) + 0.01
    ymin = float(np.nanmin(df_scored["primary_inter"].dropna())) - 0.01
    ymax = float(np.nanmax(df_scored["primary_inter"].dropna())) + 0.01
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axvline(thresh_intra, color="#222", lw=1.8, ls="--", zorder=5)
    ax.axhline(thresh_inter, color="#444", lw=1.4, ls=":", zorder=5)
    diag = np.linspace(min(xmin, ymin), max(xmax, ymax), 100)
    ax.plot(diag, diag, color="#bbb", lw=1.0, ls="-", zorder=2, label="inter=intra")

    ax.fill_between([thresh_intra, xmax], ymin, thresh_inter,
                    color=to_rgba(_CAT_COLORS["Specific"], 0.09))
    ax.fill_between([thresh_intra, xmax], thresh_inter, ymax,
                    color=to_rgba(_CAT_COLORS["Unspecific"], 0.09))
    ax.fill_betweenx([ymin, ymax], xmin, thresh_intra,
                     color=to_rgba(_CAT_COLORS["No"], 0.07))

    for lbl, x, y in [
        ("Specific", (thresh_intra + xmax) / 2, (ymin + thresh_inter) / 2),
        ("Unspecific", (thresh_intra + xmax) / 2, (thresh_inter + ymax) / 2),
        ("No", (xmin + thresh_intra) / 2, ymax * 0.8),
    ]:
        ax.text(x, y, lbl, fontsize=9, color=_CAT_COLORS[lbl],
                fontweight="bold", alpha=0.45, ha="center")

    ax.set_xlabel("Primary intra-network r", fontsize=10)
    ax.set_ylabel("Primary inter-network r", fontsize=10)
    ax.set_title("2D classification space", fontsize=10, fontweight="bold")
    ax.legend(fontsize=7, frameon=False, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    _save_or_show(fig, output_path)
    return fig

# =============================================================================
# FINGERPRINT BUILDING
# =============================================================================

def build_fingerprint_from_specific(df_scored, cats,
                                    bids_root_template=_DEFAULT_BIDS,
                                    atlas_name="EDNIxCSC", atlas_level=2,
                                    use_lr=True, n_top=10, verbose=True):
    """Build fingerprint from all Specific runs (each run = independent observation)."""
    df = df_scored.copy()
    df["_cat"] = cats.values if hasattr(cats, "values") else cats

    sp_mats = {}
    for sp in _sp_order(df):
        if "species" not in df.columns:
            break
        sp_df = df[(df["species"] == sp) & (df["_cat"] == "Specific")]
        if sp_df.empty:
            if verbose:
                print(f"  [fingerprint] {sp}: 0 Specific runs -- skipped")
            continue

        mats = []
        rois_ref = None
        for _, row in sp_df.iterrows():
            path = row.get("corr_matrix_path", None)
            if not path or not os.path.exists(str(path)):
                continue
            try:
                rois, mat = _load_matrix(path)
                if rois_ref is None:
                    rois_ref = rois
                elif len(rois) != len(rois_ref):
                    continue
                mats.append(mat)
            except Exception as e:
                warnings.warn(f"[fingerprint] {sp}: {e}")

        if mats:
            sp_mats[sp] = {"individual_mats": mats, "rois": rois_ref, "n": len(mats)}
            if verbose:
                print(f"  [fingerprint] {sp}: {len(mats)} Specific runs")

    if len(sp_mats) < 2:
        raise ValueError(f"Only {len(sp_mats)} species with Specific runs.")

    return sp_mats, build_fc_fingerprint(sp_mats, n_top=n_top, verbose=verbose)

def build_fc_fingerprint(sp_mats, n_top=10, verbose=True):
    """Build fingerprint treating all individual runs as independent observations."""
    if len(sp_mats) < 2:
        raise ValueError(f"Need >=2 species, got {len(sp_mats)}")

    first = next(iter(sp_mats.values()))
    common = set(first["rois"])
    for v in sp_mats.values():
        common &= set(v["rois"])

    rois = [r for r in first["rois"] if r in common]
    n_roi = len(rois)
    if n_roi < 4:
        raise ValueError(f"Only {n_roi} common ROIs.")

    sp_names = list(sp_mats.keys())

    # Align all individual run matrices
    all_runs = []
    for sp in sp_names:
        idx_map = [sp_mats[sp]["rois"].index(r) for r in rois]
        for mat in sp_mats[sp]["individual_mats"]:
            all_runs.append(mat[np.ix_(idx_map, idx_map)])

    t_stat = np.full((n_roi, n_roi), np.nan)
    p_val = np.full((n_roi, n_roi), np.nan)
    mean_r = np.full((n_roi, n_roi), np.nan)
    std_r = np.full((n_roi, n_roi), np.nan)

    for i in range(n_roi):
        for j in range(i + 1, n_roi):
            v = np.array([m[i, j] for m in all_runs if np.isfinite(m[i, j])])
            if len(v) >= 2:
                mean_r[i, j] = mean_r[j, i] = float(np.mean(v))
                std_r[i, j] = std_r[j, i] = float(np.std(v))
                t, p = _stats.ttest_1samp(v, 0)
                t_stat[i, j] = t_stat[j, i] = float(t)
                p_val[i, j] = p_val[j, i] = float(p)
            elif len(v) == 1:
                mean_r[i, j] = mean_r[j, i] = float(v[0])
                t_stat[i, j] = t_stat[j, i] = float(v[0] / 1e-9)
                p_val[i, j] = p_val[j, i] = 1.0

    homo_mask = np.zeros((n_roi, n_roi), bool)
    cross_mask = np.zeros((n_roi, n_roi), bool)

    for i, ri in enumerate(rois):
        for j, rj in enumerate(rois):
            if i >= j:
                continue
            l = ri.startswith(("L_", "l_")) and rj.startswith(("R_", "r_"))
            r = ri.startswith(("R_", "r_")) and rj.startswith(("L_", "l_"))
            if l or r:
                (homo_mask if _is_homotopic(ri, rj) else cross_mask)[i, j] = True
                (homo_mask if _is_homotopic(ri, rj) else cross_mask)[j, i] = True

    triu = np.triu(np.ones((n_roi, n_roi), bool), k=1)
    fp = []

    for mtype, label in [(homo_mask, "homotopic"), (cross_mask, "heterotopic"),
                         (triu & ~homo_mask & ~cross_mask, "ipsilateral")]:
        cand = triu & mtype
        idxs = np.argwhere(cand)
        if not len(idxs):
            continue
        abs_t = np.array([abs(t_stat[i, j]) if np.isfinite(t_stat[i, j]) else 0 for i, j in idxs])
        for k in np.argsort(abs_t)[::-1][:n_top]:
            i, j = idxs[k]
            fp.append((rois[i], rois[j], float(mean_r[i, j]),
                       float(t_stat[i, j]), float(p_val[i, j]), label))

    fp.sort(key=lambda x: abs(x[3]), reverse=True)
    fp_idx = [(rois.index(ri), rois.index(rj)) for ri, rj, *_ in fp]
    fp_pattern = np.array([mean_r[i, j] for i, j in fp_idx])

    if verbose:
        n_runs = len(all_runs)
        print(f"  [fingerprint] {len(sp_names)} species x {n_roi} ROIs -> "
              f"{len(fp)} connections ({n_runs} individual runs)")
        print(f"  {'Connection':48s}  {'t':>6}  {'p':>7}  type")
        print(f"  {'-' * 72}")
        for ri, rj, mr, t, p, ct in fp[:8]:
            sig = "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "ns"
            print(f"  {ri[:22]:22s} -- {rj[:22]:22s}  {t:>+6.2f}  {p:>7.4f}{sig:3s}  {ct}")

    return dict(
        rois=rois, n_species=sp_names, sp_mean_mats=None,
        stability=t_stat, mean_r_cross=mean_r, std_r_cross=std_r, pval=p_val,
        fingerprint_connections=fp, fingerprint_indices=fp_idx,
        fingerprint_pattern=fp_pattern, homo_mask=homo_mask,
        cross_mask=cross_mask, n_roi=n_roi,
        homo_in_fp=[(i, j) for (ri, rj, _, _, _, ct), (i, j) in zip(fp, fp_idx) if ct == "homotopic"],
        cross_in_fp=[(i, j) for (ri, rj, _, _, _, ct), (i, j) in zip(fp, fp_idx) if ct == "heterotopic"],
    )

# =============================================================================
# FIGURE B — FINGERPRINT
# =============================================================================
def plot_fc_fingerprint(fp, n_top_per_type=8, output_path=None):
    """
    Figure B — four panels.

    Colorbars placed under the heatmap matrix.
    """
    rois = fp["rois"]
    t_stat = fp["stability"]
    mean_r = fp["mean_r_cross"]
    fp_conn = fp["fingerprint_connections"]
    fp_idx = fp["fingerprint_indices"]
    n_roi = fp["n_roi"]
    short = [r.replace("L_", "L ").replace("R_", "R ")[:13] for r in rois]

    type_colors = {"homotopic": "#009E73", "heterotopic": "#56B4E9", "ipsilateral": "#E69F00"}
    by_type = {"homotopic": [], "heterotopic": [], "ipsilateral": []}

    for conn in fp_conn:
        ct = conn[5]
        if ct in by_type:
            by_type[ct].append(conn)

    for ct in by_type:
        by_type[ct] = sorted(by_type[ct], key=lambda c: c[2], reverse=True)[:n_top_per_type]

    type_labels = {"homotopic": "Homotopic (L↔R same region)",
                   "heterotopic": "Heterotopic (L↔R cross-region)",
                   "ipsilateral": "Ipsilateral (same hemisphere)"}

    # Figure size
    fig_height = max(9, n_top_per_type * 0.55 + 4)
    fig = plt.figure(figsize=(16, fig_height))

    # GridSpec: 4 rows x 2 columns
    # Row 0-2: heatmap (rows 0-2) and bar panels
    # Row 3: colorbars
    gs = gridspec.GridSpec(4, 2,
                           height_ratios=[1, 1, 1, 0.15],
                           width_ratios=[1.0, 0.9],
                           hspace=0.35, wspace=0.25,
                           left=0.05, right=0.96, top=0.92, bottom=0.10)

    # ========== LEFT SIDE: Heatmap (rows 0-2, column 0) ==========
    ax_heat = fig.add_subplot(gs[0:3, 0])
    t_disp = np.where(np.isfinite(t_stat), t_stat, 0)
    vabs = max(float(np.nanpercentile(np.abs(t_disp[t_disp != 0]), 97)), 1.0)
    norm_t = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
    im = ax_heat.imshow(t_disp, cmap="RdBu_r", norm=norm_t,
                        aspect="equal", interpolation="nearest")
    ax_heat.set_xticks(range(n_roi))
    ax_heat.set_yticks(range(n_roi))
    ax_heat.set_xticklabels(short, rotation=90, fontsize=3)
    ax_heat.set_yticklabels(short, fontsize=3)
    ax_heat.set_title("t-statistic per connection\n(vs 0 across all Specific runs)",
                      fontsize=9, fontweight="bold")

    # Gold stars for fingerprint connections
    for i, j in fp_idx:
        ax_heat.plot(j, i, "*", color=_GOLD, ms=4,
                     markeredgecolor=_GOLD_EDGE, markeredgewidth=0.5, zorder=5)

    # ========== COLORBAR FOR HEATMAP (row 3, column 0) ==========
    cbar_ax_heat = fig.add_subplot(gs[3, 0])
    cbar_heat = fig.colorbar(im, cax=cbar_ax_heat, orientation="horizontal", label="t-stat")
    cbar_heat.ax.tick_params(labelsize=7)

    # ========== RIGHT SIDE: Connection type panels ==========
    all_r = [c[2] for ct in by_type for c in by_type[ct] if by_type[ct]]
    r_abs = max(float(np.nanpercentile(np.abs(all_r), 97)), 0.05) if all_r else 0.5
    norm_r = TwoSlopeNorm(vmin=-r_abs, vcenter=0, vmax=r_abs)
    cmap_r = plt.cm.RdBu_r

    for panel_i, ct in enumerate(["homotopic", "heterotopic", "ipsilateral"]):
        ax = fig.add_subplot(gs[panel_i, 1])
        conns = by_type[ct]
        color = type_colors[ct]

        if not conns:
            ax.text(0.5, 0.5, f"No {ct} connections",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=9, color="#aaa")
            ax.set_title(type_labels[ct], fontsize=9, fontweight="bold", color=color)
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        n_show = len(conns)
        y_pos = np.arange(n_show)
        r_vals = [c[2] for c in conns]
        t_vals = [c[3] for c in conns]
        p_vals = [c[4] for c in conns]
        bar_colors = [cmap_r(norm_r(r)) for r in r_vals]

        bars = ax.barh(y_pos, r_vals, color=bar_colors, edgecolor="white", lw=0.5, height=0.7)

        for yi, (r, t, p) in enumerate(zip(r_vals, t_vals, p_vals)):
            sig = "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else ""
            x_txt = r + (0.02 if r >= 0 else -0.02)
            ha_txt = "left" if r >= 0 else "right"
            ax.text(x_txt, yi, f"t={t:+.1f}{sig}", va="center", ha=ha_txt, fontsize=6, color="#333")

        labels = [f"{c[0][:15]} ↔ {c[1][:15]}" for c in conns]
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=6)
        ax.axvline(0, color="#aaa", lw=0.7)
        ax.set_xlabel("Mean r", fontsize=7)
        ax.set_title(type_labels[ct], fontsize=8, fontweight="bold", color=color)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=6)

    # ========== COLORBAR FOR MEAN R (row 3, column 1) ==========
    cbar_ax_r = fig.add_subplot(gs[3, 1])
    sm = plt.cm.ScalarMappable(cmap=cmap_r, norm=norm_r)
    sm.set_array([])
    cbar_r = fig.colorbar(sm, cax=cbar_ax_r, orientation="horizontal", label="Mean r")
    cbar_r.ax.tick_params(labelsize=7)
    _save_or_show(fig, output_path)
    return fig

# =============================================================================
# FIGURE C — CLASSIFICATION MODEL
# =============================================================================

def plot_classification_model(df_scored, thresh_intra, thresh_inter, cats, output_path=None):
    sp_list = _sp_order(df_scored)
    n_sp = len(sp_list)
    if n_sp == 0:
        return None

    gc = _gcolors(sp_list)
    rng = np.random.default_rng(42)

    fig = plt.figure(figsize=(max(n_sp * 3.0, 9), 10))
    gs = gridspec.GridSpec(3, n_sp, wspace=0.06, hspace=0.40,
                           top=0.92, bottom=0.07, left=0.09, right=0.98)

    for col, sp in enumerate(sp_list):
        sp_df = df_scored[df_scored["species"] == sp] if "species" in df_scored.columns else df_scored
        sp_cats = cats[sp_df.index]
        color = gc[sp]
        c_ct = sp_cats.value_counts()

        for row, (metric, thresh, ylabel) in enumerate([
            ("primary_intra", thresh_intra, "Primary intra-network r"),
            ("primary_inter", thresh_inter, "Primary inter-network r"),
        ]):
            ax = fig.add_subplot(gs[row, col])
            av = sp_df[metric].dropna() if metric in sp_df else pd.Series(dtype=float)
            if av.empty:
                if row == 0:
                    ax.set_title(f"{sp}\n(n={len(sp_df)} runs)", fontsize=8, fontweight="bold", color=color)
                continue

            if len(av) >= 4:
                pts = ax.violinplot([av.values], [0], widths=0.65,
                                    showmedians=False, showextrema=False)
                for pc in pts["bodies"]:
                    pc.set_facecolor(color)
                    pc.set_alpha(0.15)
                    pc.set_edgecolor(color)

            jit = rng.uniform(-0.13, 0.13, len(av))
            pt_c = [_CAT_COLORS.get(sp_cats[i], "#aaa") for i in av.index]
            ax.scatter(jit, av.values, c=pt_c, s=18, alpha=0.75, zorder=3, edgecolors="none")

            if len(av) >= 2:
                q1, med, q3 = np.nanpercentile(av.values, [25, 50, 75])
                ax.plot([-0.15, 0.15], [med, med], color=color, lw=2.0, zorder=5)
                ax.plot([0, 0], [q1, q3], color=color, lw=1.2, zorder=4)

            ax.axhline(thresh, color="#333", lw=1.5, ls="--", zorder=6)
            ax.axhline(0, color="#ddd", lw=0.4)
            ax.set_xticks([])
            ax.tick_params(labelsize=7)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["bottom"].set_visible(False)

            if col == 0:
                ax.set_ylabel(ylabel, fontsize=8)
            if row == 0:
                counts = " ".join(f"{cat[0]}={c_ct.get(cat, 0)}" for cat in CATS)
                ax.set_title(f"{sp} (n={len(sp_df)} runs)\n{counts}",
                             fontsize=8, fontweight="bold", color=color)

        ax = fig.add_subplot(gs[2, col])
        xv = sp_df["primary_intra"].dropna() if "primary_intra" in sp_df else pd.Series(dtype=float)
        yv = sp_df["primary_inter"].dropna() if "primary_inter" in sp_df else pd.Series(dtype=float)
        ci = xv.index.intersection(yv.index)

        if not ci.empty:
            pt_c2 = [_CAT_COLORS.get(sp_cats[i], "#aaa") for i in ci]
            ax.scatter(xv[ci], yv[ci], c=pt_c2, s=15, alpha=0.65, zorder=3, edgecolors="none")

        ax.axvline(thresh_intra, color="#333", lw=1.1, ls="--", zorder=5)
        ax.axhline(thresh_inter, color="#555", lw=0.9, ls=":", zorder=5)

        if not ci.empty:
            lim = np.nanmax(np.abs([xv[ci].values, yv[ci].values])) * 1.1
            diag = np.linspace(-lim, lim, 50)
            ax.plot(diag, diag, color="#ccc", lw=0.8, zorder=2)

        ax.set_xlabel("Intra r", fontsize=7)
        if col == 0:
            ax.set_ylabel("Inter r", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.legend(handles=[mpatches.Patch(color=_CAT_COLORS[c], label=c) for c in CATS],
               loc="lower center", bbox_to_anchor=(0.5, 0.0),
               ncol=4, fontsize=8, frameon=False)
    plt.tight_layout(rect=[0, 0.03, 1, 0.998])
    _save_or_show(fig, output_path)
    return fig

# =============================================================================
# FIGURE D — FINGERPRINT HEATMAP
# =============================================================================

def plot_fingerprint_heatmap(df_scored, fp, cats,
                             bids_root_template=_DEFAULT_BIDS,
                             atlas_name="EDNIxCSC", atlas_level=2,
                             use_lr=True, output_path=None):
    fp_conn = fp["fingerprint_connections"]
    fp_idx = fp["fingerprint_indices"]
    fp_rois = fp["rois"]
    fp_pat = fp["fingerprint_pattern"]
    n_fp = len(fp_conn)

    if n_fp == 0:
        return None

    df = df_scored.copy()
    df["_cat"] = cats.values if hasattr(cats, "values") else cats

    if "species" not in df.columns or "bids_dir" not in df.columns:
        warnings.warn("need species and bids_dir")
        return None

    sp_list = _sp_order(df)
    groups = [(sp, bd) for sp in sp_list
              for bd in sorted(df.loc[df["species"] == sp, "bids_dir"].dropna().unique())]

    if not groups:
        return None

    n_grp = len(groups)
    sp_gc = _gcolors(sp_list)
    conn_labels = [f"{c[0][:11]}--{c[1][:11]}" for c in fp_conn]
    conn_types = [c[5] for c in fp_conn]

    group_vals = np.full((n_fp, n_grp), np.nan)
    group_frac = np.zeros(n_grp)

    print("  [fig_D] computing per-(species,BIDS) fingerprint values...")

    for gi, (sp, bd) in enumerate(groups):
        mask = (df["species"] == sp) & (df["bids_dir"] == bd)
        sp_bd = df[mask]
        n_spec = (sp_bd["_cat"] == "Specific").sum()
        group_frac[gi] = n_spec / max(len(sp_bd), 1)
        spec_df = sp_bd[sp_bd["_cat"] == "Specific"]

        conn_acc = [[] for _ in range(n_fp)]
        for _, row in spec_df.iterrows():
            path = row.get("corr_matrix_path", None)
            if not path or not os.path.exists(str(path)):
                continue
            try:
                rois2, mat = _load_matrix(path)
                for k, (fi, fj) in enumerate(fp_idx):
                    ri, rj = fp_rois[fi], fp_rois[fj]
                    if ri in rois2 and rj in rois2:
                        r = float(mat[rois2.index(ri), rois2.index(rj)])
                        if np.isfinite(r):
                            conn_acc[k].append(r)
            except Exception:
                pass

        for k in range(n_fp):
            if conn_acc[k]:
                group_vals[k, gi] = float(np.mean(conn_acc[k]))

    cell_w = max(0.55, 8.0 / n_grp)
    cell_h = max(0.28, 9.0 / n_fp)
    fig_w = float(np.clip(1.8 + n_grp * cell_w + 0.6, 8, 28))
    fig_h = float(np.clip(n_fp * cell_h + 2.0, 6, 20))

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(2, n_grp + 2, height_ratios=[n_fp, 1],
                           width_ratios=[1.2] + [1.0] * n_grp + [0.08],
                           hspace=0.12, wspace=0.04,
                           top=0.88, bottom=0.12, left=0.18, right=0.97)

    all_v = np.concatenate([fp_pat, group_vals[np.isfinite(group_vals)]])
    vabs = max(float(np.nanpercentile(np.abs(all_v), 97)), 0.05)
    norm = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
    cmap = "RdBu_r"

    ax_fp = fig.add_subplot(gs[0, 0])
    ax_fp.imshow(fp_pat.reshape(-1, 1), cmap=cmap, norm=norm,
                 aspect="auto", interpolation="nearest")
    ax_fp.set_yticks(range(n_fp))
    ax_fp.set_yticklabels(conn_labels, fontsize=max(3, 6 - n_fp // 10))
    ax_fp.set_xticks([0])
    ax_fp.set_xticklabels(["Expected"], fontsize=7, rotation=30, ha="right")
    ax_fp.set_title("Pattern", fontsize=8, fontweight="bold")

    ct_col = {"homotopic": "#009E73", "heterotopic": "#56B4E9", "ipsilateral": "#E69F00"}
    for k, ct in enumerate(conn_types):
        ax_fp.axhline(k + 0.5, color=ct_col.get(ct, "#aaa"), lw=0.4, alpha=0.5)

    sp_col_ranges = {}
    for gi, (sp, bd) in enumerate(groups):
        ax = fig.add_subplot(gs[0, gi + 1])
        ax.imshow(group_vals[:, gi].reshape(-1, 1), cmap=cmap, norm=norm,
                  aspect="auto", interpolation="nearest")
        ax.set_yticks([])
        ax.set_xticks([0])
        bd_lbl = bd.replace("BIDS_", "").replace("ds004513-download", "Hum")[:10]
        ax.set_xticklabels([bd_lbl], fontsize=max(4, 7 - n_grp // 4),
                           rotation=45, ha="right", color=sp_gc.get(sp, "#333"))
        sp_col_ranges.setdefault(sp, []).append(gi)

    for sp, col_list in sp_col_ranges.items():
        lo, hi = min(col_list), max(col_list)
        ax_lo = fig.axes[lo + 1]
        ax_hi = fig.axes[hi + 1]
        x_lo = ax_lo.get_position().x0
        x_hi = ax_hi.get_position().x1
        y_top = ax_lo.get_position().y1
        fig.text((x_lo + x_hi) / 2, y_top + 0.015, sp, ha="center", va="bottom",
                 fontsize=9, fontweight="bold", color=sp_gc.get(sp, "#333"),
                 transform=fig.transFigure)

    cbar_ax = fig.add_subplot(gs[0, -1])
    fig.colorbar(_mpl_cm.ScalarMappable(norm=norm, cmap=cmap),
                 cax=cbar_ax, label="Mean r")
    cbar_ax.yaxis.set_tick_params(labelsize=6)

    ax_bar = fig.add_subplot(gs[1, 1:n_grp + 1])
    bar_c = [sp_gc.get(sp, "#888") for sp, _ in groups]
    ax_bar.bar(np.arange(n_grp), group_frac * 100, color=bar_c, edgecolor="white", lw=0.4)

    for xi, v in enumerate(group_frac * 100):
        ax_bar.text(xi, v + 0.5, f"{v:.0f}%", ha="center", va="bottom",
                    fontsize=max(4, 7 - n_grp // 4), color="#333")

    ax_bar.set_xticks([])
    ax_bar.set_ylabel("% Specific runs", fontsize=7)
    ax_bar.set_ylim(0, 110)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)

    fig.legend(handles=[mpatches.Patch(color=ct_col[k], label=k) for k in ct_col],
               loc="lower left", bbox_to_anchor=(0.0, 0.0),
               ncol=3, fontsize=7, frameon=False)

    _save_or_show(fig, output_path)
    return fig

# =============================================================================
# FIGURE E — PROPORTION BARS
# =============================================================================

def plot_fp_category_bar(df_scored, cats, group_by="species", output_path=None):
    df = df_scored.copy()
    df["fp_category"] = cats.values if hasattr(cats, "values") else cats

    if group_by == "bids_dir" and "bids_dir" in df.columns:
        group_vals = sorted(df["bids_dir"].dropna().unique())
        col_src = df["bids_dir"]
    elif group_by == "subject" and "subject" in df.columns:
        df["_subj_lbl"] = df.get("species", "").astype(str) + "/" + df["subject"].astype(str)
        group_vals = []
        for sp in (_sp_order(df) if "species" in df.columns else []):
            for sub in sorted(df.loc[df["species"] == sp, "subject"].dropna().unique()):
                group_vals.append(f"{sp}/{sub}")
        if not group_vals:
            group_vals = sorted(df["_subj_lbl"].unique())
        col_src = df["_subj_lbl"]
    else:
        sp_vals = df["species"].dropna().unique().tolist() if "species" in df.columns else []
        group_vals = [s for s in _PHYLO if s in sp_vals] + [s for s in sp_vals if s not in _PHYLO]
        col_src = df["species"] if "species" in df.columns else pd.Series(range(len(df)), index=df.index)

    counts = {}
    for g in group_vals:
        gdf = df[col_src == g]["fp_category"].dropna()
        total = len(gdf)
        if total > 0:
            counts[g] = {c: 100.0 * (gdf == c).sum() / total for c in CATS}
        else:
            counts[g] = {c: 0.0 for c in CATS}

    n_g = len(group_vals)
    fig, ax = plt.subplots(figsize=(max(5, n_g * 1.2), 4.5))
    bottoms = np.zeros(n_g)
    x = np.arange(n_g)

    for cat in CATS:
        vals = np.array([counts[g][cat] for g in group_vals])
        ax.bar(x, vals, 0.7, bottom=bottoms, color=_CAT_COLORS[cat],
               label=cat, edgecolor="white", lw=0.5)
        for xi, (v, b) in enumerate(zip(vals, bottoms)):
            if v > 5:
                ax.text(xi, b + v / 2, f"{v:.0f}%", ha="center",
                        va="center", fontsize=7, color="white", fontweight="bold")
        bottoms += vals

    ax.set_xticks(x)
    ax.set_xticklabels(group_vals, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("% runs", fontsize=10)
    ax.set_ylim(0, 105)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22),
              ncol=len(CATS), fontsize=8, frameon=False)

    plt.tight_layout()
    _save_or_show(fig, output_path)
    return fig

# =============================================================================
# FIGURE F — MATRICES BY CATEGORY
# =============================================================================
def plot_fp_category_matrices(df_scored, cats,
                              bids_root_template=None,
                              atlas_name="EDNIxCSC", atlas_level=2,
                              use_lr=True, group_by="species",
                              output_path=None, figsize=None):
    """
    Figure F -- FC matrices by fp_category.

    FIXED: Handles varying ROI counts by finding common ROIs across all matrices.
    """
    import matplotlib.cm as _cm

    # Align cats with df_scored
    if cats is not None:
        if isinstance(cats, pd.Series):
            cats_aligned = cats.reindex(df_scored.index).fillna("No")
        else:
            cats_aligned = pd.Series(cats, index=df_scored.index, name='fp_category')

        df = df_scored.copy()
        df["fp_category"] = cats_aligned.values
    else:
        df = df_scored.copy()
        if "fp_category" not in df.columns:
            print("\n  [Fig F] No category information provided\n")
            return None

    # Define groups
    if group_by == "bids_dir" and "bids_dir" in df.columns:
        all_groups = sorted(df["bids_dir"].dropna().unique())

        def _group_mask(df_, g):
            return df_["bids_dir"] == g
    else:
        sp_vals = df["species"].dropna().unique().tolist() if "species" in df.columns else []
        all_groups = [s for s in _PHYLO if s in sp_vals] + [s for s in sp_vals if s not in _PHYLO]

        def _group_mask(df_, g):
            return df_["species"] == g if "species" in df_.columns else pd.Series(True, index=df_.index)

    cats_order = ["all", "Specific", "Unspecific", "No"]
    cat_label = {"all": "All", "Specific": "Specific",
                 "Unspecific": "Unspecific", "No": "No"}

    # FIRST PASS: Collect all matrices and their ROIs to find common ROIs
    all_matrices = {}  # (cat, g) -> list of (rois, mat)
    all_roi_sets = []  # collect all ROI sets to find intersection

    for cat in cats_order:
        all_matrices[cat] = {}
        for g in all_groups:
            g_mask = _group_mask(df, g)

            if cat == "all":
                subset = df[g_mask]
            else:
                subset = df[g_mask & (df["fp_category"] == cat)]

            matrices = []
            for idx, row in subset.iterrows():
                path = row.get("corr_matrix_path", None)
                if not path or not os.path.exists(str(path)):
                    continue
                try:
                    rois, mat = _load_matrix(path)
                    matrices.append((rois, mat))
                    all_roi_sets.append(set(rois))
                except Exception:
                    pass

            all_matrices[cat][g] = matrices
            print(f"  [Fig F] {group_by}={g:<22s} cat={cat:<12s} "
                  f"collected {len(matrices)} matrices")

    # Find common ROIs across ALL matrices (not just per category)
    if all_roi_sets:
        common_rois = list(set.intersection(*all_roi_sets) if len(all_roi_sets) > 1 else all_roi_sets[0])
        common_rois.sort()
        print(f"  [Fig F] Common ROIs across all matrices: {len(common_rois)}")
    else:
        common_rois = []

    if len(common_rois) < 4:
        print("  [Fig F] Not enough common ROIs - cannot create figure")
        return None

    # SECOND PASS: Align matrices to common ROIs and average
    results = {}
    for cat in cats_order:
        results[cat] = {}
        for g in all_groups:
            matrices = all_matrices[cat][g]

            if not matrices:
                results[cat][g] = (common_rois, None, 0, 0)
                continue

            # Align each matrix to common ROIs
            aligned_mats = []
            for rois, mat in matrices:
                # Find indices of common ROIs in this matrix
                idx_map = [rois.index(r) for r in common_rois if r in rois]
                if len(idx_map) == len(common_rois):
                    aligned_mats.append(mat[np.ix_(idx_map, idx_map)])

            n_loaded = len(aligned_mats)
            mean_mat = np.nanmean(np.stack(aligned_mats, 0), 0) if aligned_mats else None

            # Get count of runs in df (for display)
            g_mask = _group_mask(df, g)
            if cat == "all":
                n_in_df = df[g_mask].shape[0]
            else:
                n_in_df = df[g_mask & (df["fp_category"] == cat)].shape[0]

            results[cat][g] = (common_rois, mean_mat, n_loaded, n_in_df)

            if n_in_df > 0:
                pct = 100.0 * n_loaded / n_in_df if n_in_df > 0 else 0
                print(f"  [Fig F] {group_by}={g:<22s} cat={cat:<12s} "
                      f"n_df={n_in_df:3d} n_loaded={n_loaded:3d} ({pct:.0f}%)")

    # Get reference ROIs
    rois_ref = common_rois
    n_rois = len(rois_ref)
    n_grp = len(all_groups)
    n_rows = len(cats_order)

    # Create ROI labels (shortened for display)
    short = [r.replace("L_", "L ").replace("R_", "R ")[:13] for r in rois_ref]

    # Determine color limits from all loaded matrices
    all_vals = []
    for cat in cats_order:
        for g in all_groups:
            _, mean_mat, _, _ = results[cat][g]
            if mean_mat is not None:
                all_vals.extend(mean_mat.ravel())

    all_vals = np.array(all_vals)
    all_vals = all_vals[np.isfinite(all_vals)]

    if len(all_vals) > 0:
        vmax_r = max(float(np.nanpercentile(np.abs(all_vals), 97)), 0.05)
    else:
        vmax_r = 0.5

    norm_r = TwoSlopeNorm(vmin=-vmax_r, vcenter=0, vmax=vmax_r)

    # Create figure
    w, h = figsize or (max(8, n_grp * 3.5), n_rows * 3.5)
    fig = plt.figure(figsize=(w, h))
    gs = gridspec.GridSpec(n_rows, n_grp + 1,
                           width_ratios=[1.0] * n_grp + [0.04],
                           wspace=0.05, hspace=0.35)
    # Plot each category and group
    for r_i, cat in enumerate(cats_order):
        col_c = _CAT_COLORS.get(cat, "#444") if cat != "all" else "#444"

        for c_i, g in enumerate(all_groups):
            ax = fig.add_subplot(gs[r_i, c_i])
            rois_g, mean_mat, n_loaded, n_in_df = results[cat][g]

            if mean_mat is None or n_loaded == 0:
                ax.set_facecolor("#f2f2f2")
                ax.set_xticks([])
                ax.set_yticks([])
                for sp in ax.spines.values():
                    sp.set_visible(False)

                if n_in_df > 0:
                    ax.text(0.5, 0.5, f"0/{n_in_df} loaded",
                            ha="center", va="center",
                            fontsize=7, color="#999999", transform=ax.transAxes)
                else:
                    ax.text(0.5, 0.5, "no data",
                            ha="center", va="center",
                            fontsize=7, color="#bbbbbb", transform=ax.transAxes)
            else:
                im = ax.imshow(mean_mat, cmap="RdBu_r", norm=norm_r,
                               aspect="equal", interpolation="nearest")

                ax.text(0.02, 0.98, f"n={n_loaded}",
                        transform=ax.transAxes,
                        fontsize=6, va="top", ha="left",
                        color="#333333", backgroundcolor='white', alpha=0.7)

                if c_i == 0:
                    step = max(1, n_rois // 30)
                    labels_to_show = short[::step]
                    positions = range(0, n_rois, step)
                    ax.set_yticks(positions)
                    ax.set_yticklabels(labels_to_show, fontsize=3)
                else:
                    ax.set_yticks([])

                ax.set_xticks([])

            if c_i == 0:
                ax.set_ylabel(cat_label.get(cat, cat), fontsize=9,
                              color=col_c, fontweight="bold", labelpad=4)
            if r_i == 0:
                ax.set_title(str(g), fontsize=8, fontweight="bold")

    cbar_ax = fig.add_subplot(gs[:, -1])
    fig.colorbar(_cm.ScalarMappable(norm=norm_r, cmap="RdBu_r"),
                 cax=cbar_ax, label="Mean r")

    plt.tight_layout()
    _save_or_show(fig, output_path)
    return fig
# =============================================================================
# FULL PIPELINE
# =============================================================================

def run_pipeline(df_qc,
                 fit_kind='correlation',
                 bids_root_template=_DEFAULT_BIDS,
                 method="auto",
                 thresh_intra=None,
                 thresh_inter=None,
                 stringency=0.0,
                 n_top=10,
                 atlas_name="EDNIxCSC", atlas_level=2, use_lr=True,
                 fig_dir=None,
                 verbose=True):
    print("\n" + "=" * 60 + "\n  EDNiX FC Pipeline\n" + "=" * 60)

    def _p(name):
        if fig_dir is None:
            return None
        os.makedirs(fig_dir, exist_ok=True)
        return os.path.join(fig_dir, name)

    print("\n[1] Scoring all runs...")
    df_scored = score_all_runs(df_qc, fit_kind, bids_root_template=bids_root_template,
                               atlas_name=atlas_name, atlas_level=atlas_level,
                               use_lr=use_lr, verbose=verbose)

    print("\n[2] Defining thresholds (2D GMM)...")
    thresh = suggest_thresholds(df_scored, method=method,
                                thresh_intra=thresh_intra,
                                thresh_inter=thresh_inter,
                                thresh_delta=thresh_delta,      # ← new kwarg
                                stringency=stringency,
                                verbose=verbose)

    cats = classify_runs(df_scored,
                         thresh["thresh_intra"],
                         thresh["thresh_inter"],
                         thresh_delta=thresh["thresh_delta"],   # ← new kwarg
                         verbose=verbose)
    df_scored["fp_category"] = cats.values

    print("\n[4] Figure: threshold distributions...")
    plot_threshold_distributions(df_scored,
                                 thresh["thresh_intra"], thresh["thresh_inter"],
                                 output_path=_p("fingerprint_threshold_distributions.png"))

    print("\n[5] Figure C: classification model...")
    plot_classification_model(df_scored, thresh["thresh_intra"],
                              thresh["thresh_inter"], cats,
                              output_path=_p("fingerprint_C_model.png"))

    print("\n[6] Figure E: proportion bars...")
    plot_fp_category_bar(df_scored, cats, group_by="species",
                         output_path=_p("fingerprint_E_bar_species.png"))
    plot_fp_category_bar(df_scored, cats, group_by="bids_dir",
                         output_path=_p("fingerprint_E_bar_bids.png"))

    print("\n[7] Figure F: matrices by category...")
    plot_fp_category_matrices(df_scored, cats,  # Make sure cats is passed correctly
                              bids_root_template=bids_root_template,
                              atlas_name=atlas_name, atlas_level=atlas_level, use_lr=use_lr,
                              group_by="species",
                              output_path=_p("fingerprint_F_matrices_by_species.png"))
    plot_fp_category_matrices(df_scored, cats,
                              bids_root_template=bids_root_template,
                              atlas_name=atlas_name, atlas_level=atlas_level, use_lr=use_lr,
                              group_by="bids_dir",
                              output_path=_p("fingerprint_F_matrices_by_bids.png"))

    print("\n[8] Building fingerprint from Specific runs...")
    fp = None
    sp_mats = {}
    try:
        sp_mats, fp = build_fingerprint_from_specific(
            df_scored, cats, bids_root_template=bids_root_template,
            atlas_name=atlas_name, atlas_level=atlas_level,
            use_lr=use_lr, n_top=n_top, verbose=verbose)

        print("\n[9] Figure B: fingerprint...")
        plot_fc_fingerprint(fp, output_path=_p("fingerprint_B_discovery.png"))

        print("\n[10] Figure D: heatmap...")
        plot_fingerprint_heatmap(df_scored, fp, cats,
                                 bids_root_template=bids_root_template,
                                 atlas_name=atlas_name, atlas_level=atlas_level, use_lr=use_lr,
                                 output_path=_p("fingerprint_D_heatmap.png"))
    except ValueError as e:
        print(f"  [skip B+D] {e}")

    print("\n" + "=" * 60 + "\n  Pipeline complete.")
    if fig_dir:
        print(f"  Figures saved to: {fig_dir}")
    print("=" * 60 + "\n")

    return dict(df_scored=df_scored, thresh=thresh, cats=cats, fp=fp, sp_mats=sp_mats)

# =============================================================================
# CLI
# =============================================================================
def main():
    p=argparse.ArgumentParser(description="EDNiX FC pipeline")
    p.add_argument("--csv",required=True)
    p.add_argument("--bids-root",default=_DEFAULT_BIDS)
    p.add_argument("--method",default="auto",choices=["auto","manual"])
    p.add_argument("--thresh-intra",type=float,default=None)
    p.add_argument("--thresh-inter",type=float,default=None)
    p.add_argument("--stringency",type=float,default=0.0)
    p.add_argument("--n-top",type=int,default=10)
    p.add_argument("--atlas-level",type=int,default=2)
    p.add_argument("--output",default="figures/qc_recap")
    args=p.parse_args()
    df_qc=pd.read_csv(args.csv)
    run_pipeline(df_qc,bids_root_template=args.bids_root,
                  method=args.method,thresh_intra=args.thresh_intra,
                  thresh_inter=args.thresh_inter,stringency=args.stringency,
                  n_top=args.n_top,atlas_level=args.atlas_level,
                  fig_dir=args.output)

if __name__=="__main__":
    main()

