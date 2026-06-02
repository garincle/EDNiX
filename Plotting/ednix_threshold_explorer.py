"""
EDNiX — Primary Network FC Threshold Explorer & QC Classifier
==============================================================

Each run_*.csv file is treated as an independent observation at the scoring
stage. One df_qc row = one BIDS session, expanded into N run rows.

Scoring uses MEDIAN aggregation (robust to outlier ROI pairs).
Classification uses the delta rule:

    intra        : median homotopic same-network r   (L_Motor <-> R_Motor, ...)
    inter_hetero : median contralateral cross-network r
                   (L_netA <-> R_netB, R_netA <-> L_netB)
                   excludes ipsilateral cross-network (leakage risk)
    delta        : intra - inter_hetero

    Spurious   : inter_hetero > intra            (delta < 0)
    No         : intra <= thresh_intra
    Unspecific : intra > thresh_intra AND delta <= thresh_delta
    Specific   : intra > thresh_intra AND delta >  thresh_delta

Thresholds are chosen as PERCENTILES of the pooled cross-species distribution.

Fingerprint significance uses a proper hierarchical test (species_mean t-test
or a linear mixed model with species/subject random effects) + Benjamini-Hochberg
FDR, instead of the old pooled one-sample t-test (which was pseudoreplicated).
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

# Taxonomic subsets. None == all species.
_PRIMATES = ["Mouselemur", "Marmoset", "Macaque", "Human"]
_RODENTS = ["Mouse", "Rat"]
_SUBSETS = {
    "all": None,
    "primates": _PRIMATES,
    "primates_rodents": _PRIMATES + _RODENTS,
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

def _inter_col(df):
    """Pick the heterotopic inter column if present, else legacy inter."""
    return ("primary_inter_hetero"
            if "primary_inter_hetero" in df.columns else "primary_inter")

def _bilateral_pair_filter(rois):
    """
    Keep only ROIs whose L_/R_ counterpart is also present.
    Returns in order: L_a, L_b, ..., R_a, R_b, ... (all Ls then all Rs,
    sorted by base name) — guarantees a balanced bilateral matrix layout.
    """
    l_bases = {n[2:] for n in rois if n.startswith(("L_", "l_"))}
    r_bases = {n[2:] for n in rois if n.startswith(("R_", "r_"))}
    common = sorted(l_bases & r_bases)
    return [f"L_{b}" for b in common] + [f"R_{b}" for b in common]


def _common_bilateral_rois(matrices_with_rois):
    """
    Given a list of (rois, mat) tuples, return the common bilateral ROI list
    (intersection across all matrices, then bilateral-paired).
    """
    if not matrices_with_rois:
        return []
    common = set(matrices_with_rois[0][0])
    for rois, _ in matrices_with_rois[1:]:
        common &= set(rois)
    return _bilateral_pair_filter([r for r in common])


def compute_subset_common_rois(df_scored, species_subset=None, verbose=True):
    """
    Load all (or subset) corr matrices listed in df_scored and return the
    common bilateral ROI list — the same list is then used for fig 04 and
    fingerprint figures so their ROI sets are guaranteed identical.

    When verbose=True, prints a per-species and global report of which ROIs
    are dropped and why:
      - missing in species X (present elsewhere but not in any X matrix)
      - has no L/R counterpart  (homotopic bilateral pairing failed)
    """
    df = df_scored.copy()
    if species_subset is not None and "species" in df.columns:
        df = df[df["species"].isin(species_subset)]

    # collect per-species ROI sets (union of all matrices for that species)
    species_rois = {}  # species -> set of ROIs found in any of its matrices
    n_loaded_per_sp = {}
    for _, row in df.iterrows():
        path = row.get("corr_matrix_path", None)
        sp = str(row.get("species", "?"))
        if not path or not os.path.exists(str(path)):
            continue
        try:
            rois, _ = _load_matrix(path)
            species_rois.setdefault(sp, set()).update(rois)
            n_loaded_per_sp[sp] = n_loaded_per_sp.get(sp, 0) + 1
        except Exception:
            pass

    if not species_rois:
        return []

    # universe = union across all species (anything we ever saw)
    universe = set().union(*species_rois.values())

    # raw common = strict intersection of species-level unions
    raw_common = set.intersection(*species_rois.values()) if len(species_rois) > 1 \
        else set(next(iter(species_rois.values())))

    # bilateral-paired final list
    paired = _bilateral_pair_filter(list(raw_common))
    paired_set = set(paired)

    if verbose:
        print(f"  [common ROIs] species: {sorted(species_rois.keys())}")
        for sp in sorted(species_rois.keys()):
            n_sp = len(species_rois[sp])
            n_mat = n_loaded_per_sp.get(sp, 0)
            print(f"     {sp:14s} {n_sp:4d} unique ROIs across {n_mat} matrices")

        # ROIs in universe but missing from at least one species
        missing_in_species = sorted(universe - raw_common)
        # ROIs that survived intersection but lost their L/R partner
        unpaired = sorted(raw_common - paired_set)

        print(f"  [common ROIs] universe={len(universe)}  "
              f"intersection={len(raw_common)}  "
              f"bilateral-paired={len(paired)}")

        if missing_in_species:
            print(f"  [dropped: missing in >=1 species] ({len(missing_in_species)})")
            # group by which species are missing each ROI
            for roi in missing_in_species[:60]:  # cap output
                missing_from = sorted(sp for sp, rs in species_rois.items()
                                      if roi not in rs)
                print(f"     {roi:42s}  missing in: {', '.join(missing_from)}")
            if len(missing_in_species) > 60:
                print(f"     ... ({len(missing_in_species) - 60} more)")

        if unpaired:
            print(f"  [dropped: no L/R counterpart] ({len(unpaired)})")
            for roi in unpaired[:30]:
                print(f"     {roi}")
            if len(unpaired) > 30:
                print(f"     ... ({len(unpaired) - 30} more)")

        n_l = sum(1 for r in paired if r.startswith("L_"))
        n_r = sum(1 for r in paired if r.startswith("R_"))
        print(f"  [common ROIs] FINAL: {len(paired)} bilateral "
              f"({n_l} L + {n_r} R)")

    return paired


def _bh_fdr(pvals):
    """Benjamini-Hochberg FDR. NaNs are passed through as NaN."""
    p = np.asarray(pvals, float)
    out = np.full(p.shape, np.nan)
    finite = np.isfinite(p)
    if not finite.any():
        return out
    pv = p[finite]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    qf = np.empty(n)
    qf[order] = np.clip(q, 0, 1)
    out[finite] = qf
    return out

def _find_all_run_paths(bids_root, subject, session, fit_kind,
                        atlas_name="EDNIxCSC", atlas_level=2, use_lr=True):
    """
    Match files like  {atlas}{LR}_{level}_run_{N}[_flattened]_{fit_kind}_matrix.csv
    where {N} is an integer and {fit_kind} matches EXACTLY (no leading words like
    "partial " that would be inside another fit_kind). The strict end-anchor
    `_{fit_kind}_matrix.csv$` is enforced by regex, so "correlation" no longer
    accidentally matches "partial correlation".

    Files containing pval / tstat / check_fit in the name are excluded.
    "flattened" is allowed (it's a legitimate suffix produced by the pipeline).
    """
    lr = "LR" if use_lr else ""
    # Escape fit_kind for regex (it may contain spaces, e.g. "partial correlation")
    fk_re = re.escape(fit_kind)
    # Filename must end exactly with _{fit_kind}_matrix.csv
    name_re = re.compile(
        rf"^{re.escape(atlas_name + lr)}_{atlas_level}_run_(\d+)"
        rf"(?:_flattened)?_{fk_re}_matrix\.csv$"
    )
    excl = ("check_fit", "pval", "tstat")
    sub_dir = os.path.join(bids_root, f"sub-{subject}")

    ses_dirs = []
    if session:
        ses_dirs.append(os.path.join(sub_dir, f"ses-{session}"))
    ses_dirs.append(sub_dir)

    all_files = []
    for ses_d in ses_dirs:
        for mid in [
            os.path.join("func", "acpc-func", "Stats", "Correl_matrix"),
            os.path.join("func", "Stats", "Correl_matrix"),
            os.path.join("func", "acpc-func", "Correl_matrix"),
            os.path.join("func", "Correl_matrix"),
        ]:
            d = os.path.join(ses_d, mid)
            if not os.path.isdir(d):
                continue
            for fn in sorted(os.listdir(d)):
                if any(x in fn for x in excl):
                    continue
                if name_re.match(fn):
                    all_files.append(os.path.join(d, fn))

    seen, unique_files = set(), []
    for f in all_files:
        if f not in seen:
            seen.add(f)
            unique_files.append(f)
    return unique_files

def _get_all_run_paths_for_row(fit_kind, row, bids_root_template,
                                atlas_name="EDNIxCSC", atlas_level=2, use_lr=True):
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

def _classify_one(intra, inter_h, t_intra, t_delta):
    """Single-observation delta classification (shared by all classifiers)."""
    if not np.isfinite(intra) or not np.isfinite(inter_h):
        return "No"
    delta = intra - inter_h
    if inter_h > intra:
        return "Spurious"
    if intra <= t_intra:
        return "No"
    if delta <= t_delta:
        return "Unspecific"
    return "Specific"

# =============================================================================
# PRIMARY NETWORK SCORE  (median + heterotopic inter)
# =============================================================================

def compute_primary_score_from_matrix(mat, rois):
    """
    intra        : median homotopic same-network r
    inter        : median ALL cross-network pairs (reference)
    inter_hetero : median contralateral cross-network pairs only
    """
    primary = _find_primary_rois(rois)
    nets = list(primary.keys())
    if not nets:
        return dict(intra_means={}, mean_intra=np.nan, mean_inter=np.nan,
                    mean_inter_hetero=np.nan, delta=np.nan,
                    contrast_ratio=np.nan, networks_found=[])

    intra_vals = {net: float(mat[v["L"], v["R"]])
                  if np.isfinite(mat[v["L"], v["R"]]) else np.nan
                  for net, v in primary.items()}
    median_intra = float(np.nanmedian(list(intra_vals.values())))

    inter_all, inter_hetero = [], []
    for i in range(len(nets)):
        for j in range(i + 1, len(nets)):
            na, nb = nets[i], nets[j]
            ia_L, ia_R = primary[na]["L"], primary[na]["R"]
            ib_L, ib_R = primary[nb]["L"], primary[nb]["R"]
            for ia, ib in [(ia_L, ib_L), (ia_L, ib_R), (ia_R, ib_L), (ia_R, ib_R)]:
                r = float(mat[ia, ib])
                if np.isfinite(r):
                    inter_all.append(r)
            for ia, ib in [(ia_L, ib_R), (ia_R, ib_L)]:   # contralateral only
                r = float(mat[ia, ib])
                if np.isfinite(r):
                    inter_hetero.append(r)

    median_inter = float(np.nanmedian(inter_all)) if inter_all else np.nan
    median_inter_hetero = float(np.nanmedian(inter_hetero)) if inter_hetero else np.nan
    delta = (float(median_intra - median_inter_hetero)
             if np.isfinite(median_inter_hetero) else np.nan)

    denom = median_inter_hetero
    if np.isnan(median_intra) or np.isnan(denom):
        ratio = np.nan
    elif denom <= 0:
        ratio = 20.0
    else:
        ratio = median_intra / denom

    return dict(intra_means=intra_vals, mean_intra=median_intra,
                mean_inter=median_inter, mean_inter_hetero=median_inter_hetero,
                delta=delta, contrast_ratio=ratio, networks_found=nets)

# =============================================================================
# STEP 1 — SCORE ALL RUNS
# =============================================================================

def score_all_runs(df_qc, fit_kind,
                   bids_root_template=_DEFAULT_BIDS,
                   atlas_name="EDNIxCSC", atlas_level=2,
                   use_lr=True, verbose=True):
    all_runs, n_scored = [], 0
    for idx, row in df_qc.iterrows():
        run_infos = _get_all_run_paths_for_row(
            fit_kind, row, bids_root_template, atlas_name, atlas_level, use_lr)
        for run_path, sp, bd, sub, ses in run_infos:
            try:
                rois, mat = _load_matrix(run_path)
                ps = compute_primary_score_from_matrix(mat, rois)
                run_row = {
                    "species": sp, "bids_dir": bd, "subject": sub, "session": ses,
                    "run": os.path.basename(run_path),
                    "corr_matrix_path": run_path,
                    "primary_intra": ps["mean_intra"],
                    "primary_inter": ps["mean_inter"],
                    "primary_inter_hetero": ps["mean_inter_hetero"],
                    "primary_delta": ps["delta"],
                    "contrast_ratio": ps["contrast_ratio"],
                }
                for net, val in ps["intra_means"].items():
                    run_row[f"per_net_intra_{net}"] = val
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
        for c in ["primary_intra", "primary_delta"]:
            v = df_scored[c].dropna()
            if len(v):
                print(f"  [score] {c:14s}: {v.min():.3f}--{v.max():.3f} "
                      f"median={v.median():.3f}")
    return df_scored

# =============================================================================
# STEP 2 — THRESHOLDS (percentile-based, on pooled cross-species data)
# =============================================================================

def suggest_thresholds(df_scored, method="auto",
                       thresh_intra=None, thresh_delta=None,
                       stringency=0.0, verbose=True):
    """
    Percentile thresholds on the POOLED cross-species distribution.

        intra threshold : 25th percentile of intra   (stringency 0)
        delta threshold : 25th percentile of delta    (stringency 0)

    stringency > 0 moves toward the median (stricter),
    stringency < 0 moves toward the 10th percentile (more lenient).
    delta may legitimately be negative, so it is not clipped to >= 0.
    """
    inter_col = _inter_col(df_scored)
    delta = df_scored["primary_intra"] - df_scored[inter_col]
    valid = df_scored["primary_intra"].notna() & df_scored[inter_col].notna()
    intra_vals = df_scored.loc[valid, "primary_intra"].values
    delta_vals = delta.loc[valid].values

    if len(intra_vals) < 10:
        t_intra, t_delta = 0.2, 0.0
        auto_intra, auto_delta = t_intra, t_delta
    else:
        s = np.clip(stringency, -1.0, 1.0)
        # base = median (P50). symmetric anchors: P75 strict, P25 lenient.
        base_intra = np.percentile(intra_vals, 50)
        base_delta = np.percentile(delta_vals, 50)
        if s >= 0:
            auto_intra = base_intra + s * (np.percentile(intra_vals, 75) - base_intra)
            auto_delta = base_delta + s * (np.percentile(delta_vals, 75) - base_delta)
        else:
            auto_intra = base_intra + abs(s) * (np.percentile(intra_vals, 25) - base_intra)
            auto_delta = base_delta + abs(s) * (np.percentile(delta_vals, 25) - base_delta)
        auto_intra = float(np.clip(auto_intra, 0.05, 0.8))
        auto_delta = float(auto_delta)
        t_intra = thresh_intra if thresh_intra is not None else auto_intra
        t_delta = thresh_delta if thresh_delta is not None else auto_delta

    n_spec = ((intra_vals > t_intra) & (delta_vals > t_delta)).sum() if len(intra_vals) else 0
    pct = 100.0 * n_spec / len(intra_vals) if len(intra_vals) else 0

    if verbose:
        print(f"\n{'=' * 50}")
        print(f"  THRESHOLD SELECTION  (stringency={stringency:+.1f})")
        print(f"  intra > {t_intra:.3f}")
        print(f"  delta > {t_delta:.3f}   (intra - inter_hetero)")
        print(f"  -> {n_spec}/{len(intra_vals)} runs Specific ({pct:.0f}%)")
        print(f"{'=' * 50}\n")

    return dict(thresh_intra=t_intra, thresh_delta=t_delta,
                auto_intra=auto_intra, auto_delta=auto_delta,
                stringency=stringency)

# =============================================================================
# STEP 3 — CLASSIFY
# =============================================================================

def classify_runs(df_scored, thresh_intra, thresh_delta, verbose=True):
    inter_col = _inter_col(df_scored)
    cats = pd.Series("No", index=df_scored.index, dtype=object)
    intra_arr = df_scored["primary_intra"].values.astype(float)
    inter_arr = df_scored[inter_col].values.astype(float)
    for i, idx in enumerate(df_scored.index):
        cats[idx] = _classify_one(intra_arr[i], inter_arr[i], thresh_intra, thresh_delta)

    if verbose:
        print(f"  Classification  thresh_intra={thresh_intra:.4f} "
              f"thresh_delta={thresh_delta:.4f} (col={inter_col})")
        if "species" in df_scored.columns:
            for sp in _sp_order(df_scored):
                mask = df_scored["species"] == sp
                if not mask.any():
                    continue
                c = cats[mask].value_counts()
                print(f"    {sp:14s} n={mask.sum():3d}   " +
                      "  ".join(f"{cat}={c.get(cat, 0)}" for cat in CATS))
    return cats

# =============================================================================
# FIGURE 01 — THRESHOLD SELECTION DISTRIBUTIONS
# =============================================================================

def plot_threshold_selection(df_scored, thresh_intra=None, thresh_delta=None,
                             output_path=None):
    """
    Panel A : intra distribution (faded per-BIDS) + global pooled density,
              percentiles P25/P50/P75, and the chosen threshold (with its
              percentile location).
    Panel B : same, for delta = intra - inter_hetero.
    Panel C : 2D classification space (intra vs inter_hetero) with the delta
              decision boundary. Kept as before.
    """
    inter_col = _inter_col(df_scored)
    df = df_scored.copy()
    df["primary_delta"] = df["primary_intra"] - df[inter_col]

    if thresh_intra is None or thresh_delta is None:
        t = suggest_thresholds(df, verbose=False)
        thresh_intra = thresh_intra if thresh_intra is not None else t["thresh_intra"]
        thresh_delta = thresh_delta if thresh_delta is not None else t["thresh_delta"]

    bids_list = sorted(df["bids_dir"].dropna().unique()) if "bids_dir" in df else []
    sp_list = _sp_order(df)
    sp_colors = _gcolors(sp_list)
    # map each BIDS to its dominant species (for color)
    bids_species = {}
    if "species" in df.columns:
        for bd in bids_list:
            sp_counts = df.loc[df["bids_dir"] == bd, "species"].value_counts()
            bids_species[bd] = sp_counts.index[0] if len(sp_counts) else None

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(f"Threshold selection (percentile-based) — "
                 f"intra>{thresh_intra:.3f}, delta>{thresh_delta:.3f}",
                 fontsize=12, fontweight="bold")

    for ax, col, thr, xlabel, title in [
        (axes[0], "primary_intra", thresh_intra,
         "Primary intra-network r (homotopic)", "Bilateral connectivity"),
        (axes[1], "primary_delta", thresh_delta,
         "Δ = intra − inter(hetero)", "Network segregation (Δ)"),
    ]:
        # faded per-BIDS densities, colored by species
        for bd in bids_list:
            v = df.loc[df["bids_dir"] == bd, col].dropna().values
            v = v[np.isfinite(v)]
            if len(v) < 5:
                continue
            x = np.linspace(v.min(), v.max(), 300)
            try:
                k = _stats.gaussian_kde(v, bw_method="silverman")
                sp = bids_species.get(bd)
                c = sp_colors.get(sp, "#888888") if sp else "#888888"
                ax.plot(x, k(x), color=c, lw=1.0, alpha=0.30, zorder=1)
                ax.fill_between(x, k(x), color=c, alpha=0.05, zorder=1)
            except Exception:
                pass

        # global pooled density (foreground)
        g = df[col].dropna().values
        g = g[np.isfinite(g)]
        if len(g) < 5:
            continue
        xg = np.linspace(g.min(), g.max(), 500)
        kg = _stats.gaussian_kde(g, bw_method="silverman")
        ax.fill_between(xg, kg(xg), color="#0072B2", alpha=0.12, zorder=2)
        ax.plot(xg, kg(xg), color="#0072B2", lw=2.2, zorder=3,
                label="global (pooled)")

        ymax = ax.get_ylim()[1]
        # percentiles
        for p, ls in [(25, ":"), (50, "--"), (75, ":")]:
            xp = np.percentile(g, p)
            ax.axvline(xp, color="#555555", lw=1.0, ls=ls, zorder=4)
            ax.text(xp, ymax * 0.95, f"P{p}", rotation=90, fontsize=6,
                    color="#555555", ha="right", va="top")
        # chosen threshold + its percentile location
        pctl = _stats.percentileofscore(g, thr)
        ax.axvline(thr, color="#D55E00", lw=2.4, zorder=5)
        ax.text(thr, ymax * 0.99, f" thr={thr:.3f}\n (≈P{pctl:.0f})",
                color="#D55E00", fontsize=8, ha="left", va="top",
                fontweight="bold", zorder=6)

        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel("Density", fontsize=9)
        ax.set_title(title, fontsize=10, fontweight="bold")

        # Two legends: thresholds/global (upper right), species swatches (upper left)
        h0, l0 = ax.get_legend_handles_labels()
        leg1 = ax.legend(h0, l0, fontsize=7, frameon=False, loc="upper right")
        ax.add_artist(leg1)
        if "species" in df.columns and sp_list:
            sp_handles = [mpatches.Patch(color=sp_colors[sp], alpha=0.5, label=sp)
                          for sp in sp_list]
            ax.legend(handles=sp_handles, fontsize=6, frameon=False,
                      loc="upper left", ncol=1, title="species",
                      title_fontsize=6)

        ax.set_ylim(bottom=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Panel C — 2D classification space (kept)
    ax = axes[2]
    gc = sp_colors  # reuse the species palette from above
    markers = ["o", "s", "^", "D", "v", "P", "*", "h"]
    for i, sp in enumerate(sp_list):
        mask = df["species"] == sp if "species" in df.columns else pd.Series(True, index=df.index)
        x = df.loc[mask, "primary_intra"].values.astype(float)
        y = df.loc[mask, inter_col].values.astype(float)
        fin = np.isfinite(x) & np.isfinite(y)
        ax.scatter(x[fin], y[fin], color=gc[sp], marker=markers[i % len(markers)],
                   s=12, alpha=0.5, edgecolors="none", label=sp, zorder=3)

    xmin = float(np.nanmin(df["primary_intra"].dropna())) - 0.01
    xmax = float(np.nanmax(df["primary_intra"].dropna())) + 0.01
    ymin = float(np.nanmin(df[inter_col].dropna())) - 0.01
    ymax = float(np.nanmax(df[inter_col].dropna())) + 0.01
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axvline(thresh_intra, color="#222", lw=1.8, ls="--", zorder=5)
    # delta decision boundary: inter_hetero = intra - delta_thresh
    xs = np.linspace(xmin, xmax, 100)
    ax.plot(xs, xs - thresh_delta, color="#444", lw=1.4, ls=":",
            zorder=5, label=f"Δ={thresh_delta:.2f} boundary")
    ax.plot(xs, xs, color="#bbb", lw=1.0, zorder=2, label="inter=intra")

    ax.set_xlabel("Primary intra-network r", fontsize=10)
    ax.set_ylabel("Primary inter-network r (hetero)", fontsize=10)
    ax.set_title("2D classification space", fontsize=10, fontweight="bold")
    ax.legend(fontsize=7, frameon=False, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    _save_or_show(fig, output_path)
    return fig

# =============================================================================
# FIGURE 02 — CLASSIFICATION MODEL (per species; uses delta)
# =============================================================================

def plot_classification_model(df_scored, thresh_intra, thresh_delta, cats,
                              output_path=None):
    inter_col = _inter_col(df_scored)
    df = df_scored.copy()
    df["primary_delta"] = df["primary_intra"] - df[inter_col]

    sp_list = _sp_order(df)
    n_sp = len(sp_list)
    if n_sp == 0:
        return None
    gc = _gcolors(sp_list)
    rng = np.random.default_rng(42)

    fig = plt.figure(figsize=(max(n_sp * 3.0, 9), 10))
    gs = gridspec.GridSpec(3, n_sp, wspace=0.06, hspace=0.40,
                           top=0.92, bottom=0.07, left=0.09, right=0.98)

    for col, sp in enumerate(sp_list):
        sp_df = df[df["species"] == sp] if "species" in df.columns else df
        sp_cats = cats[sp_df.index]
        color = gc[sp]
        c_ct = sp_cats.value_counts()

        for row, (metric, thresh, ylabel) in enumerate([
            ("primary_intra", thresh_intra, "Primary intra-network r"),
            ("primary_delta", thresh_delta, "Δ = intra − inter(hetero)"),
        ]):
            ax = fig.add_subplot(gs[row, col])
            av = sp_df[metric].dropna() if metric in sp_df else pd.Series(dtype=float)
            if av.empty:
                if row == 0:
                    ax.set_title(f"{sp}\n(n={len(sp_df)} runs)", fontsize=8,
                                 fontweight="bold", color=color)
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
            ax.scatter(jit, av.values, c=pt_c, s=18, alpha=0.75, zorder=3,
                       edgecolors="none")
            if len(av) >= 2:
                q1, med, q3 = np.nanpercentile(av.values, [25, 50, 75])
                ax.plot([-0.15, 0.15], [med, med], color=color, lw=2.0, zorder=5)
                ax.plot([0, 0], [q1, q3], color=color, lw=1.2, zorder=4)
            ax.axhline(thresh, color="#333", lw=1.5, ls="--", zorder=6)
            ax.axhline(0, color="#ddd", lw=0.4)
            ax.set_xticks([])
            ax.tick_params(labelsize=7)
            for s in ("top", "right", "bottom"):
                ax.spines[s].set_visible(False)
            if col == 0:
                ax.set_ylabel(ylabel, fontsize=8)
            if row == 0:
                counts = " ".join(f"{cat[0]}={c_ct.get(cat, 0)}" for cat in CATS)
                ax.set_title(f"{sp} (n={len(sp_df)} runs)\n{counts}",
                             fontsize=8, fontweight="bold", color=color)

        # bottom row: intra (x) vs delta (y); Specific = top-right quadrant
        ax = fig.add_subplot(gs[2, col])
        xv = sp_df["primary_intra"].dropna() if "primary_intra" in sp_df else pd.Series(dtype=float)
        yv = sp_df["primary_delta"].dropna() if "primary_delta" in sp_df else pd.Series(dtype=float)
        ci = xv.index.intersection(yv.index)
        if not ci.empty:
            pt_c2 = [_CAT_COLORS.get(sp_cats[i], "#aaa") for i in ci]
            ax.scatter(xv[ci], yv[ci], c=pt_c2, s=15, alpha=0.65, zorder=3,
                       edgecolors="none")
        ax.axvline(thresh_intra, color="#333", lw=1.1, ls="--", zorder=5)
        ax.axhline(thresh_delta, color="#555", lw=0.9, ls=":", zorder=5)
        ax.axhline(0, color="#ddd", lw=0.4)
        ax.set_xlabel("Intra r", fontsize=7)
        if col == 0:
            ax.set_ylabel("Δ (intra − inter)", fontsize=7)
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
# FIGURE 03 — CATEGORY PROPORTION BARS
# =============================================================================

def plot_category_proportions(df_scored, cats, group_by="species", output_path=None):
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
        counts[g] = ({c: 100.0 * (gdf == c).sum() / total for c in CATS}
                     if total > 0 else {c: 0.0 for c in CATS})

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
                ax.text(xi, b + v / 2, f"{v:.0f}%", ha="center", va="center",
                        fontsize=7, color="white", fontweight="bold")
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
# FIGURE 04 — FC MATRICES BY CATEGORY (median across runs)
# =============================================================================

def plot_fc_matrices_by_category(df_scored, cats,
                                 bids_root_template=None,
                                 atlas_name="EDNIxCSC", atlas_level=2,
                                 use_lr=True, group_by="species",
                                 species_subset=None,
                                 common_rois=None,
                                 output_path=None, figsize=None):
    import matplotlib.cm as _cm

    if cats is not None:
        if isinstance(cats, pd.Series):
            cats_aligned = cats.reindex(df_scored.index).fillna("No")
        else:
            cats_aligned = pd.Series(cats, index=df_scored.index, name="fp_category")
        df = df_scored.copy()
        df["fp_category"] = cats_aligned.values
    else:
        df = df_scored.copy()
        if "fp_category" not in df.columns:
            print("\n  [Fig 04] No category information\n")
            return None

    if species_subset is not None and "species" in df.columns:
        df = df[df["species"].isin(species_subset)]
        if df.empty:
            print("  [Fig 04] subset empty — skipped")
            return None

    if group_by == "bids_dir" and "bids_dir" in df.columns:
        all_groups = sorted(df["bids_dir"].dropna().unique())
        def _gmask(d, g): return d["bids_dir"] == g
    else:
        sp_vals = df["species"].dropna().unique().tolist() if "species" in df.columns else []
        all_groups = [s for s in _PHYLO if s in sp_vals] + [s for s in sp_vals if s not in _PHYLO]
        def _gmask(d, g): return d["species"] == g if "species" in d.columns else pd.Series(True, index=d.index)

    cats_order = ["all", "Specific", "Unspecific", "No"]
    cat_label = {"all": "All", "Specific": "Specific",
                 "Unspecific": "Unspecific", "No": "No"}

    all_matrices, all_roi_sets = {}, []
    for cat in cats_order:
        all_matrices[cat] = {}
        for g in all_groups:
            gm = _gmask(df, g)
            subset = df[gm] if cat == "all" else df[gm & (df["fp_category"] == cat)]
            mats = []
            for _, row in subset.iterrows():
                path = row.get("corr_matrix_path", None)
                if not path or not os.path.exists(str(path)):
                    continue
                try:
                    rois, mat = _load_matrix(path)
                    mats.append((rois, mat))
                    all_roi_sets.append(set(rois))
                except Exception:
                    pass
            all_matrices[cat][g] = mats

    if all_roi_sets:
        if common_rois is not None:
            # use externally supplied list (already bilateral-paired)
            available = set.intersection(*all_roi_sets) if len(all_roi_sets) > 1 else all_roi_sets[0]
            common_rois = [r for r in common_rois if r in available]
        else:
            common_unfiltered = set.intersection(*all_roi_sets) if len(all_roi_sets) > 1 else all_roi_sets[0]
            common_rois = _bilateral_pair_filter(list(common_unfiltered))
    else:
        common_rois = common_rois or []
    if len(common_rois) < 4:
        print("  [Fig 04] Not enough common bilateral ROIs")
        return None

    results = {}
    for cat in cats_order:
        results[cat] = {}
        for g in all_groups:
            mats = all_matrices[cat][g]
            aligned = []
            for rois, mat in mats:
                idx_map = [rois.index(r) for r in common_rois if r in rois]
                if len(idx_map) == len(common_rois):
                    aligned.append(mat[np.ix_(idx_map, idx_map)])
            n_loaded = len(aligned)
            # MEDIAN across runs (robust)
            med_mat = np.nanmedian(np.stack(aligned, 0), 0) if aligned else None
            gm = _gmask(df, g)
            n_in_df = (df[gm].shape[0] if cat == "all"
                       else df[gm & (df["fp_category"] == cat)].shape[0])
            results[cat][g] = (common_rois, med_mat, n_loaded, n_in_df)

    n_rois = len(common_rois)
    n_grp = len(all_groups)
    n_rows = len(cats_order)
    short = [r.replace("L_", "L ").replace("R_", "R ")[:13] for r in common_rois]

    all_vals = []
    for cat in cats_order:
        for g in all_groups:
            _, m, _, _ = results[cat][g]
            if m is not None:
                all_vals.extend(m.ravel())
    all_vals = np.array(all_vals)
    all_vals = all_vals[np.isfinite(all_vals)]
    vmax_r = max(float(np.nanpercentile(np.abs(all_vals), 97)), 0.05) if len(all_vals) else 0.5
    norm_r = TwoSlopeNorm(vmin=-vmax_r, vcenter=0, vmax=vmax_r)

    w, h = figsize or (max(8, n_grp * 3.5), n_rows * 3.5)
    fig = plt.figure(figsize=(w, h))
    gs = gridspec.GridSpec(n_rows, n_grp + 1,
                           width_ratios=[1.0] * n_grp + [0.04],
                           wspace=0.05, hspace=0.35)
    for r_i, cat in enumerate(cats_order):
        col_c = _CAT_COLORS.get(cat, "#444") if cat != "all" else "#444"
        for c_i, g in enumerate(all_groups):
            ax = fig.add_subplot(gs[r_i, c_i])
            _, med_mat, n_loaded, n_in_df = results[cat][g]
            if med_mat is None or n_loaded == 0:
                ax.set_facecolor("#f2f2f2")
                ax.set_xticks([]); ax.set_yticks([])
                for sp in ax.spines.values():
                    sp.set_visible(False)
                ax.text(0.5, 0.5, f"0/{n_in_df}" if n_in_df > 0 else "no data",
                        ha="center", va="center", fontsize=7,
                        color="#999999", transform=ax.transAxes)
            else:
                ax.imshow(med_mat, cmap="RdBu_r", norm=norm_r,
                          aspect="equal", interpolation="nearest")
                ax.text(0.02, 0.98, f"n={n_loaded}", transform=ax.transAxes,
                        fontsize=6, va="top", ha="left", color="#333",
                        backgroundcolor="white", alpha=0.7)
                if c_i == 0:
                    step = max(1, n_rois // 30)
                    ax.set_yticks(range(0, n_rois, step))
                    ax.set_yticklabels(short[::step], fontsize=3)
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
                 cax=cbar_ax, label="Median r")
    plt.tight_layout()
    _save_or_show(fig, output_path)
    return fig

# =============================================================================
# FINGERPRINT — hierarchical test + BH-FDR
# =============================================================================

def _edge_test_species_mean(values, species):
    """One value per species (mean within species), one-sample t across species."""
    d = pd.DataFrame({"v": values, "sp": species}).dropna()
    if d.empty:
        return np.nan, np.nan, 0
    g = d.groupby("sp")["v"].mean().values
    g = g[np.isfinite(g)]
    if len(g) < 3:
        return (float(np.mean(g)) if len(g) else np.nan), np.nan, len(g)
    t, p = _stats.ttest_1samp(g, 0.0)
    return float(np.mean(g)), float(p), len(g)

def _edge_test_lmm(values, species, subjects):
    """LMM: val ~ 1 + (1|species) + (1|species:subject); test intercept vs 0."""
    try:
        import statsmodels.formula.api as smf
    except Exception:
        return np.nan, np.nan, 0
    d = pd.DataFrame({"v": values, "sp": species, "su": subjects}).dropna()
    if d["sp"].nunique() < 2 or len(d) < 4:
        return np.nan, np.nan, d["sp"].nunique()
    try:
        md = smf.mixedlm("v ~ 1", d, groups=d["sp"],
                         vc_formula={"su": "0 + C(su)"})
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = md.fit(reml=True, method="lbfgs", maxiter=200)
        return float(r.fe_params["Intercept"]), float(r.pvalues["Intercept"]), d["sp"].nunique()
    except Exception:
        return np.nan, np.nan, d["sp"].nunique()

def build_fingerprint_from_specific(df_scored, cats,
                                    bids_root_template=_DEFAULT_BIDS,
                                    atlas_name="EDNIxCSC", atlas_level=2,
                                    use_lr=True, n_top=10,
                                    species_subset=None,
                                    common_rois=None,
                                    fp_test="species_mean",
                                    fdr_alpha=0.05, verbose=True):
    df = df_scored.copy()
    df["_cat"] = cats.values if hasattr(cats, "values") else cats
    if species_subset is not None and "species" in df.columns:
        df = df[df["species"].isin(species_subset)]

    sp_mats = {}
    for sp in _sp_order(df):
        if "species" not in df.columns:
            break
        sp_df = df[(df["species"] == sp) & (df["_cat"] == "Specific")]
        if sp_df.empty:
            if verbose:
                print(f"  [fingerprint] {sp}: 0 Specific runs -- skipped")
            continue
        mats, subs, rois_ref = [], [], None
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
                subs.append(str(row.get("subject", "NA")))
            except Exception as e:
                warnings.warn(f"[fingerprint] {sp}: {e}")
        if mats:
            sp_mats[sp] = {"individual_mats": mats, "subjects": subs,
                           "rois": rois_ref, "n": len(mats)}
            if verbose:
                print(f"  [fingerprint] {sp}: {len(mats)} Specific runs")

    if len(sp_mats) < 2:
        raise ValueError(f"Only {len(sp_mats)} species with Specific runs.")
    return sp_mats, build_fc_fingerprint(sp_mats, n_top=n_top, fp_test=fp_test,
                                         fdr_alpha=fdr_alpha,
                                         common_rois=common_rois,
                                         verbose=verbose)

def build_fc_fingerprint(sp_mats, n_top=10, fp_test="species_mean",
                         fdr_alpha=0.05, common_rois=None, verbose=True):
    if len(sp_mats) < 2:
        raise ValueError(f"Need >=2 species, got {len(sp_mats)}")

    if common_rois is None:
        first = next(iter(sp_mats.values()))
        common = set(first["rois"])
        for v in sp_mats.values():
            common &= set(v["rois"])
        rois = _bilateral_pair_filter(list(common))
    else:
        # use externally supplied list (assumed bilateral-paired)
        rois = [r for r in common_rois
                if all(r in v["rois"] for v in sp_mats.values())]
    n_roi = len(rois)
    if n_roi < 4:
        raise ValueError(f"Only {n_roi} common bilateral ROIs.")
    sp_names = list(sp_mats.keys())

    # build aligned per-run arrays with species/subject labels
    run_arrs, run_sp, run_su = [], [], []
    for sp in sp_names:
        idx_map = [sp_mats[sp]["rois"].index(r) for r in rois]
        subs = sp_mats[sp].get("subjects", ["NA"] * len(sp_mats[sp]["individual_mats"]))
        for mat, su in zip(sp_mats[sp]["individual_mats"], subs):
            run_arrs.append(mat[np.ix_(idx_map, idx_map)])
            run_sp.append(sp)
            run_su.append(su)
    run_sp = np.array(run_sp)
    run_su = np.array(run_su)
    stack = np.stack(run_arrs, 0)   # (n_runs, n_roi, n_roi)

    effect = np.full((n_roi, n_roi), np.nan)
    pval = np.full((n_roi, n_roi), np.nan)

    triu = np.triu(np.ones((n_roi, n_roi), bool), k=1)
    for i, j in np.argwhere(triu):
        vals = stack[:, i, j]
        if fp_test == "lmm":
            e, p, _ = _edge_test_lmm(vals, run_sp, run_su)
            if not np.isfinite(p):    # fallback if LMM failed for this edge
                e, p, _ = _edge_test_species_mean(vals, run_sp)
        else:
            e, p, _ = _edge_test_species_mean(vals, run_sp)
        effect[i, j] = effect[j, i] = e
        pval[i, j] = pval[j, i] = p

    # BH-FDR across all tested edges (upper triangle)
    p_flat = pval[triu]
    q_flat = _bh_fdr(p_flat)
    qval = np.full((n_roi, n_roi), np.nan)
    qval[triu] = q_flat
    qval = qval + qval.T

    # connection type masks
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

    def _ctype(i, j):
        if homo_mask[i, j]:
            return "homotopic"
        if cross_mask[i, j]:
            return "heterotopic"
        return "ipsilateral"

    connections = []
    for i, j in np.argwhere(triu):
        q = qval[i, j]
        connections.append(dict(
            i=int(i), j=int(j), roi_i=rois[i], roi_j=rois[j],
            effect=float(effect[i, j]), pval=float(pval[i, j]),
            qval=float(q) if np.isfinite(q) else np.nan,
            type=_ctype(i, j),
            survived=bool(np.isfinite(q) and q < fdr_alpha),
        ))

    # display selection: survivors first, ranked by |effect|, capped per type
    survivors = [c for c in connections if c["survived"]]
    by_type = {"homotopic": [], "heterotopic": [], "ipsilateral": []}
    for c in survivors:
        by_type[c["type"]].append(c)
    display = []
    for ct in by_type:
        ranked = sorted(by_type[ct], key=lambda c: abs(c["effect"]), reverse=True)
        display.extend(ranked[:n_top])
    display.sort(key=lambda c: abs(c["effect"]), reverse=True)

    fp_idx = [(c["i"], c["j"]) for c in display]
    fp_pattern = np.array([c["effect"] for c in display])

    if verbose:
        n_surv = len(survivors)
        print(f"  [fingerprint] {len(sp_names)} species x {n_roi} ROIs | "
              f"test={fp_test} | {n_surv}/{int(triu.sum())} edges survive "
              f"FDR<{fdr_alpha} ({len(display)} shown)")
        print(f"  {'Connection':48s}  {'effect':>7}  {'q':>7}  type")
        print(f"  {'-' * 74}")
        for c in display[:8]:
            sig = ("***" if c["qval"] < .001 else "**" if c["qval"] < .01
                   else "*" if c["qval"] < .05 else "ns")
            print(f"  {c['roi_i'][:22]:22s} -- {c['roi_j'][:22]:22s}  "
                  f"{c['effect']:>+7.3f}  {c['qval']:>7.4f}{sig:3s}  {c['type']}")

    return dict(
        rois=rois, n_species=sp_names, n_roi=n_roi,
        effect_mat=effect, pval_mat=pval, qval_mat=qval,
        homo_mask=homo_mask, cross_mask=cross_mask,
        connections=connections, display=display,
        fingerprint_indices=fp_idx, fingerprint_pattern=fp_pattern,
        fp_test=fp_test, fdr_alpha=fdr_alpha,
    )

# =============================================================================
# FIGURE 05 — FINGERPRINT CONNECTIONS (right bars = FDR survivors only)
# =============================================================================

def plot_fingerprint_connections(fp, n_top_per_type=8, output_path=None):
    rois = fp["rois"]
    effect = fp["effect_mat"]
    qmat = fp["qval_mat"]
    n_roi = fp["n_roi"]
    alpha = fp["fdr_alpha"]
    short = [r.replace("L_", "L ").replace("R_", "R ")[:13] for r in rois]

    type_colors = {"homotopic": "#009E73", "heterotopic": "#56B4E9",
                   "ipsilateral": "#E69F00"}
    type_labels = {"homotopic": "Homotopic (L↔R same region)",
                   "heterotopic": "Heterotopic (L↔R cross-region)",
                   "ipsilateral": "Ipsilateral (same hemisphere)"}

    # survivors per type, ranked by effect
    surv = [c for c in fp["connections"] if c["survived"]]
    by_type = {"homotopic": [], "heterotopic": [], "ipsilateral": []}
    for c in surv:
        by_type[c["type"]].append(c)
    for ct in by_type:
        by_type[ct] = sorted(by_type[ct], key=lambda c: c["effect"],
                             reverse=True)[:n_top_per_type]

    fig_height = max(9, n_top_per_type * 0.55 + 4)
    fig = plt.figure(figsize=(16, fig_height))
    gs = gridspec.GridSpec(4, 2, height_ratios=[1, 1, 1, 0.15],
                           width_ratios=[1.0, 0.9], hspace=0.35, wspace=0.25,
                           left=0.05, right=0.96, top=0.92, bottom=0.10)

    # LEFT: effect heatmap, stars = FDR survivors
    ax_heat = fig.add_subplot(gs[0:3, 0])
    e_disp = np.where(np.isfinite(effect), effect, 0)
    nz = e_disp[e_disp != 0]
    vabs = max(float(np.nanpercentile(np.abs(nz), 97)) if nz.size else 0.1, 0.05)
    norm_e = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
    im = ax_heat.imshow(e_disp, cmap="RdBu_r", norm=norm_e,
                        aspect="equal", interpolation="nearest")
    ax_heat.set_xticks(range(n_roi)); ax_heat.set_yticks(range(n_roi))
    ax_heat.set_xticklabels(short, rotation=90, fontsize=3)
    ax_heat.set_yticklabels(short, fontsize=3)
    ax_heat.set_title(f"Edge effect ({fp['fp_test']}) — stars survive FDR<{alpha}",
                      fontsize=9, fontweight="bold")
    surv_ij = np.argwhere(np.isfinite(qmat) & (qmat < alpha))
    for i, j in surv_ij:
        if i < j:
            ax_heat.plot(j, i, "*", color=_GOLD, ms=4,
                         markeredgecolor=_GOLD_EDGE, markeredgewidth=0.5, zorder=5)
    cbar_ax_heat = fig.add_subplot(gs[3, 0])
    fig.colorbar(im, cax=cbar_ax_heat, orientation="horizontal",
                 label="Edge effect (estimate)").ax.tick_params(labelsize=7)

    # RIGHT: survivors per type
    all_e = [c["effect"] for ct in by_type for c in by_type[ct]]
    e_abs = max(float(np.nanpercentile(np.abs(all_e), 97)), 0.05) if all_e else 0.5
    norm_r = TwoSlopeNorm(vmin=-e_abs, vcenter=0, vmax=e_abs)
    cmap_r = plt.cm.RdBu_r

    for panel_i, ct in enumerate(["homotopic", "heterotopic", "ipsilateral"]):
        ax = fig.add_subplot(gs[panel_i, 1])
        conns = by_type[ct]
        color = type_colors[ct]
        if not conns:
            ax.text(0.5, 0.5, f"No {ct} edges survive FDR", ha="center",
                    va="center", transform=ax.transAxes, fontsize=9, color="#aaa")
            ax.set_title(type_labels[ct], fontsize=9, fontweight="bold", color=color)
            ax.set_xticks([]); ax.set_yticks([])
            continue
        y_pos = np.arange(len(conns))
        e_vals = [c["effect"] for c in conns]
        q_vals = [c["qval"] for c in conns]
        ax.barh(y_pos, e_vals, color=[cmap_r(norm_r(e)) for e in e_vals],
                edgecolor="white", lw=0.5, height=0.7)
        for yi, (e, q) in enumerate(zip(e_vals, q_vals)):
            sig = "***" if q < .001 else "**" if q < .01 else "*" if q < .05 else ""
            ax.text(e + (0.01 if e >= 0 else -0.01), yi,
                    f"q={q:.3f}{sig}", va="center",
                    ha="left" if e >= 0 else "right", fontsize=6, color="#333")
        ax.set_yticks(y_pos)
        ax.set_yticklabels([f"{c['roi_i'][:15]} ↔ {c['roi_j'][:15]}" for c in conns],
                           fontsize=6)
        ax.axvline(0, color="#aaa", lw=0.7)
        ax.set_xlabel("Edge effect", fontsize=7)
        ax.set_title(type_labels[ct], fontsize=8, fontweight="bold", color=color)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=6)

    cbar_ax_r = fig.add_subplot(gs[3, 1])
    sm = plt.cm.ScalarMappable(cmap=cmap_r, norm=norm_r); sm.set_array([])
    fig.colorbar(sm, cax=cbar_ax_r, orientation="horizontal",
                 label="Edge effect").ax.tick_params(labelsize=7)
    _save_or_show(fig, output_path)
    return fig

# =============================================================================
# FIGURE 06 — FINGERPRINT HEATMAP BY (species, BIDS)
# =============================================================================

def plot_fingerprint_heatmap(df_scored, fp, cats,
                             bids_root_template=_DEFAULT_BIDS,
                             atlas_name="EDNIxCSC", atlas_level=2,
                             use_lr=True, species_subset=None, output_path=None):
    display = fp["display"]
    fp_idx = fp["fingerprint_indices"]
    fp_rois = fp["rois"]
    fp_pat = fp["fingerprint_pattern"]
    n_fp = len(display)
    if n_fp == 0:
        print("  [Fig 06] no surviving connections to display")
        return None

    df = df_scored.copy()
    df["_cat"] = cats.values if hasattr(cats, "values") else cats
    if species_subset is not None and "species" in df.columns:
        df = df[df["species"].isin(species_subset)]
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
    conn_labels = [f"{c['roi_i'][:11]}--{c['roi_j'][:11]}" for c in display]
    conn_types = [c["type"] for c in display]

    group_vals = np.full((n_fp, n_grp), np.nan)
    group_frac = np.zeros(n_grp)

    for gi, (sp, bd) in enumerate(groups):
        mask = (df["species"] == sp) & (df["bids_dir"] == bd)
        sp_bd = df[mask]
        group_frac[gi] = (sp_bd["_cat"] == "Specific").sum() / max(len(sp_bd), 1)
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
                group_vals[k, gi] = float(np.median(conn_acc[k]))  # median across runs

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
    vabs = max(float(np.nanpercentile(np.abs(all_v), 97)), 0.05) if all_v.size else 0.1
    norm = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
    cmap = "RdBu_r"

    ax_fp = fig.add_subplot(gs[0, 0])
    ax_fp.imshow(fp_pat.reshape(-1, 1), cmap=cmap, norm=norm,
                 aspect="auto", interpolation="nearest")
    ax_fp.set_yticks(range(n_fp))
    ax_fp.set_yticklabels(conn_labels, fontsize=max(3, 6 - n_fp // 10))
    ax_fp.set_xticks([0]); ax_fp.set_xticklabels(["Expected"], fontsize=7,
                                                 rotation=30, ha="right")
    ax_fp.set_title("Pattern", fontsize=8, fontweight="bold")
    ct_col = {"homotopic": "#009E73", "heterotopic": "#56B4E9", "ipsilateral": "#E69F00"}
    for k, ct in enumerate(conn_types):
        ax_fp.axhline(k + 0.5, color=ct_col.get(ct, "#aaa"), lw=0.4, alpha=0.5)

    sp_col_ranges = {}
    for gi, (sp, bd) in enumerate(groups):
        ax = fig.add_subplot(gs[0, gi + 1])
        ax.imshow(group_vals[:, gi].reshape(-1, 1), cmap=cmap, norm=norm,
                  aspect="auto", interpolation="nearest")
        ax.set_yticks([]); ax.set_xticks([0])
        bd_lbl = bd.replace("BIDS_", "").replace("ds004513-download", "Hum")[:10]
        ax.set_xticklabels([bd_lbl], fontsize=max(4, 7 - n_grp // 4),
                           rotation=45, ha="right", color=sp_gc.get(sp, "#333"))
        sp_col_ranges.setdefault(sp, []).append(gi)

    for sp, col_list in sp_col_ranges.items():
        lo, hi = min(col_list), max(col_list)
        x_lo = fig.axes[lo + 1].get_position().x0
        x_hi = fig.axes[hi + 1].get_position().x1
        y_top = fig.axes[lo + 1].get_position().y1
        fig.text((x_lo + x_hi) / 2, y_top + 0.015, sp, ha="center", va="bottom",
                 fontsize=9, fontweight="bold", color=sp_gc.get(sp, "#333"))

    cbar_ax = fig.add_subplot(gs[0, -1])
    fig.colorbar(_mpl_cm.ScalarMappable(norm=norm, cmap=cmap),
                 cax=cbar_ax, label="Median r")
    cbar_ax.yaxis.set_tick_params(labelsize=6)

    ax_bar = fig.add_subplot(gs[1, 1:n_grp + 1])
    ax_bar.bar(np.arange(n_grp), group_frac * 100,
               color=[sp_gc.get(sp, "#888") for sp, _ in groups],
               edgecolor="white", lw=0.4)
    for xi, v in enumerate(group_frac * 100):
        ax_bar.text(xi, v + 0.5, f"{v:.0f}%", ha="center", va="bottom",
                    fontsize=max(4, 7 - n_grp // 4), color="#333")
    ax_bar.set_xticks([]); ax_bar.set_ylabel("% Specific runs", fontsize=7)
    ax_bar.set_ylim(0, 110)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)

    fig.legend(handles=[mpatches.Patch(color=ct_col[k], label=k) for k in ct_col],
               loc="lower left", bbox_to_anchor=(0.0, 0.0),
               ncol=3, fontsize=7, frameon=False)
    _save_or_show(fig, output_path)
    return fig

# =============================================================================
# FIGURE 07 — HUMAN THRESHOLD SWEEP
# =============================================================================

def _human_sweep_figure(mats, intra, inter, sweep_intra, sweep_delta,
                        thresh_intra, thresh_delta, sweep_label, species,
                        output_path=None):
    """
    Render one sweep figure for human runs.

    sweep_intra : list/array of intra thresholds (one per column)
    sweep_delta : list/array of delta thresholds (one per column)
    sweep_label : axis label describing what is varied
    """
    n_cols = len(sweep_intra)
    assert n_cols == len(sweep_delta)
    cats_order = ["all", "Specific", "Unspecific", "No"]
    n_rows = len(cats_order)

    # color norm from all matrices
    allv = np.concatenate([m.ravel() for m in mats])
    allv = allv[np.isfinite(allv)]
    vmax = max(float(np.nanpercentile(np.abs(allv), 97)), 0.05)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    fig = plt.figure(figsize=(max(8, n_cols * 3.2), n_rows * 3.2 + 1.4))
    # 2 rows of axes: matrix grid (top) + proportions strip (bottom)
    gs = gridspec.GridSpec(n_rows + 1, n_cols + 1,
                           height_ratios=[1.0] * n_rows + [0.6],
                           width_ratios=[1.0] * n_cols + [0.04],
                           wspace=0.05, hspace=0.30)

    # precompute per-column classification + proportions
    col_props = []
    col_runcats = []
    for c_i in range(n_cols):
        t_in = sweep_intra[c_i]
        t_de = sweep_delta[c_i]
        run_cats = np.array([_classify_one(intra[k], inter[k], t_in, t_de)
                             for k in range(len(mats))])
        col_runcats.append(run_cats)
        n_total = len(run_cats)
        col_props.append({c: 100.0 * (run_cats == c).sum() / max(n_total, 1)
                          for c in CATS})

    # Matrix grid
    for c_i in range(n_cols):
        t_in = sweep_intra[c_i]
        t_de = sweep_delta[c_i]
        run_cats = col_runcats[c_i]
        is_ref = (abs(t_in - thresh_intra) < 1e-9 and
                  abs(t_de - thresh_delta) < 1e-9)

        for r_i, cat in enumerate(cats_order):
            ax = fig.add_subplot(gs[r_i, c_i])
            sel = (np.ones(len(mats), bool) if cat == "all" else run_cats == cat)
            n_sel = int(sel.sum())
            if n_sel == 0:
                ax.set_facecolor("#f2f2f2")
                ax.set_xticks([]); ax.set_yticks([])
                for s in ax.spines.values():
                    s.set_visible(False)
                ax.text(0.5, 0.5, "0", ha="center", va="center",
                        fontsize=8, color="#999", transform=ax.transAxes)
            else:
                med = np.nanmedian(np.stack([m for m, s in zip(mats, sel) if s], 0), 0)
                ax.imshow(med, cmap="RdBu_r", norm=norm, aspect="equal",
                          interpolation="nearest")
                ax.text(0.02, 0.98, f"n={n_sel}", transform=ax.transAxes,
                        fontsize=6, va="top", ha="left", color="#333",
                        backgroundcolor="white", alpha=0.7)
                ax.set_xticks([]); ax.set_yticks([])
            if c_i == 0:
                ax.set_ylabel("All" if cat == "all" else cat, fontsize=9,
                              color=_CAT_COLORS.get(cat, "#444"),
                              fontweight="bold")
            if r_i == 0:
                star = "  ★" if is_ref else ""
                ax.set_title(f"intra>{t_in:.3f}\nΔ>{t_de:.3f}{star}",
                             fontsize=7, fontweight="bold")

    # Colorbar
    cbar_ax = fig.add_subplot(gs[:n_rows, -1])
    fig.colorbar(_mpl_cm.ScalarMappable(norm=norm, cmap="RdBu_r"),
                 cax=cbar_ax, label="Median r")

    # Bottom strip: stacked-bar proportions
    ax_bar = fig.add_subplot(gs[n_rows, :n_cols])
    x = np.arange(n_cols)
    bottoms = np.zeros(n_cols)
    for cat in CATS:
        vals = np.array([col_props[c_i][cat] for c_i in range(n_cols)])
        ax_bar.bar(x, vals, 0.7, bottom=bottoms, color=_CAT_COLORS[cat],
                   label=cat, edgecolor="white", lw=0.5)
        for xi, (v, b) in enumerate(zip(vals, bottoms)):
            if v > 5:
                ax_bar.text(xi, b + v / 2, f"{v:.0f}%", ha="center", va="center",
                            fontsize=7, color="white", fontweight="bold")
        bottoms += vals
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels([f"intra>{sweep_intra[c]:.3f}\nΔ>{sweep_delta[c]:.3f}"
                            for c in range(n_cols)], fontsize=7)
    ax_bar.set_ylabel("% runs", fontsize=8)
    ax_bar.set_xlabel(sweep_label, fontsize=8)
    ax_bar.set_ylim(0, 105)
    ax_bar.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30),
                  ncol=len(CATS), fontsize=7, frameon=False)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)

    fig.suptitle(f"{species}: {sweep_label}",
                 fontsize=12, fontweight="bold")
    # gridspec layout already explicit; tight_layout doesn't help with colorbar mix
    _save_or_show(fig, output_path)
    return fig


def _load_human_runs(df_scored, species="Human"):
    """Load + align human matrices once. Returns (df_filtered, mats, intra, inter)."""
    inter_col = _inter_col(df_scored)
    df = df_scored[df_scored["species"] == species].copy() \
        if "species" in df_scored.columns else df_scored.copy()
    if df.empty:
        return None, [], np.array([]), np.array([])
    mats, idx_keep, rois_ref = [], [], None
    for idx, row in df.iterrows():
        path = row.get("corr_matrix_path", None)
        if not path or not os.path.exists(str(path)):
            continue
        try:
            rois, mat = _load_matrix(path)
            if rois_ref is None:
                rois_ref = rois
            if len(rois) != len(rois_ref):
                continue
            mats.append(mat); idx_keep.append(idx)
        except Exception:
            pass
    if not mats:
        return None, [], np.array([]), np.array([])
    df = df.loc[idx_keep]
    intra = df["primary_intra"].values.astype(float)
    inter = df[inter_col].values.astype(float)
    return df, mats, intra, inter


def plot_human_threshold_sweep(df_scored, thresh_intra, thresh_delta,
                               species="Human",
                               factors=(0.8, 0.9, 1.0, 1.1, 1.2),
                               output_path=None):
    """
    Convenience wrapper: render the COMBINED sweep (both thresholds scaled
    together). For the three-panel decomposition use
    plot_human_threshold_sweeps() (plural).
    """
    df, mats, intra, inter = _load_human_runs(df_scored, species)
    if not mats:
        print(f"  [Fig 07 combined] no {species} runs")
        return None
    sweep_intra = [thresh_intra * f for f in factors]
    sweep_delta = [thresh_delta * f for f in factors]
    return _human_sweep_figure(
        mats, intra, inter, sweep_intra, sweep_delta,
        thresh_intra, thresh_delta,
        sweep_label="both thresholds × factor (intra and Δ scale together)",
        species=species, output_path=output_path)


def plot_human_threshold_sweeps(df_scored, thresh_intra, thresh_delta,
                                species="Human",
                                factors=(0.8, 0.9, 1.0, 1.1, 1.2),
                                output_paths=None):
    """
    Three separate sweep figures:
      (a) vary intra only      (Δ fixed at cross-species value)
      (b) vary Δ only          (intra fixed at cross-species value)
      (c) vary both            (intra and Δ scaled by the same factor)

    output_paths : dict with keys 'intra', 'delta', 'both' (any may be None)
                   or None to display.
    """
    df, mats, intra, inter = _load_human_runs(df_scored, species)
    if not mats:
        print(f"  [Fig 07] no {species} runs")
        return None

    op = output_paths or {}
    figs = {}

    # (a) vary intra only
    figs["intra"] = _human_sweep_figure(
        mats, intra, inter,
        sweep_intra=[thresh_intra * f for f in factors],
        sweep_delta=[thresh_delta] * len(factors),
        thresh_intra=thresh_intra, thresh_delta=thresh_delta,
        sweep_label="varying intra threshold only (Δ fixed)",
        species=species, output_path=op.get("intra"))

    # (b) vary delta only
    figs["delta"] = _human_sweep_figure(
        mats, intra, inter,
        sweep_intra=[thresh_intra] * len(factors),
        sweep_delta=[thresh_delta * f for f in factors],
        thresh_intra=thresh_intra, thresh_delta=thresh_delta,
        sweep_label="varying Δ threshold only (intra fixed)",
        species=species, output_path=op.get("delta"))

    # (c) vary both
    figs["both"] = _human_sweep_figure(
        mats, intra, inter,
        sweep_intra=[thresh_intra * f for f in factors],
        sweep_delta=[thresh_delta * f for f in factors],
        thresh_intra=thresh_intra, thresh_delta=thresh_delta,
        sweep_label="varying both thresholds together",
        species=species, output_path=op.get("both"))

    return figs

# =============================================================================
# FIGURE 08 — HUMAN SURROGATE / PERTURBATION SENSITIVITY
# =============================================================================

def _make_surrogate(mat, kind, rng, shift=0.10, alpha=0.4, sigma=0.10):
    """
    Parametric perturbations of a real matrix (NOT formal null models):
      shift  : add a constant offset to off-diagonal edges
      smooth : shrink edges toward their median (focal point), variance down
      noise  : add Gaussian noise per edge
    Off-diagonal is clipped to [-1, 1]; symmetry and diagonal are preserved.
    """
    n = mat.shape[0]
    off = ~np.eye(n, dtype=bool)
    out = mat.copy()
    if kind == "shift":
        out[off] = mat[off] + shift
    elif kind == "smooth":
        focal = np.nanmedian(mat[off])
        out[off] = focal + alpha * (mat[off] - focal)
    elif kind == "noise":
        noise = rng.normal(0, sigma, size=(n, n))
        noise = (noise + noise.T) / 2
        out[off] = mat[off] + noise[off]
    out = np.clip(out, -1, 1)
    out = (out + out.T) / 2
    np.fill_diagonal(out, np.diag(mat))
    return out

def plot_human_surrogate_sensitivity(df_scored, cats, thresh_intra, thresh_delta,
                                     species="Human", seed=42, output_path=None):
    """
    For each surrogate type (Original, Shift, Smooth, Noise):
      - apply the perturbation to EVERY human Specific run individually
      - score and classify each perturbed matrix
      - show (top) the median perturbed matrix + its classification
      - show (bottom) the stacked-bar proportion of categories across runs

    Note: shift/smooth/noise are PARAMETRIC PERTURBATIONS, not formal null
    surrogates (they don't preserve BOLD spectral / autocorr structure).
    """
    df = df_scored.copy()
    df["_cat"] = cats.values if hasattr(cats, "values") else cats
    sp_df = df[(df["species"] == species) & (df["_cat"] == "Specific")] \
        if "species" in df.columns else df
    if sp_df.empty:
        sp_df = df[df["species"] == species] if "species" in df.columns else df
    mats, rois_ref = [], None
    for _, row in sp_df.iterrows():
        path = row.get("corr_matrix_path", None)
        if not path or not os.path.exists(str(path)):
            continue
        try:
            rois, mat = _load_matrix(path)
            if rois_ref is None:
                rois_ref = rois
            if len(rois) == len(rois_ref):
                mats.append(mat)
        except Exception:
            pass
    if not mats:
        print(f"  [Fig 08] no {species} matrices")
        return None

    rng = np.random.default_rng(seed)
    surrogate_kinds = [
        ("Original",       None,      {}),
        ("Shift (+0.10)",  "shift",   dict(shift=0.10)),
        ("Smooth (α=0.4)", "smooth",  dict(alpha=0.4)),
        ("Noise (σ=0.10)", "noise",   dict(sigma=0.10)),
    ]

    # Apply each perturbation to every run, classify each, summarise.
    panel_data = []
    for label, kind, params in surrogate_kinds:
        perturbed = ([m for m in mats] if kind is None
                     else [_make_surrogate(m, kind, rng, **params) for m in mats])
        # classify each perturbed run
        run_cats = []
        for m in perturbed:
            ps = compute_primary_score_from_matrix(m, rois_ref)
            run_cats.append(_classify_one(ps["mean_intra"], ps["mean_inter_hetero"],
                                          thresh_intra, thresh_delta))
        # median across runs (for display)
        med = np.nanmedian(np.stack(perturbed, 0), 0)
        med_ps = compute_primary_score_from_matrix(med, rois_ref)
        med_cat = _classify_one(med_ps["mean_intra"], med_ps["mean_inter_hetero"],
                                thresh_intra, thresh_delta)
        # proportions
        n_total = len(run_cats)
        props = {c: 100.0 * sum(1 for x in run_cats if x == c) / max(n_total, 1)
                 for c in CATS}
        panel_data.append(dict(label=label, median_mat=med, median_cat=med_cat,
                               median_intra=med_ps["mean_intra"],
                               median_delta=med_ps["delta"],
                               props=props, n=n_total))

    # color norm from all displayed matrices
    allv = np.concatenate([d["median_mat"].ravel() for d in panel_data])
    allv = allv[np.isfinite(allv)]
    vmax = max(float(np.nanpercentile(np.abs(allv), 97)), 0.05)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    n_pan = len(panel_data)
    fig = plt.figure(figsize=(n_pan * 3.4, 7.5))
    gs = gridspec.GridSpec(2, n_pan + 1,
                           height_ratios=[1.6, 1.0],
                           width_ratios=[1.0] * n_pan + [0.04],
                           wspace=0.10, hspace=0.40)

    # top row: median perturbed matrix
    for c_i, d in enumerate(panel_data):
        ax = fig.add_subplot(gs[0, c_i])
        ax.imshow(d["median_mat"], cmap="RdBu_r", norm=norm,
                  aspect="equal", interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(d["label"], fontsize=9, fontweight="bold")
        ax.text(0.5, -0.06,
                f"median: intra={d['median_intra']:.3f}  Δ={d['median_delta']:.3f}\n"
                f"→ {d['median_cat']}",
                transform=ax.transAxes, ha="center", va="top", fontsize=7.5,
                color=_CAT_COLORS.get(d["median_cat"], "#444"), fontweight="bold")
    cbar_ax = fig.add_subplot(gs[0, -1])
    fig.colorbar(_mpl_cm.ScalarMappable(norm=norm, cmap="RdBu_r"),
                 cax=cbar_ax, label="r")

    # bottom row: stacked proportion bars across all human runs
    ax_bar = fig.add_subplot(gs[1, :n_pan])
    x = np.arange(n_pan)
    bottoms = np.zeros(n_pan)
    for cat in CATS:
        vals = np.array([d["props"][cat] for d in panel_data])
        ax_bar.bar(x, vals, 0.65, bottom=bottoms, color=_CAT_COLORS[cat],
                   label=cat, edgecolor="white", lw=0.5)
        for xi, (v, b) in enumerate(zip(vals, bottoms)):
            if v > 5:
                ax_bar.text(xi, b + v / 2, f"{v:.0f}%", ha="center", va="center",
                            fontsize=8, color="white", fontweight="bold")
        bottoms += vals
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels([f"{d['label']}\n(n={d['n']} runs)" for d in panel_data],
                           fontsize=8)
    ax_bar.set_ylabel("% of human runs", fontsize=9)
    ax_bar.set_ylim(0, 105)
    ax_bar.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18),
                  ncol=len(CATS), fontsize=8, frameon=False)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)

    fig.suptitle(f"{species}: parametric perturbation sensitivity "
                 f"(thresholds intra>{thresh_intra:.3f}, Δ>{thresh_delta:.3f})",
                 fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0, 0.02, 1, 0.94])
    _save_or_show(fig, output_path)
    return fig

# =============================================================================
# PIPELINE (single level + subset)
# =============================================================================

def run_pipeline(df_qc, fit_kind="correlation",
                 bids_root_template=_DEFAULT_BIDS,
                 thresh_intra=None, thresh_delta=None, stringency=0.0,
                 n_top=10, atlas_name="EDNIxCSC", atlas_level=2, use_lr=True,
                 subsets=("all", "primates", "primates_rodents"),
                 fp_test="species_mean", fdr_alpha=0.05,
                 make_human_figs=True, fig_dir=None, verbose=True):
    """
    Full pipeline for ONE atlas level. Within a single call:

      Steps 1-3 (score, threshold, classify) run on ALL species so the
      thresholds are a stable cross-species reference.

      Figures 04 / 05 / 06 are then re-rendered for EACH taxonomic subset
      (the subset's common bilateral ROI set is recomputed each time, so
      restricting to primates yields a richer matrix).

      Figures 07 / 08 (human-specific) are produced once.

    Output tree: fig_dir/lvl{L}/{subset}/NN_*.png
    """
    print("\n" + "=" * 60 + f"\n  EDNiX FC Pipeline — lvl{atlas_level}\n" + "=" * 60)

    def _p(name, subset_name):
        if fig_dir is None:
            return None
        d = os.path.join(fig_dir, f"lvl{atlas_level}", subset_name)
        os.makedirs(d, exist_ok=True)
        return os.path.join(d, name)

    # ── 1-3: score + threshold + classify (all species) ──────────────────────
    print("\n[1] Scoring all runs...")
    df_scored = score_all_runs(df_qc, fit_kind,
                               bids_root_template=bids_root_template,
                               atlas_name=atlas_name, atlas_level=atlas_level,
                               use_lr=use_lr, verbose=verbose)

    print("\n[2] Thresholds (pooled cross-species)...")
    thresh = suggest_thresholds(df_scored, thresh_intra=thresh_intra,
                                thresh_delta=thresh_delta, stringency=stringency,
                                verbose=verbose)
    t_in, t_de = thresh["thresh_intra"], thresh["thresh_delta"]

    print("\n[3] Classifying...")
    cats = classify_runs(df_scored, t_in, t_de, verbose=verbose)
    df_scored["fp_category"] = cats.values

    # ── threshold + class-overview figures (all species, single render) ──────
    print("\n[4] Fig 01 threshold selection...")
    plot_threshold_selection(df_scored, t_in, t_de,
                             output_path=_p("01_threshold_selection_distributions.png", "all"))
    print("\n[5] Fig 02 classification model...")
    plot_classification_model(df_scored, t_in, t_de, cats,
                              output_path=_p("02_classification_model_per_species.png", "all"))
    print("\n[6] Fig 03 category proportions...")
    plot_category_proportions(df_scored, cats, group_by="species",
                              output_path=_p("03_category_proportions_by_species.png", "all"))
    plot_category_proportions(df_scored, cats, group_by="bids_dir",
                              output_path=_p("03_category_proportions_by_bids.png", "all"))

    # ── per-subset figures 04 / 05 / 06 ──────────────────────────────────────
    fingerprints = {}
    for subset_name in subsets:
        species_subset = _SUBSETS.get(subset_name, None)
        if species_subset is not None and "species" in df_scored.columns:
            n_runs = df_scored["species"].isin(species_subset).sum()
            if n_runs == 0:
                print(f"\n  [skip subset {subset_name}] no runs found")
                continue

        print(f"\n  ─── subset: {subset_name} ───")

        # SAME common ROI list used by fig 04 and the fingerprint figures
        common_rois = compute_subset_common_rois(
            df_scored, species_subset=species_subset, verbose=verbose)

        print(f"  [7-{subset_name}] Fig 04 FC matrices by category...")
        plot_fc_matrices_by_category(
            df_scored, cats, bids_root_template=bids_root_template,
            atlas_name=atlas_name, atlas_level=atlas_level, use_lr=use_lr,
            group_by="species", species_subset=species_subset,
            common_rois=common_rois,
            output_path=_p("04_fc_matrices_by_category_species.png", subset_name))
        plot_fc_matrices_by_category(
            df_scored, cats, bids_root_template=bids_root_template,
            atlas_name=atlas_name, atlas_level=atlas_level, use_lr=use_lr,
            group_by="bids_dir", species_subset=species_subset,
            common_rois=common_rois,
            output_path=_p("04_fc_matrices_by_category_bids.png", subset_name))

        print(f"  [8-{subset_name}] Fingerprint...")
        fp, sp_mats = None, {}
        try:
            sp_mats, fp = build_fingerprint_from_specific(
                df_scored, cats, bids_root_template=bids_root_template,
                atlas_name=atlas_name, atlas_level=atlas_level, use_lr=use_lr,
                n_top=n_top, species_subset=species_subset,
                common_rois=common_rois,
                fp_test=fp_test, fdr_alpha=fdr_alpha, verbose=verbose)
            print(f"  [9-{subset_name}] Fig 05 fingerprint connections...")
            plot_fingerprint_connections(
                fp, output_path=_p("05_fingerprint_connections.png", subset_name))
            print(f"  [10-{subset_name}] Fig 06 fingerprint heatmap...")
            plot_fingerprint_heatmap(
                df_scored, fp, cats, bids_root_template=bids_root_template,
                atlas_name=atlas_name, atlas_level=atlas_level, use_lr=use_lr,
                species_subset=species_subset,
                output_path=_p("06_fingerprint_heatmap_by_bids.png", subset_name))
        except ValueError as e:
            print(f"  [skip fingerprint {subset_name}] {e}")
        fingerprints[subset_name] = dict(sp_mats=sp_mats, fp=fp,
                                          common_rois=common_rois)

    # ── human-specific figures (rendered once, under subset=all) ─────────────
    if (make_human_figs and "species" in df_scored.columns
            and (df_scored["species"] == "Human").any()):
        print("\n[11] Fig 07 human threshold sweeps (intra / Δ / both)...")
        plot_human_threshold_sweeps(
            df_scored, t_in, t_de,
            output_paths={
                "intra": _p("07a_human_sweep_intra_only.png", "all"),
                "delta": _p("07b_human_sweep_delta_only.png", "all"),
                "both":  _p("07c_human_sweep_both.png", "all"),
            })
        print("\n[12] Fig 08 human surrogate sensitivity...")
        plot_human_surrogate_sensitivity(df_scored, cats, t_in, t_de,
                                         output_path=_p("08_human_surrogate_sensitivity.png", "all"))

    print("\n" + "=" * 60 + f"\n  lvl{atlas_level} complete.")
    if fig_dir:
        print(f"  Figures: {os.path.join(fig_dir, f'lvl{atlas_level}')}")
    print("=" * 60 + "\n")
    return dict(df_scored=df_scored, thresh=thresh, cats=cats,
                fingerprints=fingerprints)

# =============================================================================
# MULTI-LEVEL / MULTI-SUBSET DRIVER
# =============================================================================

def run_pipeline_multilevel(df_qc, fit_kind="correlation",
                            bids_root_template=_DEFAULT_BIDS,
                            atlas_levels=(2, 3, 4),
                            subsets=("all", "primates", "primates_rodents"),
                            thresh_intra=None, thresh_delta=None, stringency=0.0,
                            n_top=10, atlas_name="EDNIxCSC", use_lr=True,
                            fp_test="species_mean", fdr_alpha=0.05,
                            fig_dir=None, verbose=True):
    """
    Loop the full pipeline over atlas levels. Each level call internally
    renders all subset variants (all / primates / primates_rodents).

    Output tree: fig_dir/lvl{L}/{subset}/NN_*.png
    """
    results = {}
    for lvl in atlas_levels:
        try:
            res = run_pipeline(
                df_qc, fit_kind=fit_kind, bids_root_template=bids_root_template,
                thresh_intra=thresh_intra, thresh_delta=thresh_delta,
                stringency=stringency, n_top=n_top, atlas_name=atlas_name,
                atlas_level=lvl, use_lr=use_lr, subsets=subsets,
                fp_test=fp_test, fdr_alpha=fdr_alpha,
                fig_dir=fig_dir, verbose=verbose)
            results[lvl] = res
        except Exception as e:
            print(f"  [skip lvl{lvl}] {e}")
    return results

# =============================================================================
# CLI
# =============================================================================

def main():
    p = argparse.ArgumentParser(description="EDNiX FC pipeline")
    p.add_argument("--csv", required=True)
    p.add_argument("--bids-root", default=_DEFAULT_BIDS)
    p.add_argument("--thresh-intra", type=float, default=None)
    p.add_argument("--thresh-delta", type=float, default=None)
    p.add_argument("--stringency", type=float, default=0.0)
    p.add_argument("--n-top", type=int, default=10)
    p.add_argument("--levels", type=int, nargs="+", default=[2, 3, 4])
    p.add_argument("--subsets", nargs="+",
                   default=["all", "primates", "primates_rodents"])
    p.add_argument("--fp-test", default="species_mean",
                   choices=["species_mean", "lmm"])
    p.add_argument("--fdr-alpha", type=float, default=0.05)
    p.add_argument("--output", default="figures/fc_explorer")
    args = p.parse_args()
    df_qc = pd.read_csv(args.csv)
    run_pipeline_multilevel(
        df_qc, bids_root_template=args.bids_root,
        atlas_levels=tuple(args.levels), subsets=tuple(args.subsets),
        thresh_intra=args.thresh_intra, thresh_delta=args.thresh_delta,
        stringency=args.stringency, n_top=args.n_top,
        fp_test=args.fp_test, fdr_alpha=args.fdr_alpha, fig_dir=args.output)

if __name__ == "__main__":
    main()
