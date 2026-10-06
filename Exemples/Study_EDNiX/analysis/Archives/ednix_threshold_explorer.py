"""
EDNiX — Primary Network FC Threshold Explorer & QC Classifier
==============================================================

REVISED: Uses top‑N homotopic pairs for primary_inter (matching sp_specificity_index_primary).
T1 = threshold on species_intra (all homotopic pairs, optionally common ROIs).
T2 = threshold on primary_inter (cross‑network correlation).
Classification: Specific iff species_intra > T1 AND primary_inter < T2.
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
# CONSTANTS (unchanged)
# =============================================================================

_CAT_COLORS = {
    "Specific": "#009E73",
    "Unspecific": "#F0E442",
    "No": "#56B4E9",
    "Spurious": "#D55E00",
}

_PHYLO = ["Mouse", "Rat", "Mouselemur", "Bat", "Marmoset",
          "Macaque", "Human", "Dog", "Cat", "Pig"]

_SPECIES_COLORS = {
    "Mouse": "#1A5276", "Rat": "#2E86C1", "Bat": "#00A896",
    "Mouselemur": "#58D68D", "Marmoset": "#F4D03F",
    "Macaque": "#F5B041", "Human": "#D35400",
    "Dog": "#8E44AD", "Cat": "#C0392B", "Pig": "#6E2C00",
}

_PALETTE = list(_SPECIES_COLORS.values())
_GOLD = "#FFD700"
_GOLD_EDGE = "#B8860B"
CATS = ["Specific", "Unspecific", "No", "Spurious"]
_DEFAULT_BIDS = "/scratch2/EDNiX/{species}/{bids_dir}"

_PRIMARY_PATTERNS = ['Somatosensory', 'Auditory_co', 'Visual_stri']
_PRIMATES = ["Mouselemur", "Marmoset", "Macaque", "Human"]
_RODENTS = ["Mouse", "Rat"]
_SUBSETS = {
    "all": None,
    "primates": _PRIMATES,
    "primates_rodents": _PRIMATES + _RODENTS,
}

# =============================================================================
# HELPERS (unchanged)
# =============================================================================

def _sp_order(df):
    vals = df["species"].unique().tolist() if "species" in df.columns else []
    return [s for s in _PHYLO if s in vals] + [s for s in vals if s not in _PHYLO]

def _gcolors(groups):
    return {g: _SPECIES_COLORS.get(g, _PALETTE[i % len(_PALETTE)])
            for i, g in enumerate(groups)}

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
    except ImportError:
        load_corr_matrix = None
    if load_corr_matrix is not None:
        try:
            return load_corr_matrix(path)
        except Exception:
            pass
    df = pd.read_csv(path, index_col=0)
    df.columns = df.columns.astype(str).str.strip()
    df.index = df.index.astype(str).str.strip()
    return list(df.index), df.values.astype(float)

def _bilateral_pair_filter(rois):
    l_bases = {n[2:] for n in rois if n.startswith(("L_", "l_"))}
    r_bases = {n[2:] for n in rois if n.startswith(("R_", "r_"))}
    common = sorted(l_bases & r_bases)
    return [f"L_{b}" for b in common] + [f"R_{b}" for b in common]

def compute_subset_common_rois(df_scored, species_subset=None,
                                report_csv_path=None, verbose=True):
    # unchanged
    df = df_scored.copy()
    if species_subset is not None and "species" in df.columns:
        df = df[df["species"].isin(species_subset)]

    species_union, species_inter, species_perfile = {}, {}, {}
    n_loaded_per_sp = {}

    for _, row in df.iterrows():
        path = row.get("corr_matrix_path", None)
        sp = str(row.get("species", "?"))
        if not path or not os.path.exists(str(path)):
            continue
        try:
            rois, _ = _load_matrix(path)
            roi_set = set(rois)
            species_union.setdefault(sp, set()).update(roi_set)
            if sp in species_inter:
                species_inter[sp] &= roi_set
            else:
                species_inter[sp] = set(roi_set)
            species_perfile.setdefault(sp, []).append({
                "path": path, "rois": roi_set,
                "bids_dir": str(row.get("bids_dir", "?")),
                "subject": str(row.get("subject", "?")),
                "session": str(row.get("session", "?")),
                "run": str(row.get("run", "?")),
            })
            n_loaded_per_sp[sp] = n_loaded_per_sp.get(sp, 0) + 1
        except Exception:
            pass

    if not species_inter:
        return []

    raw_common = (set.intersection(*species_inter.values())
                  if len(species_inter) > 1
                  else set(next(iter(species_inter.values()))))
    paired = _bilateral_pair_filter(list(raw_common))

    if verbose:
        print(f"  [common ROIs] FINAL: {len(paired)} bilateral pairs")
    return paired

def _bh_fdr(pvals):
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
    # unchanged
    lr = "LR" if use_lr else ""
    fk_re = re.escape(fit_kind)
    name_re = re.compile(
        rf"^{re.escape(atlas_name + lr)}_{atlas_level}_run_(\d+)"
        rf"(?:_flattened)?_{fk_re}_matrix\.csv$"
    )
    excl = ("check_fit", "pval", "tstat", "flattened")
    sub_dir = os.path.join(bids_root, f"sub-{subject}")

    ses_dirs = []
    if session:
        ses_dirs.append(os.path.join(sub_dir, f"ses-{session}"))
    ses_dirs.append(sub_dir)

    all_files = []
    for ses_d in ses_dirs:
        d = os.path.join(ses_d, "func", "acpc-func", "Stats", "Correl_matrix",
                         atlas_name + lr, fit_kind)
        if not os.path.isdir(d):
            continue
        for fn in sorted(os.listdir(d)):
            if "flattened" in fn:
                continue
            if any(x in fn for x in excl):
                continue
            if name_re.match(fn):
                all_files.append(os.path.join(d, fn))

    unique_files = sorted(set(all_files))
    return unique_files

def _get_all_run_paths_for_row(fit_kind, row, bids_root_template,
                               atlas_name="EDNIxCSC", atlas_level=2, use_lr=True):
    # unchanged
    sp = str(row.get("species", ""))
    bd = str(row.get("bids_dir", ""))
    sub = str(row.get("subject", ""))
    ses = str(row.get("session", ""))
    if not bids_root_template:
        return []
    root = bids_root_template.format(species=sp, bids_dir=bd)
    paths = _find_all_run_paths(root, sub, ses, fit_kind, atlas_name, atlas_level, use_lr)

    unique_paths = {}
    for p in paths:
        norm_p = os.path.realpath(p)
        if norm_p not in unique_paths:
            unique_paths[norm_p] = p
    return [(unique_paths[p], sp, bd, sub, ses) for p in sorted(unique_paths.keys())]

def _is_homotopic(a, b):
    def _base(n):
        return re.sub(r"^[LlRr]_", "", str(n)).lower().strip()
    return _base(a) == _base(b)

# -----------------------------------------------------------------------------
# NEW CLASSIFICATION RULE (using two thresholds)
# -----------------------------------------------------------------------------
def _classify_one(intra_primary, inter_primary, intra_species,
                  t_intra, t_inter):
    """
    Classify a single run using separate thresholds on intra and inter.

    Rules:
    - Spurious: primary_inter > primary_intra (cross > homotopic)
    - No: species_intra <= t_intra (homotopic connectivity too low)
    - Specific: species_intra > t_intra AND primary_inter < t_inter
    - Unspecific: all other cases
    """
    if not np.isfinite(intra_primary) or not np.isfinite(inter_primary) or not np.isfinite(intra_species):
        return "No"

    # Spurious: cross > homotopic
    if inter_primary > intra_primary:
        return "Spurious"

    # No: intra (all homotopic) below threshold
    if intra_species <= t_intra:
        return "No"

    # Specific: high intra AND low inter
    if intra_species > t_intra and inter_primary < t_inter:
        return "Specific"

    # Unspecific: everything else
    return "Unspecific"

# =============================================================================
# NORMALIZATION FUNCTIONS (unchanged)
# =============================================================================

def normalize_specific_matrices(all_specific_mats, target_median=None, method="translate"):
    # unchanged
    if not all_specific_mats:
        return [], 0

    all_values = []
    for m in all_specific_mats:
        triu_idx = np.triu_indices_from(m, k=1)
        vals = m[triu_idx]
        vals = vals[np.isfinite(vals)]
        all_values.extend(vals)
    all_values = np.array(all_values)

    if len(all_values) == 0:
        return all_specific_mats, 0

    ref_median = np.nanmedian(all_values)

    if target_median is None:
        target_median = ref_median
        print(f"  [normalize] Aligning all matrices to global reference median: {ref_median:.3f}")
    else:
        print(f"  [normalize] Reference median: {ref_median:.3f} → Target: {target_median:.3f}")

    normalized_mats = []
    for mat in all_specific_mats:
        mat_norm = mat.copy()
        triu_idx = np.triu_indices_from(mat, k=1)
        values = mat[triu_idx]
        values = values[np.isfinite(values)]

        if len(values) == 0:
            normalized_mats.append(mat)
            continue

        current_median = np.nanmedian(values)

        if method == "translate":
            shift = target_median - current_median
            mat_norm[triu_idx] = values + shift
            mat_norm = (mat_norm + mat_norm.T) / 2
            np.fill_diagonal(mat_norm, np.diag(mat) + shift)
            mat_norm = np.clip(mat_norm, -1, 1)
        elif method == "scale":
            q1 = np.nanpercentile(values, 25)
            q3 = np.nanpercentile(values, 75)
            iqr = q3 - q1
            if iqr < 1e-10:
                iqr = 1.0
            ref_q1 = np.nanpercentile(all_values, 25)
            ref_q3 = np.nanpercentile(all_values, 75)
            ref_iqr = ref_q3 - ref_q1
            if ref_iqr < 1e-10:
                ref_iqr = 1.0
            scale_factor = ref_iqr / iqr
            mat_norm[triu_idx] = target_median + (values - current_median) * scale_factor
            mat_norm = (mat_norm + mat_norm.T) / 2
            np.fill_diagonal(mat_norm, target_median + (np.diag(mat) - current_median) * scale_factor)
            mat_norm = np.clip(mat_norm, -1, 1)
        elif method == "robust":
            shift = target_median - current_median
            mat_norm[triu_idx] = values + shift
            mat_norm = (mat_norm + mat_norm.T) / 2
            np.fill_diagonal(mat_norm, np.diag(mat) + shift)
            mat_norm = np.clip(mat_norm, -1, 1)
            new_values = mat_norm[triu_idx]
            new_values = new_values[np.isfinite(new_values)]
            if len(new_values) > 0:
                new_median = np.nanmedian(new_values)
                if abs(new_median - target_median) > 0.01:
                    shift2 = target_median - new_median
                    mat_norm[triu_idx] = mat_norm[triu_idx] + shift2
                    mat_norm = (mat_norm + mat_norm.T) / 2
                    mat_norm = np.clip(mat_norm, -1, 1)

        normalized_mats.append(mat_norm)

    if normalized_mats:
        print(f"  [normalize] Normalized {len(normalized_mats)} matrices to median = {target_median:.3f}")

    return normalized_mats, ref_median

def normalize_fingerprint_matrix(mat, method="median_iqr"):
    # unchanged
    triu_indices = np.triu_indices_from(mat, k=1)
    values = mat[triu_indices]
    values = values[np.isfinite(values)]

    if len(values) == 0:
        return mat.copy()

    mat_norm = mat.copy()

    if method == "median_iqr":
        median = np.nanmedian(values)
        q1 = np.nanpercentile(values, 25)
        q3 = np.nanpercentile(values, 75)
        iqr = q3 - q1
        if iqr < 1e-10:
            iqr = 1.0
        mat_norm[triu_indices] = (values - median) / iqr
        mat_norm = (mat_norm + mat_norm.T) / 2
        np.fill_diagonal(mat_norm, 0)
    elif method == "zscore":
        mean = np.nanmean(values)
        std = np.nanstd(values)
        if std < 1e-10:
            std = 1.0
        mat_norm[triu_indices] = (values - mean) / std
        mat_norm = (mat_norm + mat_norm.T) / 2
        np.fill_diagonal(mat_norm, 0)
    elif method == "range":
        vmin = np.nanmin(values)
        vmax = np.nanmax(values)
        if vmax - vmin < 1e-10:
            return mat_norm
        mat_norm[triu_indices] = 2 * (values - vmin) / (vmax - vmin) - 1
        mat_norm = (mat_norm + mat_norm.T) / 2
        np.fill_diagonal(mat_norm, 0)

    return mat_norm

# =============================================================================
# REVISED DUAL SCORE COMPUTATION (top‑N primary pairs, common ROIs for species)
# =============================================================================
def compute_dual_scores_from_matrix(mat, rois,
                                    common_rois=None,
                                    n_top_homotopic=3,
                                    use_inter_all=True):
    """
    Compute primary (top‑N) and species‑wide scores.

    Parameters
    ----------
    mat : 2D ndarray
        Correlation matrix.
    rois : list of str
        ROI names (with L_/R_ prefixes).
    common_rois : list of str, optional
        Base names (no prefix) to restrict species_intra to.
    n_top_homotopic : int, default 3
        Number of strongest homotopic pairs to use for primary_inter.
    use_inter_all : bool, default True
        If True, inter = ALL cross‑hemisphere pairs between selected left/right regions;
        if False, only contralateral (non‑homotopic) pairs.

    Returns
    -------
    dict with keys:
        primary_intra, primary_inter, primary_delta (bounded specificity),
        species_intra, species_inter, species_delta,
        n_homotopic_pairs_primary, n_homotopic_pairs_all,
        networks_found
    """
    import re

    # ---- Helper functions ----
    def _base(name):
        return re.sub(r"^[LlRr]_", "", str(name)).lower().strip()

    def _get_homotopic_pairs(rois_list, mat_data, restrict_bases=None):
        l_names = [r for r in rois_list if r.startswith(("L_", "l_"))]
        r_names = [r for r in rois_list if r.startswith(("R_", "r_"))]
        restrict_norm = ({b.lower().strip() for b in restrict_bases}
                         if restrict_bases is not None else None)
        pairs = []
        for lname in l_names:
            base = _base(lname)
            if restrict_norm is not None and base not in restrict_norm:
                continue
            rname = next((r for r in r_names if _base(r) == base), None)
            if rname is not None:
                pairs.append((lname, rname))

        if not pairs:
            return [], [], []
        name_to_idx = {name: i for i, name in enumerate(rois_list)}
        l_idx = [name_to_idx[l] for l, _ in pairs]
        r_idx = [name_to_idx[r] for _, r in pairs]
        return pairs, l_idx, r_idx

    # ---- 1. Primary scores (using only primary ROIs) ----
    keep_indices = []
    for i, name in enumerate(rois):
        for pattern in _PRIMARY_PATTERNS:
            if pattern.lower() in name.lower():
                keep_indices.append(i)
                break
    if keep_indices:
        mat_primary = mat[np.ix_(keep_indices, keep_indices)]
        rois_primary = [rois[i] for i in keep_indices]
    else:
        mat_primary = mat
        rois_primary = rois

    # Get all homotopic pairs within primary ROIs
    homo_pairs_primary, l_idx_primary, r_idx_primary = _get_homotopic_pairs(rois_primary, mat_primary)

    # Compute primary_intra as median of all homotopic pairs in primary ROIs
    primary_intra_vals = []
    for i, j in zip(l_idx_primary, r_idx_primary):
        r = mat_primary[i, j]
        if np.isfinite(r):
            primary_intra_vals.append(float(r))
    primary_intra = float(np.nanmedian(primary_intra_vals)) if primary_intra_vals else np.nan

    # ---- Select top N homotopic pairs for primary_inter ----
    pair_corrs = []
    for idx, (i, j) in enumerate(zip(l_idx_primary, r_idx_primary)):
        r = mat_primary[i, j]
        if np.isfinite(r):
            pair_corrs.append((r, i, j))
    pair_corrs.sort(key=lambda x: x[0], reverse=True)
    n_top = min(n_top_homotopic, len(pair_corrs))
    top_pairs = pair_corrs[:n_top]
    top_l_idx = [i for _, i, _ in top_pairs]
    top_r_idx = [j for _, _, j in top_pairs]

    # Primary inter: cross-hemisphere correlations between selected left and right regions
    if top_l_idx and top_r_idx:
        if use_inter_all:
            inter_vals = []
            for li in top_l_idx:
                for rj in top_r_idx:
                    # skip homotopic pairs (where li and rj correspond to same base)
                    is_homotopic = any(li == i and rj == j for _, i, j in top_pairs)
                    if not is_homotopic:
                        r = mat_primary[li, rj]
                        if np.isfinite(r):
                            inter_vals.append(float(r))
        else:
            # contralateral only: between different regions
            inter_vals = []
            for i in range(len(top_l_idx)):
                for j in range(i+1, len(top_l_idx)):
                    r1 = mat_primary[top_l_idx[i], top_r_idx[j]]
                    if np.isfinite(r1):
                        inter_vals.append(float(r1))
                    r2 = mat_primary[top_r_idx[i], top_l_idx[j]]
                    if np.isfinite(r2):
                        inter_vals.append(float(r2))
        primary_inter = float(np.nanmedian(inter_vals)) if inter_vals else np.nan
    else:
        primary_inter = np.nan



    # Primary delta (bounded specificity) for reference
    primary_delta = _bounded_specificity(primary_intra, primary_inter)

    # ---- 2. Species scores (using ALL homotopic pairs, optionally restricted to common_rois) ----
    restrict_bases = set(common_rois) if common_rois is not None else None
    homo_pairs_all, l_idx_all, r_idx_all = _get_homotopic_pairs(
        rois, mat, restrict_bases=restrict_bases
    )

    species_intra_vals = []
    for i, j in zip(l_idx_all, r_idx_all):
        r = mat[i, j]
        if np.isfinite(r):
            species_intra_vals.append(float(r))
    species_intra = float(np.nanmedian(species_intra_vals)) if species_intra_vals else np.nan

    # Species inter: we reuse primary_inter (cross from top primary pairs)
    species_inter = primary_inter

    # ---- 2. Species scores (using ALL homotopic pairs, optionally restricted to common_rois) ----
    restrict_bases = set(common_rois) if common_rois is not None else None
    homo_pairs_all, l_idx_all, r_idx_all = _get_homotopic_pairs(
        rois, mat, restrict_bases=restrict_bases
    )

    species_intra_vals = []
    for i, j in zip(l_idx_all, r_idx_all):
        r = mat[i, j]
        if np.isfinite(r):
            species_intra_vals.append(float(r))
    species_intra = float(np.nanmedian(species_intra_vals)) if species_intra_vals else np.nan

    # NEW: species-wide cross-network — ALL cross-hemisphere pairs over the SAME
    # common-ROI pool used for species_intra (symmetric counterpart, not top-N primary)
    if l_idx_all and r_idx_all:
        cross_vals_all = []
        for i, li in enumerate(l_idx_all):
            for j, rj in enumerate(r_idx_all):
                if i == j:
                    continue  # skip the homotopic pair itself
                r = mat[li, rj]
                if np.isfinite(r):
                    cross_vals_all.append(float(r))
        species_inter_all = float(np.nanmedian(cross_vals_all)) if cross_vals_all else np.nan
    else:
        species_inter_all = np.nan

    # Species inter used for CLASSIFICATION stays primary_inter (top-N) — unchanged
    species_inter = primary_inter
    species_delta = (float(species_intra - species_inter)
                     if np.isfinite(species_intra) and np.isfinite(species_inter) else np.nan)

    return {
        "primary_intra": primary_intra,
        "primary_inter": primary_inter,
        "primary_delta": primary_delta,
        "n_homotopic_pairs_primary": len(homo_pairs_primary),
        "species_intra": species_intra,
        "species_inter": species_inter,
        "species_inter_all": species_inter_all,   # NEW — descriptive only, not used for classification
        "species_delta": species_delta,
        "n_homotopic_pairs_all": len(homo_pairs_all),
        "networks_found": [l for l, _ in homo_pairs_primary],
    }

def _bounded_specificity(intra, inter):
    if not (np.isfinite(intra) and np.isfinite(inter)):
        return np.nan
    z_intra = np.arctanh(np.clip(intra, -0.999, 0.999))
    z_inter = np.arctanh(np.clip(inter, -0.999, 0.999))
    return float(np.tanh(z_intra - z_inter))

# =============================================================================
# STEP 1 — SCORE ALL RUNS (with common_rois and n_top)
# =============================================================================

def score_all_runs(df_qc, fit_kind,
                   bids_root_template=_DEFAULT_BIDS,
                   atlas_name="EDNIxCSC", atlas_level=2,
                   use_lr=True, use_inter_all=True,
                   common_rois=None, n_top_homotopic=3,
                   audit_csv_path=None, verbose=True):
    """
    Score every run found under df_qc rows.
    Computes primary (top‑N) and species scores.
    """
    all_runs = []
    seen_keys = set()

    qc_key_cols = ["species", "bids_dir", "subject", "session"]
    if all(c in df_qc.columns for c in qc_key_cols):
        before = len(df_qc)
        df_qc = df_qc.drop_duplicates(subset=qc_key_cols)
        if verbose and before != len(df_qc):
            print(f"  [score] Deduplicated df_qc: {before} → {len(df_qc)} rows")

    n_scored = 0
    n_skipped_flattened = 0
    n_skipped_duplicates = 0

    for idx, row in df_qc.iterrows():
        run_infos = _get_all_run_paths_for_row(
            fit_kind, row, bids_root_template, atlas_name, atlas_level, use_lr)

        for run_path, sp, bd, sub, ses in run_infos:
            run_basename = os.path.basename(run_path)

            if "flattened" in run_basename:
                n_skipped_flattened += 1
                continue

            run_key = (sp, bd, sub, ses, run_basename)
            if run_key in seen_keys:
                n_skipped_duplicates += 1
                continue
            seen_keys.add(run_key)

            try:
                rois, mat = _load_matrix(run_path)
                ps = compute_dual_scores_from_matrix(
                    mat, rois,
                    common_rois=common_rois,
                    n_top_homotopic=n_top_homotopic,
                    use_inter_all=use_inter_all
                )
                run_row = {
                    "species": sp,
                    "bids_dir": bd,
                    "subject": sub,
                    "session": ses,
                    "run": run_basename,
                    "corr_matrix_path": run_path,
                    # Primary scores
                    "primary_intra": ps["primary_intra"],
                    "primary_inter": ps["primary_inter"],
                    "primary_delta": ps["primary_delta"],
                    # Species scores
                    "species_intra": ps["species_intra"],
                    "species_inter": ps["species_inter"],
                    "species_inter_all": ps["species_inter_all"],
                    "species_delta": ps["species_delta"],
                    # Metadata
                    "n_networks_found": len(ps["networks_found"]),
                    "n_homotopic_pairs_primary": ps.get("n_homotopic_pairs_primary", 0),
                    "n_homotopic_pairs_all": ps.get("n_homotopic_pairs_all", 0),
                }
                for col in df_qc.columns:
                    if col not in run_row:
                        run_row[col] = row.get(col, np.nan)
                all_runs.append(run_row)
                n_scored += 1
            except Exception as e:
                warnings.warn(f"[score] {run_path}: {e}")

    if verbose:
        if n_skipped_flattened > 0:
            print(f"  [score] Skipped {n_skipped_flattened} flattened matrices")
        if n_skipped_duplicates > 0:
            print(f"  [score] Skipped {n_skipped_duplicates} duplicate runs")

    if not all_runs:
        raise ValueError("No runs could be loaded!")
    df_scored = pd.DataFrame(all_runs)

    duplicate_cols = ["species", "bids_dir", "subject", "session", "run"]
    if all(c in df_scored.columns for c in duplicate_cols):
        before = len(df_scored)
        df_scored = df_scored.drop_duplicates(subset=duplicate_cols, keep='first')
        if verbose and before != len(df_scored):
            print(f"  [score] Final dedup: {before} → {len(df_scored)} rows")

    if verbose and "species" in df_scored.columns:
        print(f"\n  [score] Final runs per species:")
        for sp in sorted(df_scored["species"].unique()):
            n = len(df_scored[df_scored["species"] == sp])
            print(f"    {sp}: {n} runs")

    return df_scored

# =============================================================================
# STEP 2 — THRESHOLDS (now on species_intra and primary_inter)
# =============================================================================

def suggest_thresholds(df_scored,
                       thresh_intra=None, thresh_inter=None,
                       stringency=0.0, verbose=True,
                       intra_percentile=5, inter_percentile=75):   # 30 → 70
    """
    Percentile thresholds on pooled cross‑species distribution.

    T1 = percentile of species_intra (intra_percentile).
         Specific requires intra > T1, so a LOW percentile (30) keeps T1
         near the bottom of the distribution — lenient, ~70% of runs clear it.

    T2 = percentile of primary_inter (inter_percentile).
         Specific requires inter < T2, so a HIGH percentile (70) keeps T2
         near the top of the distribution — lenient, ~70% of runs fall
         inferior to it (pass). This mirrors T1's leniency symmetrically.
    """
    intra_vals = df_scored["species_intra"].dropna().values
    inter_vals = df_scored["primary_inter"].dropna().values

    if len(intra_vals) < 10:
        t_intra = 0.2
        t_inter = 0.2
    else:
        s = np.clip(stringency, -1.0, 1.0)

        # Intra threshold (T1)
        p_intra = np.percentile(intra_vals, intra_percentile)
        p75_intra = np.percentile(intra_vals, 75)
        p05_intra = np.percentile(intra_vals, 5)
        if s >= 0:
            t_intra = p_intra + s * (p75_intra - p_intra)
        else:
            t_intra = p_intra + abs(s) * (p_intra - p05_intra)
        t_intra = float(np.clip(t_intra, 0.05, 0.8))

        # Inter threshold (T2) – low percentile to select low inter values
        if len(inter_vals) > 0:
            p_inter = np.percentile(inter_vals, inter_percentile)
            p75_inter = np.percentile(inter_vals, 75)
            p05_inter = np.percentile(inter_vals, 5)
        else:
            p_inter = 0.2
            p75_inter = 0.4
            p05_inter = 0.0

        if s >= 0:
            t_inter = p_inter + s * (p75_inter - p_inter)
        else:
            t_inter = p_inter + abs(s) * (p_inter - p05_inter)
        t_inter = float(np.clip(t_inter, -0.1, 0.9))

    if thresh_intra is not None:
        t_intra = float(thresh_intra)
    if thresh_inter is not None:
        t_inter = float(thresh_inter)

    # Count Specific runs (species_intra > T1 AND primary_inter < T2)
    n_spec = 0
    for intra_s, inter_p in zip(df_scored["species_intra"].values,
                                df_scored["primary_inter"].values):
        if np.isfinite(intra_s) and np.isfinite(inter_p):
            if intra_s > t_intra and inter_p < t_inter:
                n_spec += 1

    if verbose:
        print(f"\n{'=' * 50}")
        print(f"  THRESHOLD SELECTION  (stringency={stringency:+.1f})")
        print(f"  T1 (intra) > {t_intra:.3f}  (P{intra_percentile} of species_intra)")
        print(f"  T2 (inter) < {t_inter:.3f}  (P{inter_percentile} of primary_inter — Specific needs inter INFERIOR to T2)")
        print(f"  -> {n_spec}/{len(intra_vals)} runs Specific")
        if len(inter_vals) > 0:
            print(f"  inter range: {inter_vals.min():.3f} to {inter_vals.max():.3f}")
            print(f"  inter P5: {np.percentile(inter_vals, 5):.3f}")
            print(f"  inter P{inter_percentile}: {np.percentile(inter_vals, inter_percentile):.3f}")
            print(f"  inter P50: {np.percentile(inter_vals, 50):.3f}")
            print(f"  inter P75: {np.percentile(inter_vals, 75):.3f}")
        print(f"{'=' * 50}\n")

    return {
        "thresh_intra": t_intra,
        "thresh_inter": t_inter,
    }

# =============================================================================
# STEP 3 — CLASSIFY (updated)
# =============================================================================

def classify_runs(df_scored, thresh_intra, thresh_inter, verbose=True):
    """Classify all runs using two thresholds."""
    cats = pd.Series("No", index=df_scored.index, dtype=object)
    intra_primary_arr = df_scored["primary_intra"].values.astype(float)
    inter_primary_arr = df_scored["primary_inter"].values.astype(float)
    intra_species_arr = df_scored["species_intra"].values.astype(float)

    for i, idx in enumerate(df_scored.index):
        cats[idx] = _classify_one(
            intra_primary_arr[i],
            inter_primary_arr[i],
            intra_species_arr[i],
            thresh_intra,
            thresh_inter
        )

    if verbose:
        print(f"\n  Classification:")
        print(f"    T1 (intra threshold) = {thresh_intra:.4f} (using species_intra)")
        print(f"    T2 (inter threshold) = {thresh_inter:.4f} (using primary_inter)")
        print(f"    Spurious: primary_inter > primary_intra")
        print(f"    No: species_intra <= T1")
        print(f"    Specific: species_intra > T1 AND primary_inter < T2")
        print(f"    Unspecific: All other cases")

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
# EXPORT CLASSIFICATION (updated column names)
# =============================================================================

def export_classification_to_excel(df_scored, cats, output_path, verbose=True):
    """Export classification with new scores."""
    if cats is not None:
        if isinstance(cats, pd.Series):
            df = df_scored.copy()
            df["fp_category"] = cats.values
        else:
            df = df_scored.copy()
            df["fp_category"] = cats
    else:
        df = df_scored.copy()
        if "fp_category" not in df.columns:
            raise ValueError("No category information found in df_scored")

    # ─── Run-level ───
    run_cols = ["species", "bids_dir", "subject", "session", "run",
                "fp_category",
                "primary_intra", "primary_inter", "primary_delta",
                "species_intra", "species_inter", "species_delta",
                "n_homotopic_pairs_primary", "n_homotopic_pairs_all"]
    run_cols = [c for c in run_cols if c in df.columns]
    df_runs = df[run_cols].copy()

    # ─── Subject summary ───
    subject_cols = ["species", "bids_dir", "subject"]
    df_subject = df.groupby(subject_cols).agg({
        "fp_category": lambda x: list(x),
        "primary_intra": "mean",
        "primary_inter": "mean",
        "primary_delta": "mean",
        "species_intra": "mean",
        "species_inter": "mean",
        "species_delta": "mean",
    }).reset_index()

    def count_categories(cat_list):
        counts = {c: 0 for c in CATS}
        for c in cat_list:
            if c in counts:
                counts[c] += 1
        return counts

    df_subject["category_counts"] = df_subject["fp_category"].apply(count_categories)
    df_subject["total_runs"] = df_subject["fp_category"].apply(len)
    df_subject["primary_category"] = df_subject["fp_category"].apply(
        lambda x: max(set(x), key=x.count) if x else "No"
    )

    for cat in CATS:
        df_subject[f"n_{cat}"] = df_subject["category_counts"].apply(lambda x: x.get(cat, 0))

    df_subject["frac_specific"] = df_subject["n_Specific"] / df_subject["total_runs"]
    df_subject = df_subject.sort_values(["species", "bids_dir", "subject"])

    # ─── Session summary ───
    if "session" in df.columns:
        session_cols = ["species", "bids_dir", "subject", "session"]
        df_session = df.groupby(session_cols).agg({
            "fp_category": lambda x: list(x),
            "primary_intra": "mean",
            "primary_inter": "mean",
            "primary_delta": "mean",
            "species_intra": "mean",
            "species_inter": "mean",
            "species_delta": "mean",
        }).reset_index()

        df_session["category_counts"] = df_session["fp_category"].apply(count_categories)
        df_session["total_runs"] = df_session["fp_category"].apply(len)
        df_session["primary_category"] = df_session["fp_category"].apply(
            lambda x: max(set(x), key=x.count) if x else "No"
        )

        for cat in CATS:
            df_session[f"n_{cat}"] = df_session["category_counts"].apply(lambda x: x.get(cat, 0))

        df_session["frac_specific"] = df_session["n_Specific"] / df_session["total_runs"]
        df_session = df_session.sort_values(["species", "bids_dir", "subject", "session"])
    else:
        df_session = pd.DataFrame({"Note": ["No session column found in data"]})

    # ─── Write to Excel ───
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        df_runs.to_excel(writer, sheet_name="Run_Classification", index=False)
        df_subject.to_excel(writer, sheet_name="Subject_Summary", index=False)
        df_session.to_excel(writer, sheet_name="Session_Summary", index=False)

        # Subjects to exclude (no Specific runs)
        non_specific_subjects = df_subject[df_subject["n_Specific"] == 0]
        if non_specific_subjects.empty:
            pd.DataFrame({"Note": ["All subjects have at least one Specific run"]}).to_excel(
                writer, sheet_name="Subjects_To_Exclude", index=False
            )
        else:
            exclude_subjects = non_specific_subjects[["species", "bids_dir", "subject", "total_runs",
                                                       "n_Specific", "n_Unspecific", "n_No", "n_Spurious",
                                                       "primary_intra", "species_intra"]].copy()
            exclude_subjects["reason"] = "No Specific runs"
            if "session" in df.columns:
                session_info = df.groupby(["species", "bids_dir", "subject"])["session"].unique().apply(list).reset_index()
                session_info.columns = ["species", "bids_dir", "subject", "sessions"]
                exclude_subjects = exclude_subjects.merge(session_info, on=["species", "bids_dir", "subject"], how="left")
            exclude_subjects.to_excel(writer, sheet_name="Subjects_To_Exclude", index=False)

        # Summary
        summary_data = {
            "Metric": ["Total Runs", "Total Subjects", "Total Sessions",
                       "Specific Runs", "Unspecific Runs", "No Runs", "Spurious Runs",
                       "Subjects with ≥1 Specific", "Subjects with No Specific"],
            "Value": [
                len(df_runs),
                df_subject["subject"].nunique(),
                df_session["session"].nunique() if "session" in df.columns else "N/A",
                df[df["fp_category"] == "Specific"].shape[0],
                df[df["fp_category"] == "Unspecific"].shape[0],
                df[df["fp_category"] == "No"].shape[0],
                df[df["fp_category"] == "Spurious"].shape[0],
                df_subject[df_subject["n_Specific"] > 0].shape[0],
                df_subject[df_subject["n_Specific"] == 0].shape[0]
            ]
        }
        pd.DataFrame(summary_data).to_excel(writer, sheet_name="Summary", index=False)

        # Species summary
        species_summary = df.groupby("species").agg({
            "fp_category": lambda x: list(x),
            "subject": "nunique",
            "primary_intra": "mean",
            "species_intra": "mean",
        }).reset_index()
        species_summary["total_runs"] = species_summary["fp_category"].apply(len)
        for cat in CATS:
            species_summary[f"n_{cat}"] = species_summary["fp_category"].apply(
                lambda x: sum(1 for c in x if c == cat)
            )
        species_summary["frac_specific"] = species_summary["n_Specific"] / species_summary["total_runs"]
        species_summary = species_summary.drop(columns=["fp_category"])
        species_summary.to_excel(writer, sheet_name="Species_Summary", index=False)

        # Formatted list_to_remove
        if not non_specific_subjects.empty:
            remove_list = []
            for _, row in non_specific_subjects.iterrows():
                if "session" in df.columns:
                    sessions = df[df["subject"] == row["subject"]]["session"].unique().tolist()
                    for ses in sessions:
                        remove_list.append(f"('{row['subject']}', '{ses}')")
                else:
                    remove_list.append(f"('{row['subject']}', '1')")
            remove_text = "[\n    " + ",\n    ".join(remove_list) + "\n]"
            pd.DataFrame({
                "list_to_remove": [remove_text],
                "instruction": ["Copy this entire cell into your species_config list_to_remove"]
            }).to_excel(writer, sheet_name="Copy_This_To_Remove", index=False)

    if verbose:
        print(f"\n  [export] Classification exported to: {output_path}")
        print(f"    {len(df_runs)} runs, {df_subject['subject'].nunique()} subjects")
        print(f"    {df_subject[df_subject['n_Specific'] == 0].shape[0]} subjects with NO Specific runs")
        if not non_specific_subjects.empty:
            print(f"\n  Subjects to exclude (no Specific runs):")
            if "session" in df.columns:
                for _, row in non_specific_subjects.iterrows():
                    sessions = df[df["subject"] == row["subject"]]["session"].unique().tolist()
                    ses_str = ", ".join([f"'{s}'" for s in sessions])
                    print(f"    ('{row['subject']}', [{ses_str}])")
            else:
                for _, row in non_specific_subjects.iterrows():
                    print(f"    ('{row['subject']}', '1')")

    return df_runs, df_subject, df_session

# =============================================================================
# FIGURES (mostly unchanged, except threshold labels and classification space)
# =============================================================================

def plot_threshold_selection(df_scored, thresh_intra=None, thresh_inter=None,
                             output_path=None):
    df = df_scored.copy()
    # Required columns
    required = ["species_intra", "primary_inter", "primary_intra"]
    if not all(c in df.columns for c in required):
        print("  [Fig 01] Missing required columns, skipping.")
        return None

    # Check if species_intra is all NaN
    if df["species_intra"].isnull().all():
        print("  [Fig 01] species_intra is all NaN — no homotopic pairs found? Skipping.")
        return None

    # Drop rows where all required are NaN (but keep rows with some data)
    df = df.dropna(subset=required, how="all")
    if df.empty:
        print("  [Fig 01] No valid data after dropping NaN, skipping.")
        return None


    if thresh_intra is None or thresh_inter is None:
        t = suggest_thresholds(df, verbose=False)
        thresh_intra = thresh_intra if thresh_intra is not None else t["thresh_intra"]
        thresh_inter = thresh_inter if thresh_inter is not None else t["thresh_inter"]

    bids_list = sorted(df["bids_dir"].dropna().unique()) if "bids_dir" in df else []
    sp_list = _sp_order(df)
    sp_colors = _gcolors(sp_list)
    bids_species = {}
    if "species" in df.columns:
        for bd in bids_list:
            sp_counts = df.loc[df["bids_dir"] == bd, "species"].value_counts()
            bids_species[bd] = sp_counts.index[0] if len(sp_counts) else None

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(f"Threshold selection — T1>{thresh_intra:.3f}, T2<{thresh_inter:.3f}",
                 fontsize=12, fontweight="bold")

    # Row 1: PRIMARY scores (primary_intra and primary_inter)
    ax = axes[0, 0]
    col = "primary_intra"
    for bd in bids_list:
        v = df.loc[df["bids_dir"] == bd, col].dropna().values
        v = v[np.isfinite(v) & (v >= -1) & (v <= 1)]
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
    g = df[col].dropna().values
    g = g[np.isfinite(g) & (g >= -1) & (g <= 1)]
    if len(g) >= 5:
        xg = np.linspace(g.min(), g.max(), 500)
        kg = _stats.gaussian_kde(g, bw_method="silverman")
        ax.fill_between(xg, kg(xg), color="#0072B2", alpha=0.12, zorder=2)
        ax.plot(xg, kg(xg), color="#0072B2", lw=2.2, zorder=3, label="global")
    ax.set_xlabel("Primary intra‑network r", fontsize=9)
    ax.set_ylabel("Density", fontsize=9)
    ax.set_title("PRIMARY: Homotopic connectivity", fontsize=10, fontweight="bold")
    ax.set_ylim(bottom=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[0, 1]
    col = "primary_inter"
    thr = thresh_inter
    for bd in bids_list:
        v = df.loc[df["bids_dir"] == bd, col].dropna().values
        v = v[np.isfinite(v) & (v >= -1) & (v <= 1)]
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
    g = df[col].dropna().values
    g = g[np.isfinite(g) & (g >= -1) & (g <= 1)]
    if len(g) >= 5:
        xg = np.linspace(g.min(), g.max(), 500)
        kg = _stats.gaussian_kde(g, bw_method="silverman")
        ax.fill_between(xg, kg(xg), color="#0072B2", alpha=0.12, zorder=2)
        ax.plot(xg, kg(xg), color="#0072B2", lw=2.2, zorder=3, label="global")
    trans = ax.get_xaxis_transform()
    if len(g) > 0:
        pctl = _stats.percentileofscore(g, thr)
        ax.axvline(thr, color="#D55E00", lw=2.4, zorder=5)
        ax.text(thr, 0.98, f" T2={thr:.3f}\n (≈P{pctl:.0f})",
                color="#D55E00", fontsize=8, ha="left", va="top",
                fontweight="bold", transform=trans)
        p25 = np.percentile(g, 25)
        p75 = np.percentile(g, 75)
        ax.axvline(p25, color="#888888", lw=1.0, ls=":", zorder=4)
        ax.axvline(p75, color="#888888", lw=1.0, ls=":", zorder=4)
        ax.text(p25, 0.92, f"P25={p25:.3f}", rotation=90, fontsize=6,
                color="#888888", ha="right", va="top", transform=trans)
        ax.text(p75, 0.92, f"P75={p75:.3f}", rotation=90, fontsize=6,
                color="#888888", ha="right", va="top", transform=trans)
    ax.set_xlabel("Primary inter‑network r (T2 threshold)", fontsize=9)
    ax.set_ylabel("Density", fontsize=9)
    ax.set_title("PRIMARY: Cross‑network connectivity", fontsize=10, fontweight="bold")
    ax.set_ylim(bottom=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # 2D classification space: primary_intra vs primary_inter
    ax = axes[0, 2]
    gc = sp_colors
    markers = ["o", "s", "^", "D", "v", "P", "*", "h"]
    for i, sp in enumerate(sp_list):
        mask = df["species"] == sp if "species" in df.columns else pd.Series(True, index=df.index)
        x = df.loc[mask, "primary_intra"].values.astype(float)
        y = df.loc[mask, "primary_inter"].values.astype(float)
        fin = np.isfinite(x) & np.isfinite(y)
        ax.scatter(x[fin], y[fin], color=gc[sp], marker=markers[i % len(markers)],
                   s=12, alpha=0.5, edgecolors="none", label=sp, zorder=3)
    xmin = float(np.nanmin(df["primary_intra"].dropna())) - 0.01
    xmax = float(np.nanmax(df["primary_intra"].dropna())) + 0.01
    ymin = float(np.nanmin(df["primary_inter"].dropna())) - 0.01
    ymax = float(np.nanmax(df["primary_inter"].dropna())) + 0.01
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    xs_fill = np.linspace(xmin, xmax, 400)
    ax.fill_between(xs_fill, xs_fill, ymax, color=to_rgba(_CAT_COLORS["Spurious"], 0.10), zorder=1)
    ax.plot(xs_fill, xs_fill, color="#bbb", lw=1.0, zorder=2)
    ax.text(xmin + 0.8*(xmax-xmin), ymin + 0.8*(ymax-ymin), "Spurious",
            color=_CAT_COLORS["Spurious"], fontsize=9, fontweight="bold", alpha=0.7)
    ax.text(xmin + 0.1*(xmax-xmin), ymin + 0.8*(ymax-ymin), "Specific\n(needs species_intra > T1 & inter < T2)",
            color=_CAT_COLORS["Specific"], fontsize=8, fontweight="bold", alpha=0.7)
    ax.text(xmin + 0.1*(xmax-xmin), ymin + 0.1*(ymax-ymin), "No / Unspecific",
            color=_CAT_COLORS["Unspecific"], fontsize=9, fontweight="bold", alpha=0.7)
    ax.set_xlabel("Primary intra‑network r", fontsize=9)
    ax.set_ylabel("Primary inter‑network r", fontsize=9)
    ax.set_title("PRIMARY: Classification space", fontsize=10, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Row 2: SPECIES scores (species_intra distribution with T1)
    ax = axes[1, 0]
    col = "species_intra"
    thr = thresh_intra
    for bd in bids_list:
        v = df.loc[df["bids_dir"] == bd, col].dropna().values
        v = v[np.isfinite(v) & (v >= -1) & (v <= 1)]
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
    g = df[col].dropna().values
    g = g[np.isfinite(g) & (g >= -1) & (g <= 1)]
    if len(g) >= 5:
        xg = np.linspace(g.min(), g.max(), 500)
        kg = _stats.gaussian_kde(g, bw_method="silverman")
        ax.fill_between(xg, kg(xg), color="#0072B2", alpha=0.12, zorder=2)
        ax.plot(xg, kg(xg), color="#0072B2", lw=2.2, zorder=3, label="global")
    trans = ax.get_xaxis_transform()
    if len(g) > 0:
        pctl = _stats.percentileofscore(g, thr)
        ax.axvline(thr, color="#D55E00", lw=2.4, zorder=5)
        ax.text(thr, 0.98, f" T1={thr:.3f}\n (≈P{pctl:.0f})",
                color="#D55E00", fontsize=8, ha="left", va="top",
                fontweight="bold", transform=trans)
        p25 = np.percentile(g, 25)
        p75 = np.percentile(g, 75)
        ax.axvline(p25, color="#888888", lw=1.0, ls=":", zorder=4)
        ax.axvline(p75, color="#888888", lw=1.0, ls=":", zorder=4)
        ax.text(p25, 0.92, f"P25={p25:.3f}", rotation=90, fontsize=6,
                color="#888888", ha="right", va="top", transform=trans)
        ax.text(p75, 0.92, f"P75={p75:.3f}", rotation=90, fontsize=6,
                color="#888888", ha="right", va="top", transform=trans)
    ax.set_xlabel("Species intra‑network r (ALL homotopic) T1", fontsize=9)
    ax.set_ylabel("Density", fontsize=9)
    ax.set_title("SPECIES: Homotopic connectivity (ALL pairs)", fontsize=10, fontweight="bold")
    ax.set_ylim(bottom=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel B2 (row2, col1): species-wide cross-network (ALL pairs) — NEW,
    # sits to the right of "SPECIES: Homotopic connectivity" (ALL pairs).
    # No threshold line here — this metric is descriptive, not used for classification.
    ax = axes[1, 1]
    col = "species_inter_all"
    for bd in bids_list:
        v = df.loc[df["bids_dir"] == bd, col].dropna().values
        v = v[np.isfinite(v) & (v >= -1) & (v <= 1)]
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
    g = df[col].dropna().values
    g = g[np.isfinite(g) & (g >= -1) & (g <= 1)]
    if len(g) >= 5:
        xg = np.linspace(g.min(), g.max(), 500)
        kg = _stats.gaussian_kde(g, bw_method="silverman")
        ax.fill_between(xg, kg(xg), color="#0072B2", alpha=0.12, zorder=2)
        ax.plot(xg, kg(xg), color="#0072B2", lw=2.2, zorder=3, label="global")
    trans = ax.get_xaxis_transform()
    if len(g) > 0:
        p25 = np.percentile(g, 25)
        p50 = np.percentile(g, 50)
        p75 = np.percentile(g, 75)
        for p, ls in [(p25, ":"), (p50, "--"), (p75, ":")]:
            ax.axvline(p, color="#888888", lw=1.0, ls=ls, zorder=4)
        ax.text(p25, 0.92, f"P25={p25:.3f}", rotation=90, fontsize=6,
                color="#888888", ha="right", va="top", transform=trans)
        ax.text(p75, 0.92, f"P75={p75:.3f}", rotation=90, fontsize=6,
                color="#888888", ha="right", va="top", transform=trans)
    ax.set_xlabel("Species inter‑network r (ALL pairs, cross‑hemisphere)", fontsize=9)
    ax.set_ylabel("Density", fontsize=9)
    ax.set_title("SPECIES: Cross‑network connectivity (ALL pairs)", fontsize=10, fontweight="bold")
    ax.set_ylim(bottom=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel C2: Species intra vs Primary inter (2D classification space for Specific)
    ax = axes[1, 2]
    for i, sp in enumerate(sp_list):
        mask = df["species"] == sp if "species" in df.columns else pd.Series(True, index=df.index)
        x = df.loc[mask, "species_intra"].values.astype(float)
        y = df.loc[mask, "primary_inter"].values.astype(float)
        fin = np.isfinite(x) & np.isfinite(y)
        ax.scatter(x[fin], y[fin], color=gc[sp], marker=markers[i % len(markers)],
                   s=12, alpha=0.5, edgecolors="none", label=sp, zorder=3)
    xmin = float(np.nanmin(df["species_intra"].dropna())) - 0.01
    xmax = float(np.nanmax(df["species_intra"].dropna())) + 0.01
    ymin = float(np.nanmin(df["primary_inter"].dropna())) - 0.01
    ymax = float(np.nanmax(df["primary_inter"].dropna())) + 0.01
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axvline(thresh_intra, color="#D55E00", lw=1.8, ls="--", zorder=5, alpha=0.6)
    ax.axhline(thresh_inter, color="#D55E00", lw=1.8, ls="--", zorder=5, alpha=0.6)
    # Specific region: intra > T1 AND inter < T2 (bottom‑right)
    ax.fill_between([thresh_intra, xmax], [ymin, ymin], [thresh_inter, thresh_inter],
                    color=to_rgba(_CAT_COLORS["Specific"], 0.15), zorder=1)
    ax.text((thresh_intra + xmax)/2, (ymin + thresh_inter)/2, "Specific",
            color=_CAT_COLORS["Specific"], fontsize=10, fontweight="bold", alpha=0.7, ha="center")
    ax.set_xlabel("Species intra‑network r (ALL homotopic)", fontsize=9)
    ax.set_ylabel("Primary inter‑network r", fontsize=9)
    ax.set_title("SPECIES: Homotopic vs Cross", fontsize=10, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if sp_list:
        sp_handles = [mpatches.Patch(color=sp_colors[sp], alpha=0.7, label=sp)
                      for sp in sp_list]
        axes[1, 2].legend(handles=sp_handles, fontsize=6, frameon=False,
                          loc="upper left", ncol=1, title="species", title_fontsize=6)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    _save_or_show(fig, output_path)
    return fig

# =============================================================================
# FIGURE 02 — CLASSIFICATION MODEL PER BIDS (UPDATED)
# =============================================================================

def plot_classification_model_by_bids(df_scored, thresh_intra, thresh_inter, cats,
                                      output_path=None):
    df = df_scored.copy()
    # Check for required columns
    if "species_intra" not in df.columns or "primary_inter" not in df.columns:
        print("  [Fig 02] Missing required columns, skipping.")
        return None

    # Check if either score is all NaN
    if df["species_intra"].isnull().all() or df["primary_inter"].isnull().all():
        print("  [Fig 02] species_intra or primary_inter all NaN — skipping.")
        return None

    # Drop rows where scores are NaN
    df = df.dropna(subset=["species_intra", "primary_inter"], how="all")
    if df.empty:
        print("  [Fig 02] No valid data, skipping.")
        return None

    bids_list = sorted(df["bids_dir"].dropna().unique()) if "bids_dir" in df else []
    if not bids_list:
        print("  [Fig 02] no BIDS directories found")
        return None

    bids_species = {}
    if "species" in df.columns:
        for bd in bids_list:
            sp_counts = df.loc[df["bids_dir"] == bd, "species"].value_counts()
            bids_species[bd] = sp_counts.index[0] if len(sp_counts) else ""

    n_bids = len(bids_list)
    n_cols = int(np.ceil(np.sqrt(n_bids)))
    n_rows = int(np.ceil(n_bids / n_cols))

    fig = plt.figure(figsize=(n_cols * 6, n_rows * 6))

    # Helper to get safe ranges
    def _safe_range(series, padding=0.01):
        vals = series.dropna().values
        if len(vals) == 0:
            return -0.1, 0.1
        return float(np.nanmin(vals)) - padding, float(np.nanmax(vals)) + padding

    xmin, xmax = _safe_range(df["species_intra"])
    ymin, ymax = _safe_range(df["primary_inter"])

    for idx, bids_name in enumerate(bids_list):
        row = idx // n_cols
        col = idx % n_cols

        gs = gridspec.GridSpecFromSubplotSpec(3, 3,
            subplot_spec=fig.add_gridspec(n_rows, n_cols)[row, col],
            width_ratios=[4, 1, 0.2], height_ratios=[1, 4, 0.2],
            wspace=0.05, hspace=0.05)

        ax_main = fig.add_subplot(gs[1, 0])
        bids_mask = df["bids_dir"] == bids_name
        bids_df = df[bids_mask]

        # Plot all points
        x = bids_df["species_intra"].values.astype(float)
        y = bids_df["primary_inter"].values.astype(float)
        fin = np.isfinite(x) & np.isfinite(y)
        if fin.sum() > 0:
            ax_main.scatter(x[fin], y[fin],
                           color="black", s=25, alpha=0.6,
                           edgecolors="white", lw=0.3, zorder=3)

        # Decision boundaries
        xs_fill = np.linspace(xmin, xmax, 400)
        ax_main.axvline(thresh_intra, color="#222", lw=1.5, ls="--", zorder=5, label=f"T1={thresh_intra:.2f}")
        ax_main.axhline(thresh_inter, color="#222", lw=1.5, ls="--", zorder=5, label=f"T2={thresh_inter:.2f}")

        # FIXED — Specific: intra > T1 AND inter < T2 (below/inferior to threshold)
        ax_main.fill_between(xs_fill, ymin, thresh_inter, where=(xs_fill > thresh_intra),
                             color=to_rgba(_CAT_COLORS["Specific"], 0.08), zorder=1)
        ax_main.fill_between(xs_fill, ymin, ymax, where=(xs_fill <= thresh_intra),
                             color=to_rgba(_CAT_COLORS["No"], 0.08), zorder=1)
        ax_main.fill_between(xs_fill, thresh_inter, ymax, where=(xs_fill > thresh_intra),
                             color=to_rgba(_CAT_COLORS["Unspecific"], 0.08), zorder=1)

        ax_main.set_xlim(xmin, xmax)
        ax_main.set_ylim(ymin, ymax)
        ax_main.tick_params(labelsize=8)
        ax_main.spines["top"].set_visible(False)
        ax_main.spines["right"].set_visible(False)

        if row == n_rows - 1:
            ax_main.set_xlabel("Species intra‑network r (ALL homotopic)", fontsize=8)
        if col == 0:
            ax_main.set_ylabel("Primary inter‑network r", fontsize=8)

        ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)
        x = bids_df["species_intra"].dropna().values
        if len(x) > 2:
            try:
                kde = _stats.gaussian_kde(x, bw_method="silverman")
                x_grid = np.linspace(x.min(), x.max(), 200)
                ax_top.plot(x_grid, kde(x_grid), color="black", lw=1.5, alpha=0.5)
            except:
                pass
        ax_top.axvline(thresh_intra, color="#D55E00", lw=1.2, ls="--", alpha=0.6)
        ax_top.set_xlim(xmin, xmax)
        ax_top.set_ylim(bottom=0)
        ax_top.set_ylabel("Density", fontsize=6, color="black")
        ax_top.tick_params(labelsize=6, colors="black")
        ax_top.spines["top"].set_visible(False)
        ax_top.spines["right"].set_visible(False)
        ax_top.spines["bottom"].set_visible(False)
        ax_top.tick_params(bottom=False)

        ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)
        y = bids_df["primary_inter"].dropna().values
        if len(y) > 2:
            try:
                kde = _stats.gaussian_kde(y, bw_method="silverman")
                y_grid = np.linspace(y.min(), y.max(), 200)
                ax_right.plot(kde(y_grid), y_grid, color="black", lw=1.5, alpha=0.5)
            except:
                pass
        ax_right.axhline(thresh_inter, color="#D55E00", lw=1.2, ls="--", alpha=0.6)
        ax_right.set_ylim(ymin, ymax)
        ax_right.set_xlim(left=0)
        ax_right.set_xlabel("Density", fontsize=6, color="black")
        ax_right.tick_params(labelsize=6, colors="black")
        ax_right.spines["top"].set_visible(False)
        ax_right.spines["right"].set_visible(False)
        ax_right.spines["bottom"].set_visible(False)
        ax_right.tick_params(bottom=False)

        ax_title = fig.add_subplot(gs[0, 1])
        ax_title.axis('off')
        species_label = f" ({bids_species.get(bids_name, '')})" if bids_species.get(bids_name) else ""
        title_text = f"{bids_name}{species_label}"
        ax_title.text(0.5, 0.5, title_text, fontsize=10, fontweight="bold",
                     color="black", ha="center", va="center")

        for gpos in [(0, 2), (1, 2), (2, 0), (2, 1), (2, 2)]:
            ax_corner = fig.add_subplot(gs[gpos[0], gpos[1]])
            ax_corner.axis('off')

    handles = []
    for cat in CATS:
        handles.append(mpatches.Patch(color=to_rgba(_CAT_COLORS[cat], 0.3),
                                      label=f"{cat} region", alpha=0.7))
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, 0.02), ncol=len(CATS),
               fontsize=8, frameon=False)

    fig.suptitle(f"Classification space per BIDS dataset with marginal distributions\n"
                 f"T1 > {thresh_intra:.3f}, T2 < {thresh_inter:.3f}",
                 fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0, 0.04, 1, 0.96])
    _save_or_show(fig, output_path)
    return fig

# =============================================================================
# FIGURE 03 — CATEGORY PROPORTION BARS (UNCHANGED)
# =============================================================================

def plot_category_proportions(df_scored, cats, group_by="species", output_path=None):
    """Plot category proportions."""
    df = df_scored.copy()
    df["fp_category"] = cats.values if hasattr(cats, "values") else cats

    if group_by == "bids_dir" and "bids_dir" in df.columns:
        group_vals = sorted(df["bids_dir"].dropna().unique())
        col_src = df["bids_dir"]
        if "species" in df.columns:
            _bids_sp_map = {}
            for bd in group_vals:
                sp_counts = df.loc[df["bids_dir"] == bd, "species"].value_counts()
                _bids_sp_map[bd] = sp_counts.index[0] if len(sp_counts) else ""
            _bids_display = {bd: f"{bd} ({_bids_sp_map[bd]})" for bd in group_vals}
        else:
            _bids_display = {bd: bd for bd in group_vals}
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

    fig, ax = plt.subplots(figsize=(max(5, n_g * 1.2), 6.5))
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
    display_labels = ([_bids_display.get(g, g) for g in group_vals]
                      if group_by == "bids_dir" and "_bids_display" in dir()
                      else group_vals)
    ax.set_xticklabels(display_labels, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("% runs", fontsize=10)
    ax.set_ylim(0, 115)

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28),
              ncol=len(CATS), fontsize=8, frameon=False)

    plt.tight_layout(rect=[0, 0, 1, 0.90])
    _save_or_show(fig, output_path)
    return fig

# =============================================================================
# FIGURE 04 — FC MATRICES BY CATEGORY (UNCHANGED)
# =============================================================================

def plot_fc_matrices_by_category(df_scored, cats,
                                 bids_root_template=None,
                                 atlas_name="EDNIxCSC", atlas_level=2,
                                 use_lr=True, group_by="species",
                                 species_subset=None,
                                 common_rois=None,
                                 output_path=None, figsize=None,
                                 normalize_method="translate"):
    """Plot FC matrices by category with normalized Specific row."""
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

    cats_order = ["all", "Specific", "Specific_norm", "Unspecific", "No"]
    cat_label = {
        "all": "All",
        "Specific": "Specific",
        "Specific_norm": "Specific\n(norm)",
        "Unspecific": "Unspecific",
        "No": "No"
    }

    all_matrices, all_roi_sets = {}, []

    for cat in cats_order:
        all_matrices[cat] = {}
        for g in all_groups:
            gm = _gmask(df, g)
            if cat == "Specific_norm":
                subset = df[gm & (df["fp_category"] == "Specific")]
            else:
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
            print(f"  [Fig 04] using externally-supplied common ROI list ({len(common_rois)} bilateral)")
        else:
            common_unfiltered = set.intersection(*all_roi_sets) if len(all_roi_sets) > 1 else all_roi_sets[0]
            common_rois = _bilateral_pair_filter(list(common_unfiltered))
            print(f"  [Fig 04] computed local common ROI list ({len(common_rois)} bilateral)")
    else:
        common_rois = common_rois or []
    if len(common_rois) < 4:
        print("  [Fig 04] Not enough common bilateral ROIs")
        return None

    results = {}
    n_rejected_total = 0

    all_specific_mats = []
    for g in all_groups:
        mats = all_matrices["Specific"][g]
        aligned = []
        for rois, mat in mats:
            idx_map = [rois.index(r) for r in common_rois if r in rois]
            if len(idx_map) == len(common_rois):
                aligned.append(mat[np.ix_(idx_map, idx_map)])
        if aligned:
            med_mat = np.nanmedian(np.stack(aligned, 0), 0)
            all_specific_mats.append(med_mat)

    normalized_specific_mats = []
    ref_median = 0
    if all_specific_mats and normalize_method is not None:
        normalized_specific_mats, ref_median = normalize_specific_matrices(
            all_specific_mats, target_median=None, method=normalize_method
        )
        print(f"  [Fig 04] Normalized {len(normalized_specific_mats)} Specific matrices")
        print(f"  [Fig 04] Reference median: {ref_median:.3f}")

    specific_idx = 0
    for cat in cats_order:
        results[cat] = {}
        for g in all_groups:
            mats = all_matrices[cat][g]
            aligned = []
            n_rejected_here = 0
            for rois, mat in mats:
                idx_map = [rois.index(r) for r in common_rois if r in rois]
                if len(idx_map) == len(common_rois):
                    aligned.append(mat[np.ix_(idx_map, idx_map)])
                else:
                    n_rejected_here += 1
            if n_rejected_here:
                print(f"    [Fig 04] {g}/{cat}: rejected {n_rejected_here}/{len(mats)} matrices")
                n_rejected_total += n_rejected_here
            n_loaded = len(aligned)

            if aligned:
                if cat == "Specific_norm" and normalized_specific_mats:
                    if specific_idx < len(normalized_specific_mats):
                        med_mat = normalized_specific_mats[specific_idx]
                        specific_idx += 1
                    else:
                        med_mat = None
                else:
                    med_mat = np.nanmedian(np.stack(aligned, 0), 0)
            else:
                med_mat = None

            gm = _gmask(df, g)
            n_in_df = (df[gm].shape[0] if cat == "all"
                       else df[gm & (df["fp_category"] == cat)].shape[0])
            results[cat][g] = (common_rois, med_mat, n_loaded, n_in_df)

    if n_rejected_total:
        print(f"  [Fig 04] total matrices rejected: {n_rejected_total}")

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
        has_data = any(results[cat][g][1] is not None for g in all_groups)
        if not has_data and cat != "all":
            continue

        col_c = _CAT_COLORS.get(cat, "#444") if cat != "all" else "#444"
        if cat == "Specific_norm":
            col_c = "#7B68EE"

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
                label = cat_label.get(cat, cat)
                ax.set_ylabel(label, fontsize=8,
                              color=col_c, fontweight="bold", labelpad=4)
            if r_i == 0:
                ax.set_title(str(g), fontsize=8, fontweight="bold")

    cbar_ax = fig.add_subplot(gs[:, -1])
    label = "Normalized r" if normalize_method is not None else "Median r"
    fig.colorbar(_mpl_cm.ScalarMappable(norm=norm_r, cmap="RdBu_r"),
                 cax=cbar_ax, label=label)
    plt.tight_layout()
    _save_or_show(fig, output_path)
    return fig
# =============================================================================
# FINGERPRINT FUNCTIONS (UNCHANGED)
# =============================================================================

def _edge_test_species_mean(values, species):
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
                                    fdr_alpha=0.05,
                                    normalize_method="median_iqr",
                                    verbose=True):
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
                                         normalize_method=normalize_method,
                                         verbose=verbose)

def build_fc_fingerprint(sp_mats, n_top=10, fp_test="species_mean",
                         fdr_alpha=0.05, common_rois=None,
                         normalize_method="median_iqr", verbose=True):
    if len(sp_mats) < 2:
        raise ValueError(f"Need >=2 species, got {len(sp_mats)}")

    if common_rois is None:
        first = next(iter(sp_mats.values()))
        common = set(first["rois"])
        for v in sp_mats.values():
            common &= set(v["rois"])
        rois = _bilateral_pair_filter(list(common))
    else:
        rois = [r for r in common_rois
                if all(r in v["rois"] for v in sp_mats.values())]
    n_roi = len(rois)
    if n_roi < 4:
        raise ValueError(f"Only {n_roi} common bilateral ROIs.")
    sp_names = list(sp_mats.keys())

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
    stack = np.stack(run_arrs, 0)

    effect = np.full((n_roi, n_roi), np.nan)
    se_mat = np.full((n_roi, n_roi), np.nan)
    pval = np.full((n_roi, n_roi), np.nan)

    triu = np.triu(np.ones((n_roi, n_roi), bool), k=1)
    for i, j in np.argwhere(triu):
        vals = stack[:, i, j]
        if fp_test == "lmm":
            e, p, _ = _edge_test_lmm(vals, run_sp, run_su)
            if not np.isfinite(p):
                e, p, _ = _edge_test_species_mean(vals, run_sp)
        else:
            e, p, _ = _edge_test_species_mean(vals, run_sp)
        effect[i, j] = effect[j, i] = e
        pval[i, j] = pval[j, i] = p
        d = pd.DataFrame({"v": vals, "sp": run_sp}).dropna()
        if not d.empty:
            g = d.groupby("sp")["v"].mean().values
            g = g[np.isfinite(g)]
            if len(g) >= 2:
                se = float(np.std(g, ddof=1) / np.sqrt(len(g)))
            else:
                se = np.nan
        else:
            se = np.nan
        se_mat[i, j] = se_mat[j, i] = se

    if normalize_method is not None:
        effect_norm = normalize_fingerprint_matrix(effect, normalize_method)
        effect = effect_norm

    p_flat = pval[triu]
    n_tested_edges = int(triu.sum())
    n_finite_p = int(np.sum(np.isfinite(p_flat)))
    expected_n = n_roi * (n_roi - 1) // 2
    assert n_tested_edges == expected_n

    q_flat = _bh_fdr(p_flat)
    qval = np.full((n_roi, n_roi), np.nan)
    ii, jj = np.argwhere(triu).T
    qval[ii, jj] = q_flat
    qval[jj, ii] = q_flat

    if verbose:
        print(f"  [fingerprint stats] matrix is {n_roi}x{n_roi}")
        print(f"                      tested {n_tested_edges} UNIQUE edges")
        print(f"                      finite p-values: {n_finite_p}/{n_tested_edges}")
        n_survive = int(np.sum(np.isfinite(q_flat) & (q_flat < fdr_alpha)))
        print(f"                      BH-FDR α={fdr_alpha}: {n_survive}/{n_finite_p} survive")
        if normalize_method is not None:
            print(f"                      Normalization: {normalize_method}")

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
            effect=float(effect[i, j]),
            se=float(se_mat[i, j]) if np.isfinite(se_mat[i, j]) else np.nan,
            pval=float(pval[i, j]),
            qval=float(q) if np.isfinite(q) else np.nan,
            type=_ctype(i, j),
            survived=bool(np.isfinite(q) and q < fdr_alpha),
        ))

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

    return dict(
        rois=rois, n_species=sp_names, n_roi=n_roi,
        effect_mat=effect, se_mat=se_mat, pval_mat=pval, qval_mat=qval,
        homo_mask=homo_mask, cross_mask=cross_mask,
        connections=connections, display=display,
        fingerprint_indices=fp_idx, fingerprint_pattern=fp_pattern,
        fp_test=fp_test, fdr_alpha=fdr_alpha,
        normalize_method=normalize_method,
    )

# =============================================================================
# PLOT 05 — FINGERPRINT CONNECTIONS (UNCHANGED)
# =============================================================================

def plot_fingerprint_connections(fp, n_top_per_type=5, output_path=None):
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

    by_type = {"homotopic": [], "heterotopic": [], "ipsilateral": []}
    for c in fp["connections"]:
        by_type[c["type"]].append(c)
    for ct in by_type:
        by_type[ct] = sorted(by_type[ct],
                             key=lambda c: abs(c["effect"]) if np.isfinite(c["effect"]) else -1,
                             reverse=True)[:n_top_per_type]

    fig_height = max(9, n_top_per_type * 0.7 + 4)
    fig = plt.figure(figsize=(16, fig_height))
    gs = gridspec.GridSpec(4, 2, height_ratios=[1, 1, 1, 0.15],
                           width_ratios=[1.0, 0.9], hspace=0.35, wspace=0.25,
                           left=0.05, right=0.96, top=0.92, bottom=0.10)

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
    n_total_surv = int(np.sum(np.triu(np.isfinite(qmat) & (qmat < alpha), k=1)))
    ax_heat.set_title(f"Edge effect ({fp['fp_test']}) — {n_total_surv} edges survive FDR<{alpha}",
                      fontsize=9, fontweight="bold")
    surv_ij = np.argwhere(np.isfinite(qmat) & (qmat < alpha))
    for i, j in surv_ij:
        if i < j:
            ax_heat.plot(j, i, "*", color=_GOLD, ms=4,
                         markeredgecolor=_GOLD_EDGE, markeredgewidth=0.5, zorder=5)
    cbar_ax_heat = fig.add_subplot(gs[3, 0])
    fig.colorbar(im, cax=cbar_ax_heat, orientation="horizontal",
                 label="Edge effect (estimate)").ax.tick_params(labelsize=7)

    all_e = [c["effect"] for ct in by_type for c in by_type[ct]
             if np.isfinite(c["effect"])]
    e_abs = max(float(np.nanpercentile(np.abs(all_e), 97)), 0.05) if all_e else 0.5
    norm_r = TwoSlopeNorm(vmin=-e_abs, vcenter=0, vmax=e_abs)
    cmap_r = plt.cm.RdBu_r

    for panel_i, ct in enumerate(["homotopic", "heterotopic", "ipsilateral"]):
        ax = fig.add_subplot(gs[panel_i, 1])
        conns = by_type[ct]
        color = type_colors[ct]
        if not conns:
            ax.text(0.5, 0.5, f"No {ct} edges", ha="center",
                    va="center", transform=ax.transAxes, fontsize=9, color="#aaa")
            ax.set_title(type_labels[ct], fontsize=9, fontweight="bold", color=color)
            ax.set_xticks([]); ax.set_yticks([])
            continue
        y_pos = np.arange(len(conns))
        e_vals = [c["effect"] for c in conns]
        se_vals = [c.get("se", np.nan) for c in conns]
        q_vals = [c["qval"] for c in conns]
        p_vals = [c["pval"] for c in conns]
        xerr = [s if np.isfinite(s) else 0 for s in se_vals]
        ax.barh(y_pos, e_vals, color=[cmap_r(norm_r(e)) if np.isfinite(e) else "#ccc"
                                       for e in e_vals],
                edgecolor="white", lw=0.5, height=0.7,
                xerr=xerr, ecolor="#444444", capsize=2, error_kw={"lw": 0.8})
        for yi, (e, q, p) in enumerate(zip(e_vals, q_vals, p_vals)):
            if np.isnan(q):
                sig = ""
                lbl = f"p={p:.3f} (no q)"
            else:
                sig = "***" if q < .001 else "**" if q < .01 else "*" if q < .05 else ""
                lbl = f"q={q:.3f}{sig}"
            x_off = 0.01 if (np.isfinite(e) and e >= 0) else -0.01
            ax.text(e + x_off if np.isfinite(e) else 0, yi,
                    lbl, va="center",
                    ha="left" if (np.isfinite(e) and e >= 0) else "right",
                    fontsize=6,
                    color="#000" if sig else "#666",
                    fontweight="bold" if sig else "normal")
        ax.set_yticks(y_pos)
        ax.set_yticklabels([f"{c['roi_i'][:15]} ↔ {c['roi_j'][:15]}" for c in conns],
                           fontsize=6)
        ax.axvline(0, color="#aaa", lw=0.7)
        ax.set_xlabel("Edge effect", fontsize=7)
        n_surv_here = sum(1 for c in conns
                          if np.isfinite(c["qval"]) and c["qval"] < alpha)
        ax.set_title(f"{type_labels[ct]}  (top {len(conns)} by |effect|; "
                     f"{n_surv_here} survive FDR)",
                     fontsize=7.5, fontweight="bold", color=color)
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
# PLOT 06 — FINGERPRINT HEATMAP (UNCHANGED)
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
                group_vals[k, gi] = float(np.median(conn_acc[k]))

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
# PIPELINES (updated to pass common_rois and n_top_homotopic)
# =============================================================================
def _to_base_names(roi_list):
    """Strip L_/R_ prefix AND lowercase, so base names match what
    _base() produces inside _get_homotopic_pairs (which also lowercases)."""
    out = []
    for r in roi_list:
        s = re.sub(r"^[LlRr]_", "", str(r))
        out.append(s.lower().strip())
    return out

def run_pipeline(df_qc, fit_kind="correlation",
                 bids_root_template=_DEFAULT_BIDS,
                 thresh_intra=None, thresh_inter=None, stringency=0.0,
                 n_top=10, atlas_name="EDNIxCSC", atlas_level=2, use_lr=True,
                 subsets=("all", "primates", "primates_rodents"),
                 fp_test="species_mean", fdr_alpha=0.05,
                 fig_dir=None, verbose=True,
                 use_inter_all=True,
                 normalize_method="translate",
                 n_top_homotopic=3):
    """Full pipeline for one atlas level, with new threshold logic."""
    print("\n" + "=" * 60 + f"\n  EDNiX FC Pipeline — lvl{atlas_level}\n" + "=" * 60)

    def _p(name, subset_name):
        if fig_dir is None:
            return None
        d = os.path.join(fig_dir, f"lvl{atlas_level}", subset_name)
        os.makedirs(d, exist_ok=True)
        return os.path.join(d, name)

    # We'll compute common ROIs first to pass to scoring
    # For the "all" subset we compute common_rois for all species
    # But we need to do this after scoring? Actually we need common_rois before scoring
    # to restrict species_intra. We'll compute common_rois from the available
    # correlation matrices. We'll do a first pass to load all matrices and compute
    # common ROIs, then rescore with that common set.

    # Step 0: gather all run paths and compute common ROIs for the full dataset
    # (this duplicates some work but is okay)
    all_run_paths = []
    for _, row in df_qc.iterrows():
        run_infos = _get_all_run_paths_for_row(
            fit_kind, row, bids_root_template, atlas_name, atlas_level, use_lr)
        for p, sp, bd, sub, ses in run_infos:
            if "flattened" not in os.path.basename(p):
                all_run_paths.append((p, sp, bd, sub, ses))

    # Build a temporary df_scored just to get species info? We'll compute common ROIs
    # by loading matrices directly.
    common_rois_all = compute_subset_common_rois(
        pd.DataFrame(all_run_paths, columns=["corr_matrix_path", "species", "bids_dir", "subject", "session"]),
        species_subset=None, verbose=verbose
    )
    common_rois_bases = _to_base_names(common_rois_all)

    print("\n[1] Scoring all runs (using common ROIs for species_intra)...")
    df_scored = score_all_runs(
        df_qc, fit_kind,
        bids_root_template=bids_root_template,
        atlas_name=atlas_name, atlas_level=atlas_level,
        use_lr=use_lr, use_inter_all=use_inter_all,
        common_rois=common_rois_bases,  # <-- use bases
        n_top_homotopic=n_top_homotopic,
        verbose=verbose
    )

    if df_scored.empty or df_scored["species_intra"].isnull().all():
        print("  [ERROR] No valid scores obtained (species_intra all NaN). "
              "Check matrix files and common ROIs. Exiting.")
        return results

    print("\n[2] Thresholds...")
    thresh = suggest_thresholds(
        df_scored,
        thresh_intra=thresh_intra,
        thresh_inter=thresh_inter,
        stringency=stringency,
        verbose=verbose
    )
    t_intra = thresh["thresh_intra"]
    t_inter = thresh["thresh_inter"]

    print("\n[3] Classifying...")
    cats = classify_runs(df_scored, t_intra, t_inter, verbose=verbose)
    df_scored["fp_category"] = cats.values

    if fig_dir:
        summary_path = _p("00_per_run_summary.csv", "all")
        df_scored.to_csv(summary_path, index=False)
        print(f"  [save] per-run summary: {summary_path}")

        excel_path = _p("classification_results.xlsx", "all")
        export_classification_to_excel(df_scored, cats, excel_path, verbose=True)

    print("\n[4] Fig 01 threshold selection...")
    plot_threshold_selection(df_scored, t_intra, t_inter,
                             output_path=_p("01_threshold_selection_distributions.png", "all"))
    print("\n[5] Fig 02 classification model...")
    plot_classification_model_by_bids(df_scored, t_intra, t_inter, cats,
                                      output_path=_p("02_classification_model_per_bids.png", "all"))
    print("\n[6] Fig 03 category proportions...")
    plot_category_proportions(df_scored, cats, group_by="species",
                              output_path=_p("03_category_proportions_by_species.png", "all"))
    plot_category_proportions(df_scored, cats, group_by="bids_dir",
                              output_path=_p("03_category_proportions_by_bids.png", "all"))

    fingerprints = {}
    for subset_name in subsets:
        species_subset = _SUBSETS.get(subset_name, None)
        if species_subset is not None and "species" in df_scored.columns:
            n_runs = df_scored["species"].isin(species_subset).sum()
            if n_runs == 0:
                print(f"\n  [skip subset {subset_name}] no runs found")
                continue

        print(f"\n  ─── subset: {subset_name} ───")
        # Compute common ROIs within this subset
        common_rois_sub = compute_subset_common_rois(
            df_scored, species_subset=species_subset,
            report_csv_path=_p("00_missing_rois_per_file.csv", subset_name),
            verbose=verbose)

        print(f"  [7-{subset_name}] Fig 04 FC matrices...")
        plot_fc_matrices_by_category(
            df_scored, cats, bids_root_template=bids_root_template,
            atlas_name=atlas_name, atlas_level=atlas_level, use_lr=use_lr,
            group_by="species", species_subset=species_subset,
            common_rois=common_rois_sub,
            normalize_method=normalize_method,
            output_path=_p("04_fc_matrices_by_category_species.png", subset_name))
        plot_fc_matrices_by_category(
            df_scored, cats, bids_root_template=bids_root_template,
            atlas_name=atlas_name, atlas_level=atlas_level, use_lr=use_lr,
            group_by="bids_dir", species_subset=species_subset,
            common_rois=common_rois_sub,
            normalize_method=normalize_method,
            output_path=_p("04_fc_matrices_by_category_bids.png", subset_name))

        print(f"  [8-{subset_name}] Fingerprint...")
        try:
            sp_mats, fp = build_fingerprint_from_specific(
                df_scored, cats, bids_root_template=bids_root_template,
                atlas_name=atlas_name, atlas_level=atlas_level, use_lr=use_lr,
                n_top=n_top, species_subset=species_subset,
                common_rois=common_rois_sub,
                fp_test=fp_test, fdr_alpha=fdr_alpha,
                normalize_method=normalize_method,
                verbose=verbose)
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
                                          common_rois=common_rois_sub)

    print("\n" + "=" * 60 + f"\n  lvl{atlas_level} complete.")
    if fig_dir:
        print(f"  Figures: {os.path.join(fig_dir, f'lvl{atlas_level}')}")
        print(f"  Classification Excel: {os.path.join(fig_dir, f'lvl{atlas_level}', 'all', 'classification_results.xlsx')}")
    print("=" * 60 + "\n")
    return dict(df_scored=df_scored, thresh=thresh, cats=cats,
                fingerprints=fingerprints)

# The _rescore_at_level and run_pipeline_multilevel functions should be updated
# to use the new threshold names (thresh_inter instead of thresh_delta) and
# pass n_top_homotopic. They are straightforward adaptations.
# For brevity, I'll include the updated versions below.

def _rescore_at_level(df_qc, df_scored_lvl2, fit_kind, bids_root_template,
                       atlas_name, atlas_level, use_lr, use_inter_all=True,
                       common_rois=None, n_top_homotopic=3, verbose=True):
    """Rescore at other atlas levels, inheriting lvl-2 categories."""
    df_l = score_all_runs(df_qc, fit_kind,
                           bids_root_template=bids_root_template,
                           atlas_name=atlas_name, atlas_level=atlas_level,
                           use_lr=use_lr, use_inter_all=use_inter_all,
                           common_rois=common_rois,
                           n_top_homotopic=n_top_homotopic,
                           verbose=False)

    def _run_idx(run_name):
        m = re.search(r"_run_(\d+)_", str(run_name))
        return int(m.group(1)) if m else -1

    df_l["_run_idx"] = df_l["run"].apply(_run_idx)
    df_l2 = df_scored_lvl2.copy()
    df_l2["_run_idx"] = df_l2["run"].apply(_run_idx)

    key_cols = ["species", "subject", "session", "_run_idx"]
    cats_only = df_l2[key_cols + ["fp_category"]].drop_duplicates(subset=key_cols)
    merged = df_l.merge(cats_only, on=key_cols, how="left",
                         suffixes=("", "_l2"))

    if verbose:
        n_unmatched = merged["fp_category"].isna().sum()
        n_total = len(merged)
        print(f"  [lvl{atlas_level}] {n_total - n_unmatched}/{n_total} runs matched")
    merged = merged[merged["fp_category"].notna()].drop(columns=["_run_idx"])
    return merged

def run_pipeline_multilevel(df_qc, fit_kind="correlation",
                            bids_root_template=_DEFAULT_BIDS,
                            atlas_levels=(2, 3, 4),
                            subsets=("all", "primates", "primates_rodents"),
                            thresh_intra=None, thresh_inter=None, stringency=0.0,
                            n_top=10, atlas_name="EDNIxCSC", use_lr=True,
                            fp_test="species_mean", fdr_alpha=0.05,
                            fig_dir=None, verbose=True,
                            use_inter_all=True,
                            normalize_method="translate",
                            n_top_homotopic=3):
    """Multi-level pipeline - level 2 defines categories."""
    results = {}

    if 2 not in atlas_levels:
        atlas_levels = (2,) + tuple(l for l in atlas_levels if l != 2)

    print("\n" + "█" * 60 + "\n  LEVEL 2 — full pipeline\n" + "█" * 60)

    if "species" in df_qc.columns:
        for sp in df_qc["species"].unique():
            n = len(df_qc[df_qc["species"] == sp])
            print(f"  [input CSV] {sp}: {n} rows")

    # Compute common ROIs for the entire dataset (for species_intra)
    all_run_paths = []
    for _, row in df_qc.iterrows():
        run_infos = _get_all_run_paths_for_row(
            fit_kind, row, bids_root_template, atlas_name, 2, use_lr)
        for p, sp, bd, sub, ses in run_infos:
            if "flattened" not in os.path.basename(p):
                all_run_paths.append((p, sp, bd, sub, ses))
    common_rois_all = compute_subset_common_rois(
        pd.DataFrame(all_run_paths, columns=["corr_matrix_path", "species", "bids_dir", "subject", "session"]),
        species_subset=None, verbose=verbose
    )

    common_rois_bases = _to_base_names(common_rois_all)

    print("\n[1] Scoring all runs (level 2)...")
    df_scored = score_all_runs(
        df_qc, fit_kind,
        bids_root_template=bids_root_template,
        atlas_name=atlas_name, atlas_level=2,
        use_lr=use_lr, use_inter_all=use_inter_all,
        common_rois=common_rois_bases,  # <-- use bases
        n_top_homotopic=n_top_homotopic,
        verbose=verbose
    )

    # Drop legacy columns if present
    legacy_cat_cols = ["sp_category", "category", "sp_category_legacy", "sp_specificity_index"]
    for col in legacy_cat_cols:
        if col in df_scored.columns:
            df_scored = df_scored.drop(columns=[col])
            if verbose:
                print(f"  [pipeline] Removed legacy column: {col}")

    print("\n[2] Thresholds...")
    thresh = suggest_thresholds(
        df_scored,
        thresh_intra=thresh_intra,
        thresh_inter=thresh_inter,
        stringency=stringency,
        verbose=verbose
    )
    t_intra = thresh["thresh_intra"]
    t_inter = thresh["thresh_inter"]

    print("\n[3] Classifying...")
    cats = classify_runs(df_scored, t_intra, t_inter, verbose=verbose)
    df_scored["fp_category"] = cats.values

    if fig_dir:
        summary_path = os.path.join(fig_dir, "lvl2", "all", "00_per_run_summary.csv")
        os.makedirs(os.path.dirname(summary_path), exist_ok=True)
        df_scored.to_csv(summary_path, index=False)
        print(f"  [save] per-run summary: {summary_path}")

        excel_path = os.path.join(fig_dir, "lvl2", "all", "classification_results.xlsx")
        os.makedirs(os.path.dirname(excel_path), exist_ok=True)
        export_classification_to_excel(df_scored, cats, excel_path, verbose=True)

    print("\n[4] Fig 01 threshold selection...")
    plot_threshold_selection(df_scored, t_intra, t_inter,
                             output_path=os.path.join(fig_dir, "lvl2", "all",
                                                      "01_threshold_selection_distributions.png") if fig_dir else None)

    print("\n[5] Fig 02 classification model...")
    plot_classification_model_by_bids(df_scored, t_intra, t_inter, cats,
                                      output_path=os.path.join(fig_dir, "lvl2", "all",
                                                               "02_classification_model_per_bids.png") if fig_dir else None)

    print("\n[6] Fig 03 category proportions...")
    plot_category_proportions(df_scored, cats, group_by="species",
                              output_path=os.path.join(fig_dir, "lvl2", "all",
                                                       "03_category_proportions_by_species.png") if fig_dir else None)
    plot_category_proportions(df_scored, cats, group_by="bids_dir",
                              output_path=os.path.join(fig_dir, "lvl2", "all",
                                                       "03_category_proportions_by_bids.png") if fig_dir else None)

    results[2] = {"df_scored": df_scored, "thresh": thresh, "cats": cats}
    df_scored_lvl2 = df_scored
    cats_lvl2 = cats

    for subset_name in subsets:
        species_subset = _SUBSETS.get(subset_name, None)
        if species_subset is not None and "species" in df_scored.columns:
            n_runs = df_scored["species"].isin(species_subset).sum()
            if n_runs == 0:
                print(f"\n  [skip subset {subset_name}] no runs found")
                continue

        print(f"\n  ─── subset: {subset_name} ───")
        common_rois_sub = compute_subset_common_rois(
            df_scored, species_subset=species_subset,
            report_csv_path=None, verbose=verbose)

        print(f"  [Fig 04] FC matrices...")
        plot_fc_matrices_by_category(
            df_scored, cats, bids_root_template=bids_root_template,
            atlas_name=atlas_name, atlas_level=2, use_lr=use_lr,
            group_by="species", species_subset=species_subset,
            common_rois=common_rois_sub,
            normalize_method=normalize_method,
            output_path=os.path.join(fig_dir, f"lvl2", subset_name,
                                     "04_fc_matrices_by_category_species.png") if fig_dir else None)
        plot_fc_matrices_by_category(
            df_scored, cats, bids_root_template=bids_root_template,
            atlas_name=atlas_name, atlas_level=2, use_lr=use_lr,
            group_by="bids_dir", species_subset=species_subset,
            common_rois=common_rois_sub,
            normalize_method=normalize_method,
            output_path=os.path.join(fig_dir, f"lvl2", subset_name,
                                     "04_fc_matrices_by_category_bids.png") if fig_dir else None)

        print(f"  [Fig 05-06] Fingerprint...")
        try:
            sp_mats, fp = build_fingerprint_from_specific(
                df_scored, cats, bids_root_template=bids_root_template,
                atlas_name=atlas_name, atlas_level=2, use_lr=use_lr,
                n_top=n_top, species_subset=species_subset,
                common_rois=common_rois_sub,
                fp_test=fp_test, fdr_alpha=fdr_alpha,
                normalize_method=normalize_method,
                verbose=verbose)
            plot_fingerprint_connections(
                fp, output_path=os.path.join(fig_dir, f"lvl2", subset_name,
                                             "05_fingerprint_connections.png") if fig_dir else None)
            plot_fingerprint_heatmap(
                df_scored, fp, cats, bids_root_template=bids_root_template,
                atlas_name=atlas_name, atlas_level=2, use_lr=use_lr,
                species_subset=species_subset,
                output_path=os.path.join(fig_dir, f"lvl2", subset_name,
                                         "06_fingerprint_heatmap_by_bids.png") if fig_dir else None)
        except ValueError as e:
            print(f"  [skip fingerprint {subset_name}] {e}")

    for lvl in atlas_levels:
        if lvl == 2:
            continue
        print("\n" + "█" * 60 + f"\n  LEVEL {lvl} — figs 04-06 only\n" + "█" * 60)

        try:
            df_scored_L = _rescore_at_level(
                df_qc, df_scored_lvl2, fit_kind=fit_kind,
                bids_root_template=bids_root_template,
                atlas_name=atlas_name, atlas_level=lvl,
                use_lr=use_lr, use_inter_all=use_inter_all,
                common_rois=common_rois_bases,
                n_top_homotopic=n_top_homotopic,
                verbose=verbose)
            if df_scored_L.empty:
                print(f"  [skip lvl{lvl}] no runs found")
                continue
            cats_L = df_scored_L["fp_category"]

            for subset_name in subsets:
                species_subset = _SUBSETS.get(subset_name, None)
                if species_subset is not None:
                    n_runs = df_scored_L["species"].isin(species_subset).sum()
                    if n_runs == 0:
                        continue

                print(f"\n  ─── lvl{lvl} subset: {subset_name} ───")
                common_rois_sub = compute_subset_common_rois(
                    df_scored_L, species_subset=species_subset,
                    report_csv_path=None, verbose=verbose)
                if not common_rois_sub:
                    print(f"  [skip lvl{lvl}/{subset_name}] no common ROIs")
                    continue

                print(f"  [Fig 04] lvl{lvl} matrices...")
                plot_fc_matrices_by_category(
                    df_scored_L, cats_L,
                    bids_root_template=bids_root_template,
                    atlas_name=atlas_name, atlas_level=lvl, use_lr=use_lr,
                    group_by="species", species_subset=species_subset,
                    common_rois=common_rois_sub,
                    output_path=os.path.join(fig_dir, f"lvl{lvl}", subset_name,
                                             "04_fc_matrices_by_category_species.png") if fig_dir else None)
                plot_fc_matrices_by_category(
                    df_scored_L, cats_L,
                    bids_root_template=bids_root_template,
                    atlas_name=atlas_name, atlas_level=lvl, use_lr=use_lr,
                    group_by="bids_dir", species_subset=species_subset,
                    common_rois=common_rois_sub,
                    output_path=os.path.join(fig_dir, f"lvl{lvl}", subset_name,
                                             "04_fc_matrices_by_category_bids.png") if fig_dir else None)

                print(f"  [Fig 05-06] lvl{lvl} fingerprint...")
                try:
                    sp_mats, fp = build_fingerprint_from_specific(
                        df_scored_L, cats_L,
                        bids_root_template=bids_root_template,
                        atlas_name=atlas_name, atlas_level=lvl, use_lr=use_lr,
                        n_top=n_top, species_subset=species_subset,
                        common_rois=common_rois_sub,
                        fp_test=fp_test, fdr_alpha=fdr_alpha,
                        normalize_method=normalize_method,
                        verbose=verbose)
                    plot_fingerprint_connections(
                        fp, output_path=os.path.join(fig_dir, f"lvl{lvl}", subset_name,
                                                     "05_fingerprint_connections.png") if fig_dir else None)
                    plot_fingerprint_heatmap(
                        df_scored_L, fp, cats_L,
                        bids_root_template=bids_root_template,
                        atlas_name=atlas_name, atlas_level=lvl, use_lr=use_lr,
                        species_subset=species_subset,
                        output_path=os.path.join(fig_dir, f"lvl{lvl}", subset_name,
                                                 "06_fingerprint_heatmap_by_bids.png") if fig_dir else None)
                except ValueError as e:
                    print(f"  [skip fingerprint lvl{lvl}/{subset_name}] {e}")
        except Exception as e:
            print(f"  [skip lvl{lvl}] {e}")

    return results

# =============================================================================
# CLI (updated)
# =============================================================================

def main():
    p = argparse.ArgumentParser(description="EDNiX FC pipeline with two-threshold classification.")
    p.add_argument("--csv", required=True)
    p.add_argument("--bids-root", default=_DEFAULT_BIDS)
    p.add_argument("--thresh-intra", type=float, default=None, help="T1 for intra (species_intra)")
    p.add_argument("--thresh-inter", type=float, default=None, help="T2 for inter (primary_inter)")
    p.add_argument("--stringency", type=float, default=0.0)
    p.add_argument("--n-top", type=int, default=10, help="Top connections for fingerprint")
    p.add_argument("--n-top-homotopic", type=int, default=3, help="Top homotopic pairs for primary_inter")
    p.add_argument("--levels", type=int, nargs="+", default=[2, 3, 4])
    p.add_argument("--subsets", nargs="+", default=["all", "primates", "primates_rodents"])
    p.add_argument("--fp-test", default="species_mean", choices=["species_mean", "lmm"])
    p.add_argument("--fdr-alpha", type=float, default=0.05)
    p.add_argument("--output", default="figures/fc_explorer")
    p.add_argument("--use-inter-all", action="store_true", default=True)
    p.add_argument("--no-inter-all", dest="use_inter_all", action="store_false")
    p.add_argument("--normalize", default="translate",
                   choices=["translate", "scale", "robust", "none"])
    args = p.parse_args()

    normalize_method = None if args.normalize == "none" else args.normalize

    df_qc = pd.read_csv(args.csv)
    run_pipeline_multilevel(
        df_qc, bids_root_template=args.bids_root,
        atlas_levels=tuple(args.levels), subsets=tuple(args.subsets),
        thresh_intra=args.thresh_intra, thresh_inter=args.thresh_inter,
        stringency=args.stringency, n_top=args.n_top,
        n_top_homotopic=args.n_top_homotopic,
        fp_test=args.fp_test, fdr_alpha=args.fdr_alpha,
        fig_dir=args.output, use_inter_all=args.use_inter_all,
        normalize_method=normalize_method)

if __name__ == "__main__":
    main()