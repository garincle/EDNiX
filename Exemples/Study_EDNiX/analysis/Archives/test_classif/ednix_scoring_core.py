"""
EDNiX — Shared Dual-Threshold Scoring Core
============================================

Single source of truth for the "species_intra / primary_inter" dual-threshold
classification, extracted from ednix_fc_threshold_explorer_simplified.py so
that:
  - the per-run QC script (ednix_qc_matrix_v2.py)
  - the batch lvl2 classifier (ednix_batch_classify_lvl2.py)
  - the threshold explorer itself

...all classify runs the exact same way. If you tune thresholds or the
common-ROI list, change DEFAULTS here (or pass overrides) and every
downstream tool picks it up.

Classification rule (unchanged from the explorer):
    Spurious   : primary_inter > primary_intra
    No         : species_intra <= T1
    Specific   : species_intra >  T1  AND  primary_inter < T2
    Unspecific : everything else
"""

import os
import re
import numpy as np
import pandas as pd

# =============================================================================
# DEFAULTS (per Clement, Aug 2026)
# =============================================================================

# Base ROI names (no L_/R_ prefix) used to restrict species_intra (the
# "ALL homotopic pairs" score). Matching is case-insensitive.
DEFAULT_COMMON_ROIS = [
    'Amygdala',
    'Auditory_cortex_(Superior_temporal)',
    'Hippocampal_formation',
    'Insula_and_others_in_lateral_sulcus',
    'Middle_Temporal_Inferior_temporal_Temporal_pole_(MIPT)',
    'Motor_and_premotor',
    'Olfactory_cortex',
    'Pallidum',
    'Periarchicortex',
    'Posterior_parietal_cortex',
    'Prefrontal_cortex_(PFC)',
    'Somatosensory_cortex',
    'Striatum',
    'Tectum',
    'Thalamus',
    'Visual_pre_and_extra_striate_cortex',
    'Visual_striate_cortex',
]

DEFAULT_THRESH_INTRA = 0.104   # T1: species_intra must be  > this  -> "Specific"/"Unspecific" vs "No"
DEFAULT_THRESH_INTER = 0.212   # T2: primary_inter must be  < this  -> "Specific" vs "Unspecific"
DEFAULT_N_TOP_HOMOTOPIC = 3    # top-N homotopic pairs used for primary_inter

_PRIMARY_PATTERNS = ['Somatosensory', 'Auditory_co', 'Visual_stri']

CATS = ["Specific", "Unspecific", "No", "Spurious"]

# =============================================================================
# LOW-LEVEL HELPERS (identical behaviour to the explorer)
# =============================================================================

def _base(name):
    return re.sub(r"^[LlRr]_", "", str(name)).lower().strip()


def load_matrix(path):
    """Load a correlation matrix CSV -> (roi_names, ndarray)."""
    df = pd.read_csv(path, index_col=0)
    df.columns = df.columns.astype(str).str.strip()
    df.index = df.index.astype(str).str.strip()
    return list(df.index), df.values.astype(float)


def bilateral_pair_filter(rois):
    l_bases = {n[2:] for n in rois if n.startswith(("L_", "l_"))}
    r_bases = {n[2:] for n in rois if n.startswith(("R_", "r_"))}
    common = sorted(l_bases & r_bases)
    return [f"L_{b}" for b in common] + [f"R_{b}" for b in common]


def _get_homotopic_pairs(rois_list, restrict_bases=None):
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


def _bounded_specificity(intra, inter):
    if not (np.isfinite(intra) and np.isfinite(inter)):
        return np.nan
    z_intra = np.arctanh(np.clip(intra, -0.999, 0.999))
    z_inter = np.arctanh(np.clip(inter, -0.999, 0.999))
    return float(np.tanh(z_intra - z_inter))


# =============================================================================
# DUAL SCORE COMPUTATION (same algorithm as ednix_fc_threshold_explorer)
# =============================================================================

def compute_dual_scores_from_matrix(mat, rois,
                                     common_rois=None,
                                     n_top_homotopic=DEFAULT_N_TOP_HOMOTOPIC,
                                     use_inter_all=True):
    """
    Returns dict with primary_intra, primary_inter, species_intra,
    species_inter (== primary_inter, used for classification),
    species_inter_all (descriptive only), plus pair counts.

    common_rois: base ROI names (no L_/R_) to restrict species_intra to.
                 If None, ALL homotopic pairs found in the matrix are used.
    """
    keep_indices = [i for i, name in enumerate(rois)
                     if any(p.lower() in name.lower() for p in _PRIMARY_PATTERNS)]
    if keep_indices:
        mat_primary = mat[np.ix_(keep_indices, keep_indices)]
        rois_primary = [rois[i] for i in keep_indices]
    else:
        mat_primary = mat
        rois_primary = rois

    homo_pairs_primary, l_idx_primary, r_idx_primary = _get_homotopic_pairs(rois_primary)

    primary_intra_vals = [float(mat_primary[i, j]) for i, j in zip(l_idx_primary, r_idx_primary)
                           if np.isfinite(mat_primary[i, j])]
    primary_intra = float(np.nanmedian(primary_intra_vals)) if primary_intra_vals else np.nan

    pair_corrs = [(mat_primary[i, j], i, j) for i, j in zip(l_idx_primary, r_idx_primary)
                  if np.isfinite(mat_primary[i, j])]
    pair_corrs.sort(key=lambda x: x[0], reverse=True)
    n_top = min(n_top_homotopic, len(pair_corrs))
    top_pairs = pair_corrs[:n_top]
    top_l_idx = [i for _, i, _ in top_pairs]
    top_r_idx = [j for _, _, j in top_pairs]

    if top_l_idx and top_r_idx:
        if use_inter_all:
            inter_vals = []
            for li in top_l_idx:
                for rj in top_r_idx:
                    is_homotopic = any(li == i and rj == j for _, i, j in top_pairs)
                    if not is_homotopic:
                        r = mat_primary[li, rj]
                        if np.isfinite(r):
                            inter_vals.append(float(r))
        else:
            inter_vals = []
            for i in range(len(top_l_idx)):
                for j in range(i + 1, len(top_l_idx)):
                    r1 = mat_primary[top_l_idx[i], top_r_idx[j]]
                    if np.isfinite(r1):
                        inter_vals.append(float(r1))
                    r2 = mat_primary[top_r_idx[i], top_l_idx[j]]
                    if np.isfinite(r2):
                        inter_vals.append(float(r2))
        primary_inter = float(np.nanmedian(inter_vals)) if inter_vals else np.nan
    else:
        primary_inter = np.nan

    primary_delta = _bounded_specificity(primary_intra, primary_inter)

    restrict_bases = set(common_rois) if common_rois is not None else None
    homo_pairs_all, l_idx_all, r_idx_all = _get_homotopic_pairs(rois, restrict_bases=restrict_bases)

    species_intra_vals = [float(mat[i, j]) for i, j in zip(l_idx_all, r_idx_all) if np.isfinite(mat[i, j])]
    species_intra = float(np.nanmedian(species_intra_vals)) if species_intra_vals else np.nan

    if l_idx_all and r_idx_all:
        cross_vals_all = []
        for i, li in enumerate(l_idx_all):
            for j, rj in enumerate(r_idx_all):
                if i == j:
                    continue
                r = mat[li, rj]
                if np.isfinite(r):
                    cross_vals_all.append(float(r))
        species_inter_all = float(np.nanmedian(cross_vals_all)) if cross_vals_all else np.nan
    else:
        species_inter_all = np.nan

    species_inter = primary_inter  # used for classification, matches explorer
    species_delta = (float(species_intra - species_inter)
                      if np.isfinite(species_intra) and np.isfinite(species_inter) else np.nan)

    return {
        "primary_intra": primary_intra,
        "primary_inter": primary_inter,
        "primary_delta": primary_delta,
        "n_homotopic_pairs_primary": len(homo_pairs_primary),
        "species_intra": species_intra,
        "species_inter": species_inter,
        "species_inter_all": species_inter_all,
        "species_delta": species_delta,
        "n_homotopic_pairs_all": len(homo_pairs_all),
    }


def classify_one(species_intra, primary_inter, primary_intra,
                  thresh_intra=DEFAULT_THRESH_INTRA, thresh_inter=DEFAULT_THRESH_INTER):
    """Same rule as ednix_fc_threshold_explorer._classify_one."""
    if not np.isfinite(primary_intra) or not np.isfinite(primary_inter) or not np.isfinite(species_intra):
        return "No"
    if primary_inter > primary_intra:
        return "Spurious"
    if species_intra <= thresh_intra:
        return "No"
    if species_intra > thresh_intra and primary_inter < thresh_inter:
        return "Specific"
    return "Unspecific"


def score_and_classify_matrix(mat, rois,
                               common_rois=None,
                               thresh_intra=DEFAULT_THRESH_INTRA,
                               thresh_inter=DEFAULT_THRESH_INTER,
                               n_top_homotopic=DEFAULT_N_TOP_HOMOTOPIC,
                               use_inter_all=True):
    """Convenience one-shot: score a single matrix and classify it."""
    if common_rois is None:
        common_rois = DEFAULT_COMMON_ROIS
    scores = compute_dual_scores_from_matrix(
        mat, rois, common_rois=common_rois,
        n_top_homotopic=n_top_homotopic, use_inter_all=use_inter_all)
    category = classify_one(
        scores["species_intra"], scores["primary_inter"], scores["primary_intra"],
        thresh_intra=thresh_intra, thresh_inter=thresh_inter)
    scores["fp_category"] = category
    scores["thresh_intra"] = thresh_intra
    scores["thresh_inter"] = thresh_inter
    return scores


def score_and_classify_path(matrix_csv_path, **kwargs):
    """Load a matrix CSV from disk, score, classify."""
    rois, mat = load_matrix(matrix_csv_path)
    return score_and_classify_matrix(mat, rois, **kwargs)
