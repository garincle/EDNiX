"""
kept_outlier_cross_species_report.py
──────────────────────────────────────
Publication-style cross-species report comparing Kept vs Outlier runs
(from SBA_Region_Profile_Index.py's classify_runs_by_* / df_class),
adapting three pieces of machinery from ednix_fc_threshold_explorer.py:

  1. The FC-matrices-by-category-by-BIDS grid (fig04 equivalent):
     rows = [all, Kept, Outlier], columns = one per BIDS, common
     cross-species ROIs, shared colorbar.
  2. The "fingerprint" of edges that robustly differ between Kept and
     Outlier ACROSS EVERY BIDS (fig05/fig06 equivalent) — same edge-test
     + FDR + display-selection machinery as the original fingerprint,
     but testing the KEPT-MINUS-OUTLIER difference per BIDS instead of
     the raw per-species mean the original tested.
  3. A combined connectivity + QC variable ranking against the
     Kept/Outlier label (Cohen's d + AUC), reusing
     compute_dual_scores_from_matrix (connectivity, copied from the
     explorer) and build_candidate_qc_table (QC, from SBA region
     profiles — qc_index_kept_outlier_relevance.py).

⚠ UNVERIFIED matrix path convention — same caveat as
compare_kept_outlier_matrices.py: the exact lvl3 matrix path/filename
convention hasn't been confirmed against the current pipeline's real
disk layout. find_lvl3_matrix_path (imported from
compare_kept_outlier_matrices) tries several candidate patterns and
reports which one worked, or reports what it DID find if nothing
matched the run token — check that output on the first real run.
"""

import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.cm as _mpl_cm
from matplotlib.colors import TwoSlopeNorm
from scipy import stats as _stats

from Exemples.Study_EDNiX.analysis.Archives.Fig5.compare_kept_outlier_matrices import find_lvl3_matrix_path, load_matrix, debug_lvl3_search
from Exemples.Study_EDNiX.analysis.Archives.Fig5.qc_index_kept_outlier_relevance import build_candidate_qc_table, rank_qc_variables

opj = os.path.join
opb = os.path.basename
ope = os.path.exists

_CAT_COLORS = {"Kept": "#009E73", "Outlier": "#D55E00", "all": "#444444"}
_GOLD, _GOLD_EDGE = "#FFD700", "#B8860B"


# =============================================================================
# 1) Cross-species matrix collection + common ROIs
# =============================================================================

def _bilateral_pair_filter(rois):
    """Same logic as ednix_fc_threshold_explorer._bilateral_pair_filter —
    keep only L_/R_ ROI names that have BOTH hemispheres present, paired
    L_ then R_ per base name."""
    l_bases = {n[2:] for n in rois if n.startswith(("L_", "l_"))}
    r_bases = {n[2:] for n in rois if n.startswith(("R_", "r_"))}
    common = sorted(l_bases & r_bases)
    return [f"L_{b}" for b in common] + [f"R_{b}" for b in common]


def collect_all_bids_matrices(all_results, bids_dirs_by_label, level=3,
                              atlas_name="EDNIxCSCLR"):
    """
    For every (label, seed) in all_results (e.g. from
    launch_SBA_MultiConfig.launch()), load every run's lvl<level> matrix
    and split kept/outlier — one entry per BIDS label.

    all_results: {(label, seed): {"df_class": df_class, ...}} — same
    structure launch() returns.
    bids_dirs_by_label: {label: bids_dir} — needed since all_results
    itself doesn't carry the BIDS path.

    Returns {label: {"kept": [(roi_names, ndarray), ...],
                     "outlier": [(roi_names, ndarray), ...]}} — each run
    keeps its OWN roi_names; NO uniformity is assumed or enforced within
    a BIDS. Confirmed real symptom this fixes: individual subjects' lvl3
    matrices can have genuinely different region coverage/ordering even
    within the SAME BIDS (the same per-subject atlas-coverage
    variability already confirmed elsewhere in this pipeline for the
    SBA parcellation) — requiring an EXACT match against the
    first-loaded matrix discarded almost every run for some BIDS (e.g.
    Macaque_BenHamed and Macaque_ZhuGarin both dropped to a group of 0
    usable runs), rather than aligning each matrix individually to
    whatever common ROI set ends up being used downstream (see
    compute_common_rois_across_bids / _align_to_common below).

    A BIDS with zero usable matrices in EITHER group is omitted
    (printed why), since a group comparison needs both.
    """
    out = {}
    for (label, seed), result in all_results.items():
        bids_dir = bids_dirs_by_label.get(label)
        if bids_dir is None:
            print(f"  [WARN] {label}: no bids_dir given — skipped")
            continue
        df_class = result["df_class"]

        kept, outlier = [], []
        n_missing = 0
        for row in df_class.itertuples():
            path = find_lvl3_matrix_path(bids_dir, row.subject, row.session, row.run,
                                         level=level, atlas_name=atlas_name)
            if path is None:
                n_missing += 1
                continue
            roi_names, mat = load_matrix(path)
            (outlier if bool(row.is_outlier) else kept).append((roi_names, mat))

        if not kept or not outlier:
            print(f"  [SKIP] {label}: kept={len(kept)}, outlier={len(outlier)} "
                 f"(need >=1 each) -- {n_missing} run(s) had no matrix")
            if n_missing == len(df_class) and len(df_class) > 0:
                # every single run found nothing at all -- find_lvl3_matrix_path
                # stays silent for that case (only warns when it found candidate
                # files but none matched the token), so run one diagnostic on the
                # first row to show exactly what was searched.
                first_row = next(df_class.itertuples())
                debug_lvl3_search(bids_dir, first_row.subject, first_row.session,
                                  level=level, atlas_name=atlas_name)
            continue

        out[label] = dict(kept=kept, outlier=outlier)
        print(f"  {label}: {len(kept)} kept + {len(outlier)} outlier lvl{level} matrices "
             f"({n_missing} skipped)")
    return out


def compute_common_rois_across_bids(matrices_by_bids):
    """
    Base ROI intersection across EVERY INDIVIDUAL matrix in
    matrices_by_bids (not just one representative set per BIDS, since
    region coverage can vary per run/subject even within the same
    BIDS — see collect_all_bids_matrices), filtered to bilateral
    (L_/R_ both present) pairs — same two-step logic as
    ednix_fc_threshold_explorer.compute_subset_common_rois +
    _bilateral_pair_filter, just operating on already-loaded roi_names
    instead of re-reading files from disk.
    """
    roi_sets = []
    for v in matrices_by_bids.values():
        for group in ("kept", "outlier"):
            for roi_names, _ in v[group]:
                roi_sets.append(set(roi_names))
    if not roi_sets:
        return []
    common = set.intersection(*roi_sets) if len(roi_sets) > 1 else roi_sets[0]
    paired = _bilateral_pair_filter(sorted(common))
    print(f"  {len(paired)} common bilateral ROI(s) across {len(matrices_by_bids)} BIDS "
         f"({len(roi_sets)} individual matrices)")
    return paired


def _align_to_common(mats_with_rois, common_rois):
    """
    Reindex EACH (roi_names, matrix) pair individually onto common_rois
    (same order across every result, needed to stack/compare). A
    matrix that doesn't contain every common ROI is DROPPED
    individually (with a count printed) rather than discarding the
    whole group — this is the per-run/per-subject region-coverage-
    variability fix described in collect_all_bids_matrices.
    """
    if not mats_with_rois:
        return []
    aligned = []
    n_dropped = 0
    for roi_names, mat in mats_with_rois:
        idx_map = [roi_names.index(r) for r in common_rois if r in roi_names]
        if len(idx_map) != len(common_rois):
            n_dropped += 1
            continue
        aligned.append(mat[np.ix_(idx_map, idx_map)])
    if n_dropped:
        print(f"  [INFO] {n_dropped}/{len(mats_with_rois)} matrix(es) dropped -- "
             f"missing one or more common ROIs")
    return aligned


# =============================================================================
# 2) FC matrices by category x BIDS grid (fig04 equivalent)
# =============================================================================

def plot_fc_matrices_kept_outlier_by_bids(matrices_by_bids, common_rois,
                                          output_path=None, figsize=None):
    """
    Same visual design as ednix_fc_threshold_explorer.plot_fc_matrices_by_category
    (rows=category, columns=BIDS, shared RdBu_r colorbar, "n=" annotation
    per cell, common-ROI y-tick labels on the leftmost column) — adapted
    to TWO categories (Kept / Outlier) plus "all" (kept+outlier pooled),
    instead of the original Specific/Unspecific/No/Spurious four.
    """
    labels = sorted(matrices_by_bids.keys())
    if not labels:
        print("  [fig04] no BIDS with usable matrices -- nothing to plot")
        return None
    n_roi = len(common_rois)
    if n_roi < 4:
        print(f"  [fig04] only {n_roi} common ROI(s) -- too few to plot")
        return None

    cats_order = ["all", "Kept", "Outlier"]
    aligned = {cat: {} for cat in cats_order}
    for label in labels:
        v = matrices_by_bids[label]
        kept_a = _align_to_common(v["kept"], common_rois)
        outlier_a = _align_to_common(v["outlier"], common_rois)
        aligned["Kept"][label] = kept_a
        aligned["Outlier"][label] = outlier_a
        aligned["all"][label] = kept_a + outlier_a

    results = {}
    for cat in cats_order:
        results[cat] = {}
        for label in labels:
            mats = aligned[cat][label]
            med_mat = np.nanmedian(np.stack(mats), axis=0) if mats else None
            results[cat][label] = (med_mat, len(mats))

    short = [r.replace("L_", "L ").replace("R_", "R ")[:13] for r in common_rois]

    all_vals = []
    for cat in cats_order:
        for label in labels:
            m, _ = results[cat][label]
            if m is not None:
                all_vals.extend(m.ravel())
    all_vals = np.array(all_vals)
    all_vals = all_vals[np.isfinite(all_vals)]
    vmax_r = max(float(np.nanpercentile(np.abs(all_vals), 97)), 0.05) if len(all_vals) else 0.5
    norm_r = TwoSlopeNorm(vmin=-vmax_r, vcenter=0, vmax=vmax_r)

    n_grp = len(labels)
    n_rows = len(cats_order)
    w, h = figsize or (max(8, n_grp * 3.5), n_rows * 3.5)
    fig = plt.figure(figsize=(w, h))
    gs = gridspec.GridSpec(n_rows, n_grp + 1,
                           width_ratios=[1.0] * n_grp + [0.04],
                           wspace=0.05, hspace=0.35)

    for r_i, cat in enumerate(cats_order):
        col_c = _CAT_COLORS.get(cat, "#444")
        for c_i, label in enumerate(labels):
            ax = fig.add_subplot(gs[r_i, c_i])
            med_mat, n_loaded = results[cat][label]
            if med_mat is None or n_loaded == 0:
                ax.set_facecolor("#f2f2f2")
                ax.set_xticks([]); ax.set_yticks([])
                for sp in ax.spines.values():
                    sp.set_visible(False)
                ax.text(0.5, 0.5, "no data", ha="center", va="center",
                       fontsize=7, color="#999999", transform=ax.transAxes)
            else:
                ax.imshow(med_mat, cmap="RdBu_r", norm=norm_r,
                         aspect="equal", interpolation="nearest")
                ax.text(0.02, 0.98, f"n={n_loaded}", transform=ax.transAxes,
                       fontsize=6, va="top", ha="left", color="#333",
                       backgroundcolor="white", alpha=0.7)
                if c_i == 0:
                    step = max(1, n_roi // 30)
                    ax.set_yticks(range(0, n_roi, step))
                    ax.set_yticklabels(short[::step], fontsize=3)
                else:
                    ax.set_yticks([])
                ax.set_xticks([])
            if c_i == 0:
                ax.set_ylabel(cat, fontsize=8, color=col_c, fontweight="bold", labelpad=4)
            if r_i == 0:
                ax.set_title(str(label), fontsize=8, fontweight="bold")

    cbar_ax = fig.add_subplot(gs[:, -1])
    fig.colorbar(_mpl_cm.ScalarMappable(norm=norm_r, cmap="RdBu_r"),
                cax=cbar_ax, label="Median r")
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


# =============================================================================
# 3) Kept-vs-Outlier "fingerprint" (fig05/fig06 equivalent)
# =============================================================================

def _bh_fdr(pvals):
    """Benjamini-Hochberg FDR — identical to ednix_fc_threshold_explorer._bh_fdr."""
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


def _is_homotopic(a, b):
    def _base(n):
        return re.sub(r"^[LlRr]_", "", str(n)).lower().strip()
    return _base(a) == _base(b)


def build_kept_outlier_fingerprint(matrices_by_bids, common_rois, fdr_alpha=0.05, n_top=10):
    """
    Edges that robustly differ between Kept and Outlier ACROSS EVERY BIDS —
    same edge-test + FDR + top-N-per-type display-selection machinery as
    ednix_fc_threshold_explorer.build_fc_fingerprint, adapted: instead of
    testing whether each edge's per-SPECIES mean differs from 0 (the
    original fingerprint), this tests whether each edge's KEPT-MINUS-
    OUTLIER difference, computed separately PER BIDS, is consistently
    non-zero ACROSS BIDS (one-sample t-test across BIDS-level
    differences) — "this edge changes the same way between kept/outlier
    in every species/BIDS we have, not just by chance in one of them."

    Returns the SAME dict shape as build_fc_fingerprint (rois, n_roi,
    effect_mat, se_mat, pval_mat, qval_mat, homo_mask, cross_mask,
    connections, display, fingerprint_indices, fingerprint_pattern,
    fdr_alpha) so plot_kept_outlier_fingerprint_connections below can
    reuse plot_fingerprint_connections's exact visual grammar.
    """
    labels = sorted(matrices_by_bids.keys())
    n_roi = len(common_rois)
    if n_roi < 4:
        raise ValueError(f"Only {n_roi} common ROI(s) -- need >=4.")

    per_bids_diff = []  # list of (n_roi, n_roi) kept_median - outlier_median matrices
    for label in labels:
        v = matrices_by_bids[label]
        kept_a = _align_to_common(v["kept"], common_rois)
        outlier_a = _align_to_common(v["outlier"], common_rois)
        if not kept_a or not outlier_a:
            continue
        kept_med = np.nanmedian(np.stack(kept_a), axis=0)
        outlier_med = np.nanmedian(np.stack(outlier_a), axis=0)
        per_bids_diff.append(kept_med - outlier_med)

    if len(per_bids_diff) < 2:
        raise ValueError(
            f"Need >=2 BIDS with both kept and outlier data to test "
            f"cross-BIDS consistency (found {len(per_bids_diff)})."
        )

    diff_stack = np.stack(per_bids_diff)  # (n_bids, n_roi, n_roi)

    effect = np.full((n_roi, n_roi), np.nan)
    se_mat = np.full((n_roi, n_roi), np.nan)
    pval = np.full((n_roi, n_roi), np.nan)

    triu = np.triu(np.ones((n_roi, n_roi), bool), k=1)
    for i, j in np.argwhere(triu):
        vals = diff_stack[:, i, j]
        vals = vals[np.isfinite(vals)]
        if len(vals) < 2:
            continue
        e = float(np.mean(vals))
        se = float(np.std(vals, ddof=1) / np.sqrt(len(vals)))
        t, p = _stats.ttest_1samp(vals, 0.0)
        effect[i, j] = effect[j, i] = e
        se_mat[i, j] = se_mat[j, i] = se
        pval[i, j] = pval[j, i] = float(p)

    p_flat = pval[triu]
    q_flat = _bh_fdr(p_flat)
    qval = np.full((n_roi, n_roi), np.nan)
    ii, jj = np.argwhere(triu).T
    qval[ii, jj] = q_flat
    qval[jj, ii] = q_flat

    n_finite_p = int(np.sum(np.isfinite(p_flat)))
    n_survive = int(np.sum(np.isfinite(q_flat) & (q_flat < fdr_alpha)))
    print(f"  [fingerprint] {len(labels)} BIDS x {n_roi} ROIs | "
         f"{n_finite_p} edges tested | {n_survive} survive FDR<{fdr_alpha}")

    homo_mask = np.zeros((n_roi, n_roi), bool)
    cross_mask = np.zeros((n_roi, n_roi), bool)
    for i, ri in enumerate(common_rois):
        for j, rj in enumerate(common_rois):
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
            i=int(i), j=int(j), roi_i=common_rois[i], roi_j=common_rois[j],
            effect=float(effect[i, j]) if np.isfinite(effect[i, j]) else np.nan,
            se=float(se_mat[i, j]) if np.isfinite(se_mat[i, j]) else np.nan,
            pval=float(pval[i, j]) if np.isfinite(pval[i, j]) else np.nan,
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

    return dict(
        rois=common_rois, n_roi=n_roi,
        effect_mat=effect, se_mat=se_mat, pval_mat=pval, qval_mat=qval,
        homo_mask=homo_mask, cross_mask=cross_mask,
        connections=connections, display=display,
        fingerprint_indices=fp_idx, fingerprint_pattern=fp_pattern,
        fp_test="kept_minus_outlier_across_bids", fdr_alpha=fdr_alpha,
    )


# =============================================================================
# 3b) Global permutation test on the fingerprint
# =============================================================================

def _fingerprint_stats_from_assignment(pooled, common_rois, fdr_alpha):
    """
    Given a specific kept/outlier assignment PER BIDS (the observed one,
    or a permuted one), compute the same per-BIDS-diff -> cross-BIDS
    t-test -> FDR pipeline as build_kept_outlier_fingerprint, but return
    only the two summary statistics the permutation test needs:
      - n_survive: edges surviving FDR<fdr_alpha (same count
        build_kept_outlier_fingerprint reports)
      - global_stat: mean(effect^2) across every tested edge -- a
        continuous, threshold-free summary of connectome-wide signal
        strength. Unlike n_survive (which degenerates to ~0 under most
        permutations and gives a coarse null), this stays informative
        even when no single edge individually survives FDR, giving a
        better-powered permutation p-value.

    pooled: {label: (mats_with_roi, kept_index_set)} -- mats_with_roi is
    the FULL list of (roi_names, matrix) tuples for that BIDS (kept+
    outlier pooled, in a fixed order), kept_index_set is the set of
    positions in that list currently labeled "kept" (the true labeling,
    or a permuted one with the same group size).

    Returns None if fewer than 2 BIDS have both groups non-empty under
    this assignment (mirrors build_kept_outlier_fingerprint's own
    ValueError condition, but returns None instead of raising since this
    is called thousands of times during permutation).
    """
    per_bids_diff = []
    for label, (mats, kept_idx) in pooled.items():
        kept_mats = [mats[i] for i in kept_idx]
        outlier_mats = [mats[i] for i in range(len(mats)) if i not in kept_idx]
        kept_a = _align_to_common(kept_mats, common_rois)
        outlier_a = _align_to_common(outlier_mats, common_rois)
        if not kept_a or not outlier_a:
            continue
        kept_med = np.nanmedian(np.stack(kept_a), axis=0)
        outlier_med = np.nanmedian(np.stack(outlier_a), axis=0)
        per_bids_diff.append(kept_med - outlier_med)

    if len(per_bids_diff) < 2:
        return None

    diff_stack = np.stack(per_bids_diff)
    n_roi = len(common_rois)
    triu = np.triu(np.ones((n_roi, n_roi), bool), k=1)
    effect = np.full((n_roi, n_roi), np.nan)
    pval = np.full((n_roi, n_roi), np.nan)
    for i, j in np.argwhere(triu):
        vals = diff_stack[:, i, j]
        vals = vals[np.isfinite(vals)]
        if len(vals) < 2:
            continue
        e = float(np.mean(vals))
        _, p = _stats.ttest_1samp(vals, 0.0)
        effect[i, j] = effect[j, i] = e
        pval[i, j] = pval[j, i] = float(p)

    p_flat = pval[triu]
    q_flat = _bh_fdr(p_flat)
    n_survive = int(np.sum(np.isfinite(q_flat) & (q_flat < fdr_alpha)))

    e_flat = effect[triu]
    finite_e = e_flat[np.isfinite(e_flat)]
    global_stat = float(np.mean(finite_e ** 2)) if finite_e.size else np.nan

    return dict(n_survive=n_survive, global_stat=global_stat)


def permutation_test_fingerprint(matrices_by_bids, common_rois, n_perm=1000,
                                 fdr_alpha=0.05, random_seed=42, verbose=True):
    """
    Global permutation test for the Kept-Outlier fingerprint: is the
    connectome-wide pattern of kept/outlier differences stronger than
    would arise from randomly relabeling kept/outlier WITHIN each BIDS
    (same per-BIDS group sizes preserved as observed)? This is the
    non-circular complement to the per-edge FDR test already in
    build_kept_outlier_fingerprint -- FDR tells you WHICH edges differ;
    this tells you whether the OVERALL pattern is stronger than chance,
    with a single, clean p-value suitable for a methods section.

    Two observed statistics are tested against their null distributions:
      - n_survive: edges surviving FDR<fdr_alpha (matches
        build_kept_outlier_fingerprint's own count)
      - global_stat: mean(effect^2) across every tested edge -- stays
        informative even when the null n_survive distribution collapses
        to mostly 0 (a common outcome with few BIDS/edges), which alone
        would only support "p < 1/n_perm" without further resolution.

    Null distribution: n_perm random within-BIDS relabelings, same
    kept/outlier group SIZES preserved per BIDS as the real data (e.g.
    if a BIDS has 10 kept + 30 outlier, every permutation for that BIDS
    also draws 10 "kept" at random from its 40 runs) -- this tests
    specifically whether the LABELING carries information, not whether
    having any two unequal-sized groups would look different by chance.

    Returns dict(observed_n_survive, observed_global_stat, p_n_survive,
    p_global_stat, null_n_survive, null_global_stat, n_perm) -- the two
    null_* arrays are the full null distributions, for
    plot_permutation_null_distribution below or your own figure.
    """
    rng = np.random.default_rng(random_seed)
    pooled = {label: (v["kept"] + v["outlier"], set(range(len(v["kept"]))))
             for label, v in matrices_by_bids.items()}

    observed = _fingerprint_stats_from_assignment(pooled, common_rois, fdr_alpha)
    if observed is None:
        raise ValueError(
            "Need >=2 BIDS with both kept and outlier data to run this test "
            "(same requirement as build_kept_outlier_fingerprint)."
        )

    null_n_survive = np.zeros(n_perm, dtype=int)
    null_global_stat = np.full(n_perm, np.nan)

    for p_i in range(n_perm):
        perm_pooled = {}
        for label, (mats, kept_idx) in pooled.items():
            n = len(mats)
            n_kept = len(kept_idx)
            perm_kept = set(rng.choice(n, size=n_kept, replace=False).tolist())
            perm_pooled[label] = (mats, perm_kept)
        stats_p = _fingerprint_stats_from_assignment(perm_pooled, common_rois, fdr_alpha)
        if stats_p is None:
            continue
        null_n_survive[p_i] = stats_p["n_survive"]
        null_global_stat[p_i] = stats_p["global_stat"]

    # +1/+1 (add-one) correction: a permutation p-value of exactly 0 is never
    # reportable -- the true null probability is at least 1/(n_perm+1), since
    # the observed assignment IS itself one member of the permutation
    # distribution (Phipson & Smyth, 2010) -- this is the standard correction
    # for exact permutation p-values.
    p_n_survive = (1 + int(np.sum(null_n_survive >= observed["n_survive"]))) / (n_perm + 1)
    valid_null_global = null_global_stat[np.isfinite(null_global_stat)]
    if valid_null_global.size:
        p_global_stat = (1 + int(np.sum(valid_null_global >= observed["global_stat"]))) / \
                        (valid_null_global.size + 1)
    else:
        p_global_stat = np.nan

    if verbose:
        print(f"  [permutation test] {n_perm} within-BIDS relabelings")
        print(f"    n_survive (FDR<{fdr_alpha}):  observed={observed['n_survive']}, "
             f"null mean={null_n_survive.mean():.2f}, p={p_n_survive:.4f}")
        print(f"    global_stat (mean effect^2): observed={observed['global_stat']:.6f}, "
             f"null mean={np.nanmean(null_global_stat):.6f}, p={p_global_stat:.4f}")

    return dict(
        observed_n_survive=observed["n_survive"],
        observed_global_stat=observed["global_stat"],
        p_n_survive=float(p_n_survive),
        p_global_stat=float(p_global_stat),
        null_n_survive=null_n_survive,
        null_global_stat=null_global_stat,
        n_perm=n_perm,
    )


def plot_permutation_null_distribution(perm_result, output_path=None):
    """
    Publication-style figure: null distributions (histograms) from
    permutation_test_fingerprint, with the observed value marked and the
    p-value annotated -- one panel per statistic (global_stat, n_survive).
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    ax = axes[0]
    valid = perm_result["null_global_stat"][np.isfinite(perm_result["null_global_stat"])]
    ax.hist(valid, bins=40, color="#888888", alpha=0.7, edgecolor="white")
    ax.axvline(perm_result["observed_global_stat"], color="#D55E00", lw=2.5,
              label=f"observed (p={perm_result['p_global_stat']:.4f})")
    ax.set_xlabel("mean(effect²) across all edges")
    ax.set_ylabel("permutations")
    ax.set_title("Global connectome-wide signal strength")
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    max_n = int(max(perm_result["null_n_survive"].max(), perm_result["observed_n_survive"]))
    ax.hist(perm_result["null_n_survive"], bins=range(0, max_n + 3),
           color="#888888", alpha=0.7, edgecolor="white", align="left")
    ax.axvline(perm_result["observed_n_survive"], color="#D55E00", lw=2.5,
              label=f"observed (p={perm_result['p_n_survive']:.4f})")
    ax.set_xlabel("edges surviving FDR<0.05")
    ax.set_ylabel("permutations")
    ax.set_title("Number of significant edges")
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.suptitle(f"Permutation test: {perm_result['n_perm']} within-BIDS relabelings", fontsize=11)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_kept_outlier_fingerprint_connections(fp, n_top_per_type=5, output_path=None):
    """
    Same figure as ednix_fc_threshold_explorer.plot_fingerprint_connections
    — build_kept_outlier_fingerprint returns the identical dict shape, so
    this is a direct visual re-implementation with only the axis/title
    text changed (Kept-Outlier difference, not a per-species mean).
    """
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
    ax_heat.set_title(f"Kept − Outlier edge difference (consistent across BIDS) — "
                      f"{n_total_surv} edges survive FDR<{alpha}",
                      fontsize=9, fontweight="bold")
    surv_ij = np.argwhere(np.isfinite(qmat) & (qmat < alpha))
    for i, j in surv_ij:
        if i < j:
            ax_heat.plot(j, i, "*", color=_GOLD, ms=4,
                        markeredgecolor=_GOLD_EDGE, markeredgewidth=0.5, zorder=5)
    cbar_ax_heat = fig.add_subplot(gs[3, 0])
    fig.colorbar(im, cax=cbar_ax_heat, orientation="horizontal",
                label="Kept − Outlier (r)").ax.tick_params(labelsize=7)

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
                sig, lbl = "", f"p={p:.3f} (no q)"
            else:
                sig = "***" if q < .001 else "**" if q < .01 else "*" if q < .05 else ""
                lbl = f"q={q:.3f}{sig}"
            x_off = 0.01 if (np.isfinite(e) and e >= 0) else -0.01
            ax.text(e + x_off if np.isfinite(e) else 0, yi, lbl, va="center",
                   ha="left" if (np.isfinite(e) and e >= 0) else "right",
                   fontsize=6, color="#000" if sig else "#666",
                   fontweight="bold" if sig else "normal")
        ax.set_yticks(y_pos)
        ax.set_yticklabels([f"{c['roi_i'][:15]} ↔ {c['roi_j'][:15]}" for c in conns], fontsize=6)
        ax.axvline(0, color="#aaa", lw=0.7)
        ax.set_xlabel("Kept − Outlier (r)", fontsize=7)
        n_surv_here = sum(1 for c in conns if np.isfinite(c["qval"]) and c["qval"] < alpha)
        ax.set_title(f"{type_labels[ct]}  (top {len(conns)} by |diff|; {n_surv_here} survive FDR)",
                    fontsize=7.5, fontweight="bold", color=color)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=6)

    cbar_ax_r = fig.add_subplot(gs[3, 1])
    sm = plt.cm.ScalarMappable(cmap=cmap_r, norm=norm_r); sm.set_array([])
    fig.colorbar(sm, cax=cbar_ax_r, orientation="horizontal",
                label="Kept − Outlier (r)").ax.tick_params(labelsize=7)
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_kept_outlier_fingerprint_heatmap(matrices_by_bids, fp, output_path=None):
    """
    Fig06 equivalent: for every BIDS, show each fingerprint edge's KEPT
    value and OUTLIER value side by side (two columns per BIDS instead of
    the original's one column per species/BIDS), plus the overall
    Kept-Outlier "pattern" column on the left — same visual grammar
    (shared colorbar, per-BIDS group labels) as
    ednix_fc_threshold_explorer.plot_fingerprint_heatmap.
    """
    display = fp["display"]
    fp_idx = fp["fingerprint_indices"]
    fp_rois = fp["rois"]
    fp_pat = fp["fingerprint_pattern"]
    n_fp = len(display)
    if n_fp == 0:
        print("  [fig06] no surviving connections to display")
        return None

    labels = sorted(matrices_by_bids.keys())
    n_grp = len(labels) * 2  # kept + outlier column per BIDS
    conn_labels = [f"{c['roi_i'][:11]}--{c['roi_j'][:11]}" for c in display]
    conn_types = [c["type"] for c in display]

    group_vals = np.full((n_fp, n_grp), np.nan)
    col_labels = []
    for label in labels:
        v = matrices_by_bids[label]
        kept_a = _align_to_common(v["kept"], fp_rois)
        outlier_a = _align_to_common(v["outlier"], fp_rois)
        kept_med = np.nanmedian(np.stack(kept_a), axis=0) if kept_a else None
        outlier_med = np.nanmedian(np.stack(outlier_a), axis=0) if outlier_a else None
        for arr, tag in ((kept_med, "K"), (outlier_med, "O")):
            col = len(col_labels)
            if arr is not None:
                for k, (fi, fj) in enumerate(fp_idx):
                    group_vals[k, col] = arr[fi, fj]
            col_labels.append(f"{label[:8]}\n{tag}")

    cell_w = max(0.4, 8.0 / n_grp)
    cell_h = max(0.28, 9.0 / n_fp)
    fig_w = float(np.clip(1.8 + n_grp * cell_w + 0.6, 8, 30))
    fig_h = float(np.clip(n_fp * cell_h + 2.0, 6, 20))
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(1, n_grp + 2,
                           width_ratios=[1.2] + [1.0] * n_grp + [0.08],
                           wspace=0.04, top=0.88, bottom=0.20, left=0.18, right=0.97)

    all_v = np.concatenate([fp_pat, group_vals[np.isfinite(group_vals)]])
    vabs = max(float(np.nanpercentile(np.abs(all_v), 97)), 0.05) if all_v.size else 0.1
    norm = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
    cmap = "RdBu_r"

    ax_fp = fig.add_subplot(gs[0, 0])
    ax_fp.imshow(fp_pat.reshape(-1, 1), cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
    ax_fp.set_yticks(range(n_fp))
    ax_fp.set_yticklabels(conn_labels, fontsize=max(3, 6 - n_fp // 10))
    ax_fp.set_xticks([0]); ax_fp.set_xticklabels(["Kept-Outlier"], fontsize=7, rotation=30, ha="right")
    ax_fp.set_title("Pattern", fontsize=8, fontweight="bold")
    ct_col = {"homotopic": "#009E73", "heterotopic": "#56B4E9", "ipsilateral": "#E69F00"}
    for k, ct in enumerate(conn_types):
        ax_fp.axhline(k + 0.5, color=ct_col.get(ct, "#aaa"), lw=0.4, alpha=0.5)

    for gi in range(n_grp):
        ax = fig.add_subplot(gs[0, gi + 1])
        ax.imshow(group_vals[:, gi].reshape(-1, 1), cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
        ax.set_yticks([]); ax.set_xticks([0])
        ax.set_xticklabels([col_labels[gi]], fontsize=max(4, 7 - n_grp // 6), rotation=45, ha="right")

    cbar_ax = fig.add_subplot(gs[0, -1])
    fig.colorbar(_mpl_cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax, label="Median r")
    cbar_ax.yaxis.set_tick_params(labelsize=6)

    fig.legend(handles=[mpatches.Patch(color=ct_col[k], label=k) for k in ct_col],
              loc="lower left", bbox_to_anchor=(0.0, 0.0), ncol=3, fontsize=7, frameon=False)
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


# =============================================================================
# 4) Connectivity + QC combined variable ranking
# =============================================================================

def compute_dual_scores_from_matrix(mat, rois, n_top_homotopic=3):
    """
    Same algorithm as ednix_fc_threshold_explorer.compute_dual_scores_from_matrix
    (primary_intra/inter from top-N homotopic pairs among primary sensory
    ROIs, species_intra from ALL homotopic pairs) — copied here so this
    report script doesn't depend on importing the explorer module
    directly. See that module for the full algorithm rationale;
    unchanged here, trimmed to the three fields this report actually uses.
    """
    _PRIMARY_PATTERNS = ['Somatosensory', 'Auditory_co', 'Visual_stri']

    def _base(name):
        return re.sub(r"^[LlRr]_", "", str(name)).lower().strip()

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

    keep_indices = [i for i, name in enumerate(rois)
                    if any(p.lower() in name.lower() for p in _PRIMARY_PATTERNS)]
    if keep_indices:
        mat_primary = mat[np.ix_(keep_indices, keep_indices)]
        rois_primary = [rois[i] for i in keep_indices]
    else:
        mat_primary, rois_primary = mat, rois

    homo_pairs_primary, l_idx_p, r_idx_p = _get_homotopic_pairs(rois_primary)
    primary_intra_vals = [float(mat_primary[i, j]) for i, j in zip(l_idx_p, r_idx_p)
                          if np.isfinite(mat_primary[i, j])]
    primary_intra = float(np.nanmedian(primary_intra_vals)) if primary_intra_vals else np.nan

    pair_corrs = [(mat_primary[i, j], i, j) for i, j in zip(l_idx_p, r_idx_p)
                 if np.isfinite(mat_primary[i, j])]
    pair_corrs.sort(key=lambda x: x[0], reverse=True)
    n_top = min(n_top_homotopic, len(pair_corrs))
    top_pairs = pair_corrs[:n_top]
    top_l = [i for _, i, _ in top_pairs]
    top_r = [j for _, _, j in top_pairs]

    if top_l and top_r:
        inter_vals = []
        for li in top_l:
            for rj in top_r:
                is_homo = any(li == i and rj == j for _, i, j in top_pairs)
                if not is_homo and np.isfinite(mat_primary[li, rj]):
                    inter_vals.append(float(mat_primary[li, rj]))
        primary_inter = float(np.nanmedian(inter_vals)) if inter_vals else np.nan
    else:
        primary_inter = np.nan

    homo_pairs_all, l_idx_all, r_idx_all = _get_homotopic_pairs(rois)
    species_intra_vals = [float(mat[i, j]) for i, j in zip(l_idx_all, r_idx_all)
                          if np.isfinite(mat[i, j])]
    species_intra = float(np.nanmedian(species_intra_vals)) if species_intra_vals else np.nan

    return dict(primary_intra=primary_intra, primary_inter=primary_inter,
               species_intra=species_intra)


def build_connectivity_table(all_results, bids_dirs_by_label, level=3,
                             atlas_name="EDNIxCSCLR", n_top_homotopic=3):
    """One row per run: label/subject/session/run/seed/is_outlier, plus
    primary_intra/primary_inter/species_intra from that run's OWN
    lvl<level> matrix (not the group median — this is the per-run
    connectivity fingerprint used to test whether it predicts
    kept/outlier)."""
    rows = []
    for (label, seed), result in all_results.items():
        bids_dir = bids_dirs_by_label.get(label)
        if bids_dir is None:
            continue
        for row in result["df_class"].itertuples():
            path = find_lvl3_matrix_path(bids_dir, row.subject, row.session, row.run,
                                         level=level, atlas_name=atlas_name)
            if path is None:
                continue
            roi_names, mat = load_matrix(path)
            scores = compute_dual_scores_from_matrix(mat, roi_names, n_top_homotopic=n_top_homotopic)
            rows.append(dict(label=label, subject=row.subject, session=row.session,
                            run=row.run, seed=row.seed, is_outlier=bool(row.is_outlier),
                            **scores))
    return pd.DataFrame(rows)


def build_combined_variable_table(all_results, bids_dirs_by_label, profiles_by_label=None,
                                  level=3, atlas_name="EDNIxCSCLR", n_top_homotopic=3,
                                  seed_region_hint=None):
    """
    Combine connectivity indices (from lvl<level> matrices) with QC
    indices (from SBA region profiles, if profiles_by_label given) into
    ONE table, one row per run, ready for rank_qc_variables-style
    analysis against is_outlier.

    profiles_by_label: {label: profiles} — the `profiles` dict from
    compute_all_run_profiles for that BIDS, if you want QC indices
    included. None (default) returns connectivity indices only.
    """
    conn_df = build_connectivity_table(all_results, bids_dirs_by_label, level=level,
                                       atlas_name=atlas_name, n_top_homotopic=n_top_homotopic)
    if profiles_by_label is None or conn_df.empty:
        return conn_df

    qc_frames = []
    for (label, seed), result in all_results.items():
        profiles = profiles_by_label.get(label)
        if profiles is None:
            continue
        qc_t = build_candidate_qc_table(result["df_class"], profiles, seed_region_hint=seed_region_hint)
        if not qc_t.empty:
            qc_t["label"] = label
            qc_frames.append(qc_t)
    if not qc_frames:
        return conn_df

    qc_df = pd.concat(qc_frames, ignore_index=True)
    merge_cols = ["label", "subject", "session", "run", "seed"]
    for col in merge_cols:
        conn_df[col] = conn_df[col].astype(str)
        qc_df[col] = qc_df[col].astype(str)

    combined = conn_df.merge(qc_df, on=merge_cols, how="outer", suffixes=("", "_qc"))
    if "is_outlier_qc" in combined.columns:
        combined["is_outlier"] = combined["is_outlier"].fillna(combined["is_outlier_qc"])
        combined = combined.drop(columns=["is_outlier_qc"])
    return combined


def plot_combined_ranking_publication(ranked, connectivity_cols=None, top_n=15, output_path=None):
    """
    Publication-style horizontal bar chart of the combined
    connectivity+QC variable ranking — connectivity-derived variables in
    one color pair, QC-derived variables in another, so a reader
    immediately sees WHICH KIND of variable best explains the
    kept/outlier split.
    """
    if ranked.empty:
        print("  [ranking fig] nothing to plot")
        return None
    connectivity_cols = connectivity_cols or {"primary_intra", "primary_inter", "species_intra"}

    top = ranked.head(top_n).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9, max(4.5, 0.45 * len(top))))
    colors = []
    for _, r in top.iterrows():
        is_conn = r["variable"] in connectivity_cols
        if r["cohens_d"] > 0:
            colors.append("#0072B2" if is_conn else "#009E73")
        else:
            colors.append("#D55E00" if is_conn else "#CC79A7")
    ax.barh(top["variable"], top["cohens_d"].abs(), color=colors, edgecolor="white", lw=0.5)
    ax.set_xlabel("|Cohen's d|  (Kept vs Outlier)", fontsize=10)
    ax.set_title("What best explains the Kept / Outlier split?", fontsize=12, fontweight="bold")
    for i, (d, auc) in enumerate(zip(top["cohens_d"], top["auc"])):
        ax.text(abs(d) + 0.02, i, f"AUC={auc:.2f}", va="center", fontsize=8)
    handles = [
        mpatches.Patch(color="#0072B2", label="Connectivity (kept higher)"),
        mpatches.Patch(color="#D55E00", label="Connectivity (outlier higher)"),
        mpatches.Patch(color="#009E73", label="QC index (kept higher)"),
        mpatches.Patch(color="#CC79A7", label="QC index (outlier higher)"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=7, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


# =============================================================================
# 5) Orchestration
# =============================================================================

def run_cross_species_report(all_results, bids_dirs_by_label, output_dir,
                             profiles_by_label=None, level=3, atlas_name="EDNIxCSCLR",
                             fdr_alpha=0.05, n_top=10, seed_region_hint=None):
    """
    One-call pipeline: fig04 (FC matrices by category x BIDS), fig05/06
    (kept-outlier fingerprint), and the combined connectivity+QC
    variable ranking. Returns a dict with every intermediate result for
    further inspection.

    all_results: {(label, seed): {"df_class": ...}} from
    launch_SBA_MultiConfig.launch().
    bids_dirs_by_label: {label: bids_dir} — e.g.
    {cfg["label"]: cfg["bids_dir"] for cfg in BIDS_CONFIGS}.
    profiles_by_label: optional {label: profiles} to include QC indices
    in the ranking alongside connectivity indices — profiles is
    compute_all_run_profiles's own return value for that BIDS.
    """
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n[1] Loading lvl{level} matrices for every BIDS...")
    matrices_by_bids = collect_all_bids_matrices(all_results, bids_dirs_by_label,
                                                 level=level, atlas_name=atlas_name)
    if not matrices_by_bids:
        raise ValueError("No BIDS produced usable matrices -- nothing to report.")

    print("\n[2] Common ROIs across BIDS...")
    common_rois = compute_common_rois_across_bids(matrices_by_bids)

    print("\n[3] Fig 04 — FC matrices by category x BIDS...")
    plot_fc_matrices_kept_outlier_by_bids(
        matrices_by_bids, common_rois,
        output_path=opj(output_dir, "04_fc_matrices_kept_outlier_by_bids.png"))

    print("\n[4] Fig 05/06 — Kept-Outlier fingerprint...")
    fp = None
    try:
        fp = build_kept_outlier_fingerprint(matrices_by_bids, common_rois,
                                            fdr_alpha=fdr_alpha, n_top=n_top)
        plot_kept_outlier_fingerprint_connections(
            fp, output_path=opj(output_dir, "05_kept_outlier_fingerprint_connections.png"))
        plot_kept_outlier_fingerprint_heatmap(
            matrices_by_bids, fp, output_path=opj(output_dir, "06_kept_outlier_fingerprint_heatmap.png"))
    except ValueError as e:
        print(f"  [skip fingerprint] {e}")

    print("\n[5] Connectivity + QC variable ranking...")
    combined = build_combined_variable_table(
        all_results, bids_dirs_by_label, profiles_by_label=profiles_by_label,
        level=level, atlas_name=atlas_name, seed_region_hint=seed_region_hint)

    ranked = rank_qc_variables(combined, exclude_cols=("subject", "session", "run", "seed", "label"))
    if not ranked.empty:
        excel_path = opj(output_dir, "07_combined_variable_ranking.xlsx")
        with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
            ranked.to_excel(writer, sheet_name="ranking", index=False)
            combined.to_excel(writer, sheet_name="raw_values", index=False)
        print(f"  [saved] {excel_path}")
        plot_combined_ranking_publication(
            ranked, output_path=opj(output_dir, "07_combined_variable_ranking.png"))

    return dict(matrices_by_bids=matrices_by_bids, common_rois=common_rois,
               fingerprint=fp, combined_table=combined, ranked=ranked)


# =============================================================================
# 6) Per-BIDS reference-correlation distributions
# =============================================================================

def plot_reference_correlation_distributions(all_results, percentiles=(5, 25, 50, 75, 95),
                                              min_correlation_by_label=None,
                                              output_path=None, ncols=3):
    """
    For each BIDS/seed in all_results, plot the distribution of
    df_class's own "pearson_r" column — the run's correlation to the
    reference profile, i.e. the exact metric
    classify_runs_by_absolute_correlation used to produce Kept/Outlier —
    as a histogram with the requested percentiles marked. If given, the
    actual min_correlation threshold used for that BIDS is drawn
    separately (red, solid) so you can see where the classification cut
    sits relative to the run-to-run spread, per species.

    min_correlation_by_label: optional {label: threshold}, e.g.
    {cfg["label"]: cfg["min_correlation"] for cfg in launcher.BIDS_CONFIGS}
    if you want the real per-BIDS cutoff drawn (each BIDS_CONFIGS entry
    can use a different min_correlation). None (default) skips the
    threshold line — the percentiles still get drawn either way.
    """
    labels = sorted({label for (label, seed) in all_results.keys()})
    n = len(labels)
    if n == 0:
        print("  [plot] no BIDS in all_results -- nothing to plot")
        return None
    ncols = min(ncols, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4.2, nrows * 3.2), squeeze=False)

    pct_colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(percentiles)))

    for idx, label in enumerate(labels):
        ax = axes[idx // ncols][idx % ncols]
        vals = []
        for (lbl, seed), result in all_results.items():
            if lbl != label:
                continue
            df_class = result["df_class"]
            if "pearson_r" in df_class.columns:
                vals.extend(pd.to_numeric(df_class["pearson_r"], errors="coerce").dropna().tolist())
        vals = np.asarray(vals, dtype=float)
        if vals.size == 0:
            ax.axis("off")
            ax.set_title(f"{label} (no pearson_r data)", fontsize=8, color="#999999")
            continue

        ax.hist(vals, bins=30, color="#4472c4", alpha=0.7, edgecolor="white")
        ymax = ax.get_ylim()[1]
        for p, color in zip(percentiles, pct_colors):
            v = float(np.percentile(vals, p))
            ax.axvline(v, color=color, lw=1.3, ls="--", zorder=3)
            ax.text(v, ymax * 0.97, f"P{p}", rotation=90, fontsize=6,
                   color=color, va="top", ha="right", zorder=4)

        thr = (min_correlation_by_label or {}).get(label)
        if thr is not None:
            ax.axvline(thr, color="#D55E00", lw=2.2, zorder=5)
            n_out = int(np.sum(vals < thr))
            ax.text(thr, ymax * 0.55, f" thr={thr:.2f}\n {n_out}/{len(vals)} outlier",
                   color="#D55E00", fontsize=7, va="center", zorder=6)

        ax.set_title(f"{label} (n={len(vals)})", fontsize=9, fontweight="bold")
        ax.set_xlabel("correlation to reference (pearson_r)", fontsize=7)
        ax.set_ylabel("runs", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig