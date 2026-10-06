"""
cluster_based_selection.py
--------------------------------
Unsupervised alternative to human-anchored run selection: instead of
picking runs by similarity to a hand-chosen reference (individual or
group), CLUSTER every run of every BIDS by profile similarity, let the
data itself decide how many groups exist (bounded, so it can't degenerate
to "every run is its own cluster" or "one cluster with everyone"), run
the real group statistics (3dttest++) on each cluster, and only THEN let
the researcher pick -- per BIDS -- which cluster looks biologically
sensible, informed by the cluster's own group map, profile, and jpg
previews.

This keeps human judgment in the loop (still needed to say "cluster 2 is
the clean DMN pattern"), but moves it to a coarser, more defensible level
than picking individual runs or hand-tuning a similarity threshold: you
are choosing among a SMALL number of DATA-DRIVEN, already-averaged
groups, not eyeballing single noisy runs.

Pipeline:
  1. Per-run region profiles (median_r) -- already computed by
     compute_all_run_profiles.
  2. Pairwise correlation between every run's profile within a BIDS
     (compute_pairwise_similarity, from SBA_Region_Profile_Index).
  3. Cluster count auto-selected by silhouette score on the 1-r distance
     matrix (agglomerative, average linkage), bounded to [k_min, k_max]
     (default 2-5) -- NEVER k=1 (no split) or k=n_runs (every run its
     own cluster).
  4. Per BIDS, copy each cluster's SBA preview jpgs into k1/k2/...
     subfolders, and run the real group 3dttest++ (optionally with
     surface reconstruction) on each cluster.
  5. FIGURE 1: every (BIDS, cluster) profile overlaid, color-coded by
     species, on the cross-species common regions -- for visual triage.
  6. FIGURE 2 (after YOU pick, per BIDS, which cluster looks right):
     distance-to-human-reference vs divergence time, exactly like
     human_anchored_phylogeny's own phylogenetic figure, but computed
     from the cluster group maps you chose rather than a
     similarity-selected run subset.

Depends on scikit-learn (AgglomerativeClustering, silhouette_score) --
confirmed available in this environment; if it isn't in yours, install it
or ask for a dependency-free fallback.
"""

import os
from pathlib import Path
import glob
import shutil

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score

from Exemples.Study_EDNiX.analysis.Archives.Fig5.human_anchored_phylogeny import (
    SPECIES_COLORS, restrict_profile,
    _ordered_common_values, )

opj = os.path.join
opb = os.path.basename


# =============================================================================
# 1) Clustering
# =============================================================================

def determine_optimal_k(sim_matrix, k_min=2, k_max=5, min_cluster_size=None,
                        linkage="complete", verbose=True):
    """
    Pick the cluster count k in [k_min, k_max] that maximizes the
    silhouette score of agglomerative clustering on the 1-r distance
    derived from sim_matrix. Silhouette is bounded [-1, 1] and doesn't
    need a "true" k to compare against, unlike inertia/elbow methods --
    it directly answers "how well-separated are the clusters at this k",
    which is exactly what deciding k without a ground truth needs.

    ⚠ Confirmed real failure mode this guards against: silhouette can be
    maximized by carving 1-2 extreme outliers into their OWN tiny
    cluster (trivially "cohesive" since it's near-empty, and far from
    everything else) rather than finding genuinely balanced structure --
    a k that isolates outliers can score HIGHER than the k that actually
    splits the data into real groups. Two changes address this:
      - min_cluster_size: any k where the SMALLEST resulting cluster has
        fewer than this many runs is excluded from consideration
        entirely (not just penalized) -- outlier-isolating solutions
        never even compete for "best k". Defaults to
        max(2, round(0.1 * n)) if not given, i.e. a cluster must hold at
        least ~10% of that BIDS's runs (and never fewer than 2) to count.
      - linkage="complete" (default, was "average"): complete linkage
        merges based on the FARTHEST pair between two clusters, so it
        only merges when EVERY point in both clusters is close -- this
        resists the "chain a couple of outliers off the side" pattern
        average linkage is prone to. Pass linkage="average" to get the
        old behavior if you want to compare.

    k_max is automatically capped at n_runs - 1 (can't have more clusters
    than runs, and k=n_runs would put every run alone, one of the
    degenerate cases this function exists to avoid).

    Returns (best_k, labels, scores) -- labels is an array of cluster
    indices (0-based) the same length/order as sim_matrix's rows; scores
    is {k: silhouette_score} for every k that PASSED the size floor (so
    plot_cluster_diagnostics only shows candidates that were actually
    eligible). If n_runs is too small to support k_min clusters at all,
    returns (1, all-zeros, {}) with a printed reason -- the caller
    decides how to handle a BIDS that can't be meaningfully clustered.
    """
    n = sim_matrix.shape[0]
    if n < k_min + 1:
        if verbose:
            print(f"  [WARN] only {n} run(s) -- too few to cluster into >= {k_min} "
                 f"groups. Returning a single cluster.")
        return 1, np.zeros(n, dtype=int), {}

    if min_cluster_size is None:
        min_cluster_size = max(2, round(0.1 * n))

    k_max_eff = min(k_max, n - 1)
    dist = 1 - np.nan_to_num(sim_matrix, nan=0.0)
    dist = (dist + dist.T) / 2
    np.fill_diagonal(dist, 0.0)
    dist = np.clip(dist, 0, None)  # guard against tiny negative floating-point noise

    best_k, best_score, best_labels = None, -np.inf, None
    scores = {}
    rejected = []
    for k in range(k_min, k_max_eff + 1):
        model = AgglomerativeClustering(n_clusters=k, metric="precomputed", linkage=linkage)
        labels = model.fit_predict(dist)
        if len(set(labels)) < 2:
            continue
        sizes = np.bincount(labels)
        if sizes.min() < min_cluster_size:
            rejected.append((k, list(sizes)))
            continue
        try:
            score = silhouette_score(dist, labels, metric="precomputed")
        except Exception:
            continue
        scores[k] = score
        if score > best_score:
            best_k, best_score, best_labels = k, score, labels

    if best_k is None:
        if verbose:
            print(f"  [WARN] no k in [{k_min}, {k_max_eff}] produced >=2 clusters "
                 f"all with >= {min_cluster_size} runs -- falling back to a single "
                 f"cluster (rejected: {rejected})")
        return 1, np.zeros(n, dtype=int), {}

    if verbose and rejected:
        print(f"  [INFO] rejected for having a cluster smaller than "
             f"min_cluster_size={min_cluster_size}: "
             f"{[(k, sizes) for k, sizes in rejected]}")

    if verbose:
        score_str = ", ".join(f"k={k}:{s:.3f}" for k, s in sorted(scores.items()))
        sizes = np.bincount(best_labels)
        print(f"  Silhouette scores: {score_str}")
        print(f"  Chosen k={best_k}, cluster sizes: {list(sizes)}")
    return best_k, best_labels, scores


def plot_cluster_diagnostics(sim_matrix, key_list, k, labels, scores,
                             title=None, output_path=None, dpi=300):
    """
    Explains WHY k was chosen, in two panels:
      left  -- silhouette score at every k tried, chosen k highlighted.
      right -- the pairwise correlation matrix, rows/columns REORDERED by
               cluster assignment, with white grid lines at cluster
               boundaries -- the block-diagonal pattern (high correlation
               within a block, lower between blocks) IS the evidence the
               clustering acted on; if the blocks aren't visually
               distinct, the chosen k is not capturing much real
               structure, and that's exactly the case this panel is meant
               to expose rather than hide behind a bare k number.

    key_list: run keys in the SAME order as sim_matrix's rows/columns
    (as returned by compute_pairwise_similarity / cluster_bids_runs).
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw=dict(width_ratios=[1, 1.3]))

    ax = axes[0]
    if scores:
        ks = sorted(scores)
        vals = [scores[kk] for kk in ks]
        colors = ["#D55E00" if kk == k else "#4472C4" for kk in ks]
        ax.bar(ks, vals, color=colors, edgecolor="white")
        for kk, v in zip(ks, vals):
            ax.text(kk, v + 0.01, f"{v:.3f}", ha="center", fontsize=8)
    ax.set_xlabel("k (number of clusters)", fontsize=9)
    ax.set_ylabel("silhouette score", fontsize=9)
    ax.set_title(f"Chosen k={k} (highest silhouette)", fontsize=10, fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    ax2 = axes[1]
    order = np.argsort(labels)
    sim_reordered = sim_matrix[np.ix_(order, order)]
    im = ax2.imshow(sim_reordered, cmap="RdBu_r", vmin=-1, vmax=1)
    plt.colorbar(im, ax=ax2, fraction=0.046, label="Pearson r")

    boundaries = np.where(np.diff(labels[order]) != 0)[0] + 0.5
    for b in boundaries:
        ax2.axhline(b, color="white", lw=1.6)
        ax2.axvline(b, color="white", lw=1.6)

    sizes = np.bincount(labels)
    tick_pos, pos = [], 0
    for s in sizes:
        tick_pos.append(pos + s / 2 - 0.5)
        pos += s
    ax2.set_xticks(tick_pos)
    ax2.set_xticklabels([f"k{i}\n(n={s})" for i, s in enumerate(sizes)], fontsize=8)
    ax2.set_yticks([])
    ax2.set_title("Pairwise correlation, reordered by cluster", fontsize=10, fontweight="bold")

    fig.suptitle(title or "Cluster choice diagnostics", fontsize=12, fontweight="bold")
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def cluster_bids_runs(sri, cfg, launcher, k_min=2, k_max=5, min_cluster_size=None,
                      linkage="complete", value_col=None, metric=None,
                      verbose=True):
    """
    Compute profiles for ONE BIDS, cluster them, and return everything
    needed downstream.

    min_cluster_size / linkage: passed straight through to
    determine_optimal_k -- see its docstring for why the defaults changed
    from plain silhouette-only selection (outlier-isolating k's used to
    win).

    Returns dict(profiles=, key_list=, sim_matrix=, k=, labels=, scores=,
    clusters={cluster_id: [key, ...]}) or None if no profiles could be
    computed. `scores` is {k: silhouette_score} for every k tried, for
    plot_cluster_diagnostics.
    """
    value_col = value_col or launcher.VALUE_COL
    metric = metric or launcher.METRIC
    label, seed = cfg["label"], cfg["seeds"][0]

    label_df = sri.parse_itk_snap_label_file(cfg["label_path"])
    profiles, _ = sri.compute_all_run_profiles(
        cfg["species"], cfg["bids_dir"], label_df, launcher.CLASSIFICATION_PATH,
        category=launcher.CATEGORY, seeds=[seed], atlas_level=launcher.ATLAS_LEVEL)
    if not profiles:
        print(f"  [SKIP] {label}: no run profiles computed")
        return None

    sim_matrix, key_list = sri.compute_pairwise_similarity(
        profiles, value_col=value_col, metric=metric)
    if verbose:
        print(f"  {label}: {len(key_list)} run(s), clustering...")
    k, labels, scores = determine_optimal_k(
        sim_matrix, k_min=k_min, k_max=k_max, min_cluster_size=min_cluster_size,
        linkage=linkage, verbose=verbose)

    clusters = {}
    for key, lab in zip(key_list, labels):
        clusters.setdefault(int(lab), []).append(key)

    return dict(profiles=profiles, key_list=key_list, sim_matrix=sim_matrix,
               k=k, labels=labels, scores=scores, clusters=clusters)


# =============================================================================
# 2) Per-cluster jpg copy + group 3dttest
# =============================================================================

def copy_cluster_jpgs(bids_dir, seed, clusters, output_dir, verbose=True):
    """
    Copy each cluster's SBA preview jpg(s) into <output_dir>/k<id>/,
    matching the confirmed real jpg location convention:
        sub-<subj>/ses-<ses>/func/templates/EDNiX/Stats/SBA/<seed>/
            sub-<subj>_ses-<ses>_..._run-<run>_bold_*.jpg

    clusters: {cluster_id: [(subject, session, run, seed), ...]}.
    """
    for cluster_id, keys in clusters.items():
        cluster_dir = opj(output_dir, f"k{cluster_id}")
        os.makedirs(cluster_dir, exist_ok=True)
        n_copied = 0
        for subject, session, run, _seed in keys:
            sba_dir = opj(bids_dir, f"sub-{subject}", f"ses-{session}",
                          "func", "templates", "EDNiX", "Stats", "SBA", str(seed))
            if not os.path.isdir(sba_dir):
                continue
            hits = glob.glob(opj(sba_dir, f"sub-{subject}_ses-{session}_*run-{run}*bold_*.jpg"))
            if not hits:
                hits = glob.glob(opj(sba_dir, f"sub-{subject}_ses-{session}_*bold_*.jpg"))
            for h in sorted(hits):
                shutil.copy2(h, opj(cluster_dir, opb(h)))
                n_copied += 1
        if verbose:
            print(f"    k{cluster_id}: {n_copied} jpg(s) -> {cluster_dir}")


def run_3dttest_per_cluster(sri, cfg, seed, clusters, output_suffix_base,
                            render_surface=False, verbose=True):
    """
    Run the real group 3dttest++ (via rerun_balanced_by_run) SEPARATELY
    for each cluster: cluster members are "kept", everyone else in that
    BIDS is treated as "outlier" for that call (the resulting outlier
    group map is not used -- rerun_balanced_by_run always computes both
    sides, this just reuses that proven code path once per cluster
    rather than writing a new single-group 3dttest function).

    Returns {cluster_id: output_suffix} for the successfully computed
    clusters, so their group maps can be located afterward with
    find_group_map_path(..., output_suffix=that_suffix).
    """
    suffixes = {}
    for cluster_id, keys in clusters.items():
        cluster_set = {tuple(str(v) for v in k) for k in keys}
        rows = []
        for key in [k for c in clusters.values() for k in c]:
            str_key = tuple(str(v) for v in key)
            rows.append(dict(subject=key[0], session=key[1], run=key[2], seed=key[3],
                            is_outlier=str_key not in cluster_set))
        df_class = pd.DataFrame(rows)

        output_suffix = f"{output_suffix_base}_k{cluster_id}"
        try:
            sri.rerun_balanced_by_run(
                df_class, cfg["species"], cfg["bids_dir"], seed=seed,
                output_suffix_kept=output_suffix,
                output_suffix_outlier=f"{output_suffix}_NotInCluster",
                exclude_lateralized=cfg.get("exclude_lateralized", True),
                render_surface=render_surface,
            )
            suffixes[cluster_id] = output_suffix
        except Exception as e:
            print(f"  [WARN] {cfg['label']} k{cluster_id}: group 3dttest failed: {e}")
    return suffixes


# =============================================================================
# 3) Driver: cluster + copy + 3dttest for every BIDS
# =============================================================================

def run_cluster_pipeline(sri, launcher, output_dir, bids_configs=None,
                         k_min=2, k_max=5, min_cluster_size=None, linkage="complete",
                         run_group_analysis=True,
                         render_surface=False, copy_jpgs=True, value_col=None,
                         metric=None):
    """
    Steps 1-4 for every BIDS_CONFIGS entry: cluster, copy jpgs, run
    per-cluster 3dttest.

    min_cluster_size / linkage: see determine_optimal_k's docstring --
    guards against silhouette selecting a k that isolates 1-2 outliers
    into their own tiny cluster instead of finding real structure.

    Returns {label: dict(profiles=, key_list=, sim_matrix=, k=, labels=,
    clusters=, group_suffixes=)} -- group_suffixes is {cluster_id:
    output_suffix} for use with find_group_map_path afterward.
    """
    bids_configs = bids_configs if bids_configs is not None else launcher.BIDS_CONFIGS
    value_col = value_col or launcher.VALUE_COL
    metric = metric or launcher.METRIC
    os.makedirs(output_dir, exist_ok=True)

    results = {}
    for cfg in bids_configs:
        label, seed = cfg["label"], cfg["seeds"][0]
        print(f"\n{'=' * 70}\n  {label} / {seed}\n{'=' * 70}")
        try:
            cluster_result = cluster_bids_runs(
                sri, cfg, launcher, k_min=k_min, k_max=k_max,
                min_cluster_size=min_cluster_size, linkage=linkage,
                value_col=value_col, metric=metric)
            if cluster_result is None:
                continue

            diag_dir = opj(output_dir, "cluster_diagnostics")
            os.makedirs(diag_dir, exist_ok=True)
            plot_cluster_diagnostics(
                cluster_result["sim_matrix"], cluster_result["key_list"],
                cluster_result["k"], cluster_result["labels"], cluster_result["scores"],
                title=f"{label} -- cluster choice",
                output_path=opj(diag_dir, f"{label}_cluster_diagnostics.png"))

            if copy_jpgs:
                jpg_dir = opj(output_dir, "jpgs", label)
                copy_cluster_jpgs(cfg["bids_dir"], seed, cluster_result["clusters"], jpg_dir)

            if run_group_analysis:
                suffixes = run_3dttest_per_cluster(
                    sri, cfg, seed, cluster_result["clusters"],
                    output_suffix_base=str(seed), render_surface=render_surface)
                cluster_result["group_suffixes"] = suffixes
            else:
                cluster_result["group_suffixes"] = {}

            results[label] = cluster_result
        except Exception as e:
            print(f"  [ERROR] {label}: {e}")

    return results


# =============================================================================
# 4) Figure 1: every (BIDS, cluster) profile, color-coded by species
# =============================================================================

def plot_all_cluster_profiles(sri, launcher, cluster_results, common_base_names,
                              value_col=None, use_group_maps=True, title=None,
                              output_path=None, dpi=300):
    """
    FIGURE 1 -- one line per (BIDS, cluster), colored by species, all on
    the cross-species common-region axis. Cluster index shown by line
    style (solid/dashed/dotted/dash-dot, cycling).

    use_group_maps: True (default) plots each cluster's REAL 3dttest
    group map profile (needs run_group_analysis=True to have been run,
    and reads it back via find_group_map_path) -- the statistically
    real result. False plots the cluster's simple profile average
    instead (compute_group_median_profile over the cluster's own runs)
    -- faster, no AFNI dependency, useful for a first look before
    committing to the group analysis.
    """
    value_col = value_col or launcher.VALUE_COL
    linestyles = ["-", "--", ":", "-."]

    fig, ax = plt.subplots(figsize=(max(11, len(common_base_names) * 0.22), 6.5))
    x = np.arange(len(common_base_names))

    cfg_by_label = {c["label"]: c for c in launcher.BIDS_CONFIGS}
    for label, res in sorted(cluster_results.items()):
        cfg = cfg_by_label.get(label)
        color = SPECIES_COLORS.get(cfg["species"] if cfg else None, "#888888")
        for cluster_id, keys in sorted(res["clusters"].items()):
            prof = None
            if use_group_maps and cluster_id in res.get("group_suffixes", {}):
                cfg_ = cfg_by_label[label]
                seed = cfg_["seeds"][0]
                bids_name = Path(cfg_["bids_dir"]).name
                map_path = sri.find_group_map_path(
                    launcher.RESULTS_BASE_DIR, cfg_["species"], bids_name, seed,
                    stat_suffix=launcher.STAT_SUFFIX, category=launcher.CATEGORY,
                    output_suffix=res["group_suffixes"][cluster_id])
                if map_path and os.path.exists(map_path):
                    rep_atlas = sri.pick_representative_atlas_volume(cfg_["bids_dir"])
                    atlas_vol = launcher._find_matching_atlas_volume(
                        cfg_["bids_dir"], map_path, rep_atlas)
                    if atlas_vol:
                        label_df = sri.parse_itk_snap_label_file(cfg_["label_path"])
                        prof = sri.compute_region_profile(
                            map_path, atlas_vol, label_df, atlas_level=launcher.ATLAS_LEVEL)
            if prof is None:
                prof = sri.compute_group_median_profile(res["profiles"], set(keys))

            ls = linestyles[cluster_id % len(linestyles)]
            ax.plot(x, _ordered_common_values(prof, common_base_names, value_col),
                   color=color, ls=ls, lw=1.8, alpha=0.9,
                   label=f"{label} k{cluster_id} (n={len(keys)})")

    ax.axhline(0, color="#999999", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(common_base_names, rotation=90, fontsize=4)
    ax.set_ylabel(value_col, fontsize=10)
    ax.set_title(title or "Every (BIDS, cluster) profile on common regions",
                fontsize=11, fontweight="bold")
    ax.legend(fontsize=6, ncol=3, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


# =============================================================================
# 5) Figure 2: distance-to-human from the researcher-chosen "good" clusters
# =============================================================================

def build_distance_from_chosen_clusters(sri, launcher, cluster_results, chosen,
                                        human_label, common_base_names,
                                        value_col=None, metric=None, verbose=True):
    """
    chosen: {label: cluster_id} -- the cluster YOU judged as the clean/
    biologically-sensible one for each BIDS, including human_label
    itself (its own chosen cluster becomes the human reference).

    Uses each chosen cluster's REAL group map (from run_3dttest_per_cluster)
    -- raises a clear error naming the BIDS if that cluster's group map
    isn't available (group analysis wasn't run, or failed for it).

    Returns distance_df (label, species, distance_to_human, pearson_r)
    exactly like human_anchored_phylogeny's own distance_df, ready for
    test_phylogenetic_signal / plot_distance_vs_divergence.
    """
    value_col = value_col or launcher.VALUE_COL
    metric = metric or launcher.METRIC
    cfg_by_label = {c["label"]: c for c in launcher.BIDS_CONFIGS}

    def _cluster_group_profile(label, cluster_id):
        cfg = cfg_by_label[label]
        res = cluster_results[label]
        seed = cfg["seeds"][0]
        suffix = res.get("group_suffixes", {}).get(cluster_id)
        if suffix is None:
            raise ValueError(f"{label} k{cluster_id}: no group map available "
                            f"(run_group_analysis must be True for this cluster)")
        bids_name = Path(cfg["bids_dir"]).name
        map_path = sri.find_group_map_path(
            launcher.RESULTS_BASE_DIR, cfg["species"], bids_name, seed,
            stat_suffix=launcher.STAT_SUFFIX, category=launcher.CATEGORY,
            output_suffix=suffix)
        if not map_path or not os.path.exists(map_path):
            raise ValueError(f"{label} k{cluster_id}: group map not found at {map_path}")
        rep_atlas = sri.pick_representative_atlas_volume(cfg["bids_dir"])
        atlas_vol = launcher._find_matching_atlas_volume(cfg["bids_dir"], map_path, rep_atlas)
        if atlas_vol is None:
            raise ValueError(f"{label} k{cluster_id}: no atlas volume matching the group map's grid")
        label_df = sri.parse_itk_snap_label_file(cfg["label_path"])
        return sri.compute_region_profile(map_path, atlas_vol, label_df,
                                          atlas_level=launcher.ATLAS_LEVEL)

    if human_label not in chosen:
        raise ValueError(f"{human_label} must have a chosen cluster too (it becomes "
                         f"the human reference).")
    human_prof = restrict_profile(_cluster_group_profile(human_label, chosen[human_label]),
                                  common_base_names)

    rows = [dict(label=human_label, species="Human", distance_to_human=0.0, pearson_r=1.0)]
    for label, cluster_id in chosen.items():
        if label == human_label:
            continue
        cfg = cfg_by_label[label]
        try:
            prof = restrict_profile(_cluster_group_profile(label, cluster_id), common_base_names)
            dist, r, n = sri.profile_distance(prof, human_prof, value_col=value_col, metric=metric)
            rows.append(dict(label=label, species=cfg["species"],
                            distance_to_human=dist, pearson_r=r))
            if verbose:
                print(f"  {label} k{cluster_id}: distance_to_human={dist:.4f} (r={r:.4f})")
        except Exception as e:
            print(f"  [WARN] {label}: {e}")

    return pd.DataFrame(rows)