"""
cluster_based_selection_statmaps.py
--------------------------------------
Voxel-wise variant of cluster_based_selection: instead of reducing each
run to a per-region median_r profile, work directly on the 3D correlation
map and compute similarity between runs from the FULL voxelwise pattern
within a common analysis mask.

Why this can be better:
  - Uses all spatial information. Two runs can share identical region
    medians while looking completely different inside each region
    (focal vs diffuse); voxelwise sees that, medians don't.
  - Atlas-agnostic: no dependence on atlas level choice, per-region
    voxel counts, or region boundary definitions.

Why this can be worse:
  - High-dimensional → small registration errors and edge/CSF voxels
    contribute directly to similarity.
  - No free noise-averaging across voxels within a region (region
    medians do that implicitly).

Both approaches plug into cluster_based_selection's downstream unchanged
(same k selection, same jpg copy, same per-cluster 3dttest, same plots) —
only the run × run similarity matrix comes from a different source.

⚠ Backend: uses nibabel + scipy only. nilearn not required.

⚠ Common-mask policy: intersection of non-zero voxels across every run.
  A voxel only counts if EVERY run has a real value there — the safest
  choice when comparing maps, since including voxels that are zero in
  some runs (AFNI's masked-out fill) deflates similarity for those runs
  without contributing real signal. The intersection/union ratio is
  reported; if it drops below ~5%, that is a real data problem to
  surface, not something to silently paper over with a union mask.
"""

import os
import numpy as np
import nibabel as nib
from nibabel.processing import resample_from_to
from Exemples.Study_EDNiX.analysis.Archives.Fig5.cluster_based_selection import (
    determine_optimal_k,
    plot_cluster_diagnostics,
    copy_cluster_jpgs,
    run_3dttest_per_cluster,
)

opj = os.path.join
opb = os.path.basename
ope = os.path.exists


# =============================================================================
# 1) Locate the statmap file for each run key
# =============================================================================

def _find_statmap_for_key(sri, bids_dir, subject, session, run, seed,
                           use_fisher=True):
    """
    Given a (subject, session, run, seed) key, return the path to that
    run's correlation map. Prefers the Fisher-z variant
    (_correlations_fish.nii.gz) when use_fisher=True, since Fisher-z is
    the right space for Pearson correlation ACROSS voxels (its values
    are closer to Gaussian and its variance is roughly stabilized,
    unlike raw r near ±1). Falls back to the raw map if the Fisher
    variant doesn't exist on disk.

    Uses sri.find_raw_corr_maps_for_run + sri._extract_run_id — the SAME
    two functions compute_all_run_profiles itself uses, so a key that
    survived that call will match here by construction.
    """
    files = sri.find_raw_corr_maps_for_run(bids_dir, subject, session, seed)
    for f in files:
        if sri._extract_run_id(opb(f)) == str(run):
            if use_fisher:
                fish = f.replace("_correlations.nii.gz", "_correlations_fish.nii.gz")
                if ope(fish):
                    return fish
            return f
    return None


# =============================================================================
# 2) Load + align + mask + vectorize + similarity
# =============================================================================


def _intersection_mask(stack, verbose=True, min_intersection_frac=0.05):
    """
    Build a mask as the intersection of NON-ZERO voxels across all runs
    in `stack`. AFNI writes literal 0.0 for voxels outside the analysis
    mask used by 3dttest++ — including those is not "real data", it's
    fill, and it would deflate similarity between a run with a small
    coverage and a run with a large one.

    Reports the intersection/union ratio; if it drops below
    min_intersection_frac, prints a WARN naming the specific shortfall
    rather than silently falling back to a permissive union.
    """
    n = stack.shape[0]
    inter = np.ones(stack.shape[1:], dtype=bool)
    union = np.zeros(stack.shape[1:], dtype=bool)
    for i in range(n):
        nz = stack[i] != 0
        inter &= nz
        union |= nz
    n_inter = int(inter.sum())
    n_union = int(union.sum())
    frac = n_inter / max(n_union, 1)
    if verbose:
        print(f"  Mask: intersection={n_inter} vox, union={n_union} vox, "
             f"ratio={frac:.3f}")
    if frac < min_intersection_frac:
        print(f"  [WARN] intersection/union = {frac:.3f} < "
             f"{min_intersection_frac:.2f} — coverage is very heterogeneous "
             f"across runs. Clustering on the intersection may be dominated "
             f"by whichever subjects happened to share the smallest mask; "
             f"consider comparing with a union or a group mask instead.")
    return inter

def _squeeze_trailing_singleton(arr):
    """
    Drop any trailing singleton dimensions — nibabel/AFNI sometimes writes
    (X, Y, Z, 1) for a scalar volume. Every other load in this codebase
    squeezes this for the same reason (qc._squeeze_trailing_singleton);
    this local copy avoids adding an import dependency.
    """
    while arr.ndim > 3 and arr.shape[-1] == 1:
        arr = arr[..., 0]
    return arr


def _load_aligned(path, master_shape, master_affine, verbose=False):
    """
    Load a NIfTI, squeeze any trailing singleton, and resample onto the
    master grid ONLY if shape or affine actually differ after squeezing.

    ⚠ Confirmed real bug fixed here (Macaque_BenHamed, subject 'Elak'):
    the previous version tested `img.shape[:3] == master_shape` for the
    "already matches" case — which SLICES OFF a trailing singleton and so
    reports a (51, 62, 36, 1) volume as matching a (51, 62, 36) grid. It
    then returned the UNSQUEEZED 4D array, the caller's own shape check
    caught the mismatch, and it raised a "resample=True would fix this"
    error that was both wrong (resample was never wired up) and
    misleading (no resample was needed — just a squeeze).
    """
    img = nib.load(path)
    data = _squeeze_trailing_singleton(np.asarray(img.dataobj, dtype=np.float32))

    grid_matches = (data.shape == master_shape
                    and np.allclose(img.affine, master_affine, atol=1e-3))
    if grid_matches:
        return data

    # Real grid mismatch (shape or affine) → resample. Log it, since a
    # resample firing means the source map is genuinely off-grid and may
    # deserve its own investigation (it usually shouldn't happen within
    # one BIDS session).
    if verbose:
        print(f"  [resample] {opb(path)}: {data.shape} -> {master_shape}")
    img3d = nib.Nifti1Image(data, img.affine, img.header)
    resampled = resample_from_to(img3d, (master_shape, master_affine))
    out = np.asarray(resampled.dataobj, dtype=np.float32)
    return _squeeze_trailing_singleton(out)

def compute_statmap_similarity_matrix(sri, bids_dir, keys, seed,
                                       mask_path=None, use_fisher=True,
                                       metric="correlation",
                                       resample=True, verbose=True):
    """
    Build a run × run similarity matrix from the actual 3D statmaps.

    Parameters
    ----------
    sri : module reference — the SBA_Region_Profile_Index module, needed
          for find_raw_corr_maps_for_run / _extract_run_id.
    bids_dir : the BIDS root, to locate the maps.
    keys : every (subject, session, run, seed) tuple to include — the
           keys of the profiles dict from compute_all_run_profiles.
    seed : the seed folder under Stats/SBA/ to look in.
    mask_path : optional path to an explicit mask NIfTI. None (default)
           uses the intersection-of-nonzero mask (see _intersection_mask).
    use_fisher : True (default) — prefer _correlations_fish.nii.gz.
    metric : "correlation" (Pearson across voxels) or "cosine".
    resample : True — resample every map onto the FIRST run's grid if the
           grids differ. False (default) resamples only if necessary AND
           reports every resample. Set True to always force the first
           run's grid (useful when grids differ by tiny affine noise).

    Returns (sim_matrix, valid_keys, info) where:
      - sim_matrix : (n_valid, n_valid) ndarray, NaN-safe (diagonal 1.0).
      - valid_keys : sorted list of keys actually included (a key whose
                     statmap couldn't be located is dropped, with a warn).
      - info       : dict(mask, shape, n_voxels, n_dropped, mask_path).
    """
    # 1. Locate a statmap file for every key
    map_paths = {}
    dropped = []
    for key in keys:
        subj, ses, run, _ = key
        path = _find_statmap_for_key(sri, bids_dir, subj, ses, run, seed,
                                      use_fisher=use_fisher)
        if path is not None:
            map_paths[key] = path
        else:
            dropped.append(key)
            if verbose:
                print(f"  [WARN] no statmap for {key}")

    if len(map_paths) < 3:
        raise ValueError(f"Too few statmaps located ({len(map_paths)}) to cluster.")

    valid_keys = sorted(map_paths.keys(), key=lambda k: tuple(str(x) for x in k))
    master_img = nib.load(map_paths[valid_keys[0]])
    master_shape = master_img.shape[:3]
    master_affine = master_img.affine

    # 2. Load + align into one stack
    stack = np.empty((len(valid_keys),) + master_shape, dtype=np.float32)
    for i, key in enumerate(valid_keys):
        data = _load_aligned(map_paths[key], master_shape, master_affine,
                             verbose=(verbose and i > 0))
        data = _squeeze_trailing_singleton(data)  # <-- ADD
        if data.shape != master_shape:
            raise ValueError(
                f"{key}: shape {data.shape} after loading/squeezing does not match "
                f"master {master_shape} — this is a real grid mismatch, not a "
                f"trailing-singleton artifact.")
        stack[i] = data

    # 3. Build the common mask
    if mask_path is not None:
        mask_data = _load_aligned(mask_path, master_shape, master_affine)
        mask = mask_data != 0
        if verbose:
            print(f"  Mask from {opb(mask_path)}: {int(mask.sum())} vox")
    else:
        mask = _intersection_mask(stack, verbose=verbose)

    if mask.sum() < 10:
        raise ValueError(f"Common mask has only {int(mask.sum())} voxel(s) — "
                         f"cannot compute a meaningful similarity.")

    # 4. Vectorize
    mat = stack[:, mask]  # (n_runs, n_voxels)

    # 5. Similarity
    if metric == "correlation":
        # Pearson across voxels — center each row before normalizing
        mat_c = mat - mat.mean(axis=1, keepdims=True)
        norms = np.linalg.norm(mat_c, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        mat_c = mat_c / norms
        sim = mat_c @ mat_c.T
    elif metric == "cosine":
        norms = np.linalg.norm(mat, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        mat_n = mat / norms
        sim = mat_n @ mat_n.T
    else:
        raise ValueError(f"Unknown metric {metric!r}")

    sim = np.clip(sim, -1.0, 1.0)
    np.fill_diagonal(sim, 1.0)

    if verbose:
        off = sim[~np.eye(len(valid_keys), dtype=bool)]
        print(f"  Statmap similarity: n={len(valid_keys)} runs, "
             f"mean off-diag r={off.mean():.4f}, "
             f"median={np.median(off):.4f}, min={off.min():.4f}, max={off.max():.4f}")

    info = dict(mask=mask, shape=master_shape, n_voxels=int(mask.sum()),
                n_dropped=len(dropped), dropped=dropped, mask_path=mask_path)
    return sim, valid_keys, info


# =============================================================================
# 3) Cluster ONE BIDS from statmaps — same shape as cluster_bids_runs
# =============================================================================

def cluster_bids_runs_from_statmaps(sri, cfg, launcher, k_min=2, k_max=5,
                                     min_cluster_size=None, linkage="complete",
                                     mask_path=None, use_fisher=True,
                                     metric="correlation",
                                     value_col=None, verbose=True):
    """
    Drop-in alternative to cluster_based_selection.cluster_bids_runs:
    same return shape, same downstream expectations — the ONLY change is
    that the similarity matrix comes from the full 3D statmaps instead of
    the region-median profiles.

    Region profiles are STILL computed (they're needed by the plotting
    helpers in cluster_based_selection, e.g. compute_group_median_profile
    for use_group_maps=False), but they are no longer what drives the
    clustering. Reusing those helpers keeps the two methods directly
    comparable: same cluster count logic, same figure code, only the
    similarity source differs.

    Returns dict(profiles=, key_list=, sim_matrix=, k=, labels=, scores=,
    clusters=, mask_info=), or None if no profiles could be computed.
    """
    value_col = value_col or launcher.VALUE_COL
    label, seed = cfg["label"], cfg["seeds"][0]

    # Profiles are computed only so downstream plotting has them available.
    # The clustering itself does not use them.
    label_df = sri.parse_itk_snap_label_file(cfg["label_path"])
    profiles, _ = sri.compute_all_run_profiles(
        cfg["species"], cfg["bids_dir"], label_df, launcher.CLASSIFICATION_PATH,
        category=launcher.CATEGORY, seeds=[seed], atlas_level=launcher.ATLAS_LEVEL)
    if not profiles:
        print(f"  [SKIP] {label}: no run profiles computed")
        return None

    if verbose:
        print(f"  {label}: {len(profiles)} run(s), computing voxel-wise similarity...")
    sim_matrix, key_list, mask_info = compute_statmap_similarity_from_profiles(
        sri, cfg, profiles, mask_path=mask_path, use_fisher=use_fisher,
        metric=metric, verbose=verbose)

    # Restrict profiles to the keys that survived statmap location (a key
    # whose map file is missing must not silently participate downstream).
    profiles = {k: profiles[k] for k in key_list}

    k, labels, scores = determine_optimal_k(
        sim_matrix, k_min=k_min, k_max=k_max, min_cluster_size=min_cluster_size,
        linkage=linkage, verbose=verbose)

    clusters = {}
    for key, lab in zip(key_list, labels):
        clusters.setdefault(int(lab), []).append(key)

    return dict(profiles=profiles, key_list=key_list, sim_matrix=sim_matrix,
                k=k, labels=labels, scores=scores, clusters=clusters,
                mask_info=mask_info)


def compute_statmap_similarity_from_profiles(sri, cfg, profiles, mask_path=None,
                                              use_fisher=True, metric="correlation",
                                              verbose=True):
    """Convenience wrapper: turn a profiles dict into its key list and hand
    off to compute_statmap_similarity_matrix with the right seed."""
    keys = list(profiles.keys())
    seed = cfg["seeds"][0]
    return compute_statmap_similarity_matrix(
        sri, cfg["bids_dir"], keys, seed,
        mask_path=mask_path, use_fisher=use_fisher, metric=metric,
        verbose=verbose)


# =============================================================================
# 4) Driver — same shape as cluster_based_selection.run_cluster_pipeline
# =============================================================================

def run_cluster_pipeline_from_statmaps(sri, launcher, output_dir,
                                        bids_configs=None,
                                        k_min=2, k_max=5, min_cluster_size=None,
                                        linkage="complete",
                                        use_fisher=True, metric="correlation",
                                        mask_path=None,
                                        run_group_analysis=True,
                                        render_surface=False, copy_jpgs=True,
                                        value_col=None):
    """
    Same flow as cluster_based_selection.run_cluster_pipeline, using the
    statmap-based cluster function. Identical downstream steps (k
    selection diagnostics, jpg copy, per-cluster 3dttest) — the only
    difference is what the run × run similarity matrix is built from.

    Outputs land under <output_dir>/cluster_diagnostics and
    <output_dir>/jpgs/<label>/k<N>/ exactly like the profile version, so
    running both into DIFFERENT output_dir values lets you diff them
    directly.

    Returns {label: cluster_result_dict} — same keys as
    run_cluster_pipeline, plus mask_info.
    """
    bids_configs = bids_configs if bids_configs is not None else launcher.BIDS_CONFIGS
    os.makedirs(output_dir, exist_ok=True)

    results = {}
    for cfg in bids_configs:
        label, seed = cfg["label"], cfg["seeds"][0]
        print(f"\n{'=' * 70}\n  {label} / {seed}  [statmaps]\n{'=' * 70}")
        try:
            cluster_result = cluster_bids_runs_from_statmaps(
                sri, cfg, launcher, k_min=k_min, k_max=k_max,
                min_cluster_size=min_cluster_size, linkage=linkage,
                mask_path=mask_path, use_fisher=use_fisher, metric=metric,
                value_col=value_col)
            if cluster_result is None:
                continue

            diag_dir = opj(output_dir, "cluster_diagnostics")
            os.makedirs(diag_dir, exist_ok=True)
            plot_cluster_diagnostics(
                cluster_result["sim_matrix"], cluster_result["key_list"],
                cluster_result["k"], cluster_result["labels"],
                cluster_result["scores"],
                title=f"{label} -- cluster choice [statmaps]",
                output_path=opj(diag_dir, f"{label}_cluster_diagnostics.png"))

            if copy_jpgs:
                jpg_dir = opj(output_dir, "jpgs", label)
                copy_cluster_jpgs(cfg["bids_dir"], seed,
                                  cluster_result["clusters"], jpg_dir)

            if run_group_analysis:
                suffixes = run_3dttest_per_cluster(
                    sri, cfg, seed, cluster_result["clusters"],
                    output_suffix_base=f"{seed}_statmap",
                    render_surface=render_surface)
                cluster_result["group_suffixes"] = suffixes
            else:
                cluster_result["group_suffixes"] = {}

            results[label] = cluster_result
        except Exception as e:
            print(f"  [ERROR] {label}: {e}")

    return results