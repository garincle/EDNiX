import nilearn
from nilearn import plotting, image
import glob
import subprocess
import os
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from scipy.stats import norm
from Tools import Load_EDNiX_requirement
from Statistics.Group_fMRI._mask_utils import build_group_mask
from Statistics.Group_fMRI._Group_anal_3dLMEr_SBA import (
    _clustsim_acf, _extract_cluster_threshold, _cluster_mask)

opj = os.path.join
opb = os.path.basename
ope = os.path.exists
spgo = subprocess.getoutput

BOLD_SUFFIX = "_space-template_desc-fMRI_residual.nii.gz"


def _format_seed(name):
    for ch in " ()/,:;.-":
        name = name.replace(ch, "_")
    return "".join(c for c in name if c.isalnum() or c == "_").strip("_")

# =============================================================================
# Helper functions for plotting and normalisation
# =============================================================================
def _bold_root(bold_path):
    bn = opb(bold_path)
    for suffix in (BOLD_SUFFIX, ".nii.gz", ".nii"):
        if bn.endswith(suffix):
            return bn[: -len(suffix)]
    return bn


def _sba_fish_path(bold_path, seed_name):
    return opj(
        os.path.dirname(os.path.dirname(bold_path)),
        "Stats", "SBA", seed_name,
        f"{_bold_root(bold_path)}_correlations_fish.nii.gz",
    )

def _safe_plot(map_path, output_folder, template, title, cmap, cut_coords,
               vmax=None, threshold=None, dpi=300):
    """
    Robustly plot a NIfTI map and save as PNG.
    Returns True if saved, False otherwise.
    """
    if not ope(map_path):
        print(f"  [WARN] {map_path} not found, skipping plot")
        return False
    try:
        img = nib.load(map_path)
        data = img.get_fdata()
        if np.all(data == 0):
            print(f"  [WARN] {map_path} contains only zeros, skipping plot")
            return False
        display = plotting.plot_stat_map(
            map_path, bg_img=template, display_mode="y", cut_coords=cut_coords,
            title=title, cmap=cmap, colorbar=True,
            vmax=vmax, threshold=threshold,
        )
        out_png = opj(output_folder, f"{os.path.basename(map_path).replace('.nii.gz','')}.png")
        display.savefig(out_png, dpi=dpi)
        display.close()
        print(f"  [saved] {out_png}")
        return True
    except Exception as e:
        print(f"  [WARN] Could not save plot for {map_path}: {e}")
        return False

def _normalize_map(data, target_percentile=95):
    """
    Per-map robust scale normalization for cross-species comparability.
    Rescales non-zero voxels so the `target_percentile`-th percentile of
    |data| maps to ±1, then clips to [-1, 1]. Background stays zero.
    """
    nonzero_mask = np.isfinite(data) & (data != 0)
    if not np.any(nonzero_mask):
        return None
    vals = np.abs(data[nonzero_mask])
    scale = np.percentile(vals, target_percentile)
    if scale <= 0:
        return None
    norm_data = np.zeros_like(data)
    norm_data[nonzero_mask] = data[nonzero_mask] / scale
    norm_data[nonzero_mask] = np.clip(norm_data[nonzero_mask], -1, 1)
    return norm_data


# =============================================================================
# Main 3dttest++ function
# =============================================================================

def _3dttest_EDNiX(
    bids_dir,
    templatehigh,
    templatelow,
    oversample_map,
    mask_func,
    cut_coords,
    label_df,
    lower_cutoff,
    upper_cutoff,
    MAIN_PATH,
    alpha,
    bold_paths,
    mean_imgs,
    regions_of_interest=None,
    cluster_size_default=10,
    opening=1,
    csim_pthr=0.001,
    vmax_percentile=99,
    method_mask_func='mask_func_over_Gray',
    output_suffix=None,
    normalize_target=95,
    exclude_lateralized=False,  # NEW: if True, drop any region whose name starts with L_/R_
):
    """
    Seed-based group t-test using AFNI 3dttest++.

    Parameters
    ----------
    bold_paths       : dict, output of extract_bold_paths()
    label_df         : pd.DataFrame, output of parse_label_file()
    regions_of_interest : list of region_name substrings or None (all regions)
    cluster_size_default : kept for backward compatibility, no longer used (no arbitrary fallback)
    opening          : opening radius for compute_epi_mask
    csim_pthr        : per-voxel p-threshold for ClustSim table lookup
    vmax_percentile  : percentile [0-100] for z-stat colorscale clipping
    method_mask_func : 'mask_func_over_Gray' | 'mask_func_minus_White' | 'onlyprovidedmask'
    output_suffix    : optional string appended to output file/folder names
    normalize_target : float or None; if provided, scale median(abs(Fisher‑z)) to this value
                       for an extra normalized non‑thresholded correlation map.
    """
    sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _, sing_synstrip, Unetpath = \
        Load_EDNiX_requirement.load_requirement(MAIN_PATH, templatehigh, bids_dir, 'yes')

    output_results1    = opj(bids_dir, "Results")
    os.makedirs(output_results1, exist_ok=True)
    studytemplatebrain = templatehigh if oversample_map else templatelow

    # Build group mask
    mask_overlap = build_group_mask(
        mean_imgs        = mean_imgs,
        mask_func        = mask_func,
        output_results1  = output_results1,
        sing_afni        = sing_afni,
        studytemplatebrain = studytemplatebrain,
        lower_cutoff     = lower_cutoff,
        upper_cutoff     = upper_cutoff,
        opening          = opening,
        method_mask_func = method_mask_func,
    )

    all_subjects  = bold_paths["subject"]
    all_bold_list = bold_paths["bold_path"]

    folder_name    = "Grp_SBA_3dTTEST" + (f"_{output_suffix}" if output_suffix else "")
    output_results = opj(output_results1, folder_name)
    os.makedirs(output_results, exist_ok=True)

    # Filter seeds
    seeds = label_df.copy()
    if regions_of_interest:
        seeds = seeds[seeds["region_name"].apply(
            lambda b: any(r in _format_seed(b) for r in regions_of_interest))]

    if exclude_lateralized:
        seeds = seeds[~seeds["region_name"].apply(
            lambda b: _format_seed(b).startswith(("L_", "R_")))]

    tag = f"_{output_suffix}" if output_suffix else ""

    for _, seed_row in seeds.iterrows():
        seed_name     = _format_seed(seed_row["region_name"])
        print(f"  Processing seed: {seed_name}")
        output_folder = opj(output_results, seed_name)
        os.makedirs(output_folder, exist_ok=True)

        # ---- 1) Average Fisher‑z maps per subject ----
        mean_per_subject = []
        master_file      = None

        for unique_id in sorted(set(all_subjects)):
            fish_maps = [
                _sba_fish_path(bold_p, seed_name)
                for subj, bold_p in zip(all_subjects, all_bold_list)
                if subj == unique_id and ope(_sba_fish_path(bold_p, seed_name))
            ]
            if not fish_maps:
                print(f"  [WARN] No Fisher‑z maps: subject={unique_id}  seed={seed_name}")
                continue

            master_run     = fish_maps[0]
            resampled_maps = [master_run]
            for fm in fish_maps[1:]:
                fm_rs = fm.replace(".nii.gz", "_rs.nii.gz")
                print(spgo(
                    f"{sing_afni} 3dresample -master {master_run} -input {fm} "
                    f"-prefix {fm_rs} -overwrite -bound_type SLAB"
                ))
                resampled_maps.append(fm_rs if ope(fm_rs) else fm)

            avg_file = opj(output_folder, f"{unique_id}_avg_fisher_map.nii.gz")
            if ope(avg_file):
                os.remove(avg_file)
            print(spgo(
                f"{sing_afni} 3dMean -prefix {avg_file} "
                f"{subprocess.list2cmdline(resampled_maps)}"
            ))
            if ope(avg_file):
                mean_per_subject.append(avg_file)
                if master_file is None:
                    master_file = avg_file

        if len(mean_per_subject) < 2:
            print(f"  [SKIP] {seed_name}: only {len(mean_per_subject)} subject(s)")
            continue

        # ---- 2) Resample all avg maps to master ----
        resampled_per_subject = [master_file]
        for avg in mean_per_subject[1:]:
            avg_rs = avg.replace(".nii.gz", "_rs.nii.gz")
            print(spgo(
                f"{sing_afni} 3dresample -master {master_file} -input {avg} "
                f"-prefix {avg_rs} -overwrite -bound_type SLAB"
            ))
            resampled_per_subject.append(avg_rs if ope(avg_rs) else avg)

        # ---- 3) 3dttest++ ----
        for f in glob.glob(opj(output_folder, "TTnew*")):
            os.remove(f)
        os.chdir(output_folder)

        use_clustsim = len(resampled_per_subject) > 14
        tt_resid = opj(output_folder, "TTresid.nii.gz")
        if ope(tt_resid):
            os.remove(tt_resid)
        print(spgo(
            f"{sing_afni} "
            f"3dttest++ -setA {subprocess.list2cmdline(resampled_per_subject)} "
            f"-toz {'-Clustsim ' if use_clustsim else ''}"
            f"-resid {tt_resid} "
            f"-mask {mask_overlap}"
        ))

        # ---- 4) Save AFNI outputs as NIfTI ----
        zmap        = opj(output_folder, f"{seed_name}{tag}_ttest-stat_fisher_zmap.nii.gz")
        correlation = opj(output_folder, f"{seed_name}{tag}_ttest-stat_fisher_cmap.nii.gz")

        print(spgo(
            f"{sing_afni} 3dcalc -overwrite "
            f"-a {opj(output_folder, 'TTnew+orig.HEAD[1]')} "
            f"-expr a -prefix {zmap}"
        ))
        print(spgo(
            f"{sing_afni} 3dcalc -overwrite "
            f"-a {opj(output_folder, 'TTnew+orig.HEAD[0]')} "
            f"-expr a -prefix {correlation}"
        ))

        # ---- 5) Cluster threshold ----
        # n > 14: 3dttest++ -Clustsim (randomised residuals); otherwise ACF of the
        # 3dttest++ residuals + 3dClustSim. No arbitrary fallback value.
        if use_clustsim:
            cluster_table = opj(output_folder, "TTnew.CSimA.NN1_2sided.1D")
        else:
            cluster_table, acf_abc = _clustsim_acf(
                sing_afni, mask_overlap, tt_resid, opj(output_folder, "Clust_"))
            print(f"  ACF (a, b, c) = {acf_abc}")
        cluster_size = _extract_cluster_threshold(sing_afni, cluster_table, csim_pthr, alpha)

        # ---- 6) Threshold and prepare images ----
        zmap_img  = nib.load(zmap)
        zmap_data = zmap_img.get_fdata()
        vmax_z    = np.percentile(np.abs(zmap_data), vmax_percentile) if np.any(zmap_data) else 5.0
        z_score   = norm.ppf(1 - csim_pthr / 2)
        print(f"  vmax={vmax_z:.3f}  voxel p={csim_pthr}  z={z_score:.3f}  k>={cluster_size:.0f}  FWE alpha={alpha}")

        # Voxelwise p < csim_pthr AND cluster extent >= cluster_size
        stat_mask = _cluster_mask(zmap_data, csim_pthr, cluster_size)

        # Thresholded z-stat map
        thr_data = zmap_data.copy()
        thr_data[~stat_mask] = 0
        thr_img  = nib.Nifti1Image(thr_data, zmap_img.affine, zmap_img.header)
        thr_path = opj(output_folder, f"{seed_name}{tag}_thresholded_ttest-stat_fisher.nii.gz")
        thr_img.to_filename(thr_path)

        # Thresholded correlation map
        corr_img      = nib.load(correlation)
        corr_data     = corr_img.get_fdata()
        thr_corr_data = corr_data.copy()
        thr_corr_data[~stat_mask] = 0
        thr_cmap      = nib.Nifti1Image(thr_corr_data, corr_img.affine, corr_img.header)
        thr_cmap_path = opj(output_folder, f"{seed_name}{tag}_thresholded_ttest-stat_fisher_cmap.nii.gz")
        thr_cmap.to_filename(thr_cmap_path)

        vmax_corr = np.max(np.abs(corr_data)) if np.any(corr_data) else 1.0

        # ---- 7) Plots using safe wrapper ----

        # 7a) Unthresholded correlation map
        _safe_plot(
            correlation, output_folder, studytemplatebrain,
            f"{seed_name} mean Fisher‑z (unthresholded)",
            "cold_hot" if np.any(corr_data < 0) and np.any(corr_data > 0) else "hot",
            cut_coords, vmax=vmax_corr
        )

        # 7b) Thresholded correlation map
        _safe_plot(
            thr_cmap_path, output_folder, studytemplatebrain,
            f"{seed_name} Fisher‑z (p<{csim_pthr}, k≥{cluster_size:.0f}, FWE {alpha})",
            "cold_hot" if np.any(thr_corr_data < 0) and np.any(thr_corr_data > 0) else "hot",
            cut_coords, vmax=vmax_corr, threshold=1e-6
        )

        # 7c) Thresholded z‑stat map
        _safe_plot(
            thr_path, output_folder, studytemplatebrain,
            f"{seed_name} t‑stat z‑map (p<{csim_pthr}, k≥{cluster_size:.0f}, FWE {alpha})",
            "cold_hot" if np.any(thr_data < 0) and np.any(thr_data > 0) else "hot",
            cut_coords, vmax=vmax_z, threshold=z_score
        )

        if normalize_target is not None:
            norm_data = _normalize_map(corr_data, normalize_target)  # or corr_data for ttest
            if norm_data is not None:
                norm_path = opj(output_folder, f"{seed_name}{tag}_zstat_norm.nii.gz")
                nib.save(nib.Nifti1Image(norm_data, corr_img.affine), norm_path)
                nz = norm_data[norm_data != 0]
                cmap_norm = "cold_hot" if (np.any(nz < 0) and np.any(nz > 0)) else "hot"
                vmax_norm = min(1.0, np.percentile(np.abs(nz), vmax_percentile)) if nz.size else 1.0
                _safe_plot(norm_path, output_folder, studytemplatebrain,
                           f"{seed_name} (Normalized, med={normalize_target:.2f})",
                           cmap_norm, cut_coords, vmax=vmax_norm)
        # 7d) Normalized non-thresholded correlation map (for cross-species comparison)
        if normalize_target is not None:
            norm_data = _normalize_map(corr_data, normalize_target)
            if norm_data is not None:
                norm_path = opj(output_folder, f"{seed_name}{tag}_correlation_norm.nii.gz")
                nib.save(nib.Nifti1Image(norm_data, corr_img.affine), norm_path)
                nz = norm_data[norm_data != 0]
                cmap_norm = "cold_hot" if (np.any(nz < 0) and np.any(nz > 0)) else "hot"
                vmax_norm = min(1.0, np.percentile(np.abs(nz), vmax_percentile)) if nz.size else 1.0
                _safe_plot(
                    norm_path, output_folder, studytemplatebrain,
                    f"{seed_name} Normalized Fisher-z corr (p{normalize_target:.0f})",
                    cmap_norm, cut_coords, vmax=vmax_norm
                )

        plt.close("all")