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

opj = os.path.join
opb = os.path.basename
ope = os.path.exists
spgo = subprocess.getoutput

BOLD_SUFFIX = "_space-template_desc-fMRI_residual.nii.gz"


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


def _format_seed(name):
    for ch in " ()/,:;.-":
        name = name.replace(ch, "_")
    return "".join(c for c in name if c.isalnum() or c == "_").strip("_")


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
    cluster_size_default=10,        # fallback cluster size when ClustSim unavailable
    opening=1,                      # morphological opening for compute_epi_mask
    csim_pthr=0.05,                 # per-voxel p-threshold for ClustSim lookup
    vmax_percentile=99,             # percentile for z-stat colorscale clipping
    method_mask_func='mask_func_over_Gray',  # mask strategy
    output_suffix=None,             # optional suffix appended to all output filenames
):
    """
    Seed-based group t-test using AFNI 3dttest++.

    Parameters
    ----------
    bold_paths       : dict, output of extract_bold_paths()
    label_df         : pd.DataFrame, output of parse_label_file()
    regions_of_interest : list of region_name substrings or None (all regions)
    cluster_size_default : fallback minimum cluster size (voxels) when ClustSim fails
    opening          : opening radius for compute_epi_mask
    csim_pthr        : per-voxel p-threshold for ClustSim table lookup
    vmax_percentile  : percentile [0-100] for z-stat colorscale clipping
    method_mask_func : 'mask_func_over_Gray' | 'mask_func_minus_White' | 'onlyprovidedmask'
    output_suffix    : optional string appended to output file/folder names
    """
    sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _, sing_synstrip, Unetpath = \
        Load_EDNiX_requirement.load_requirement(MAIN_PATH, templatehigh, bids_dir, 'yes')

    output_results1    = opj(bids_dir, "Results")
    os.makedirs(output_results1, exist_ok=True)
    studytemplatebrain = templatehigh if oversample_map else templatelow

    # ?? Group brain mask (shared utility) ????????????????????????????????????
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

    # ?? Filter seeds ??????????????????????????????????????????????????????????
    seeds = label_df.copy()
    if regions_of_interest:
        seeds = seeds[seeds["region_name"].apply(
            lambda b: any(r in _format_seed(b) for r in regions_of_interest))]

    for _, seed_row in seeds.iterrows():
        seed_name     = _format_seed(seed_row["region_name"])
        tag           = f"_{output_suffix}" if output_suffix else ""
        print(f"  Processing seed: {seed_name}")
        output_folder = opj(output_results, seed_name)
        os.makedirs(output_folder, exist_ok=True)

        # ?? Per-subject average Fisher-z maps ?????????????????????????????????
        mean_per_subject = []
        master_file      = None

        for unique_id in sorted(set(all_subjects)):
            fish_maps = [
                _sba_fish_path(bold_p, seed_name)
                for subj, bold_p in zip(all_subjects, all_bold_list)
                if subj == unique_id and ope(_sba_fish_path(bold_p, seed_name))
            ]
            if not fish_maps:
                print(f"  [WARN] No Fisher-z maps: subject={unique_id}  seed={seed_name}")
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

        # ?? Resample all avg maps to master ???????????????????????????????????
        resampled_per_subject = [master_file]
        for avg in mean_per_subject[1:]:
            avg_rs = avg.replace(".nii.gz", "_rs.nii.gz")
            print(spgo(
                f"{sing_afni} 3dresample -master {master_file} -input {avg} "
                f"-prefix {avg_rs} -overwrite -bound_type SLAB"
            ))
            resampled_per_subject.append(avg_rs if ope(avg_rs) else avg)

        # ?? 3dttest++ ?????????????????????????????????????????????????????????
        for f in glob.glob(opj(output_folder, "TTnew*")):
            os.remove(f)
        os.chdir(output_folder)

        use_clustsim = len(resampled_per_subject) > 14
        print(spgo(
            f"{sing_afni} "
            f"3dttest++ -setA {subprocess.list2cmdline(resampled_per_subject)} "
            f"-toz {'-Clustsim ' if use_clustsim else ''}"
            f"-mask {mask_overlap}"
        ))

        # ?? Save AFNI outputs as NIfTI ????????????????????????????????????????
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

        # ?? Cluster threshold ?????????????????????????????????????????????????
        if use_clustsim:
            csim_out = spgo(
                f"{sing_afni} "
                f"1d_tool.py -infile {output_folder}/TTnew.CSimA.NN1_2sided.1D "
                f"-csim_show_clustsize -verb 0 -csim_pthr {csim_pthr} -csim_alpha {alpha}"
            )
            try:
                cluster_size = next(
                    int(l.strip()) for l in csim_out.splitlines() if l.strip().isdigit()
                )
            except StopIteration:
                print(f"  [WARN] ClustSim parse failed, using default={cluster_size_default}")
                cluster_size = cluster_size_default
        else:
            cluster_size = cluster_size_default

        # ?? Threshold & visualise ?????????????????????????????????????????????
        zmap_img  = nib.load(zmap)
        zmap_data = zmap_img.get_fdata()
        vmax      = np.percentile(np.abs(zmap_data), vmax_percentile) if np.any(zmap_data) else 5.0
        z_score   = norm.ppf(1 - alpha / 2)
        print(f"  vmax={vmax:.3f}  alpha={alpha}  z_score={z_score:.3f}")

        stat_mask = np.abs(zmap_data) >= z_score

        # Thresholded z-stat map
        thr_data = zmap_data.copy()
        thr_data[~stat_mask] = 0
        thr_img  = nib.Nifti1Image(thr_data, zmap_img.affine, zmap_img.header)
        thr_path = opj(output_folder, f"{seed_name}{tag}_thresholded_ttest-stat_fisher.nii.gz")
        thr_img.to_filename(thr_path)

        # Correlation map masked by stat threshold
        corr_img      = nib.load(correlation)
        corr_data     = corr_img.get_fdata()
        thr_corr_data = corr_data.copy()
        thr_corr_data[~stat_mask] = 0
        thr_cmap      = nib.Nifti1Image(thr_corr_data, corr_img.affine, corr_img.header)
        thr_cmap_path = opj(output_folder, f"{seed_name}{tag}_thresholded_ttest-stat_fisher_cmap.nii.gz")
        thr_cmap.to_filename(thr_cmap_path)

        vmax_corr = np.max(np.abs(corr_data)) if np.any(corr_data) else 1.0

        # ?? Plot 1: unthresholded correlation map ?????????????????????????????
        display = plotting.plot_stat_map(
            correlation, dim=0, vmax=vmax_corr,
            colorbar=True, bg_img=studytemplatebrain,
            display_mode="y", cut_coords=cut_coords,
            title=f"{seed_name} ? mean Fisher-z (unthresholded)",
        )
        display.savefig(opj(output_folder, f"{seed_name}{tag}_correlation_mosaic.jpg"))
        display.close()

        # ?? Plot 2: thresholded correlation map ???????????????????????????????
        display = plotting.plot_stat_map(
            thr_cmap, dim=0, threshold=1e-6, vmax=vmax_corr,
            colorbar=True, bg_img=studytemplatebrain,
            display_mode="y", cut_coords=cut_coords,
            title=f"{seed_name} ? Fisher-z (p<{alpha}, cluster?{cluster_size})",
        )
        display.savefig(opj(output_folder, f"{seed_name}{tag}_correlation_thresholded_mosaic.jpg"))
        display.close()

        # ?? Plot 3: thresholded z-stat map ????????????????????????????????????
        display = plotting.plot_stat_map(
            thr_img, dim=0, threshold=z_score, vmax=vmax,
            colorbar=True, bg_img=studytemplatebrain,
            display_mode="y", cut_coords=cut_coords,
            title=f"{seed_name} ? t-stat z-map (p<{alpha}, cluster?{cluster_size})",
        )
        display.savefig(opj(output_folder, f"{seed_name}{tag}_thresholded_stat_mosaic.jpg"))
        display.close()
        plt.close("all")