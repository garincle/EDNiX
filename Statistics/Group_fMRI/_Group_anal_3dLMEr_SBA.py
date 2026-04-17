import nilearn
from nilearn import plotting
from Tools import Load_EDNiX_requirement
from Statistics.Group_fMRI._mask_utils import build_group_mask
import subprocess
import os
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt
from scipy.stats import norm

opj = os.path.join
opb = os.path.basename
ope = os.path.exists
spco = subprocess.check_output
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


def _extract_cluster_threshold(cluster_file, pthr, alpha):
    try:
        with open(cluster_file) as f:
            lines = f.readlines()
        alpha_values = [float(a) for a in lines[6].split("|")[1].split()]
        if alpha not in alpha_values:
            return 10
        col_idx  = alpha_values.index(alpha) + 1
        data     = np.loadtxt(cluster_file, comments="#")
        row_idx  = np.where(np.isclose(data[:, 0], pthr, atol=1e-6))[0]
        if not len(row_idx):
            return 10
        cluster_size = data[row_idx[0], col_idx]
        print(f"  Cluster threshold pthr={pthr} alpha={alpha} ? {cluster_size:.1f} voxels")
        return cluster_size
    except Exception as e:
        print(f"  [WARN] cluster file unreadable: {e}; using 10")
        return 10


def _3dLMEr_EDNiX(
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
    model,
    gltCode,
    stat_img,
    regions_of_interest=None,
    opening=1,                      # morphological opening for compute_epi_mask
    csim_pthr=0.05,                 # per-voxel p-threshold for ClustSim lookup
    vmax_percentile=99,             # percentile for colorscale clipping
    method_mask_func='mask_func_over_Gray',
    output_suffix=None,             # optional suffix for output names
):
    """
    Seed-based group analysis using AFNI 3dLMEr (mixed-effects).

    Parameters
    ----------
    model      : 3dLMEr model string, e.g. '"Sess*(1|Subj)"'
    gltCode    : full -gltCode argument string(s) passed verbatim to 3dLMEr
    stat_img   : list of contrast labels matching gltCode order
    opening    : opening radius for compute_epi_mask
    csim_pthr  : per-voxel p-threshold for ClustSim table lookup
    vmax_percentile : percentile [0-100] for colorscale clipping
    method_mask_func : mask strategy (see _mask_utils.build_group_mask)
    output_suffix    : optional string appended to output file/folder names
    """
    sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _, sing_synstrip, Unetpath = \
        Load_EDNiX_requirement.load_requirement(MAIN_PATH, templatehigh, bids_dir, 'yes')

    output_results1    = opj(bids_dir, "Results")
    os.makedirs(output_results1, exist_ok=True)
    studytemplatebrain = templatehigh if oversample_map else templatelow

    # ?? Group brain mask ??????????????????????????????????????????????????????
    mask_overlap = build_group_mask(
        mean_imgs          = mean_imgs,
        mask_func          = mask_func,
        output_results1    = output_results1,
        sing_afni          = sing_afni,
        studytemplatebrain = studytemplatebrain,
        lower_cutoff       = lower_cutoff,
        upper_cutoff       = upper_cutoff,
        opening            = opening,
        method_mask_func   = method_mask_func,
    )

    all_subjects  = bold_paths["subject"]
    all_sessions  = bold_paths["session"]
    all_runs      = bold_paths["run"]
    all_bold_list = bold_paths["bold_path"]

    folder_name    = "Grp_SBA_3dLME_network" + (f"_{output_suffix}" if output_suffix else "")
    output_results = opj(output_results1, folder_name)
    os.makedirs(output_results, exist_ok=True)

    seeds = label_df.copy()
    if regions_of_interest:
        seeds = seeds[seeds["region_name"].apply(
            lambda b: any(r in _format_seed(b) for r in regions_of_interest))]

    for _, seed_row in seeds.iterrows():
        seed_name     = _format_seed(seed_row["region_name"])
        tag           = f"_{output_suffix}" if output_suffix else ""
        output_folder = opj(output_results, seed_name)
        os.makedirs(output_folder, exist_ok=True)

        # ?? Collect Fisher-z maps, build design matrix ????????????????????????
        all_images, all_IDs, all_Sess, all_Runs = [], [], [], []
        Resample_master = None

        for subj, ses, run, bold_p in zip(all_subjects, all_sessions, all_runs, all_bold_list):
            fish = _sba_fish_path(bold_p, seed_name)
            if not ope(fish):
                print(f"  [SKIP] {fish}")
                continue

            if Resample_master is None:
                Resample_master = fish
                final_fish = fish
            else:
                check = subprocess.run(
                    f"{sing_afni} 3dMatch -quiet -source {Resample_master} -input {fish}",
                    shell=True, capture_output=True, text=True,
                )
                if check.returncode != 0:
                    resampled = fish.replace(".nii.gz", "_resampled.nii.gz")
                    subprocess.run(
                        f"{sing_afni} 3dresample -overwrite -master {Resample_master} "
                        f"-input {fish} -prefix {resampled}",
                        shell=True, check=True,
                    )
                    final_fish = resampled
                else:
                    final_fish = fish

            all_images.append(final_fish)
            all_IDs.append(str(subj))
            all_Sess.append(f"Sess_{ses}")
            all_Runs.append(f"run_{run}" if run is not None else "run_01")

        if len(all_images) < 2:
            print(f"  [SKIP] {seed_name}: only {len(all_images)} image(s)")
            continue

        design_df = (
            pd.DataFrame({"Subj": all_IDs, "Sess": all_Sess, "run": all_Runs, "InputFile": all_images})
            .drop_duplicates().sort_values(by=["Subj", "run"])
        )
        dm_path = opj(output_folder, "design_matrix.txt")
        if ope(dm_path):
            os.remove(dm_path)
        design_df.to_csv(dm_path, index=False, sep="\t")

        # ?? 3dLMEr ????????????????????????????????????????????????????????????
        stat_maps  = opj(output_folder, f"{seed_name}{tag}_3dLME_glt.nii.gz")
        resid_path = opj(output_folder, f"{seed_name}{tag}_resid.nii.gz")
        for f in [stat_maps, resid_path]:
            if ope(f):
                os.remove(f)

        os.chdir(output_results)
        lmer_cmd = (
            f"{sing_afni} "
            f"3dLMEr -prefix {stat_maps} -jobs 20 -mask {mask_overlap} "
            f"-model {model} {gltCode} "
            f"-dataTable @{dm_path} -resid {resid_path}"
        )
        print(lmer_cmd)
        spco(lmer_cmd, shell=True)

        spco(
            f"{sing_afni} 3dClustSim -mask {mask_overlap} -LOTS "
            f"-prefix {output_folder}/Clust_",
            shell=True,
        )

        # ?? Sub-brick extraction ? NIfTI + threshold + visualise ??????????????
        z_score   = norm.ppf(1 - alpha / 2)
        full_img  = nib.load(stat_maps)
        full_data = full_img.get_fdata()
        affine    = full_img.affine

        for i, glt_label in enumerate(stat_img):
            coef_brick  = 1 + 2 * i
            zstat_brick = 2 + 2 * i

            coef_data  = full_data[..., coef_brick]
            zstat_data = full_data[..., zstat_brick]

            coef_path  = opj(output_folder, f"{seed_name}{tag}_{glt_label}_coef.nii.gz")
            zstat_path = opj(output_folder, f"{seed_name}{tag}_{glt_label}_zstat.nii.gz")
            nib.save(nib.Nifti1Image(coef_data,  affine), coef_path)
            nib.save(nib.Nifti1Image(zstat_data, affine), zstat_path)

            vmax_coef  = np.percentile(np.abs(coef_data),  vmax_percentile) if np.any(coef_data)  else 1.0
            vmax_zstat = np.percentile(np.abs(zstat_data), vmax_percentile) if np.any(zstat_data) else 5.0

            # Plot 1: unthresholded coefficient map
            display = plotting.plot_stat_map(
                coef_path, dim=0, vmax=vmax_coef,
                colorbar=True, bg_img=studytemplatebrain,
                display_mode="y", cut_coords=cut_coords,
                title=f"{seed_name} ? {glt_label} coef (unthresholded)",
            )
            display.savefig(opj(output_folder, f"{seed_name}{tag}_{glt_label}_correlation.jpg"))
            display.close()
            plt.close("all")

            if glt_label == "groupeffect":
                # No cluster thresholding for group mean ? just plot the Z-stat
                display = plotting.plot_stat_map(
                    zstat_path, dim=0, vmax=vmax_zstat,
                    colorbar=True, bg_img=studytemplatebrain,
                    display_mode="y", cut_coords=cut_coords,
                    title=f"{seed_name} ? {glt_label} Z-stat (unthresholded)",
                )
                display.savefig(opj(output_folder, f"{seed_name}{tag}_{glt_label}_zstat.jpg"))
                display.close()
                plt.close("all")
            else:
                cluster_size = _extract_cluster_threshold(
                    opj(output_folder, "Clust_.NN1_2sided.1D"), csim_pthr, alpha
                )

                # Build stat mask from z-stat map
                stat_mask = np.abs(zstat_data) >= z_score

                # Thresholded z-stat
                thr_zstat_data = zstat_data.copy()
                thr_zstat_data[~stat_mask] = 0
                thr_zstat_img  = nib.Nifti1Image(thr_zstat_data, affine)
                thr_zstat_path = opj(output_folder, f"{seed_name}{tag}_{glt_label}_zstat_thr.nii.gz")
                thr_zstat_img.to_filename(thr_zstat_path)

                # Coefficient map masked by stat threshold
                thr_coef_data = coef_data.copy()
                thr_coef_data[~stat_mask] = 0
                thr_coef_img  = nib.Nifti1Image(thr_coef_data, affine)
                thr_coef_path = opj(output_folder, f"{seed_name}{tag}_{glt_label}_coef_thr.nii.gz")
                thr_coef_img.to_filename(thr_coef_path)

                # Plot 2: thresholded coefficient map
                display = plotting.plot_stat_map(
                    thr_coef_path, dim=0, threshold=1e-6, vmax=vmax_coef,
                    colorbar=True, bg_img=studytemplatebrain,
                    display_mode="y", cut_coords=cut_coords,
                    title=f"{seed_name} ? {glt_label} coef (p<{alpha}, cluster?{cluster_size})",
                )
                display.savefig(opj(output_folder, f"{seed_name}{tag}_{glt_label}_correlation_thr.jpg"))
                display.close()

                # Plot 3: thresholded z-stat map
                display = plotting.plot_stat_map(
                    thr_zstat_path, dim=0, threshold=z_score, vmax=vmax_zstat,
                    colorbar=True, bg_img=studytemplatebrain,
                    display_mode="y", cut_coords=cut_coords,
                    title=f"{seed_name} ? {glt_label} Z-stat (p<{alpha}, cluster?{cluster_size})",
                )
                display.savefig(opj(output_folder, f"{seed_name}{tag}_{glt_label}_zstat_thr.jpg"))
                display.close()
                plt.close("all")