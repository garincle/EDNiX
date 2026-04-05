import nilearn
from nilearn import plotting
from Tools import Load_EDNiX_requirement, check_nii, getpath
import subprocess
import os
import numpy as np
import pandas as pd
import nibabel as nib
from nilearn.masking import compute_epi_mask
import matplotlib.pyplot as plt
from scipy.stats import norm
import Tools.Load_EDNiX_requirement

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
        os.path.dirname(
            os.path.dirname(bold_path)),
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
        col_idx = alpha_values.index(alpha) + 1
        data    = np.loadtxt(cluster_file, comments="#")
        row_idx = np.where(np.isclose(data[:, 0], pthr, atol=1e-6))[0]
        if not len(row_idx):
            return 10
        cluster_size = data[row_idx[0], col_idx]
        print(f"  Cluster threshold pthr={pthr} alpha={alpha} → {cluster_size:.1f} voxels")
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
    regions_of_interest=None):
    """
    Seed-based group analysis using AFNI 3dLMEr (mixed-effects).

    Parameters
    ----------
    bold_paths : dict — output of extract_bold_paths()
                 keys: 'subject', 'session', 'run', 'bold_path'

    label_df   : pd.DataFrame — output of parse_label_file()
                 columns: region_name, base_region, label_id, hemisphere, R, G, B, A
                 Each row defines one seed; region_name (e.g. 'L_Isocortex') is used
                 as the SBA seed folder name.

    model      : 3dLMEr model string, e.g. '~1+(1|Subj)'
    gltCode    : full -gltCode argument string(s) passed verbatim to 3dLMEr
    stat_img   : list of contrast labels matching gltCode order,
                 e.g. ['groupeffect', 'contrast1'].
                 'groupeffect' skips cluster-based thresholding.

    regions_of_interest : list of base_region substrings to restrict the
                          analysis.  None = all regions in label_df.
    """
    sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _,sing_synstrip,Unetpath =  Load_EDNiX_requirement.load_requirement(MAIN_PATH,templatehigh,bids_dir,'yes')

    output_results1 = opj(bids_dir, "Results")
    os.makedirs(output_results1, exist_ok=True)
    studytemplatebrain = templatehigh if oversample_map else templatelow

    # ── Group brain mask ─────────────────────────────────────────────────────
    mean_imgs_rs = nilearn.image.concat_imgs(mean_imgs, auto_resample=True, verbose=0)
    mask_img     = compute_epi_mask(
        mean_imgs_rs,
        lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff,
        connected=True, opening=1, exclude_zeros=True, ensure_finite=True,
    )
    mask_path    = opj(output_results1, "mask_mean_func.nii.gz")
    mask_orig    = opj(output_results1, "mask_mean_func_orig.nii.gz")
    mask_overlap = opj(output_results1, "mask_mean_func_overlapp.nii.gz")
    mask_img.to_filename(mask_path)

    for cmd in [
        f"3dmask_tool -overwrite -input {mask_path} -prefix {mask_path} -fill_holes",
        f"3dresample -master {mask_path} -input {mask_func} -prefix {mask_orig} -overwrite -bound_type SLAB",
        f"3dcalc -a {mask_path} -b {mask_orig} -expr 'a*b' -prefix {mask_overlap} -overwrite",
    ]:
        print(spgo(f"{sing_afni} {cmd}"))

    all_subjects  = bold_paths["subject"]
    all_sessions  = bold_paths["session"]
    all_runs      = bold_paths["run"]
    all_bold_list = bold_paths["bold_path"]

    output_results = opj(output_results1, "Grp_SBA_3dLME_network")
    os.makedirs(output_results, exist_ok=True)

    seeds = label_df.copy()
    if regions_of_interest:
        seeds = seeds[seeds["region_name"].apply(
            lambda b: any(r in _format_seed(b) for r in regions_of_interest))]

    for _, seed_row in seeds.iterrows():
        seed_name     = _format_seed(seed_row["region_name"])
        output_folder = opj(output_results, seed_name)
        os.makedirs(output_folder, exist_ok=True)

        # ── Collect Fisher-z maps, build design matrix ────────────────────
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
                    f"{sing_afni} "
                    f"3dMatch -quiet -source {Resample_master} -input {fish}",
                    shell=True, capture_output=True, text=True,
                )
                if check.returncode != 0:
                    resampled = fish.replace(".nii.gz", "_resampled.nii.gz")
                    subprocess.run(
                        f"{sing_afni} "
                        f"3dresample -overwrite -master {Resample_master} "
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
        if ope(dm_path): os.remove(dm_path)
        design_df.to_csv(dm_path, index=False, sep="\t")

        # ── 3dLMEr ────────────────────────────────────────────────────────
        stat_maps = opj(output_folder, "3dLME_glt.nii.gz")
        resid_path = opj(output_folder, "resid.nii.gz")
        for f in [stat_maps, resid_path]:
            if ope(f): os.remove(f)

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
            f"{sing_afni} "
            f"3dClustSim -mask {mask_overlap} -LOTS -prefix {output_folder}/Clust_",
            shell=True,
        )

        # ── Sub-brick extraction, threshold, visualise ────────────────────
        # 3dLME sub-brick layout per glt:
        #   brick 0          : Chi-sq (global F-test)
        #   brick 1 + 2*i    : coefficient (effect / correlation) for glt i
        #   brick 2 + 2*i    : Z-stat for glt i
        os.chdir(output_results)
        z_score = norm.ppf(1 - alpha / 2)
        full_data = nib.load(stat_maps).get_fdata()
        affine = nib.load(stat_maps).affine

        for i, glt_label in enumerate(stat_img):
            coef_brick = 1 + 2 * i  # effect / correlation map
            zstat_brick = 2 + 2 * i  # Z-stat map

            coef_path = opj(output_folder, f"{seed_name}_{glt_label}_coef.nii.gz")
            zstat_path = opj(output_folder, f"{seed_name}_{glt_label}_zstat.nii.gz")

            nib.save(nib.Nifti1Image(full_data[..., coef_brick], affine), coef_path)
            nib.save(nib.Nifti1Image(full_data[..., zstat_brick], affine), zstat_path)

            coef_data = full_data[..., coef_brick]
            zstat_data = full_data[..., zstat_brick]
            vmax_coef = np.max(np.abs(coef_data)) if np.any(coef_data) else 1.0
            vmax_zstat = np.max(np.abs(zstat_data)) if np.any(zstat_data) else 5.0

            # ── Plot 1: unthresholded correlation / coefficient map ───────
            display = plotting.plot_stat_map(
                coef_path, dim=0, vmax=vmax_coef,
                colorbar=True, bg_img=studytemplatebrain,
                display_mode="y", cut_coords=cut_coords,
                title=f"{seed_name} — {glt_label} coef (unthresholded)",
            )
            display.savefig(opj(output_folder, f"{seed_name}_{glt_label}_correlation.jpg"))
            display.close()
            plt.close("all")

            if glt_label == "groupeffect":
                # groupeffect: no cluster thresholding — just plot the Z-stat
                display = plotting.plot_stat_map(
                    zstat_path, dim=0, vmax=vmax_zstat,
                    colorbar=True, bg_img=studytemplatebrain,
                    display_mode="y", cut_coords=cut_coords,
                    title=f"{seed_name} — {glt_label} Z-stat (unthresholded)",
                )
                display.savefig(opj(output_folder, f"{seed_name}_{glt_label}_zstat.jpg"))
                display.close()
                plt.close("all")
            else:
                cluster_size = _extract_cluster_threshold(
                    opj(output_folder, "Clust_.NN1_2sided.1D"), alpha, 0.05
                )

                # Threshold Z-stat map
                thr_zstat = nilearn.image.threshold_img(zstat_path, z_score, cluster_threshold=float(cluster_size))
                thr_zstat_path = opj(output_folder, f"{seed_name}_{glt_label}_zstat_thr.nii.gz")
                thr_zstat.to_filename(thr_zstat_path)

                # Apply same spatial mask to coefficient map
                thr_coef = nilearn.image.threshold_img(coef_path, z_score, cluster_threshold=float(cluster_size))
                thr_coef_path = opj(output_folder, f"{seed_name}_{glt_label}_coef_thr.nii.gz")
                thr_coef.to_filename(thr_coef_path)

                # ── Plot 2: coefficient map masked by stat threshold ──────
                display = plotting.plot_stat_map(
                    thr_coef_path, dim=0, threshold=0, vmax=vmax_coef,
                    colorbar=True, bg_img=studytemplatebrain,
                    display_mode="y", cut_coords=cut_coords,
                    title=f"{seed_name} — {glt_label} coef (p<{alpha}, cluster≥{cluster_size})",
                )
                display.savefig(opj(output_folder, f"{seed_name}_{glt_label}_correlation_thr.jpg"))
                display.close()

                # ── Plot 3: thresholded Z-stat map ────────────────────────
                display = plotting.plot_stat_map(
                    thr_zstat_path, dim=0, threshold=z_score, vmax=vmax_zstat,
                    colorbar=True, bg_img=studytemplatebrain,
                    display_mode="y", cut_coords=cut_coords,
                    title=f"{seed_name} — {glt_label} Z-stat (p<{alpha}, cluster≥{cluster_size})",
                )
                display.savefig(opj(output_folder, f"{seed_name}_{glt_label}_zstat_thr.jpg"))
                display.close()
                plt.close("all")