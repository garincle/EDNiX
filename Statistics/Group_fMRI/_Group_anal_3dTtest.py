import nilearn
from nilearn import plotting, image
import glob
import subprocess
import os
import numpy as np
import nibabel as nib
from nilearn.masking import compute_epi_mask
import matplotlib.pyplot as plt
from scipy.stats import norm
from Tools import Load_EDNiX_requirement, check_nii, getpath

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
        os.path.dirname(
            os.path.dirname(bold_path)),
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
    label_df,           # ← output of parse_label_file()
    lower_cutoff,
    upper_cutoff,
    MAIN_PATH,
    alpha,
    bold_paths,         # ← output of extract_bold_paths()
    mean_imgs,
    regions_of_interest=None,  # optional list of base_region substrings to process
):
    """
    Seed-based group t-test using AFNI 3dttest++.

    Parameters
    ----------
    bold_paths : dict — output of extract_bold_paths()
                 keys: 'subject', 'session', 'run', 'bold_path'

    label_df   : pd.DataFrame — output of parse_label_file()
                 columns: region_name, base_region, label_id, hemisphere, R, G, B, A
                 Each row is one seed; region_name (e.g. 'L_Isocortex') is used as
                 the SBA seed folder name.

    regions_of_interest : list of base_region substrings to restrict the analysis,
                          e.g. ['Isocortex', 'Frontal'].  None = all regions.
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

    all_subjects   = bold_paths["subject"]
    all_bold_list  = bold_paths["bold_path"]

    output_results = opj(output_results1, "Grp_SBA_3dTTEST")
    os.makedirs(output_results, exist_ok=True)

    # Filter label_df by regions_of_interest if requested
    seeds = label_df.copy()
    if regions_of_interest:
        seeds = seeds[seeds["region_name"].apply(
            lambda b: any(r in _format_seed(b) for r in regions_of_interest))]

    for _, seed_row in seeds.iterrows():
        seed_name     = _format_seed(seed_row["region_name"])
        print('proccessing ' + str(seed_name) + ' seed')
        output_folder = opj(output_results, seed_name)
        os.makedirs(output_folder, exist_ok=True)

        # ── Per-subject average Fisher-z maps ─────────────────────────────
        mean_per_subject = []
        master_file = None   # first valid avg map — used as resample target for ttest

        for unique_id in sorted(set(all_subjects)):
            fish_maps = [
                _sba_fish_path(bold_p, seed_name)
                for subj, bold_p in zip(all_subjects, all_bold_list)
                if subj == unique_id and ope(_sba_fish_path(bold_p, seed_name))
            ]
            if not fish_maps:
                print(f"  [WARN] No Fisher-z maps: subject={unique_id}  seed={seed_name}")
                continue

            # ── Resample all runs to the first run before 3dMean ──────────
            # Runs from different sessions may have different grid sizes
            master_run     = fish_maps[0]
            resampled_maps = [master_run]
            for fm in fish_maps[1:]:
                fm_rs = fm.replace(".nii.gz", "_rs.nii.gz")
                print(spgo(
                    f"{sing_afni} "
                    f"3dresample -master {master_run} -input {fm} "
                    f"-prefix {fm_rs} -overwrite -bound_type SLAB"
                ))
                resampled_maps.append(fm_rs if ope(fm_rs) else fm)

            avg_file = opj(output_folder, f"{unique_id}_avg_fisher_map.nii.gz")
            if ope(avg_file):
                os.remove(avg_file)
            print(spgo(
                f"{sing_afni} "
                f"3dMean -prefix {avg_file} "
                f"{subprocess.list2cmdline(resampled_maps)}"
            ))
            if ope(avg_file):
                mean_per_subject.append(avg_file)
                if master_file is None:
                    master_file = avg_file   # first valid avg = master for ttest

        if len(mean_per_subject) < 2:
            print(f"  [SKIP] {seed_name}: only {len(mean_per_subject)} subject(s)")
            continue

        # ── Resample all avg maps to master before 3dttest++ ──────────────
        # Subjects processed on different sessions may have different grids
        resampled_per_subject = [master_file]
        for avg in mean_per_subject[1:]:
            avg_rs = avg.replace(".nii.gz", "_rs.nii.gz")
            print(spgo(
                f"{sing_afni} "
                f"3dresample -master {master_file} -input {avg} "
                f"-prefix {avg_rs} -overwrite -bound_type SLAB"
            ))
            resampled_per_subject.append(avg_rs if ope(avg_rs) else avg)

        # ── 3dttest++ ─────────────────────────────────────────────────────
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

        zmap = opj(output_folder, f"{seed_name}_ttest-stat_fisher_zmap.nii.gz")
        print(spgo(
            f"{sing_afni} "
            f"3dcalc -overwrite "
            f"-a {opj(output_folder, 'TTnew+orig.HEAD[1]')} "
            f"-expr a -prefix {zmap}"
        ))

        correlation = opj(output_folder, f"{seed_name}_ttest-stat_fisher_cmap.nii.gz")
        print(spgo(
            f"{sing_afni} "
            f"3dcalc -overwrite "
            f"-a {opj(output_folder, 'TTnew+orig.HEAD[0]')} "
            f"-expr a -prefix {correlation}"))

        # ── Cluster threshold ─────────────────────────────────────────────
        if use_clustsim:
            csim_out = spgo(
                f"{sing_afni} "
                f"1d_tool.py -infile {output_folder}/TTnew.CSimA.NN1_2sided.1D "
                f"-csim_show_clustsize -verb 0 -csim_pthr 0.05 -csim_alpha {alpha}"
            )
            try:
                cluster_size = next(
                    int(l.strip()) for l in csim_out.splitlines() if l.strip().isdigit()
                )
            except StopIteration:
                cluster_size = 10
        else:
            cluster_size = 10

        # ── Threshold & visualise ─────────────────────────────────────────
        loadimg = nib.load(zmap).get_fdata()
        vmax = np.percentile(np.abs(loadimg), 99) if np.any(loadimg) else 5.0
        z_score = norm.ppf(1 - alpha / 2)

        # Threshold the z-stat map (for stat plot)
        thr_img = nilearn.image.threshold_img(zmap, z_score, cluster_threshold=int(cluster_size))
        thr_path = opj(output_folder, f"{seed_name}_thresholded_ttest-stat_fisher.nii.gz")
        thr_img.to_filename(thr_path)

        # Threshold the correlation map using the same z-score mask
        # (show only voxels that survive the stat threshold)
        thr_cmap = nilearn.image.threshold_img(correlation, z_score, cluster_threshold=int(cluster_size))
        thr_cmap_path = opj(output_folder, f"{seed_name}_thresholded_ttest-stat_fisher_cmap.nii.gz")
        thr_cmap.to_filename(thr_cmap_path)

        corr_data = nib.load(correlation).get_fdata()
        vmax_corr = np.max(np.abs(corr_data)) if np.any(corr_data) else 1.0

        # ── Plot 1: unthresholded group-mean Fisher-z correlation map ────
        display = plotting.plot_stat_map(
            correlation, dim=0, vmax=vmax_corr,
            colorbar=True, bg_img=studytemplatebrain,
            display_mode="y", cut_coords=cut_coords,
            title=f"{seed_name} — mean Fisher-z (unthresholded)",
        )
        display.savefig(opj(output_folder, f"{seed_name}_correlation_mosaic.jpg"))
        display.close()

        # ── Plot 2: correlation map masked by stat threshold ──────────────
        display = plotting.plot_stat_map(
            thr_cmap, dim=0, threshold=0, vmax=vmax_corr,
            colorbar=True, bg_img=studytemplatebrain,
            display_mode="y", cut_coords=cut_coords,
            title=f"{seed_name} — Fisher-z (p<{alpha}, cluster≥{cluster_size})",
        )
        display.savefig(opj(output_folder, f"{seed_name}_correlation_thresholded_mosaic.jpg"))
        display.close()

        # ── Plot 3: thresholded t-stat z-score map ────────────────────────
        display = plotting.plot_stat_map(
            thr_img, dim=0, threshold=z_score, vmax=vmax,
            colorbar=True, bg_img=studytemplatebrain,
            display_mode="y", cut_coords=cut_coords,
            title=f"{seed_name} — t-stat z-map (p<{alpha}, cluster≥{cluster_size})",
        )
        display.savefig(opj(output_folder, f"{seed_name}_thresholded_stat_mosaic.jpg"))
        display.close()
        plt.close("all")