import os
import subprocess
import numpy as np
import pandas as pd
import nibabel as nib
import nilearn
from nilearn import plotting
from nilearn.masking import compute_epi_mask
import matplotlib.pyplot as plt
from scipy.stats import norm
from Tools import Load_EDNiX_requirement
from Statistics.Group_fMRI._mask_utils import build_group_mask

opj = os.path.join
opb = os.path.basename
ope = os.path.exists

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


def _format_seed_group(name):
    name = name.replace("L_", "").replace("R_", "")
    return _format_seed(name)


def check_grid_match(reference, target):
    ref, tgt = nib.load(reference), nib.load(target)
    return np.allclose(ref.affine, tgt.affine) and ref.shape == tgt.shape


def mirror_hemisphere(img_path, output_dir, midline_x=0.0, is_left_seed=True):
    """Mirror one hemisphere onto the other across midline_x (mm)."""
    img    = nib.load(img_path)
    data   = img.get_fdata()
    affine = img.affine
    mirror = np.copy(data)
    nx     = data.shape[0]

    j0, k0 = data.shape[1] // 2, data.shape[2] // 2
    xs          = np.arange(nx)
    vox_coords  = np.stack([xs, np.full(nx, j0), np.full(nx, k0), np.ones(nx)], axis=0)
    world_xs    = (affine @ vox_coords)[0]
    midline_vox = int(np.argmin(np.abs(world_xs - midline_x)))
    max_copy    = min(midline_vox, nx - midline_vox - 1)

    if is_left_seed:
        for i in range(max_copy):
            lx, rx = midline_vox - i - 1, midline_vox + i + 1
            if lx >= 0 and rx < nx:
                mirror[rx] = data[lx]
    else:
        for i in range(max_copy):
            rx, lx = midline_vox + i + 1, midline_vox - i - 1
            if rx < nx and lx >= 0:
                mirror[lx] = data[rx]

    out_path = opj(output_dir, opb(img_path).replace(".nii.gz", "_mirrored.nii.gz"))
    nib.save(nib.Nifti1Image(mirror, affine, img.header), out_path)
    return out_path


def run_3dLMEr_analysis(output_folder, mask_path, design_matrix_path, model, glt_spec, sing_afni, tag=""):
    stat_maps  = opj(output_folder, f"3dLME_results{tag}.nii.gz")
    resid_path = opj(output_folder, f"resid{tag}.nii.gz")
    log_path   = opj(output_folder, f"3dLME_results{tag}_log.txt")
    for f in [stat_maps, resid_path, log_path]:
        if ope(f):
            os.remove(f)

    glt_args = " ".join(
        f"-gltCode {code} '{name}' '{contrast}'" for code, name, contrast in glt_spec
    )
    cmd = (
        f"{sing_afni} "
        f"3dLMEr -prefix {stat_maps} -jobs 20 -mask {mask_path} "
        f"-model \"{model}\" {glt_args} "
        f"-dataTable @{design_matrix_path} -resid {resid_path}"
    )
    print(f"  Running 3dLMEr  model={model}")
    try:
        result = subprocess.run(
            cmd, shell=True, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        with open(log_path, "w") as f:
            f.write(f"=== Command ===\n{cmd}\n\n=== STDOUT ===\n{result.stdout}\n=== STDERR ===\n{result.stderr}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"3dLMEr failed (exit {e.returncode}):\n{e.stderr}") from e
    return stat_maps


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
    glt_spec,
    contrast_names,
    midline_x,
    visualize,
    percent,
    regions_of_interest=None,
    opening=1,
    csim_pthr=0.05,
    vmax_percentile=99,
    method_mask_func='mask_func_over_Gray',
    output_suffix=None,
):
    """
    Seed-based group analysis with hemisphere mirroring using AFNI 3dLMEr.

    Parameters
    ----------
    glt_spec       : list of (code, name, contrast) tuples for 3dLMEr -gltCode
    contrast_names : list of output labels matching glt_spec order
    midline_x      : MNI x-coordinate of brain midline in mm
    visualize      : 'threshold' | 'percentile' | None
    percent        : percentile cutoff when visualize='percentile'
    opening        : opening radius for compute_epi_mask
    csim_pthr      : per-voxel p-threshold for ClustSim table lookup
    vmax_percentile : percentile [0-100] for colorscale clipping
    method_mask_func : mask strategy (see _mask_utils.build_group_mask)
    output_suffix  : optional string appended to output file/folder names
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

    tag            = f"_{output_suffix}" if output_suffix else ""
    folder_name    = "Grp_SBA_3dLME_network" + tag
    output_results = opj(output_results1, folder_name)
    os.makedirs(output_results, exist_ok=True)

    # ?? Build bilateral seed pairs ????????????????????????????????????????????
    region_names = set(label_df["region_name"].values)
    seed_pairs   = {}
    for _, row in label_df.iterrows():
        name = row["region_name"]
        if not name.startswith("L_"):
            continue
        if regions_of_interest and not any(r in _format_seed(name) for r in regions_of_interest):
            continue
        right = "R_" + name[2:]
        if right in region_names:
            seed_pairs[name] = right

    for left_raw, right_raw in seed_pairs.items():
        left_seed  = _format_seed(left_raw)
        right_seed = _format_seed(right_raw)
        base_name  = _format_seed_group(left_raw)
        print(f"\n  Seed pair: {left_seed} / {right_seed}  ?  {base_name}")

        output_folder = opj(output_results, base_name)
        os.makedirs(output_folder, exist_ok=True)

        design_data = []
        master_file = None

        for current_seed, is_left in [(left_seed, True), (right_seed, False)]:
            for subj, ses, run, bold_p in zip(all_subjects, all_sessions, all_runs, all_bold_list):
                fish = _sba_fish_path(bold_p, current_seed)
                if not ope(fish):
                    print(f"  [SKIP] {fish}")
                    continue

                if master_file is None:
                    master_file = fish
                    final_fish  = fish
                elif not check_grid_match(master_file, fish):
                    out_rs = fish.replace(".nii.gz", "_resampled.nii.gz")
                    subprocess.run(
                        f"{sing_afni} 3dresample -master {master_file} "
                        f"-input {fish} -prefix {out_rs} -overwrite",
                        shell=True, check=True,
                    )
                    final_fish = out_rs if ope(out_rs) else fish
                else:
                    final_fish = fish

                sba_dir       = os.path.dirname(fish)
                mirrored_fish = mirror_hemisphere(final_fish, sba_dir, midline_x, is_left)
                run_label     = f"run_{run}" if run is not None else "run_01"

                design_data.append({
                    "Subj":       str(subj),
                    "Sess":       f"Sess_{ses}",
                    "run":        run_label,
                    "Hemisphere": "L" if is_left else "R",
                    "InputFile":  mirrored_fish,
                })

        if not design_data:
            print(f"  [SKIP] No valid data for {base_name}")
            continue

        dm_df   = pd.DataFrame(design_data)[["Subj", "Sess", "run", "Hemisphere", "InputFile"]]
        dm_path = opj(output_folder, "design_matrix.txt")
        dm_df.to_csv(dm_path, sep="\t", index=False)

        stat_maps = run_3dLMEr_analysis(
            output_folder, mask_overlap, dm_path, model, glt_spec, sing_afni, tag=tag
        )

        subprocess.run(
            f"{sing_afni} 3dClustSim -mask {mask_overlap} -LOTS "
            f"-prefix {output_folder}/Clust_",
            shell=True,
        )

        # ?? Extract sub-bricks ? NIfTI ????????????????????????????????????????
        full_img  = nib.load(stat_maps)
        full_data = full_img.get_fdata()
        affine    = full_img.affine
        z_score   = norm.ppf(1 - alpha / 2)

        # Chi-sq map (brick 0)
        fmap_path = opj(output_folder, f"Hemisphere_ChiSq{tag}.nii.gz")
        nib.save(nib.Nifti1Image(full_data[..., 0], affine), fmap_path)

        for i, name in enumerate(contrast_names):
            coef_data  = full_data[..., 1 + 2 * i]
            zstat_data = full_data[..., 2 + 2 * i]

            coef_path  = opj(output_folder, f"{name}{tag}_coef.nii.gz")
            zstat_path = opj(output_folder, f"{name}{tag}_zstat.nii.gz")
            nib.save(nib.Nifti1Image(coef_data,  affine), coef_path)
            nib.save(nib.Nifti1Image(zstat_data, affine), zstat_path)

            vmax_coef  = np.percentile(np.abs(coef_data),  vmax_percentile) if np.any(coef_data)  else 1.0
            vmax_zstat = np.percentile(np.abs(zstat_data), vmax_percentile) if np.any(zstat_data) else 5.0

        if visualize == "threshold":
            _visualize_threshold(
                full_data, affine, output_folder, contrast_names,
                studytemplatebrain, cut_coords, alpha, csim_pthr,
                vmax_percentile, tag,
            )
        elif visualize == "percentile":
            _visualize_percentile(
                full_data, affine, output_folder, contrast_names,
                studytemplatebrain, cut_coords, percent,
                vmax_percentile, tag,
            )


def _visualize_threshold(
    full_data, affine, output_folder, contrast_names,
    template, cut_coords, alpha, csim_pthr, vmax_percentile, tag,
):
    cluster_file = opj(output_folder, "Clust_.NN1_2sided.1D")
    try:
        cluster_thresh = float(open(cluster_file).readlines()[8].split()[1])
    except Exception as e:
        print(f"  [WARN] {e}; cluster_thresh=10")
        cluster_thresh = 10

    z_score = norm.ppf(1 - alpha / 2)

    # Chi-sq plot
    fmap_path = opj(output_folder, f"Hemisphere_ChiSq{tag}.nii.gz")
    plotting.plot_stat_map(
        fmap_path, bg_img=template, display_mode="y", cut_coords=cut_coords,
        title="Hemisphere Effect (?²)", cmap="hot", colorbar=True,
    ).savefig(opj(output_folder, f"Hemisphere_ChiSq{tag}.png"), dpi=300)
    plt.close()

    for i, name in enumerate(contrast_names):
        coef_data  = full_data[..., 1 + 2 * i]
        zstat_data = full_data[..., 2 + 2 * i]

        coef_path  = opj(output_folder, f"{name}{tag}_coef.nii.gz")
        zstat_path = opj(output_folder, f"{name}{tag}_zstat.nii.gz")

        vmax_coef  = np.percentile(np.abs(coef_data),  vmax_percentile) if np.any(coef_data)  else 1.0
        vmax_zstat = np.percentile(np.abs(zstat_data), vmax_percentile) if np.any(zstat_data) else 5.0

        # Unthresholded coef plot
        cmap_coef = "cold_hot" if (np.any(coef_data < 0) and np.any(coef_data > 0)) else "hot"
        plotting.plot_stat_map(
            coef_path, bg_img=template, display_mode="y", cut_coords=cut_coords,
            title=f"{name} (Coefficients)", cmap=cmap_coef, colorbar=True,
        ).savefig(opj(output_folder, f"{name}{tag}_coef.png"), dpi=300)
        plt.close()

        # Build stat mask and threshold both maps
        stat_mask = np.abs(zstat_data) >= z_score

        thr_zstat_data = zstat_data.copy()
        thr_zstat_data[~stat_mask] = 0
        thr_coef_data  = coef_data.copy()
        thr_coef_data[~stat_mask]  = 0

        thr_zstat_path = opj(output_folder, f"{name}{tag}_zstat_thr.nii.gz")
        thr_coef_path  = opj(output_folder, f"{name}{tag}_coef_thr.nii.gz")
        nib.save(nib.Nifti1Image(thr_zstat_data, affine), thr_zstat_path)
        nib.save(nib.Nifti1Image(thr_coef_data,  affine), thr_coef_path)

        has_pos = np.any(thr_zstat_data > 0)
        has_neg = np.any(thr_zstat_data < 0)
        cmap    = "cold_hot" if (has_pos and has_neg) else "hot" if has_pos else "winter"

        plotting.plot_stat_map(
            thr_zstat_path, bg_img=template, display_mode="y", cut_coords=cut_coords,
            title=f"{name} (|Z|?{z_score:.2f})", cmap=cmap, colorbar=True,
            threshold=z_score, vmax=vmax_zstat,
        ).savefig(opj(output_folder, f"{name}{tag}_thr.png"), dpi=300)
        plt.close()


def _visualize_percentile(
    full_data, affine, output_folder, contrast_names,
    template, cut_coords, percent, vmax_percentile, tag,
):
    # Chi-sq plot
    fmap_path = opj(output_folder, f"Hemisphere_ChiSq{tag}.nii.gz")
    plotting.plot_stat_map(
        fmap_path, bg_img=template, display_mode="y", cut_coords=cut_coords,
        title="Hemisphere Effect (?²)", cmap="hot", colorbar=True,
    ).savefig(opj(output_folder, f"Hemisphere_ChiSq{tag}.png"), dpi=300)
    plt.close()

    for i, name in enumerate(contrast_names):
        coef_data  = full_data[..., 1 + 2 * i]
        zstat_data = full_data[..., 2 + 2 * i]

        coef_path  = opj(output_folder, f"{name}{tag}_coef.nii.gz")
        zstat_path = opj(output_folder, f"{name}{tag}_zstat.nii.gz")

        pos_vals = zstat_data[zstat_data > 0]
        neg_vals = zstat_data[zstat_data < 0]
        pos_thr  = np.percentile(pos_vals, 100 - percent) if len(pos_vals) else np.inf
        neg_thr  = np.percentile(neg_vals, percent)       if len(neg_vals) else -np.inf

        thr_data       = np.zeros_like(zstat_data)
        mask           = (zstat_data >= pos_thr) | (zstat_data <= neg_thr)
        thr_data[mask] = zstat_data[mask]

        thr_path = opj(output_folder, f"{name}{tag}_thr_pct.nii.gz")
        nib.save(nib.Nifti1Image(thr_data, affine), thr_path)

        if   len(pos_vals) > 100 and len(neg_vals) > 100: title, cmap = f"{name} (Top/Bottom {percent}%)", "cold_hot"
        elif len(pos_vals) > 100:                          title, cmap = f"{name} (Top {percent}%)",        "hot"
        elif len(neg_vals) > 100:                          title, cmap = f"{name} (Bottom {percent}%)",     "winter"
        else:
            print(f"  [SKIP] {name} ? no non-zero values")
            continue

        vmax_thr = np.percentile(np.abs(thr_data[mask]), vmax_percentile) if mask.any() else 5.0
        plotting.plot_stat_map(
            thr_path, bg_img=template, display_mode="y", cut_coords=cut_coords,
            title=title, cmap=cmap, colorbar=True, vmax=vmax_thr,
        ).savefig(opj(output_folder, f"{name}{tag}_thr_pct.png"), dpi=300)
        plt.close()