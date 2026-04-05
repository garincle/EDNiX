import os
from Tools import Load_EDNiX_requirement, check_nii, getpath
import subprocess
import numpy as np
import pandas as pd
import nibabel as nib
import nilearn
from nilearn import plotting
from nilearn.masking import compute_epi_mask
import matplotlib.pyplot as plt
from scipy.stats import norm
import Tools.Load_EDNiX_requirement

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
        os.path.dirname(
            os.path.dirname(bold_path)),
        "Stats", "SBA", seed_name,
        f"{_bold_root(bold_path)}_correlations_fish.nii.gz",
    )


def _format_seed(name):
    for ch in " ()/,:;.-":
        name = name.replace(ch, "_")
    return "".join(c for c in name if c.isalnum() or c == "_").strip("_")


def _format_seed_group(name):
    """Remove L_/R_ hemisphere prefix then sanitise."""
    name = name.replace("L_", "").replace("R_", "")
    return _format_seed(name)


def check_grid_match(reference, target):
    ref, tgt = nib.load(reference), nib.load(target)
    return np.allclose(ref.affine, tgt.affine) and ref.shape == tgt.shape


def resample_to_match(reference, target, sing_afni, sing_wb):
    out = target.replace(".nii.gz", "_resampled.nii.gz")
    subprocess.run(
        f"{sing_afni} "
        f"3dresample -master {reference} -input {target} -prefix {out} -overwrite",
        shell=True, check=True,
    )
    return out


def mirror_hemisphere(img_path, output_dir, midline_x=0.0, is_left_seed=True):
    """
    Mirror one hemisphere onto the other across the midline (x = midline_x mm).

    Strategy
    --------
    For a LEFT seed  : copy left-of-midline voxels → symmetric right positions.
    For a RIGHT seed : copy right-of-midline voxels → symmetric left positions.
    The midline voxel is found using the full affine (handles oblique volumes).
    The affine and header are preserved — the image stays in the same space.

    Parameters
    ----------
    midline_x : x-coordinate of the midline in world (mm) space.
                Default 0.0 (standard MNI/ACPC midline).
    """
    img = nib.load(img_path)
    data = img.get_fdata()
    affine = img.affine
    mirror = np.copy(data)
    nx = data.shape[0]

    # World x-coordinate for every voxel index along axis 0
    # Use full affine: world_x = affine[0,0]*i + affine[0,1]*j + affine[0,2]*k + affine[0,3]
    # Here we want the x-coordinate along the first axis at j=k=0 (voxel col/slice centre)
    # More robustly: compute world x at each i with j,k fixed at image centre
    j0, k0 = data.shape[1] // 2, data.shape[2] // 2
    xs = np.arange(nx)
    # homogeneous coords for all i at (j0, k0)
    vox_coords = np.stack([xs,
                           np.full(nx, j0),
                           np.full(nx, k0),
                           np.ones(nx)], axis=0)  # 4 × nx
    world_coords = affine @ vox_coords  # 4 × nx
    world_xs = world_coords[0]  # x in mm

    midline_vox = int(np.argmin(np.abs(world_xs - midline_x)))
    max_copy = min(midline_vox, nx - midline_vox - 1)

    if is_left_seed:
        # Copy each left voxel to its symmetric right position
        for i in range(max_copy):
            lx = midline_vox - i - 1
            rx = midline_vox + i + 1
            if lx >= 0 and rx < nx:
                mirror[rx] = data[lx]
    else:
        # Copy each right voxel to its symmetric left position
        for i in range(max_copy):
            rx = midline_vox + i + 1
            lx = midline_vox - i - 1
            if rx < nx and lx >= 0:
                mirror[lx] = data[rx]
    # NOTE: no mirror[::-1] — that would flip the entire volume and break the affine

    out_path = opj(output_dir, opb(img_path).replace(".nii.gz", "_mirrored.nii.gz"))
    nib.save(nib.Nifti1Image(mirror, affine, img.header), out_path)
    return out_path



def run_3dLMEr_analysis(output_folder, mask_path, design_matrix_path, model, glt_spec, sing_afni, sing_wb):
    stat_maps  = opj(output_folder, "3dLME_results.nii.gz")
    resid_path = opj(output_folder, "resid.nii.gz")
    log_path   = opj(output_folder, "3dLME_results_log.txt")
    for f in [stat_maps, resid_path, log_path]:
        if ope(f): os.remove(f)

    glt_args = " ".join(
        f"-gltCode {code} '{name}' '{contrast}'" for code, name, contrast in glt_spec
    )
    cmd = (
        f"{sing_afni} "
        f"3dLMEr -prefix {stat_maps} -jobs 20 -mask {mask_path} "
        f"-model \"{model}\" {glt_args} "
        f"-dataTable @{design_matrix_path} -resid {resid_path}"
    )
    print(f"Running 3dLMEr  model={model}  GLTs={glt_spec}")
    try:
        result = subprocess.run(cmd, shell=True, check=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True)
        with open(log_path, "w") as f:
            f.write(f"=== Command ===\n{cmd}\n\n=== STDOUT ===\n{result.stdout}\n=== STDERR ===\n{result.stderr}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"3dLMEr failed (exit {e.returncode}):\n{e.stderr}") from e
    return stat_maps


def visualize_results(stat_maps, output_folder, contrast_names, template, cut_coords, alpha):
    cluster_file = opj(output_folder, "Clust_.NN1_2sided.1D")
    try:
        cluster_thresh = float(open(cluster_file).readlines()[8].split()[1])
    except Exception as e:
        print(f"  [WARN] {e}; cluster_thresh=10")
        cluster_thresh = 10

    z_score      = norm.ppf(1 - alpha / 2)
    img          = nib.load(stat_maps)
    data, affine = img.get_fdata(), img.affine

    fmap = opj(output_folder, "Hemisphere_ChiSq.nii.gz")
    nib.save(nib.Nifti1Image(data[..., 0], affine), fmap)
    plotting.plot_stat_map(fmap, bg_img=template, display_mode="y", cut_coords=cut_coords,
                           title="Hemisphere Effect (χ²)", cmap="hot", colorbar=True
                           ).savefig(opj(output_folder, "Hemisphere_ChiSq.png"), dpi=300)
    plt.close()

    for i, name in enumerate(contrast_names):
        coef_data  = data[..., 1 + 2 * i]
        zstat_data = data[..., 2 + 2 * i]
        nib.save(nib.Nifti1Image(coef_data,  affine), opj(output_folder, f"{name}_coef.nii.gz"))
        nib.save(nib.Nifti1Image(zstat_data, affine), opj(output_folder, f"{name}_zstat.nii.gz"))

        cmap_coef = "cold_hot" if (np.any(coef_data < 0) and np.any(coef_data > 0)) else "hot"
        plotting.plot_stat_map(nib.Nifti1Image(coef_data, affine), bg_img=template,
                               display_mode="y", cut_coords=cut_coords,
                               title=f"{name} (Coefficients)", cmap=cmap_coef, colorbar=True
                               ).savefig(opj(output_folder, f"{name}_coef.png"), dpi=300)
        plt.close()

        # threshold positive tail: keep voxels >= +z_score
        pos_thr  = nilearn.image.threshold_img(nib.Nifti1Image(zstat_data, affine), z_score, cluster_threshold=cluster_thresh)
        # threshold negative tail: negate, threshold at +z_score, negate back
        neg_thr  = -nilearn.image.threshold_img(nib.Nifti1Image(-zstat_data, affine), z_score, cluster_threshold=cluster_thresh).get_fdata()
        combined = pos_thr.get_fdata() + neg_thr
        comb_img = nib.Nifti1Image(combined, affine)
        nib.save(comb_img, opj(output_folder, f"{name}_thr.nii.gz"))
        has_pos, has_neg = np.any(combined > 0), np.any(combined < 0)
        cmap = "cold_hot" if (has_pos and has_neg) else "hot" if has_pos else "winter"
        plotting.plot_stat_map(comb_img, bg_img=template, display_mode="y", cut_coords=cut_coords,
                               title=f"{name} (|Z|≥{z_score:.2f})", cmap=cmap, colorbar=True
                               ).savefig(opj(output_folder, f"{name}_thr.png"), dpi=300)
        plt.close()


def visualize_results_percentile(stat_maps, output_folder, contrast_names, template, cut_coords, percent):
    img          = nib.load(stat_maps)
    data, affine = img.get_fdata(), img.affine

    fmap = opj(output_folder, "Hemisphere_ChiSq.nii.gz")
    nib.save(nib.Nifti1Image(data[..., 0], affine), fmap)
    plotting.plot_stat_map(fmap, bg_img=template, display_mode="y", cut_coords=cut_coords,
                           title="Hemisphere Effect (χ²)", cmap="hot", colorbar=True
                           ).savefig(opj(output_folder, "Hemisphere_ChiSq.png"), dpi=300)
    plt.close()

    for i, name in enumerate(contrast_names):
        coef_data  = data[..., 1 + 2 * i]
        zstat_data = data[..., 2 + 2 * i]
        nib.save(nib.Nifti1Image(coef_data,  affine), opj(output_folder, f"{name}_coef.nii.gz"))
        nib.save(nib.Nifti1Image(zstat_data, affine), opj(output_folder, f"{name}_zstat.nii.gz"))

        pos_vals = zstat_data[zstat_data > 0]
        neg_vals = zstat_data[zstat_data < 0]
        pos_thr  = np.percentile(pos_vals, 100 - percent) if len(pos_vals) else np.inf
        neg_thr  = np.percentile(neg_vals, percent)       if len(neg_vals) else -np.inf
        thr_data = np.zeros_like(zstat_data)
        mask     = (zstat_data >= pos_thr) | (zstat_data <= neg_thr)
        thr_data[mask] = zstat_data[mask]

        if   len(pos_vals) > 100 and len(neg_vals) > 100: title, cmap = f"{name} (Top/Bottom {percent}%)", "cold_hot"
        elif len(pos_vals) > 100:                          title, cmap = f"{name} (Top {percent}%)", "hot"
        elif len(neg_vals) > 100:                          title, cmap = f"{name} (Bottom {percent}%)", "winter"
        else:
            print(f"  [SKIP] {name} — no non-zero values"); continue

        plotting.plot_stat_map(nib.Nifti1Image(thr_data, affine), bg_img=template,
                               display_mode="y", cut_coords=cut_coords,
                               title=title, cmap=cmap, colorbar=True
                               ).savefig(opj(output_folder, f"{name}_thr_pct.png"), dpi=300)
        plt.close()


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
):
    """
    Seed-based group analysis with hemisphere mirroring using AFNI 3dLMEr.

    Parameters
    ----------
    bold_paths : dict — output of extract_bold_paths()
                 keys: 'subject', 'session', 'run', 'bold_path'

    label_df   : pd.DataFrame — output of parse_label_file()
                 columns: region_name, base_region, label_id, hemisphere, R, G, B, A
                 Bilateral seed pairs are detected automatically via the L_/R_ prefix.

    glt_spec   : list of (code, name, contrast) tuples for 3dLMEr -gltCode
    contrast_names : list of output labels matching glt_spec order
    midline_x  : MNI x-coordinate of brain midline in mm
    visualize  : 'threshold' | 'percentile' | None
    percent    : percentile cutoff when visualize='percentile'

    regions_of_interest : list of base_region substrings to restrict the
                          analysis.  None = all bilateral pairs in label_df.
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
        f"3dresample -overwrite -master {mask_path} -input {mask_func} -prefix {mask_orig}",
        f"3dcalc -overwrite -a {mask_path} -b {mask_orig} -expr 'a*b' -prefix {mask_overlap}",
    ]:
        subprocess.run(f"{sing_afni} {cmd}", shell=True, check=True)

    all_subjects  = bold_paths["subject"]
    all_sessions  = bold_paths["session"]
    all_runs      = bold_paths["run"]
    all_bold_list = bold_paths["bold_path"]

    output_results = opj(output_results1, "Grp_SBA_3dLME_network")
    os.makedirs(output_results, exist_ok=True)

    # ── Build bilateral seed pairs from label_df ──────────────────────────
    region_names = set(label_df["region_name"].values)
    seed_pairs = {}
    for _, row in label_df.iterrows():
        name = row["region_name"]
        if not name.startswith("L_"):
            continue
        if regions_of_interest and not any(r in _format_seed(row["region_name"]) for r in regions_of_interest):
            continue
        right = "R_" + name[2:]
        if right in region_names:
            seed_pairs[name] = right

    for left_raw, right_raw in seed_pairs.items():
        left_seed  = _format_seed(left_raw)
        right_seed = _format_seed(right_raw)
        base_name  = _format_seed_group(left_raw)
        print(f"\n  Seed pair: {left_seed} / {right_seed}  →  {base_name}")

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
                    final_fish  = resample_to_match(master_file, fish, sing_afni, sing_wb)
                else:
                    final_fish  = fish

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
            output_folder, mask_overlap, dm_path, model, glt_spec, sing_afni, sing_wb
        )

        subprocess.run(
            f"{sing_afni} "
            f"3dClustSim -mask {mask_overlap} -LOTS -prefix {output_folder}/Clust_",
            shell=True,
        )

        if visualize == "threshold":
            visualize_results(stat_maps, output_folder, contrast_names, studytemplatebrain, cut_coords, alpha)
        elif visualize == "percentile":
            visualize_results_percentile(stat_maps, output_folder, contrast_names, studytemplatebrain, cut_coords, percent)