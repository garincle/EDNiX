import os
import glob
import subprocess
import numpy as np
import pandas as pd
import nibabel as nib
import nilearn
from nilearn import plotting
from nilearn.masking import compute_epi_mask
import matplotlib.pyplot as plt
from scipy.stats import norm

from nilearn.masking import compute_epi_mask
import matplotlib.pyplot as plt
from scipy.stats import norm
import Tools.Load_EDNiX_requirement
from fonctions.extract_filename import extract_filename
# Path utilities
opj = os.path.join
opb = os.path.basename
ope = os.path.exists


def check_grid_match(reference, target):
    """Check if two images have matching grids"""
    ref_img = nib.load(reference)
    target_img = nib.load(target)
    return np.allclose(ref_img.affine, target_img.affine) and ref_img.shape == target_img.shape


def resample_to_match(reference, target, afni_sif, s_bind):
    """Resample target to match reference grid"""
    output_path = target.replace('.nii.gz', '_resampled.nii.gz')
    command = f"singularity run {s_bind} {afni_sif} 3dresample -master {reference} -input {target} -prefix {output_path} -overwrite"
    subprocess.run(command, shell=True, check=True)
    return output_path


def format_seed_name(seed_name):
    replace_chars = [' ', '(', ')', ',', '/', ':', ';', '.', '-']
    for char in replace_chars:
        seed_name = seed_name.replace(char, '_')
    formatted_name = ''.join(char for char in seed_name if char.isalnum() or char == '_')
    formatted_name = formatted_name.strip('_')
    return formatted_name

def format_seed_name_Group(seed_name):
    """Format seed name by removing special characters and hemisphere prefixes"""
    seed_name = seed_name.replace('L_', '').replace('R_', '')
    replace_chars = [' ', '(', ')', ',', '/', ':', ';', '.', '-']
    for char in replace_chars:
        seed_name = seed_name.replace(char, '_')
    return ''.join(char for char in seed_name if char.isalnum() or char == '_').strip('_')


def mirror_hemisphere(img_path, output_dir, midline_x=1, is_left_seed=True):
    """Mirror image by copying values from one hemisphere to the other,
    skipping incomplete edge slices.

    Parameters:
    -----------
    img_path : str
        Path to input NIfTI image
    output_dir : str
        Directory to save mirrored image
    midline_x : float
        MNI x-coordinate of midline (default: 0.05166)
    is_left_seed : bool
        Whether the seed is in the left hemisphere

    Returns:
    --------
    str: Path to mirrored image
    """
    # Load the image
    img = nib.load(img_path)
    data = img.get_fdata()
    affine = img.affine

    # Create output array
    mirrored_data = np.copy(data)

    # Get precise midline voxel index
    coords = np.array([[i, 0, 0, 1] for i in range(data.shape[0])])
    world_coords = np.dot(affine, coords.T).T[:, 0]
    midline_vox = np.argmin(np.abs(world_coords - midline_x))

    if is_left_seed:
        # Left seed: copy left to right
        max_voxels_to_copy = min(midline_vox, data.shape[0] - midline_vox - 1)

        for i in range(max_voxels_to_copy):
            # Left source voxel
            left_x = midline_vox - i - 1
            # Right target voxel
            right_x = midline_vox + i + 1

            if left_x >= 0 and right_x < data.shape[0]:
                mirrored_data[right_x] = data[left_x]
    else:
        # Right seed: copy right to left
        max_voxels_to_copy = min(midline_vox, data.shape[0] - midline_vox - 1)

        for i in range(max_voxels_to_copy):
            # Right source voxel
            right_x = midline_vox + i + 1
            # Left target voxel
            left_x = midline_vox - i - 1

            if right_x < data.shape[0] and left_x >= 0:
                mirrored_data[left_x] = data[right_x]

        # Flip to match left-hemisphere orientation
        mirrored_data = mirrored_data[::-1]

    # Save result
    mirrored_img = nib.Nifti1Image(mirrored_data, affine, img.header)
    mirrored_path = opj(output_dir, f"{opb(img_path).replace('.nii.gz', '_mirrored.nii.gz')}")
    mirrored_img.to_filename(mirrored_path)

    return mirrored_path

def run_3dLMEr_analysis(output_folder, mask_path, design_matrix_path, model, glt_spec, afni_sif, s_bind):
    """
    Run 3dLMEr analysis with configurable model and GLT specs

    Parameters:
    -----------
    output_folder : str
        Path to output directory
    mask_path : str
        Path to mask file
    design_matrix_path : str
        Path to design matrix file
    model : str
        Mixed effects model formula (e.g., "(1|Subj)*Hemisphere")
    glt_spec : list of tuples
        GLT specifications as (code, name, contrast) tuples
    afni_sif : str
        Path to AFNI Singularity image
    s_bind : str
        Singularity bind paths

    Returns:
    --------
    str: Path to output statistical map
    """
    import os
    import subprocess
    from os.path import join as opj

    # Set output paths
    stat_maps = opj(output_folder, '3dLME_results.nii.gz')
    resid_path = opj(output_folder, 'resid.nii.gz')
    log_path = opj(output_folder, '3dLME_results_log.txt')

    # Clean previous results
    for f in [stat_maps, resid_path, log_path]:
        if os.path.exists(f):
            os.remove(f)

    # Build GLT arguments
    glt_args = []
    for code, name, contrast in glt_spec:
        glt_args.append(f"-gltCode {code} '{name}' '{contrast}'")
    glt_cmd = " ".join(glt_args)

    # Build the full command
    cmd = (
        f"singularity run {s_bind} {afni_sif} "
        f"3dLMEr "
        f"-prefix {stat_maps} "
        f"-jobs 20 "
        f"-mask {mask_path} "
        f"-model \"{model}\" "
        f"{glt_cmd} "
        f"-dataTable @{design_matrix_path} "
        f"-resid {resid_path}"
    )

    # Run with error handling
    try:
        print(f"Running 3dLMEr with:\nModel: {model}\nGLTs: {glt_spec}")
        result = subprocess.run(
            cmd,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )

        # Save logs
        with open(log_path, 'w') as f:
            f.write(f"=== Command ===\n{cmd}\n\n")
            f.write("=== STDOUT ===\n")
            f.write(result.stdout)
            f.write("\n=== STDERR ===\n")
            f.write(result.stderr)

        print("3dLMEr completed successfully!")
        return stat_maps

    except subprocess.CalledProcessError as e:
        error_msg = f"""
        3dLMEr failed with exit code {e.returncode}
        Command: {cmd}
        Error output:
        {e.stderr}
        """
        raise RuntimeError(error_msg) from e


def visualize_results(stat_maps, output_folder, contrast_names, studytemplatebrain, cut_coords, alpha):
    """Visualize AFNI 3dLMEr results with proper contrast handling"""
    # Load cluster thresholds
    cluster_file = opj(output_folder, 'Clust_.NN1_2sided.1D')
    try:
        with open(cluster_file) as f:
            lines = f.readlines()
        cluster_thresh = float(lines[8].split()[1])
    except Exception as e:
        print(f"Could not read cluster threshold file: {e}")
        cluster_thresh = 10  # Fallback value

    z_score = norm.ppf(1 - alpha / 2)

    # Load AFNI output
    img = nib.load(stat_maps)
    data = img.get_fdata()
    affine = img.affine

    # ==================================================================
    # 1. First handle the overall F-test (sub-brick 0 - Hemisphere Chi-sq)
    # ==================================================================
    fstat_data = data[..., 0]
    fstat_map = opj(output_folder, 'Hemisphere_ChiSq.nii.gz')
    nib.save(nib.Nifti1Image(fstat_data, affine), fstat_map)

    plotting.plot_stat_map(
        fstat_map,
        bg_img=studytemplatebrain,
        display_mode='mosaic',
        cut_coords=cut_coords,
        title='Hemisphere Effect (χ²)',
        cmap='hot',
        colorbar=True
    ).savefig(opj(output_folder, 'Hemisphere_ChiSq.png'), dpi=300)
    plt.close()

    # ==================================================================
    # 2. Handle each contrast (coefficients and Z-stats)
    # ==================================================================
    # Verify we have enough sub-bricks for all contrasts
    n_contrasts = len(contrast_names)
    expected_bricks = 1 + 2 * n_contrasts  # 1 F-stat + 2 per contrast (coef + zstat)

    if data.shape[-1] < expected_bricks:
        raise ValueError(f"Expected {expected_bricks} sub-bricks but got {data.shape[-1]}")

    # Process each contrast
    for i, contrast_name in enumerate(contrast_names):
        # Calculate indices for this contrast
        coef_idx = 1 + 2 * i
        zstat_idx = 2 + 2 * i

        # Save coefficient map
        coef_data = data[..., coef_idx]
        coef_map = opj(output_folder, f'{contrast_name}_coef.nii.gz')
        nib.save(nib.Nifti1Image(coef_data, affine), coef_map)

        # Save Z-stat map
        zstat_data = data[..., zstat_idx]
        zstat_map = opj(output_folder, f'{contrast_name}_zstat.nii.gz')
        nib.save(nib.Nifti1Image(zstat_data, affine), zstat_map)

        # Plot coefficients
        plotting.plot_stat_map(
            nib.Nifti1Image(coef_data, affine),
            bg_img=studytemplatebrain,
            display_mode='mosaic',
            cut_coords=cut_coords,
            title=f'{contrast_name} (Coefficients)',
            cmap='cold_hot' if (np.any(coef_data < 0) and np.any(coef_data > 0)) else 'hot',
            colorbar=True
        ).savefig(opj(output_folder, f'{contrast_name}_coef.png'), dpi=300)
        plt.close()

        # Threshold and plot Z-stats
        pos_thr = nilearn.image.threshold_img(
            nib.Nifti1Image(zstat_data, affine),
            threshold=z_score,
            cluster_threshold=cluster_thresh
        )
        neg_thr = nilearn.image.threshold_img(
            nib.Nifti1Image(zstat_data, affine),
            threshold=-z_score,
            cluster_threshold=cluster_thresh
        )

        # Combine positive and negative thresholds
        combined_data = pos_thr.get_fdata() - neg_thr.get_fdata()
        combined_map = nib.Nifti1Image(combined_data, affine)
        combined_path = opj(output_folder, f'{contrast_name}_thr.nii.gz')
        nib.save(combined_map, combined_path)

        # Determine colormap
        has_pos = np.any(combined_data > 0)
        has_neg = np.any(combined_data < 0)

        plotting.plot_stat_map(
            combined_map,
            bg_img=studytemplatebrain,
            display_mode='mosaic',
            cut_coords=cut_coords,
            title=f'{contrast_name} (|Z|>={z_score:.2f})',
            cmap='cold_hot' if (has_pos and has_neg) else 'hot' if has_pos else 'winter',
            colorbar=True
        ).savefig(opj(output_folder, f'{contrast_name}_thr.png'), dpi=300)
        plt.close()


def visualize_results_percentile(stat_maps, output_folder, contrast_names, studytemplatebrain, cut_coords, percent):
    """Visualize AFNI results using percentile-based thresholding with proper contrast handling"""
    # Load AFNI output
    img = nib.load(stat_maps)
    data = img.get_fdata()
    affine = img.affine

    # Verify we have enough sub-bricks for all contrasts
    n_contrasts = len(contrast_names)
    expected_bricks = 1 + 2 * n_contrasts  # 1 F-stat + 2 per contrast (coef + zstat)

    if data.shape[-1] < expected_bricks:
        raise ValueError(f"Expected {expected_bricks} sub-bricks but got {data.shape[-1]}")

    # ==================================================================
    # 1. First handle the overall F-test (sub-brick 0 - Hemisphere Chi-sq)
    # ==================================================================
    fstat_data = data[..., 0]
    fstat_map = opj(output_folder, 'Hemisphere_ChiSq.nii.gz')
    nib.save(nib.Nifti1Image(fstat_data, affine), fstat_map)

    # Plot F-stat with hot colormap (only positive values)
    plotting.plot_stat_map(
        fstat_map,
        bg_img=studytemplatebrain,
        display_mode='mosaic',
        cut_coords=cut_coords,
        title='Hemisphere Effect (χ²)',
        cmap='hot',
        colorbar=True
    ).savefig(opj(output_folder, 'Hemisphere_ChiSq.png'), dpi=300)
    plt.close()

    # ==================================================================
    # 2. Handle each contrast (coefficients and Z-stats)
    # ==================================================================
    for i, contrast_name in enumerate(contrast_names):
        # Calculate indices for this contrast
        coef_idx = 1 + 2 * i
        zstat_idx = 2 + 2 * i

        # Save coefficient map
        coef_data = data[..., coef_idx]
        coef_map = opj(output_folder, f'{contrast_name}_coef.nii.gz')
        print(coef_map)
        nib.save(nib.Nifti1Image(coef_data, affine), coef_map)

        # Save Z-stat map
        zstat_data = data[..., zstat_idx]
        zstat_map = opj(output_folder, f'{contrast_name}_zstat.nii.gz')
        nib.save(nib.Nifti1Image(zstat_data, affine), zstat_map)


    # Process each contrast
    for i, contrast_name in enumerate(contrast_names):
        # Calculate indices for this contrast
        coef_idx = 1 + 2 * i
        zstat_idx = 2 + 2 * i

        # Save coefficient map
        coef_data = data[..., coef_idx]
        coef_map = opj(output_folder, f'{contrast_name}_coef.nii.gz')
        print(coef_map)
        nib.save(nib.Nifti1Image(coef_data, affine), coef_map)

        # Save Z-stat map
        zstat_data = data[..., zstat_idx]
        zstat_map = opj(output_folder, f'{contrast_name}_zstat.nii.gz')
        nib.save(nib.Nifti1Image(zstat_data, affine), zstat_map)

        zstat_data = data[..., zstat_idx]
        # Separate positive and negative values
        pos_values = zstat_data[zstat_data > 0]
        neg_values = zstat_data[zstat_data < 0]

        # Initialize thresholds to extremes
        pos_threshold = np.inf  # Will exclude everything if no pos values
        neg_threshold = -np.inf  # Will exclude everything if no neg values

        # Calculate percentiles only if values exist
        if len(pos_values) > 0:
            pos_threshold = np.percentile(pos_values, 100 - percent)

        if len(neg_values) > 0:
            neg_threshold = np.percentile(neg_values, percent)

        # Apply thresholds
        thr_data = np.zeros_like(zstat_data)
        thr_data[(zstat_data >= pos_threshold) | (zstat_data <= neg_threshold)] = zstat_data[
            (zstat_data >= pos_threshold) | (zstat_data <= neg_threshold)]

        # Generate appropriate title based on what exists
        if len(pos_values) > 100 and len(neg_values) > 100:
            title = f"{contrast_name} (Top/Bottom {percent}%)"
            cmap = 'cold_hot'
        elif len(pos_values) > 100:
            title = f"{contrast_name} (Top {percent}%)"
            cmap = 'hot'
        elif len(neg_values) > 100:
            title = f"{contrast_name} (Bottom {percent}%)"
            cmap = 'winter'
        else:
            print(f"Skipping {contrast_name} - no non-zero values")
            continue

        # Visualization
        plotting.plot_stat_map(
            nib.Nifti1Image(thr_data, affine),
            bg_img=studytemplatebrain,
            display_mode='mosaic',
            cut_coords=cut_coords,
            title=title,
            cmap=cmap,
            colorbar=True
        ).savefig(opj(output_folder, f'{contrast_name}_thr_pct.png'), dpi=300)
        plt.close()

def _3dLMEr_EDNiX(bids_dir, templatehigh, templatelow, oversample_map, mask_func, cut_coords,
                  panda_files, selected_atlases, lower_cutoff, upper_cutoff, MAIN_PATH, FS_dir, alpha,
                  all_ID, all_Session, all_data_path, endfmri, mean_imgs,
                  ntimepoint_treshold, model, glt_spec, contrast_names, midline_x, visualize, percent):
    """
    Complete seed-based analysis with proper hemisphere mirroring

    For left seeds: Uses right-hemisphere mirrored version
    For right seeds: Uses left-hemisphere mirrored version
    """

    # Load requirements and setup output
    s_path, afni_sif, fsl_sif, fs_sif, itk_sif, wb_sif, strip_sif, s_bind = Tools.Load_EDNiX_requirement.load_requirement(
        MAIN_PATH, bids_dir, FS_dir)
    output_results1 = opj(bids_dir, 'Results')
    os.makedirs(output_results1, exist_ok=True)

    # Template selection
    studytemplatebrain = templatehigh if oversample_map else templatelow

    # Create group mask
    mean_imgs_rs = nilearn.image.concat_imgs(mean_imgs)
    mask_img = compute_epi_mask(mean_imgs_rs,
                                lower_cutoff=lower_cutoff,
                                upper_cutoff=upper_cutoff)
    mask_img.to_filename(opj(output_results1, 'mask_mean_func.nii.gz'))

    # Process mask with AFNI
    commands = [
        f"3dmask_tool -overwrite -input {opj(output_results1, 'mask_mean_func.nii.gz')} "
        f"-prefix {opj(output_results1, 'mask_mean_func.nii.gz')} -fill_holes",

        f"3dresample -overwrite -master {opj(output_results1, 'mask_mean_func.nii.gz')} "
        f"-input {mask_func} -prefix {opj(output_results1, 'mask_mean_func_orig.nii.gz')}",

        f"3dcalc -overwrite -a {opj(output_results1, 'mask_mean_func.nii.gz')} "
        f"-b {opj(output_results1, 'mask_mean_func_orig.nii.gz')} "
        f"-expr 'a*b' -prefix {opj(output_results1, 'mask_mean_func_overlapp.nii.gz')}"
    ]

    for cmd in commands:
        subprocess.run(f"singularity run {s_bind} {afni_sif} {cmd}", shell=True, check=True)

    # Main analysis loop
    for panda_file, atlas in zip(panda_files, selected_atlases):
        output_results = opj(output_results1, 'Grp_SBA_3dLME_network')
        os.makedirs(output_results, exist_ok=True)

        # Create seed pairs dictionary (only store one direction)
        seed_pairs = {}
        for _, row in panda_file.iterrows():
            seed_name = row['region']
            if seed_name.startswith('L_'):
                right_counterpart = 'R_' + seed_name[2:]
                if right_counterpart in panda_file['region'].values:
                    seed_pairs[seed_name] = right_counterpart  # Only store L->R mapping

        # Process each unique seed pair only once
        for key, value in seed_pairs.items():
            key = format_seed_name(key)
            value = format_seed_name(value)
            base_seed_name = format_seed_name_Group(key)
            print(base_seed_name)

            output_folder = opj(output_results, base_seed_name)
            os.makedirs(output_folder, exist_ok=True)

            # Build design matrix
            design_data = []
            master_file = None

            # Process both original seed and its counterpart if exists
            seeds_to_process = [key, value]

            for current_seed in seeds_to_process:
                current_is_left = current_seed.startswith('L_')
                base_seed_name = format_seed_name(current_seed)  # e.g. "OB" for both L_OB and R_OB
                for ID, Session, data_path in zip(all_ID, all_Session, all_data_path):
                    # Get input directory for current seed (either L or R version)
                    input_dir = opj(data_path, 'func', '01_prepro', '03_atlas_space', '10_Results', 'SBA', current_seed)

                    if not ope(input_dir):
                        continue

                    # Process each run
                    for run_idx, run_path in enumerate(sorted(glob.glob(opj(data_path, 'func', endfmri)))):
                        img = nib.load(run_path)
                    if len(img.shape) == 4 and img.shape[3] < ntimepoint_treshold:
                        continue

                    root_RS = extract_filename(os.path.basename(run_path))
                    original_file = opj(input_dir, f"{root_RS}_correlations_fish.nii.gz")

                    if not ope(original_file):
                        continue

                    # Handle resampling if needed
                    if master_file is None:
                        master_file = original_file
                    elif not check_grid_match(master_file, original_file):
                        original_file = resample_to_match(master_file, original_file, afni_sif, s_bind)

                    # Create mirrored version
                    mirrored_file = mirror_hemisphere(original_file, input_dir, midline_x, current_is_left)

                    # Determine which hemisphere we're analyzing
                    analyzed_hemisphere = 'L' if current_is_left else 'R'

                    # Add to design matrix
                    design_data.append({
                        'Subj': str(ID),
                        'Sess': f'Sess_{Session}',
                        'run': f'run_{run_idx}',
                        'Hemisphere': analyzed_hemisphere,
                        'InputFile': mirrored_file
                    })

            # Skip if no data
            if not design_data:
                print(f"No valid data for seed {base_seed_name}")
                continue

            # Save design matrix with proper column order
            df = pd.DataFrame(design_data)[['Subj', 'Sess', 'run', 'Hemisphere', 'InputFile']]
            design_matrix_path = opj(output_folder, 'design_matrix.txt')
            df.to_csv(design_matrix_path, sep='\t', index=False)

            # Run analysis
            mask_path = opj(output_results1, 'mask_mean_func_overlapp.nii.gz')
            stat_maps = run_3dLMEr_analysis(output_folder, mask_path, design_matrix_path, model, glt_spec, afni_sif, s_bind)

            # Cluster simulation
            command = (
                f"singularity run {s_bind} {afni_sif} 3dClustSim -mask "
                f"{opj(output_results1, 'mask_mean_func_overlapp.nii.gz')} "
                f"-LOTS -prefix {output_folder}/Clust_")
            subprocess.run(command, shell=True)

            if visualize == 'threshold':
                # Visualization
                visualize_results(stat_maps, output_folder, contrast_names, studytemplatebrain, cut_coords, alpha)

            elif visualize == 'percentile':
                visualize_results_percentile(stat_maps, output_folder, contrast_names, studytemplatebrain, cut_coords,
                                             percent)