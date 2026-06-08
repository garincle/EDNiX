import os
import ants
import json
import numpy as np
import nibabel as nib
import subprocess
from pathlib import Path

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile
spgo = subprocess.getoutput

from Tools import run_cmd, diaryfile, check_nii
from modalities.fMRI.extract_filename import extract_filename
from modalities.fMRI import _2b_fix_orient


# =============================================================================
# Module-level helpers
# =============================================================================

def _ants_affine_to_motion_params(affine):
    """
    Convert ANTs affine parameters to (dx, dy, dz, rot_x_deg, rot_y_deg, rot_z_deg).

    ANTsPy BOLDRigid stores a 12-element parameter vector:
      [r00, r10, r20, r01, r11, r21, r02, r12, r22, tx, ty, tz]
    (column-major rotation matrix followed by translation).

    Returns
    -------
    tuple : (dx, dy, dz, rot_x_deg, rot_y_deg, rot_z_deg)
        Translations in mm, rotations in degrees.
    """
    dx, dy, dz = affine[9], affine[10], affine[11]

    rot_x = np.arcsin(np.clip(affine[6], -1.0, 1.0))
    cos_rx = np.cos(rot_x)
    if abs(cos_rx) > 1e-6:
        rot_y = np.arctan2(affine[7] / cos_rx, affine[8] / cos_rx)
        rot_z = np.arctan2(affine[3] / cos_rx, affine[0] / cos_rx)
    else:
        # Gimbal lock fallback
        rot_y = np.arctan2(-affine[2], affine[1])
        rot_z = 0.0

    return dx, dy, dz, np.degrees(rot_x), np.degrees(rot_y), np.degrees(rot_z)


def _extract_motion_params(matrices_dir, output_1D, output_matrix=None):
    """
    Extract motion parameters from ANTs *_0GenericAffine.mat files.

    Parameters
    ----------
    matrices_dir : str
        Directory containing *_0GenericAffine.mat files (one per volume).
    output_1D : str
        Output path for motion parameters (N_vols × 6: dx dy dz rx ry rz).
        Translations in mm, rotations in degrees.
    output_matrix : str, optional
        Output path for flattened 4×3 affine matrices (N_vols × 12), AFNI format.

    Returns
    -------
    motion_params : np.ndarray, shape (N_vols, 6)
    """
    mat_files = sorted(Path(matrices_dir).glob('*_0GenericAffine.mat'))
    if not mat_files:
        raise FileNotFoundError(f"No *_0GenericAffine.mat files found in {matrices_dir}")

    motion_params = []
    all_matrices = []

    for mat in mat_files:
        transform = ants.read_transform(str(mat))
        affine = transform.parameters
        dx, dy, dz, rx, ry, rz = _ants_affine_to_motion_params(affine)
        motion_params.append([dx, dy, dz, rx, ry, rz])

        if output_matrix is not None:
            ants_aff = affine.reshape(4, 3).T
            all_matrices.append(ants_aff.reshape(-1))

    motion_params = np.array(motion_params)
    np.savetxt(output_1D, motion_params, fmt="%.6f", delimiter=" ")
    print(f"[INFO] {len(motion_params)} volumes — motion params saved (dx dy dz rx_deg ry_deg rz_deg)")

    if output_matrix is not None:
        np.savetxt(output_matrix, np.array(all_matrices), fmt="%.6f", delimiter=" ")
        print(f"[INFO] Affine matrices saved: {output_matrix}")

    return motion_params


def _compute_adaptive_fd(motion_params, brain_volume_mm3):
    """
    Compute framewise displacement with brain-size-adapted rotation scaling.

    Rotations (degrees) are converted to mm of cortical surface displacement
    using the brain radius estimated from brain volume (sphere approximation),
    following the logic of Power et al. 2012 but with a species-specific radius
    instead of the fixed 50 mm human value used by fMRIPrep.

    Parameters
    ----------
    motion_params : np.ndarray, shape (N_vols, 6)
        Columns: dx dy dz rx_deg ry_deg rz_deg
    brain_volume_mm3 : float
        Brain volume in mm³ (from adaptive N4 block).

    Returns
    -------
    fd : np.ndarray, shape (N_vols - 1,)
        Framewise displacement in mm.
    r_mm : float
        Brain radius (mm) used for rotation scaling.
    """
    # Brain radius from volume assuming sphere
    r_mm = (3.0 * brain_volume_mm3 / (4.0 * np.pi)) ** (1.0 / 3.0)
    print(f"[INFO] Brain radius for FD scaling: {r_mm:.1f} mm "
          f"(volume={brain_volume_mm3:.0f} mm³)")

    d = np.diff(motion_params, axis=0)           # shape (N-1, 6)
    trans_mm  = d[:, :3]                         # already in mm
    rot_mm    = np.deg2rad(d[:, 3:]) * r_mm      # degrees → mm at surface

    fd = np.sqrt(np.sum(trans_mm ** 2 + rot_mm ** 2, axis=1))
    return fd, r_mm


def _write_censor_files(fd, threshold, n_vols, censore1D, censoretxt,
                        censor_prev=True):
    """
    Write AFNI-compatible censor files from precomputed FD vector.

    AFNI convention: 1 = keep, 0 = censor.

    Parameters
    ----------
    fd : np.ndarray, shape (N_vols - 1,)
    threshold : float
        FD threshold in mm.
    n_vols : int
        Total number of volumes (fd has n_vols - 1 values).
    censore1D : str
        Output path for per-volume censor column (.1D).
    censoretxt : str
        Output path for CENSORTR string (.txt).
    censor_prev : bool
        Also censor the TR preceding each bad volume.

    Returns
    -------
    censor : np.ndarray of int, shape (n_vols,)
    n_censored : int
    """
    censor = np.ones(n_vols, dtype=int)

    bad = np.where(fd > threshold)[0] + 1   # derivative shifts index by 1
    censor[bad] = 0
    if censor_prev:
        prev = bad - 1
        censor[prev[prev >= 0]] = 0

    np.savetxt(censore1D, censor, fmt='%d')

    censored_trs = np.where(censor == 0)[0]
    censor_str = ('0:' + ','.join(map(str, censored_trs))
                  if len(censored_trs) > 0 else '')
    with open(censoretxt, 'w') as f:
        f.write(censor_str + '\n')

    n_censored = int(np.sum(censor == 0))
    print(f"[INFO] Censoring: {n_censored}/{n_vols} volumes removed "
          f"(FD threshold = {threshold:.3f} mm)")
    return censor, n_censored


def _adaptive_n4_correction(input_image, output_path, mask=None,
                             auto_tune=True, species_info=None):
    """
    Adaptive N4 bias field correction for multi-species data.

    Automatically adjusts spline parameter, shrink factor, and convergence
    iterations based on brain volume, covering the ~1000× size range from
    mouse (~500 mm³) to human (~1,200,000 mm³).

    Parameters
    ----------
    input_image : str or ANTsImage
    output_path : str
    mask : str or ANTsImage, optional
        If None, an Otsu-based mask is generated automatically.
    auto_tune : bool
        Adjust N4 parameters from brain volume when True.
    species_info : dict, optional
        Extra fields merged into the JSON sidecar.

    Returns
    -------
    N4 : ANTsImage
        Bias-field-corrected image.
    brain_volume_mm3 : float
        Brain volume used for parameter selection (useful downstream for FD).
    metadata : dict
    """
    IMG = ants.image_read(input_image) if isinstance(input_image, str) else input_image

    voxel_volume   = float(np.prod(IMG.spacing))
    image_shape    = list(IMG.shape)

    # ---- Mask ---------------------------------------------------------------
    if mask is None:
        img_array = IMG.numpy()
        threshold = np.mean(img_array[img_array > 0]) * 0.15
        mask_img  = ants.get_mask(IMG, low_thresh=threshold)
    elif isinstance(mask, str):
        mask_img = ants.image_read(mask)
    else:
        mask_img = mask

    brain_voxels      = int(np.sum(mask_img.numpy() > 0))
    brain_volume_mm3  = brain_voxels * voxel_volume
    print(f"[INFO] N4 — brain volume: {brain_volume_mm3:.0f} mm³  "
          f"voxel spacing: {IMG.spacing}")

    # ---- Parameter selection ------------------------------------------------
    if auto_tune:
        if brain_volume_mm3 < 1_000:
            spline_param, shrink_factor = 50,  2
            convergence_iters = [50, 50, 30, 20]
            species_category  = "small (<1 000 mm³, e.g. mouse)"
        elif brain_volume_mm3 < 5_000:
            spline_param, shrink_factor = 100, 3
            convergence_iters = [50, 50, 40, 30]
            species_category  = "medium (<5 000 mm³, e.g. rat/marmoset)"
        elif brain_volume_mm3 < 100_000:
            spline_param, shrink_factor = 150, 3
            convergence_iters = [50, 50, 50, 40]
            species_category  = "large (<100 000 mm³, e.g. macaque/cat)"
        else:
            spline_param, shrink_factor = 200, 4
            convergence_iters = [50, 50, 50, 50]
            species_category  = "very large (≥100 000 mm³, e.g. human)"
        print(f"[INFO] N4 auto-tuned: {species_category}  "
              f"spline={spline_param}  shrink={shrink_factor}")
    else:
        spline_param, shrink_factor = 100, 3
        convergence_iters = [50, 50, 40, 30]
        species_category  = "default"

    # ---- N4 -----------------------------------------------------------------
    N4 = ants.n4_bias_field_correction(
        IMG,
        mask=mask_img,
        shrink_factor=shrink_factor,
        convergence={'iters': convergence_iters, 'tol': 1e-07},
        spline_param=spline_param,
        verbose=True)

    ants.image_write(N4, output_path, ri=False)

    # Save bias field for QC
    bias_field_path = output_path.replace('.nii.gz', '_biasfield.nii.gz')
    ants.image_write(IMG / (N4 + 1e-10), bias_field_path, ri=False)

    # ---- Sidecar ------------------------------------------------------------
    metadata = {
        "Sources": input_image if isinstance(input_image, str) else "ANTsImage",
        "Description": "Bias field correction (N4) with adaptive parameters",
        "BrainVolume_mm3": brain_volume_mm3,
        "BrainVolume_voxels": brain_voxels,
        "VoxelSpacing": list(IMG.spacing),
        "ImageShape": image_shape,
        "SpeciesCategory": species_category,
        "Parameters": {
            "shrink_factor": shrink_factor,
            "convergence_iters": convergence_iters,
            "convergence_tol": 1e-07,
            "spline_param": spline_param,
            "mask_used": True,
            "auto_tuned": auto_tune},
        "BiasField": bias_field_path}
    if species_info:
        metadata.update(species_info)
    with open(output_path.replace('.nii.gz', '.json'), 'w') as f:
        json.dump(metadata, f, indent=3)

    return N4, brain_volume_mm3, metadata


def _write_json(path, sources, description, command):
    """Write a minimal BIDS-style JSON sidecar."""
    with open(path.replace('.nii.gz', '.json'), 'w') as f:
        json.dump({"Sources": sources,
                   "Description": description,
                   "Command": command}, f, indent=3)


# =============================================================================
# Main preprocessing function
# =============================================================================

def preprocess_data(dir_prepro_raw_process, RS, list_RS, nb_run, T1_eq, TR,
                    Slice_timing_info, dir_prepro_raw_matrices, n_for_ANTS,
                    overwrite, sing_afni, diary_file, animalPosition,
                    humanPosition, orientation, doWARPonfunc, diary_WARNING):

    run_cmd.msg('##  Working on step 1 (function: _1_fMRI_preTTT_in_fMRIspace).  ##',
                diary_file, 'HEADER')

    for i in range(int(nb_run)):

        run_cmd.msg(f'Working on {dir_prepro_raw_process}  run {i + 1}',
                    diary_file, 'OKGREEN')

        root_RS  = extract_filename(RS[i])
        raw_func = list_RS[i]

        # ---- Output paths ---------------------------------------------------
        base_fMRI_targeted  = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-vol_rmv_from_txt.nii.gz')
        base_fMRI           = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-vol_rmv.nii.gz')
        fMRI_despike        = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-despiked.nii.gz')
        fMRI_SliceT         = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-SliceTfixed.nii.gz')
        fMRI_runMean        = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-runMean.nii.gz')
        fMRI_stc            = opj(dir_prepro_raw_process, 'stc.txt')
        fMRI_outcount       = opj(dir_prepro_raw_process, root_RS + f'_space-func_desc-outcount_run{i}.1D')

        outpuprefix_motion_folder = opj(dir_prepro_raw_matrices, root_RS + '_space-func_desc-motion_correction')
        outpuprefix_motion        = opj(outpuprefix_motion_folder, 'img')
        mat_files_pattern         = opj(outpuprefix_motion_folder, 'img*.mat')

        fMRI_run_motion_corrected = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-motion_corrected.nii.gz')
        file_motion_correction    = opj(dir_prepro_raw_matrices, root_RS + '_space-func_desc-motion_correction.1D')
        matrix_motion_params      = opj(dir_prepro_raw_matrices, root_RS + '_space-func_desc-matrix_motion_params.1D')

        censore1D    = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-censor.1D')
        censoretxt   = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-censor.txt')
        demean       = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-demean.1D')
        deriv        = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-deriv.1D')
        motion_enorm = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-motion_enorm.1D')

        fMRI_runMean_align  = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-runMean_align.nii.gz')
        fMRI_runMean_n4Bias = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-runMean_n4Bias.nii.gz')
        fMRI_BASE           = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_BASE.nii.gz')
        fMRI_BASE_Mean      = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_BASE_Mean.nii.gz')

        # ---- 1. Remove bad volumes from text file (optional) ----------------
        txt_path = opj(opd(list_RS[i]), root_RS + '.txt')
        if ope(txt_path):
            with open(txt_path, 'r') as fh:
                cut_low  = int(fh.readline().strip())
                cut_high = int(fh.readline().strip())
            command = (sing_afni + '3dTcat -prefix ' + base_fMRI_targeted +
                       ' ' + raw_func + f'[{cut_low}-{cut_high - 1}]' + overwrite)
            _write_json(base_fMRI_targeted,
                        [base_fMRI_targeted, txt_path],
                        'Remove volumes listed in text file.', command)
            raw_func = base_fMRI_targeted
            run_cmd.do(command, diary_file)

        # ---- 2. Remove T1-equilibration volumes -----------------------------
        nb_vol = nib.load(raw_func).shape[-1]
        command = (sing_afni + '3dTcat -prefix ' + base_fMRI +
                   ' ' + raw_func + f'[{T1_eq}-{nb_vol - 1}]' + overwrite)
        _write_json(base_fMRI, base_fMRI, 'Remove first T1-equilibration volumes.', command)
        run_cmd.do(command, diary_file)

        # ---- 3. Despiking ---------------------------------------------------
        command = (sing_afni + '3dDespike -NEW -nomask' + overwrite +
                   ' -prefix ' + fMRI_despike + ' ' + base_fMRI)
        _write_json(fMRI_despike, base_fMRI, 'Despiking.', command)
        run_cmd.run(command, diary_file)

        # ---- 4. Slice-timing correction -------------------------------------
        stc_applied = False

        if Slice_timing_info == 'Auto':
            if opi(fMRI_stc):
                command = (sing_afni + '3dTshift -wsinc9' + overwrite +
                           f' -TR {TR} -tpattern @{fMRI_stc}' +
                           ' -prefix ' + fMRI_SliceT + ' ' + fMRI_despike)
                stc_applied = True
            else:
                nl  = ('export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";' +
                       sing_afni + '3dinfo -slice_timing ' + list_RS[0])
                STC = list(map(float, spgo(nl).split('\n')[-1].split('|')))
                if np.sum(STC) > 0:
                    run_cmd.msg(f'INFO: SliceTiming = {STC}', diary_file, 'OKGREEN')
                    command = (sing_afni + '3dTshift -wsinc9' + overwrite +
                               ' -prefix ' + fMRI_SliceT + ' ' + fMRI_despike)
                    stc_applied = True
                else:
                    run_cmd.msg(
                        'WARNING: Slice timing not found — copying without STC. '
                        'You SHOULD provide slice timing manually!',
                        diary_file, 'WARNING')
                    diaryfile.create(diary_WARNING,
                                     'WARNING: Slice timing missing, STC skipped.')
                    command = (sing_afni + '3dcalc -a ' + fMRI_despike +
                               ' -prefix ' + fMRI_SliceT +
                               ' -expr "a"' + overwrite)
                    run_cmd.do(command, diary_file)

        elif (not isinstance(Slice_timing_info, list) and
              Slice_timing_info.split(' ')[0] == '-tpattern'):
            command = (sing_afni + '3dTshift -wsinc9' + overwrite +
                       f' -TR {TR} {Slice_timing_info}' +
                       ' -prefix ' + fMRI_SliceT + ' ' + fMRI_despike)
            stc_applied = True

        elif isinstance(Slice_timing_info, list):
            command = (sing_afni + '3dTshift -wsinc9' + overwrite +
                       f' -TR {TR} -tpattern @{fMRI_stc}' +
                       ' -prefix ' + fMRI_SliceT + ' ' + fMRI_despike)
            stc_applied = True

        else:
            raise ValueError(run_cmd.error(
                'ERROR: Slice_timing_info is not defined as expected.', diary_file))

        if stc_applied:
            _write_json(fMRI_SliceT, fMRI_despike, 'Slice-timing correction.', command)
            run_cmd.run(command, diary_file)

        # ---- 5. Run mean (pre-motion) ---------------------------------------
        command = (sing_afni + '3dTstat' + overwrite +
                   ' -mean -prefix ' + fMRI_runMean + ' ' + fMRI_SliceT)
        _write_json(fMRI_runMean, fMRI_SliceT, 'Temporal mean image.', command)
        run_cmd.run(command, diary_file)

        # ---- 6. Outlier count per volume ------------------------------------
        command = (sing_afni + '3dToutcount' + overwrite +
                   ' -automask -fraction -polort 4 -legendre ' +
                   fMRI_SliceT + ' > ' + fMRI_outcount)
        subprocess.run(command, shell=True, check=True)

        # ---- 7. Orientation fix (pre-motion, sets fMRI_BASE) ----------------
        print(animalPosition)
        _2b_fix_orient.fix_orient(fMRI_BASE, fMRI_SliceT, list_RS,
                                  animalPosition, humanPosition, orientation,
                                  doWARPonfunc, sing_afni, diary_file)

        # ---- 8. Motion correction (ANTs BOLDRigid) --------------------------
        os.makedirs(outpuprefix_motion_folder, exist_ok=True)

        motion_result = ants.motion_correction(
            image=ants.image_read(fMRI_SliceT),
            fixed=ants.image_read(fMRI_runMean),
            verbose=True,
            type_of_transform='BOLDRigid',
            interpolator=n_for_ANTS,
            outprefix=outpuprefix_motion)

        motion_result['motion_corrected'].to_filename(fMRI_run_motion_corrected)
        check_nii.keep_header(fMRI_run_motion_corrected, fMRI_SliceT)

        _write_json(fMRI_run_motion_corrected,
                    [fMRI_SliceT, fMRI_runMean],
                    'Rigid realignment (ANTs BOLDRigid).',
                    'ants.motion_correction — BOLDRigid, fixed=runMean')

        # ---- 9. Extract motion parameters from ANTs .mat files --------------
        motion_params = _extract_motion_params(
            matrices_dir=outpuprefix_motion_folder,
            output_1D=file_motion_correction,
            output_matrix=matrix_motion_params)

        # ---- 10. Adaptive FD censoring with brain-size-scaled rotations -----
        #
        # Unlike fMRIPrep (fixed r=50 mm, humans only), rotations are converted
        # to mm of cortical surface displacement using the species-specific brain
        # radius estimated from brain volume. This is critical for small-brained
        # species: a 1° rotation in a mouse (r≈5 mm) produces ~0.09 mm
        # displacement vs ~0.87 mm in a human (r≈50 mm).
        #
        # FD threshold = 2 × mean voxel size (species-adaptive).
        # NOTE: brain_volume_mm3 is retrieved from the N4 step (step 12 below).
        # We first compute mean voxel size from the motion-corrected image,
        # then defer brain_volume_mm3 until after N4 has run.
        # To keep the pipeline linear we do a lightweight mask here instead.

        _img_tmp     = nib.load(fMRI_run_motion_corrected)
        _zooms       = _img_tmp.header.get_zooms()[:3]
        _mean_voxel  = float(np.mean(_zooms))
        _fd_threshold = round(2.0 * _mean_voxel, 3)

        run_cmd.msg(
            f'INFO: mean voxel size = {_mean_voxel:.3f} mm  →  '
            f'FD threshold = {_fd_threshold:.3f} mm',
            diary_file, 'OKGREEN')

        # Brain volume estimate for rotation scaling (lightweight Otsu mask)
        _ants_tmp       = ants.image_read(fMRI_runMean)   # use mean for stability
        _arr            = _ants_tmp.numpy()
        _thr_quick      = np.mean(_arr[_arr > 0]) * 0.15
        _mask_quick     = (_arr > _thr_quick).astype(float)
        _brain_vol_quick = float(np.sum(_mask_quick)) * float(np.prod(_ants_tmp.spacing))

        fd, r_mm = _compute_adaptive_fd(motion_params, _brain_vol_quick)

        n_vols = motion_params.shape[0]
        censor, n_censored = _write_censor_files(
            fd, _fd_threshold, n_vols, censore1D, censoretxt, censor_prev=True)

        run_cmd.msg(
            f'INFO: FD censoring — brain radius {r_mm:.1f} mm, '
            f'{n_censored}/{n_vols} volumes censored',
            diary_file, 'OKGREEN')

        # ---- 11. Motion regressors (demean, derivative, enorm) --------------
        # enorm: raw Euclidean norm for QC visualisation only (units not critical)
        command = (sing_afni + '1d_tool.py' + overwrite +
                   ' -infile ' + file_motion_correction +
                   ' -set_nruns 1 -derivative -collapse_cols euclidean_norm' +
                   ' -write ' + motion_enorm)
        run_cmd.run(command, diary_file)

        command = (sing_afni + '1d_tool.py' + overwrite +
                   ' -infile ' + file_motion_correction +
                   ' -derivative -write ' + deriv)
        run_cmd.run(command, diary_file)

        command = (sing_afni + '1d_tool.py' + overwrite +
                   ' -infile ' + file_motion_correction +
                   ' -demean -write ' + demean)
        run_cmd.run(command, diary_file)

        # ---- 12. Run mean of motion-corrected data --------------------------
        command = (sing_afni + '3dTstat' + overwrite +
                   ' -mean -prefix ' + fMRI_runMean_align +
                   ' ' + fMRI_run_motion_corrected)
        _write_json(fMRI_runMean_align, fMRI_run_motion_corrected,
                    'Temporal mean of motion-corrected data.', command)
        run_cmd.run(command, diary_file)

        # ---- 13. Adaptive N4 bias field correction --------------------------
        _, brain_volume_mm3, _ = _adaptive_n4_correction(
            fMRI_runMean_align,
            fMRI_runMean_n4Bias,
            mask=None,
            auto_tune=True,
            species_info=None)

        # ---- 14. Orientation fix (post-N4, sets fMRI_BASE_Mean) -------------
        _2b_fix_orient.fix_orient(fMRI_BASE_Mean, fMRI_runMean_n4Bias, list_RS,
                                  animalPosition, humanPosition, orientation,
                                  doWARPonfunc, sing_afni, diary_file)
