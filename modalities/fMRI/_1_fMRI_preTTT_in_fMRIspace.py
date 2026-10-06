import os
import ants
import json
import numpy as np
import nibabel as nib
import subprocess
import shutil
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
from modalities.fMRI import motion


# =============================================================================
# Module-level helpers
# =============================================================================


def _same_physical_frame(img_a, img_b, atol=1e-4):
    """True if both images have the same voxel grid and voxel-to-world affine."""
    a, b = nib.load(img_a), nib.load(img_b)
    return a.shape[:3] == b.shape[:3] and np.allclose(a.affine, b.affine, atol=atol)

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
            # ITK parameters: 3x3 matrix row-major, then translation. Row-major 3x4
            # [R | t] (AFNI aff12 layout). The previous .T transposed R (i.e. R^-1
            # combined with t), an inconsistent affine.
            R = np.asarray(affine[:9]).reshape(3, 3)
            t = np.asarray(affine[9:12]).reshape(3, 1)
            all_matrices.append(np.hstack([R, t]).reshape(-1))

    motion_params = np.array(motion_params)
    np.savetxt(output_1D, motion_params, fmt="%.6f", delimiter=" ")
    print(f"[INFO] {len(motion_params)} volumes — motion params saved (dx dy dz rx_deg ry_deg rz_deg)")

    if output_matrix is not None:
        np.savetxt(output_matrix, np.array(all_matrices), fmt="%.6f", delimiter=" ")
        print(f"[INFO] Affine matrices saved: {output_matrix}")

    return motion_params


def _compute_adaptive_fd(motion_params, brain_volume_mm3,
                         fd_ref_mm=0.5,
                         fd_floor_mm=0.10,
                         fd_cap_mm=0.50,
                         rotation_exponent=0.5):
    """
    Compute framewise displacement with brain-size-adapted rotation scaling
    and a physiologically bounded censoring threshold.

    FD formula (Power et al. 2012):
        FD_t = sqrt( sum_i (trans_i_t)^2 + sum_j (rot_j_t * r_mm)^2 )
    where rot_j_t is in radians and r_mm is the brain radius (sphere
    approximation from brain volume). Translations are already in mm.

    Threshold selection (multi-species):
        A purely linear scaling (threshold ∝ r_mm) is too strict for small
        brains: a rat (r ≈ 7 mm) would receive a threshold of ~0.05 mm,
        below the physiological noise floor (~0.1 mm). We therefore use
        a power-law scaling with a physiologically motivated floor and cap:

            threshold = clip( fd_ref_mm * (r_mm / r_ref_mm)^rotation_exponent,
                              fd_floor_mm,
                              fd_cap_mm )

    with r_ref_mm = radius of a 1200 cm³ human brain (≈ 65.9 mm).

    Parameters
    ----------
    motion_params : np.ndarray, shape (N_vols, 6)
        Columns: dx dy dz rx_deg ry_deg rz_deg
    brain_volume_mm3 : float
        Brain volume in mm³ (from adaptive N4 block).
    fd_ref_mm : float
        Reference FD threshold for a human brain (default 0.5, fMRIPrep).
    fd_floor_mm : float
        Absolute minimum threshold (default 0.10 mm, physiological floor).
    fd_cap_mm : float
        Absolute maximum threshold (default 0.50 mm, human default).
    rotation_exponent : float
        Power-law exponent for species scaling (default 0.5, softer than
        linear; set to 1.0 to recover the old linear behaviour).

    Returns
    -------
    fd : np.ndarray, shape (N_vols - 1,)
        Framewise displacement in mm.
    r_mm : float
        Brain radius (mm) used for rotation scaling.
    fd_threshold : float
        Adaptive censoring threshold in mm, physiologically bounded.
    """
    # ---- Brain radius (sphere approximation) -------------------------------
    r_mm = (3.0 * brain_volume_mm3 / (4.0 * np.pi)) ** (1.0 / 3.0)
    print(f"[INFO] Brain radius for FD scaling: {r_mm:.1f} mm "
          f"(volume={brain_volume_mm3:.0f} mm³)")

    # ---- Framewise displacement (Power et al. 2012) ------------------------
    d = np.diff(motion_params, axis=0)           # shape (N-1, 6)
    trans_mm = d[:, :3]                          # already in mm
    rot_mm   = np.deg2rad(d[:, 3:]) * r_mm       # degrees → mm at surface

    fd = np.sqrt(np.sum(trans_mm ** 2 + rot_mm ** 2, axis=1))

    # ---- Adaptive threshold with physiological bounds ----------------------
    r_ref_mm = (3.0 * 1.2e6 / (4.0 * np.pi)) ** (1.0 / 3.0)   # ≈ 65.9 mm

    fd_threshold = fd_ref_mm * (r_mm / r_ref_mm) ** rotation_exponent
    fd_threshold = float(np.clip(fd_threshold, fd_floor_mm, fd_cap_mm))

    print(f"[INFO] FD threshold: {fd_threshold:.4f} mm "
          f"(ref={fd_ref_mm} mm, exponent={rotation_exponent}, "
          f"floor={fd_floor_mm}, cap={fd_cap_mm})")

    return fd, r_mm, fd_threshold


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
        #print(f"[INFO] N4 auto-tuned: {species_category}  "
        #      f"spline={spline_param}  shrink={shrink_factor}")
    else:
        spline_param, shrink_factor = 100, 3
        convergence_iters = [50, 50, 40, 30]
        species_category  = "default"

    # ---- B-spline distance proportional to brain size ----------------------
    # A scalar spline_param is the B-spline mesh SPACING in mm (ANTsPy/N4).
    # The previous values (50 mm mouse ... 200 mm human) were coarser relative
    # to the brain in small species (50 mm on a ~10 mm mouse brain = no local
    # correction). Keep the ANTs human default (200 mm for ~1200 cm3) and
    # scale it with the cube root of the brain volume, never below 4 shrunken voxels.
    if auto_tune:
        spline_param = 200.0 * (brain_volume_mm3 / 1.2e6) ** (1.0 / 3.0)
        spline_param = float(max(spline_param, 4.0 * max(IMG.spacing[:3]) * shrink_factor))

    # ---- N4 -----------------------------------------------------------------
    N4 = ants.n4_bias_field_correction(
        IMG,
        mask=mask_img,
        shrink_factor=shrink_factor,
        convergence={'iters': convergence_iters, 'tol': 1e-07},
        spline_param=spline_param,
        verbose=False)

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
                    humanPosition, orientation, doWARPonfunc, diary_WARNING,
                    motion_mask_path=None, motion_mask_erode=3, motion_method='ants', sing_fsl=''):

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
        fMRI_outcount       = opj(dir_prepro_raw_process, root_RS + f'_space-func_desc-outcount.1D')

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
                                  doWARPonfunc, sing_afni, diary_file, run_idx=i)

        # ---- 8. Motion correction (ANTs BOLDRigid) --------------------------
        # remove transforms from previous executions (stale img_XXXXX files would be globbed)
        shutil.rmtree(outpuprefix_motion_folder, ignore_errors=True)
        os.makedirs(outpuprefix_motion_folder, exist_ok=True)

        if motion_method == 'ants':
            motion_result = ants.motion_correction(
                image=ants.image_read(fMRI_SliceT),
                fixed=ants.image_read(fMRI_runMean),
                verbose=False,
                type_of_transform='BOLDRigid',
                outprefix=outpuprefix_motion)        # (an `interpolator` argument was ignored by ants.registration)

            motion_result['motion_corrected'].to_filename(fMRI_run_motion_corrected)
            check_nii.keep_header(fMRI_run_motion_corrected, fMRI_SliceT)

            _write_json(fMRI_run_motion_corrected,
                        [fMRI_SliceT, fMRI_runMean],
                        'Rigid realignment (ANTs BOLDRigid).',
                        'ants.motion_correction — BOLDRigid, fixed=runMean')
        else:
            # 3dvolreg / MCFLIRT: their transforms are converted into the ANTs matrices read by step 5 (seems to give worth results but can be usefull if ants failed)
            try:
                _ext = motion.estimate_external(
                    motion_method, fMRI_SliceT, fMRI_runMean, outpuprefix_motion_folder, fMRI_run_motion_corrected,
                    lambda cmd: run_cmd.run(cmd, diary_file), sing_afni, sing_fsl)
            except (ValueError, ImportError, FileNotFoundError) as exc:
                raise Exception(run_cmd.error('ERROR: motion_method=' + motion_method + ' failed for ' + root_RS + ': ' + str(exc) +
                                              " (motion_method='ants' is the default)", diary_file))
            check_nii.keep_header(fMRI_run_motion_corrected, fMRI_SliceT)
            _write_json(fMRI_run_motion_corrected, [fMRI_SliceT, fMRI_runMean],
                        'Rigid realignment (' + motion_method + ').', motion_method + ', base=runMean')
            run_cmd.msg(
                f"INFO: motion estimated by {motion_method}; its matrices are converted for ANTs (convention '{_ext['convention']}'"
                + (', checked against its registered volumes: agreement '
                   + ' / '.join(f'{k} {v:.4f}' for k, v in _ext['scores'].items()) + f" (uncorrected {_ext['uncorrected']:.4f})"
                   if _ext['conclusive'] else
                   f", NOT checkable: the largest motion is {_ext['max_displacement_vox']:.2f} voxel") + ')',
                diary_file, 'OKGREEN' if _ext['conclusive'] else 'WARNING')

        # ---- 8b. Motion matrices in the reoriented frame (fMRI_BASE) ---------
        # The matrices above are expressed in the physical frame of fMRI_SliceT.
        # _5 applies them to fMRI_BASE; when _2b_fix_orient changed the frame
        # (3drefit -orient / -deoblique, 3dWarp), they must be re-estimated in
        # the fMRI_BASE frame, otherwise rotations/translations are applied
        # about the wrong axes.
        outpuprefix_motion_folder_reor = outpuprefix_motion_folder + '_reoriented'
        shutil.rmtree(outpuprefix_motion_folder_reor, ignore_errors=True)
        if _same_physical_frame(fMRI_SliceT, fMRI_BASE):
            run_cmd.msg('INFO: fMRI_BASE shares the fMRI_SliceT frame — motion matrices reused.',
                        diary_file, 'OKGREEN')
        else:
            os.makedirs(outpuprefix_motion_folder_reor, exist_ok=True)
            _base_img = ants.image_read(fMRI_BASE)
            if motion_method == 'ants':
                ants.motion_correction(
                    image=_base_img,
                    fixed=ants.get_average_of_timeseries(_base_img),
                    verbose=False,
                    type_of_transform='BOLDRigid',
                    outprefix=opj(outpuprefix_motion_folder_reor, 'img'))
            else:
                _ref_base = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-runMean_BASE.nii.gz')
                _tmp_corrected = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-motion_corrected_BASE_tmp.nii.gz')
                ants.get_average_of_timeseries(_base_img).to_filename(_ref_base)
                try:
                    motion.estimate_external(motion_method, fMRI_BASE, _ref_base, outpuprefix_motion_folder_reor, _tmp_corrected,
                                             lambda cmd: run_cmd.run(cmd, diary_file), sing_afni, sing_fsl)
                except (ValueError, ImportError, FileNotFoundError) as exc:
                    raise Exception(run_cmd.error('ERROR: motion_method=' + motion_method + ' failed on the reoriented series of ' +
                                                  root_RS + ': ' + str(exc), diary_file))
                for _f in (_ref_base, _tmp_corrected):
                    if os.path.exists(_f):
                        os.remove(_f)
            del _base_img
            run_cmd.msg('INFO: frame changed by the orientation fix — motion matrices '
                        're-estimated on fMRI_BASE in ' + outpuprefix_motion_folder_reor,
                        diary_file, 'OKGREEN')

        # ---- 8c. Optional: re-estimate the motion with a brain mask ----------
        # Whole-image registration follows the bright tissue that does not move with the brain
        # (skull, muscle, skin); a mask INSIDE the brain follows the brain (see motion.py). The
        # matrices of the frame used by step 5 are replaced, and the motion parameters, FD and
        # censoring below are computed from them.
        motion_dir_for_params = outpuprefix_motion_folder
        brain_volume_mm3 = None
        if motion_mask_path and motion_method != 'ants':
            run_cmd.msg('WARNING: motion_mask is only used with motion_method=ants: ignored with ' + motion_method,
                        diary_file, 'WARNING')
        elif motion_mask_path:
            _same_frame = _same_physical_frame(fMRI_SliceT, fMRI_BASE)
            _est_dir = outpuprefix_motion_folder if _same_frame else outpuprefix_motion_folder_reor
            try:
                refined = motion.refine_with_mask(fMRI_BASE, _est_dir, motion_mask_path, motion_mask_erode)
            except ValueError as exc:
                run_cmd.msg('WARNING: motion_mask not used for ' + root_RS + ' (motion estimated without a mask): ' + str(exc),
                            diary_file, 'WARNING')
            else:
                motion_dir_for_params = _est_dir
                brain_volume_mm3 = refined['volume_mm3_full']
                if _same_frame:     # the realigned series gives the mean used by the next steps
                    refined['motion_corrected'].to_filename(fMRI_run_motion_corrected)
                    check_nii.keep_header(fMRI_run_motion_corrected, fMRI_SliceT)
                run_cmd.msg(
                    f"INFO: motion re-estimated inside the brain mask ({refined['voxels_eroded']} voxels, "
                    f"eroded by {motion_mask_erode}; brain {refined['volume_mm3_full']:.0f} mm3). Difference with the "
                    f"estimation without mask, between volumes: median {refined['median']:.3f} mm, 90th percentile "
                    f"{refined['p90']:.3f} mm, maximum {refined['max']:.3f} mm. A difference well below the FD threshold "
                    f"means the mask makes no difference for this run; a large one means that tissue outside the brain "
                    f"was driving the estimation. Unmasked matrices kept in {_est_dir}_unmasked.",
                    diary_file, 'OKGREEN' if refined['median'] < 0.02 else 'WARNING')

        # ---- 9. Extract motion parameters from ANTs .mat files --------------
        motion_params = _extract_motion_params(
            matrices_dir=motion_dir_for_params,
            output_1D=file_motion_correction,
            output_matrix=matrix_motion_params)

        # ---- 10. FD censoring, displacement AND threshold scaled to brain size ----
        #
        # Rotations are converted to mm at the surface of a sphere of the same
        # volume as the brain (radius r). The threshold is scaled by the same
        # radius so that the criterion is the same *relative* displacement in
        # every species: FD_thr = FD_REF_MM * r / R_REF_MM, with FD_REF_MM = 0.5 mm
        # (fMRIPrep default) for a 1200 cm3 human brain (R_REF_MM ~ 65.9 mm).
        _ants_tmp        = ants.image_read(fMRI_runMean)
        _arr             = _ants_tmp.numpy()
        _thr_quick       = np.mean(_arr[_arr > 0]) * 0.15
        _brain_vol_quick = float(np.sum(_arr > _thr_quick)) * float(np.prod(_ants_tmp.spacing))
        if brain_volume_mm3 is not None:       # a brain mask is a better volume than a threshold on the mean image
            _brain_vol_quick = brain_volume_mm3

        fd, r_mm, _fd_threshold = _compute_adaptive_fd(
            motion_params, _brain_vol_quick,
            fd_ref_mm=0.5,  # référence humaine
            fd_floor_mm=0.10,  # plancher physiologique
            fd_cap_mm=0.50,  # plafond humain
            rotation_exponent=0.5)  # loi de puissance douce

        n_vols = motion_params.shape[0]
        censor, n_censored = _write_censor_files(
            fd, _fd_threshold, n_vols, censore1D, censoretxt, censor_prev=True)

        run_cmd.msg(
            f'INFO: FD censoring — brain radius {r_mm:.1f} mm, threshold {_fd_threshold:.4f} mm, '
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
                                  doWARPonfunc, sing_afni, diary_file, run_idx=i)
