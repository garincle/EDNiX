import os
import nibabel as nib
import json
from math import pi

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname

from Tools import run_cmd
from modalities.fMRI.extract_filename import extract_filename


# =============================================================================
# Module-level helper  (same pattern as _2_coregistration_to_norm.py)
# =============================================================================

def _write_json(nii_path, sources, description, command=''):
    """Write a minimal BIDS-style JSON sidecar next to a NIfTI file."""
    meta = {"Sources": sources, "Description": description}
    if command:
        meta["Command"] = command
    with open(nii_path.replace('.nii.gz', '.json'), 'w') as f:
        json.dump(meta, f, indent=3)


# =============================================================================
# Main function
#
# Called from _2_coregistration_to_norm.correct_img() as:
#
#   correct_img(dir_prepro_orig_process, dir_prepro_fmap,
#               fMRI_runMean_n4Bias, RS, list_map, RS_map,
#               map_idx, run_idx, recordings,
#               overwrite, sing_afni, sing_fsl, topup_file, diary_file)
#
# Outputs consumed by the caller:
#   • {root}_space-func_desc-runMean_fieldmap_rads.nii.gz  (read by FUGUE)
#
# The function is NOT called for 'very_old' recordings (handled upstream).
# =============================================================================

def correct_img(dir_prepro_orig_process, dir_prepro_fmap,
                fMRI_runMean_n4Bias, RS, list_map, RS_map,
                i, r, recordings,
                overwrite, sing_afni, sing_fsl, topup_file, diary_file):

    run_cmd.msg(
        '##  Working on step 2a (function: _2a_correct_img).  ##',
        diary_file, 'HEADER')

    # =========================================================================
    # BRANCH A — '2_mapdir': two separate phase-encoding directories
    #
    # Two SE-EPI fieldmap runs (z=0: one PE direction, z=1: opposite).
    # Each is stabilised (mean → volreg → mean), then concatenated in the
    # order [z=1, z=0] = [opposite, forward] as expected by topup's datain.
    # NOTE: the concatenation order must match the rows in topup_file[0].
    # =========================================================================

    if recordings == '2_mapdir':

        aligned_means = []   # will hold [fmri_to_fmap_align_Mean_z1,
                             #            fmri_to_fmap_align_Mean_z0]
                             # (filled in reverse below to match topup order)

        for z in [0, 1]:
            root = extract_filename(RS_map[z])

            fmap_mean               = opj(dir_prepro_fmap, root + f'_fmap_mean{z}.nii.gz')
            fmri_to_fmap_align      = opj(dir_prepro_fmap, root + f'_fmri_to_fmap_align{z}.nii.gz')
            fmri_to_fmap_align_Mean = opj(dir_prepro_fmap, root + f'_fmri_to_fmap_align_Mean{z}.nii.gz')

            # Step 1 — temporal mean of the fieldmap volume
            command = (sing_afni + '3dTstat' + overwrite +
                       ' -mean -prefix ' + fmap_mean +
                       ' ' + list_map[z])
            run_cmd.run(command, diary_file)
            _write_json(fmap_mean,
                        list_map[z],
                        '4D temporal mean (3dTstat, AFNI).',
                        command)
            # NOTE: sidecar saved next to fmap_mean, not in dir_prepro_orig_process
            # (original wrote to orig_process which is inconsistent with the file
            # living in dir_prepro_fmap — fixed here)

            # Step 2 — stabilise volumes by registering each to the mean
            # 3dvolreg is used here for consistency with the original; for
            # small-animal data with sub-mm voxels, ANTs would be more accurate
            # (cf. motion correction in step 1 which uses ANTs BOLDRigid).
            command = (sing_afni + '3dvolreg' + overwrite +
                       ' -verbose -zpad 1' +
                       ' -base ' + fmap_mean +
                       ' -prefix ' + fmri_to_fmap_align +
                       ' -cubic ' + opj(dir_prepro_orig_process, RS_map[z]))
            run_cmd.run(command, diary_file)
            _write_json(fmri_to_fmap_align,
                        [opj(dir_prepro_orig_process, RS_map[z]), fmap_mean],
                        'Fieldmap stabilisation — rigid realignment to mean (3dvolreg, AFNI).',
                        command)

            # Step 3 — mean of the stabilised volumes
            command = (sing_afni + '3dTstat' + overwrite +
                       ' -mean -prefix ' + fmri_to_fmap_align_Mean +
                       ' ' + fmri_to_fmap_align)
            run_cmd.run(command, diary_file)
            _write_json(fmri_to_fmap_align_Mean,
                        fmri_to_fmap_align,
                        '4D temporal mean of stabilised fieldmap (3dTstat, AFNI).',
                        command)

            aligned_means.append(fmri_to_fmap_align_Mean)

        # Concatenate: [z=1 (opposite PE), z=0 (forward PE)] — must match
        # the row order in topup_file[0] (datain acqparams.txt).
        root_concat      = extract_filename(RS_map[0])
        fmri_to_fmap_concat = opj(dir_prepro_fmap,
                                  root_concat + '_fmri_to_fmap_concat.nii.gz')
        command = (sing_afni + '3dTcat' + overwrite +
                   ' -prefix ' + fmri_to_fmap_concat +
                   ' ' + aligned_means[1] +   # z=1 first (opposite PE)
                   ' ' + aligned_means[0])    # z=0 second (forward PE)
        run_cmd.run(command, diary_file)
        _write_json(fmri_to_fmap_concat,
                    [aligned_means[1], aligned_means[0]],
                    'Concatenation of opposite-PE fieldmap means for topup '
                    '(order: opposite-PE first, forward-PE second — must match datain). '
                    '(3dTcat, AFNI).',
                    command)

        # root for downstream outputs (fieldmap files named after RS_map[0])
        root = root_concat

    # =========================================================================
    # BRANCH B — 'new' / 'old': single fieldmap directory
    #
    # One SE-EPI fieldmap run per call (indexed by i).
    # Stabilise, then concatenate with the functional mean as the second volume
    # so topup sees [fieldmap_mean, func_mean] → estimates the warp field
    # between the two PE conditions.
    # =========================================================================

    else:
        root    = extract_filename(RS_map[i])
        root_rs = extract_filename(RS[r])   # available for logging if needed

        fmap_mean               = opj(dir_prepro_fmap, root + '_fmap_mean.nii.gz')
        fmri_to_fmap_align      = opj(dir_prepro_fmap, root + '_fmri_to_fmap_align.nii.gz')
        fmri_to_fmap_align_Mean = opj(dir_prepro_fmap, root + '_fmri_to_fmap_align_Mean.nii.gz')
        fmri_to_fmap_concat     = opj(dir_prepro_fmap, root + '_fmri_to_fmap_concat.nii.gz')

        # Step 1 — temporal mean of the fieldmap volume
        command = (sing_afni + '3dTstat' + overwrite +
                   ' -mean -prefix ' + fmap_mean +
                   ' ' + opj(dir_prepro_orig_process, RS_map[i]))
        run_cmd.run(command, diary_file)
        _write_json(fmap_mean,
                    opj(dir_prepro_orig_process, RS_map[i]),
                    '4D temporal mean (3dTstat, AFNI).',
                    command)

        # Step 2 — stabilise fieldmap volumes
        command = (sing_afni + '3dvolreg' + overwrite +
                   ' -verbose -zpad 1' +
                   ' -base ' + fmap_mean +
                   ' -prefix ' + fmri_to_fmap_align +
                   ' -cubic ' + opj(dir_prepro_orig_process, RS_map[i]))
        run_cmd.run(command, diary_file)
        _write_json(fmri_to_fmap_align,
                    [opj(dir_prepro_orig_process, RS_map[i]), fmap_mean],
                    'Fieldmap stabilisation — rigid realignment to mean (3dvolreg, AFNI).',
                    command)

        # Step 3 — mean of stabilised fieldmap
        command = (sing_afni + '3dTstat' + overwrite +
                   ' -mean -prefix ' + fmri_to_fmap_align_Mean +
                   ' ' + fmri_to_fmap_align)
        run_cmd.run(command, diary_file)
        _write_json(fmri_to_fmap_align_Mean,
                    fmri_to_fmap_align,
                    '4D temporal mean of stabilised fieldmap (3dTstat, AFNI).',
                    command)

        # Step 4 — concatenate [fieldmap_mean, func_mean] for topup
        # Order: fieldmap (opposite PE) first, functional mean (forward PE)
        # second — must match datain acqparams rows.
        command = (sing_afni + '3dTcat' + overwrite +
                   ' -prefix ' + fmri_to_fmap_concat +
                   ' ' + fmri_to_fmap_align_Mean +
                   ' ' + fMRI_runMean_n4Bias)
        run_cmd.run(command, diary_file)
        _write_json(fmri_to_fmap_concat,
                    [fmri_to_fmap_align_Mean, fMRI_runMean_n4Bias],
                    'Concatenation of fieldmap mean + functional mean for topup '
                    '(order must match datain). (3dTcat, AFNI).',
                    command)

    # =========================================================================
    # COMMON — even-dimension enforcement, topup, Hz→rad/s conversion
    # =========================================================================

    fmri_to_fmap_even          = opj(dir_prepro_fmap, root + '_fmri_to_fmap_even.nii.gz')
    fMRI_runMean_fieldmap      = opj(dir_prepro_fmap, root + '_space-func_desc-runMean_fieldmap.nii.gz')
    fMRI_runMean_fieldmap_rads = opj(dir_prepro_fmap, root + '_space-func_desc-runMean_fieldmap_rads.nii.gz')
    fMRI_runMean_fieldmap_mag  = opj(dir_prepro_fmap, root + '_space-func_desc-runMean_fieldmap_mag.nii.gz')
    fMRI_b0_distortion_corrected = opj(dir_prepro_fmap, root + '_b0_distortion_corrected.nii.gz')

    # ---- Even-dimension enforcement -----------------------------------------
    # topup requires even voxel counts in all spatial dimensions.
    # If any dimension is odd, the last slice along that axis is dropped.
    # Header shape is updated AFTER trimming so it reflects the actual data.

    im     = nib.load(fmri_to_fmap_concat)
    imdata = im.get_fdata()

    for ax, size in enumerate(imdata.shape[:3]):
        if size % 2 == 0:
            run_cmd.msg(
                f'INFO: axis {ax} size={size} — even, no trimming needed.',
                diary_file, 'OKGREEN')
        else:
            run_cmd.msg(
                f'INFO: axis {ax} size={size} — odd, dropping last slice.',
                diary_file, 'OKGREEN')
            imdata = imdata.take(range(size - 1), axis=ax)

    # Build header AFTER trimming so shape is correct
    hdr = im.header.copy()
    hdr.set_data_shape(imdata.shape)
    nib.Nifti1Image(imdata, im.affine, hdr).to_filename(fmri_to_fmap_even)

    _write_json(fmri_to_fmap_even,
                fmri_to_fmap_concat,
                'Even-dimension enforcement for topup '
                '(last slice dropped on odd axes). (nibabel)',
                'nibabel.Nifti1Image — axis trimming in numpy')

    # ---- topup --------------------------------------------------------------
    # Estimates the susceptibility-induced off-resonance field from the
    # concatenated opposite-PE pair. Outputs:
    #   fMRI_runMean_fieldmap      — field in Hz
    #   fMRI_b0_distortion_corrected — corrected b0 image for QC
    command = (sing_fsl + 'topup' +
               ' --imain='   + fmri_to_fmap_even +
               ' --datain='  + topup_file[0] +
               ' --config='  + topup_file[1] +
               ' --fout='    + fMRI_runMean_fieldmap +
               ' --iout='    + fMRI_b0_distortion_corrected)
    run_cmd.run(command, diary_file)
    _write_json(fMRI_runMean_fieldmap,
                [fmri_to_fmap_even, topup_file[0], topup_file[1]],
                'Susceptibility fieldmap estimation (FSL topup). Output in Hz.',
                command)

    # ---- Hz → rad/s conversion ----------------------------------------------
    # FUGUE expects the fieldmap in rad/s.
    # fslmaths multiplies voxel-wise by 2π.
    command = (sing_fsl + 'fslmaths ' + fMRI_runMean_fieldmap +
               ' -mul ' + str(2 * pi) +
               ' ' + fMRI_runMean_fieldmap_rads)
    run_cmd.run(command, diary_file)
    _write_json(fMRI_runMean_fieldmap_rads,
                fMRI_runMean_fieldmap,
                'Fieldmap unit conversion Hz → rad/s (×2π) for FUGUE. (fslmaths, FSL).',
                command)

    # ---- Magnitude image (for QC / brain masking) ---------------------------
    command = (sing_fsl + 'fslmaths ' + fMRI_b0_distortion_corrected +
               ' -Tmean ' + fMRI_runMean_fieldmap_mag)
    run_cmd.run(command, diary_file)
    _write_json(fMRI_runMean_fieldmap_mag,
                fMRI_b0_distortion_corrected,
                'Temporal mean of topup-corrected b0 image — magnitude for QC. '
                '(fslmaths, FSL).',
                command)
