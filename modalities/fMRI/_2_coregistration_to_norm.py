import os
import ants
import json

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd
from modalities.fMRI.extract_filename import extract_filename
from modalities.fMRI import _2b_fix_orient, _2a_correct_img
from modalities.fMRI import plot_QC_func


# =============================================================================
# Module-level helper
# =============================================================================

def _write_json(nii_path, sources, description, command=''):
    """Write a minimal BIDS-style JSON sidecar next to a NIfTI file."""
    meta = {"Sources": sources,
            "Description": description}
    if command:
        meta["Command"] = command
    with open(nii_path.replace('.nii.gz', '.json'), 'w') as f:
        json.dump(meta, f, indent=3)


# =============================================================================
# Main coregistration function
# =============================================================================

def coregist_to_norm(correction_direction, list_RS, dir_prepro_fmap,
                     dir_prepro_raw_process, dir_prepro_orig_process,
                     RS, RS_map, nb_run, recordings, bids_dir,
                     REF_int, list_map, animalPosition, humanPosition,
                     doWARPonfunc, dir_prepro_raw_matrices, orientation,
                     DwellT, n_for_ANTS, overwrite, sing_afni, sing_fsl,
                     dmap, dbold, config_f, diary_file):

    run_cmd.msg(
        '##  Working on step 2 (function: _2_coregistration_to_norm).  ##',
        diary_file, 'HEADER')

    # ---- Topup config file --------------------------------------------------
    topup_file_path = opj(dir_prepro_raw_process, '4topup.txt')
    with open(topup_file_path, 'w') as f:
        f.write(dmap + ' \n')
        f.write(dbold + ' \n')
    topup_file = [topup_file_path, config_f]

    if recordings != 'very_old':
        run_cmd.msg(f'INFO: DwellT = {DwellT}  — please verify.',
                    diary_file, 'OKGREEN')

    # =========================================================================
    # LOOP 1 — Distortion correction + orientation fix of the mean image
    #
    # For 'new' recordings: one fieldmap per run → zip(list_map, runs).
    # For 'old' / '2_mapdir': one shared fieldmap → process REF run only.
    # For 'very_old': no fieldmap available → simple copy.
    # =========================================================================

    if recordings == 'new':
        items_to_process = list(zip(range(len(list_map)), range(nb_run)))
    else:
        items_to_process = [(0, REF_int)]

    for map_idx, run_idx in items_to_process:

        run_cmd.msg(
            f'INFO: distortion correction — run {run_idx + 1}/{nb_run}',
            diary_file, 'OKGREEN')

        root_RS = extract_filename(RS[run_idx])

        fMRI_runMean_n4Bias         = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-runMean_n4Bias.nii.gz')
        fMRI_runMean_unwarped       = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-runMean_unwarped.nii.gz')
        runMean_unwarped_reoriented = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-runMean_unwarped_reoriented.nii.gz')

        if recordings == 'very_old':
            # No fieldmap — pass through unchanged
            command = (sing_afni + '3dcopy ' + fMRI_runMean_n4Bias +
                       ' ' + fMRI_runMean_unwarped + overwrite)
            run_cmd.run(command, diary_file)
            _write_json(fMRI_runMean_unwarped,
                        fMRI_runMean_n4Bias,
                        'Copy (no fieldmap available for very_old recordings).',
                        command)

        else:
            # Correct image geometry / intensity before unwarping
            _2a_correct_img.correct_img(
                dir_prepro_orig_process, dir_prepro_fmap,
                fMRI_runMean_n4Bias, RS, list_map, RS_map,
                map_idx, run_idx, recordings,
                overwrite, sing_afni, sing_fsl, topup_file, diary_file)

            root_map = extract_filename(RS_map[map_idx])
            fMRI_runMean_fieldmap_rads = opj(
                dir_prepro_fmap,
                root_map + '_space-func_desc-runMean_fieldmap_rads.nii.gz')

            # FUGUE EPI distortion correction
            # NOTE: FUGUE operates in voxel space; correction_direction must be
            # one of {x, x-, y, y-, z, z-} matching the phase-encoding axis.
            # DwellT (effective echo spacing, seconds) must be verified per
            # acquisition — a wrong value will over- or under-correct.
            command = (sing_fsl + 'fugue'
                       ' -i '          + fMRI_runMean_n4Bias +
                       ' --dwell='     + str(DwellT) +
                       ' --loadfmap='  + fMRI_runMean_fieldmap_rads +
                       ' --unwarpdir=' + correction_direction +
                       ' -u '          + fMRI_runMean_unwarped)
            run_cmd.run(command, diary_file)
            _write_json(fMRI_runMean_unwarped,
                        [fMRI_runMean_n4Bias, fMRI_runMean_fieldmap_rads],
                        'EPI distortion correction (FSL FUGUE).',
                        command)

        # Orientation fix on the unwarped mean
        doWARPonfunc_spe = doWARPonfunc if doWARPonfunc in ['WARP', 'header'] else 'WARP'
        _2b_fix_orient.fix_orient(
            runMean_unwarped_reoriented, fMRI_runMean_unwarped, list_RS,
            animalPosition, humanPosition, orientation,
            doWARPonfunc_spe, sing_afni, diary_file)

    # =========================================================================
    # LOOP 2 — Reorient full 4D functional + mean for every run
    #
    # Kept as a separate loop because orientation fix on the 4D timeseries
    # is independent of distortion correction (which only applies to the mean).
    # =========================================================================

    for run_idx in range(int(nb_run)):

        run_cmd.msg(
            f'INFO: reorientation — run {run_idx + 1}/{nb_run}',
            diary_file, 'OKGREEN')

        root_RS = extract_filename(RS[run_idx])

        fMRI_run_motion_corrected = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-motion_corrected.nii.gz')
        fMRI_reoriented           = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_reoriented.nii.gz')
        fMRI_runMean_n4Bias       = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-runMean_n4Bias.nii.gz')
        fMRI_runMean_reoriented   = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_runMean_reoriented.nii.gz')

        _2b_fix_orient.fix_orient(
            fMRI_reoriented, fMRI_run_motion_corrected, list_RS,
            animalPosition, humanPosition, orientation,
            doWARPonfunc, sing_afni, diary_file)

        _2b_fix_orient.fix_orient(
            fMRI_runMean_reoriented, fMRI_runMean_n4Bias, list_RS,
            animalPosition, humanPosition, orientation,
            doWARPonfunc, sing_afni, diary_file)

    # =========================================================================
    # LOOP 3 — Co-registration of each run mean to the reference run mean
    #
    # Reference = REF_int run, unwarped and reoriented.
    # For non-reference runs: ANTs BOLDRigid registration of the run mean to
    # the reference mean, then the resulting transform is applied to the full
    # 4D timeseries (imagetype=3 for time series).
    # For the reference run itself under very_old + WARP/header: identity copy
    # (the reference is already in its own space by definition).
    # =========================================================================

    root_RS_ref               = extract_filename(RS[REF_int])
    runMean_unwarped_reoriented_Ref = opj(
        dir_prepro_raw_process,
        root_RS_ref + '_space-func_desc-runMean_unwarped_reoriented.nii.gz')
    REF = ants.image_read(runMean_unwarped_reoriented_Ref)

    os.makedirs(opj(bids_dir, 'QC', 'fMRI_runMean_in_REF'), exist_ok=True)

    for run_idx in range(int(nb_run)):

        run_cmd.msg(
            f'INFO: coregistration to reference — run {run_idx + 1}/{nb_run}',
            diary_file, 'OKGREEN')

        root_RS = extract_filename(RS[run_idx])
        is_ref  = (root_RS == root_RS_ref)

        fMRI_reoriented         = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_reoriented.nii.gz')
        fMRI_runMean_reoriented = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_runMean_reoriented.nii.gz')
        fMRI_runMean_inRef      = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_runMean_inRef.nii.gz')
        fMRI_run_inRef_mat      = opj(dir_prepro_raw_matrices, root_RS + '_space-func_desc-fMRI_run_inRef')
        fMRI_run_inRef          = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_run_inRef.nii.gz')

        if is_ref and recordings == 'very_old' and doWARPonfunc in ['WARP', 'header']:
            # Reference run under very_old — already in reference space, copy only.
            for src, dst in [(fMRI_reoriented,         fMRI_run_inRef),
                             (fMRI_runMean_reoriented,  fMRI_runMean_inRef)]:
                command = sing_afni + '3dcopy ' + src + ' ' + dst + overwrite
                run_cmd.run(command, diary_file)
                _write_json(dst, src,
                            'Identity copy — reference run already in reference space.',
                            command)

        else:
            # -----------------------------------------------------------------
            # Step A: register mean functional to reference mean
            # BOLDRigid = rigid body (6 DOF), appropriate for within-session
            # inter-run alignment where no anatomical deformation is expected.
            # -----------------------------------------------------------------
            IMG = ants.image_read(fMRI_runMean_reoriented)
            mTx = ants.registration(
                fixed=REF,
                moving=IMG,
                type_of_transform='BOLDRigid',
                initial_transform=None,
                outprefix=fMRI_run_inRef_mat)

            moved_mean = ants.apply_transforms(
                fixed=REF,
                moving=IMG,
                transformlist=mTx['fwdtransforms'],
                interpolator=n_for_ANTS)
            ants.image_write(moved_mean, fMRI_runMean_inRef, ri=False)
            _write_json(fMRI_runMean_inRef,
                        [fMRI_runMean_reoriented, runMean_unwarped_reoriented_Ref],
                        'Run-mean coregistration to reference run (ANTs BOLDRigid).',
                        f'transform: {fMRI_run_inRef_mat}')

            # -----------------------------------------------------------------
            # Step B: apply the same rigid transform to the full 4D timeseries.
            # imagetype=3 is required by ANTsPy for 4D time series — it treats
            # the last dimension as time and applies the spatial transform to
            # each volume independently.
            # -----------------------------------------------------------------
            FUNC = ants.image_read(fMRI_reoriented)
            moved_func = ants.apply_transforms(
                fixed=REF,
                moving=FUNC,
                transformlist=mTx['fwdtransforms'],
                interpolator=n_for_ANTS,
                imagetype=3)
            ants.image_write(moved_func, fMRI_run_inRef, ri=False)
            _write_json(fMRI_run_inRef,
                        [fMRI_reoriented, runMean_unwarped_reoriented_Ref],
                        'Full 4D timeseries coregistered to reference run (ANTs BOLDRigid).',
                        f'transform: {fMRI_run_inRef_mat}')

        # QC plot — always, for every run including reference
        plot_QC_func.plot_qc(
            runMean_unwarped_reoriented_Ref,
            fMRI_runMean_reoriented,
            opj(bids_dir, 'QC', 'fMRI_runMean_in_REF',
                root_RS + '_fMRI_runMean_in_REF.png'))
