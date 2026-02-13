#import
import os
import ants
import json

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile
from Tools import run_cmd
from fMRI.extract_filename import extract_filename
from fMRI import _2b_fix_orient, _2a_correct_img
from fMRI import plot_QC_func

def coregist_to_norm(correction_direction, list_RS, dir_prepro_fmap, dir_prepro_raw_process, dir_prepro_orig_process, RS, RS_map, nb_run, recordings, bids_dir,
                     REF_int, list_map, animalPosition, humanPosition, doWARPonfunc, dir_prepro_raw_matrices, orientation, DwellT, n_for_ANTS,
                     overwrite, sing_afni, sing_fsl, dmap, dbold, config_f, diary_file):
    nl = '##  Working on step ' + str(2) + '(function: _2_coregistration_to_norm).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    topup_f = open(opj(dir_prepro_raw_process, '4topup.txt'), "w")
    topup_f.write(dmap + ' \n')
    topup_f.write(dbold + ' \n')
    topup_f.close()
    topup_file = [opj(dir_prepro_raw_process, '4topup.txt'), config_f]

    if recordings != 'very_old':
        nl = 'INFO: DwellT is equal to ' + str(DwellT) + ' please check!!'
        run_cmd.msg(nl, diary_file, 'OKGREEN')

    # Determine what to process
    if recordings == 'new':
        items_to_process = zip(range(len(list_map)), range(nb_run))
    else:  # 'old', '2_mapdir', 'very_old'
        items_to_process = [(0, REF_int)]

    for i, r in items_to_process:
        # Get root names
        if recordings != 'very_old':
            root = extract_filename(RS_map[i])  # For fieldmap
            fMRI_runMean_fieldmap_rads = opj(dir_prepro_fmap, root + '_space-func_desc-runMean_fieldmap_rads.nii.gz')

        root_RS = extract_filename(RS[r])  # For functional data
        # Define file paths with correct naming
        fMRI_runMean_n4Bias = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-runMean_n4Bias.nii.gz')
        fMRI_runMean_unwarped = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-runMean_unwarped.nii.gz')
        runMean_unwarped_reoriented = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-runMean_unwarped_reoriented.nii.gz')

        if recordings == 'very_old':
            # Simple copy for very_old
            command = (sing_afni + '3dcopy ' + fMRI_runMean_n4Bias + ' ' + fMRI_runMean_unwarped + overwrite)
            run_cmd.run(command, diary_file)

            dictionary = {
                "Sources": fMRI_runMean_n4Bias,
                "Description": 'Copy.',
                "Command": command,}
        else:
            # Standard processing for other types
            # Correct image
            _2a_correct_img.correct_img(dir_prepro_orig_process, dir_prepro_fmap, fMRI_runMean_n4Bias, RS, list_map, RS_map, i, r, recordings,
				overwrite,sing_afni,sing_fsl,topup_file,diary_file)

            # Apply FUGUE correction
            command = (sing_fsl + 'fugue -i ' + fMRI_runMean_n4Bias +
                       ' --dwell=' + str(DwellT) + ' --loadfmap=' + fMRI_runMean_fieldmap_rads +
                       ' --unwarpdir=' + correction_direction + ' -u ' + fMRI_runMean_unwarped)
            run_cmd.run(command, diary_file)

            dictionary = {
                "Sources": [fMRI_runMean_n4Bias, fMRI_runMean_fieldmap_rads],
                "Description": 'Distortion correction (fugue from FSL).',
                "Command": command,}

        # Common steps for all types
        json_object = json.dumps(dictionary, indent=3)
        with open(fMRI_runMean_unwarped.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)

        if doWARPonfunc in ['WARP', 'header']:
            doWARPonfunc_spe = str(doWARPonfunc)
        else:
            doWARPonfunc_spe = 'WARP'
        # Fix orientation
        _2b_fix_orient.fix_orient(runMean_unwarped_reoriented, fMRI_runMean_unwarped, list_RS,
                                  animalPosition, humanPosition, orientation, doWARPonfunc_spe, sing_afni, diary_file)

    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ###
    ###                                              fix header problems                                                                ###
    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ###

    # 1.0 first you need to fix the header problems and potentially fit anat images
    for r in range(0, int(nb_run)):
        nl = 'work on ' + str(dir_prepro_raw_process) + ' run ' + str(i +1)
        run_cmd.msg(nl, diary_file, 'OKGREEN')
        root_RS = extract_filename(RS[r])

        fMRI_run_motion_corrected = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-motion_corrected.nii.gz')
        fMRI_reoriented = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_reoriented.nii.gz')
        fMRI_runMean_n4Bias = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-runMean_n4Bias.nii.gz')
        fMRI_runMean_reoriented = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_runMean_reoriented.nii.gz')

        ### 2.0 Start fix_orient
        _2b_fix_orient.fix_orient(fMRI_reoriented, fMRI_run_motion_corrected, list_RS,
                                  animalPosition, humanPosition, orientation, doWARPonfunc, sing_afni, diary_file)
        ### 2.0 Start fix_orient
        _2b_fix_orient.fix_orient(fMRI_runMean_reoriented, fMRI_runMean_n4Bias, list_RS,
                                  animalPosition, humanPosition, orientation, doWARPonfunc, sing_afni, diary_file)

    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ###
    ###                                          co-registration of each run to the norm                                                ###
    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ###

    for r in range(0, int(nb_run)):
        root_RS = extract_filename(RS[r])
        root_RS_ref = extract_filename(RS[REF_int])

        fMRI_reoriented = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_reoriented.nii.gz')
        #### run mean REF, in new orient, and unwarped
        runMean_unwarped_reoriented_Ref = opj(dir_prepro_raw_process, root_RS_ref + '_space-func_desc-runMean_unwarped_reoriented.nii.gz')
        REF = ants.image_read(runMean_unwarped_reoriented_Ref)

        #### run mean REF, in new orient (3dWARPed if needed)
        fMRI_runMean_reoriented = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_runMean_reoriented.nii.gz')

        #### result of the co-registration
        fMRI_runMean_inRef = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_runMean_inRef.nii.gz')
        fMRI_run_inRef_mat = opj(dir_prepro_raw_matrices, root_RS + '_space-func_desc-fMRI_run_inRef')
        fMRI_run_inRef = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_run_inRef.nii.gz')

        if root_RS == root_RS_ref and recordings == 'very_old' and doWARPonfunc in ['WARP', 'header']:  # do not process ref not corrected...

            command = (sing_afni + '3dcopy ' + fMRI_reoriented +
                       ' ' + fMRI_run_inRef + overwrite)
            run_cmd.run(command, diary_file)
            dictionary = {"Sources": fMRI_reoriented,
                          "Description": 'Copy.', "Command": command,}
            # Common steps for all types
            json_object = json.dumps(dictionary, indent=3)
            with open(fMRI_run_inRef.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)

            command = (sing_afni + '3dcopy ' + fMRI_runMean_reoriented +
                       ' ' + fMRI_runMean_inRef + overwrite)
            run_cmd.run(command, diary_file)
            dictionary = {"Sources": fMRI_runMean_reoriented,
                          "Description": 'Copy.', "Command": command,}
            # Common steps for all types
            json_object = json.dumps(dictionary, indent=3)
            with open(fMRI_runMean_inRef.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)

        else:
            # 1.0 calculate co-registration mean image ref to mean image func
            IMG = ants.image_read(fMRI_runMean_reoriented)
            mTx = ants.registration(fixed=REF, moving=IMG,
                                    type_of_transform='BOLDRigid',
                                    initial_transform=None,
                                    outprefix=fMRI_run_inRef_mat)

            ##  Apply the transformation to the mean image
            moved = ants.apply_transforms(fixed=REF, moving=IMG,
                                          transformlist=mTx['fwdtransforms'],
                                          interpolator=n_for_ANTS)
            ants.image_write(moved, fMRI_runMean_inRef, ri=False)

            dictionary = {"Sources": [fMRI_runMean_reoriented,
                                      runMean_unwarped_reoriented_Ref],
                          "Description": 'Coregistration (ANTspy).', }
            json_object = json.dumps(dictionary, indent=3)
            with open(fMRI_runMean_inRef.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)

            # 2.0 apply to all the volume in the func
            FUNC = ants.image_read(fMRI_reoriented)
            moved = ants.apply_transforms(fixed=REF, moving=FUNC,
                                          transformlist=mTx['fwdtransforms'],
                                          interpolator=n_for_ANTS, imagetype=3)
            ants.image_write(moved, fMRI_run_inRef, ri=False)

            dictionary = {"Sources": [runMean_unwarped_reoriented_Ref,
                                      fMRI_reoriented],
                          "Description": 'Coregistration (ANTspy).',}
            json_object = json.dumps(dictionary, indent=3)
            with open(fMRI_run_inRef.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)

            if not ope(opj(bids_dir, 'QC')):
                os.mkdir(opj(bids_dir, 'QC'))
            if not ope(opj(bids_dir, 'QC', 'fMRI_runMean_in_REF')):
                os.mkdir(opj(bids_dir, 'QC', 'fMRI_runMean_in_REF'))

        ####plot the QC
        plot_QC_func.plot_qc(runMean_unwarped_reoriented_Ref,
                             fMRI_runMean_reoriented,
                             opj(bids_dir, 'QC', 'fMRI_runMean_in_REF', root_RS + '_fMRI_runMean_in_REF.png'))

