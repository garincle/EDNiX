# import
import os
import subprocess
import glob
import json
import nibabel as nib
import numpy as np
import pandas as pd

opj = os.path.join
opb = os.path.basename
opr = os.path.relpath
opd = os.path.dirname
ope = os.path.exists
ops = os.path.splitext
opi = os.path.isfile

spgo = subprocess.getoutput

from Tools import run_cmd, getpath, diaryfile, Load_EDNiX_requirement, check_nii

from fonctions import chooseanat, _1_fMRI_preTTT_in_fMRIspace, _2_coregistration_to_norm, _3_mask_fMRI, _4_check_mask, \
    _5_anat_to_fMRI, _6_Melodic, _7_post_TTT, _8_fMRI_to_anat
from fonctions import _9_coregistration_to_template_space, _10_Correl_matrix, _11_Seed_base_many_regionsatlas, \
    _12_fMRI_QC, _14_fMRI_QC_matrix, _100_Data_Clean, _200_Data_QC
from fonctions.extract_filename import extract_filename


def preprocess_data(species, all_ID, all_Session, all_data_path, all_Session_max, stdy_template, stdy_template_mask,
                    BASE_SS, BASE_mask, T1_eq, Slice_timing_info, anat_func_same_space, use_master_for_Allineate,
                    correction_direction, REF_int, SBAspace, erod_seed, smoothSBA, deoblique, orientation,
                    TfMRI, GM_mask_studyT, GM_mask, creat_study_template, type_norm, coregistration_longitudinal,
                    dilate_mask, overwrite_option, nb_ICA_run, blur, ICA_cleaning, extract_exterior_CSF, extract_WM,
                    n_for_ANTS, aff_metric_ants, aff_metric_ants_Transl, list_atlases, selected_atlases, panda_files,
                    endfmri, endjson, endmap,
                    oversample_map, use_cortical_mask_func, cut_coordsX, cut_coordsY, cut_coordsZ, threshold_val,
                    Skip_step,
                    bids_dir, costAllin, use_erode_WM_func_masks, do_not_correct_signal, use_erode_V_func_masks,
                    folderforTemplate_Anat, IhaveanANAT, do_anat_to_func, Method_mask_func, segmentation_name_list,
                    band,
                    extract_Vc, selected_atlases_matrix, specific_roi_tresh, delta_thresh, extract_GS, MAIN_PATH,
                    DwellT, SED, TR, TRT, type_of_transform, ntimepoint_treshold, registration_fast, normalize,
                    reftemplate_path, reference, function_is_rest):
    sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _, sing_synstrip, Unetpath = Load_EDNiX_requirement.load_requirement(
        MAIN_PATH, reftemplate_path, bids_dir, 'yes')

    if species in ['Human', 'Chimpanzee']:
        config_f = opj(MAIN_PATH, 'Tool_library', 'config', 'b02b0Human.cnf')
    else:
        config_f = opj(MAIN_PATH, 'Tool_library', 'config', 'b02b0Macaque.cnf')

    # Set the environment variable for the current process
    os.environ["AFNI_NIFTI_TYPE_WARN"] = "NO"
    overwrite = ''
    if overwrite_option == True:
        overwrite = ' -overwrite'
    # Usage
    type_norm = check_nii.normalize_anat_type(type_norm)

    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, all_Session_max):
        nl = 'INFO: Work on ' + str(ID) + ' session ' + str(Session)
        run_cmd.printcolor(nl, 'HEADER')

        # Resting data organization folders ###########################################################################
        (path_func, dir_fmap, dir_prepro_orig, dir_prepro_orig_labels, dir_prepro_orig_masks,
         dir_prepro_orig_process, dir_prepro_orig_rs, dir_prepro_orig_task,
         dir_prepro_acpc, dir_prepro_acpc_labels, dir_prepro_acpc_masks,
         dir_prepro_acpc_process, dir_prepro_acpc_rs, dir_prepro_acpc_task,
         dir_prepro_template, dir_prepro_template_labels, dir_prepro_template_masks,
         dir_prepro_template_process, dir_prepro_template_rs, dir_prepro_template_task) = getpath.func(data_path,
                                                                                                       reference)

        for path in [path_func, dir_fmap, dir_prepro_fmap, dir_prepro_orig, dir_prepro_orig_labels, dir_prepro_orig_masks,
                     dir_prepro_orig_process, dir_prepro_orig_rs, dir_prepro_orig_task,
                     dir_prepro_acpc, dir_prepro_acpc_labels, dir_prepro_acpc_masks,
                     dir_prepro_acpc_process, dir_prepro_acpc_rs, dir_prepro_acpc_task,
                     dir_prepro_template, dir_prepro_template_labels, dir_prepro_template_masks,
                     dir_prepro_template_process, dir_prepro_template_rs, dir_prepro_template_task]:

            if not os.path.exists(path):
                if function_is_rest and path in (dir_prepro_orig_rs, dir_prepro_acpc_rs, dir_prepro_template_rs):
                    os.makedirs(path)
                elif not function_is_rest and path in (
                dir_prepro_orig_task, dir_prepro_acpc_task, dir_prepro_template_task):
                    os.makedirs(path)
                else:
                    os.makedirs(path)

        diary_file = diaryfile.create(opj(path_func, str(ID) + ' session ' + str(Session)), nl)
        diary_WARNING = opj(path_func, 'MAJOR_WARNING.txt')

        # link with the individual anatomical template ################################################################
        if IhaveanANAT == False:
            (anat_subject, brainmask, G_mask, V_mask, W_mask, dir_transfo, FS_dir,
             dir_prepro, volumes_dir, labels_dir, masks_dir) = chooseanat.create(folderforTemplate_Anat, diary_file)

        else:
            (anat_subject, brainmask, G_mask, V_mask, W_mask, dir_transfo, FS_dir,
             dir_prepro, volumes_dir, labels_dir, masks_dir) = chooseanat.retrieve(ID, data_path,
                                                                                   Session, anat_func_same_space,
                                                                                   use_erode_V_func_masks,
                                                                                   use_erode_WM_func_masks,
                                                                                   TfMRI, diary_file)

        ################# longitudinal co-registration  Yes or No ######################################################
        if coregistration_longitudinal == True:
            if creat_study_template == True:
                BASE_SS_coregistr = stdy_template
                BASE_SS_mask = stdy_template_mask
            else:
                BASE_SS_coregistr = BASE_SS
                BASE_SS_mask = BASE_mask

            if Session == max_ses:
                transfo_concat_Anat = [
                    opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat'),
                    opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_max_1InverseWarp.nii.gz')]
                w2inv_Anat = [True, False]
            else:
                data_path_max = opj(bids_dir, 'sub-' + ID, 'ses-' + str(max_ses))
                path_anat_max = opj(data_path_max, 'anat')
                dir_transfo_max = opj(path_anat_max, 'matrices')

                transfo_concat_Anat = [
                    opj(dir_transfo_max, 'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat'),
                    opj(dir_transfo_max, 'template_to_' + type_norm + '_SyN_final_max_1InverseWarp.nii.gz'),
                    opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
                    opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_1InverseWarp.nii.gz')]
                w2inv_Anat = [True, False, True, False]

        else:  # No
            if creat_study_template == True:
                BASE_SS_coregistr = stdy_template
                BASE_SS_mask = stdy_template_mask
            else:
                BASE_SS_coregistr = BASE_SS
                BASE_SS_mask = BASE_mask

            transfo_concat_Anat = [opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
                                   opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_1InverseWarp.nii.gz')]
            w2inv_Anat = [True, False]

        # Check the func runs
        list_RS = sorted(glob.glob(opj(path_func, endfmri)))
        RS = [opb(i) for i in list_RS]

        if len(list_RS) == 0:
            nl = ('ERROR : No func image found, we are look for an image define such as opj(dir_fMRI_Refth_RS, endfmri) and here it is ' +
                        str(opj(path_func, endfmri)) + ' I would check how you define "endfmri"')

            raise ValueError(run_cmd.error(nl, diary_file))

        nl = ("INFO: now let's check that this is a real a 4D fMRI image with enough time point as in define with the variable ntimepoint_treshold=" +
                    str(ntimepoint_treshold))
        run_cmd.msg(nl, diary_file, 'OKGREEN')

        list_RS_list = list_RS.copy()
        list_pop_index = []

        for imageF in list_RS_list:
            # Load the fMRI NIfTI image
            fmri_image = nib.load(imageF)
            # Get the shape of the image (x, y, z, t)
            image_shape = fmri_image.shape
            # Check the number of time points (4th dimension)
            if len(image_shape) == 4:
                ntimepoint = image_shape[3]  # The 4th dimension represents time
                nl = "INFO: " + f"Number of time points: {ntimepoint}"
                run_cmd.msg(nl, diary_file, 'OKGREEN')

                if int(ntimepoint) < ntimepoint_treshold:
                    index_of_imageF = list_RS.index(imageF)
                    nl = "INFO: We will not analyze " + str(imageF) + " because there are not enough time points"
                    run_cmd.msg(nl, diary_file, 'WARNING')

                    diary_WARNING_file = diaryfile.create(diary_WARNING, nl)

                    list_RS.pop(index_of_imageF)
                    list_pop_index.append(index_of_imageF)
            else:
                nl = "INFO: " + str(imageF) + " is not a 4D fMRI image"
                run_cmd.msg(nl, diary_file, 'WARNING')

                diary_WARNING_file = diaryfile.create(diary_WARNING, nl)

        nb_run = len(list_RS)
        # Setup for distortion correction using Fieldmaps

        #### find and correct the confound files
        for imageF in list_RS:
            # Load the fMRI NIfTI image
            fmri_image = nib.load(imageF)
            # Get the shape of the image (x, y, z, t)
            image_shape = fmri_image.shape
            # Check the number of time points (4th dimension)
            ntimepoint = image_shape[3]  # The 4th dimension represents time

            if ope(opj(path_func, extract_filename(imageF) + '_confounds.tsv')):
                confounds_df = pd.read_csv(opj(path_func, extract_filename(imageF) + '_confounds.tsv'), sep='\t')

                # Check if the number of rows in the confounds file matches the number of time points
                if len(confounds_df) != ntimepoint:
                    nl = ('Mismatch in the number of time points: ' + len(confounds_df) +
                          ' in confounds file vs ' + ntimepoint + ' in fMRI image.')

                    diary_WARNING_file = diaryfile.create(diary_WARNING)
                    run_cmd.msg(nl, diary_WARNING_file, 'WARNING')

                # Remove the first T1_eq rows from the confounds DataFrame
                if T1_eq > 0:
                    # Preserve the first column (assuming it's a label or name)
                    first_column = confounds_df.iloc[:, 0]
                    remaining_columns = confounds_df.iloc[:, 1:]

                    # Remove the first T1_eq rows from the remaining columns
                    remaining_columns = remaining_columns.iloc[T1_eq:]

                    # Concatenate the first column with the remaining columns
                    confounds_df = pd.concat([first_column, remaining_columns], axis=1)

                    # Save the modified confounds DataFrame back to the TSV file
                    confounds_df.to_csv(opj(dir_prepro_orig, extract_filename(imageF) + '_confounds_correct.tsv'), sep='\t', index=False)
                    print(f"Removed the first {T1_eq} TRs from the confounds file.")
                else:
                    print("T1_eq is not positive; no rows removed.")
            else:
                print(f"Confounds file not found: {opj(opd(imageF), extract_filename(imageF) + '_confounds.tsv')}")

        # Find the fmap images
        list_map = sorted(glob.glob(opj(dir_fmap, endmap)))
        nl = "looking for fmap image with the command glob.glob(" + str(opj(dir_fmap, endmap))
        run_cmd.msg(nl, diary_file, 'OKGREEN')
        if len(list_map) == 0:
            list_map = sorted(glob.glob(opj(path_func, endmap)))
            nl = "No fmap found in fmap folder, looking now in func folder with the command glob.glob(" + str(
                opj(path_func, endmap))
            run_cmd.msg(nl, diary_file, 'OKGREEN')
            if len(list_map) > 0:
                nl = "WARNING: We found some fieldmap images in the func folder, be sure that this is what you want"
                run_cmd.msg(nl, diary_file, 'WARNING')
            else:
                nl = "No fmap found in func folder either"
                run_cmd.msg(nl, diary_file, 'WARNING')
        else:
            nl = "We found " + str(list_map)
            run_cmd.msg(nl, diary_file, 'OKGREEN')

        if len(list_RS) > 0:
            RS_map = [opb(i) for i in list_map]

            ######### choose TOPUP strategy ##########
            if len(list_map) == 0:
                recordings = 'very_old'
                nl = 'WARNING: There is no image available for building a fieldmaps'
                run_cmd.msg(nl, diary_file, 'WARNING')

            elif len(list_map) == 1:
                comp = check_nii.comphd(opj(path_func, RS_map[0]), opj(path_func, RS[int(REF_int) - 1]))

                if comp == True:
                    recordings = 'old'  # there is only one AP recording to correct for field distortion
                else:
                    recordings = 'very_old'  # the only one AP recording to correct for field distortion is useless !!!
                    nl = 'WARNING : Before moving on, check the quality of the AP image, you may decide to NOT use it for correction'
                    run_cmd.msg(nl, diary_file, 'WARNING')

            elif len(list_map) == 2:
                recordings = '2_mapdir'  # there one AP per PA recordings in total
            elif len(list_map) > 2:
                recordings = 'new'  # there is one AP per PA recordings to correct for field distortion
                if len(list_pop_index) > 0:
                    list_pop_index.sort(reverse=True)
                    for index in list_pop_index:
                        list_map.pop(index)
                if not len(list_map) == len(list_RS):
                    nl = 'ERROR: Check the runs. There is probably one or two broken files that has been repeated and that you should remove !'
                    raise NameError(run_cmd.error(nl, diary_file))

            nl = 'INFO: recordings type detected: ' + str(recordings)
            run_cmd.msg(nl, diary_file, 'OKGREEN')

            list_json = sorted(glob.glob(opj(path_func, endjson)))

            # get useful information about the func images #############################################################
            ### check if we found some .json file
            if not list_json:
                nl = 'WARNING: no .json found!!, you will need to at least provide the TR. Beware that no TOPUP correction can be applied if you do not provide the DwellT as well.'
                run_cmd.msg(nl, diary_file, 'WARNING')

            else:
                f = open(list_json[0])
                info_RS = json.load(f)

            ## find metrics in header (for fun)
            nl = 'Get infos about the func image'
            run_cmd.msg(nl, diary_file, 'OKGREEN')

            ## TR
            if TR == 'Auto':
                try:
                    TR_val = info_RS["RepetitionTime"]
                    nl = 'TR = ' + str(TR_val)
                    run_cmd.msg(nl, diary_file, 'OKGREEN')

                except:
                    try:
                        # Calculate the time difference
                        slice_timing = info_RS["SliceTiming"]
                        nslice = int(len(slice_timing))
                        slice_timing.sort()
                        slice_intervals = [slice_timing[i + 1] - slice_timing[i] for i in range(len(slice_timing) - 1)]
                        # Calculate the TR
                        TR_val = ((sum(slice_intervals) / (nslice - 1)) * 1000) * nslice
                        nl = 'WARNING: TR not found in Header file!!!!! Repetition  Time (TR) calculated: ' + str(
                            TR_val) + ' seconds. YOU ABSOLUTELY NEED TO DOUBLE CHECK THAT!'
                        run_cmd.msg(nl, diary_file, 'WARNING')

                    except:
                        nl = (
                            "ERROR: TR was set to auto, but we were unable to find it inside the .json file, I know that's crazy but something might be wrong with it. "
                            "Either it was not available in this file or our automatic technic didn't work properly. Restart and provide the TR value as a string. It should solve this issue")
                        raise Exception(run_cmd.error(nl, diary_file))
            else:
                TR_val = TR
                nl = 'TR = ' + str(TR_val)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            if Slice_timing_info == 'Auto':
                try:
                    slice_timing = info_RS["SliceTiming"]
                    nl = "INFO: SliceTiming = " + str(slice_timing)
                    run_cmd.msg(nl, diary_file, 'OKGREEN')

                    STC = map(str, slice_timing)
                    stc = ' '.join(STC)
                    with open(opj(dir_prepro_orig, 'stc.txt'), 'w') as f:
                        f.write(stc)
                    f.close()
                except:
                    cmd = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + sing_afni + '3dinfo -slice_timing ' + \
                          list_RS[0]
                    nl = spgo(cmd).split('\n')
                    STC = nl[-1].split('|')
                    STC = list(map(float, STC))
                    if np.sum(STC) > 0:
                        nl = "INFO: SliceTiming = " + str(STC)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')

                    else:
                        nl = "WARNING: Slice Timing not found, this will be particularly DANGEROUS, you SHOULD PROVIDE MANUALLY ONE!"
                        run_cmd.msg(nl, diary_file, 'WARNING')

                        diary_WARNING_file = diaryfile.create(diary_WARNING, nl)

            elif isinstance(Slice_timing_info, list) == True:
                nl = "INFO: SliceTiming = " + str(Slice_timing_info)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

                STC = map(str, Slice_timing_info)
                stc = ' '.join(STC)
                with open(opj(dir_prepro_orig, 'stc.txt'), 'w') as f:
                    f.write(stc)
                f.close()

            try:
                TE = info_RS["EchoTime"]
                nl = "INFO: EchoTime = " + str(TE)
                run_cmd.msg(nl, diary_file, 'OKGREEN')
            except:
                nl = "INFO: EchoTime not found in header"
                run_cmd.msg(nl, diary_file, 'WARNING')

            try:
                EES = info_RS["EffectiveEchoSpacing"]
                nl = "INFO: Effective Echo Spacing = " + str(EES)
                run_cmd.msg(nl, diary_file, 'OKGREEN')
            except:
                nl = "INFO: Effective Echo Spacing not found in header"
                run_cmd.msg(nl, diary_file, 'WARNING')

            if TRT == 'Auto':
                try:
                    TRT_val = info_RS['TotalReadoutTime']
                    nl = "INFO: Total Readout Time = " + str(TRT_val)
                    run_cmd.msg(nl, diary_file, 'OKGREEN')
                except:
                    nl = "INFO: Total Readout Time not found in header"
                    run_cmd.msg(nl, diary_file, 'WARNING')

                    TRT_val = 'None'
            else:
                TRT_val = TRT
                nl = 'TRT = ' + str(TRT_val)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            # Find correction_direction
            if correction_direction == 'Auto':
                nl = 'INFO: input correction_direction was empty, let s try to find what is with the header'
                run_cmd.msg(nl, diary_file, 'OKGREEN')

                try:
                    PE_d2 = info_RS['PhaseEncodingDirection']
                except:
                    nl = "INFO: Phase Encoding Direction not found in header"
                    run_cmd.msg(nl, diary_file, 'WARNING')

                    recordings = 'very_old'
                    nl = 'WARNING : No distortion correction will be applied with fugue'
                    run_cmd.msg(nl, diary_file, 'WARNING')

                    correction_direction_val = 'None'
                    dmap = ''
                    dbold = ''
                    PE_d2 = 'None'
            else:
                nl = 'INFO: input correction_direction is the launcher was determined as' + str(correction_direction)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

                PE_d2 = 'None'

            if TRT_val != 'None':
                if PE_d2 == 'j' or correction_direction == 'y-':
                    dmap = '0 1 0 ' + str(TRT_val)
                    dbold = '0 -1 0 ' + str(TRT_val)
                    correction_direction_val = 'y-'
                elif PE_d2 == 'j-' or correction_direction == 'y':
                    dmap = '0 -1 0 ' + str(TRT_val)
                    dbold = '0 1 0 ' + str(TRT_val)
                    correction_direction_val = 'y'
                elif PE_d2 == 'i' or correction_direction == 'x-':
                    dmap = '1 0 0 ' + str(TRT_val)
                    dbold = '-1 0 0 ' + str(TRT_val)
                    correction_direction_val = 'x-'
                elif PE_d2 == 'i-' or correction_direction == 'x':
                    dmap = '-1 0 0 ' + str(TRT_val)
                    dbold = '1 0 0 ' + str(TRT_val)
                    correction_direction_val = 'x'
                else:
                    recordings = 'very_old'
                    nl = 'WARNING : No distortion correction will be applied with fugue'
                    run_cmd.msg(nl, diary_file, 'WARNING')
                    dmap = ''
                    dbold = ''
                    correction_direction_val = 'None'

            else:
                nl = 'WARNING: TRT not found'
                run_cmd.msg(nl, diary_file, 'WARNING')

                nl = 'WARNING : No distortion correction will be applied with fugue'
                run_cmd.msg(nl, diary_file, 'WARNING')

                recordings = 'very_old'
                dmap = ''
                dbold = ''
                correction_direction_val = ''

            # Find SliceEncodingDirection
            if SED == 'Auto':
                try:
                    SED_val = info_RS["SliceEncodingDirection"]
                except KeyError:
                    try:
                        if info_RS["ImageOrientationPatientDICOM"][0] == 1:
                            SED_val = "i"
                        elif info_RS["ImageOrientationPatientDICOM"][0] == -1:
                            SED_val = "i-"
                        elif info_RS["ImageOrientationPatientDICOM"][1] == 1:
                            SED_val = "j"
                        elif info_RS["ImageOrientationPatientDICOM"][1] == -1:
                            SED_val = "j-"
                        elif info_RS["ImageOrientationPatientDICOM"][2] == 1:
                            SED_val = "k"
                        elif info_RS["ImageOrientationPatientDICOM"][2] == -1:
                            SED_val = "k-"
                        else:
                            SED_val = 'None'
                            nl = 'WARNING !!!! Can not find the Slice Encoding Direction. No restriction of deformation can be applied (not a big deal)'
                            run_cmd.msg(nl, diary_file, 'WARNING')

                            nl = 'SED = ' + str(SED_val)
                            run_cmd.msg(nl, diary_file, 'WARNING')

                    except KeyError:
                        SED_val = 'None'
                        nl = 'WARNING !!!! Can not find the Slice Encoding Direction. No restriction of deformation can be applied (not a big deal)'
                        run_cmd.msg(nl, diary_file, 'WARNING')

                        nl = 'SED = ' + str(SED_val)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')
            else:
                SED_val = SED
                nl = 'SED = ' + str(SED_val)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            ### #Find Dwell time (to double check)
            if DwellT == 'Auto':
                try:
                    DwellT_val = "%.16f" % (float(info_RS["DwellTime"]))
                except:
                    try:
                        nslice = int(len(info_RS["SliceTiming"]))
                        DwellT_val = "%.16f" % (float(info_RS["TotalReadOutTimeEPI"] / nslice))
                    except:
                        DwellT_val = 'None'
                        nl = 'WARNING: could not find the Dwell time. Beware that no TOPUP correction can be applied! If you still want to do it provide a DwellT value manually as a string.'
                        run_cmd.msg(nl, diary_file, 'WARNING')

            else:
                DwellT_val = DwellT
                nl = 'DwellT = ' + str(DwellT_val)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            DIR = os.getcwd()
            nl = 'Working path : ' + DIR
            run_cmd.msg(nl, diary_file, 'OKGREEN')

            # Run the preprocessing ####################################################################################

            if 1 in Skip_step:
                nl = 'skip step ' + str(1)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            else:
                _1_fMRI_preTTT_in_fMRIspace.preprocess_data(dir_prepro_orig_process,
                                                            RS,
                                                            list_RS,
                                                            nb_run,
                                                            T1_eq,
                                                            TR_val,
                                                            Slice_timing_info,
                                                            overwrite,
                                                            sing_afni, diary_file, diary_WARNING)

            if 2 in Skip_step:
                nl = 'skip step ' + str(2)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            else:
                _2_coregistration_to_norm.coregist_to_norm(correction_direction_val,
                                                           dir_prepro_orig_process,
                                                           RS, RS_map, nb_run, recordings, REF_int, list_map, deoblique,
                                                           orientation, DwellT_val, n_for_ANTS,
                                                           overwrite, sing_afni, sing_fsl, dmap, dbold, config_f,
                                                           diary_file)

            if 3 in Skip_step:
                nl = 'skip step ' + str(3)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            else:
                _3_mask_fMRI.Refimg_to_meanfMRI(anat_func_same_space, BASE_SS_coregistr, TfMRI, dir_prepro_orig,
                                                dir_prepro_acpc, use_master_for_Allineate,
                                                dir_prepro_template, RS, nb_run, REF_int, ID, dir_prepro, brainmask,
                                                V_mask, W_mask, G_mask, dilate_mask,
                                                costAllin, anat_subject, Method_mask_func, overwrite, type_of_transform,
                                                aff_metric_ants,
                                                sing_afni, sing_fs, sing_fsl, sing_itk, diary_file)

            if 4 in Skip_step:
                nl = 'skip step ' + str(4)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            else:
                _4_check_mask._itk_check_masks(dir_prepro_orig, sing_itk, diary_file, sing_afni, overwrite)

            if 5 in Skip_step:
                nl = 'skip step ' + str(5)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            else:
                _5_anat_to_fMRI.Refimg_to_meanfMRI(REF_int, SED_val, anat_func_same_space,
                                                   TfMRI, dir_prepro_orig, dir_prepro_acpc, RS, nb_run,
                                                   ID, dir_prepro, n_for_ANTS, aff_metric_ants, aff_metric_ants_Transl,
                                                   list_atlases, labels_dir, anat_subject,
                                                   IhaveanANAT, do_anat_to_func, type_of_transform, registration_fast,
                                                   overwrite, sing_afni, diary_file)

            if 6 in Skip_step or ICA_cleaning == 'Skip':
                nl = 'skip step ' + str(6)
                run_cmd.msg(nl, diary_file, 'OKGREEN')


            else:
                _6_Melodic.Melodic_correct(dir_residuals_template, dir_residuals_acpc, dir_prepro_orig, dir_prepro_acpc,
                                           nb_ICA_run, nb_run, RS, ICA_cleaning, sing_fsl, sing_itk, TR_val, diary_file)

            if 7 in Skip_step:
                nl = 'skip step ' + str(7)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            else:
                _7_post_TTT.signal_regression(dir_prepro_orig, dir_residuals_acpc,
                                              nb_run, RS, blur, TR_val, ICA_cleaning, extract_exterior_CSF, extract_WM,
                                              normalize,
                                              do_not_correct_signal, band, extract_Vc, extract_GS, overwrite, sing_afni,
                                              diary_file)

            if 8 in Skip_step:
                nl = 'skip step ' + str(8)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            else:
                _8_fMRI_to_anat.to_anat_space(dir_prepro_orig, dir_prepro_acpc,
                                              nb_run, RS, n_for_ANTS, do_anat_to_func, anat_func_same_space, diary_file)

            if 9 in Skip_step:
                nl = 'skip step ' + str(9)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            else:
                _9_coregistration_to_template_space.to_common_template_space(deoblique, dir_prepro_orig,
                                                                             dir_prepro_acpc, dir_prepro_template,
                                                                             nb_run, RS, transfo_concat_Anat,
                                                                             w2inv_Anat, do_anat_to_func, list_atlases,
                                                                             BASE_SS_mask, GM_mask, GM_mask_studyT,
                                                                             creat_study_template,
                                                                             anat_func_same_space, orientation, REF_int,
                                                                             IhaveanANAT,
                                                                             overwrite, sing_afni, diary_file)

            if 10 in Skip_step:
                nl = 'skip step ' + str(10)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            else:
                _10_Correl_matrix.correl_matrix(dir_prepro_orig, RS, nb_run,
                                                selected_atlases_matrix, segmentation_name_list, ID, Session, TR_val,
                                                bids_dir, sing_afni, diary_file)

            if 11 in Skip_step:
                nl = 'skip step ' + str(11)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            else:
                _11_Seed_base_many_regionsatlas.SBA(SBAspace, BASE_SS_coregistr, erod_seed, dir_prepro_orig,
                                                    dir_prepro_acpc,
                                                    dir_prepro_template, RS, nb_run, selected_atlases, panda_files,
                                                    oversample_map,
                                                    use_cortical_mask_func, cut_coordsX, cut_coordsY, cut_coordsZ,
                                                    threshold_val,
                                                    sing_afni, diary_file, smoothSBA, TR_val)

            if 12 in Skip_step:
                nl = 'skip step ' + str(12)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            else:
                _12_fMRI_QC.fMRI_QC(correction_direction, dir_prepro_orig, dir_prepro_template, RS, nb_run, sing_afni,
                                    diary_file)

            if 14 in Skip_step:
                nl = 'skip step ' + str(14)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            else:
                _14_fMRI_QC_matrix.fMRI_QC_matrix(dir_prepro_orig, dir_prepro_acpc, dir_prepro_template,
                                                  specific_roi_tresh, delta_thresh, RS, nb_run, diary_file)

            if 100 in Skip_step:
                nl = 'skip step ' + str(100)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            else:
                _100_Data_Clean.clean(dir_prepro_orig, dir_prepro_acpc,
                                      dir_prepro_template, RS, nb_run, diary_file)

            if 200 in Skip_step:
                nl = 'skip step ' + str(200)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            else:
                _200_Data_QC._itk_check_coregistr(dir_prepro_template, BASE_SS_coregistr, sing_itk, diary_file)