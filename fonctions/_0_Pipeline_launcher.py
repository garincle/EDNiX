#import
import os
import subprocess
import glob
import json
import sys
import nibabel as nib
import datetime
import numpy as np
import pandas as pd
#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
ops = os.path.splitext
opi = os.path.isfile
spco = subprocess.check_output
spgo = subprocess.getoutput

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

import fonctions._1_fMRI_preTTT_in_fMRIspace
import fonctions._2_coregistration_to_norm
import fonctions._3_mask_fMRI
import fonctions._4_check_mask
import fonctions._5_anat_to_fMRI
import fonctions._6_Melodic
import fonctions._7_post_TTT
import fonctions._8_fMRI_to_anat
import fonctions._9_coregistration_to_template_space
import fonctions._10_Correl_matrix
import fonctions._11_Seed_base_many_regionsatlas
import fonctions._12_fMRI_QC
import fonctions._14_fMRI_QC_matrix
import fonctions._100_Data_Clean
import fonctions._200_Data_QC
from fonctions.extract_filename import extract_filename
import Tools.Load_EDNiX_requirement
def preprocess_data(all_ID, all_Session, all_data_path, all_Session_max, stdy_template, stdy_template_mask,
                    BASE_SS, BASE_mask, T1_eq, Slice_timing_info, anat_func_same_space,
                    correction_direction, REF_int, SBAspace, erod_seed, smoothSBA, deoblique, orientation,
                    TfMRI, GM_mask_studyT, GM_mask, creat_study_template, type_norm, coregistration_longitudinal,
                    dilate_mask, overwrite_option, nb_ICA_run, blur, ICA_cleaning, extract_exterior_CSF, extract_WM,
                    n_for_ANTS, aff_metric_ants, aff_metric_ants_Transl, list_atlases, selected_atlases, panda_files, endfmri, endjson, endmap,
                    oversample_map, use_cortical_mask_func, cut_coordsX, cut_coordsY, cut_coordsZ, threshold_val, Skip_step,
                    bids_dir, costAllin, use_erode_WM_func_masks, do_not_correct_signal, use_erode_V_func_masks,
                    folderforTemplate_Anat, IhaveanANAT, do_anat_to_func, Method_mask_func, segmentation_name_list, band,
                    extract_Vc, selected_atlases_matrix, specific_roi_tresh, unspecific_ROI_thresh, extract_GS, MAIN_PATH,
                    DwellT, SED, TR, TRT, type_of_transform, ntimepoint_treshold, registration_fast, FS_dir, normalize):

    s_path, afni_sif, fsl_sif, fs_sif, itk_sif, wb_sif, strip_sif, s_bind =  Tools.Load_EDNiX_requirement.load_requirement(MAIN_PATH, bids_dir, FS_dir)
    config_f = opj(MAIN_PATH, 'code', 'config', 'b02b0.cnf')

    if overwrite_option == True:
        overwrite = ' -overwrite'
    else:
        overwrite = ''


    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, all_Session_max):
        nl = 'INFO: Work on ' + str(ID) + ' session ' + str(Session)
        print(bcolors.OKGREEN + nl + bcolors.ENDC)

        # Resting data organization folders ###########################################################################

        dir_fMRI_Refth_RS          = opj(data_path, 'func')
        dir_fMRI_Refth_RS_prepro   = opj(dir_fMRI_Refth_RS, '01_prepro')
        dir_fMRI_Refth_RS_prepro1  = opj(dir_fMRI_Refth_RS_prepro, '01_funcspace')
        dir_fMRI_Refth_RS_prepro2  = opj(dir_fMRI_Refth_RS_prepro, '02_anatspace')
        dir_fMRI_Refth_RS_prepro3  = opj(dir_fMRI_Refth_RS_prepro, '03_atlas_space')
        dir_fMRI_Refth_RS_residual = opj(dir_fMRI_Refth_RS, '02_residual')
        dir_RS_ICA_native          = opj(dir_fMRI_Refth_RS_residual, '01_ICA_native')
        dir_RS_ICA_native_PreTT    = opj(dir_fMRI_Refth_RS_residual, '01_ICA_native_PreTT')
        dir_fMRI_Refth_map         = opj(data_path, 'fmap')

        if not ope(dir_fMRI_Refth_RS):
            os.makedirs(dir_fMRI_Refth_RS)
        if not ope(dir_fMRI_Refth_RS_prepro):
            os.makedirs(dir_fMRI_Refth_RS_prepro)
        if not ope(dir_fMRI_Refth_RS_prepro1):
            os.makedirs(dir_fMRI_Refth_RS_prepro1)

        date_file  = datetime.date.today()
        ct         = datetime.datetime.now()
        diary_name = str(ID) + ' session ' + str(Session) + str(date_file) + '.txt'
        diary_file = opj(dir_fMRI_Refth_RS, diary_name)
        if not opi(diary_file):
            diary = open(diary_file, "w")
            diary.write(f'\n{ct}')
            diary.write(f'\n{nl}')
        else:
            diary = open(diary_file, "a")
            diary.write(f'\n{ct}')
            diary.write(f'\n{nl}')

        # link with the individual anatomical template ################################################################

        if IhaveanANAT==False:

            # The anatomy organization folders
            dir_transfo  = ''
            dir_prepro   = ''
            volumes_dir  = ''
            labels_dir   = ''
            masks_dir    = ''

            anat_subject = opj(folderforTemplate_Anat,'template.nii.gz')
            brainmask    = opj(folderforTemplate_Anat,'brainmask.nii.gz')
            V_mask       = opj(folderforTemplate_Anat,'Vmask.nii.gz')
            W_mask       = opj(folderforTemplate_Anat,'Wmask.nii.gz')
            G_mask       = opj(folderforTemplate_Anat,'Gmask.nii.gz')

            nl = ' No individual anatomical template will be used for this animal'
            diary.write(f'\n{nl}')
        else:
            if anat_func_same_space == True:

                # The anatomy organization folders
                path_anat     = opj(data_path,'anat')
                dir_transfo   = opj(path_anat,'matrices')
                dir_native    = opj(path_anat,'native')
                dir_prepro    = opj(dir_native,'01_preprocess')
                wb_native_dir = opj(dir_native,'02_Wb')
                volumes_dir   = opj(wb_native_dir,'volumes')
                labels_dir    = opj(volumes_dir,'labels')
                masks_dir     = opj(volumes_dir,'masks')

            else:
                template_anat_for_fmri = glob.glob(opj(opd(data_path),'ses-' + str(Session),'anat','native','02_Wb', 'volumes', '*' + TfMRI + '_brain.nii*'))

                if len(template_anat_for_fmri) > 1:
                    nl = "WARNING: we found multiple anat template for this animal, please choose one!"
                    print(bcolors.WARNING + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    data_path_anat = input("Please enter manually a data_path_anat for preprocessing:")
                elif len(template_anat_for_fmri) == 0:
                    nl = "ERROR: We haven't found any anat template for this animal! We can't continue ! please provide a valid link for at least one anat image! current link is :" + str(template_anat_for_fmri)
                    diary.write(f'\n{nl}')
                    raise Exception(bcolors.FAIL + nl + bcolors.ENDC)
                else :
                    data_path_anat = opd(opd(opd(opd(opd(template_anat_for_fmri[0])))))
                    nl = 'INFO: We found this image as template: ' + str(template_anat_for_fmri)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                # The anatomy organization folders
                path_anat     = opj(data_path_anat,'anat')
                dir_transfo   = opj(path_anat,'matrices')
                dir_native    = opj(path_anat,'native')
                dir_prepro    = opj(dir_native,'01_preprocess')
                wb_native_dir = opj(dir_native,'02_Wb')
                volumes_dir   = opj(wb_native_dir,'volumes')
                labels_dir    = opj(volumes_dir,'labels')
                masks_dir     = opj(volumes_dir,'masks')


            anat_subject =  opj(volumes_dir,ID + type_norm + '_brain.nii.gz')
            brainmask    = opj(masks_dir,'brain_mask_in_anat_DC.nii.gz')
            G_mask       = opj(labels_dir, type_norm + 'Gmask.nii.gz')
            if use_erode_V_func_masks == True:
                V_mask = opj(masks_dir, type_norm + 'Vmask_erod.nii.gz')
            else:
                V_mask = opj(masks_dir, type_norm + 'Vmask.nii.gz')
            if use_erode_WM_func_masks == True:
                W_mask = opj(masks_dir, type_norm + 'Wmask_erod.nii.gz')
            else:
                W_mask = opj(masks_dir, type_norm + 'Wmask.nii.gz')



        ################# longitudinal co-registration  Yes or No ######################################################

        if coregistration_longitudinal == True:                   # Yes
            if creat_study_template == True:
                BASE_SS_coregistr = stdy_template
                BASE_SS_mask      = stdy_template_mask
            else:
                BASE_SS_coregistr = BASE_SS
                BASE_SS_mask      = BASE_mask

            if Session == max_ses:
                transfo_concat_Anat = [opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat'),
                                       opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_max_1InverseWarp.nii.gz')]
                w2inv_Anat = [True, False]
            else:
                data_path_max   = opj(bids_dir, 'sub-' + ID, 'ses-' + str(max_ses))
                path_anat_max   = opj(data_path_max, 'anat')
                dir_transfo_max = opj(path_anat_max, 'matrices')

                transfo_concat_Anat = [opj(dir_transfo_max, 'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat'),
                                       opj(dir_transfo_max, 'template_to_' + type_norm + '_SyN_final_max_1InverseWarp.nii.gz'),
                                       opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
                                       opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_1InverseWarp.nii.gz')]
                w2inv_Anat = [True, False, True, False]

        else:                                                    # No
            if creat_study_template == True:
                BASE_SS_coregistr = stdy_template
                BASE_SS_mask      = stdy_template_mask
            else:
                BASE_SS_coregistr = BASE_SS
                BASE_SS_mask      = BASE_mask

            transfo_concat_Anat = [opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
                                   opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_1InverseWarp.nii.gz')]
            w2inv_Anat = [True, False]


        # Check the func runs

        list_RS = sorted(glob.glob(opj(dir_fMRI_Refth_RS, endfmri)))
        RS      = [os.path.basename(i) for i in list_RS]

        if len(list_RS) == 0:
            nl = 'ERROR : No func image found, we are look for an image define such as opj(dir_fMRI_Refth_RS, endfmri) and here it is ' + str(opj(dir_fMRI_Refth_RS, endfmri)) + ' I would check how you define "endfmri"'
            diary.write(f'\n{nl}')
            raise ValueError(bcolors.FAIL + nl + bcolors.ENDC)

        nl = "INFO: now let's check that this is a real a 4D fMRI image with enough time point as in define with the variable ntimepoint_treshold=" + str(ntimepoint_treshold)
        print(bcolors.OKGREEN + nl + bcolors.ENDC)
        diary.write(f'\n{nl}')

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
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

                if int(ntimepoint) < ntimepoint_treshold:
                    index_of_imageF = list_RS.index(imageF)
                    nl = "INFO: We will not analyze " + str(imageF) + " because there are not enough time points"
                    print(bcolors.WARNING +nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    diary_file_WARNING = opj(opd(opd(dir_fMRI_Refth_RS_prepro1)), 'MAJOR_WARNING.txt')
                    if not opi(diary_file_WARNING):
                        diary_file_WARNING_file = open(diary_file_WARNING, "w")
                        diary_file_WARNING_file.write(f'\n{nl}')
                    else:
                        diary_file_WARNING_file = open(diary_file_WARNING, "a")
                        diary_file_WARNING_file.write(f'\n{nl}')
                    diary_file_WARNING_file.close()

                    list_RS.pop(index_of_imageF)
                    list_pop_index.append(index_of_imageF)
            else:
                nl = "INFO: " + str(imageF) + " is not a 4D fMRI image"
                print(bcolors.WARNING + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

                diary_file_WARNING = opj(opd(opd(dir_fMRI_Refth_RS_prepro1)), 'MAJOR_WARNING.txt')
                if not opi(diary_file_WARNING):
                    diary_file_WARNING_file = open(diary_file_WARNING, "w")
                    diary_file_WARNING_file.write(f'\n{nl}')
                else:
                    diary_file_WARNING_file = open(diary_file_WARNING, "a")
                    diary_file_WARNING_file.write(f'\n{nl}')
                diary_file_WARNING_file.close()

        nb_run  = len(list_RS)
        # Setup for distortion correction using Fieldmaps

        #### find and correct the counfound files
        for imageF in list_RS:
            # Load the fMRI NIfTI image
            fmri_image = nib.load(imageF)
            # Get the shape of the image (x, y, z, t)
            image_shape = fmri_image.shape
            # Check the number of time points (4th dimension)
            ntimepoint = image_shape[3]  # The 4th dimension represents time

            if ope(opj(opd(imageF), extract_filename(imageF) + '_confounds.tsv')):
                confounds_df = pd.read_csv(opj(opd(imageF), extract_filename(imageF) + '_confounds.tsv'), sep='\t')

                # Check if the number of rows in the confounds file matches the number of time points
                if len(confounds_df) != ntimepoint:
                    print(f"Mismatch in the number of time points: {len(confounds_df)} in confounds file vs {ntimepoint} in fMRI image.")

                    diary_file_WARNING = opj(opd(opd(dir_fMRI_Refth_RS_prepro1)), 'MAJOR_WARNING.txt')
                    if not opi(diary_file_WARNING):
                        diary_file_WARNING_file = open(diary_file_WARNING, "w")
                        diary_file_WARNING_file.write(f'\n{f"Mismatch in the number of time points: {len(confounds_df)} in confounds file vs {ntimepoint} in fMRI image."}')
                    else:
                        diary_file_WARNING_file = open(diary_file_WARNING, "a")
                        diary_file_WARNING_file.write(f'\n{f"Mismatch in the number of time points: {len(confounds_df)} in confounds file vs {ntimepoint} in fMRI image."}')
                    diary_file_WARNING_file.close()

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
                    confounds_df.to_csv(opj(dir_fMRI_Refth_RS_prepro1, extract_filename(imageF) + '_confounds_correct.tsv'), sep='\t', index=False)
                    print(f"Removed the first {T1_eq} TRs from the confounds file.")
                else:
                    print("T1_eq is not positive; no rows removed.")
            else:
                print(f"Confounds file not found: {opj(opd(imageF), extract_filename(imageF) + '_confounds.tsv')}")

        #### find and correct the fmap files
        list_map = sorted(glob.glob(opj(dir_fMRI_Refth_map, endmap)))
        nl = "looking for fmap image with the command glob.glob(" + str(opj(dir_fMRI_Refth_map, endmap))
        print(bcolors.OKGREEN + nl + bcolors.ENDC)
        diary.write(f'\n{nl}')
        nl = "We found " + str(list_map)
        print(bcolors.OKGREEN + nl + bcolors.ENDC)
        diary.write(f'\n{nl}')
        if len(list_RS)>0:
            RS_map   = [os.path.basename(i) for i in list_map]

            ######### choose TOPUP strategy ##########
            if len(list_map) == 0:
                recordings = 'very_old'
                nl = 'WARNING: There is no image available for building a fieldmaps'
                print(bcolors.WARNING + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
            elif len(list_map) == 1:
                cmd = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + \
                      '3dinfo -same_grid ' + opj(dir_fMRI_Refth_map, RS_map[0]) + ' ' + opj(dir_fMRI_Refth_RS, RS[int(REF_int)-1])
                dummy = spgo(cmd).split('\n')
                grid = int(dummy[-1])

                if int(grid) == 0:
                    recordings = 'very_old' # the only one AP recording to correct for field distortion is useless !!!
                    nl = 'WARNING : Before moving on, check the quality of the AP image, you may decide to NOT use it for correction'
                    print(bcolors.WARNING + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                elif int(grid[-3:][0]) == 1:
                    recordings = 'old'    # there is only one AP recording to correct for field distortion
            elif len(list_map) == 2:
                recordings = '2_mapdir'    # there one AP per PA recordings in total
            elif len(list_map) > 2:
                recordings = 'new' # there is one AP per PA recordings to correct for field distortion
                if len(list_pop_index) > 0:
                    list_pop_index.sort(reverse=True)
                    for index in list_pop_index:
                        list_map.pop(index)
                if not len(list_map) == len(list_RS):
                    nl = 'ERROR: Check the runs. There is probably one or two broken files that has been repeated and that you should remove !'
                    raise NameError(bcolors.FAIL + nl + bcolors.ENDC)

            nl = 'INFO: recordings type detected: ' + str(recordings)
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')
            list_json = sorted(glob.glob(opj(dir_fMRI_Refth_RS, endjson)))


            # get useful information about the func images #############################################################
            ### check if we found some .json file
            if not list_json:
                nl = 'WARNING: no .json found!!, you will need to at least provide the TR. Beware that no TOPUP correction can be applied if you do not provide the DwellT as well.'
                print(bcolors.WARNING + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
            else:
                f = open(list_json[0])
                info_RS = json.load(f)

            ## find metrics in header (for fun)
            nl = 'Get infos about the func image'
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

            ## TR
            if TR == 'Auto':
                try:
                    TR_val     = info_RS["RepetitionTime"]
                    nl = 'TR = ' + str(TR_val)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                except:
                    try:
                        # Calculate the time difference
                        slice_timing = info_RS["SliceTiming"]
                        nslice = int(len(slice_timing))
                        slice_timing.sort()
                        slice_intervals = [slice_timing[i + 1] - slice_timing[i] for i in range(len(slice_timing) - 1)]
                        # Calculate the TR
                        TR_val = ((sum(slice_intervals)/(nslice-1))*1000)*nslice
                        nl = 'WARNING: TR not found in Header file!!!!! Repetition  Time (TR) calculated: ' + str(TR_val) + ' seconds. YOU ABSOLUTELY NEED TO DOUBLE CHECK THAT!'
                        print(bcolors.WARNING + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')
                    except:
                        nl = ("ERROR: TR was set to auto, but we were unable to find it inside the .json file, I know that's crazy but something might be wrong with it. "
                              "Either it was not available in this file or our automatic technic didn't work properly. Restart and provide the TR value as a string. It should solve this issue")
                        diary.write(f'\n{nl}')
                        raise Exception(bcolors.FAIL + nl + bcolors.ENDC)
            else:
                TR_val = TR
                nl = 'TR = ' + str(TR_val)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

            if Slice_timing_info == 'Auto':
                try:
                    slice_timing = info_RS["SliceTiming"]
                    nl = "INFO: SliceTiming = " + str(slice_timing)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    STC = map(str,slice_timing)
                    stc = ' '.join(STC)
                    with open(opj(dir_fMRI_Refth_RS_prepro1,'stc.txt'),'w') as f:
                        f.write(stc)
                    f.close()
                except:
                    cmd = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -slice_timing ' + \
                          list_RS[0]
                    nl = spgo(cmd).split('\n')
                    STC = nl[-1].split('|')
                    STC = list(map(float, STC))
                    if np.sum(STC) > 0:
                        nl = "INFO: SliceTiming = " + str(STC)
                        print(bcolors.OKGREEN + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')
                    else :
                        nl = "WARNING: Slice Timing not found, this will be particularly DANGEROUS, you SHOULD PROVIDE MANUALLY ONE!"
                        print(bcolors.WARNING + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')

                        diary_file_WARNING = opj(opd(opd(dir_fMRI_Refth_RS_prepro1)), 'MAJOR_WARNING.txt')
                        if not opi(diary_file_WARNING):
                            diary_file_WARNING_file = open(diary_file_WARNING, "w")
                            diary_file_WARNING_file.write(f'\n{nl}')
                        else:
                            diary_file_WARNING_file = open(diary_file_WARNING, "a")
                            diary_file_WARNING_file.write(f'\n{nl}')
                        diary_file_WARNING_file.close()

            elif isinstance(Slice_timing_info, list) == True:
                nl = "INFO: SliceTiming = " + str(Slice_timing_info)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                STC = map(str, Slice_timing_info)
                stc = ' '.join(STC)
                with open(opj(dir_fMRI_Refth_RS_prepro1, 'stc.txt'), 'w') as f:
                    f.write(stc)
                f.close()

            try:
                TE = info_RS["EchoTime"]
                nl = "INFO: EchoTime = " + str(TE)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
            except:
                nl = "INFO: EchoTime not found in header"
                print(bcolors.WARNING + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
            try:
                EES    = info_RS["EffectiveEchoSpacing"]
                nl = "INFO: Effective Echo Spacing = " + str(EES)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
            except:
                nl = "INFO: Effective Echo Spacing not found in header"
                print(bcolors.WARNING + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

            if TRT == 'Auto':
                try:
                    TRT_val = info_RS['TotalReadoutTime']
                    nl = "INFO: Total Readout Time = " + str(TRT_val)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                except:
                    nl = "INFO: Total Readout Time not found in header"
                    print(bcolors.WARNING + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    TRT_val = 'None'
            else:
                TRT_val = TRT
                nl = 'TRT = ' + str(TRT_val)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

            # Find correction_direction
            if correction_direction == 'Auto':
                nl = 'INFO: input correction_direction was empty, let s try to find what is with the header'
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                try:
                    PE_d2 = info_RS['PhaseEncodingDirection']
                except:
                    nl = "INFO: Phase Encoding Direction not found in header"
                    print(bcolors.WARNING + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    recordings = 'very_old'
                    nl = 'WARNING : No distortion correction will be applied with fugue'
                    print(bcolors.WARNING + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    correction_direction_val = 'None'
                    dmap = ''
                    dbold = ''
                    PE_d2 = 'None'
            else :
                nl = 'INFO: input correction_direction is the launcher was determined as' + str(correction_direction)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
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
                    print(bcolors.WARNING + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    dmap = ''
                    dbold = ''
                    correction_direction_val = 'None'

            else:
                nl =  'WARNING: TRT not found'
                print(bcolors.WARNING + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                nl = 'WARNING : No distortion correction will be applied with fugue'
                print(bcolors.WARNING + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
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
                        if   info_RS["ImageOrientationPatientDICOM"][0] == 1:  SED_val = "i"
                        elif info_RS["ImageOrientationPatientDICOM"][0] == -1: SED_val = "i-"
                        elif info_RS["ImageOrientationPatientDICOM"][1] == 1:  SED_val = "j"
                        elif info_RS["ImageOrientationPatientDICOM"][1] == -1: SED_val = "j-"
                        elif info_RS["ImageOrientationPatientDICOM"][2] == 1:  SED_val = "k"
                        elif info_RS["ImageOrientationPatientDICOM"][2] == -1: SED_val = "k-"
                        else:
                            SED_val = 'None'
                            nl = 'WARNING !!!! Can not find the Slice Encoding Direction. No restriction of deformation can be applied (not a big deal)'
                            print(bcolors.WARNING + nl + bcolors.ENDC)
                            diary.write(f'\n{nl}')
                            nl = 'SED = ' + str(SED_val)
                            print(bcolors.OKGREEN + nl + bcolors.ENDC)
                            diary.write(f'\n{nl}')
                    except KeyError:
                        SED_val = 'None'
                        nl = 'WARNING !!!! Can not find the Slice Encoding Direction. No restriction of deformation can be applied (not a big deal)'
                        print(bcolors.WARNING + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')
                        nl = 'SED = ' + str(SED_val)
                        print(bcolors.OKGREEN + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')
            else:
                SED_val=SED
                nl = 'SED = ' + str(SED_val)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

            ### #Find Dwell time (to double check)
            if DwellT=='Auto':
                try:
                    DwellT_val    = "%.16f" % (float(info_RS["DwellTime"]))
                except:
                    try:
                        nslice = int(len(info_RS["SliceTiming"]))
                        DwellT_val   = "%.16f" % (float(info_RS["TotalReadOutTimeEPI"] / nslice))
                    except:
                        DwellT_val = 'None'
                        nl = 'WARNING: could not find the Dwell time. Beware that no TOPUP correction can be applied! If you still want to do it provide a DwellT value manually as a string.'
                        print(bcolors.WARNING + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')

            else:
                DwellT_val=DwellT
                nl = 'DwellT = ' + str(DwellT_val)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')


            DIR = os.getcwd()
            nl = 'Working path : ' + DIR
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

            diary.write(f'\n')
            diary.close()



            # Run the preprocessing ####################################################################################

            if 1 in Skip_step:
                ct = datetime.datetime.now()
                diary = open(diary_file, "a")
                diary.write(f'\n{ct}')
                nl = 'skip step ' + str(1)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                diary.write(f'\n')
                diary.close()

            else:
               fonctions._1_fMRI_preTTT_in_fMRIspace.preprocess_data(dir_fMRI_Refth_RS_prepro1,
                                                                      RS,
                                                                      list_RS,
                                                                      nb_run,
                                                                      T1_eq,
                                                                      TR_val,
                                                                      Slice_timing_info,
                                                                      overwrite,
                                                                      s_bind,afni_sif,diary_file)

            if 2 in Skip_step:
                ct = datetime.datetime.now()
                diary = open(diary_file, "a")
                diary.write(f'\n{ct}')
                nl = 'skip step ' + str(2)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

            else:
                fonctions._2_coregistration_to_norm.coregist_to_norm(correction_direction_val,
                                                                     dir_fMRI_Refth_RS_prepro1,
                                                                     RS, RS_map, nb_run, recordings, REF_int, list_map, deoblique,
                                                                     orientation, DwellT_val, n_for_ANTS,
                                                                     overwrite,s_bind,afni_sif,fsl_sif,dmap,dbold,config_f,diary_file)

            if 3 in Skip_step:
                ct = datetime.datetime.now()
                diary = open(diary_file, "a")
                diary.write(f'\n{ct}')
                nl = 'skip step ' + str(3)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                diary.write(f'\n')
                diary.close()

            else:
                fonctions._3_mask_fMRI.Refimg_to_meanfMRI(anat_func_same_space, BASE_SS_coregistr,TfMRI , dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                       dir_fMRI_Refth_RS_prepro3, RS, nb_run, REF_int, ID, dir_prepro, brainmask, V_mask, W_mask, G_mask, dilate_mask,
                       costAllin, anat_subject, Method_mask_func, overwrite, type_of_transform, aff_metric_ants,
                       s_bind,afni_sif,fs_sif, fsl_sif, itk_sif,diary_file)

            if 4 in Skip_step:
                ct = datetime.datetime.now()
                diary = open(diary_file, "a")
                diary.write(f'\n{ct}')
                nl = 'skip step ' + str(4)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                diary.write(f'\n')
                diary.close()

            else:
                fonctions._4_check_mask._itk_check_masks(dir_fMRI_Refth_RS_prepro1,
                                                         s_bind,itk_sif,diary_file, afni_sif, overwrite)

            if 5 in Skip_step:
                ct = datetime.datetime.now()
                diary = open(diary_file, "a")
                diary.write(f'\n{ct}')
                nl = 'skip step ' + str(5)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                diary.write(f'\n')
                diary.close()

            else:
                fonctions._5_anat_to_fMRI.Refimg_to_meanfMRI(REF_int, SED_val, anat_func_same_space,
                                                             TfMRI, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, RS, nb_run,
                                                             ID, dir_prepro, n_for_ANTS, aff_metric_ants, aff_metric_ants_Transl, list_atlases, labels_dir, anat_subject,
                                                             IhaveanANAT, do_anat_to_func, type_of_transform, registration_fast,
                                                             overwrite, s_bind, afni_sif,diary_file)

            if 6 in Skip_step or ICA_cleaning == 'Skip':
                nl = 'skip step ' + str(6)
                ct = datetime.datetime.now()
                diary = open(diary_file, "a")
                diary.write(f'\n{ct}')
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                diary.write(f'\n')
                diary.close()

            else:
                fonctions._6_Melodic.Melodic_correct(dir_RS_ICA_native_PreTT, dir_RS_ICA_native, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                    nb_ICA_run, nb_run, RS, ICA_cleaning, MAIN_PATH,s_bind,fsl_sif,itk_sif,TR_val, diary_file)

            if 7 in Skip_step:
                ct = datetime.datetime.now()
                diary = open(diary_file, "a")
                diary.write(f'\n{ct}')
                nl = 'skip step ' + str(7)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                diary.write(f'\n')
                diary.close()

            else:
                fonctions._7_post_TTT.signal_regression(dir_fMRI_Refth_RS_prepro1, dir_RS_ICA_native,
                    nb_run, RS, blur, TR_val, ICA_cleaning, extract_exterior_CSF, extract_WM, normalize,
                    do_not_correct_signal, band, extract_Vc, extract_GS, overwrite, s_bind,afni_sif,diary_file)

            if 8 in Skip_step:
                ct = datetime.datetime.now()
                diary = open(diary_file, "a")
                diary.write(f'\n{ct}')
                nl = 'skip step ' + str(8)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                diary.write(f'\n')
                diary.close()

            else:
                fonctions._8_fMRI_to_anat.to_anat_space(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                                                        nb_run, RS, n_for_ANTS, do_anat_to_func, anat_func_same_space,diary_file)

            if 9 in Skip_step:
                ct = datetime.datetime.now()
                diary = open(diary_file, "a")
                diary.write(f'\n{ct}')
                nl = 'skip step ' + str(9)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                diary.write(f'\n')
                diary.close()

            else:
                fonctions._9_coregistration_to_template_space.to_common_template_space(deoblique, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3,
                                                                                       nb_run, RS, transfo_concat_Anat,w2inv_Anat,do_anat_to_func, list_atlases,
                                                                                       BASE_SS_mask, GM_mask, GM_mask_studyT, creat_study_template,
                                                                                       anat_func_same_space, orientation, REF_int, IhaveanANAT,
                                                                                       overwrite,s_bind,afni_sif,diary_file)

            if 10 in Skip_step:
                ct = datetime.datetime.now()
                diary = open(diary_file, "a")
                diary.write(f'\n{ct}')
                nl = 'skip step ' + str(10)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                diary.write(f'\n')
                diary.close()

            else:
                fonctions._10_Correl_matrix.correl_matrix(dir_fMRI_Refth_RS_prepro1, RS, nb_run,
                                                          selected_atlases_matrix, segmentation_name_list, ID, Session, TR_val,
                                                          bids_dir,s_bind,afni_sif,diary_file)

            if 11 in Skip_step:
                ct = datetime.datetime.now()
                diary = open(diary_file, "a")
                diary.write(f'\n{ct}')
                nl = 'skip step ' + str(11)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                diary.write(f'\n')
                diary.close()

            else:
                fonctions._11_Seed_base_many_regionsatlas.SBA(SBAspace, BASE_SS_coregistr, erod_seed, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                                                              dir_fMRI_Refth_RS_prepro3, RS, nb_run, selected_atlases, panda_files, oversample_map,
                                                              use_cortical_mask_func,cut_coordsX, cut_coordsY, cut_coordsZ, threshold_val,
                                                              s_bind, afni_sif,diary_file, smoothSBA, TR_val)

            if 12 in Skip_step:
                ct = datetime.datetime.now()
                diary = open(diary_file, "a")
                diary.write(f'\n{ct}')
                nl = 'skip step ' + str(12)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                diary.write(f'\n')
                diary.close()

            else:
                fonctions._12_fMRI_QC.fMRI_QC(correction_direction, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro3, RS, nb_run, s_bind, afni_sif,diary_file)

            if 14 in Skip_step:
                ct = datetime.datetime.now()
                diary = open(diary_file, "a")
                diary.write(f'\n{ct}')
                nl = 'skip step ' + str(14)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                diary.write(f'\n')
                diary.close()

            else:
                fonctions._14_fMRI_QC_matrix.fMRI_QC_matrix(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3,
                   specific_roi_tresh, unspecific_ROI_thresh, RS, nb_run, diary_file)

            if 100 in Skip_step:
                ct = datetime.datetime.now()
                diary = open(diary_file, "a")
                diary.write(f'\n{ct}')
                nl = 'skip step ' + str(100)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                diary.write(f'\n')
                diary.close()

            else:
                fonctions._100_Data_Clean.clean(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                                                dir_fMRI_Refth_RS_prepro3, RS, nb_run,diary_file)

            if 200 in Skip_step:
                ct = datetime.datetime.now()
                diary = open(diary_file, "a")
                diary.write(f'\n{ct}')
                nl = 'skip step ' + str(200)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                diary.write(f'\n')
                diary.close()

            else:
                fonctions._200_Data_QC._itk_check_coregistr(dir_fMRI_Refth_RS_prepro3, BASE_SS_coregistr,
                                                            s_bind,itk_sif,diary_file)