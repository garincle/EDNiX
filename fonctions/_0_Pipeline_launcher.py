#import
import os
import subprocess
import glob
import json
import sys
import nibabel as nib

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
ops = os.path.splitext
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
import fonctions._13_fMRI_QC_SBA
import fonctions._14_fMRI_QC_matrix
import fonctions._100_Data_Clean
import fonctions._200_Data_QC

def preprocess_data(all_ID, all_Session, all_data_path, max_sessionlist, stdy_template, stdy_template_mask, BASE_SS, BASE_mask, T1_eq, anat_func_same_space, 
    correction_direction, REF_int, study_fMRI_Refth, SBAspace, erod_seed, deoblique, orientation,
    TfMRI, GM_mask_studyT, GM_mask, creat_study_template, type_norm, coregistration_longitudinal, dilate_mask, overwrite_option, nb_ICA_run, blur, melodic_prior_post_TTT,
    extract_exterior_CSF, extract_WM, n_for_ANTS, aff_metric_ants, list_atlases, selected_atlases, panda_files, endfmri, endjson, endmap, oversample_map, use_cortical_mask_func,
    cut_coordsX, cut_coordsY, cut_coordsZ, threshold_val, Skip_step, bids_dir, costAllin, use_erode_WM_func_masks, do_not_correct_signal, use_erode_V_func_masks,
    folderforTemplate_Anat, IhaveanANAT, doMaskingfMRI, do_anat_to_func, Method_mask_func, segmentation_name_list, band, extract_Vc, lower_cutoff, upper_cutoff, selected_atlases_matrix,
    specific_roi_tresh, unspecific_ROI_thresh, Seed_name, extract_GS, MAIN_PATH, DwellT, SED, TR, TRT, type_of_transform, ntimepoint_treshold, s_bind, s_path):
    sys.path.append(opj(MAIN_PATH + 'Code', 'EasyMRI_brain-master'))
    ### singularity set up
    afni_sif = ' ' + opj(s_path, 'afni_make_build_24_2_01.sif') + ' '
    fsl_sif = ' ' + opj(s_path, 'fsl_6.0.5.1-cuda9.1.sif') + ' '
    fs_sif = ' ' + opj(s_path, 'freesurfer_NHP.sif') + ' '
    itk_sif = ' ' + opj(s_path, 'itksnap_5.0.9.sif') + ' '
    config_f = opj(MAIN_PATH, 'code', 'config', 'b02b0.cnf')
    if overwrite_option == True:
        overwrite = ' -overwrite'
    else:
        overwrite = ''

    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, max_sessionlist):
        print(bcolors.OKGREEN + 'INFO: Work on ' + str(ID) + ' session ' + str(Session) + bcolors.ENDC)
        if IhaveanANAT==False:
            # The anatomy
            dir_transfo  = ''
            dir_prepro    = ''
            volumes_dir   = ''
            labels_dir    = ''
            masks_dir     = ''
        else:
            if anat_func_same_space == True:
                # The anatomy
                path_anat    = opj(data_path,'anat/')
                dir_transfo  = opj(path_anat,'matrices')
                dir_native    = opj(path_anat,'native')
                dir_prepro    = opj(dir_native,'01_preprocess')
                wb_native_dir = opj(dir_native,'02_Wb')
                volumes_dir   = opj(wb_native_dir,'volumes')
                labels_dir    = opj(volumes_dir,'labels')
                masks_dir     = opj(volumes_dir,'masks')
            else:
                template_anat_for_fmri = glob.glob(opj(opd(data_path),'ses-' + str(Session),'anat','native','02_Wb', 'volumes', '*' + TfMRI + '_brain.nii*'))
                print(bcolors.OKGREEN + 'INFO: We found this image as template: ' + str(template_anat_for_fmri) + bcolors.ENDC)
                if len(template_anat_for_fmri) == 1:
                    data_path_anat = opd(opd(opd(opd(opd(template_anat_for_fmri[0])))))
                elif len(template_anat_for_fmri) > 1:
                    print(bcolors.WARNING + "WARNING: we found multiple anat template for this animal, please choose one!" + bcolors.ENDC)
                    data_path_anat = input("Please enter manually a data_path_anat for preprocessing:")
                elif len(template_anat_for_fmri) == 0:
                    raise Exception(bcolors.FAIL + "ERROR: We havn't found any anat template for this animal! We can't continue ! please provid a valide link for at least one anat image! "
                                            "current link is :" + str(template_anat_for_fmri)+ bcolors.ENDC)
                # The anatomy
                path_anat    = opj(data_path_anat,'anat')
                dir_transfo  = opj(path_anat,'matrices')
                dir_native    = opj(path_anat,'native')
                dir_prepro    = opj(dir_native,'01_preprocess')
                wb_native_dir = opj(dir_native,'02_Wb')
                volumes_dir   = opj(wb_native_dir,'volumes')
                labels_dir    = opj(volumes_dir,'labels')
                masks_dir     = opj(volumes_dir,'masks')

        ################# coregistration longitudinal ???? #################
        if coregistration_longitudinal == True:
            if creat_study_template == True:
                BASE_SS_coregistr = stdy_template
                BASE_SS_mask = stdy_template_mask
            else:
                BASE_SS_coregistr = BASE_SS
                BASE_SS_mask = BASE_mask
            if Session == max_ses:
                transfo_concat_Anat = \
                    [opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_1Warp.nii.gz'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat')]
                w2inv_fwd = [False, False]
            else:
                data_path_max = opj(bids_dir, 'sub-' + ID, 'ses-' + str(max_ses))
                path_anat_max = opj(data_path_max, 'anat/')
                dir_transfo_max = opj(path_anat_max, 'matrices')

                transfo_concat_Anat = \
                    [opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_1Warp.nii.gz'),
                     opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
                     opj(dir_transfo_max, 'template_to_' + type_norm + '_SyN_final_max_1Warp.nii.gz'),
                     opj(dir_transfo_max, 'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat')]
                w2inv_Anat = [False, False, False, False]

        ################# coregistration non longitudinal #################
        else:
            if creat_study_template == True:
                BASE_SS_coregistr = stdy_template
                BASE_SS_mask = stdy_template_mask
            else:
                BASE_SS_coregistr = BASE_SS
                BASE_SS_mask = BASE_mask

            transfo_concat_Anat = [opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
                                   opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_1InverseWarp.nii.gz')]
            w2inv_Anat = [True, False]

        if IhaveanANAT == True:
            anat_subject =  opj(dir_prepro, ID + '_acpc_cropped' + TfMRI + '.nii.gz')
            brainmask     = opj(masks_dir,'brain_mask_in_anat_DC.nii.gz')
            if use_erode_V_func_masks == True:
                V_mask        = opj(masks_dir, type_norm + 'Vmask_erod.nii.gz')
            else:
                V_mask = opj(masks_dir, type_norm + 'Vmask.nii.gz')
            if use_erode_WM_func_masks == True:
                W_mask        = opj(masks_dir, type_norm + 'Wmask_erod.nii.gz')
            else:
                W_mask = opj(masks_dir, type_norm + 'Wmask.nii.gz')
            G_mask        = opj(labels_dir, type_norm + 'Gmask.nii.gz')
        else:
            anat_subject = opj(folderforTemplate_Anat,'template.nii.gz')
            brainmask     = opj(folderforTemplate_Anat,'brainmask.nii.gz')
            V_mask        = opj(folderforTemplate_Anat,'Vmask.nii.gz')
            W_mask = opj(folderforTemplate_Anat,'Wmask.nii.gz')
            G_mask = opj(folderforTemplate_Anat,'Gmask.nii.gz')

        # Resting data;
        dir_fMRI_Refth_RS          = opj(data_path,'func')
        dir_fMRI_Refth_RS_prepro   = opj(dir_fMRI_Refth_RS,'01_prepro')
        dir_fMRI_Refth_RS_prepro1  = opj(dir_fMRI_Refth_RS_prepro,'01_funcspace')
        dir_fMRI_Refth_RS_prepro2  = opj(dir_fMRI_Refth_RS_prepro,'02_anatspace')
        dir_fMRI_Refth_RS_prepro3  = opj(dir_fMRI_Refth_RS_prepro,'03_atlas_space')
        dir_fMRI_Refth_RS_residual = opj(dir_fMRI_Refth_RS,'02_residual')
        dir_RS_ICA_native    = opj(dir_fMRI_Refth_RS_residual,'01_ICA_native')
        dir_RS_ICA_native_PreTT    = opj(dir_fMRI_Refth_RS_residual,'01_ICA_native_PreTT')
        dir_fMRI_Refth_map          = opj(data_path,'fmap')

        # get useful informations
        list_RS = sorted(glob.glob(opj(dir_fMRI_Refth_RS, endfmri)))

        if len(list_RS) == 0:
            raise ValueError(bcolors.FAIL + 'ERROR : No func image found, we are look for an image define such as opj(dir_fMRI_Refth_RS, endfmri) '
                                            'and here it is ' + str(opj(dir_fMRI_Refth_RS, endfmri)) + ' I would check how you define "endfmri"' + bcolors.ENDC)
        print(bcolors.OKGREEN + "INFO: now let's check that this is a real a 4D fMRI image with enough time point as in define with the variable ntimepoint_treshold="
              + str(ntimepoint_treshold) + bcolors.ENDC)
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
                print(bcolors.OKGREEN + "INFO: " + f"Number of time points: {ntimepoint}" + bcolors.ENDC)
                if int(ntimepoint) < ntimepoint_treshold:
                    index_of_imageF = list_RS.index(imageF)
                    print(bcolors.OKGREEN + "INFO: We will not analyze " + str(imageF) + " because not enough time point" + bcolors.ENDC)
                    list_RS.pop(index_of_imageF)
                    list_pop_index.append(index_of_imageF)
            else:
                print(bcolors.OKGREEN + "INFO: This is not a 4D fMRI image" + bcolors.ENDC)

        nb_run = len(list_RS)
        RS = [os.path.basename(i) for i in list_RS]
        print(bcolors.OKGREEN + "INFO: We will analyse run " + str(list_RS) + bcolors.ENDC)

        list_map = sorted(glob.glob(opj(dir_fMRI_Refth_map, endmap))) #[]
        print(bcolors.OKGREEN + "looking for fmap image with the command glob.glob(" + str(opj(dir_fMRI_Refth_map, endmap)) + bcolors.ENDC)
        print(bcolors.OKGREEN + "We found " + str(list_map) + bcolors.ENDC)
        if len(list_RS)>0:
            RS_map   = [os.path.basename(i) for i in list_map]
            ############################################
            #### choose TOPUP stragtegy ################
            ############################################
            if len(list_map) == 0:
                recordings = 'very_old'
                print(bcolors.WARNING + 'WARNING: Before moving on, chech the quality of the AP image, you may decide to NOT use it for correction' + bcolors.ENDC)
            elif len(list_map) == 1:
                cmd = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + \
                      '3dinfo -same_grid ' + opj(dir_fMRI_Refth_map, RS_map[0]) + ' ' + opj(dir_fMRI_Refth_RS, RS[int(REF_int)-1])
                dummy = spgo(cmd).split('\n')
                grid = int(dummy[-1])
                if int(grid) == 0:
                    recordings = 'very_old' # the only one AP recording to correct for field distorsion is useless !!!
                    print(bcolors.WARNING + 'WARNING : Before moving on, chech the quality of the AP image, you may decide to NOT use it for correction' + bcolors.ENDC)

                elif int(grid[-3:][0]) == 1:
                    recordings = 'old'    # there is only one AP recording to correct for field distorsion
            elif len(list_map) == 2:
                recordings = '2_mapdir'    # there one AP per PA recordings in total
            elif len(list_map) > 2:
                recordings = 'new' # there is one AP per PA recordings to correct for field distorsion
                if len(list_pop_index) > 0:
                    list_pop_index.sort(reverse=True)
                    for index in list_pop_index:
                        list_map.pop(index)
                if not len(list_map) == len(list_RS):
                   raise NameError(bcolors.FAIL + 'ERROR: Check the runs there is probably one broken file that has been repeated and that you should remove !' + bcolors.ENDC)

            print(bcolors.OKGREEN + 'INFO: recordings type detected: ' + str(recordings) + bcolors.ENDC)
            list_json = sorted(glob.glob(opj(dir_fMRI_Refth_RS, endjson)))

            ### check if we found some .json file
            if not list_json:
                print(bcolors.WARNING + 'WARNING: no .json found!!, you will need to at least provide the TR, and not TOPUP correction will be applied if you do not provide DwellT' + bcolors.ENDC)
            else:
                f = open(list_json[0])
                info_RS = json.load(f)
            if TR == 'Auto':
                ## TR
                try:
                    TR_val     = info_RS["RepetitionTime"]
                    print(bcolors.OKGREEN + 'TR = ' + str(TR_val) + bcolors.ENDC)
                except:
                    try:
                        # Calculate the time difference
                        slice_timing = info_RS["SliceTiming"]
                        nslice = int(len(info_RS["SliceTiming"]))
                        slice_timing.sort()
                        slice_intervals = [slice_timing[i + 1] - slice_timing[i] for i in range(len(slice_timing) - 1)]
                        # Calculate the TR
                        TR_val = ((sum(slice_intervals)/(nslice-1))*1000)*nslice
                        print(bcolors.WARNING + 'WARNING: TR not found in Header file!!!!! Repetition  Time (TR) calculated: ' + str(TR_val) + ' seconds.'
                                                'YOU ABSOLUTELY NEED TO DOUBLE CHECK THAT!' + bcolors.ENDC)
                    except:
                        raise Exception(bcolors.FAIL + "ERROR: TR was set to auto, but we were unable to find it with the .json file, "
                                                       "I know that's crazy but something might be wrong with it. Either it was not "
                                                       "avaialble in this file or our auto technic didn't work. Restart and provide the "
                                                       "TR value as str should solve this issue" + bcolors.ENDC)
            else:
                TR_val = TR
                print(bcolors.OKGREEN + 'TR = ' + str(TR_val) + bcolors.ENDC)

            ## find metrics in header (for fun)
            print('find metrics in header (for fun)')
            try:
                slice_timing = info_RS["SliceTiming"]
                print(bcolors.OKGREEN + "INFO: SliceTiming = " + str(slice_timing) + bcolors.ENDC)
            except:
                print(bcolors.OKGREEN + "INFO: Slice Timing not found" + bcolors.ENDC)
            try:
                TE     = info_RS["EchoTime"]
                print(bcolors.OKGREEN + "INFO: EchoTime = " + str(TE) + bcolors.ENDC)
            except:
                print(bcolors.OKGREEN + "INFO: EchoTime not found in header" + bcolors.ENDC)
            try:
                EES    = info_RS["EffectiveEchoSpacing"]
                print(bcolors.OKGREEN + "INFO: Effective Echo Spacing = " + str(EES) + bcolors.ENDC)
            except:
                print(bcolors.OKGREEN + "INFO: Effective Echo Spacing not found in header" + bcolors.ENDC)

            if TRT == 'Auto':
                try:
                    TRT_val = info_RS['TotalReadoutTime']
                    print(bcolors.OKGREEN + "INFO: Total Readout Time = " + str(TRT_val) + bcolors.ENDC)
                except:
                    print(bcolors.OKGREEN + "INFO: Total Readout Time not found in header" + bcolors.ENDC)
                    TRT_val = 'None'
            else:
                TRT_val = TRT
                print(bcolors.OKGREEN + 'TRT = ' + str(TRT_val) + bcolors.ENDC)

            #Find correction_direction
            if correction_direction == 'Auto':
                print(bcolors.OKGREEN + 'INFO: input correction_direction was empty, let s try to find what is with the header' + bcolors.ENDC)
                try:
                    PE_d2 = info_RS['PhaseEncodingDirection']
                except:
                    print(bcolors.OKGREEN + "INFO: Phase Encoding Direction not found in header" + bcolors.ENDC)
                    dmap=''
                    dbold=''
                    print(bcolors.WARNING + "WARNING :Phase Encoding Direction not found in header and you didn't provided any" + bcolors.ENDC)
                    recordings = 'very_old'
                    print(bcolors.WARNING + 'WARNING : recordings = very_old no distortion correction will be applied with fugue' + bcolors.ENDC)
                    correction_direction_val = 'None'
                    PE_d2 = 'None'
            else :
                (bcolors.OKGREEN + 'INFO: input correction_direction is the launcher was determined as' + str(correction_direction) + bcolors.ENDC)
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
                    print(bcolors.WARNING + 'WARNING : recordings = very_old no distortion correction will be applied with fugue' + bcolors.ENDC)
                    dmap = ''
                    dbold = ''
            else:
                print(bcolors.WARNING + 'WARNING: TRT not found' + bcolors.ENDC)
                recordings = 'very_old'
                print(bcolors.WARNING + 'WARNING : recordings = very_old no distortion correction will be applied with fugue' + bcolors.ENDC)
                dmap = ''
                dbold = ''
                correction_direction_val = ''

            #Find SliceEncodingDirection
            if SED == 'Auto':
                try:
                    SED_val = info_RS["SliceEncodingDirection"]
                except KeyError:
                    try:
                        if info_RS["ImageOrientationPatientDICOM"][0]   == 1:  SED_val = "i"
                        elif info_RS["ImageOrientationPatientDICOM"][0] == -1: SED_val = "i-"
                        elif info_RS["ImageOrientationPatientDICOM"][1] == 1:  SED_val = "j"
                        elif info_RS["ImageOrientationPatientDICOM"][1] == -1: SED_val = "j-"
                        elif info_RS["ImageOrientationPatientDICOM"][2] == 1:  SED_val = "k"
                        elif info_RS["ImageOrientationPatientDICOM"][2] == -1: SED_val = "k-"
                        else:
                            SED_val = 'None'
                            print(bcolors.WARNING + 'WARNING !!!! Can not find SliceEncodingDirection in the DICOM, no restriction of '
                                                    'deformation will be applied (not a big deal)' + bcolors.ENDC)
                            print(bcolors.OKGREEN + 'SED = ' + str(SED_val))
                    except KeyError:
                        SED_val = 'None'
                        print(bcolors.WARNING + 'WARNING !!!! Can not find SliceEncodingDirection in the DICOM, no restriction of '
                                                'deformation will be applied (not a big deal)' + bcolors.ENDC)
                        print(bcolors.OKGREEN + 'SED = ' + str(SED_val))
            else:
                SED_val=SED
                print(bcolors.OKGREEN + 'SED = ' + str(SED_val) + bcolors.ENDC)

            ####Find Dwell time (to double check)
            if DwellT=='Auto':
                try:
                    DwellT_val    = "%.16f" % (float(info_RS["DwellTime"]))
                except:
                    try:
                        DwellT_val   = "%.16f" % (float(info_RS["TotalReadOutTimeEPI"] / nslice))
                    except:
                        DwellT_val = 'None'
                        print(bcolors.WARNING + 'WARNING: couldn t find Dwell time, will not apply TOPUP correction !!! if you want to do it, '
                                                'then provide a DwellT value manually as str!!!!' + bcolors.ENDC)
            else:
                DwellT_val=DwellT
                print(bcolors.OKGREEN + 'DwellT = ' + str(DwellT_val) + bcolors.ENDC)
            DIR = os.getcwd()
            print(bcolors.OKGREEN + 'Working path : ' + DIR + bcolors.ENDC)

            if 1 in Skip_step:
                print(bcolors.OKGREEN + 'skip step ' + str(1) + bcolors.ENDC)
            else:
                print(bcolors.OKGREEN + '##########   Working on step ' + str(1) + ' _1_fMRI_preTTT_in_fMRIspace  ###############' + bcolors.ENDC)
                print(bcolors.OKGREEN + str(ID) + ' Session ' + str(Session) + bcolors.ENDC)
                fonctions._1_fMRI_preTTT_in_fMRIspace.preprocess_data(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, RS, list_RS, nb_run, T1_eq, overwrite,
                                                                      s_bind,afni_sif)
            if 2 in Skip_step:
                print('skip step ' + str(2))
            else:
                print(bcolors.OKGREEN + '##########   Working on step ' + str(2) + ' _2_coregistration_to_norm  ###############' + bcolors.ENDC)
                print(bcolors.OKGREEN + str(ID) + ' Session ' + str(Session) + bcolors.ENDC)
                fonctions._2_coregistration_to_norm.coregist_to_norm(correction_direction_val, dir_fMRI_Refth_RS_prepro1, RS, RS_map, nb_run, recordings, REF_int, list_map, study_fMRI_Refth, deoblique,
                                                                     orientation, DwellT_val, n_for_ANTS, overwrite,s_bind,afni_sif,fsl_sif,dmap,dbold,config_f)

            if 3 in Skip_step:
                print(bcolors.OKGREEN + 'skip step ' + str(3) + bcolors.ENDC)
            else:
                print(bcolors.OKGREEN + '##########   Working on step ' + str(3) + ' _3_mask_fMRI  ###############' + bcolors.ENDC)
                print(bcolors.OKGREEN + str(ID) + ' Session ' + str(Session))
                fonctions._3_mask_fMRI.Refimg_to_meanfMRI(anat_func_same_space, BASE_SS_coregistr, TfMRI, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                       dir_fMRI_Refth_RS_prepro3, RS, nb_run, REF_int, ID, dir_prepro, brainmask, V_mask, W_mask, G_mask, dilate_mask,
                       costAllin, anat_subject, doMaskingfMRI, Method_mask_func, lower_cutoff, upper_cutoff, overwrite, type_of_transform, aff_metric_ants,
                       s_bind,afni_sif,fs_sif,  fsl_sif, itk_sif)
            if 4 in Skip_step:
                print(bcolors.OKGREEN + 'skip step ' + str(4) + bcolors.ENDC)
            else:
                print(bcolors.OKGREEN + '##########   Working on step ' + str(4) + ' _itk_check_masks  ###############' + bcolors.ENDC)
                print(bcolors.OKGREEN + str(ID) + ' Session ' + str(Session) + bcolors.ENDC)
                fonctions._4_check_mask._itk_check_masks(dir_fMRI_Refth_RS_prepro1,s_bind,itk_sif)

            if 5 in Skip_step:
                print(bcolors.OKGREEN + 'skip step ' + str(5) + bcolors.ENDC)
            else:
                print(bcolors.OKGREEN + '##########   Working on step ' + str(5) + ' _5_anat_to_fMRI  ###############' + bcolors.ENDC)
                print(bcolors.OKGREEN + str(ID) + ' Session ' + str(Session) + bcolors.ENDC)
                fonctions._5_anat_to_fMRI.Refimg_to_meanfMRI(REF_int, SED_val, anat_func_same_space, TfMRI, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, RS, nb_run, ID, dir_prepro, n_for_ANTS, aff_metric_ants, list_atlases, labels_dir, anat_subject, IhaveanANAT, do_anat_to_func, type_of_transform, overwrite, s_bind, afni_sif)

            if 6 in Skip_step:
                print(bcolors.OKGREEN + 'skip step ' + str(6) + bcolors.ENDC)
            else:
                if melodic_prior_post_TTT == True:
                    print(bcolors.OKGREEN + '##########   Working on step ' + str(6) + ' Melodic_correct  ###############' + bcolors.ENDC)
                    print(bcolors.OKGREEN + str(ID) + ' Session ' + str(Session) + bcolors.ENDC)
                    fonctions._6_Melodic.Melodic_correct(dir_RS_ICA_native_PreTT, dir_RS_ICA_native, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                nb_ICA_run, nb_run, RS, TfMRI, overwrite,s_bind,fsl_sif,itk_sif,TR_val)

            if 7 in Skip_step:
                print(bcolors.OKGREEN + 'skip step ' + str(7) + bcolors.ENDC)
            else:
                print(bcolors.OKGREEN + '##########   Working on step ' + str(7) + ' _7_post_TTT  ###############' + bcolors.ENDC)
                print(bcolors.OKGREEN + str(ID) + ' Session ' + str(Session))
                fonctions._7_post_TTT.signal_regression(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_RS_ICA_native,
                nb_run, RS, blur, TR_val, melodic_prior_post_TTT, extract_exterior_CSF, extract_WM, do_not_correct_signal, band, extract_Vc, extract_GS, overwrite,s_bind,afni_sif)

            if 8 in Skip_step:
                print(bcolors.OKGREEN + 'skip step ' + str(8) + bcolors.ENDC)
            else:
                print(bcolors.OKGREEN + '##########   Working on step ' + str(8) + ' _8_fMRI_to_anat  ###############' + bcolors.ENDC)
                print(bcolors.OKGREEN + str(ID) + ' Session ' + str(Session))
                fonctions._8_fMRI_to_anat.to_anat_space(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                nb_run, RS, n_for_ANTS, do_anat_to_func, anat_func_same_space)

            if 9 in Skip_step:
                print(bcolors.OKGREEN + 'skip step ' + str(9) + bcolors.ENDC)
            else:
                print(bcolors.OKGREEN + '##########   Working on step ' + str(9) + ' _9_coregistration_to_template_space  ###############' + bcolors.ENDC)
                print(bcolors.OKGREEN + str(ID) + ' Session ' + str(Session))
                fonctions._9_coregistration_to_template_space.to_common_template_space(deoblique, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3,
                nb_run, RS, transfo_concat_Anat,w2inv_Anat,do_anat_to_func, n_for_ANTS, list_atlases, TfMRI, BASE_SS_mask, GM_mask, GM_mask_studyT, creat_study_template,
                anat_func_same_space, orientation, dir_prepro, ID, REF_int, IhaveanANAT, overwrite,s_bind,afni_sif)

            if 10 in Skip_step:
                print(bcolors.OKGREEN + 'skip step ' + str(10) + bcolors.ENDC)
            else:
                print(bcolors.OKGREEN + '##########   Working on step ' + str(10) + ' _10_Correl_matrix  ###############' + bcolors.ENDC)
                print(bcolors.OKGREEN + str(ID) + ' Session ' + str(Session))
                fonctions._10_Correl_matrix.correl_matrix(dir_fMRI_Refth_RS_prepro1, RS, nb_run, selected_atlases_matrix, segmentation_name_list, ID, Session, bids_dir,s_bind,afni_sif)

            if 11 in Skip_step:
                print(bcolors.OKGREEN + 'skip step ' + str(11) + bcolors.ENDC)
            else:
                print(bcolors.OKGREEN + '##########   Working on step ' + str(11) + ' _11_Seed_base_many_regionsatlas  ###############' + bcolors.ENDC)
                print(bcolors.OKGREEN + str(ID) + ' Session ' + str(Session) + bcolors.ENDC)
                fonctions._11_Seed_base_many_regionsatlas.SBA(SBAspace, BASE_SS_coregistr, erod_seed, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                dir_fMRI_Refth_RS_prepro3, RS, nb_run, selected_atlases, panda_files, oversample_map, use_cortical_mask_func,
                cut_coordsX, cut_coordsY, cut_coordsZ, threshold_val, s_bind, afni_sif)

            if 12 in Skip_step:
                print(bcolors.OKGREEN + 'skip step ' + str(12) + bcolors.ENDC)
            else:
                print(bcolors.OKGREEN + '##########   Working on step ' + str(12) + ' _12_fMRI_QC  ###############' + bcolors.ENDC)
                print(bcolors.OKGREEN + str(ID) + ' Session ' + str(Session) + bcolors.ENDC)
                fonctions._12_fMRI_QC.fMRI_QC(correction_direction, ID, Session, segmentation_name_list, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3, specific_roi_tresh, unspecific_ROI_thresh, RS, nb_run, bids_dir,s_bind,afni_sif)

            if 13 in Skip_step:
                print(bcolors.OKGREEN + 'skip step ' + str(13) + bcolors.ENDC)
            else:
                print(bcolors.OKGREEN + '##########   Working on step ' + str(12) + ' _13_fMRI_QC_SBA  ###############' + bcolors.ENDC)
                print(bcolors.OKGREEN + str(ID) + ' Session ' + str(Session) + bcolors.ENDC)
                fonctions._13_fMRI_QC_SBA.fMRI_QC_SBA(Seed_name, BASE_SS_coregistr, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                            dir_fMRI_Refth_RS_prepro3, RS, nb_run, selected_atlases, panda_files, oversample_map,
                            use_cortical_mask_func)

            if 14 in Skip_step:
                print(bcolors.OKGREEN + 'skip step ' + str(14) + bcolors.ENDC)
            else:
                print(bcolors.OKGREEN + '##########   Working on step ' + str(12) + ' _14_fMRI_QC_matrix  ###############' + bcolors.ENDC)
                print(bcolors.OKGREEN + str(ID) + ' Session ' + str(Session) + bcolors.ENDC)
                fonctions._14_fMRI_QC_matrix.fMRI_QC_matrix(ID, Session, segmentation_name_list, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3, specific_roi_tresh, unspecific_ROI_thresh, RS, nb_run, bids_dir,s_bind,afni_sif)

            if 100 in Skip_step:
                print(bcolors.OKGREEN + 'skip step ' + str(100) + bcolors.ENDC)
            else:
                print(bcolors.OKGREEN + '##########   Working on step ' + str(12) + ' _100_Data_Clean  ###############' + bcolors.ENDC)
                print(bcolors.OKGREEN + str(ID) + ' Session ' + str(Session) + bcolors.ENDC)
                fonctions._100_Data_Clean.clean(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3, RS, nb_run)

            if 200 in Skip_step:
                print(bcolors.OKGREEN + 'skip step ' + str(200) + bcolors.ENDC)
            else:
                print(bcolors.OKGREEN + '##########   Working on step ' + str(12) + ' _200_Data_QC  ###############' + bcolors.ENDC)
                print(bcolors.OKGREEN + str(ID) + ' Session ' + str(Session) + bcolors.ENDC)
                fonctions._200_Data_QC._itk_check_coregistr(dir_fMRI_Refth_RS_prepro3, BASE_SS_coregistr,s_bind,itk_sif)