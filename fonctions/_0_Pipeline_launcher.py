#import
import os
import subprocess
import glob
import json
import sys

################################### if re-use this script auto: ####################################################
##### Trinity Session 6 (T1 or T2) and Unity Session take the other T1. Otherwise: anat image not in same space ####
########################################v###########################################################################

##########################################
########### Subject loader################
##########################################

#https://bids-standard.github.io/pybids/reports/index.html
#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
ops = os.path.splitext
spco = subprocess.check_output
spgo = subprocess.getoutput



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
import fonctions._13_spatial_QC
import fonctions._14_fMRI_QC_SBA
import fonctions._100_Data_Clean



##### to ask or find solution ####
##### if change the orientation should we change the  phase encoding direction as well? ####
def preprocess_data(all_ID, all_Session, all_data_path, max_sessionlist, stdy_template, stdy_template_mask, BASE_SS, BASE_mask, T1_eq, anat_func_same_space, 
    correction_direction, REF_int, study_fMRI_Refth, IgotbothT1T2, otheranat, deoblique_exeption1, deoblique_exeption2, deoblique, orientation,
    TfMRI, GM_mask_studyT, GM_mask, creat_study_template, type_norm, coregistration_longitudinal, dilate_mask, overwrite_option, nb_ICA_run, blur, melodic_prior_post_TTT,
    extract_exterior_CSF, extract_WM, n_for_ANTS, list_atlases, selected_atlases, panda_files, useT1T2_for_coregis, endfmri, endjson, endmap, oversample_map, use_cortical_mask_func,
    cut_coordsX, cut_coordsY, cut_coordsZ, threshold_val, Skip_step, bids_dir, costAllin, use_erode_WM_func_masks, do_not_correct_signal, use_erode_V_func_masks,
    folderforTemplate_Anat, IhaveanANAT, doMaskingfMRI, do_anat_to_func, Method_mask_func, segmentation_name_list, band, extract_Vc, lower_cutoff, upper_cutoff, selected_atlases_matrix, specific_roi_tresh, unspecific_ROI_thresh, Seed_name, extract_GS, MAIN_PATH):

    sys.path.append(opj(MAIN_PATH + 'Code', 'EasyMRI_brain-master'))

    ### singularity set up
    s_bind = ' --bind ' + opj('/', 'scratch', 'in_Process/') + ',' + MAIN_PATH
    s_path = opj(MAIN_PATH, 'code', 'singularity')
    afni_sif = ' ' + opj(s_path, 'afni_make_build_24_2_01.sif') + ' '
    fsl_sif = ' ' + opj(s_path, 'fsl_6.0.5.1-cuda9.1.sif') + ' '
    fs_sif = ' ' + opj(s_path, 'freesurfer_NHP.sif') + ' '
    itk_sif = ' ' + opj(s_path, 'itksnap_5.0.9.sif') + ' '
    wb_sif = ' ' + opj(s_path, 'connectome_workbench_1.5.0-freesurfer-update.sif') + ' '

    config_f = opj(MAIN_PATH, 'code', 'config', 'b02b0.cnf')

    ###########################################################################################################################################################
    ############################################################## start the proces ###########################################################################
    ###########################################################################################################################################################
    ######### define other usefull paramater automatically (do no touch)#########

    if overwrite_option == True:
        overwrite = ' -overwrite'
    else:
        overwrite = ''

    listTimage = []
    if IgotbothT1T2 == True:
        listTimage = [otheranat, type_norm]
    else: 
        listTimage = [type_norm]

    ####################################################################################
    ########################## Start the pipeline !!!!!!!!!!!!!!!!!!!!!!   #############
    ####################################################################################

    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, max_sessionlist):
        print('work on ' + str(ID) + ' Session ' + str(Session))

        if IhaveanANAT==False:
            # The anatomy
            path_anat    = ''
            dir_transfo  = ''

            dir_native    = ''
            dir_prepro    = ''
            wb_native_dir = ''
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

                max_ses_anat = max_ses
                Session_anat = Session


            else:
                template_anat_for_fmri = glob.glob(opj(opd(data_path),'**','anat','native','02_Wb', 'volumes', '*' + TfMRI + '_brain.nii*'))
                print(opj(opd(data_path),'**','anat','native','02_Wb', 'volumes', '*' + TfMRI + '_brain.nii*'))
                print(template_anat_for_fmri)

                if len(template_anat_for_fmri) == 1:
                    data_path_anat = opd(opd(opd(opd(opd(template_anat_for_fmri[0])))))
                elif len(template_anat_for_fmri) > 1:
                    print("we found multiple anat template for this animal, please choose one!")
                    data_path_anat = input("Please enter manually a data_path_anat for preprocessing:")
                elif len(template_anat_for_fmri) == 0:
                    print("ERROR: We havn't found any anat template for this animal!!!! We can't continue !!! please provid at least one anat image!")
                    print(template_anat_for_fmri + 'not found!!!')

                # The anatomy
                path_anat    = opj(data_path_anat,'anat')
                dir_transfo  = opj(path_anat,'matrices')

                dir_native    = opj(path_anat,'native')
                dir_prepro    = opj(dir_native,'01_preprocess')
                wb_native_dir = opj(dir_native,'02_Wb')
                volumes_dir   = opj(wb_native_dir,'volumes')
                labels_dir    = opj(volumes_dir,'labels')
                masks_dir     = opj(volumes_dir,'masks')


        ####################################################################################
        ########################## Coregistration template to anat #########################
        ####################################################################################

        ################# coregistration longitudinal ???? #################

        if coregistration_longitudinal==True:
            if Session_anat == max_ses_anat:

                transfo_concat = \
                    [opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
                    opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_1Warp.nii.gz'),
                    opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat'),
                    opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_1Warp.nii.gz')]
                w2inv_fwd = [False,False,False,False]

                ####
                transfo_concat_inv = \
                    [opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_1InverseWarp.nii.gz'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_1InverseWarp.nii.gz'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat')]
                w2inv_inv = [False, True, False, True]

            else:
                data_path_max = opj(bids_dir,'sub-' + ID,'ses-' + str(max_ses_anat))
                path_anat_max    = opj(data_path_max,'anat')
                dir_transfo_max  = opj(path_anat_max,'matrices')
                dir_native_max    = opj(path_anat_max,'native')
                wb_native_dir_max = opj(dir_native_max,'02_Wb')
                volumes_dir_max   = opj(wb_native_dir_max,'volumes')
                masks_dir_max     = opj(volumes_dir_max,'masks')

                transfo_concat = \
                    [opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
                     opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_1Warp.nii.gz'),
                     opj(dir_transfo_max, 'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat'),
                     opj(dir_transfo_max, 'template_to_' + type_norm + '_SyN_final_max_1Warp.nii.gz')]

                w2inv_fwd = [False,False,False,False]

                #### if doesn't work change the order
                transfo_concat_inv = \
                    [opj(dir_transfo_max,'template_to_' + type_norm + '_SyN_final_max_1InverseWarp.nii.gz'),
                     opj(dir_transfo_max,'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_1InverseWarp.nii.gz'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat')]
                w2inv_inv = [False, True, False, True]


            if creat_study_template == True:
                BASE_SS_coregistr     = stdy_template
                BASE_SS_mask = stdy_template_mask
            else:
                BASE_SS_coregistr     = BASE_SS
                BASE_SS_mask = BASE_mask


        ################# coregistration non longitudinal #################

        else:
            if creat_study_template == True:
                BASE_SS_coregistr     = stdy_template
                BASE_SS_mask = stdy_template_mask
            else:
                BASE_SS_coregistr     = BASE_SS
                BASE_SS_mask = BASE_mask

            transfo_concat = \
                [opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
                 opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_1Warp.nii.gz')]
            w2inv_fwd = [False, False]

            transfo_concat_inv = \
                [opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_1InverseWarp.nii.gz'),
                 opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat')]
            w2inv_inv = [False, True]

        ####################################
        ########### choose Ref_file ########
        ####################################
        if useT1T2_for_coregis == True:
            Ref_file = opj(volumes_dir,ID + type_norm + '_' + otheranat + '_brain.nii.gz')
            ###creat brain of the subject template
            command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + opj(volumes_dir, ID + '_' + type_norm + '_' + otheranat + '_template.nii.gz') + \
            ' -b ' + output_for_mask + \
            ' -prefix ' + Ref_file + ' -expr "a*b"'
            spco([command], shell=True)
        else:
        # Organization of the folders
            Ref_file = opj(volumes_dir,ID + type_norm + '_brain.nii.gz')

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
            G_mask        = opj(labels_dir, type_norm + 'GM_select.nii.gz')
        else:
            anat_subject = opj(folderforTemplate_Anat,'template.nii.gz')
            brainmask     = opj(folderforTemplate_Anat,'brainmask.nii.gz')
            V_mask        = opj(folderforTemplate_Anat,'V_mask.nii.gz')
            W_mask = opj(folderforTemplate_Anat,'W_mask.nii.gz')
            G_mask = opj(folderforTemplate_Anat,'G_mask.nii.gz')

        # Resting data;
        dir_fMRI_Refth_RS          = opj(data_path,'func')
        dir_fMRI_Refth_RS_prepro   = opj(dir_fMRI_Refth_RS,'01_prepro')
        dir_fMRI_Refth_RS_prepro1  = opj(dir_fMRI_Refth_RS_prepro,'01_funcspace')
        dir_fMRI_Refth_RS_prepro2  = opj(dir_fMRI_Refth_RS_prepro,'02_anatspace')
        dir_fMRI_Refth_RS_prepro3  = opj(dir_fMRI_Refth_RS_prepro,'03_atlas_space')

        dir_fMRI_Refth_RS_residual = opj(dir_fMRI_Refth_RS,'02_residual')
        dir_RS_ICA_native    = opj(dir_fMRI_Refth_RS_residual,'01_ICA_native')
        dir_RS_ICA_native_PreTT    = opj(dir_fMRI_Refth_RS_residual,'01_ICA_native_PreTT')
        dir_RS_norm          = opj(dir_fMRI_Refth_RS_residual,'02_ICA_Norm')

        dir_fMRI_Refth_map          = opj(data_path,'fmap')


        ############################################
        ####find an auto way to defin it XXX #######
        ############################################

        # get useful informations
        list_RS = sorted(glob.glob(opj(dir_fMRI_Refth_RS, endfmri)))
        list_json = sorted(glob.glob(opj(dir_fMRI_Refth_RS, endjson)))

        if len(list_RS) == 0:
            print('ERROR : No func image found')
            print(opj(dir_fMRI_Refth_RS, endfmri))
            print(list_RS)

        nb_run  = len(list_RS)
        RS      = [os.path.basename(i) for i in list_RS]

        list_map = sorted(glob.glob(opj(dir_fMRI_Refth_map, endmap))) #[]

        if len(list_RS)>0:
            RS_map   = [os.path.basename(i) for i in list_map]

            ref_nom_IRMf = list_RS[REF_int]

            ############################################
            #### choose TOPUP stragtegy ################
            ############################################

            if len(list_map) == 0:
                recordings = 'very_old'
                print('WARNING : Before moving on, chech the quality of the AP image, you may decide to NOT use it for correction')

            elif len(list_map) == 1:
                cmd = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + \
                      '3dinfo -same_grid ' + opj(dir_fMRI_Refth_map, RS_map[0]) + ' ' + opj(dir_fMRI_Refth_RS, RS[int(REF_int)-1])
                dummy = spgo(cmd).split('\n')
                grid = int(dummy[-1])

                if int(grid) == 0:
                    recordings = 'very_old' # the only one AP recording to correct for field distorsion is useless !!!
                    print('WARNING : Before moving on, chech the quality of the AP image, you may decide to NOT use it for correction')

                elif int(grid[-3:][0]) == 1:
                    recordings = 'old'    # there is only one AP recording to correct for field distorsion

            elif len(list_map) == 2:
                recordings = '2_mapdir'    # there one AP per PA recordings in total

            elif len(list_map) > 2:
                recordings = 'new' # there is one AP per PA recordings to correct for field distorsion
                if not len(list_map) == len(list_RS):
                   raise NameError('WARNING : Check the runs there is probably one broken file that has been repeated and that you should remove !')


            print('######################################################## recordings type detected: ' + str(recordings) + '##############################################')
            f = open(list_json[0])
            info_RS = json.load(f)

            #nslice?
            try:
                nslice = int(len(info_RS["SliceTiming"]))
            except:
                nslice = int(len(info_RS["SlicePosition"]))

            ## TR
            try:
                TR     = info_RS["RepetitionTime"]
            except:
                # Calculate the time difference
                slice_timing = info_RS["SliceTiming"]
                slice_timing.sort()
                slice_intervals = [slice_timing[i + 1] - slice_timing[i] for i in range(len(slice_timing) - 1)]

                # Calculate the TR
                TR = ((sum(slice_intervals)/(nslice-1))*1000)*nslice
                print(f"######################################################################################## WARNING !!!! TR not found in Header file!!!!! Repetition  Time (TR) calculated: {TR} seconds. YOU ABSOLUTELY NEED TO DOUBLE CHECK THAT!!!!!!!!!!!!!!!!!########################################################################################")

            ## slice_timing
            try:
                slice_timing = info_RS["SliceTiming"]
            except:
                print("Slice Timing not found!!!")
            try:
                TE     = info_RS["EchoTime"]
            except:
                print("EchoTime not found in header")
            try:
                EES    = info_RS["EffectiveEchoSpacing"]
            except:
                print("Effective Echo Spacing not found in header")
            try:
                TRT = info_RS['TotalReadoutTime']
            except:
                print("Total Readout Time not found in header")

            if correction_direction:
                ('INFO: input correction_direction is the launcher was determined as' + correction_direction)
            else:
                ('INFO: input correction_direction was empty, let s try to find what is with the header')
                try:
                    PE_d2 = info_RS['PhaseEncodingDirection']
                    if PE_d2 == 'j':
                        dmap = '0 1 0 ' + str(TRT)
                        dbold = '0 -1 0 ' + str(TRT)
                        correction_direction = 'y-'
                    elif PE_d2 == 'j-':
                        dmap = '0 -1 0 ' + str(TRT)
                        dbold = '0 1 0 ' + str(TRT)
                        correction_direction = 'y'
                    elif PE_d2 == 'i':
                        dmap = '1 0 0 ' + str(TRT)
                        dbold = '-1 0 0 ' + str(TRT)
                        correction_direction = 'x-'
                    elif PE_d2 == 'i-':
                        dmap = '-1 0 0 ' + str(TRT)
                        dbold = '1 0 0 ' + str(TRT)
                        correction_direction = 'x'
                except:
                    print("INFO: Phase Encoding Direction not found in header")
                    dmap=''
                    dbold=''

                    print("WARNING :Phase Encoding Direction not found in header and you didn't provided any")
                    recordings = 'very_old'
                    print('WARNING : recordings = very_old no distortion correction will be applied with fugue')


            SED = 'i'  # or initialize it with some default value if needed
            try:
                SED = info_RS["SliceEncodingDirection"]
            except KeyError:
                try:
                    IOPD = info_RS["ImageOrientationPatientDICOM"]
                    if info_RS["ImageOrientationPatientDICOM"][0]   == 1:  SED = "i"
                    elif info_RS["ImageOrientationPatientDICOM"][0] == -1: SED = "i-"
                    elif info_RS["ImageOrientationPatientDICOM"][1] == 1:  SED = "j"
                    elif info_RS["ImageOrientationPatientDICOM"][1] == -1: SED = "j-"
                    elif info_RS["ImageOrientationPatientDICOM"][2] == 1:  SED = "k"
                    elif info_RS["ImageOrientationPatientDICOM"][2] == -1: SED = "k-"
                except KeyError:
                    print('WARNING !!!! Can not find SliceEncodingDirection in the DICOM, applied i as default values!!!!')

                    '''
                    while True:
                        SED = input("Please define a variable (i, j, or k) for SliceEncodingDirection "
                                    "since it was not found in the .json file: ")
                        if SED in ["i", "i-", "j", "j-", "k", "k-"]:
                            break
                        else:
                            print("Invalid input. Please enter either i, j, or k.")
                    '''

            #### to check
            try:
                DwellT    = "%.16f" % (float(info_RS["DwellTime"]))
            except:
                DwellT   = "%.16f" % (float(info_RS["TotalReadOutTimeEPI"] / nslice))

            ############################################
            ####marche pas XXX #######
            ############################################
            '''
            if info_RS["PatientPosition"]   == 'FFS':
                    orientation = 'RSA'
            elif info_RS["PatientPosition"] == 'HFS':
                    orientation = 'RSP'
            print('Orientation : ' + orientation)
            '''


            DIR = os.getcwd()
            print('Working path : '+ DIR)

            if 1 in Skip_step:
                print('skip step ' + str(1))
            else:
                print('##########   Working on step ' + str(1) + ' _1_fMRI_preTTT_in_fMRIspace  ###############')
                print(str(ID) + ' Session ' + str(Session))
                fonctions._1_fMRI_preTTT_in_fMRIspace.preprocess_data(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, RS, list_RS, nb_run, T1_eq, overwrite,
                                                                      s_bind,afni_sif)


            if 2 in Skip_step:
                print('skip step ' + str(2))
            else:
                print('##########   Working on step ' + str(2) + ' _2_coregistration_to_norm  ###############')
                print(str(ID) + ' Session ' + str(Session))

                fonctions._2_coregistration_to_norm.coregist_to_norm(TR, anat_func_same_space, dir_prepro, correction_direction, dir_fMRI_Refth_RS_prepro1, RS, RS_map, nb_run, recordings,
                    REF_int, list_map, study_fMRI_Refth, IgotbothT1T2, path_anat, otheranat, ID, Session, deoblique_exeption1, deoblique_exeption2, deoblique, orientation, DwellT, n_for_ANTS, overwrite,
                                                                     s_bind,afni_sif,fsl_sif,dmap,dbold,config_f)

            if 3 in Skip_step:
                print('skip step ' + str(3))
            else:
                print('##########   Working on step ' + str(3) + ' _3_mask_fMRI  ###############')
                print(str(ID) + ' Session ' + str(Session))
                fonctions._3_mask_fMRI.Refimg_to_meanfMRI(SED, anat_func_same_space, BASE_SS_coregistr, TfMRI, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                    dir_fMRI_Refth_RS_prepro3, RS, nb_run, REF_int, ID, dir_prepro, n_for_ANTS, brainmask, V_mask, W_mask, G_mask, dilate_mask,
                                                             list_atlases, labels_dir, costAllin, anat_subject, IhaveanANAT, doMaskingfMRI, do_anat_to_func, Method_mask_func, lower_cutoff, upper_cutoff, overwrite,
                                                          s_bind,afni_sif,fs_sif)
            if 4 in Skip_step:
                print('skip step ' + str(4))
            else:
                print('##########   Working on step ' + str(4) + ' _itk_check_masks  ###############')
                print(str(ID) + ' Session ' + str(Session))
                fonctions._4_check_mask._itk_check_masks(dir_fMRI_Refth_RS_prepro1,s_bind,itk_sif)

            if 5 in Skip_step:
                print('skip step ' + str(5))
            else:
                print('##########   Working on step ' + str(5) + ' _5_anat_to_fMRI  ###############')
                print(str(ID) + ' Session ' + str(Session))
                fonctions._5_anat_to_fMRI.Refimg_to_meanfMRI(SED, anat_func_same_space, BASE_SS_coregistr, TfMRI, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                    dir_fMRI_Refth_RS_prepro3, RS, nb_run, REF_int, ID, dir_prepro, n_for_ANTS, brainmask, V_mask, W_mask, G_mask, dilate_mask,
                                                             list_atlases, labels_dir, costAllin, anat_subject, IhaveanANAT, doMaskingfMRI, do_anat_to_func, Method_mask_func, lower_cutoff, upper_cutoff, overwrite,
                                                             s_bind,afni_sif)
            if 6 in Skip_step:
                print('skip step ' + str(6))
            else:
                if melodic_prior_post_TTT == True:
                    print('##########   Working on step ' + str(6) + ' Melodic_correct  ###############')
                    print(str(ID) + ' Session ' + str(Session))
                    fonctions._6_Melodic.Melodic_correct(dir_RS_ICA_native_PreTT, dir_RS_ICA_native, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                nb_ICA_run, nb_run, RS, TfMRI, overwrite,s_bind,fsl_sif,itk_sif,TR)

            if 7 in Skip_step:
                print('skip step ' + str(7))
            else:
                print('##########   Working on step ' + str(7) + ' _7_post_TTT  ###############')
                print(str(ID) + ' Session ' + str(Session))
                fonctions._7_post_TTT.signal_regression(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_RS_ICA_native,
                nb_run, RS, blur, TR, melodic_prior_post_TTT, extract_exterior_CSF, extract_WM, do_not_correct_signal, band, extract_Vc, extract_GS, overwrite,s_bind,afni_sif)

            if 8 in Skip_step:
                print('skip step ' + str(8))
            else:
                print('##########   Working on step ' + str(8) + ' _8_fMRI_to_anat  ###############')
                print(str(ID) + ' Session ' + str(Session))
                fonctions._8_fMRI_to_anat.to_anat_space(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                nb_run, RS, n_for_ANTS, TR, overwrite)

            if 9 in Skip_step:
                print('skip step ' + str(9))
            else:
                print('##########   Working on step ' + str(9) + ' _9_coregistration_to_template_space  ###############')
                print(str(ID) + ' Session ' + str(Session))
                fonctions._9_coregistration_to_template_space.to_common_template_space(Session, deoblique_exeption1, deoblique_exeption2, deoblique, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3, BASE_SS_coregistr,
                    nb_run, RS, transfo_concat,w2inv_fwd,do_anat_to_func, TR, n_for_ANTS, list_atlases, TfMRI, BASE_SS_mask, GM_mask, GM_mask_studyT, creat_study_template, anat_func_same_space, orientation, path_anat, ID, REF_int, dir_prepro, IhaveanANAT, overwrite,s_bind,afni_sif)

            if 10 in Skip_step:
                print('skip step ' + str(10))
            else:
                print('##########   Working on step ' + str(10) + ' _10_Correl_matrix  ###############')
                print(str(ID) + ' Session ' + str(Session))
                fonctions._10_Correl_matrix.correl_matrix(dir_fMRI_Refth_RS_prepro1, RS, nb_run, selected_atlases_matrix, segmentation_name_list, ID, Session, bids_dir,s_bind,afni_sif)

            if 11 in Skip_step:
                print('skip step ' + str(11))
            else:
                print('##########   Working on step ' + str(11) + ' _11_Seed_base_many_regionsatlas  ###############')
                print(str(ID) + ' Session ' + str(Session))
                fonctions._11_Seed_base_many_regionsatlas.SBA(volumes_dir, BASE_SS_coregistr, TfMRI, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                dir_fMRI_Refth_RS_prepro3, RS, nb_run, ID, selected_atlases, panda_files, oversample_map, use_cortical_mask_func, cut_coordsX, cut_coordsY, cut_coordsZ, threshold_val, overwrite,s_bind,afni_sif)

            if 12 in Skip_step:
                print('skip step ' + str(12))
            else:
                print('##########   Working on step ' + str(12) + ' _12_fMRI_QC  ###############')
                print(str(ID) + ' Session ' + str(Session))
                fonctions._12_fMRI_QC.fMRI_QC(ID, Session, segmentation_name_list, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3, specific_roi_tresh, unspecific_ROI_thresh, RS, nb_run, bids_dir,s_bind,afni_sif)

            if 13 in Skip_step:
                print('skip step ' + str(13))
            else:
                print('##########   Working on step ' + str(13) + ' _itk_check_masks  ###############')
                print(str(ID) + ' Session ' + str(Session))
                fonctions._13_spatial_QC._itk_check_spatial_co(dir_fMRI_Refth_RS_prepro3,s_bind,itk_sif)

            if 14 in Skip_step:
                print('skip step ' + str(14))
            else:
                fonctions._14_fMRI_QC_SBA.fMRI_QC_SBA(Seed_name, BASE_SS_coregistr, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                            dir_fMRI_Refth_RS_prepro3, RS, nb_run, selected_atlases, panda_files, oversample_map,
                            use_cortical_mask_func)

            if 100 in Skip_step:
                print('skip step ' + str(100))
            else:
                fonctions._100_Data_Clean.clean(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3, RS, nb_run)