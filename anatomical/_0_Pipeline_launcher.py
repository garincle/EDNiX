#import
import os
import subprocess
import sys

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

######### to do
### add a verbose feature
### add an option to not redow if it already exists
### improve overwrite option
### QC
### add script to say what have been done per animal

##########################################
########### Subject loader################
##########################################

#https://bids-standard.github.io/pybids/reports/index.html


import anatomical._0_Pipeline_launcher
import anatomical._1_correct_orient
import anatomical._2_clean_anat
import anatomical._3_make_template
import anatomical._4_create_template_brain
import anatomical._5_brainT_to_stdyT
import anatomical._5_brainT_to_stdyT_max
import anatomical._6_stdyT_to_AtlasT
import anatomical._7_prepar_aseg
import anatomical._8_do_fMRImasks
import anatomical._9_nii_to_mgz
import anatomical._10_FS_1_white
import anatomical._12_make_pial
import anatomical._13_FS_freeview
import anatomical._14_Finalise
import anatomical._15_to_WB
import anatomical._100_Data_Clean
import anatomical._200_Data_QC



def preprocess_anat(BIDStype, deoblique_exeption1, deoblique_exeption2, deoblique, BASE_mask, coregistration_longitudinal, creat_study_template,
    orientation, masking_img, brain_skullstrip_1, brain_skullstrip_2, n_for_ANTS, Skip_step, check_visualy_each_img, do_manual_crop, do_fMRImasks,
    BASE_SS,BASE_bet, which_on, all_ID_max, max_session, all_data_path_max, all_ID, all_Session, all_data_path, study_template_atlas_forlder, template_skullstrip,
    IgotbothT1T2, list_atlases, Aseg_ref, Aseg_refLR, dir_out, FS_dir, do_surfacewith, Atemplate_to_Stemplate,
    FS_buckner40_TIF,FS_buckner40_GCS, Hmin, Lut_file, otheranat, type_norm, max_sessionlist, bids_dir, check_visualy_final_mask, useT1T2_for_coregis, FreeSlabel_ctab_list, list_atlases_2, cost3dAllineate, Align_img_to_template,
    species, overwrite_option,MAIN_PATH):

    sys.path.append(opj(MAIN_PATH + 'Code', 'EasyMRI_brain-master'))

    ### singularity set up
    s_bind = ' --bind ' + opj('/', 'scratch', 'in_Process/') + ',' + MAIN_PATH
    s_path = opj(MAIN_PATH, 'code', 'singularity')
    afni_sif = ' ' + opj(s_path, 'afni_make_build_24_2_01.sif') + ' '
    fsl_sif = ' ' + opj(s_path, 'fsl_6.0.5.1-cuda9.1.sif') + ' '
    fs_sif = ' ' + opj(s_path, 'freesurfer_NHP.sif') + ' '
    itk_sif = ' ' + opj(s_path, 'itksnap_5.0.9.sif') + ' '
    wb_sif = ' ' + opj(s_path, 'connectome_workbench_1.5.0-freesurfer-update.sif') + ' '

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
        print( '###################################################### work on subject' + str(ID) + ' Session ' + str(Session) + '###############################################################' +
        '############################################################################################################################################################################'+
        '############################################################################################################################################################################')
        animal_folder =   'sub-' + ID + '_ses-' + str(Session)


        # The anatomy
        path_anat     = opj(data_path,'anat')
        dir_transfo   = opj(path_anat,'matrices')
        dir_native    = opj(path_anat,'native')
        dir_prepro    = opj(dir_native,'01_preprocess')
        wb_native_dir = opj(dir_native,'02_Wb')
        volumes_dir   = opj(wb_native_dir,'volumes')
        labels_dir    = opj(volumes_dir,'labels')
        masks_dir     = opj(volumes_dir,'masks')

        BASE_SS_coregistr     = BASE_SS
        BASE_SS_mask = BASE_mask

        command = 'AFNI_NIFTI_TYPE_WARM=NO'
        spco([command], shell=True)

        DIR = os.getcwd()
        print('Working path : '+ DIR)

        ###########define orientation#############
        ##########################################

        #creat path
        if ope(dir_prepro) == False:
            os.makedirs(dir_prepro)

        if ope(dir_transfo) == False:
            os.makedirs(dir_transfo)

        #creat path
        if ope(wb_native_dir) == False:
            os.makedirs(labels_dir)
            os.makedirs(masks_dir)

        #creat path
        if ope(opj(bids_dir, 'QC')) == False:
            os.makedirs(opj(bids_dir, 'QC'))

        print(str(ID) + ' ' +  str(Session) + ' ' + str(type_norm))

        print('ok')

        if 1 in Skip_step:
            print('skip step ' + str(1))

        else:
            anatomical._1_correct_orient.correct_orient(BIDStype, listTimage, path_anat, ID, Session, otheranat, type_norm, deoblique_exeption1, deoblique_exeption2, deoblique, orientation, dir_prepro, IgotbothT1T2, overwrite,s_bind,afni_sif,fs_sif)

        if 2 in Skip_step:
            print('skip step ' + str(2))

        else:
            anatomical._2_clean_anat.clean_anat(Align_img_to_template, cost3dAllineate, bids_dir, listTimage, path_anat, ID, Session, otheranat, type_norm, deoblique_exeption1, deoblique_exeption2, deoblique, orientation, dir_prepro, masking_img, do_manual_crop,
            brain_skullstrip_1, brain_skullstrip_2, masks_dir, volumes_dir, n_for_ANTS, dir_transfo, BASE_SS_coregistr, BASE_SS_mask, BASE_SS,BASE_bet, IgotbothT1T2, check_visualy_each_img, check_visualy_final_mask, overwrite,
                                                s_bind,afni_sif,fsl_sif,fs_sif, itk_sif)

    ###### STOP THE LOOP ######

    if creat_study_template==True:

        if 3 in Skip_step:
            print('skip step ' + str(3))

        else:
            anatomical._3_make_template.make_template(which_on, all_ID_max, max_session, all_data_path_max, all_ID, all_Session, all_data_path, type_norm, study_template_atlas_forlder, template_skullstrip, BASE_SS, BASE_mask, overwrite,s_bind,afni_sif,fsl_sif,fs_sif, itk_sif)

    ###### LOOP AGAIN ######
    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, max_sessionlist):
        print('###################################################### work on subject' + str(ID) + ' Session ' + str(Session) + '###############################################################' +
            '############################################################################################################################################################################'+
            '############################################################################################################################################################################')
        animal_folder =   'sub-' + ID + '_ses-' + str(Session)


        # The anatomy
        path_anat     = opj(data_path,'anat/')
        dir_transfo   = opj(path_anat,'matrices')
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
            data_path_max     = opj(bids_dir,'sub-' + ID,'ses-' + str(max_ses))
            path_anat_max     = opj(data_path_max,'anat')
            dir_transfo_max   = opj(path_anat_max,'matrices')
            dir_native_max    = opj(path_anat_max,'native')
            wb_native_dir_max = opj(dir_native_max,'02_Wb')
            volumes_dir_max   = opj(wb_native_dir_max,'volumes')
            masks_dir_max     = opj(volumes_dir_max,'masks')

            #####!!! use the last image as template
            BASE_SS_coregistr     = opj(volumes_dir_max,ID + type_norm + '_brain.nii.gz')
            masking_img = type_norm
            BASE_SS_mask = opj(masks_dir_max, ID + masking_img + '_mask_2.nii.gz')

            transfo_concat = \
            [opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_1Warp.nii.gz'),
             opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
             opj(dir_transfo_max,'template_to_' + type_norm + '_SyN_final_max_1Warp.nii.gz'),
             opj(dir_transfo_max,'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat')]
            w2inv_fwd = [False,False,False,False]

            ####
            transfo_concat_inv = \
            [opj(dir_transfo_max,'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat'),
             opj(dir_transfo_max,'template_to_' + type_norm + '_SyN_final_max_1InverseWarp.nii.gz'),
             opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
             opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_1InverseWarp.nii.gz')]
            w2inv_inv = [True, False, True, False]


            if coregistration_longitudinal == True:
                if Session == max_ses:
                    transfo_concat = \
                        [opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_1Warp.nii.gz'),
                         opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat')]
                    w2inv_fwd = [False, False]

                    transfo_concat_inv = \
                        [opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_1InverseWarp.nii.gz'),
                         opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat')]
                    w2inv_inv = [False, True]

        ################# coregistration non longitudinal #################

        else:
            if creat_study_template == True:
                stdy_template_mask = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm,
                                         'study_template_mask.nii.gz')
                stdy_template = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm,
                                    'study_template.nii.gz')

                BASE_SS_coregistr     = stdy_template
                BASE_SS_mask = stdy_template_mask
            else:
                BASE_SS_coregistr     = BASE_SS
                BASE_SS_mask = BASE_mask

            transfo_concat = \
                [opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_1Warp.nii.gz'),
                 opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat')]
            w2inv_fwd = [False, False]

            transfo_concat_inv = \
                [opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_1InverseWarp.nii.gz'),
                 opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat')]
            w2inv_inv = [False, True]

        if useT1T2_for_coregis == True:
            Ref_file = opj(volumes_dir, ID + type_norm + '_' + otheranat + '_brain.nii.gz')
        else:
            # Organization of the folders
            Ref_file = opj(volumes_dir, ID + type_norm + '_brain.nii.gz')

        # creat a template for each individual (nice skulltrip corrected of each mean of anat img)
        if 4 in Skip_step:
            print('skip step ' + str(4))

        else:
            anatomical._4_create_template_brain.create_indiv_template_brain(dir_prepro, ID, Session, listTimage, volumes_dir, masking_img, brain_skullstrip_1, brain_skullstrip_2, 
                masks_dir, type_norm, n_for_ANTS, dir_transfo, BASE_SS_coregistr, BASE_SS_mask, otheranat, check_visualy_final_mask, useT1T2_for_coregis, bids_dir, overwrite, BASE_bet, s_bind,afni_sif,fsl_sif,fs_sif, itk_sif)


        ###### coregistration of each indiv template to the selected template (sty, atlas)
        if 5 in Skip_step:
            print('skip step ' + str(5))

        else:
            anatomical._5_brainT_to_stdyT.brainT_to_T(dir_prepro, ID, Session, listTimage, n_for_ANTS, dir_transfo, type_norm, BASE_SS_coregistr, Ref_file, volumes_dir, transfo_concat_inv,w2inv_inv,creat_study_template, which_on, all_ID_max, max_session, all_data_path_max, all_ID, all_Session, all_data_path, study_template_atlas_forlder, otheranat, bids_dir, overwrite)

        if coregistration_longitudinal == True:
            if Session == max_ses:
                if creat_study_template == True:
                    stdy_template_mask = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm,
                                             'study_template_mask.nii.gz')
                    stdy_template = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm,
                                        'study_template.nii.gz')
                    BASE_SS_coregistr     = stdy_template
                    BASE_SS_mask = stdy_template_mask
                else:
                    BASE_SS_coregistr     = BASE_SS
                    BASE_SS_mask = BASE_mask

                transfo_concat = \
                    [opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_1Warp.nii.gz'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_1Warp.nii.gz'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat')]
                w2inv_fwd = [False, False,False,False]

                ####
                transfo_concat_inv = \
                    [opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_1InverseWarp.nii.gz'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_1InverseWarp.nii.gz')]
                w2inv_inv = [True, False, True, False]

                ###### coregistration of each indiv template to the selected template (sty, atlas)
                if 5 in Skip_step:
                    print('skip step ' + str(5))

                else:
                    anatomical._5_brainT_to_stdyT_max.brainT_to_T_max(dir_prepro, ID, Session, listTimage, n_for_ANTS,
                                                                      dir_transfo, type_norm, BASE_SS_coregistr, Ref_file,
                                                                      volumes_dir, transfo_concat_inv,w2inv_inv, creat_study_template,
                                                                      which_on, all_ID_max, max_session, all_data_path_max,
                                                                      all_ID, all_Session, all_data_path,
                                                                      study_template_atlas_forlder, otheranat, bids_dir,
                                                                      overwrite)
    ###### STOP THE LOOP ######

    
    if creat_study_template==True:
        stdy_template_mask = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm,
                                 'study_template_mask.nii.gz')
        stdy_template = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm,
                            'study_template.nii.gz')

        ###################redefine new atlases variable!!!

        if 6 in Skip_step:
            print('skip step ' + str(6))
        else:
            anatomical._6_stdyT_to_AtlasT.stdyT_to_AtlasT(list_atlases, Aseg_ref, Aseg_refLR, BASE_SS, dir_out, n_for_ANTS, study_template_atlas_forlder, Atemplate_to_Stemplate, overwrite,
                                                          s_bind,afni_sif)

        ###### re-define the variable: study template atlas, etc should now be the new template !!
        list_atlases3 = list_atlases
        list_atlases = []
        for atlas_new in list_atlases3:
            atlas2 = opj(dir_out, opb(atlas_new))
            list_atlases.append(atlas2)

        BASE_SS2 = BASE_SS
        BASE_mask2 = BASE_mask
        Aseg_ref2 = Aseg_ref
        Aseg_refLR2 = Aseg_refLR
        BASE_SS     = opj(dir_out, opb(BASE_SS2))
        BASE_mask   = opj(dir_out, opb(BASE_mask2))

        ####atlases files
        if ope(Aseg_ref2):
            Aseg_ref = opj(dir_out, opb(Aseg_ref2))
        else:
            Aseg_ref = ''
        if ope(Aseg_ref2):
            Aseg_refLR = opj(dir_out, opb(Aseg_refLR2))
        else:
            Aseg_refLR = ''

    ###### LOOP AGAIN ######

    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, max_sessionlist):
        print( '###################################################### work on subject' + str(ID) + ' Session ' + str(Session) + '###############################################################' +
        '############################################################################################################################################################################'+
        '############################################################################################################################################################################' )
        animal_folder =   'sub-' + ID + '_ses-' + str(Session)

        # The anatomy
        path_anat    = opj(data_path,'anat/')
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
            if creat_study_template == True:
                stdy_template_mask = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm,
                                         'study_template_mask.nii.gz')
                stdy_template = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm,
                                    'study_template.nii.gz')
                BASE_SS_coregistr     = stdy_template
                BASE_SS_mask = stdy_template_mask
            else:
                BASE_SS_coregistr     = BASE_SS
                BASE_SS_mask = BASE_mask
                
            if Session == max_ses:
                transfo_concat = \
                    [opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_1Warp.nii.gz'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_1Warp.nii.gz'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat')]
                w2inv_fwd = [False, False, False, False]

            else:
                data_path_max = opj(bids_dir,'sub-' + ID,'ses-' + str(max_ses))
                path_anat_max    = opj(data_path_max,'anat/')
                dir_transfo_max  = opj(path_anat_max,'matrices')
                dir_native_max    = opj(path_anat_max,'native')
                wb_native_dir_max = opj(dir_native_max,'02_Wb')
                volumes_dir_max   = opj(wb_native_dir_max,'volumes')
                masks_dir_max     = opj(volumes_dir_max,'masks')

                transfo_concat = \
                    [opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_1Warp.nii.gz'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
                     opj(dir_transfo_max,'template_to_' + type_norm + '_SyN_final_max_1Warp.nii.gz'),
                     opj(dir_transfo_max,'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat')]
                w2inv_fwd = [False, False, False, False]

        ################# coregistration non longitudinal #################

        else:
            if creat_study_template == True:
                stdy_template_mask = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm,
                                         'study_template_mask.nii.gz')
                stdy_template = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm,
                                    'study_template.nii.gz')
                BASE_SS_coregistr     = stdy_template
                BASE_SS_mask = stdy_template_mask
            else:
                BASE_SS_coregistr     = BASE_SS
                BASE_SS_mask = BASE_mask

            transfo_concat = \
                [opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_1Warp.nii.gz'),
                 opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat')]
            w2inv_fwd = [False, False]

        if useT1T2_for_coregis == True:
            Ref_file = opj(volumes_dir,ID + type_norm + '_' + otheranat + '_brain.nii.gz')
        else:
        # Organization of the folders
            Ref_file = opj(volumes_dir,ID + type_norm + '_brain.nii.gz')


        if 7 in Skip_step:
            print('skip step ' + str(7))
        else:
            anatomical._7_prepar_aseg.prepar_aseg(Ref_file, labels_dir, volumes_dir, masks_dir, dir_transfo, BASE_SS_mask, BASE_SS_coregistr, Aseg_refLR, Aseg_ref, type_norm, ID, transfo_concat,w2inv_fwd, dir_prepro, list_atlases, check_visualy_each_img, n_for_ANTS, overwrite,
                                                  s_bind,afni_sif,itk_sif)
        
        if do_fMRImasks == True:

            if 8 in Skip_step:
                print('skip step ' + str(8))

            else:
                ########################## Building fMRI masks for EPI analysis ##############################
                anatomical._8_do_fMRImasks.do_fMRImasks(masks_dir, labels_dir, type_norm, overwrite,
                                                        s_bind,afni_sif,)

        if 9 in Skip_step:
            print('skip step ' + str(9))

        else:
            ########################## White Surface construction ##############################
            # You can go grab a cup of coffe, it can take more than an hour...
            anatomical._9_nii_to_mgz.nii_to_mgz(ID, Session, FS_dir, Ref_file, labels_dir, volumes_dir, otheranat, IgotbothT1T2, type_norm, overwrite,
                                                s_bind,fs_sif)

        if 10 in Skip_step:
            print('skip step ' + str(10))

        else:
            ########################## White Surface construction ##############################
            anatomical._10_FS_1_white.White_create(FS_dir, animal_folder,s_bind,fs_sif)

        if 11 in Skip_step:
            print('skip step ' + str(11))

        else:
            ########################## White Surface construction ##############################
            # You can go grab a cup of coffe, it can take more than an hour...
            anatomical._10_FS_1_white.White_more(FS_dir, animal_folder, FS_buckner40_TIF,FS_buckner40_GCS,s_bind,fs_sif)

        if 12 in Skip_step:
            print('skip step ' + str(12))

        else:
            anatomical._12_make_pial.make_pial(FS_dir, animal_folder, type_norm, otheranat, Hmin, Ref_file, do_surfacewith, overwrite,
                                               s_bind,fs_sif)


        if check_visualy_each_img == True:
            if 13 in Skip_step:
                print('skip step ' + str(13))

            else:
                anatomical._13_FS_freeview.FS_Freeview(FS_dir, animal_folder, 'pial', Lut_file,
                                                       s_bind,fs_sif)

        if 14 in Skip_step:
            print('skip step ' + str(14))

        else:
            anatomical._14_Finalise.FS_finalise(FS_dir, animal_folder, FreeSlabel_ctab_list, list_atlases_2, labels_dir, type_norm, Ref_file,
                                                s_bind,fs_sif)
        if 15 in Skip_step:
            print('skip step ' + str(15))

        else:
            anatomical._15_to_WB.WB_prep(FS_dir, dir_native, animal_folder, Ref_file, species, list_atlases_2,s_bind,afni_sif,fsl_sif,fs_sif,wb_sif)

        if 100 in Skip_step:
                print('skip step ' + str(100))
        else:
            anatomical._100_Data_Clean.clean(all_ID, all_Session, all_data_path)

        if 200 in Skip_step:
            print('skip step ' + str(200))

        else:
            anatomical._200_Data_QC._itk_check_masks(dir_prepro, masks_dir, ID, type_norm,s_bind, itk_sif)