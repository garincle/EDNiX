#import
import os
import subprocess
import datetime

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

import anatomical._1_correct_orient
import anatomical._2_clean_anat
import anatomical._3_make_template
import anatomical._4_skullstrip_template
import anatomical._5_create_template_brain
import anatomical._6_brainT_to_stdyT
import anatomical._6_brainT_to_stdyT_max
import anatomical._7_stdyT_to_AtlasT
import anatomical._8_prepar_aseg
import anatomical._9_do_fMRImasks
import anatomical._10_nii_to_mgz
import anatomical._11_FS_1_white
import anatomical._12_make_pial
import anatomical._13_FS_freeview
import anatomical._14_Finalise
import anatomical._15_to_WB
import anatomical._16_anat_QC_SNR
import anatomical._100_Data_Clean
import anatomical._200_Data_QC
import Tools.Load_EDNiX_requirement

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

def preprocess_anat(BIDStype, deoblique, BASE_mask, coregistration_longitudinal, creat_study_template,
    orientation, masking_img, brain_skullstrip_1, brain_skullstrip_2, n_for_ANTS, aff_metric_ants, Skip_step,
    check_visualy_each_img, do_fMRImasks, BASE_SS, which_on, all_ID_max, all_data_path_max, all_ID,
    all_Session, all_data_path, template_skullstrip, list_atlases, Aseg_ref, Aseg_refLR, FS_dir,
    do_surfacewith, Atemplate_to_Stemplate, FS_buckner40_TIF,FS_buckner40_GCS, Lut_file, otheranat,
    type_norm, all_Session_max, bids_dir, check_visualy_final_mask, FreeSlabel_ctab_list,
    list_atlases_2, cost3dAllineate, Align_img_to_template, species, type_of_transform,
    type_of_transform_stdyT, fMRImasks, overwrite_option, MAIN_PATH, aff_metric_ants_Transl, aff_metric_ants_Transl_template):

    ### singularity set up
    Hmin = ['l', 'r']
    s_path, afni_sif, fsl_sif, fs_sif, itk_sif, wb_sif, strip_sif, s_bind =  Tools.Load_EDNiX_requirement.load_requirement(MAIN_PATH, bids_dir, FS_dir)

    ###########################################################################################################################################################
    ############################################################## start the proces ###########################################################################
    ###########################################################################################################################################################
    ######### define other usefull paramater automatically (do no touch)#########

    if overwrite_option == True:
        overwrite = ' -overwrite'
    else:
        overwrite = ''

    IgotbothT1T2 = bool(otheranat.strip())
    if IgotbothT1T2 == True:
        listTimage = [otheranat, type_norm]
    else: 
        listTimage = [type_norm]

    ####################################################################################
    ########################## Start the pipeline !!!!!!!!!!!!!!!!!!!!!!   #############
    ####################################################################################

    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, all_Session_max):
        NL1 = '###################################################### work on subject: ' + str(ID) + ' Session ' + str(Session) + ' BLOCK 1 ###############################################################'
        print(bcolors.HEADER + NL1 + bcolors.ENDC)
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
        # folder where you want to store the study template
        study_template_atlas_folder = bids_dir + '/sty_template'
        # then where do you want your atlases in sty template to be
        dir_out = bids_dir + '/sty_template/atlases'

        # Set the environment variable for the current process
        os.environ["AFNI_NIFTI_TYPE_WARN"] = "NO"

        DIR = os.getcwd()
        print(bcolors.OKGREEN + 'INFO: Working path : ' + DIR + bcolors.ENDC)

        # create the folder's architecture
        if ope(dir_prepro) == False:
            os.makedirs(dir_prepro)
        if ope(dir_transfo) == False:
            os.makedirs(dir_transfo)
        if ope(labels_dir) == False:
            os.makedirs(labels_dir)
        if ope(masks_dir) == False:
                os.makedirs(masks_dir)
        if ope(opj(bids_dir, 'QC')) == False:
            os.makedirs(opj(bids_dir, 'QC'))

        date_file = datetime.date.today()
        ct = datetime.datetime.now()
        diary_name = str(ID) + ' session ' + str(Session) + str(date_file) + '.txt'
        diary_file = opj(path_anat, diary_name)
        if not opi(diary_file):
            diary = open(diary_file, "w")

            diary.write(NL1)
        else :
            diary = open(diary_file, "a")
            diary.write(NL1)
        diary.write(f'\n{ct}')
        diary.write(f'\n')
        diary.close()

        if 1 in Skip_step:
            print(bcolors.OKGREEN + 'INFO: skip step ' + str(1) + bcolors.ENDC)
        else:
            anatomical._1_correct_orient.correct_orient(BIDStype, listTimage, path_anat, ID, Session, otheranat, type_norm, deoblique, orientation, dir_prepro, IgotbothT1T2, overwrite,s_bind,afni_sif,fs_sif,diary_file)

        if 2 in Skip_step:
            print(bcolors.OKGREEN + 'INFO: skip step ' + str(2) + bcolors.ENDC)
        else:
            anatomical._2_clean_anat.clean_anat(Align_img_to_template, cost3dAllineate, bids_dir, listTimage, type_of_transform, ID, aff_metric_ants, Session, otheranat,
               type_norm, dir_prepro, masking_img, brain_skullstrip_1, brain_skullstrip_2, masks_dir, volumes_dir, BASE_SS_coregistr,
               BASE_SS_mask, BASE_SS, IgotbothT1T2, check_visualy_each_img, check_visualy_final_mask, template_skullstrip,
               study_template_atlas_folder, overwrite,
               s_bind,afni_sif,fsl_sif,fs_sif, itk_sif, strip_sif,diary_file)

    ###### STOP THE LOOP ######
    if creat_study_template==True:
        if not ope(study_template_atlas_folder): os.mkdir(study_template_atlas_folder)
        if not ope(dir_out): os.mkdir(dir_out)
        ct = datetime.datetime.now()
        diary_name = 'Study_template_BLOCK1.txt'
        diary_file = opj(study_template_atlas_folder, diary_name)
        if not opi(diary_file):
            diary = open(diary_file, "w")
            diary.write('Create the study template')
        else:
            diary = open(diary_file, "a")
        diary.write(f'\n{ct}')
        diary.write(f'\n')
        diary.close()

        if 3 in Skip_step:
            print(bcolors.OKGREEN + 'INFO: skip step ' + str(3) + bcolors.ENDC)
        else:
            anatomical._3_make_template.make_template(which_on, all_ID_max, all_Session_max, all_data_path_max, all_ID, all_Session, all_data_path, type_norm, study_template_atlas_folder,
                  s_bind, afni_sif,diary_file)

        if 4 in Skip_step:
            print(bcolors.OKGREEN + 'INFO: skip step ' + str(4) + bcolors.ENDC)
        else:
            anatomical._4_skullstrip_template.skullstrip_T(BASE_SS, BASE_mask, dir_prepro, ID, Session,
                 dir_transfo, type_norm, volumes_dir, BASE_SS_coregistr, BASE_SS_mask, type_of_transform, aff_metric_ants,
                 study_template_atlas_folder, otheranat, template_skullstrip,
                 masking_img, brain_skullstrip_1, brain_skullstrip_2, masks_dir,
                 check_visualy_final_mask, s_bind, afni_sif, fsl_sif, fs_sif, itk_sif, strip_sif, overwrite,diary_file)

    ###### LOOP AGAIN ######
    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, all_Session_max):

        NL2 = '###################################################### work on subject: ' + str(ID) + ' Session ' + str(Session) + ' BLOCK 2 ###############################################################'
        print(bcolors.HEADER + NL2 + bcolors.ENDC)

        # The anatomy
        path_anat     = opj(data_path,'anat/')
        dir_transfo   = opj(path_anat,'matrices')
        dir_native    = opj(path_anat,'native')
        dir_prepro    = opj(dir_native,'01_preprocess')
        wb_native_dir = opj(dir_native,'02_Wb')
        volumes_dir   = opj(wb_native_dir,'volumes')
        masks_dir     = opj(volumes_dir,'masks')
        Ref_file = opj(volumes_dir, ID + type_norm + '_brain_step_1.nii.gz')

        date_file = datetime.date.today()
        ct = datetime.datetime.now()
        diary_name = str(ID) + ' session ' + str(Session) + str(date_file) + '.txt'
        diary_file = opj(path_anat, diary_name)
        if not opi(diary_file):
            diary = open(diary_file, "w")
            diary.write(NL2)
        else:
            diary = open(diary_file, "a")
            diary.write(NL2)
        diary.write(f'\n{ct}')
        diary.write(f'\n')
        diary.close()

        ####################################################################################
        ########################## Coregistration template to anat #########################
        ####################################################################################
        if (coregistration_longitudinal and Session != max_ses) or coregistration_longitudinal==False:
            if coregistration_longitudinal:
                data_path_max     = opj(bids_dir,'sub-' + ID,'ses-' + str(max_ses))
                path_anat_max     = opj(data_path_max,'anat')
                dir_transfo_max   = opj(path_anat_max,'matrices')
                dir_native_max    = opj(path_anat_max,'native')
                wb_native_dir_max = opj(dir_native_max,'02_Wb')
                volumes_dir_max   = opj(wb_native_dir_max,'volumes')
                masks_dir_max     = opj(volumes_dir_max,'masks')

                #####!!! use the last image as template
                BASE_SS_coregistr     = opj(volumes_dir_max,ID + type_norm + '_brain_step_1.nii.gz')
                masking_img = type_norm
                BASE_SS_mask = opj(masks_dir_max, ID + masking_img + '_mask_2.nii.gz')

                transfo_concat_inv = \
                [opj(dir_transfo_max,'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat'),
                 opj(dir_transfo_max,'template_to_' + type_norm + '_SyN_final_max_1InverseWarp.nii.gz'),
                 opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
                 opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_1InverseWarp.nii.gz')]
                w2inv_inv = [True, False, True, False]

            ################# coregistration non longitudinal #################
            else:
                if creat_study_template == True:
                    stdy_template_mask = opj(study_template_atlas_folder, 'studytemplate2_' + type_norm,
                                             'study_template_mask.nii.gz')
                    stdy_template = opj(study_template_atlas_folder, 'studytemplate2_' + type_norm,
                                        'study_template.nii.gz')
                    BASE_SS_coregistr     = stdy_template
                    BASE_SS_mask = stdy_template_mask
                else:
                    BASE_SS_coregistr     = BASE_SS
                    BASE_SS_mask = BASE_mask

                transfo_concat_inv = \
                    [opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_1InverseWarp.nii.gz')]
                w2inv_inv = [True, False]

            # create a template for each individual (nice skulltrip corrected of each mean of anat img)
            if 5 in Skip_step:
                print(bcolors.OKGREEN + 'INFO: skip step ' + str(5) + bcolors.ENDC)
            else:
                anatomical._5_create_template_brain.create_indiv_template_brain(dir_prepro, type_of_transform, ID, aff_metric_ants, Session, listTimage, volumes_dir, masking_img, brain_skullstrip_1, brain_skullstrip_2, masks_dir, type_norm, BASE_SS_coregistr, BASE_SS_mask, otheranat,
                check_visualy_final_mask, template_skullstrip, study_template_atlas_folder, bids_dir, s_bind,afni_sif,fsl_sif,fs_sif, itk_sif, strip_sif,diary_file)

            ######## co-registration to template
            if 6 in Skip_step:
                print(bcolors.OKGREEN + 'INFO: skip step ' + str(6) + bcolors.ENDC)
            else:
                anatomical._6_brainT_to_stdyT.brainT_to_T(aff_metric_ants_Transl, dir_prepro, ID, Session, listTimage, n_for_ANTS, dir_transfo, type_norm, BASE_SS_coregistr, Ref_file, volumes_dir, transfo_concat_inv,w2inv_inv,bids_dir, type_of_transform, aff_metric_ants,diary_file)

        ######## ADD max img co-registration to template
        if coregistration_longitudinal == True and Session == max_ses:
            if creat_study_template == True:
                stdy_template_mask = opj(study_template_atlas_folder, 'studytemplate2_' + type_norm,
                                         'study_template_mask.nii.gz')
                stdy_template = opj(study_template_atlas_folder, 'studytemplate2_' + type_norm,
                                    'study_template.nii.gz')
                BASE_SS_coregistr     = stdy_template
                BASE_SS_mask = stdy_template_mask
            else:
                BASE_SS_coregistr     = BASE_SS
                BASE_SS_mask = BASE_mask

            transfo_concat_inv = \
                [opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat'),
                 opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_1InverseWarp.nii.gz')]
            w2inv_inv = [True, False]

            # create a template for each individual (nice skulltrip corrected of each mean of anat img)
            if 5 in Skip_step:
                print(bcolors.OKGREEN + 'INFO: skip step ' + str(5) + bcolors.ENDC)
            else:
                anatomical._5_create_template_brain.create_indiv_template_brain(dir_prepro, type_of_transform, ID, aff_metric_ants, Session, listTimage, volumes_dir, masking_img, brain_skullstrip_1, brain_skullstrip_2, masks_dir, type_norm, BASE_SS_coregistr, BASE_SS_mask, otheranat,
                check_visualy_final_mask, template_skullstrip, study_template_atlas_folder, bids_dir, s_bind,afni_sif,fsl_sif,fs_sif, itk_sif, strip_sif,diary_file)

            ###### co-registration of each indiv template to the selected template (sty, atlas)
            if 6 in Skip_step:
                print(bcolors.OKGREEN + 'INFO: skip step ' + str(6) + bcolors.ENDC)
            else:
                anatomical._6_brainT_to_stdyT_max.brainT_to_T_max(aff_metric_ants_Transl, aff_metric_ants, creat_study_template, dir_prepro, ID, Session, listTimage, n_for_ANTS,
                    dir_transfo, type_norm, BASE_SS_coregistr, Ref_file, volumes_dir, transfo_concat_inv,w2inv_inv,
                    study_template_atlas_folder, otheranat, bids_dir, type_of_transform,
                    which_on, all_data_path_max, IgotbothT1T2, all_data_path,
                    s_bind, afni_sif,diary_file)

    if creat_study_template==True:

        ct = datetime.datetime.now()
        diary_name = 'Study_template_BLOCK2.txt'
        diary_file = opj(study_template_atlas_folder, diary_name)
        if not opi(diary_file):
            diary = open(diary_file, "w")
            diary.write('ADD the atlases')
        else:
            diary = open(diary_file, "a")
        diary.write(f'\n{ct}')
        diary.write(f'\n')
        diary.close()

        ################### redefine new atlases variable!!!
        if 7 in Skip_step:
            print(bcolors.OKGREEN + 'INFO: skip step ' + str(7) + bcolors.ENDC)
        else:
            anatomical._7_stdyT_to_AtlasT.stdyT_to_AtlasT(aff_metric_ants_Transl_template, list_atlases, Aseg_ref, Aseg_refLR, BASE_SS, dir_out, n_for_ANTS, aff_metric_ants, study_template_atlas_folder, Atemplate_to_Stemplate, type_of_transform_stdyT, overwrite,
                                                          s_bind,afni_sif,diary_file)

        ###### re-define the variable: study template atlas, etc should now be the new template !!
        list_atlases3 = list_atlases
        list_atlases = []
        for atlas_new in list_atlases3:
            atlas2 = opj(dir_out, opb(atlas_new))
            list_atlases.append(atlas2)

        BASE_SS2    = BASE_SS
        BASE_mask2  = BASE_mask
        Aseg_ref2   = Aseg_ref
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
    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, all_Session_max):
        NL3 ='###################################################### work on subject: ' + str(ID) + ' Session ' + str(Session) + ' BLOCK 3 ###############################################################'
        print(bcolors.HEADER + NL3 + bcolors.ENDC)

        animal_folder = 'sub-' + ID + '_ses-' + str(Session)

        # The anatomy
        path_anat     = opj(data_path,'anat')
        dir_transfo   = opj(path_anat,'matrices')
        dir_native    = opj(path_anat,'native')
        dir_prepro    = opj(dir_native,'01_preprocess')
        wb_native_dir = opj(dir_native,'02_Wb')
        volumes_dir   = opj(wb_native_dir,'volumes')
        labels_dir    = opj(volumes_dir,'labels')
        masks_dir     = opj(volumes_dir,'masks')

        date_file = datetime.date.today()
        ct = datetime.datetime.now()
        diary_name = str(ID) + ' session ' + str(Session) + str(date_file) + '.txt'
        diary_file = opj(path_anat, diary_name)
        if not opi(diary_file):
            diary = open(diary_file, "w")
            diary.write(NL3)
        else:
            diary = open(diary_file, "a")
            diary.write(NL3)
        diary.write(f'\n{ct}')
        diary.write(f'\n')
        diary.close()

        ####################################################################################
        ########################## Coregistration template to anat #########################
        ####################################################################################
        if coregistration_longitudinal==True:
            if creat_study_template == True:
                stdy_template_mask = opj(study_template_atlas_folder, 'studytemplate2_' + type_norm,
                                         'study_template_mask.nii.gz')
                stdy_template = opj(study_template_atlas_folder, 'studytemplate2_' + type_norm,
                                    'study_template.nii.gz')
                BASE_SS_coregistr     = stdy_template
                BASE_SS_mask = stdy_template_mask
            else:
                BASE_SS_coregistr     = BASE_SS
                BASE_SS_mask = BASE_mask
                
            if Session == max_ses:
                transfo_concat = \
                    [opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_1Warp.nii.gz'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat')]
                w2inv_fwd = [False, False]

            else:
                data_path_max = opj(bids_dir,'sub-' + ID,'ses-' + str(max_ses))
                path_anat_max    = opj(data_path_max,'anat/')
                dir_transfo_max  = opj(path_anat_max,'matrices')

                transfo_concat = \
                    [opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_1Warp.nii.gz'),
                     opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_0GenericAffine.mat'),
                     opj(dir_transfo_max,'template_to_' + type_norm + '_SyN_final_max_1Warp.nii.gz'),
                     opj(dir_transfo_max,'template_to_' + type_norm + '_SyN_final_max_0GenericAffine.mat')]
                w2inv_fwd = [False, False, False, False]

        ################# coregistration non longitudinal #################
        else:
            if creat_study_template == True:
                stdy_template_mask = opj(study_template_atlas_folder, 'studytemplate2_' + type_norm,
                                         'study_template_mask.nii.gz')
                stdy_template = opj(study_template_atlas_folder, 'studytemplate2_' + type_norm,
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

        Ref_file = opj(volumes_dir,ID + type_norm + '_brain.nii.gz')

        if 8 in Skip_step:
            print(bcolors.OKGREEN + 'INFO: skip step ' + str(8) + bcolors.ENDC)
        else:
            anatomical._8_prepar_aseg.prepar_aseg(IgotbothT1T2, Ref_file, labels_dir, volumes_dir, masks_dir, dir_transfo, BASE_SS_mask, BASE_SS_coregistr, Aseg_refLR, Aseg_ref,
                type_norm, ID, transfo_concat,w2inv_fwd, dir_prepro, list_atlases, check_visualy_each_img, n_for_ANTS, otheranat, overwrite, bids_dir, Session,
                s_bind,afni_sif,itk_sif,diary_file)
        
        if do_fMRImasks == True:
            if 9 in Skip_step:
                print(bcolors.OKGREEN + 'INFO: skip step ' + str(9) + bcolors.ENDC)
            else:
                ########################## Building fMRI masks for EPI analysis ##############################
                anatomical._9_do_fMRImasks.do_fMRImasks(masks_dir, labels_dir, type_norm, fMRImasks, overwrite,s_bind,afni_sif,diary_file)

        if 10 in Skip_step:
            print(bcolors.OKGREEN + 'INFO: skip step ' + str(10) + bcolors.ENDC)
        else:
            ########################## White Surface construction ##############################
            # You can go grab a cup of coffe, it can take more than an hour...
            anatomical._10_nii_to_mgz.nii_to_mgz(ID, Session, FS_dir, Ref_file, labels_dir, volumes_dir, otheranat, IgotbothT1T2, type_norm, overwrite,
                                                s_bind,fs_sif, diary_file)

        if 11 in Skip_step:
            print(bcolors.OKGREEN + 'INFO: skip step ' + str(11) + bcolors.ENDC)
        else:
            ########################## White Surface construction ##############################
            anatomical._11_FS_1_white.White_create(FS_dir, animal_folder,s_bind,fs_sif,diary_file)

            # You can go grab a cup of coffe, it can take more than an hour...
            anatomical._11_FS_1_white.White_more(FS_dir, animal_folder, FS_buckner40_TIF,FS_buckner40_GCS,s_bind,fs_sif,diary_file)

        if 12 in Skip_step:
            print(bcolors.OKGREEN + 'INFO: skip step ' + str(12) + bcolors.ENDC)
        else:
            anatomical._12_make_pial.make_pial(FS_dir, animal_folder, type_norm, otheranat, Hmin, Ref_file, do_surfacewith, overwrite,
                                               s_bind,fs_sif,diary_file)

        if check_visualy_each_img == True:
            if 13 in Skip_step:
                print(bcolors.OKGREEN + 'INFO: skip step ' + str(13) + bcolors.ENDC)
            else:
                anatomical._13_FS_freeview.FS_Freeview(FS_dir, animal_folder, 'pial', Lut_file,
                                                       s_bind,fs_sif)

        if 14 in Skip_step:
            print(bcolors.OKGREEN + 'INFO: skip step ' + str(14) + bcolors.ENDC)
        else:
            anatomical._14_Finalise.FS_finalise(FS_dir, animal_folder, FreeSlabel_ctab_list, list_atlases_2, labels_dir, type_norm, Ref_file,
                                                s_bind,fs_sif,diary_file)
        if 15 in Skip_step:
            print(bcolors.OKGREEN + 'INFO: skip step ' + str(15) + bcolors.ENDC)
        else:
            anatomical._15_to_WB.WB_prep(FS_dir, dir_native, animal_folder, Ref_file, species, list_atlases_2,s_bind,afni_sif,fsl_sif,fs_sif,wb_sif,diary_file)

        if 16 in Skip_step:
            print(bcolors.OKGREEN + 'INFO: skip step ' + str(16) + bcolors.ENDC)
        else:
            anatomical._16_anat_QC_SNR.anat_QC(BASE_SS_coregistr, Ref_file, type_norm, labels_dir, dir_prepro, ID, listTimage, masks_dir,
            s_bind, afni_sif, diary_file)

        if 100 in Skip_step:
                print(bcolors.OKGREEN + 'INFO: skip step ' + str(100) + bcolors.ENDC)
        else:
            anatomical._100_Data_Clean.clean(all_ID, all_Session, all_data_path,diary_file)

        if 200 in Skip_step:
            print(bcolors.OKGREEN + 'INFO: skip step ' + str(200) + bcolors.ENDC)

        else:
            anatomical._200_Data_QC._itk_check_masks(dir_prepro, masks_dir, ID, type_norm,s_bind, itk_sif,diary_file)