import os
opj = os.path.join
from Tools import diaryfile, getpath, run_cmd
from anat.loop2 import _mask2_Data_QC, _5_create_template_brain, _6_stdyTmax


def run(ID, Session, data_path, max_ses,type_norm,ref_suffix,coregistration_longitudinal,bids_dir, mask_suffix,
        creat_study_template,stdy_template,stdy_template_mask,BASE_SS,BASE_mask,Skip_step,listTimage,list_transfo,
        brain_skullstrip_2,check_visualy_final_mask, MNIBcorrect_indiv,sing_afni, sing_fsl, sing_fs, sing_itk,
        sing_synstrip,sing_wb,preftool,reference,masking_img,Unetpath):

    nl = '####################### work on subject: ' + str(ID) + ' Session ' + str(Session) + ' BLOCK 2 ###############################'
    run_cmd.printcolor(nl, 'HEADER')

    # The anatomy
    path_anat, _,_, dir_prepro, _, volumes_dir, _, masks_dir =getpath.anat(data_path, reference,
                                                                     '', False,
                                                                     True, 'native')


    diary_file = diaryfile.create(opj(path_anat, str(ID) + ' session ' + str(Session)), nl)

    ####################################################################################
    ########################## Coregistration template to anat #########################
    ####################################################################################


    if coregistration_longitudinal == False:

        if creat_study_template == True:
            BASE_SS_coregistr = stdy_template
            BASE_SS_mask      = stdy_template_mask
            sourcename        = 'studyTemplate'
        else:
            BASE_SS_coregistr = BASE_SS
            BASE_SS_mask      = BASE_mask
            sourcename        = 'EDNiX'

    else:
        if Session != max_ses:
            data_path_max = opj(bids_dir, 'sub-' + ID, 'ses-' + str(max_ses))
            path_anat_max, dir_transfo_max, _, _, dir_native_max, volumes_dir_max, _,masks_dir_max = getpath.anat(data_path_max, reference,
                                                                     '', False,
                                                                     True, 'native')
            #####!!! use the last image as template
            BASE_SS_coregistr = opj(volumes_dir_max, ID + ref_suffix  + type_norm + '.nii.gz')
            BASE_SS_mask      = opj(masks_dir_max,   ID + mask_suffix + '_mask.nii.gz')
            sourcename        = 'refSession'

            masking_img = type_norm

        else : #Session == max_ses
            if creat_study_template == True:
                BASE_SS_coregistr = stdy_template
                BASE_SS_mask      = stdy_template_mask
                sourcename        = 'studyTemplate'
            else:
                BASE_SS_coregistr = BASE_SS
                BASE_SS_mask      = BASE_mask
                sourcename        = 'EDNiX'


    # create a template for each individual (nice skulltrip corrected of each mean of anat img)
    if 5 in Skip_step:
        run_cmd.msg('INFO: skip step ' + str(5), diary_file, 'OKGREEN')
    else:

        _5_create_template_brain.create_indiv_template_brain(dir_prepro, list_transfo, ID,
                                                             Session, listTimage, volumes_dir, masking_img,
                                                             brain_skullstrip_2, masks_dir, BASE_SS_coregistr,
                                                             BASE_SS_mask,check_visualy_final_mask,
                                                             bids_dir, MNIBcorrect_indiv,sing_afni, sing_fsl, sing_fs, sing_itk, sing_synstrip,Unetpath,
                                                             Skip_step, diary_file, preftool)

    if 'itk_1' in Skip_step:
        run_cmd.msg('INFO: skip step ' + str('itk_1'), diary_file, 'OKGREEN')
    else:
        output4mask = opj(masks_dir, ID + '_desc-step2_mask.nii.gz')
        end_maskname = '_'.join([ID, 'final', 'mask', '2.nii.gz'])
        input4msk = opj(volumes_dir, ID + '_space-acpc_desc-SS-step1_' + masking_img + '.nii.gz')
        _mask2_Data_QC._itk_check_masks(output4mask, input4msk, end_maskname, masks_dir, sing_itk, diary_file)

    print(sourcename)
    ###### co-registration of each indiv template to the selected template (sty, atlas)
    if 6 in Skip_step:
        run_cmd.msg('INFO: skip step ' + str(6), diary_file, 'OKGREEN')
    else:
        _6_stdyTmax.nativetoTemplate(sourcename,ID, listTimage, data_path,
        bids_dir, creat_study_template, coregistration_longitudinal, reference, ref_suffix, BASE_SS_coregistr,
        list_transfo, type_norm,sing_wb, diary_file)













