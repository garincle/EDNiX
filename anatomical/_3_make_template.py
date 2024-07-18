##############################################################################
####### CREATE THE STUDY TEMPLATE (IF YOU WANT ON) ###########################
##############################################################################
import os
import subprocess
import shutil
import anatomical.anat_to_common_EMB
from nilearn.masking import compute_epi_mask

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput


############# use #############
#######Name: nilearn
#######Version: 0.9.0
################################


def make_template(which_on, all_ID_max, max_session, all_data_path_max, all_ID, all_Session, all_data_path, type_norm, study_template_atlas_forlder, template_skullstrip, BASE_SS, BASE_mask, overwrite,
                  s_bind, afni_sif, fsl_sif, fs_sif):

    if which_on == 'max': # all or max
        all_ID_temp = all_ID_max
        all_Session_temp = max_session
        all_data_path_temp = all_data_path_max

    elif which_on == 'all':
        all_ID_temp = all_ID
        all_Session_temp = all_Session
        all_data_path_temp = all_data_path

    else:
        print("error with which_on name, need to be all or max")

    template_list = []
    for ID, Session, data_path in zip(all_ID_temp, all_Session_temp, all_data_path_temp):

        # The anatomy
        path_anat    = opj(data_path,'anat/')
        dir_transfo  = opj(path_anat,'matrices')

        dir_native    = opj(path_anat,'native')
        dir_prepro    = opj(dir_native,'01_preprocess')
        wb_native_dir = opj(dir_native,'02_Wb')
        volumes_dir   = opj(wb_native_dir,'volumes')
        labels_dir    = opj(volumes_dir,'labels')
        masks_dir     = opj(volumes_dir,'masks')

        ############load the image for the template
        #creat a large image for co-registration
        if not ope(opj(opj(all_data_path[0], 'anat', 'native', '01_preprocess'), all_ID[0] + '_acpc_cropped' + type_norm + 'Zp.nii.gz')):
            command = 'singularity run' + s_bind + afni_sif + '3dZeropad -overwrite -I 20 -S 20 -A 20 -P 20 -L 20 -R 20 -S 20 -prefix ' + opj(opj(all_data_path[0], 'anat', 'native', '01_preprocess'), all_ID[0] + '_acpc_cropped' + type_norm + 'Zp.nii.gz') + \
                      ' ' + opj(opj(all_data_path[0], 'anat', 'native', '01_preprocess'), all_ID[0] + '_acpc_cropped' + type_norm + '.nii.gz')
            spco(command, shell=True)

        #resemple the other anat to this large img
        command = 'singularity run' + s_bind + afni_sif + '3dresample -overwrite -master ' + opj(opj(all_data_path[0], 'anat', 'native', '01_preprocess'), all_ID[0] + '_acpc_cropped' + type_norm + 'Zp.nii.gz') + ' -overwrite' + \
        ' -input ' + opj(dir_prepro, ID + '_acpc_cropped' + type_norm + '.nii.gz') + ' -prefix ' + opj(dir_prepro, ID + '_acpc_cropped' + type_norm + 'Zp.nii.gz') + ' -bound_type SLAB'
        spco(command, shell=True)
        template_list.append(opj(dir_prepro, ID + '_acpc_cropped' + type_norm + 'Zp.nii.gz'))

            ###########################template

    templatedir2 = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm)
    if not os.path.exists(templatedir2): os.mkdir(templatedir2)
    anatomical.anat_to_common_EMB.anats_to_common(
        template_list,
        templatedir2,
        brain_volume=300000,
        blur_radius_coarse=11,
        convergence=0,
        registration_kind='nonlinear',
        nonlinear_levels=[4,7],
        nonlinear_minimal_patches=[13,9])

    warp_adj = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'warped_3_adjusted_mean.nii.gz')
    #shutil.copyfile(warp_adj, opj(study_template_atlas_forlder, 'studytemplate_' + type_norm, 'study_template.nii.gz'))
    stdy_template_mask = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template_mask.nii.gz')
    stdy_template = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template.nii.gz')

    if template_skullstrip == 'Custum_Macaque':
        #convert to float
        command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + warp_adj + ' ' + opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template_float.nii.gz')
        spco([command], shell=True)
        mask_img = compute_epi_mask(opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template_float.nii.gz'), lower_cutoff=0.75, upper_cutoff=0.85, connected=True, opening=3,
            exclude_zeros=False, ensure_finite=True)
        mask_img.to_filename(stdy_template_mask)
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + stdy_template_mask + \
        ' -input ' + stdy_template_mask + ' -fill_holes'
        spco(command, shell=True)


    #elif template_skullstrip == 'MachinL':
    #    spco(['python3', '/home/cgarin/Documents/0_Clement/CODE/CODE/NHP-BrainExtraction-master/UNet_Model/muSkullStrip.py',
    #    '-in', warp_adj, '-model', '/home/cgarin/Documents/1_Macaque_MRI/4_SSwarper_muSkull/result/model-20-epoch', '-out',
    #     opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm)])
    #    shutil.copyfile(opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'warped_3_adjusted_mean_pre_mask.nii.gz'), stdy_template_mask)

    elif template_skullstrip == 'Custum_Macaque2':
        #convert to float
        mask_img = compute_epi_mask(warp_adj, lower_cutoff=0.75, upper_cutoff=0.90, connected=True, opening=3,
            exclude_zeros=False, ensure_finite=True)
        mask_img.to_filename(stdy_template_mask)

    elif template_skullstrip == 'Manual':
        while True:
            Your_choice = input('Enter "OK" if you have finished doing/corrected yourself the template mask img to continu and save it as "study_template_mask" in the template folder')
            if Your_choice == '' or not Your_choice in ['C','OK']: 
                print('Please answer with OK!') 
            else: 
                break 
        if Your_choice=="OK":
            print("CONTINUE!!")

    elif template_skullstrip =='Custum_QWARPT2':

        command = 'singularity run' + s_bind + afni_sif + '3dQwarp -overwrite -iwarp' + \
        ' -base ' + BASE_SS + \
        ' -prefix ' + stdy_template_mask + \
        ' -source ' + warp_adj + ' -maxlev 3 -resample'
        spco(command, shell=True)

        command = 'singularity run' + s_bind + afni_sif + '3dNwarpApply -nwarp ' + opj(study_template_atlas_forlder,'study_template_mask_WARPINV.nii.gz') + \
        ' -source ' + BASE_mask + ' -master ' + warp_adj + ' -interp NN' + \
        ' -prefix ' + stdy_template_mask + ' -overwrite'
        spco(command, shell=True)


    elif template_skullstrip == '3dSkullStrip':
        command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -overwrite -prefix ' + stdy_template_mask + \
            ' -input ' + warp_adj + ' -blur_fwhm 2 -orig_vol -mask_vol -use_skull -monkey'
        spco(command, shell=True)


    elif template_skullstrip =='bet2':
        #####creat an approximate brain mask
        command = 'singularity run' + s_bind + fsl_sif + 'bet2 ' + warp_adj + ' ' + stdy_template_mask + \
        ' -f 0.70'
        spco([command], shell=True)

        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' +stdy_template_mask + \
        ' -input ' + stdy_template_mask + ' -fill_holes'
        spco(command, shell=True)


    #elif brain_skullstrip == 'antsBrainExtraction':
    #    spco(['antsBrainExtraction.sh', '-d', '3',
    #        '-a', warp_adj,
    #        '-e', BASE_SS, '-m', BASE_mask,
    #        '-o', study_template_atlas_forlder + '/studytemplate2_' + type_norm + '/'])
    #    shutil.copyfile(opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'BrainExtractionMask.nii.gz'), stdy_template_mask)


    '''
    'singularity run' + s_bind + afni_sif + 3dresample -master /media/cgarin/Clement_1/1_Macaques/1_PFC_study/10_2023_PP_U/Study_template/studytemplate2_T1/warped_3_adjusted_mean.nii.gz -input /media/cgarin/Clement_1/1_Macaques/1_PFC_study/10_2023_PP_U/Study_template/studytemplate2_T1/study_template_maskcp.nii.gz -prefix /media/cgarin/Clement_1/1_Macaques/1_PFC_study/10_2023_PP_U/Study_template/studytemplate2_T1/study_template_mask.nii.gz
    '''

    ##extract brain
    command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + stdy_template_mask + ' -b ' + warp_adj + ' -expr "(a*b)" -prefix ' + opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template_not_align.nii.gz')
    spco(command, shell=True)

    ### to move in brainT to styT??? XXX
    ##align template to BASE_SS (atlas template)
    command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -warp shift_rotate -cmass -overwrite -base ' + BASE_SS + \
    ' -nomask -onepass' + \
    ' -master ' + opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template_not_align.nii.gz') + \
    ' -prefix ' + stdy_template + \
    ' -source ' + opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template_not_align.nii.gz') + ' -1Dmatrix_save ' + \
    opj(study_template_atlas_forlder,ID + '_allign_template_to_stdy_template.1D')
    spco(command, shell=True)

    command = 'singularity run' + s_bind + afni_sif + '3dAllineate -overwrite -interp NN -1Dmatrix_apply ' + opj(study_template_atlas_forlder,ID + '_allign_template_to_stdy_template.1D') + \
    ' -prefix ' + stdy_template_mask + \
    ' -master ' + stdy_template + \
    ' -input  ' + stdy_template_mask
    spco([command], shell=True)  
