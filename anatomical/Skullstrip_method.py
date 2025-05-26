###################################
###      Skullstrip method      ###
###################################
import subprocess
import os
import shutil
import nibabel as nib
import ants
import math
import nilearn
from nilearn.masking import compute_epi_mask
import numpy as np
from anatomical import Histrogram_mask_EMB
from fonctions.extract_filename import extract_filename
import datetime
import json

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

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput


def Skullstrip_method(step_skullstrip, template_skullstrip, study_template_atlas_folder, masking_img, brain_skullstrip_1, brain_skullstrip_2, masks_dir, volumes_dir, dir_prepro, type_norm, BASE_SS_coregistr, BASE_SS_mask,
    type_of_transform, ID, aff_metric_ants, check_visualy_final_mask, s_bind, afni_sif, fsl_sif, fs_sif, itk_sif, strip_sif,diary_file):
    ct = datetime.datetime.now()
    nl = 'Run anatomical.Skullstrip_method.Skullstrip_method'
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{nl}')

    nl = 'INFO: If you can not find a good solution for Skullstriping due to bad image quality, you can always modify it by hands and save it as: ' + \
         opj(masks_dir, ID + masking_img + 'final_mask.nii.gz') + ' for step 1' + \
         opj(masks_dir, ID + masking_img + 'final_mask_2.nii.gz') + ' for step 2' + \
         opj(study_template_atlas_folder, 'studytemplate2_' + type_norm, 'studytemplate_final_mask.nii.gz') + ' for step "sty template"'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    nl = 'WARNING: Note that any manual segmentation (saved with the correct name) it will AUTOMATICALLY be selected for Skullstriping'
    print(bcolors.WARNING + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    nl1=nl2=''
    if step_skullstrip == 1:
        masking_img = masking_img
        input_for_msk = opj(dir_prepro, ID + '_anat_reorient_NU' + masking_img + '.nii.gz')
        output_for_mask = opj(masks_dir, ID + masking_img + '_mask_1.nii.gz')
        brain_skullstrip = brain_skullstrip_1
        process_dir = masks_dir

        nl1 = "INFO: brain_skullstrip method is " + brain_skullstrip
        nl2 = 'INFO: looking for manual segmentation named:' + opj(masks_dir, ID + masking_img + 'final_mask.nii.gz') + '...'


    elif step_skullstrip == 2:
        masking_img = masking_img
        input_for_msk = opj(volumes_dir, ID + '_' + masking_img + '_template.nii.gz')
        output_for_mask = opj(masks_dir, ID + masking_img + '_mask_2.nii.gz')
        brain_skullstrip = brain_skullstrip_2
        process_dir = masks_dir

        nl1 = "INFO: brain_skullstrip method is " + brain_skullstrip
        nl2 = 'INFO: looking for manual segmentation named:' + opj(masks_dir, ID + masking_img + 'final_mask_2.nii.gz') + '...'

    elif step_skullstrip == 3:
        masking_img = type_norm
        input_for_msk = opj(study_template_atlas_folder, 'studytemplate2_' + type_norm, 'warped_3_adjusted_mean.nii.gz')
        output_for_mask = opj(study_template_atlas_folder, 'studytemplate2_' + type_norm, 'study_template_mask.nii.gz')
        brain_skullstrip = template_skullstrip
        process_dir = study_template_atlas_folder + '/studytemplate2_' + type_norm + '/'
        nl1 = "INFO: brain_skullstrip method is " + brain_skullstrip
        nl2 = 'INFO: looking for manual segmentation named:' + opj(study_template_atlas_folder, 'studytemplate2_' + type_norm, 'studytemplate_final_mask.nii.gz') + '...'

    else:
        nl1 ='No step_skullstrip ?'
        diary.write(f'\n{nl1}')
        nl2 = 'Nothing will be done.'
        print(nl2)
        diary.write(f'\n{nl2}')
        raise Exception(bcolors.FAIL + nl1 + ': '+ nl2 + bcolors.ENDC)

    if not nl1 == '':
        print(bcolors.OKGREEN + nl1 + bcolors.ENDC)
        diary.write(f'\n{nl1}')
        print(bcolors.OKGREEN + nl2 + bcolors.ENDC)
        diary.write(f'\n{nl2}')

    ###########################################################################################################################################################################
    ######################################################### Brain Skullstrip Library (do your own if you need to!!!) ########################################################
    ###########################################################################################################################################################################

    # This Brain Skullstrip Library will be evolving and can/should be personalized, don't hesitate to propose other efficient solution on GIT!!!
    nl1 = nl2 =''
    if os.path.exists(opj(masks_dir, ID + masking_img + 'final_mask.nii.gz')) and step_skullstrip == 1:
        nl1 = 'WARNING: We found a final mask for Skullstrip 1!!! no Skullstrip will be calculated!'
        nl2 = 'INFO: please delete' + opj(masks_dir, ID + masking_img + 'final_mask.nii.gz') + ' if you want retry to create a skulstripp images'
        shutil.copyfile(opj(masks_dir, ID + masking_img + 'final_mask.nii.gz'), output_for_mask)
        dictionary = {"Sources": opj(masks_dir, ID + masking_img + 'final_mask.nii.gz'),
                      "Description": 'Copy.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)
        print(bcolors.WARNING + nl1 + bcolors.ENDC)
        diary.write(f'\n{nl1}')
        print(bcolors.OKGREEN + nl2 + bcolors.ENDC)
        diary.write(f'\n{nl2}')

    elif os.path.exists(opj(masks_dir, ID + masking_img + 'final_mask_2.nii.gz')) and step_skullstrip == 2:
        nl1 = 'WARNING: We found a final mask for Skullstrip 2!!! no Skullstrip will be calculated!'
        nl2 = 'INFO: please delete' + opj(masks_dir, ID + masking_img + 'final_mask_2.nii.gz') + ' if you want retry to create a skulstripp images'
        shutil.copyfile(opj(masks_dir, ID + masking_img + 'final_mask_2.nii.gz'), output_for_mask)
        dictionary = {"Sources": opj(masks_dir, ID + masking_img + 'final_mask_2.nii.gz'),
                      "Description": 'Copy.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

        print(bcolors.WARNING + nl1 + bcolors.ENDC)
        diary.write(f'\n{nl1}')
        print(bcolors.OKGREEN + nl2 + bcolors.ENDC)
        diary.write(f'\n{nl2}')

    elif os.path.exists(opj(study_template_atlas_folder, 'studytemplate2_' + type_norm, 'studytemplate_final_mask.nii.gz')) and step_skullstrip == 3:
        nl1 = 'WARNING: We found a final mask for Skullstriping the study template, no Skullstrip will be calculated!'
        nl2 = 'INFO: please delete' + opj(study_template_atlas_folder, 'studytemplate2_' + type_norm, 'studytemplate_final_mask.nii.gz') + ' if you want retry to create a skulstripp images'
        shutil.copyfile(opj(study_template_atlas_folder, 'studytemplate2_' + type_norm, 'studytemplate_final_mask.nii.gz'), output_for_mask)
        dictionary = {"Sources": opj(study_template_atlas_folder, 'studytemplate2_' + type_norm, 'studytemplate_final_mask.nii.gz'),
                      "Description": 'Copy.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

        print(bcolors.WARNING + nl1 + bcolors.ENDC)
        diary.write(f'\n{nl1}')
        print(bcolors.OKGREEN + nl2 + bcolors.ENDC)
        diary.write(f'\n{nl2}')

    else:
        nl1 = 'WARNING: Manual segmentation not found'
        nl2 = 'INFO: continuing with the selected skullstriping method'

        print(bcolors.WARNING + nl1 + bcolors.ENDC)
        diary.write(f'\n{nl1}')
        print(bcolors.OKGREEN + nl2 + bcolors.ENDC)
        diary.write(f'\n{nl2}')

        #### species specific ######################################## MACAQUE ####################################

        if brain_skullstrip == '3dSkullStrip_macaque':
            command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
                '-input ' + input_for_msk + ' -blur_fwhm 1 -orig_vol -mask_vol -monkey'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes -dilate_input 4'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (3dSkullStrip from AFNI).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)


        #### species specific ######################################## MARMOSET ####################################

        elif brain_skullstrip == '3dSkullStrip_marmoset':
            command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
                '-input ' + input_for_msk + ' -blur_fwhm 1 -shrink_fac 0.1 -fac 100 -norm_vol -mask_vol -marmoset'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (3dSkullStrip from AFNI).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip == '3dSkullStrip_dog_macaque':
            command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
                '-input ' + input_for_msk + ' -blur_fwhm 2 -orig_vol -mask_vol -monkey'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes -dilate_input 4 -4'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (3dSkullStrip from AFNI).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        #################################### BABOON ####################################

        elif brain_skullstrip == 'Custum_Baboon':
            #convert to float
            command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            mask_img = compute_epi_mask(input_for_msk[:-7] + '_float.nii.gz', lower_cutoff=0.85, upper_cutoff=0.95, connected=True, opening=3,
                exclude_zeros=False, ensure_finite=True)
            mask_img.to_filename(output_for_mask)
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes -dilate_input 15'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (compute_epi_mask from Nilearn).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        #################################### MACAQUE ####################################

        elif brain_skullstrip == 'Custum_Macaque':
            #convert to float
            command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            mask_img = compute_epi_mask(input_for_msk[:-7] + '_float.nii.gz', lower_cutoff=0.2, upper_cutoff=0.90, connected=True, opening=3,
                exclude_zeros=False, ensure_finite=True)
            mask_img.to_filename(output_for_mask)
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes -dilate_input 8 -2'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (compute_epi_mask from Nilearn).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        #################################### CHIMPANZEE ####################################

        elif brain_skullstrip == 'Custum_Chimp':
            loadimg = nib.load(input_for_msk).get_fdata()
            loadimgsort85 =  np.percentile(np.abs(loadimg)[np.abs(loadimg)>0], 40)
            mask_imag = nilearn.image.threshold_img(input_for_msk, loadimgsort85, cluster_threshold=10)
            mask_imag.to_filename(output_for_mask)
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes -dilate_input 3'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            command = 'singularity run' + s_bind + afni_sif + '3dresample -master ' + input_for_msk + ' -prefix ' + output_for_mask + ' -input ' + output_for_mask + ' -overwrite -bound_type SLAB'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (image threshold from Nilearn).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        #################################### Mouse ####################################

        elif brain_skullstrip == 'Custum_mouse':
            #convert to float
            command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            mask_img = compute_epi_mask(input_for_msk[:-7] + '_float.nii.gz', lower_cutoff=0.2, upper_cutoff=0.60, connected=True, opening=3,
                exclude_zeros=False, ensure_finite=True)
            mask_img.to_filename(output_for_mask)
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes -dilate_input 3 -1'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (compute_epi_mask from Nilearn).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)


        #################################### BAT ####################################

        elif brain_skullstrip == 'Custum_bat':
            #convert to float
            command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            mask_img = compute_epi_mask(input_for_msk[:-7] + '_float.nii.gz', lower_cutoff=0.7, upper_cutoff=0.80, connected=True, opening=3,
                exclude_zeros=False, ensure_finite=True)
            mask_img.to_filename(output_for_mask)
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes -dilate_input 8 -2'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (compute_epi_mask from Nilearn).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        #################################### DOG ####################################

        elif brain_skullstrip == 'Custum_dog':
            #convert to float
            command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            mask_img = compute_epi_mask(input_for_msk[:-7] + '_float.nii.gz', lower_cutoff=0.85, upper_cutoff=0.95, connected=True, opening=3,
                exclude_zeros=False, ensure_finite=True)
            mask_img.to_filename(output_for_mask)
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes -dilate_input 23'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (compute_epi_mask from Nilearn).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        #################################### RAT ####################################
        ### specific for GD study
        elif brain_skullstrip == '3dSkullStrip_Rat':
            if ID in ['301502', '302101', '302105', '302106', '301603', '300908', '301500', '301501', '301503',
                      '301504', '301505', '301508', '301509', '302107', '302108', '300600']:
                command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite -input ' + input_for_msk + ' -orig_vol -mask_vol -rat'
                nl = spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
                command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
                          ' -input ' + output_for_mask + ' -fill_holes -dilate_input 15'
                nl= spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
            else:
                command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite -input ' + input_for_msk + ' -orig_vol -mask_vol -surface_coil -rat'
                nl = spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
                command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
                          ' -input ' + output_for_mask + ' -fill_holes -dilate_input 4'
                nl= spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (3dSkullStrip from AFNI).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip == '3dSkullStrip_rat':
            command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
                '-input ' + input_for_msk + ' -blur_fwhm 1 -orig_vol -mask_vol -rat'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes -dilate_input 2'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (3dSkullStrip from AFNI).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip == '3dSkullStrip_rat_no_dil':
            command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
                '-input ' + input_for_msk + ' -blur_fwhm 1 -orig_vol -mask_vol -rat'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes -dilate_input 2'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (3dSkullStrip from AFNI).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip == 'sammba_rat':
            command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            nichols_masker = Histrogram_mask_EMB.HistogramMask()
            nichols_masker.inputs.in_file = input_for_msk[:-7] + '_float.nii.gz'
            nichols_masker.inputs.volume_threshold = 2500
            # nichols_masker.inputs.upper_cutoff = 0.2
            # nichols_masker.inputs.lower_cutoff = 0.8
            # nichols_masker.inputs.intensity_threshold = 500
            # nichols_masker.inputs.opening = 2
            # nichols_masker.inputs.closing = 10
            nichols_masker.inputs.dilation_size = (1, 2, 3)
            nichols_masker.inputs.connected = True
            nichols_masker.inputs.out_file = output_for_mask
            res = nichols_masker.run()  # doctest: +SKIP
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (Histrogram_mask_EMB from sammba).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip == 'sammba_dog':
            command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            nichols_masker = Histrogram_mask_EMB.HistogramMask()
            nichols_masker.inputs.in_file = input_for_msk[:-7] + '_float.nii.gz'
            #nichols_masker.inputs.volume_threshold = 2500
            nichols_masker.inputs.upper_cutoff = 0.85
            nichols_masker.inputs.lower_cutoff = 0.7
            # nichols_masker.inputs.intensity_threshold = 500
            nichols_masker.inputs.opening = 2
            nichols_masker.inputs.closing = 10
            nichols_masker.inputs.dilation_size = (1, 2, 3)
            nichols_masker.inputs.connected = True
            nichols_masker.inputs.out_file = output_for_mask
            res = nichols_masker.run()  # doctest: +SKIP
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (Histrogram_mask_EMB from sammba).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
                      ' -input ' + output_for_mask + ' -fill_holes -dilate_input -5 7'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

        elif brain_skullstrip == 'custum_rat':
            # convert to float
            command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + input_for_msk + ' ' + opj(input_for_msk[:-7] + '_float.nii.gz')
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            mask_img = compute_epi_mask(input_for_msk[:-7] + '_float.nii.gz', lower_cutoff=0.8, upper_cutoff=0.85,
                                        connected=True, opening=3,
                                        exclude_zeros=True, ensure_finite=True)
            mask_img.to_filename(output_for_mask)
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
                      ' -input ' + output_for_mask + ' -fill_holes -dilate_input 3'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (compute_epi_mask from Nilearn).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip == 'sammba_mouse':
            command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
            spco([command], shell=True)
            nichols_masker = Histrogram_mask_EMB.HistogramMask()
            nichols_masker.inputs.in_file = input_for_msk[:-7] + '_float.nii.gz'
            nichols_masker.inputs.volume_threshold = 360
            # nichols_masker.inputs.upper_cutoff = 0.2
            # nichols_masker.inputs.lower_cutoff = 0.8
            # nichols_masker.inputs.intensity_threshold = 500
            # nichols_masker.inputs.opening = 2
            # nichols_masker.inputs.closing = 10
            nichols_masker.inputs.dilation_size = (1, 2, 3)
            nichols_masker.inputs.connected = True
            nichols_masker.inputs.out_file = output_for_mask
            res = nichols_masker.run()  # doctest: +SKIP
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (Histrogram_mask_EMB from sammba).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        #################################### General function ####################################
        elif brain_skullstrip == 'bet2_ANTS':
            #####creat an approximate brain mask

            hd_IMG = ants.image_header_info(input_for_msk)
            command = 'singularity run' + s_bind + fsl_sif + 'bet2 ' + input_for_msk + \
                      ' ' + opj(opj(process_dir, extract_filename(input_for_msk)) + '_bet' + masking_img) + \
                      ' -f 0.40 -c ' + str(math.ceil(int(hd_IMG['dimensions'][0]) / 2)) + \
                      ' ' + str(math.ceil(int(hd_IMG['dimensions'][1]) / 2)) + \
                      ' ' + str(math.ceil(int(hd_IMG['dimensions'][2]) / 2)) + ' -m'
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            IMG = ants.image_read(input_for_msk)
            REF_BET = ants.image_read(BASE_SS_coregistr)
            tmp_bet = ants.image_read(opj(opj(process_dir, extract_filename(input_for_msk)) + '_bet' + masking_img + '.nii.gz'))
            REF_MASK = ants.image_read(BASE_SS_mask)

            mtx1 = ants.registration(fixed=tmp_bet, moving=REF_BET, type_of_transform='Translation',
                                     outprefix=opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_SyN_shift_'))
            MEAN_tr = ants.apply_transforms(fixed=tmp_bet, moving=REF_BET, transformlist=mtx1['fwdtransforms'],
                                            interpolator='nearestNeighbor')
            ants.image_write(MEAN_tr, opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_SyN_shift.nii.gz'),
                             ri=False)
            mTx = ants.registration(fixed=tmp_bet, moving=REF_BET,
                                    outprefix=opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_SyN_'),
                                    type_of_transform=type_of_transform,
                                    initial_transform=mtx1['fwdtransforms'],
                                    aff_metric=aff_metric_ants,
                                    grad_step=0.1,
                                    flow_sigma=3,
                                    total_sigma=0,
                                    aff_sampling=32,
                                    aff_random_sampling_rate=0.2,
                                    syn_sampling=32,
                                    aff_iterations=(1000, 500, 250, 100),
                                    aff_shrink_factors=(8, 4, 2, 1),
                                    aff_smoothing_sigmas=(3, 2, 1, 0),
                                    reg_iterations=(1000, 500, 250, 100),
                                    reg_smoothing_sigmas=(3, 2, 1, 0),
                                    reg_shrink_factors=(8, 4, 2, 1),
                                    verbose=True)

            transfo_concat = [opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_SyN_1Warp.nii.gz'),
                              opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_SyN_0GenericAffine.mat')]

            tmp_mask1 = ants.apply_transforms(fixed=tmp_bet, moving=REF_MASK,
                                              transformlist=transfo_concat, interpolator='nearestNeighbor')
            tmp_mask1 = ants.threshold_image(tmp_mask1, 0.5, 1, 1, 0, True)
            tmp_mask1 = ants.morphology(tmp_mask1, operation='dilate', radius=2, mtype='binary', shape='ball')
            tmp_mask1 = ants.iMath(tmp_mask1, operation='GetLargestComponent')
            seg_tmp = ants.atropos(a=IMG, m='[0.1,1x1x1]', c='[3,0]', i='kmeans[3]', x=tmp_mask1)
            #  Clean up:
            seg_tmp2 = ants.iMath(seg_tmp['segmentation'], 'Pad', 10)
            W = ants.threshold_image(seg_tmp2, 3, 3, 1, 0, True)
            W = ants.iMath(W, operation='GetLargestComponent')
            W = W * 3
            G = ants.threshold_image(seg_tmp2, 2, 2, 1, 0, True)
            G = ants.iMath(G, operation='GetLargestComponent')
            TMP1 = ants.iMath(G, 'FillHoles', 2)
            G = G * TMP1
            C = ants.threshold_image(seg_tmp2, 1, 1, 1, 0, True)
            TMP2 = ants.morphology(C, operation='erode', radius=10, mtype='binary', shape='ball')
            G[G == 0] = TMP2[G == 0]
            G = G * 2
            seg_tmp2 = W
            seg_tmp2[W == 0] = G[W == 0]
            #  clean the Brainmask
            tmp_mask3 = ants.threshold_image(seg_tmp2, 3, 3, 1, 0, True)
            TMP3 = ants.threshold_image(seg_tmp2, 2, 2, 1, 0, True)
            tmp_mask3[tmp_mask3 == 0] = TMP3[tmp_mask3 == 0]
            tmp_mask3 = ants.morphology(tmp_mask3, operation='erode', radius=2, mtype='binary', shape='ball')
            tmp_mask3 = ants.iMath(tmp_mask3, operation='GetLargestComponent')
            tmp_mask3 = ants.morphology(tmp_mask3, operation='dilate', radius=4, mtype='binary', shape='ball')
            tmp_mask3 = ants.iMath(tmp_mask3, 'FillHoles', 2)
            tmp_mask3 = ants.morphology(tmp_mask3, operation='dilate', radius=5, mtype='binary', shape='ball')
            tmp_mask3 = ants.morphology(tmp_mask3, operation='erode', radius=5, mtype='binary', shape='ball')
            tmp_mask3 = ants.iMath(tmp_mask3, 'Pad', -10)

            ants.image_write(tmp_mask3, output_for_mask, ri=False)
            dictionary = {"Sources": [input_for_msk,
                                      BASE_SS_coregistr,
                                      BASE_SS_mask],
                          "Description": 'Brain mask (bet2 from FSL and atropos from ANTspy).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip == 'Custum_ANTS_NL':
            IMG = ants.image_read(input_for_msk)
            REF_IMG = ants.image_read(BASE_SS_coregistr)
            mtx1 = ants.registration(fixed=IMG, moving=REF_IMG, type_of_transform='Translation',
                                     outprefix=opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_SyN_'))
            REF_MASK = ants.image_read(BASE_SS_mask)
            tmp_mask1 = ants.apply_transforms(fixed=IMG, moving=REF_MASK,
                                              transformlist=mtx1['fwdtransforms'], interpolator='nearestNeighbor')
            tmp_mask1 = ants.threshold_image(tmp_mask1, 0.5, 1, 1, 0, True)
            tmp_mask1 = ants.morphology(tmp_mask1, operation='dilate', radius=2, mtype='binary', shape='ball')
            tmp_mask1 = ants.iMath(tmp_mask1, operation='GetLargestComponent')
            seg_tmp = ants.atropos(a=IMG, m='[0.1,1x1x1]', c='[3,0]', i='kmeans[3]', x=tmp_mask1)
            #  Clean up:
            seg_tmp2 = ants.iMath(seg_tmp['segmentation'], 'Pad', 10)
            W = ants.threshold_image(seg_tmp2, 3, 3, 1, 0, True)
            W = ants.iMath(W, operation='GetLargestComponent')
            W = W * 3
            G = ants.threshold_image(seg_tmp2, 2, 2, 1, 0, True)
            G = ants.iMath(G, operation='GetLargestComponent')
            TMP1 = ants.iMath(G, 'FillHoles', 2)
            G = G * TMP1
            C = ants.threshold_image(seg_tmp2, 1, 1, 1, 0, True)
            TMP2 = ants.morphology(C, operation='erode', radius=10, mtype='binary', shape='ball')
            G[G == 0] = TMP2[G == 0]
            G = G * 2
            seg_tmp2 = W
            seg_tmp2[W == 0] = G[W == 0]
            #  clean the Brainmask
            tmp_mask3 = ants.threshold_image(seg_tmp2, 3, 3, 1, 0, True)
            TMP3 = ants.threshold_image(seg_tmp2, 2, 2, 1, 0, True)
            tmp_mask3[tmp_mask3 == 0] = TMP3[tmp_mask3 == 0]
            tmp_mask3 = ants.morphology(tmp_mask3, operation='erode', radius=2, mtype='binary', shape='ball')
            tmp_mask3 = ants.iMath(tmp_mask3, operation='GetLargestComponent')
            tmp_mask3 = ants.morphology(tmp_mask3, operation='dilate', radius=4, mtype='binary', shape='ball')
            tmp_mask3 = ants.iMath(tmp_mask3, 'FillHoles', 2)
            # tmp_mask1 = ants.iMath(tmp_mask1,'Pad',10)
            # tmp_mask3[tmp_mask3==0]=tmp_mask1[tmp_mask3==0]
            tmp_mask3 = ants.morphology(tmp_mask3, operation='dilate', radius=5, mtype='binary', shape='ball')
            tmp_mask3 = ants.morphology(tmp_mask3, operation='erode', radius=5, mtype='binary', shape='ball')
            tmp_mask3 = ants.iMath(tmp_mask3, 'Pad', -10)
            ants.image_write(tmp_mask3, output_for_mask, ri=False)

            dictionary = {"Sources": [input_for_msk,
                                      BASE_SS_coregistr,
                                      BASE_SS_mask],
                          "Description": 'Brain mask (Atropos from ANTspy).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip == 'Custum_ANTS_Garin':
            nl = 'INFO: type_of_transform is ' + str(type_of_transform)
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

            IMG = ants.image_read(input_for_msk)
            REF_IMG = ants.image_read(BASE_SS_coregistr)

            mtx1 = ants.registration(fixed=IMG, moving=REF_IMG, type_of_transform='Translation',
                                     outprefix=opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_SyN_shift_'))
            MEAN_tr = ants.apply_transforms(fixed=IMG, moving=REF_IMG, transformlist=mtx1['fwdtransforms'],
                                            interpolator='nearestNeighbor')
            ants.image_write(MEAN_tr, opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_SyN_shift.nii.gz'),
                             ri=False)

            nl = 'INFO: type_of_transform is ' + str(type_of_transform)
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

            mTx = ants.registration(fixed=IMG, moving=REF_IMG,
                                    outprefix=opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_SyN_'),
                                    type_of_transform= type_of_transform,
                                    initial_transform=mtx1['fwdtransforms'],
                                    aff_metric=aff_metric_ants,
                                    grad_step=0.1,
                                    flow_sigma=3,
                                    total_sigma=0,
                                    aff_sampling=32,
                                    aff_random_sampling_rate=0.2,
                                    syn_sampling=32,
                                    aff_iterations=(1000, 500, 250, 100),
                                    aff_shrink_factors=(8, 4, 2, 1),
                                    aff_smoothing_sigmas=(3, 2, 1, 0),
                                    reg_iterations=(1000, 500, 250, 100),
                                    reg_smoothing_sigmas=(3, 2, 1, 0),
                                    reg_shrink_factors=(8, 4, 2, 1),
                                    verbose=True)
            nl = mTx
            print(nl)
            diary.write(f'\n{nl}')

            transfo_concat = [opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_SyN_1Warp.nii.gz'),
                              opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_SyN_0GenericAffine.mat')]

            REF_MASK = ants.image_read(BASE_SS_mask)
            tmp_mask1 = ants.apply_transforms(fixed=IMG, moving=REF_MASK,
                                              transformlist=transfo_concat, interpolator='nearestNeighbor')
            spacing = tmp_mask1.spacing  # This will give you the voxel size in x, y, z (e.g., (1.0, 1.0, 1.2) mm)
            # Use voxel size as sigma for Gaussian smoothing
            sigma = spacing  # Set sigma to voxel size for each dimension
            # Apply Gaussian smoothing (sigma controls the amount of smoothing)
            smoothed_mask = ants.smooth_image(tmp_mask1, sigma=sigma)  # Adjust sigma for more or less smoothing
            # Threshold to return to binary mask
            binary_smoothed_mask = ants.threshold_image(smoothed_mask, low_thresh=0.5, high_thresh=1)
            binary_smoothed_mask = ants.iMath(binary_smoothed_mask, operation='GetLargestComponent')
            ants.image_write(binary_smoothed_mask, output_for_mask, ri=False)

            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
                      ' -input ' + output_for_mask + ' -fill_holes'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": [input_for_msk,
                                      BASE_SS_coregistr,
                                      BASE_SS_mask],
                          "Description": 'Brain mask (Atropos from ANTspy).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip.startswith('CustumANTSGarin_'):
            # Extract the percentile from the string
            dilate_b = int(brain_skullstrip.split('_')[1])

            nl = 'INFO: type_of_transform is ' + str(type_of_transform)
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

            IMG = ants.image_read(input_for_msk)
            REF_IMG = ants.image_read(BASE_SS_coregistr)

            mtx1 = ants.registration(fixed=IMG, moving=REF_IMG, type_of_transform='Translation',
                                     outprefix=opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_SyN_shift_'))
            MEAN_tr = ants.apply_transforms(fixed=IMG, moving=REF_IMG, transformlist=mtx1['fwdtransforms'],
                                            interpolator='nearestNeighbor')
            ants.image_write(MEAN_tr, opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_SyN_shift.nii.gz'),
                             ri=False)

            nl = 'INFO: type_of_transform is ' + str(type_of_transform)
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

            mTx = ants.registration(fixed=IMG, moving=REF_IMG,
                                    outprefix=opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_SyN_'),
                                    type_of_transform= type_of_transform,
                                    initial_transform=mtx1['fwdtransforms'],
                                    aff_metric=aff_metric_ants,
                                    grad_step=0.1,
                                    flow_sigma=3,
                                    total_sigma=0,
                                    aff_sampling=32,
                                    aff_random_sampling_rate=0.2,
                                    syn_sampling=32,
                                    aff_iterations=(1000, 500, 250, 100),
                                    aff_shrink_factors=(8, 4, 2, 1),
                                    aff_smoothing_sigmas=(3, 2, 1, 0),
                                    reg_iterations=(1000, 500, 250, 100),
                                    reg_smoothing_sigmas=(3, 2, 1, 0),
                                    reg_shrink_factors=(8, 4, 2, 1),
                                    verbose=True)
            nl = mTx
            print(nl)
            diary.write(f'\n{nl}')

            transfo_concat = [opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_SyN_1Warp.nii.gz'),
                              opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_SyN_0GenericAffine.mat')]

            REF_MASK = ants.image_read(BASE_SS_mask)
            tmp_mask1 = ants.apply_transforms(fixed=IMG, moving=REF_MASK,
                                              transformlist=transfo_concat, interpolator='nearestNeighbor')
            spacing = tmp_mask1.spacing  # This will give you the voxel size in x, y, z (e.g., (1.0, 1.0, 1.2) mm)
            # Use voxel size as sigma for Gaussian smoothing
            sigma = spacing  # Set sigma to voxel size for each dimension
            # Apply Gaussian smoothing (sigma controls the amount of smoothing)
            smoothed_mask = ants.smooth_image(tmp_mask1, sigma=sigma)  # Adjust sigma for more or less smoothing
            # Threshold to return to binary mask
            binary_smoothed_mask = ants.threshold_image(smoothed_mask, low_thresh=0.5, high_thresh=1)
            binary_smoothed_mask = ants.iMath(binary_smoothed_mask, operation='GetLargestComponent')
            ants.image_write(binary_smoothed_mask, output_for_mask, ri=False)

            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
                      ' -input ' + output_for_mask + ' -fill_holes -dilate_input ' + str(dilate_b) + ' -' + str(dilate_b)
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": [input_for_msk,
                                      BASE_SS_coregistr,
                                      BASE_SS_mask],
                          "Description": 'Brain mask (Atropos from ANTspy).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip == '3dSkullStrip_noDilate':
            command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
                      '-input ' + input_for_msk + ' -blur_fwhm 2 -orig_vol -mask_vol -use_skull'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
                      ' -input ' + output_for_mask + ' -fill_holes'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (3dSkullStrip from AFNI).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip == '3dSkullStrip':
            command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
                      '-input ' + input_for_msk + ' -blur_fwhm 2 -orig_vol -mask_vol -use_skull'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
                      ' -input ' + output_for_mask + ' -fill_holes -dilate_input 2'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (3dSkullStrip from AFNI).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip.startswith('_bet'):
            # Extract the last two digits to use as the -f value
            f_value = brain_skullstrip[-4:]
            # Create the approximate brain mask using bet2
            command = f'singularity run {s_bind}{fsl_sif} bet2 {input_for_msk} ' \
                      f'{opj(process_dir, extract_filename(input_for_msk) + "_bet" + masking_img + ".nii.gz")} -f {f_value}'
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            # Run the AFNI 3dcalc command to create the final mask
            command = f'singularity run {s_bind}{afni_sif} 3dcalc -a ' \
                      f'{opj(process_dir, extract_filename(input_for_msk) + "_bet" + masking_img + ".nii.gz")} ' \
                      f'-expr "step(a)" -prefix {output_for_mask} -overwrite'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (bet from FSL).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip.startswith('CustumNilearn_'):
            # Extract the cutoff values from the string
            _, lower_cutoff, upper_cutoff = brain_skullstrip.split('_')
            lower_cutoff = float(lower_cutoff)
            upper_cutoff = float(upper_cutoff)

            # Convert to float
            command = f'singularity run {s_bind}{fs_sif} mri_convert -odt float {input_for_msk} {opj(process_dir, extract_filename(input_for_msk))}_float.nii.gz'
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            # Compute the EPI mask using nilearn with the given cutoff values
            mask_img = compute_epi_mask(f"{opj(process_dir, extract_filename(input_for_msk))}_float.nii.gz", lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff, connected=True, opening=3,
                                        exclude_zeros=False, ensure_finite=True)
            mask_img.to_filename(output_for_mask)
            # Use AFNI to process the mask
            command = f'singularity run {s_bind}{afni_sif} 3dmask_tool -overwrite -prefix {output_for_mask} -input {output_for_mask} -fill_holes -dilate_input 1'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

            tmp_mask1 = ants.image_read(output_for_mask)
            spacing = tmp_mask1.spacing  # This will give you the voxel size in x, y, z (e.g., (1.0, 1.0, 1.2) mm)
            # Use voxel size as sigma for Gaussian smoothing
            sigma = spacing  # Set sigma to voxel size for each dimension
            # Apply Gaussian smoothing (sigma controls the amount of smoothing)
            smoothed_mask = ants.smooth_image(tmp_mask1, sigma=sigma)  # Adjust sigma for more or less smoothing
            # Threshold to return to binary mask
            binary_smoothed_mask = ants.threshold_image(smoothed_mask, low_thresh=0.5, high_thresh=1)
            binary_smoothed_mask = ants.iMath(binary_smoothed_mask, operation='GetLargestComponent')
            ants.image_write(binary_smoothed_mask, output_for_mask, ri=False)
            # Resample the mask image
            command = f'singularity run {s_bind}{afni_sif} 3dresample -master {input_for_msk} -prefix {output_for_mask} -input {output_for_mask} -overwrite -bound_type SLAB'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (compute_epi_mask from Nilearn).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip.startswith('CustumNilearnExcludeZeros_'):
            # Extract the cutoff values from the string
            _, lower_cutoff, upper_cutoff = brain_skullstrip.split('_')
            lower_cutoff = float(lower_cutoff)
            upper_cutoff = float(upper_cutoff)
            # Convert to float
            command = f'singularity run {s_bind}{fs_sif} mri_convert -odt float {input_for_msk} {opj(process_dir, extract_filename(input_for_msk))}_float.nii.gz'
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            # Compute the EPI mask using nilearn with the given cutoff values
            mask_img = compute_epi_mask(f"{opj(process_dir, extract_filename(input_for_msk))}_float.nii.gz", lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff, connected=True, opening=3,
                                        exclude_zeros=True, ensure_finite=True)
            mask_img.to_filename(output_for_mask)
            # Use AFNI to process the mask
            command = f'singularity run {s_bind}{afni_sif} 3dmask_tool -overwrite -prefix {output_for_mask} -input {output_for_mask} -fill_holes -dilate_input 1'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            tmp_mask1 = ants.image_read(output_for_mask)
            spacing = tmp_mask1.spacing  # This will give you the voxel size in x, y, z (e.g., (1.0, 1.0, 1.2) mm)
            # Use voxel size as sigma for Gaussian smoothing
            sigma = spacing  # Set sigma to voxel size for each dimension
            # Apply Gaussian smoothing (sigma controls the amount of smoothing)
            smoothed_mask = ants.smooth_image(tmp_mask1, sigma=sigma)  # Adjust sigma for more or less smoothing
            # Threshold to return to binary mask
            binary_smoothed_mask = ants.threshold_image(smoothed_mask, low_thresh=0.5, high_thresh=1)
            binary_smoothed_mask = ants.iMath(binary_smoothed_mask, operation='GetLargestComponent')
            ants.image_write(binary_smoothed_mask, output_for_mask, ri=False)
            # Resample the mask image
            command = f'singularity run {s_bind}{afni_sif} 3dresample -master {input_for_msk} -prefix {output_for_mask} -input {output_for_mask} -overwrite -bound_type SLAB'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (compute_epi_mask from Nilearn).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        # Handle 'Custum' type skullstripping dynamically for percentile
        elif brain_skullstrip.startswith('CustumThreshold_'):
            # Extract the percentile from the string
            percentile = int(brain_skullstrip.split('_')[1])
            # Load the image and calculate the threshold at the given percentile
            loadimg = nib.load(input_for_msk).get_fdata()
            loadimgsort = np.percentile(np.abs(loadimg)[np.abs(loadimg) > 0], percentile)
            # Threshold the image using nilearn
            mask_imag = nilearn.image.threshold_img(input_for_msk, threshold=loadimgsort, cluster_threshold=10)
            mask_imag.to_filename(output_for_mask)
            # Use AFNI to process the mask
            command = f'singularity run {s_bind}{afni_sif} 3dmask_tool -overwrite -prefix {output_for_mask} -input {output_for_mask} -fill_holes -dilate_input 1'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            tmp_mask1 = ants.image_read(output_for_mask)
            spacing = tmp_mask1.spacing  # This will give you the voxel size in x, y, z (e.g., (1.0, 1.0, 1.2) mm)
            # Use voxel size as sigma for Gaussian smoothing
            sigma = spacing  # Set sigma to voxel size for each dimension
            # Apply Gaussian smoothing (sigma controls the amount of smoothing)
            smoothed_mask = ants.smooth_image(tmp_mask1, sigma=sigma)  # Adjust sigma for more or less smoothing
            # Threshold to return to binary mask
            binary_smoothed_mask = ants.threshold_image(smoothed_mask, low_thresh=0.5, high_thresh=1)
            binary_smoothed_mask = ants.iMath(binary_smoothed_mask, operation='GetLargestComponent')
            ants.image_write(binary_smoothed_mask, output_for_mask, ri=False)
            # Resample the mask image
            command = f'singularity run {s_bind}{afni_sif} 3dresample -master {input_for_msk} -prefix {output_for_mask} -input {output_for_mask} -overwrite -bound_type SLAB'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (threshold_img from Nilearn).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip.startswith('Vol_sammba_'):
            volume = int(brain_skullstrip.split('_')[2])
            command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            nichols_masker = Histrogram_mask_EMB.HistogramMask()
            nichols_masker.inputs.in_file = input_for_msk[:-7] + '_float.nii.gz'
            nichols_masker.inputs.volume_threshold = int(volume)
            # nichols_masker.inputs.upper_cutoff = 0.2
            # nichols_masker.inputs.lower_cutoff = 0.8
            # nichols_masker.inputs.intensity_threshold = 500
            # nichols_masker.inputs.opening = 2
            # nichols_masker.inputs.closing = 10
            nichols_masker.inputs.dilation_size = (1, 2, 3)
            nichols_masker.inputs.connected = True
            nichols_masker.inputs.out_file = output_for_mask
            res = nichols_masker.run()  # doctest: +SKIP
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (Histrogram_mask_EMB from sammba).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip =='QWARP':
            command = 'singularity run' + s_bind + afni_sif + '3dQwarp -overwrite -iwarp' + \
            ' -base ' + BASE_SS_coregistr + \
            ' -prefix ' + opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_AFNIQ.nii.gz') + \
            ' -source ' + input_for_msk + ' -maxlev 5 -resample'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            command = 'singularity run' + s_bind + afni_sif + '3dNwarpApply -nwarp ' + opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_AFNIQ_WARPINV.nii.gz') + \
            ' -source ' + BASE_SS_mask + ' -master ' + input_for_msk + ' -interp NN' + \
            ' -prefix ' + opj(opj(process_dir, extract_filename(input_for_msk)) + masking_img + 'template_brainmask.nii.gz') + ' -overwrite'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            Ex_Mask    = opj(opj(process_dir, extract_filename(input_for_msk)) + 'mask_tmp' + masking_img + '.nii.gz')
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + Ex_Mask + \
            ' -input ' + opj(opj(process_dir, extract_filename(input_for_msk)) + masking_img + 'template_brainmask.nii.gz') + ' -fill_holes' # -dilate_input 2'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            shutil.copyfile(Ex_Mask, output_for_mask)
            dictionary = {"Sources": [input_for_msk,
                                      BASE_SS_coregistr,
                                      BASE_SS_mask],
                          "Description": 'Brain mask (3dQwarp from AFNI).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip =='synthstrip':
                command = 'singularity run' + s_bind + strip_sif + \
                ' -o ' + opj(opj(process_dir, extract_filename(input_for_msk)) + 'skullstriped.nii.gz') + \
                ' -m ' + output_for_mask + \
                          ' -i ' + input_for_msk
                nl= spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
                dictionary = {"Sources": input_for_msk,
                              "Description": 'Brain mask (synthstrip).', }
                json_object = json.dumps(dictionary, indent=2)
                with open(output_for_mask[:-7] + '.json', "w") as outfile:
                    outfile.write(json_object)

        elif brain_skullstrip == 'muSkullStrip_Human':
            command = 'python3 ' + opj(opd(afni_sif), 'NHP-BrainExtraction', 'UNet_Model', 'muSkullStrip.py') + \
                      ' -in ' + input_for_msk + \
                      ' -model ' + opj(opd(afni_sif), 'NHP-BrainExtraction', 'UNet_Model', 'models','Site-Human-epoch_08.model') + \
                      ' -out ' + process_dir
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            shutil.copyfile(opj(opd(output_for_mask), extract_filename(input_for_msk) + '_pre_mask.nii.gz'), output_for_mask)
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (U-Net).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip == 'Custum_Macaque2':
            import re
            # Use a regular expression to find the session number (ses-XX)
            match = re.search(r'/ses-(\d+)/', process_dir)
            if match:
                Session = match.group(1)
                nl = f"Session number: {Session}"
            else:
                nl = "No session number found."
            print(nl)
            diary.write(f'\n{nl}')

            if ID == 'Pickle' and Session == 5:
                command = 'python3 ' + opj(opd(afni_sif),'NHP-BrainExtraction','UNet_Model','muSkullStrip.py') + \
                          ' -in ' + input_for_msk + \
                          ' -model ' + opj(opd(afni_sif),'NHP-BrainExtraction','UNet_Model','models','model-20_macaque-epoch') + \
                          ' -out ' + process_dir
                nl= spgo(command)
                diary.write(f'\n{nl}')
                print(nl)

                shutil.copyfile(opj(opd(output_for_mask), extract_filename(input_for_msk) + '_pre_mask.nii.gz'), output_for_mask)
                command = f'singularity run {s_bind}{afni_sif} 3dmask_tool -overwrite -prefix {output_for_mask} -input {output_for_mask} -fill_holes -dilate_input 10'
                nl= spgo(command)
                diary.write(f'\n{nl}')
                print(nl)

                dictionary = {"Sources": input_for_msk,
                              "Description": 'Brain mask (U-Net).', }
                json_object = json.dumps(dictionary, indent=2)
                with open(output_for_mask[:-7] + '.json', "w") as outfile:
                    outfile.write(json_object)

            elif ID == 'Trinity' and Session == 6:
                command = 'python3 ' + opj(opd(afni_sif),'NHP-BrainExtraction','UNet_Model','muSkullStrip.py') + \
                          ' -in ' + input_for_msk + \
                          ' -model ' + opj(opd(afni_sif),'NHP-BrainExtraction','UNet_Model','models','model-20_macaque-epoch') + \
                          ' -out ' + process_dir
                nl= spgo(command)
                diary.write(f'\n{nl}')
                print(nl)

                shutil.copyfile(opj(opd(output_for_mask), extract_filename(input_for_msk) + '_pre_mask.nii.gz'), output_for_mask)
                command = f'singularity run {s_bind}{afni_sif} 3dmask_tool -overwrite -prefix {output_for_mask} -input {output_for_mask} -fill_holes -dilate_input 10'
                nl= spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
                dictionary = {"Sources": input_for_msk,
                              "Description": 'Brain mask (U-Net).', }
                json_object = json.dumps(dictionary, indent=2)
                with open(output_for_mask[:-7] + '.json', "w") as outfile:
                    outfile.write(json_object)

            else:
                command = 'python3 ' + opj(opd(afni_sif),'NHP-BrainExtraction','UNet_Model','muSkullStrip.py') + \
                          ' -in ' + input_for_msk + \
                          ' -model ' + opj(opd(afni_sif),'NHP-BrainExtraction','UNet_Model','models','model-20_macaque-epoch') + \
                          ' -out ' + process_dir
                nl= spgo(command)
                diary.write(f'\n{nl}')
                print(nl)

                shutil.copyfile(opj(opd(output_for_mask), extract_filename(input_for_msk) + '_pre_mask.nii.gz'), output_for_mask)
                command = f'singularity run {s_bind}{afni_sif} 3dmask_tool -overwrite -prefix {output_for_mask} -input {output_for_mask} -fill_holes -dilate_input 3'
                nl= spgo(command)
                diary.write(f'\n{nl}')
                print(nl)

                dictionary = {"Sources": input_for_msk,
                              "Description": 'Brain mask (U-Net).', }
                json_object = json.dumps(dictionary, indent=2)
                with open(output_for_mask[:-7] + '.json', "w") as outfile:
                    outfile.write(json_object)


        elif brain_skullstrip =='Custum_QWARP':
                command = 'singularity run' + s_bind + afni_sif + '3dQwarp -overwrite -iwarp' + \
                ' -base ' + BASE_SS_coregistr + \
                ' -prefix ' + opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_AFNIQ.nii.gz') + \
                ' -source ' + input_for_msk + ' -maxlev 5 -resample'
                nl= spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
                command = 'singularity run' + s_bind + afni_sif + '3dNwarpApply -nwarp ' + opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_AFNIQ_WARPINV.nii.gz') + \
                ' -source ' + BASE_SS_mask + ' -master ' + input_for_msk + ' -interp NN' + \
                ' -prefix ' + opj(opj(process_dir, extract_filename(input_for_msk)) + masking_img + 'template_brainmask.nii.gz') + ' -overwrite'
                nl= spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
                Ex_Mask    = opj(opj(process_dir, extract_filename(input_for_msk)) + 'mask_tmp' + masking_img + '.nii.gz')
                command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + Ex_Mask + \
                ' -input ' + opj(opj(process_dir, extract_filename(input_for_msk)) + masking_img + 'template_brainmask.nii.gz') + ' -fill_holes' # -dilate_input 2'
                nl= spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
                shutil.copyfile(Ex_Mask,output_for_mask)
                dictionary = {"Sources": [input_for_msk,
                                          BASE_SS_coregistr,
                                          BASE_SS_mask],
                              "Description": 'Brain mask (3dQwarp from AFNI).', }
                json_object = json.dumps(dictionary, indent=2)
                with open(output_for_mask[:-7] + '.json', "w") as outfile:
                    outfile.write(json_object)

        elif brain_skullstrip =='Custum_QWARPT2':
            command = 'singularity run' + s_bind + afni_sif + '3dQwarp -overwrite -lpa -iwarp' + \
            ' -base ' + BASE_SS_coregistr + \
            ' -prefix ' + opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_AFNIQ.nii.gz') + \
            ' -source ' + input_for_msk + ' -maxlev 3 -resample'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            command = 'singularity run' + s_bind + afni_sif + '3dNwarpApply -nwarp ' + opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_AFNIQ_WARPINV.nii.gz') + \
            ' -source ' + BASE_SS_mask + ' -master ' + input_for_msk + ' -interp NN' + \
            ' -prefix ' + opj(opj(process_dir, extract_filename(input_for_msk)) + masking_img + 'template_brainmask.nii.gz') + ' -overwrite'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            Ex_Mask    = opj(opj(process_dir, extract_filename(input_for_msk)) + 'mask_tmp' + masking_img + '.nii.gz')
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + Ex_Mask + \
            ' -input ' + opj(opj(process_dir, extract_filename(input_for_msk)) + masking_img + 'template_brainmask.nii.gz') + ' -fill_holes'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            shutil.copyfile(Ex_Mask, output_for_mask)
            dictionary = {"Sources": [input_for_msk,
                                      BASE_SS_coregistr,
                                      BASE_SS_mask],
                          "Description": 'Brain mask (3dQwarp from AFNI).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip =='Custum_QWARPT2_dil':
            command = 'singularity run' + s_bind + afni_sif + '3dQwarp -overwrite -lpa -iwarp' + \
            ' -base ' + BASE_SS_coregistr + \
            ' -prefix ' + opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_AFNIQ.nii.gz') + \
            ' -source ' + input_for_msk + ' -maxlev 3 -resample'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            command = 'singularity run' + s_bind + afni_sif + '3dNwarpApply -nwarp ' + opj(opj(process_dir, extract_filename(input_for_msk)) + 'template_to_' + masking_img + '_AFNIQ_WARPINV.nii.gz') + \
            ' -source ' + BASE_SS_mask + ' -master ' + input_for_msk + ' -interp NN' + \
            ' -prefix ' + opj(opj(process_dir, extract_filename(input_for_msk)) + masking_img + 'template_brainmask.nii.gz') + ' -overwrite'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            Ex_Mask    = opj(opj(process_dir, extract_filename(input_for_msk)) + 'mask_tmp' + masking_img + '.nii.gz')
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + Ex_Mask + \
            ' -input ' + opj(opj(process_dir, extract_filename(input_for_msk)) + masking_img + 'template_brainmask.nii.gz') + ' -fill_holes -dilate_input 2'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": [input_for_msk,
                                      BASE_SS_coregistr,
                                      BASE_SS_mask],
                          "Description": 'Brain mask (3dQwarp from AFNI).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

            shutil.copyfile(Ex_Mask,output_for_mask)

        elif brain_skullstrip == 'NoSkullStrip':
            command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + input_for_msk + \
                      ' -expr "step(a)" -prefix ' + output_for_mask + ' -overwrite'
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (positive voxels).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif brain_skullstrip == 'Manual':
            def run_command_and_wait(command):
                nl = 'INFO: Running command:' + command
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                result = subprocess.run(command, shell=True)
                if result.returncode == 0:
                    nl_final = 'INFO: Command completed successfully.'
                    print(bcolors.OKGREEN + nl_final + bcolors.ENDC)
                else:
                    nl_final = 'WARNING: Command failed with return code:' + result.returncode
                    print(bcolors.WARNING + nl_final + bcolors.ENDC)
                return nl_final

            if not os.path.exists(output_for_mask):
                command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + input_for_msk + ' -expr "step(a)" -prefix ' + output_for_mask
                nl= spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
            command = ('singularity run' + s_bind + itk_sif + 'itksnap -g ' + input_for_msk + ' -s ' + output_for_mask)
            nl_final=run_command_and_wait(command)
            diary.write(f'\n{nl_final}')
            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (manual drawing with itksnap).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        else:
            nl = "ERROR: brain_skullstrip not recognized, check that brain_skullstrip_1 or brain_skullstrip_2 are correctly written!!"
            diary.write(f'\n{nl}')
            raise Exception(bcolors.FAIL + nl + bcolors.ENDC)

        if check_visualy_final_mask == True:
            if step_skullstrip == 1:
                nl = 'WARNING1: any modifications should be save as : ' + opj(opj(process_dir, extract_filename(input_for_msk)) + masking_img + 'final_mask.nii.gz')
                print(bcolors.WARNING + nl + bcolors.ENDC)
            elif step_skullstrip == 2:
                nl = 'WARNING1: any modifications should be save as : ' + opj(opj(process_dir, extract_filename(input_for_msk)) + masking_img + 'final_mask_2.nii.gz')
                print(nl)
            diary.write(f'\n{nl}')
            nl = 'WARNING2: These modifications will automatically be taken as "final mask", delete the this file if you want to use a Skulltrip method!!'
            print(bcolors.WARNING + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

            def run_command_and_wait(command):
                nl = "Running command:" + command
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                result = subprocess.run(command, shell=True)
                if result.returncode == 0:
                    nl_final = 'INFO: Command completed successfully.'
                    print(bcolors.OKGREEN + nl_final + bcolors.ENDC)
                else:
                    nl_final = 'WARNING: Command failed with return code:' + result.returncode
                    print(bcolors.WARNING + nl_final + bcolors.ENDC)
                return nl_final

            # Example usage
            command = ('singularity run' + s_bind + itk_sif + 'itksnap -g ' + input_for_msk + ' -s ' + output_for_mask)
            nl_final = run_command_and_wait(command)
            diary.write(f'\n{nl_final}')

            if os.path.exists(opj(masks_dir, ID + masking_img + 'final_mask.nii.gz')):
                shutil.copyfile(opj(masks_dir, ID + masking_img + 'final_mask.nii.gz'), output_for_mask)

            elif os.path.exists(opj(masks_dir, ID + masking_img + 'final_mask_2.nii.gz')):
                shutil.copyfile(opj(masks_dir, ID + masking_img + 'final_mask.nii.gz'), output_for_mask)

            elif os.path.exists(opj(study_template_atlas_folder, 'studytemplate2_' + type_norm,
                              'studytemplate_final_mask_study_template.nii.gz')):
                shutil.copyfile(opj(study_template_atlas_folder, 'studytemplate2_' + type_norm,
                              'studytemplate_final_mask_study_template.nii.gz'), output_for_mask)

            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (manual drawing with itksnap).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

    return(output_for_mask)

    diary.write(f'\n')
    diary.close()



