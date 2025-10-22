import os
import subprocess
import json
import re
import nibabel as nib

from nilearn.image import resample_to_img
from nilearn.image.image import mean_img

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile



from Tools import run_cmd,get_orientation
from fonctions.extract_filename import extract_filename
from fonctions import Skullstrip_func

def Refimg_to_meanfMRI(anat_func_same_space, BASE_SS_coregistr,TfMRI , dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, use_master_for_Allineate,
                       dir_fMRI_Refth_RS_prepro3, RS, nb_run, REF_int, ID, dir_prepro, brainmask, V_mask, W_mask, G_mask, dilate_mask,
                       costAllin, anat_subject, Method_mask_func, overwrite, type_of_transform, aff_metric_ants,
                       sing_afni,fs_sif, sing_fsl, sing_itk,diary_file):

    nl = '##  Working on step ' + str(3) + '(function: _3_mask_fMRI).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    ##### create new variable for template space (we will need to store and downsample template img to func resolution)
    if ope(dir_fMRI_Refth_RS_prepro3) == False:

        os.makedirs(dir_fMRI_Refth_RS_prepro3)
        os.makedirs(opj(dir_fMRI_Refth_RS_prepro3,'tmp'))

    ### create a list of the image to be corrected

    root_RS_ref = extract_filename(RS[REF_int])

    MEAN_im_list = ' ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS_ref + '_xdtr_mean_deob.nii.gz') #image "reference" (to be created)
    MEAN_im_list_1    = [opj(dir_fMRI_Refth_RS_prepro1,root_RS_ref + '_xdtr_mean_deob.nii.gz')] #image "reference" (to be created)

    for r in range(int(nb_run)):
        root_RS = extract_filename(RS[r])
        if not RS[r] == RS[REF_int]: # do not process ref...
            # add the co-registered mean image to the list
            MEAN_im_list = MEAN_im_list + ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref.nii.gz') #add images in the same space
            MEAN_im_list_1.append(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref.nii.gz'))

    ################################################################# create a mean image to use for the anat to func and recenter

    ###### average all func data and clean the image #####

    command = (sing_afni + '3dTcat' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'residual_motion.nii.gz') + MEAN_im_list)
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": MEAN_im_list,
                  "Description": '4D concatenation (3dTcat,AFNI).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro1,'residual_motion.json'), "w") as outfile:
        outfile.write(json_object)

    #################################### production of Mean image ####################################

    mean_haxby = mean_img(MEAN_im_list_1)
    mean_haxby.to_filename(opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz'))
    dictionary = {"Sources": MEAN_im_list_1,
                  "Description": 'mean image (mean_image, nilearn).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.json'), "w") as outfile:
        outfile.write(json_object)

    mean_haxby.to_filename(opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_test.nii.gz'))
    dictionary = {"Sources": MEAN_im_list_1,
                  "Description": 'mean image (mean_image, nilearn).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_test.json'), "w") as outfile:
        outfile.write(json_object)

    # Load the image directly
    mean_img_path = opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz')
    img = nib.load(mean_img_path)
    # Get voxel sizes
    delta_x, delta_y, delta_z = [str(round(abs(x), 10)) for x in img.header.get_zooms()[:3]]


    # Usage
    orient_meanimg = get_orientation.get_orientation_nibabel(mean_img_path)
    nl = 'Orientation: ' + orient_meanimg
    run_cmd.msg(nl, diary_file, 'OKGREEN')

    command = (sing_afni + '3dcalc' + overwrite + ' -a ' + BASE_SS_coregistr +
               ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + ' -expr "a"')
    run_cmd.do(command, diary_file)

    command = (sing_afni + '3dresample' + overwrite +
               ' -orient ' + orient_meanimg +
               ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') +
               ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z +
               ' -rmode Cu -input ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz'))
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": [BASE_SS_coregistr,
                              opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')],
                  "Description": 'Resampling (3dresample, AFNI).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro3, 'BASE_SS_fMRI.json'), "w") as outfile:
        outfile.write(json_object)

    # create anat space dir
    if ope(dir_fMRI_Refth_RS_prepro2) == False:
        os.makedirs(dir_fMRI_Refth_RS_prepro2)

    ############################### ############################### ############################### 
    ##                         resample the masks for signal extraction                         ###
    ############################### ############################### ###############################

    ##### first you need to re-create a dilate ref anat image (works better for co-registration anat to fMRI)
    # dilate a little bit MORE the "maskDilat"

    if dilate_mask != 0:
        command = (sing_afni + '3dmask_tool' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') +
                   ' -input ' + brainmask + ' -fill_holes -dilate_input ' + str(dilate_mask))
    else:
        command = (sing_afni + '3dmask_tool' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') +
                   ' -input ' + brainmask + ' -fill_holes')
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": brainmask,
                  "Description": 'dilation of a factor ' + str(dilate_mask) + ' (3dmask_tool, nilearn).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro2, 'maskDilatanat.json'), "w") as outfile:
        outfile.write(json_object)


    for input1, output2 in zip([anat_subject, brainmask, opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz'),
                                V_mask, W_mask, G_mask],
                               [opj(dir_fMRI_Refth_RS_prepro2, ('_').join([ID, 'res-func', TfMRI + '.nii.gz'])),
                                opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz'),
                                opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz'),
                                opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz'),
                                opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz'),
                                opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz')]):

        if ope(input1):
            if input1 == anat_subject:
                command = (sing_afni + '3dresample' + overwrite +
                           ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') +
                           ' -master ' + anat_subject +
                           ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz'))
                run_cmd.run(command, diary_file)

                # skullstrip the anat
                command = (sing_afni + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') +
                           ' -b ' + anat_subject +
                           ' -prefix ' + output2 + ' -expr "a*b"')
                run_cmd.do(command, diary_file)

                command = (sing_afni + '3dresample' + overwrite +
                           ' -prefix ' + output2 +
                           ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' +
                           ' -input  ' + output2)
                run_cmd.run(command, diary_file)

                dictionary = {"Sources": [anat_subject,
                                          opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz'),
                                          opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')],
                              "Description": 'Skull stripping and resampling (3dcalc and 3dresample, AFNI).', }
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.json'), "w") as outfile:
                    outfile.write(json_object)

            else:
                command = (sing_afni + '3dresample' + overwrite +
                           ' -prefix ' + output2 +
                           ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' +
                           ' -input  ' + input1)
                run_cmd.run(command, diary_file)

                dictionary = {"Sources": [input1,
                                          opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')],
                              "Description": 'resampling (3dresample, AFNI).', }
                json_object = json.dumps(dictionary, indent=2)
                with open(output2[:-7] + '.json', "w") as outfile:
                    outfile.write(json_object)
        else:
            nl =  ('WARNING:' + str(input1) + ' not found!!! this may be because you have not provided an aseg file, then no '
                                  'extraction of WM or Ventricles or GM will be possible... please check that!')
            run_cmd.msg(nl, diary_file, 'WARNING')

    nl = "INFO: brain_skullstrip method is " + Method_mask_func
    run_cmd.msg(nl, diary_file, 'OKGREEN')
    nl = 'INFO: looking for manual segmentation named:' + opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz') + '...'
    run_cmd.msg(nl, diary_file, 'OKGREEN')

    #### explore if manual_mask.nii.gz exists?
    if opi(opj(dir_fMRI_Refth_RS_prepro1,'manual_mask.nii.gz')):

        nl = 'WARNING: We found a final mask to skullstrip the functional image !!! no Skullstrip will be calculated!'
        run_cmd.msg(nl, diary_file, 'WARNING')
        nl = 'INFO: please delete' +opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz') + ' if you want retry to create a skulstripp images'
        run_cmd.msg(nl, diary_file, 'OKGREEN')


        command = (sing_afni + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1,'manual_mask.nii.gz') +
                   ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1,'maskDilat_Allineate_in_func.nii.gz') + ' -expr "a"')
        run_cmd.do(command, diary_file)

        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1,'manual_mask.nii.gz'),
                      "Description": 'Copy.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1,'maskDilat_Allineate_in_func.json'), "w") as outfile:
            outfile.write(json_object)

    #### if not, create an fMRI mask
    else:
        nl = "INFO: no manual mask found "
        run_cmd.msg(nl, diary_file, 'OKGREEN')
        if anat_func_same_space:
            command = (sing_afni + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro2, 'maskDilat.nii.gz') +
                       ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.nii.gz') + ' -expr "a"')
            run_cmd.do(command, diary_file)

            dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro2, 'maskDilat.nii.gz'),
                          "Description": 'Copy.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.json'), "w") as outfile:
                outfile.write(json_object)

            nl = 'you are using the mask from the anat img'
            run_cmd.msg(nl, diary_file, 'OKGREEN')
        else:

            Skullstrip_func.Skullstrip_func(Method_mask_func, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                                                      overwrite, costAllin, type_of_transform,
                                                      aff_metric_ants, sing_afni, sing_fsl, fs_sif, sing_itk, diary_file)
