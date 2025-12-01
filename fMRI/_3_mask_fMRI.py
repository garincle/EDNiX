import os
import json
import ants
import nibabel as nib

from nilearn.image import resample_to_img
from nilearn.image.image import mean_img

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile



from Tools import run_cmd,get_orientation
from fMRI.extract_filename import extract_filename
from fMRI import Skullstrip_func
from fMRI import plot_QC_func


def Refimg_to_meanfMRI(anat_func_same_space, BASE_SS_coregistr,TfMRI , dir_prepro_raw_process, dir_prepro_raw_masks, dir_prepro_acpc_masks, dir_prepro_acpc_process,
                       dir_prepro_template_process, RS, nb_run, REF_int, ID, dir_transfo, brainmask, V_mask, W_mask, G_mask, dilate_mask, n_for_ANTS, bids_dir,
                       costAllin, anat_subject, Method_mask_func, overwrite, type_of_transform, aff_metric_ants,
                       sing_afni, sing_fs, sing_fsl, sing_itk, diary_file):

    nl = '##  Working on step ' + str(3) + '(function: _3_mask_fMRI).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')
    ### create a list of the image to be corrected

    root_RS_ref = extract_filename(RS[REF_int])
    ffMRI_runMean_inRef_list = ' ' + opj(dir_prepro_raw_process, root_RS_ref + '_space-func_desc-fMRI_run_inRef.nii.gz')
    ffMRI_runMean_inRef_list1 = [opj(dir_prepro_raw_process, root_RS_ref + '_space-func_desc-fMRI_runMean_inRef.nii.gz')]
    residual_motion = opj(dir_prepro_raw_process, 'all_runs_space-func_desc-fMRI_residual_motion.nii.gz')
    Mean_Image = opj(dir_prepro_raw_process, 'all_runs_space-func_desc-fMRI_Mean_Image.nii.gz')
    BASE_SS_fMRI = opj(dir_prepro_template_process, 'BASE_SS_fMRI.nii.gz')
    BASE_SS = opj(dir_prepro_template_process, 'BASE_SS.nii.gz')
    maskDilatanat = opj(dir_prepro_acpc_masks, ID + '_space-acpc_mask_dilated.nii.gz')
    maskDilatfunc = opj(dir_prepro_acpc_masks, ID + '_space-acpc_mask_dilated_res_func.nii.gz')
    anat_res_func = opj(dir_prepro_acpc_process, ('_').join(['anat_space-acpc_res-func', TfMRI + '.nii.gz']))

    ### ref of the manual mask
    final_mask = opj(dir_prepro_raw_masks, ID + '_final_mask.nii.gz')
    Prepro_fMRI_mask = opj(dir_prepro_raw_masks, ID + '_fMRI_mask.nii.gz')

    for r in range(int(nb_run)):
        root_RS = extract_filename(RS[r])
        if not RS[r] == RS[REF_int]: # do not process ref...
            # add the co-registered mean image to the list
            fMRI_run_inRef = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_run_inRef.nii.gz')
            ffMRI_runMean_inRef_list = ffMRI_runMean_inRef_list + ' ' + fMRI_run_inRef #add images in the same space
            ffMRI_runMean_inRef_list1.append(fMRI_run_inRef)

    ################################################################# create a mean image to use for the anat to func and recenter
    ###### average all func data and clean the image #####

    command = (sing_afni + '3dTcat' + overwrite + ' -prefix ' + residual_motion + ffMRI_runMean_inRef_list)
    dictionary = {"Sources": ffMRI_runMean_inRef_list,
                  "Description": '4D concatenation (3dTcat,AFNI).',
                  "Command": command,}
    json_object = json.dumps(dictionary, indent=3)
    with open(residual_motion.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)
    run_cmd.run(command, diary_file)

    #################################### production of Mean image ####################################

    mean_haxby = mean_img(ffMRI_runMean_inRef_list1)
    mean_haxby.to_filename(Mean_Image)
    dictionary = {"Sources": ffMRI_runMean_inRef_list1,
                  "Description": 'mean image (mean_image, nilearn).',}
    json_object = json.dumps(dictionary, indent=3)
    with open(Mean_Image.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)

    # Load the image directly
    img = nib.load(Mean_Image)
    # Get voxel sizes
    delta_x, delta_y, delta_z = [str(round(abs(x), 10)) for x in img.header.get_zooms()[:3]]

    # Usage
    orient_meanimg = get_orientation.get_orientation_nibabel(Mean_Image)
    nl = 'Orientation: ' + orient_meanimg
    run_cmd.msg(nl, diary_file, 'OKGREEN')

    command = (sing_afni + '3dcalc -a ' + BASE_SS_coregistr +
               ' -prefix ' + BASE_SS +
               ' -expr "a"' + overwrite)
    dictionary = {"Sources": BASE_SS_coregistr,
                  "Description": 'copy.',
                  "Command": command, }
    json_object = json.dumps(dictionary, indent=3)
    with open(BASE_SS.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)
    run_cmd.do(command, diary_file)

    command = (sing_afni + '3dresample' + overwrite +
               ' -orient ' + orient_meanimg +
               ' -prefix ' + BASE_SS_fMRI +
               ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z +
               ' -rmode Cu -input ' + BASE_SS_coregistr)
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": [BASE_SS_coregistr,
                              Mean_Image],
                  "Description": 'Resampling (3dresample, AFNI).',
                  "Command": command, }
    json_object = json.dumps(dictionary, indent=3)
    with open(BASE_SS_fMRI.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)
    run_cmd.run(command, diary_file)

    ############################### ############################### ############################### 
    ##                         resample the masks for signal extraction                         ###
    ############################### ############################### ###############################

    ##### first you need to re-create a dilate ref anat image (works better for co-registration anat to fMRI)
    # dilate a little bit MORE the "maskDilat"

    if dilate_mask != 0:
        command = (sing_afni + '3dmask_tool' + overwrite + ' -prefix ' + maskDilatanat +
                   ' -input ' + brainmask + ' -fill_holes -dilate_input ' + str(dilate_mask))
    else:
        command = (sing_afni + '3dmask_tool' + overwrite + ' -prefix ' + maskDilatanat +
                   ' -input ' + brainmask + ' -fill_holes')
    dictionary = {"Sources": brainmask,
                  "Description": 'dilation of a factor ' + str(dilate_mask) + ' (3dmask_tool, nilearn).',
                  "Command": command, }
    json_object = json.dumps(dictionary, indent=3)
    with open(maskDilatanat.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)
    run_cmd.run(command, diary_file)

    for input1, output2 in zip([anat_subject, brainmask, maskDilatanat,
                                V_mask, W_mask, G_mask],
                               [anat_res_func,
                                opj(dir_prepro_acpc_masks,'mask_ref.nii.gz'),
                                maskDilatfunc,
                                opj(dir_prepro_acpc_masks,'Vmask.nii.gz'),
                                opj(dir_prepro_acpc_masks,'Wmask.nii.gz'),
                                opj(dir_prepro_acpc_masks,'Gmask.nii.gz')]):
        if ope(input1):
            if input1 == anat_subject:
                command = (sing_afni + '3dresample' + overwrite +
                           ' -prefix ' + maskDilatanat +
                           ' -master ' + anat_subject +
                           ' -input ' + maskDilatanat)
                run_cmd.run(command, diary_file)

                # skullstrip the anat
                command = (sing_afni + '3dcalc' + overwrite + ' -a ' + maskDilatanat +
                           ' -b ' + anat_subject +
                           ' -prefix ' + output2 + ' -expr "a*b"')
                run_cmd.do(command, diary_file)

                command = (sing_afni + '3dresample' + overwrite +
                           ' -prefix ' + output2 +
                           ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z +
                           ' -input ' + output2)
                run_cmd.run(command, diary_file)

                dictionary = {"Sources": [anat_subject,
                                          maskDilatanat,
                                          Mean_Image],
                              "Description": 'Skull stripping and resampling (3dcalc and 3dresample, AFNI).',
                              "Command": command, }
                json_object = json.dumps(dictionary, indent=3)
                with open(anat_res_func.replace('.nii.gz', '.json'), "w") as outfile:
                    outfile.write(json_object)
                run_cmd.run(command, diary_file)
            else:
                command = (sing_afni + '3dresample' + overwrite +
                           ' -prefix ' + output2 +
                           ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z +
                           ' -input ' + input1)
                run_cmd.run(command, diary_file)

                dictionary = {"Sources": [input1,
                                          Mean_Image],
                              "Description": 'resampling (3dresample, AFNI).',
                              "Command": command, }
                json_object = json.dumps(dictionary, indent=3)
                with open(output2.replace('.nii.gz', '.json'), "w") as outfile:
                    outfile.write(json_object)
                run_cmd.run(command, diary_file)
        else:
            nl =  ('WARNING:' + str(input1) + ' not found!!! this may be because you have not provided an aseg file, then no '
                                  'extraction of WM or Ventricles or GM will be possible... please check that!')
            run_cmd.msg(nl, diary_file, 'WARNING')

    nl = "INFO: brain_skullstrip method is " + Method_mask_func
    run_cmd.msg(nl, diary_file, 'OKGREEN')
    nl = 'INFO: looking for manual segmentation named:' + final_mask + '...'
    run_cmd.msg(nl, diary_file, 'OKGREEN')

    #### explore if manual_mask.nii.gz exists?
    if opi(final_mask):
        nl = 'WARNING: We found a final mask to skullstrip the functional image !!! no Skullstrip will be calculated!'

        nl = 'INFO: please delete' + final_mask + ' if you want retry to create a skulstripp images'
        run_cmd.msg(nl, diary_file, 'OKGREEN')

        command = (sing_afni + '3dresample' + overwrite + ' -input ' + final_mask + ' -master ' + Mean_Image +
                   ' -prefix ' + final_mask + ' -overwrite')
        run_cmd.do(command, diary_file)

        command = (sing_afni + '3dcalc' + overwrite + ' -a ' + final_mask +
                   ' -prefix ' + Prepro_fMRI_mask + ' -expr "a"')
        run_cmd.do(command, diary_file)

        dictionary = {"Sources": final_mask,
                      "Description": 'Copy.',
                      "Command": command, }
        json_object = json.dumps(dictionary, indent=3)
        with open(Prepro_fMRI_mask.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)
        run_cmd.run(command, diary_file)

    #### if not, create an fMRI mask
    else:
        nl = "INFO: no manual mask found "
        run_cmd.msg(nl, diary_file, 'OKGREEN')
        if anat_func_same_space:
            MEAN = ants.image_read(Mean_Image)
            IMG = ants.image_read(maskDilatanat)
            matrix = opj(dir_transfo, 'acpc_0GenericAffine.mat')
            MEAN_tr = ants.apply_transforms(fixed=MEAN, moving=IMG, whichtoinvert=[True],
                                            transformlist=matrix, interpolator='nearestNeighbor')
            ants.image_write(MEAN_tr, Prepro_fMRI_mask, ri=False)

            nl = 'you are using the mask from the anat img'
            run_cmd.msg(nl, diary_file, 'OKGREEN')

        else:
            Skullstrip_func.Skullstrip_func(Method_mask_func, Mean_Image, Prepro_fMRI_mask, anat_res_func, maskDilatfunc, dir_prepro_raw_process,
                                                      overwrite, costAllin, type_of_transform,
                                                      aff_metric_ants, sing_afni, sing_fsl, sing_fs, sing_itk, diary_file)

    if not ope(opj(bids_dir, 'QC')):
        os.mkdir(opj(bids_dir, 'QC'))
    if not ope(opj(bids_dir, 'QC', 'check_mask_fMRI')):
        os.mkdir(opj(bids_dir, 'QC', 'check_mask_fMRI'))

    ####plot the QC
    plot_QC_func.plot_qc(Mean_Image,
                         Prepro_fMRI_mask,
                         opj(bids_dir, 'QC', 'check_fmri_mask', root_RS_ref + '_check_mask_fMRI.png'))