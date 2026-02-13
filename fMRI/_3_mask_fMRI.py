import os
import json
import ants
import nibabel as nib
from scipy import ndimage
import numpy as np
from nilearn.image import mean_img
from Tools import run_cmd,get_orientation, check_nii
from fMRI.extract_filename import extract_filename
from fMRI import Skullstrip_func, plot_QC_func

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile



def Refimg_to_meanfMRI(MAIN_PATH, anat_func_same_space, BASE_SS_coregistr,TfMRI , dir_prepro_raw_process, dir_prepro_raw_masks, dir_prepro_acpc_masks, dir_prepro_acpc_process,
                       dir_prepro_template_process, RS, nb_run, REF_int, ID, dir_transfo, brainmask, V_mask, W_mask, G_mask, WBG_mask, dilate_mask, n_for_ANTS, bids_dir,
                       costAllin, anat_subject, Method_mask_func, overwrite, type_of_transform, aff_metric_ants, extra_erode,
                       sing_afni, sing_fs, sing_fsl, sing_itk, diary_file):

    nl = '##  Working on step ' + str(3) + '(function: _3_mask_fMRI).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')
    ### create a list of the image to be corrected

    root_RS_ref = extract_filename(RS[REF_int])
    Mean_Image = opj(dir_prepro_raw_process, 'all_runs_space-func_desc-fMRI_Mean_Image.nii.gz')
    BASE_SS_fMRI = opj(dir_prepro_template_process, 'BASE_SS_fMRI.nii.gz')
    BASE_SS = opj(dir_prepro_template_process, 'BASE_SS.nii.gz')
    maskDilatanat = opj(dir_prepro_acpc_masks, ID + '_space-acpc_mask_dilated.nii.gz')
    maskDilatfunc = opj(dir_prepro_acpc_masks, ID + '_space-acpc_mask_dilated_res_func.nii.gz')
    anat_res_func = opj(dir_prepro_acpc_process, ('_').join(['anat_space-acpc_res-func', TfMRI + '.nii.gz']))

    ### ref of the manual mask
    final_mask = opj(dir_prepro_raw_masks, ID + '_final_mask.nii.gz')
    Prepro_fMRI_mask = opj(dir_prepro_raw_masks, ID + '_fMRI_mask.nii.gz')

    ffMRI_runMean_inRef_list1 = []
    for r in range(int(nb_run)):
        root_RS = extract_filename(RS[r])
        fMRI_run_inRef = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_runMean_inRef.nii.gz')
        ffMRI_runMean_inRef_list1.append(fMRI_run_inRef)

    ################################################################# create a mean image to use for the anat to func and recenter
    ###### average all func data and clean the image #####
    '''
    command = (sing_afni + '3dTcat' + overwrite + ' -prefix ' + residual_motion + ffMRI_runMean_inRef_list)
    dictionary = {"Sources": ffMRI_runMean_inRef_list,
                  "Description": '4D concatenation (3dTcat,AFNI).',
                  "Command": command,}
    json_object = json.dumps(dictionary, indent=3)
    with open(residual_motion.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)
    run_cmd.run(command, diary_file)
    '''
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
                                V_mask, W_mask, G_mask, WBG_mask],
                               [anat_res_func,
                                opj(dir_prepro_acpc_masks,'mask_ref.nii.gz'),
                                maskDilatfunc,
                                opj(dir_prepro_acpc_masks,'Vmask.nii.gz'),
                                opj(dir_prepro_acpc_masks,'Wmask.nii.gz'),
                                opj(dir_prepro_acpc_masks,'Gmask.nii.gz'),
                                opj(dir_prepro_acpc_masks,'WBGmask.nii.gz')]):
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
        run_cmd.msg(nl, diary_file, 'OKGREEN')

        nl = 'INFO: please delete' + final_mask + ' if you want retry to create a skulstripp images'
        run_cmd.msg(nl, diary_file, 'OKGREEN')

        command = (sing_afni + '3dresample' + overwrite + ' -input ' + final_mask + ' -master ' + Mean_Image +
                   ' -prefix ' + Prepro_fMRI_mask + ' -overwrite')
        run_cmd.do(command, diary_file)
        check_nii.keep_header(Prepro_fMRI_mask, Mean_Image)

        command = (sing_afni + '3dcalc' + overwrite + ' -a ' + Prepro_fMRI_mask +
                   ' -prefix ' + Prepro_fMRI_mask + ' -expr "a"')
        run_cmd.do(command, diary_file)

        dictionary = {"Sources": Prepro_fMRI_mask,
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
            Skullstrip_func.Skullstrip_func(MAIN_PATH, Method_mask_func, Mean_Image, Prepro_fMRI_mask, anat_res_func, maskDilatfunc, dir_prepro_raw_process,
                                                      overwrite, costAllin, type_of_transform,
                                                      aff_metric_ants, sing_afni, sing_fsl, sing_fs, sing_itk, diary_file)

    if extra_erode > 0:
        # 1. Éroder final_mask basé sur l'intensité de Mean_Image
        nl = 'Eroding final_mask based on Mean_Image intensity at edges'
        run_cmd.msg(nl, diary_file, 'OKGREEN')

        # Charger les images
        mask_img = nib.load(Prepro_fMRI_mask)
        mask_data = mask_img.get_fdata()
        mean_img_nifti = nib.load(Mean_Image)
        mean_data = mean_img_nifti.get_fdata()

        # Créer un masque érodé basé sur l'intensité
        # Détecter les contours du masque
        eroded_mask = ndimage.binary_erosion(mask_data, iterations=extra_erode)
        edge_mask = mask_data.astype(bool) & ~eroded_mask

        # Calculer l'intensité moyenne dans le masque et sur les bords
        mean_intensity_inside = np.mean(mean_data[eroded_mask])
        std_intensity_inside = np.std(mean_data[eroded_mask])

        # Seuil adaptatif : retirer les voxels de bord avec intensité < moyenne - 1*std
        threshold = mean_intensity_inside - 1.0 * std_intensity_inside

        # Nouveau masque : garder les voxels internes + voxels de bord avec intensité suffisante
        new_mask = eroded_mask.copy()
        edge_voxels_to_keep = edge_mask & (mean_data >= threshold)
        new_mask = new_mask | edge_voxels_to_keep

        # Sauvegarder le nouveau masque
        new_mask_img = nib.Nifti1Image(new_mask.astype(np.float32), mask_img.affine, mask_img.header)
        nib.save(new_mask_img, Prepro_fMRI_mask)

        nl = f'Mask eroded: threshold={threshold:.2f}, voxels removed from edges: {np.sum(edge_mask) - np.sum(edge_voxels_to_keep)}'
        run_cmd.msg(nl, diary_file, 'OKGREEN')

        dictionary = {"Sources": [Prepro_fMRI_mask, Mean_Image],
                      "Description": f'Intensity-based erosion at edges. Threshold: {threshold:.2f}',
                      "Command": 'Python erosion based on Mean_Image intensity'}
        json_object = json.dumps(dictionary, indent=3)
        with open(Prepro_fMRI_mask.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)

    if not ope(opj(bids_dir, 'QC')):
        os.mkdir(opj(bids_dir, 'QC'))
    if not ope(opj(bids_dir, 'QC', 'check_mask_fMRI')):
        os.mkdir(opj(bids_dir, 'QC', 'check_mask_fMRI'))

    ####plot the QC
    plot_QC_func.plot_qc(Mean_Image,
                         Prepro_fMRI_mask,
                         opj(bids_dir, 'QC', 'check_fmri_mask', root_RS_ref + '_check_mask_fMRI.png'))