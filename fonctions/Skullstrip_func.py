###################################
###      Skullstrip method      ###
###################################
import subprocess
import os
import nibabel as nib
import ants
import nilearn
import numpy as np
from nilearn.masking import compute_epi_mask
from anatomical.skullstrip import Histrogram_mask_EMB
from fonctions.extract_filename import extract_filename
import datetime
import json
import shutil


opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

from Tools import run_cmd

def Skullstrip_func(Method_mask_func, input_for_msk, output_for_mask, master, maskDilatfunc, dir_prepro_orig_process,
                                                      overwrite, costAllin, type_of_transform,
                                                      aff_metric_ants, s_bind, afni_sif, fsl_sif, sing_fs, itk_sif,diary_file):

    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  Working on step ' + str(3) + '(function: Skullstrip_func).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    Mean_Image_RcT_for_mask = opj(dir_prepro_orig_process, 'Mean_Image_RcT_for_mask.nii.gz')
    Mean_Image_RcT_for_mask_1D = opj(dir_prepro_orig_process, 'Mean_Image_RcT_for_mask.1D')
    Mean_Image_RcT_for_mask_1D_INV = opj(dir_prepro_orig_process, 'Mean_Image_RcT_for_mask_INV.1D')
    Mean_Image_RcT_for_mask_desc = opj(dir_prepro_orig_process, 'Mean_Image_RcT_for_mask_')
    Mean_Image_RcT_for_mask_shift = opj(dir_prepro_orig_process, 'Mean_Image_RcT_for_mask_shift.nii.gz')
    transo = [opj(dir_prepro_orig_process, 'Mean_Image_RcT_for_mask_1Warp.nii.gz'),
    opj(dir_prepro_orig_process, 'Mean_Image_RcT_for_mask_0GenericAffine.mat')]
    brain_skullstrip = Method_mask_func

    if brain_skullstrip == "3dAllineate":
        ##### mask the func img
        command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -cmass -EPI -final NN -float -twobest 5 -fineblur 0 -nomask -base ' + \
                  master + ' -prefix ' + Mean_Image_RcT_for_mask + \
                  ' -source ' + input_for_msk + ' -' + costAllin + ' -1Dmatrix_save ' + \
                  Mean_Image_RcT_for_mask_1D + \
                  ' -master ' + master
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": [input_for_msk,
                                  master],
                      "Description": 'Co-registration (3dAllineate,AFNI).',
                      "Command": command, }
        json_object = json.dumps(dictionary, indent=3)
        with open(Mean_Image_RcT_for_mask.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)          

        command = 'singularity run' + s_bind + afni_sif + 'cat_matvec ' + Mean_Image_RcT_for_mask_1D + \
                  ' -I | tail -n +3 > ' + Mean_Image_RcT_for_mask_1D_INV
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -final NN -1Dmatrix_apply ' + Mean_Image_RcT_for_mask_1D_INV + \
                  ' -prefix ' + output_for_mask + \
                  ' -master ' + input_for_msk + \
                  ' -input  ' + maskDilatfunc
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": [maskDilatfunc,
                                  input_for_msk,
                                  Mean_Image_RcT_for_mask_1D_INV],
                      "Description": 'Co-registration (3dAllineate,AFNI).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

    elif brain_skullstrip == '3dSkullStrip':
        command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
            '-input ' + input_for_msk + ' -blur_fwhm 2 -orig_vol -mask_vol'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input -1 1'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": input_for_msk,
                      "Description": ['binary brain mask (3dSkullStrip, AFNI)',
                                      'fill holes and dilation (3dmask_tool,AFNI).'], }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

    elif brain_skullstrip == '3dSkullStrip_monkey':

        command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
            '-input ' + input_for_msk + ' -blur_fwhm 2 -orig_vol -mask_vol -monkey'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 2'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": input_for_msk,
                      "Description": ['binary brain mask (3dSkullStrip, AFNI)',
                                      'fill holes and dilation (3dmask_tool,AFNI).'], }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)


    elif brain_skullstrip == '3dSkullStrip_monkeynodil':
        command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
            '-input ' + input_for_msk + ' -blur_fwhm 2 -orig_vol -mask_vol -monkey'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        dictionary = {"Sources": input_for_msk,
                      "Description": ['binary brain mask (3dSkullStrip, AFNI)',
                                      'fill holes  (3dmask_tool,AFNI).'], }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

    elif brain_skullstrip == '3dSkullStrip_dog':
        command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
            '-input ' + input_for_msk + ' -orig_vol -mask_vol -monkey -use_skull'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + output_for_mask + ' -expr "step(a-4)" -prefix ' + output_for_mask + ' -overwrite'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input -1 2'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        dictionary = {"Sources": input_for_msk,
                      "Description": ['binary brain mask (3dSkullStrip, AFNI)',
                                      'fill holes  (3dmask_tool,AFNI).'], }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

    elif brain_skullstrip == '3dSkullStrip_marmoset':
        command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
            '-input ' + input_for_msk + ' -blur_fwhm 1 -orig_vol -mask_vol -marmoset'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 1'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": input_for_msk,
                      "Description": ['binary brain mask (3dSkullStrip, AFNI)',
                                      'fill holes and dilation (3dmask_tool,AFNI).'], }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

    elif brain_skullstrip == 'muSkullStrip_cross_species':
        command = 'python3 ' + opj(opd(afni_sif), 'NHP-BrainExtraction', 'UNet_Model', 'muSkullStrip.py') + \
                  ' -in ' + input_for_msk + \
                  ' -model ' + opj(opd(afni_sif), 'NHP-BrainExtraction', 'UNet_Model', 'models',
                                   'model-02-_cross_species-epoch') + \
                  ' -out ' + opd(input_for_msk)
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        shutil.copyfile(opj(opd(output_for_mask), extract_filename(input_for_msk) + '_pre_mask.nii.gz'),
                        output_for_mask)
        command = f'singularity run {s_bind}{afni_sif} 3dmask_tool -overwrite -prefix {output_for_mask} -input {output_for_mask} -fill_holes -dilate_input'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": input_for_msk,
                      "Description": 'Brain mask (U-Net).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

    elif brain_skullstrip =='bet2':
        ##### create an approximate brain mask
        command = 'singularity run' + s_bind + fsl_sif + 'bet2 ' + input_for_msk + ' ' + Mean_Image_RcT_for_mask + \
        ' -f 0.70'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": input_for_msk,
                      "Description": 'skull stripping (bet2,FSL).',
                      "Command": command, }
        json_object = json.dumps(dictionary, indent=3)
        with open(Mean_Image_RcT_for_mask.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)

        command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + Mean_Image_RcT_for_mask + ' -expr "step(a)" -prefix ' + output_for_mask + ' -overwrite'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": Mean_Image_RcT_for_mask,
                      "Description": 'binary mask (3dcalc "step(x)", AFNI)', }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)


    elif brain_skullstrip.startswith('_bet'):
        # only difference with bet2 : the "-f 70" parameter

        # Extract the last two digits to use as the -f value
        f_value = brain_skullstrip[-4:]
        # Create the approximate brain mask using bet2
        command = f'singularity run {s_bind}{fsl_sif} bet2 {input_for_msk} ' \
                  f'{Mean_Image_RcT_for_mask} -f {f_value}'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": input_for_msk,
                      "Description": 'skull stripping (bet2,FSL).',
                      "Command": command, }
        json_object = json.dumps(dictionary, indent=3)
        with open(Mean_Image_RcT_for_mask.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)

        # Run the AFNI 3dcalc command to create the final mask
        command = f'singularity run {s_bind}{afni_sif} 3dcalc -a ' \
                  f'{Mean_Image_RcT_for_mask} ' \
                  f'-expr "step(a)" -prefix {output_for_mask} -overwrite'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": Mean_Image_RcT_for_mask,
                      "Description": 'binary mask (3dcalc "step(x)", AFNI)', }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

    elif brain_skullstrip.startswith('CustomNilearn_'):
        # Extract the cutoff values from the string
        _, lower_cutoff, upper_cutoff = brain_skullstrip.split('_')
        lower_cutoff = float(lower_cutoff)
        upper_cutoff = float(upper_cutoff)

        # Convert to float
        command = f'singularity run {s_bind}{sing_fs} mri_convert -odt float {input_for_msk} {Mean_Image_RcT_for_mask}'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        # Compute the EPI mask using nilearn with the given cutoff values
        mask_img = compute_epi_mask(f'{Mean_Image_RcT_for_mask}', lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff, connected=True, opening=3,
                                    exclude_zeros=False, ensure_finite=True)
        mask_img.to_filename(output_for_mask)


        # Use AFNI to process the mask
        command = f'singularity run {s_bind}{afni_sif} 3dmask_tool -overwrite -prefix {output_for_mask} -input {output_for_mask} -fill_holes -dilate_input 1'
        nl = spgo(command)
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
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": [input_for_msk,
                                  Mean_Image_RcT_for_mask],
                      "Description": ['convertion to float (mri_convert, Freesurfer)',
                                      'binary brain mask (compute_epi_mask,nilearn).',
                                      'fill holes and dilation (3dmask_tool, AFNI'],}
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)


    elif brain_skullstrip.startswith('CustomNilearnExcludeZeros_'):
        # Extract the cutoff values from the string
        _, lower_cutoff, upper_cutoff = brain_skullstrip.split('_')
        lower_cutoff = float(lower_cutoff)
        upper_cutoff = float(upper_cutoff)
        # Convert to float
        command = f'singularity run {s_bind}{sing_fs} mri_convert -odt float {input_for_msk} {Mean_Image_RcT_for_mask}'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        # Compute the EPI mask using nilearn with the given cutoff values
        mask_img = compute_epi_mask(f'{Mean_Image_RcT_for_mask}', lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff, connected=True, opening=3,
                                    exclude_zeros=True, ensure_finite=True)
        mask_img.to_filename(output_for_mask)
        # Use AFNI to process the mask
        command = f'singularity run {s_bind}{afni_sif} 3dmask_tool -overwrite -prefix {output_for_mask} -input {output_for_mask} -fill_holes -dilate_input 1'
        nl = spgo(command)
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
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": [input_for_msk,
                                  Mean_Image_RcT_for_mask],
                      "Description": ['convertion to float (mri_convert, Freesurfer)',
                                      'binary brain mask (compute_epi_mask,nilearn).',
                                      'fill holes and dilation (3dmask_tool, AFNI',
                                      'smooth and threshold (ANTSpy)'], }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

    # Handle 'Custom' type skullstripping dynamically for percentile
    elif brain_skullstrip.startswith('CustomThreshold_'):
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
        nl = spgo(command)
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
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": input_for_msk,
                      "Description": ['binary brain mask (threshold_img,nilearn).',
                                      'fill holes and dilation (3dmask_tool, AFNI',
                                      'smooth and threshold (ANTSpy)'], }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

    elif brain_skullstrip.startswith('Vol_sammba_'):
        volume = int(brain_skullstrip.split('_')[2])

        command = 'singularity run' + s_bind + sing_fs + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
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
                      "Description": ['convertion to float (mri_convert, Freesurfer)',
                                      'binary brain mask (HistogramMask,Sammba).'], }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

    elif brain_skullstrip == 'Custom_ANTS_NL':
        IMG = ants.image_read(input_for_msk)
        REF_IMG = ants.image_read(master)
        mtx1 = ants.registration(fixed=IMG, moving=REF_IMG, type_of_transform='Translation',
                                 outprefix=Mean_Image_RcT_for_mask_desc)
        REF_MASK = ants.image_read(maskDilatfunc)
        tmp_mask1 = ants.apply_transforms(fixed=IMG, moving=REF_MASK,
                                          transformlist=mtx1['fwdtransforms'], interpolator='nearestNeighbor')
        tmp_mask1 = ants.threshold_image(tmp_mask1, 0.5, 1, 1, 0, True)
        tmp_mask1 = ants.iMath(tmp_mask1, operation='GetLargestComponent')
        ants.image_write(tmp_mask1, output_for_mask, ri=False)
        dictionary = {"Sources": [input_for_msk,
                                  master,
                                  maskDilatfunc],
                      "Description": ['binary brain mask (ANTSpy).'], }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)


    elif brain_skullstrip == 'Custom_ANTS_Garin':
        IMG = ants.image_read(input_for_msk)
        REF_IMG = ants.image_read(master)

        mtx1 = ants.registration(fixed=IMG, moving=REF_IMG, type_of_transform='Translation',
                                 outprefix=Mean_Image_RcT_for_mask_shift)
        MEAN_tr = ants.apply_transforms(fixed=IMG, moving=REF_IMG, transformlist=mtx1['fwdtransforms'],
                                        interpolator='nearestNeighbor')
        ants.image_write(MEAN_tr, Mean_Image_RcT_for_mask_shift,
                         ri=False)
        mTx = ants.registration(fixed=IMG, moving=REF_IMG,
                                outprefix=Mean_Image_RcT_for_mask_desc,
                                initial_transform=mtx1['fwdtransforms'],
                                type_of_transform=type_of_transform,
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

        transfo_concat = transo

        REF_MASK = ants.image_read(maskDilatfunc)
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
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        command = f'singularity run {s_bind}{afni_sif} 3dresample -master {input_for_msk} -prefix {output_for_mask} -input {output_for_mask} -overwrite -bound_type SLAB'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": [input_for_msk,
                                  master,
                                  maskDilatfunc],
                      "Description": ['binary brain mask (ANTSpy).',
                                      'fill holes (3dmasktool,AFNI)',
                                      'resampling (3dresample, AFNI)'], }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)


    elif brain_skullstrip == 'NoSkullStrip':
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + input_for_msk + \
                  ' -expr "step(a)" -prefix ' + output_for_mask + ' -overwrite'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": input_for_msk,
                      "Description": 'binary image (3dcalc "step(x)", AFNI).',}
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

    elif brain_skullstrip == 'Manual':
        def run_command_and_wait(cmd):
            nl ='INFO: Running command:' +  cmd
            run_cmd.msg(nl, diary_file, 'OKGREEN')
            result = subprocess.run(cmd, shell=True)
            if result.returncode == 0:
                nl = 'INFO: Command completed successfully.'
                run_cmd.msg(nl, diary_file, 'OKGREEN')
            else:
                nl = 'WARNING: Command failed with return code: ' +  result.returncode
                run_cmd.msg(nl, diary_file, 'WARNING')

        if not ope(output_for_mask):
            command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + input_for_msk + ' -expr "step(a)" -prefix ' + output_for_mask
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

        command = ('singularity run' + s_bind + itk_sif + 'itksnap -g ' + input_for_msk + ' -s ' + output_for_mask)
        run_command_and_wait(command)
        nl = command
        diary.write(f'\n{nl}')

        dictionary = {"Sources": input_for_msk,
                      "Description": 'Manual binary brain image (itksnap).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

    else:
        nl = "ERROR: brain_skullstrip not recognized, check that Method_mask_func is correctly written!!"
        diary.write(f'\n{nl}')
        raise Exception(run_cmd.error(nl, diary_file))

    diary.write(f'\n')
    diary.close()

    return output_for_mask


