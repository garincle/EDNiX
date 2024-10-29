###################################
###      Skullstrip method      ###
###################################
import subprocess
import os
import sys
import shutil
import nibabel as nib
import ants
import math
import nilearn
from nilearn.image import math_img
import numpy as np
from nilearn.masking import compute_epi_mask
from anatomical import Histrogram_mask_EMB
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


def Skullstrip_func(Method_mask_func, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, overwrite, costAllin, lower_cutoff, upper_cutoff, type_of_transform, aff_metric_ants, s_bind, afni_sif, fsl_sif, fs_sif, itk_sif):

    input_for_msk = opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz')
    output_for_mask = opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.nii.gz')
    master = opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.nii.gz')
    brain_skullstrip = Method_mask_func

    if brain_skullstrip == "3dAllineate":
        ##### mask the func img
        command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -cmass -EPI -final NN -float -twobest 5 -fineblur 0 -nomask -base ' + \
                  master + ' -prefix ' + opj(
            dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_for_mask.nii.gz') + \
                  ' -source ' + opj(dir_fMRI_Refth_RS_prepro1,
                                    'Mean_Image.nii.gz') + ' -' + costAllin + ' -1Dmatrix_save ' + \
                  opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_for_mask.1D') + \
                  ' -master ' + master
        spco(command, shell=True)

        command = 'singularity run' + s_bind + afni_sif + 'cat_matvec ' + opj(dir_fMRI_Refth_RS_prepro1,
                                                                              'Mean_Image_RcT_for_mask.1D') + \
                  ' -I | tail -n +3 > ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_for_mask_INV.1D')
        spco([command], shell=True)

        command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -final NN -1Dmatrix_apply ' + opj(
            dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_for_mask_INV.1D') + \
                  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.nii.gz') + \
                  ' -master ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz') + \
                  ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2, 'maskDilat.nii.gz')
        spco([command], shell=True)

    elif brain_skullstrip == "nilearn":
        # convert to float
        command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + opj(dir_fMRI_Refth_RS_prepro1,
                                                                                        'Mean_Image.nii.gz') + ' ' + opj(
            dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz')[:-7] + '_float.nii.gz'
        spco([command], shell=True)
        mask_img = compute_epi_mask(opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz')[:-7] + '_float.nii.gz',
                                    lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff,
                                    connected=True, opening=1,
                                    exclude_zeros=False, ensure_finite=True)
        mask_img.to_filename(opj(dir_fMRI_Refth_RS_prepro2, 'maskDilat_nilearn.nii.gz'))
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + opj(
            dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.nii.gz') + \
                  ' -input ' + opj(dir_fMRI_Refth_RS_prepro2,
                                   'maskDilat_nilearn.nii.gz') + ' -fill_holes -dilate_input -1 1'
        spco(command, shell=True)

    elif brain_skullstrip == '3dSkullStrip':
        command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
            '-input ' + input_for_msk + ' -blur_fwhm 2 -orig_vol -mask_vol'
        spco(command, shell=True)
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input -1 1'
        spco(command, shell=True)

    elif brain_skullstrip == '3dSkullStrip_monkey':
        command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
            '-input ' + input_for_msk + ' -blur_fwhm 2 -orig_vol -mask_vol -monkey'
        spco(command, shell=True)
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 2'
        spco(command, shell=True)

    elif brain_skullstrip == '3dSkullStrip_monkeynodil':
        command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
            '-input ' + input_for_msk + ' -blur_fwhm 2 -orig_vol -mask_vol -monkey'
        spco(command, shell=True)
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes'
        spco(command, shell=True)

    elif brain_skullstrip == '3dSkullStrip_marmoset':
        command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
            '-input ' + input_for_msk + ' -blur_fwhm 1 -orig_vol -mask_vol -marmoset'
        spco(command, shell=True)
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 1'
        spco(command, shell=True)

    elif brain_skullstrip =='bet2':
        #####creat an approximate brain mask
        command = 'singularity run' + s_bind + fsl_sif + 'bet2 ' + input_for_msk + ' ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_for_mask.nii.gz') + \
        ' -f 0.70'
        spco([command], shell=True)
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_for_mask.nii.gz') + ' -expr "step(a)" -prefix ' + output_for_mask + ' -overwrite'
        spco(command, shell=True)

    elif brain_skullstrip.startswith('_bet'):
        # Extract the last two digits to use as the -f value
        f_value = brain_skullstrip[-4:]
        # Create the approximate brain mask using bet2
        command = f'singularity run {s_bind}{fsl_sif} bet2 {input_for_msk} ' \
                  f'{opj(dir_fMRI_Refth_RS_prepro1, "Mean_Image_RcT_for_mask.nii.gz")} -f {f_value}'
        spco([command], shell=True)
        # Run the AFNI 3dcalc command to create the final mask
        command = f'singularity run {s_bind}{afni_sif} 3dcalc -a ' \
                  f'{opj(dir_fMRI_Refth_RS_prepro1, "Mean_Image_RcT_for_mask.nii.gz")} ' \
                  f'-expr "step(a)" -prefix {output_for_mask} -overwrite'
        spco(command, shell=True)

    elif brain_skullstrip.startswith('CustumNilearn_'):
        # Extract the cutoff values from the string
        _, lower_cutoff, upper_cutoff = brain_skullstrip.split('_')
        lower_cutoff = float(lower_cutoff)
        upper_cutoff = float(upper_cutoff)

        # Convert to float
        command = f'singularity run {s_bind}{fs_sif} mri_convert -odt float {input_for_msk} {opj(dir_fMRI_Refth_RS_prepro1, "Mean_Image_RcT_for_mask.nii.gz")}'
        spco([command], shell=True)
        # Compute the EPI mask using nilearn with the given cutoff values
        mask_img = compute_epi_mask(f'{opj(dir_fMRI_Refth_RS_prepro1, "Mean_Image_RcT_for_mask.nii.gz")}', lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff, connected=True, opening=3,
                                    exclude_zeros=False, ensure_finite=True)
        mask_img.to_filename(output_for_mask)
        # Use AFNI to process the mask
        command = f'singularity run {s_bind}{afni_sif} 3dmask_tool -overwrite -prefix {output_for_mask} -input {output_for_mask} -fill_holes -dilate_input 1'
        spco(command, shell=True)
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
        spco(command, shell=True)

    elif brain_skullstrip.startswith('CustumNilearnExcludeZeros_'):
        # Extract the cutoff values from the string
        _, lower_cutoff, upper_cutoff = brain_skullstrip.split('_')
        lower_cutoff = float(lower_cutoff)
        upper_cutoff = float(upper_cutoff)
        # Convert to float
        command = f'singularity run {s_bind}{fs_sif} mri_convert -odt float {input_for_msk} {opj(dir_fMRI_Refth_RS_prepro1, "Mean_Image_RcT_for_mask.nii.gz")}'
        spco([command], shell=True)

        # Compute the EPI mask using nilearn with the given cutoff values
        mask_img = compute_epi_mask(f'{opj(dir_fMRI_Refth_RS_prepro1, "Mean_Image_RcT_for_mask.nii.gz")}', lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff, connected=True, opening=3,
                                    exclude_zeros=True, ensure_finite=True)
        mask_img.to_filename(output_for_mask)
        # Use AFNI to process the mask
        command = f'singularity run {s_bind}{afni_sif} 3dmask_tool -overwrite -prefix {output_for_mask} -input {output_for_mask} -fill_holes -dilate_input 1'
        spco(command, shell=True)
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
        spco(command, shell=True)

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
        spco(command, shell=True)
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
        spco(command, shell=True)

    elif brain_skullstrip.startswith('Vol_sammba_'):
        volume = int(brain_skullstrip.split('_')[2])

        command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
        spco([command], shell=True)
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

    elif brain_skullstrip == 'Custum_ANTS_NL':
        IMG = ants.image_read(input_for_msk)
        REF_IMG = ants.image_read(master)
        mtx1 = ants.registration(fixed=IMG, moving=REF_IMG, type_of_transform='Translation',
                                 outprefix=opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_for_mask_'))
        REF_MASK = ants.image_read(opj(dir_fMRI_Refth_RS_prepro2, 'maskDilat.nii.gz'))
        tmp_mask1 = ants.apply_transforms(fixed=IMG, moving=REF_MASK,
                                          transformlist=mtx1['fwdtransforms'], interpolator='nearestNeighbor')
        tmp_mask1 = ants.threshold_image(tmp_mask1, 0.5, 1, 1, 0, True)
        tmp_mask1 = ants.iMath(tmp_mask1, operation='GetLargestComponent')

    elif brain_skullstrip == 'Custum_ANTS_Garin':
        IMG = ants.image_read(input_for_msk)
        REF_IMG = ants.image_read(master)

        mtx1 = ants.registration(fixed=IMG, moving=REF_IMG, type_of_transform='Translation',
                                 outprefix=opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_for_mask_shift'))
        MEAN_tr = ants.apply_transforms(fixed=IMG, moving=REF_IMG, transformlist=mtx1['fwdtransforms'],
                                        interpolator='nearestNeighbor')
        ants.image_write(MEAN_tr, opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_for_mask_shift.nii.gz'),
                         ri=False)
        mTx = ants.registration(fixed=IMG, moving=REF_IMG,
                                outprefix=opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_for_mask_'),
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

        transfo_concat = \
            [opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_for_mask_1Warp.nii.gz'),
             opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_for_mask_0GenericAffine.mat')]

        REF_MASK = ants.image_read(opj(dir_fMRI_Refth_RS_prepro2, 'maskDilat.nii.gz'))
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
        spco(command, shell=True)
        command = f'singularity run {s_bind}{afni_sif} 3dresample -master {input_for_msk} -prefix {output_for_mask} -input {output_for_mask} -overwrite -bound_type SLAB'
        spco(command, shell=True)

    elif brain_skullstrip == 'NoSkullStrip':
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + input_for_msk + \
                  ' -expr "step(a)" -prefix ' + output_for_mask + ' -overwrite'
        spco(command, shell=True)

    elif brain_skullstrip == 'Manual':
        def run_command_and_wait(command):
            print(bcolors.OKGREEN + 'INFO: Running command:', command + bcolors.ENDC)
            result = subprocess.run(command, shell=True)
            if result.returncode == 0:
                print(bcolors.OKGREEN + 'INFO: Command completed successfully.' + bcolors.ENDC)
            else:
                print(bcolors.WARNING + 'WARNING: Command failed with return code:', result.returncode, bcolors.ENDC)

        if not os.path.exists(output_for_mask):
            command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + input_for_msk + ' -expr "step(a)" -prefix ' + output_for_mask
            spco(command, shell=True)
        command = ('singularity run' + s_bind + itk_sif + 'itksnap -g ' + input_for_msk + ' -s ' + output_for_mask)
        run_command_and_wait(command)

    else:
        raise Exception(bcolors.FAIL + "ERROR: brain_skullstrip not recognized, check that brain_skullstrip_1 or brain_skullstrip_2 are correctly written!!" + bcolors.ENDC)


    return(output_for_mask)


