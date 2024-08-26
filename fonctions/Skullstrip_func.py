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

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput


def Skullstrip_func(Method_mask_func, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, overwrite, costAllin, lower_cutoff, upper_cutoff, s_bind, afni_sif, fsl_sif, fs_sif, itk_sif):


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
                                   'maskDilat_nilearn.nii.gz') + ' -fill_holes -dilate_input 1'
        spco(command, shell=True)

    elif brain_skullstrip == '3dSkullStrip':
        command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
            '-input ' + input_for_msk + ' -blur_fwhm 2 -orig_vol -mask_vol -monkey'
        spco(command, shell=True)
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 2'
        spco(command, shell=True)

    elif brain_skullstrip == '3dSkullStrip_marmoset':
        command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
            '-input ' + input_for_msk + ' -blur_fwhm 2 -orig_vol -mask_vol -marmoset'
        spco(command, shell=True)
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 2'
        spco(command, shell=True)

    elif brain_skullstrip == 'Custum_1':
        loadimg = nib.load(input_for_msk).get_fdata()
        loadimgsort85 =  np.percentile(np.abs(loadimg)[np.abs(loadimg)>0], 40)
        mask_imag = nilearn.image.threshold_img(input_for_msk, loadimgsort85, cluster_threshold=10)
        mask_imag.to_filename(output_for_mask)
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 1'
        spco(command, shell=True)
        command = 'singularity run' + s_bind + afni_sif + '3dresample -master ' + input_for_msk + ' -prefix' + output_for_mask + ' -input ' + output_for_mask + ' -overwrite -bound_type SLAB'

    elif brain_skullstrip =='bet2':
        #####creat an approximate brain mask
        command = 'singularity run' + s_bind + fsl_sif + 'bet2 ' + input_for_msk + ' ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_for_mask.nii.gz') + \
        ' -f 0.70'
        spco([command], shell=True)
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_for_mask.nii.gz') + ' -expr "step(a)" -prefix ' + output_for_mask + ' -overwrite'
        spco(command, shell=True)

    elif brain_skullstrip =='bet2_high':
        print('brain_skullstrip: applying ' + brain_skullstrip + ' method')
        #####creat an approximate brain mask
        command = 'singularity run' + s_bind + fsl_sif + 'bet2 ' + input_for_msk + ' ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_for_mask.nii.gz') + \
        ' -f 0.20'
        spco([command], shell=True)
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_for_mask.nii.gz') + ' -expr "step(a)" -prefix ' + output_for_mask + ' -overwrite'
        spco(command, shell=True)

    return(output_for_mask)


