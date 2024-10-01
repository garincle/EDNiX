#################################################
########    creat brain image of animal  ########
#################################################
#co-register linear to the new sty template...
##############################################################################
####### CREATE THE STUDY TEMPLATE (IF YOU WANT ON) ###########################
##############################################################################
import os
import nibabel as nib
import subprocess
from nilearn import plotting
import anatomical.Skullstrip_method
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

def skullstrip_T(BASE_SS, BASE_mask, dir_prepro, ID, Session,
                 dir_transfo, type_norm, volumes_dir, BASE_SS_coregistr, BASE_SS_mask, type_of_transform, aff_metric_ants,
                 study_template_atlas_forlder, otheranat, template_skullstrip,
                 masking_img, brain_skullstrip_1, brain_skullstrip_2, masks_dir,
                 check_visualy_final_mask, s_bind, afni_sif, fsl_sif, fs_sif, itk_sif, strip_sif, overwrite):

    warp_adj = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'warped_3_adjusted_mean.nii.gz')
    stdy_template = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template.nii.gz')

    step_skullstrip = 3
    stdy_template_mask = anatomical.Skullstrip_method.Skullstrip_method(step_skullstrip, template_skullstrip, study_template_atlas_forlder, masking_img, brain_skullstrip_1, brain_skullstrip_2, masks_dir, volumes_dir, dir_prepro, type_norm, BASE_SS_coregistr, BASE_SS_mask,
    type_of_transform, ID, aff_metric_ants, check_visualy_final_mask, s_bind, afni_sif, fsl_sif, fs_sif, itk_sif, strip_sif)

    ##extract brain
    command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + stdy_template_mask + ' -b ' + warp_adj + ' -expr "(a*b)" -prefix ' + opj \
        (study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template_not_align.nii.gz')
    spco(command, shell=True)

    ##align template to BASE_SS (atlas template)
    command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -warp shift_rotate -cmass -overwrite -base ' + BASE_SS + \
              ' -nomask -onepass -final NN' + \
              ' -master ' + opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template_not_align.nii.gz') + \
              ' -prefix ' + stdy_template + \
              ' -source ' + opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template_not_align.nii.gz') + ' -1Dmatrix_save ' + \
              opj(study_template_atlas_forlder , 'studytemplate2_' + type_norm, 'align_template_to_stdy_template.1D')
    spco(command, shell=True)
    '''
    current_working_directory = os.getcwd()
    os.chdir(study_template_atlas_forlder + '/studytemplate2_' + type_norm)
    command = 'singularity run' + s_bind + afni_sif + '@Align_Centers -base ' + BASE_SS + \
    ' -dset ' + opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template_not_align.nii.gz') + ' -cm -prefix ' + 'align_template_to_stdy_template' + overwrite
    spco([command], shell=True)
    os.chdir(str(current_working_directory))
    '''
    command = 'singularity run' + s_bind + afni_sif + '3dAllineate -overwrite -final NN -1Dmatrix_apply ' + \
        opj(study_template_atlas_forlder , 'studytemplate2_' + type_norm, 'align_template_to_stdy_template.1D') + \
              ' -prefix ' + stdy_template_mask + \
              ' -master ' + stdy_template + \
              ' -input  ' + stdy_template_mask
    spco([command], shell=True)