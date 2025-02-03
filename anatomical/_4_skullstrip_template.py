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
import json
import datetime

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
                 study_template_atlas_folder, otheranat, template_skullstrip,
                 masking_img, brain_skullstrip_1, brain_skullstrip_2, masks_dir,
                 check_visualy_final_mask, s_bind, afni_sif, fsl_sif, fs_sif, itk_sif, strip_sif, overwrite,diary_file):

    ct = datetime.datetime.now()
    nl = 'Run anatomical._4_skullstrip_template.skullstrip_T'
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{nl}')
    diary.write(f'\n')
    diary.close()

    warp_adj      = opj(study_template_atlas_folder, 'studytemplate2_' + type_norm, 'warped_3_adjusted_mean.nii.gz')
    stdy_template = opj(study_template_atlas_folder, 'studytemplate2_' + type_norm, 'study_template.nii.gz')

    step_skullstrip = 3
    stdy_template_mask = anatomical.Skullstrip_method.Skullstrip_method(step_skullstrip, template_skullstrip, study_template_atlas_folder, masking_img, brain_skullstrip_1, brain_skullstrip_2, masks_dir, volumes_dir, dir_prepro, type_norm, BASE_SS_coregistr, BASE_SS_mask,
    type_of_transform, ID, aff_metric_ants, check_visualy_final_mask, s_bind, afni_sif, fsl_sif, fs_sif, itk_sif, strip_sif,diary_file)

    ct = datetime.datetime.now()
    nl = 'Run anatomical._4_skullstrip_template.skullstrip_T, after Skullstrip_method'
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{nl}')

    ## extract brain
    command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + stdy_template_mask + ' -b ' + warp_adj + \
              ' -expr "(a*b)" -prefix ' + stdy_template
    nl = spgo(command)
    diary.write(f'\n{nl}')
    print(nl)
    dictionary = {"Sources": [stdy_template_mask,
                              warp_adj],
                  "Description": 'Skull stripping.', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(study_template_atlas_folder, 'studytemplate2_' + type_norm, 'study_template.json'), "w") as outfile:
        outfile.write(json_object)

    diary.close()