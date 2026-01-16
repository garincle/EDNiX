import os
import json
import pathlib

opj = os.path.join
opi = os.path.isfile

from Tools import run_cmd

def _itk_check_func_in_template(dir_prepro_template_postprocessed, dir_prepro_template_masks, dir_prepro_template_process, sing_itk,diary_file, sing_afni, overwrite):

    nl = '##  Working on step ' + str(100) + '(function: _itk_check_func_in_template).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    Prepro_fMRI_mask = opj(dir_prepro_template_masks, 'Gmask.nii.gz')
    Mean_Image = opj(dir_prepro_template_process, 'BASE_SS_fMRI.nii.gz')
    Mean_Image_template = opj(dir_prepro_template_postprocessed, 'all_runs_space-template-func_desc-fMRI_Mean_Image_SS.nii.gz')

    command = (sing_itk + 'itksnap -g ' + Mean_Image +
               ' -o ' + Mean_Image_template)
    run_cmd.wait(command, diary_file)
