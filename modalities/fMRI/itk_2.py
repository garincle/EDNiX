import os
import json
import pathlib

opj = os.path.join
opi = os.path.isfile

from Tools import run_cmd

def _itk_check_func_in_template(TfMRI, dir_prepro_acpc_process, sing_itk,diary_file):

    nl = '##  Working on step ' + str(100) + '(function: _itk_check_func_in_template).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    Mean_Image = opj(dir_prepro_acpc_process, 'all_runs_space-anat_desc-fMRI_Mean_Image_unwarped.nii.gz')
    Mean_Image_template = opj(dir_prepro_acpc_process, 'anat_space-acpc_res-func_'+ TfMRI + '.nii.gz')

    command = (sing_itk + 'itksnap -g ' + Mean_Image +
               ' -o ' + Mean_Image_template)
    run_cmd.wait(command, diary_file)
