import os
import json
import pathlib

opj = os.path.join
opi = os.path.isfile

from Tools import run_cmd

def _itk_check_func_in_template(TfMRI, dir_prepro_raw_process, dir_prepro_raw_masks, sing_itk,diary_file):

    nl = '##  Working on step ' + str(100) + '(function: _itk_check_func_in_template).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    Mean_Image = opj(dir_prepro_raw_process, 'all_runs_space-func_desc-fMRI_Mean_Image_SS.nii.gz')
    Mean_Image_template = opj(dir_prepro_raw_process, ('_').join(['anat_space-raw-func', TfMRI + '.nii.gz']))
    mask = opj(dir_prepro_raw_masks,'mask_ref.nii.gz')

    command = (sing_itk + 'itksnap -g ' + Mean_Image +
               ' -o ' + Mean_Image_template +
               ' -s ' + mask)
    run_cmd.wait(command, diary_file)
