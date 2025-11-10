import os
import json
import pathlib

opj = os.path.join
opi = os.path.isfile

from Tools import run_cmd

def _itk_check_func_in_template(dir_prepro_template_postprocessed, dir_prepro_template_masks, dir_prepro_template_process, sing_itk,diary_file, sing_afni, overwrite):

    nl = '##  Working on step ' + str(100) + '(function: _itk_check_func_in_template).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    final_mask = opj(dir_prepro_template_postprocessed, 'all_runs_space-template-func_desc-fMRI_Mean_Image_SS.nii.gz')
    Prepro_fMRI_mask = opj(dir_prepro_template_masks, 'mask_ref.nii.gz')
    Mean_Image = opj(dir_prepro_template_process, 'BASE_SS_fMRI.nii.gz')

    File1 = pathlib.Path(Prepro_fMRI_mask)
    Orig = File1.stat().st_ctime

    if not opi(final_mask):
        COND = 0
    else:
        File2 = pathlib.Path(final_mask)
        COND = 1

    nl = 'WARNING: if you modify the mask, please save it as ' + str(final_mask)
    run_cmd.msg(nl, diary_file, 'WARNING')

    command = (sing_itk + 'itksnap -g ' + Mean_Image +
               ' -s ' + Prepro_fMRI_mask)
    run_cmd.wait(command, diary_file)

    if COND == 0:
        if opi(final_mask):
            dictionary = {"Sources": [Mean_Image,
                                      Prepro_fMRI_mask],
                          "Description": 'manual modification (itksnap).',
                          "Command": command, }
            json_object = json.dumps(dictionary, indent=3)
            with open(final_mask.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)
    else:
        New = File2.stat().st_mtime

        if New - Orig >0:
            dictionary = {"Modif_f rom_Sources": [Mean_Image,
                                                 Prepro_fMRI_mask],
                          "Modif_Description": 'manual modification (itksnap).', "Command": command, }
            json_object = json.dumps(dictionary, indent=3)
            with open(final_mask.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)

            nl = 'the file: ' + final_mask + 'has been manually modified'
            run_cmd.msg(nl, diary_file, 'OKGREEN')

            command = (sing_afni + '3dcalc' + overwrite + ' -a ' + final_mask +
                       ' -prefix ' + Prepro_fMRI_mask + ' -expr "a"')
            run_cmd.do(command, diary_file)

            dictionary = {"Sources": final_mask,
                          "Description": 'Copy.', "Command": command, }
            json_object = json.dumps(dictionary, indent=3)
            with open(Prepro_fMRI_mask.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)