import os
import json
import pathlib

opj = os.path.join
opi = os.path.isfile

from Tools import run_cmd

def _itk_check_masks(dir_fMRI_Refth_RS_prepro1,sing_itk,diary_file, sing_afni, overwrite):

    nl = '##  Working on step ' + str(4) + '(function: _itk_check_masks).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    File1 = pathlib.Path(opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.nii.gz'))
    Orig = File1.stat().st_ctime


    if not opi(opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz')):
        COND = 0
    else:
        File2 = pathlib.Path(opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz'))
        COND = 1

    nl = 'WARNING: if you modify the mask, please save it as ' + str(opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz'))
    run_cmd.msg(nl, diary_file, 'WARNING')

    command = (sing_itk + 'itksnap -g ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_test.nii.gz') +
               ' -s ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat_Allineate_in_func.nii.gz'))
    run_cmd.wait(command, diary_file)

    if COND == 0:
        if opi(opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz')):
            dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_test.nii.gz'),
                                      opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.nii.gz')],
                          "Description": 'manual modification (itksnap).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.json'), "w") as outfile:
                outfile.write(json_object)
    else:
        New = File2.stat().st_mtime

        if New - Orig >0:
            dictionary = {"Modif_f rom_Sources": [opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_test.nii.gz'),
                                                 opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.nii.gz')],
                          "Modif_Description": 'manual modification (itksnap).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.json'), "a") as outfile:
                outfile.write(json_object)

            nl = 'the file: ' + opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz') + 'has been manually modified'
            run_cmd.msg(nl, diary_file, 'OKGREEN')

            command = (sing_afni + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz') +
                       ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.nii.gz') + ' -expr "a"')
            run_cmd.do(command, diary_file)

            dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz'),
                          "Description": 'Copy.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.json'), "w") as outfile:
                outfile.write(json_object)