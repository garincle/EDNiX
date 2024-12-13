import os
import subprocess
import datetime
import json
import pathlib

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
opi = os.path.isfile

spco = subprocess.check_output
spgo = subprocess.getoutput

from fonctions.extract_filename import extract_filename

def _itk_check_masks(dir_fMRI_Refth_RS_prepro1,s_bind,itk_sif,diary_file):
    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  Working on step ' + str(4) + '(function: _itk_check_masks).  ##'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    File1 = pathlib.Path(opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.nii.gz'))
    Orig = File1.stat().st_ctime


    if not opi(opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz')):
        COND = 0
    else:
        File2 = pathlib.Path(opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz'))
        COND = 1


    nl = 'WARNING: if you modify the mask, please save it as ' + str(opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz'))
    print(bcolors.WARNING + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    command = 'gnome-terminal -- singularity run ' + s_bind + itk_sif + 'itksnap -g ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_test.nii.gz') + \
              ' -s ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat_Allineate_in_func.nii.gz')
    subprocess.run(f"gnome-terminal --wait -- bash -c '{command}; exec bash'", shell=True)

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
            dictionary = {"Modif_from_Sources": [opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_test.nii.gz'),
                                                 opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.nii.gz')],
                          "Modif_Description": 'manual modification (itksnap).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.json'), "a") as outfile:
                outfile.write(json_object)

            nl = 'the file: ' + opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz') + 'has been manually modified'
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

    diary.write(f'\n')
    diary.close()