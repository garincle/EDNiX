import os
import subprocess

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

from fonctions.extract_filename import extract_filename

def _itk_check_masks(dir_fMRI_Refth_RS_prepro1,s_bind,itk_sif):

    print(bcolors.WARNING + 'WARNING: if you modify the mask, please save it as ' + str(opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz')) + bcolors.ENDC)
    command = 'gnome-terminal -- singularity run ' + s_bind + itk_sif + 'itksnap -g ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_test.nii.gz') + \
              ' -s ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat_Allineate_in_func.nii.gz')
    subprocess.run(f"gnome-terminal --wait -- bash -c '{command}; exec bash'", shell=True)