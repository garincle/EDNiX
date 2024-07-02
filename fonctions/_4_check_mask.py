import os
import subprocess


#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

from fonctions.extract_filename import extract_filename

def _itk_check_masks(dir_fMRI_Refth_RS_prepro1):

    print('if you modify the mask, please save it as ' + str(opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz')))
    command = 'gnome-terminal -- itksnap -g ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_test.nii.gz') + \
              ' -s ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat_Allineate_in_func.nii.gz')
    subprocess.run(f"gnome-terminal --wait -- bash -c '{command}; exec bash'", shell=True)