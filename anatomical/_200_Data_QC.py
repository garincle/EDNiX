#coef dice??
##see fMRI prep²


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

def _itk_check_masks(dir_prepro, masks_dir, ID, type_norm):

    print('if you modify the mask, please save it as ' + str(opj(masks_dir, ID + '_finalmask.nii.gz')))
    command = 'gnome-terminal -- itksnap -g ' + opj(dir_prepro, ID + '_acpc_cropped' + type_norm + '.nii.gz') + \
              ' -s ' +opj(masks_dir, 'brain_mask_in_anat_DC.nii.gz')
    subprocess.run(f"gnome-terminal --wait -- bash -c '{command}; exec bash'", shell=True)