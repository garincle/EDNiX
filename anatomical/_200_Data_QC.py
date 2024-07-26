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

def _itk_check_masks(dir_prepro, masks_dir, ID, type_norm,s_bind,itk_sif):

    def run_command_and_wait(command):
        print("Running command:", command)
        result = subprocess.run(command, shell=True)
        if result.returncode == 0:
            print("Command completed successfully.")
        else:
            print("Command failed with return code:", result.returncode)

    # Example usage
    command = ('singularity run' + s_bind + itk_sif + 'itksnap -g ' + opj(dir_prepro, ID + '_acpc_cropped' + type_norm + '.nii.gz') + \
              ' -s ' + opj(masks_dir, 'brain_mask_in_anat_DC.nii.gz'))
    run_command_and_wait(command)

    # Continue with the rest of your script
    print("Continuing with the rest of the script...")


