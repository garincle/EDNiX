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

def _itk_check_spatial_co(dir_fMRI_Refth_RS_prepro3,s_bind,itk_sif):

    def run_command_and_wait(command):
        print("Running command:", command)
        result = subprocess.run(command, shell=True)
        if result.returncode == 0:
            print("Command completed successfully.")
        else:
            print("Command failed with return code:", result.returncode)

    # Example usage
    command = ('singularity run' + s_bind + itk_sif + 'itksnap -g ' + opj(dir_fMRI_Refth_RS_prepro3,'Mean_Image_RcT_SS_in_template.nii.gz') + \
              ' -o ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz'))
    run_command_and_wait(command)

    # Continue with the rest of your script
    print("Continuing with the rest of the script...")
