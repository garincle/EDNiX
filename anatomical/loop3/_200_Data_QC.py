#coef dice??
##see fMRI prep²
import os
import subprocess
import datetime


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


def _itk_check_masks(dir_prepro, masks_dir, ID, type_norm,s_bind,itk_sif,diary_file):
    ct = datetime.datetime.now()
    nl = 'Run anatomical._200_Data_QC.clean._itk_check_masks'
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{nl}')

    def run_command_and_wait(command):
        print(bcolors.OKGREEN + "INFO: Running command:" + bcolors.ENDC, command)
        result = subprocess.run(command, shell=True)
        if result.returncode == 0:
            print(bcolors.OKGREEN + "INFO: Command completed successfully." + bcolors.ENDC)
        else:
            print(bcolors.OKGREEN + "INFO: Command failed with return code:" + bcolors.ENDC, result.returncode)

    # Example usage
    command = ('singularity run' + s_bind + itk_sif + 'itksnap -g ' + opj(dir_prepro, ID + '_acpc_cropped' + type_norm + '.nii.gz') + \
              ' -s ' + opj(masks_dir, 'brain_mask_in_anat_DC.nii.gz'))
    run_command_and_wait(command)

    diary.write(f'\n')
    diary.close()

