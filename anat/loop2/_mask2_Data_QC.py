#coef dice??
##see fMRI prep²
import os
import subprocess
import datetime
from Tools import run_cmd

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


def _itk_check_masks(output4mask, input4msk, end_maskname, masks_dir, itk_sif, diary_file):
    ct = datetime.datetime.now()
    nl = 'Run anat._mask2_Data_QC.clean._itk_check_masks'
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{nl}')


    nl = ('INFO: If you can not find a good solution for Skullstriping due to bad image quality, you can always modify it by hands and save it as: ' +
          opj(masks_dir,end_maskname))
    run_cmd.msg(nl, diary_file,'OKGREEN')


    def run_command_and_wait(command):
        print(bcolors.OKGREEN + "INFO: Running command:" + bcolors.ENDC, command)
        result = subprocess.run(command, shell=True)
        if result.returncode == 0:
            print(bcolors.OKGREEN + "INFO: Command completed successfully." + bcolors.ENDC)
        else:
            print(bcolors.OKGREEN + "INFO: Command failed with return code:" + bcolors.ENDC, result.returncode)

    command = (itk_sif + 'itksnap -g ' + input4msk +
               ' -s ' + output4mask)
    run_command_and_wait(command)


    diary.write(f'\n')
    diary.close()

