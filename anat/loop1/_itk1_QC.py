import os
from Tools import getpath
import subprocess
import datetime
from Tools import run_cmd
import shutil
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


def _itk_check_masks(data_path, ID, type_norm, itk_sif,diary_file):
    ct = datetime.datetime.now()
    nl = 'Run anat._itk1_Data_QC.clean._itk_check_masks'
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{nl}')

    _,  dir_transfo,_, dir_prepro, _, volumes_dir, _, masks_dir = getpath.anat(data_path,
                                                                     '', '', False, False, 'native')
    anat_input1 = opj(dir_prepro, ID + '_space-raw_desc-n4Bias_')
    output4mask  = opj(masks_dir, ID + '_desc-step1_mask.nii.gz')

    end_maskname = '_'.join([ID, 'final', 'mask.nii.gz'])
    input4msk = anat_input1 + type_norm + '.nii.gz'

    def run_command_and_wait(command):
        print(bcolors.OKGREEN + "INFO: Running command:" + bcolors.ENDC, command)
        result = subprocess.run(command, shell=True)
        if result.returncode == 0:
            print(bcolors.OKGREEN + "INFO: Command completed successfully." + bcolors.ENDC)
        else:
            print(bcolors.OKGREEN + "INFO: Command failed with return code:" + bcolors.ENDC, result.returncode)

    nl = ('INFO: If you can not find a good solution for Skullstriping due to bad image quality, you can always modify it by hands and save it as: ' +
          opj(masks_dir,end_maskname))
    run_cmd.msg(nl, diary_file,'OKGREEN')

    command = (itk_sif + 'itksnap -g ' + input4msk +
               ' -s ' + output4mask)
    run_command_and_wait(command)

    if ope(opj(masks_dir,end_maskname)):
        nl = 'WARNING: We found an already existing final mask !!! it will be used! but restart step 2 to re-skullstrip your image !'
        run_cmd.msg(nl, diary_file, 'WARNING')
        output4mask = opj(masks_dir, ID + '_desc-step1_mask.nii.gz')
        shutil.copyfile(opj(masks_dir,end_maskname), output4mask)

    diary.write(f'\n')
    diary.close()
