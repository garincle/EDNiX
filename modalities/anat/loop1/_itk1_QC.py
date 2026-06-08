import os
from Tools import getpath
import subprocess
import datetime
from Tools import run_cmd
import shutil
from Tools import check_nii
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

# FIX — snapshot mtime before, compare after, return bool
def _itk_check_masks(data_path, ID, type_norm, itk_sif, diary_file):
    ct = datetime.datetime.now()
    nl = 'Run anat._itk1_Data_QC.clean._itk_check_masks'

    diary = open(diary_file, 'a')
    diary.write(f'\n{ct}')
    diary.write(f'\n{nl}')

    _, dir_transfo, _, dir_prepro, _, volumes_dir, _, masks_dir = getpath.anat(
        data_path, '', '', False, False, 'native')

    anat_input1  = opj(dir_prepro,  ID + '_space-raw_desc-Bias_')
    output4mask  = opj(masks_dir,   ID + '_desc-step1_mask.nii.gz')
    end_maskname = '_'.join([ID, 'final', 'mask.nii.gz'])
    input4msk    = anat_input1 + type_norm + '.nii.gz'
    final_mask = opj(masks_dir, end_maskname)

    if os.path.exists(final_mask):
        nl = 'WARNING: We found an already existing final mask !!! No Skullstrip will be performed!'
        run_cmd.msg(nl, diary_file, 'WARNING')
        shutil.copyfile(final_mask, output4mask)
        check_nii.resamp(output4mask, input4msk, 'msk', '', '',
                         diary_file,
                         '')

    nl = ('INFO: If you cannot find a good skullstrip solution, you can always '
          'modify the mask by hand and save it as: ' + final_mask)
    run_cmd.msg(nl, diary_file, 'OKGREEN')

    # --- Snapshot mtime BEFORE opening ITK-SNAP ----------------------------
    mtime_before = os.path.getmtime(final_mask) if os.path.exists(final_mask) else None

    command = (itk_sif + 'itksnap -g ' + input4msk + ' -s ' + output4mask)

    def run_command_and_wait(command):
        print(bcolors.OKGREEN + 'INFO: Running command:' + bcolors.ENDC, command)
        result = subprocess.run(command, shell=True)
        if result.returncode == 0:
            print(bcolors.OKGREEN + 'INFO: Command completed successfully.' + bcolors.ENDC)
        else:
            print(bcolors.OKGREEN + 'INFO: Command failed with return code:'
                  + bcolors.ENDC, result.returncode)

    run_command_and_wait(command)

    # --- Check whether final_mask is new or was modified -------------------
    mask_was_updated = False
    if os.path.exists(final_mask):
        mtime_after = os.path.getmtime(final_mask)
        if mtime_before is None:
            # mask did not exist before → user created it
            mask_was_updated = True
            run_cmd.msg(
                'INFO: new final mask saved by user — '
                '_2_clean_anat will be re-run.',
                diary_file, 'OKGREEN')
        elif mtime_after > mtime_before:
            # mask existed but was overwritten
            mask_was_updated = True
            before_str = datetime.datetime.fromtimestamp(mtime_before).strftime('%Y-%m-%d %H:%M:%S')
            after_str  = datetime.datetime.fromtimestamp(mtime_after).strftime('%Y-%m-%d %H:%M:%S')
            run_cmd.msg(
                f'INFO: final mask was updated ({before_str} → {after_str}) — '
                f'_2_clean_anat will be re-run.',
                diary_file, 'OKGREEN')
        else:
            run_cmd.msg(
                'INFO: final mask unchanged — _2_clean_anat will not be re-run.',
                diary_file, 'OKGREEN')
    else:
        run_cmd.msg(
            'INFO: no final mask saved — _2_clean_anat will not be re-run.',
            diary_file, 'OKGREEN')

    diary.write(f'\n')
    diary.close()

    return mask_was_updated   # ← caller uses this to decide re-run