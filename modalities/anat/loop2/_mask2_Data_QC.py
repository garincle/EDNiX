import os
import shutil
import subprocess
import datetime
import numpy as np
import nibabel as nib

from Tools import run_cmd

opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile


class bcolors:
    HEADER    = '\033[95m'
    OKBLUE    = '\033[94m'
    OKCYAN    = '\033[96m'
    OKGREEN   = '\033[92m'
    WARNING   = '\033[93m'
    FAIL      = '\033[91m'
    ENDC      = '\033[0m'
    BOLD      = '\033[1m'
    UNDERLINE = '\033[4m'


def _dice(mask_before_path, mask_after_path, diary_file):
    """
    Compute Dice similarity coefficient between two binary masks.

    Dice = 2|A∩B| / (|A|+|B|)

    A value of 1.0 means identical masks. Below ~0.95 typically indicates
    a substantial manual correction for a skull-strip mask.

    Returns float in [0,1] or None if computation failed.
    """
    try:
        a = np.rint(nib.load(mask_before_path).get_fdata()).astype(bool)
        b = np.rint(nib.load(mask_after_path).get_fdata()).astype(bool)
        union = np.sum(a) + np.sum(b)
        if union == 0:
            run_cmd.msg('WARNING: both masks empty — Dice undefined.',
                        diary_file, 'WARNING')
            return None
        return float(2.0 * np.sum(a & b) / union)
    except Exception as exc:
        run_cmd.msg(f'WARNING: Dice computation failed: {exc}',
                    diary_file, 'WARNING')
        return None


def _itk_check_masks(output4mask, input4msk, end_maskname, masks_dir,
                     itk_sif, diary_file):
    """
    Open ITK-SNAP for manual inspection and editing of the step-2 skull-strip
    mask (itk_2 step).

    The ITK-SNAP session shows:
      - Greyscale:    input4msk    (space-acpc_desc-SS-step1 anatomy)
      - Segmentation: output4mask  (current step-2 mask)

    If the user saves a corrected mask it should be saved as:
      {masks_dir}/{end_maskname}  (e.g. {ID}_final_mask_2.nii.gz)

    If that file is created or updated, _5_create_template_brain should be
    rerun by the caller to propagate the correction.

    Parameters
    ----------
    output4mask : str
        Path to the current step-2 mask shown as segmentation layer.
    input4msk : str
        Path to the anatomy shown as the greyscale layer.
    end_maskname : str
        Filename (not full path) of the final manual mask the user saves to.
    masks_dir : str
        Directory where end_maskname lives.
    itk_sif : str
        ITK-SNAP Singularity prefix.
    diary_file : str
        Diary log path.

    Returns
    -------
    mask_was_updated : bool
        True if the user created or updated the final mask during this session.
    """
    ct = datetime.datetime.now()
    nl = 'Run anat._mask2_Data_QC._itk_check_masks (itk_2)'
    with open(diary_file, 'a') as diary:
        diary.write(f'\n{ct}')
        diary.write(f'\n{nl}\n')

    final_mask = opj(masks_dir, end_maskname)

    run_cmd.msg(
        'INFO: If you cannot find a good skullstrip solution, modify the '
        'mask by hand and save it as: ' + final_mask,
        diary_file, 'OKGREEN')

    # ---- Snapshot mtime BEFORE opening ITK-SNAP ----------------------------
    mtime_before = os.path.getmtime(final_mask) if opi(final_mask) else None

    # Save a backup for Dice comparison (only if mask already exists)
    backup_path = None
    if mtime_before is not None:
        backup_path = final_mask.replace('.nii.gz', '_pre_itk2_backup.nii.gz')
        try:
            shutil.copy2(final_mask, backup_path)
        except Exception as exc:
            run_cmd.msg(f'WARNING: could not save mask backup for Dice: {exc}',
                        diary_file, 'WARNING')
            backup_path = None

    # ---- Open ITK-SNAP (blocking) ------------------------------------------
    def _run_and_wait(command):
        print(bcolors.OKGREEN + 'INFO: Running command:' + bcolors.ENDC,
              command)
        result = subprocess.run(command, shell=True)
        if result.returncode == 0:
            print(bcolors.OKGREEN + 'INFO: Completed successfully.'
                  + bcolors.ENDC)
        else:
            print(bcolors.WARNING + 'WARNING: Failed with return code: '
                  + bcolors.ENDC + str(result.returncode))

    command = itk_sif + 'itksnap -g ' + input4msk + ' -s ' + output4mask
    _run_and_wait(command)

    # ---- Detect whether final_mask was created or updated ------------------
    mask_was_updated = False

    if opi(final_mask):
        mtime_after = os.path.getmtime(final_mask)

        if mtime_before is None:
            mask_was_updated = True
            run_cmd.msg(
                'INFO: new final mask saved — _5_create_template_brain '
                'will be re-run.',
                diary_file, 'OKGREEN')

        elif mtime_after > mtime_before:
            mask_was_updated = True
            before_str = datetime.datetime.fromtimestamp(
                mtime_before).strftime('%Y-%m-%d %H:%M:%S')
            after_str = datetime.datetime.fromtimestamp(
                mtime_after).strftime('%Y-%m-%d %H:%M:%S')
            run_cmd.msg(
                f'INFO: final mask updated ({before_str} → {after_str}) — '
                f'_5_create_template_brain will be re-run.',
                diary_file, 'OKGREEN')

            # ---- Dice coefficient ------------------------------------------
            if backup_path is not None and opi(backup_path):
                dice = _dice(backup_path, final_mask, diary_file)
                if dice is not None:
                    level = ('minor'    if dice >= 0.98 else
                             'moderate' if dice >= 0.90 else
                             'major')
                    run_cmd.msg(
                        f'INFO: Dice (old vs new mask) = {dice:.4f} '
                        f'— {level} correction.',
                        diary_file, 'OKGREEN')
                try:
                    os.remove(backup_path)
                except Exception:
                    pass

        else:
            run_cmd.msg(
                'INFO: final mask unchanged — _5_create_template_brain '
                'will not be re-run.',
                diary_file, 'OKGREEN')

    else:
        run_cmd.msg(
            'INFO: no final mask saved — _5_create_template_brain '
            'will not be re-run.',
            diary_file, 'OKGREEN')

    return mask_was_updated
