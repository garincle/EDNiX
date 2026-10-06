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
    UNDERLINE = '\033[4e'


def _backtonative_mask_path(masks_dir, ID):
    """
    Canonical path for the backtonative mask.

    Convention: {masks_dir}/{ID}_backtonative_mask.nii.gz

    This is the manually edited final mask produced via ITK-SNAP (step itk_3).
    When it exists it takes priority over the automatically generated ACPC mask
    for the final skull-stripped anatomy in backtonative.apply().

    This function is the single source of truth for this path — both
    _200_Data_QC and backtonative.py import and call this function so the
    naming is always consistent.
    """
    return opj(masks_dir, ID + '_backtonative_mask.nii.gz')


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


def _itk_check_masks(info, brain_mask, volumes_dir, masks_dir,
                     ID, type_norm, itk_sif, diary_file):
    """
    Open ITK-SNAP for manual inspection and editing of the final ACPC mask
    (itk_3 — the last manual correction step before surface generation).

    The ITK-SNAP session shows:
      - Greyscale:    space-acpc_desc-template_{type_norm}  (anatomical ref)
      - Overlay:      template-space SS anatomy             (registration check)
      - Segmentation: brain_mask                           (current ACPC mask)

    If the user wants to apply a correction they must save the edited mask as
    the backtonative mask:
        {masks_dir}/{ID}_backtonative_mask.nii.gz

    This file is the canonical manual override consumed by backtonative.apply()
    to produce the final skull-stripped anatomy (space-acpc_desc-SS_T1w).
    Because the mask feeds into the registration (step 6) and the final SS
    anatomy (step 8), if it is updated the caller (_loop3) must rerun:
        [0] regenerate SS T1w from new mask (backtonative.apply)
        [1] re-run step 6 (_6_stdyTmax.nativetoTemplate)
        [2] re-run step 8 (backtonative.apply with fresh transforms)

    Parameters
    ----------
    info : list
        Transformation package from backtonative.get().
    brain_mask : str
        Path to the current ACPC brain mask shown as the segmentation layer.
    volumes_dir : str
        Path to the volumes directory for this subject/session.
    masks_dir : str
        Path to the masks directory for this subject/session.
    ID : str
        Subject ID.
    type_norm : str
        Anatomical image type used as reference (e.g. 'T1w').
    itk_sif : str
        ITK-SNAP Singularity prefix.
    diary_file : str
        Diary log path.

    Returns
    -------
    mask_was_updated : bool
        True if the user created or updated the backtonative mask.
        The caller is responsible for triggering the re-registration cascade.
    """
    ct = datetime.datetime.now()
    nl = 'Run anat._200_Data_QC._itk_check_masks (itk_3)'
    with open(diary_file, 'a') as diary:
        diary.write(f'\n{ct}')
        diary.write(f'\n{nl}\n')

    # Canonical backtonative mask path
    backtonative_mask = _backtonative_mask_path(masks_dir, ID)

    run_cmd.msg(
        'INFO: ITK-SNAP (itk_3) — inspect the final ACPC mask.\n'
        'To apply a manual correction, save the edited mask as:\n  '
        + backtonative_mask + '\n'
        'If saved, EDNiX will automatically rerun step 6 (registration) '
        'and step 8 (backtonative) to propagate the correction.',
        diary_file, 'OKGREEN')

    # ---- Snapshot mtime BEFORE opening ITK-SNAP ----------------------------
    mtime_before = (os.path.getmtime(backtonative_mask)
                    if opi(backtonative_mask) else None)

    # Save a backup for Dice comparison (only if mask already exists)
    backup_path = None
    if mtime_before is not None:
        backup_path = backtonative_mask.replace('.nii.gz',
                                                '_pre_itk3_backup.nii.gz')
        try:
            shutil.copy2(backtonative_mask, backup_path)
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

    command = (
        itk_sif + 'itksnap'
        ' -g ' + opj(volumes_dir,
                     ID + '_space-acpc_desc-template_' + type_norm + '.nii.gz') +
        ' -o ' + opj(info[0][4],
                     '_'.join([info[0][0], 'space-acpc_desc-SS', type_norm])
                     + '.nii.gz') +
        ' -s ' + brain_mask)
    _run_and_wait(command)

    # ---- Detect whether backtonative mask was created or updated -----------
    mask_was_updated = False

    if opi(backtonative_mask):
        mtime_after = os.path.getmtime(backtonative_mask)

        if mtime_before is None:
            mask_was_updated = True
            run_cmd.msg(
                'INFO: backtonative mask created — re-registration cascade '
                '(steps 6 + 8) will be triggered by caller.',
                diary_file, 'OKGREEN')

        elif mtime_after > mtime_before:
            mask_was_updated = True
            before_str = datetime.datetime.fromtimestamp(
                mtime_before).strftime('%Y-%m-%d %H:%M:%S')
            after_str = datetime.datetime.fromtimestamp(
                mtime_after).strftime('%Y-%m-%d %H:%M:%S')
            run_cmd.msg(
                f'INFO: backtonative mask updated ({before_str} → {after_str})'
                f' — re-registration cascade (steps 6 + 8) will be triggered '
                f'by caller.',
                diary_file, 'OKGREEN')

            # ---- Dice coefficient ------------------------------------------
            if backup_path is not None and opi(backup_path):
                dice = _dice(backup_path, backtonative_mask, diary_file)
                if dice is not None:
                    level = ('minor'    if dice >= 0.98 else
                             'moderate' if dice >= 0.90 else
                             'major')
                    run_cmd.msg(
                        f'INFO: Dice (old vs new backtonative mask) = '
                        f'{dice:.4f} — {level} correction.',
                        diary_file, 'OKGREEN')
                try:
                    os.remove(backup_path)
                except Exception:
                    pass

        else:
            run_cmd.msg(
                'INFO: backtonative mask unchanged — no cascade needed.',
                diary_file, 'OKGREEN')

    else:
        run_cmd.msg(
            'INFO: no backtonative mask saved — automatic ACPC mask will '
            'continue to be used. No cascade triggered.',
            diary_file, 'OKGREEN')

    return mask_was_updated
