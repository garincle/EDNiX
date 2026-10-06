import os
import datetime

opj = os.path.join
opi = os.path.isfile
ope = os.path.exists

from Tools import run_cmd


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


def clean(dir_prepro, volumes_dir, masks_dir, ID, diary_file):
    """
    Remove intermediate anatomical preprocessing files that are no longer
    needed after the full pipeline has completed.

    Files preserved (never removed):
      - {ID}_final_mask.nii.gz  — manual mask, irreplaceable
      - {ID}_space-acpc_desc-template_*  — used by step 5 completion check
      - {ID}_space-acpc_desc-Bias_*      — N4-corrected, used by step 5 check
      - {ID}_space-acpc_desc-SS_*        — final skull-strip, used in fMRI
      - {ID}_space-acpc_mask.nii.gz      — final ACPC mask, used in fMRI
      - {ID}_space-raw_desc-*            — orientation-corrected raw images
      - {ID}_desc-norm_*                 — used by FreeSurfer (step 10)

    Files removed (intermediate, large, reproducible):
      - {ID}_space-acpc_desc-SS-step1_*  — replaced by final SS
      - {ID}_space-acpc_desc-SS-step2_*  — replaced by final SS
      - {ID}_space-acpc_desc-64_*        — cropping intermediate
      - {ID}_space-acpc_desc-cropped_*   — cropping intermediate
      - {ID}_brain_for_Align_Center_*    — ACPC alignment intermediate
      - {ID}_desc-step1_mask             — replaced by final mask
      - {ID}_desc-step2_mask             — replaced by final mask
      - {ID}_space-acpc_desc-Zp_*        — zero-padding intermediate (small brain)
      - {ID}_space-acpc_desc-64_Zp_*     — zero-padding intermediate
    """

    ct = datetime.datetime.now()
    with open(diary_file, 'a') as diary:
        diary.write(f'\n{ct}\n')
        diary.write('Run anat._100_Data_Clean.clean\n')

    run_cmd.msg('Run anat._100_Data_Clean.clean', diary_file, 'HEADER')

    # ---- Files to preserve explicitly ------------------------------------
    final_mask = opj(masks_dir, ID + '_final_mask.nii.gz')

    # ---- Build removal list ----------------------------------------------
    list_to_remove = []

    for T in ['T1w', 'T2w']:
        list_to_remove += [
            opj(volumes_dir, ID + '_space-acpc_desc-SS-step1_' + T),
            opj(volumes_dir, ID + '_space-acpc_desc-SS-step2_' + T),
            opj(dir_prepro,  ID + '_space-acpc_desc-64_' + T),
            opj(dir_prepro,  ID + '_space-acpc_desc-cropped_' + T),
            opj(dir_prepro,  ID + '_brain_for_Align_Center_' + T),
            # Zero-padding intermediates (small brain pipeline)
            opj(dir_prepro,  ID + '_space-acpc_desc-Zp_' + T),
            opj(dir_prepro,  ID + '_space-acpc_desc-64_Zp_' + T),
        ]

    list_to_remove += [
        opj(masks_dir, ID + '_desc-step1_mask'),
        opj(masks_dir, ID + '_desc-step2_mask'),
    ]

    # ---- Show what will be deleted ---------------------------------------
    files_to_delete = []
    for stem in list_to_remove:
        for ext in ['.nii.gz', '.json']:
            path = stem + ext
            # Never remove final_mask even if it somehow ends up in the list
            if path in (final_mask, final_mask.replace('.nii.gz', '.json')):
                continue
            if opi(path):
                files_to_delete.append(path)

    if not files_to_delete:
        run_cmd.msg(f'INFO: Nothing to clean for {ID} — all intermediate files '
                    f'already removed or never created.', diary_file, 'OKGREEN')
        return

    print(bcolors.WARNING + f'\nWARNING: The following intermediate files for {ID} '
          f'will be permanently deleted:' + bcolors.ENDC)
    total_mb = 0
    for path in files_to_delete:
        try:
            size_mb = os.path.getsize(path) / (1024 * 1024)
            total_mb += size_mb
            print(f'  🗑️  {os.path.basename(path)}  ({size_mb:.1f} MB)')
        except Exception:
            print(f'  🗑️  {os.path.basename(path)}')

    print(bcolors.WARNING + f'\n  Total to free: ~{total_mb:.0f} MB' + bcolors.ENDC)

    # ---- Confirm preservation -------------------------------------------
    if opi(final_mask):
        print(bcolors.OKGREEN + f'\n  PRESERVING: {os.path.basename(final_mask)} '
              f'(manual mask — never deleted)' + bcolors.ENDC)
        run_cmd.msg(f'PRESERVING: {final_mask}', diary_file, 'OKGREEN')

    # ---- Confirmation dialog --------------------------------------------
    print(bcolors.FAIL + f'\n🚨 DELETE {len(files_to_delete)} files for {ID}?' + bcolors.ENDC)
    print(bcolors.FAIL + 'This action cannot be undone!' + bcolors.ENDC)
    response = input(bcolors.FAIL + "Type 'YES' to confirm, anything else to cancel: "
                     + bcolors.ENDC)

    if response.strip() != 'YES':
        msg = f'CANCELLED: Deletion cancelled for {ID}.'
        print(bcolors.OKBLUE + msg + bcolors.ENDC)
        run_cmd.msg(msg, diary_file, 'OKBLUE')
        return

    # ---- Delete ---------------------------------------------------------
    run_cmd.msg(f'PROCEEDING: Deleting {len(files_to_delete)} files for {ID}.',
                diary_file, 'OKGREEN')

    for path in files_to_delete:
        try:
            os.remove(path)
            run_cmd.msg(f'INFO: Removed {os.path.basename(path)}', diary_file, 'OKGREEN')
        except Exception as e:
            run_cmd.msg(f'ERROR: Cannot remove {path}: {e}', diary_file, 'FAIL')

    run_cmd.msg(f'INFO: Clean complete for {ID}. '
                f'Freed ~{total_mb:.0f} MB.', diary_file, 'OKGREEN')
