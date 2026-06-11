import os
import shutil
import datetime

opj = os.path.join
opb = os.path.basename
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


def clean(dir_prepro_raw_process, dir_prepro_fmap,
          dir_prepro_acpc_process, dir_prepro_orig_process,
          dir_prepro_template_process, diary_file):
    """
    Remove intermediate fMRI preprocessing directories.

    ╔══════════════════════════════════════════════════════════════════════╗
    ║  WARNING — ONLY RUN AFTER THE FULL PIPELINE IS COMPLETE             ║
    ║                                                                      ║
    ║  The following directories are read by downstream steps:            ║
    ║                                                                      ║
    ║  dir_prepro_raw_process  (raw/preprocessing/)                       ║
    ║    → Step 7 reads: censor.1D, demean.1D, deriv.1D                  ║
    ║    → Step 12 reads: motion_enorm.1D, deriv.1D                      ║
    ║                                                                      ║
    ║  dir_prepro_orig_process (acpc-func/preprocessing/)                 ║
    ║    → Step 7 reads: run_inRef.nii.gz, run_inRef_SS.nii.gz           ║
    ║    → Step 8 reads: Mean_Image_SS.nii.gz                             ║
    ║    → Step 9 reads: Mean_Image_SS.nii.gz                             ║
    ║    → Step 11 reads: anat_space-acpc-func_{T}.nii.gz                ║
    ║    → Step 12 reads: run_inRef.nii.gz, tSNR files                   ║
    ║                                                                      ║
    ║  dir_prepro_acpc_process (acpc-anat/preprocessing/)                 ║
    ║    → Step 8 reads: Mean_Image_unwarped, Mean_Image_apply            ║
    ║    → Step 9 reads: Mean_Image_unwarped                              ║
    ║    → Step 11 reads: anat_space-acpc_res-func_{T}.nii.gz            ║
    ║                                                                      ║
    ║  dir_prepro_template_process (templates/EDNiX/preprocessing/)       ║
    ║    → Step 9  reads: BASE_SS_fMRI.nii.gz                             ║
    ║    → Step 11 reads: BASE_SS_fMRI.nii.gz                             ║
    ║    → Step 12 reads: BASE_SS_fMRI.nii.gz                             ║
    ║                                                                      ║
    ║  Deleting these means steps 7-12 CANNOT be rerun.                   ║
    ║  Only proceed if all analysis is complete and verified.              ║
    ╚══════════════════════════════════════════════════════════════════════╝

    Directories PRESERVED (never touched):
      - postprocessed_rs/      — residuals (FC, SBA inputs)
      - Stats/                 — correlation matrices, SBA results
      - raw/masks/             — brain masks
      - raw/matrices/          — motion parameters and transforms
      - func/QC/               — QC reports
    """

    ct = datetime.datetime.now()
    with open(diary_file, 'a') as f:
        f.write(f'\n{ct}\nRun fMRI._100_Data_Clean.clean\n')

    # ---- Directories to remove ------------------------------------------
    dirs_to_remove = [
        dir_prepro_raw_process,
        dir_prepro_fmap,
        dir_prepro_acpc_process,
        dir_prepro_orig_process,
        dir_prepro_template_process,
    ]

    # ---- Summarise -------------------------------------------------------
    print(bcolors.FAIL
          + '\n╔══════════════════════════════════════════════════════════════╗'
          + bcolors.ENDC)
    print(bcolors.FAIL
          + '║  ⚠  IMPORTANT: These directories are needed to rerun        ║'
          + bcolors.ENDC)
    print(bcolors.FAIL
          + '║     steps 7–12 (signal regression, SBA, QC).                ║'
          + bcolors.ENDC)
    print(bcolors.FAIL
          + '║     Only proceed if the full pipeline is COMPLETE.           ║'
          + bcolors.ENDC)
    print(bcolors.FAIL
          + '╚══════════════════════════════════════════════════════════════╝'
          + bcolors.ENDC)

    print(bcolors.WARNING
          + '\nDirectories to be PERMANENTLY DELETED:' + bcolors.ENDC)

    total_mb = 0
    found = []
    for d in dirs_to_remove:
        if ope(d) and os.path.isdir(d):
            try:
                size_mb = sum(
                    os.path.getsize(opj(root, fn))
                    for root, _, files in os.walk(d)
                    for fn in files) / (1024 * 1024)
                n = sum(len(files) for _, _, files in os.walk(d))
                total_mb += size_mb
                found.append((d, n, size_mb))
                print(f'  🗑️  {d}')
                print(f'       {n} files  ({size_mb:.0f} MB)')
            except Exception:
                found.append((d, '?', '?'))
                print(f'  🗑️  {d}')
        else:
            print(f'  —   {d}  (not found — already deleted)')

    if not found:
        msg = 'INFO: Nothing to clean — all directories already removed.'
        print(bcolors.OKGREEN + msg + bcolors.ENDC)
        with open(diary_file, 'a') as f:
            f.write(f'\n{msg}\n')
        return

    print(bcolors.WARNING
          + f'\n  Total to free: ~{total_mb:.0f} MB' + bcolors.ENDC)

    print(bcolors.OKGREEN + '\nDirectories PRESERVED:' + bcolors.ENDC)
    func_root = os.path.dirname(os.path.dirname(dir_prepro_raw_process))
    for label, path in [
        ('postprocessed_rs (residuals)',  'acpc-func/postprocessed_rs + acpc-anat/postprocessed_rs'),
        ('Stats (FC matrices, SBA)',      'acpc-func/Stats + templates/EDNiX/Stats'),
        ('raw/masks (brain masks)',       opj(func_root, 'raw', 'masks')),
        ('raw/matrices (motion params)', opj(func_root, 'raw', 'matrices')),
        ('func/QC (QC reports)',          opj(func_root, 'QC')),
    ]:
        print(f'  ✓  {label}')

    # ---- Confirmation ---------------------------------------------------
    print(bcolors.FAIL
          + '\n🚨 ARE YOU SURE? Deleting these means steps 7-12 cannot be rerun.'
          + bcolors.ENDC)
    response = input(bcolors.FAIL
                     + "Type 'YES' to confirm deletion, anything else to cancel: "
                     + bcolors.ENDC)

    if response.strip() != 'YES':
        msg = 'CANCELLED: User cancelled deletion.'
        print(bcolors.OKBLUE + msg + bcolors.ENDC)
        with open(diary_file, 'a') as f:
            f.write(f'\n{msg}\n')
        return

    # ---- Delete ---------------------------------------------------------
    msg = 'PROCEEDING: User confirmed deletion.'
    print(bcolors.OKGREEN + msg + bcolors.ENDC)
    with open(diary_file, 'a') as f:
        f.write(f'\n{msg}\n')

    freed_mb = 0
    for d in dirs_to_remove:
        if ope(d) and os.path.isdir(d):
            try:
                size_mb = sum(
                    os.path.getsize(opj(root, fn))
                    for root, _, files in os.walk(d)
                    for fn in files) / (1024 * 1024)
                shutil.rmtree(d)
                freed_mb += size_mb
                msg = f'SUCCESS: {d} deleted ({size_mb:.0f} MB)'
                print(bcolors.OKGREEN + msg + bcolors.ENDC)
                with open(diary_file, 'a') as f:
                    f.write(f'\n{msg}\n')
            except Exception as e:
                msg = f'ERROR: Cannot delete {d}: {e}'
                print(bcolors.FAIL + msg + bcolors.ENDC)
                with open(diary_file, 'a') as f:
                    f.write(f'\n{msg}\n')
        else:
            with open(diary_file, 'a') as f:
                f.write(f'\nINFO: {d} not found — skipped\n')

    msg = f'INFO: fMRI Clean complete. Freed ~{freed_mb:.0f} MB.'
    print(bcolors.OKGREEN + msg + bcolors.ENDC)
    with open(diary_file, 'a') as f:
        f.write(f'\n{msg}\n')
