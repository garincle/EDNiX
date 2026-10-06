import os
import shutil
import datetime
import glob

opj = os.path.join
opb = os.path.basename
ope = os.path.exists
opi = os.path.isfile


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


def clean(dir_prepro_raw_process, dir_prepro_fmap,
          dir_prepro_acpc_process, dir_prepro_orig_process,
          dir_prepro_template_process, diary_file):
    """
    Modified cleaning: Selective removal in specific directories.

    - func/templates/          → KEPT entirely
    - func/anat/preprocessing/ → KEPT entirely
    - func/acpc-func/preprocessing/ → Selective removal
    - func/raw/preprocessing/  → Selective removal
    """

    ct = datetime.datetime.now()
    with open(diary_file, 'a') as f:
        f.write(f'\n{ct}\nRun fMRI._100_Data_Clean.clean (modified)\n')

    # ---- Summarise -------------------------------------------------------
    print(bcolors.FAIL
          + '\n╔══════════════════════════════════════════════════════════════╗'
          + bcolors.ENDC)
    print(bcolors.FAIL
          + '║  ⚠  MODIFIED CLEANING: Selective file removal                ║'
          + bcolors.ENDC)
    print(bcolors.FAIL
          + '║     - func/templates/ kept entirely                          ║'
          + bcolors.ENDC)
    print(bcolors.FAIL
          + '║     - func/anat/preprocessing/ kept entirely                 ║'
          + bcolors.ENDC)
    print(bcolors.FAIL
          + '║     - func/acpc-func/preprocessing/ selectively cleaned      ║'
          + bcolors.ENDC)
    print(bcolors.FAIL
          + '║     - func/raw/preprocessing/ selectively cleaned            ║'
          + bcolors.ENDC)
    print(bcolors.FAIL
          + '╚══════════════════════════════════════════════════════════════╝'
          + bcolors.ENDC)

    # ---- Show what will be removed --------------------------------------
    print(bcolors.WARNING + '\nDirectories with SELECTIVE removal:' + bcolors.ENDC)
    print(f'  📁 func/acpc-func/preprocessing/')
    print(f'     → Keeping: 7 session-level files + per-run fMRI_run_inRef/runMean_inRef files (needed to rerun step 7)')
    print(f'     → Removing: all other files')
    print(f'  📁 func/raw/preprocessing/')
    print(f'     → Keeping: *.1D, *.txt, and 2 specific files')
    print(f'     → Removing: all other files')
    print(bcolors.OKGREEN + '\nDirectories KEPT entirely:' + bcolors.ENDC)
    print(f'  ✅ func/templates/')
    print(f'  ✅ func/anat/preprocessing/')

    # ---- Confirmation ---------------------------------------------------
    print(bcolors.FAIL
          + '\n🚨 ARE YOU SURE? This will permanently delete files.'
          + bcolors.ENDC)
    response = input(bcolors.FAIL
                     + "Type 'YES' to confirm selective deletion: "
                     + bcolors.ENDC)

    if response.strip() != 'YES':
        msg = 'CANCELLED: User cancelled deletion.'
        print(bcolors.OKBLUE + msg + bcolors.ENDC)
        with open(diary_file, 'a') as f:
            f.write(f'\n{msg}\n')
        return

    msg = 'PROCEEDING: User confirmed deletion.'
    print(bcolors.OKGREEN + msg + bcolors.ENDC)
    with open(diary_file, 'a') as f:
        f.write(f'\n{msg}\n')

    # ---- 1. Clean func/acpc-func/preprocessing/ -------------------------
    acpc_func_preproc = dir_prepro_orig_process  # Should point to func/acpc-func/preprocessing/

    if ope(acpc_func_preproc) and os.path.isdir(acpc_func_preproc):
        # Files to KEEP — exact session-level filenames
        keep_files = [
            'all_runs_space-acpc-func_desc-fMRI_Mean_Image_SS.nii.gz',
            'all_runs_space-acpc-func_desc-fMRI_Mean_Image_SS.json',
            'all_runs_space-acpc-func_desc-fMRI_Mean_Image_SS_0GenericAffine.mat',
            'anat_space-acpc-func_T1w.nii.gz',
            'anat_space-acpc-func_T1w.json',
            'all_runs_space-acpc-func_desc-fMRI_Mean_Image_res-anat.json',
            'all_runs_space-acpc-func_desc-fMRI_Mean_Image_res-anat.nii.gz'
        ]

        # Per-run files to KEEP — these are the actual inputs step 7
        # (_7_post_TTT.py) reads to (re)compute residuals. Without these,
        # rerunning step 7 alone is impossible and forces a reprocess from
        # step 5 onward. Matched by suffix since the root_RS prefix varies
        # per run.
        keep_suffixes = [
            '_space-acpc-func_desc-fMRI_run_inRef.nii.gz',
            '_space-acpc-func_desc-fMRI_run_inRef.json',
            '_space-acpc-func_desc-fMRI_run_inRef_SS.nii.gz',
            '_space-acpc-func_desc-fMRI_run_inRef_SS.json',
            '_space-acpc-func_desc-fMRI_runMean_inRef.nii.gz',
            '_space-acpc-func_desc-fMRI_runMean_inRef.json',
        ]

        print(bcolors.WARNING + f'\nCleaning: {acpc_func_preproc}' + bcolors.ENDC)
        removed_count = 0
        kept_count = 0

        for root, dirs, files in os.walk(acpc_func_preproc):
            for file in files:
                file_path = opj(root, file)
                if file in keep_files or any(file.endswith(suf) for suf in keep_suffixes):
                    kept_count += 1
                    print(f'  ✅ KEPT: {file}')
                else:
                    try:
                        os.remove(file_path)
                        removed_count += 1
                        print(f'  🗑️  REMOVED: {file}')
                        with open(diary_file, 'a') as f:
                            f.write(f'\nRemoved: {file_path}\n')
                    except Exception as e:
                        msg = f'ERROR removing {file_path}: {e}'
                        print(bcolors.FAIL + msg + bcolors.ENDC)
                        with open(diary_file, 'a') as f:
                            f.write(f'\n{msg}\n')

        msg = f'SUCCESS: {acpc_func_preproc} cleaned. Kept {kept_count} files, removed {removed_count} files.'
        print(bcolors.OKGREEN + msg + bcolors.ENDC)
        with open(diary_file, 'a') as f:
            f.write(f'\n{msg}\n')
    else:
        msg = f'INFO: {acpc_func_preproc} not found — skipping'
        print(bcolors.WARNING + msg + bcolors.ENDC)
        with open(diary_file, 'a') as f:
            f.write(f'\n{msg}\n')

    # ---- 2. Clean func/raw/preprocessing/ -------------------------------
    raw_preproc = dir_prepro_raw_process  # Should point to func/raw/preprocessing/

    if ope(raw_preproc) and os.path.isdir(raw_preproc):
        # Patterns to KEEP
        keep_patterns = [
            '.1D',  # All .1D files
            '.txt',  # All .txt files
            'all_runs_space-func_desc-fMRI_Mean_Image_SS.nii.gz',
            'all_runs_space-func_desc-fMRI_Mean_Image.nii.gz'
        ]

        print(bcolors.WARNING + f'\nCleaning: {raw_preproc}' + bcolors.ENDC)
        removed_count = 0
        kept_count = 0

        for root, dirs, files in os.walk(raw_preproc):
            for file in files:
                file_path = opj(root, file)

                # Check if file should be kept
                keep = False
                for pattern in keep_patterns:
                    if file.endswith(pattern) or file == pattern:
                        keep = True
                        break

                if keep:
                    kept_count += 1
                    print(f'  ✅ KEPT: {file}')
                else:
                    try:
                        os.remove(file_path)
                        removed_count += 1
                        print(f'  🗑️  REMOVED: {file}')
                        with open(diary_file, 'a') as f:
                            f.write(f'\nRemoved: {file_path}\n')
                    except Exception as e:
                        msg = f'ERROR removing {file_path}: {e}'
                        print(bcolors.FAIL + msg + bcolors.ENDC)
                        with open(diary_file, 'a') as f:
                            f.write(f'\n{msg}\n')

        msg = f'SUCCESS: {raw_preproc} cleaned. Kept {kept_count} files, removed {removed_count} files.'
        print(bcolors.OKGREEN + msg + bcolors.ENDC)
        with open(diary_file, 'a') as f:
            f.write(f'\n{msg}\n')
    else:
        msg = f'INFO: {raw_preproc} not found — skipping'
        print(bcolors.WARNING + msg + bcolors.ENDC)
        with open(diary_file, 'a') as f:
            f.write(f'\n{msg}\n')

    # ---- 3. Note about directories kept entirely ------------------------
    msg = 'INFO: func/templates/ and func/anat/preprocessing/ kept entirely (no changes).'
    print(bcolors.OKGREEN + msg + bcolors.ENDC)
    with open(diary_file, 'a') as f:
        f.write(f'\n{msg}\n')

    msg = 'INFO: fMRI Clean complete (modified selective cleaning).'
    print(bcolors.OKGREEN + msg + bcolors.ENDC)
    with open(diary_file, 'a') as f:
        f.write(f'\n{msg}\n')