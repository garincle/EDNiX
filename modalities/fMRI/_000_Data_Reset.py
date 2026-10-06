import os
import shutil
import datetime
from Tools import Load_subject_with_BIDS

opj = os.path.join
ope = os.path.exists


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


def _log(diary_file, msg):
    with open(diary_file, 'a') as f:
        f.write(f'\n{msg}\n')


def _remove_dir(path, diary_file, dry_run):
    """Remove a directory tree, logging what happened."""
    if not ope(path):
        return
    if dry_run:
        print(f'  🔎 WOULD REMOVE: {path}')
        _log(diary_file, f'[DRY RUN] Would remove: {path}')
    else:
        try:
            shutil.rmtree(path)
            print(f'  🗑️  REMOVED: {path}')
            _log(diary_file, f'Removed: {path}')
        except Exception as e:
            msg = f'ERROR removing {path}: {e}'
            print(bcolors.FAIL + msg + bcolors.ENDC)
            _log(diary_file, msg)


def reset_func(data_path, diary_file, dry_run=False, skip_confirmation=False):
    """
    Full reset of a subject/session's func/ (and sibling fmap/) folders,
    to relaunch EDNiX.func with no leftover files from a previous attempt.

    KEPT:
      - Any loose file directly inside func/ whose name contains '_task-'
        (the raw BIDS scans, their .json sidecars, and any per-run .txt
        volume-trim files) — never touched.
      - func/raw/masks/ entirely — this is where the manually-verified
        final_mask lives (_4_check_mask.py), so skull-stripping QC does
        not need to be redone.
      - Any loose file directly inside fmap/ (raw fieldmap acquisitions,
        if present) — never touched.

    REMOVED:
      - Any loose file directly inside func/ whose name does NOT contain
        '_task-' — e.g. the diary and launcher_parameters .txt files
        written by diaryfile.create(), which are pipeline output, not
        raw data.
      - func/raw/preprocessing/ and func/raw/matrices/ (step 1 outputs —
        motion correction will be recomputed from scratch).
      - func/acpc-func/, func/acpc-anat/, func/templates/ entirely.
      - func/QC/ entirely, and any other directory found directly under
        func/ that isn't 'raw' (this makes the reset robust to output
        folders not explicitly listed here).
      - fmap/preprocessing/.

    This mirrors the confirm-before-delete pattern used by
    _100_Data_Clean.py, but is a full reset rather than a selective
    space-saving cleanup — use this BEFORE relaunching the pipeline on
    a subject, not after finishing it.

    Parameters
    ----------
    data_path : str
        The subject/session directory that directly contains 'func' and
        'fmap' subfolders (i.e. what's normally passed to getpath.func).
    diary_file : str
        Path to the run's diary/log file.
    dry_run : bool
        If True, print/log what WOULD be removed without deleting
        anything. Recommended as a first pass before actually running.
    skip_confirmation : bool
        If True, skip the per-subject 'YES' prompt and proceed directly.
        Used by reset_all() to ask for confirmation once for the whole
        batch instead of once per subject/session.
    """

    path_func = opj(data_path, 'func')
    path_fmap = opj(data_path, 'fmap')

    ct = datetime.datetime.now()
    _log(diary_file, f'{ct}\nRun fMRI._000_Data_Reset.reset_func'
                      f'{" [DRY RUN]" if dry_run else ""}\n')

    print(bcolors.FAIL
          + '\n╔══════════════════════════════════════════════════════════════╗'
          + bcolors.ENDC)
    print(bcolors.FAIL
          + '║  ⚠  FULL RESET: remove all preprocessing output              ║'
          + bcolors.ENDC)
    print(bcolors.FAIL
          + '║     - Raw BIDS scans/sidecars/trim-files kept (name has "_task-") ║'
          + bcolors.ENDC)
    print(bcolors.FAIL
          + '║     - Other loose files in func/ (diary, launcher_parameters) removed ║'
          + bcolors.ENDC)
    print(bcolors.FAIL
          + '║     - func/raw/masks/ (manual mask) kept                     ║'
          + bcolors.ENDC)
    print(bcolors.FAIL
          + '║     - Raw fieldmap files (loose files in fmap/) kept          ║'
          + bcolors.ENDC)
    print(bcolors.FAIL
          + '║     - Everything else under func/ and fmap/ removed          ║'
          + bcolors.ENDC)
    print(bcolors.FAIL
          + '╚══════════════════════════════════════════════════════════════╝'
          + bcolors.ENDC)

    if not ope(path_func):
        msg = f'INFO: {path_func} not found — nothing to reset.'
        print(bcolors.WARNING + msg + bcolors.ENDC)
        _log(diary_file, msg)
        return

    # ---- Figure out exactly what will be removed, before touching anything
    to_remove = []      # directories
    files_to_remove = []  # loose files directly in func/ that aren't raw BIDS data

    raw_dir = opj(path_func, 'raw')
    if ope(raw_dir):
        for sub in os.listdir(raw_dir):
            if sub == 'masks':
                continue  # kept — manual mask lives here
            sub_path = opj(raw_dir, sub)
            if os.path.isdir(sub_path):
                to_remove.append(sub_path)

    for entry in os.listdir(path_func) if ope(path_func) else []:
        entry_path = opj(path_func, entry)
        if entry == 'raw':
            continue  # handled above (selectively)
        if os.path.isdir(entry_path):
            to_remove.append(entry_path)
        elif os.path.isfile(entry_path):
            # Raw BIDS scans/sidecars/trim-files all contain '_task-' in
            # their name (matches this pipeline's own endfmri/endjson
            # patterns). Anything else loose in func/ — e.g. the diary and
            # launcher_parameters .txt files written by diaryfile.create()
            # — is pipeline output, not raw data, and gets removed.
            if '_task-' not in entry:
                files_to_remove.append(entry_path)

    fmap_prepro = opj(path_fmap, 'preprocessing')
    if ope(fmap_prepro):
        to_remove.append(fmap_prepro)

    if not to_remove and not files_to_remove:
        msg = 'INFO: nothing found to remove — folder already clean.'
        print(bcolors.OKGREEN + msg + bcolors.ENDC)
        _log(diary_file, msg)
        return

    print(bcolors.WARNING + f'\n{len(to_remove)} folder(s) will be removed:' + bcolors.ENDC)
    for p in to_remove:
        print(f'  📁 {p}')

    print(bcolors.WARNING + f'\n{len(files_to_remove)} loose file(s) will be removed:' + bcolors.ENDC)
    for p in files_to_remove:
        print(f'  📄 {p}')

    print(bcolors.OKGREEN + '\nKept (not touched):' + bcolors.ENDC)
    print(f'  ✅ Raw BIDS scans/sidecars/trim-files in {path_func}/ (files containing "_task-")')
    print(f'  ✅ {raw_dir}/masks/  (manual mask)')
    print(f'  ✅ Loose files directly in {path_fmap}/ (raw fieldmap data)')

    if dry_run:
        print(bcolors.OKBLUE + '\nDRY RUN — nothing was deleted. Re-run with dry_run=False to actually remove these.' + bcolors.ENDC)
        _log(diary_file, 'DRY RUN complete — nothing deleted.')
        return

    # ---- Confirmation ---------------------------------------------------
    if skip_confirmation:
        msg = 'PROCEEDING: confirmation given for the whole batch (reset_all).'
        print(bcolors.OKGREEN + msg + bcolors.ENDC)
        _log(diary_file, msg)
    else:
        print(bcolors.FAIL
              + '\n🚨 ARE YOU SURE? This will permanently delete the folders listed above.'
              + bcolors.ENDC)
        response = input(bcolors.FAIL
                         + "Type 'YES' to confirm full reset: "
                         + bcolors.ENDC)

        if response.strip() != 'YES':
            msg = 'CANCELLED: User cancelled reset.'
            print(bcolors.OKBLUE + msg + bcolors.ENDC)
            _log(diary_file, msg)
            return

        msg = 'PROCEEDING: User confirmed full reset.'
        print(bcolors.OKGREEN + msg + bcolors.ENDC)
        _log(diary_file, msg)

    for p in to_remove:
        _remove_dir(p, diary_file, dry_run=False)

    for p in files_to_remove:
        try:
            os.remove(p)
            print(f'  🗑️  REMOVED: {p}')
            _log(diary_file, f'Removed: {p}')
        except Exception as e:
            msg = f'ERROR removing {p}: {e}'
            print(bcolors.FAIL + msg + bcolors.ENDC)
            _log(diary_file, msg)

    msg = 'SUCCESS: full reset complete. Ready to relaunch EDNiX.func.'
    print(bcolors.OKGREEN + msg + bcolors.ENDC)
    _log(diary_file, msg)


def reset_all(allinfo_study_c, bids_dir, list_to_keep=None, list_to_remove=None, dry_run=True):
    """
    Run reset_func on every subject/session found in allinfo_study_c.

    Pass the SAME allinfo_study_c (already built by your own script via
    load_bids.Load_BIDS_to_pandas with whatever modality/suffix fits
    your study — e.g. modalities=['anat'], suffixes=['T2w'] for some
    studies, modalities=['func'], suffixes=['bold'] for others) and the
    SAME list_to_keep/list_to_remove you use for your actual relaunch
    call, so the reset and the relaunch are always scoped to the same
    subjects. This avoids re-deriving subjects with a hardcoded
    modality/suffix internally, which can silently return the wrong
    (or an empty) subject list for studies that discover subjects
    differently.

    Parameters
    ----------
    allinfo_study_c : pandas.DataFrame
        Output of load_bids.Load_BIDS_to_pandas(...) from your own script.
    bids_dir : str
        Root of the BIDS dataset (same value passed to Load_BIDS_to_pandas).
    list_to_keep : list of (ID, Session) tuples, optional
        If non-empty, ONLY these subject/sessions are reset.
    list_to_remove : list of (ID, Session) tuples, optional
        These subject/sessions are excluded from the reset.
    dry_run : bool
        If True (default), only preview what would be removed — nothing
        is deleted and no confirmation is asked. Always run this first.
        If False, previews every subject/session, then asks for ONE
        'YES' confirmation covering the whole batch (not one prompt per
        subject).
    """
    list_to_keep = list_to_keep or []
    list_to_remove = list_to_remove or []

    all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max = \
        Load_subject_with_BIDS.load_data_bids(allinfo_study_c, bids_dir, list_to_keep, list_to_remove)

    if dry_run:
        print('\n' + '=' * 70)
        print('DRY RUN — previewing reset for all selected subjects/sessions')
        print('=' * 70)

        for ID, Session, data_path in zip(all_ID, all_Session, all_data_path):
            print(f'\n--- {ID} session {Session} ---')
            diary_file = opj(data_path, str(ID) + '_session_' + str(Session) + '_reset_diary.txt')
            reset_func(data_path, diary_file, dry_run=True)
        return

    # dry_run=False: preview every subject/session first, then ask for ONE
    # confirmation covering the whole batch instead of prompting per-subject.
    print('\n' + '=' * 70)
    print(f'PREVIEW — {len(all_ID)} subject/session(s) selected for reset')
    print('=' * 70)

    for ID, Session, data_path in zip(all_ID, all_Session, all_data_path):
        print(f'\n--- {ID} session {Session} ---')
        diary_file = opj(data_path, str(ID) + '_session_' + str(Session) + '_reset_diary.txt')
        reset_func(data_path, diary_file, dry_run=True)

    print('\n' + '=' * 70)
    print(bcolors.FAIL
          + f"🚨 ARE YOU SURE? This will permanently delete files for all "
            f"{len(all_ID)} subject/session(s) listed above."
          + bcolors.ENDC)
    response = input(bcolors.FAIL
                     + "Type 'YES' to confirm reset for ALL of the above: "
                     + bcolors.ENDC)

    if response.strip() != 'YES':
        print(bcolors.OKBLUE + 'CANCELLED — nothing was deleted.' + bcolors.ENDC)
        return

    print(bcolors.OKGREEN + f'\nProceeding with reset for all {len(all_ID)} subject/session(s)...' + bcolors.ENDC)
    for ID, Session, data_path in zip(all_ID, all_Session, all_data_path):
        print(f'\n--- {ID} session {Session} ---')
        diary_file = opj(data_path, str(ID) + '_session_' + str(Session) + '_reset_diary.txt')
        reset_func(data_path, diary_file, dry_run=False, skip_confirmation=True)