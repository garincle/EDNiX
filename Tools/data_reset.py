"""
Generic engine of the EDNiX "reset before relaunch" cleaners (anat, PET, CT).

The fMRI cleaner (modalities/fMRI/_000_Data_Reset.py) set the rules followed here:
  - raw BIDS files are never touched,
  - files drawn or corrected by hand (masks, crops, manual segmentations) are never touched,
  - everything else written by the pipeline under the modality folder is removed,
  - dry_run=True previews without deleting; a single 'YES' confirms a whole batch.

A modality cleaner only provides a plan: the folder, the loose files to keep and a `keep(path)` predicate for the
files inside its sub-folders. Files are removed one by one (never a blind rmtree), so a manual file stays exactly where
the pipeline looks for it; empty folders are pruned afterwards.
"""
import os
import fnmatch
import datetime

opj = os.path.join
ope = os.path.exists
opb = os.path.basename


class bcolors:
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'


def log(diary_file, msg):
    if diary_file:
        os.makedirs(os.path.dirname(diary_file) or '.', exist_ok=True)
        with open(diary_file, 'a') as f:
            f.write(f'\n{msg}\n')


def matches(name, patterns):
    """True if the file name matches one of the glob patterns (case-sensitive, basename only)."""
    return any(fnmatch.fnmatchcase(name, p) for p in patterns)


def is_raw_bids(name):
    """Raw BIDS files of a modality folder start with 'sub-'; EDNiX diaries and parameter files do not
    ('<ID> session <S><date>.txt', '<ID>launcher_parameters session ...')."""
    return name.startswith('sub-')


def collect(root, keep_loose, keep_file, skip_dirs=()):
    """
    Files to delete under `root`.

    keep_loose(name) -> bool : for files directly in root
    keep_file(path) -> bool  : for files inside the sub-folders of root
    skip_dirs                : names of sub-folders of root left entirely untouched

    Returns (to_delete, kept): two sorted lists of paths.
    """
    to_delete, kept = [], []
    if not os.path.isdir(root):
        return to_delete, kept
    for entry in sorted(os.listdir(root)):
        path = opj(root, entry)
        if os.path.isfile(path) or os.path.islink(path):
            (kept if keep_loose(entry) else to_delete).append(path)
        elif os.path.isdir(path):
            if entry in skip_dirs:
                kept.append(path + os.sep)
                continue
            for dirpath, _, files in os.walk(path):
                for f in files:
                    p = opj(dirpath, f)
                    (kept if keep_file(p) else to_delete).append(p)
    return sorted(to_delete), sorted(kept)


def prune_empty_dirs(root):
    """Remove the empty folders below root (root itself is kept)."""
    for dirpath, dirnames, filenames in os.walk(root, topdown=False):
        if dirpath == root:
            continue
        try:
            if not os.listdir(dirpath):
                os.rmdir(dirpath)
        except OSError:
            pass


def _size_mb(paths):
    tot = 0
    for p in paths:
        try:
            tot += os.path.getsize(p)
        except OSError:
            pass
    return tot / 1024 ** 2


def execute(label, root, to_delete, kept, diary_file, dry_run=True, skip_confirmation=False, verbose=True):
    """
    Preview, confirm and delete. Returns the number of files removed (0 for a dry run or a cancel).
    """
    ct = datetime.datetime.now()
    log(diary_file, f'{ct}\nRun {label}{" [DRY RUN]" if dry_run else ""}')

    if not to_delete:
        msg = f'INFO: {root}: nothing to remove.'
        print(bcolors.OKGREEN + msg + bcolors.ENDC)
        log(diary_file, msg)
        return 0

    if verbose:
        print(bcolors.WARNING + f'\n{len(to_delete)} file(s) to remove in {root} (~{_size_mb(to_delete):.0f} MB)'
              + bcolors.ENDC)
        for p in to_delete[:40]:
            print(f'  🗑️  {os.path.relpath(p, root)}')
        if len(to_delete) > 40:
            print(f'  ... and {len(to_delete) - 40} more')
        manual = [p for p in kept if not is_raw_bids(opb(p.rstrip(os.sep)))]
        print(bcolors.OKGREEN + f'\nKept: {len(kept)} file(s)/folder(s), of which manual / protected:' + bcolors.ENDC)
        for p in manual:
            print(f'  ✅ {os.path.relpath(p, root)}')
        if not manual:
            print('  (none)')

    if dry_run:
        print(bcolors.OKBLUE + 'DRY RUN — nothing was deleted.' + bcolors.ENDC)
        log(diary_file, f'DRY RUN: {len(to_delete)} file(s) would be removed.')
        return 0

    if not skip_confirmation:
        response = input(bcolors.FAIL + f"Type 'YES' to delete {len(to_delete)} file(s) in {root}: " + bcolors.ENDC)
        if response.strip() != 'YES':
            print(bcolors.OKBLUE + 'CANCELLED — nothing was deleted.' + bcolors.ENDC)
            log(diary_file, 'CANCELLED by the user.')
            return 0

    n = 0
    for p in to_delete:
        try:
            os.remove(p)
            n += 1
        except OSError as e:
            msg = f'ERROR removing {p}: {e}'
            print(bcolors.FAIL + msg + bcolors.ENDC)
            log(diary_file, msg)
    prune_empty_dirs(root)
    msg = f'SUCCESS: {n} file(s) removed in {root}; {len(kept)} kept.'
    print(bcolors.OKGREEN + msg + bcolors.ENDC)
    log(diary_file, msg)
    return n


def run_batch(plans, dry_run=True):
    """
    plans: list of dicts (label, root, to_delete, kept, diary_file), one per subject / session.
    dry_run=True: preview all. dry_run=False: preview all, ONE confirmation, then delete all.
    Returns the total number of files removed.
    EDNIX_SKIP_RESET=1 in the environment (set by the run_all_ednix launcher for its workers, which have no
    terminal to answer 'YES') turns every reset into a no-op.
    """
    if os.environ.get('EDNIX_SKIP_RESET'):
        print(bcolors.OKBLUE + 'INFO: EDNIX_SKIP_RESET is set: reset skipped.' + bcolors.ENDC)
        return 0
    for pl in plans:
        print('\n' + '-' * 70 + f'\n{pl["label"]}')
        execute(pl['label'], pl['root'], pl['to_delete'], pl['kept'], pl['diary_file'], dry_run=True)
    n_files = sum(len(pl['to_delete']) for pl in plans)
    if dry_run or n_files == 0:
        print(bcolors.OKBLUE + f'\n{len(plans)} subject/session(s), {n_files} file(s) would be removed. '
              f'Run again with dry_run=False to delete.' + bcolors.ENDC)
        return 0
    print('\n' + '=' * 70)
    response = input(bcolors.FAIL + f"🚨 Delete {n_files} file(s) for the {len(plans)} subject/session(s) above? "
                                    f"Type 'YES': " + bcolors.ENDC)
    if response.strip() != 'YES':
        print(bcolors.OKBLUE + 'CANCELLED — nothing was deleted.' + bcolors.ENDC)
        return 0
    total = 0
    for pl in plans:
        total += execute(pl['label'], pl['root'], pl['to_delete'], pl['kept'], pl['diary_file'],
                         dry_run=False, skip_confirmation=True, verbose=False)
    return total


def subjects(allinfo_study_c, bids_dir, list_to_keep=None, list_to_remove=None):
    """(ID, Session, data_path) of the selected subjects / sessions, as the launchers select them."""
    from Tools import Load_subject_with_BIDS
    all_ID, all_Session, all_data_path, *_ = Load_subject_with_BIDS.load_data_bids(
        allinfo_study_c, bids_dir, list_to_keep or [], list_to_remove or [])
    return list(zip(all_ID, all_Session, all_data_path))
