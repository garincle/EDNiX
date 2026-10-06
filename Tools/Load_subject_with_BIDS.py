import os
os.sep = '/'  # Override to force forward slashes globally
import numpy as np
import pandas as pd
from pathlib import Path

def linux_path(path_input):
    # Ensure it's a Path object first, then convert it to POSIX (Linux) style
    if not isinstance(path_input, Path):
        path_input = Path(path_input)  # Convert to Path if it's a string

    # Convert to Linux-style path
    linux_path = path_input.as_posix()

    # Check if it's a Windows path with a drive letter (e.g., 'C:/' or 'D:/')
    if len(linux_path) > 1 and linux_path[1] == ":":
        # Extract the drive letter and map it to /mnt/ in Linux (e.g., 'C:' -> '/mnt/c')
        drive_letter = linux_path[0].lower()  # Convert drive letter to lowercase
        linux_path = f"/mnt/{drive_letter}" + linux_path[2:].replace("\\", "/")  # Replace backslashes with forward slashes
    return linux_path


def print_included_tuples(allinfo_study_c):
    included_tuples = []

    for ID in pd.unique(allinfo_study_c.subject.dropna()):
        # Get sessions for this ID
        sessions = allinfo_study_c.loc[allinfo_study_c['subject'] == ID, 'session'].dropna().astype(str)

        # Sort the same way as load_data_bids
        try:
            sorted_sessions = sorted(sessions, key=lambda x: int(x), reverse=True)
        except (ValueError, TypeError):
            sorted_sessions = sorted(sessions, reverse=True)

        # Add to list
        for session in sorted_sessions:
            included_tuples.append((ID, session))

    print("Included tuples (ID, session):")
    for i, tup in enumerate(included_tuples):
        print(tup, end=', ' if i < len(included_tuples) - 1 else '')

def _session_sort(sessions, reverse=True):
    """Numeric order when every session is an integer, else alphabetic."""
    sessions = [str(x) for x in sessions]
    try:
        return sorted(sessions, key=lambda x: int(x), reverse=reverse)
    except (ValueError, TypeError):
        return sorted(sessions, reverse=reverse)


def load_data_bids(allinfo_study_c, bids_dir, list_to_keep=[], list_to_remove=[]):
    """
    Build the (subject, session, path) lists, then apply list_to_keep /
    list_to_remove.

    The "max session" lists (study template with which_on='max', longitudinal
    reference) now exclude list_to_remove entries: previously excluded
    sessions still entered the study template. They deliberately do NOT
    follow list_to_keep, which is also used to process subjects in batches
    (ednix_checkpoint): the template keeps using every non-excluded subject.
    """
    # all (ID, session) entries, sessions in descending order per subject
    entries = []
    seen = set()
    for ID in pd.unique(allinfo_study_c.subject.dropna()):
        sessions = allinfo_study_c.loc[allinfo_study_c['subject'] == ID].session.dropna()
        for session in _session_sort(sessions, reverse=True):
            if (ID, session) not in seen:
                seen.add((ID, session))
                entries.append((ID, session))

    # filters (IDs and sessions compared as strings; unmatched filter entries are reported)
    keep = {(str(i), str(s)) for i, s in list_to_keep}
    remove = {(str(i), str(s)) for i, s in list_to_remove}
    present = {(str(i), str(s)) for i, s in entries}
    for name, filt in (('list_to_keep', keep), ('list_to_remove', remove)):
        unmatched = sorted(filt - present)
        if unmatched:
            print(f'WARNING: {len(unmatched)} {name} entries match no (subject, session): {unmatched[:10]}')
    not_removed = [e for e in entries if (str(e[0]), str(e[1])) not in remove]
    if keep:
        entries = [e for e in entries if (str(e[0]), str(e[1])) in keep]
    elif remove:
        entries = not_removed

    # max session per subject among the non-excluded sessions (independent of list_to_keep)
    max_ses = {}
    for ID, session in not_removed:
        max_ses.setdefault(ID, []).append(session)
    max_ses = {ID: _session_sort(ss, reverse=True)[0] for ID, ss in max_ses.items()}

    all_ID = [e[0] for e in entries]
    all_Session = [e[1] for e in entries]
    all_data_path = [os.path.join(bids_dir, f'sub-{ID}', f'ses-{ses}') for ID, ses in entries]
    all_Session_max = [max_ses[ID] for ID in all_ID]
    all_ID_max = list(max_ses.keys())
    all_data_path_max = [os.path.join(bids_dir, f'sub-{ID}', f'ses-{max_ses[ID]}') for ID in all_ID_max]

    print("Remaining all_data_path:", all_data_path)
    return all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max
