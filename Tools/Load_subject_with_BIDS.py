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

def load_data_bids(allinfo_study_c, bids_dir, list_to_keep=[], list_to_remove=[]):
    all_data_path = []
    all_ID = []
    all_Session = []
    all_data_path_max = []
    all_ID_max = []
    all_Session_max = []
    animal_ID = []

    seen_entries = set()  # track unique combinations

    # First, compute max session per subject
    max_session_dict = {}
    for ID in pd.unique(allinfo_study_c.subject.dropna()):
        list_session = allinfo_study_c.loc[allinfo_study_c['subject'] == ID].session.dropna()
        try:
            # Numeric sort
            listereverse = sorted(list_session.astype(str), key=lambda x: int(x), reverse=True)
            max_ses = max(listereverse, key=lambda x: int(x))
        except (ValueError, TypeError):
            # Alphabetic sort fallback
            listereverse = sorted(list_session.astype(str), reverse=True)
            max_ses = max(listereverse)
        max_session_dict[ID] = (listereverse, max_ses)

    # Collect all entries
    for ID, (sessions, max_ses) in max_session_dict.items():
        for session in sessions:
            data_path = os.path.join(bids_dir, f'sub-{ID}', f'ses-{session}')
            entry = (ID, session, data_path, max_ses)
            if entry not in seen_entries:
                seen_entries.add(entry)
                all_ID.append(ID)
                all_Session.append(session)
                all_data_path.append(data_path)
                all_Session_max.append(max_ses)
                animal_ID.append(f'{ID}ses-{session}')

        # Max session entries
        data_path_max = os.path.join(bids_dir, f'sub-{ID}', f'ses-{max_ses}')
        entry_max = (ID, max_ses, data_path_max)
        if entry_max not in seen_entries:
            seen_entries.add(entry_max)
            all_ID_max.append(ID)
            all_data_path_max.append(data_path_max)

    # Apply list_to_keep / list_to_remove filters
    def filter_lists(IDs, Sessions, Paths, Sessions_max):
        if len(list_to_keep) > 0:
            mask = [(ID, Session) in list_to_keep for ID, Session in zip(IDs, Sessions)]
        elif len(list_to_remove) > 0:
            mask = [(ID, Session) not in list_to_remove for ID, Session in zip(IDs, Sessions)]
        else:
            mask = [True] * len(IDs)
        IDs = [id_ for id_, m in zip(IDs, mask) if m]
        Sessions = [s for s, m in zip(Sessions, mask) if m]
        Paths = [p for p, m in zip(Paths, mask) if m]
        Sessions_max = [sm for sm, m in zip(Sessions_max, mask) if m]
        return IDs, Sessions, Paths, Sessions_max

    all_ID, all_Session, all_data_path, all_Session_max = filter_lists(
        all_ID, all_Session, all_data_path, all_Session_max
    )

    print("Remaining all_data_path:", all_data_path)
    return all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max