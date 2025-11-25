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
    included_tuples = allinfo_study_c[['subject', 'session']].dropna().drop_duplicates()
    included_tuples.sort_values(
        by=['subject', 'session'],
        ascending=[True, False],
        inplace=True)
    included_tuples.reset_index(drop=True, inplace=True)
    print("Included tuples (ID, session):")
    for _, row in included_tuples.iterrows():
        print((row['subject'], row['session']), end=', ')

def load_data_bids(allinfo_study_c, bids_dir, list_to_keep, list_to_remove):
    all_data_path = []
    all_ID = []
    all_Session = []
    all_data_path_max = []
    all_ID_max = []
    max_session = []
    all_Session_max = []
    animal_ID = []

    for ID in pd.unique(allinfo_study_c.subject.dropna()):
        list_session = allinfo_study_c.loc[allinfo_study_c['subject'] == ID].session.dropna()

        # Step 1: Extract and reverse the sessions while keeping leading zeros
        try:
            # Essayer de traiter comme des nombres
            listereverse = sorted(list_session.astype(str), key=lambda x: int(x), reverse=True)
            # Calculer la session max (comme string, en gardant les zéros initiaux)
            max_ses = max(listereverse, key=lambda x: int(x))
        except (ValueError, TypeError):
            # Si échec (valeurs non numériques), trier par ordre alphabétique
            listereverse = sorted(list_session.astype(str), reverse=True)
            max_ses = max(listereverse)
        max_session.append(max_ses)

        for session in listereverse:
            data_path = os.path.join(bids_dir, 'sub-' + ID, 'ses-' + str(session))
            all_data_path.append(data_path)
            all_Session.append(session)
            all_ID.append(ID)
            all_Session_max.append(max_ses)
            animal_ID.append(ID + 'ses-' + str(session))

    for ID, Session in zip(pd.unique(allinfo_study_c.subject.dropna()), max_session):
        data_path = os.path.join(bids_dir, 'sub-' + ID, 'ses-' + str(Session))
        all_data_path_max.append(data_path)
        all_ID_max.append(ID)

    if len(list_to_keep) > 0 and len(list_to_remove) == 0:
        removelist = []
        for num, (ID, Session, data_path, max_ses) in enumerate(zip(all_ID, all_Session, all_data_path, all_Session_max)):
            if (ID, Session) in list_to_keep:
                removelist.append(num)

        all_ID = [item for i, item in enumerate(all_ID) if i in removelist]
        all_Session = [item for i, item in enumerate(all_Session) if i in removelist]
        all_data_path = [item for i, item in enumerate(all_data_path) if i in removelist]
        all_Session_max = [item for i, item in enumerate(all_Session_max) if i in removelist]

    elif len(list_to_keep) == 0 and len(list_to_remove) > 0:
        removelist = []
        for num, (ID, Session, data_path, max_ses) in enumerate(zip(all_ID, all_Session, all_data_path, all_Session_max)):
            if (ID, Session) in list_to_remove:
                removelist.append(num)

        all_ID = [item for i, item in enumerate(all_ID) if i not in removelist]
        all_Session = [item for i, item in enumerate(all_Session) if i not in removelist]
        all_data_path = [item for i, item in enumerate(all_data_path) if i not in removelist]
        all_Session_max = [item for i, item in enumerate(all_Session_max) if i not in removelist]

    elif len(list_to_keep) == 0 and len(list_to_remove) == 0:
        print('All subjects are included for the analysis')

    else:
        raise ValueError('One of list_to_keep or list_to_remove needs to be empty')

    print("Remaining all_data_path:", all_data_path)
    return all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max