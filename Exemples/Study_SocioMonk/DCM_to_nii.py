import os
import glob
import subprocess
import pandas as pd
import re
import os
# Precompiled binary linux_ubuntu_16_64: May  4 2020 (Version AFNI_20.1.06 'Otho')
# Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
spco = subprocess.check_output
opd = os.path.dirname

############# variables
out_analysis = '/home/cgarin/Documents/macaque/'
DICOMdir = '/home/common/benhalab/CASCAD/Macaque/MRI/Sociomonk/Test1_Litchy/sub-Litchi_ses-02/dicom'

in_DICOM = DICOMdir
out_code = out_analysis + '/code'
config_file = out_code + '/dcm2bids_config.json'
if not os.path.exists(out_analysis): os.mkdir((out_analysis))
if not os.path.exists(in_DICOM): os.mkdir((in_DICOM))

##creat  your scaffold
spco(['dcm2bids_scaffold', '-o', out_analysis, '--force'])

##creat configuration file
spco(['dcm2bids_helper', '-d', in_DICOM, '-o', out_analysis + '/sourcedata/', '--force'])

def extract_sub_and_session_id(text):
    # Define the regex pattern to match the subject and session IDs
    pattern = r'MD(\d+)_Mc([A-Z0-9]+)_P(\d+)_1_1'

    # Search for the pattern in the given text
    match = re.search(pattern, text)

    # If a match is found, format and return the extracted subject ID and session ID
    if match:
        subject_id = f'{match.group(2)}'
        session_id = f'{match.group(3)}'
        return subject_id, session_id
    else:
        return None, None

def list_folders(DICOMdir):
    """
    List all folders in the specified directory.

    Parameters:
    directory (str): The path to the directory to search.

    Returns:
    list: A list of folder names.
    """
    try:
        # Get a list of all entries in the directory
        entries = os.listdir(DICOMdir)

        # Filter out the entries that are directories
        folders = [entry for entry in entries if os.path.isdir(os.path.join(DICOMdir, entry))]
        paths =    [os.path.join(DICOMdir, folder) for folder in folders]
        return folders, paths
    except FileNotFoundError:
        print(f"The directory {DICOMdir} does not exist.")
        return []
    except PermissionError:
        print(f"Permission denied to access the directory {DICOMdir}.")
        return []

# List all folders in the specified directory
folders, paths = list_folders(DICOMdir)

# Print the list of folders
print("Folders in the directory:", out_analysis)
for folder in folders:
    print(folder)

# Méthode alternative pour supprimer les fichiers problématiques
import shutil

def remove_denoising_files(subject_dir):
    """Supprime tous les fichiers contenant 'Denoising' dans le dossier sujet"""
    for root, dirs, files in os.walk(subject_dir):
        for file in files:
            if 'Denoising' in file:
                file_path = os.path.join(root, file)
                print(f"Deleting: {file_path}")
                os.remove(file_path)


import os
import glob
import subprocess
import re


def renumber_runs(bids_subject_dir):
    """Renumérote les runs pour combler les trous"""
    func_dir = os.path.join(bids_subject_dir, "func")
    if not os.path.exists(func_dir):
        return

    # Trouver tous les fichiers bold restants
    bold_files = glob.glob(os.path.join(func_dir, "*_task-rest_*_bold.nii.gz"))
    bold_files.sort()  # Trier par ordre naturel

    # Extraire les numéros de run actuels
    runs_info = []
    for file in bold_files:
        match = re.search(r'_run-(\d+)_', file)
        if match:
            current_run = int(match.group(1))
            runs_info.append((current_run, file))

    # Trier par numéro de run original
    runs_info.sort(key=lambda x: x[0])

    # Renumérotter séquentiellement
    for new_run_num, (old_run_num, file_path) in enumerate(runs_info, 1):
        if new_run_num != old_run_num:
            # Nouveaux noms
            new_file = re.sub(r'_run-\d+_', f'_run-{new_run_num:02d}_', file_path)
            new_json = file_path.replace('.nii.gz', '.json')
            new_json_renamed = new_file.replace('.nii.gz', '.json')

            # Renommer .nii.gz
            print(f"📁 RENAMING: {os.path.basename(file_path)} -> {os.path.basename(new_file)}")
            os.rename(file_path, new_file)

            # Renommer .json
            if os.path.exists(new_json):
                os.rename(new_json, new_json_renamed)


for img, path in zip(folders, paths):
    sub_id = img[2:]
    ses_id = '1'
    print(f"Processing: {img} -> sub-{sub_id} ses-{ses_id}")

    try:
        # 1. CONVERTIR TOUT
        subprocess.check_output([
            'dcm2bids', '-d', path, '-p', sub_id, '-s', ses_id,
            '-c', config_file, '-o', out_analysis,
            '--auto_extract_entities', '--bids_validate', '--force_dcm2bids'
        ])

        # 2. SUPPRIMER LES FICHIERS DENOISING
        bids_subject_dir = os.path.join(out_analysis, f"sub-{sub_id}", f"ses-{ses_id}")

        if os.path.exists(bids_subject_dir):
            deleted_files = []

            # Supprimer les fichiers Denoising
            for root, dirs, files in os.walk(bids_subject_dir):
                for file in files:
                    file_path = os.path.join(root, file)

                    if 'Denoising' in file:
                        print(f"🗑️ DELETING: {file_path}")
                        os.remove(file_path)
                        deleted_files.append(file_path)

                    elif file.endswith('.json'):
                        try:
                            with open(file_path, 'r') as f:
                                content = f.read()
                                if 'Denoising' in content:
                                    print(f"🗑️ DELETING JSON: {file_path}")
                                    os.remove(file_path)
                                    deleted_files.append(file_path)

                                    # Supprimer le .nii.gz associé
                                    nii_file = file_path.replace('.json', '.nii.gz')
                                    if os.path.exists(nii_file):
                                        print(f"🗑️ DELETING NIfTI: {nii_file}")
                                        os.remove(nii_file)
                                        deleted_files.append(nii_file)
                        except:
                            pass

            # 3. RENUMÉROTER LES RUNS SI DES FICHIERS ONT ÉTÉ SUPPRIMÉS
            if deleted_files:
                print(f"🔄 Renumbering runs after deleting {len(deleted_files)} files")
                renumber_runs(bids_subject_dir)
            else:
                print("✅ No Denoising files found")

    except subprocess.CalledProcessError as e:
        print(f"❌ Error processing {img}: {e}")