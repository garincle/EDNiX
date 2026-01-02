import os
import glob
import shutil
from tqdm import tqdm

opj = os.path.join


def cleanBIDS(BIDS_folder):
    # Collecting the files you want to keep
    keep_files = set()

    # Patterns plus spécifiques pour ne garder que les fichiers BIDS de base
    patterns = [
        # Anatomical files in BIDS root structure
        '/sub-*/ses-*/anat/*.nii.gz',
        '/sub-*/ses-*/anat/*.json',
        '/sub-*/anat/*.nii.gz',  # For no session structure
        '/sub-*/anat/*.json',  # For no session structure

        # Functional files in BIDS root structure
        '/sub-*/ses-*/func/*.nii.gz',
        '/sub-*/ses-*/func/*.json',
        '/sub-*/func/*.nii.gz',  # For no session structure
        '/sub-*/func/*.json',  # For no session structure

        # Specific mask files you want to keep
        '/sub-*/ses-*/anat/native/volumes/masks/*final_mask_orig.nii.gz',
        '/sub-*/ses-*/anat/native/volumes/masks/*final_mask.nii.gz',
        '/sub-*/ses-*/anat/native/volumes/masks/*final_mask2.nii.gz',
        '/sub-*/ses-*/func/preprocessing/space-func/masks/*_final_mask.nii.gz',
        '/sub-*/ses-*/func/preprocessing/space-func/masks/*_final_mask_orig.nii.gz',

        # Study template files
        '/sty_template/derivatives/acpc/volumes/masks/studyTemplate_final_mask.nii.gz',

        # Top-level BIDS files
        '/code/*',
        '/sourcedata/*',
        '/derivatives/*',

        # Surface exemples
        '/sub-*/ses-*/anat/native/surfaces/Native_resol/Exemple1.scene',

        # BIDS required files
        '/dataset_description.json',
        '/participants.tsv',
        '/participants.json',
        '/CHANGES',
        '/README'
    ]

    for pattern in patterns:
        full_pattern = BIDS_folder + pattern
        keep_files.update(glob.glob(full_pattern, recursive=True))

    # Also keep top-level files in BIDS_folder
    keep_files.update(f for f in glob.glob(BIDS_folder + '/*') if os.path.isfile(f))

    print(f"\nFound {len(keep_files)} files to KEEP")
    print("Sample of files to keep:")
    for f in sorted(keep_files)[:20]:  # Show only first 20 to avoid flooding
        print(f"  {f}")

    # Ask for user confirmation before proceeding
    confirm = 'yes'
    if confirm not in ['yes', 'y']:
        print("Operation aborted. No files were deleted.")
        return

    # Collect all files, not directories
    all_files = set()
    for root, dirs, files in os.walk(BIDS_folder):
        for file in files:
            full_path = opj(root, file)
            all_files.add(full_path)

    # Determine files to delete
    to_delete = all_files - keep_files

    print(f"\nDeleting {len(to_delete)} files...")

    # Deleting files with progress bar
    for path in tqdm(to_delete, desc="Deleting files"):
        try:
            os.remove(path)
        except Exception as e:
            print(f"Error deleting {path}: {e}")

    # Collect all directories, recursively using os.walk for better control
    all_dirs = set()
    for root, dirs, files in os.walk(BIDS_folder):
        for dir_name in dirs:
            full_path = opj(root, dir_name)
            all_dirs.add(full_path)

    print(f"\nCleaning empty directories...")

    # Deleting only directories that are truly empty
    deleted_dirs = 0
    for dir_path in sorted(all_dirs, reverse=True):  # reverse so children are processed before parents
        try:
            if not os.listdir(dir_path):
                shutil.rmtree(dir_path)
                deleted_dirs += 1
                print(f"Deleted empty directory: {dir_path}")
        except Exception as e:
            print(f"Error checking/deleting directory {dir_path}: {e}")

    print(f"\nCleanup complete!")
    print(f"Deleted {len(to_delete)} files and {deleted_dirs} empty directories")
    print(f"Kept {len(keep_files)} files")

list_paths = ['/scratch2/EDNiX/Rat/BIDS_Gd/',
              '/scratch2/EDNiX/Dog/BIDS_k9/',
              '/scratch2/EDNiX/Marmoset/BIDS_NIH_MBM/',
              '/scratch2/EDNiX/Human/BIDS_ds004513-raw-data/',
              '/scratch2/EDNiX/Human/ds004856/',
              '/scratch2/EDNiX/Macaque/BIDS_Cdt_Garin/',
              '/scratch2/EDNiX/Macaque/BIDS_BenHamed/',
              '/scratch2/EDNiX/Mouse_lemur/BIDS_Garin/',
              '/scratch2/EDNiX/Bat/BIDS_bat/',
              '/scratch2/EDNiX/Mouse/BIDS_Gd/']

#list_paths = glob.glob('/scratch2/EDNiX/*/*/')

for path in list_paths:
    print(f"\n{'=' * 80}")
    print(f"Processing: {path}")
    print(f"{'=' * 80}")
    cleanBIDS(path)
