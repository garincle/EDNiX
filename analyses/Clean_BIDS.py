import os
import glob
import shutil
from tqdm import tqdm

def cleanBIDS(BIDS_folder):
    # Collecting the files you want to keep
    keep_files = set()
    patterns = [
        '/**/**/**/sub-*',
        '/**/**/anat/native/02_Wb/volumes/masks/*final_mask.nii.gz',
        '/**/**/anat/native/02_Wb/volumes/masks/*final_mask_2.nii.gz',
        '/**/**/**/**/**/manual_mask.nii.gz',
        '/sty_template/studytemplate2_T2w/study_template_mask.nii.gz']

    for pattern in patterns:
        keep_files.update(glob.glob(BIDS_folder + pattern))

    # Also keep top-level files in BIDS_folder
    keep_files.update(f for f in glob.glob(BIDS_folder + '/*') if os.path.isfile(f))

    # Show the list of files to keep
    print("\nList of files to KEEP:")
    for f in sorted(keep_files):
        print(f)

    # Ask for user confirmation before proceeding
    confirm = input("\nDo you want to continue with deletion of all other files? (yes/no): ").strip().lower()
    if confirm not in ['yes', 'y']:
        print("Operation aborted. No files were deleted.")
        return

    # Collect all files, not directories
    all_files = [f for f in glob.glob(BIDS_folder + '/**', recursive=True) if os.path.isfile(f)]

    # Determine files to delete
    to_delete = set(all_files) - keep_files

    # Deleting files
    for path in to_delete:
        try:
            os.remove(path)
            print(f"Deleted file: {path}")
        except Exception as e:
            print(f"Error deleting {path}: {e}")

    # Collect all directories, recursively
    all_dirs = [d for d in glob.glob(BIDS_folder + '/**', recursive=True) if os.path.isdir(d)]

    # Keep running the cleanup until no more empty directories remain
    while True:
        empty_dirs = [dir for dir in all_dirs if not os.listdir(dir)]  # Find empty directories

        if not empty_dirs:
            break  # Stop if no empty directories are found

        # Deleting empty directories
        for dir in empty_dirs:
            try:
                shutil.rmtree(dir)
                print(f"Deleted empty directory: {dir}")
            except Exception as e:
                print(f"Error deleting directory {dir}: {e}")

        # Refresh directory list
        all_dirs = [d for d in glob.glob(BIDS_folder + '/**', recursive=True) if os.path.isdir(d)]

    print("Cleanup complete!")

def clean_func_preserve_manual_and_topfiles(BIDS_folder):
    print("🔒 Starting cleanup of `func/` — keeping only `manual_mask.nii.gz` and top-level files")

    # Find manual masks anywhere under func/
    manual_masks = set(
        os.path.abspath(f)
        for f in glob.glob(os.path.join(BIDS_folder, '**/func/**/manual_mask.nii.gz'), recursive=True)
    )

    # Find all files directly in func/ (not in subdirectories)
    top_func_files = set()
    for dirpath in glob.glob(os.path.join(BIDS_folder, '**/func/'), recursive=True):
        if os.path.isdir(dirpath):
            files_here = glob.glob(os.path.join(dirpath, '*'))
            top_func_files.update(os.path.abspath(f) for f in files_here if os.path.isfile(f))

    keep_files = manual_masks.union(top_func_files)

    print(f"🛡️ Total files to KEEP: {len(keep_files)}")
    for f in list(keep_files)[:5]:
        print(f"  ✓ {f}")
    if len(keep_files) > 5:
        print("  ...")

    # List all files under any func/ folder
    all_func_files = set(
        os.path.abspath(f)
        for f in glob.glob(os.path.join(BIDS_folder, '**/func/**/*'), recursive=True)
        if os.path.isfile(f)
    )

    to_delete = all_func_files - keep_files

    print(f"🗑️ Total files to DELETE: {len(to_delete)}")
    for f in tqdm(to_delete, desc="Deleting unwanted func files"):
        try:
            os.remove(f)
        except Exception as e:
            print(f"⚠️ Error deleting {f}: {e}")

    print("✅ Cleanup done.")

# Example usage:
# clean_func_preserve_manual_and_topfiles("/path/to/BIDS")
