import os
import glob
import shutil

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