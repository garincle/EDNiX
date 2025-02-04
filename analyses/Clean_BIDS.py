#### clean a folder
import os
import glob
import shutil

def cleanBIDS(BIDS_folder):
    # Collecting the files you want to keep
    caca1 = glob.glob(BIDS_folder + '/**/**/anat/sub-**')
    caca2 = glob.glob(BIDS_folder + '/**/**/fmap/sub-')
    caca3 = glob.glob(BIDS_folder + '/**/**/func/sub-*')
    caca5 = glob.glob(BIDS_folder + '/**/**/anat/native/02_Wb/volumes/masks/*final_mask.nii.gz')
    caca6 = glob.glob(BIDS_folder + '/**/**/anat/native/02_Wb/volumes/masks/*final_mask_2.nii.gz')
    caca = glob.glob(BIDS_folder + '/**/**/**/**/**/manual_mask.nii.gz')
    caca4 = glob.glob(BIDS_folder + '/sty_template/studytemplate2_T2w/study_template_mask.nii.gz')
    caca7 = [f for f in glob.glob(BIDS_folder + '/**') if os.path.isfile(f)]

    caca_all = set(caca + caca1 + caca2 + caca3 + caca4 + caca5 + caca6 + caca7)

    # Collect all files, not directories
    all_files = [f for f in glob.glob(BIDS_folder + '/**', recursive=True) if os.path.isfile(f)]

    # Determine files to delete
    to_delete = set(all_files) - caca_all

    # Deleting files
    for path in to_delete:
        os.remove(path)  # Remove file
        print(f"Deleted file: {path}")

    # Collect all directories, recursively, under the BIDS_file
    all_dirs = [d for d in glob.glob(BIDS_folder + '/**', recursive=True) if os.path.isdir(d)]

    # Keep running the cleanup until no more empty directories remain
    while True:
        # Collect empty directories
        empty_dirs = [dir for dir in all_dirs if not os.listdir(dir)]  # Check if the directory is empty

        if not empty_dirs:
            break  # Stop if no empty directories are found

        # Deleting empty directories
        for dir in empty_dirs:
            shutil.rmtree(dir)  # Remove empty directory
            print(f"Deleted empty directory: {dir}")

        # Rebuild the list of all directories as some have been deleted
        all_dirs = [d for d in glob.glob(BIDS_folder + '/**', recursive=True) if os.path.isdir(d)]

    print("Empty directory cleanup complete!")