import os
import subprocess

# Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput


def extract_filename(file_path):
    # Use os.path.basename to get the base name of the file
    base_name = opb(file_path)
    # Check if the file has a ".gz" extension and remove it if present
    if base_name.endswith(".nii.gz"):
        base_name = base_name[:-7]
    # Return the desired part of the file path
    elif base_name.endswith(".nii"):
        base_name = base_name[:-4]
    # Return the desired part of the file path
    return base_name