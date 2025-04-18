################################################
### convert nifi file in mgz for freesurfer ####
################################################
import os
import subprocess
import datetime
import re

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

def nii_to_mgz(ID, Session, FS_dir, Ref_file, labels_dir, volumes_dir, otheranat, IgotbothT1T2, type_norm, overwrite, s_bind,fs_sif,diary_file):
    ct = datetime.datetime.now()
    nl = 'Run anatomical._10_nii_to_mgz.nii_to_mgz'
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{nl}')

    animal_folder =   'sub-' + ID + '_ses-' + str(Session)

    def get_reorient(orient_code):
        """Generate reorientation parameters based on the given orientation code."""
        mapping = {'R': -1, 'L': 1, 'A': 3, 'P': -3, 'S': -2, 'I': 2}
        # Convert orientation code to reorient parameters
        try:
            return f" -r {' '.join(str(mapping[axis]) for axis in orient_code)} "
        except KeyError:
            return None  # Invalid orientation code

    # Construct the command
    cmd = f'singularity run {s_bind}{fs_sif}mri_info --orientation {Ref_file}'
    orient_raw = subprocess.getoutput(cmd)

    # Extract the orientation code from the output
    match = re.search(r'([RLAPSI]{3})\b', orient_raw)  # Extract last valid orientation code
    orient_code = match.group(1) if match else None

    # Get reorientation parameters
    reorient = get_reorient(orient_code) if orient_code else None

    # Handle invalid or missing orientation
    if reorient is None:
        error_msg = f'ERROR: No valid orientation code found in the output: {orient_raw}'
        raise ValueError(f"{error_msg}")
    else:
        print(f"Valid original rientation: {orient_raw}")
        print(f"New fs orientation: {reorient}")

    fwdFS_cmd = ' --in_orientation ' + orient_raw + reorient

    if ope(opj(FS_dir, animal_folder)) == False:
        os.makedirs(opj(FS_dir, animal_folder))
        os.makedirs(opj(FS_dir, animal_folder, 'mri'))
        os.makedirs(opj(FS_dir, animal_folder, 'surf'))
        os.makedirs(opj(FS_dir, animal_folder, 'stats'))
        os.makedirs(opj(FS_dir, animal_folder, 'label'))
        os.makedirs(opj(FS_dir, animal_folder, 'scripts'))

    command = 'singularity run' + s_bind + fs_sif + 'mri_convert ' + fwdFS_cmd + ' ' + Ref_file + ' ' + opj(FS_dir,animal_folder, 'mri','orig.mgz') + \
    ';singularity run' + s_bind + fs_sif + 'mri_convert ' + fwdFS_cmd + ' ' + Ref_file.replace('.nii.gz', '_norm_' + type_norm + '.nii.gz') + ' ' + opj(FS_dir,animal_folder,'mri','brain.mgz') + \
    ';singularity run' + s_bind + fs_sif + 'mri_convert ' + fwdFS_cmd + ' -odt uchar ' + opj(labels_dir, type_norm + 'wm.nii.gz') + ' ' + opj(FS_dir,animal_folder,'mri','wm.seg.mgz') + \
    ';singularity run' + s_bind + fs_sif + 'mri_convert ' + fwdFS_cmd + ' ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' ' + opj(FS_dir,animal_folder,'mri','aseg.mgz') + \
    ';singularity run' + s_bind + fs_sif + 'mri_convert ' + fwdFS_cmd + ' ' + opj(labels_dir, type_norm + 'filled.nii.gz') + ' '  + opj(FS_dir,animal_folder,'mri','filled.mgz')
    nl = spgo(command)
    diary.write(f'\n{nl}')
    print(nl)


    if IgotbothT1T2 ==True:
        otheranatimg  = opj(volumes_dir, ID + otheranat + '_brain.nii.gz')
        if ope(otheranatimg):
            ####### attention!! change LPS 
            command = 'singularity run' + s_bind + fs_sif + 'mri_convert ' + fwdFS_cmd + ' ' + otheranatimg + ' ' + opj(FS_dir,animal_folder, 'mri', otheranat + 'brain.mgz')
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
        else:
            IgotbothT1T2==False
            nl = 'Warning No ' + otheranat + ' found !!!!!!!!!!!! (you said yes!!)'
            diary.write(f'\n{nl}')
            raise Exception(bcolors.FAIL + nl + bcolors.ENDC)

    diary.write(f'\n')
    diary.close()