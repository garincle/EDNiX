################################################
### convert nifi file in mgz for freesurfer ####
################################################
import os
import subprocess
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

def nii_to_mgz(ID, Session, FS_dir, Ref_file, labels_dir, volumes_dir, otheranat, IgotbothT1T2, type_norm, overwrite,s_bind,fs_sif):
    animal_folder =   'sub-' + ID + '_ses-' + str(Session)

    cmd = 'singularity run' + s_bind + fs_sif + 'mri_info --orientation ' + Ref_file
    orient_raw = spgo(cmd)

    # purpouse for FS : LIA
    if orient_raw   == 'RAS': reorient = ' -r -1 3 -2 '
    elif orient_raw == 'LAS': reorient = ' -r 1 3 -2 '
    elif orient_raw == 'RSA': reorient = ' -r -1 -2 3'
    elif orient_raw == 'LSA': reorient = ' -r 1 -2 3 '
    elif orient_raw == 'RIA': reorient = ' -r -1 2 3 '
    elif orient_raw == 'LIA': reorient = ' -r 1 2 3 '
    elif orient_raw == 'LAI': reorient = ' -r 1 3 2 '
    elif orient_raw == 'RAI': reorient = ' -r -1 3 2 '
    elif orient_raw == 'RSP': reorient = ' -r -1 -2 -3 '
    elif orient_raw == 'LSP': reorient = ' -r 1 -2 -3 '
    elif orient_raw == 'RPS': reorient = ' -r -1 -3 -2 '
    elif orient_raw == 'LPS': reorient = ' -r 1 -3 -2 '
    else:
        print(bcolors.FAIL + 'ERROR: ' + orient_raw + ' not found in the list' + bcolors.ENDC)

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
    spco([command], shell=True)


    if IgotbothT1T2 ==True:
        otheranatimg  = opj(volumes_dir, ID + otheranat + '_brain.nii.gz')
        if os.path.exists(otheranatimg):
            ####### attention!! change LPS 
            command = 'singularity run' + s_bind + fs_sif + 'mri_convert ' + fwdFS_cmd + ' ' + otheranatimg + ' ' + opj(FS_dir,animal_folder, 'mri', otheranat + 'brain.mgz')
            spco([command], shell=True)
        else:
            IgotbothT1T2==False
            raise Exception(bcolors.FAIL + 'Warning No ' + otheranat + ' found !!!!!!!!!!!! (you said yes!!)' + bcolors.ENDC)