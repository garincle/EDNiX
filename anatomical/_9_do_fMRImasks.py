import os
import subprocess
import ants

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

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

def do_fMRImasks(masks_dir, labels_dir, type_norm, fMRImasks, overwrite,s_bind,afni_sif):
    # for RS analysis
    Ventri_mask = ',4,5,43,44,14'
    White_mask  = ',2,41'

    # 2) MASKS for denoising Resting State Data
    img = ants.image_read(opj(masks_dir,'brain_mask_in_anat_DC.nii.gz'))
    dilate = ants.morphology(img, operation='dilate', radius=2, mtype='binary', shape='ball')
    ants.image_write(dilate, opj(masks_dir, type_norm + 'brainmask_dilat.nii.gz'), ri=False)

    if fMRImasks=='aseg':
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "amongst(a' + Ventri_mask + ')" -prefix ' + opj(masks_dir, type_norm + 'Vmask.nii.gz')
        spco([command], shell=True)
        img = ants.image_read(opj(masks_dir, type_norm + 'Vmask.nii.gz'))
        eroded = ants.morphology(img, operation='erode', radius=1, mtype='binary', shape='ball')
        ants.image_write(eroded, opj(masks_dir, type_norm + 'Vmask_erod.nii.gz'), ri=False)
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "amongst(a' + White_mask + ')" -prefix ' + opj(masks_dir, type_norm + 'Wmask.nii.gz')
        spco([command], shell=True)
        img = ants.image_read(opj(masks_dir, type_norm + 'Wmask.nii.gz'))
        eroded = ants.morphology(img, operation='erode', radius=1, mtype='binary', shape='ball')
        ants.image_write(eroded, opj(masks_dir, type_norm + 'Wmask_erod.nii.gz'), ri=False)

    elif fMRImasks == 'custom':
        if ope(opj(masks_dir, type_norm + 'Vmask.nii.gz')):
            img = ants.image_read(opj(masks_dir, type_norm + 'Vmask.nii.gz'))
            eroded = ants.morphology(img, operation='erode', radius=1, mtype='binary', shape='ball')
            ants.image_write(eroded, opj(masks_dir, type_norm + 'Vmask_erod.nii.gz'), ri=False)
        else:
            print(bcolors.WARNING + 'We have not found the Ventricular mask, it might cause problem if you want to extract/regress the ventricular signal' + bcolors.ENDC)
        if ope(opj(masks_dir, type_norm + 'Wmask.nii.gz')):
            img = ants.image_read(opj(masks_dir, type_norm + 'Wmask.nii.gz'))
            eroded = ants.morphology(img, operation='erode', radius=1, mtype='binary', shape='ball')
            ants.image_write(eroded, opj(masks_dir, type_norm + 'Wmask_erod.nii.gz'), ri=False)
        else:
            print(bcolors.WARNING + 'We have not found the white matter mask, it might cause problem if you want to extract/regress the white matter signal' + bcolors.ENDC)
        if ope(opj(masks_dir, type_norm + 'Gmask.nii.gz')):
            print(bcolors.OKGREEN + 'We have foud a Gray mask' + bcolors.ENDC)
        else:
            print(bcolors.WARNING + 'We have not found the Gray mask, it might cause problem if you want to extract the Gray signal specifically' + bcolors.ENDC)
    else:
        raise Exception(bcolors.FAIL + 'fMRImasks must be custom or aseg' + bcolors.ENDC)