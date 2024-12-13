import os
import subprocess
import nibabel as nb
import numpy as np
from math import pi
from fonctions.extract_filename import extract_filename
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
# Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

# here we need a fmap_dir
#

# list_fmap =

def correct_img_HCP(dir_fMRI_Refth_RS_prepro1, RS, list_map, RS_map, study_fMRI_Refth, r, overwrite,s_bind,afni_sif,fsl_sif,topup_file):
    for i in [0, 1]:
        root = extract_filename(RS_map[i])
        root_RS = extract_filename(RS[r])

        # copy map imag in new location
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + list_map[i] + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,
                                                                                  RS_map[i]) + ' -expr "a"'
        spco([command], shell=True)

        # mean of the map img
        command = 'singularity run' + s_bind + afni_sif + '3dTstat' + overwrite + ' -mean -prefix ' + \
                  opj(dir_fMRI_Refth_RS_prepro1,root + '_map_mean_pre' + str(i) + '.nii.gz') + ' ' + \
                  opj(dir_fMRI_Refth_RS_prepro1, RS_map[i])
        spco([command], shell=True)

        # register each volume to the base image
        command = 'singularity run' + s_bind + afni_sif + '3dvolreg' + overwrite + ' -verbose -zpad 1 -base ' + \
                  opj(dir_fMRI_Refth_RS_prepro1,root + '_map_mean_pre' + str(i) + '.nii.gz') + \
                  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_map_align' + str(i) + '.nii.gz') + \
                  ' -cubic ' + \
                  opj(dir_fMRI_Refth_RS_prepro1, RS_map[i])
        spco([command], shell=True)

        # mean of the map img to ref img
        command = 'singularity run' + s_bind + afni_sif + '3dTstat' + overwrite + ' -mean -prefix ' + \
                  opj(dir_fMRI_Refth_RS_prepro1,root + '_map_mean' + str(i) + '.nii.gz') + ' ' + \
                  opj(dir_fMRI_Refth_RS_prepro1, root + '_map_align' + str(i) + '.nii.gz')
        spco([command], shell=True)

    root = extract_filename(RS_map[0])
    root1 = extract_filename(RS_map[1])
    command = 'singularity run' + s_bind + afni_sif + '3dTcat' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_se.nii.gz') + \
              ' ' + opj(dir_fMRI_Refth_RS_prepro1, root1 + '_map_mean' + str(1) + '.nii.gz') + \
                ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean' + str(0) + '.nii.gz')
    spco([command], shell=True)

    #### zeropad?? add a slice instate of removing!!!
    im = nb.load(opj(dir_fMRI_Refth_RS_prepro1, root + '_se.nii.gz'))
    imdata = im.get_fdata()
    s = imdata.shape
    dests = np.array(s)
    hdr = im.header.copy()
    hdr.set_data_shape(imdata.shape)
    for b, d in enumerate(s):
        if b < 3:
            if (d % 2) == 0:
                print(bcolors.OKGREEN + "{0} est paire, no need to remove a slice" + bcolors.ENDC)
            else:
                print(bcolors.OKGREEN + "{0} est impaire, we will have to remove a slice" + bcolors.ENDC)
                imdata = imdata.take(range(d - 1), axis=b)
    nb.Nifti1Image(imdata, im.affine, hdr).to_filename(opj(dir_fMRI_Refth_RS_prepro1, root + '_se1.nii.gz'))

    ### se_map don't change but 1 -1
    ### b02b0 don't change

    command = 'singularity run' + s_bind + fsl_sif + 'topup --imain=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_se1.nii.gz') + \
              ' --datain=' + topup_file[0] + \
              ' --config=' + topup_file[1] + \
              ' --fout=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap.nii.gz') + \
              ' --iout=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_unwarped.nii.gz')
    spco([command], shell=True)

    ##### for fugue
    command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap.nii.gz') + ' -mul ' + str(2 * pi) + \
              ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_rads.nii.gz')
    spco([command], shell=True)

    command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_unwarped.nii.gz') + \
              ' -Tmean ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_mag.nii.gz')
    spco([command], shell=True)