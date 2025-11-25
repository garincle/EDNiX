import os
import ants
from Tools import run_cmd

def read(file,diary_name,sing_fs):

    cmd = (sing_fs + 'mri_convert ' + file + ' ' + file.replace('.mgz','.nii.gz'))
    run_cmd.run(cmd, diary_name)
    var = ants.image_read(file.replace('.mgz','.nii.gz'))
    os.remove(file.replace('.mgz','.nii.gz'))
    return var

def write(var,file,diary_name,sing_fs):

    ants.image_write(var,file + '.nii.gz')
    cmd = (sing_fs + 'mri_convert ' + file + '.nii.gz ' + file + '.mgz')
    run_cmd.run(cmd, diary_name)
    os.remove(file +'.nii.gz')

def readHD(file,diary_name,sing_fs):

    cmd = (sing_fs + 'mri_convert ' + file + ' ' + file.replace('.mgz','.nii.gz'))
    run_cmd.run(cmd, diary_name)
    var = ants.image_header_info(file.replace('.mgz','.nii.gz'))
    os.remove(file.replace('.mgz','.nii.gz'))
    return var
