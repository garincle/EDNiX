###################################################################################################
######################################  fix header orient #########################################
###################################################################################################
import os
import subprocess
import shutil


#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

def fix_orient(imgO, imgI, dir_fMRI_Refth_RS_prepro1, root_RS, deoblique, orientation, overwrite, s_bind,afni_sif):

    if deoblique == 'no_deoblique':  # do nothing
        print('no_deoblique')
        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))

    elif deoblique == 'WARP_without_3drefit':  # do nothing
        print('WARP_without_3drefit')
        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))

    elif deoblique == 'WARP_Gridset':  # do nothing
        print('WARP_without_3drefit')
        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))

    elif deoblique == 'WARP':
        command = 'singularity run' + s_bind + afni_sif + '3drefit ' + overwrite + ' -orient ' + orientation + ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI)
        spco([command], shell=True)
        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))

    elif deoblique == 'header':
        command = 'singularity run' + s_bind + afni_sif + '3drefit ' + overwrite + ' -deoblique -orient ' + orientation + ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI)
        spco([command], shell=True)
        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))

    elif deoblique == 'header_WO_deob':
        command = 'singularity run' + s_bind + afni_sif + '3drefit ' + overwrite + ' -orient ' + orientation + ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI)
        spco([command], shell=True)
        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))


