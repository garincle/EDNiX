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

def fix_orient(anat_func_same_space, imgO, imgI, dir_prepro, IgotbothT1T2, path_anat, otheranat, dir_fMRI_Refth_RS_prepro1, root_RS, ID, Session, deoblique_exeption1, deoblique_exeption2, deoblique, orientation, overwrite,
               s_bind,afni_sif):
    debolique_spe = ID + 'ses-' + str(Session)
    if debolique_spe in deoblique_exeption1:
        deoblique_1='exeption1'
        print(deoblique_1)
    elif debolique_spe in deoblique_exeption2:
        deoblique_1='exeption2'
    else:
        deoblique_1=deoblique
        print(deoblique_1)

    if deoblique_1 == 'exeption1':  # do nothing
        print('exeption1')
        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))

    elif deoblique_1 == 'exeption2':  # do nothing
        print('exeption2')
        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))

    elif deoblique_1 == 'WARP':
        command = 'singularity run' + s_bind + afni_sif + '3drefit ' + overwrite + ' -orient ' + orientation + ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI)
        spco([command], shell=True)
        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))

    elif deoblique_1 == 'header':
        command = 'singularity run' + s_bind + afni_sif + '3drefit ' + overwrite + ' -orient ' + orientation + ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI)
        spco([command], shell=True)
        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))

    elif deoblique_1 == 'header_deob':
        command = 'singularity run' + s_bind + afni_sif + '3drefit ' + overwrite + ' -deoblique -orient ' + orientation + ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI)
        spco([command], shell=True)
        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))


