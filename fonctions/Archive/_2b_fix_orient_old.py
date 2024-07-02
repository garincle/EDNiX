###################################################################################################
######################################  fix header orient #########################################
###################################################################################################
import os
import subprocess
import glob
import shutil
import sys


#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

def fix_orient(anat_func_same_space, imgO, imgI, dir_prepro, IgotbothT1T2, path_anat, otheranat, dir_fMRI_Refth_RS_prepro1, root_RS, ID, Session, deoblique_exeption1, deoblique_exeption2, deoblique, orientation, overwrite):

    if IgotbothT1T2 == True:
        command = '3dinfo -same_obl ' + opj(path_anat, ID + 'template_indiv' + otheranat + '.nii.gz') + ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI)
        nx = spco([command], shell=True)
    else:
        nx="b'1\\n1\\n'" #code to say "same oblique"...

    debolique_spe = ID + 'ses-' + str(Session)
    if debolique_spe in deoblique_exeption1:
        deoblique_1='exeption1'
    elif debolique_spe in deoblique_exeption2:
        deoblique_1='exeption2'
    else:
        deoblique_1=deoblique

    ##do not work !!!!!!
    if deoblique_1=='header':
        if str(nx)=="b'1\\n1\\n'":
            print('same oblique')
            command = '3drefit -deoblique ' + overwrite + ' -orient ' + orientation + ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI)
            spco([command], shell=True)
            # reorient the fiedls according to the json file

        else:
            print('oblique different')
            # reorient the fiedls according to the json file
            command = '3drefit -deoblique ' + overwrite + ' -orient ' + orientation + ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI)
            spco([command], shell=True)

        shutil.copyfile( opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))
        
    #work
    elif deoblique_1=='WARP':
        if str(nx)=="b'1\\n1\\n'":
            print('same oblique')
            command = '3drefit' + overwrite + ' -orient ' + orientation + ' ' +  opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgI)
            spco([command], shell=True)
            # reorient the fiedls according to the json file
            command = '3dWarp' + overwrite + ' -deoblique -wsinc5 -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO) + \
            ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI)
            spco([command], shell=True)

        else:
            print('oblique different')
            # reorient the fiedls according to the json file
            command = '3dWarp' + overwrite + ' -deoblique -wsinc5 -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO) + \
            ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI)
            spco([command], shell=True)
            command = '3drefit' + overwrite + ' -orient ' + orientation + ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO)
            spco([command], shell=True)

    elif deoblique_1=='exeption1': #do nothing
        print('exeption1')
        shutil.copyfile( opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))

    elif deoblique_1=='exeption2': #re-alineate
        print('exeption2')
        command = '3drefit' + overwrite + ' -orient ' + orientation + ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI)
        spco([command], shell=True)
        # reorient the fiedls according to the json file
        command = '3dWarp' + overwrite + ' -deoblique -wsinc5 -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO) + \
        ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI)
        spco([command], shell=True)
        print('need to realinate or use just one')


