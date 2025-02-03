###################################################################################################
######################################  fix header orient #########################################
###################################################################################################
import os
import subprocess
import shutil
import datetime
import json

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

def fix_orient(imgO, imgI, dir_fMRI_Refth_RS_prepro1, root_RS, deoblique, orientation, overwrite, s_bind,afni_sif,diary_file):

    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  Working on step ' + str(2) + '(function: _2a_fix_orient).  ##'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    if deoblique == 'WARP_without_3drefit' or deoblique == 'no_deoblique' or deoblique == 'deob_WO_orient':  # do nothing
        nl = 'do not reorient with 3drefit, just copy'

        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))
        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgI),
                      "Description": ' Copy', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO[:-7] + '.json'), "w") as outfile:
            outfile.write(json_object)

    elif deoblique == 'WARP' or deoblique == 'header_WO_deob' or deoblique == 'WARP_Gridset':
        command = 'singularity run' + s_bind + afni_sif + '3drefit ' + overwrite + ' -orient ' + orientation + ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI)
        nl = spgo(command)

        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))
        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgI),
                      "Description": ' Change header orientation (3drefit,AFNI', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO[:-7] + '.json'), "w") as outfile:
            outfile.write(json_object)

    elif deoblique == 'header':
        command = 'singularity run' + s_bind + afni_sif + '3drefit ' + overwrite + ' -deoblique -orient ' + orientation + ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI)
        nl = spgo(command)

        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))
        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgI),
                      "Description": ' deoblique and Change header orientation (3drefit,AFNI', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO[:-7] + '.json'), "w") as outfile:
            outfile.write(json_object)

    diary.write(f'\n{nl}')
    print(nl)
    diary.write(f'\n')
    diary.close()


