###################################################################################################
######################################  fix header orient #########################################
###################################################################################################
import os
import shutil
import json

opj = os.path.join

from Tools import run_cmd

def fix_orient(imgO, imgI, dir_fMRI_Refth_RS_prepro1, root_RS, deoblique, orientation, overwrite, sing_afni,diary_file):

    nl = '##  Working on step ' + str(2) + '(function: _2a_fix_orient).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    if deoblique == 'WARP_without_3drefit' or deoblique == 'no_deoblique' or deoblique == 'deob_WO_orient':  # do nothing
        nl = 'do not reorient with 3drefit, just copy'
        run_cmd.msg(nl, diary_file, 'OKGREEN')

        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))
        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgI),
                      "Description": ' Copy', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO[:-7] + '.json'), "w") as outfile:
            outfile.write(json_object)

    elif deoblique == 'WARP' or deoblique == 'header_WO_deob' or deoblique == 'WARP_Gridset':
        command = (sing_afni + '3drefit ' + overwrite + ' -orient ' + orientation + ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI))
        run_cmd.run(command, diary_file)

        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))
        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgI),
                      "Description": ' Change header orientation (3drefit,AFNI', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO[:-7] + '.json'), "w") as outfile:
            outfile.write(json_object)

    elif deoblique == 'header':
        command = (sing_afni + '3drefit ' + overwrite + ' -deoblique -orient ' + orientation + ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI))
        run_cmd.run(command, diary_file)

        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))
        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgI),
                      "Description": ' deoblique and Change header orientation (3drefit,AFNI', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO[:-7] + '.json'), "w") as outfile:
            outfile.write(json_object)

    elif deoblique == 'WARPbaboon':
        command = (sing_afni + '3drefit ' + overwrite + ' -deoblique -duporigin ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI) +
                   '  -orient ' + orientation + ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI))
        run_cmd.run(command, diary_file)

        shutil.copyfile(opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI), opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO))
        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgI),
                      "Description": ' deoblique and Change header orientation (3drefit,AFNI', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + imgO[:-7] + '.json'), "w") as outfile:
            outfile.write(json_object)



