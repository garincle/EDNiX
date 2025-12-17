###################################################################################################
######################################  fix header orient #########################################
###################################################################################################
import os
import shutil
import json

opj = os.path.join
opi = os.path.isfile
from Tools import run_cmd
from Tools import get_orientation
def fix_orient(runMean_reorient, fMRI_runMean_unwarpped, list_RS, animalP, humanP, orientation, doWARPonfunc, sing_afni, diary_file):

    nl = '##  Working on step ' + str(2) + '(function: _2b_fix_orient).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    if len(list_RS)>1:
        if len(humanP) != len(list_RS):
            if len(humanP) == 1:
                humanP = humanP*len(list_RS)
                run_cmd.msg('it is assumed that the subject position in the scanner remained the same in every image ', diary_file, 'WARNING')
            else:
                nl = 'check the parameters: ther must be a position set for each type of images. If you do not know, leave it empty'
                run_cmd.msg(nl, diary_file, 'WARNING')

        if len(animalP) != len(list_RS):
            if len(animalP) == 1:
                animalP = animalP * len(list_RS)
                run_cmd.msg('it is assumed that the subject position in the scanner remained the same in every image ',
                            diary_file, 'WARNING')
            else:
                nl = 'check the parameters: ther must be a position set for each type of images. If you do not know, leave it empty'
                run_cmd.msg(nl, diary_file, 'WARNING')

    if humanP[0] == '':
        if opi(list_RS[0].replace('.nii.gz', 'json')):
            f = open(list_RS[0].replace('.nii.gz', 'json'))
            info = json.load(f)
            try:
                humanP[0] = info['PatientPosition']
                run_cmd.msg('the subject was scanned with the following position parameter : ' + humanP[0],
                            diary_file, 'OKGREEN')
            except:
                humanP[0] = humanP[0]

    if animalP[0] == '':
        animalP[0] = 'humanlike'

    cmd = sing_afni + '3dinfo -orient ' + fMRI_runMean_unwarpped
    msg, _ = run_cmd.get(cmd, diary_file)
    orient = msg.decode("utf-8").split('\n')[-2]

    cmd = sing_afni + '3dinfo -is_oblique ' + fMRI_runMean_unwarpped
    msg, _ = run_cmd.get(cmd, diary_file)
    obli = msg.decode("utf-8").split('\n')[-2]

    print('obli = ' + obli)
    print('orient = ' + orient)

    if obli == '1' and doWARPonfunc=='WARP':
        cmd = sing_afni + '3dWarp -overwrite -deoblique -prefix ' + runMean_reorient + ' ' + fMRI_runMean_unwarpped
        print(cmd)
        run_cmd.run(cmd, diary_file)
        desc = 'Correction of the obliquity.'
        run_cmd.do(cmd, diary_file)
        cmd = sing_afni + '3dinfo -orient ' + runMean_reorient
        msg, _ = run_cmd.get(cmd, diary_file)

        orient = msg.decode("utf-8").split('\n')[-2]

    if obli == '1' and doWARPonfunc=='header':
        cmd = sing_afni + '3dcalc -overwrite -a ' + fMRI_runMean_unwarpped + ' -prefix ' + runMean_reorient + ' -expr "a"'
        run_cmd.do(cmd, diary_file)
        desc = 'Copy .'

        cmd = sing_afni + '3drefit -overwrite -deoblique -orient ' + orientation + ' ' +  runMean_reorient
        run_cmd.run(cmd, diary_file)
        desc = 'Correction of the obliquity.'

        orient = msg.decode("utf-8").split('\n')[-2]

    else:
        cmd = sing_afni + '3dcalc -overwrite -a ' + fMRI_runMean_unwarpped + ' -prefix ' + runMean_reorient + ' -expr "a"'
        run_cmd.do(cmd, diary_file)
        desc = 'Copy .'

    if not orientation == '':
        cmd = sing_afni + '3drefit -overwrite -orient ' + orientation + ' ' + runMean_reorient
        run_cmd.run(cmd, diary_file)

    else:
        if not animalP == 'humanlike':
            neworient = get_orientation.getreal(humanP, animalP, orient)
            run_cmd.msg('the new orientation will be : ' + neworient, diary_file, 'OKGREEN')

            cmd = sing_afni + '3drefit -overwrite -orient ' + neworient + ' ' + runMean_reorient
            run_cmd.run(cmd, diary_file)
            desc = 'Correction of the fields orientation.'

    dictionary = {"Sources": fMRI_runMean_unwarpped,
                  "Description": desc, }
    json_object = json.dumps(dictionary, indent=2)
    with open(runMean_reorient.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)

