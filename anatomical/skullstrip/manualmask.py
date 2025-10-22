import os
import shutil
import json

ope = os.path.exists

from Tools import run_cmd

def do(input_for_msk,output_for_mask,sing_afni,preftool,sing_itk,sing_fs,diary_file):

    if not ope(output_for_mask):
        command = (sing_afni + '3dcalc -a ' + input_for_msk + ' -expr "step(a)" -prefix ' + output_for_mask)
        run_cmd.do(command, diary_file)

    manualmask = output_for_mask.replace('.nii.gz','_desc-manual.nii.gz')

    nl = 'WARNING1: any modifications should be save as : ' + manualmask
    run_cmd.msg(nl, diary_file, 'WARNING')

    nl = 'WARNING2: These modifications will automatically be taken as "final mask", delete this file if you want to use a Skulltrip method!!'
    run_cmd.msg(nl, diary_file, 'WARNING')

    if preftool == 'itksnap':
        command = (sing_itk + 'itksnap -g ' + input_for_msk + ' -s ' + output_for_mask)
        nl_final = run_cmd.wait(command, diary_file)

    elif preftool == 'freeview':
        cmd = (sing_fs + 'freeview -v ' + input_for_msk + ' ' + output_for_mask +
               ':colormap=heat:opacity=0.5:visible=1')
        run_cmd.do(cmd, diary_file)

    if ope(manualmask):
        shutil.copyfile(manualmask, output_for_mask)

    dictionary = {"Sources": input_for_msk,
                  "Description": 'Brain mask (manual drawing with ' + preftool + ').', }
    json_object = json.dumps(dictionary, indent=2)
    with open(output_for_mask[:-7] + '.json', "w") as outfile:
        outfile.write(json_object)