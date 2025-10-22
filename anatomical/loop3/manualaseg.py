import os
import json
import ants
import math

ope = os.path.exists
opj = os.path.join

from Tools import run_cmd
from atlases import correctaseg


def do(input,output,anatfile,preftool,sing_itk,sing_fs,diary_file,FS_refs):

    Aseg_lut = opj(FS_refs, 'FreeSurferAllLut.txt')
    itk_lut  = opj(FS_refs, 'FreeSurferAllLut_ITK.txt')

    manualaseg = output.replace('.nii.gz','_desc-manual.nii.gz')

    nl = 'WARNING1: any modifications should be save as : ' + manualaseg
    run_cmd.msg(nl, diary_file, 'WARNING')

    nl = 'WARNING2: These modifications will automatically be taken as "final mask", delete this file if you want to use a Skulltrip method!!'
    run_cmd.msg(nl, diary_file, 'WARNING')

    if preftool == 'itksnap':
        command = (sing_itk + 'itksnap -g ' + anatfile + ' -s ' + input + '-l' + itk_lut)
        nl_final = run_cmd.wait(command, diary_file)

    elif preftool == 'freeview':
        cmd = (sing_fs + 'freeview -v ' + anatfile +
               ' ' + input +':visible=0:colormap=lut:lut=' + Aseg_lut + ':opacity=0.5')
        run_cmd.do(cmd, diary_file)

    if ope(manualaseg):
        hd_BRAIN = ants.image_header_info(anatfile)
        midsection = math.ceil(int(hd_BRAIN['dimensions'][0]) / 2)
        correctaseg.smooth(manualaseg, output, 'lr', midsection, diary_file)

    dictionary = {"Sources": input,
                  "Description": 'Brain aseg (manual drawing with ' + preftool + ').', }
    json_object = json.dumps(dictionary, indent=2)
    with open(output.replace('.nii.gz' + '.json'), "w") as outfile:
        outfile.write(json_object)