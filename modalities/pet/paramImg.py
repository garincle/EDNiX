#import
import os
import sys
import ants
import numpy as np
import json

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname


sys.path.insert(1, 'tools')
import run_cmd

Palette = (' MODE_USER_SCALE -thresholding THRESHOLD_TYPE_OFF THRESHOLD_TEST_SHOW_INSIDE 0 2.5' +
              ' -pos-user 0 2.5 -interpolate true -palette-name ROY-BIG -disp-pos true -disp-neg false -disp-zero false')


def bratio(img1,imgSS,diary_name,sing_wb):
    
    nl =  'normalize the voxel signal by the mean brain signal'
    run_cmd.msg(nl,diary_name)

    SS = ants.image_read(imgSS)
    SUV = ants.image_clone(SS)
    meanS = SS[SS>0].mean()
    print(meanS)

    Dir_pet = opd(img1)
    N = opb(img1).split('_')
    E = [i for i, elem in enumerate(N) if 'desc' in elem]
    Name1 = '_'.join(N[0:E[0]+1])
    
    N = opb(imgSS).split('_')
    E = [i for i, elem in enumerate(N) if 'desc' in elem]
    Name2 = '_'.join(N[0:E[0]+1])

    SUV = SS/meanS
    ants.image_write(SUV, opj(Dir_pet, Name2 + '-bratio_pet.nii.gz'))
    cmd = (sing_wb + 'wb_command -volume-palette' + ' ' + opj(Dir_pet, Name2 + '-bratio_pet.nii.gz') + Palette)
    run_cmd.run(cmd,diary_name)

    dictionary = {"Sources": imgSS,
                      "Description": 'normilisation of the signal by the global mean (ANTspy)', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(Dir_pet, Name2 + '-bratio_pet.json'), "w") as outfile:
        outfile.write(json_object)

    IMG = ants.image_read(img1)
    SUV = ants.image_clone(IMG)
    SUV = IMG/meanS
    ants.image_write(SUV, opj(Dir_pet, Name1 + '-bratio_pet.nii.gz'))
    cmd = (sing_wb + 'wb_command -volume-palette' + ' ' + opj(Dir_pet, Name1 + '-bratio_pet.nii.gz') + Palette)
    run_cmd.run(cmd,diary_name)
    dictionary = {"Sources": [img1,
    imgSS],
    "Description": 'normilisation of the signal by the global mean in the brain(ANTspy)', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(Dir_pet, Name1 + '-bratio_pet.json'), "w") as outfile:
        outfile.write(json_object)

def suv(file,ratio):

    Dir_pet = opd(file)
    N = opb(file).split('_')
    E = [i for i, elem in enumerate(N) if 'desc' in elem]
    Name = '_'.join(N[0:E[0]])
    suffix = N[E[0]]
    newname = '_'.join([Name,suffix + '-SUV','pet'])
    
    IMG = ants.image_read(file)
    new_img = ants.image_clone(IMG)
    header = ants.image_header_info

    if len(headedr['dimensions'])==4:
        img_unmerged = ants.ndimage_to_list(IMG)
        img_list = list()
        for i in range(len(img_unmerged)):
            temp = img_unmerged[i] / ratio
            img_list.append(temp)

        new_img = ants.list_to_ndimage(img, img_list)
    else:
        new_img = IMG/ ratio
    
    ants.image_write(new_img, opj(Dir_pet, newname + '.nii.gz'))
    dictionary = {"Sources": [file,
    str(ratio)],
    "Description": 'Mean image ponderate by subject weight and tracor qte (3dTstat, AFNI and ANtspy).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(Dir_pet, newname + '.json'), "w") as outfile:
        outfile.write(json_object)
    

    



