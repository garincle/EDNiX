# import
import os
import numpy as np
import ants
import json
import nibabel as nib

opj = os.path.join
opi = os.path.isfile

from Tools import run_cmd, get_orientation

volume_ref = 107000

def get(animal,dir,file,brain_mask,diary_name,sing_fs):

    nl = 'Check if the brain is not too small to fit the freesurfer requirement'
    run_cmd.msg(nl, diary_name, 'HEADER')

    datfile = ''

    orient_raw, _, _, _ ,_= get_orientation.use_ants(file, sing_fs)
    
    img_hd     = ants.image_header_info(file)
    voxX = img_hd['spacing'][0]

    if brain_mask == '':
        img = ants.image_read(file)
        img[img>0]=1
        img_size = img.sum()
    else:
        msk = ants.image_read(brain_mask)
        img_size = msk.sum()
    volume = img_size * np.power(voxX, 3)

    change_hd = 0
    spacing = [0]
    resamp = 0

    # Load the image directly
    img = nib.load(file)
    # Get voxel sizes
    new_size, delta_y, delta_z = [str(round(abs(x), 10)) for x in img.header.get_zooms()[:3]]
    if volume < 3000:
        nl = ('That brain is too small. The header of the nifty file to be used by freessurfer will have the size and the resolution modified.'
            'hence a dat file is created to be able to reconstruct the surfaces later on')
        run_cmd.msg(nl, diary_name, 'OKGREEN')

        datfile = opj(dir, animal + '_rescale.dat')

        change_hd = 1
        new_voxsize = int(pow(volume_ref / img_size, 1 / 3) * 100) / 100
        spacing = [new_voxsize, new_voxsize, new_voxsize]
        resamp = 0
        pos = [0, 0, 0]

        if img_size > 2000:
            resamp = 1
            new_size = int(new_voxsize * 20) / 10
            if orient_raw[0] == 'R':
                pos[0] = new_size
            if orient_raw[1] == 'A':
                pos[1] = 0 #  or new_size : need to check
            elif orient_raw[1] == 'P':
                pos[1] = new_size * -1
            if orient_raw[2] == 'S':
                pos[2] = new_size

        nl = 'Creation of the ' + animal + ' dat file'
        run_cmd.msg(nl, diary_name, 'OKGREEN')

        dat = open(datfile, 'w')
        dat.write(animal + '\n')
        dat.write(str(voxX) + '\n')
        dat.write(str(voxX) + '\n')
        dat.write('0.15\n')
        dat.write(str(spacing[0] / voxX) + ' 0 0 ' + str(pos[0]) + ' \n')
        dat.write('0 ' + str(spacing[0] / voxX) + ' 0 ' + str(pos[1]) + ' \n')
        dat.write('0 0 ' + str(spacing[0] / voxX) + ' ' + str(pos[2]) + ' \n')
        dat.write('0 0 0 1\n')
        dat.write('round\n')
        dat.close()

        dictionary = {"change_hd": change_hd,
                      "spacing": spacing,
                      "resamp": resamp,
                      "new_size": new_size,
                      "datfile": datfile}

        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir, animal + '_rescale.json'), "w") as outfile:
            outfile.write(json_object)

    return change_hd,spacing,resamp,new_size,datfile

def check(jsonfile):

    change_hd = 0
    spacing = ''
    resamp = ''
    new_size = ''
    datfile = ''

    if opi(jsonfile):
        f = open(jsonfile)
        info = json.load(f)
        change_hd = info["change_hd"]
        spacing   = info["spacing"]
        resamp    = info["resamp"]
        new_size  = info["new_size"]
        datfile   = info["datfile"]

    return change_hd,spacing,resamp,new_size,datfile

