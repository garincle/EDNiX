import os
import json
import ants
import numpy as np
import pandas as pd

opj = os.path.join
ope = os.path.exists
opb = os.path.basename
opd = os.path.dirname


def check_atlasjson(labelname,atlas_name,atlas_level,atlas_file,info):

    hd = ants.image_header_info(atlas_file)
    if hd['nDimensions']==3:
        atlas_level = 0

    newline=pd.DataFrame(columns=['ROIs','labels'])
    newfile=''
    levelnb  = atlas_level
    labelnb = []
    toget = info['atlas_definitions'][atlas_name + '[' + str(levelnb) + ']']['regions']
    for i in toget:
        if labelname in i['region']:
            labelnb.append(i['label'])

    if not labelnb:
        if hd['nDimensions']==3:
            print('sorry the label is not found! maybe change the name')
        else:

            nb = int(hd['dimensions'][3])
            nextvol = list(np.arange(nb))
            nextvol.remove(levelnb)
            for levelnb in nextvol:
                toget = info['atlas_definitions'][atlas_name + '[' + str(levelnb) + ']']['regions']
                for i in toget:
                    if labelname in i['region']:
                        labelnb.append(i['label'])
                        newline.loc[len(newline)] = {'ROIs': i['label'],'labels': i['region']}
                if labelnb:

                    forpetdir = opj(opd(atlas_file),'for_pet')
                    os.makedirs(forpetdir,exist_ok=True)
                    nii = ants.image_read(atlas_file)
                    img = ants.ndimage_to_list(nii)
                    for z in labelnb:
                        img[atlas_level][img[levelnb]==z]=z
                    ants.image_write(img[atlas_level], opj(forpetdir,opb(atlas_file).replace('.nii.gz','_test.nii.gz')))
                    newfile = opj(forpetdir,opb(atlas_file).replace('.nii.gz','_desc-modif4pet.nii.gz'))
                    break
            if not labelnb:
                print('sorry the label was not found! maybe change the name')
                labelnb=None
                levelnb=None

    return labelnb,levelnb,newline,newfile


