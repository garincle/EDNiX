#import
import os
import sys
import json
import shutil
import ants
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statistics import mean

opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

sys.path.insert(1, '../tools')
import run_cmd




def check2(nii,listdata,sep1,listroi,sep2,diary_name):
    
    # need to control if all the ROIs are present in the PET resolution (small nuclei such as the habenula could be missing)
    
    data   = pd.read_csv(listdata,sep=sep1, header=None)
    ROIS   = pd.read_csv(listroi, sep=sep2)
    
    list_rois = list(ROIS['labels'])
    list_nb   = list(ROIS['ROIs'])

    D = data.columns[data.isna().any()].tolist()
    print(len(D))
    
    if D:
        if D[0] == 0:
            print('need to remove the first column')
            data = data.drop(0, axis=1)

    # need to control if all the ROIs are present in the PET resolution (small nuclei such as the habenula could be missing)
    
    img = ants.image_read(nii)
    rr = np.setdiff1d(list_nb, img.unique()[1:])
    
    print(rr)
    
    if len(rr) > 0:
        for z in rr:
            e = list_rois.index(z)
            nl = 'beware : the column nb :' + str(e) + ' is missing in the segmented file'
            run_cmd.msg(nl, diary_name)
            
            data.insert(loc=e, column='new_' + str(e), value=0)

    data.columns = list_rois

    return data,list_rois

def check3(seg,atlas_ctab,sep1,Roi_vol,sep2,Roi_kin,sep3,diary_name):

    ROIS      = pd.read_csv(atlas_ctab,sep=sep1)
    list_rois = list(ROIS['labels'])
    list_nb   = list(ROIS['ROIs'])


    data   = pd.read_csv(Roi_vol, sep=sep2, header=None)
    nb_vox = pd.read_csv(Roi_kin, sep=sep3, header=None)
    
    D = data.columns[data.isna().any()].tolist()
    
    if D:
        if D[0] == 0:
            nl = 'need to remove the first column'
            run_cmd.msg(nl, diary_name)
            data   = data.drop(0, axis=1)
            nb_vox = nb_vox.drop(0, axis=1)

    # need to control if all the ROIs are present in the PET resolution (small nuclei such as the habenula could be missing)
    img = ants.image_read(seg)
    rr = np.setdiff1d(list_nb, img.unique()[1:])

    if len(rr) > 0:
        for z in rr:
            print(z)
            e  = list_nb.index(z)
            s = list_rois[e]
            print(s)
            nl = 'beware : the column nb :' + str(e) + ' is missing in the *_seg- file'
            run_cmd.msg(nl, diary_name)
            data.insert(loc=e, column='new_' + s, value=0)
            nb_vox.insert(loc=e, column='new_' + s, value=0)

    data.columns   = list_rois
    nb_vox.columns = list_rois

    return data,nb_vox,list_rois,list_nb


