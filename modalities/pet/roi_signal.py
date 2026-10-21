#import
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ants
import json

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd,roiStats
import dftFile
from turkutools import turkubfbp,turkustrm,turkudv

dftname = 'Kinetic'

def read(img,seg,atlas_ctab,Do_cut,duration,diary_name,sing_afni):

    # check up 
    for fneeded in [img,seg]:
        if not opi(fneeded):
            raise FileNotFoundError(f"Input file not found: {fneeded}")

    nl = 'extract the mean ROIs value '
    run_cmd.msg(nl, diary_name,'OKGREEN')

    Dir_pet  = opd(img)
    Roi_dir  = opj(Dir_pet,'ROI_analysis')
    Pict_dir = opj(Roi_dir,'pictures')
    
    pictname = '_'.join(['time','activity','curve'])
    meanpict = '_'.join(['mean',pictname])

    # Extract the ROIs data ..............................................................................................................................

    if ope(Roi_dir) == False:
        os.makedirs(Roi_dir)

    for option,source,suffix in zip(['-nzvolume -nomeanout ',''],[seg,img],['vol',dftname]):
        roiStats.get(source,Roi_dir,seg,atlas_ctab,option,suffix,sing_afni,diary_name)

    im_head = ants.image_header_info(img)

    if len(im_head['dimensions'])>3:

        if ope(Pict_dir) == False:
            os.makedirs(Pict_dir)

        jsonfile = opj(Dir_pet ,'_'.join(opb(img).split('_')[0:2] + ['desc-frametiming','pet.json']))
        f = open(jsonfile)
        info_pet = json.load(f)

        data = pd.read_csv(opj(Roi_dir, opb(img.replace('.nii.gz','_'.join(['',dftname,'ROIs.txt'])))), sep='\t', header=None)

        nii = ants.image_read(seg)
        dummy = nii.unique()
        dummy = list(map(int, dummy[1:]))

        ROIS         = pd.read_csv(atlas_ctab, sep='\t')
        ROIS_present = ROIS[ROIS['ROIs'].isin(dummy)]
        list_rois    = list(ROIS_present['labels'])

        M = data.mean(axis=1)
        St = data.std(axis=1)
        max_val = M.max()+St.max()

        x = np.zeros(len(info_pet["FrameReferenceTime"]))
        for l in range(len(info_pet["FrameReferenceTime"])):
            x[l] = info_pet["FrameReferenceTime"][l]/60 # in min

        fig, ax = plt.subplots()

        for i in range(len(list_rois)):
            y=data[i]
            ax.plot(x, y, 'o-', linewidth=1, label=list_rois[i])
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        # fig.legend(loc='upper right')
        
        for p in ['png','svg']:
            plt.savefig(opj(Pict_dir, '.'.join([pictname,p])), format=p)
        plt.close()

        plt.figure()
        plt.errorbar(x,M,yerr=St,fmt='o-',linewidth=3,color='red')

        if Do_cut==1:
            for j in duration:
                x_val = j/60
                plt.plot([x_val,x_val],[0,max_val],color='black')

        for p in ['png', 'svg']:
            plt.savefig(opj(Pict_dir, '.'.join([meanpict, p])), format=p)
        plt.close()


def dft(img,seg,Roi_vol,Roi_kin,atlas_ctab,mskfile,
        REF,DV_time,Study_n,isotope,
        diary_name,sing_turku):

    # 1 # set up :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # check up ........................................................................................................

    for fneeded in [img,seg]:
        if not opi(fneeded):
            raise FileNotFoundError(f"Input file not found: {fneeded}")

    # get variables ...................................................................................................

    Dir_pet = opd(img)
    #fname= filename.get_typename(opb(img))

    # new folders .....................................................................................................
    Roi_dir  = opj(Dir_pet, 'ROI_analysis')
    pict_dir = opj(Roi_dir, 'pictures')
    if not ope(Roi_dir):
        os.makedirs(Roi_dir)
    if ope(pict_dir) == False:
        os.makedirs(pict_dir)

    jsonfile = opj(Dir_pet, '_'.join(opb(img).split('_')[0:2] + ['desc-frametiming', 'pet.json']))
    f = open(jsonfile)
    info_pet = json.load(f)

    ## files to be created .................................................................................................
    # configuration files
    ROIdft     = opj(Roi_dir, '_'.join([opb(img).split('.')[0], dftname,'ROIs','ref.dft']))

    ROIdft_tmp2 = ROIdft.replace('_ref.dft','.dft')
    ROIdft_tmp1 = ROIdft_tmp2.replace('.dft','.txt')

    Refdft = opj(Roi_dir, '_'.join([opb(img).split('.')[0], REF['name'], dftname + '.dft']))


    nl = 'extract the mean ROIs value'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    #data,nb_vox,list_rois,list_nb = niiRoi.check3(seg,atlas_ctab,'\t',Roi_vol,'\t',Roi_kin,'\t',diary_name)

    nb_vox = pd.read_csv(Roi_vol, sep='\t',header=None)
    nb_vox = nb_vox.dropna(axis=1, how="all")
    data   = pd.read_csv(Roi_kin, sep='\t',header=None)
    data = data.dropna(axis=1, how="all")


    nii = ants.image_read(seg)
    dummy = nii.unique()
    dummy = list(map(int, dummy[1:]))

    ROIS = pd.read_csv(atlas_ctab, sep='\t')
    ROIS_present = ROIS[ROIS['ROIs'].isin(dummy)]
    list_rois    = list(ROIS_present['labels'])
    list_nb      = list(ROIS_present['ROIs'])

    # create the ROIs.dtf ..............................................................................................
    dftFile.roi(ROIdft_tmp1,data,nb_vox,list_rois,list_nb,Study_n,info_pet["FrameTimesStart"],info_pet["FrameDuration"],diary_name)

    # create the ref.dtf
    if REF['type'] == 'atlas':
        L_index = []
        for z in REF["rois"]:
            L_index.append(list_nb.index(z))
        dftFile.ref(Refdft.replace('.dft', '.txt'),data,nb_vox,REF['name'],L_index,Study_n,info_pet["FrameTimesStart"],info_pet["FrameDuration"],diary_name)
    elif REF['type'] == 'mask':
        dftFile.mask(Refdft.replace('.dft', '.txt'),img,mskfile,REF['name'],Study_n,info_pet["FrameTimesStart"],info_pet["FrameDuration"],diary_name)

    # combine the info from both files
    lst = []
    for i in range(len(list_nb) + 2):
        lst.append(i)

    file1 = pd.read_csv(ROIdft_tmp2, sep='\t', header=None, names=lst)
    file2 = pd.read_csv(Refdft, sep='\t', header=None, names=['0', '1', '2', '3'])

    file1[len(list_nb) + 3] = file2['2']
    col = len(list_nb) + 1
    file1[col] = file1[col].astype(object)
    file1.loc[0:3, col] = [REF['name'],'.','All',file2['1'][3]]

    file1.to_csv(ROIdft, sep='\t', header=None, index=None)


    # 2 execution ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # 2.1 Estimate the binding potential from simplified reference tissue model'.
    nl = '# Estimate the binding potential from simplified reference tissue model'
    run_cmd.msg(nl, diary_name,'OKGREEN')
    turkustrm.bfm(ROIdft, img.replace('.nii.gz', '.sif'), list_rois, isotope, REF['name'], 1, diary_name, sing_turku)

    # 2.2 Computation of Parametric Images
    nl = '# Computation of Parametric Images:'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    # 2.2.1 Binding potential (BP) using simplified reference tissue model (SRTM)
    # 2.2.1.1 Basis function approach
    turkubfbp.run(img, Refdft, diary_name, sing_turku)
    # 2.2.1.2 Linear approach'
    turkustrm.img(img, Refdft, diary_name, sing_turku)

    # 2.2.2 Distribution Volume Ratio (DVR)
    #  applying the multiple time graphical analysis (MTGA) approach for reversible tracer uptake
    if not DV_time == 0:
        turkudv.run(img, Refdft, DV_time, diary_name, sing_turku)

