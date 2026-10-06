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
import niiRoi

import dftFile

def read(img,seg,pict,Do_cut,diary_name,sing_afni):

    Dir_pet = opd(img)
    N = opb(img).split('_')
    E = [i for i, elem in enumerate(N) if 'desc' in elem]
    Name = '_'.join(N[0:E[0]])

    nl = 'extract the mean ROIs value '
    run_cmd.msg(nl, diary_name)

    # Extract the ROIs data ..............................................................................................................................

    if ope(opj(Dir_pet, 'ROI_analysis')) == False:
        os.makedirs(opj(Dir_pet, 'ROI_analysis'))

    for option,source,suffix in zip(['-nzvolume -nomeanout -mask ','-mask '],[seg,img],['vol','Kinetic']):
        cmd = (sing_afni + '3dROIstats -quiet ' + option + seg + ' ' + source)
        Bp, _ = run_cmd.get(cmd, diary_name)
        list1D = Bp.decode("utf-8").split('\n')
        Bp1d = open(opj(Dir_pet, 'ROI_analysis', Name + '_ROIs_'+ suffix + '.txt'), 'w')
        for i in range(len(list1D)):
            line = '\t'.join(list1D[i].split()) + '\n'
            Bp1d.write(line)
        Bp1d.close()

    if pict[0] == 1:
        if ope(opj(Dir_pet, 'ROI_analysis', 'pictures')) == False:
            os.makedirs(opj(Dir_pet, 'ROI_analysis', 'pictures'))
    
        data,list_rois = niiRoi.check2(seg,opj(Dir_pet, 'ROI_analysis', Name  + '_ROIs_Kinetic.txt'),'\t',
        pict[1],'\t',diary_name)
    
        M = data.mean(axis=1)
        St = data.std(axis=1)
        max_val = M.max()+St.max()

        x = np.zeros(len(pict[2]))
        for l in range(len(pict[2])):
            x[l] = pict[2][l]/60

        fig, ax = plt.subplots()

        for i in range(len(list_rois)):
            y = data[list_rois[i]]
            ax.plot(x, y, 'o-', linewidth=1, label=list_rois[i])
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        # fig.legend(loc='upper right')

        plt.savefig(opj(Dir_pet, 'ROI_analysis', 'pictures', 'time_activity_curve.png'), format='png')
        plt.savefig(opj(Dir_pet, 'ROI_analysis', 'pictures', 'time_activity_curve.svg'), format='svg')
        plt.close()

        plt.figure()
        plt.errorbar(x,M,yerr=St,fmt='o-',linewidth=3,color='red')
        if Do_cut[0]==1:
            for j in Do_cut[1]:
                x_val = j/60
                plt.plot([x_val,x_val],[0,max_val],color='black')

        plt.savefig(opj(Dir_pet, 'ROI_analysis', 'pictures', 'mean_time_activity_curve.png'), format='png')
        plt.savefig(opj(Dir_pet, 'ROI_analysis', 'pictures', 'mean_time_activity_curve.svg'), format='svg')
        plt.close()


def dft(img,seg,Roi_vol,Roi_kin,atlas_ctab,
        REF,DV_time,Study_n,REF_name,isotope,FT,FD,
        pet_raw_dir,
        diary_name,sing_turku):

    Dir_pet = opd(img)
    N = opb(img).split('_')
    E = [i for i, elem in enumerate(N) if 'desc' in elem]
    Name = '_'.join(N[0:E[0]])

    nl = 'extract the mean ROIs value '
    run_cmd.msg(nl, diary_name)

    data,nb_vox,list_rois,list_nb = niiRoi.check3(seg,atlas_ctab,'\t',Roi_vol,'\t',Roi_kin,'\t',diary_name)

    L_index   = [list_nb.index(REF[0]), list_nb.index(REF[1])]

    B = data.to_numpy()
    A = nb_vox.to_numpy()


    # create the ROIs.dtf ...............................................................................................................................
    dft1_name = Name + '_ROIs_Kinetic'
    dftFile.roi(opj(Dir_pet, 'ROI_analysis', dft1_name + '.txt'),data,nb_vox,list_rois,list_nb,Study_n,FT,FD,diary_name)
    
    # create the ref.dtf
    dft2_name = Name + '_' + REF_name + '_Kinetic'
    dftFile.ref(opj(Dir_pet, 'ROI_analysis', dft2_name + '.txt'),data,nb_vox,REF_name,L_index,Study_n,FT,FD,diary_name)



    nl = '# Estimate the binding potential from simplified reference tissue model'
    run_cmd.msg(nl, diary_name)

    lst = []
    for i in range(len(list_nb) + 2):
        lst.append(i)
    file1 = pd.read_csv(opj(Dir_pet, 'ROI_analysis', dft1_name + '.dft'), sep='\t', header=None, names=lst)
    file2 = pd.read_csv(opj(Dir_pet, 'ROI_analysis', dft2_name + '.dft'), sep='\t', header=None,
                        names=['0', '1', '2', '3'])
    file1[len(list_nb) + 3] = file2['2']
    file1[len(list_nb) + 1][0] = REF_name
    file1[len(list_nb) + 1][1] = '.'
    file1[len(list_nb) + 1][2] = 'All'
    file1[len(list_nb) + 1][3] = int(A[0][L_index[0]]) + int(A[0][L_index[1]])
    file1.to_csv(opj(Dir_pet, 'ROI_analysis', dft1_name + '_ref.dft'), sep='\t', header=None, index=None)

    cmd = (sing_turku + 'dftweigh ' + opj(Dir_pet, 'ROI_analysis', dft1_name + '_ref.dft') + ' ' + img.replace('.nii.gz','.sif'))
    run_cmd.run(cmd,diary_name)

    cmd = (sing_turku + 'bfmsrtm -fit=' + opj(Dir_pet, 'ROI_analysis', dft1_name + '_fit.dft') +
           ' -i=' + isotope +
           ' ' + opj(Dir_pet, 'ROI_analysis', dft1_name + '_ref.dft') +
           ' ' + REF_name +
           ' ' + opj(Dir_pet, 'ROI_analysis', dft1_name + '_bp.txt'))
    run_cmd.run(cmd,diary_name)

    if ope(opj(Dir_pet, 'ROI_analysis', 'pictures')) == False:
        os.makedirs(opj(Dir_pet, 'ROI_analysis', 'pictures'))

    file0 = pd.read_csv(opj(Dir_pet, 'ROI_analysis', dft1_name + '_ref.dft'), sep='\t', header=None, skiprows=4)

    x = file0[0]
    y1 = file0[len(list_rois) + 2]
    fig, ax = plt.subplots()

    for i in range(16):
        y = file0[2 + i]
        ax.plot(x, y, 'o-', linewidth=1, label=list_rois[i])

    ax.plot(x, y1, 'ro-', linewidth=3, label='reference: ' + REF_name)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    # fig.legend(loc='upper right')

    plt.savefig(opj(Dir_pet, 'ROI_analysis', 'pictures', 'Before_fit.png'), format='png')
    plt.savefig(opj(Dir_pet, 'ROI_analysis', 'pictures', 'Before_fit.svg'), format='svg')
    plt.close()

    file1 = pd.read_csv(opj(Dir_pet, 'ROI_analysis', dft1_name + '_fit.dft'), sep='\t', header=None, skiprows=4,
                        skipfooter=2,engine='python')
    x = file1[0]
    y1 = file1[len(list_rois) + 2]

    
    for i in range(len(list_rois)):
        fig, ax = plt.subplots()
        y = file1[2 + i]
        ax.plot(x, y, 'o-', linewidth=2)
        ax.plot(x, y1, 'o-', linewidth=2, )
        ax.set_title(list_rois[i])
        plt.savefig(opj(Dir_pet, 'ROI_analysis', 'pictures', list_rois[i] + '_fit.png'), format='png')
        plt.close()

    nl = '# Computation of Parametric Images:\n'
    run_cmd.msg(nl, diary_name)

    cmd = 'gzip -d ' + img
    run_cmd.run(cmd,diary_name)

    nl = '# Computation of parametric image with binding potential (BP) using simplified reference tissue model (SRTM) with basis function approach\n'
    run_cmd.msg(nl, diary_name)
    # http://www.turkupetcentre.net/petanalysis/tpcclib/doc/imgbfbp.html

    cmd = (sing_turku + 'imgbfbp -err=' + opj(Dir_pet, Name + '_desc-BP-err_pet.nii') +
           ' ' + img.replace('.nii.gz','.nii') +
           ' ' + opj(Dir_pet, 'ROI_analysis', dft2_name + '.dft') +
           ' ' + opj(Dir_pet, Name + '_desc-BP_pet.nii') +
           ' ' + img.replace('.nii.gz','.sif'))
    run_cmd.run(cmd,diary_name)

    cmd = 'gzip ' + opj(Dir_pet, Name + '_desc-BP_pet.nii')
    run_cmd.run(cmd,diary_name)

    dictionary = {"Sources": [img.replace('.nii.gz', '.nii'),
                              opj(Dir_pet, 'ROI_analysis', dft2_name + '.dft'),
                              img.replace('.nii.gz','.sif')],
                  "Description": 'Computation of parametric image with binding potential (BP) using simplified reference tissue model (SRTM) with basis function approach (turku,tpcclib: imgbfbp).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(Dir_pet, Name + '_desc-BP_pet.json'), "w") as outfile:
        outfile.write(json_object)

    cmd = 'gzip ' + opj(Dir_pet, Name + '_desc-BP-err_pet.nii')
    run_cmd.run(cmd,diary_name)



    nl = '# Computation of parametric image with binding potential (BP) using simplified reference tissue model (SRTM) with linear approach\n'
    run_cmd.msg(nl, diary_name)
    # http://www.turkupetcentre.net/petanalysis/tpcclib/doc/bfmsrtm.html

    cmd = (sing_turku + 'imgsrtm' + ' ' + img.replace('.nii.gz','.nii') +
          ' ' + opj(Dir_pet, 'ROI_analysis', dft2_name + '.dft') +
          ' ' + opj(Dir_pet, Name + '_desc-STRM_pet.nii'))
    run_cmd.run(cmd,diary_name)

    cmd = 'gzip ' + opj(Dir_pet, Name + '_desc-STRM_pet.nii')
    run_cmd.run(cmd,diary_name)

    dictionary = {"Sources": [img.replace('.nii.gz', '.nii'),
                              opj(Dir_pet, 'ROI_analysis', dft2_name + '.dft'),
                              img.replace('.nii.gz','.sif')],
                  "Description": 'Computation of parametric image with binding potential (BP) using simplified reference tissue model (SRTM) with linear approach (turku,tpcclib: imgsrtm).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(Dir_pet, Name + '_desc-STRM_pet.json'), "w") as outfile:
        outfile.write(json_object)


    if not DV_time == '':

        nl = '# Computation of parametric image with Distribution Volume Ratio (DVR) applying the multiple time graphical analysis (MTGA) approach for reversible tracer uptake\n'
        run_cmd.msg(nl, diary_name)
        # http://www.turkupetcentre.net/petanalysis/tpcclib/doc/imgdv.html
        
        cmd = (sing_turku + 'imgdv -n=' + opj(Dir_pet, Name + '_desc-DVR_plot_pet') +
           ' ' + opj(Dir_pet, 'ROI_analysis', dft2_name + '.dft') +
           ' ' + img.replace('.nii.gz','.nii') +
           ' ' + DV_time + ' ' + opj(Dir_pet, Name + '_desc-DVR_pet.nii'))
        run_cmd.run(cmd,diary_name)

        cmd = 'gzip ' + opj(Dir_pet, Name + '_desc-DVR_plot_pet.nii')
        run_cmd.run(cmd,diary_name)

        cmd = 'gzip ' + opj(Dir_pet, Name + '_desc-DVR_pet.nii')
        run_cmd.run(cmd,diary_name)

        dictionary = {"Sources": [img.replace('.nii.gz', '.nii'),
                                opj(Dir_pet, 'ROI_analysis', dft2_name + '.dft')],
                    "Description": 'Computation of parametric image with Distribution Volume Ratio (DVR) applying the multiple time graphical analysis (MTGA) approach for reversible tracer uptake (turku,tpcclib: imgdv).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(Dir_pet, Name + '_desc-DVR_pet.json'), "w") as outfile:
            outfile.write(json_object)



    cmd = 'gzip ' + img.replace('.nii.gz','.nii')
    run_cmd.run(cmd,diary_name)
