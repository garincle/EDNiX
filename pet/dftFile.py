#import
import sys
import pandas
import shutil
import numpy as np
from statistics import mean

sys.path.insert(1, '../tools')
import run_cmd

def roi(file,data,nb_vox,list_rois,list_nb,Study_n,FT,FD,diary_name):
    
    nl = 'create the ROIs.dtf'
    run_cmd.msg(nl, diary_name)

    B = data.to_numpy()
    A = nb_vox.to_numpy()

    dft = open(file, "w")
    dft.write('DFT')

    for i in range(len(list_rois)):
        dft.write(f'\t{list_nb[i]}')
    ct = Study_n
    dft.write(f'\n{ct}')
    for i in range(len(list_rois)):
        ct = '.'
        dft.write(f'\t{ct}')
    ct = 'kBq/mL'
    dft.write(f'\n{ct}')
    ct = 'All'
    for i in range(len(list_rois)):
        dft.write(f'\t{ct}')
    ct = 'Times (sec)'
    dft.write(f'\n{ct}')
    for i in range(len(list_rois)):
        ct = str(B[0][i])
        dft.write(f'\t{ct}')
    for j in range(len(FT)):
        dft.write('\n' + str(FT[j]) + '\t' + str(int(FT[j]) + int(FD[j])))
        for i in range(len(list_rois)):
            ct = str(A[j][i])
            dft.write(f'\t{ct}')
    dft.write(f'\n')
    dft.close()
    shutil.move(file, file.replace('.txt','.dft'))


def ref(file,data,nb_vox,REF_name,L_index,Study_n,FT,FD,diary_name):
    
    nl = 'create the "ref".dtf'
    run_cmd.msg(nl, diary_name)

    B = data.to_numpy()
    A = nb_vox.to_numpy()

    dft = open(file, "w")
    dft.write('DFT\t' + REF_name)
    ct = Study_n
    dft.write(f'\n{ct}')
    ct = '.'
    dft.write(f'\t{ct}')
    ct = 'kBq/mL'
    dft.write(f'\n{ct}')
    ct = '.'
    dft.write(f'\t{ct}')
    ct = 'Times (sec)'
    dft.write(f'\n{ct}')
    ct = str(int(B[0][L_index[0]]) + int(B[0][L_index[1]]))
    dft.write(f'\t{ct}')
    for j in range(len(FT)):
        dft.write('\n' + str(FT[j]) + '\t' + str(int(FT[j]) + int(FD[j])))
        ct = str(mean([float(A[j][L_index[0]]), float(A[j][L_index[1]])]))
        dft.write(f'\t{ct}')
    dft.write(f'\n')
    dft.close()
    shutil.move(file, file.replace('.txt','.dft'))

