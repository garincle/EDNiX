#import
import os
import shutil
import ants
import numpy as np
opi = os.path.isfile

from Tools import run_cmd


def roi(file,data,nb_vox,list_rois,list_nb,Study_n,FS,FD,diary_name):
    
    nl = 'create the ROIs.dtf'
    run_cmd.msg(nl, diary_name,'OKGREEN')


    A = nb_vox.to_numpy()
    B = data.to_numpy()

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
        ct = str(A[0][i])
        dft.write(f'\t{ct}')
    for j in range(len(FS)):
        dft.write('\n' + str(FS[j]) + '\t' + str(int(FS[j]) + int(FD[j])))
        for i in range(len(list_rois)):
            ct = str(B[j][i])
            dft.write(f'\t{ct}')
    dft.write(f'\n')
    dft.close()
    shutil.move(file, file.replace('.txt','.dft'))


def ref(file,data,nb_vox,REF_name,L_index,Study_n,FS,FD,diary_name):

    nl = 'create the "ref".dtf'
    run_cmd.msg(nl, diary_name,'OKGREEN')


    A = nb_vox.to_numpy()
    tt=0
    for z in L_index:
        tt+=int(A[0][z])

    B = data.to_numpy()

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
    ct = str(tt)
    dft.write(f'\t{ct}')
    for j in range(len(FS)):
        dft.write('\n' + str(FS[j]) + '\t' + str(int(FS[j]) + int(FD[j])))
        aa = 0
        for z in L_index:
            aa += float(B[j][z])
        ct = str(aa/len(L_index))
        dft.write(f'\t{ct}')
    dft.write(f'\n')
    dft.close()
    shutil.move(file, file.replace('.txt','.dft'))

def mask(file, img,mskfile,REF_name, Study_n, FS, FD, diary_name):

        msk = ants.image_read(mskfile)
        vol = msk.sum() * np.prod(msk.spacing)

        IMG = ants.image_read(img)
        data = ants.ndimage_to_list(IMG)

        nl = 'create the "ref".dtf'
        run_cmd.msg(nl, diary_name,'OKGREEN')

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
        ct = str(vol)
        dft.write(f'\t{ct}')
        for j in range(len(FS)):
            dft.write('\n' + str(FS[j]) + '\t' + str(int(FS[j]) + int(FD[j])))
            dummy = data[j]*msk
            ct=str(dummy.sum())
            dft.write(f'\t{ct}')
        dft.write(f'\n')
        dft.close()
        shutil.move(file, file.replace('.txt', '.dft'))

