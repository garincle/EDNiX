#import
import os
import shutil

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
opi = os.path.isfile

import run_cmd

def create(file,START,FT,FS,FD,isotope,diary_name,sing_turku):
    
    nl ='create the sif file'
    run_cmd.msg(nl,diary_name)
    
    # Weighting the Radioactivity concentration according to signal intensity, physical decay, frame duration  and voxel volume
    Study_n = 'study_01234'

    dir  = opd(file)

    if opi(file.replace('.nii.gz','.sif')) == False:
    
        SIF = open(opj(dir,'sif.txt'),"w")
        SIF.write(START + ' ' + str(len(FT)) + ' 4 1 ' + Study_n + ' ' + isotope)
        for j in range(len(FT)):
            SIF.write('\n' + str(FS[j]) + ' ' + str(int(FS[j]) + int(FD[j])) + ' 0 0')
        SIF.close()
        shutil.move(opj(dir,'sif.txt'), file.replace('.nii.gz','.sif'))
        
        cmd = 'gzip -d ' + file
        run_cmd.run(cmd,diary_name)
        
        cmd = (sing_turku + 'imgweigh ' + file.replace('.nii.gz','.nii') + ' ' + file.replace('.nii.gz','.sif'))
        run_cmd.run(cmd,diary_name)
        
        cmd = 'gzip ' + file.replace('.nii.gz','.nii')
        run_cmd.run(cmd,diary_name)
