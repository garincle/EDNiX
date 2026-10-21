#import
import os
import shutil
import glob
import json


opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
opi = os.path.isfile
ope = os.path.exist


from Tools import run_cmd
from motilities.pet.turkutools import turkusif
from motilities.pet import dcmFile

def create(file,START,FT,FS,FD,isotope,diary_name,sing_turku):
    
    nl ='create the sif file'
    run_cmd.msg(nl,diary_name,'OKGREEN')

    if not opi(file):
        raise FileNotFoundError(f"Input file not found: {file}")

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

        turkusif.weight(file, diary_name, sing_turku)

def read(file):

    SIF = open(file, "r")
    lignes = SIF.readlines()
    SIF.close()
    listn = lignes[0].split(' ')
    START = listn[0]
    isotope = listn[-1].split('\n')[0]

    return START,isotope


def check(file,dir_dcm,diary_file, endfile,sing_turku):
    # create the sif file if needed:
    if opi(file.replace('.nii.gz', '.sif')):
        START, isotope = read(file.replace('.nii.gz', '.sif'))
    else:
        jsonfile = file.replace('.nii.gz', '.json')

        # From the dicom if exist
        if ope(dir_dcm) and opi(jsonfile):
            list_dcm = sorted(glob.glob(opj(dir_dcm, '*.' + endfile)))
            isotope, START = dcmFile.extract(opj(dir_dcm, list_dcm[0]))
            f = open(jsonfile)
            info_pet = json.load(f)
            if "FrameTimesStart" in info_pet and "FrameReferenceTime" in info_pet and "FrameDuration" in info_pet:
                create(file, START, info_pet["FrameReferenceTime"], info_pet["FrameTimesStart"],
                               info_pet["FrameDuration"], isotope, diary_file, sing_turku)
        else:
            print('you must create a sif file')

