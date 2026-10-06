#import
import os
import sys
import subprocess
import glob
import json
import shutil
import datetime
import pydicom as pdcm
import numpy as np
import fnmatch

opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

spgo = subprocess.getoutput


def path(Major_path,study,Subject_Dir,Sess_to_process,tracor,nbInj):

    N = Subject_Dir.split('-')
    Subname = '-'.join([N[0], N[1]])
    
    Sess      = Sess_to_process.split('-')
    print(Sess)

    type_sess = Sess[2].lower()

    if len(Sess)==5:
        loca_sess = Sess[3].lower()
    
    pet_raw_dcm   =''
    cond          =''
    data_path_raw = ''
    dir_dcm       = ''


    if study == 'eatingcircuit':
    
        data_path_raw   = opj(Major_path, study,'BIDS','sourcedata', Subject_Dir, Sess_to_process)
        
        # get the correct path for the data i nDICOM format
        
        level1 = '_'.join([N[1],N[2]])
        
        if not tracor.lower() =='raclo':
            level3 = '_'.join([tracor.lower(),'new','pet'])
        else:
            level3 = '_'.join([tracor.lower()+'pride','new','pet'])

        if type_sess == 'control':
            level2 = type_sess
            name   = '_'.join([N[1],Sess[3],Sess[1],'CTL'])

        elif type_sess == 'flx':
            level2 = '_'.join([Sess[2],Sess[3]])
            name   = '_'.join([N[1],Sess[4],Sess[1],type_sess.upper()])

        elif type_sess == 'way':
            if loca_sess == 'sys':
                level2 = '_'.join([Sess[2],Sess[3]])
                name = '_'.join([N[1],Sess[4],Sess[1],type_sess.upper()])
            else:
                level2 = '_'.join([Sess[2],'loc',Sess[3]+nbInj])
                name = '_'.join([N[1],Sess[4],Sess[1],'-'.join([type_sess.upper(),Sess[3].upper(), 'INJ'+nbInj])])

        elif type_sess == 'sb':
            if loca_sess == 'sys':
                level2 = '_'.join([Sess[2],Sess[3]])
                name = '_'.join([N[1],Sess[4],Sess[1],type_sess.upper()])
            else:
                level2 = '_'.join([Sess[2],'loc',Sess[3]+nbInj])
                name = '_'.join([N[1],Sess[4],Sess[1],'-'.join([type_sess.upper(),Sess[3].upper(), 'INJ'+nbInj])])

        elif type_sess == 'ari':
    
            level2 = '_'.join([Sess[2],Sess[3]])
            name   = '_'.join([N[1],Sess[4],Sess[1],'ARP'])
        
        elif type_sess == 'bicu':
            level2 = '_'.join([Sess[2],'loc',Sess[3]+nbInj])  
            
            if Sess[3] == 'cdn':
                inj_site = 'vCdN'
            elif Sess[3] == 'vs':
                inj_site = 'vSV'
            name = '_'.join([N[1],Sess[4],Sess[1],'-'.join(['Bicu', inj_site+ nbInj])])
        
        elif type_sess == 'methyl':
            level2 = '_'.join([Sess[2],'loc',Sess[3]+nbInj])  
            name = '_'.join([N[1],Sess[4],Sess[1],'MPH',Sess[3].upper()+ nbInj])
        
        elif type_sess == 'ppx':
            level2 = '_'.join([Sess[2],'loc',Sess[3]+nbInj])  
            name = '_'.join([N[1],Sess[4],Sess[1],'MPH',Sess[3].upper()+ nbInj])

        
        pet_raw_dcm = opj(Major_path,study ,'dicom',level1,level2,level3,name)
        #print(pet_raw_dcm)
        dir_dcm = sorted(os.listdir(pet_raw_dcm))
        dir_dcm = opj(pet_raw_dcm,dir_dcm[-1])
    
    elif study == 'imagina':
        data_path_raw   = opj(Major_path, study,'BIDS', 'rawdata', Subject_Dir, Sess_to_process)
    
        if type_sess[:-1] == 'control':
            cond = Sess[0] + '-baseline_' + type_sess[-1]  
        else: 
            if len(loca_sess)>1:
                cond = Sess[0] + '-' + loca_sess[0] + '_' + loca_sess[1][:-1] + '_' + loca_sess[1][-1] 
        
        pet_raw_dcm = opj(Major_path,study,'BIDS','sourcedata', Subname, cond ,'DICOM','PET_FDG')
        dir_dcm = opj(pet_raw_dcm,'PET_DATA_dyn')


    pet_raw_dir = opj(data_path_raw,'pet')
    print(dir_dcm)
    print(pet_raw_dir)
    
    return pet_raw_dir,dir_dcm