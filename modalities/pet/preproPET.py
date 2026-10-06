#import
import os
import sys
import glob
import ants

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile


from Tools import run_cmd

from pet import roi_signal
from pet import cut_image
from pet import norm2template
from pet import paramImg


def FDG(file,dseg,template,T1wacpcmask,T1wacpc,templateT1w,templatemask,maxval,dir_transfo,atlas_label,
FT,do_cut,steps,frame_stop,duration,ratio,diary_name,sing_afni,sing_wb):

    im_head = ants.image_header_info(file)
    Dir_pet = opd(file)
    N = opb(file).split('_')
    E = [i for i, elem in enumerate(N) if 'desc' in elem]
    Name = '_'.join(N[0:E[0]])
    suffix = N[E[0]]

    if not ratio == '':
        paramImg.suv(file,ratio)
        file = opj(Dir_pet,'_'.join([Name,suffix + '-SUV','pet.nii.gz']))

    if len(im_head['dimensions'])>3:
        option1 = [1,atlas_label,FT]
    else : 
        option1 = [0]
    
    if do_cut[0]==0:

        nl =  'get the mean amplitude signal per each Roi'
        run_cmd.msg(nl,diary_name)
        roi_signal.read(file,dseg,option1,[0],[0],diary_name,sing_afni)

        nl =  'normalize the image to the anat and to the reference template and skull-strip it'
        run_cmd.msg(nl,diary_name)
        norm2template.pet(file,template,T1wacpcmask,T1wacpc,templateT1w,templatemask,maxval,dir_transfo,diary_name,sing_wb)
        
        for x,z in zip(['acpc',opj('templates', template)],['acpc',template]):
            img1 = opj(Dir_pet,x,'_'.join([Name, 'space-' + z, suffix      ,'pet.nii.gz']))
            imgSS= opj(Dir_pet,x,'_'.join([Name, 'space-' + z, suffix+'-SS','pet.nii.gz']))
            paramImg.bratio(img1,imgSS,diary_name,sing_wb)

    else:
        nl =  'get the mean amplitude signal per each Roi'
        run_cmd.msg(nl,diary_name)
        roi_signal.read(file,dseg,option1,[do_cut[0],steps[0]],diary_name,sing_afni)
        nl =  'normalize the image to the anat and to the reference template and skull-strip it'
        run_cmd.msg(nl,diary_name)
        norm2template.pet(file,template,T1wacpcmask,T1wacpc,templateT1w,templatemask,maxval,dir_transfo,diary_name,sing_wb)
    
        nl =  'mean by each time frames'
        run_cmd.msg(nl,diary_name)
        cut_image.doit(file,steps[0],frame_stop[0],duration,ratio,diary_name,sing_afni)
        E = -1
        N = opb(file).split('_')
        E = [i for i, elem in enumerate(N) if 'desc' in elem]
        print(E)
        Name = '_'.join(N[0:E[0]])
        
        for i in range(len(steps[0])):
            
            Suffix = 'desc-mean+' + str(int(steps[0][i]/60))
            img = opj(Dir_pet,'_'.join([Name,Suffix,'pet.nii.gz']))
            
            nl =  'get the mean amplitude signal per each Roi'
            run_cmd.msg(nl,diary_name)
            roi_signal.read(img,dseg,[0],[0],diary_name,sing_afni)
            
            nl =  'normalize the image to the anat and to the reference template and skull-strip it'
            run_cmd.msg(nl,diary_name)
            norm2template.pet(img,template,T1wacpcmask,T1wacpc,templateT1w,templatemask,maxval,dir_transfo,diary_name,sing_wb)
        
        
            for x,z in zip(['acpc',opj('templates', template)],['acpc',template]):
                img1 = opj(Dir_pet,x,'_'.join([Name, 'space-' + z, Suffix      ,'pet.nii.gz']))
                imgSS= opj(Dir_pet,x,'_'.join([Name, 'space-' + z, Suffix+'-SS','pet.nii.gz']))
                paramImg.bratio(img1,imgSS,diary_name,sing_wb)


def srtm(file,dseg,ROIs_vol,ROIs_Kin,template,T1wacpcmask,T1wacpc,templateT1w,templatemask,maxval,dir_transfo,
atlas_label,REF,REF_name,pet_raw_dir,isotope,FT,FD,DV_time,diary_name,sing_afni,sing_turku,sing_wb):
    
    nl =  'get the mean amplitude signal per each Roi'
    run_cmd.msg(nl,diary_name)
    roi_signal.read(file,dseg,[0],[0], diary_name,sing_afni)
    
    roi_signal.dft(file,dseg,ROIs_vol,ROIs_Kin,
                      atlas_label,
                      REF,
                      DV_time,
                      'study_01234',
                      REF_name,
                      isotope,
                      FT,
                      FD,
                      pet_raw_dir,
                      diary_name,sing_turku)
    

    IMGs=['BP','BP-err','STRM','DVR']
    E = -1
    N = opb(file).split('_')
    dir = opd(file)
    E = [i for i, elem in enumerate(N) if 'desc' in elem]
    Name = '_'.join(N[0:E[0]])    
    for i in range(len(IMGs)):
        img = opj(dir,Name + '_desc-' + IMGs[i] + '_pet.nii.gz')
        if opi(img):
            nl =  'normalize the image to the anat and to the reference template and skull-strip it'
            run_cmd.msg(nl,diary_name)
            norm2template.pet(img,template,T1wacpcmask,T1wacpc,templateT1w,templatemask,maxval,dir_transfo,diary_name,sing_wb)
    

def suv(ratio,file,dseg,template,T1wacpcmask,T1wacpc,templateT1w,templatemask,maxval,dir_transfo,diary_name,sing_afni,sing_wb):
    
    Dir_pet = opd(file)
    N = opb(file).split('_')
    E = [i for i, elem in enumerate(N) if 'desc' in elem]
    Name = '_'.join(N[0:E[0]])
    suffix = N[E[0]]
    if not ratio == '':
        paramImg.suv(file,ratio)
        file = opj(Dir_pet,'_'.join([Name,suffix + '-SUV','pet.nii.gz']))

    nl =  'get the mean amplitude signal per each Roi'
    run_cmd.msg(nl,diary_name)
    roi_signal.read(file,dseg,[0],[0],diary_name,sing_afni)
    nl =  'normalize the image to the anat and to the reference template and skull-strip it'
    run_cmd.msg(nl,diary_name)
    norm2template.pet(file,template,T1wacpcmask,T1wacpc,templateT1w,templatemask,maxval,dir_transfo,diary_name,sing_wb)
    
    for i,j in zip(['acpc',opj('templates', template)],['acpc',template]):
        img1= opj(Dir_pet,i,'_'.join([Name, 'space-' + j, suffix,'pet.nii.gz']))
        imgSS= opj(Dir_pet,i,'_'.join([Name, 'space-'+ j, suffix + '-SS','pet.nii.gz']))
        paramImg.bratio(img1,imgSS,diary_name,sing_wb)
