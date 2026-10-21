#import
import os
import ants
import json
import glob

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd,extract_filename
import roi_signal,cut_image,norm2template,paramImg

def suv(file,dseg,template,T1wacpcmask,T1wacpc,templateT1w,templatemask,maxval,dir_transfo,atlas_label,
        do_cut,duration,ratio,bratio, diary_name,sing_afni,sing_wb,same_space):

    Dir_pet = opd(file)
    fname,suffix = extract_filename.get_desc(opb(file),'desc')
    hd           = ants.image_header_info(file)

    patterns = ['', '-mean', '-SS', '-SS-mean']

    if not any(x in opb(file).split('_') for x in ("space-acpc", "space-" + template)):
        maskfile = sorted(glob.glob(opj(Dir_pet, '*desc-coreg_mask.nii.gz')))[0]

    msk = ants.image_read(maskfile)

    if not ratio == 1:
        # Standardized Uptake Value
        root = 'SUV'
        newSuffix = [f'desc-{root}{p}'for p in patterns]
        print(newSuffix)

        output_1a = opj(Dir_pet, '_'.join([fname, newSuffix[0], 'pet.nii.gz']))
        output_1b = opj(Dir_pet, '_'.join([fname, newSuffix[1], 'pet.nii.gz']))
        output_1c = opj(Dir_pet, '_'.join([fname, newSuffix[2], 'pet.nii.gz']))
        output_1d = opj(Dir_pet, '_'.join([fname, newSuffix[3], 'pet.nii.gz']))

        paramImg.suv(file,ratio,diary_name)

        nl = 'get the mean amplitude signal per each Roi'
        run_cmd.msg(nl, diary_name,'OKGREEN')
        if len(hd['dimensions']) == 4:
            # take the result mean image instead
            roi_signal.read(output_1b, dseg, atlas_label, 0, '', diary_name, sing_afni)
        else:
            # take the result image
            roi_signal.read(output_1a, dseg, atlas_label, 0, '', diary_name, sing_afni)

        # skullstripping
        IMG = ants.image_read(output_1a)

        if len(hd['dimensions']) == 4:
            img_unmerged = ants.ndimage_to_list(IMG)
            nii = img_unmerged[0] * 0
            img_list = list()
            for i in range(len(img_unmerged)):
                temp = img_unmerged[i] * msk
                nii += temp
                img_list.append(temp)
                nii = nii / len(img_unmerged)
            new_img = ants.list_to_ndimage(IMG, img_list)
            ants.image_write(nii, output_1d)
            dictionary = {"Sources": output_1c,
                          "Description": 'Mean image (3dTstat, AFNI).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_1d.replace('.nii.gz','.json'), "w") as outfile:
                outfile.write(json_object)
        else:
            new_img=IMG*msk

        ants.image_write(new_img, output_1c)
        dictionary = {"Sources": [output_1a,
                                  maskfile],
                      "Description": 'Skull stripping (with ants).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_1c.replace('.nii.gz','.json'), "w") as outfile:
            outfile.write(json_object)

        maxval=2.5

    else:
        output_1a = file
        newSuffix=[]
        [newSuffix.append(suffix+i) for i in patterns]

    if bratio == 1:
        root = 'bratio'
        newnewSuffix = [f'desc-{root}{p}' for p in patterns]
        print(newnewSuffix)

        output_2a = opj(Dir_pet, '_'.join([fname, newnewSuffix[0], 'pet.nii.gz']))
        output_2b = opj(Dir_pet, '_'.join([fname, newnewSuffix[1], 'pet.nii.gz']))

        imgSS     = opj(Dir_pet, '_'.join([fname, newSuffix[2], 'pet.nii.gz']))
        meanimgSS = opj(Dir_pet, '_'.join([fname, newSuffix[3], 'pet.nii.gz']))
        print(imgSS)
        if not opi(imgSS):
            # skullstripping
            IMG = ants.image_read(output_1a)

            if len(hd['dimensions']) == 4:
                img_unmerged = ants.ndimage_to_list(IMG)
                nii = img_unmerged[0] * 0
                img_list = list()
                for i in range(len(img_unmerged)):
                    temp = img_unmerged[i] * msk
                    nii += temp
                    img_list.append(temp)
                    nii = nii / len(img_unmerged)
                new_img = ants.list_to_ndimage(IMG, img_list)
                ants.image_write(nii, meanimgSS)
                dictionary = {"Sources": imgSS,
                              "Description": 'Mean image (3dTstat, AFNI).', }
                json_object = json.dumps(dictionary, indent=2)
                with open(meanimgSS.replace('.nii.gz', '.json'), "w") as outfile:
                    outfile.write(json_object)
            else:
                new_img = IMG * msk

            ants.image_write(new_img, imgSS)
            dictionary = {"Sources": [output_1a,
                                      maskfile],
                          "Description": 'Skull stripping (with ants).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(imgSS.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)

        nl = 'get the mean amplitude signal per each Roi'
        run_cmd.msg(nl, diary_name)
        if len(hd['dimensions']) == 4:
            roi_signal.read(meanimgSS, dseg, atlas_label, 0, '', diary_name, sing_afni)
        else:
            roi_signal.read(imgSS, dseg, atlas_label, 0, '', diary_name, sing_afni)

        paramImg.bratio(output_1a, imgSS, diary_name, sing_wb)

        nl = 'get the mean amplitude signal per each Roi'
        run_cmd.msg(nl, diary_name)
        if len(hd['dimensions']) == 4:
            roi_signal.read(output_2b, dseg, atlas_label, 0, '', diary_name, sing_afni)
        else:
            roi_signal.read(output_2a, dseg, atlas_label, 0, '', diary_name, sing_afni)

        maxval=1.5
    else:
        output_2a = output_1a

    nl = 'normalize the image to the anatomical and to the reference template and skull-strip it'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    unique = list(dict.fromkeys([file, output_1a,output_2a]))
    print(unique)
    for files in unique:
        norm2template.pet(files, template, T1wacpcmask, T1wacpc, templateT1w, templatemask, maxval, dir_transfo, diary_name,sing_wb,same_space)
        if do_cut==1:
            nl =  'mean by each time frames'
            run_cmd.msg(nl,diary_name)
            steps, frame_stop = cut_image.setTime(file, duration)
            cutlist = cut_image.doit(files,maskfile,steps[0],frame_stop[0],duration,ratio,diary_name,sing_afni)

            for i in range(len(steps[0])):

                nl =  'get the mean amplitude signal per each Roi'
                run_cmd.msg(nl,diary_name,'OKGREEN')
                roi_signal.read(cutlist[i],dseg,atlas_label,0,'',diary_name,sing_afni)

                nl =  'normalize the image to the anatomical and to the reference template and skull-strip it'
                run_cmd.msg(nl,diary_name)
                norm2template.pet(cutlist[i],template,T1wacpcmask,T1wacpc,templateT1w,templatemask,maxval,dir_transfo,diary_name,sing_wb,same_space)


def srtm(file,dseg,ROIs_vol,ROIs_Kin,template,T1wacpcmask,whiteacpcmask,csfacpcmask,grayacpcmask,T1wacpc,templateT1w,templatemask,maxval,dir_transfo,
atlas_label,REF,isotope,DV_time,diary_name,sing_afni,sing_turku,sing_wb,same_space):

    for fneeded in [file,dseg,T1wacpcmask,whiteacpcmask,csfacpcmask,T1wacpc,templateT1w,templatemask]:
        if not opi(fneeded):
            raise FileNotFoundError(f"Input file not found: {fneeded}")

    nl =  'get the mean amplitude signal per each Roi'
    run_cmd.msg(nl,diary_name,'OKGREEN')

    IMGs = ['BP', 'BP-err', 'STRM', 'DVR']

    roi_signal.read(file,dseg,atlas_label,0,'', diary_name,sing_afni)

    mskfile = ''
    if REF['type'] == 'mask':
        selectmsk = {'white':whiteacpcmask,
                     'csf':csfacpcmask,
                     'white+csf':opj(opd(whiteacpcmask)),
                     'gray':grayacpcmask,
                     'brain':T1wacpcmask}
        if REF['name']=='white+csf':
            cmd = (sing_afni + '3dcalc -a ' + whiteacpcmask + ' -b '+ csfacpcmask +
                   ' -prefix ' + opj(opd(whiteacpcmask),'white+csf_mask.nii.gz')
                   + '-expr "a*b"')
            run_cmd(cmd,diary_name)
        mskfile = selectmsk[REF['name']]

    roi_signal.dft(file,dseg,ROIs_vol,ROIs_Kin,atlas_label,mskfile,REF,DV_time,'study_01234',isotope,
                   diary_name,sing_turku)

    dir = opd(file)
    fname,_ = extract_filename.get_desc(opb(file),'desc')
    for i in range(len(IMGs)):
        img = opj(dir,'_'.join([fname,'desc-' + IMGs[i],'pet.nii.gz']))
        if opi(img):
            nl =  'normalize the image to the anatomical and to the reference template and skull-strip it'
            run_cmd.msg(nl,diary_name,'OKGREEN')
            norm2template.pet(img,template,T1wacpcmask,T1wacpc,templateT1w,templatemask,maxval,dir_transfo,diary_name,sing_wb,same_space)

