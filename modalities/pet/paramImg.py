#import
import os
import ants
import json

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
opi = os.path.isfile

from Tools import run_cmd, filename

Palette = (' MODE_USER_SCALE -thresholding THRESHOLD_TYPE_OFF THRESHOLD_TEST_SHOW_INSIDE 0 2.5' +
              ' -pos-user 0 2.5 -interpolate true -palette-name ROY-BIG -disp-pos true -disp-neg false -disp-zero false')

def bratio(img1,imgSS,diary_name,sing_wb):
    
    nl =  'normalize the voxel signal by the mean brain signal'
    run_cmd.msg(nl,diary_name,'OKGREEN')

    for fneeded in [img1,imgSS]:
        if not opi(fneeded):
            raise FileNotFoundError(f"Input file not found: {fneeded}")

    # get the mean value from the skull-stripped image
    SS    = ants.image_read(imgSS)
    meanS = SS[SS>0].mean()
    nl = 'mean value of ' + imgSS + ': ' + str(meanS)
    run_cmd.msg(nl, diary_name,'OKGREEN')

    for img in[imgSS,img1]:
        fname, suffix = filename.get_desc(opb(img), 'desc')
        Dir_pet = opd(img)

        IMG     = ants.image_read(img)
        header  = ants.image_header_info(img)

        newsuffix = ['desc','bratio']
        if 'SS' in suffix:
            newsuffix.append('SS')

        newfile = opj(Dir_pet, '_'.join([fname, '-'.join(newsuffix), 'pet.nii.gz']))

        if len(header['dimensions'])==4:
            addsuffix = ['mean']
            meanfile = opj(Dir_pet, '_'.join([fname, '-'.join(newsuffix + addsuffix), 'pet.nii.gz']))

            img_unmerged = ants.ndimage_to_list(IMG)
            nii = img_unmerged[0] * 0
            img_list = list()
            for i in range(len(img_unmerged)):
                temp = img_unmerged[i] / meanS
                nii+=temp
                img_list.append(temp)
                nii=nii/len(img_unmerged)
            new_img = ants.list_to_ndimage(IMG, img_list)

            ants.image_write(nii, meanfile)

            dictionary = {"Sources": newfile,
                          "Description": 'Mean image (3dTstat, AFNI).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(meanfile.replace('.nii.gz','.json'), "w") as outfile:
                outfile.write(json_object)

        else:
            new_img = ants.image_clone(IMG)
            new_img = new_img/meanS

        ants.image_write(new_img, newfile)
        cmd = (sing_wb + 'wb_command -volume-palette' + ' ' + newfile + Palette)
        run_cmd.run(cmd,diary_name)

        dictionary = {"Sources": img,
        "Description": 'normalisation of the signal by the global mean (ANTspy)', }
        json_object = json.dumps(dictionary, indent=2)
        with open(newfile.replace('.nii.gz','.json'), "w") as outfile:
            outfile.write(json_object)



def suv(file,ratio,diary_name):

    ''' SUV stands for Standardized Uptake Value. '''

    nl =  'normalize the voxel signal by the ratio'
    run_cmd.msg(nl,diary_name)
    fname, _ = filename.get_desc(opb(file), 'desc')
    Dir_pet = opd(file)
    newname = '_'.join([fname,'desc-SUV','pet'])
    
    IMG     = ants.image_read(file)
    header  = ants.image_header_info(file)

    if len(header['dimensions'])==4:
        img_unmerged = ants.ndimage_to_list(IMG)
        nii = img_unmerged[0] * 0
        img_list = list()
        for i in range(len(img_unmerged)):
            temp = img_unmerged[i] / ratio
            img_list.append(temp)
            nii += temp
        nii=nii/len(img_unmerged)
        ants.image_write(nii, opj(Dir_pet, '_'.join([fname,'desc-SUV-mean','pet.nii.gz'])))
        new_img = ants.list_to_ndimage(IMG, img_list)
    else:
        new_img = ants.image_clone(IMG)
        new_img = new_img / ratio
    ants.image_write(new_img, opj(Dir_pet, newname + '.nii.gz'))

    dictionary = {"Sources": [file,
                              str(ratio)],
                  "Description": 'Mean image weighted by subject weight and tracor Qte (3dTstat, AFNI and ANtspy).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(Dir_pet, newname + '.json'), "w") as outfile:
        outfile.write(json_object)
    

    



