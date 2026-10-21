#import
import os
import glob
from os.path import samefile

import ants
import json

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd,getpath
from modalities.pet import reorient

def save_json(jsonfile,FS,FT,FD):
  with open(jsonfile, "w", encoding="utf-8") as f:
    json.dump({
            "FrameTimesStart": FS,
            "FrameReferenceTime": FT,
            "FrameDuration": FD,
        }, f, indent=4)

def run(ID,tracor,data_path, reference,humanPosition,animalPosition, anatraw, Need_check_orient, IhaveCT,ref_CT, ct_pet_same_space,pethd,IhaveAnat,anat_pet_same_space,Skip_step,diary_file,sing_afni, sing_fs):

    (dir_prepro_raw, dir_prepro_raw_process, _, _, path_pet,
     dir_prepro_orig, _, _,dir_prepro_orig_process, _,
     _, _, _, _, _, _, _, _,_, _, _, _) = getpath.pet(data_path, reference)

    list_img = sorted(glob.glob(opj(dir_prepro_raw, ID + '_trc-' + tracor + '*_pet.nii.gz')))
    print(opj(dir_prepro_raw, ID + '_trc-' + tracor + '*_pet.nii.gz'))
    print(list_img)
    # mean image

    os.makedirs(path_pet,exist_ok=True)
    os.makedirs(dir_prepro_raw_process,exist_ok=True)
    os.makedirs(dir_prepro_orig, exist_ok=True)

    for img in list_img:
        imghd = ants.image_header_info(img)
        if imghd['nDimensions'] == 4:
            mean_img = opj(dir_prepro_raw_process,opb(img.replace('_pet.nii.gz','_desc-mean_pet.nii.gz')))
            if not opi(mean_img):
                cmd = (sing_afni + '3dTstat -mean -prefix ' + mean_img + ' ' + img)
                run_cmd.run(cmd, diary_file)
                dictionary = {"Sources": img,
                              "Description": 'mean image (3dTstat from AFNI).', }
                json_object = json.dumps(dictionary, indent=2)
                with open(mean_img.replace('.nii.gz', '.json'), "w") as outfile:
                    outfile.write(json_object)

            jsonfile = img.replace('.nii.gz', '.json')
            if opi(jsonfile):
                f = open(jsonfile)
                info_pet = json.load(f)
            if "FrameTimesStart" in info_pet and "FrameReferenceTime" in info_pet and "FrameDuration" in info_pet:
                timejson = opj(dir_prepro_orig, img.replace('_pet.nii.gz', '_desc-frametiming_pet.json'))
                if not opi(timejson):
                    save_json(timejson, info_pet["FrameTimesStart"], info_pet["FrameReferenceTime"],
                              info_pet["FrameDuration"])

    if 1 in Skip_step:
        run_cmd.msg('INFO: skip step ' + str(1), diary_file, 'OKGREEN')
    else:
        for img in list_img:
            jsonfile= img.replace('.nii.gz', '.json')
            if opi(jsonfile):
                f = open(jsonfile)
                info_pet = json.load(f)
                if humanPosition =='unknown':
                    if "PatientPosition" in info_pet:
                        humanPosition = info_pet["PatientPosition"]

                if not pethd :
                    # you know that the reconstruction used a wrong header format
                    if IhaveAnat and anat_pet_same_space :
                        REF_image = anatraw
                    else:
                        if IhaveCT and ct_pet_same_space:
                            REF_image = ref_CT
                        reorient.doit(dir_prepro_orig_process, img,humanPosition,animalPosition, REF_image, Need_check_orient, True,
                              diary_file, sing_afni, sing_fs)