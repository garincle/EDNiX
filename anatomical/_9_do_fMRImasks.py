import os
import subprocess
import ants
import json
import datetime

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def do_fMRImasks(masks_dir, labels_dir, type_norm, fMRImasks, overwrite,s_bind,afni_sif,diary_file):

    ct = datetime.datetime.now()
    nl = 'Run anatomical._9_do_fMRImasks.do_fMRImasks'
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{nl}')


    # for RS analysis
    Ventri_mask = ',4,5,43,44,14'
    White_mask  = ',2,41'
    Gmask_mask = ',3,42'
    # 2) MASKS for denoising Resting State Data
    img = ants.image_read(opj(masks_dir,'brain_mask_in_anat_DC.nii.gz'))
    dilate = ants.morphology(img, operation='dilate', radius=2, mtype='binary', shape='ball')
    ants.image_write(dilate, opj(masks_dir, type_norm + 'brainmask_dilat.nii.gz'), ri=False)
    dictionary = {"Sources": opj(masks_dir,'brain_mask_in_anat_DC.nii.gz'),
                  "Description": 'Dilationg .', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(masks_dir, type_norm + 'brainmask_dilat.json'), "w") as outfile:
        outfile.write(json_object)

    if fMRImasks=='aseg':

        ##### Creat Vmask
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "amongst(a' + Ventri_mask + ')" -prefix ' + opj(masks_dir, type_norm + 'Vmask.nii.gz')
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": opj(labels_dir, type_norm + 'aseg.nii.gz'),
                      "Description": 'ventricular mask.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(masks_dir, type_norm + 'Vmask.json'), "w") as outfile:
            outfile.write(json_object)
        img = ants.image_read(opj(masks_dir, type_norm + 'Vmask.nii.gz'))
        eroded = ants.morphology(img, operation='erode', radius=1, mtype='binary', shape='ball')
        ants.image_write(eroded, opj(masks_dir, type_norm + 'Vmask_erod.nii.gz'), ri=False)
        dictionary = {"Sources": opj(masks_dir, type_norm + 'Vmask.nii.gz'),
                      "Description": 'eroded.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(masks_dir, type_norm + 'Vmask_erod.json'), "w") as outfile:
            outfile.write(json_object)

        ##### Creat Gmask
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "amongst(a' + Gmask_mask + ')" -prefix ' + opj(masks_dir, type_norm + 'Gmask.nii.gz')
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": opj(labels_dir, type_norm + 'aseg.nii.gz'),
                      "Description": 'Gray matterr mask.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(masks_dir, type_norm + 'Gmask.json'), "w") as outfile:
            outfile.write(json_object)

        img = ants.image_read(opj(masks_dir, type_norm + 'Gmask.nii.gz'))
        eroded = ants.morphology(img, operation='erode', radius=1, mtype='binary', shape='ball')
        ants.image_write(eroded, opj(masks_dir, type_norm + 'Gmask_erod.nii.gz'), ri=False)

        dictionary = {"Sources": opj(masks_dir, type_norm + 'Gmask.nii.gz'),
                      "Description": 'eroded.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(masks_dir, type_norm + 'Gmask_erod.json'), "w") as outfile:
            outfile.write(json_object)

        ##### Creat Wmask
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "amongst(a' + White_mask + ')" -prefix ' + opj(masks_dir, type_norm + 'Wmask.nii.gz')
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": opj(labels_dir, type_norm + 'aseg.nii.gz'),
                      "Description": 'white matterr mask.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(masks_dir, type_norm + 'Wmask.json'), "w") as outfile:
            outfile.write(json_object)

        img = ants.image_read(opj(masks_dir, type_norm + 'Wmask.nii.gz'))
        eroded = ants.morphology(img, operation='erode', radius=1, mtype='binary', shape='ball')
        ants.image_write(eroded, opj(masks_dir, type_norm + 'Wmask_erod.nii.gz'), ri=False)

        dictionary = {"Sources": opj(masks_dir, type_norm + 'Wmask.nii.gz'),
                      "Description": 'eroded.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(masks_dir, type_norm + 'Wmask_erod.json'), "w") as outfile:
            outfile.write(json_object)

    elif fMRImasks == 'custom':
        if ope(opj(masks_dir, type_norm + 'Vmask.nii.gz')):
            img = ants.image_read(opj(masks_dir, type_norm + 'Vmask.nii.gz'))
            eroded = ants.morphology(img, operation='erode', radius=1, mtype='binary', shape='ball')
            ants.image_write(eroded, opj(masks_dir, type_norm + 'Vmask_erod.nii.gz'), ri=False)
            dictionary = {"Sources": opj(masks_dir, type_norm + 'Vmask.nii.gz'),
                          "Description": 'eroded.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(masks_dir, type_norm + 'Vmask_erod.json'), "w") as outfile:
                outfile.write(json_object)


        else:
            nl = 'We have not found the Ventricular mask, it might cause problem if you want to extract/regress the ventricular signal'
            print(bcolors.WARNING + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')
        if ope(opj(masks_dir, type_norm + 'Wmask.nii.gz')):
            img = ants.image_read(opj(masks_dir, type_norm + 'Wmask.nii.gz'))
            eroded = ants.morphology(img, operation='erode', radius=1, mtype='binary', shape='ball')
            ants.image_write(eroded, opj(masks_dir, type_norm + 'Wmask_erod.nii.gz'), ri=False)
            dictionary = {"Sources": opj(masks_dir, type_norm + 'Wmask.nii.gz'),
                          "Description": 'eroded.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(masks_dir, type_norm + 'Wmask_erod.json'), "w") as outfile:
                outfile.write(json_object)
        else:
            nl = 'We have not found the white matter mask, it might cause problem if you want to extract/regress the white matter signal'
            print(bcolors.WARNING + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')
        if ope(opj(masks_dir, type_norm + 'Gmask.nii.gz')):
            nl = 'We have found a Gray mask'
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')
        else:
            nl = 'We have not found the Gray mask, it might cause problem if you want to extract the Gray signal specifically'
            print(bcolors.WARNING + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')
    else:
        nl = 'fMRImasks must be custom or aseg'
        diary.write(f'\n{nl}')
        raise Exception(bcolors.FAIL + nl + bcolors.ENDC)

    diary.write(f'\n')
    diary.close()
