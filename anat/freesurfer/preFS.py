#import
import os
import json
import math
import ants
import numpy as np

opj = os.path.join
ope = os.path.exists
opb = os.path.basename
opd = os.path.dirname
opi = os.path.isfile

from Tools import run_cmd, get_orientation
from Tools import getpath

# for wm.mgz
aseg_edit    = [4,43,10,49,11,50,12,51,26,58,28,60,31,63,77,85]
wm_aseg_edit = [2,41,251,252,253,254,255,13,52,16]


filled_edit  = [[4,2,10,11,12,13,26,28,31],
                [43,41,49,50,51,52,58,60,63],
                [77,85]]
CSF_edit = [4,5,43,44,14,15,24,72,30,62,31,63]
WM_edit  = [2,41,85,7,46,77,251,252,253,254,255]
GM_edit  = [3,42,8,47,10,49,11,12,13,16,17,18,26,28,50,51,52,53,54,58,60]

# CSF, WM, GM min and max for T1 and T2
min_max = [[1,35,80,110,40,90],
           [10,50,80,110,120,190]]

# for RS analysis
Ventri_mask = [4,43,5,44,14,15]
White_mask  = [2,41,251,252,253,254,255]

def prepa_img(Subname,Ref_file,dir_prepro,aseg_file,labels_dir,midsection,diary_name):

    nl = ' step 10 : Prepare for the left/rigth hemisphere separation (get the filled.nii.gz file)'
    run_cmd.msg(nl, diary_name,'HEADER')

    # species_check = ['human', 'chimpanzee', 'macaque', 'baboon', 'vervet', 'marmoset', 'pig']
    # if species in species_check:
    #    # for building relevant surfaces in Freesurfer the olfactif bulb has been labeled as optic tract in some species
    #    WM_edit.append(85)

    BRAIN = ants.image_read(Ref_file)
    aseg = ants.image_read(aseg_file)
    hd_BRAIN = ants.image_header_info(Ref_file)

    # 1) Strip off extensions (.nii.gz or .nii)
    fname = opb(Ref_file)
    if fname.lower().endswith('.nii.gz'):
        stem = fname[:-7]
    else:
        stem = os.path.splitext(fname)[0]

    # 2) Get the last “chunk” after the final underscore
    #    e.g. “…_desc-SS_T1”    → “T1”
    #         “…_desc-SS_T1w”   → “T1w”
    #         “…_desc-SS_customT2extra” → “customT2extra”
    suffix = stem.rsplit('_', 1)[-1]

    # 3) Determine Type_norm by looking for the most specific tokens first
    if 'T1w' in suffix:
        Type_norm = 'T1w'
    elif 'T2w' in suffix:
        Type_norm = 'T2w'
    elif 'T1' in suffix:
        Type_norm = 'T1w'
        run_cmd.msg(
            f'WARNING: Found "{suffix}", mapping to T1w for surface purposes.',
            diary_name, 'WARNING')
    elif 'T2' in suffix:
        Type_norm = 'T2w'
        run_cmd.msg(
            f'WARNING: Found "{suffix}", mapping to T2w for surface purposes.',
            diary_name, 'WARNING'
        )
    else:
        raise ValueError(
            f"Cannot infer T1/T2 suffix from filename '{Ref_file}'. "
            "Expected a chunk containing 'T1' or 'T2'."
        )

    # Now Type_norm is one of 'T1w' or 'T2w'
    print("Using Type_norm =", Type_norm)

    Sig_max = BRAIN.max()
    norm = (BRAIN / Sig_max) * 110

    if Type_norm == 'T2w':
        new_min_max = min_max[1]
    else:
        # consider that it is a T1w image
        new_min_max = min_max[0]

    new_img = ants.image_clone(BRAIN)*0

    for select,min_im,max_im in zip([CSF_edit,GM_edit],
                                    [new_min_max[0],new_min_max[4]],
                                    [new_min_max[1],new_min_max[5]]):
        img = ants.mask_image(norm,aseg,select)
        mask_img = ants.image_clone(img)
        mask_img[mask_img>0]=1
        img1 = (img - img.min()) / int(img.max() - img.min())
        img2 = img1 * int(max_im - min_im) + int(min_im)
        img3 = ants.mask_image(img2, mask_img, 1)
        new_img = new_img + img3

    WM = ants.mask_image(norm, aseg, WM_edit)
    mask_img = ants.image_clone(WM)
    mask_img[mask_img > 0] = 1
    WM1 = (WM - WM.min()) / int(WM.max() - WM.min())
    WM2 = WM1 * int(new_min_max[3] - new_min_max[2]) + int(new_min_max[2])
    WM3 = ants.mask_image(WM2, mask_img, 1)
    new_img = new_img + WM3

    ants.image_write(new_img, opj(dir_prepro, Subname + '_desc-norm_' + Type_norm + '.nii.gz'), ri=False)

    img_aseg = ants.mask_image(aseg, aseg, aseg_edit)
    img_aseg[img_aseg > 0] = 250

    wm_aseg = ants.mask_image(aseg, aseg, wm_aseg_edit)
    wm_aseg[wm_aseg > 0] = 1

    WM_seg = ants.mask_image(WM2, wm_aseg, 1)
    WM_FS  = img_aseg + WM_seg
    ants.image_write(WM_FS, opj(labels_dir, Subname + '_seg-wm_dseg.nii.gz'), ri=False)

    # filled image (define left and right hemispheres)
    filled = ants.image_clone(BRAIN) * 0
    if midsection =='':
        midsection = math.ceil(int(hd_BRAIN['dimensions'][0])/2)

    thalamus = ants.mask_image(aseg, aseg, 10)
    cent = ants.get_centroids(thalamus)
    pos  = cent[0][0] * hd_BRAIN['direction'][0][0]
    orig = hd_BRAIN['origin'][0] * hd_BRAIN['direction'][0][0]
    dist = np.abs(pos - orig)/ hd_BRAIN['spacing'][0]

    for i,value in zip([0,1],[127,255]):
        img_filled = ants.mask_image(aseg, aseg, filled_edit[i])
        filled[img_filled > 0] = value

    # deal with bilateral regions
    img_filled_left = ants.mask_image(aseg, aseg, filled_edit[2])
    img_filled_right = ants.image_clone(img_filled_left)
    img_filled_left[midsection:,:,:]   = 0
    img_filled_right[0:midsection,:,:] = 0

    if dist < midsection:
        filled[img_filled_left > 0] = 127
        filled[img_filled_right > 0] = 255
    else:
        filled[img_filled_left > 0] = 255
        filled[img_filled_right > 0] = 127

    ants.image_write(filled, opj(labels_dir, Subname + '_seg-filled_dseg.nii.gz'), ri=False)
    nl = 'Done'
    run_cmd.msg(nl, diary_name,'OKGREEN')


def prepa_T2img(Subname,Ref_file,dir_prepro,asegfile,diary_name):

    nl = 'Prepare for the T2w image for the pial surface'
    run_cmd.msg(nl, diary_name,'HEADER')

    BRAIN   = ants.image_read(Ref_file)
    aseg = ants.image_read(asegfile)
    new_img = ants.image_clone(BRAIN) * 0

    Sig_max = BRAIN.max()
    norm    = (BRAIN / Sig_max) * 110
    new_min_max = min_max[1]

    for select, min_im, max_im in zip([CSF_edit, GM_edit],
                                      [new_min_max[0], new_min_max[4]],
                                      [new_min_max[1], new_min_max[5]]):
        img = ants.mask_image(norm, aseg, select)
        mask_img = ants.image_clone(img)
        mask_img[mask_img > 0] = 1
        img1 = (img - img.min()) / int(img.max() - img.min())
        img2 = img1 * int(max_im - min_im) + int(min_im)
        img3 = ants.mask_image(img2, mask_img, 1)
        new_img = new_img + img3

    WM = ants.mask_image(norm, aseg, WM_edit)
    mask_img = ants.image_clone(WM)
    mask_img[mask_img > 0] = 1
    WM1 = (WM - WM.min()) / int(WM.max() - WM.min())
    WM2 = WM1 * int(new_min_max[3] - new_min_max[2]) + int(new_min_max[2])
    WM3 = ants.mask_image(WM2, mask_img, 1)
    new_img = new_img + WM3

    ants.image_write(new_img, opj(dir_prepro, Subname + '_desc-norm_T2w.nii.gz'), ri=False)

    nl ='Done'
    run_cmd.msg(nl, diary_name,'OKGREEN')


def toFS(list1, Subname,change_hd,scaling,data_path,reference,BALSAname,diary_name,sing_fs):

    nl = 'Create the folders and the files for Freesurfer'
    run_cmd.msg(nl, diary_name,'HEADER')

    path_anat, dir_transfo, FS_dir, prepro_dir, dir_native, volumes_dir, labels_dir, masks_dir = getpath.anat(data_path,reference,BALSAname,False,False,'native')

    if ope(opj(FS_dir, Subname)) == False:
        os.makedirs(opj(FS_dir, Subname))
        os.makedirs(opj(FS_dir, Subname, 'mri'))
        os.makedirs(opj(FS_dir, Subname, 'surf'))
        os.makedirs(opj(FS_dir, Subname, 'stats'))
        os.makedirs(opj(FS_dir, Subname, 'label'))
        os.makedirs(opj(FS_dir, Subname, 'scripts'))

    if change_hd == 0:
        Ref = ''
    elif change_hd == 1:
        Ref = '_resamp-' + str(scaling)
        list1[0] = opj(prepro_dir, opb(list1[0].replace('.nii.gz', Ref + '.nii.gz')))
        labels_dir  = prepro_dir


    list1[1] = opj(prepro_dir,list1[1] + Ref + '.nii.gz')
    for i in range(2,5):
        list1[i] = opj(labels_dir,list1[i] + Ref + '.nii.gz')
    
    [_, _, fwdFS_cmd, _,_] = get_orientation.use_ants(list1[0], sing_fs)

    list2 = ['orig.mgz',
             'brain.mgz',
             'wm.mgz',
             'aseg.mgz',
             'filled.mgz']

    for i in range(len(list1)):
        if opi(list1[i])==True:
            cmd = sing_fs + 'mri_convert ' + fwdFS_cmd + list1[i] + ' ' + opj(FS_dir, Subname, 'mri', list2[i])
            run_cmd.run(cmd,diary_name)

    nl = 'Done'
    run_cmd.msg(nl, diary_name,'OKGREEN')


def msk_RS(Subname,brain_msk,aseg_img,diary_name,mskversion):

    nl = 'MASKS for denoising Resting State Data'
    run_cmd.msg(nl, diary_name,'HEADER')

    dir_msk = opd(brain_msk)
    msk = ants.image_read(brain_msk)

    dilate = ants.morphology(msk, operation='dilate', radius=2, mtype='binary', shape='ball')
    ants.image_write(dilate, opj(dir_msk, Subname + '_desc-dilat_mask.nii.gz'), ri=False)

    dictionary = {"Sources": opj(dir_msk, Subname + '_mask.nii.gz'),
                  "Description": 'Dilated version of the mask.', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_msk, Subname + '_desc-dilat_mask.json'), "w") as outfile:
        outfile.write(json_object)

    if mskversion == 'aseg':
        aseg = ants.image_read(aseg_img)
        img = ants.mask_image(aseg, aseg, GM_edit)
        img[img > 0] = 1
        dilate = ants.morphology(img, operation='dilate', radius=1, mtype='binary', shape='ball')

        ants.image_write(dilate, opj(dir_msk, Subname + '_desc-Gray_mask.nii.gz'), ri=False)

        dictionary = {"Sources": aseg_img,
                      "Description": 'Gray mask.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_msk, Subname + '_desc-Gray_mask.json'), "w") as outfile:
            outfile.write(json_object)


    for select,name,desc in zip([White_mask,Ventri_mask],['White','Vent'],['white matter','ventricles']):

        if mskversion == 'aseg':
            img = ants.mask_image(aseg, aseg, select)
            img[img > 0] = 1
            ants.image_write(img, opj(dir_msk, Subname + '_desc-' + name + '_mask.nii.gz'), ri=False)
            dictionary = {"Sources": aseg_img,
                          "Description": desc + ' mask.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_msk, Subname + '_desc-' + name + '_mask.json'), "w") as outfile:
                outfile.write(json_object)
            eroded = ants.morphology(img, operation='erode', radius=1, mtype='binary', shape='ball')
            ants.image_write(eroded, opj(dir_msk, Subname + '_desc-erod-' + name + '_mask.nii.gz'), ri=False)

            dictionary = {"Sources": opj(dir_msk, Subname + '_desc-' + name + '_mask.nii.gz'),
                          "Description": 'Eroded version of the ' + desc + ' mask.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_msk, Subname + '_desc-erod-' + name + '_mask.json'), "w") as outfile:
                outfile.write(json_object)

        elif mskversion == 'custom':
            if opi(opj(dir_msk, Subname + '_desc-' + name + '_mask.nii.gz')):
                img = ants.image_read(opj(dir_msk, Subname + '_desc-' + name + '_mask.nii.gz'))
                eroded = ants.morphology(img, operation='erode', radius=1, mtype='binary', shape='ball')
                ants.image_write(eroded, opj(dir_msk, Subname + '_desc-erod-' + name + '_mask.nii.gz'), ri=False)

                dictionary = {"Sources": opj(dir_msk, Subname + '_desc-' + name + '_mask.nii.gz'),
                              "Description": 'Eroded version of the ' + desc + ' mask.', }
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(dir_msk, Subname + '_desc-erod-' + name + '_mask.json'), "w") as outfile:
                    outfile.write(json_object)
            else:
                nl = 'We have not found the desc mask, it might cause problem if you want to extract/regress the desc signal'
                run_cmd.msg(nl, diary_name, 'WARNING')
        else:
            nl = 'fMRImasks must be custom or aseg'
            raise Exception(run_cmd.error(nl, diary_name))
    nl = 'Done'
    run_cmd.msg(nl, diary_name,'OKGREEN')



def resamp(list,data_path,reference,BALSAname,scaling,new_size,resamp, diary_name):

    nl = '## Change the resolution of the MRI images to ' + str(scaling)
    run_cmd.msg(nl, diary_name,'HEADER')

    _, _, _, dir_prepro, _, volumes_dir, labels_dir, _ = getpath.anat(data_path, reference,BALSAname, '', '', 'native')

    # set images to resolution 1mm
    name = opb(list[0]).split('.')[0].split('_')[0]
    #print(name)
    brain_img = ants.image_read(list[0])

    X = brain_img.origin
    S = brain_img.spacing

    new_orig = np.zeros(3)
    new_orig[0] = X[0] * scaling[0] / S[0]
    new_orig[1] = X[1] * scaling[1] / S[1]
    new_orig[2] = X[2] * scaling[2] / S[2]

    for new_img,dir_file,interp in zip(list,
                                        ['',
                                        dir_prepro,
                                        labels_dir,
                                        labels_dir,
                                        labels_dir],
                                        [3,3,1,1,1]):
        print(new_img)
        if new_img == list[0]:
            resamp_img = ants.image_read(new_img)
        else:
            resamp_img = ants.image_read(opj(dir_file,new_img + '.nii.gz'))
        ants.set_spacing(resamp_img, scaling)
        ants.set_origin(resamp_img, [new_orig[0], new_orig[1], new_orig[2]])
        if resamp == 1:
            resamp_img = ants.resample_image(resamp_img, (new_size, new_size, new_size), False, interp)
            resamp_img[resamp_img < 0] = 0
            if 'desc-norm' in new_img:
                resamp_img = ants.smooth_image(resamp_img, new_size)
            descript = 'reformat image with a voxel size of ' + str(new_size) + 'and a resampling of ' + str(scaling)
        else:
            descript = 'reformat image with a resampling of ' + str(scaling)
        if new_img == list[0]:
            newname = opj(dir_prepro, opb(new_img).replace('.nii.gz','_resamp-' + str(scaling[0]) + '.nii.gz'))
            print(newname)
        else:
            newname = opj(dir_prepro, new_img + '_resamp-' + str(scaling[0]) + '.nii.gz')
            print(newname)
        ants.image_write(resamp_img, newname)
        dictionary = {"Sources": [opj(dir_file,new_img),
                                  list[0]],
                      "Description": descript + ' (Antspy)' }
        json_object = json.dumps(dictionary, indent=2)
        with open(newname.replace('.nii.gz','.json'), "w") as outfile:
            outfile.write(json_object)

    nl = 'Done'
    run_cmd.msg(nl, diary_name,'OKGREEN')


def toannot(file,ref, change_hd,scaling,new_size,FS_dir, diary_name,sing_fs):

    nl = 'Convert the volume for FS to get a annot surface file'
    run_cmd.msg(nl, diary_name,'HEADER')

    [_, _, fwdFS_cmd, _, _] = get_orientation.use_ants(ref, sing_fs)

    name   = opb(file).split('.')[0]
    animal = name.split('_')[0]
    label  = name.split('_')[-2].split('-')[-1]

    if change_hd == 1:

        nl = '## Change the resolution of the MRI images to ' + str(scaling)
        run_cmd.msg(nl, diary_name,'OKGREEN')

        brain_img = ants.image_read(ref)
        dir_file = opd(ref)

        new_orig = brain_img.origin
        file_name = opj(dir_file, name + '_resamp-' + str(scaling[0]) + '.nii.gz')
        resamp_img = ants.image_read(file)
        ants.set_spacing(resamp_img, scaling)
        ants.set_origin(resamp_img, [new_orig[0], new_orig[1], new_orig[2]])
        resamp_img = ants.resample_image(resamp_img, (new_size, new_size, new_size), False, 1)
        resamp_img[resamp_img < 0] = 0

        ants.image_write(resamp_img, file_name)

        dictionary = {"Sources": [file,
                                  ref],
                      "Description": 'reformat image with a voxel size of ' + str(
                          new_size) + 'and a resampling of ' + str(scaling) + ' (Antspy)'}
        json_object = json.dumps(dictionary, indent=2)
        with open(file_name.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)

    else :
        file_name = file
    
    cmd = sing_fs + 'mri_convert ' + fwdFS_cmd + file_name + ' ' + opj(FS_dir, animal, 'mri', animal + '_' + label + '.mgz')
    
    run_cmd.run(cmd,diary_name)







