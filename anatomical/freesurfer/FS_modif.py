#import
import os
import shutil
import sys
import math
import ants
import numpy as np

opj = os.path.join
ope = os.path.exists
opb = os.path.basename
opd = os.path.dirname
opi = os.path.isfile

from Tools import run_cmd

from anatomical.freesurfer import mgz2ants

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
#min_max = [[1,35,80,110,40,90],[10,50,80,110,120,190]]
min_max = [[1,35,80,110,82,90],[10,50,80,110,120,190]]

# for RS analysis
Ventri_mask = [4,43,5,44,14,15]
White_mask  = [2,41,251,252,253,254,255]


def prepa(Dir,midsection,new_size,bckup,diary_name,sing_fs):

    nl = '2.0 : Prepare for the left/rigth hemisphere separation (get the filled.nii.gz file)'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    # species_check = ['human', 'chimpanzee', 'macaque', 'baboon', 'vervet', 'marmoset', 'pig']
    # if species in species_check:
    #    # for building relevant surfaces in Freesurfer the olfactif bulb has been labeled as optic tract in some species
    #    WM_edit.append(85)
    if bckup ==1:
        for data in ['brain','filled','wm']:
            shutil.copy(opj(Dir, data + '.mgz'),opj(Dir, data + '_old.mgz'))

    BRAIN = mgz2ants.read(opj(Dir,'orig.mgz'),diary_name,sing_fs)
    aseg  = mgz2ants.read(opj(Dir,'aseg.mgz'),diary_name,sing_fs)

    Sig_max = BRAIN.max()
    norm = (BRAIN / Sig_max) * 110
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
    if not new_size == '':
        new_img = ants.smooth_image(new_img, new_size)

    ants.image_write(new_img, opj(Dir, 'brain.nii.gz'), ri=False)
    cmd = (sing_fs + 'mri_convert ' + opj(Dir, 'brain.nii.gz') + ' '  + opj(Dir, 'brain.mgz'))
    run_cmd.run(cmd, diary_name)

    img_aseg = ants.mask_image(aseg, aseg, aseg_edit)
    img_aseg[img_aseg > 0] = 250

    wm_aseg = ants.mask_image(aseg, aseg, wm_aseg_edit)
    wm_aseg[wm_aseg > 0] = 1

    WM_seg = ants.mask_image(WM2, wm_aseg, 1)
    WM_FS  = img_aseg + WM_seg
    ants.image_write(WM_FS,  opj(Dir, 'wm.nii.gz'), ri=False)
    cmd = (sing_fs + 'mri_convert ' + opj(Dir, 'wm.nii.gz') + ' ' + opj(Dir, 'wm.mgz'))
    run_cmd.run(cmd, diary_name)

    # filled image (define left and right hemispheres)
    filled = ants.image_clone(BRAIN) * 0
    hd_BRAIN = mgz2ants.readHD(opj(Dir, 'orig.mgz'), diary_name, sing_fs)
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

    ants.image_write(filled, opj(Dir, 'filled.nii.gz'), ri=False)
    cmd = (sing_fs + 'mri_convert ' + opj(Dir, 'filled.nii.gz') + ' ' + opj(Dir, 'filled.mgz'))
    run_cmd.run(cmd, diary_name)

    nl = 'Done but better check the result'
    run_cmd.msg(nl, diary_name,'OKGREEN')