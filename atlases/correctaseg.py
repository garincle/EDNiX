#import
import ants
import numpy as np

from Tools import run_cmd

csf        = 24
cortex     = [3,42]
white      = [2,78,41,79,14,15,5,4,44,43,10,11,12,13,17,18,26,28,49,50,51,52,53,54,58,60]
others     = [10,11,12,13,17,18,26,28,49,50,51,52,53,54,58,60,85]
not_smooth = [16,14,15,5,4,7,8,44,43,46,47]

def smooth(file1,file2,LR,midsection,diary_name):

    nl = 'gentle smooth of the aseg file afer manual correction to reduce the 3d noises'
    run_cmd.msg(nl, diary_name,'HEADER')

    img = ants.image_read(file1)
    imghd   = ants.image_header_info(file1)
    spacing = np.min(imghd['spacing'])

    run_cmd.msg(spacing, diary_name,'OKGREEN')

    new_img = ants.image_clone(img) * 0
    if LR == 'l':
        img[midsection:,:,:] = 0
    elif LR == 'r':
        img[0:midsection,:,:] = 0

    img0 = ants.get_mask(img, low_thresh=1)
    img0 = ants.smooth_image(img0,spacing)
    img0 = ants.threshold_image(img0, 0.5, 1, 1, 0, True)
    img0 = ants.morphology(img0, operation='dilate', radius=spacing*2, mtype='binary', shape='ball')
    img0=img0*csf
    new_img[img0==csf]=csf

    img1 = ants.image_clone(img) * 0
    img2 = ants.mask_image(img, img, white)
    img1[img2 > 0] = 1
    img1 = ants.iMath(img1, 'FillHoles', 2)
    img1 = ants.smooth_image(img1, spacing)
    img1 = ants.threshold_image(img1, 0.5, 1, 1, 0, True)

    if LR == 'l':
        new_img[img1==1]=2
    elif LR == 'r':
        new_img[img1==1]=41

    for i in others:
        img1 = ants.get_mask(img, low_thresh=i,high_thresh=i)
        img1 = ants.smooth_image(img1, spacing)
        img1 = ants.threshold_image(img1, 0.3, 1, 1, 0, True) * i
        new_img[img1==i] = i

    for i in not_smooth:
        new_img[img == i] = i

    for i in cortex:
        img1 = ants.get_mask(img, low_thresh=i,high_thresh=i)
        img1 = ants.smooth_image(img1, spacing)
        img1 = ants.threshold_image(img1, 0.2, 1, 1, 0, True) * i
        new_img[img1==i]=i

    if LR == 'l':
        new_img[midsection:,:,:] = 0
    elif LR == 'r':
        new_img[0:midsection,:,:] = 0

    ants.image_write(new_img,file2)













