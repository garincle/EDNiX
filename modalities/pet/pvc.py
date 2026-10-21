import os
import ants
import numpy as np
from nipype.interfaces.petpvc import PETPVC

opj = os.path.join
opd = os.path.dirname
opi = os.path.isfile

def correct(img,seg,algo='RBV'):

    '''
    algo = MG, GTM, LAB,MTC, RBV or IY
    see Thomas et al 2016 Phys. Med. Biol. 61 7975
    '''

    for fneeded in [img,seg]:
        if not opi(fneeded):
            raise FileNotFoundError(f"Input file not found: {fneeded}")

    hd1 = ants.image_header_info(img)
    vox = list(hd1['spacing'])
    print(vox)

    hd = ants.image_header_info(seg)
    print(hd['nDimensions'])
    if hd['nDimensions']==3:
        print('need to do something')
        new_seg = seg.replace('.nii.gz','-4D.nii.gz')
        SEG     = ants.image_read(seg)
        spacing   = SEG.spacing + (1,)
        orig      = SEG.origin + (0,)
        direction = SEG.direction
        new_direction = np.eye(4)
        new_direction[:3, :3] = direction

        GM  = ants.get_mask(SEG,low_thresh=2,high_thresh=2,cleanup=0)
        WM  = ants.get_mask(SEG, low_thresh=3, high_thresh=3,cleanup=0)
        CSF = ants.get_mask(SEG, low_thresh=1, high_thresh=1,cleanup=0)
        imageTar = ants.make_image((*GM.shape, 3))
        imageTar.set_spacing(spacing)
        imageTar.set_origin(orig)
        imageTar.set_direction(new_direction)
        res = ants.list_to_ndimage(imageTar, [GM, WM, CSF])
        ants.image_write(res,new_seg)
    else:
        new_seg = seg

    pvcf = PETPVC()
    pvcf.inputs.in_file = img
    pvcf.inputs.mask_file = new_seg
    pvcf.inputs.out_file = img.replace('.nii.gz','_desc-PVC.nii.gz')
    pvcf.inputs.pvc = algo
    pvcf.inputs.fwhm_x = vox[0]
    pvcf.inputs.fwhm_y = vox[1]
    pvcf.inputs.fwhm_z = vox[2]
    pvcf.run()
