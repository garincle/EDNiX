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
WBGM_edit  = [3,42,8,47,10,49,11,12,13,16,17,18,26,28,50,51,52,53,54,58,60]
Cereb_GM_edit  = [3,42,10,49,11,12,13,17,18,26,28,50,51,52,53,54,58,60]

# CSF, WM, GM min and max for T1 and T2
min_max = [[1,35,80,110,40,90],
           [10,50,80,110,120,190]]

# for RS analysis
Ventri_mask = [4,43,5,44,14,15]
White_mask  = [2,41,251,252,253,254,255]

# Surface reconstruction method (step 10 brain/wm/filled images, step 12 pial search distance).
#   'V2h' (default, chosen Nov 2026 from the surface tests of PR #84 on 5 library templates + 10 BIDS subjects):
#         synthetic brain.mgz with non-overlapping classes (CSF 10-35, GM 55-85 keeping the within-class texture,
#         WM 110, Gaussian sigma 0.5 voxel), wm.mgz WM 110 / subcortex 250, hippocampus/amygdala kept in the GM class
#         and out of wm/filled; pial search capped at MAX_REAL_THICKNESS real mm (see pial_max).
#   'legacy': EDNiX v1 (overlapping class ranges of min_max, wm.mgz with the rescaled WM intensities, FreeSurfer
#         default 5 FS-mm pial cap whatever the smallbrain rescaling).
SURFACE_METHOD = 'V2h'
MAX_REAL_THICKNESS = 3.0     # mm, real units: above the isocortex thickness of every species of the library

def pial_max(real_mm_per_fs_mm, method=None):
    """-max of mris_make_surfaces (FS mm) for a MAX_REAL_THICKNESS real-mm cap, or None for FreeSurfer's default.
    real_mm_per_fs_mm: real voxel size / FreeSurfer voxel size (1 without the smallbrain header rescale)."""
    if (method or SURFACE_METHOD) == 'legacy':
        return None
    return max(5.0, math.ceil(MAX_REAL_THICKNESS / float(real_mm_per_fs_mm)))


def real_mm_per_fs(ref_file, rescale_json):
    """Real mm per FreeSurfer mm from the smallbrain rescale json (1.0 when the header was not rescaled)."""
    if not opi(rescale_json):
        return 1.0
    with open(rescale_json) as f:
        info = json.load(f)
    if not info.get('change_hd'):
        return 1.0
    return float(ants.image_header_info(ref_file)['spacing'][0]) / float(info['spacing'][0])


def prepa_img(Subname,Ref_file,dir_prepro,aseg_file,labels_dir,midsection,diary_name,method=None):
    """brain (desc-norm), wm (seg-wm) and filled (seg-filled) images for FreeSurfer, by SURFACE_METHOD."""
    if (method or SURFACE_METHOD) == 'V2h':
        return prepa_img_v2h(Subname, Ref_file, dir_prepro, aseg_file, labels_dir, midsection, diary_name)
    return prepa_img_legacy(Subname, Ref_file, dir_prepro, aseg_file, labels_dir, midsection, diary_name)


def _type_norm(Ref_file, diary_name):
    fname = opb(Ref_file)
    stem = fname[:-7] if fname.lower().endswith('.nii.gz') else os.path.splitext(fname)[0]
    suffix = stem.rsplit('_', 1)[-1]
    for token, out in (('T1w', 'T1w'), ('T2w', 'T2w'), ('T1', 'T1w'), ('T2', 'T2w')):
        if token in suffix:
            if token != out:
                run_cmd.msg(f'WARNING: Found "{suffix}", mapping to {out} for surface purposes.', diary_name, 'WARNING')
            return out
    raise ValueError(f"Cannot infer T1/T2 suffix from filename '{Ref_file}'. Expected a chunk containing 'T1' or 'T2'.")


def _filled(aseg, a, hd_BRAIN, midsection, extra_left=(), extra_right=()):
    """filled image: 127 left / 255 right; bilateral labels split at midsection, side decided by the thalamus
    centroid (same rule as the legacy prepa_img)."""
    if midsection == '':
        midsection = math.ceil(int(hd_BRAIN['dimensions'][0]) / 2)
    cent = ants.get_centroids(ants.mask_image(aseg, aseg, 10))
    pos = cent[0][0] * hd_BRAIN['direction'][0][0]
    orig = hd_BRAIN['origin'][0] * hd_BRAIN['direction'][0][0]
    dist = np.abs(pos - orig) / hd_BRAIN['spacing'][0]
    filled = np.zeros(a.shape, np.float32)
    filled[np.isin(a, filled_edit[0] + list(extra_left))] = 127
    filled[np.isin(a, filled_edit[1] + list(extra_right))] = 255
    bil = np.isin(a, filled_edit[2])
    left, right = bil.copy(), bil.copy()
    left[midsection:, :, :] = False
    right[:midsection, :, :] = False
    filled[left] = 127 if dist < midsection else 255
    filled[right] = 255 if dist < midsection else 127
    return filled


def prepa_img_v2h(Subname, Ref_file, dir_prepro, aseg_file, labels_dir, midsection, diary_name):
    """V2h images (see SURFACE_METHOD). Same file names and folders as the legacy prepa_img."""
    run_cmd.msg(' step 10 : synthetic brain / wm / filled images for FreeSurfer (surface method V2h)',
                diary_name, 'HEADER')
    Type_norm = _type_norm(Ref_file, diary_name)
    BRAIN = ants.image_read(Ref_file)
    aseg = ants.image_read(aseg_file)
    hd_BRAIN = ants.image_header_info(Ref_file)
    b = BRAIN.numpy().astype(np.float32)
    a = np.round(aseg.numpy()).astype(np.int32)

    def band(codes, lo, hi):
        """class rescaled into [lo, hi] (1st-99th percentile), keeping its within-class texture"""
        m = np.isin(a, codes)
        o = np.zeros(a.shape, np.float32)
        if m.sum() < 10:
            o[m] = (lo + hi) / 2.
            return o
        v = b[m].astype(np.float64)
        p1, p99 = np.percentile(v, [1, 99])
        o[m] = lo + np.clip((v - p1) / max(p99 - p1, 1e-6), 0, 1) * (hi - lo)
        return o

    norm = band(CSF_edit, 10., 35.) + band(Cereb_GM_edit, 55., 85.)    # hippocampus/amygdala stay in the GM class
    norm[np.isin(a, WM_edit)] = 110.
    norm_img = ants.smooth_image(BRAIN.new_image_like(norm), 0.5, sigma_in_physical_coordinates=False)
    ants.image_write(norm_img, opj(dir_prepro, Subname + '_desc-norm_' + Type_norm + '.nii.gz'), ri=False)

    wm = np.zeros(a.shape, np.float32)
    wm[np.isin(a, wm_aseg_edit)] = 110.
    wm[np.isin(a, aseg_edit)] = 250.
    ants.image_write(BRAIN.new_image_like(wm), opj(labels_dir, Subname + '_seg-wm_dseg.nii.gz'), ri=False)

    filled = _filled(aseg, a, hd_BRAIN, midsection)
    ants.image_write(BRAIN.new_image_like(filled), opj(labels_dir, Subname + '_seg-filled_dseg.nii.gz'), ri=False)
    run_cmd.msg('Done', diary_name, 'OKGREEN')


def prepa_img_legacy(Subname,Ref_file,dir_prepro,aseg_file,labels_dir,midsection,diary_name):

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

    for select,min_im,max_im in zip([CSF_edit,Cereb_GM_edit],
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

    for select, min_im, max_im in zip([CSF_edit, Cereb_GM_edit],
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
    freesurferlist = ['mri', 'surf', 'stats', 'label', 'scripts']
    for folder in freesurferlist:
        if ope(opj(FS_dir, Subname, folder)) == False:
            os.makedirs(opj(FS_dir, Subname, folder))
    if change_hd == 0:
        Ref = ''
    elif change_hd == 1:
        Ref = '_resamp-' + str(scaling)
        list1[0] = opj(prepro_dir, opb(list1[0].replace('.nii.gz', Ref + '.nii.gz')))
        labels_dir  = prepro_dir


    list1[1] = opj(prepro_dir,list1[1] + Ref + '.nii.gz')
    for i in range(2,5):
        list1[i] = opj(labels_dir,list1[i] + Ref + '.nii.gz')
    
    [_, _, fwdFS_cmd, _] = get_orientation.use_ants(list1[0])

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

    for select,name,desc in zip([White_mask,Ventri_mask, Cereb_GM_edit, WBGM_edit],['White','Vent', 'Cerebral_Gray', 'Whole_Brain_Gray'],['white matter','ventricles', 'cerebral gray', 'Whole brain gray']):
        if mskversion == 'aseg':
            aseg = ants.image_read(aseg_img)
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

    # --in_orientation tells mri_convert how the voxels of the file it reads are stored: it must come from that
    # file. It was taken from `ref` (desc-norm T1w); when the label volume is stored with another axis order
    # (e.g. the res-iso copy written by 3dresample -orient in surfFS), the .mgz was mirrored relative to the
    # FreeSurfer volumes and every surface label was sampled at the mirrored position.
    [_, _, fwdFS_cmd, _] = get_orientation.use_ants(file_name)

    cmd = sing_fs + 'mri_convert ' + fwdFS_cmd + file_name + ' ' + opj(FS_dir, animal, 'mri', animal + '_' + label + '.mgz')

    run_cmd.run(cmd,diary_name)







