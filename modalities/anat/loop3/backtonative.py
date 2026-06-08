import os
import shutil
import glob
import json

opj = os.path.join
ope = os.path.exists
opi = os.path.isfile
opd = os.path.dirname
opb = os.path.basename

from Tools import run_cmd
from Tools import getpath
from Tools import QC_plot
from Tools import check_nii

from modalities.anat import norm2template

# Single source of truth for the backtonative mask path — defined in
# _200_Data_QC so both the ITK-SNAP step and the apply() function always
# refer to the same filename convention.
from modalities.anat.loop3._200_Data_QC import _backtonative_mask_path

'''          ┌────────────┐
             │  Reference │  (e.g. EDNiX)
             └────┬───────┘
                  │
          (native → reference)
                  │
        ┌─────────┴─────────┐
        │                   │
   [Longitudinal]      [Study Template]
  (refSession)       (studyTemplate)
        │                   │
        └─────────┬─────────┘
                  │
            Native space
'''

def get(ID,datapath,bids_dir,Session,max_ses,suffix,type_norm,BASE_SS, BASE_mask,BASE_atlas_folder,study_template_atlas_folder,
        creat_study_template,coregistration_longitudinal,reference,transfoS,species,diary_file):
    """
    TRANSFORMATION PACKAGE STRUCTURE:

    Each package contains 11 elements organized as follows:

    package = [
        ref_name,           # 1 # str: Reference space name ('template', 'refSession', 'studyTemplate')
        source_image,       # 2 # str: Path to source image to transform FROM
        source_mask,        # 3 # str: Path to source mask/brain mask
        label_base_path,    # 4 # str: Base path for segmentation labels (append region names)
        template_volume,    # 5 # str: Path to target template volume
        template_masks,     # 6 # str: Directory for template-space masks
        template_labels,    # 7 # str: Directory for template-space labels
        transfoinv,         # 8 # list: Inverse transformation files [mat, nii.gz...]
        w2inv,              # 9 # list: Boolean flags for using warps in inverse [True, False...]
        transfofwd,         # 10 # list: Forward transformation files [nii.gz, mat...]
        w2fwd               # 11 # list: Boolean flags for using warps in forward [False, True...]]

    SPECIFIC EXAMPLES:

    1. NATIVE → TEMPLATE (package1 - forward transformation):
       ref_name = 'template'
       transfoinv = ['native_to_template_0GenericAffine.mat', 'native_to_template_1InverseWarp.nii.gz']
       w2inv = [True, False]  # Use affine, then apply inverse warp
       transfofwd = ['native_to_template_1Warp.nii.gz', 'native_to_template_0GenericAffine.mat']
       w2fwd = [False, False] # Apply warp, then affine

    2. NATIVE → REFSESSION (package2 - longitudinal):
       ref_name = 'refSession'

    3. NATIVE → STUDYTEMPLATE (package3):
       ref_name = 'studyTemplate'

    USAGE PATTERNS:
    Forward Transformation (Native → Template):
    result = apply_transforms(
        fixed_image=package[4],      # template_volume
        moving_image=package[1],     # source_image
        transform_list=package[8],   # transfofwd
        warp_flags=package[10]       # w2fwd)
    """

    nl = 'set every step of transformation'
    run_cmd.msg(nl, diary_file, 'HEADER')

    (path_anat,dir_transfo,FS_dir,dir_prepro,dir_native,volumes_dir,labels_dir,masks_dir,
     wb_template_dir,wb_template_vol,wb_template_labels,wb_template_masks,
     _,_,_,_,wb_studytemplate_dir,wb_studytemplate_vol,wb_studytemplate_labels,wb_studytemplate_masks,
     wb_refsession_dir,wb_refsession_vol,wb_refsession_labels,wb_refsession_masks) = getpath.anat(datapath,reference,'',
                                                                                                  creat_study_template,
                                                                                                  coregistration_longitudinal,
                                                                                                  'all')
    TDir = dir_transfo
    transfoinv2 = ''
    transfofwd2 = ''
    w2inv2 = ''
    w2fwd2 = ''
    transfoinv3 = ''
    transfofwd3 = ''
    w2inv3 = ''
    w2fwd3 = ''
    transfoinv4 = ''
    transfofwd4 = ''
    w2inv4 = ''
    w2fwd4 = ''
    ref1          = reference
    sourceimage1  = BASE_SS
    sourcemask1   = BASE_mask
    labelname1    = opj(BASE_atlas_folder, species + '_seg-')
    transforname1 = ['native', 'to', ref1, transfoS]
    package1 = [ref1, sourceimage1, sourcemask1, labelname1, wb_template_vol, wb_template_masks, wb_template_labels]

    if coregistration_longitudinal == True:
        data_path_max = opj(bids_dir, 'sub-' + ID, 'ses-' + str(max_ses))
        (_, transfo_refSession, _, _, _, volumes_refSession, label_refSession, mask_refSession) = getpath.anat(data_path_max, reference,
                                                                     '', False,
                                                                     True, 'native')
        TDir = transfo_refSession
        ref2          = 'refSession'
        sourceimage2  = opj(volumes_refSession, '_'.join([ID, suffix, type_norm]) + '.nii.gz')
        sourcemask2   = opj(mask_refSession, ID + '_mask.nii.gz')
        labelname2    = opj(label_refSession, ID + '_seg-')
        transforname2 = ['native', 'to', ref2, transfoS]

        package2 = [ref2, sourceimage2, sourcemask2, labelname2, wb_refsession_vol, wb_refsession_masks,
                    wb_refsession_labels]

    if creat_study_template == True:
        studyacpc_dir, studyprepro_dir, studytransfo_dir, studyvolume_dir, studylabels_dir, studymasks_dir,_,_,_,_,_,_,_,_ = getpath.stytemplate(
            study_template_atlas_folder, reference, '')
        print(studyvolume_dir)
        ref3          = 'studyTemplate'
        sourceimage3 = opj(studyvolume_dir, '_'.join([ref3, type_norm]) + '.nii.gz')
        print(sourceimage3)
        sourcemask3   = opj(studymasks_dir,  ref3 + '_mask.nii.gz')
        labelname3    = opj(studylabels_dir, ref3 + '_seg-')
        transforname3 = ['native', 'to', ref3, transfoS]
        transforname4 = ['studyTemplate', 'to', ref1, transfoS]

        package3      = [ref3, sourceimage3, sourcemask3, labelname3, wb_studytemplate_vol, wb_studytemplate_masks,
                         wb_studytemplate_labels]
    # set each step of transformation ..................................................................................
    # from native to template
    if opi(opj(TDir, '_'.join(transforname1 + ['1InverseWarp.nii.gz']))):
        transfoinv1 = [opj(TDir, '_'.join(transforname1 + ['0GenericAffine.mat'])),
                       opj(TDir, '_'.join(transforname1 + ['1InverseWarp.nii.gz']))]
        w2inv1      = [True, False]

        transfofwd1 = [opj(TDir, '_'.join(transforname1 + ['1Warp.nii.gz'])),
                       opj(TDir, '_'.join(transforname1 + ['0GenericAffine.mat'])),]
        w2fwd1      = [False, False]
    else:
        transfoinv1 = [opj(TDir, '_'.join(transforname1 + ['0GenericAffine.mat']))]
        w2inv1      = [True]
        transfofwd1 = transfoinv1
        w2fwd1      = [False]

    if coregistration_longitudinal == True and Session != max_ses:
        if opi(opj(dir_transfo, '_'.join(transforname2 + ['1InverseWarp.nii.gz']))):
            transfoinv2 = [opj(dir_transfo, '_'.join(transforname2 + ['0GenericAffine.mat'])),
                           opj(dir_transfo, '_'.join(transforname2 + ['1InverseWarp.nii.gz']))]
            w2inv2      = [True, False]
            transfofwd2 = [opj(dir_transfo, '_'.join(transforname2 + ['1Warp.nii.gz'])),
                           opj(dir_transfo, '_'.join(transforname2 + ['0GenericAffine.mat'])), ]
            w2fwd2      = [False, False]
        else:
            transfoinv2 = [opj(dir_transfo, '_'.join(transforname1 + ['0GenericAffine.mat']))]
            w2inv2      = [True]
            transfofwd2 = transfoinv2
            w2fwd2      = [False]

    # from native to study template
    if creat_study_template == True:
        if opi(opj(TDir, '_'.join(transforname3 + ['1InverseWarp.nii.gz']))):
            transfoinv3 = [opj(TDir, '_'.join(transforname3 + ['0GenericAffine.mat'])),
                           opj(TDir, '_'.join(transforname3 + ['1InverseWarp.nii.gz']))]
            w2inv3      = [True, False]

            transfofwd3 = [opj(TDir, '_'.join(transforname3 + ['1Warp.nii.gz'])),
                           opj(TDir, '_'.join(transforname3 + ['0GenericAffine.mat'])),]
            w2fwd3      = [False, False]
        else:
            transfoinv3 = [opj(TDir, '_'.join(transforname3 + ['0GenericAffine.mat']))]
            w2inv3      = [True]
            transfofwd3 = transfoinv3
            w2fwd3      = [False]

        if opi(opj(studytransfo_dir, '_'.join(transforname4 + ['1InverseWarp.nii.gz']))):
            transfoinv4 = [opj(studytransfo_dir, '_'.join(transforname4 + ['0GenericAffine.mat'])),
                           opj(studytransfo_dir, '_'.join(transforname4 + ['1InverseWarp.nii.gz']))]
            w2inv4      = [True, False]
            transfofwd4 = [opj(studytransfo_dir, '_'.join(transforname4 + ['1Warp.nii.gz'])),
                           opj(studytransfo_dir, '_'.join(transforname4 + ['0GenericAffine.mat'])), ]
            w2fwd4      = [False, False]
        else:
            transfoinv4 = [opj(studytransfo_dir, '_'.join(transforname4 + ['0GenericAffine.mat']))]
            w2inv4      = [True]
            transfofwd4 = transfoinv4
            w2fwd4      = [False]


    ####################################################################################################################
    # adapt for each situation
    if creat_study_template == False:
        if coregistration_longitudinal == False:
            # classical situation
            transfoinv = transfoinv1
            w2inv      = w2inv1
            transfofwd = transfofwd1
            w2fwd      = w2fwd1

            package1.extend([transfoinv,w2inv,transfofwd,w2fwd])
            package = [package1]

        else:
            # longitudinal study
            if Session == max_ses:
                transfoinv = transfoinv1
                w2inv      = w2inv1
                transfofwd = transfofwd1
                w2fwd      = w2fwd1
                package1.extend([transfoinv, w2inv, transfofwd, w2fwd])
                package = [package1]
            else:
                transfoinv = transfoinv2 + transfoinv1
                w2inv      = w2inv1 + w2inv2
                transfofwd = transfofwd1 + transfofwd2
                w2fwd      = w2fwd1 + w2fwd2

                package1.extend([transfoinv, w2inv, transfofwd, w2fwd])
                package2.extend([transfoinv2, w2inv2, transfofwd2, w2fwd2])
                package = [package1,package2]

    else:
        if coregistration_longitudinal == False:
            transfoinv = transfoinv3 + transfoinv4
            w2inv      = w2inv3 + w2inv4
            transfofwd = transfofwd4 + transfofwd3
            print(w2fwd4)
            print(w2fwd3)
            w2fwd      =  w2fwd4 + w2fwd3
            package1.extend([transfoinv, w2inv, transfofwd, w2fwd])
            package3.extend([transfoinv3, w2inv3, transfofwd3, w2fwd3])
            package = [package1, package3]

        else:
            # longitudinal study
            if Session == max_ses:
                transfoinv = transfoinv3 + transfoinv4
                w2inv      = w2inv3 + w2inv4
                transfofwd = transfofwd4 + transfofwd3
                w2fwd      = w2fwd4 + w2fwd3
                package1.extend([transfoinv, w2inv, transfofwd, w2fwd])
                package3.extend([transfoinv3, w2inv3, transfofwd3, w2fwd3])
                package = [package1,package3]

            else:
                transfoinv = transfoinv2 + transfoinv3 + transfoinv4
                w2inv      = w2inv2 + w2inv3 + w2inv4
                transfofwd = transfofwd4 + transfofwd3 + transfofwd2
                w2fwd      = w2fwd4 + w2fwd3 + w2fwd2
                package1.extend([transfoinv, w2inv, transfofwd, w2fwd])

                transfoinv3new = transfoinv2 + transfoinv3
                w2inv3new      = w2inv2 + w2inv3
                transfofwd3new = transfofwd3 + transfofwd2
                w2fwd3new      = w2fwd3 + w2fwd2

                package3.extend([transfoinv3new, w2inv3new, transfofwd3new, w2fwd3new])
                package2.extend([transfoinv2, w2inv2, transfofwd2, w2fwd2])
                package = [package1,package3,package2]
    return package


# =============================================================================
# Backtonative mask helpers
# =============================================================================

# _backtonative_mask_path is imported from anat.loop3._200_Data_QC above.
# Do not redefine it here — _200_Data_QC is the single source of truth.


def _final_ss_path(volumes_dir, ID, Timage):
    """Canonical path for the final skull-stripped anatomy."""
    return opj(volumes_dir, '_'.join([ID, 'space-acpc', 'desc-SS', Timage]) + '.nii.gz')


def _needs_regeneration(mask_path, ss_path, diary_file):
    """
    Return True if the final skull-stripped image needs to be regenerated.

    Regeneration is required when:
      - The final SS image does not exist yet, OR
      - The backtonative mask is newer than the final SS image (the mask was
        edited after the last run).

    Parameters
    ----------
    mask_path : str
        Path to the mask that will be used (backtonative or automatic).
    ss_path : str
        Path to the final skull-stripped anatomy.
    diary_file : str
        Diary log path.

    Returns
    -------
    bool
    """
    if not opi(ss_path):
        run_cmd.msg(
            f'INFO: final SS image not found — will generate: {opb(ss_path)}',
            diary_file, 'OKGREEN')
        return True

    mask_mtime = os.path.getmtime(mask_path)
    ss_mtime   = os.path.getmtime(ss_path)

    if mask_mtime > ss_mtime:
        import datetime
        mask_dt = datetime.datetime.fromtimestamp(mask_mtime).strftime('%Y-%m-%d %H:%M:%S')
        ss_dt   = datetime.datetime.fromtimestamp(ss_mtime).strftime('%Y-%m-%d %H:%M:%S')
        run_cmd.msg(
            f'INFO: backtonative mask ({mask_dt}) is newer than final SS image '
            f'({ss_dt}) — regenerating {opb(ss_path)}.',
            diary_file, 'OKGREEN')
        return True

    run_cmd.msg(
        f'INFO: final SS image is up to date — skipping regeneration of '
        f'{opb(ss_path)}.',
        diary_file, 'OKGREEN')
    return False


def _apply_final_skull_strip(imgtemplate, mask_path, ss_path,
                              sing_afni, diary_file, mask_source_label):
    """
    Apply a brain mask to the template image to produce the final SS anatomy.

    Parameters
    ----------
    imgtemplate : str
        Path to the template image (space-acpc_desc-template_T1w).
    mask_path : str
        Path to the mask to apply.
    ss_path : str
        Output path for the final skull-stripped image.
    sing_afni : str
        AFNI Singularity prefix.
    diary_file : str
        Diary log path.
    mask_source_label : str
        Human-readable label describing where the mask came from, used in the
        JSON sidecar ('backtonative mask' or 'automatic ACPC mask').
    """
    cmd = (sing_afni + '3dcalc -overwrite'
           ' -a ' + imgtemplate +
           ' -b ' + mask_path +
           ' -prefix ' + ss_path +
           ' -expr "a*b"')
    run_cmd.do(cmd, diary_file)

    sidecar = {
        'Sources':     [imgtemplate, mask_path],
        'Description': f'Final skull-stripped anatomy (mask source: {mask_source_label}).',
        'Command':     cmd,
    }
    with open(ss_path.replace('.nii.gz', '.json'), 'w') as f:
        json.dump(sidecar, f, indent=3)


# =============================================================================
# apply()
# =============================================================================

def apply(ID, volumes_dir, masks_dir, labels_dir, bids_dir, info, listTimage,
          targetsuffix, n_for_ANTS, list_atlas, path_label_code, type_norm,
          diary_file, sing_afni, sing_wb):
    '''
    to debug if necessary: structure of info should be (check the function just above ...):
    info[i] = [sourcename, sourcefile, sourcemask, sourcelabeldir, target_voldir,
               target_maskdir, target_labeldir, transfoinv, w2inv, transfofwd, w2fwd]
    '''
    print(info)
    for i in range(len(info)):
        for Timage in listTimage:
            img_ref    = opj(volumes_dir, '_'.join([ID, targetsuffix, Timage]) + '.nii.gz')
            imgtemplate = opj(volumes_dir, '_'.join([ID, 'space-acpc', 'desc-template', Timage]) + '.nii.gz')

            if opi(img_ref):
                if Timage == type_norm:

                    # ------------------------------------------------------------------
                    # Mask selection: backtonative mask takes priority over automatic
                    # ACPC mask when it exists.
                    #
                    # The backtonative mask ({ID}_space-acpc_desc-backtonative_mask.nii.gz)
                    # is created by the user via ITK-SNAP (step itk_3). It represents
                    # the manually corrected final mask for this subject/session.
                    #
                    # Priority logic:
                    #   1. backtonative mask exists  → use it unconditionally
                    #   2. backtonative mask absent  → use the automatic ACPC mask
                    #      propagated back from the template by norm2template.apply()
                    #
                    # Regeneration logic (checked via _needs_regeneration):
                    #   - If the final SS image does not yet exist → generate it.
                    #   - If the chosen mask is newer than the final SS image
                    #     (i.e. edited after the last pipeline run) → regenerate.
                    #   - Otherwise → skip to avoid unnecessary recomputation.
                    # ------------------------------------------------------------------

                    backtonative_mask = _backtonative_mask_path(masks_dir, ID)
                    acpc_mask         = opj(masks_dir, ID + '_space-acpc_mask.nii.gz')
                    ss_path           = _final_ss_path(volumes_dir, ID, Timage)

                    if opi(backtonative_mask):
                        # ---- Use manually edited backtonative mask ----------------
                        run_cmd.msg(
                            f'INFO: backtonative mask found for {ID} — '
                            f'using it for final skull stripping: '
                            f'{opb(backtonative_mask)}',
                            diary_file, 'OKGREEN')

                        if _needs_regeneration(backtonative_mask, ss_path, diary_file):
                            # Resample backtonative mask to match imgtemplate grid
                            # (safety check: user may have saved at different resolution)
                            check_nii.resamp(backtonative_mask, imgtemplate,
                                             'msk', '', '', diary_file, '')

                            _apply_final_skull_strip(
                                imgtemplate, backtonative_mask, ss_path,
                                sing_afni, diary_file,
                                mask_source_label='backtonative mask (manual)')

                    else:
                        # ---- Use automatic ACPC mask --------------------------------
                        run_cmd.msg(
                            f'INFO: no backtonative mask found for {ID} — '
                            f'using automatic ACPC mask.',
                            diary_file, 'OKGREEN')

                        if opi(info[0][2]):
                            # Propagate mask from template space back to ACPC space
                            norm2template.apply(
                                'acpc', img_ref, masks_dir, ID, info[0][2],
                                info[0][7], info[0][8], '', '', '', Timage,
                                n_for_ANTS, diary_file, sing_wb)

                            check_nii.resamp(acpc_mask, img_ref,
                                             'msk', '', '', diary_file, sing_wb)

                            if _needs_regeneration(acpc_mask, ss_path, diary_file):
                                _apply_final_skull_strip(
                                    imgtemplate, acpc_mask, ss_path,
                                    sing_afni, diary_file,
                                    mask_source_label='automatic ACPC mask')

                # forward normalisation to reference space
                norm2template.apply(
                    info[i][0], info[i][1], info[i][4], ID, img_ref,
                    info[i][9], info[i][10], 'SS', '', '', Timage,
                    n_for_ANTS, diary_file, sing_wb)

                # backward normalisation (template → ACPC)
                oldname = opb(info[i][1]).split('_')
                print('oldname=' + str(oldname))

                newname = opj(opd(info[i][1]), '_'.join(oldname[0:-1] + [Timage]) + '.nii.gz')
                print('newname=' + str(newname))
                print('ope newname = ' + str(ope(newname)))

                if ope(newname):
                    norm2template.apply(
                        'acpc', img_ref, info[i][4], info[i][0], newname,
                        info[i][7], info[i][8], 'SS', '', '', Timage,
                        n_for_ANTS, diary_file, sing_wb)

                    QC_plot.mosaic(
                        info[i][1],
                        opj(info[i][4], '_'.join([ID, 'space-' + info[i][0], 'desc-SS', Timage]) + '.nii.gz'),
                        opj(bids_dir, 'QC',
                            '_'.join(['backtonative_native_in_space-' + info[i][0], Timage]),
                            '_'.join([ID, 'space-' + info[i][0], Timage]) + '.png'))

                    QC_plot.mosaic(
                        img_ref,
                        opj(info[i][4], '_'.join([info[i][0], 'space-acpc', 'desc-SS', Timage]) + '.nii.gz'),
                        opj(bids_dir, 'QC',
                            '_'.join([info[i][0], 'space-acpc', Timage]),
                            '_'.join([info[i][0], 'space-acpc', Timage]) + '.png'))

                elif opi(img_ref) == False and Timage == type_norm:
                    ValueError("Can't find the template name in the library, check naming and path!")
                else:
                    print("no " + str(Timage) + " template found, not a big deal")

        # myelin image
        if opi(opj(volumes_dir, ID + targetsuffix + '_T1wdividedbyT2w.nii.gz')):
            norm2template.apply(
                info[i][0], info[i][1], info[i][4], ID,
                opj(volumes_dir, ID + targetsuffix + '_T1wdividedbyT2w.nii.gz'),
                info[i][9], info[i][10], 'desc-SS', '', '', 'T1wdividedbyT2w',
                n_for_ANTS, diary_file, sing_wb)
            check_nii.resamp(
                opj(masks_dir, ID + '_space-' + info[i][0] + '_T1wdividedbyT2w.nii.gz'),
                info[i][1], 'T1w', '', '', diary_file, sing_wb)

    # the aseg image
    if opi(info[0][3] + '4FS_dseg.nii.gz'):
        norm2template.apply(
            'acpc', img_ref, labels_dir, ID, info[0][3] + '4FS_dseg.nii.gz',
            info[0][7], info[0][8], '4FS', path_label_code, 'FreeSurfer', '', '',
            diary_file, sing_wb)
        check_nii.resamp(
            opj(labels_dir, ID + '_space-acpc_seg-4FS_dseg.nii.gz'),
            img_ref, 'label', path_label_code, 'FreeSurfer', diary_file, sing_wb)
    else:
        masks = glob.glob(opd(info[0][2]) + '/*desc-*_mask.nii.gz')
        descriptors = []
        for mask_path in masks:
            filename   = os.path.basename(mask_path)
            desc_part  = filename.split('desc-')[1].split('_mask')[0]
            descriptors.append(desc_part)
            norm2template.apply(
                'acpc', img_ref, masks_dir, ID, mask_path,
                info[0][7], info[0][8], desc_part, '', '', Timage,
                n_for_ANTS, diary_file, sing_wb)
            filename = opj(masks_dir, ID + '_space-acpc' + desc_part + '_mask.nii.gz')
            new_filename = filename.replace('space-acpc', 'desc-')
            os.rename(filename, new_filename)
            filename = opj(masks_dir, ID + '_space-acpc' + desc_part + '_mask.json')
            new_filename = filename.replace('space-acpc', 'desc-')
            os.rename(filename, new_filename)

    # the atlas images
    for i in range(len(list_atlas[0])):
        if opi(info[0][3] + list_atlas[0][i] + '_dseg.nii.gz'):
            norm2template.apply(
                'acpc', img_ref, labels_dir, ID,
                info[0][3] + list_atlas[0][i] + '_dseg.nii.gz',
                info[0][7], info[0][8], list_atlas[0][i],
                path_label_code, list_atlas[0][i], '', '',
                diary_file, sing_wb)
            check_nii.resamp_no_check(
                opj(labels_dir, ID + '_space-acpc_seg-' + list_atlas[0][i] + '_dseg.nii.gz'),
                img_ref, 'label', path_label_code, list_atlas[0][i],
                diary_file, sing_wb)

    # remove the _space-acpc_ section from mask and label filenames
    for directory, suffix in zip([masks_dir, labels_dir], ['mask', 'dseg']):
        listimg = sorted(glob.glob(opj(directory, ID + '_space-acpc_*_' + suffix + '*')))
        for i in range(len(listimg)):
            name    = opb(listimg[i]).split('_')
            newname = '_'.join([name[0]] + name[2:])
            shutil.move(listimg[i], opj(opd(listimg[i]), newname))
