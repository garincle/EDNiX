import os
import shutil
import glob

opj = os.path.join
ope = os.path.exists
opi = os.path.isfile
opd = os.path.dirname
opb = os.path.basename

from Tools import run_cmd
from Tools import getpath
from Tools import QC_plot
from Tools import check_nii

from anatomical import norm2template

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
            w2fwd2      = [False]


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

def apply(ID,volumes_dir,masks_dir,labels_dir,bids_dir,info,listTimage,targetsuffix,
          n_for_ANTS,list_atlas,path_label_code,type_norm,diary_file,sing_afni,sing_wb):
    '''
    to debug if necessary :  structure of info should be (check the function just above ....):
    info[i] = [sourcename,sourcefile,sourcemask,sourcelabeldir,target_voldir,target_maskdir,target_labeldir,
    transfoinv,w2inv,transfofwd,w2fwd]
    '''
    print(info)
    for i in range(len(info)):
        for Timage in listTimage:
            img_ref = opj(volumes_dir, '_'.join([ID,targetsuffix,Timage]) + '.nii.gz')
            imgtemplate = opj(volumes_dir, '_'.join([ID,'space-acpc','desc-template',Timage]) + '.nii.gz')
            if opi(img_ref):
                if Timage == type_norm:
                    # the brainmask
                    if opi(info[0][2]):
                        # backward
                        norm2template.apply('acpc', img_ref, masks_dir, ID, info[0][2],
                                            info[i][7], info[i][8], '', '', '', Timage, n_for_ANTS,
                                            diary_file, sing_wb)

                        check_nii.resamp(opj(masks_dir, ID + '_space-acpc_mask.nii.gz'), img_ref, 'msk', '', '', diary_file,
                                         sing_wb)
                        #last skull stripping
                        cmd = (sing_afni + '3dcalc -overwrite -a ' + imgtemplate + ' -b ' + opj(masks_dir, ID + '_space-acpc_mask.nii.gz') +
                               ' -prefix ' + opj(volumes_dir, '_'.join([ID, 'space-acpc', 'desc-SS', Timage]) + '.nii.gz') + ' -expr "a*b"')
                        run_cmd.do(cmd, diary_file)


                # forward
                norm2template.apply(info[i][0], info[i][1], info[i][4], ID, img_ref,
                                    info[i][9], info[i][10], 'SS', '', '', Timage, n_for_ANTS,
                                    diary_file, sing_wb)

                # backward
                oldname = opb(info[i][1]).split('_')
                print('oldname=' + str(oldname))

                newname = opj(opd(info[i][1]), '_'.join(oldname[0:-1] + [Timage]) + '.nii.gz')
                print('newname=' + str(newname))
                print('ope newname = ' + str(ope(newname)))
                norm2template.apply('acpc', img_ref, info[i][4], info[i][0], newname,
                                    info[i][7], info[i][8], 'SS', '', '', Timage, n_for_ANTS,
                                    diary_file, sing_wb)

                QC_plot.mosaic(info[i][1],
                               opj(info[i][4], '_'.join([ID, 'space-' + info[i][0], 'desc-SS', Timage]) + '.nii.gz'),
                               opj(bids_dir, 'QC', '_'.join(['backtonative_native_in_space-' + info[i][0], Timage]),
                               '_'.join([ID, 'space-' + info[i][0], Timage]) + '.png'))

                QC_plot.mosaic(img_ref,
                               opj(info[i][4], '_'.join([info[i][0], 'space-acpc', 'desc-SS', Timage]) + '.nii.gz'),
                               opj(bids_dir, 'QC', '_'.join([info[i][0], 'space-acpc', Timage]),
                               '_'.join([info[i][0], 'space-acpc', Timage]) + '.png'))

        # myelin image
        if opi(opj(volumes_dir, ID + targetsuffix + '_T1wdividedbyT2w.nii.gz')):
            # forward
            norm2template.apply(info[i][0], info[i][1], info[i][4], ID, opj(volumes_dir, ID + targetsuffix + '_T1wdividedbyT2w.nii.gz'),
                                info[i][9], info[i][10], 'desc-SS', '', '', 'T1wdividedbyT2w', n_for_ANTS,
                                diary_file, sing_wb)
            check_nii.resamp(opj(masks_dir, ID + '_space-' + info[i][0]+ '_T1wdividedbyT2w.nii.gz'), info[i][1], 'T1w', '', '', diary_file,
                             sing_wb)

    # the aseg image
    if opi(info[0][3] + '4FS_dseg.nii.gz'):
        # backward
        norm2template.apply('acpc', img_ref, labels_dir,  ID, info[0][3] + '4FS_dseg.nii.gz',
                            info[0][7], info[0][8], '4FS', path_label_code, 'FreeSurfer', '', '',
                            diary_file, sing_wb)
        check_nii.resamp(opj(labels_dir, ID + '_space-acpc_seg-4FS_dseg.nii.gz'), img_ref, 'label',
                         path_label_code, 'FreeSurfer', diary_file, sing_wb)
    # the atlas images
    for i in range(len(list_atlas[0])):
        # backward
        if opi(info[0][3] + list_atlas[0][i] + '_dseg.nii.gz'):
            norm2template.apply('acpc', img_ref, labels_dir,ID, info[0][3] + list_atlas[0][i] + '_dseg.nii.gz',
                                 info[0][7], info[0][8], list_atlas[0][i], path_label_code, list_atlas[0][i], '', '',
                                diary_file, sing_wb)
            check_nii.resamp_no_check(opj(labels_dir, ID + '_space-acpc_seg-' + list_atlas[0][i] + '_dseg.nii.gz'), img_ref, 'label',
                             path_label_code, list_atlas[0][i], diary_file, sing_wb)


# remove the _space_acpc_ section of the name
    for dir,suffix in zip([masks_dir,labels_dir],['mask','dseg']):
        listimg = sorted(glob.glob(opj(dir, ID + '_space-acpc_*_' + suffix + '*')))

        for i in range(len(listimg)):
            name = opb(listimg[i]).split('_')
            newname = '_'.join([name[0]] + name[2:])
            shutil.move(listimg[i], opj(opd(listimg[i]), newname))


