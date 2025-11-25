import os

opj = os.path.join
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd
from Tools import QC_plot
from Tools import getpath
from anat import norm2template

def nativetoTemplate(sourcename,ID,listTimage, data_path,
                    bids_dir,creat_study_template,coregistration_longitudinal,reference,targetsuffix,sourceimage,
                    list_transfo,type_norm,sing_wb,diary_file):

    nl = 'Run anat._6_stdyTmax.nativetostudyT'
    run_cmd.msg(nl, diary_file,'HEADER')

    refnb = 0
    for i, j in enumerate(list_transfo):
        if list_transfo[i]["name"] == 'coreg':
            refnb = i


    # important variables : check if it is ok .........................................................................
    targetname  = 'native'
    transfoname = [targetname, 'to', sourcename]
    transfoT    = 'Shift'
    transfoS    = 'Final'
    suffix      = 'SS'
    #..................................................................................................................


    if any(sourcename in word for word in ['studyTemplate','refSession']):
        arg = sourcename
        (_, dir_transfo, _, _, _, volumes_dir, _, _,
         _, template_vol, _, _) = getpath.anat(data_path, reference, '', creat_study_template,
                                                coregistration_longitudinal, arg)
    else:
        arg = 'template'
        (_, dir_transfo, _, _, _, volumes_dir, _, _,
         _, template_vol, _, _, _, _, _, _) = getpath.anat(data_path, reference, '', creat_study_template,
                                                           coregistration_longitudinal, arg)

    targefile    = opj(volumes_dir, ID + targetsuffix + type_norm + '.nii.gz')

    transfonameR = list(reversed(transfoname))
    transfonameT = '_'.join(transfoname + [transfoT])
    transfonameS = '_'.join(transfoname + [transfoS])

    if not ope(dir_transfo):
        os.makedirs(dir_transfo)
    if not ope(template_vol):
        os.makedirs(template_vol)

    #######################################################################
    ############### coregistration                      ###################
    #######################################################################

    nl = 'INFO: Co-registration ready'
    run_cmd.msg(nl, diary_file, 'OKGREEN')
    nl = 'INFO: WE WILL COREGISTER: ' + targefile + '(the anat img)' + ' to ' + sourceimage + '(the reference)'
    run_cmd.msg(nl, diary_file, 'OKGREEN')
    nl = 'INFO: type_of_transform=' + list_transfo[refnb]["type_of_transform"]
    run_cmd.msg(nl, diary_file, 'OKGREEN')

    mTx = norm2template.norm(targetname,targefile,'',template_vol,
                             sourcename,sourceimage,'',dir_transfo,
                             list_transfo[refnb]["type_of_transform"],transfonameT,transfonameS,
                             list_transfo[refnb]["affmetricT"], list_transfo[refnb]["affmetric"], list_transfo[refnb]["interpol"],
                             diary_file,sing_wb,'_' + suffix + '_' + type_norm,0)

    if opi(opj(dir_transfo,transfonameS + '_1Warp.nii.gz')):
        trlisti = [opj(dir_transfo,transfonameS + '_0GenericAffine.mat'),
                   opj(dir_transfo,transfonameS + '_1InverseWarp.nii.gz')]
        w2i   = [True,False]

        w2fwd = [False,False]
        trlistf = [opj(dir_transfo, transfonameS + '_1Warp.nii.gz'),
            opj(dir_transfo, transfonameS + '_0GenericAffine.mat')]

    else:
        trlisti = trlistf = [opj(dir_transfo, transfonameS + '_0GenericAffine.mat')]
        w2i   = [True]
        w2fwd = [False]

    for Timage in listTimage:

        img_ref = opj(volumes_dir, ID + targetsuffix + Timage + '.nii.gz')

        norm2template.apply(sourcename, sourceimage, template_vol, targetname, img_ref,
                            trlistf, w2fwd,suffix, '', '', Timage,
                            list_transfo[refnb]["interpol"], diary_file, sing_wb)

        norm2template.apply(targetname, img_ref, template_vol, sourcename, sourceimage,
                            trlisti, w2i, suffix, '', '', Timage,
                            list_transfo[refnb]["interpol"], diary_file, sing_wb)

        ### plot the QC
        QC_plot.mosaic(sourceimage,
                       opj(template_vol, '_'.join([targetname, 'space-' + sourcename, 'desc-' + suffix,type_norm]) + '.nii.gz'),
                       opj(bids_dir, 'QC', '_'.join(transfoname + [Timage]),
                       '_'.join([ID] + transfoname + [Timage]) + '.png'))

        QC_plot.mosaic(img_ref,
                       opj(template_vol, '_'.join([sourcename, 'space-' + targetname, 'desc-' + suffix,type_norm]) + '.nii.gz'),
                       opj(bids_dir, 'QC', '_'.join(transfonameR + [Timage]),
                       '_'.join([ID] + transfonameR + [Timage]) + '.png'))







