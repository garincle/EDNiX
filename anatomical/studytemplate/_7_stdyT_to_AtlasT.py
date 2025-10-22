##############################################################################
####### Co-register your study template to the atlas template ################
##############################################################################

import os

opj = os.path.join
opb = os.path.basename
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd
from Tools import getpath
from anatomical import norm2template

def stdyT_to_AtlasT(aff_metric_ants_Transl_template, list_atlases, BASE_SS, BASE_atlas_folder,BASE_mask,species,
                    n_for_ANTS, aff_metric_ants, study_template_atlas_folder, stdy_template, type_of_transform_stdyT,path_label_code,
                    diary_file,reference,sing_wb):

    nl = 'Run anatomical._7_stdyT_to_AtlasT.stdyT_to_AtlasT'
    run_cmd.msg(nl, diary_file,'HEADER')

    targetname = 'studyTemplate'
    sourcename = 'EDNiX'

    transfoname = [targetname, 'to', sourcename]
    transfoT    = 'Shift'
    transfoS    = 'Final'

    (_, _, dir_transfo, _, labels_dir, masks_dir,_, wb_template_vol,
     _, _,_, _, _, _) = getpath.stytemplate(study_template_atlas_folder,reference, '')

    transfonameR = list(reversed(transfoname))
    transfonameT = '_'.join(transfoname + [transfoT])
    transfonameS = '_'.join(transfoname + [transfoS])


    if not ope(dir_transfo):
        os.makedirs(dir_transfo)
    if not ope(wb_template_vol):
        os.makedirs(wb_template_vol)
    if not ope(labels_dir):
        os.makedirs(labels_dir)

    #######################################################################
    ############### coregistration                      ###################
    #######################################################################

    nl = 'INFO: Co-registration ready'
    run_cmd.msg(nl, diary_file, 'OKGREEN')
    nl = 'INFO: WE WILL COREGISTER: ' + stdy_template + '(the anat img)' + ' to ' + BASE_SS + '(the reference)'
    run_cmd.msg(nl, diary_file, 'OKGREEN')
    nl = 'INFO: type_of_transform=' + type_of_transform_stdyT
    run_cmd.msg(nl, diary_file, 'OKGREEN')


    mTx = norm2template.norm(targetname,stdy_template,'',wb_template_vol,
                             sourcename,BASE_SS,'',dir_transfo,
                             type_of_transform_stdyT,transfonameT,transfonameS,
                             aff_metric_ants_Transl_template,aff_metric_ants,n_for_ANTS,
                             diary_file,sing_wb,'',0)

    if opi(opj(dir_transfo,transfonameS + '_1Warp.nii.gz')):
        w2i = [True, False]
        w2fwd = [False, False]
    else:
        w2i = [True]
        w2fwd = [False]

    ####################################################################################
    ########################## seg into the native space ###################
    sourcename = 'acpc'
    for atlas in list_atlases[0]:
        atlas = str(atlas)
        atlasfile = opj(BASE_atlas_folder, species + '_seg-' + atlas + '_dseg.nii.gz')

        if opi(atlasfile):

            nl = 'INFO: Working in sending ' + atlas + ' into the ' + targetname + 'space'
            run_cmd.msg(nl, diary_file, 'OKGREEN')

            norm2template.apply(sourcename,stdy_template,labels_dir,targetname,atlasfile,
                                mTx['invtransforms'],w2i,atlas,path_label_code,atlas,'','', diary_file,sing_wb)
        else:
            nl = 'WARNING: ' + atlasfile  + ' not found in list_atlases, we can continue but it might restrict several outcome of the script'
            run_cmd.msg(nl, diary_file, 'WARNING')


    segfile = opj(BASE_atlas_folder, species + '_seg-4FS_dseg.nii.gz')
    if opi(segfile):
        nl = 'INFO: Working in sending ' + segfile + ' in ' + targetname + ' space'
        run_cmd.msg(nl, diary_file, 'OKGREEN')
        print(mTx['invtransforms'])
        norm2template.apply(sourcename, stdy_template, labels_dir, targetname, segfile,
                                mTx['invtransforms'], w2i,'4FS',path_label_code, 'FreeSurfer', '','',
                                diary_file, sing_wb)
    else:
        nl = 'WARNING: ' + segfile + ' not found in list_atlases, we can continue but it might restrict several outcome of the script'
        run_cmd.msg(nl, diary_file, 'WARNING')

    if opi(BASE_mask):
        nl = 'INFO: Working in sending ' + BASE_mask + ' in ' + targetname + ' space'
        run_cmd.msg(nl, diary_file, 'OKGREEN')

        norm2template.apply(sourcename, stdy_template, masks_dir, targetname, BASE_mask,
                            mTx['invtransforms'], w2i,'','', '', '','',
                            diary_file, sing_wb)
    else:
        nl = 'WARNING: ' + segfile + ' not found in list_atlases, we can continue but it might restrict several outcome of the script'
        run_cmd.msg(nl, diary_file, 'WARNING')

