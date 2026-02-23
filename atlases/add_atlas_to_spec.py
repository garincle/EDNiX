#import
import os
import glob
import ants



opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd
from anat.connectomeWB import WB_label
from atlases import templatefeat

T1interp    = 'hammingWindowedSinc'
labelinterp = 'genericLabel'

antsaff     = '0GenericAffine.mat'
antswarp    = '1Warp.nii.gz'
antsinvwarp = '1InverseWarp.nii.gz'


def add(specie,animal,newAtlasfile,newAtlasAnatfile,labelname,sameAnimal,doT,wFS,LR,diary_name,path_ATLAS, fs_tools,cmd_mris,sing_afni,sing_fs,export_fs,sing_wb):

    '''
    add_atlas_to_spec takes a labelled volume file fr oma an atlas and create a new volume file readable for WB and the corresponding labelled surfaces
    Parameters
    ----------
    file : the labelled volume
    source : the T1/T2w volume sharing the same space as the labelled volume
    template : the destination of the new atlas
    name
    T_transfo: do you need to change the T1 into a T2
    doT :[yes-or-no,source,exist-or-not,T_transfo,msk,type] type could 'Translation','Rigid','SyNCC'
    label the name of the atlas
    diary_name : where to write the processings steps
    Returns
    -------
    '''

    nl = 'add atlas to specfile'
    run_cmd.msg(nl, diary_name)

    [refDir,refile ] = templatefeat.get(specie, path_ATLAS, fs_tools,
                                         animal, '', '', 'wb_atlas', 'anat')


    # Targetparams:

    matrice_dir   = opj(refDir,'matrices')
    refFS_dir     = opj(refDir, 'freesurfer')
    reflabeldir = opj(refDir,'volumes','labels')
    check         = glob.glob(opj(refDir, 'preprocessing', '*_resamp-*'))
    refWB_dir     = opj(refDir,'surfaces','Native_resol')
    surftype = '.native'
    roi      = 'roi'
    spec = opj(refWB_dir, animal + '_native_LR.wb.spec')


    Tmode         = 'SyNCC' # use 'SyN' with syn_metric='mattes' for T2/T1

    #newAtlaslabel   is the file readable by workbench' :blabla_label.txt'
    #labelname : the name you want to give to the new atlas ex 'BAL2026'


    newAtlas2ref  = '-'.join(['from',labelname])
    newAtlas2refT = {'Tname' :opj(matrice_dir, '_'.join([newAtlas2ref,''])),
                     'affine': opj(matrice_dir, '_'.join([newAtlas2ref,antsaff])),
                     'warp'  :   opj(matrice_dir, '_'.join([newAtlas2ref,antswarp])),
                     'invwarp':opj(matrice_dir, '_'.join([newAtlas2ref,antsinvwarp])),
                     'type':Tmode}

    if doT == 1:
        if not opi(newAtlas2refT['affine']):
            if not ope(matrice_dir):
                os.makedirs(matrice_dir)

            source = ants.image_read(newAtlasAnatfile)
            target = ants.image_read(refile)

            if not sameAnimal:

                ants.registration(fixed=source, moving=target,
                                  type_of_transform='Translation',outprefix=newAtlas2refT['Tname'])
                ants.registration(fixed=source, moving=target,
                                    type_of_transform=Tmode,
                                    initial_transform=newAtlas2refT['affine'],
                                    outprefix=newAtlas2refT['Tname'])

                if opi(newAtlas2refT['invwarp']):
                    transfo = [newAtlas2refT['affine'],newAtlas2refT['invwarp']]
                    w2i     = [True,False]
                else:
                    transfo = newAtlas2refT['affine']
                    w2i     = [True]

                WB_label.vol(newAtlasfile, refile, reflabeldir,spec, transfo, w2i, opj(path_code_label,labelname + '_label.txt'), diary_name,sing_wb)

            else:
                ants.registration(fixed=source, moving=target,
                                  type_of_transform='Translation', outprefix=newAtlas2refT['Tname'])
                ants.registration(fixed=source, moving=target,
                                  type_of_transform='Rigid',
                                  initial_transform=newAtlas2refT['affine'],
                                  outprefix=newAtlas2refT['Tname'])

                transfo = newAtlas2refT['affine']
                w2i = [True]
                WB_label.vol(newAtlasfile, refile, reflabeldir,spec, transfo, w2i, opj(path_code_label,labelname + '_label.txt'), diary_name,sing_wb)

        else :
            if opi(newAtlas2refT['invwarp']):
                transfo = [newAtlas2refT['affine'],newAtlas2refT['invwarp']]
                w2i     = [True,False]
            else:
                transfo = newAtlas2refT['affine']
                w2i     = [True]

        WB_label.vol(newAtlasfile, refile, reflabeldir,spec, transfo, w2i, opj(path_code_label,labelname + '_label.txt'), diary_name,sing_wb)

    if wFS == 1:
        if not check:
            change_hd = 0
        else:
            change_hd = 1

        if LR == 0:
            suffix = 'StatsLUT.txt'
        else:
            suffix = 'l.ctab'
        WB_label.surfFS(opj(reflabeldir, animal + '_seg-' + labelname + '_dseg.nii.gz'),
                        opj(path_code_label,'_'.join([labelname,suffix])),
                        change_hd, diary_name,sing_fs,path_code_label,export_fs,sing_afni,sing_wb)

        WB_label.surfWB(cmd_mris , refFS_dir, animal, labelname, refWB_dir, surftype, diary_name)
    else:
        for h in [0,1]:
            WB_label.surfvolWB(newAtlasfile, animal, labelname, refWB_dir, surftype, h, diary_name, sing_wb)

    WB_label.cifti(animal, labelname, refWB_dir, surftype,roi, spec, diary_name)













