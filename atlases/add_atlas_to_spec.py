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
from Tools import getpath
from anat import transfo_T2toT1w
from anat.connectomeWB import WB_label


def add(cmd_mris,atlas,template,doT,label,name,LR,diary_name):
    '''
    add_atlas_to_spec takes a labelled volume file fr oma an atlas and create a new volume file readable for WB and the corresponding labelled surfaces
    Parameters
    ----------
    file : the labelled volume
    source : the T1/T2w volume sharing the same space as the labelled volume
    template : the destination of the new atlas
    name
    T_transfo: do you need to change the T1 into a T2
    doT :[yesorno,source,existornot,T_transfo,msk,type] type could 'Translation','Rigid','SyNCC'
    label the name of the atlas
    diary_name : where to write the processings steps
    Returns
    -------
    '''

    nl = 'add atlas to specfile'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    main_path = opd(opd(opd(template)))
    dir_path  = opd(template)
    animal    = opb(template).split('_')[0]

    check = glob.glob(opj(main_path,'preprocessing','*_resamp-*'))
    if check == []:
        change_hd = 0
    else:
        change_hd = 1

    FS_dir   = opj(main_path,'freesurfer')
    WB_dir   = opj(main_path,'native','surfaces','Native_resol')
    surftype = '.native'
    roi = 'roi'
    spec     = opj(main_path,'native', 'surfaces', 'Native_resol', animal + '_native_LR.wb.spec')

    if LR ==0:
        suffix = '_StatsLUT.txt'
    else:
        suffix = '_l.ctab'


    if doT[0] == 1:

        if doT[2] == 0:
            source = doT[1]
            if not ope(opj(main_path, 'matrices')):
                os.makedirs(opj(main_path, 'matrices'))
            if doT[3] == 1:
                transfo_T2toT1w.prepa(source, doT[4], '', diary_name)
                source = source.replace('_T2w.nii.gz','_T1w.nii.gz')

            orig = ants.image_read(source)
            new  = ants.image_read(template)
            ants.registration(fixed=orig, moving=new,
                                type_of_transform=doT[5],
                                initial_transform=None,
                                outprefix=opj(main_path, 'matrices', 'from_' + name + '_'))

            if type == 'SyNCC':
                transfo = [opj(main_path, 'matrices', 'from_' + name + '_1InverseWarp.nii.gz'),
                           opj(main_path, 'matrices', 'from_' + name + '_0GenericAffine.mat')]
                w2i     = [False,True]
            else:
                transfo = opj(main_path, 'matrices', 'from_' + name + '_0GenericAffine.mat')
                w2i     = [True]

            WB_label.vol(atlas, template, spec, transfo, w2i, label + '_label.txt', diary_name)

        elif doT[2] == 1 and name == 'to_ac':
            transfo = opj(main_path, 'matrices', 'to_ac_0GenericAffine.mat')
            w2i = [True]
            WB_label.vol(atlas, template, spec, transfo, w2i, label + '_label.txt', diary_name)

        else :
            if opi(opj(main_path, 'matrices', 'from_' + name + '_0GenericAffine.mat'))==True:
                if opi(opj(main_path, 'matrices', 'from_' + name + '_1InverseWarp.nii.gz')):
                    transfo = [opj(main_path, 'matrices', 'from_' + name + '_1InverseWarp.nii.gz'),
                               opj(main_path, 'matrices', 'from_' + name + '_0GenericAffine.mat')]
                    w2i = [False,True]

                    #transfo = [opj(main_path, 'matrices', 'from_' + name + '_0GenericAffine.mat'),
                    #           opj(main_path, 'matrices', 'from_' + name + '_1Warp.nii.gz')]
                    #w2i = [False, False]
                else :
                    transfo = opj(main_path, 'matrices', 'from_' + name + '_0GenericAffine.mat')
                    w2i     = [True]

            WB_label.vol(atlas, template, spec, transfo, w2i, label + '_label.txt', diary_name)

    else :
        WB_label.vol(atlas, template, spec, '', '', label + '_label.txt', diary_name)

    WB_label.surfFS(opj(dir_path, 'labels', animal + '_seg-' + opb(label).split('_')[0] + '_dseg.nii.gz'),
                    label + suffix,
                    change_hd, diary_name)

    WB_label.surfWB(cmd_mris, FS_dir, animal, opb(label).split('_')[0], WB_dir, surftype, diary_name)
    WB_label.cifti(animal, opb(label).split('_')[0], WB_dir, surftype, roi, spec, diary_name)









