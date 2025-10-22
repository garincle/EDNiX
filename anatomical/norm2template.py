#import
import os
import json
import numpy as np
import ants

opj = os.path.join
ope = os.path.exists
opb = os.path.basename
opd = os.path.dirname
opi = os.path.isfile

from Tools import run_cmd

def norm(targetname,target_anatfile,masktarget,target_vol,
         sourcename,source_anatfile,masksource,
         target_transfo,type_of_transform,transfonameT,transfonameS,
         affmetricTrans,affmetric,interpol,
         diary_name,sing_wb,suffix,for_surf):

    nl = 'Normalization to ' + sourcename
    run_cmd.msg(nl, diary_name,'OKGREEN')

    if ope(target_transfo) == False:
        os.makedirs(target_transfo)

    TARGET = ants.image_read(target_anatfile)
    SOURCE = ants.image_read(source_anatfile)
    if masktarget== '':
        MASKTARGET = None
    else:
        MASKTARGET = ants.image_read(masktarget)
    if masksource == '':
        MASKSOURCE = None
    else:
        MASKSOURCE = ants.image_read(masksource)

    # always start with translation (very helpful if the two images are far from each other)
    mtx1 = ants.registration(fixed=SOURCE, moving=TARGET, type_of_transform='Translation',
                             aff_metric=affmetricTrans,
                             mask=MASKSOURCE,moving_mask=MASKTARGET,
                             outprefix=opj(target_transfo, transfonameT + '_'))

    translat_img = ants.apply_transforms(fixed=SOURCE, moving=TARGET,transformlist=mtx1['fwdtransforms'],
                                         interpolator=interpol)
    ants.image_write(translat_img, opj(target_transfo, transfonameT + '.nii.gz'), ri=False)

    dictionary = {"Sources": [target_anatfile,
                              source_anatfile],
                  "Description": 'Co-registration (translation,ANTspy).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(target_transfo, transfonameT + '.json'), "w") as outfile:
        outfile.write(json_object)


    if type_of_transform == 'SyNCC':
        mTx2 = ants.registration(fixed=SOURCE, moving=TARGET,
                                 type_of_transform=type_of_transform,
                                 initial_transform=mtx1['fwdtransforms'],
                                 mask=MASKSOURCE,moving_mask=MASKTARGET,
                                 outprefix=opj(target_transfo, transfonameS +'_'))

    elif type_of_transform == 'Rigid' or type_of_transform == 'Affine':
        mTx2 = ants.registration(fixed=SOURCE, moving=TARGET,initial_transform=mtx1['fwdtransforms'],
                                 mask=MASKSOURCE,moving_mask=MASKTARGET,
                                 type_of_transform=type_of_transform, aff_metric=affmetric,
                                 outprefix=opj(target_transfo, transfonameS +'_'))


    else:
        mTx2 = ants.registration(fixed=SOURCE, moving=TARGET,
                                 type_of_transform=type_of_transform,
                                 aff_metric=affmetric,
                                 initial_transform=mtx1['fwdtransforms'],
                                 mask=MASKSOURCE,moving_mask=MASKTARGET,
                                 outprefix=opj(target_transfo, transfonameS +'_'),
                                 grad_step=0.1,flow_sigma=3,total_sigma=0,aff_sampling=32,aff_random_sampling_rate=0.2,
                                 syn_sampling=32,
                                 aff_iterations=(1000, 500, 250, 100),aff_shrink_factors=(8, 4, 2, 1),
                                 aff_smoothing_sigmas=(3, 2, 1, 0),
                                 reg_iterations=(1000, 500, 250, 100),reg_smoothing_sigmas=(3, 2, 1, 0),
                                 reg_shrink_factors=(8, 4, 2, 1),
                                 verbose=True)

    Norm = ants.apply_transforms(fixed=SOURCE, moving=TARGET, transformlist=mTx2['fwdtransforms'],
                                 interpolator=interpol)
    ants.image_write(Norm, opj(target_vol, targetname + '_space-' + sourcename + suffix + '.nii.gz'), ri=False)

    dictionary = {"Sources": [target_anatfile,
                              source_anatfile],
                  "Description": 'Coregistration to the ' + sourcename + ' template.', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(target_vol, targetname + '_space-' + sourcename + suffix + '.json'), "w") as outfile:
        outfile.write(json_object)
    
    # for surface transformation
    if for_surf == 1:
        Affine = ants.read_transform(mTx2['fwdtransforms'][1])
        PTX = Affine.apply_to_point((0, 0, 0))
        test = [None] * 4, [None] * 4, [None] * 4, [None] * 4
    
        test[0][:-1] = Affine.parameters[0:3]
        test[0][2] = test[0][2] * -1
        test[0][3] = PTX[0] * -1
    
        test[1][:-1] = Affine.parameters[3:6]
        test[1][2] = test[1][2] * -1
        test[1][3] = PTX[1] * -1
    
        test[2][:-1] = Affine.parameters[6:9]
        test[2][0] = test[2][0] * -1
        test[2][1] = test[2][1] * -1
        test[2][3] = PTX[2]
    
        test[3][:] = 0, 0, 0, 1
    
        np.savetxt(opj(target_transfo, targetname + '_to_' + sourcename + '_for_surface.mat'), test, fmt='%10.8f', delimiter=' ')
    
        cmd = (sing_wb + 'wb_command -convert-affine -from-world ' + opj(target_transfo, targetname + '_to_' + sourcename + '_for_surface.mat') +
               ' -inverse -to-world ' + opj(target_transfo, 'affine_' + sourcename + '.nii.gz'))
        run_cmd.run(cmd,diary_name)
    
        cmd = (sing_wb + 'wb_command -convert-warpfield -from-itk ' + opj(target_transfo, transfonameS + '_1InverseWarp.nii.gz') +
               ' -to-world ' + opj(target_transfo, 'standard_' + sourcename + '.nii.gz'))
        run_cmd.run(cmd,diary_name)

    return mTx2


def myelin(targetname,sourcename,ref_template,labels_dir,atlas_label,transfolist,w2inv_inv,interpol,diary_name):

    nl = 'Normalization to the ' + sourcename + ' template'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    myelin = ants.image_read(opj(labels_dir, targetname + '_space-acpc_desc-template_T1wDividedByT2w.nii.gz'))
    SOURCE    = ants.image_read(ref_template)

    Norm = ants.apply_transforms(fixed=SOURCE, moving=myelin,
                                 transformlist=transfolist,
                                 interpolator=interpol,whichtoinvert=w2inv_inv)
    ants.image_write(Norm, opj(atlas_label, targetname + '_space-' + sourcename + '_T1wDividedByT2w.nii.gz'), ri=False)


def apply(template_name,template_file,destination,targetname,file,transfolist,w2inv_inv,descript,path_code,labelname, Timage,
          interpol,diary_name,sing_wb):

    nl = 'Normalization to the ' + template_name + ' template'
    run_cmd.msg(nl, diary_name,'OKGREEN')
    N = opb(file).split('.')[0].split('_')
    type = N[-1]

    if type == 'dseg':
        interp = 'genericLabel'
        name = targetname + '_space-' + template_name + '_seg-' + descript + '_dseg.nii.gz'
    elif type == 'mask':
        interp = 'genericLabel'
        name = targetname + '_space-' + template_name + descript + '_mask.nii.gz'
    elif 'T1' in type:
        interp = interpol
        type = Timage
        name = targetname + '_space-' + template_name + '_desc-' + descript + '_' + type + '.nii.gz'
    elif 'T2' in type:
        type = Timage
        interp = interpol
        name = targetname + '_space-' + template_name + '_desc-' + descript + '_' + type + '.nii.gz'
    elif type == 'bold':
        interp = interpol
    elif type == 'pet':
        interp = 'linear'
        name = targetname + '_space-' + template_name + '_desc-' + descript + '_T1w.nii.gz'
    print('file = '  + file)
    print('ope = ' + str(ope(file)))
    header = ants.image_header_info(file)
    if len(header['dimensions']) > 3:
        dim4 = 3
    else:
        dim4 = 0

    IMG    = ants.image_read(file)
    SOURCE = ants.image_read(template_file)

    A = ants.apply_transforms(fixed=SOURCE, moving=IMG,
                              transformlist=transfolist,
                              interpolator=interp,imagetype=dim4,whichtoinvert=w2inv_inv)
    os.makedirs(destination, exist_ok=True)
    ants.image_write(A, opj(destination, name), ri=False)

    if type == 'dseg':
        cmd = (sing_wb + 'wb_command -volume-label-import' + ' ' + opj(destination, name) +
        ' ' + opj(path_code,labelname + '_label.txt') + ' ' + opj(destination, name) + ' -drop-unused-labels')
        run_cmd.run(cmd,diary_name)
    dictionary = {"Sources": [file,
                              template_file],
          "Description": 'normalization to the ' + template_name + ' template.', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(destination, name).replace('.nii.gz','.json'), "w") as outfile:
        outfile.write(json_object)

    return opj(destination, name)

