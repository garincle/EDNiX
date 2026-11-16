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


    if type_of_transform == 'SyNadaptive':
        # SyNCC, refined with a compliant masked SyN only when the overlap checks fail (see registration_adaptive.py)
        from modalities.anat import registration_adaptive
        mTx2 = registration_adaptive.register(SOURCE, TARGET, mtx1['fwdtransforms'],
                                              opj(target_transfo, transfonameS + '_'),
                                              fixed_mask=MASKSOURCE, moving_mask=MASKTARGET,
                                              log=lambda m: run_cmd.msg(m, diary_name, 'OKGREEN'))

    elif type_of_transform == 'SyNCC':
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
                                 verbose=False)

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


# acpc rigid alignment (loop1). antsRegistration with MI samples 20 % of the voxels at random (no fixed seed in
# ANTsPy): on some sessions a run lands in a wrong local minimum that pushes the brain out of the template grid
# (Rat sub-301909: one third of the brain kept), while a re-run works. The default run is kept when it passes the
# check below; otherwise other starting points / seeds are tried and the best one is kept.
ACPC_MIN_KEPT = 0.90          # fraction of the brain (in mm3) that must stay inside the template grid
ACPC_SEEDS = (1, 2, 3)


def _brain_kept_and_dice(fixed, moving, tx):
    """fraction of the moving brain (voxels > 0, in mm3) kept inside the fixed grid after tx, and its Dice with the
    fixed brain (voxels > 0)"""
    mb = ants.threshold_image(moving, 1e-6, 1e30)
    w = ants.apply_transforms(fixed=fixed, moving=mb, transformlist=tx, interpolator='nearestNeighbor').numpy() > 0.5
    fb = fixed.numpy() > 0
    kept = w.sum() * np.prod(fixed.spacing) / max(float((mb.numpy() > 0.5).sum() * np.prod(moving.spacing)), 1e-9)
    dice = 2 * (w & fb).sum() / max(float(w.sum() + fb.sum()), 1e-9)
    return float(kept), float(dice)


def _acpc_candidate(kind, seed, fixed, moving, moving_mask, tmp_prefix, type_of_transform, affmetricT, affmetric):
    """kind 'default' = translation then rigid (the historical EDNiX alignment);
    'rotsearch' = ants.affine_initializer (search over rotations from the image centres) then rigid"""
    os.environ['ANTS_RANDOM_SEED'] = str(seed)
    if kind == 'default':
        t1 = ants.registration(fixed=fixed, moving=moving, type_of_transform='Translation', aff_metric=affmetricT,
                               moving_mask=moving_mask, outprefix=tmp_prefix + 'T_')
        init = t1['fwdtransforms']
    else:
        init = [ants.affine_initializer(fixed, moving, search_factor=20, radian_fraction=0.15,
                                        use_principal_axis=False, local_search_iterations=10,
                                        txfn=tmp_prefix + 'init.mat')]
    r = ants.registration(fixed=fixed, moving=moving, initial_transform=init, moving_mask=moving_mask,
                          type_of_transform=type_of_transform, aff_metric=affmetric, outprefix=tmp_prefix + 'R_')
    return r['fwdtransforms']


def acpc(targetname, target_anatfile, masktarget, target_vol, source_anatfile, target_transfo,
         type_of_transform, affmetricTrans, affmetric, interpol, diary_name, suffix):
    """rigid alignment of the skull-stripped anat to the template (step acpc of loop1), with a check that the brain
    stays in the template grid. Writes, as norm(): <target_transfo>/acpc_0GenericAffine.mat and
    <target_vol>/<targetname>_space-acpc<suffix>.nii.gz"""
    import glob, shutil
    os.makedirs(target_transfo, exist_ok=True)
    fixed = ants.image_read(source_anatfile)
    moving = ants.image_read(target_anatfile)
    mmask = ants.image_read(masktarget) if masktarget else None
    seed0 = os.environ.get('ANTS_RANDOM_SEED')
    tried = []
    try:
        plan = [('default', ACPC_SEEDS[0])] + \
               [('rotsearch', ACPC_SEEDS[0])] + [(k, s) for s in ACPC_SEEDS[1:] for k in ('default', 'rotsearch')]
        for n, (kind, seed) in enumerate(plan):
            prefix = opj(target_transfo, 'acpc_try' + str(n) + '_')
            tx = _acpc_candidate(kind, seed, fixed, moving, mmask, prefix, type_of_transform, affmetricTrans, affmetric)
            kept, dice = _brain_kept_and_dice(fixed, moving, tx)
            tried.append(dict(kind=kind, seed=seed, kept=kept, dice=dice, tx=tx[0]))
            run_cmd.msg('INFO: acpc ' + kind + ' seed ' + str(seed) + ': brain kept in the template grid ' +
                        str(round(100 * kept, 1)) + ' %, Dice with the template ' + str(round(dice, 3)),
                        diary_name, 'OKGREEN')
            if n == 0 and kept >= ACPC_MIN_KEPT:
                break
            if n == 1:
                best = max(tried, key=lambda t: t['dice'])
                if best['kept'] >= ACPC_MIN_KEPT:
                    break
        best = tried[0] if len(tried) == 1 else max(tried, key=lambda t: t['dice'])
        if len(tried) > 1:
            run_cmd.msg('WARNING: acpc of ' + targetname + ': the default alignment kept only ' +
                        str(round(100 * tried[0]['kept'], 1)) + ' % of the brain in the template grid; kept ' +
                        best['kind'] + ' seed ' + str(best['seed']) + ' (' + str(round(100 * best['kept'], 1)) +
                        ' %, Dice ' + str(round(best['dice'], 3)) + ')', diary_name, 'WARNING')
        if best['kept'] < ACPC_MIN_KEPT:
            run_cmd.msg('WARNING: acpc of ' + targetname + ': no alignment keeps ' + str(int(100 * ACPC_MIN_KEPT)) +
                        ' % of the brain in the template grid: check the skull stripping or put a manual ' +
                        targetname + '_acpc_tmp.nii.gz in preprocessing', diary_name, 'WARNING')
        out_mat = opj(target_transfo, 'acpc_0GenericAffine.mat')
        shutil.copyfile(best['tx'], out_mat)
    finally:
        if seed0 is None:
            os.environ.pop('ANTS_RANDOM_SEED', None)
        else:
            os.environ['ANTS_RANDOM_SEED'] = seed0
        for f in glob.glob(opj(target_transfo, 'acpc_try*')):
            os.remove(f)

    Norm = ants.apply_transforms(fixed=fixed, moving=moving, transformlist=[out_mat], interpolator=interpol)
    out_img = opj(target_vol, targetname + '_space-acpc' + suffix + '.nii.gz')
    ants.image_write(Norm, out_img, ri=False)
    with open(out_img.replace('.nii.gz', '.json'), 'w') as f:
        json.dump({"Sources": [target_anatfile, source_anatfile],
                   "Description": 'Rigid alignment to the template (acpc), ' + best['kind'] + ' seed ' + str(best['seed']) +
                                  ', brain kept in grid ' + str(round(100 * best['kept'], 1)) + ' %.'}, f, indent=2)
    return [out_mat]


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
    print('type = ' + type)
    print('interp = ' + interp)
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
    print(A)
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

