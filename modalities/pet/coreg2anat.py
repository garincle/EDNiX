"""
Step 6 of the PET pipeline: registration of the PET session to the anatomy, the CT or a template
(scenarios of Illustration/pet4.svg), and atlases / masks / anatomy brought into the PET space.

The registration is estimated once per session, on the mean image of the session after the
between-runs motion correction (all the runs share its space), and stored as a registration record
(transfo.py) with the chains to every target space:

1v1  anatomy, same session (hybrid PET-MR): PET -> acpc = the AC-PC matrix of the anatomy
     (anat/matrices/acpc_0GenericAffine.mat), as for the functional images acquired in the space of
     the anatomy (anat_func_same_space): the PET was reoriented like the anatomy at step 1.
1v2a anatomy, other session: rigid registration PET -> anatomy (AC-PC), mutual information, restricted
     to the dilated brain of the anatomy (<session>_from-pet-to-<T>_0GenericAffine.mat).
1v2b anatomy and CT: PET -> CT (same space: nothing; otherwise rigid registration) then CT -> anatomy
     (record of the CT pipeline).
2v1  CT and CT template, no anatomy: PET -> CT -> template (CT pipeline).
2v2  CT, no CT template: PET -> CT -> CT study template (CT pipeline); no atlas.
3v1  PET template in the space of the reference template: PET -> PET template (affine + SyN).
3v2  no anatomy, no template: PET -> PET study template (step 2); no atlas.
With an anatomy, the AC-PC chain is composed with the anatomy -> template chains of
modalities/anat/loop3/backtonative.get (template, study template, reference session), the bridge
used by the functional pipeline.
"""
import os
import shutil

import ants

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd, QC_plot
from modalities.pet import petutils as pu
from modalities.pet import transfo
from modalities.anat.freesurfer import FS_freeview

labelinterp = 'genericLabel'
petinterp = 'linear'
anatinterp = 'linear'


def scenario(IhaveanANAT, anat_pet_same_space, IhaveaCT, ct_targets, pet_template, study_template):
    """Scenario code of pet4.svg from the launcher options (ct_targets: target names of the CT record)."""
    if IhaveanANAT:
        if anat_pet_same_space:
            return '1v1'
        return '1v2b' if IhaveaCT else '1v2a'
    if IhaveaCT:
        if ct_targets is None:
            raise FileNotFoundError('IhaveaCT=True but no CT registration record: run the CT pipeline first')
        return '2v2' if ct_targets == ['studyTemplate'] else '2v1'
    if pet_template:
        return '3v1'
    if study_template:
        return '3v2'
    raise ValueError('no anatomy, no CT, no PET template and no PET study template: nothing to register the PET '
                     'to (set IhaveanANAT, IhaveaCT, template_PET or creat_study_template_pet)')


def _dilated(mask_file, mm):
    m = ants.image_read(mask_file)
    r = max(1, int(round(mm / min(m.spacing[:3]))))
    return ants.iMath(ants.threshold_image(m, 0.5, 1e9, 1, 0), 'MD', r)


def rigid_to(fixed_file, moving_file, prefix, metric, mask_file=None, dilate_mm=5.0, initial=None):
    fixed = ants.image_read(fixed_file)
    moving = ants.image_read(moving_file)
    kw = {}
    if mask_file and opi(mask_file):
        kw['mask'] = _dilated(mask_file, dilate_mm)
    if initial:
        kw['initial_transform'] = initial
    ants.registration(fixed=fixed, moving=moving, type_of_transform='Rigid', aff_metric=metric,
                      outprefix=prefix, **kw)
    return transfo.registration(prefix)


def manual_correction(pet_mean, fixed_file, chain, prefix, pet_in_fixed, imgtype, diary_name, sing_fs):
    """
    PET engineer's manual correction: the PET in the anatomical space is opened in freeview with its
    'Transform Volume' tool; when the user saves it (same name), the registration is estimated again from
    the edited image, starting from the current one (PET -> PET, mean squares).
    """
    transfo.apply(fixed_file, pet_mean, chain, 'to', petinterp, out=pet_in_fixed)
    before = os.path.getmtime(pet_in_fixed)
    run_cmd.msg('MANUAL CORRECTION: move the PET (Tools > Transform Volume) and save it under the same name '
                '(File > Save Volume); close freeview to continue', diary_name, 'WARNING')
    FS_freeview.check_coreg(pet_in_fixed, fixed_file, imgtype, sing_fs, transform_volume=True)
    if os.path.getmtime(pet_in_fixed) == before:
        run_cmd.msg('INFO: the PET was not modified: registration kept', diary_name, 'OKGREEN')
        return chain
    edited = ants.image_read(pet_in_fixed)
    ants.registration(fixed=edited, moving=ants.image_read(pet_mean), type_of_transform='Rigid',
                      aff_metric='meansquares', initial_transform=[f for f, _ in chain['to']],
                      outprefix=prefix)
    run_cmd.msg('INFO: registration re-estimated from the manual correction', diary_name, 'OKGREEN')
    return transfo.registration(prefix)


def register(sc, session_mean, session, paths, anat, packages, ct_record, pet_template, reference,
             study_template, T_PET, ct_pet_same_space, coreg_metric, type_of_transform_template, do_manual,
             Need_check_coreg, bids_dir, diary_name, sing_fs):
    """
    Registration of the session (scenario `sc`); writes and returns the registration record.
    anat: dict of the anatomical files in AC-PC space (see _0_Pipeline_launcher.anat_files) or None.
    packages: output of backtonative.get (or None without anatomy).
    """
    nl = 'Registration of the PET session ' + session + ' (scenario ' + sc + ' of pet4.svg)'
    run_cmd.msg(nl, diary_name, 'HEADER')
    mat_dir = paths['orig_matrices']
    os.makedirs(mat_dir, exist_ok=True)
    targets = {}

    def pet_to_ct():
        """PET -> CT chain (empty when they share the space)."""
        ct_img = ct_record['reference']
        if ct_pet_same_space:
            return {'to': [], 'from': []}
        return rigid_to(ct_img, session_mean, opj(mat_dir, session + '_from-pet-to-ct_'), coreg_metric)

    if sc in ('1v1', '1v2a', '1v2b'):
        if sc == '1v1':
            acpc = transfo.check(transfo.affine(opj(anat['dir_transfo'], 'acpc_0GenericAffine.mat')))
        elif sc == '1v2a':
            prefix = opj(mat_dir, session + '_from-pet-to-' + T_PET + '_')
            acpc = rigid_to(anat['acpc'], session_mean, prefix, coreg_metric, mask_file=anat['mask'])
            if do_manual:
                acpc = manual_correction(session_mean, anat['acpc'], acpc, prefix,
                                         opj(paths['orig_formc'], 'pet2' + T_PET + '.nii.gz'), 'pet',
                                         diary_name, sing_fs)
        else:
            if 'acpc' not in ct_record['targets']:
                raise ValueError('IhaveanANAT and IhaveaCT: the CT record has no AC-PC target (run the CT pipeline '
                                 'with IhaveanANAT=True)')
            ct_acpc = ct_record['targets']['acpc']
            acpc = transfo.compose(pet_to_ct(), {'to': ct_acpc['to'], 'from': ct_acpc['from']})
        targets['acpc'] = dict(image=anat['acpc'], brain=anat['brain'], mask=anat['mask'], **acpc)
        for pkg in packages or []:
            name = pkg[0]
            chain = transfo.compose(acpc, transfo.from_backtonative(pkg))
            targets[name] = dict(image=pkg[1], brain=pkg[1], mask=pkg[2], **chain)
    elif sc in ('2v1', '2v2'):
        to_ct = pet_to_ct()
        for name, t in ct_record['targets'].items():
            chain = transfo.compose(to_ct, {'to': t['to'], 'from': t['from']})
            targets[name] = dict(image=t['image'], brain=t.get('brain', t['image']), mask=t.get('mask', ''), **chain)
    elif sc == '3v1':
        prefix = opj(mat_dir, session + '_to_' + reference + '_')
        ants.registration(fixed=ants.image_read(pet_template), moving=ants.image_read(session_mean),
                          type_of_transform=type_of_transform_template, aff_metric=coreg_metric, outprefix=prefix)
        targets[reference] = dict(image=pet_template, brain=pet_template,
                                  mask=pet_template.replace('.nii.gz', '_mask.nii.gz')
                                  if opi(pet_template.replace('.nii.gz', '_mask.nii.gz')) else '',
                                  **transfo.registration(prefix))
    elif sc == '3v2':
        prefix = opj(mat_dir, session + '_to_studyTemplate_')
        targets['studyTemplate'] = dict(image=study_template, brain=study_template, mask='',
                                        **transfo.registration(prefix))
    else:
        raise ValueError(f'unknown scenario {sc}')

    for name, t in targets.items():
        transfo.check(t)
    record = dict(modality='pet', scenario=sc, session=session, reference=session_mean, targets=targets,
                  created=pu.now())
    rec_file = record_path(paths, session)
    transfo.write_record(rec_file, record)

    # quality control: PET in every target space, brain mask of the target as contour
    os.makedirs(paths['qc'], exist_ok=True)
    for name, t in targets.items():
        pet_in = opj(paths['orig_formc'], session + '_space-' + name + '_desc-mean_pet.nii.gz')
        transfo.apply(t['image'], session_mean, t, 'to', petinterp, out=pet_in)
        if t.get('mask') and opi(t['mask']):
            png = opj(paths['qc'], session + '_space-' + name + '_pet.png')
            QC_plot.mosaic(pet_in, t['mask'], png)
            os.makedirs(opj(bids_dir, 'QC', 'pet_in_' + ('acpc' if name == 'acpc' else 'template')), exist_ok=True)
            shutil.copyfile(png, opj(bids_dir, 'QC', 'pet_in_' + ('acpc' if name == 'acpc' else 'template'), opb(png)))
        if Need_check_coreg:
            FS_freeview.check_coreg(pet_in, t['image'], 'pet', sing_fs, transform_volume=False)
    return rec_file


def record_path(paths, session):
    return opj(paths['orig_matrices'], session + '_desc-registration_pet.json')


def to_pet_space(record, session, paths, anat, template_files, list_atlas, diary_name, T_PET='T1w'):
    """
    Atlases, masks and anatomy (or template) in the PET space of the session (orig/labels, orig/masks,
    orig/preprocessing). Atlases: nearest-label interpolation (genericLabel); masks: linear interpolation
    then threshold 0.5 (fraction of the PET voxel covered); 4D atlases (levels) volume by volume.
    template_files: dict(name, atlas_dir, species, mask, gray, white, csf) used without anatomy.
    Returns dict(atlases={name: file}, masks={name: file}, anat=file).
    """
    nl = 'Atlases and masks into the PET space'
    run_cmd.msg(nl, diary_name, 'OKGREEN')
    for d in (paths['orig_labels'], paths['orig_masks'], paths['orig_process']):
        os.makedirs(d, exist_ok=True)
    ref = record['reference']
    out = dict(atlases={}, masks={}, anat='')
    if 'acpc' in record['targets'] and anat:
        t = record['targets']['acpc']
        atlas_src = {a: opj(anat['labels_dir'], anat['ID'] + '_seg-' + a + '_dseg.nii.gz') for a in list_atlas}
        mask_src = dict(brain=anat['mask'], gray=anat['gray'], white=anat['white'], csf=anat['csf'])
        anat_src, anat_out = anat['acpc'], pu.out_name(paths['orig_process'], session, space='orig', suffix=T_PET)
    elif template_files and template_files['name'] in record['targets'] and template_files.get('atlas_dir'):
        t = record['targets'][template_files['name']]
        atlas_src = {a: opj(template_files['atlas_dir'], template_files['species'] + '_seg-' + a + '_dseg.nii.gz')
                     for a in list_atlas}
        mask_src = {k: template_files.get(k, '') for k in ('brain', 'gray', 'white', 'csf')}
        anat_src = t['image']
        anat_out = pu.out_name(paths['orig_process'], session, space='orig', desc='template', suffix=T_PET)
    else:
        run_cmd.msg('WARNING: no atlas for this scenario (' + record['scenario'] + '): the regional analysis, the '
                    'brain ratio and the reference-region models need an atlas or masks', diary_name, 'WARNING')
        return out

    for name, src in atlas_src.items():
        if not opi(src):
            run_cmd.msg('WARNING: atlas not found: ' + src, diary_name, 'WARNING')
            continue
        dst = opj(paths['orig_labels'], session + '_space-orig_seg-' + name + '_dseg.nii.gz')
        transfo.apply(ref, src, t, 'from', labelinterp, out=dst)
        pu.write_json(dst, [src, ref], 'Atlas in the PET space (genericLabel interpolation).')
        out['atlases'][name] = dst
    for name, src in mask_src.items():
        if not src or not opi(src):
            run_cmd.msg(f'WARNING: {name} mask not found: {src}', diary_name, 'WARNING')
            continue
        dst = opj(paths['orig_masks'], session + '_space-orig' + ('' if name == 'brain' else '_desc-' + name)
                  + '_mask.nii.gz')
        m = transfo.apply(ref, src, t, 'from', 'linear')
        ants.image_write(ants.threshold_image(m, 0.5, 1e9, 1, 0), dst, ri=False)
        pu.write_json(dst, [src, ref], 'Mask in the PET space (linear interpolation, threshold 0.5).')
        out['masks'][name] = dst
    if anat_src and opi(anat_src):
        transfo.apply(ref, anat_src, t, 'from', anatinterp, out=anat_out)
        pu.write_json(anat_out, [anat_src, ref], 'Anatomical image in the PET space (visualisation).')
        out['anat'] = anat_out
    return out
