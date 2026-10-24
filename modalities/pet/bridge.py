"""
Group step 'bridge' of the PET pipeline: link of the CT study template to the reference (MRI) template through a PET
study template (scenario 2v2 of Illustration/pet4.svg: CT, no CT template, no anatomy).

Without an anatomy, the PET of every session reaches the CT study template through the CT (PET -> CT -> CT study
template, record of the CT pipeline), but the CT study template is not linked to the reference template: no atlas,
no brain mask. A CT has little contrast inside the brain, so the CT template is a poor target for the T1 template.
The FDG images have it (grey / white matter):

1. PET study template in the space of the CT study template: mean of the session means brought there through
   their own chains (each scaled by its 99th percentile, so that the doses do not weigh differently). The mean, not
   ants.build_template: the template must be in the space the sessions are registered to; the warps of a template
   built again from these images would not be part of the session chains.
2. Reference template (skull-stripped T1, moving) -> PET study template (fixed): multi-start similarity
   (coreg2anat.rigid_to) then affine, brain mask of the reference template as moving mask (metric in the brain of the
   T1); optional SyN stage (bridge_transform='SyNRA') in the bright (FDG) mask of the PET template.
3. Brain mask of the reference template in the PET study template, ratio image (template / its brain mean), PET study
   template in the reference space (QC).
4. Every 2v2 session: its record receives the reference target (session -> CT study template -> reference) and the
   atlases / masks of the reference template are brought into the PET space (coreg2anat.to_pet_space).

Same chain as the dog notebook of the PET engineer (CT study template, PET through the CT warps, PET study template,
T1 template affinely on the PET template, brain mask brought back).
"""
import os
import shutil

import numpy as np
import ants

opj = os.path.join
opb = os.path.basename
opi = os.path.isfile

from Tools import run_cmd, QC_plot
from modalities.pet import petutils as pu
from modalities.pet import transfo, coreg2anat

BRIDGE_TRANSFORMS = ('Affine', 'SyNRA')


def record_path(study_dir, reference):
    return opj(study_dir, 'matrices', 'studyTemplate_to_' + reference + '_desc-bridge_pet.json')


def template_in_ct_space(records, out, diary_file):
    """Mean of the session means (scaled by their 99th percentile) in the space of the CT study template."""
    imgs, ref, sources = [], None, []
    for rec in records:
        t = rec['targets']['studyTemplate']
        ref = ref or ants.image_read(t['image'])
        mv = transfo.apply(ref, rec['reference'], t, 'to', 'linear').numpy()
        p = np.percentile(mv[np.isfinite(mv) & (mv > 0)], 99) if np.any(mv > 0) else 1.0
        imgs.append(mv / float(p if p > 0 else 1.0))
        sources.append(rec['reference'])
    if len(imgs) < 2:
        raise ValueError(f'the PET study template needs at least 2 sessions registered to the CT study template; '
                         f'found {len(imgs)}')
    ants.image_write(ref.new_image_like(np.mean(imgs, axis=0).astype(np.float32)), out, ri=False)
    pu.write_json(out, sources, f'PET study template: mean of {len(imgs)} session means in the space of the CT study '
                                f'template (through their PET -> CT -> CT study template chains), each scaled by its '
                                f'99th percentile.')
    run_cmd.msg('INFO: PET study template in the CT study template space: ' + opb(out), diary_file, 'OKGREEN')
    return out


def register_reference(pet_template, ref_brain, ref_mask, study_dir, reference, metric, bridge_transform, diary_file,
                       percentile=50.0):
    """Reference template (skull-stripped T1) -> PET study template; returns the chains studyTemplate <-> reference
    ({'to': resample a study-template image in the reference space, 'from': the reverse})."""
    if bridge_transform not in BRIDGE_TRANSFORMS:
        raise ValueError(f'bridge_transform={bridge_transform!r}: expected one of {BRIDGE_TRANSFORMS}')
    mat = opj(study_dir, 'matrices')
    os.makedirs(mat, exist_ok=True)
    fixed_mask = coreg2anat.threshold_mask(pet_template, percentile, 3.0)
    ants.image_write(fixed_mask, opj(mat, 'studyTemplate_desc-coreg_mask.nii.gz'), ri=False)
    moving_mask = coreg2anat._dilated(ref_mask, 2.0)
    # similarity (multi-start), then affine from it: an affine started from a translation stops in a sheared local
    # optimum (6 mm error on the phantom against 1.2 mm with the cascade, as the Translation / Rigid / Similarity /
    # Affine stages of ANTs' TRSAA). The skull-stripped T1 covers only the brain of the head-wide PET template: no
    # fixed mask, no coverage condition (the moving mask restricts the metric to the brain of the T1).
    sim = coreg2anat.rigid_to(pet_template, ref_brain, opj(mat, reference + '_to_studyTemplate_desc-similarity_'),
                              metric, coreg_transform='Similarity', moving_mask=moving_mask, min_coverage=0.0,
                              diary_name=diary_file, same_subject=False)
    prefix = opj(mat, reference + '_to_studyTemplate_')
    chain = transfo.register(pet_template, ref_brain, prefix, 'Affine', aff_metric=metric,
                             initial_transform=[f for f, _ in sim['to']], moving_mask=moving_mask,
                             **coreg2anat._adaptive_reg_params(ants.image_read(pet_template),
                                                               ants.image_read(ref_brain)))
    if bridge_transform == 'SyNRA':
        nl = opj(mat, reference + '_to_studyTemplate_desc-SyN_')
        chain = transfo.register(pet_template, ref_brain, nl, 'SyNOnly', aff_metric=metric, syn_metric='mattes',
                                 initial_transform=[f for f, _ in chain['to']], mask=coreg2anat._dilated(fixed_mask, 3))
    # chain: M = reference, T = study template; swapped: M = study template, T = reference
    return {'to': chain['from'], 'from': chain['to']}


def run(records, study_dir, tracer, reference, ref_brain, ref_mask, metric, bridge_transform, diary_file):
    """Steps 1-3 of the module docstring; writes and returns the bridge record."""
    nl = 'Bridge of the CT study template to the ' + reference + ' template (PET study template)'
    run_cmd.msg(nl, diary_file, 'HEADER')
    for f in (ref_brain, ref_mask):
        if not opi(f):
            raise FileNotFoundError(f'bridge: reference template file not found: {f}')
    tag = ('trc-' + tracer + '_') if tracer else ''
    pet_tpl = template_in_ct_space(records, opj(study_dir, 'studyTemplate_' + tag + 'space-CTstudyTemplate_pet.nii.gz'),
                                   diary_file)
    chain = register_reference(pet_tpl, ref_brain, ref_mask, study_dir, reference, metric, bridge_transform,
                               diary_file)
    mask = pet_tpl.replace('_pet.nii.gz', '_mask.nii.gz')
    m = transfo.apply(pet_tpl, ref_mask, chain, 'from', 'linear')
    ants.image_write(ants.threshold_image(m, 0.5, 1e30, 1, 0), mask, ri=False)
    pu.write_json(mask, [ref_mask, pet_tpl], f'Brain mask of the {reference} template in the PET study template '
                                             f'(linear interpolation, threshold 0.5).')
    tpl, msk = ants.image_read(pet_tpl), ants.image_read(mask)
    if float(msk.sum()) == 0:
        raise ValueError('bridge: the brain mask of the reference template does not overlap the PET study template')
    ratio = pet_tpl.replace('_pet.nii.gz', '_desc-bratio_pet.nii.gz')
    mean = float(tpl.numpy()[msk.numpy() > 0].mean())
    ants.image_write((tpl * msk) / (mean if mean > 0 else 1.0), ratio, ri=False)
    pu.write_json(ratio, [pet_tpl, mask], 'PET study template / its mean in the brain mask.')
    in_ref = opj(study_dir, 'studyTemplate_' + tag + 'space-' + reference + '_pet.nii.gz')
    transfo.apply(ref_brain, pet_tpl, chain, 'to', 'linear', out=in_ref)
    pu.write_json(in_ref, [pet_tpl, ref_brain], f'PET study template in the {reference} space.')
    png = opj(study_dir, 'studyTemplate_' + tag + 'space-' + reference + '_pet.png')
    QC_plot.mosaic(in_ref, ref_mask, png)
    rec = dict(modality='pet', kind='bridge', reference=reference, image=pet_tpl, mask=mask,
               target=dict(image=ref_brain, brain=ref_brain, mask=ref_mask, **chain),
               bridge_transform=bridge_transform, sessions=[r['session'] for r in records], created=pu.now())
    rf = record_path(study_dir, reference)
    transfo.write_record(rf, rec)
    run_cmd.msg('INFO: bridge record: ' + rf + ' (QC: ' + opb(png) + ')', diary_file, 'OKGREEN')
    return rf


def extend_session(rec_file, bridge_file, paths, template_files, list_atlas, bids_dir, diary_file):
    """Step 4: reference target and atlases / masks in the PET space for one 2v2 session; returns the record."""
    rec = transfo.read_record(rec_file)
    br = transfo.read_record(bridge_file)
    name = br['reference']
    if 'studyTemplate' not in rec['targets']:
        raise ValueError(f'{rec_file}: no studyTemplate target (scenario {rec.get("scenario")}): not a 2v2 session')
    st = rec['targets']['studyTemplate']
    chain = transfo.compose({'to': st['to'], 'from': st['from']}, {'to': br['target']['to'], 'from': br['target']['from']})
    rec['targets'][name] = dict(image=br['target']['image'], brain=br['target']['brain'], mask=br['target']['mask'],
                                **transfo.check(chain))
    rec['bridge'] = bridge_file
    transfo.write_record(rec_file, rec)
    rec['pet_space'] = coreg2anat.to_pet_space(rec, rec['session'], paths, None, template_files, list_atlas,
                                               diary_file)
    transfo.write_record(rec_file, rec)
    # QC: PET of the session in the reference space
    os.makedirs(paths['qc'], exist_ok=True)
    t = rec['targets'][name]
    pet_in = opj(paths['orig_formc'], rec['session'] + '_space-' + name + '_desc-mean_pet.nii.gz')
    transfo.apply(t['image'], rec['reference'], t, 'to', 'linear', out=pet_in)
    png = opj(paths['qc'], rec['session'] + '_space-' + name + '_pet.png')
    QC_plot.mosaic(pet_in, t['mask'], png)
    qc_dir = opj(bids_dir, 'QC', 'pet_in_template')
    os.makedirs(qc_dir, exist_ok=True)
    shutil.copyfile(png, opj(qc_dir, opb(png)))
    run_cmd.msg('INFO: ' + rec['session'] + ': ' + name + ' target added through the PET study template',
                diary_file, 'OKGREEN')
    return rec

