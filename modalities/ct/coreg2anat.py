"""
Step 3 of the CT pipeline: registration of the CT (scenarios 1v2b, 2v1, 2v2 of Illustration/pet4.svg) and
registration record read by the PET pipeline (same format as the PET record, modalities/pet/transfo.py).

with an anatomy (1v2b)     rigid CT -> anatomy (AC-PC head image), mutual information on the CT clipped to
                           [-1000, 2000] HU, restricted to the brain mask dilated by 10 mm (skull included);
                           anat_ct_same_space=True: the AC-PC matrix of the anatomy (same scanner space).
CT template (2v1)          affine + SyN CT -> CT template given in the space of the reference template.
CT study template (2v2)    registration to the study template of the CT images (step 2).
The PET engineer registered the CT to the raw anatomy and composed the result with 'shift' and 'acpc'
matrices (the 'shift' matrix is already contained in the AC-PC matrix of the anatomical pipeline).
"""
import os
import shutil

import numpy as np
import ants

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
opi = os.path.isfile

from Tools import run_cmd, QC_plot
from modalities.pet import petutils as pu
from modalities.pet import transfo
from modalities.pet.coreg2anat import (manual_correction, qc_contour, rigid_to, ct_head_mask, refine_in_brain,
                                       ct_intracranial_mask, _validate_coreg_transform)
from modalities.anat.freesurfer import FS_freeview

CT_CLIP = (-1000.0, 2000.0)


def clipped(ct_file, out):
    img = ants.image_read(ct_file)
    arr = np.clip(img.numpy(), *CT_CLIP)
    ants.image_write(ants.from_numpy(arr, origin=img.origin, spacing=img.spacing, direction=img.direction), out,
                     ri=False)
    return out


CT_BRAIN_WINDOW = (0.0, 80.0)
CT_BONE_HU = 200.0


def brain_window(ct_file, out, window=CT_BRAIN_WINDOW, smooth_mm=1.0):
    """
    CT in the brain window for the second pass of the CT -> anatomy registration (ct_brain_refine): Gaussian
    smoothing of sigma smooth_mm (noise of a low-dose attenuation CT), then 0 - 1 over `window` (HU: CSF ~5, white
    matter ~30, grey matter ~40), 0 below it (air, fat), 1 up to CT_BONE_HU (partial volume) and 0 above it: the
    bone is dark, as the cortical bone of a T1w. On [-1000, 2000] HU the whole parenchyma falls in one or two of
    the 32 bins of the Mattes mutual information and only the skull drives the registration; in the brain window
    the grey / white matter and the ventricles do (brain window of head CT: Muschelli et al. 2015).
    """
    ct = ants.image_read(ct_file)
    sm = ants.smooth_image(ct, smooth_mm, sigma_in_physical_coordinates=True) if smooth_mm else ct
    v = sm.numpy()
    lo, hi = window
    w = np.clip((v - lo) / (hi - lo), 0.0, 1.0)
    w[v > CT_BONE_HU] = 0.0
    ants.image_write(ct.new_image_like(w.astype(np.float32)), out, ri=False)
    pu.write_json(out, [ct_file], f'CT in the brain window {window} HU (smoothed {smooth_mm} mm, 0 above '
                                  f'{CT_BONE_HU} HU): moving image of the second CT -> anatomy pass.')
    return out


def head_image(ct_file, out, pad_mm=5.0):
    """CT restricted to the head for the CT study template (the manual crop of the PET engineer's notebook, made
    automatic): voxels outside the head mask (ct_head_mask: > -500 HU, largest component, holes filled, dilated by
    pad_mm) set to -1000 HU (air), image cropped to the box of the mask. The bed and the head holder no longer drive
    ants.build_template; world coordinates are kept, so the transforms apply to the full CT."""
    ct = ants.image_read(ct_file)
    head = ants.iMath(ct_head_mask(ct), 'MD', max(1, int(round(pad_mm / min(ct.spacing[:3])))))
    arr = np.where(head.numpy() > 0, np.clip(ct.numpy(), *CT_CLIP), CT_CLIP[0]).astype(np.float32)
    img = ants.crop_image(ct.new_image_like(arr), label_image=head, label=1)
    ants.image_write(img, out, ri=False)
    pu.write_json(out, [ct_file], f'CT restricted to the head (> -500 HU, largest component, dilated {pad_mm} mm; '
                                  f'outside -1000 HU), cropped: input of the CT study template.')
    return out


def record_path(paths, base):
    return opj(paths['orig_matrices'], base + '_desc-registration_ct.json')


def register(sc, ct_img, raw, base, paths, anat, packages, template_CT, reference, study_template, T_CT, metric,
             type_of_transform_template, do_manual, Need_check_coreg, bids_dir, diary_name, sing_fs,
             template_mask='', coreg_transform='Rigid', coreg_init='auto', ct_brain_refine=False):
    """Registration of the CT (scenario sc); writes and returns the registration record.

    Parameters
    ----------
    sc : str
        'anat_same_space' (same scanner space as the anatomy, the AC-PC matrix is copied),
        '1v2b' (rigid CT -> anatomy AC-PC), '2v1' (CT -> CT template), '2v2' (CT -> CT study template).
    ct_img : str
        Raw CT (already reoriented), used as the moving image of the registrations and as the reference
        of the record.
    raw : str
        Raw CT before reorientation (kept in the record for traceability).
    base : str
        Base name of the CT (BIDS entities without suffix).
    paths : dict
        Output paths (getpath.ct).
    anat : dict | None
        Anatomical files in AC-PC space (see modalities.pet._0_Pipeline_launcher._anat_files) or None
        when there is no anatomy.
    packages : list | None
        Output of backtonative.get (or None without anatomy).
    template_CT : str
        CT template in the reference template space (scenario 2v1) or ''.
    reference : str
        Reference template name.
    study_template : str
        CT study template image (scenario 2v2) or ''.
    T_CT : str
        Anatomical contrast used as the registration target ('T1w' or 'T2w').
    metric : str
        Metric of the CT -> anatomy registration (typically 'mattes').
    type_of_transform_template : str
        ANTs transform used for scenario 2v1 (CT -> CT template). Default 'SyNRA' (rigid + affine + SyN).
    do_manual : bool
        Manual correction after the automatic CT -> anatomy registration (freeview, scenario 1v2b).
    Need_check_coreg : bool
        Open freeview to visually check every target registration.
    bids_dir : str
        BIDS dataset root (for the QC copies).
    diary_name : str
        Diary file for the log messages.
    sing_fs : str
        Singularity prefix for FreeSurfer (freeview).
    template_mask : str
        Brain mask of the reference template (target mask of scenario 2v1, also read by the PET pipeline).
    coreg_transform : str
        Any ANTsPy `type_of_transform` (see modalities.pet.coreg2anat.ANTS_TRANSFORMS). Used for the
        CT -> anatomy registration (scenario 1v2b). 'Rigid' (default) is the standard for CT -> anatomy;
        'Affine' compensates a size difference (atrophy, scanner FOV); 'SyNRA' is rare (cross-scanner,
        template building). Scenario 2v1 uses `type_of_transform_template`; 'anat_same_space' and
        '2v2' do not use it.
    coreg_init : str
        Initialisation of the CT -> anatomy registration (modalities.pet.coreg2anat.init_candidates): 'auto'
        (centres of mass of the CT head without the bed and of the anatomy, centres of the fields of view, AC-PC
        matrix of the anatomy, identity, search over the rotations, centroids of the intracranial cavity of the CT
        and of the brain of the anatomy; the best start is kept), 'com', 'fov', 'acpc', 'identity', 'search' or
        'brain'.
    ct_brain_refine : bool
        Scenario 1v2b: second pass in the brain after the head registration: CT in the brain window (brain_window)
        within the brain of the anatomy brought into the CT and dilated by 3 mm (refine_in_brain), from the first
        pass. The grey / white matter and the ventricles then drive the registration instead of the skull.
    """
    _validate_coreg_transform(coreg_transform, diary_name)

    nl = 'Registration of the CT ' + base + ' (scenario ' + sc + ')'
    run_cmd.msg(nl, diary_name, 'HEADER')
    mat_dir = paths['orig_matrices']
    os.makedirs(mat_dir, exist_ok=True)
    os.makedirs(paths['orig_formc'], exist_ok=True)
    # clip to [-1000, 2000] HU: keeps the bone, removes the metallic artefacts and the air outside
    work = clipped(ct_img, opj(paths['orig_formc'], base + '_desc-clipped_ct.nii.gz'))
    targets = {}

    if sc == 'anat_same_space':
        # CT and anatomy in the same scanner space: the CT follows the AC-PC matrix of the anatomy
        acpc = transfo.check(transfo.affine(opj(anat['dir_transfo'], 'acpc_0GenericAffine.mat')))

    elif sc == '1v2b':
        # CT and anatomy in different spaces: rigid CT -> anatomy, metric in the brain mask of the anatomy dilated
        # by 10 mm (skull included: the bone edges constrain the registration), CT restricted to the head (the bed
        # and the head holder removed, > -500 HU, largest component), multi-start initialisation. rigid_to shifts
        # the CT to >= 0 for the centre of mass: in HU the air (-1000) is a negative mass larger than the head, and
        # the centre of mass of the raw CT lies away from the head (ANTs' default initialisation).
        prefix = opj(mat_dir, base + '_from-ct-to-' + T_CT + '_')
        head = ct_head_mask(work)
        ants.image_write(head, opj(paths['orig_formc'], base + '_desc-head_mask.nii.gz'), ri=False)
        # intracranial cavity (skull ring filled slice by slice): 'brain' start and centre of the rotation search,
        # independent of the face, neck and shoulders in the field of view
        cavity = ct_intracranial_mask(work)
        if cavity is not None:
            ants.image_write(cavity, opj(paths['orig_formc'], base + '_desc-intracranial_mask.nii.gz'), ri=False)
        else:
            run_cmd.msg('WARNING: no closed skull found in the CT (bone > 300 HU): no intracranial start',
                        diary_name, 'WARNING')
        acpc_mat = opj(anat['dir_transfo'], 'acpc_0GenericAffine.mat')
        acpc = rigid_to(anat['acpc'], work, prefix, metric, mask_file=anat['mask'], dilate_mm=10.0,
                        coreg_transform=coreg_transform, moving_mask=head, init=coreg_init,
                        acpc=acpc_mat if opi(acpc_mat) else None, diary_name=diary_name, moving_brain=cavity)
        if ct_brain_refine:
            win = brain_window(work, opj(paths['orig_formc'], base + '_desc-brainwindow_ct.nii.gz'))
            acpc = refine_in_brain(anat['acpc'], anat['mask'], win, acpc, prefix + 'desc-brain_', metric, 3.0,
                                   coreg_transform, opj(paths['orig_formc'], base + '_desc-coregbrain_mask.nii.gz'),
                                   diary_name)
        if do_manual:
            # manual correction done on the full CT so the engineer sees the skull; the chain is
            # world-to-world, valid either way, and the re-registration uses initial_transform
            acpc = manual_correction(work, anat['acpc'], acpc, prefix,
                                     opj(paths['orig_formc'], 'ct2' + T_CT + '.nii.gz'),
                                     'ct', diary_name, sing_fs)

    if sc in ('anat_same_space', '1v2b'):
        targets['acpc'] = dict(image=anat['acpc'], brain=anat['brain'], mask=anat['mask'], **acpc)
        for pkg in packages or []:
            targets[pkg[0]] = dict(image=pkg[1], brain=pkg[1], mask=pkg[2],
                                   **transfo.compose(acpc, transfo.from_backtonative(pkg)))

    elif sc == '2v1':
        # CT -> CT template (already in the reference template space): the transform level is
        # controlled by type_of_transform_template, not by coreg_transform (template-building case)
        prefix = opj(mat_dir, base + '_to_' + reference + '_')
        own = template_CT.replace('.nii.gz', '_mask.nii.gz')
        transfo.register(template_CT, work, prefix, type_of_transform_template, aff_metric=metric)
        targets[reference] = dict(image=template_CT, brain=template_CT,
                                  mask=own if opi(own) else (template_mask if template_mask and opi(template_mask) else ''),
                                  **transfo.registration(prefix))

    elif sc == '2v2':
        # CT -> CT study template (registration already done in step 2: studyT.make_template)
        prefix = opj(mat_dir, base + '_to_studyTemplate_')
        targets['studyTemplate'] = dict(image=study_template, brain=study_template, mask='',
                                        **transfo.registration(prefix))

    else:
        raise ValueError(f'unknown CT scenario {sc}')

    for t in targets.values():
        transfo.check(t)
    record = dict(modality='ct', scenario=sc, session=base, reference=ct_img, raw=raw, targets=targets,
                  coreg_transform=coreg_transform, coreg_init=coreg_init, ct_brain_refine=bool(ct_brain_refine),
                  created=pu.now())
    transfo.write_record(record_path(paths, base), record)

    # quality control: CT in every target space, brain mask of the target (or head of a study template) as contour
    os.makedirs(paths['qc'], exist_ok=True)
    for name, t in targets.items():
        out_dir = paths['acpc'] if name == 'acpc' else opj(paths['templates'], name)
        os.makedirs(out_dir, exist_ok=True)
        dst = pu.out_name(out_dir, base, space=name, suffix='ct')
        transfo.apply(t['image'], ct_img, t, 'to', 'linear', out=dst)
        pu.write_json(dst, [ct_img, t['image']], f'CT in the {name} space (ANTsPy, linear).', Transforms=t['to'])
        png = opj(paths['qc'], base + '_space-' + name + '_ct.png')
        QC_plot.mosaic(dst, qc_contour(t, paths['qc'], name), png)
        qc_dir = opj(bids_dir, 'QC', 'ct_in_' + ('acpc' if name == 'acpc' else 'template'))
        os.makedirs(qc_dir, exist_ok=True)
        shutil.copyfile(png, opj(qc_dir, opb(png)))
        if Need_check_coreg:
            FS_freeview.check_coreg(dst, t['image'], 'ct', sing_fs, transform_volume=False)

    # anatomy in the CT space (visualisation)
    if 'acpc' in targets:
        anat_in_ct = pu.out_name(paths['orig_process'], base, space='orig', suffix=T_CT)
        transfo.apply(ct_img, anat['acpc'], targets['acpc'], 'from', 'linear', out=anat_in_ct)
        pu.write_json(anat_in_ct, [anat['acpc'], ct_img], 'Anatomy in the CT space (visualisation).')
    return record_path(paths, base)