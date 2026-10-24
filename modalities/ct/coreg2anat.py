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
from modalities.pet.coreg2anat import (manual_correction, qc_contour, rigid_to,
                                       _validate_coreg_transform)
from modalities.anat.freesurfer import FS_freeview

CT_CLIP = (-1000.0, 2000.0)


def clipped(ct_file, out):
    img = ants.image_read(ct_file)
    arr = np.clip(img.numpy(), *CT_CLIP)
    ants.image_write(ants.from_numpy(arr, origin=img.origin, spacing=img.spacing, direction=img.direction), out,
                     ri=False)
    return out


def record_path(paths, base):
    return opj(paths['orig_matrices'], base + '_desc-registration_ct.json')


def register(sc, ct_img, raw, base, paths, anat, packages, template_CT, reference, study_template, T_CT, metric,
             type_of_transform_template, do_manual, Need_check_coreg, bids_dir, diary_name, sing_fs,
             template_mask='', coreg_transform='Rigid'):
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
    """
    _validate_coreg_transform(coreg_transform)

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
        # CT and anatomy in different spaces: rigid (or higher) CT -> anatomy, with a brain mask dilated
        # by 10 mm (skull included: the bone edges constrain the registration)
        prefix = opj(mat_dir, base + '_from-ct-to-' + T_CT + '_')
        acpc = rigid_to(anat['acpc'], work, prefix, metric,
                        mask_file=anat['mask'], dilate_mm=10.0,
                        coreg_transform=coreg_transform)
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
                  coreg_transform=coreg_transform, created=pu.now())
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