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

Registration of 1v2a and pet_to_ct (1v2b, 2v1, 2v2) goes through `rigid_to`:
    - species-adaptive pyramid (`_adaptive_reg_params`): the coarsest shrink level always keeps at
      least 32 voxels per dimension, so the same code works on a mouse (50^3), a macaque cropped to
      the brain (73x54x39) or a human (256^3);
    - `coreg_transform` option (any ANTs `type_of_transform`: 'Rigid', 'Affine', 'TRSAA', 'SyNRA',
      'BsplineSyN', 'SyNCC'...) to escalate the transform when a rigid alignment is not enough
      (atrophy between sessions, cross-scanner PET/MR, template building);
    - automatic quality check (`_registration_quality`) + retry with a fallback metric (meansquares):
      the moving image must cover the target mask, otherwise the registration is retried and the
      best result is kept;
    - **no double composition**: ANTs already includes `initial_transform` in the affine it writes at
      `outprefix`; the chain returned by `transfo.registration(prefix)` is the full transformation.
      Composing the translation with it would apply the translation twice (bug from an earlier
      implementation).
"""
import os
import shutil

import numpy as np
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


# ---------------------------------------------------------------------------
# Scenario selection
# ---------------------------------------------------------------------------

def scenario(IhaveanANAT, anat_pet_same_space, IhaveaCT, ct_targets, pet_template, study_template):
    """Scenario code of pet4.svg from the launcher options (ct_targets: target names of the CT record)."""
    if IhaveanANAT:
        if anat_pet_same_space:
            return '1v1'
        if IhaveaCT and ct_targets is None:
            raise FileNotFoundError('IhaveanANAT=True and IhaveaCT=True but no CT registration record: '
                                    'run the CT pipeline first')
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


# ---------------------------------------------------------------------------
# Helpers: mask dilation, adaptive pyramid, quality, transform validation
# ---------------------------------------------------------------------------

def _dilated(mask_file, mm):
    m = ants.image_read(mask_file)
    r = max(1, int(round(mm / min(m.spacing[:3]))))
    return ants.iMath(ants.threshold_image(m, 0.5, 1e9, 1, 0), 'MD', r)


def _adaptive_reg_params(fixed_img, moving_img):
    """Registration parameters adapted to the smaller of the two images: the shrink factors are
    chosen so that the coarsest pyramid level keeps at least 32 voxels per dimension. This is what
    makes the registration robust across species (human 256^3, macaque 128^3, rodent 50^3) and for
    a PET cropped to the brain only (73x54x39): the default ANTsPy pyramid (8, 4, 2, 1) would shrink
    a 54-voxel axis to 7 voxels at the first level, where a Mattes metric has too few samples to
    constrain six parameters."""
    min_dim = min(min(fixed_img.shape[:3]), min(moving_img.shape[:3]))
    if min_dim >= 128:
        return dict(aff_shrink_factors=(8, 4, 2, 1), aff_smoothing_sigmas=(4, 2, 1, 0),
                    aff_iterations=(1000, 500, 500, 250), aff_random_sampling_rate=1.0)
    if min_dim >= 64:
        return dict(aff_shrink_factors=(4, 2, 1), aff_smoothing_sigmas=(2, 1, 0),
                    aff_iterations=(1000, 500, 250), aff_random_sampling_rate=1.0)
    if min_dim >= 32:
        return dict(aff_shrink_factors=(2, 1), aff_smoothing_sigmas=(1, 0),
                    aff_iterations=(500, 250), aff_random_sampling_rate=1.0)
    return dict(aff_shrink_factors=(1,), aff_smoothing_sigmas=(0,),
                aff_iterations=(250,), aff_random_sampling_rate=1.0)


def _registration_quality(warped_moving, brain_mask):
    """Quality = mean intensity of the resampled moving inside the target mask, normalised by the
    mean of the moving on its positive voxels. Close to 1 when the PET covers the brain (good
    registration), close to 0 when it is placed outside (failed registration)."""
    w = warped_moving.numpy()
    m = brain_mask.numpy() > 0.5
    if m.sum() == 0:
        return 0.0
    inside = w[m]
    inside = inside[np.isfinite(inside)]
    positive = w[np.isfinite(w) & (w > 0)]
    if inside.size == 0 or positive.size == 0:
        return 0.0
    denom = float(np.mean(positive))
    return float(np.mean(inside)) / denom if denom > 0 else 0.0


# Any value accepted by ants.registration(type_of_transform=...). Composite types cascade
# internally in ANTs: e.g. 'SyNRA' = Rigid + Affine + SyN, 'TRSAA' = Translation + Rigid + Affine.
ANTS_TRANSFORMS = (
    'Translation', 'Rigid', 'Similarity', 'Affine', 'AffineFast',
    'QuickRigid', 'DenseRigid', 'BOLDRigid', 'TRSAA',
    'Bspline', 'BsplineSyN', 'SyN', 'SyNRA', 'SyNOnly', 'SyNCC',
    'Exponential', 'Polynomial', 'TVMSQ', 'TVMSQC', 'SyNBoldAff', 'SyNBoldAff', "SyNRA"
)

# transforms that are slow: warn if chosen
SLOW_TRANSFORMS = ('SyN', 'SyNRA', 'SyNOnly', 'SyNCC', 'BsplineSyN', 'TVMSQ', 'TVMSQC', 'Exponential', 'Polynomial')


def _validate_coreg_transform(t):
    if t not in ANTS_TRANSFORMS:
        raise ValueError(f'coreg_transform={t!r}: expected one of ANTs type_of_transform {ANTS_TRANSFORMS}')
    if t in SLOW_TRANSFORMS:
        run_cmd.msg(f'WARNING: coreg_transform={t!r} is a non-linear ANTs transform: expect 5-10x the '
                    f'rigid time. Inspect the QC carefully (PET has little structural contrast in the '
                    f'cortex; deformations can be non-physiological).', '', 'WARNING')
    return t


# ---------------------------------------------------------------------------
# Debug: log every ANTs call
# ---------------------------------------------------------------------------

import functools


def _log_ants_calls(enable=True):
    """Patch ants.registration / ants.apply_transforms to print their arguments.
    Idempotent: enabling twice does not stack the patches."""
    if getattr(_log_ants_calls, '_patched', False) == enable:
        return
    if enable:
        _log_ants_calls._orig_reg = ants.registration
        _log_ants_calls._orig_apply = ants.apply_transforms

        @functools.wraps(_log_ants_calls._orig_reg)
        def reg_log(*args, **kwargs):
            print('[ants.registration] args=', args, 'kwargs=', kwargs)
            return _log_ants_calls._orig_reg(*args, **kwargs)

        @functools.wraps(_log_ants_calls._orig_apply)
        def apply_log(*args, **kwargs):
            print('[ants.apply_transforms] args=', args, 'kwargs=', kwargs)
            return _log_ants_calls._orig_apply(*args, **kwargs)

        ants.registration = reg_log
        ants.apply_transforms = apply_log
        _log_ants_calls._patched = True
    else:
        ants.registration = _log_ants_calls._orig_reg
        ants.apply_transforms = _log_ants_calls._orig_apply
        _log_ants_calls._patched = False


# ---------------------------------------------------------------------------
# rigid_to: translation + coreg_transform, adaptive pyramid, quality check
# ---------------------------------------------------------------------------

def rigid_to(fixed_file, moving_file, prefix, metric, mask_file=None, dilate_mm=5.0, initial=None,
             save_translated=True, translated_interp='linear', verbose=True, min_quality=0.3,
             coreg_transform='Rigid'):
    """
    PET -> target registration with species-adaptive parameters and automatic metric retry.

    1. Translation (coarse pre-alignment), then `coreg_transform`, with a pyramid whose coarsest
       level keeps >= 32 voxels per dimension (see `_adaptive_reg_params`): same code for a mouse,
       a macaque, a human, or a PET cropped to the brain. ANTs composes `initial_transform` with
       the found transform internally; the chain returned by `transfo.registration(prefix)` is
       already the full translation + transform (no explicit compose).

    2. Quality check: the resampled moving image must cover the target mask. On failure, retry
       with the same transform but the `meansquares` metric (the Mattes mutual information is
       undefined when the overlap is empty). The best of the two is kept.

    Parameters
    ----------
    coreg_transform : str
        Any ANTsPy `type_of_transform` (see ANTS_TRANSFORMS above). Common choices for PET -> anatomy:
          'Rigid'  (default) rigid 6 DOF, standard for PET -> anatomy on the same subject;
          'Affine'          12 DOF, when the two brains differ in size (atrophy between sessions, scanner FOV);
          'SyNRA'           rigid + affine + SyN in one call (ANTs cascade): rare, template building
                            or cross-scanner PET/MR; slow, inspect the QC;
          'TRSAA'           translation + rigid + affine in one call: linear alternative to SyNRA.
        Any other ANTs transform (BsplineSyN, SyNCC, TVMSQ...) is passed verbatim.
    """
    _validate_coreg_transform(coreg_transform)

    fixed_img = ants.image_read(fixed_file) if isinstance(fixed_file, str) else fixed_file
    moving_img = ants.image_read(moving_file) if isinstance(moving_file, str) else moving_file

    kw = {}
    brain_mask_quality = None
    if isinstance(mask_file, ants.ANTsImage):
        kw['mask'] = mask_file
        brain_mask_quality = mask_file
    elif mask_file and opi(mask_file):
        brain_mask_quality = ants.image_read(mask_file)
        kw['mask'] = _dilated(mask_file, dilate_mm)

    adaptive = _adaptive_reg_params(fixed_img, moving_img)
    if verbose:
        _log_ants_calls(True)
        print(f'[rigid_to] fixed       shape  = {fixed_img.shape}')
        print(f'[rigid_to] moving      shape  = {moving_img.shape}')
        print(f'[rigid_to] metric             = {metric}')
        print(f'[rigid_to] coreg_transform    = {coreg_transform}')
        print(f'[rigid_to] pyramid            = {adaptive["aff_shrink_factors"]}')

    # 1. Translation (pre-alignment)
    trans_prefix = prefix + 'desc-translation_'
    trans_chain = transfo.register(fixed_file, moving_file, trans_prefix, 'Translation',
                                    aff_metric=metric, **adaptive, **kw)
    trans_files = [f for f, _ in trans_chain['to']]
    if verbose:
        print(f'[rigid_to] translation        -> {trans_files[0]}')

    if save_translated:
        translated_out = prefix + 'desc-translation.nii.gz'
        transfo.apply(fixed_file, moving_file, trans_chain, 'to', translated_interp, out=translated_out)
        if verbose:
            print(f'[rigid_to] translated image   -> {translated_out}  (exists={opi(translated_out)})')

    # 2. Two attempts with the user-chosen transform: metric, then meansquares fallback
    init_list = ([initial] if isinstance(initial, str) else list(initial or [])) + trans_files
    attempts = [('metric', metric), ('meansquares', 'meansquares')]

    best = None
    best_quality = -1.0
    for name, met in attempts:
        attempt_prefix = prefix + 'attempt-' + coreg_transform.lower() + '-' + name + '_'
        try:
            if verbose:
                print(f'[rigid_to] attempt {coreg_transform}/{name} (metric={met})')
            # ANTs already includes initial_transform in the affine it writes: chain is the full
            # translation + transform, no explicit compose needed.
            chain = transfo.register(fixed_file, moving_file, attempt_prefix, coreg_transform,
                                      aff_metric=met, initial_transform=init_list, **adaptive, **kw)
            if brain_mask_quality is not None:
                warped = transfo.apply(fixed_file, moving_file, chain, 'to', translated_interp)
                q = _registration_quality(warped, brain_mask_quality)
            else:
                q = 1.0
            if verbose:
                print(f'[rigid_to] attempt {coreg_transform}/{name}: quality = {q:.3f}')
            if q > best_quality:
                best_quality = q
                best = (coreg_transform.lower() + '-' + name, chain)
            if q >= min_quality:
                break
        except Exception as e:
            if verbose:
                print(f'[rigid_to] attempt {coreg_transform}/{name} failed: {e}')

    if best is None:
        raise RuntimeError(f'registration failed for {moving_file}: all attempts raised')

    name, _ = best
    if verbose:
        print(f'[rigid_to] best = {name} (quality {best_quality:.3f})')
    if best_quality < min_quality:
        run_cmd.msg(f'WARNING: registration quality {best_quality:.3f} < {min_quality} '
                    f'({opb(moving_file)} -> {opb(fixed_file)}): result may be wrong; '
                    f'consider a stronger coreg_transform (Affine, SyNRA) or do_manual=True', '', 'WARNING')

    # 3. Copy the winning attempt files under the standard prefix
    transfo.clear(prefix)
    for suffix in (transfo.ANTS_AFF, transfo.ANTS_WARP, transfo.ANTS_INVWARP):
        src = prefix + 'attempt-' + name + '_' + suffix
        if opi(src):
            shutil.copyfile(src, prefix + suffix)

    if verbose:
        _log_ants_calls(False)
    # ANTs includes initial_transform in the affine at prefix: the chain is the full PET -> target
    # transformation. DO NOT compose with trans_chain again (would apply the translation twice).
    return transfo.registration(prefix)


# ---------------------------------------------------------------------------
# Manual correction (scenario 1v2a, do_manual=True)
# ---------------------------------------------------------------------------

def manual_correction(pet_mean, fixed_file, chain, prefix, pet_in_fixed, imgtype, diary_name, sing_fs):
    """
    PET engineer's manual correction: the PET in the anatomical space is opened in freeview with its
    'Transform Volume' tool; when the user saves it (same name), the registration is estimated again
    from the edited image, starting from the current one (PET -> PET, mean squares).
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


# ---------------------------------------------------------------------------
# Main entry point: register()
# ---------------------------------------------------------------------------

def register(sc, session_mean, session, paths, anat, packages, ct_record, pet_template, reference,
             study_template, T_PET, ct_pet_same_space, coreg_metric, type_of_transform_template, do_manual,
             Need_check_coreg, bids_dir, diary_name, sing_fs, template_mask='', anchor=None,
             moving_for_estimate=None, coreg_transform='Rigid'):
    """
    Registration of the session (scenario `sc`); writes and returns the registration record.

    Parameters
    ----------
    sc : str
        Scenario code from `scenario()`.
    session_mean : str
        Path of the session mean image (after between-runs MoCo). Kept as the reference of the record.
    paths, session : str
        Output paths (getpath.pet) and session name.
    anat : dict | None
        Anatomical files in AC-PC space (see _0_Pipeline_launcher._anat_files) or None.
    packages : list | None
        Output of backtonative.get (or None without anatomy).
    ct_record : dict | None
        Registration record of the CT pipeline (modality 'ct'), or None.
    pet_template : str
        PET template image (scenario 3v1) or ''.
    reference : str
        Reference template name.
    study_template : str
        PET study template image (scenario 3v2) or ''.
    T_PET : str
        Anatomical contrast used as the registration target ('T1w' or 'T2w').
    ct_pet_same_space : bool
        PET and CT share the scanner space.
    coreg_metric : str
        Metric of the rigid PET -> anatomy / CT (typically 'mattes').
    type_of_transform_template : str
        Transform used for scenario 3v1 (PET -> PET template).
    do_manual : bool
        Manual correction after the automatic registration (freeview, scenario 1v2a).
    Need_check_coreg : bool
        Open freeview to visually check every target registration.
    bids_dir : str
        BIDS dataset root (for the QC copies).
    diary_name : str
        Diary file for the log messages.
    sing_fs : str
        Singularity prefix for FreeSurfer (freeview).
    template_mask : str
        Brain mask of the reference template (target mask of scenarios 2v1 and 3v1).
    anchor : dict | None
        Chains session -> scanner space (MoCo.scanner_anchor), used when PET and anatomy / CT
        share the scanner space (1v1, ct_pet_same_space). None otherwise.
    moving_for_estimate : str | None
        Image used as the moving image of the registration (typically the session mean cropped to
        the MoCo crop box). ANTs initialises a rigid registration from the centre of mass of the
        *whole* moving image; a full-FOV PET whose bladder, liver or injection site dominates the
        signal is therefore placed at a biased starting position even when the metric is restricted
        to a brain mask. The transforms produced are world-to-world and remain valid for the full
        session mean, kept as the reference of the record. Defaults to session_mean when not given.
    coreg_transform : str
        Any ANTsPy `type_of_transform` (see `ANTS_TRANSFORMS`). Used for the PET -> anatomy / CT
        registration (scenarios 1v2a, 1v2b, 2v1, 2v2). 'Rigid' (default) is the standard for
        PET -> anatomy on the same subject; 'Affine' compensates a size difference; 'SyNRA' is a
        composite (rigid + affine + SyN) for rare cross-scanner cases (slow). Scenario 3v1 uses
        `type_of_transform_template` instead.
    """
    _validate_coreg_transform(coreg_transform)

    nl = 'Registration of the PET session ' + session + ' (scenario ' + sc + ' of pet4.svg)'
    run_cmd.msg(nl, diary_name, 'HEADER')
    mat_dir = paths['orig_matrices']
    os.makedirs(mat_dir, exist_ok=True)
    targets = {}

    # moving image of every registration: the cropped session mean when given, full session mean otherwise
    moving = moving_for_estimate or session_mean
    if moving_for_estimate and os.path.normpath(moving_for_estimate) != os.path.normpath(session_mean):
        run_cmd.msg('INFO: registration estimated on ' + opb(moving_for_estimate) + ' (cropped session mean); '
                    'the record keeps the full session mean as its reference', diary_name, 'OKGREEN')
    run_cmd.msg('INFO: coreg_transform = ' + coreg_transform + ' (PET -> anatomy / CT)', diary_name, 'OKGREEN')

    def pet_to_ct():
        """PET -> CT chain (scanner position of the CT acquisition when they share the space)."""
        ct_img = ct_record['reference']
        if ct_pet_same_space:
            return anchor or {'to': [], 'from': []}
        ct = ants.image_read(ct_img)
        head = ants.iMath(ants.threshold_image(ct, -500.0, 1e9, 1, 0), 'MD', 3)   # head of the CT (> -500 HU)
        return rigid_to(ct_img, moving, opj(mat_dir, session + '_from-pet-to-ct_'), coreg_metric,
                        mask_file=head, coreg_transform=coreg_transform)

    # ------------------------------------------------------------------
    # Branch 1: anatomy present (1v1, 1v2a, 1v2b)
    # ------------------------------------------------------------------
    if sc in ('1v1', '1v2a', '1v2b'):
        if sc == '1v1':
            # hybrid PET-MR: the AC-PC matrix of the anatomy is used as-is (the PET was reoriented
            # like the anatomy at step 1); no registration is estimated
            acpc = transfo.check(transfo.affine(opj(anat['dir_transfo'], 'acpc_0GenericAffine.mat')))
            if anchor:
                acpc = transfo.compose(anchor, acpc)

        elif sc == '1v2a':
            # anatomy from another session: rigid (or higher) PET -> anatomy, then optional manual correction
            prefix = opj(mat_dir, session + '_from-pet-to-' + T_PET + '_')
            acpc = rigid_to(anat['acpc'], moving, prefix, coreg_metric,
                            mask_file=anat['mask'], coreg_transform=coreg_transform)
            if do_manual:
                # manual correction done on the FULL session mean so the engineer sees the whole field
                # of view; the chain is world-to-world, valid either way, and the re-registration uses
                # initial_transform so ANTs' centre-of-mass initialisation is bypassed
                acpc = manual_correction(session_mean, anat['acpc'], acpc, prefix,
                                         opj(paths['orig_formc'], 'pet2' + T_PET + '.nii.gz'), 'pet',
                                         diary_name, sing_fs)

        else:  # sc == '1v2b'
            # anatomy + CT: PET -> CT (same space: nothing; else rigid) then CT -> anatomy (record of the CT pipeline)
            if 'acpc' not in ct_record['targets']:
                raise ValueError('IhaveanANAT and IhaveaCT: the CT record has no AC-PC target (run the CT pipeline '
                                 'with IhaveanANAT=True)')
            ct_acpc = ct_record['targets']['acpc']
            if os.path.normpath(ct_acpc['image']) != os.path.normpath(anat['acpc']):
                raise ValueError(f"the CT was registered to {ct_acpc['image']} but the anatomy of this PET session is "
                                 f"{anat['acpc']}: run the CT pipeline with the anatomy of the PET session")
            acpc = transfo.compose(pet_to_ct(), {'to': ct_acpc['to'], 'from': ct_acpc['from']})

        # targets: AC-PC of the anatomy, then each template package from backtonative
        targets['acpc'] = dict(image=anat['acpc'], brain=anat['brain'], mask=anat['mask'], **acpc)
        for pkg in packages or []:
            name = pkg[0]
            chain = transfo.compose(acpc, transfo.from_backtonative(pkg))
            targets[name] = dict(image=pkg[1], brain=pkg[1], mask=pkg[2], **chain)

    # ------------------------------------------------------------------
    # Branch 2: CT only (2v1, 2v2)
    # ------------------------------------------------------------------
    elif sc in ('2v1', '2v2'):
        to_ct = pet_to_ct()
        for name, t in ct_record['targets'].items():
            chain = transfo.compose(to_ct, {'to': t['to'], 'from': t['from']})
            mask = t.get('mask') or (template_mask if name == reference else '')
            targets[name] = dict(image=t['image'], brain=t.get('brain', t['image']), mask=mask, **chain)

    # ------------------------------------------------------------------
    # Branch 3: PET template (3v1)
    # ------------------------------------------------------------------
    elif sc == '3v1':
        prefix = opj(mat_dir, session + '_to_' + reference + '_')
        own = pet_template.replace('.nii.gz', '_mask.nii.gz')
        # whole head (no brain mask: the scalp and skull signal constrains the affine stages; a brain mask
        # dilated by 5-10 mm gave 1.4-4.6 mm errors against < 1 mm on the test phantom)
        chain = transfo.register(pet_template, moving, prefix, type_of_transform_template,
                                 aff_metric=coreg_metric)
        targets[reference] = dict(image=pet_template, brain=pet_template,
                                  mask=own if opi(own) else template_mask, **chain)

    # ------------------------------------------------------------------
    # Branch 4: PET study template (3v2, no atlas)
    # ------------------------------------------------------------------
    elif sc == '3v2':
        # registration of the motion-corrected session mean to the PET study template (done in
        # loop2.run before this function is called; the step-2 matrices
        # <session>_to_studyTemplate_ of the first reoriented run are kept for the crop of the
        # motion correction)
        prefix = study_prefix(paths, session)
        targets['studyTemplate'] = dict(image=study_template, brain=study_template, mask='',
                                        **transfo.registration(prefix))

    else:
        raise ValueError(f'unknown scenario {sc}')

    # sanity check: every transform file referenced by the record must exist
    for name, t in targets.items():
        transfo.check(t)

    record = dict(modality='pet', scenario=sc, session=session, reference=session_mean, targets=targets,
                  coreg_transform=coreg_transform, created=pu.now())
    rec_file = record_path(paths, session)
    transfo.write_record(rec_file, record)

    # quality control: PET in every target space, brain mask of the target (or head of a study
    # template) as contour
    os.makedirs(paths['qc'], exist_ok=True)
    for name, t in targets.items():
        pet_in = opj(paths['orig_formc'], session + '_space-' + name + '_desc-mean_pet.nii.gz')
        transfo.apply(t['image'], session_mean, t, 'to', petinterp, out=pet_in)
        png = opj(paths['qc'], session + '_space-' + name + '_pet.png')
        QC_plot.mosaic(pet_in, qc_contour(t, paths['qc'], name), png)
        qc_dir = opj(bids_dir, 'QC', 'pet_in_' + ('acpc' if name == 'acpc' else 'template'))
        os.makedirs(qc_dir, exist_ok=True)
        shutil.copyfile(png, opj(qc_dir, opb(png)))
        if Need_check_coreg:
            FS_freeview.check_coreg(pet_in, t['image'], 'pet', sing_fs, transform_volume=False)
    return rec_file


# ---------------------------------------------------------------------------
# Helpers: QC contour, study prefix, record path, tissue fractions, atlases
# ---------------------------------------------------------------------------

def qc_contour(target, qc_dir, name):
    """Contour of the QC images: brain mask of the target space, else a head mask of the target image
    (study template without mask, ants.get_mask)."""
    if target.get('mask') and opi(target['mask']):
        return target['mask']
    out = opj(qc_dir, name + '_desc-qccontour_mask.nii.gz')
    if not opi(out):
        ants.image_write(ants.get_mask(ants.image_read(target['image'])), out, ri=False)
    return out


def study_prefix(paths, session):
    """Prefix of the registration of the motion-corrected session mean to the PET study template (3v2)."""
    return opj(paths['orig_matrices'], session + '_desc-MoCo_to_studyTemplate_')


def record_path(paths, session):
    return opj(paths['orig_matrices'], session + '_desc-registration_pet.json')


def tissue_fraction(src, ref, chain, pet_spacing, out):
    """
    Fraction of every PET voxel covered by a binary mask of the target space: the mask is averaged over
    the extent of a PET voxel (Gaussian of the variance of a box of the PET voxel size, sigma = voxel /
    sqrt(12)) before the linear resampling into the PET space (a plain linear resampling samples the
    mask at the PET voxel centres only). Used by the partial-volume correction.
    """
    m = ants.threshold_image(ants.image_read(src), 0.5, 1e9, 1, 0)
    sigma = [float(v) / np.sqrt(12.0) for v in pet_spacing[:3]]
    sm = ants.smooth_image(m, sigma, sigma_in_physical_coordinates=True)
    f = transfo.apply(ref, sm, chain, 'from', 'linear')
    a = np.clip(f.numpy(), 0.0, 1.0)
    ants.image_write(f.new_image_like(a), out, ri=False)
    return out


def to_pet_space(record, session, paths, anat, template_files, list_atlas, diary_name, T_PET='T1w'):
    """
    Atlases, masks and anatomy (or template) in the PET space of the session (orig/labels, orig/masks,
    orig/preprocessing). Atlases: nearest-label interpolation (genericLabel); masks: linear interpolation
    then threshold 0.5 (fraction of the PET voxel covered); 4D atlases (levels) volume by volume.
    Tissue fractions (label-GM/WM/brain _probseg, tissue_fraction) of the whole-brain grey matter, the
    non-eroded white matter and the brain mask, for the partial-volume correction.
    template_files: dict(name, atlas_dir, species, brain, gray, gray_wb, white, csf, white_full) without anatomy.
    Returns dict(atlases={name: file}, masks={name: file}, probseg={GM, WM, brain: file}, anat=file).
    """
    nl = 'Atlases and masks into the PET space'
    run_cmd.msg(nl, diary_name, 'OKGREEN')
    for d in (paths['orig_labels'], paths['orig_masks'], paths['orig_process']):
        os.makedirs(d, exist_ok=True)
    ref = record['reference']
    out = dict(atlases={}, masks={}, probseg={}, anat='')
    if 'acpc' in record['targets'] and anat:
        t = record['targets']['acpc']
        atlas_src = {a: opj(anat['labels_dir'], anat['ID'] + '_seg-' + a + '_dseg.nii.gz') for a in list_atlas}
        mask_src = dict(brain=anat['mask'], gray=anat['gray'], white=anat['white'], csf=anat['csf'])
        frac_src = dict(GM=anat.get('gray_wb', ''), WM=anat.get('white_full', ''), brain=anat['mask'])
        anat_src, anat_out = anat['acpc'], pu.out_name(paths['orig_process'], session, space='orig', suffix=T_PET)
    elif template_files and template_files['name'] in record['targets'] and template_files.get('atlas_dir'):
        t = record['targets'][template_files['name']]
        atlas_src = {a: opj(template_files['atlas_dir'], template_files['species'] + '_seg-' + a + '_dseg.nii.gz')
                     for a in list_atlas}
        mask_src = {k: template_files.get(k, '') for k in ('brain', 'gray', 'white', 'csf')}
        frac_src = dict(GM=template_files.get('gray_wb', ''), WM=template_files.get('white_full', ''),
                        brain=template_files.get('brain', ''))
        anat_src = t['image']
        anat_out = pu.out_name(paths['orig_process'], session, space='orig', desc='template', suffix=T_PET)
    else:
        run_cmd.msg('WARNING: no atlas for this scenario (' + record['scenario'] + '): the regional analysis, the '
                    'brain ratio and the reference-region models need an atlas or masks', diary_name, 'WARNING')
        return out
    pet_spacing = ants.image_read(ref).spacing

    # tissue fractions for the partial-volume correction (GM / WM / brain)
    for name, src in frac_src.items():
        if src and opi(src):
            dst = opj(paths['orig_masks'], session + '_space-orig_label-' + name + '_probseg.nii.gz')
            tissue_fraction(src, ref, t, pet_spacing, dst)
            pu.write_json(dst, [src, ref], 'Fraction of the PET voxel in the mask (box-average, then linear '
                                           'resampling), for the partial-volume correction.')
            out['probseg'][name] = dst

    # atlases (nearest-label interpolation)
    for name, src in atlas_src.items():
        if not opi(src):
            run_cmd.msg('WARNING: atlas not found: ' + src, diary_name, 'WARNING')
            continue
        dst = opj(paths['orig_labels'], session + '_space-orig_seg-' + name + '_dseg.nii.gz')
        transfo.apply(ref, src, t, 'from', labelinterp, out=dst)
        pu.write_json(dst, [src, ref], 'Atlas in the PET space (genericLabel interpolation).')
        out['atlases'][name] = dst

    # masks (linear + threshold 0.5)
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

    # anatomical image in PET space (visualisation)
    if anat_src and opi(anat_src):
        transfo.apply(ref, anat_src, t, 'from', anatinterp, out=anat_out)
        pu.write_json(anat_out, [anat_src, ref], 'Anatomical image in the PET space (visualisation).')
        out['anat'] = anat_out
    return out