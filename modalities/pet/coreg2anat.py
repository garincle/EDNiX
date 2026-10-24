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
    - multi-start initialisation (`coreg_init`): centres of mass of the masked images, centres of the fields of
      view, AC-PC matrix of the anatomy (the PET is near the native anatomy), identity; the start with the highest
      mutual information that covers the brain is kept (`<prefix>candidates.json`);
    - optional mask of the PET (`pet_coreg_mask`): percentile threshold + largest cluster, a mask file, or the
      brain of the anatomy back-projected after a first pass ('anat', second pass in the brain);
    - optional FDG-like target built from the tissue masks of the anatomy (`coreg_target='synthetic'`);
    - species-adaptive pyramid (`_adaptive_reg_params`): the coarsest shrink level keeps at least 32 voxels per
      dimension (mouse 50^3, macaque cropped to the brain, human 256^3);
    - ANTs includes `initial_transform` in the affine it writes at `outprefix`: the chain returned by
      `transfo.registration(prefix)` is the full transformation (no explicit composition).
"""
import os
import json
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
# Helpers: masks, adaptive pyramid, transform validation
# ---------------------------------------------------------------------------

def _read(img):
    return ants.image_read(img) if isinstance(img, str) else img


def _binary(img):
    return ants.threshold_image(img, 0.5, 1e30, 1, 0)


def _dilated(mask, mm):
    """Binary mask (path or ANTsImage) dilated by mm (at least one voxel)."""
    m = _binary(_read(mask))
    r = max(1, int(round(mm / min(m.spacing[:3]))))
    return ants.iMath(m, 'MD', r)


def _on_grid(mask, ref):
    """Binary mask (path or ANTsImage, any grid) resampled on the grid of `ref` by world coordinates."""
    return _binary(ants.resample_image_to_target(_read(mask), ref, interp_type='nearestNeighbor'))


def _nonnegative(img):
    """Image shifted so that its minimum is 0. The centre of mass used to initialise a registration is
    intensity-weighted: in a CT in HU the air (-1000 HU) is a negative mass larger than the head, and the centre
    of mass of the whole image ends up on the far side of the field of view, away from the head."""
    mn = float(img.min())
    return img - mn if mn < 0 else img


def threshold_mask(img, percentile=30.0, dilate_mm=3.0):
    """
    Mask of the brightest connected object of an image (PET, or CT shifted to >= 0): voxels above the
    `percentile` of the positive voxels, largest connected component, holes filled, dilated by dilate_mm. Same
    recipe as the previous PET scripts (np.percentile + nilearn threshold_img(cluster_threshold) + 3dmask_tool
    -fill_holes). For an FDG image of the head the brain is the brightest large object; with the torso in the
    field of view, crop first (crop_for_MoCo) or give the mask as a file.
    """
    img = _read(img)
    a = img.numpy()
    pos = a[np.isfinite(a) & (a > 0)]
    if pos.size == 0:
        raise ValueError('threshold_mask: the image has no positive voxel')
    m = ants.threshold_image(img, float(np.percentile(pos, percentile)), 1e30, 1, 0)
    m = ants.iMath(ants.iMath(m, 'GetLargestComponent'), 'FillHoles')
    return _dilated(m, dilate_mm) if dilate_mm > 0 else m


def ct_head_mask(ct, hu=-500.0):
    """Head (and body) of a CT in HU: > -500 HU, largest connected component (removes the bed and the head
    holder, which are separated from the head by air), holes filled (sinuses, airways)."""
    ct = _read(ct)
    m = ants.threshold_image(ct, hu, 1e30, 1, 0)
    return ants.iMath(ants.iMath(m, 'GetLargestComponent'), 'FillHoles')


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


# Values accepted by ants.registration(type_of_transform=...) in ANTsPy (allowable_tx of ants.registration).
# Composite types cascade internally in ANTs: e.g. 'SyNRA' = Rigid + Affine + SyN, 'TRSAA' = Translation + Rigid +
# Similarity + Affine + Affine.
ANTS_TRANSFORMS = (
    'Translation', 'Rigid', 'Similarity', 'Affine', 'AffineFast', 'BOLDAffine',
    'QuickRigid', 'DenseRigid', 'BOLDRigid', 'TRSAA',
    'SyN', 'SyNRA', 'SyNOnly', 'SyNCC', 'SyNAggro', 'SyNLessAggro', 'SyNabp', 'SyNBold', 'SyNBoldAff',
    'ElasticSyN', 'Elastic', 'TVMSQ', 'TVMSQC',
)
LINEAR_TRANSFORMS = ('Translation', 'Rigid', 'Similarity', 'Affine', 'AffineFast', 'BOLDAffine', 'QuickRigid',
                     'DenseRigid', 'BOLDRigid', 'TRSAA')


def _msg(text, diary_name, code='ENDC'):
    """run_cmd.msg, printed only when there is no diary (functions called outside a pipeline run)."""
    if diary_name:
        run_cmd.msg(text, diary_name, code)
    else:
        run_cmd.printcolor(text, code)


def _validate_coreg_transform(t, diary_name=''):
    if t not in ANTS_TRANSFORMS:
        raise ValueError(f'coreg_transform={t!r}: expected one of the ANTsPy type_of_transform {ANTS_TRANSFORMS}')
    if t not in LINEAR_TRANSFORMS:
        _msg(f'WARNING: coreg_transform={t!r} is non-linear: a PET or a CT of the same subject differs from '
                    f'its MRI by a rigid motion (plus small scanner scale errors); a deformation fitted on the low '
                    f'cortical contrast of a PET is not physiological. Inspect the QC.', diary_name, 'WARNING')
    elif t not in ('Translation', 'Rigid', 'QuickRigid', 'DenseRigid', 'BOLDRigid'):
        _msg(f'WARNING: coreg_transform={t!r} has scaling parameters: for the same subject they can absorb a '
                    f'misregistration (the PET brain shrunk or enlarged into the MRI brain); prefer Rigid unless a '
                    f'scanner scale error is documented.', diary_name, 'WARNING')
    return t


# ---------------------------------------------------------------------------
# Initialisation candidates and multi-start registration
# ---------------------------------------------------------------------------

INIT_MODES = ('auto', 'com', 'fov', 'acpc', 'identity', 'search')


def _grid_center(img):
    """World (ANTs, LPS) coordinates of the centre of the grid of img."""
    sp = np.asarray(img.spacing[:3], float)
    n = np.asarray(img.shape[:3], float) - 1.0
    return np.asarray(img.origin[:3], float) + np.asarray(img.direction, float)[:3, :3] @ (sp * n / 2.0)


def _translation_file(path, t):
    """ANTs matrix file of the translation x -> x + t (fixed points -> moving points, LPS)."""
    tr = ants.create_ants_transform(transform_type='AffineTransform', dimension=3,
                                    translation=tuple(float(v) for v in t))
    ants.write_transform(tr, path)
    return path


def _rigid_part(tr):
    """Closest rotation (SVD projection) of the linear part of an ANTs affine transform, same centre and
    translation: the starts of a rigid registration must not carry a scale or a shear (ANTs writes the initial
    transform into the result)."""
    p = np.asarray(tr.parameters, float)
    u, _, vt = np.linalg.svd(p[:9].reshape(3, 3))
    r = u @ np.diag([1.0, 1.0, np.sign(np.linalg.det(u @ vt))]) @ vt
    return ants.create_ants_transform(transform_type='AffineTransform', dimension=3, matrix=r,
                                      translation=tuple(p[9:12]), center=tuple(np.asarray(tr.fixed_parameters, float)))


def rotation_search(fixed, fixed_mask, moving, path, step_deg=15.0, arc_fraction=0.5, voxel_mm=None):
    """
    Start found by a search over the rotations (ANTs antsAffineInitializer: every rotation of the sphere by
    `step_deg` within arc_fraction x 180 degrees around the centres of mass, a short local optimisation of the mutual
    information from each, the best kept), for a head rotated far from the anatomy (a CT
    head flexed on the head holder, a PET of an animal in another position). Images resampled to `voxel_mm`
    (default 3 voxels of the fixed image, at least the coarsest voxel of the moving image) to keep it fast. The
    result is reduced to its rotation (_rigid_part) and written to `path`.
    """
    vox = voxel_mm or max(3.0 * min(fixed.spacing[:3]), max(moving.spacing[:3]))
    f = ants.resample_image(_nonnegative(fixed), (vox,) * 3, use_voxels=False, interp_type=0)
    m = ants.resample_image(moving, (vox,) * 3, use_voxels=False, interp_type=0)
    mask = None
    if fixed_mask is not None:
        mask = _binary(ants.resample_image_to_target(fixed_mask, f, interp_type='nearestNeighbor'))
    tmp = path.replace('.mat', '_affine.mat')
    ants.affine_initializer(f, m, search_factor=step_deg, radian_fraction=arc_fraction, use_principal_axis=False,
                            local_search_iterations=10, mask=mask, txfn=tmp)
    ants.write_transform(_rigid_part(ants.read_transform(tmp)), path)
    os.remove(tmp)
    return path


def init_candidates(fixed, fixed_mask, moving, prefix, init='auto', acpc=None, initial=None):
    """
    Initial transforms (fixed points -> moving points) of rigid_to, as a list of (name, [files]):
      com       centres of mass of the fixed image in its mask and of the (masked) moving image: moments
                initialisation (Alpert et al. 1990, J Nucl Med 31:1717; AFNI 3dAllineate -cmass);
      fov       centres of the two fields of view (the moving image cropped to its mask): 'FOV to FOV', for images
                whose intensity distributions are too different for the moments;
      acpc      AC-PC matrix of the anatomy: the PET / CT is assumed to be near the native anatomy (consistent scanner
                coordinates, similar head position), the matrix the anatomical pipeline already estimated is reused;
      identity  world coordinates of the headers as they are;
      search    best start of a search over the rotations (rotation_search): large rotations between the images;
      initial   a given transform (manual correction, first pass).
    init: 'auto' (every available one), one of the names above, or 'initial' (only `initial`).
    """
    if init not in INIT_MODES + ('initial',):
        raise ValueError(f'coreg_init={init!r}: expected one of {INIT_MODES}')
    names = ('com', 'fov', 'acpc', 'identity', 'search') if init == 'auto' else (init,)
    out = []
    if initial:
        out.append(('initial', [initial] if isinstance(initial, str) else list(initial)))
    if init == 'initial':
        if not out:
            raise ValueError('init_candidates: init="initial" without an initial transform')
        return out
    fx = _nonnegative(fixed)
    for name in names:
        if name == 'com':
            f = fx * fixed_mask if fixed_mask is not None else fx
            t = np.asarray(ants.get_center_of_mass(moving)) - np.asarray(ants.get_center_of_mass(f))
            out.append((name, [_translation_file(prefix + 'init-com.mat', t)]))
        elif name == 'fov':
            t = _grid_center(moving) - _grid_center(fixed)
            out.append((name, [_translation_file(prefix + 'init-fov.mat', t)]))
        elif name == 'acpc':
            if acpc and opi(acpc):
                out.append((name, [acpc]))
            elif init == 'acpc':
                raise FileNotFoundError(f'coreg_init="acpc" but no AC-PC matrix: {acpc}')
        elif name == 'identity':
            out.append((name, [_translation_file(prefix + 'init-identity.mat', (0.0, 0.0, 0.0))]))
        elif name == 'search':
            out.append((name, [rotation_search(fixed, fixed_mask, moving, prefix + 'init-search.mat')]))
    return out


def _score(fixed, moving, chain, fixed_mask, coverage_mask, moving_support):
    """(Mattes mutual information, coverage) of a registration: the ANTs metric value of the fixed image and the
    resampled moving image in the fixed mask (negative: lower is better), and the fraction of `coverage_mask`
    (brain of the fixed image) inside the resampled support of the moving image (its mask / field of view).
    MI can be increased by sliding the object out of the mask, hence the coverage condition."""
    warped = transfo.apply(fixed, moving, chain, 'to', 'linear')
    mi = float(ants.image_similarity(fixed, warped, metric_type='MattesMutualInformation', fixed_mask=fixed_mask))
    sup = transfo.apply(fixed, moving_support, chain, 'to', 'linear').numpy() > 0.5
    cm = coverage_mask.numpy() > 0.5
    cov = float(np.logical_and(sup, cm).sum()) / max(1.0, float(cm.sum()))
    return mi, cov


def rigid_to(fixed_file, moving_file, prefix, metric, mask_file=None, dilate_mm=5.0, initial=None,
             coreg_transform='Rigid', moving_mask=None, init='auto', acpc=None, min_coverage=0.9, diary_name='',
             verbose=False, same_subject=True):
    """
    Registration moving (PET, CT) -> fixed (anatomy, CT), multi-start; returns the chains of `prefix`.

    mask_file   brain mask of the fixed image: path (dilated by dilate_mm) or ANTsImage (used as is). The metric is
                computed in it.
    moving_mask mask of the moving image (path or ANTsImage, any grid; 3dAllineate -source_mask): the moving image
                is multiplied by it and cropped to its box, and the metric ignores the voxels outside it. Removes
                the torso, the bed and the extracerebral uptake from both the metric and the initialisation.
    init        initialisation candidates (init_candidates): 'auto' = centres of mass, centres of the fields of
                view, AC-PC matrix (`acpc`), identity, search over the rotations, plus `initial` if given.

    For every candidate: Translation, then `coreg_transform`, with the adaptive pyramid; ANTs includes the initial
    transform in the matrix it writes. The candidate kept has the lowest Mattes metric value (highest mutual
    information, Wells et al. 1996; Maes et al. 1997) among those whose moving image covers at least
    `min_coverage` of the fixed brain; a multi-start search against the local optima of mutual information
    (Jenkinson & Smith 2001, FLIRT; AFNI 3dAllineate -twobest). Scores in <prefix>candidates.json.
    same_subject=False (template to template): no warning on the scale / non-linear transforms.
    """
    if same_subject:
        _validate_coreg_transform(coreg_transform, diary_name)
    elif coreg_transform not in ANTS_TRANSFORMS:
        raise ValueError(f'coreg_transform={coreg_transform!r}: expected one of {ANTS_TRANSFORMS}')
    fixed = _read(fixed_file)
    moving = _nonnegative(_read(moving_file))
    if moving.dimension == 4:
        raise ValueError(f'rigid_to: the moving image must be 3D ({moving_file})')

    # mask of the fixed image (metric) and brain of the fixed image (coverage)
    fmask, cover = None, None
    if isinstance(mask_file, ants.ANTsImage):
        fmask = _binary(mask_file)
        cover = fmask
    elif mask_file and opi(mask_file):
        cover = _binary(ants.image_read(mask_file))
        fmask = _dilated(cover, dilate_mm)

    # moving image masked and cropped to its mask (the transforms are world-to-world: valid for the full image)
    mmask = None
    if moving_mask is not None and not (isinstance(moving_mask, str) and not moving_mask):
        mmask = _on_grid(moving_mask, moving)
        if float(mmask.sum()) == 0:
            raise ValueError('rigid_to: the moving mask does not overlap the moving image')
        moving = ants.crop_image(moving * mmask, label_image=mmask, label=1)
        mmask = ants.crop_image(mmask, label_image=mmask, label=1)
    support = mmask if mmask is not None else moving * 0 + 1
    if cover is None:
        cover = _binary(_nonnegative(fixed) * 0 + 1)

    kw = {}
    if fmask is not None:
        kw['mask'] = fmask
    if mmask is not None:
        kw['moving_mask'] = mmask
    adaptive = _adaptive_reg_params(fixed, moving)
    cands = init_candidates(fixed, fmask, moving, prefix, init, acpc, initial)

    results = []
    for name, init_files in cands:
        cp = prefix + 'init-' + name + '_'
        try:
            transfo.register(fixed, moving, cp + 'desc-translation_', 'Translation', aff_metric=metric,
                             initial_transform=init_files, **adaptive, **kw)
            chain = transfo.register(fixed, moving, cp, coreg_transform, aff_metric=metric,
                                     initial_transform=[cp + 'desc-translation_' + transfo.ANTS_AFF], **adaptive, **kw)
            mi, cov = _score(fixed, moving, chain, fmask, cover, support)
        except Exception as e:                           # a start that diverges (empty overlap) is skipped
            _msg(f'WARNING: registration start {name} failed: {e}', diary_name, 'WARNING')
            continue
        results.append(dict(init=name, init_files=init_files, metric=mi, coverage=cov, prefix=cp))
        if verbose:
            print(f'[rigid_to] {opb(prefix)} start {name}: MI metric {mi:.4f}, coverage {cov:.3f}')
    if not results:
        raise RuntimeError(f'registration failed for every start: {moving_file} -> {fixed_file}')

    valid = [r for r in results if r['coverage'] >= min_coverage]
    best = min(valid, key=lambda r: r['metric']) if valid else max(results, key=lambda r: r['coverage'])
    msg = ', '.join(f"{r['init']}: MI {r['metric']:.4f} / coverage {r['coverage']:.2f}" for r in results)
    _msg(f"INFO: {opb(prefix)}: starts [{msg}] -> kept {best['init']}", diary_name, 'OKGREEN')
    if not valid:
        _msg(f"WARNING: no start covers {min_coverage:.0%} of the brain of {opb(str(fixed_file))}: kept the "
                    f"best coverage ({best['coverage']:.2f}); check the orientation of the images, the masks and the "
                    f"QC (do_manual=True to correct by hand)", diary_name, 'WARNING')

    transfo.clear(prefix)
    for suffix in (transfo.ANTS_AFF, transfo.ANTS_WARP, transfo.ANTS_INVWARP):
        if opi(best['prefix'] + suffix):
            shutil.copyfile(best['prefix'] + suffix, prefix + suffix)
    for r in results:                                    # the starts that were not kept
        transfo.clear(r['prefix'])
        transfo.clear(r['prefix'] + 'desc-translation_')
    with open(prefix + 'candidates.json', 'w') as f:
        json.dump(dict(fixed=str(fixed_file), moving=str(moving_file), transform=coreg_transform, metric=metric,
                       kept=best['init'], starts=[{k: v for k, v in r.items() if k != 'prefix'} for r in results]),
                  f, indent=2)
    return transfo.registration(prefix)


# ---------------------------------------------------------------------------
# PET mask, second pass in the brain, synthetic FDG target
# ---------------------------------------------------------------------------

PET_MASKS = ('none', 'threshold', 'anat')


def pet_mask(mode, pet, out, percentile=30.0, dilate_mm=3.0):
    """
    Mask of the PET used by the registration (moving mask, `pet_coreg_mask`):
      'none'       no mask (the crop box of the motion correction is still applied by loop2);
      'threshold'  percentile threshold + largest connected component + filled holes (threshold_mask);
      'anat'       no mask for the first pass, then the brain of the anatomy back-projected (refine_in_brain);
      <file>       a mask drawn on the PET (any grid in the PET world space).
    Returns the mask file or None.
    """
    if not mode or mode in ('none', 'anat'):
        return None
    if mode == 'threshold':
        ants.image_write(threshold_mask(pet, percentile, dilate_mm), out, ri=False)
        pu.write_json(out, [pet], f'PET mask for the registration: > {percentile} percentile of the positive voxels, '
                                  f'largest component, holes filled, dilated {dilate_mm} mm.')
        return out
    if opi(mode):
        return mode
    raise ValueError(f'pet_coreg_mask={mode!r}: expected {PET_MASKS} or a mask file')


def refine_in_brain(fixed, brain_mask, moving, chain, prefix, metric, dilate_mm, coreg_transform, out_mask,
                    diary_name=''):
    """
    Second pass in the brain: the brain mask of the anatomy is brought into the PET space with the first
    registration (chain), dilated by dilate_mm (PET resolution, residual error of the first pass), and the
    registration is estimated again from the first one on the PET brain only (moving mask) within the same dilated
    brain of the anatomy (fixed mask). The extracerebral uptake (temporal muscles, Harderian / lacrimal glands, eye
    muscles, nasal mucosa) and the scalp no longer pull the cortex of the PET onto the CSF or the skull of the MRI.
    Same idea as the previous PET scripts (template mask back-projected to the subject, 3dAllineate -source_mask).
    """
    mv = _read(moving)
    bm = _dilated(_binary(transfo.apply(mv, brain_mask, chain, 'from', 'linear')), dilate_mm)
    ants.image_write(bm, out_mask, ri=False)
    pu.write_json(out_mask, [brain_mask, moving], f'Brain of the anatomy in the PET space (first registration), '
                                                  f'dilated {dilate_mm} mm: moving mask of the second pass.')
    return rigid_to(fixed, moving, prefix, metric, mask_file=_dilated(brain_mask, dilate_mm),
                    initial=[f for f, _ in chain['to']], init='initial', moving_mask=bm,
                    coreg_transform=coreg_transform, diary_name=diary_name)


COREG_TARGETS = ('anat', 'brain', 'synthetic')
FDG_WEIGHTS = dict(GM=1.0, WM=0.35, CSF=0.05)


def synthetic_pet(anat, pet_spacing, out, fwhm_mm=None):
    """
    FDG-like image of the anatomy (coreg_target='synthetic'): grey matter 1, white matter 0.35, CSF and the rest of
    the brain 0.05, outside the brain 0 (grey / white FDG uptake ratio about 3, Phelps et al. 1979; Heiss et al.
    1984), smoothed to the PET resolution (Gaussian of FWHM fwhm_mm, default twice the largest PET voxel). The PET is
    registered to an image of its own contrast built from the segmentation of the MRI (Ashburner & Friston 1997,
    NeuroImage 6:209): the position of the cortical ribbon drives the alignment, not the T1 contrast of the skull
    and scalp. FDG only (other tracers have other grey / white contrasts).
    """
    brain = _binary(ants.image_read(anat['mask']))
    gm_f = next((f for f in (anat.get('gray_wb'), anat.get('gray')) if f and opi(f)), '')
    wm_f = next((f for f in (anat.get('white_full'), anat.get('white')) if f and opi(f)), '')
    if not gm_f or not wm_f:
        raise FileNotFoundError("coreg_target='synthetic' needs the grey and white matter masks of the anatomy")
    gm = _on_grid(gm_f, brain).numpy() > 0
    wm = _on_grid(wm_f, brain).numpy() > 0
    s = np.where(brain.numpy() > 0, FDG_WEIGHTS['CSF'], 0.0).astype(np.float32)
    s[wm] = FDG_WEIGHTS['WM']
    s[gm] = FDG_WEIGHTS['GM']
    fwhm = fwhm_mm or 2.0 * float(max(pet_spacing[:3]))
    img = ants.smooth_image(brain.new_image_like(s), fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0))),
                            sigma_in_physical_coordinates=True)
    ants.image_write(img, out, ri=False)
    pu.write_json(out, [anat['mask'], gm_f, wm_f], f'Synthetic FDG image of the anatomy (GM {FDG_WEIGHTS["GM"]}, '
                                                   f'WM {FDG_WEIGHTS["WM"]}, CSF {FDG_WEIGHTS["CSF"]}), FWHM {fwhm} mm.')
    return out


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
             moving_for_estimate=None, coreg_transform='Rigid', coreg_init='auto', pet_coreg_mask='none',
             pet_mask_percentile=30.0, pet_mask_dilate_mm=3.0, coreg_target='anat', synthetic_fwhm=None):
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
    coreg_init : str
        Initialisation of the PET -> anatomy / CT registration (init_candidates): 'auto' (centres of mass, centres
        of the fields of view, AC-PC matrix of the anatomy, identity; the best start is kept), 'com', 'fov', 'acpc'
        or 'identity'.
    pet_coreg_mask : str
        Mask of the PET for the registration (pet_mask): 'none', 'threshold' (pet_mask_percentile, largest cluster,
        filled holes, dilated pet_mask_dilate_mm), 'anat' (second pass on the brain of the anatomy back-projected
        into the PET and dilated pet_mask_dilate_mm), or a mask file in the PET space.
    coreg_target : str
        Fixed image of the 1v2a registration: 'anat' (AC-PC head image), 'brain' (skull-stripped anatomy) or
        'synthetic' (FDG-like image of the tissue masks, synthetic_pet; FDG only, smoothed to synthetic_fwhm mm,
        default twice the largest PET voxel).
    """
    _validate_coreg_transform(coreg_transform, diary_name)
    if coreg_target not in COREG_TARGETS:
        raise ValueError(f'coreg_target={coreg_target!r}: expected one of {COREG_TARGETS}')

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
        head = _dilated(ct_head_mask(ct_img), 3.0)          # head of the CT (> -500 HU, without the bed)
        return rigid_to(ct_img, moving, opj(mat_dir, session + '_from-pet-to-ct_'), coreg_metric,
                        mask_file=head, coreg_transform=coreg_transform, init=coreg_init,
                        moving_mask=pet_mask(pet_coreg_mask if pet_coreg_mask != 'anat' else 'none', moving,
                                             opj(paths['orig_formc'], session + '_desc-coreg_mask.nii.gz'),
                                             pet_mask_percentile, pet_mask_dilate_mm),
                        diary_name=diary_name)

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
            # anatomy from another session: rigid PET -> anatomy (multi-start, optional PET mask / second pass in the
            # brain / FDG-like target), then optional manual correction
            prefix = opj(mat_dir, session + '_from-pet-to-' + T_PET + '_')
            os.makedirs(paths['orig_formc'], exist_ok=True)
            if coreg_target == 'synthetic':
                fixed = synthetic_pet(anat, ants.image_read(moving).spacing,
                                      opj(paths['orig_formc'], session + '_desc-synthetic_' + T_PET + '.nii.gz'),
                                      synthetic_fwhm)
            else:
                fixed = anat['brain'] if coreg_target == 'brain' else anat['acpc']
            run_cmd.msg(f'INFO: PET -> anatomy: target {coreg_target} ({opb(fixed)}), init {coreg_init}, PET mask '
                        f'{pet_coreg_mask}', diary_name, 'OKGREEN')
            acpc_mat = opj(anat['dir_transfo'], 'acpc_0GenericAffine.mat')
            two_pass = pet_coreg_mask == 'anat'
            acpc = rigid_to(fixed, moving, prefix + ('pass-1_' if two_pass else ''), coreg_metric,
                            mask_file=anat['mask'], coreg_transform=coreg_transform, init=coreg_init,
                            acpc=acpc_mat if opi(acpc_mat) else None,
                            moving_mask=pet_mask(pet_coreg_mask, moving,
                                                 opj(paths['orig_formc'], session + '_desc-coreg_mask.nii.gz'),
                                                 pet_mask_percentile, pet_mask_dilate_mm),
                            diary_name=diary_name)
            if two_pass:
                acpc = refine_in_brain(fixed, anat['mask'], moving, acpc, prefix, coreg_metric, pet_mask_dilate_mm,
                                       coreg_transform,
                                       opj(paths['orig_formc'], session + '_desc-coregbrain_mask.nii.gz'), diary_name)
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
                  coreg_transform=coreg_transform, coreg_init=coreg_init, pet_coreg_mask=pet_coreg_mask,
                  coreg_target=coreg_target, created=pu.now())
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