"""
Boundary check of a PET -> anatomy registration: is the PET brain where the MRI brain is, or larger / smaller /
shifted?

A dim rim of FDG in the outer cortex, with the signal left mostly in the sulci, has two causes that look alike on a
QC image:
  1. partial volume: the macaque cortex (1.5-2.5 mm) is thinner than the PET resolution, the uptake of the gyral
     crowns spreads into the CSF and the skull, two facing banks of a sulcus add up (Hoffman et al. 1979, J Comput
     Assist Tomogr 3:299; Rousset et al. 1998, J Nucl Med 39:904). With a perfect registration the half maximum of
     the PET lies on the pial boundary and the outer cortex is dim. No registration fixes it: partial-volume
     correction (step 'pvc') or a deeper surface sampling does.
  2. registration: a PET brain larger than the MRI brain (scale: PET voxel size of the header, gradient
     non-linearity of the MRI; or a fit that slid the PET outward) puts its dim rim on the MRI cortex, a translation
     does it on one side only.

The check separates them. The PET resampled on the anatomy is averaged in shells of signed distance to the brain
boundary (negative inside), in six sectors (+-x, +-y, +-z of the AC-PC space, by the direction of the outward normal
of the boundary). The expected profile is the FDG-like image of the anatomy (grey 1, white 0.35, CSF 0.05;
coreg2anat.synthetic_pet) smoothed to the PET resolution, or the smoothed brain mask when no tissue masks exist. In
both profiles the edge is the half-maximum point of the outer fall (between the maximum inside and the minimum just
outside): smoothing a step leaves its half maximum on the step, whatever the width, so the difference of the two
edges is the displacement of the PET boundary in that sector, positive when the PET brain extends beyond the MRI
brain. Median of the six = dilation (scale-like error); half difference of opposite sectors = translation.

Accuracy on a folded-cortex phantom (2 mm cortex, PET FWHM 3-4.5 mm, tests/unit/test_pet.py): within 0.1 mm of 0 for
a correct registration in the sectors without extracerebral uptake; 70-85 % of the true mean normal displacement for
a 4 % scale or a 1.5 mm shift (the estimate is conservative); an assumed FWHM of 2 to 6 mm for a true 4.5 mm changes
the shifts by < 0.1 mm. An uptake 2 mm outside the brain (temporal muscles) biases its sectors inward by 0.3-0.4 mm:
the dorsal sector (+z) is the most reliable.
"""
import json
import os

import ants
import numpy as np
from scipy import ndimage

opj = os.path.join
opi = os.path.isfile

from modalities.pet import transfo

AXES = ('x', 'y', 'z')
SECTORS = tuple(s + a for a in AXES for s in '+-')
STEP_MM = 0.5
RANGE_MM = (-6.0, 5.0)          # profile: 6 mm inside the brain, 5 mm outside
EDGE_MM = (-4.0, 4.5)           # edge search window around the boundary of the MRI brain


def _gauss(img, fwhm_mm):
    return ants.smooth_image(img, fwhm_mm / (2.0 * np.sqrt(2.0 * np.log(2.0))), sigma_in_physical_coordinates=True)


def _first_file(*names):
    """First existing file among `names` ('' if none)."""
    return next((f for f in names if f and opi(f)), '')


def tissue_image(anat, brain):
    """FDG-like image of the anatomy before smoothing (grey / white masks), else the brain mask; with its model."""
    from modalities.pet import coreg2anat
    gm_f = _first_file(anat.get('gray_wb'), anat.get('gray'))
    wm_f = _first_file(anat.get('white_full'), anat.get('white'))
    if not (gm_f and wm_f):
        return brain.clone('float'), 'brain mask'
    w = coreg2anat.FDG_WEIGHTS
    s = np.where(brain.numpy() > 0, w['CSF'], 0.0).astype(np.float32)
    s[coreg2anat._on_grid(wm_f, brain).numpy() > 0] = w['WM']
    s[coreg2anat._on_grid(gm_f, brain).numpy() > 0] = w['GM']
    return brain.new_image_like(s), 'grey / white matter'


def expected_image(anat, brain, fwhm_mm):
    """FDG-like image of the anatomy smoothed to the PET resolution, with its model (tissue_image)."""
    img, model = tissue_image(anat, brain)
    return _gauss(img, fwhm_mm), model


def estimate_fwhm(tissue, brain, pet, chain, dilate_mm=2.0, grid=None):
    """
    Effective resolution of the PET (FWHM, mm), estimated from the data: the FDG-like image of the anatomy is smoothed
    at every FWHM of `grid` (1-10 mm by 0.5), brought on the PET grid with the registration, and the FWHM kept is the
    one whose image correlates best with the PET in the brain dilated by dilate_mm (the extracerebral uptake further
    out is not modelled). The resolution of a PET-MRI pair matched by smoothing the segmented MRI to the PET (as the
    effective-resolution estimate of Joshi et al. 2009 on a phantom). An upper bound: a residual misregistration and a
    grey / white contrast other than 1 / 0.35 are compensated by more smoothing. Returns (fwhm, {fwhm: r}).
    """
    from modalities.pet import coreg2anat
    pet_img = ants.image_read(pet) if isinstance(pet, str) else pet
    if pet_img.dimension == 4:
        pet_img = ants.slice_image(pet_img, axis=3, idx=0)
    roi = coreg2anat._dilated(coreg2anat._binary(transfo.apply(pet_img, brain, chain, 'from', 'linear')), dilate_mm)
    sel = roi.numpy() > 0
    y = pet_img.numpy()[sel]
    grid = np.arange(1.0, 10.01, 0.5) if grid is None else grid
    curve = {}
    for f in grid:
        x = transfo.apply(pet_img, _gauss(tissue, float(f)), chain, 'from', 'linear').numpy()[sel]
        curve[round(float(f), 2)] = round(float(np.corrcoef(x, y)[0, 1]), 4)
    best = max(curve, key=curve.get)
    return best, curve


def signed_distance(brain):
    """Signed distance (mm) of the voxel centres to the boundary of the brain mask, negative inside: the boundary
    lies half a voxel beyond the last voxel of the mask."""
    m = brain.numpy() > 0.5
    sp = brain.spacing[:3]
    h = 0.5 * float(np.mean(sp))
    return np.where(m, h - ndimage.distance_transform_edt(m, sampling=sp),
                    ndimage.distance_transform_edt(~m, sampling=sp) - h)


def sectors(dist, brain):
    """Sector label (index of SECTORS) of every voxel: dominant world axis of the outward normal of the boundary (the
    gradient of the signed distance), so that a shift along an axis moves the boundary of its two sectors."""
    g = np.stack(np.gradient(dist.astype(np.float32), *brain.spacing[:3]))
    g = np.tensordot(np.asarray(brain.direction)[:3, :3], g, axes=1)
    ax = np.abs(g).argmax(axis=0)
    comp = np.take_along_axis(g, ax[None], axis=0)[0]
    return 2 * ax + (comp < 0).astype(int)


def profiles(values, dist, sec):
    """Mean of `values` per distance bin (centres, STEP_MM) for every sector; NaN for empty bins."""
    edges = np.arange(RANGE_MM[0], RANGE_MM[1] + STEP_MM, STEP_MM)
    centres = 0.5 * (edges[:-1] + edges[1:])
    b = np.digitize(dist, edges) - 1
    ok = (b >= 0) & (b < len(centres))
    out = {}
    for k, name in enumerate(SECTORS):
        sel = ok & (sec == k)
        s = np.bincount(b[sel], weights=values[sel], minlength=len(centres))
        n = np.bincount(b[sel], minlength=len(centres))
        out[name] = np.where(n > 0, s / np.maximum(n, 1), np.nan)
    return centres, out


def edge(d, p):
    """Position (mm) of the outer edge of a profile: where it falls, going outward from its maximum within
    EDGE_MM[0]..0, to half way between that maximum and its minimum within 0..EDGE_MM[1] (half maximum of a
    smoothed step; a shift of the profile moves this point by the same distance). NaN when not found."""
    ok = np.isfinite(p)
    if ok.sum() < 4:
        return np.nan
    dd = np.linspace(EDGE_MM[0], EDGE_MM[1], 401)
    pp = np.interp(dd, d[ok], p[ok])
    inner, outer = dd <= 0, dd >= 0
    i0 = int(np.argmax(np.where(inner, pp, -np.inf)))
    half = 0.5 * (pp[i0] + pp[outer].min())
    below = np.nonzero((np.arange(len(dd)) > i0) & (pp <= half))[0]
    if not below.size:
        return np.nan
    j = below[0]
    return float(np.interp(half, [pp[j], pp[j - 1]], [dd[j], dd[j - 1]]))


def edge_shift(d, observed, expected):
    """Displacement (mm) of the edge of the PET profile relative to the expected one, positive outward."""
    s = edge(d, observed) - edge(d, expected)
    return round(s, 2) if np.isfinite(s) else 0.0


def summary(shifts):
    """Dilation (median of the six shifts: two sectors biased by an extracerebral uptake, e.g. the temporal
    muscles, do not move it) and translation per axis (half difference of opposite sectors) from the shifts."""
    dil = float(np.median([shifts[s] for s in SECTORS]))
    tr = {a: round(0.5 * (shifts['+' + a] - shifts['-' + a]), 2) for a in AXES}
    return round(dil, 2), tr


def _verdict(dil, tr, tol_mm):
    """Text of the result: scale error, translation, or partial volume."""
    worst_tr = max(abs(v) for v in tr.values())
    if abs(dil) > tol_mm and abs(dil) >= worst_tr:
        return (f'PET brain {"larger" if dil > 0 else "smaller"} than the MRI brain by {abs(dil):.1f} mm per side: '
                f'scale error (check the PET voxel size; coreg_transform="Similarity")')
    if worst_tr > tol_mm:
        ax = max(tr, key=lambda a: abs(tr[a]))
        return f'PET shifted by {tr[ax]:+.1f} mm along {ax}: registration error (QC, coreg_init, pet_coreg_mask)'
    return (f'PET boundary on the MRI boundary (within {tol_mm:.1f} mm): a dim outer cortex is partial volume, '
            f'not registration (step pvc with pvc_fwhm)')


def check(anat_image, brain_mask, pet, chain, out_prefix, fwhm_mm=None, anat=None, tol_mm=0.5):
    """
    Boundary check of the registration `chain` (PET -> anatomy, 'to' direction) of the PET image `pet`.
    Writes <out_prefix>boundary.json and <out_prefix>boundary.png; returns the result dict with
    shifts (mm per sector, positive outward), dilation_mm, translation_mm, verdict, and the resolution estimated from
    the data (fwhm_estimated_mm, estimate_fwhm). fwhm_mm: PET resolution of the expected profile (pvc_fwhm); None:
    the estimated one.
    """
    brain = ants.image_read(brain_mask) if isinstance(brain_mask, str) else brain_mask
    fixed = ants.image_read(anat_image) if isinstance(anat_image, str) else anat_image
    brain = ants.resample_image_to_target(brain, fixed, interp_type='nearestNeighbor')
    tissue, model = tissue_image(anat or {}, brain)
    est, curve = estimate_fwhm(tissue, brain, pet, chain)
    fwhm = float(fwhm_mm) if fwhm_mm else est
    pet_in = transfo.apply(fixed, pet, chain, 'to', 'linear').numpy()
    dist = signed_distance(brain)
    sec = sectors(dist, brain)
    d, obs = profiles(pet_in, dist, sec)
    _, exp = profiles(_gauss(tissue, fwhm).numpy(), dist, sec)
    shifts = {s: edge_shift(d, obs[s], exp[s]) for s in SECTORS}
    dil, tr = summary(shifts)
    verdict = _verdict(dil, tr, tol_mm)
    res = dict(fwhm_mm=fwhm, fwhm_estimated_mm=est, fwhm_correlation=curve, model=model, shifts_mm=shifts,
               dilation_mm=dil, translation_mm=tr, tolerance_mm=tol_mm, verdict=verdict,
               note='positive shift: the PET brain extends beyond the MRI brain in that sector; fwhm_estimated_mm: '
                    'upper bound of the PET resolution (best correlation with the smoothed tissue image)')
    with open(out_prefix + 'boundary.json', 'w') as f:
        json.dump(res, f, indent=2)
    _plot(d, obs, exp, shifts, verdict, out_prefix + 'boundary.png')
    return res


def _unit(v, ok):
    """Profile scaled to 0-1 over the bins `ok`."""
    return (v - np.nanmin(v[ok])) / max(1e-12, np.nanmax(v[ok]) - np.nanmin(v[ok]))


def _plot(d, obs, exp, shifts, title, png):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 3, figsize=(11, 6), sharex=True)
    for ax, s in zip(axs.T.ravel(), SECTORS):
        o, e = obs[s], exp[s]
        ok = np.isfinite(o) & np.isfinite(e)
        if ok.sum() > 1:
            ax.plot(d, _unit(o, ok), 'k-', label='PET')
            ax.plot(d, _unit(e, ok), 'C1--', label='expected')
        ax.axvline(0, color='0.6', lw=0.8)
        ax.set_title(f'{s}: {shifts[s]:+.1f} mm', fontsize=9)
    axs[0, 0].legend(fontsize=8)
    for ax in axs[1]:
        ax.set_xlabel('distance to the brain boundary (mm)')
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(png, dpi=100)
    plt.close(fig)


def check_record(rec_file, fwhm_mm=None, bids_dir=None):
    """
    Boundary check of an existing registration record (step 6 output, <session>_desc-registration_pet.json): the
    session mean of the record in the AC-PC space of the anatomy, the brain mask of the record and, next to it, the
    grey / white matter masks of the anatomical pipeline. Writes <session>_desc-registration_boundary.json / .png next
    to the record, copies the png to pet/QC and, with bids_dir, to <bids_dir>/QC/pet_boundary. Returns the result.
    """
    import shutil
    rec = transfo.read_record(rec_file)
    t = rec['targets'].get('acpc')
    if t is None:
        raise ValueError(f'{rec_file}: no AC-PC target (registration without anatomy): no boundary check')
    mask = t['mask']
    anat = dict(gray_wb=mask.replace('_space-acpc_mask', '_desc-Whole_Brain_Gray_mask'),
                white_full=mask.replace('_space-acpc_mask', '_desc-White_mask'))
    prefix = rec_file.replace('_desc-registration_pet.json', '_desc-registration_')
    res = check(t['image'], mask, rec['reference'], t, prefix, fwhm_mm, anat=anat)
    qc = opj(os.path.dirname(os.path.dirname(os.path.dirname(rec_file))), 'QC')
    dests = [qc] + ([opj(bids_dir, 'QC', 'pet_boundary')] if bids_dir else [])
    for d in dests:
        os.makedirs(d, exist_ok=True)
        shutil.copyfile(prefix + 'boundary.png', opj(d, os.path.basename(prefix) + 'boundary.png'))
    return res
