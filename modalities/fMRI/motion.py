"""
Motion estimation helpers for step 1.

Rigid registration (6 degrees of freedom) is the right model of head motion (fMRIPrep also estimates
6 rigid parameters per volume), but what it follows depends on the voxels it sees. On the whole
image, bright tissue that does not move with the brain (skull, muscle, skin, specially with a
surface coil) dominates the metric and the estimated motion of the brain is under-estimated. A brain
mask restricts the metric to the brain, but only works if

  * the mask lies INSIDE the brain, further from its border than the head moves (a mask that touches
    the border lets the bright tissue enter it as soon as the brain moves), and
  * the registration has a proper multi-resolution schedule: the 'BOLDRigid' preset of ANTsPy has two
    levels and 100 + 20 iterations, which is not enough once a mask is used; correlation is also
    much more stable than mutual information on a small mask.

See docs/functional_pipeline.md ("Motion") for the simulations behind these choices.
"""
import glob
import os
import shutil

import ants
import nibabel as nib
import numpy as np
from scipy import ndimage as ndi

opj = os.path.join

# 4 levels, global correlation: used with a brain mask (mutual information diverged on small masks).
# All the voxels of the mask are used (the default samples 20 % at random): in simulations the error
# fell from 0.29 to 0.19 of the uncorrected misalignment, with less dispersion between repeats (and one
# repeat out of ten had failed at 20 %). The mask is small, so the cost stays moderate (~0.6 s per volume).
SCHEDULE = dict(type_of_transform='Rigid', aff_metric='GC', aff_iterations=(200, 100, 50, 20),
                aff_shrink_factors=(4, 2, 1, 1), aff_smoothing_sigmas=(2, 1, 0, 0),
                aff_random_sampling_rate=1.0)


def prepare_mask(mask_path, reference_path, erode_voxels=3, min_retained=0.8, min_voxels=200):
    """
    Brain mask on the grid of the series, eroded by `erode_voxels`.
    Returns (ANTs mask image, dict(volume_mm3_full, volume_mm3_eroded, voxels_eroded)).
    Raises ValueError when the mask cannot be used: other orientation frame (the final functional
    mask is in the frame of the REORIENTED series), not in the field of view, or too small once eroded.
    """
    ref_nib, mask_nib = nib.load(reference_path), nib.load(mask_path)
    ax_ref, ax_mask = ''.join(nib.aff2axcodes(ref_nib.affine)), ''.join(nib.aff2axcodes(mask_nib.affine))
    if ax_ref != ax_mask:
        raise ValueError(f'the mask is in the orientation frame {ax_mask} and the series in {ax_ref}: give a mask in the '
                         f'frame of the reoriented series (the one written by step 3)')
    series = ants.image_read(reference_path)
    grid = ants.slice_image(series, axis=3, idx=0) if series.dimension == 4 else series
    mask = ants.image_read(mask_path)
    vol_before = float((mask.numpy() > 0.5).sum()) * float(np.prod(mask.spacing))
    resampled = ants.resample_image_to_target(mask, grid, interp_type='nearestNeighbor')
    arr = resampled.numpy() > 0.5
    vol_voxel = float(np.prod(grid.spacing))
    vol_after = float(arr.sum()) * vol_voxel
    if vol_before <= 0 or not (min_retained <= vol_after / vol_before <= 1.0 / min_retained):
        raise ValueError(f'the mask has {vol_before:.1f} mm3 and {vol_after:.1f} mm3 once on the grid of the series: '
                         f'it does not overlap the field of view (another frame?)')
    eroded = ndi.binary_erosion(arr, iterations=int(erode_voxels)) if erode_voxels > 0 else arr
    if int(eroded.sum()) < min_voxels:
        raise ValueError(f'the mask keeps {int(eroded.sum())} voxels once eroded by {erode_voxels} voxels '
                         f'(minimum {min_voxels}): lower motion_mask_erode')
    out = ants.from_numpy(eroded.astype('float32'), origin=grid.origin, spacing=grid.spacing, direction=grid.direction)
    return out, dict(volume_mm3_full=vol_after, volume_mm3_eroded=float(eroded.sum()) * vol_voxel,
                     voxels_eroded=int(eroded.sum()))


def estimate(series, reference, outprefix, mask=None):
    """Motion of every volume of `series` to `reference` (ANTs rigid). Without a mask: the 'BOLDRigid' preset."""
    if mask is None:
        return ants.motion_correction(image=series, fixed=reference, type_of_transform='BOLDRigid',
                                      outprefix=outprefix, verbose=False)
    return ants.motion_correction(image=series, fixed=reference, mask=mask, mask_all_stages=True,
                                  outprefix=outprefix, verbose=False, **SCHEDULE)


def reference_points(mask):
    """Eight points (physical, mm) spread inside the mask, to compare transforms where the brain is."""
    idx = np.argwhere(mask.numpy() > 0.5)
    centre, spread = idx.mean(axis=0), idx.std(axis=0) * 1.2
    pts = []
    for sx in (-1, 1):
        for sy in (-1, 1):
            for sz in (-1, 1):
                vox = centre + np.array([sx, sy, sz]) * spread
                pts.append(np.array(mask.origin) + (vox * np.array(mask.spacing)) @ np.array(mask.direction).T)
    return np.array(pts)


def discrepancy(dir_a, dir_b, points):
    """
    How differently two estimations of the same series place the brain, volume by volume, in mm
    (median, 90th percentile, maximum). A constant difference between the two reference poses is removed:
    only the relative alignment of the volumes matters.
    """
    fa, fb = sorted(glob.glob(opj(dir_a, 'img*.mat'))), sorted(glob.glob(opj(dir_b, 'img*.mat')))
    if not fa or len(fa) != len(fb):
        raise ValueError(f'cannot compare {len(fa)} and {len(fb)} motion matrices')
    d = np.array([[np.array(ants.read_transform(a).apply_to_point(tuple(p))) -
                   np.array(ants.read_transform(b).apply_to_point(tuple(p))) for p in points] for a, b in zip(fa, fb)])
    d = d - d.mean(axis=0, keepdims=True)
    e = np.sqrt((d ** 2).sum(-1).mean(-1))
    return float(np.median(e)), float(np.percentile(e, 90)), float(e.max())


def refine_with_mask(series_path, est_dir, mask_path, erode_voxels=3):
    """
    Re-estimate the motion of `series_path` with a brain mask and replace the matrices of `est_dir` by the
    result (the unmasked ones are kept in <est_dir>_unmasked). Returns dict(median, p90, max: difference
    between the two estimations in mm, volume_mm3_full, volume_mm3_eroded, voxels_eroded, motion_corrected).
    Raises ValueError (nothing changed) when the mask cannot be used.
    """
    mask, info = prepare_mask(mask_path, series_path, erode_voxels)
    series = ants.image_read(series_path)
    masked_dir, unmasked_dir = est_dir + '_masked', est_dir + '_unmasked'
    shutil.rmtree(masked_dir, ignore_errors=True)
    os.makedirs(masked_dir)
    result = estimate(series, ants.get_average_of_timeseries(series), opj(masked_dir, 'img'), mask)
    median, p90, worst = discrepancy(est_dir, masked_dir, reference_points(mask))
    shutil.rmtree(unmasked_dir, ignore_errors=True)
    os.rename(est_dir, unmasked_dir)
    os.rename(masked_dir, est_dir)
    return dict(info, median=median, p90=p90, max=worst, motion_corrected=result['motion_corrected'])


# =============================================================================
# Motion estimated by AFNI (3dvolreg) or FSL (MCFLIRT), converted to ANTs matrices
# =============================================================================
# Step 5 applies the motion of every volume in one interpolation with ANTs, from img_XXXXX_0GenericAffine.mat
# files. 3dvolreg and MCFLIRT write their transforms in their own conventions; the conversion is done by
# `nitransforms` (whose AFNI / FSL conventions are tested against the programs by its authors), then checked
# against the registered volumes the program itself wrote: the run stops if the converted matrices do not
# reproduce them, rather than silently applying a wrong or inverted motion.
EXTERNAL_METHODS = ('3dvolreg', 'mcflirt')
_LPS_RAS = np.diag([-1.0, -1.0, 1.0, 1.0])


def _nitransforms():
    try:
        import nitransforms as nt
    except ImportError as exc:
        raise ImportError("motion_method '3dvolreg' / 'mcflirt' needs the python package nitransforms "
                          "(pip install nitransforms)") from exc
    return nt


def external_command(method, series, reference, est_dir, corrected, sing_afni='', sing_fsl=''):
    """(command, matrices file or folder written by the program, registered series written by the program)."""
    if method == '3dvolreg':
        mats = opj(est_dir, 'volreg.aff12.1D')
        cmd = (f'{sing_afni}3dvolreg -prefix {corrected} -base {reference} -1Dfile {opj(est_dir, "volreg_params.1D")} '
               f'-1Dmatrix_save {mats} -twopass -zpad 4 -overwrite {series}')
        return cmd, mats, corrected
    if method == 'mcflirt':
        base = opj(est_dir, 'mcflirt')
        cmd = f'{sing_fsl}mcflirt -in {series} -reffile {reference} -out {base} -mats -plots -cost normcorr -dof 6'
        return cmd, base + '.mat', base + '.nii.gz'
    raise ValueError(f"unknown motion_method {method!r}: expected 'ants', '3dvolreg' or 'mcflirt'")


def read_ras_matrices(method, mats, reference_path, series_path):
    """(n_volumes, 4, 4) RAS matrices (reference to volume, nitransforms convention) read from the program's files."""
    nt = _nitransforms()
    if method == '3dvolreg':
        return np.asarray(nt.linear.load(mats, fmt='afni', reference=reference_path, moving=series_path).matrix)
    files = sorted(glob.glob(opj(mats, 'MAT_*')))
    if not files:
        raise FileNotFoundError(f'no MAT_* file in {mats}')
    return np.array([np.asarray(nt.linear.load(f, fmt='fsl', reference=reference_path, moving=series_path).matrix)
                     for f in files])


def fov_centre(grid):
    """Physical (LPS) centre of the field of view of an ANTs image."""
    return np.array(grid.origin) + ((np.array(grid.shape) - 1) / 2.0 * np.array(grid.spacing)) @ np.array(grid.direction).T


def write_itk_matrix(lps, path, centre):
    """
    Write a 4x4 LPS point mapping p_volume = A p_reference + b as an ANTs affine file, with the rotation
    centre `centre` (translation re-expressed about it, so that tx, ty, tz are the translation of the
    centre of the field of view, as in the ANTs estimation: FD takes them as such).
    """
    A, b = lps[:3, :3], lps[:3, 3]
    tx = ants.create_ants_transform(transform_type='AffineTransform', dimension=3)
    tx.set_parameters(list(A.ravel()) + list((b + (A - np.eye(3)) @ centre).ravel()))
    tx.set_fixed_parameters(list(centre))
    ants.write_transform(tx, path)


def _displacement_vox(lps, grid):
    """Largest displacement (voxels) of the corners of the central half of the field of view."""
    c, half = fov_centre(grid), (np.array(grid.shape) * np.array(grid.spacing)) / 4.0
    pts = np.array([c + np.array([sx, sy, sz]) * half for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)])
    moved = (lps[:3, :3] @ pts.T).T + lps[:3, 3]
    return float(np.linalg.norm(moved - pts, axis=1).max() / np.min(grid.spacing))


def choose_convention(itk_mats, series_path, registered_path, min_motion_vox=0.5, n_check=3, min_gap=2e-3, work_dir=None):
    """
    Whether the matrices of the program map the reference to the volume ('as_is') or the volume to the
    reference ('inverse'), decided on the volumes that moved most: the convention for which ANTs resampling
    of the original volume is closest to the registered volume written by the program.
    Returns (name, info). Raises ValueError when neither reproduces the program's volume.
    When no volume moved by `min_motion_vox` voxels the two conventions cannot be told apart and differ by
    little: 'as_is' is returned with info['conclusive'] = False.
    """
    import tempfile
    series, registered = ants.image_read(series_path), ants.image_read(registered_path)
    grid = ants.slice_image(registered, axis=3, idx=0)
    disp = np.array([_displacement_vox(m, grid) for m in itk_mats])
    order = np.argsort(disp)[::-1][:n_check]
    info = dict(max_displacement_vox=float(disp.max()), conclusive=False, scores={})
    if disp.max() < min_motion_vox:
        return 'as_is', info
    centre = fov_centre(grid)
    tmp = work_dir or tempfile.mkdtemp()
    candidates = {'as_is': lambda m: m, 'inverse': np.linalg.inv}

    def corr(a, b, keep):
        return float(np.corrcoef(a[keep], b[keep])[0, 1])

    scores, baseline = {name: [] for name in candidates}, []
    for k in order:
        vol, prog = ants.slice_image(series, axis=3, idx=int(k)), ants.slice_image(registered, axis=3, idx=int(k))
        keep = prog.numpy() > np.percentile(prog.numpy(), 30)
        baseline.append(corr(vol.numpy(), prog.numpy(), keep))
        for name, fn in candidates.items():
            f = opj(tmp, f'check_{name}_{k}.mat')
            write_itk_matrix(fn(itk_mats[k]), f, centre)
            moved = ants.apply_transforms(fixed=grid, moving=vol, transformlist=[f], interpolator='linear')
            scores[name].append(corr(moved.numpy(), prog.numpy(), keep))
    mean = {name: float(np.mean(v)) for name, v in scores.items()}
    info['scores'], info['uncorrected'] = mean, float(np.mean(baseline))
    best = max(mean, key=mean.get)
    other = min(mean, key=mean.get)
    if mean[best] <= info['uncorrected'] + 1e-3:
        raise ValueError(f"neither the matrices of the program nor their inverse reproduce its registered volumes "
                         f"(agreement {mean['as_is']:.4f} / {mean['inverse']:.4f}, uncorrected {info['uncorrected']:.4f}): "
                         f"the conversion cannot be trusted, use motion_method='ants'")
    info['conclusive'] = (mean[best] - mean[other]) >= min_gap
    return (best if info['conclusive'] else 'as_is'), info


def estimate_external(method, series_path, reference_path, est_dir, corrected_path, runner, sing_afni='', sing_fsl=''):
    """
    Estimate the motion of `series_path` with 3dvolreg or MCFLIRT, `runner(command)` executing the command, and
    leave in `est_dir` the img_XXXXX_0GenericAffine.mat files the rest of the pipeline reads. The registered
    series of the program is left in `corrected_path`. Returns a dict describing the convention found.
    """
    shutil.rmtree(est_dir, ignore_errors=True)
    os.makedirs(est_dir)
    cmd, mats, registered = external_command(method, series_path, reference_path, est_dir, corrected_path, sing_afni, sing_fsl)
    runner(cmd)
    if not os.path.exists(mats):
        raise FileNotFoundError(f'{method} did not write {mats}: see the log of the command')
    if registered != corrected_path:
        shutil.move(registered, corrected_path)
        registered = corrected_path
    ras = read_ras_matrices(method, mats, reference_path, series_path)
    itk = np.array([_LPS_RAS @ m @ _LPS_RAS for m in ras])            # RAS (nitransforms) -> LPS (ITK / ANTs)
    n_vol = ants.image_header_info(series_path)['dimensions'][3]
    if len(itk) != n_vol:
        raise ValueError(f'{method} gave {len(itk)} transforms for {n_vol} volumes')
    name, info = choose_convention(itk, series_path, registered, work_dir=est_dir)
    grid = ants.slice_image(ants.image_read(registered), axis=3, idx=0)
    centre = fov_centre(grid)
    for k, m in enumerate(itk):
        write_itk_matrix(m if name == 'as_is' else np.linalg.inv(m), opj(est_dir, f'img_{k:05d}_0GenericAffine.mat'), centre)
    for f in glob.glob(opj(est_dir, 'check_*.mat')):
        os.remove(f)
    return dict(info, convention=name, n_volumes=int(n_vol), method=method)
