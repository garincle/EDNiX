"""
Susceptibility distortion correction (SDC) of the functional time series.

FUGUE unwarps the run mean and saves its voxel shift map (--saveshift, in
voxels along the phase-encoding axis). This module turns that shift map into
an ANTs displacement field (physical coordinates, mm) so that the correction
can be inserted in the one-shot resampling chain of every volume:

    point in reference space -> run_inRef -> SDC field -> motion_j -> volume j

The sign convention of the shift map is not assumed: it is determined by
reproducing FUGUE's own unwarped image with both signs and keeping the one
that matches (determine_sign). If neither matches, SDC is not applied.
"""
import json
import os
import tempfile
import numpy as np
import ants

AXES = {'x': 0, 'y': 1, 'z': 2}


def pe_axis(correction_direction):
    """Voxel axis (0, 1, 2) of an FSL unwarp direction such as 'y-'."""
    d = str(correction_direction).strip().lower()
    if not d or d[0] not in AXES:
        raise ValueError(f"correction_direction must be x, x-, y, y-, z or z- (got {correction_direction!r})")
    return AXES[d[0]]


def shift_to_field(shift_img, axis, sign):
    """
    ANTs displacement field from a voxel shift map: u(x) = sign * shift(x) * e_axis,
    with e_axis the physical vector (mm, ANTs LPS frame) of one voxel step along `axis`.
    `shift_img` is an ants image (or a path).
    """
    if isinstance(shift_img, str):
        shift_img = ants.image_read(shift_img)
    shift = shift_img.numpy().astype(np.float32)
    if shift.ndim == 4:
        shift = shift[..., 0]
    step = np.asarray(shift_img.direction)[:3, axis] * float(shift_img.spacing[axis])
    vec = (sign * shift)[..., None] * step[None, None, None, :].astype(np.float32)
    return ants.from_numpy(vec, origin=shift_img.origin[:3], spacing=shift_img.spacing[:3],
                           direction=np.asarray(shift_img.direction)[:3, :3], has_components=True)


def write_field(field, out_path):
    ants.image_write(field, out_path)
    return out_path


def _corr(a, b, mask):
    a, b = a[mask], b[mask]
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def determine_sign(distorted, fugue_unwarped, shift, axis, min_shift_vox=0.1, min_gain=0.002, min_corr=0.95):
    """
    Sign of the shift map such that warping `distorted` with the field reproduces
    FUGUE's unwarped image. Returns dict(sign, corr_plus, corr_minus, corr_none, max_shift, status)
    with status 'ok', 'negligible' (shift below min_shift_vox everywhere) or 'ambiguous'.

    The right sign moves the image towards the FUGUE output and the wrong one away
    from it: 'ok' when the best correlation is >= min_corr and either beats the other
    by min_gain, or the two signs fall on both sides of the uncorrected correlation
    (small distortions, where both effects are small but opposite).
    """
    dist = ants.image_read(distorted) if isinstance(distorted, str) else distorted
    ref = ants.image_read(fugue_unwarped) if isinstance(fugue_unwarped, str) else fugue_unwarped
    sh = ants.image_read(shift) if isinstance(shift, str) else shift
    max_shift = float(np.nanmax(np.abs(sh.numpy())))
    out = dict(sign=None, corr_plus=np.nan, corr_minus=np.nan, corr_none=np.nan, max_shift=max_shift, status='ambiguous')
    if max_shift < min_shift_vox:
        out['status'] = 'negligible'
        return out
    r = ref.numpy()
    pos = r[r > 0]
    mask = r > (0.1 * np.percentile(pos, 99)) if pos.size else r > 0
    out['corr_none'] = _corr(dist.numpy(), r, mask)
    fd, tmp = tempfile.mkstemp(suffix='_sdc_test_field.nii.gz')   # unique: runs may be processed in parallel
    os.close(fd)
    try:
        for sign, key in ((1, 'corr_plus'), (-1, 'corr_minus')):
            write_field(shift_to_field(sh, axis, sign), tmp)
            moved = ants.apply_transforms(fixed=dist, moving=dist, transformlist=[tmp], interpolator='linear')
            out[key] = _corr(moved.numpy(), r, mask)
    finally:
        os.remove(tmp)
    best = 1 if out['corr_plus'] >= out['corr_minus'] else -1
    hi, lo = max(out['corr_plus'], out['corr_minus']), min(out['corr_plus'], out['corr_minus'])
    opposite = np.isfinite(out['corr_none']) and lo < out['corr_none'] < hi
    if np.isfinite(hi) and hi >= min_corr and (hi - lo >= min_gain or opposite):
        out.update(sign=best, status='ok')
    return out


def run_field(matrices_dir, root_RS, root_RS_ref, recordings):
    """
    SDC field to apply to a run, or None.
      'very_old'          : no field map, no correction
      'new'               : one field map per run -> only the run's own field (a run whose
                            field was rejected is NOT corrected with another run's field)
      'old', '2_mapdir'   : one shared field map -> the reference run's field for every run
    """
    if recordings == 'very_old':
        return None
    own = os.path.join(matrices_dir, root_RS + '_space-func_desc-sdc_1Warp.nii.gz')
    if os.path.exists(own):
        return own
    if recordings in ('old', '2_mapdir'):
        ref = os.path.join(matrices_dir, root_RS_ref + '_space-func_desc-sdc_1Warp.nii.gz')
        return ref if os.path.exists(ref) else None
    return None


def save_report(path, info, sources):
    with open(path, 'w') as fh:
        json.dump(dict(info, Sources=sources,
                       Description='SDC displacement field (mm) from the FUGUE shift map; sign checked '
                                   'against the FUGUE-unwarped mean.'), fh, indent=2, default=float)


def pe_restriction(pe_code, direction):
    """
    Weights of restrict_transformation for a non-linear registration of an EPI image:
    1 along the physical axis closest to its phase-encoding voxel axis (EPI distortions
    happen along it), 0.1 on the others. ANTs weights the components of the displacement
    in physical space, hence the use of the image direction matrix.
    pe_code: 'i', 'j', 'k' (BIDS PhaseEncodingDirection, sign ignored) or 'x', 'y', 'z'
    (FSL unwarp direction). Anything else gives no restriction (1, 1, 1).
    """
    axis = {'i': 0, 'x': 0, 'j': 1, 'y': 1, 'k': 2, 'z': 2}.get(str(pe_code).strip().lower()[:1])
    if axis is None:
        return (1, 1, 1)
    phys = int(np.argmax(np.abs(np.asarray(direction, dtype=float)[:3, axis])))
    return tuple(1 if a == phys else 0.1 for a in range(3))

