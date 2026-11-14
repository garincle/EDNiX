"""
Brain mask drawn by hand on the raw PET of a session (same convention as the functional pipeline,
func/raw/masks/<ID>_final_mask.nii.gz):

    <BIDS>/sub-<ID>/ses-<S>/pet/raw/masks/<session>_desc-manual_mask.nii.gz
    (<session> = name of the runs without run / desc entities, e.g. sub-01_ses-1_trc-FDG)

The mask is on the grid of the raw (BIDS) PET, before the reorientation of step 1. It is brought into the
reoriented space (pet/orig/masks/<session>_space-orig_desc-manual_mask.nii.gz):
  - reorientation that only rewrote the header (3drefit -orient, header copy, world transform of a reference):
    the voxel array of the mask is kept and receives the header of the reoriented PET (exact);
  - reorientation that resampled the PET (oblique image, deoblique='WARP'): the world transform of the
    reorientation (sidecar 'WorldTransform') is applied to the header of the mask, then the mask is resampled
    on the reoriented grid (nearest neighbour).

Uses (modalities/pet/loop2.py):
  - brain mask of the session in the PET space (bratio, skull-stripped dynamic image, SUV-window ROI): it
    replaces the brain of the atlas when both exist (scenario 3v2 has no atlas, hence no brain mask otherwise);
  - crop_for_MoCo='mask' (or 'auto' without crop of the study template): head box of steps 4-6 = bounding box of
    the dilated mask, instead of a box drawn in freeview;
  - pet_coreg_mask='manual': dilated mask as moving mask of the PET -> anatomy / CT registration;
  - PET study template (step 2) and scenario 3v2 (step 6): when every session has a manual mask, the template is built
    from the brain of the session means (masked_image: mean x dilated mask, cropped to its bounding box) and each
    motion-corrected session mean is registered to it the same way (the extracerebral uptake, the bed and the FOV
    edges do not drive the registrations).
"""
import glob
import os

import numpy as np
import nibabel as nib
from nibabel.processing import resample_from_to
from scipy import ndimage

opj = os.path.join
opb = os.path.basename
opi = os.path.isfile

from Tools import run_cmd
from modalities.pet import petutils as pu


def raw_file(paths, session):
    """Manual mask of the session in pet/raw/masks ('' when there is none). A single *_desc-manual_mask file of
    another name is accepted (several: ambiguous, not used)."""
    d = opj(paths['raw'], 'masks')
    f = opj(d, session + '_desc-manual_mask.nii.gz')
    if opi(f):
        return f
    found = sorted(glob.glob(opj(d, '*_desc-manual_mask.nii.gz')))
    return found[0] if len(found) == 1 else ''


def orig_file(paths, session):
    return opj(paths['orig_masks'], session + '_space-orig_desc-manual_mask.nii.gz')


def _first3d(img):
    """3D nibabel image (first volume of a 4D image)."""
    if len(img.shape) > 3:
        return nib.Nifti1Image(np.asarray(img.dataobj[..., 0], dtype=np.float32), img.affine)
    return img


def _same_grid(a, b, atol=1e-3):
    return tuple(a.shape[:3]) == tuple(b.shape[:3]) and np.allclose(a.affine, b.affine, atol=atol)


def _binary(img):
    d = np.asarray(img.dataobj, dtype=np.float32)
    if d.ndim > 3:
        d = d[..., 0]
    return d > 0.5


def _write(mask_bool, affine, out, sources, description):
    os.makedirs(os.path.dirname(out), exist_ok=True)
    img = nib.Nifti1Image(mask_bool.astype(np.uint8), affine)
    img.set_sform(affine, code=1)
    img.set_qform(affine, code=1)
    nib.save(img, out)
    pu.write_json(out, sources, description)
    return out


def _nearest(mask_bool, affine, ref3d):
    """Mask resampled on the grid of ref3d (nearest neighbour, world coordinates)."""
    src = nib.Nifti1Image(mask_bool.astype(np.float32), affine)
    return np.asarray(resample_from_to(src, ref3d, order=0).dataobj) > 0.5


def to_orig(raw_mask, raw_img, reoriented, out, diary_file):
    """
    Manual mask of the raw PET `raw_img` brought into the space of `reoriented` (step 1 output of the same run).
    Returns `out`.
    """
    m_img = nib.load(raw_mask)
    raw3 = _first3d(nib.load(raw_img))
    reo = nib.load(reoriented)
    reo3 = _first3d(reo)
    m = _binary(m_img)
    if not m.any():
        raise ValueError(run_cmd.error(f'ERROR: manual mask {raw_mask} is empty', diary_file))
    note = ''
    if not _same_grid(m_img, raw3):
        # drawn on another grid of the same scanner space: resampled on the raw PET first
        m = _nearest(m, m_img.affine, raw3)
        note = 'mask resampled on the grid of the raw PET (nearest neighbour); '
        run_cmd.msg(f'WARNING: {opb(raw_mask)} is not on the grid of {opb(raw_img)}: resampled (nearest neighbour)',
                    diary_file, 'WARNING')
    side = pu.read_json(pu.sidecar(reoriented))
    desc = str(side.get('Description', ''))
    if tuple(reo3.shape[:3]) == tuple(raw3.shape[:3]) and '3dWarp' not in desc:
        out_mask, how = m, 'voxel array kept, header of the reoriented PET (reorientation of the header only)'
    else:
        W = side.get('WorldTransform')
        if W is None:
            raise ValueError(run_cmd.error(f'ERROR: {reoriented} was resampled by step 1 but its sidecar has no '
                                           'WorldTransform: run step 1 again', diary_file))
        A = np.asarray(W, dtype=float) @ raw3.affine
        out_mask = _nearest(m, A, reo3)
        how = 'world transform of step 1 applied, then nearest-neighbour resampling on the reoriented grid'
    if not out_mask.any():
        raise ValueError(run_cmd.error(f'ERROR: manual mask {raw_mask} is empty after the reorientation', diary_file))
    _write(out_mask, reo3.affine, out, [raw_mask, reoriented], 'Brain mask drawn by hand on the raw PET, in the '
           'reoriented space (' + note + how + ').')
    run_cmd.msg(f'INFO: manual brain mask {opb(raw_mask)} -> {opb(out)} ({how})', diary_file, 'OKGREEN')
    return out


def prepare(paths, session, raw_img, reoriented, diary_file):
    """
    Manual mask of the session in the reoriented space ('' without manual mask). Computed again when the raw mask
    or the reoriented image is newer than it (step 1 run again, mask redrawn).
    """
    raw_mask = raw_file(paths, session)
    if not raw_mask:
        return ''
    out = orig_file(paths, session)
    if opi(out) and os.path.getmtime(out) >= max(os.path.getmtime(raw_mask), os.path.getmtime(reoriented)):
        return out
    return to_orig(raw_mask, raw_img, reoriented, out, diary_file)


def on_grid(mask, ref, out):
    """Mask on the grid of the image `ref` (nearest neighbour); `mask` itself when it is already on it."""
    m_img, r3 = nib.load(mask), _first3d(nib.load(ref))
    if _same_grid(m_img, r3):
        return mask
    _write(_nearest(_binary(m_img), m_img.affine, r3), r3.affine, out, [mask, ref],
           'Manual brain mask on the grid of the PET (nearest neighbour).')
    return out


def dilated(mask, dilate_mm, out):
    """Mask dilated by dilate_mm (ball, voxel size of the mask)."""
    img = nib.load(mask)
    m = _binary(img)
    m = _ball_dilation(m, np.asarray(img.header.get_zooms()[:3], dtype=float), dilate_mm)
    return _write(m, img.affine, out, [mask], f'Manual brain mask dilated by {dilate_mm} mm.')


def box(mask, image, out, margin_mm=10.0):
    """Head box of steps 4-6 (crop_for_MoCo='mask'): bounding box of the manual mask enlarged by margin_mm, on the
    grid of `image` (mean image of the first run of the session)."""
    r3 = _first3d(nib.load(image))
    m = _binary(nib.load(on_grid(mask, image, out.replace('.nii.gz', '_tmp.nii.gz'))))
    tmp = out.replace('.nii.gz', '_tmp.nii.gz')
    if opi(tmp):
        os.remove(tmp)
        js = pu.sidecar(tmp)
        if opi(js):
            os.remove(js)
    idx = np.argwhere(m)
    zooms = np.asarray(r3.header.get_zooms()[:3], dtype=float)
    pad = np.ceil(margin_mm / zooms).astype(int)
    lo = np.maximum(idx.min(0) - pad, 0)
    hi = np.minimum(idx.max(0) + pad + 1, np.asarray(m.shape))
    b = np.zeros(m.shape, bool)
    b[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] = True
    return _write(b, r3.affine, out, [mask, image], f'Head box of the motion correction: bounding box of the manual '
                                                    f'brain mask enlarged by {margin_mm} mm.')


def _ball_dilation(m, zooms, dilate_mm):
    if dilate_mm <= 0:
        return m
    r = np.ceil(dilate_mm / zooms).astype(int)
    g = np.mgrid[-r[0]:r[0] + 1, -r[1]:r[1] + 1, -r[2]:r[2] + 1].astype(float)
    ball = ((g[0] * zooms[0]) ** 2 + (g[1] * zooms[1]) ** 2 + (g[2] * zooms[2]) ** 2) <= dilate_mm ** 2
    return ndimage.binary_dilation(m, structure=ball)


def masked_image(mask, image, out, dilate_mm=3.0, margin_mm=5.0):
    """
    Brain of `image` (3D, or first volume of a 4D image): image x (manual mask on its grid, dilated by dilate_mm),
    cropped to the bounding box of the dilated mask enlarged by margin_mm. Same world coordinates as `image` (the
    affine of the crop is shifted accordingly), so a registration estimated on it is valid for `image`.
    """
    img = _first3d(nib.load(image))
    d = np.asarray(img.dataobj, dtype=np.float32)
    m_img = nib.load(mask)
    m = _binary(m_img) if _same_grid(m_img, img) else _nearest(_binary(m_img), m_img.affine, img)
    if not m.any():
        raise ValueError(f'manual mask {mask} does not overlap {image}')
    zooms = np.asarray(img.header.get_zooms()[:3], dtype=float)
    m = _ball_dilation(m, zooms, dilate_mm)
    idx = np.argwhere(m)
    pad = np.ceil(margin_mm / zooms).astype(int)
    lo = np.maximum(idx.min(0) - pad, 0)
    hi = np.minimum(idx.max(0) + pad + 1, np.asarray(m.shape))
    sl = tuple(slice(a, b) for a, b in zip(lo, hi))
    data = np.where(m, d, 0.0)[sl].astype(np.float32)
    aff = img.affine.copy()
    aff[:3, 3] = (img.affine @ np.r_[lo, 1.0])[:3]
    os.makedirs(os.path.dirname(out), exist_ok=True)
    o = nib.Nifti1Image(data, aff)
    o.set_sform(aff, code=1)
    o.set_qform(aff, code=1)
    nib.save(o, out)
    pu.write_json(out, [image, mask], f'Brain of the image: manual brain mask dilated by {dilate_mm} mm, cropped to its '
                                      f'bounding box + {margin_mm} mm (same world coordinates).',
                  BrainMasked=True)
    return out
