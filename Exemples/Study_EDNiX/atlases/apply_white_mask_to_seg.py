#!/usr/bin/env python3
"""
apply_white_mask_to_seg.py: carry a hand-corrected <Species>_desc-White_mask.nii.gz of a library template into its
<Species>_seg-4FS_dseg.nii.gz, the file the surfaces (preFS.prepa_img) and the subjects (warped seg-4FS) are built from.
The White mask itself is not read by the surface reconstruction: it is derived from seg-4FS (preFS.msk_RS).

    python Exemples/Study_EDNiX/atlases/apply_white_mask_to_seg.py                      # Marmoset, V2 library
    python Exemples/Study_EDNiX/atlases/apply_white_mask_to_seg.py --species Marmoset --library <path> --plan

Rules (labels FreeSurfer, seg-4FS):
  - inside the mask, cortex (3/42), CSF (24) or background (0) -> white matter of the hemisphere of the nearest white
    matter voxel (2 or 41); other labels (subcortical, hippocampus, ventricles...) are kept;
  - outside the mask, white matter (2/41) -> label of the nearest non-white-matter voxel among cortex 3/42
    (white matter removed by the correction becomes cortex);
  - the eroded mask <Species>_desc-erod-White_mask.nii.gz is rewritten from the corrected mask (radius 1, as msk_RS).
The original seg-4FS is kept as <Species>_seg-4FS_desc-beforeWMfix_dseg.nii.gz (never overwritten).
"""
import argparse
import os
import shutil
import sys
from os.path import join as opj

import numpy as np
import nibabel as nib
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
import build_library_V2 as L  # noqa: E402

WM = (2, 41)
CORTEX = (3, 42)
TO_WM = (0, 3, 42, 24)


def relabel(seg, mask):
    """(new seg, n voxels -> WM, n voxels WM -> cortex)."""
    seg = seg.copy()
    is_wm = np.isin(seg, WM)
    # hemisphere of the nearest white matter voxel
    _, idx = ndimage.distance_transform_edt(~is_wm, return_indices=True)
    nearest_wm = seg[tuple(idx)]
    add = mask & ~is_wm & np.isin(seg, TO_WM)
    seg[add] = nearest_wm[add]
    remove = ~mask & np.isin(seg, WM)
    if remove.any():
        is_cx = np.isin(seg, CORTEX)
        _, idx = ndimage.distance_transform_edt(~is_cx, return_indices=True)
        seg[remove] = seg[tuple(idx)][remove]
    return seg, int(add.sum()), int(remove.sum())


def apply(species, library, plan=False):
    for sp, ref, tdir in L.find_templates(library, species, L.REFERENCES):
        vol = opj(tdir, ref, 'volumes')
        seg_f = opj(vol, 'labels', f'{sp}_seg-4FS_dseg.nii.gz')
        msk_f = opj(vol, 'masks', f'{sp}_desc-White_mask.nii.gz')
        if not (os.path.isfile(seg_f) and os.path.isfile(msk_f)):
            print(f'[{sp}] missing {seg_f if not os.path.isfile(seg_f) else msk_f}: skipped')
            continue
        backup = seg_f.replace('_seg-4FS_dseg', '_seg-4FS_desc-beforeWMfix_dseg')
        src = backup if os.path.isfile(backup) else seg_f          # rerunnable: always start from the original
        seg_img, msk_img = nib.load(src), nib.load(msk_f)
        if seg_img.shape != msk_img.shape or not np.allclose(seg_img.affine, msk_img.affine, atol=1e-4):
            sys.exit(f'[{sp}] White mask and seg-4FS are not on the same grid: {msk_f}')
        seg = np.asarray(seg_img.dataobj).round().astype(np.int32)
        mask = np.asarray(msk_img.dataobj) > 0.5
        new, n_add, n_rm = relabel(seg, mask)
        vox_mm3 = float(np.prod(seg_img.header.get_zooms()[:3]))
        print(f'[{sp}] {seg_f}\n   -> white matter: +{n_add} voxels ({n_add * vox_mm3:.1f} mm3), '
              f'-{n_rm} voxels ({n_rm * vox_mm3:.1f} mm3); WM total {np.isin(seg, WM).sum()} -> {np.isin(new, WM).sum()}')
        if plan:
            continue
        if not os.path.isfile(backup):
            shutil.copy2(seg_f, backup)
        hdr = seg_img.header.copy()
        hdr.set_data_dtype(np.int32)
        nib.save(nib.Nifti1Image(new, seg_img.affine, hdr), seg_f)
        ero = ndimage.binary_erosion(np.isin(new, WM), structure=ndimage.generate_binary_structure(3, 1))
        nib.save(nib.Nifti1Image(ero.astype(np.uint8), msk_img.affine), opj(vol, 'masks', f'{sp}_desc-erod-White_mask.nii.gz'))
        print(f'   written (original kept: {backup})')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--species', nargs='+', default=['Marmoset'])
    ap.add_argument('--library', default=L.LIBRARY_V2 if os.path.isdir(L.LIBRARY_V2) else L.LIBRARY_V1)
    ap.add_argument('--plan', action='store_true', help='report the changes, write nothing')
    a = ap.parse_args()
    apply(a.species, a.library, a.plan)
