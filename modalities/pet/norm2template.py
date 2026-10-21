"""
Step 9 of the PET pipeline: parametric images of the PET space exported to the AC-PC space of the
anatomy (pet/acpc) and to every template space (pet/templates/<name>), with their skull-stripped
version (-SS), through the chains of the registration record (one resampling per image).

Grid of the exported images (export_resolution):
    'pet'    the target space resampled at the PET voxel size (as the template of the functional
             pipeline is resampled at the functional resolution): no false precision, small files;
    'target' the grid of the target image (anatomy / template resolution).
The atlases and masks of the template spaces are written on the same grid (templates/<name>/labels,
masks), for group analyses.
"""
import os

import numpy as np
import ants

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd, getpath
from modalities.pet import petutils as pu
from modalities.pet import transfo

petinterp = 'linear'

PALETTE = (' MODE_USER_SCALE -thresholding THRESHOLD_TYPE_OFF THRESHOLD_TEST_SHOW_INSIDE 0 {v}'
           ' -pos-user 0 {v} -interpolate true -palette-name ROY-BIG -disp-pos true -disp-neg false -disp-zero false')


def target_grid(target, pet_ref, out_dir, name, export_resolution):
    """Reference image of the exported images in one target space (created once)."""
    if export_resolution == 'target':
        return target['image']
    grid = opj(out_dir, name + '_desc-grid-pet-resolution.nii.gz')
    if not opi(grid):
        os.makedirs(out_dir, exist_ok=True)
        img = ants.image_read(target['image'])
        vox = float(min(ants.image_read(pet_ref).spacing[:3]))
        if min(img.spacing[:3]) < vox:
            img = ants.resample_image(img, (vox, vox, vox), use_voxels=False, interp_type=0)
        ants.image_write(img, grid, ri=False)
    return grid


def _palette(img, maxval, diary_name, sing_wb):
    if sing_wb is None or maxval is None:
        return
    run_cmd.run(sing_wb + 'wb_command -volume-palette ' + img + PALETTE.format(v=maxval), diary_name)


def export(images, record, session, paths, bids_dir, export_resolution, diary_name, sing_wb,
           maxvals=None, target_files=None):
    """
    images: parametric images of orig/postprocessing (3D, or 4D when asked). For every target of the record
    (acpc, templates): <dir>/<base>_space-<name>_desc-<desc>_pet.nii.gz and its -SS version.
    target_files: {target name: dict(label_base, list_atlas, masks)} to also write the atlases
    (<label_base><atlas>_dseg.nii.gz, as in backtonative.get) and masks of a template space on the export grid.
    """
    maxvals = maxvals or {}
    out = []
    for name, t in record['targets'].items():
        if name == 'acpc':
            out_dir, lab_dir, msk_dir = paths['acpc'], None, None
        else:
            out_dir, lab_dir, msk_dir = getpath.emission_template_dirs(paths, name)
        os.makedirs(out_dir, exist_ok=True)
        grid = target_grid(t, record['reference'], out_dir, name, export_resolution)
        mask = None
        if t.get('mask') and opi(t['mask']):
            mask = ants.threshold_image(ants.resample_image_to_target(ants.image_read(t['mask']), ants.image_read(grid),
                                                                      interp_type='nearestNeighbor'), 0.5, 1e9, 1, 0)
        for img in images:
            base, desc = _split(img)
            dst = pu.out_name(out_dir, base, desc=desc, space=name)
            res = transfo.apply(grid, img, t, 'to', petinterp)
            ants.image_write(res, dst, ri=False)
            pu.write_json(dst, [img, grid], f'PET image in the {name} space (ANTsPy, linear, one resampling).',
                          Transforms=t['to'])
            _palette(dst, maxvals.get(desc), diary_name, sing_wb)
            out.append(dst)
            if mask is not None and res.dimension == 3:
                ss = pu.out_name(out_dir, base, desc=desc + '-SS', space=name)
                ants.image_write(res * mask, ss, ri=False)
                pu.write_json(ss, [dst, t['mask']], 'Skull-stripped (brain mask of the target space).')
                _palette(ss, maxvals.get(desc), diary_name, sing_wb)
                out.append(ss)
        if lab_dir and target_files and name in target_files:
            _template_labels(target_files[name], grid, lab_dir, msk_dir, diary_name)
    return out


def _split(img):
    """(base, desc) of an image of orig/postprocessing."""
    b = pu.strip_ext(img).split('_')
    desc = next((p[5:] for p in b if p.startswith('desc-')), 'pet')
    return pu.base_name(img), desc


def _template_labels(tf, grid, lab_dir, msk_dir, diary_name):
    os.makedirs(lab_dir, exist_ok=True)
    os.makedirs(msk_dir, exist_ok=True)
    g = ants.image_read(grid)
    for atlas in tf.get('list_atlas', []):
        src = tf['label_base'] + atlas + '_dseg.nii.gz'
        dst = opj(lab_dir, opb(src))
        if opi(src) and not opi(dst):
            a = ants.image_read(src)
            if a.dimension == 4:
                vols = [ants.resample_image_to_target(v, g, interp_type='genericLabel') for v in ants.ndimage_to_list(a)]
                ants.image_write(ants.list_to_ndimage(_stack_like(a, g), vols), dst, ri=False)
            else:
                ants.image_write(ants.resample_image_to_target(a, g, interp_type='genericLabel'), dst, ri=False)
            pu.write_json(dst, [src, grid], 'Atlas on the export grid (genericLabel).')
    for src in tf.get('masks', []):
        dst = opj(msk_dir, opb(src))
        if src and opi(src) and not opi(dst):
            m = ants.resample_image_to_target(ants.image_read(src), g, interp_type='linear')
            ants.image_write(ants.threshold_image(m, 0.5, 1e9, 1, 0), dst, ri=False)
            pu.write_json(dst, [src, grid], 'Mask on the export grid (linear, threshold 0.5).')


def _stack_like(a4, g3):
    """4D image with the spatial geometry of g3 and the number of volumes of a4 (template for list_to_ndimage)."""
    n = a4.shape[3]
    arr = np.zeros(tuple(g3.shape) + (n,), dtype=np.float32)
    sp = tuple(g3.spacing) + (1.0,)
    org = tuple(g3.origin) + (0.0,)
    d = np.eye(4)
    d[:3, :3] = g3.direction
    return ants.from_numpy(arr, origin=org, spacing=sp, direction=d)
