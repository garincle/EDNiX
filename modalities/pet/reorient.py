"""
Step 1 of the PET pipeline (also used by the CT pipeline): orientation of the header.

Independent image (no other modality in the same scanner space, or one registered later):
the image is handled as the anatomical images (modalities/anat/loop1/_1_correct_orient) and the
functional images (modalities/fMRI/_2b_fix_orient) of EDNiX:
    1. an oblique image is resampled on a cardinal grid (3dWarp -deoblique, world coordinates kept),
       or its obliquity is dropped from the header (deoblique='header', 3drefit -deoblique);
    2. the axes are relabelled (3drefit -orient) with the `orientation` code when it is given,
       otherwise with get_orientation.getreal(humanPosition, animalPosition, current code).

Image in the scanner space of a reference already reoriented (hybrid PET-MR: the anatomy; PET/CT: the
CT), `follow` given: after the obliquity correction, the image receives the world transform W that the
relabelling applied to the reference, W = A_after . A_before^-1 (voxel-to-world matrices of the
reference after and before its relabelling). 3drefit -orient keeps the magnitude of each origin and
sets its sign from the new orientation code (3drefit.c: xxorg = (ORIENT_sign == '+') ? -|xxorg| :
|xxorg|): two images of one scanner space with different fields of view (PET and MR or CT grids) are
relabelled consistently only when their first voxels lie on the same side of the coordinate origin.
Copying W makes the PET follow the reference exactly, whatever the grids.

With header_wrong=True (reconstruction that wrote a wrong orientation in the header, PET engineer's
case), the orientation code of a raw reference image acquired in the same space is copied first.
Every reoriented image records its own world transform in its sidecar ('WorldTransform'), read when it
is the reference of another modality.
"""
import os
import shutil

import numpy as np
import nibabel as nib

opj = os.path.join
opd = os.path.dirname
opb = os.path.basename
opi = os.path.isfile

from Tools import run_cmd, get_orientation
from modalities.pet import petutils as pu

# linear interpolation: no over/undershoot of the activity concentrations (the engineer used wsinc5)
deob_interp = '-linear'


def _afni_info(opt, img, sing_afni, diary_name):
    out, _ = run_cmd.get(sing_afni + '3dinfo ' + opt + ' ' + img, diary_name)
    lines = [l for l in out.decode('utf-8', errors='replace').split('\n') if l.strip()]
    if not lines:
        raise RuntimeError(f'3dinfo {opt} gave no output for {img}')
    return lines[-1].strip()


def no_reorientation(animalPosition):
    return animalPosition in ('', 'na', 'humanlike', None)


def deoblique_apply(input, output, deoblique, diary_name, sing_afni):
    """Obliquity correction (see module docstring). Returns the descriptions of what was done."""
    if not opi(input):
        raise FileNotFoundError(f"Input file not found: {input}")
    obli = _afni_info('-is_oblique', input, sing_afni, diary_name)
    if obli == '1' and deoblique == 'WARP':
        run_cmd.run(sing_afni + '3dWarp -overwrite -deoblique ' + deob_interp + ' -prefix ' + output + ' ' + input,
                    diary_name)
        return ['obliquity corrected by resampling (3dWarp -deoblique -linear)']
    shutil.copyfile(input, output)
    if obli == '1':
        run_cmd.run(sing_afni + '3drefit -deoblique ' + output, diary_name)
        return ['obliquity removed from the header (3drefit -deoblique)']
    return []


def world_map(before, after):
    """World transform between two headers of the same voxel array: W = A_after . A_before^-1 (RAS, mm)."""
    b, a = nib.load(before), nib.load(after)
    if tuple(b.shape[:3]) != tuple(a.shape[:3]):
        raise ValueError(f'{after} ({a.shape[:3]}) and {before} ({b.shape[:3]}) do not share their voxel array: the '
                         f'world transform of the reorientation cannot be derived')
    return a.affine @ np.linalg.inv(b.affine)


def apply_world_map(input, W, output):
    """Header of `input` moved by the world transform W (data untouched): A_new = W . A."""
    im = nib.load(input)
    A = np.asarray(W, dtype=float) @ im.affine
    hdr = im.header.copy()
    hdr.set_data_dtype(np.float32)          # scaled values kept exactly (int16 x slope fits in float32)
    out = nib.Nifti1Image(np.asarray(im.dataobj, dtype=np.float32), A, hdr)
    out.header['scl_slope'] = 1
    out.header['scl_inter'] = 0
    code = int(hdr['sform_code']) or int(hdr['qform_code']) or 1
    out.set_sform(A, code=code)
    out.set_qform(A, code=code)
    nib.save(out, output)
    return output


def reference_world_map(ref_raw, ref_reoriented, work_dir, diary_name, sing_afni):
    """
    World transform applied by the reorientation of a reference image (raw -> reoriented): read from its
    sidecar when it was reoriented by this module (CT), else derived from the two images (anatomical
    pipeline: copy of the raw image `_space-raw_desc-template_`, and `_space-raw_desc-reorient_`), the
    deobliqued grid being recomputed with 3dWarp when the raw image is oblique (world coordinates kept).
    """
    side = pu.read_json(pu.sidecar(ref_reoriented))
    if side.get('WorldTransform') is not None:
        return np.asarray(side['WorldTransform'], dtype=float)
    for f in (ref_raw, ref_reoriented):
        if not f or not opi(f):
            raise FileNotFoundError(f'reference image of the same scanner space not found: {f!r}')
    before = ref_raw
    if _afni_info('-is_oblique', ref_raw, sing_afni, diary_name) == '1':
        os.makedirs(work_dir, exist_ok=True)
        before = opj(work_dir, 'reference_deoblique_tmp.nii.gz')
        run_cmd.run(sing_afni + '3dWarp -overwrite -deoblique -prefix ' + before + ' ' + ref_raw, diary_name)
    W = world_map(before, ref_reoriented)
    if before != ref_raw:
        os.remove(before)
    return W


def reorient_apply(input, output, humanPosition, animalPosition, orientation, deoblique, diary_name, sing_afni):
    """Deoblique + relabel the axes (see module docstring). Returns (description, world transform W)."""
    desc = deoblique_apply(input, output, deoblique, diary_name, sing_afni)
    A_before = nib.load(output).affine
    orient_img = _afni_info('-orient', output, sing_afni, diary_name)
    if orientation:
        new_orient = orientation
    elif not no_reorientation(animalPosition):
        new_orient = get_orientation.getreal(humanPosition, animalPosition, orient_img)
    else:
        new_orient = orient_img
    if new_orient != orient_img:
        run_cmd.run(sing_afni + '3drefit -orient ' + new_orient + ' ' + output, diary_name)
        desc.append(f'axes relabelled {orient_img} -> {new_orient} (3drefit -orient; humanPosition={humanPosition}, '
                    f'animalPosition={animalPosition}, orientation={orientation!r})')
    run_cmd.msg(f'INFO: orientation {orient_img} -> {new_orient}', diary_name, 'OKGREEN')
    W = nib.load(output).affine @ np.linalg.inv(A_before)
    return '; '.join(desc) or 'no change of orientation needed (copy)', W


def copy_orientation(input, ref, output, diary_name, sing_afni):
    """Header written with a wrong orientation: copy the orientation code of `ref` (same scanner space)."""
    shutil.copyfile(input, output)
    orient_ref = _afni_info('-orient', ref, sing_afni, diary_name)
    run_cmd.run(sing_afni + '3drefit -orient ' + orient_ref + ' ' + output, diary_name)
    return f'orientation code copied from {ref} ({orient_ref}) (3drefit -orient)'


def doit(new_dir, img, humanPosition, animalPosition, orientation, Ref_native_raw, Need_check, header_wrong,
         timing, diary_name, sing_afni, sing_fs, suffix='pet', deoblique='WARP', sif='', header_ref='',
         follow=None):
    """
    Reorient one raw image into new_dir: <base>_desc-reorient_<suffix>.nii.gz (+ its duration-weighted
    mean <base>_desc-reorient-mean_<suffix>.nii.gz for a dynamic image, + the SIF of the raw image).
    Ref_native_raw: image displayed under the result for the visual check; may be '' when there is none.
    header_ref: raw image of the same scanner space whose orientation code is copied when header_wrong=True.
    sif: SIF of the raw image (Turku format, step 1), copied next to the reoriented image.
    follow: dict(W=4x4 world transform, source=description) of the reference of the same scanner space
    (reference_world_map); the positions and the orientation code are then not used.
    Returns (reoriented image, image to look at).
    """
    nl = 'REORIENT THE FILE ' + opb(img)
    run_cmd.msg(nl, diary_name, 'OKGREEN')
    if not opi(img):
        raise FileNotFoundError(f"Input file not found: {img}")
    os.makedirs(new_dir, exist_ok=True)

    base = pu.base_name(img, suffix)
    step3 = pu.out_name(new_dir, base, desc='reorient', suffix=suffix)
    step2 = pu.out_name(new_dir, base, desc='reorient-mean', suffix=suffix)

    if header_wrong:
        # the orientation code of the reference (same scanner space, raw image) is copied, then the image is
        # reoriented like the reference was (the engineer copied the code but skipped the reorientation)
        if not header_ref or not opi(header_ref):
            raise FileNotFoundError('header_wrong=True needs the raw image of a modality acquired in the same space '
                                    f'(anatomy or CT); not found: {header_ref!r}')
        tmp = step3.replace('.nii.gz', '_tmp.nii.gz')
        desc0 = copy_orientation(img, header_ref, tmp, diary_name, sing_afni)
        desc, W = reorient_apply(tmp, step3, humanPosition, animalPosition, orientation, deoblique, diary_name,
                                 sing_afni)
        desc = desc0 + '; ' + desc
        os.remove(tmp)
    elif follow is not None:
        tmp = step3.replace('.nii.gz', '_tmp.nii.gz')
        d = deoblique_apply(img, tmp, deoblique, diary_name, sing_afni)
        W = np.asarray(follow['W'], dtype=float)
        apply_world_map(tmp, W, step3)
        os.remove(tmp)
        d.append('world transform of the reorientation of the reference of the same scanner space applied to the '
                 'header (' + str(follow.get('source', '')) + ')')
        desc = '; '.join(d)
        run_cmd.msg('INFO: same scanner space: the reorientation of ' + str(follow.get('source', 'the reference')) +
                    ' is applied (positions and orientation code of this modality not used)', diary_name, 'OKGREEN')
    else:
        desc, W = reorient_apply(img, step3, humanPosition, animalPosition, orientation, deoblique,
                                 diary_name, sing_afni)
    im3 = nib.load(step3)
    if len(im3.shape) == 4 and im3.shape[3] == 1:
        # static image stored with one frame: written as 3D (one code path for the static images)
        pu.save_like(np.asarray(im3.dataobj, dtype=np.float32)[..., 0], im3, step3)
        desc += '; single frame written as a 3D image'
    pu.write_json(step3, img, desc, WorldTransform=np.asarray(W).tolist())

    if sif and opi(sif):
        shutil.copyfile(sif, step3.replace('.nii.gz', '.sif'))

    to_see = step3
    if pu.is_dynamic(step3):
        pu.write_mean(step3, step2, timing, None, 'Time average of the reoriented frames (weights = frame durations).')
        to_see = step2

    if Need_check:
        run_cmd.msg('INFO: check the orientation (close freeview to continue)', diary_name, 'OKGREEN')
        ref = (' ' + Ref_native_raw + ':colormap=grayscale:opacity=1') if Ref_native_raw and opi(Ref_native_raw) else ''
        run_cmd.do(sing_fs + 'freeview -v' + ref + ' ' + to_see + ':colormap=heat:opacity=0.5:visible=1', diary_name)
    return step3, to_see
