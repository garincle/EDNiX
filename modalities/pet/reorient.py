"""
Step 1 of the PET pipeline (also used by the CT pipeline): orientation of the header.

The image is handled exactly as the anatomical images (modalities/anat/loop1/_1_correct_orient)
and the functional images (modalities/fMRI/_2b_fix_orient) of EDNiX, so that images acquired in
the same scanner space stay in the same space after the correction:
    1. an oblique image is resampled on a cardinal grid (3dWarp -deoblique), or its obliquity is
       dropped from the header (deoblique='header', 3drefit -deoblique, like doWARPonfunc='header');
    2. the axes are relabelled (3drefit -orient) with the `orientation` code when it is given,
       otherwise with get_orientation.getreal(humanPosition, animalPosition, current code).
With header_wrong=True (reconstruction that wrote a wrong orientation in the header, PET
engineer's case), the orientation code of a reference image acquired in the same space (anatomy or
CT) is copied instead.
"""
import os
import shutil

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


def reorient_apply(input, output, humanPosition, animalPosition, orientation, deoblique, diary_name, sing_afni):
    """Deoblique + relabel the axes (see module docstring). Returns the description written in the sidecar."""
    if not opi(input):
        raise FileNotFoundError(f"Input file not found: {input}")

    obli = _afni_info('-is_oblique', input, sing_afni, diary_name)
    desc = []
    if obli == '1' and deoblique == 'WARP':
        run_cmd.run(sing_afni + '3dWarp -overwrite -deoblique ' + deob_interp + ' -prefix ' + output + ' ' + input,
                    diary_name)
        desc.append('obliquity corrected by resampling (3dWarp -deoblique -linear)')
    else:
        shutil.copyfile(input, output)
        if obli == '1':
            run_cmd.run(sing_afni + '3drefit -deoblique ' + output, diary_name)
            desc.append('obliquity removed from the header (3drefit -deoblique)')

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
    return '; '.join(desc) or 'no change of orientation needed (copy)'


def copy_orientation(input, ref, output, diary_name, sing_afni):
    """Header written with a wrong orientation: copy the orientation code of `ref` (same scanner space)."""
    shutil.copyfile(input, output)
    orient_ref = _afni_info('-orient', ref, sing_afni, diary_name)
    run_cmd.run(sing_afni + '3drefit -orient ' + orient_ref + ' ' + output, diary_name)
    return f'orientation code copied from {ref} ({orient_ref}) (3drefit -orient)'


def doit(new_dir, img, humanPosition, animalPosition, orientation, Ref_native_raw, Need_check, header_wrong,
         timing, diary_name, sing_afni, sing_fs, suffix='pet', deoblique='WARP', sif='', header_ref=''):
    """
    Reorient one raw image into new_dir: <base>_desc-reorient_<suffix>.nii.gz (+ its duration-weighted
    mean <base>_desc-reorient-mean_<suffix>.nii.gz for a dynamic image, + the SIF of the raw image).
    Ref_native_raw: image displayed under the result for the visual check; may be '' when there is none.
    header_ref: raw image of the same scanner space whose orientation code is copied when header_wrong=True.
    sif: SIF of the raw image (Turku format, step 1), copied next to the reoriented image.
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
        desc = desc0 + '; ' + reorient_apply(tmp, step3, humanPosition, animalPosition, orientation, deoblique,
                                             diary_name, sing_afni)
        os.remove(tmp)
    else:
        desc = reorient_apply(img, step3, humanPosition, animalPosition, orientation, deoblique,
                              diary_name, sing_afni)
    pu.write_json(step3, img, desc)

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
