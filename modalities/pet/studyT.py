"""
Study template of a PET (or CT) modality (steps 2 and 3 of the PET pipeline).

- make_template: template built from the mean image of every session (ants.build_template), then
  every session registered to it (<session>_to_studyTemplate_* in orig/matrices). It is used to crop the
  images before the motion correction (the crop drawn once on the template is brought back to every
  session) and, without anatomy and without modality template (scenarios 2v2 / 3v2 of pet4.svg), as the
  common space of the study.
- manual_crop: crop drawn by hand on the template with freeview (kept interactive, PET engineer).
- crop_mask: box of the cropped template brought into the space of one session.
"""
import os
import shutil

import numpy as np
import nibabel as nib
import ants

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd
from modalities.pet import petutils as pu
from modalities.pet import transfo


def template_name(study_dir, tracer, suffix='pet'):
    tag = ('trc-' + tracer + '_') if tracer else ''
    return opj(study_dir, 'studyTemplate_' + tag + suffix + '.nii.gz')


def crop_name(template):
    return template.replace('.nii.gz', '').rsplit('_', 1)[0] + '_desc-crop_' + template.rsplit('_', 1)[1]


def _normalised(path):
    """Mean image scaled by its 99th percentile (different doses must not weigh differently in the template)."""
    img = ants.image_read(path)
    a = img.numpy()
    p = np.percentile(a[np.isfinite(a) & (a > 0)], 99) if np.any(a > 0) else 1.0
    return img / float(p if p > 0 else 1.0)


def is_brain_masked(template):
    """True when the study template was built from the brain of the session means (manual masks)."""
    side = pu.sidecar(template)
    return bool(opi(side) and pu.read_json(side).get('BrainMasked', False))


def initial_template(images, pad_mm=15.0):
    """
    Domain and starting image of ants.build_template. Without them ANTsPy starts from the sum of the images resampled
    in world coordinates on the grid of the first one: sessions whose head sits elsewhere in the scanner add shifted
    copies of the brain (ghost), the template drifts towards the middle of them and leaves the grid of the first image
    (brain cut at the edges of the template, 12 mm on 3 faces on a macaque FDG study of 9 sessions; reproduced on the
    phantom with 15 mm offsets). Same order as sammba-mri anats_to_common (Celestine et al. 2020): centring, rigid,
    then the affine / non-linear iterations.
    - grid: spacing and axes of the first image, size = largest extent of the images on each axis + 2 x pad_mm,
      centred on the centre of mass of the first image;
    - start: mean of every image registered rigidly to the first one (centre-of-mass start of ants.registration).
    """
    from modalities.pet import coreg2anat
    # CT (air < 0 after scaling): shifted to >= 0 so that the centres of mass and the zero padding mean air
    low = min(float(im.min()) for im in images)
    shift = -low if low < 0 else 0.0
    images = [im + shift for im in images] if shift else images
    ref = images[0]
    sp = np.asarray(ref.spacing[:3], float)
    D = np.asarray(ref.direction)[:3, :3]
    ext = np.max([np.asarray(im.shape[:3], float) * np.asarray(im.spacing[:3], float) for im in images], axis=0)
    size = np.ceil((ext + 2.0 * pad_mm) / sp).astype(int)
    com = np.asarray(ants.get_center_of_mass(ref), float)
    origin = com - D @ ((size - 1) / 2.0 * sp)
    dom = ants.make_image(tuple(int(x) for x in size), 0.0, spacing=tuple(sp), origin=tuple(origin), direction=D)
    fixed = ants.resample_image_to_target(ref, dom, interp_type='linear')
    acc = fixed.clone()
    for im in images[1:]:
        reg = ants.registration(fixed, im, type_of_transform='Rigid', **coreg2anat._adaptive_reg_params(fixed, im))
        acc = acc + reg['warpedmovout']
    return acc / float(len(images)) - shift


def make_template(study_dir, sessions, tracer, transform, diary_file, suffix='pet', iterations=4, masked=False):
    """
    sessions: list of (session mean image, orig/matrices folder, session base name).
    masked: the images are the brain of the session means (manualmask.masked_image); written in the sidecar
    (BrainMasked), read by step 6 of scenario 3v2 to register the brain of each motion-corrected mean.
    Writes the template and, for every session, <matrices>/<base>_to_studyTemplate_* (fixed = template).
    """
    nl = 'Run make_template (study template of the ' + suffix.upper() + ' images)'
    run_cmd.msg(nl, diary_file, 'HEADER')
    sessions = [s for s in sessions if opi(s[0])]
    if len(sessions) < 2:
        raise ValueError(f'a study template needs at least 2 sessions with a mean image; found {len(sessions)}')
    os.makedirs(opj(study_dir, 'matrices'), exist_ok=True)
    images = [_normalised(s[0]) for s in sessions]
    tpl = ants.build_template(initial_template=initial_template(images), image_list=images, iterations=iterations,
                              gradient_step=0.2, blending_weight=0.75, type_of_transform=transform)
    fname = template_name(study_dir, tracer, suffix)
    ants.image_write(tpl, fname, ri=False)
    pu.write_json(fname, [s[0] for s in sessions],
                  f'Study template ({len(sessions)} sessions, ants.build_template started from the rigid mean on a '
                  f'padded grid centred on the first session, {transform}, '
                  f'{iterations} iterations; images scaled by their 99th percentile' +
                  ('; brain of the session means (manual masks dilated, cropped)' if masked else '') + ').',
                  BrainMasked=bool(masked))

    for (img, matdir, base), im in zip(sessions, images):
        os.makedirs(matdir, exist_ok=True)
        transfo.register(tpl, im, opj(matdir, base + '_to_studyTemplate_'), transform)
        run_cmd.msg('INFO: ' + base + ' registered to the study template', diary_file, 'OKGREEN')
    return fname


def manual_crop(template, diary_name, sing_fs, redo=False):
    """
    Crop drawn by hand (freeview: Tools > Crop Volume, then File > Save Volume under the proposed name).
    The crop file is created only once: the engineer's version overwrote a previous crop at every run.
    """
    crop = crop_name(template)
    if opi(crop) and not redo:
        run_cmd.msg('INFO: crop already drawn: ' + crop, diary_name, 'OKGREEN')
        return crop
    shutil.copyfile(template, crop)
    print('CROP: open Tools / Crop Volume')
    print('For a brain image: keep the skin around the skull and the eyeballs (if visible); remove most of '
          'the external soft tissues (ears, tip of the nose) and the lower jaw.')
    print('Once satisfied: Apply, then File / Save Volume (keep the name ' + opb(crop) + ')')
    run_cmd.do(sing_fs + 'freeview -v ' + crop + ' -subtitle "CROP-THE-IMAGE!"', diary_name)
    if nib.load(crop).shape == nib.load(template).shape:
        run_cmd.msg('WARNING: the crop has the size of the template: no crop was drawn (the whole field of '
                    'view will be used for the motion correction)', diary_name, 'WARNING')
    pu.write_json(crop, template, 'Manually cropped image (freeview, FreeSurfer).')
    return crop


def manual_crop_session(image, crop, diary_name, sing_fs, redo=False):
    """Same crop, drawn on the mean image of one session (no study template); written to `crop`."""
    os.makedirs(opd(crop), exist_ok=True)
    if opi(crop) and not redo:
        return crop
    shutil.copyfile(image, crop)
    print('CROP the mean image of this session: Tools / Crop Volume, Apply, File / Save Volume (' + opb(crop) + ')')
    run_cmd.do(sing_fs + 'freeview -v ' + crop + ' -subtitle "CROP-THE-IMAGE!"', diary_name)
    pu.write_json(crop, image, 'Manually cropped image (freeview, FreeSurfer).')
    return crop


def box(crop):
    """Image of ones on the grid of a cropped image (the box drawn by the user)."""
    img = ants.image_read(crop)
    return img * 0 + 1


def crop_mask(crop, session_image, out_mask, chain=None):
    """
    Box of `crop` in the space of `session_image`: through the session -> template chain (chain['from'],
    fixed = session image) when the crop was drawn on the study template, by world coordinates otherwise.
    The engineer applied [affine, inverse warp] without inverting the affine, and thresholded the
    template intensities instead of the box.
    """
    ref = ants.image_read(session_image)
    if ref.dimension == 4:
        ref = ants.slice_image(ref, axis=3, idx=0)
    b = box(crop)
    if chain is None:
        moved = ants.resample_image_to_target(b, ref, interp_type='nearestNeighbor')
    else:
        moved = transfo.apply(ref, b, chain, 'from', interpolator='nearestNeighbor')
    moved = ants.threshold_image(moved, 0.5, 1e9, 1, 0)
    if moved.sum() == 0:
        raise ValueError(f'the crop box of {crop} does not overlap {session_image}')
    ants.image_write(moved, out_mask, ri=False)
    pu.write_json(out_mask, [crop, session_image], 'Box of the manual crop in the space of this session.')
    return out_mask


def cropped(img, mask):
    """img cropped to the bounding box of mask (ANTsImage, 3D)."""
    return ants.crop_image(img, label_image=mask, label=1)
