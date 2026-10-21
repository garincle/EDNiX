"""
Parametric images of the PET space (step 7): SUV, glucose-corrected SUV, ratio to the brain mean,
skull-stripped versions. numpy / nibabel; every image keeps the geometry of its input.

SUV (body weight)   = C (Bq/mL) x body weight (g) / injected radioactivity (Bq), with C decay-corrected
                      to the injection time (or the dose decayed to the decay-correction time); 1 g/mL
                      (Boellaard et al. 2009, 2015; Kinahan & Fletcher 2010). The PET engineer's ratio
                      (MBq / kg x 1000) is the same quantity when C is in Bq/mL; the image units are now
                      read from the sidecar.
SUVglu              = SUV x glucose / 5.0 mmol/L (EANM; the engineer multiplied by the glycaemia in g/L
                      without reference value).
bratio              = C / mean of C over the brain mask (global normalisation, dimensionless). The
                      engineer's mean included the zeros outside the brain in the time blocks.
"""
import numpy as np

from modalities.pet import petutils as pu

GLUCOSE_REF_MMOL_L = 5.0


def suv(file, factor, out, report):
    img, d = pu.load(file)
    pu.save_like(d * float(factor), img, out)
    pu.write_json(out, file, 'Standardized uptake value (body weight), g/mL.', SUV=report)
    return out


def suvglu(suv_file, glucose_mmol_L, out):
    img, d = pu.load(suv_file)
    f = float(glucose_mmol_L) / GLUCOSE_REF_MMOL_L
    pu.save_like(d * f, img, out)
    pu.write_json(out, suv_file, f'SUV corrected for the blood glucose: x {glucose_mmol_L:.3f} / '
                  f'{GLUCOSE_REF_MMOL_L} mmol/L.')
    return out


def brain_mean(data3d, brain):
    v = data3d[brain & np.isfinite(data3d)]
    if v.size == 0:
        raise ValueError('empty brain mask in the PET space')
    m = float(v.mean())
    if not np.isfinite(m) or m == 0:
        raise ValueError('mean of the image over the brain mask is zero or not finite')
    return m


def bratio(file, brain_mask, out):
    """Image divided by its mean over the brain mask (3D image)."""
    img, d = pu.load(file)
    if d.ndim > 3:
        raise ValueError(f'bratio: {file} is 4D; use its time average')
    brain = pu.mask_bool(brain_mask, file)
    m = brain_mean(d, brain)
    pu.save_like(d / m, img, out)
    pu.write_json(out, [file, brain_mask], 'Ratio to the mean over the brain mask (global normalisation).',
                  BrainMean=m)
    return out, m


def skullstrip(file, brain_mask, out):
    img, d = pu.load(file)
    brain = pu.mask_bool(brain_mask, file)
    if d.ndim > 3:
        d = d * brain[..., None]
    else:
        d = d * brain
    pu.save_like(d, img, out)
    pu.write_json(out, [file, brain_mask], 'Skull-stripped (brain mask in the PET space).')
    return out
