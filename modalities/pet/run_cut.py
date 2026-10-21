"""
PET signal during an fMRI run of a simultaneous PET-MR session (step 10, optional).

The run occupies [AcquisitionTime + T1_eq x TR, AcquisitionTime + n_volumes x TR] (clock time of the BOLD
sidecar); the PET frames are placed on the same clock with the BIDS TimeZero of the PET sidecar. The mean
concentration during the run is the average of the frames weighted by their overlap with the run.
The PET engineer matched the clock times with time.mktime on strings (date and sub-second parts lost) and
took whole frames on both sides; frames are now weighted by their overlap.
"""
import datetime

import numpy as np

from modalities.pet import petutils as pu


def _clock(s):
    s = str(s)
    for fmt in ('%H:%M:%S.%f', '%H:%M:%S'):
        try:
            t = datetime.datetime.strptime(s, fmt)
            return t.hour * 3600 + t.minute * 60 + t.second + t.microsecond / 1e6
        except ValueError:
            continue
    raise ValueError(f'clock time not understood: {s!r} (expected hh:mm:ss[.ffffff])')


def run_window(func_info, n_volumes, T1_eq):
    """(start, end) of the analysed part of the run, seconds of clock time."""
    tr = float(func_info['RepetitionTime'])
    t0 = _clock(func_info['AcquisitionTime'])
    return t0 + T1_eq * tr, t0 + n_volumes * tr


def getframes(pet_info, timing, func_info, n_volumes, T1_eq):
    """Overlap (s) of every PET frame with the run, and the run window relative to the PET TimeZero."""
    tz = _clock(pet_info.get('TimeZero') or pet_info['AcquisitionTime'])
    a, b = run_window(func_info, n_volumes, T1_eq)
    a, b = a - tz, b - tz
    if b < a:
        raise ValueError('fMRI run ends before it starts (check AcquisitionTime / RepetitionTime)')
    ov = np.clip(np.minimum(timing['end'], b) - np.maximum(timing['start'], a), 0, None)
    return ov, (a, b)


def cut(pet_file, timing, overlap, out, brain_mask=None):
    """Mean concentration during the run (frames weighted by their overlap); ratio to the brain mean if a mask."""
    from modalities.pet import paramImg
    if overlap.sum() <= 0:
        raise ValueError('no PET frame overlaps the fMRI run')
    img, d = pu.load(pet_file)
    w = overlap / overlap.sum()
    m = np.tensordot(d, w, axes=([3], [0]))
    descr = 'PET during the fMRI run (frames weighted by their overlap with the run).'
    if brain_mask:
        bm = paramImg.brain_mean(m, pu.mask_bool(brain_mask, pet_file))
        m = m / bm
        descr += f' Divided by its brain mean ({bm:.4g}).'
    pu.save_like(m, img, out)
    pu.write_json(out, pet_file, descr, FrameWeights=w.tolist())
    return out
