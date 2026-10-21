"""
Time blocks of a dynamic PET image (Do_cut, step 7): duration-weighted mean of the frames of every block
of `duration` minutes, divided by its mean over the brain mask.

The PET engineer assigned the frames to the blocks with frame indices shifted by one and averaged the
frames without weights; here a frame belongs to the block that contains its middle, and the frames are
weighted by their duration.
"""
import numpy as np

from modalities.pet import petutils as pu
from modalities.pet import paramImg


def setTime(timing, duration):
    """[(t0, t1, frames)] of the blocks of `duration` minutes from the start of the first frame (seconds)."""
    if timing is None:
        raise ValueError('time blocks need the frame timing')
    step = 60.0 * float(duration)
    t0 = float(timing['start'][0])
    t_end = float(timing['end'][-1])
    blocks = []
    while t0 < t_end - 1e-6:
        t1 = t0 + step
        frames = np.where((timing['mid'] >= t0) & (timing['mid'] < t1))[0]
        if frames.size:
            blocks.append((t0, t1, frames))
        t0 = t1
    return blocks


def doit(img, timing, duration, brain_mask, out_dir, base, desc):
    """Writes <base>_desc-<desc>-period<t0>to<t1>min_pet.nii.gz for every block; returns the files."""
    nii, d = pu.load(img)
    brain = pu.mask_bool(brain_mask, img) if brain_mask else None
    files = []
    for t0, t1, frames in setTime(timing, duration):
        m = pu.weighted_mean(d, timing, frames)
        tag = f'{desc}-period{t0 / 60:g}to{t1 / 60:g}min'
        out = pu.out_name(out_dir, base, desc=tag)
        if brain is not None:
            bm = paramImg.brain_mean(m, brain)
            m = m / bm
            descr = f'Duration-weighted mean of frames {frames.tolist()} divided by its brain mean ({bm:.4g}).'
        else:
            descr = f'Duration-weighted mean of frames {frames.tolist()} (no brain mask: not normalised).'
        pu.save_like(m, nii, out)
        pu.write_json(out, img, descr, Block_min=[t0 / 60, t1 / 60])
        files.append(out)
    return files
