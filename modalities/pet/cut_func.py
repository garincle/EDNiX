"""PET signal during every fMRI run of the session (step 10, simultaneous PET-MR): see run_cut."""
import os
import glob

import nibabel as nib

opj = os.path.join
opb = os.path.basename
opi = os.path.isfile

from Tools import run_cmd
from modalities.pet import petutils as pu
from modalities.pet import run_cut


def cut(data_path, endfmri, pet_file, timing, pet_info, T1_eq, brain_mask, out_dir, diary_name):
    """For every raw BOLD run of <data_path>/func matching endfmri: <pet base>_desc-<task><run>fMRI_pet.nii.gz."""
    runs = sorted(glob.glob(opj(data_path, 'func', endfmri)))
    runs = [r for r in runs if r.endswith(('.nii.gz', '.nii'))]
    if not runs:
        run_cmd.msg('WARNING: no fMRI run found with ' + opj(data_path, 'func', endfmri), diary_name, 'WARNING')
        return []
    base = pu.base_name(pet_file)
    out = []
    for r in runs:
        info = pu.raw_sidecar(r)
        if 'AcquisitionTime' not in info or 'RepetitionTime' not in info:
            run_cmd.msg('WARNING: ' + opb(r) + ': no AcquisitionTime / RepetitionTime in its sidecar', diary_name,
                        'WARNING')
            continue
        nvol = nib.load(r).shape[3]
        ov, (a, b) = run_cut.getframes(pet_info, timing, info, nvol, T1_eq)
        if ov.sum() <= 0:
            run_cmd.msg(f'WARNING: {opb(r)} [{a:.0f}, {b:.0f}] s after the PET TimeZero: outside the PET acquisition',
                        diary_name, 'WARNING')
            continue
        ents = {p.split('-')[0]: p.split('-', 1)[1] for p in pu.strip_ext(r).split('_') if '-' in p}
        tag = ''.join(ents.get(k, '') for k in ('task', 'run')) + 'fMRI'
        dst = pu.out_name(out_dir, base, desc=tag)
        run_cut.cut(pet_file, timing, ov, dst, brain_mask)
        run_cmd.msg(f'INFO: {opb(r)}: PET frames weighted {[round(x, 1) for x in ov[ov > 0]]} s', diary_name, 'OKGREEN')
        out.append(dst)
    return out
