"""
Partial-volume correction of a static PET image (optional, PETPVC: Thomas et al. 2016, Phys Med Biol
61:7975), with the grey matter, white matter and CSF masks in the PET space.

The PET engineer passed the voxel size as the point-spread function: the PSF is the spatial resolution
of the scanner and reconstruction (FWHM, mm), a property of the acquisition that must be given
(pvc_fwhm in the launcher); the correction is skipped without it. The segmentation given to PETPVC is a
4D image (grey, white, CSF, rest) whose volumes sum to one in every voxel.
"""
import os

import numpy as np

opi = os.path.isfile

from Tools import run_cmd
from modalities.pet import petutils as pu


def correct(img, masks, fwhm_mm, out, diary_name, algo='RBV'):
    """
    algo: MG, GTM, LABBE, MTC, RBV, IY... (PETPVC names). masks: dict gray / white / csf (PET space).
    Needs nipype and the petpvc program.
    """
    try:
        from nipype.interfaces.petpvc import PETPVC
    except ImportError as e:
        raise ImportError('partial-volume correction needs nipype and PETPVC (pip install nipype; petpvc on the '
                          'PATH)') from e
    if fwhm_mm is None:
        raise ValueError('pvc_fwhm (scanner resolution, FWHM in mm) is needed for the partial-volume correction')
    fw = np.atleast_1d(np.asarray(fwhm_mm, dtype=float))
    fw = np.repeat(fw, 3) if fw.size == 1 else fw
    for k in ('gray', 'white', 'csf'):
        if not masks.get(k) or not opi(masks[k]):
            raise FileNotFoundError(f'{k} mask in the PET space needed for the partial-volume correction')
    nii, d = pu.load(img)
    seg = [pu.mask_bool(masks[k], img).astype(np.float32) for k in ('gray', 'white', 'csf')]
    total = np.clip(seg[0] + seg[1] + seg[2], 0, None)
    seg = [s / np.where(total > 1, total, 1) for s in seg]
    seg.append(np.clip(1 - sum(seg), 0, 1))
    seg_file = out.replace('.nii.gz', '_seg4D.nii.gz')
    pu.save_like(np.stack(seg, axis=3), nii, seg_file)
    run_cmd.msg(f'Partial-volume correction ({algo}, PSF FWHM {fw.tolist()} mm)', diary_name, 'OKGREEN')
    p = PETPVC()
    p.inputs.in_file = img
    p.inputs.mask_file = seg_file
    p.inputs.out_file = out
    p.inputs.pvc = algo
    p.inputs.fwhm_x, p.inputs.fwhm_y, p.inputs.fwhm_z = [float(x) for x in fw]
    p.run()
    pu.write_json(out, [img, seg_file], f'Partial-volume correction (PETPVC {algo}, PSF FWHM {fw.tolist()} mm).')
    return out
