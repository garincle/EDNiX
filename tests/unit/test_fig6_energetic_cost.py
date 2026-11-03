import os
import numpy as np
import nibabel as nib
from scipy import ndimage

from Exemples.Study_EDNiX.analysis import fig6_energetic_cost as F
from Exemples.Study_EDNiX.analysis.library import multimodal_figure as mmf


def test_coupling_for_species_end_to_end(tmp_path, monkeypatch):
    rng = np.random.default_rng(0)
    shape, nt = (14, 14, 14), 120
    aff = np.eye(4)
    bids = tmp_path / 'BIDS_X'
    pet_dir = tmp_path / 'BIDS_PET'
    # a smooth "degree" field built from a spatially structured shared signal
    field = ndimage.gaussian_filter(rng.standard_normal(shape), 2.5)
    field = (field - field.min()) / np.ptp(field)
    shared = rng.standard_normal(nt)
    for sub in ('01', '02', '03'):
        d = bids / f'sub-{sub}' / 'ses-1' / 'func' / 'templates' / 'EDNiX'
        (d / 'postprocessed_rs').mkdir(parents=True)
        (d / 'labels').mkdir()
        data = rng.standard_normal(shape + (nt,)).astype(np.float32) + 3 * field[..., None] * shared
        nib.save(nib.Nifti1Image(data, aff), str(d / 'postprocessed_rs' / f'sub-{sub}_task-rest_bold_space-template_desc-fMRI_residual.nii.gz'))
        lab = np.zeros(shape + (3,), np.int16)
        lab[2:12, 2:12, 2:12, 2] = 1
        nib.save(nib.Nifti1Image(lab, aff), str(d / 'labels' / f'{sub}_seg-EDNIxCSCLR_dseg.nii.gz'))
    (pet_dir / 'Results' / 'Grp_PET_mean_desc-bratio').mkdir(parents=True)
    pet = field + .3 * ndimage.gaussian_filter(rng.standard_normal(shape), 2.5)
    nib.save(nib.Nifti1Image(pet.astype(np.float32), aff), str(pet_dir / F.PET_MAP))
    lab_txt = tmp_path / 'labels.txt'
    lab_txt.write_text('0 0 0 0 0 0 0 "Clear Label"\n1 255 0 0 1 1 1 "Prefrontal_cortex_(PFC)"\n')
    c = dict(bids_configs=[dict(label='DS', species='Test', bids_dir=str(bids), label_path=str(lab_txt))],
             kept_runs_final={'DS': {('01', '1'), ('02', '1'), ('03', '1')}})
    monkeypatch.setitem(mmf.SOURCES['metabolism'], 'Test', [str(pet_dir)])
    monkeypatch.setitem(F.FWHM_MM, 'Test', 2.0)
    monkeypatch.setattr(F, 'OUT', str(tmp_path / 'out'))
    monkeypatch.setattr(F, 'N_SURR', 100)
    os.makedirs(F.OUT)
    st, reg, used = F.coupling_for_species(c, 'Test')
    assert st['n_runs'] == 3 and st['r'] > .4 and st['pet_overlap'] > .99
    assert os.path.exists(os.path.join(F.OUT, 'Test_energetic_cost.nii.gz'))
    assert len(used) == 3 and 'residual_bratio' in reg.columns
