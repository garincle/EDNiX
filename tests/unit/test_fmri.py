import os
import numpy as np
import nibabel as nib
import ants

from modalities.fMRI._7_post_TTT import _apply_normalisation
from modalities.fMRI._10_Correl_matrix import erode_atlas_labels
from modalities.fMRI.transfo_for_func import apply_motion_correction_and_transforms


def test_psc_uses_the_pre_regression_mean(tmp_path):
    rng = np.random.default_rng(0)
    base = 1000 + rng.normal(0, 10, (4, 4, 4, 60))
    res = rng.normal(0, 10, (4, 4, 4, 60))
    keep = np.ones(60, int)
    keep[[3, 4]] = 0
    res[..., keep == 0] = 0
    np.savetxt(tmp_path / 'c.1D', keep, fmt='%d')
    for n, a in (('f', base), ('r', res)):
        nib.save(nib.Nifti1Image(a.astype('f4'), np.eye(4)), tmp_path / f'{n}.nii.gz')
    _apply_normalisation(str(tmp_path / 'r.nii.gz'), str(tmp_path / 'm.nii.gz'), str(tmp_path / 's.nii.gz'), 'psc', '',
                         str(tmp_path / 'd.txt'), str(tmp_path / 'f.nii.gz'), str(tmp_path / 'c.1D'))
    out = nib.load(tmp_path / 'r.nii.gz').get_fdata()
    assert 0.8 < out[..., keep == 1].std() < 1.2          # 10 / 1000 * 100 = 1 %
    assert np.all(out[..., keep == 0] == 0)


def test_erosion_really_erodes_two_voxels():
    a = np.zeros((20, 20, 20))
    a[2:12, 2:12, 2:12] = 1
    out, report = erode_atlas_labels(a, 2)
    assert report[1] == 2 and (out == 1).sum() == 6 ** 3


def test_one_shot_resampling_recovers_motion(tmp_path):
    x, y, z = np.mgrid[0:30, 0:30, 0:30].astype(float)
    ref = np.exp(-((x - 14) ** 2 / 30 + (y - 16) ** 2 / 50 + (z - 15) ** 2 / 20)) * 100
    A = np.diag([2.0, 2.0, 2.0, 1.0])
    nib.save(nib.Nifti1Image(ref.astype('f4'), A), tmp_path / 'ref.nii.gz')
    R = ants.image_read(str(tmp_path / 'ref.nii.gz'))
    vols = []
    for k in range(4):
        t = ants.create_ants_transform(transform_type='Euler3DTransform', dimension=3)
        t.set_parameters([0.02 * k, 0, 0, 0.8 * k, 0, 0.4 * k])
        t.set_fixed_parameters([0, 0, 0])
        vols.append(t.apply_to_image(R, reference=R, interpolation='linear').numpy())
    nib.save(nib.Nifti1Image(np.stack(vols, -1).astype('f4'), A), tmp_path / 'run.nii.gz')
    os.makedirs(tmp_path / 'mc')
    ants.motion_correction(image=ants.image_read(str(tmp_path / 'run.nii.gz')), fixed=R,
                           type_of_transform='BOLDRigid', outprefix=str(tmp_path / 'mc' / 'img'))
    apply_motion_correction_and_transforms(str(tmp_path / 'run.nii.gz'), str(tmp_path / 'ref.nii.gz'),
                                           str(tmp_path / 'out'), str(tmp_path / 'mc' / 'img*.mat'), [], 2.0, '')
    out = nib.load(tmp_path / 'out.nii.gz').get_fdata()
    for k in range(4):
        assert np.corrcoef(out[..., k].ravel(), ref.ravel())[0, 1] > 0.99


def test_tr_is_always_a_float_in_seconds():
    import pytest
    from modalities.fMRI._0_Pipeline_launcher import resolve_TR
    assert resolve_TR('2.0', None)[0] == 2.0                        # given as text in a launcher
    assert resolve_TR(1.5, None)[0] == 1.5
    assert resolve_TR('Auto', {'RepetitionTime': 2.4})[0] == 2.4
    # single-band: 30 slices, 0.1 s apart -> 3.0 s
    st = [i * 0.1 for i in (list(range(0, 30, 2)) + list(range(1, 30, 2)))]
    assert abs(resolve_TR('Auto', {'SliceTiming': st})[0] - 3.0) < 1e-6
    # multiband 3: 30 slices, 10 distinct times 0.1 s apart -> 1.0 s
    mb = [(i % 10) * 0.1 for i in range(30)]
    assert abs(resolve_TR('Auto', {'SliceTiming': mb})[0] - 1.0) < 1e-6
    with pytest.raises(ValueError):
        resolve_TR('Auto', {})
    with pytest.raises(ValueError):
        resolve_TR('two seconds', None)



def test_dwell_time_is_the_effective_echo_spacing():
    from modalities.fMRI._0_Pipeline_launcher import resolve_dwell
    assert float(resolve_dwell('Auto', {'EffectiveEchoSpacing': 0.00058})[0]) == 0.00058
    v = float(resolve_dwell('Auto', {'TotalReadoutTime': 0.0551, 'ReconMatrixPE': 96})[0])
    assert abs(v - 0.0551 / 95) < 1e-12
    # BIDS DwellTime (receiver sample dwell, ~2.6 us) must not be used as an echo spacing
    assert resolve_dwell('Auto', {'DwellTime': 2.6e-6})[0] == 'None'
    assert resolve_dwell('Auto', {})[0] == 'None'
    assert resolve_dwell('0.0005', None)[0] == '0.0005'
