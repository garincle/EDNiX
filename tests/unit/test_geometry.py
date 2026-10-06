import numpy as np
import nibabel as nib
import pytest

from Tools import get_orientation, check_nii

AXES = {'R': (0, 1), 'L': (0, -1), 'A': (1, 1), 'P': (1, -1), 'S': (2, 1), 'I': (2, -1)}


def _matrix(code):
    m = np.zeros((3, 3))
    for i, c in enumerate(code):
        j, s = AXES[c]
        m[j, i] = s
    return m


@pytest.mark.parametrize('animal', ['AHF', 'AFF', 'AFS', 'AHS'])
@pytest.mark.parametrize('human', ['HFS', 'FFS', 'HFP', 'FFP'])
def test_orientation_is_a_proper_rotation(animal, human):
    # a reflection would silently swap left and right
    assert round(np.linalg.det(_matrix(get_orientation.getreal(human, animal, 'RAS')))) == 1


def test_orientation_humanlike_is_identity():
    assert get_orientation.getreal('HFS', 'humanlike', 'LPI') == 'LPI'


def test_sphinx_head_first_supine_scanner_labels():
    # sphinx (prone, head first) declared as HFS: scanner R = animal L, A = S, S = A
    assert get_orientation.getreal('HFS', 'AHF', 'RAS') == 'LSA'


def _pair(tmp_path, shift_mm):
    A = np.diag([0.1, 0.1, 0.1, 1.0])
    B = A.copy()
    B[0, 3] = shift_mm
    x = np.zeros((20, 20, 20), 'f4')
    x[5:10, 5:10, 5:10] = 1
    a, b = tmp_path / 'a.nii.gz', tmp_path / 'b.nii.gz'
    nib.save(nib.Nifti1Image(x, A), a)
    nib.save(nib.Nifti1Image(x, B), b)
    return str(a), str(b)


def test_comphd_detects_a_5_voxel_shift(tmp_path):
    a, b = _pair(tmp_path, 0.5)
    assert check_nii.comphd(a, a)
    assert not check_nii.comphd(a, b)


def test_resamp_resamples_when_only_the_origin_differs(tmp_path):
    a, b = _pair(tmp_path, 0.5)
    check_nii.resamp(b, a, 'msk', '', '', str(tmp_path / 'diary.txt'), '')
    assert np.allclose(nib.load(b).affine, nib.load(a).affine)


def test_keep_header_refuses_mismatched_grids(tmp_path):
    a, _ = _pair(tmp_path, 0.0)
    c = tmp_path / 'c.nii.gz'
    nib.save(nib.Nifti1Image(np.zeros((10, 10, 10), 'f4'), np.eye(4)), c)
    with pytest.raises(ValueError):
        check_nii.keep_header(str(c), a)
