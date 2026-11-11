"""V2h surface method (preFS.SURFACE_METHOD): images identical to the validated reference implementation
(test_surface_variants_library.prepa_v2, hippo_wm=False, PR #84), pial cap in real mm, Pcreate command."""
import importlib.util
import os

import ants
import numpy as np
import pytest

from modalities.anat.freesurfer import preFS, FS_surf

HERE = os.path.dirname(__file__)


def _reference():
    spec = importlib.util.spec_from_file_location('ref', os.path.join(HERE, '_surface_variants_reference.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _synthetic(tmp_path):
    rng = np.random.default_rng(0)
    n = 40
    z, y, x = np.mgrid[:n, :n, :n]
    d = np.sqrt((z - 20) ** 2 + (y - 20) ** 2 + (x - 20) ** 2)
    left = z < 20                      # first axis = left/right, as in the filled rule
    a = np.zeros((n, n, n), np.int32)
    a[d < 16] = np.where(left, 3, 42)[d < 16]
    a[d < 12] = np.where(left, 2, 41)[d < 12]
    a[d < 4] = np.where(left, 4, 43)[d < 4]
    a[(d >= 5) & (d < 7) & (y < 20)] = np.where(left, 10, 49)[(d >= 5) & (d < 7) & (y < 20)]
    a[(d >= 8) & (d < 10) & (x < 12)] = np.where(left, 17, 53)[(d >= 8) & (d < 10) & (x < 12)]
    a[(d >= 16) & (d < 18)] = 24
    t1 = np.select([np.isin(a, [4, 43, 24]), np.isin(a, [3, 42, 10, 49, 17, 53]), np.isin(a, [2, 41])],
                   [20., 60., 100.], 0.) + rng.normal(0, 3, a.shape) * (a > 0)
    sp = (0.2, 0.2, 0.2)
    t1_f, a_f = str(tmp_path / 'sub_space-acpc_desc-SS_T1w.nii.gz'), str(tmp_path / 'aseg.nii.gz')
    ants.image_write(ants.from_numpy(t1.astype('float32'), spacing=sp), t1_f)
    ants.image_write(ants.from_numpy(a.astype('float32'), spacing=sp), a_f)
    return t1_f, a_f, a


def test_v2h_matches_reference(tmp_path):
    t1_f, a_f, a = _synthetic(tmp_path)
    mine, ref = tmp_path / 'mine', tmp_path / 'ref'
    mine.mkdir()
    ref.mkdir()
    preFS.prepa_img('sub', t1_f, str(mine), a_f, str(mine), '', str(tmp_path / 'diary.txt'), method='V2h')
    _reference().prepa_v2({'preFS': preFS}, 'sub', t1_f, a_f, str(ref), const=False, hippo_wm=False)
    for name in ('sub_desc-norm_T1w', 'sub_seg-wm_dseg', 'sub_seg-filled_dseg'):
        m = ants.image_read(str(mine / (name + '.nii.gz'))).numpy()
        r = ants.image_read(str(ref / (name + '.nii.gz'))).numpy()
        assert np.allclose(m, r, atol=1e-4), name
    norm = ants.image_read(str(mine / 'sub_desc-norm_T1w.nii.gz')).numpy()
    wm = ants.image_read(str(mine / 'sub_seg-wm_dseg.nii.gz')).numpy()
    filled = ants.image_read(str(mine / 'sub_seg-filled_dseg.nii.gz')).numpy()
    hippo = np.isin(a, [17, 53])
    assert (wm[hippo] == 0).all() and (filled[hippo] == 0).all()     # hippocampus out of wm / filled
    assert 50 < np.median(norm[hippo]) < 90                          # and in the GM class of brain.mgz
    assert set(np.unique(filled)) <= {0, 127, 255}


def test_pial_max():
    assert preFS.pial_max(1.0, 'V2h') == 5            # human / macaque: FreeSurfer default kept
    assert preFS.pial_max(0.283, 'V2h') == 11         # rat: 3 real mm
    assert preFS.pial_max(0.174, 'V2h') == 18         # mouse
    assert preFS.pial_max(0.283, 'legacy') is None


@pytest.mark.parametrize('mx, expected', [(None, '-mgz -T1 brain'), (11.0, '-mgz -max 11 -T1 brain')])
def test_pcreate_command(monkeypatch, mx, expected):
    cmds = []
    monkeypatch.setattr(FS_surf.run_cmd, 'do', lambda c, d: cmds.append(c))
    monkeypatch.setattr(FS_surf.run_cmd, 'msg', lambda *a, **k: None)
    FS_surf.Pcreate('sub', 'lr', 'diary', '', max_fs_mm=mx)
    assert len(cmds) == 2 and all(expected in c for c in cmds)


def test_fs_scale_factor_from_rescale_dat(tmp_path):
    """FS_2_WB_short unscales small-brain surfaces by 1/s read from <ID>_rescale.dat (smallbrain.get layout)."""
    from modalities.anat.freesurfer import FS_2_WB_short as W
    dat = tmp_path / 'Marmoset_rescale.dat'
    dat.write_text('Marmoset\n0.15\n0.15\n0.15\n3.2 0 0 0.9 \n0 3.2 0 0 \n0 0 3.2 0.9 \n0 0 0 1\nround\n')
    assert W.fs_scale_factor(str(dat)) == 3.2
    assert W.fs_scale_factor('') == 1.0
