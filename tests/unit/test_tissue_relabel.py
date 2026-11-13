"""tissue_relabel: recovers the GM/WM boundary of a synthetic brain, never grows ventricles through sulci,
does nothing without contrast, and the 'coreg' option reaches the pipeline."""
import numpy as np

from modalities.anat import tissue_relabel as T
from modalities.anat import transfoparams


def _brain(n=48):
    rng = np.random.default_rng(0)
    z, y, x = np.mgrid[:n, :n, :n]
    c = n / 2
    d = np.sqrt((z - c) ** 2 + (y - c) ** 2 + (x - c) ** 2)
    ang = np.arctan2(y - c, x - c)
    left = x < c
    true = np.zeros(d.shape, int)
    true[d < 20] = np.where(left, 3, 42)[d < 20]
    wm = d < 14 + 2 * np.sin(6 * ang)
    true[wm] = np.where(left, 2, 41)[wm]
    true[d < 5] = np.where(left, 4, 43)[d < 5]
    mask = d < 22
    csf = mask & ((d >= 20) | np.isin(true, (4, 43)))
    t1 = np.select([csf, np.isin(true, (3, 42)), np.isin(true, (2, 41))], [20, 60, 100], 0).astype(float)
    t1 = (t1 + rng.normal(0, 5, t1.shape)) * mask
    reg = np.zeros(d.shape, int)
    reg[d < 20] = np.where(left, 3, 42)[d < 20]
    reg[d < 14] = np.where(left, 2, 41)[d < 14]
    reg[d < 4] = np.where(left, 4, 43)[d < 4]
    return true, reg, t1, mask, d


def _dice(a, b, codes):
    A, B = np.isin(a, codes), np.isin(b, codes)
    return 2 * (A & B).sum() / (A.sum() + B.sum())


def test_recovers_boundary():
    true, reg, t1, mask, _ = _brain()
    out, rep = T.compute(reg, [t1], mask, (1, 1, 1))
    assert rep['applied']
    for codes in (T.CORTEX, T.WM, T.VENT):
        assert _dice(true, out, codes) > _dice(true, reg, codes)
    assert rep['boundary_gradient_after'] > rep['boundary_gradient_before']


def test_no_ventricle_in_sulci():
    true, reg, t1, mask, d = _brain()
    out, _ = T.compute(reg, [t1], mask, (1, 1, 1))
    assert not np.isin(out, T.VENT)[d >= 10].any()


def test_no_contrast_keeps_registration():
    _, reg, t1, mask, _ = _brain()
    noise = np.random.default_rng(1).normal(50, 5, t1.shape) * mask
    out, rep = T.compute(reg, [noise], mask, (1, 1, 1))
    assert not rep['applied'] and np.array_equal(out, reg)


def test_coreg_option_passthrough():
    blocks = [dict(name='coreg', type_of_transform='SyNCC', tissue_relabel=True)]
    _, out = transfoparams.set('do_as_I_said', True, False, '', '', '', '', blocks)
    coreg = next(b for b in out if b['name'] == 'coreg')
    assert coreg['tissue_relabel'] is True
    assert all('tissue_relabel' not in b for b in out if b['name'] != 'coreg')
    _, out = transfoparams.set('do_it_for_me', True, False, '', '', '', '', [])
    assert next(b for b in out if b['name'] == 'coreg')['tissue_relabel'] is False
