"""Step 7 (denoising): helpers, and the three methods end to end with AFNI / FSL emulated in numpy."""
import os
import numpy as np
import pytest

import modalities.fMRI._7_post_TTT as N


def test_band_parsing():
    assert N._parse_band('0.01 0.1', 2.0)[4] == ' -passband 0.01 0.1'
    assert N._parse_band('0.01', 2.0)[4] == ' -stopband 0 0.01'
    assert N._parse_band('0.01 0.3', 2.0)[1] is None          # low-pass at or above Nyquist: high-pass only
    for bad in ('0.1 0.01', '0', 'abc', '0.01 0.05 0.1'):
        with pytest.raises(ValueError):
            N._parse_band(bad, 2.0)


def test_filter_degrees_of_freedom():
    # 0.01-0.1 Hz removes most of the spectrum at short TR
    assert N._n_filter_regressors(300, 1.0, 0.01, 0.1) == 243
    assert N._n_filter_regressors(300, 2.0, 0.01, 0.1) == 189
    with pytest.raises(ValueError):
        N._check_dof(117, 28, 97, 'run', os.devnull)


@pytest.fixture
def emu(monkeypatch):
    import _afni_emulator as E
    E.install(monkeypatch)
    return E


@pytest.mark.parametrize('method', ['AFNI', 'Grandjean', 'fMRIPrep_style'])
@pytest.mark.parametrize('normalize', ['zscore', 'Skip'])
def test_methods_end_to_end(emu, method, normalize):
    out, keep, warn = emu.run(method, normalize)
    assert out is not None, warn
    if normalize == 'zscore':
        assert np.all(out[..., keep == 0] == 0)
        assert abs(out[6:10, 6:10, 4][..., keep == 1].std() - 1) < 0.05


def test_too_many_regressors_is_an_error(emu, monkeypatch):
    out, keep, warn = emu.run('fMRIPrep_style', 'zscore', T=120, tr=1.0)
    assert out is None and any('degrees of freedom' in w for w in warn)
    monkeypatch.setenv('EDNIX_STRICT', '1')
    with pytest.raises(RuntimeError):
        emu.run('fMRIPrep_style', 'zscore', T=120, tr=1.0)


def test_exterior_csf_is_not_a_parameter():
    # no step produces its mask and the signals come from skull-stripped data: the option was removed from
    # the denoising step and from preprocess_fmri (a launcher passing it gets the unknown-argument warning)
    import inspect
    from modalities.fMRI import _0_Pipeline_launcher as L
    from modalities.fMRI import _7_post_TTT as N
    assert 'extract_exterior_CSF' not in inspect.signature(N.signal_regression).parameters
    assert 'extract_exterior_CSF' not in inspect.signature(L.preprocess_fmri).parameters


def _out_of_band_fraction(out, tr, lo, hi):
    """Share of the (non-constant) spectral power of the residual outside [lo, hi] Hz."""
    ts = out[6:10, 6:10, 4].reshape(-1, out.shape[-1]).astype(np.float64)
    power = np.abs(np.fft.rfft(ts - ts.mean(axis=1, keepdims=True), axis=1)) ** 2
    f = np.fft.rfftfreq(out.shape[-1], tr)
    outside = (f < lo) | (f > hi)
    return power[:, outside].sum() / power.sum()


@pytest.mark.parametrize('method', ['AFNI', 'Grandjean', 'fMRIPrep_style'])
def test_bandpass_really_filters_both_edges(emu, method):
    """High-pass AND low-pass: the residual has no power outside the band (3dTproject / 1dBport emulated)."""
    out, keep, warn = emu.run(method, 'Skip', censor=False, GS=False, Vc=False)
    assert out is not None, warn
    assert _out_of_band_fraction(out, 2.0, 0.01, 0.1) < 1e-8
    # control: without nuisance filtering the data have most of their power outside the band
    raw = emu.make_subject(censor=False)
    import nibabel as nib
    d = np.asarray(nib.load(os.path.join(raw[0]['proc'], raw[1] + '_space-acpc-func_desc-fMRI_run_inRef.nii.gz')).dataobj)
    assert _out_of_band_fraction(d, 2.0, 0.01, 0.1) > 0.5


def test_high_pass_only_keeps_the_high_frequencies(emu):
    out, keep, warn = emu.run('Grandjean', 'Skip', band='0.01', censor=False, GS=False, Vc=False)
    assert out is not None, warn
    assert _out_of_band_fraction(out, 2.0, 0.01, 0.25) < 1e-8          # nothing below 0.01 Hz
    assert _out_of_band_fraction(out, 2.0, 0.0, 0.1) > 0.2             # but the power above 0.1 Hz is kept


@pytest.mark.parametrize('method', ['Grandjean', 'fMRIPrep_style'])
def test_header_tr_different_from_the_analysis_tr_is_refused(emu, method):
    out, keep, warn = emu.run(method, 'Skip', header_tr=1.0)
    assert out is None and any('wrong frequencies' in w for w in warn)


def test_header_tr_units(tmp_path):
    import nibabel as nib
    img = nib.Nifti1Image(np.zeros((2, 2, 2, 5), 'f4'), np.eye(4))
    img.header.set_zooms((1, 1, 1, 2000.0)); img.header.set_xyzt_units('mm', 'msec')
    nib.save(img, tmp_path / 'a.nii.gz')
    assert N._header_tr(str(tmp_path / 'a.nii.gz')) == 2.0
    N._check_header_tr(str(tmp_path / 'a.nii.gz'), 2.0)
    with pytest.raises(ValueError, match='3drefit'):
        N._check_header_tr(str(tmp_path / 'a.nii.gz'), 1.5)


@pytest.mark.parametrize('mode', ['fixed', 'adaptive', 'parallel'])
def test_acompcor_modes_run_and_report_the_parallel_analysis(emu, mode):
    out, keep, warn = emu.run('fMRIPrep_style', 'Skip', GS=False, Vc=False, acompcor_mode=mode)
    assert out is not None, warn
    assert any('parallel analysis' in m for c, m in emu.MSG)
    with pytest.raises(ValueError, match='acompcor_mode'):
        emu.run('fMRIPrep_style', 'Skip', acompcor_mode='variance')


def test_parallel_analysis_counts_shared_sources_not_noise():
    rng = np.random.default_rng(0)
    T, V = 200, 400
    t_basis = np.column_stack([np.ones(T), N._dct_basis(T, 2.0, 0.01)])

    def prepare(X):
        X = X - t_basis @ np.linalg.lstsq(t_basis, X, rcond=None)[0]
        return X / X.std(0)
    noise = prepare(rng.normal(size=(T, V)))
    assert N._n_components_parallel(noise, seed=1)[0] <= 1
    S = rng.normal(size=(T, 3)); S = np.cumsum(S, axis=0); S /= S.std(0)
    A = rng.normal(size=(3, V)); A /= np.sqrt((A ** 2).sum(0))
    mixed = prepare(0.6 * (S @ A) / (S @ A).std(0).mean() + rng.normal(size=(T, V)))
    assert N._n_components_parallel(mixed, seed=1)[0] == 3
    # the variance criterion of aCompCor needs far more components on the same data
    s = np.linalg.svd(mixed, compute_uv=False)
    assert np.searchsorted(np.cumsum(s ** 2 / (s ** 2).sum()), 0.5) + 1 > 15


def test_confounds_table_with_nan_in_the_first_row(tmp_path):
    import numpy as np
    import pandas as pd
    import pytest
    from modalities.fMRI._7_post_TTT import _read_confounds
    n = 50
    # fMRIPrep layout: derivative columns start with NaN; a text column may be present
    pd.DataFrame({'trans_x': np.arange(n, dtype=float),
                  'trans_x_derivative1': [np.nan] + [1.0] * (n - 1),
                  'framewise_displacement': [np.nan] + [0.1] * (n - 1),
                  'label': ['a'] * n}).to_csv(tmp_path / 'c.tsv', sep='\t', index=False)
    df = _read_confounds(str(tmp_path / 'c.tsv'), n)
    assert list(df.columns) == ['trans_x', 'trans_x_derivative1', 'framewise_displacement']
    assert not df.isna().any().any() and df['trans_x_derivative1'].iloc[0] == 0.0
    with pytest.raises(ValueError):
        _read_confounds(str(tmp_path / 'c.tsv'), n + 1)


def _shrink_mask(monkeypatch, name, n_voxels):
    """make_subject writes Wmask / Vmask; keep only the first n_voxels of one of them (small brain, eroded mask)."""
    import _afni_emulator as E
    import nibabel as nib
    original = E.make_subject

    def make(*a, **k):
        dirs, rs, keep = original(*a, **k)
        path = os.path.join(dirs['masks'], name + '.nii.gz')
        img = nib.load(path); m = np.asarray(img.dataobj) > 0
        idx = np.argwhere(m)[:n_voxels]
        small = np.zeros_like(m, dtype='u1'); small[tuple(idx.T)] = 1
        nib.save(nib.Nifti1Image(small, img.affine), path)
        return dirs, rs, keep
    monkeypatch.setattr(E, 'make_subject', make)


@pytest.mark.parametrize('method', ['AFNI', 'Grandjean', 'fMRIPrep_style'])
def test_tiny_ventricle_mask_is_skipped_and_reported(emu, monkeypatch, method):
    _shrink_mask(monkeypatch, 'Vmask', 3)
    out, keep, warn = emu.run(method, 'Skip', GS=False, Vc=True)
    assert out is not None, warn
    assert any('CSF mask has 3 voxel' in w and 'NOT used' in w for w in warn)
    assert not np.isnan(out).any()


def test_tissue_report_is_written_for_each_run(emu, monkeypatch):
    import glob
    import json
    _shrink_mask(monkeypatch, 'Vmask', 3)
    seen = {}
    original = N._tissue_too_small

    def spy(label, mask_path, root_RS, diary_file, report):
        r = original(label, mask_path, root_RS, diary_file, report); seen.update(report); return r
    monkeypatch.setattr(N, '_tissue_too_small', spy)
    emu.run('Grandjean', 'Skip', GS=False, Vc=False)
    assert seen['CSF']['voxels'] == 3 and seen['CSF']['used'] is False
    assert seen['WM']['used'] is True and seen['WM']['voxels'] > N.MIN_TISSUE_VOXELS
