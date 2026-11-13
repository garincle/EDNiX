"""Supp. Fig. 8 (library/longitudinal_ablation.py) on synthetic longitudinal data."""
import numpy as np
import pandas as pd
import pytest

pytest.importorskip('statsmodels')
from Exemples.Study_EDNiX.analysis.library import longitudinal_ablation as LA
from Statistics.Group_fMRI.EDNiX_longitudinal import fit_lme


def _morph(noise, rng, n_sub=10, n_ses=4, slope=-5.0, ref_noise=0.0, base=None):
    """Left-isocortex rows in the format of longitudinal.load_morphometry. The latest
    session (reference) gets `ref_noise`, the others `noise`."""
    rows, ages = [], []
    base = rng.normal(1000, 50, n_sub) if base is None else base
    for i in range(n_sub):
        a0 = 70 + 2 * i
        for s in range(1, n_ses + 1):
            age = a0 + 1.5 * (s - 1)
            e = rng.normal(0, ref_noise if s == n_ses else noise)
            rows.append(dict(subject=str(i), session=str(s), region='Isocortex', hemisphere='left',
                             atlas_level=1, surface_area_mm2=base[i] + slope * (age - 70) + e))
            ages.append(dict(subject=str(i), session=str(s), age=age))
    return pd.DataFrame(rows), pd.DataFrame(ages).drop_duplicates()


def test_lmm_matches_statsmodels():
    rng = np.random.default_rng(1)
    df, age = _morph(3.0, rng)
    m = df.merge(age, on=['subject', 'session'])
    ours = LA.lmm_random_intercept(m['surface_area_mm2'], m['age'], m['subject'])
    ref = fit_lme(m, 'surface_area_mm2', 'age')
    assert ours['slope'] == pytest.approx(ref['slope_per_unit'], rel=1e-4)
    assert ours['se'] == pytest.approx(ref['se_per_unit'], rel=2e-2)


def test_compare_modes_detects_noise_reduction_and_no_bias():
    rng = np.random.default_rng(2)
    base = rng.normal(1000, 50, 12)
    long_df, age = _morph(1.0, rng, n_sub=12, base=base)
    cross_df, _ = _morph(6.0, rng, n_sub=12, base=base)
    # reference session identical in both modes, as in the pipeline
    last = long_df['session'] == '4'
    cross_df.loc[last, 'surface_area_mm2'] = long_df.loc[last, 'surface_area_mm2'].to_numpy()
    morph = {LA.MODE_LONG: {'surface': long_df}, LA.MODE_CROSS: {'surface': cross_df}}
    m = LA.matched_sessions(morph, age, 'age', 'surface', 'surface_area_mm2')
    assert m['is_reference'].sum() == 12
    r = LA.compare_modes(m, 'age', n_boot=300)
    assert r['ref_abs_diff_pct_median'] == pytest.approx(0.0, abs=1e-9)
    assert r['sigma_e_ratio'] > 2 and r['sigma_e_ratio_ci_low'] > 1 and r['sigma_e_ratio_p'] < 0.05
    assert r['se_ratio'] > 1
    assert r['delta_slope_ci_low'] < 0 < r['delta_slope_ci_high']      # same true slope
    assert r['n_subjects_less_noisy_long'] >= 9


def test_find_folder_requires_a_single_match(tmp_path):
    (tmp_path / 'BIDS_Park_nolong').mkdir()
    assert LA.find_folder(str(tmp_path / 'BIDS_Park*no*long*')).endswith('BIDS_Park_nolong')
    (tmp_path / 'BIDS_Park_no_long_old').mkdir()
    with pytest.raises(FileNotFoundError):
        LA.find_folder(str(tmp_path / 'BIDS_Park*no*long*'))
