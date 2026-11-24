"""Fig. 3 a-f composites (library/longitudinal.py): primary / association cortex per species."""
import numpy as np
import pandas as pd
import pytest

pytest.importorskip('statsmodels')
from Paper.analysis.library import longitudinal as LG


def _morph():
    rows = []
    for sub in ('a', 'b'):
        for ses in ('1', '2'):
            for i, reg in enumerate(LG.PRIMARY_REGIONS + LG.ASSOCIATION_REGIONS + ['Extra_striate_Cx']):
                for hemi in ('left', 'right'):
                    rows.append(dict(subject=sub, session=ses, region=reg, hemisphere=hemi, atlas_level=3,
                                     surface_area_mm2=100.0 * (i + 1), volume_mm3=10.0 * (i + 1),
                                     thickness_mm=1.0 + i))
            rows.append(dict(subject=sub, session=ses, region='Isocortex', hemisphere='left',
                             atlas_level=1, surface_area_mm2=1e4, volume_mm3=1e3, thickness_mm=2.0))
    df = pd.DataFrame(rows)
    return {k: df[['subject', 'session', 'region', 'hemisphere', 'atlas_level', c]]
            for k, c in (('surface', 'surface_area_mm2'), ('volume', 'volume_mm3'),
                         ('thickness', 'thickness_mm'))}


def test_sum_and_weighted_thickness():
    m = _morph()
    s, info = LG.composite_region_rows(m, 'surface', 'Primary')
    assert info['level'] == 3 and info['hemisphere'] == 'left' and not info['missing']
    assert info['matched'] == LG.PRIMARY_REGIONS                 # not Extra_striate_Cx
    assert np.allclose(s['surface_area_mm2'], 100 + 200 + 300 + 400)
    t, info = LG.composite_region_rows(m, 'thickness', 'Association')
    w = np.array([500., 600.]); x = np.array([5., 6.])
    assert info['weighting'] == 'surface-weighted mean'
    assert np.allclose(t['thickness_mm'], (w * x).sum() / w.sum())
    assert set(t['region']) == {'Association'} and len(t) == 4


def test_incomplete_session_dropped():
    m = _morph()
    m['volume'] = m['volume'][~((m['volume'].subject == 'a') & (m['volume'].session == '2')
                                & (m['volume'].region == 'PFC'))]
    v, _ = LG.composite_region_rows(m, 'volume', 'Association')
    assert len(v) == 3


def test_main_and_supp_sets_are_complementary():
    for sp in ('Human', 'Macaque'):
        assert {LG.FIG3_MAIN_REGIONS[sp], LG.FIG3_SUPP_REGIONS[sp]} == {'Primary', 'Association'}


def _long_morph(driver=None, n_sub=12, n_ses=4, seed=0):
    """Left level-3 regions, linear % change per year: association -2 %/yr, primary 0,
    except `driver` (a primary region) at -6 %/yr."""
    rng = np.random.default_rng(seed)
    rows, ages = [], []
    for i in range(n_sub):
        base = rng.normal(1000, 50)
        for s in range(1, n_ses + 1):
            age = 60 + 2 * i + s
            ages.append(dict(subject=str(i), session=str(s), age=age))
            for reg in LG.PRIMARY_REGIONS + LG.ASSOCIATION_REGIONS:
                r = -0.02 if reg in LG.ASSOCIATION_REGIONS else (-0.06 if reg == driver else 0.0)
                v = base * np.exp(r * (age - 60)) * rng.normal(1, 0.002)
                rows.append(dict(subject=str(i), session=str(s), region=reg, hemisphere='left',
                                 atlas_level=3, surface_area_mm2=v))
    return {'surface': pd.DataFrame(rows)}, pd.DataFrame(ages)


def test_region_stats_annual_change_and_loro():
    from Paper.analysis.library import longitudinal_regions as LR
    morph, age = _long_morph(driver='Striate_Cx')
    st = LR.region_stats([('Human', morph, age, 'age', 'Human')],
                         modalities=[LG.FIG3_MODALITIES[0]])
    t = st['table'].set_index('estimate')
    assert t.loc['PFC', 'pct_per_year'] == pytest.approx(100 * np.expm1(-0.02), abs=0.1)
    assert t.loc['Striate_Cx', 'pct_per_year'] == pytest.approx(100 * np.expm1(-0.06), abs=0.1)
    assert abs(t.loc['Primary without Striate_Cx', 'pct_per_year']) < 0.1      # driver removed
    assert t.loc['Primary', 'pct_per_year'] < -0.5                        # driven by V1
    assert any('OVERLAP' in l for l in st['lines'])
    assert any('sign change without' in l and 'Striate (V1)' in l for l in st['lines'])
    assert 'Association without PFC' not in t.index      # 1 region left = PPC alone
