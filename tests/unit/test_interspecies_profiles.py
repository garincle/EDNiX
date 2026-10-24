"""Cross-species regional profiles (Plotting/Plot_interspecies_regional_profiles.py) on synthetic ROI tables."""
import os

import numpy as np
import pandas as pd
import pytest

from Plotting import Plot_interspecies_regional_profiles as P

opj = os.path.join

REGIONS = ['Frontal', 'Motor', 'Somatosensory', 'Parietal', 'Temporal', 'Visual', 'Cingulate', 'Insula']


def _write(root, sub, ses, run, values, nvox, desc, extra=None):
    d = opj(root, 'derivatives', f'sub-{sub}', f'ses-{ses}', P.ROI_SUBDIR)
    os.makedirs(d, exist_ok=True)
    rows = []
    for i, r in enumerate(REGIONS):
        for h, off in (('L', 0), ('R', 1000)):
            rows.append(dict(label=10 + i + off, region=f'{h}_{r}', n_voxels=nvox[r], volume_mm3=nvox[r] * 8.0,
                             n_valid=nvox[r], mean=values[r] + (0.01 if h == 'R' else -0.01), sd=0.1,
                             median=values[r]))
    for lab, name, v in (extra or []):
        rows.append(dict(label=lab, region=name, n_voxels=200, volume_mm3=1600.0, n_valid=200, mean=v, sd=0.1,
                         median=v))
    run_ent = f'_run-{run}' if run else ''
    pd.DataFrame(rows).to_csv(opj(d, f'sub-{sub}_ses-{ses}_trc-FDG{run_ent}_seg-EDNIxCSCLR-lvl1_desc-{desc}_ROIs.tsv'),
                              sep='\t', index=False)


@pytest.fixture
def study(tmp_path):
    rng = np.random.default_rng(1)
    base = dict(zip(REGIONS, [1.10, 1.05, 1.00, 1.08, 0.95, 1.20, 1.02, 0.98]))
    hum, mac = tmp_path / 'Human', tmp_path / 'Macaque'
    for s in range(8):
        v = {r: base[r] + rng.normal(0, 0.02) for r in REGIONS}
        _write(hum, f'H{s:02d}', '1', None, v, {r: 500 for r in REGIONS}, 'bratio',
               extra=[(99, 'HumanOnly', 1.0)])
    for s in range(6):
        nv = {r: 200 for r in REGIONS}
        if s == 0:
            nv['Insula'] = 5                       # too small in one macaque -> region dropped
        for run in (1, 2) if s == 0 else (None,):
            v = {r: base[r] + rng.normal(0, 0.02) for r in REGIONS}
            v['Visual'] += 0.25                   # true species difference
            _write(mac, f'M{s:02d}', '1', run, v, nv, 'bratio-mean')
    return {'Human': str(hum), 'Macaque': str(mac)}


def test_collect_finds_each_species_desc(study):
    df = P.collect(study, level=1)
    assert set(df.groupby('species')['desc'].first().items()) == {('Human', 'bratio'), ('Macaque', 'bratio-mean')}
    assert df.loc[df.species == 'Macaque', 'file'].nunique() == 7


def test_subject_level_merges_runs_and_hemispheres(study):
    sub = P.subject_level(P.collect(study, level=1))
    m0 = sub[(sub.subject == 'M00')]
    assert m0['n_runs'].eq(2).all()
    assert set(m0['region_key']) == set(REGIONS)          # L_/R_ merged
    assert sub.groupby('species')['subject'].nunique().to_dict() == {'Human': 8, 'Macaque': 6}


def test_filter_keeps_common_large_regions_in_atlas_order(study):
    sub = P.subject_level(P.collect(study, level=1))
    kept, rep = P.filter_regions(sub, min_voxels=10)
    assert kept == [r for r in REGIONS if r != 'Insula']   # Insula: 5 voxels in M00; HumanOnly: one species
    assert not rep.loc[rep.region_key == 'HumanOnly', 'present_in_all_species'].any()


def test_size_filter_uses_fwhm(study):
    sub = P.subject_level(P.collect(study, level=1))
    # macaque regions: 2 x 200 voxels x 8 mm3 = 3200 mm3 -> side 14.7 mm; FWHM 8 mm -> 1.8 FWHM < 2
    kept, rep = P.filter_regions(sub, 10, fwhm_mm={'Human': 4.0, 'Macaque': 8.0}, min_size_fwhm=2)
    assert kept == []
    assert rep['size_in_fwhm'].notna().all()


def test_stats_detect_the_planted_difference(study):
    sub = P.subject_level(P.collect(study, level=1))
    kept, _ = P.filter_regions(sub)
    st = P.compare_species(P.profile_matrix(sub, kept), 'Human', 'permutation', 2000)
    sig = set(st.loc[st.q_fdr < 0.05, 'region'])
    assert 'Visual' in sig and len(sig) <= 2
    v = st.set_index('region').loc['Visual']
    assert abs(v.difference - 0.25) < 0.05 and v.hedges_g > 3


def test_normalisation_zscore_and_reference(study):
    sub = P.subject_level(P.collect(study, level=1))
    z = P.normalise(sub, 'zscore')
    assert np.allclose(z.groupby(['species', 'subject'])['value'].mean(), 0, atol=1e-12)
    r = P.normalise(sub, 'reference', ['Motor'])
    assert np.allclose(r.loc[r.region_key == 'Motor', 'value'], 1.0)


def test_bh_fdr_matches_statsmodels():
    from statsmodels.stats.multitest import multipletests
    p = np.array([0.001, 0.02, 0.03, 0.2, 0.5, np.nan])
    q = P.bh_fdr(p)
    assert np.allclose(q[:5], multipletests(p[:5], method='fdr_bh')[1]) and np.isnan(q[5])


def test_run_writes_tables_figures_and_legends(study, tmp_path):
    out = tmp_path / 'out'
    res = P.run(study, str(out), dict(N_PERM=500, N_BOOT=200))
    tag = 'EDNIxCSCLR-lvl1'
    for f in (f'profiles_{tag}.xlsx', f'methods_{tag}.txt'):
        assert (out / f).is_file()
    for fig in ('F1_fingerprint', 'F2_bars', 'F3_scatter'):
        for ext in ('png', 'svg', 'txt'):
            assert (out / f'{fig}_{tag}.{ext}').is_file()
    sim = res[1]['similarity'].iloc[0]
    assert sim.pearson_ci_low <= sim.pearson_r <= sim.pearson_ci_high
    assert 0 < sim.reliability_ref <= 1


def test_run_rejects_unknown_config(study, tmp_path):
    with pytest.raises(ValueError, match='unknown config'):
        P.run(study, str(tmp_path), {'NPERM': 10})
