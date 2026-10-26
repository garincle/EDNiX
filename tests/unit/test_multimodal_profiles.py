"""Figure XX data path (Exemples/Study_EDNiX/analysis/multimodal_profiles.py) on synthetic BIDS trees."""
import os

import numpy as np
import pandas as pd
import nibabel as nib
import pytest

pytest.importorskip('seaborn')
from Exemples.Study_EDNiX.analysis import multimodal_profiles as MP
from Plotting import Plot_interspecies_regional_profiles as PIRP

opj = os.path.join
REGIONS = ['Frontal_lobe', 'Temporal_lobe', 'Parietal_lobe', 'Occipital_lobe', 'Insular_lobe']
EXTRA = ['Lateral_ventricules', 'Brainstem']          # not cortical; ventricles excluded from the brain volume
NAMES = REGIONS + EXTRA


def _labels(atlas_lib, sp):
    d = opj(atlas_lib, 'atlas', 'mammals', 'primates', sp, 'label_code')
    os.makedirs(d, exist_ok=True)
    with open(opj(d, 'EDNIxCSCLR_StatsLUT.txt'), 'w') as f:
        for i, n in enumerate(NAMES):
            f.write(f'{10 + i} L_{n} 10 10 10 0\n{1010 + i} R_{n} 10 10 10 0\n')


def _dseg(bids, sub, sizes):
    """4D dseg: level 1 (index 0) one label; level 2 (index 1) the regions with `sizes` voxels per hemisphere."""
    d = opj(bids, f'sub-{sub}', 'ses-1', 'anat', 'native', 'volumes', 'labels')
    os.makedirs(d, exist_ok=True)
    n = 2 * sum(sizes)
    lab = np.zeros((n, 1, 1, 2), np.int32)
    lab[..., 0] = 1
    k = 0
    for i, s in enumerate(sizes):
        for off in (0, 1000):
            lab[k:k + s, 0, 0, 1] = 10 + i + off
            k += s
    nib.save(nib.Nifti1Image(lab, np.eye(4)), opj(d, f'sub-{sub}_seg-EDNIxCSCLR_dseg.nii.gz'))


def _matrix(bids, sub, strength, rng):
    d = opj(bids, f'sub-{sub}', 'ses-1', 'func', 'acpc-func', 'Stats', 'Correl_matrix', 'EDNIxCSCLR', 'correlation')
    os.makedirs(d, exist_ok=True)
    names = [f'{h}_{n}' for h in 'LR' for n in NAMES]
    s = np.array([strength[n] for h in 'LR' for n in NAMES])
    r = np.clip(np.sqrt(np.outer(s, s)) + rng.normal(0, 0.01, (len(s), len(s))), -0.95, 0.95)
    r = (r + r.T) / 2
    np.fill_diagonal(r, 0)
    for run in (0, 1):
        pd.DataFrame(r, index=names, columns=names).to_csv(opj(d, f'EDNIxCSCLR_1_run_{run}_correlation_matrix.csv'))


def _pet(bids, sub, values):
    d = opj(bids, f'sub-{sub}', 'ses-1', PIRP.ROI_SUBDIR)
    os.makedirs(d, exist_ok=True)
    rows = [dict(label=10 + i + off, region=f'{h}_{n}', n_voxels=100, n_valid=100, volume_mm3=800.0,
                 mean=values[n], sd=0.1, median=values[n])
            for i, n in enumerate(NAMES) for h, off in (('L', 0), ('R', 1000))]
    pd.DataFrame(rows).to_csv(opj(d, f'sub-{sub}_ses-1_trc-FDG_seg-EDNIxCSCLR-lvl1_desc-bratio_ROIs.tsv'),
                              sep='\t', index=False)


@pytest.fixture
def cfg(tmp_path):
    rng = np.random.default_rng(0)
    lib = str(tmp_path / 'lib')
    src = {'morphometry': {}, 'FC': {}, 'metabolism': {}}
    for sp, nsub, occ in (('Human', 6, 1.10), ('Macaque', 4, 1.40)):
        _labels(lib, sp)
        bids = str(tmp_path / sp)
        for s in range(nsub):
            _dseg(bids, f'{sp[0]}{s}', [40, 30, 20, 10, 5, 50, 60])
            _matrix(bids, f'{sp[0]}{s}', dict(zip(NAMES, [0.3, 0.2, 0.25, 0.4, 0.1, 0.05, 0.05])), rng)
            v = dict(zip(NAMES, [1.1, 0.95, 1.05, occ, 1.0, 0.5, 0.8]))
            _pet(bids, f'{sp[0]}{s}', {k: x + rng.normal(0, 0.02) for k, x in v.items()})
        for m in src:
            src[m][sp] = [bids]
    return dict(SOURCES=src, ATLAS='EDNIxCSCLR', LEVEL=2, ATLAS_LIB=lib, DESC=('bratio', 'bratio-mean'))


def test_modalities_are_read_at_level_2_with_canonical_names(cfg):
    t = MP.collect_modalities(cfg)
    assert set(t) == {'morphometry', 'FC', 'metabolism'}
    for df in t.values():
        assert set(REGIONS) <= set(df['region'])
    m = t['morphometry']
    h0 = m[(m.subject == 'H0')].set_index('region')['value']
    # ventricles excluded from the brain: 2 x (40+30+20+10+5+60) = 330 voxels
    assert 'Lateral_ventricules' not in h0.index
    assert np.isclose(h0['Frontal_lobe'], 80 / 330 * 100)
    assert m.groupby('species')['subject'].nunique().to_dict() == {'Human': 6, 'Macaque': 4}


def test_node_strength_is_mean_fisher_z_without_diagonal(tmp_path):
    r = np.array([[0, 0.5, 0.2], [0.5, 0, 0.1], [0.2, 0.1, 0]])
    p = tmp_path / 'm.csv'
    pd.DataFrame(r, index=list('abc'), columns=list('abc')).to_csv(p)
    s = MP.node_strength(str(p))
    assert np.isclose(s['a'], np.mean(np.arctanh([0.5, 0.2])))


def test_profiles_and_figure(cfg, tmp_path):
    from Exemples.Study_EDNiX.analysis import Figure_XX_multimodal as F
    t = MP.collect_modalities(cfg)
    F.REGIONS, F.LEVEL = None, 2
    assert F.select_regions(t) == ['Frontal_lobe', 'Insular_lobe', 'Temporal_lobe', 'Parietal_lobe',
                                   'Occipital_lobe']                      # Atlas_nomenclature order, no Brainstem
    prof_b, stats_b, cov_b = F.metabolism_stats(t['metabolism'], REGIONS)
    assert (cov_b.n_missing == 0).all()
    assert list(prof_b.columns) == REGIONS
    occ = stats_b.set_index('region').loc['Occipital_lobe']
    assert occ.q_fdr < 0.05 and occ.difference > 0.2
    prof_c, regions_c = MP.multimodal_profiles(t, REGIONS, F.SPECIES)
    assert regions_c == REGIONS
    z = prof_c['Human']['FC']
    assert np.allclose(z.mean(axis=1), 0) and np.allclose(z.std(axis=1, ddof=1), 1)
    cmc = MP.cross_modal_correlations(prof_c)
    assert len(cmc) == 6
    fig = F.build_figure({}, (0.5, 1.5), prof_b, stats_b, prof_c, regions_c)
    fig.savefig(tmp_path / 'fig.png')
    assert 'Figure XX' in F.legend_text(prof_b, stats_b, prof_c, regions_c, (0.5, 1.5), cmc)


def test_region_missing_in_one_species_is_dropped_not_the_species(cfg):
    """Regression: a region absent from every human table used to drop all humans (listwise deletion)."""
    from Exemples.Study_EDNiX.analysis import Figure_XX_multimodal as F
    met = MP.collect_modalities(cfg)['metabolism']
    met = met[~((met.species == 'Human') & (met.region == 'Insular_lobe'))]          # macaque-only region
    met.loc[(met.subject == 'M0') & (met.region == 'Parietal_lobe'), 'value'] = np.nan  # one NaN subject
    prof_b, stats_b, cov_b = F.metabolism_stats(met, REGIONS)
    assert set(prof_b.index.get_level_values(0)) == {'Human', 'Macaque'}
    assert list(prof_b.columns) == ['Frontal_lobe', 'Temporal_lobe', 'Occipital_lobe']
    miss = cov_b.set_index(['region', 'species'])['n_missing']
    assert miss[('Insular_lobe', 'Human')] == 6 and miss[('Parietal_lobe', 'Macaque')] == 1
    with pytest.raises(ValueError, match='species found'):
        PIRP.compare_species(prof_b.loc[['Macaque']], 'Human', 'permutation', 10)


def test_fingerprint_names_are_kept_and_classified_at_level_3():
    names = ['L_Orbital_frontal_cortex_(oFC)', 'R_Prefrontal_cortex_(PFC)', 'L_Thalamus', 'R_Striate_Cx']
    keys = [MP.region_key(n) for n in names]
    assert keys[0] == 'Orbital_frontal_cortex_(oFC)'                     # base name, as the Fig 5 profiles
    assert MP.cortical_regions(keys, 3) == ['Orbital_frontal_cortex_(oFC)', 'Prefrontal_cortex_(PFC)',
                                            'Striate_Cx']
