"""Figure 6 a-c: FDG datasets (library/pet_datasets.py) on synthetic ROI tables."""
import os

import numpy as np
import pandas as pd
import nibabel as nib
import pytest

pytest.importorskip('seaborn')
from Exemples.Study_EDNiX.analysis.library import pet_datasets as PD
from Plotting import Plot_interspecies_regional_profiles as PIRP

opj = os.path.join
NAMES = ['Frontal_lobe', 'Temporal_lobe', 'Parietal_lobe', 'Occipital_lobe', 'Insular_lobe']


def _write(bids, sub, ses, values, desc):
    d = opj(bids, f'sub-{sub}', f'ses-{ses}', PIRP.ROI_SUBDIR)
    os.makedirs(d, exist_ok=True)
    rows = [dict(label=10 + i + off, region=f'{h}_{n}', n_voxels=100, n_valid=100, volume_mm3=800.0,
                 mean=values[n], sd=0.1, median=values[n])
            for i, n in enumerate(NAMES) for h, off in (('L', 0), ('R', 1000))]
    pd.DataFrame(rows).to_csv(opj(d, f'sub-{sub}_ses-{ses}_trc-FDG_seg-EDNIxCSCLR-lvl1_desc-{desc}_ROIs.tsv'),
                              sep='\t', index=False)


@pytest.fixture
def datasets(tmp_path):
    rng = np.random.default_rng(0)
    base = dict(zip(NAMES, [1.10, 0.95, 1.05, 1.15, 1.00]))
    out = []
    for label, sp, nsub, sessions, desc, occ in (('Human_Merida', 'Human', 8, ['1'], 'bratio', 0.0),
                                                 ('Human_Castrillon', 'Human', 6, ['open', 'closed'], 'bratio-mean', 0.0),
                                                 ('Macaque_Tremblay', 'Macaque', 3, ['1'], 'bratio-mean', 0.3)):
        bids = str(tmp_path / label)
        for s in range(nsub):
            for ses in sessions:
                v = {k: x + (occ if k == 'Occipital_lobe' else 0) + (0.05 if ses == 'open' and k == 'Occipital_lobe'
                                                                     else 0) + rng.normal(0, 0.02)
                     for k, x in base.items()}
                _write(bids, f'{label[:2]}{s}', ses, v, desc)
        out.append((label, sp, bids))
    return out


def test_sessions_are_kept_or_averaged(datasets):
    df = PD.dataset_table(datasets, 'EDNIxCSCLR', 2)
    ds = PD.dataset_table(datasets, 'EDNIxCSCLR', 2, sessions=True)
    n = df.groupby('dataset')['subject'].nunique().to_dict()
    assert n == {'Human_Merida': 8, 'Human_Castrillon': 6, 'Macaque_Tremblay': 3}
    c = ds[ds.dataset == 'Human_Castrillon']
    assert set(c.session) == {'open', 'closed'} and c.subject.nunique() == 6
    # the subject mean is the mean of its sessions
    m = df[(df.dataset == 'Human_Castrillon') & (df.region == 'Occipital_lobe')].set_index('subject')['value']
    s = c[c.region == 'Occipital_lobe'].groupby('subject')['value'].mean()
    assert np.allclose(m.sort_index(), s.sort_index())


def test_reference_difference_and_agreement(datasets):
    df = PD.dataset_table(datasets, 'EDNIxCSCLR', 2)
    w = PD.dataset_matrix(df, NAMES, [d[0] for d in datasets])
    st = PD.compare_to_reference(w, 'Macaque_Tremblay', 'permutation', 200)
    assert set(st.dataset) == {'Human_Merida', 'Human_Castrillon'}
    occ = st[st.region == 'Occipital_lobe']
    assert (occ.difference < -0.2).all()                     # human - macaque
    ag = PD.agreement(w, 'Human_Merida', 'Human_Castrillon', n_boot=200).iloc[0]
    assert ag.pearson_r > 0.9 and 0 < ag.lin_ccc <= ag.pearson_r + 1e-9
    assert ag.pearson_r_ci_low <= ag.pearson_r <= ag.pearson_r_ci_high


def test_lin_ccc_penalises_offset():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    assert np.isclose(PD.lin_ccc(x, x), 1.0)
    assert np.isclose(np.corrcoef(x, x + 1)[0, 1], 1.0) and PD.lin_ccc(x, x + 1) < 0.75


def test_pair_classes_and_permutation(datasets):
    ds = PD.dataset_table(datasets, 'EDNIxCSCLR', 2, sessions=True)
    labels = [d[0] for d in datasets]
    ws = PD.session_matrix(ds, NAMES, labels)
    one, chosen = PD.one_session(ws, {'Human_Castrillon': 'open'})
    assert (chosen[chosen.dataset == 'Human_Castrillon'].session == 'open').all()
    assert one.shape[0] == 8 + 6 + 3
    multi = ws.loc[['Human_Castrillon']]
    pairs, test = PD.pair_similarity(one, multi, humans=('Human_Merida', 'Human_Castrillon'),
                                     reference='Macaque_Tremblay', n_perm=200)
    n = pairs.groupby('pair_class').size().to_dict()
    assert n['same subject, other session'] == 6
    assert n['within Human_Merida'] == 28 and n['within Human_Castrillon'] == 15
    assert n['between human datasets'] == 48 and n['human vs macaque'] == (8 + 6) * 3
    assert len(test) == 1 and 0 < test.p_permutation.iat[0] <= 1


def test_paint_regions(tmp_path):
    lab = np.zeros((4, 1, 1, 2), np.int32)
    lab[:2, 0, 0, 1], lab[2:, 0, 0, 1] = 10, 1013          # L_Frontal, R_Occipital at level 2 (index 1)
    f = str(tmp_path / 'a.nii.gz')
    nib.save(nib.Nifti1Image(lab, np.eye(4)), f)
    ref = nib.Nifti1Image(np.ones((4, 1, 1), np.float32), np.eye(4))
    ldf = pd.DataFrame({'label_id': [10, 1013], 'region_name': ['L_Frontal_lobe', 'R_Occipital_lobe']})
    img = PD.paint_regions(pd.Series({'frontal_lobe': 0.2, 'Occipital_lobe': -0.1}), f, 1, ref, ldf,
                           mask=np.array([True, True, True, False]).reshape(4, 1, 1))
    v = np.asarray(img.dataobj).ravel()
    assert np.allclose(v, [0.2, 0.2, -0.1, 0.0])
