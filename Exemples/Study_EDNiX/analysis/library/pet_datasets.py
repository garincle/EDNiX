"""
FDG-PET at the level of the DATASET (Figure 6 a-c): two human datasets (CERMEP / Merida, static 10-min frame; Castrillon,
dynamic simultaneous PET/MR, late frames 20-45 min) and one macaque dataset (Tremblay), on the regions of the EDNiX atlas.

  dataset_table()          long table dataset, species, subject, session, region, value (bratio), one row per session
                           (session=None: sessions of a subject averaged, the unit of the group maps of panel a)
  dataset_matrix()         wide table index (dataset, subject), columns regions
  compare_to_reference()   region-wise difference of every dataset with the reference dataset (macaque), permutation
                           test, BH-FDR within each dataset (PIRP.compare_species)
  agreement()              agreement of the regional group means of two datasets: Pearson / Spearman r, Lin's
                           concordance (CCC, absolute agreement: both are bratio), 95% CI by bootstrap of the subjects,
                           split-half reliability of each dataset and disattenuated r (noise ceiling)
  pair_similarity()        Pearson r between the regional profiles of every pair of subjects, in classes: same subject /
                           other session, between subjects within a dataset, between the two human datasets, human vs
                           macaque; permutation test of the dataset labels (within - between)
  paint_regions()          regional values painted into the template atlas (surface of the difference map of panel a)

Region values are those of Plot_interspecies_regional_profiles (step 7 of the PET pipeline): left and right labels
merged by valid-voxel weighting; a region with < min_voxels valid voxels in one hemisphere label of a session is NaN.
"""
import os

import numpy as np
import pandas as pd
import nibabel as nib
from scipy import stats as sstats

from Exemples.Study_EDNiX.analysis.library import multimodal_profiles as MP
from Plotting import Plot_interspecies_regional_profiles as PIRP

opj = os.path.join
opb = os.path.basename

SEP = '|'                                       # subject|session key used to keep sessions apart in PIRP.subject_level


def dataset_table(datasets, atlas, level, desc=('bratio', 'bratio-mean'), min_voxels=10, sessions=False):
    """datasets: [(label, species, bids_dir), ...]. Long table dataset, species, subject, session, region, value.
    sessions=False: one row per subject (sessions averaged, as pet_group.run_group_mean); True: one row per session."""
    rows = []
    for label, sp, bd in datasets:
        df = PIRP.collect({label: bd}, atlas, level - 1, desc)
        df['session'] = df['session'].fillna('1')
        if sessions:
            df['subject'] = df['subject'].astype(str) + SEP + df['session'].astype(str)
        sub = PIRP.subject_level(df, 'merge')
        sub.loc[sub['min_label_voxels'] < min_voxels, 'value'] = np.nan
        sub['region'] = sub['region_key'].map(MP.region_key)
        if sessions:
            parts = sub['subject'].str.split(SEP, n=1, expand=True)
            sub['subject'], sub['session'] = parts[0], parts[1]
        else:
            sub['session'] = 'mean'
        sub['dataset'], sub['species'], sub['bids_dir'] = label, sp, opb(os.path.normpath(bd))
        rows.append(sub[['dataset', 'species', 'bids_dir', 'subject', 'session', 'region', 'value']])
    return pd.concat(rows, ignore_index=True)


def complete_regions(df, regions, label='metabolism'):
    """Regions with a value in every subject of every dataset (MP.complete_regions with dataset as the group)."""
    d = df.rename(columns={'species': 'species_name', 'dataset': 'species'})
    return MP.complete_regions(d, regions, list(dict.fromkeys(d['species'])), label)


def dataset_matrix(df, regions, datasets=None):
    """Wide table index (dataset, subject), columns regions (in that order), datasets in the order given."""
    w = df[df['region'].isin(regions)].pivot_table(index=['dataset', 'subject'], columns='region', values='value')
    w = w.reindex(columns=list(regions)).dropna(axis=0, how='any')
    if datasets:
        w = w.reindex([d for d in datasets if d in w.index.get_level_values(0)], level=0)
    return w


def compare_to_reference(w, reference, test='permutation', n_perm=10000, seed=0):
    """Region-wise test of every dataset against `reference` (difference = dataset - reference; BH-FDR over the
    regions, within each dataset). With n = 2 macaques a permutation test has C(n + 2, 2) distinct relabellings
    (e.g. 231 for 19 humans: smallest P ~ 0.004): the difference map is descriptive."""
    return PIRP.compare_species(w, reference, test, n_perm, seed).rename(
        columns={'species': 'dataset', 'n_species': 'n_dataset', 'mean_species': 'mean_dataset',
                 'sem_species': 'sem_dataset'})


def lin_ccc(x, y):
    """Lin's concordance correlation coefficient (Lin 1989, Biometrics 45:255): agreement on the identity line."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    mx, my = x.mean(), y.mean()
    vx, vy = x.var(), y.var()
    cov = ((x - mx) * (y - my)).mean()
    return 2 * cov / (vx + vy + (mx - my) ** 2)


def agreement(w, a, b, n_boot=5000, seed=0):
    """Agreement of the regional group-mean profiles of datasets a and b (subjects x regions in w): Pearson, Spearman,
    Lin's CCC, mean absolute difference; 95% CI by resampling the subjects of each dataset; split-half (Spearman-Brown)
    reliability of each dataset and disattenuated Pearson r = r / sqrt(rel_a rel_b) (noise ceiling)."""
    rng = np.random.default_rng(seed)
    A, B = w.loc[a].values, w.loc[b].values
    ma, mb = A.mean(0), B.mean(0)

    def all_(x, y):
        return (np.corrcoef(x, y)[0, 1], sstats.spearmanr(x, y).statistic, lin_ccc(x, y), np.mean(np.abs(x - y)))

    est = all_(ma, mb)
    boot = np.array([all_(A[rng.integers(0, len(A), len(A))].mean(0), B[rng.integers(0, len(B), len(B))].mean(0))
                     for _ in range(n_boot)])
    lo, hi = np.nanpercentile(boot, 2.5, axis=0), np.nanpercentile(boot, 97.5, axis=0)
    rel_a = PIRP.split_half_reliability(A, rng=rng)
    rel_b = PIRP.split_half_reliability(B, rng=rng)
    ceil = np.sqrt(rel_a * rel_b) if np.isfinite(rel_a) and np.isfinite(rel_b) and rel_a > 0 and rel_b > 0 else np.nan
    out = dict(dataset_a=a, dataset_b=b, n_a=len(A), n_b=len(B), n_regions=w.shape[1])
    for k, name in enumerate(('pearson_r', 'spearman_r', 'lin_ccc', 'mean_abs_difference')):
        out.update({name: est[k], f'{name}_ci_low': lo[k], f'{name}_ci_high': hi[k]})
    out.update(reliability_a=rel_a, reliability_b=rel_b, noise_ceiling=ceil,
               pearson_r_disattenuated=est[0] / ceil if np.isfinite(ceil) else np.nan)
    return pd.DataFrame([out])


def _r(x, y):
    return float(np.corrcoef(x, y)[0, 1])


def pair_similarity(subj, sess=None, humans=(), reference=None, n_perm=10000, seed=0):
    """Pearson r between the regional profiles of pairs of subjects.
    subj : wide table (dataset, subject) x regions, ONE session per subject (so that every subject profile has the same
           noise; a mean of two sessions is more reliable than a single session and would inflate its class);
    sess : wide table (dataset, subject|session) x regions of the datasets with repeated sessions (same subject, other
           session class); None skips it.
    humans: the two human datasets (classes 'within <dataset>' and 'between human datasets'); reference: the macaque
           dataset (class 'human vs macaque').
    Returns (pairs, test): pairs = long table class, a, b, r; test = permutation of the dataset labels of the human
    subjects (sizes kept): statistic = mean r within datasets - mean r between datasets, two-sided P. The pairs share
    subjects (not independent): only the label permutation gives a valid P (Mantel-like)."""
    rng = np.random.default_rng(seed)
    rows = []
    if sess is not None and len(sess):
        idx = sess.index.to_frame(index=False)
        idx[['subject', 'session']] = idx['subject'].str.split(SEP, n=1, expand=True)
        for (ds, sb), g in idx.groupby(['dataset', 'subject'], sort=False):
            if len(g) < 2:
                continue
            v = [sess.iloc[i].values for i in g.index]
            for i in range(len(v)):
                for j in range(i + 1, len(v)):
                    rows.append(dict(pair_class='same subject, other session', dataset_a=ds, dataset_b=ds,
                                     a=f'{sb}:{g.session.iloc[i]}', b=f'{sb}:{g.session.iloc[j]}', r=_r(v[i], v[j])))
    keys = list(subj.index)
    X = subj.values
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            da, db = keys[i][0], keys[j][0]
            if da == db:
                cls = f'within {da}'
            elif da in humans and db in humans:
                cls = 'between human datasets'
            elif reference in (da, db) and (da in humans or db in humans):
                cls = 'human vs macaque'
            else:
                continue
            rows.append(dict(pair_class=cls, dataset_a=da, dataset_b=db, a=keys[i][1], b=keys[j][1], r=_r(X[i], X[j])))
    pairs = pd.DataFrame(rows)
    # permutation of the human dataset labels
    hm = [k for k, (d, _) in enumerate(keys) if d in humans]
    test = pd.DataFrame()
    if len(set(keys[k][0] for k in hm)) == 2:
        C = np.corrcoef(X[hm])
        lab = np.array([keys[k][0] for k in hm])
        iu = np.triu_indices(len(hm), 1)
        cv = C[iu]

        def stat(l_):
            same = l_[iu[0]] == l_[iu[1]]
            return cv[same].mean() - cv[~same].mean()

        obs = stat(lab)
        null = np.array([stat(rng.permutation(lab)) for _ in range(n_perm)])
        p = (np.sum(np.abs(null) >= abs(obs) - 1e-12) + 1) / (n_perm + 1)
        test = pd.DataFrame([dict(test='within - between human datasets (mean pairwise r)', statistic=obs,
                                  null_mean=null.mean(), null_sd=null.std(ddof=1), p_permutation=p, n_perm=n_perm,
                                  n_subjects=len(hm))])
    return pairs, test


def one_session(sess_w, prefer=None):
    """From a wide table (dataset, subject|session) x regions, one row per subject: the session `prefer`
    ({dataset: session}) when it exists, otherwise the first session (alphabetical)."""
    prefer = prefer or {}
    idx = sess_w.index.to_frame(index=False)
    idx[['sub', 'ses']] = idx['subject'].str.split(SEP, n=1, expand=True)
    keep = []
    for (ds, sb), g in idx.groupby(['dataset', 'sub'], sort=False):
        want = prefer.get(ds)
        hit = g[g.ses == want] if want is not None else g.iloc[:0]
        keep.append((hit if len(hit) else g.sort_values('ses')).index[0])
    out = sess_w.iloc[keep].copy()
    out.index = pd.MultiIndex.from_arrays([idx.loc[keep, 'dataset'].values, idx.loc[keep, 'sub'].values],
                                          names=['dataset', 'subject'])
    chosen = idx.loc[keep, ['dataset', 'sub', 'ses']].rename(columns={'sub': 'subject', 'ses': 'session'})
    return out, chosen.reset_index(drop=True)


def session_matrix(df_sess, regions, datasets=None):
    """Wide table (dataset, subject|session) x regions from dataset_table(..., sessions=True)."""
    d = df_sess.copy()
    d['subject'] = d['subject'].astype(str) + SEP + d['session'].astype(str)
    return dataset_matrix(d, regions, datasets)


def paint_regions(values, atlas_file, level_index, ref_img, label_df, mask=None):
    """Volume on the grid of ref_img (nibabel image) where every voxel of a label takes the value of its region
    (values: Series region -> value, region = base name, compared case-insensitively); 0 elsewhere and outside mask
    (bool volume on the same grid). level_index: 0-based volume of a 4D atlas."""
    from Statistics.Group_PET import pet_group
    lab = pet_group.labels_on_grid(atlas_file, level_index, ref_img)
    names = dict(zip(label_df['label_id'].astype(int), label_df['region_name']))
    v = {str(k).lower().strip(): float(x) for k, x in values.items() if np.isfinite(x)}
    out = np.zeros(lab.shape, np.float32)
    for i in np.unique(lab[lab > 0]):
        key = MP.region_key(names.get(int(i), '')).lower().strip()
        if key in v:
            out[lab == i] = v[key]
    if mask is not None:
        out[~mask] = 0
    return nib.Nifti1Image(out, ref_img.affine)
