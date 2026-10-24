"""
Regional group analysis of PET images: the regional values computed by step 7 of the PET pipeline in the
PET space (on the measured voxels, atlases brought into the PET space), for every subject.

    <BIDS>/sub-<ID>/ses-<ses>/pet/orig/postprocessing/ROI_analysis/
        <base>_seg-<atlas>-lvl<k>_desc-<desc>_ROIs.tsv     label, region, n_voxels, n_valid, mean, sd, median

This is the EDNiX version of the earlier ROI scripts (3_correct_atlas_for_small_V.py, 4_signal_to_excel.py):
 1. one table subjects x regions (mean of every region);
 2. regions with fewer than `min_voxels` PET voxels in at least one subject are dropped (the earlier rule:
    more than 10 voxels in 100 % of the subjects): a region of a few PET voxels is dominated by the partial
    volume effect and by the registration error. The criterion is in PET voxels, so it follows the
    resolution of each scanner; it is also reported in mm3;
 3. normalisation (SUVR): division by the mean of a reference region, i.e. the voxel-weighted mean of its
    labels (union of the voxels, as the reference TAC of the kinetic models; the earlier script averaged
    the left and right means without weights), or by the voxel-weighted mean of all the regions of the
    atlas level ('global'); 'none' keeps SUV;
 4. per region: group means and SD, GLM (same designs as the voxel-wise analysis), two-sided p and
    Benjamini-Hochberg FDR over the regions (Benjamini & Hochberg 1995);
 5. tables (TSV and Excel, like all_animals.xlsx / all_animals_norm.xlsx) and a figure of the effects.
"""
import os
import glob

import numpy as np
import pandas as pd

from Statistics.Group_PET import _glm
from Statistics.Group_PET.pet_group import _entities, _filter, subject_table, _log

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname


def collect_roi_tables(bids_dir, atlas, level, desc, list_to_keep=None, list_to_remove=None):
    """Long table: subject, session, run, base, label, region, n_voxels, n_valid, mean (one row per image and
    region). desc: desc of the regional table ('SUV', 'SUV-mean', 'SUV-PVC', 'bratio'...)."""
    name = f'*_seg-{atlas}-lvl{level}_desc-{desc}_ROIs.tsv'
    pats = [opj(bids_dir, 'sub-*', 'ses-*', 'pet', 'orig', 'postprocessing', 'ROI_analysis', name),
            opj(bids_dir, 'sub-*', 'pet', 'orig', 'postprocessing', 'ROI_analysis', name)]
    rows = []
    for f in sorted(set(sum((glob.glob(p) for p in pats), []))):
        ent = _entities(f)
        t = pd.read_csv(f, sep='\t')
        base = opb(f).split('_seg-')[0]
        t.insert(0, 'base', base)
        t.insert(0, 'run', ent.get('run', ''))
        t.insert(0, 'session', ent.get('ses', ''))
        t.insert(0, 'subject', ent.get('sub'))
        rows.append(t)
    if not rows:
        raise FileNotFoundError(f'no regional table {name} in {bids_dir}: run step 7 of the PET pipeline '
                                f'(atlas, level and desc as in the file names of ROI_analysis/)')
    long = pd.concat(rows, ignore_index=True)
    keys = _filter(long[['subject', 'session']].drop_duplicates().reset_index(drop=True), list_to_keep,
                   list_to_remove)
    long = long.merge(keys, on=['subject', 'session'])
    if 'n_valid' not in long.columns:
        long['n_valid'] = long['n_voxels']
    return long


def session_values(long):
    """Runs of one session combined (voxel-weighted mean of the regional means)."""
    def comb(g):
        w = g['n_valid'].to_numpy(float)
        m = g['mean'].to_numpy(float)
        ok = np.isfinite(m) & (w > 0)
        return pd.Series(dict(region=g['region'].iloc[0], n_voxels=g['n_voxels'].min(), n_valid=g['n_valid'].min(),
                              volume_mm3=g['volume_mm3'].min() if 'volume_mm3' in g else np.nan,
                              mean=float(np.sum(w[ok] * m[ok]) / np.sum(w[ok])) if ok.any() else np.nan))
    return long.groupby(['subject', 'session', 'label'], sort=True).apply(comb, include_groups=False).reset_index()


def filter_regions(ses, min_voxels=10):
    """Regions with at least min_voxels valid PET voxels in every session. Returns (kept table, report)."""
    by = ses.groupby('label').agg(region=('region', 'first'), min_voxels=('n_valid', 'min'),
                                  min_volume_mm3=('volume_mm3', 'min'), n_sessions=('n_valid', 'size'))
    n_ses = ses[['subject', 'session']].drop_duplicates().shape[0]
    by['kept'] = (by['min_voxels'] >= min_voxels) & (by['n_sessions'] == n_ses)
    kept = set(by.index[by['kept']])
    return ses[ses['label'].isin(kept)].reset_index(drop=True), by.reset_index()


def normalise_regions(ses, normalization, all_regions=None):
    """
    Divides the regional means of every session by its reference value. normalization:
    'none'; 'global' (voxel-weighted mean of all the regions of `all_regions`, i.e. before the filter: the
    atlas coverage of the brain); a list of labels or {'labels': [...]} (reference region, union of the
    labels, voxel-weighted). Returns (table with 'value' = normalised mean, reference values per session).
    """
    out = ses.copy()
    if normalization == 'none':
        out['value'] = out['mean']
        return out, None
    src = all_regions if all_regions is not None else ses
    if normalization == 'global':
        sel = src
    else:
        labels = normalization.get('labels') if isinstance(normalization, dict) else normalization
        labels = [int(l) for l in np.atleast_1d(labels)]
        sel = src[src['label'].isin(labels)]
        missing = set(labels) - set(sel['label'])
        if missing:
            raise ValueError(f'reference label(s) {sorted(missing)} absent from the regional tables')
    ok = sel[np.isfinite(sel['mean']) & (sel['n_valid'] > 0)]
    ref = ok.groupby(['subject', 'session']).apply(
        lambda g: float(np.sum(g['n_valid'] * g['mean']) / np.sum(g['n_valid'])), include_groups=False)
    ref = ref.rename('normalising_value').reset_index()
    out = out.merge(ref, on=['subject', 'session'], how='left')
    out['value'] = out['mean'] / out['normalising_value']
    return out, ref


def run_roi_analysis(bids_dir, atlas, level, desc, design, group_col=None, groups=None, covariates=None,
                     test_covariate=None, paired_col='session', paired_levels=None, normalization='global',
                     min_voxels=10, aggregate='subject', list_to_keep=None, list_to_remove=None, output_dir=None,
                     output_suffix=None, fdr_q=0.05):
    """
    Regional group analysis (module docstring). Parameters as pet_group.run_group_analysis, plus:
    atlas, level : atlas and level of the regional tables (file names of ROI_analysis/, e.g. 'EDNIxCSCLR', 0)
    desc         : desc of the tables ('SUV', 'SUV-mean', 'SUV-PVC'...)
    normalization: 'global' | 'none' | list of labels of this atlas level (e.g. [7, 1007]: cerebellum)
    min_voxels   : minimum number of valid PET voxels of a region in every session (10: earlier rule)
    Returns dict(output_dir, values (subjects x regions), stats (one row per region), report).
    """
    covariates = list(covariates or [])
    long = collect_roi_tables(bids_dir, atlas, level, desc, list_to_keep, list_to_remove)
    ses = session_values(long)
    kept, report = filter_regions(ses, min_voxels)
    dropped = report[~report['kept']]
    _log(f'ROI analysis {atlas} level {level} desc-{desc}: {report["kept"].sum()} of {len(report)} regions kept '
         f'(>= {min_voxels} PET voxels in every session)'
         + (': dropped ' + ', '.join(f'{r} ({v})' for r, v in zip(dropped.region, dropped.min_voxels))
            if len(dropped) else ''))
    norm, ref = normalise_regions(kept, normalization, all_regions=ses)
    norm_tag = normalization if isinstance(normalization, str) else 'ref'
    out = output_dir or opj(bids_dir, 'Results', f'Grp_PET_ROI_{design}_seg-{atlas}-lvl{level}_desc-{desc}'
                            f'_norm-{norm_tag}' + (f'_{output_suffix}' if output_suffix else ''))
    os.makedirs(out, exist_ok=True)

    # wide tables (one row per session), as the earlier Excel files
    wide_raw = kept.pivot_table(index=['subject', 'session'], columns='label', values='mean')
    wide = norm.pivot_table(index=['subject', 'session'], columns='label', values='value')
    names = kept.groupby('label')['region'].first()
    for w, f in ((wide_raw, 'values_raw'), (wide, 'values_normalised')):
        w2 = w.copy()
        w2.columns = [f'{names[c]} ({c})' for c in w2.columns]
        w2.to_csv(opj(out, f + '.tsv'), sep='\t')
        try:
            w2.to_excel(opj(out, f + '.xlsx'))
        except Exception as e:  # openpyxl missing
            _log(f'  WARNING: {f}.xlsx not written ({e})')
    report.to_csv(opj(out, 'regions_voxel_counts.tsv'), sep='\t', index=False)

    # observations
    obs = wide.reset_index()
    needed = covariates + ([group_col] if design == 'two-sample' else []) + \
        ([paired_col] if design == 'paired' and paired_col not in ('session', 'subject') else [])
    info = subject_table(obs[['subject', 'session']].assign(image='', run='', base='', desc=desc, data_path=''),
                         bids_dir, needed)
    obs = info.drop(columns=['image', 'run', 'base', 'data_path']).merge(obs, on=['subject', 'session'])
    if ref is not None:
        obs = obs.merge(ref, on=['subject', 'session'], how='left')
    labels = list(wide.columns)
    Y = obs[labels].to_numpy(float)
    if design == 'paired':
        if not paired_levels or len(paired_levels) != 2:
            raise ValueError("design='paired' needs paired_levels=[a, b]")
        a, b = (str(x) for x in paired_levels)
        ka = obs[paired_col].astype(str)
        ia = {s: i for i, s in enumerate(obs.subject) if ka.iloc[i] == a}
        ib = {s: i for i, s in enumerate(obs.subject) if ka.iloc[i] == b}
        subs = sorted(set(ia) & set(ib))
        Y = np.stack([Y[ia[s]] - Y[ib[s]] for s in subs])
        obs = obs.loc[[ia[s] for s in subs]].reset_index(drop=True)
    elif aggregate == 'subject' and obs.subject.duplicated().any():
        num = obs.select_dtypes('number').columns
        g = obs.groupby('subject', sort=True)
        first = g.first()
        first[num] = g[num].mean()
        obs = first.reset_index()
        Y = obs[labels].to_numpy(float)
        _log(f'  sessions of the same subject averaged: {len(obs)} subjects')
    design_glm = 'one-sample' if design == 'paired' else design
    X, c, cols, rows = _glm.design_matrix(obs, design_glm, group_col, groups, covariates, test_covariate)
    Y, obs = Y[rows], obs[rows].reset_index(drop=True)
    cname = _glm.contrast_name(design, c, cols, groups)
    glm = _glm.GLM(X, c)
    cope, t, resid = glm.fit(Y)
    p = _glm.parametric_p(t, glm.df)
    se = np.sqrt((resid ** 2).sum(0) / glm.df * glm.var_c)
    tcrit = float(_glm.stats.t.isf(0.025, glm.df))
    res = pd.DataFrame(dict(label=labels, region=[names[l] for l in labels], cope=cope, se=se,
                            ci95_low=cope - tcrit * se, ci95_high=cope + tcrit * se, t=t, df=glm.df, p=p,
                            p_FDR=_glm.bh_fdr(p)))
    res['significant_FDR'] = res['p_FDR'] < fdr_q
    if design == 'two-sample':
        for gname in groups:
            sel = (obs[group_col].astype(str) == str(gname)).to_numpy()
            res[f'mean_{gname}'] = Y[sel].mean(0)
            res[f'sd_{gname}'] = Y[sel].std(0, ddof=1)
        res['percent_difference'] = 100 * (res[f'mean_{groups[0]}'] - res[f'mean_{groups[1]}']) / res[f'mean_{groups[1]}']
    else:
        res['mean'] = Y.mean(0)
        res['sd'] = Y.std(0, ddof=1)
    res = res.sort_values('p').reset_index(drop=True)
    res.to_csv(opj(out, f'{cname}_regions.tsv'), sep='\t', index=False, float_format='%.6g')
    obs.drop(columns=labels).to_csv(opj(out, 'inputs.tsv'), sep='\t', index=False)
    _plot_effects(res, opj(out, f'{cname}_regions.png'), cname, fdr_q)
    _log(f'  {cname}: {int(res["significant_FDR"].sum())} region(s) with FDR q < {fdr_q} '
         f'({X.shape[0]} observations, df = {glm.df})')
    return dict(output_dir=out, values=wide, stats=res, report=report)


def _plot_effects(res, png, title, q, n_max=40):
    """Effect (with 95 % CI) of the regions with the smallest p."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        r = res.head(n_max).iloc[::-1]
        fig, ax = plt.subplots(figsize=(6, 0.22 * len(r) + 1.2))
        col = np.where(r['significant_FDR'], '#c0392b', '#7f8c8d')
        ax.errorbar(r['cope'], np.arange(len(r)), xerr=[r['cope'] - r['ci95_low'], r['ci95_high'] - r['cope']],
                    fmt='none', ecolor='#555555', lw=0.8)
        ax.scatter(r['cope'], np.arange(len(r)), c=col, s=14, zorder=3)
        ax.axvline(0, color='k', lw=0.6)
        ax.set_yticks(np.arange(len(r)))
        ax.set_yticklabels(r['region'], fontsize=6)
        ax.set_xlabel('effect (normalised units), 95 % CI')
        ax.set_title(f'{title} (red: FDR q < {q})', fontsize=8)
        fig.tight_layout()
        fig.savefig(png, dpi=200)
        fig.savefig(png.replace('.png', '.svg'))
        plt.close(fig)
    except Exception as e:
        _log(f'  [plot] WARNING: {e}')
