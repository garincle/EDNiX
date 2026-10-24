"""
Cross-species regional profiles: fingerprints and bar plots (e.g. FDG-PET, human vs macaque).

The PET pipeline writes one table per run, atlas and level (modalities/pet/roi_signal.py):
    <bids>/derivatives/sub-*/ses-*/pet/orig/postprocessing/ROI_analysis/
        <base>_seg-<atlas>-lvl<k>_desc-<desc>_ROIs.tsv     label, region, n_voxels, n_valid, mean, sd, median
This module gathers those tables for several species and compares them region by region.

HOMOLOGY
--------
The EDNiX atlases (EDNIxCSCLR / EDNIxCSC) share ONE label table across species (atlases/label_translate.py
copies the same LUT into every species' label_code folder): a label id, or its name, designates the same
region in every species. Regions are therefore matched by name (hemisphere prefix L_/R_ removed when the
hemispheres are merged), and only regions present in every species enter the comparison. A non-EDNiX atlas
(different tables per species) must not be used here.

WHAT IS COMPARED
----------------
Absolute uptake (SUV) is not comparable across species: injected dose per kg, body-mass scaling of the
plasma input, scanner, reconstruction and anaesthesia all differ. What is compared is the within-subject
RELATIVE profile: by default the ratio of each region to the whole-brain mean ('bratio' image of the
pipeline, i.e. proportional scaling; Fox et al. 1988; Dukart et al. 2010 for the choice of reference),
or, optionally, a ratio to reference regions or a within-subject z-score across regions.

UNIT OF ANALYSIS
----------------
One subject = one sample (runs and sessions of a subject are averaged first), as in the Figure 5 pipeline.

PARTIAL VOLUME
--------------
A macaque brain is ~15x smaller than a human brain, scanned at a comparable resolution: small macaque
regions lose more signal to partial volume, which lowers cortical ratios systematically, independently of
metabolism (Hoffman et al. 1979; Soret, Bacharach & Buvat 2007). Two guards:
  - a region is kept only if it has >= MIN_VOXELS valid voxels in EVERY subject (old script
    3_correct_atlas_for_small_V.py rule, 10 voxels);
  - with FWHM_MM given per species, regions whose equivalent cube side is below MIN_SIZE_FWHM x FWHM are
    dropped, and size_in_fwhm is written for every region and species either way.
Use the partial-volume corrected images ('<desc>-PVC') when they exist.

STATISTICS
----------
Per region, every species vs REF_SPECIES: difference of means, Hedges' g with its 95% CI (Hedges & Olkin
1985), p from a two-sided permutation test of the species labels (exact when the number of relabelings is
small; Welch's t-test or Mann-Whitney on request), Benjamini-Hochberg FDR across regions within each
comparison. Profile similarity between species: Pearson and Spearman r between the species mean profiles,
95% CI by bootstrap of the subjects within species (regions are not exchangeable -- spatial
autocorrelation -- so permuting regions would give an anti-conservative p, and none is reported), and the
split-half reliability of each species' profile (Spearman-Brown) with the disattenuated r
r / sqrt(rel_a * rel_b) (Spearman 1904): the between-species similarity a perfectly reliable profile would
reach.

OUTPUT (OUT_DIR)
----------------
    profiles_<tag>.xlsx      data (long), region filter, subject profiles, species stats, similarity
    F1_fingerprint_<tag>     polar profile per species (mean +- SEM) + polar effect sizes (g) vs reference
    F2_bars_<tag>            grouped bars (mean +- SEM) with every subject, FDR stars
    F3_scatter_<tag>         region means of each species vs reference, identity line
    methods_<tag>.txt        parameters and method text
Every figure: .png, .svg and a .txt legend.

References: Passingham, Stephan & Kotter 2002 Nat Rev Neurosci (regional fingerprints); Mars et al. 2018
eLife (comparative fingerprints human-macaque); Weissgerber et al. 2015 PLoS Biol (show every data point).
"""
import os
import re
import glob
import textwrap

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats as sstats

opj = os.path.join
opb = os.path.basename

# Okabe-Ito palette, assigned to the sorted species names exactly as fig5_analysis.SPECIES_COLORS, so that a
# species keeps the same colour in every EDNiX figure.
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#000000"]
_FIG5_SPECIES = ["Human", "Macaque", "Marmoset", "Mouselemur", "Rat", "Mouse", "Dog"]
SPECIES_COLORS = {n: PALETTE[i % len(PALETTE)] for i, n in enumerate(sorted(set(_FIG5_SPECIES)))}

ROI_SUBDIR = opj('pet', 'orig', 'postprocessing', 'ROI_analysis')
_HEMI = re.compile(r'^([LlRr])_')


def species_color(sp, i=0):
    return SPECIES_COLORS.get(sp, PALETTE[(i + len(SPECIES_COLORS)) % len(PALETTE)])


# ======================================================================================================== data
def find_roi_tables(bids_dir, atlas, level, desc):
    """ROI tables of every run under <bids_dir>/derivatives (or bids_dir itself when it is the derivatives
    folder) for one atlas, level and desc (exact desc: 'bratio' does not match 'bratio-mean')."""
    root = opj(bids_dir, 'derivatives') if os.path.isdir(opj(bids_dir, 'derivatives')) else bids_dir
    pat = opj(root, '**', ROI_SUBDIR, f'*_seg-{atlas}-lvl{level}_desc-{desc}_ROIs.tsv')
    return sorted(glob.glob(pat, recursive=True))


def _entity(name, key):
    m = re.search(rf'(?:^|_){key}-([A-Za-z0-9]+)', name)
    return m.group(1) if m else None


def read_roi_table(path, value_col='mean'):
    tab = pd.read_csv(path, sep='\t')
    missing = {'label', 'region', value_col} - set(tab.columns)
    if missing:
        raise ValueError(f'{path}: columns {sorted(missing)} missing')
    n = tab['n_valid'] if 'n_valid' in tab.columns else tab['n_voxels']
    name = opb(path)
    hemi = tab['region'].astype(str).str.extract(_HEMI, expand=False).str.upper()
    return pd.DataFrame({
        'subject': _entity(name, 'sub'), 'session': _entity(name, 'ses'), 'run': _entity(name, 'run'),
        'file': path, 'label': tab['label'].astype(int), 'region': tab['region'].astype(str),
        'base_region': tab['region'].astype(str).str.replace(_HEMI, '', regex=True),
        'hemisphere': hemi.map({'L': 'left', 'R': 'right'}).fillna('bilateral'),
        'n_valid': n.astype(float),
        'volume_mm3': tab['volume_mm3'].astype(float) if 'volume_mm3' in tab.columns else np.nan,
        'value': tab[value_col].astype(float)})


def collect(sources, atlas='EDNIxCSCLR', level=1, desc=('bratio', 'bratio-mean'), value_col='mean'):
    """
    Long table of every run of every species.
    sources: {species: bids_dir} or {species: {'bids_dir': ..., 'desc': ..., 'files': [...]}}; 'files' (explicit
    list of tables) overrides the search. desc: one desc or candidates tried in order (the first one with
    tables is used for that species: static human 'bratio', dynamic macaque 'bratio-mean').
    """
    rows = []
    for sp, src in sources.items():
        src = {'bids_dir': src} if isinstance(src, str) else dict(src)
        descs = src.get('desc', desc)
        descs = [descs] if isinstance(descs, str) else list(descs)
        files, used = list(src.get('files') or []), src.get('desc') if src.get('files') else None
        if not files:
            for d in descs:
                files = find_roi_tables(src['bids_dir'], atlas, level, d)
                if files:
                    used = d
                    break
        if not files:
            raise FileNotFoundError(f'{sp}: no {atlas} lvl{level} ROI table with desc in {descs} under '
                                    f'{src.get("bids_dir")} (run step 7 of the PET pipeline with atlases)')
        print(f'  {sp}: {len(files)} table(s), desc-{used}')
        for f in files:
            t = read_roi_table(f, value_col)
            t.insert(0, 'species', sp)
            t['desc'] = used
            if t['subject'].isna().any():
                raise ValueError(f'{f}: no sub-<id> entity in the file name')
            rows.append(t)
    return pd.concat(rows, ignore_index=True)


def subject_level(df, hemispheres='merge'):
    """Mean over the runs / sessions of a subject, then (hemispheres='merge') volume-weighted mean of the left
    and right parts of each region (weights: valid voxels). Region key: base name (merge) or full name."""
    if hemispheres not in ('merge', 'separate'):
        raise ValueError(f'hemispheres {hemispheres!r}: expected "merge" or "separate"')
    key = 'base_region' if hemispheres == 'merge' else 'region'
    d = df.copy()
    ok = d['value'].notna() & (d['n_valid'] > 0)
    d['_wv'] = np.where(ok, d['value'] * d['n_valid'], 0.0)
    d['_w'] = np.where(ok, d['n_valid'], 0.0)
    run = (d.groupby(['species', 'subject', 'file', key], sort=False)
             .agg(_wv=('_wv', 'sum'), _w=('_w', 'sum'), n_valid=('n_valid', 'sum'),
                  min_label_voxels=('n_valid', 'min'),
                  volume_mm3=('volume_mm3', lambda v: v.sum(min_count=1)), label=('label', 'min'))
             .reset_index())
    run['value'] = np.where(run['_w'] > 0, run['_wv'] / run['_w'].where(run['_w'] > 0), np.nan)
    sub = (run.groupby(['species', 'subject', key], sort=False)
              .agg(value=('value', 'mean'), n_valid=('n_valid', 'min'), min_label_voxels=('min_label_voxels', 'min'),
                   volume_mm3=('volume_mm3', 'mean'),
                   label=('label', 'min'), n_runs=('file', 'nunique'))
              .reset_index().rename(columns={key: 'region_key'}))
    return sub


def normalise(sub, method=None, ref_regions=None):
    """None: values as they are (already a ratio, e.g. bratio); 'reference': ratio to the volume-weighted mean
    of ref_regions (list of region keys, or a regular expression); 'zscore': within-subject z-score across
    regions (keeps only the shape of the profile)."""
    if method in (None, 'none'):
        return sub
    out = []
    for _, g in sub.groupby(['species', 'subject'], sort=False):
        g = g.copy()
        if method == 'reference':
            if isinstance(ref_regions, str):
                m = g['region_key'].str.contains(ref_regions, flags=re.IGNORECASE, regex=True)
            else:
                m = g['region_key'].isin(ref_regions or [])
            if not m.any():
                raise ValueError(f'reference regions {ref_regions!r} not found for {g.species.iat[0]} '
                                 f'sub-{g.subject.iat[0]}')
            g['value'] = g['value'] / np.average(g.loc[m, 'value'], weights=g.loc[m, 'n_valid'])
        elif method == 'zscore':
            g['value'] = (g['value'] - g['value'].mean()) / g['value'].std(ddof=1)
        else:
            raise ValueError(f'normalisation {method!r}: expected None, "reference" or "zscore"')
        out.append(g)
    return pd.concat(out, ignore_index=True)


def filter_regions(sub, min_voxels=10, fwhm_mm=None, min_size_fwhm=None, exclude_regex=None):
    """
    Regions kept for the comparison: present in every species, >= min_voxels valid voxels in every subject (in
    each hemisphere label when the hemispheres are merged, as the old per-label rule),
    not matched by exclude_regex, and (fwhm_mm {species: FWHM in mm} and min_size_fwhm given) equivalent cube
    side (median volume ** 1/3) >= min_size_fwhm x FWHM in every species. Returns (kept keys in atlas order,
    report table).
    """
    species = list(dict.fromkeys(sub['species']))
    rep = (sub.groupby(['region_key', 'species'])
              .agg(label=('label', 'min'), min_voxels=('min_label_voxels', 'min'), median_volume_mm3=('volume_mm3', 'median'),
                   n_subjects=('subject', 'nunique'), n_missing=('value', lambda v: int(v.isna().sum())))
              .reset_index())
    fwhm_mm = fwhm_mm or {}
    rep['fwhm_mm'] = rep['species'].map(fwhm_mm)
    rep['size_in_fwhm'] = np.cbrt(rep['median_volume_mm3']) / rep['fwhm_mm']
    n_sub = sub.groupby('species')['subject'].nunique()
    rep['ok_voxels'] = (rep['min_voxels'] >= min_voxels) & (rep['n_missing'] == 0) \
        & (rep['n_subjects'] == rep['species'].map(n_sub))
    rep['ok_size'] = True
    if min_size_fwhm:
        rep['ok_size'] = rep['size_in_fwhm'].isna() | (rep['size_in_fwhm'] >= min_size_fwhm)
    wide = rep.pivot(index='region_key', columns='species', values='ok_voxels')
    size = rep.pivot(index='region_key', columns='species', values='ok_size')
    present = wide.reindex(columns=species).notna().all(axis=1)
    keep = present & wide.reindex(columns=species).fillna(False).astype(bool).all(axis=1) \
        & size.reindex(columns=species).fillna(False).astype(bool).all(axis=1)
    if exclude_regex:
        keep &= ~keep.index.to_series().str.contains(exclude_regex, flags=re.IGNORECASE, regex=True)
    order = rep.groupby('region_key')['label'].min().sort_values()
    kept = [k for k in order.index if keep.get(k, False)]
    rep['kept'] = rep['region_key'].isin(kept)
    rep['present_in_all_species'] = rep['region_key'].map(present)
    return kept, rep.sort_values(['label', 'species']).reset_index(drop=True)


def profile_matrix(sub, regions):
    """Wide table: index (species, subject), columns regions (atlas order)."""
    return (sub[sub['region_key'].isin(regions)]
            .pivot_table(index=['species', 'subject'], columns='region_key', values='value')
            .reindex(columns=regions))


# ================================================================================================== statistics
def hedges_g(a, b):
    """Hedges' g of b - a (pooled SD, small-sample correction) and its approximate 95% CI."""
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return np.nan, np.nan, np.nan
    sp = np.sqrt(((n1 - 1) * np.var(a, ddof=1) + (n2 - 1) * np.var(b, ddof=1)) / (n1 + n2 - 2))
    if sp == 0:
        return np.nan, np.nan, np.nan
    j = 1 - 3 / (4 * (n1 + n2) - 9)
    g = j * (np.mean(b) - np.mean(a)) / sp
    se = np.sqrt((n1 + n2) / (n1 * n2) + g ** 2 / (2 * (n1 + n2)))
    return g, g - 1.96 * se, g + 1.96 * se


def _pvalue(a, b, test, n_perm, rng):
    if test == 'welch':
        return sstats.ttest_ind(b, a, equal_var=False).pvalue
    if test == 'mannwhitney':
        return sstats.mannwhitneyu(b, a, alternative='two-sided').pvalue
    if test == 'permutation':
        res = sstats.permutation_test((b, a), lambda x, y: np.mean(x) - np.mean(y), permutation_type='independent',
                                      n_resamples=n_perm, alternative='two-sided', random_state=rng)
        return res.pvalue
    raise ValueError(f'test {test!r}: expected "permutation", "welch" or "mannwhitney"')


def bh_fdr(p):
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        pv = p[ok]
        o = np.argsort(pv)
        r = pv[o] * len(pv) / np.arange(1, len(pv) + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]
        q_ok = np.empty_like(pv)
        q_ok[o] = np.minimum(r, 1)
        q[ok] = q_ok
    return q


def compare_species(prof, ref_species, test='permutation', n_perm=10000, seed=0):
    """Region-wise comparison of every species with ref_species (difference = species - reference)."""
    rng = np.random.default_rng(seed)
    species = [s for s in prof.index.get_level_values(0).unique() if s != ref_species]
    if ref_species not in prof.index.get_level_values(0):
        raise ValueError(f'reference species {ref_species} not in the data')
    a_all = prof.loc[ref_species]
    out = []
    for sp in species:
        b_all = prof.loc[sp]
        rows = []
        for reg in prof.columns:
            a, b = a_all[reg].dropna().values, b_all[reg].dropna().values
            g, lo, hi = hedges_g(a, b)
            p = _pvalue(a, b, test, n_perm, rng) if len(a) >= 2 and len(b) >= 2 else np.nan
            rows.append(dict(region=reg, species=sp, reference=ref_species, n_ref=len(a), n_species=len(b),
                             mean_ref=a.mean(), sem_ref=sstats.sem(a) if len(a) > 1 else np.nan,
                             mean_species=b.mean(), sem_species=sstats.sem(b) if len(b) > 1 else np.nan,
                             difference=b.mean() - a.mean(), hedges_g=g, g_ci_low=lo, g_ci_high=hi, p=p))
        r = pd.DataFrame(rows)
        r['q_fdr'] = bh_fdr(r['p'])
        r['test'] = test
        out.append(r)
    return pd.concat(out, ignore_index=True)


def split_half_reliability(m, n_split=1000, rng=None):
    """Spearman-Brown corrected split-half reliability of the mean profile (subjects x regions); NaN for < 4
    subjects."""
    rng = rng or np.random.default_rng(0)
    n = m.shape[0]
    if n < 4:
        return np.nan
    rs = []
    for _ in range(n_split):
        idx = rng.permutation(n)
        h1, h2 = m[idx[:n // 2]].mean(0), m[idx[n // 2:]].mean(0)
        rs.append(np.corrcoef(h1, h2)[0, 1])
    r = np.nanmean(rs)
    return 2 * r / (1 + r)


def profile_similarity(prof, ref_species, n_boot=5000, seed=0):
    """Pearson / Spearman r of every species' mean profile with the reference's, bootstrap CI over subjects,
    split-half reliabilities and disattenuated r."""
    rng = np.random.default_rng(seed)
    a = prof.loc[ref_species].values
    rel_a = split_half_reliability(a, rng=rng)
    rows = []
    for sp in prof.index.get_level_values(0).unique():
        if sp == ref_species:
            continue
        b = prof.loc[sp].values
        r_p = np.corrcoef(a.mean(0), b.mean(0))[0, 1]
        r_s = sstats.spearmanr(a.mean(0), b.mean(0)).statistic
        boot_p, boot_s = [], []
        for _ in range(n_boot):
            ma = a[rng.integers(0, len(a), len(a))].mean(0)
            mb = b[rng.integers(0, len(b), len(b))].mean(0)
            boot_p.append(np.corrcoef(ma, mb)[0, 1])
            boot_s.append(sstats.spearmanr(ma, mb).statistic)
        rel_b = split_half_reliability(b, rng=rng)
        ceil = np.sqrt(rel_a * rel_b) if np.isfinite(rel_a) and np.isfinite(rel_b) and rel_a > 0 and rel_b > 0 \
            else np.nan
        rows.append(dict(species=sp, reference=ref_species, n_regions=a.shape[1], n_ref=len(a), n_species=len(b),
                         pearson_r=r_p, pearson_ci_low=np.nanpercentile(boot_p, 2.5),
                         pearson_ci_high=np.nanpercentile(boot_p, 97.5),
                         spearman_r=r_s, spearman_ci_low=np.nanpercentile(boot_s, 2.5),
                         spearman_ci_high=np.nanpercentile(boot_s, 97.5),
                         reliability_ref=rel_a, reliability_species=rel_b, noise_ceiling=ceil,
                         pearson_r_disattenuated=r_p / ceil if np.isfinite(ceil) else np.nan))
    return pd.DataFrame(rows)


# ======================================================================================================= plots
def _stars(q, alpha):
    if not np.isfinite(q) or q >= alpha:
        return ''
    return '***' if q < alpha / 50 else '**' if q < alpha / 5 else '*'


def _short(name, n=22):
    name = name.replace('_', ' ')
    return name if len(name) <= n else name[:n - 1] + '.'


def _save(fig, pre, legend):
    for ext in ('png', 'svg'):
        fig.savefig(f'{pre}.{ext}', dpi=300 if ext == 'png' else None, bbox_inches='tight')
    plt.close(fig)
    with open(pre + '.txt', 'w') as f:
        f.write(textwrap.fill(legend, 100) + '\n')
    print(f'  [FIG] {pre}.png')


def _baseline(ylabel):
    return 0.0 if 'z-score' in ylabel else 1.0


def plot_fingerprint(prof, stats, ref_species, pre, ylabel, alpha=0.05, legend_extra=''):
    """(a) polar profile of every species, mean +- SEM, regions in atlas order, dashed circle at the baseline
    (1 = whole-brain / reference mean, 0 for z-scores), stars at q < alpha vs reference;
    (b) polar bars of Hedges' g (species - reference), filled when q < alpha."""
    regions = list(prof.columns)
    n = len(regions)
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    angc = np.r_[ang, ang[:1]]
    species = list(prof.index.get_level_values(0).unique())
    others = [s for s in species if s != ref_species]
    fig = plt.figure(figsize=(7.2 + 3.6 * len(others), 7.2))
    ax = fig.add_subplot(1, 1 + len(others), 1, projection='polar')
    mx = []
    for i, sp in enumerate(species):
        m, s = prof.loc[sp].mean().values, prof.loc[sp].sem().values
        s = np.nan_to_num(s)
        c = species_color(sp, i)
        ax.plot(angc, np.r_[m, m[:1]], color=c, lw=1.8, label=f'{sp} (n={prof.loc[sp].shape[0]})')
        ax.fill_between(angc, np.r_[m - s, m[:1] - s[:1]], np.r_[m + s, m[:1] + s[:1]], color=c, alpha=0.25, lw=0)
        mx.append(np.nanmax(m + s))
        mx.append(np.nanmin(m - s))
    lo, hi = min(mx), max(mx)
    pad = 0.12 * (hi - lo if hi > lo else 1)
    ax.set_ylim(lo - pad, hi + 2 * pad)
    ax.plot(np.linspace(0, 2 * np.pi, 200), np.full(200, _baseline(ylabel)), 'k--', lw=0.7)
    st = stats.set_index(['species', 'region'])
    for sp in others:
        for a_, reg in zip(ang, regions):
            t = _stars(st.loc[(sp, reg), 'q_fdr'], alpha)
            if t:
                ax.text(a_, hi + 1.4 * pad, t, ha='center', va='center', fontsize=7,
                        color=species_color(sp, species.index(sp)))
    ax.set_xticks(ang)
    ax.set_xticklabels([_short(r) for r in regions], fontsize=6 if n > 20 else 8)
    ax.tick_params(axis='y', labelsize=6)
    ax.set_title(f'a  regional profile ({ylabel})', loc='left', fontsize=9)
    ax.legend(loc='upper right', bbox_to_anchor=(1.25, 1.12), fontsize=7, frameon=False)
    for k, sp in enumerate(others):
        bx = fig.add_subplot(1, 1 + len(others), 2 + k, projection='polar')
        s = st.loc[sp].reindex(regions)
        g = s['hedges_g'].values
        sig = s['q_fdr'].values < alpha
        width = 2 * np.pi / n * 0.85
        base = np.nanmax(np.abs(g)) if np.isfinite(g).any() else 1
        c = species_color(sp, species.index(sp))
        cref = species_color(ref_species, species.index(ref_species))
        cols = [c if v > 0 else cref for v in np.nan_to_num(g)]
        bx.bar(ang, np.abs(np.nan_to_num(g)), width=width, color=cols, alpha=0.9,
               edgecolor=cols, linewidth=0.6, fill=False)
        bx.bar(ang[sig], np.abs(g[sig]), width=width, color=np.array(cols, dtype=object)[sig], alpha=0.9)
        bx.set_ylim(0, base * 1.1)
        bx.set_xticks(ang)
        bx.set_xticklabels([_short(r) for r in regions], fontsize=6 if n > 20 else 8)
        bx.tick_params(axis='y', labelsize=6)
        bx.set_title(f'{"bcdefg"[k]}  |Hedges g|, {sp} vs {ref_species}\n'
                     f'colour = higher species; filled = q < {alpha}', loc='left', fontsize=9)
    fig.tight_layout()
    legend = (f'Fingerprint of the regional {ylabel} in {", ".join(species)}. (a) Mean across subjects '
              f'(line) +- SEM (band) of each region, regions in atlas order; dashed circle: baseline '
              f'({_baseline(ylabel):g}); stars: Benjamini-Hochberg q < {alpha} (*), < {alpha / 5:g} (**), '
              f'< {alpha / 50:g} (***), test {stats["test"].iat[0]}, each species vs {ref_species}. '
              f'(b-) Absolute Hedges g per region (species - {ref_species}); bar colour = the species with the '
              f'higher value; filled bars: q < {alpha}. One subject = one sample (runs averaged). '
              + legend_extra)
    _save(fig, pre, legend)


def plot_bars(prof, stats, ref_species, pre, ylabel, alpha=0.05, per_row=24, legend_extra=''):
    """Grouped bars (mean +- SEM) with every subject, regions in atlas order, wrapped every per_row regions."""
    regions = list(prof.columns)
    species = list(prof.index.get_level_values(0).unique())
    rows = [regions[i:i + per_row] for i in range(0, len(regions), per_row)]
    st = stats.set_index(['species', 'region'])
    w = 0.8 / len(species)
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(len(rows), 1, figsize=(max(6.0, 0.42 * min(per_row, len(regions)) * len(species)),
                                                     3.2 * len(rows)), squeeze=False)
    vals = prof.values[np.isfinite(prof.values)]
    ylo, yhi = np.min(vals), np.max(vals)
    base = _baseline(ylabel)
    ylo, yhi = min(ylo, base), max(yhi, base)
    span = yhi - ylo if yhi > ylo else 1
    for ax, regs in zip(axes[:, 0], rows):
        x = np.arange(len(regs))
        for i, sp in enumerate(species):
            d = prof.loc[sp, regs]
            c = species_color(sp, i)
            xs = x - 0.4 + w * (i + 0.5)
            # bars start at the baseline (1 = whole-brain mean, 0 for z-scores): their length is the
            # deviation from it, so the axis can be zoomed without truncating a bar drawn from 0
            ax.bar(xs, d.mean().values - base, w * 0.9, bottom=base, yerr=d.sem().values, color=c, alpha=0.55,
                   edgecolor=c,
                   error_kw=dict(lw=0.8, capsize=1.5), label=f'{sp} (n={d.shape[0]})' if regs is rows[0] else None)
            for j, reg in enumerate(regs):
                v = d[reg].dropna().values
                ax.scatter(xs[j] + rng.uniform(-w * 0.25, w * 0.25, len(v)), v, s=6, color=c,
                           edgecolors='white', linewidths=0.3, zorder=3)
        for j, reg in enumerate(regs):
            t = ' '.join(_stars(st.loc[(sp, reg), 'q_fdr'], alpha) for sp in species if sp != ref_species).strip()
            if t:
                ax.text(x[j], yhi + 0.04 * span, t, ha='center', va='bottom', fontsize=7)
        ax.axhline(_baseline(ylabel), color='k', lw=0.6, ls='--')
        ax.set_xticks(x)
        ax.set_xticklabels([_short(r, 28) for r in regs], rotation=60, ha='right', fontsize=7)
        ax.set_xlim(-0.6, len(regs) - 0.4)
        ax.set_ylim(ylo - 0.05 * span, yhi + 0.15 * span)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(axis='y', labelsize=7)
        ax.spines[['top', 'right']].set_visible(False)
    axes[0, 0].legend(fontsize=7, frameon=False, ncol=len(species), loc='upper left', bbox_to_anchor=(0, 1.25))
    fig.tight_layout()
    legend = (f'Regional {ylabel} in {", ".join(species)}: bars, mean across subjects +- SEM, drawn from the '
              f'baseline ({_baseline(ylabel):g}) so that their length is the deviation from it; dots, one subject '
              f'(runs averaged); stars, Benjamini-Hochberg q < '
              f'{alpha} (*), < {alpha / 5:g} (**), < {alpha / 50:g} (***) of the {stats["test"].iat[0]} test of each '
              f'species vs {ref_species}, corrected across the {len(regions)} regions. ' + legend_extra)
    _save(fig, pre, legend)


def plot_scatter(prof, sim, ref_species, pre, ylabel, legend_extra=''):
    """Region means of every species vs the reference (+- SEM both axes), identity line, r [95% CI]."""
    species = [s for s in prof.index.get_level_values(0).unique() if s != ref_species]
    fig, axes = plt.subplots(1, len(species), figsize=(3.6 * len(species), 3.6), squeeze=False)
    a = prof.loc[ref_species]
    sim = sim.set_index('species')
    for ax, sp in zip(axes[0], species):
        b = prof.loc[sp]
        c = species_color(sp, 1)
        ax.errorbar(a.mean(), b.mean(), xerr=a.sem(), yerr=b.sem(), fmt='o', ms=3, color=c, ecolor=c,
                    elinewidth=0.5, alpha=0.85)
        lim = [np.nanmin([a.mean().min(), b.mean().min()]), np.nanmax([a.mean().max(), b.mean().max()])]
        pad = 0.05 * (lim[1] - lim[0] if lim[1] > lim[0] else 1)
        lim = [lim[0] - pad, lim[1] + pad]
        ax.plot(lim, lim, 'k--', lw=0.6)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_aspect('equal')
        s = sim.loc[sp]
        txt = (f'r = {s.pearson_r:.2f} [{s.pearson_ci_low:.2f}, {s.pearson_ci_high:.2f}]\n'
               f'rho = {s.spearman_r:.2f} [{s.spearman_ci_low:.2f}, {s.spearman_ci_high:.2f}]')
        if np.isfinite(s.noise_ceiling):
            txt += f'\nceiling = {s.noise_ceiling:.2f}'
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va='top', fontsize=7)
        ax.set_xlabel(f'{ref_species} {ylabel}', fontsize=8)
        ax.set_ylabel(f'{sp} {ylabel}', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    legend = (f'Region means ({ylabel}) of each species vs {ref_species}; error bars, SEM across subjects; dashed, '
              f'identity. r (Pearson) and rho (Spearman) between species mean profiles, 95% CI by bootstrap of the '
              f'subjects within species (no p-value: regions are spatially autocorrelated, not exchangeable). '
              f'ceiling: sqrt of the product of the Spearman-Brown split-half reliabilities of the two profiles, '
              f'the r two perfectly matching but noisy profiles would reach. ' + legend_extra)
    _save(fig, pre, legend)


# ========================================================================================================= run
DEFAULTS = dict(
    ATLAS='EDNIxCSCLR', LEVELS=[1], DESC=('bratio', 'bratio-mean'), VALUE_COL='mean',
    HEMISPHERES='merge',            # 'merge' (volume-weighted L+R) or 'separate'
    NORMALISATION=None,             # None (bratio is already a ratio), 'reference', 'zscore'
    REF_REGIONS=None,               # region keys or regex for NORMALISATION='reference'
    REF_SPECIES='Human',
    MIN_VOXELS=10, FWHM_MM=None, MIN_SIZE_FWHM=None, EXCLUDE_REGEX=None,
    TEST='permutation', N_PERM=10000, N_BOOT=5000, ALPHA=0.05, SEED=0, BARS_PER_ROW=24,
)


def _ylabel(cfg, descs):
    if cfg['NORMALISATION'] == 'zscore':
        return 'within-subject z-score'
    if cfg['NORMALISATION'] == 'reference':
        return 'ratio to reference regions'
    if all(d and d.startswith('bratio') for d in descs):
        return 'ratio to whole-brain mean'
    return '/'.join(sorted(set(d for d in descs if d)))


def _methods(cfg, level, kept, rep, sub, sim, descs):
    n = sub.groupby('species')['subject'].nunique().to_dict()
    lines = [
        f'Cross-species regional profiles -- atlas {cfg["ATLAS"]}, level {level}.',
        f'Images: desc {descs}; regional statistic: {cfg["VALUE_COL"]} over the valid voxels of each label '
        f'(Tools/roiStats.py), in the native PET space.',
        f'Subjects per species: {n}. Runs and sessions of a subject averaged (one subject = one sample).',
        f'Hemispheres: {cfg["HEMISPHERES"]} (merge = valid-voxel-weighted mean of L_ and R_ labels).',
        f'Normalisation: {cfg["NORMALISATION"] or "none (values used as written)"}'
        + (f', reference regions {cfg["REF_REGIONS"]}' if cfg['NORMALISATION'] == 'reference' else '') + '.',
        f'Regions: matched by name across species (shared EDNiX label table); kept when present in every species, '
        f'>= {cfg["MIN_VOXELS"]} valid voxels in every subject'
        + (f', equivalent cube side >= {cfg["MIN_SIZE_FWHM"]} x FWHM {cfg["FWHM_MM"]}' if cfg['MIN_SIZE_FWHM'] else '')
        + (f', not matching {cfg["EXCLUDE_REGEX"]!r}' if cfg['EXCLUDE_REGEX'] else '')
        + f': {len(kept)} of {rep["region_key"].nunique()} regions.',
        f'Statistics: each species vs {cfg["REF_SPECIES"]}, region by region, {cfg["TEST"]} test'
        + (f' ({cfg["N_PERM"]} relabelings, exact when fewer exist)' if cfg['TEST'] == 'permutation' else '')
        + f', two-sided; Hedges g with approximate 95% CI; Benjamini-Hochberg FDR across regions, alpha '
        f'{cfg["ALPHA"]}.',
        f'Profile similarity: Pearson and Spearman r of the species mean profiles, 95% CI from {cfg["N_BOOT"]} '
        f'bootstrap resamples of subjects within species; split-half reliability (1000 random splits, '
        f'Spearman-Brown) and disattenuated r = r / sqrt(rel_ref x rel_species).',
    ]
    for _, s in sim.iterrows():
        lines.append(f'  {s.species} vs {s.reference}: r = {s.pearson_r:.3f} [{s.pearson_ci_low:.3f}, '
                     f'{s.pearson_ci_high:.3f}], rho = {s.spearman_r:.3f}, noise ceiling = {s.noise_ceiling:.3f}')
    small = rep[(rep['size_in_fwhm'] < 2) & rep['kept']]
    if len(small):
        lines.append(f'WARNING: {small["region_key"].nunique()} kept region(s) are smaller than 2 FWHM in at least '
                     f'one species (size_in_fwhm in the region sheet): their values are dominated by partial volume.')
    return '\n'.join(lines) + '\n'


def run(sources, out_dir, config=None):
    """
    sources: {species: bids_dir or {'bids_dir', 'desc', 'files'}}; out_dir: output folder; config: overrides of
    DEFAULTS. Returns {level: dict of tables}.
    """
    cfg = dict(DEFAULTS, **(config or {}))
    unknown = set(config or {}) - set(DEFAULTS)
    if unknown:
        raise ValueError(f'unknown config key(s) {sorted(unknown)}; expected {sorted(DEFAULTS)}')
    if cfg['REF_SPECIES'] not in sources:
        raise ValueError(f'REF_SPECIES {cfg["REF_SPECIES"]} not in sources {list(sources)}')
    os.makedirs(out_dir, exist_ok=True)
    results = {}
    for level in cfg['LEVELS']:
        tag = f'{cfg["ATLAS"]}-lvl{level}' + ('_LR' if cfg['HEMISPHERES'] != 'merge' else '') \
            + (f'_{cfg["NORMALISATION"]}' if cfg['NORMALISATION'] else '')
        print(f'\n=== {tag} ===')
        df = collect(sources, cfg['ATLAS'], level, cfg['DESC'], cfg['VALUE_COL'])
        descs = sorted(df.groupby('species')['desc'].first().to_dict().items())
        sub = subject_level(df, cfg['HEMISPHERES'])
        sub = normalise(sub, cfg['NORMALISATION'], cfg['REF_REGIONS'])
        kept, rep = filter_regions(sub, cfg['MIN_VOXELS'], cfg['FWHM_MM'], cfg['MIN_SIZE_FWHM'], cfg['EXCLUDE_REGEX'])
        print(f'  {len(kept)} region(s) kept of {rep["region_key"].nunique()}')
        if len(kept) < 3:
            raise ValueError(f'{tag}: only {len(kept)} region(s) common to every species pass the filters '
                             f'(see the region table); nothing to compare')
        prof = profile_matrix(sub, kept)
        order = [cfg['REF_SPECIES']] + [s for s in sources if s != cfg['REF_SPECIES']]
        prof = prof.reindex(order, level=0)
        stats = compare_species(prof, cfg['REF_SPECIES'], cfg['TEST'], cfg['N_PERM'], cfg['SEED'])
        sim = profile_similarity(prof, cfg['REF_SPECIES'], cfg['N_BOOT'], cfg['SEED'])
        ylabel = _ylabel(cfg, [d for _, d in descs])
        extra = (f'Atlas {cfg["ATLAS"]} level {level}, {len(kept)} regions common to all species '
                 f'(>= {cfg["MIN_VOXELS"]} voxels in every subject); images: '
                 + ', '.join(f'{s} desc-{d}' for s, d in descs) + '.')
        with pd.ExcelWriter(opj(out_dir, f'profiles_{tag}.xlsx')) as xw:
            df.to_excel(xw, sheet_name='runs_long', index=False)
            rep.to_excel(xw, sheet_name='regions', index=False)
            prof.reset_index().to_excel(xw, sheet_name='subject_profiles', index=False)
            stats.to_excel(xw, sheet_name='species_stats', index=False)
            sim.to_excel(xw, sheet_name='similarity', index=False)
        plot_fingerprint(prof, stats, cfg['REF_SPECIES'], opj(out_dir, f'F1_fingerprint_{tag}'), ylabel,
                         cfg['ALPHA'], extra)
        plot_bars(prof, stats, cfg['REF_SPECIES'], opj(out_dir, f'F2_bars_{tag}'), ylabel, cfg['ALPHA'],
                  cfg['BARS_PER_ROW'], extra)
        plot_scatter(prof, sim, cfg['REF_SPECIES'], opj(out_dir, f'F3_scatter_{tag}'), ylabel, extra)
        with open(opj(out_dir, f'methods_{tag}.txt'), 'w') as f:
            f.write(_methods(cfg, level, kept, rep, sub, sim, descs))
        results[level] = dict(runs=df, regions=rep, profiles=prof, stats=stats, similarity=sim)
    return results
