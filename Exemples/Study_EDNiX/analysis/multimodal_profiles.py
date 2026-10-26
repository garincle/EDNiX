"""
Multimodal, multi-species regional profiles: morphometry, functional connectivity and glucose metabolism on the
same EDNiX regions (Figure XX, Figure_XX_multimodal.py).

One value per subject and region, for each modality, from the outputs of the EDNiX pipelines:
  morphometry  anat/native/volumes/labels/<sub>_seg-<atlas>_dseg.nii.gz (4D, one level per volume): volume of the
               region (L + R) in % of the total brain volume of that level (every label except ventricles / CSF).
  FC           func/acpc-func/Stats/Correl_matrix/<atlas>/<fit>/<atlas>_<k>_run_<i>_<fit>_matrix.csv
               (modalities/fMRI/_10_Correl_matrix.py): average connectivity of a region = mean Fisher z of its
               correlations with every other region of the matrix (node strength / (N - 1)); runs averaged in z,
               then L and R averaged.
  metabolism   pet/orig/postprocessing/ROI_analysis/*_seg-<atlas>-lvl<k>_desc-bratio[-mean]_ROIs.tsv (step 7 of
               the PET pipeline): brain ratio (activity / mean activity of the brain mask), L and R merged by
               valid-voxel weighting (Plotting/Plot_interspecies_regional_profiles.py).
Sessions of a subject are averaged: one subject = one sample.

LEVELS. `level` is the hierarchy level of Atlas_nomenclature.py (1 = coarsest ... 4), as ATLAS_LEVEL of
Figure_2_anat.py. The pipelines index the 4D atlases from 0: level 2 is volume 1 of the dseg, `lvl1` of the PET
tables and `<atlas>_1_run_*` of the FC matrices.

REGION NAMES. A region is keyed by the base name of its label (label name without the L_ / R_ prefix), the
names of the Figure 5 connectivity fingerprints (base_name of SBA_Region_Profile_Index profiles). The canonical
names of Atlas_nomenclature.DLABEL_NAME_ALIASES serve only to sort the regions and to tell cortical regions
(Atlas_nomenclature level<k>_cortical) from the others.
"""
import os
import re
import glob

import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Tools.atlas_labels import parse_label_file, get_atlas_label_path, find_species_path
from Exemples.Study_EDNiX.atlases import Atlas_nomenclature as AN
from Plotting import Plot_interspecies_regional_profiles as PIRP

opj = os.path.join
opb = os.path.basename

NON_BRAIN = re.compile(r'ventric|csf', re.IGNORECASE)
MODALITIES = ('morphometry', 'FC', 'metabolism')
MODALITY_LABEL = {'morphometry': 'volume (% brain)', 'FC': 'average connectivity (z)',
                  'metabolism': 'glucose metabolism (bratio)'}
# Okabe-Ito colours, distinct from the species colours
MODALITY_COLORS = {'morphometry': '#0072B2', 'FC': '#009E73', 'metabolism': '#D55E00'}


def region_key(name):
    """Base name of a label (hemisphere prefix removed): the region names of the Figure 5 fingerprints."""
    return re.sub(r'^[lLrR]_', '', str(name)).strip()


def canonical(name):
    """Canonical EDNiX region name of a label name (Atlas_nomenclature), for sorting and classification."""
    n = AN._norm(name)
    return AN.DLABEL_NAME_ALIASES_NORM.get(n, region_key(name))


def cortical_regions(names, level):
    """Region keys of `names` whose canonical name is a cortical region of `level` (Atlas_nomenclature), sorted
    by the order of Atlas_nomenclature (lobes grouped), then by name."""
    order = AN.get_level_regions(level, cortical_only=True)
    keep = [n for n in dict.fromkeys(names) if canonical(n) in order]
    return sorted(keep, key=lambda n: (order.index(canonical(n)), n))


def _entity(path, key):
    m = re.search(rf'(?:^|[/_]){key}-([A-Za-z0-9]+)', path)
    return m.group(1) if m else None


def _as_list(x):
    return [x] if isinstance(x, str) else list(x or [])


# ================================================================================================ morphometry
def find_dseg(bids_dir, atlas):
    root = opj(bids_dir, 'derivatives') if os.path.isdir(opj(bids_dir, 'derivatives')) else bids_dir
    pat = opj(root, 'sub-*', '**', 'anat', 'native', 'volumes', 'labels', f'*_seg-{atlas}_dseg.nii.gz')
    return sorted(glob.glob(pat, recursive=True))


def volume_percent(dseg, label_df, level):
    """% of brain volume of every region (L + R) of one 4D dseg at `level` (1-based)."""
    img = nib.load(dseg)
    lab = np.asarray(img.dataobj)
    if lab.ndim == 4:
        if level - 1 >= lab.shape[3]:
            raise ValueError(f'{dseg}: {lab.shape[3]} levels, level {level} requested')
        lab = lab[..., level - 1]
    lab = np.rint(lab).astype(np.int64)
    ids, counts = np.unique(lab[lab > 0], return_counts=True)
    names = dict(zip(label_df['label_id'].astype(int), label_df['region_name']))
    df = pd.DataFrame({'label': ids, 'n': counts})
    df['region'] = [region_key(names.get(int(i), f'label_{i}')) for i in df['label']]
    df = df[~df['region'].map(canonical).str.contains(NON_BRAIN) & ~df['region'].str.contains(NON_BRAIN)]
    tot = df['n'].sum()
    return (df.groupby('region')['n'].sum() / tot * 100.0).rename('value')


def morphometry_table(sources, atlas, level, atlas_lib):
    """sources: {species: [bids_dir, ...]}. Long table species, bids_dir, subject, region, value (% brain)."""
    rows = []
    for sp, dirs in sources.items():
        label_df = parse_label_file(get_atlas_label_path(atlas_lib, find_species_path(atlas_lib, sp), atlas))
        for bd in _as_list(dirs):
            files = find_dseg(bd, atlas)
            print(f'  [morphometry] {sp}: {len(files)} dseg in {bd}')
            for f in files:
                v = volume_percent(f, label_df, level).reset_index()
                v['species'], v['bids_dir'] = sp, opb(os.path.normpath(bd))
                v['subject'], v['session'] = _entity(f, 'sub'), _entity(f, 'ses') or '1'
                rows.append(v)
    if not rows:
        raise FileNotFoundError(f'no {atlas} dseg found for {list(sources)}')
    df = pd.concat(rows, ignore_index=True)
    return _subject_mean(df)


# ========================================================================================================= FC
def find_matrices(bids_dir, atlas, level, fit_kind='correlation'):
    root = opj(bids_dir, 'derivatives') if os.path.isdir(opj(bids_dir, 'derivatives')) else bids_dir
    pat = opj(root, 'sub-*', '**', 'func', 'acpc-func', 'Stats', 'Correl_matrix', atlas, fit_kind,
              f'{atlas}_{level - 1}_run_*_{fit_kind}_matrix.csv')
    return [f for f in sorted(glob.glob(pat, recursive=True))
            if not any(x in opb(f) for x in ('check_fit', 'flattened', 'pval', 'tstat'))]


def node_strength(path):
    """Mean Fisher z of each region with every other region of one matrix (diagonal excluded)."""
    m = pd.read_csv(path, index_col=0)
    r = m.values.astype(float)
    np.fill_diagonal(r, np.nan)
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    return pd.Series(np.nanmean(z, axis=1), index=[str(i) for i in m.index])


def fc_table(sources, atlas, level, fit_kind='correlation'):
    """sources: {species: [bids_dir, ...]}. Long table species, bids_dir, subject, region, value (mean z)."""
    rows = []
    for sp, dirs in sources.items():
        for bd in _as_list(dirs):
            files = find_matrices(bd, atlas, level, fit_kind)
            print(f'  [FC] {sp}: {len(files)} matrices in {bd}')
            for f in files:
                s = node_strength(f)
                v = pd.DataFrame({'region': [region_key(n) for n in s.index], 'value': s.values})
                v = v.groupby('region', as_index=False)['value'].mean()          # L and R averaged
                v['species'], v['bids_dir'] = sp, opb(os.path.normpath(bd))
                v['subject'], v['session'] = _entity(f, 'sub'), _entity(f, 'ses') or '1'
                v['run'] = f
                rows.append(v)
    if not rows:
        raise FileNotFoundError(f'no {atlas} level-{level} {fit_kind} matrix found for {list(sources)}')
    df = pd.concat(rows, ignore_index=True)
    # runs then sessions averaged in z
    df = df.groupby(['species', 'bids_dir', 'subject', 'session', 'region'], as_index=False)['value'].mean()
    return _subject_mean(df)


# ================================================================================================= metabolism
def metabolism_table(sources, atlas, level, desc=('bratio', 'bratio-mean'), min_voxels=10):
    """sources: {species: bids_dir or [bids_dir, ...]}. Long table species, subject, region, value (bratio);
    regions with < min_voxels valid voxels in one hemisphere label of a subject are set to NaN for that subject."""
    rows = []
    for sp, dirs in sources.items():
        for bd in _as_list(dirs):
            df = PIRP.collect({sp: bd}, atlas, level - 1, desc)
            sub = PIRP.subject_level(df, 'merge')
            sub.loc[sub['min_label_voxels'] < min_voxels, 'value'] = np.nan
            sub['region'] = sub['region_key'].map(region_key)
            sub['bids_dir'] = opb(os.path.normpath(bd))
            rows.append(sub[['species', 'bids_dir', 'subject', 'region', 'value']])
    return pd.concat(rows, ignore_index=True)


def _subject_mean(df):
    return df.groupby(['species', 'bids_dir', 'subject', 'region'], as_index=False)['value'].mean()


def collect_modalities(cfg):
    """{modality: long table} for the modalities whose sources are given in cfg['SOURCES'][modality]."""
    out = {}
    src = cfg['SOURCES']
    if src.get('morphometry'):
        out['morphometry'] = morphometry_table(src['morphometry'], cfg['ATLAS'], cfg['LEVEL'], cfg['ATLAS_LIB'])
    if src.get('FC'):
        out['FC'] = fc_table(src['FC'], cfg['ATLAS'], cfg['LEVEL'], cfg.get('FIT_KIND', 'correlation'))
    if src.get('metabolism'):
        out['metabolism'] = metabolism_table(src['metabolism'], cfg['ATLAS'], cfg['LEVEL'], cfg['DESC'],
                                             cfg.get('MIN_VOXELS', 10))
    for k, v in out.items():
        v.insert(0, 'modality', k)
    return out


# ================================================================================================== profiles
def subject_matrix(df, regions):
    """Wide table index (species, subject), columns regions; subjects missing a region are dropped (listwise)."""
    w = df[df['region'].isin(regions)].pivot_table(index=['species', 'subject'], columns='region', values='value')
    w = w.reindex(columns=regions)
    return w.dropna(axis=0, how='any')


def zscore_rows(w):
    """Within-subject z-score across the regions (keeps the shape of each subject's profile)."""
    return w.sub(w.mean(axis=1), axis=0).div(w.std(axis=1, ddof=1), axis=0)


def multimodal_profiles(tables, regions, species):
    """{species: {modality: subjects x regions z-scored within subject}} on the regions present in every
    modality of every species; returns (profiles, regions kept)."""
    keep = list(regions)
    for mod, df in tables.items():
        for sp in species:
            present = set(df.loc[(df.species == sp) & df.value.notna(), 'region'])
            keep = [r for r in keep if r in present]
    prof = {sp: {} for sp in species}
    for mod, df in tables.items():
        for sp in species:
            w = subject_matrix(df[df.species == sp], keep)
            if len(w):
                prof[sp][mod] = zscore_rows(w).droplevel(0)
    return prof, keep


def cross_modal_correlations(prof):
    """Spearman r between the mean (z) profiles of every pair of modalities, per species (descriptive: the
    regions are few and spatially dependent, no p-value)."""
    from scipy.stats import spearmanr
    rows = []
    for sp, mods in prof.items():
        names = [m for m in MODALITIES if m in mods]
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                rows.append(dict(species=sp, modality_a=a, modality_b=b, n_regions=mods[a].shape[1],
                                 spearman_r=spearmanr(mods[a].mean(), mods[b].mean()).statistic))
    return pd.DataFrame(rows)


# ===================================================================================================== panels
def short_name(region, width=14):
    """Display name: the fingerprint name with spaces, wrapped on lines of `width` characters."""
    import textwrap
    return '\n'.join(textwrap.wrap(region.replace('_', ' '), width)) or region


def draw_radar(ax, mods, regions, title, color_title='black', fontsize=7, rlim=None):
    """One species: one closed line per modality (mean within-subject z across subjects +- SEM band) on radial
    axes = regions."""
    n = len(regions)
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    angc = np.r_[ang, ang[:1]]
    for mod in MODALITIES:
        if mod not in mods:
            continue
        w = mods[mod].reindex(columns=regions)
        m, s = w.mean().values, np.nan_to_num(w.sem().values) if len(w) > 1 else np.zeros(n)
        c = MODALITY_COLORS[mod]
        ax.plot(angc, np.r_[m, m[:1]], color=c, lw=1.5, label=f'{mod} (n={len(w)})')
        ax.fill_between(angc, np.r_[m - s, m[:1] - s[:1]], np.r_[m + s, m[:1] + s[:1]], color=c, alpha=0.2, lw=0)
    ax.plot(np.linspace(0, 2 * np.pi, 200), np.zeros(200), 'k--', lw=0.6)
    ax.set_theta_offset(np.pi / 2)          # first region at the top, clockwise (as the fingerprints)
    ax.set_theta_direction(-1)
    if rlim:
        ax.set_ylim(*rlim)
    ax.set_xticks(ang)
    ax.set_xticklabels([short_name(r) for r in regions], fontsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize - 1)
    ax.set_title(title, fontsize=fontsize + 2, color=color_title, fontweight='bold', pad=6)


def autocrop(img, tol=0.02):
    """Image without its uniform border (background colour of the corner pixel), so that the brain fills the panel."""
    rgb = img[..., :3].astype(float)
    if rgb.max() > 1:
        rgb = rgb / 255.0
    diff = np.abs(rgb - rgb[0, 0]).max(axis=2) > tol
    if img.shape[-1] == 4:
        diff &= img[..., 3] > 0
    rows, cols = np.where(diff.any(axis=1))[0], np.where(diff.any(axis=0))[0]
    if not len(rows) or not len(cols):
        return img
    return img[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1]


def draw_surface(ax, png, title, fontsize=8):
    ax.axis('off')
    if png and os.path.isfile(png):
        ax.imshow(autocrop(plt.imread(png)))
    else:
        ax.text(0.5, 0.5, f'surface not rendered\n{opb(png) if png else ""}', ha='center', va='center',
                fontsize=fontsize, color='0.4', transform=ax.transAxes)
    ax.set_title(title, fontsize=fontsize + 1, fontweight='bold')
