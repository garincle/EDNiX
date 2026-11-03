"""
Multimodal, multi-species regional profiles: morphometry, functional connectivity and glucose metabolism on the
same EDNiX regions (Figure 6, library/multimodal_figure.py).

One value per subject and region, for each modality, from the outputs of the EDNiX pipelines:
  morphometry  anat/native/volumes/labels/<sub>_seg-<atlas>_dseg.nii.gz (4D, one level per volume): volume of the
               region (L + R) in % of the total brain volume of that level (every label except ventricles / CSF).
  FC           func/acpc-func/Stats/Correl_matrix/<atlas>/<fit>/<atlas>_<k>_run_<i>_<fit>_matrix.csv
               (modalities/fMRI/_10_Correl_matrix.py): connectivity strength of a region = mean over the other
               regions of the matrix (N - 1) of the Fisher z of its correlation, STRENGTHS below: 'positive'
               (default; negative z counted as 0, weight-conserving positive strength of Rubinov & Sporns 2011,
               as the hub/degree maps of Buckner 2009 and Tomasi & Volkow 2010), 'absolute' (|z|), 'signed'
               (plain mean, which lets positive and negative correlations cancel); runs averaged, then L and R
               averaged.
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
MODALITIES = ('morphometry', 'FC', 'Q80', 'Q20', 'metabolism')     # FC (matrix strength) and Q80/Q20 (SBA) are alternatives
MODALITY_LABEL = {'morphometry': 'volume (% brain)', 'FC': 'connectivity strength (z)',
                  'Q80': 'retrosplenial SBA, 80th percentile (z)', 'Q20': 'retrosplenial SBA, 20th percentile (z)',
                  'metabolism': 'glucose metabolism (bratio)'}
# Okabe-Ito blue / green + a purple: none close to the Human orange or the Macaque magenta (Fig 6)
MODALITY_LINESTYLE = {'morphometry': '-', 'FC': (0, (4, 1.5)), 'Q80': '-', 'Q20': (0, (4, 1.5)), 'metabolism': (0, (1, 1.2))}
MODALITY_COLORS = {'morphometry': '#0072B2', 'FC': '#009E73', 'Q80': '#009E73', 'Q20': '#7FCBB5', 'metabolism': '#6A3D9A'}


def region_key(name):
    """Base name of a label (hemisphere prefix removed): the region names of the Figure 5 fingerprints."""
    return re.sub(r'^[lLrR]_', '', str(name)).strip()


def canonical(name):
    """Canonical EDNiX region name of a label name (Atlas_nomenclature), for sorting and classification."""
    return AN.canonical_region(name)


def cortical_regions(names, level):
    """Region keys of `names` whose canonical name is a cortical region of `level` (Atlas_nomenclature), grouped
    by lobe in the order of the atlas figure (Atlas_nomenclature.atlas_order), then by name."""
    cortical = set(AN.get_level_regions(level, cortical_only=True))
    keep = sorted(n for n in dict.fromkeys(names) if canonical(n) in cortical)
    return AN.atlas_order(keep)


def brain_regions(names, level):
    """Region keys of `names` whose canonical name is any region of `level` (cortical and subcortical, e.g. striatum,
    thalamus, hippocampus; ventricles and CSF excluded), in the order of the atlas figure."""
    allowed = set(AN.get_level_regions(level))
    keep = sorted(n for n in dict.fromkeys(names)
                  if canonical(n) in allowed and not NON_BRAIN.search(str(canonical(n))) and not NON_BRAIN.search(str(n)))
    return AN.atlas_order(keep)


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


STRENGTHS = ('positive', 'absolute', 'signed')


def node_strength(path, kind='positive'):
    """Connectivity strength of each region in one matrix: mean over the other regions (diagonal excluded) of the
    Fisher z of the correlation, with kind = 'positive' (negative z set to 0, positive strength), 'absolute' (|z|)
    or 'signed' (plain mean: positive and negative correlations cancel)."""
    if kind not in STRENGTHS:
        raise ValueError(f'kind must be one of {STRENGTHS}, got {kind!r}')
    m = pd.read_csv(path, index_col=0)
    r = m.values.astype(float)
    np.fill_diagonal(r, np.nan)
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    if kind == 'positive':
        z = np.where(np.isnan(z), np.nan, np.clip(z, 0, None))
    elif kind == 'absolute':
        z = np.abs(z)
    return pd.Series(np.nanmean(z, axis=1), index=[str(i) for i in m.index])


def fc_table(sources, atlas, level, fit_kind='correlation', strength='positive'):
    """sources: {species: [bids_dir, ...]}. Long table species, bids_dir, subject, region, value (mean z)."""
    rows = []
    for sp, dirs in sources.items():
        for bd in _as_list(dirs):
            files = find_matrices(bd, atlas, level, fit_kind)
            print(f'  [FC] {sp}: {len(files)} matrices in {bd}')
            for f in files:
                s = node_strength(f, strength)
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


# ======================================================================================== SBA (Fig. 4 selection)
def sba_tables(c, species, which='selected'):
    """Connectivity from the retrosplenial seed-based maps of Figure 4 instead of the correlation matrices.
    c = Inputs().fmri() (fmri_fingerprints.context + fmri_group). For every run of the datasets of `species`
    (which = 'selected': the runs kept by the final round of the run selection, 10 per dataset; 'all': every run), the
    regional 80th and 20th percentile of the seed correlation (the 'region|q80' / 'region|q20' values of the
    fingerprints, scaled as in Figure 4). Returns {'Q80': long table, 'Q20': long table} (species, bids_dir,
    subject = dataset:run, region, value); one 'subject' is one run."""
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import profile_stats as fa
    cfg, common, sol = c['cfg'], c['common'], c['species_of_label']
    chan = [str(x).rpartition('|')[2].lower() for x in common]
    base = [region_key(str(x).rpartition('|')[0] or str(x)) for x in common]
    rows = {'Q80': [], 'Q20': []}
    for lab, profs in c['profiles_by_label'].items():
        sp = sol.get(lab)
        if sp not in species:
            continue
        units = c['selected_units'].get(lab, {}) if which == 'selected' else profs
        for key, prof in units.items():
            v = np.asarray(fa._normalize_values(fa._ordered_common_values(prof, common, cfg['VALUE_COL']),
                                                cfg['ZSCORE']), float)
            for i, ch in enumerate(chan):
                m = {'q80': 'Q80', 'q20': 'Q20'}.get(ch)
                if m and np.isfinite(v[i]):
                    rows[m].append(dict(species=sp, bids_dir=lab, subject=f'{lab}:{key}', region=base[i],
                                        value=float(v[i])))
    out = {}
    for m, r in rows.items():
        if not r:
            raise ValueError(f'no {m} value in the selected runs of {species}')
        out[m] = pd.DataFrame(r).groupby(['species', 'bids_dir', 'subject', 'region'], as_index=False)['value'].mean()
    return out


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
        out['FC'] = fc_table(src['FC'], cfg['ATLAS'], cfg['LEVEL'], cfg.get('FIT_KIND', 'correlation'),
                             cfg.get('FC_STRENGTH', 'positive'))
    if src.get('metabolism'):
        out['metabolism'] = metabolism_table(src['metabolism'], cfg['ATLAS'], cfg['LEVEL'], cfg['DESC'],
                                             cfg.get('MIN_VOXELS', 10))
    for k, v in out.items():
        v.insert(0, 'modality', k)
    return out


# ================================================================================================== profiles
def region_coverage(df, regions, species):
    """Number of subjects of each species with a value for each region (long table region, species, n_subjects,
    n_with_value, n_missing). A region absent from a species' atlas, or below the voxel threshold, has n_missing > 0."""
    rows = []
    for sp in species:
        d = df[df.species == sp]
        n = d['subject'].nunique()
        ok = d[d['value'].notna()].groupby('region')['subject'].nunique()
        for r in regions:
            k = int(ok.get(r, 0))
            rows.append(dict(region=r, species=sp, n_subjects=n, n_with_value=k, n_missing=n - k))
    return pd.DataFrame(rows)


def complete_regions(df, regions, species, label=''):
    """Regions of `regions` with a value in EVERY subject of EVERY species (a region missing in one subject is
    dropped, not the subject: dropping subjects listwise empties a species as soon as one region is absent from
    its atlas). Prints the regions dropped and why. Returns (kept regions, coverage table)."""
    cov = region_coverage(df, regions, species)
    bad = cov[cov.n_missing > 0]
    kept = [r for r in regions if r not in set(bad['region'])]
    for sp in species:
        if not (df.species == sp).any():
            raise ValueError(f'{label}: no data for species {sp} (species found: {sorted(df.species.unique())})')
    if len(bad):
        print(f'  [{label}] {bad.region.nunique()} region(s) dropped (missing in some subjects):')
        for r, g in bad.groupby('region', sort=False):
            print('    ' + r + ': ' + ', '.join(f'{x.species} {x.n_missing}/{x.n_subjects}' for x in g.itertuples()))
    print(f'  [{label}] {len(kept)} of {len(regions)} region(s) complete in every subject of {species}')
    return kept, cov


def subject_matrix(df, regions):
    """Wide table index (species, subject), columns regions; subjects still missing a region are dropped (none
    after complete_regions)."""
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
        keep, _ = complete_regions(df, keep, species, mod)
    if len(keep) < 3:
        raise ValueError(f'only {len(keep)} region(s) complete in every modality and species: no radar possible')
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


RADAR_NAMES = {   # shorter spoke names of the radars (full names on the bars of panel b and in the tables)
    'Middle_Temporal_Inferior_temporal_Temporal_pole_(MIPT)': 'Mid./inf. temporal (MIPT)',
    'Ventral_areas_of_the_temporal_lobe_(vent_Temp)': 'Ventral temporal',
    'Insula_and_others_in_lateral_sulcus': 'Insula and lateral sulcus',
    'Auditory_cortex_(Superior_temporal)': 'Auditory cortex',
    'Visual_pre_and_extra_striate_cortex': 'Visual extrastriate',
    'Posterior_medial_cortex_(PMC)': 'Posterior medial (PMC)',
}


def radar_name(region, width=14):
    import textwrap
    return '\n'.join(textwrap.wrap(RADAR_NAMES.get(region, region.replace('_', ' ')), width)) or region


def draw_radar(ax, mods, regions, title, color_title='black', fontsize=7, rlim=None, width=14):
    """One species: one closed line per modality (mean within-subject z across subjects +- SEM band) on radial
    axes = regions (atlas order, first region at the top, counter-clockwise, lobe ring: PIRP.setup_radar); region
    names wrapped every `width` characters."""
    n = len(regions)
    if rlim is None:
        v = [np.abs(mods[m].reindex(columns=regions).mean().values) for m in MODALITIES if m in mods]
        r = float(np.ceil(np.nanmax(np.concatenate(v)) * 1.2 * 2) / 2) if v else 1.0
        rlim = (-r, r)
    ang = PIRP.setup_radar(ax, regions, rlim)          # top, counter-clockwise, lobe ring (atlas figure)
    angc = np.r_[ang, ang[:1]]
    for mod in MODALITIES:
        if mod not in mods:
            continue
        w = mods[mod].reindex(columns=regions)
        m, s = w.mean().values, np.nan_to_num(w.sem().values) if len(w) > 1 else np.zeros(n)
        c = MODALITY_COLORS[mod]
        ax.plot(angc, np.r_[m, m[:1]], color=c, lw=1.4, ls=MODALITY_LINESTYLE[mod], label=f'{mod} (n={len(w)})')
        ax.fill_between(angc, np.r_[m - s, m[:1] - s[:1]], np.r_[m + s, m[:1] + s[:1]], color=c, alpha=0.2, lw=0)
    ax.plot(np.linspace(0, 2 * np.pi, 200), np.zeros(200), 'k--', lw=0.6)
    labels = [radar_name(r, width) for r in regions]
    PIRP.radar_labels(ax, ang, labels, fontsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize - 1)
    top = labels[0].count('\n') + 1 if labels else 0           # lines of the name at the top, above the circle
    ax.set_title(title, fontsize=fontsize + 2, color=color_title, fontweight='bold', pad=4 + 1.25 * fontsize * top)


def draw_radar_overlay(ax, prof, species, regions, colors, fontsize=6, rlim=None, width=18):
    """All species on ONE radar (regions = axes, atlas order, first at the top, counter-clockwise, lobe ring):
    colour = species, line style = modality (volume solid, connectivity dashed, metabolism dotted); thin line =
    mean within-subject z across subjects."""
    n = len(regions)
    mods = {sp: prof[sp] for sp in species if sp in prof}
    if rlim is None:
        v = [np.abs(m.reindex(columns=regions).mean().values) for d in mods.values() for m in d.values()]
        r = float(np.ceil(np.nanmax(np.concatenate(v)) * 1.2 * 2) / 2)
        rlim = (-r, r)
    ang = PIRP.setup_radar(ax, regions, rlim)
    angc = np.r_[ang, ang[:1]]
    for sp, d in mods.items():
        for mod in MODALITIES:
            if mod in d:
                m = d[mod].reindex(columns=regions).mean().values
                ax.plot(angc, np.r_[m, m[:1]], color=colors.get(sp, 'black'), lw=1.1,
                        ls=MODALITY_LINESTYLE[mod], label=f'{sp} {mod}')
    ax.plot(np.linspace(0, 2 * np.pi, 200), np.zeros(200), color='0.3', ls='-', lw=0.4)
    labels = [radar_name(r, width) for r in regions]
    PIRP.radar_labels(ax, ang, labels, fontsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize - 1)
    return rlim


def consensus_profiles(prof, species, regions):
    """{modality: (mean, lo, hi)} across species: mean of the species mean profiles (z) and, as 'variance', the
    range (min-max) of the species means (two species: the two values; standard deviation is not defined well
    with so few species)."""
    out = {}
    for mod in MODALITIES:
        arr = [prof[sp][mod].reindex(columns=regions).mean().values for sp in species if mod in prof.get(sp, {})]
        if len(arr) == len(species):
            a = np.vstack(arr)
            out[mod] = (a.mean(axis=0), a.min(axis=0), a.max(axis=0))
    return out


def draw_radar_consensus(ax, cons, regions, fontsize=6, rlim=None, width=14):
    """Consensus across species, one line per modality (mean of the species means) with the range across species
    shaded; same axes as draw_radar."""
    n = len(regions)
    if rlim is None:
        r = float(np.ceil(max(np.nanmax(np.abs(np.r_[v[1], v[2]])) for v in cons.values()) * 1.2 * 2) / 2)
        rlim = (-r, r)
    ang = PIRP.setup_radar(ax, regions, rlim)
    angc = np.r_[ang, ang[:1]]
    for mod, (m, lo, hi) in cons.items():
        c = MODALITY_COLORS[mod]
        ax.plot(angc, np.r_[m, m[:1]], color=c, lw=1.3, ls=MODALITY_LINESTYLE[mod])
        ax.fill_between(angc, np.r_[lo, lo[:1]], np.r_[hi, hi[:1]], color=c, alpha=0.18, lw=0)
    ax.plot(np.linspace(0, 2 * np.pi, 200), np.zeros(200), 'k--', lw=0.6)
    PIRP.radar_labels(ax, ang, [radar_name(r, width) for r in regions], fontsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize - 1)
    return rlim


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
