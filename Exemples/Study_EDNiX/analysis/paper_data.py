"""
paper_data.py -- datasets and data loading shared by the paper scripts.

One place for the BIDS folders, the atlas, the output root and the
"one value per subject" rule, so that Figure 2, Supplementary Figs 4-8 and the
sensitivity analyses read exactly the same table.

Loading is the slow part (it parses every xlsx / json of every subject): the
collected tables are cached in <PAPER_DIR>/_cache/. Set refresh=True (or delete
the cache) after re-running EDNiX on a dataset.
"""
import os
import sys
import pickle

import numpy as np
import pandas as pd

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

opj = os.path.join

ATLAS_LIB = '/home/cgarin/PycharmProjects/EDNiX/Atlases_library'
ATLAS_NAME = 'EDNIxCSC'
NEWICK_FILE = '/scratch2/EDNiX/Liste des espèces.nwk'
HERCULANO_XLSX = '/scratch2/EDNiX/Herculano-Houzel_mammalian_brain_dataset_complete.xlsx'

FIGURES_ROOT = '/scratch2/EDNiX/results/EDNiX_Figures'
PAPER_DIR = opj(FIGURES_ROOT, 'Paper')                  # every paper figure lands here
SENSITIVITY_DIR = opj(FIGURES_ROOT, 'sensitivity_analysis')
PHYLO_DIR = opj(FIGURES_ROOT, 'phylo_preliminary')
CACHE_DIR = opj(PAPER_DIR, '_cache')

# Output of the Fig 5 pipeline (Launcher_Fig5.BASE_CONFIG["OUTPUT_ROOT"], percentiles q=20)
FIG5_ROOT = opj(FIGURES_ROOT, 'Fig_5_V15__quant20')
FIG5_CONTROL_ROOT = FIG5_ROOT + '_control_common_target'
FIG3_DIR = opj(FIGURES_ROOT, 'Figure_3_longit')

# ── datasets ────────────────────────────────────────────────────────────────
# Anatomy (Figure 2, Supp. Figs 4, 6-8): every dataset with a T1/T2 processed by EDNiX.
ANAT_DATASETS = {
    'Mouse': ['/scratch2/EDNiX/Mouse/BIDS_Grandjean2'],
    'Rat': ['/scratch2/EDNiX/Rat/BIDS_Grandjean'],
    'Mouselemur': ['/scratch2/EDNiX/Mouselemur/BIDS_Garin'],
    'Marmoset': ['/scratch2/EDNiX/Marmoset/BIDS_Tian',
                 '/scratch2/EDNiX/Marmoset/BIDS_Hori'],
    'Macaque': ['/scratch2/EDNiX/Macaque/BIDS_Tremblay',
                '/scratch2/EDNiX/Macaque/BIDS_BenHamed',
                '/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin'],
    'Human': ['/scratch2/EDNiX/Human/BIDS_Merida',
              '/scratch2/EDNiX/Human/BIDS_Park',
              '/scratch2/EDNiX/Human/BIDS_Castrillon',
              '/scratch2/EDNiX/Human/BIDS_Babayan'],
    'Dog': ['/scratch2/EDNiX/Dog/BIDS_Boch_K9'],
}
# Resting-state fMRI (Supp. Fig. 5, Figs 4-5): datasets with BOLD runs.
FUNC_DATASETS = {
    'Mouse': ['/scratch2/EDNiX/Mouse/BIDS_Grandjean2'],
    'Rat': ['/scratch2/EDNiX/Rat/BIDS_Grandjean'],
    'Mouselemur': ['/scratch2/EDNiX/Mouselemur/BIDS_Garin'],
    'Marmoset': ['/scratch2/EDNiX/Marmoset/BIDS_Tian',
                 '/scratch2/EDNiX/Marmoset/BIDS_Hori'],
    'Macaque': ['/scratch2/EDNiX/Macaque/BIDS_Tremblay',
                '/scratch2/EDNiX/Macaque/BIDS_BenHamed',
                '/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin'],
    'Human': ['/scratch2/EDNiX/Human/BIDS_Castrillon',
              '/scratch2/EDNiX/Human/BIDS_Babayan'],
    'Dog': ['/scratch2/EDNiX/Dog/BIDS_Boch_K9'],
}

# Cortical regions drawn per atlas level (regions with vertices on the pial
# surface in every species; see the exclusion notes in Figure_2_anat history).
CORTICAL_REGIONS = {
    1: ['Isocortex'],
    2: ['Frontal_lobe', 'Insular_lobe', 'Temporal_lobe', 'Parietal_lobe', 'Occipital_lobe'],
    3: ['PFC', 'Motor_premotor', 'Insula', 'Auditory_Cx', 'Olfactory_Cx',
        'Somatosensory_Cx', 'Post_parietal_Cx', 'Striate_Cx', 'Extra_striate_Cx'],
    4: ['BA_6', 'BA_4_M1', 'Insular_Cx', 'Belt_A2', 'BA_27', 'BA_35_36',
        'BA_1_3_S1', 'S2', 'BA_18_V2'],
}

MORPHO_METRICS = [  # (modality, column, axis label)
    ('surface', 'surface_area_mm2', 'Surface area (mm$^2$)'),
    ('volume', 'volume_mm3', 'Volume (mm$^3$)'),
    ('thickness', 'thickness_mm', 'Thickness (mm)'),
]


def species_config(datasets):
    return {sp: {'bids_dirs': list(bds), 'list_to_keep': [], 'list_to_remove': []}
            for sp, bds in datasets.items()}


def _atlas_paths(datasets):
    from Exemples.Study_EDNiX.analysis.ednix_bids_tools import (
        get_atlas_label_path, find_species_path)
    paths, frags = {}, {}
    for sp in datasets:
        try:
            frags[sp] = find_species_path(ATLAS_LIB, sp)
            paths[sp] = get_atlas_label_path(ATLAS_LIB, frags[sp], ATLAS_NAME,
                                             prefer_statslut=False)
        except (FileNotFoundError, ValueError) as e:
            print(f"  [WARN] {sp}: {e}")
    return paths, frags


def collect(datasets, extract, tag, refresh=False):
    """collect_multi_species() on `datasets`, cached as <CACHE_DIR>/<tag>.pkl."""
    cache = opj(CACHE_DIR, f'{tag}.pkl')
    if not refresh and os.path.isfile(cache):
        with open(cache, 'rb') as fh:
            print(f"  [cache] {cache}")
            return pickle.load(fh)
    from Exemples.Study_EDNiX.analysis.ednix_bids_tools import collect_multi_species
    paths, frags = _atlas_paths(datasets)
    data = collect_multi_species(
        species_config(datasets), regions_of_interest=None, extract=tuple(extract),
        atlas_name=ATLAS_NAME, atlas_label_paths=paths, fit_kind='correlation',
        atlas_library_root=ATLAS_LIB, species_atlas_fragments=frags)
    data = {k: v for k, v in data.items() if k in extract}
    os.makedirs(CACHE_DIR, exist_ok=True)
    with open(cache, 'wb') as fh:
        pickle.dump(data, fh)
    return data


def bids_name(path):
    return os.path.basename(str(path).rstrip('/'))


def per_subject(df, value_col):
    """One row per (species, dataset, subject, region, hemisphere, atlas level):
    mean over sessions. Sessions of one animal are not independent samples."""
    if df is None or df.empty or value_col not in df.columns:
        return df
    keys = [c for c in ('species', 'bids_dir', 'subject', 'region', 'hemisphere', 'atlas_level')
            if c in df.columns]
    out = df.groupby(keys, as_index=False).agg(**{value_col: (value_col, 'mean')})
    if 'session' in df.columns:
        n = df.groupby(keys)['session'].nunique().reset_index(name='n_sessions')
        out = out.merge(n, on=keys, how='left')
    return out


def bilateral(df, value_col):
    """Bilateral value per subject-region: the 'bilateral' row when it exists,
    otherwise left + right. Hemispheres present only on one side are dropped
    (a single hemisphere is not a bilateral value)."""
    if df is None or df.empty:
        return df
    if 'hemisphere' not in df.columns:
        return df.assign(hemisphere='bilateral')
    keys = [c for c in ('species', 'bids_dir', 'subject', 'session', 'region', 'atlas_level')
            if c in df.columns]
    bil = df[df['hemisphere'] == 'bilateral']
    lr = df[df['hemisphere'].isin(['left', 'right'])]
    lr = lr[~lr.set_index(keys).index.isin(bil.set_index(keys).index)]
    nside = lr.groupby(keys)['hemisphere'].transform('nunique')
    lr = lr[nside == 2]
    if value_col == 'thickness_mm':
        lr_sum = lr.groupby(keys, as_index=False)[value_col].mean()   # thickness: L/R mean
    else:
        lr_sum = lr.groupby(keys, as_index=False)[value_col].sum()    # area, volume: L + R
    lr_sum['hemisphere'] = 'bilateral'
    bil = bil.groupby(keys, as_index=False)[value_col].mean().assign(hemisphere='bilateral')
    return pd.concat([bil, lr_sum], ignore_index=True)


def load_morphometry(level=None, refresh=False):
    """{'surface'|'volume'|'thickness': bilateral, one-row-per-subject DataFrame}.
    bids_dir is reduced to the folder name (BIDS_xxx)."""
    raw = collect(ANAT_DATASETS, ('surface', 'volume', 'thickness'), 'morphometry', refresh)
    out = {}
    for mod, col, _ in MORPHO_METRICS:
        df = raw.get(mod)
        if df is None or df.empty or col not in df.columns:
            out[mod] = pd.DataFrame()
            continue
        df = df.copy()
        df['bids_dir'] = df['bids_dir'].map(bids_name)
        if level is not None and 'atlas_level' in df.columns:
            df = df[df['atlas_level'] == level]
        out[mod] = per_subject(bilateral(df, col), col)
    return out


def select_regions(df, regions):
    """Rows of `df` whose region matches one of `regions` (EDNiX name matching),
    with the 'region' column set to the requested name."""
    from Exemples.Study_EDNiX.analysis.ednix_bids_tools import _df_region_mask
    parts = []
    for r in regions:
        m = _df_region_mask(df['region'], [r])
        if m.any():
            parts.append(df[m].assign(region=r))
    return pd.concat(parts, ignore_index=True) if parts else df.iloc[0:0]


def load_qc(datasets, tag, refresh=False):
    """QC table, one row per (subject, session[, run]) -- QC is not averaged:
    a bad session must stay visible."""
    raw = collect(datasets, ('qc',), tag, refresh)
    df = raw.get('qc', pd.DataFrame()).copy()
    if not df.empty:
        df['bids_dir'] = df['bids_dir'].map(bids_name)
    return df


def robust_outliers(values, threshold=2.5):
    """Boolean mask of |x - median| / (1.4826 * MAD) > threshold (Leys et al. 2013)."""
    v = np.asarray(values, float)
    med = np.nanmedian(v)
    mad = np.nanmedian(np.abs(v - med))
    if not np.isfinite(mad) or mad == 0:
        return np.zeros(v.shape, bool)
    return np.abs(v - med) / (1.4826 * mad) > threshold


def fig5_panel(name, root=None):
    """Path of a panel written by the Fig 5 pipeline (06_figures/Paper, Fig5, Fig6 or
    supplementary), '<name>.png'; None when it does not exist (yet)."""
    root = root or FIG5_ROOT
    for sub in ('Paper', 'Fig5', 'Fig6', 'supplementary', ''):
        p = opj(root, '06_figures', sub, f'{name}.png')
        if os.path.isfile(p):
            return p
    print(f"  [WARN] panel {name}.png not found under {root}/06_figures -- run "
          f"Fig5/Launcher_Fig5.py (RUN_PART='figures' is enough)")
    return None
