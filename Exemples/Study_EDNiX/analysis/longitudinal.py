"""
EDNiX Longitudinal — Run Script
=================================
Dataset-specific formatting + launch of longitudinal analysis.

Edit only this file for each new dataset.
All analysis logic lives in ednix_longitudinal_tools.py.

Output structure
----------------
<OUT_DIR>/
  Macaque/surface/    <region>_trajectory.png  +  _pct_change.png  +  stats.xlsx
  Macaque/thickness/  ...
  Macaque/volume/     ...
  Human/surface/      ...
  ...
  longitudinal_summary.xlsx
"""

import os
import sys
import warnings
import pandas as pd
import numpy as np

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')
from Plotting.ednix_bids_tools import (
    collect_multi_species,
    get_atlas_label_path,
    find_species_path,
)
from Statistics.Group_fMRI.EDNiX_longitudinal import run_longitudinal

opj = os.path.join

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit here
# ══════════════════════════════════════════════════════════════════════════════

ATLAS_LIB   = '/home/cgarin/PycharmProjects/EDNiX/Atlases_library'
ATLAS_NAME  = 'EDNIxCSC'
ATLAS_LEVEL = 1
OUT_DIR     = '/scratch2/EDNiX/results/longitudinal_analysis'
REGIONS     = ['Isocortex', 'Allocortex', 'Periallocortex']
MIN_SESSIONS = 3

MACAQUE_BIDS_DIRS = [
    '/scratch2/EDNiX/Macaque/BIDS_Cdt_Garin',
]
MACAQUE_EXCEL = '/scratch2/EDNiX/Macaque/BIDS_Cdt_Garin/Garin_macaque.xlsx'

HUMAN_BIDS_DIRS = [
    '/scratch2/EDNiX/Human/ds004856/',
]

# participants.tsv column that contains age at each wave
# The table has columns: AgeMRI_W1, AgeMRI_W2, etc.
HUMAN_AGE_PREFIX = 'AgeMRI_W'   # prefix before wave number


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Load morphometry via ednix_bids_tools
# ══════════════════════════════════════════════════════════════════════════════

def _build_atlas_paths(species_list):
    paths, frags = {}, {}
    for sp in species_list:
        try:
            frag        = find_species_path(ATLAS_LIB, sp)
            frags[sp]   = frag
            paths[sp]   = get_atlas_label_path(
                ATLAS_LIB, frag, ATLAS_NAME, prefer_statslut=False)
            print(f"  [{sp}] atlas: {paths[sp]}")
        except Exception as e:
            print(f"  [WARN] {sp}: {e}")
    return paths, frags


def load_morphometry(species_bids: dict, atlas_label_paths: dict) -> dict:
    """
    Call collect_multi_species and return the raw data dict.
    Keys: 'surface', 'volume', 'thickness' (DataFrames or None).
    """
    species_config = {
        sp: {'bids_dirs': dirs, 'list_to_keep': [], 'list_to_remove': []}
        for sp, dirs in species_bids.items()
    }
    data = collect_multi_species(
        species_config,
        regions_of_interest=REGIONS,
        extract=('surface', 'volume', 'thickness'),
        atlas_name=ATLAS_NAME,
        atlas_label_paths=atlas_label_paths,
        atlas_library_root=ATLAS_LIB,
    )
    # collect_multi_species returns a dict; values can be list of DFs or a DF
    # Normalise to single DataFrame per modality
    out = {}
    for key in ('surface', 'volume', 'thickness'):
        val = data.get(key)
        if val is None:
            out[key] = pd.DataFrame()
        elif isinstance(val, list):
            out[key] = pd.concat(val, ignore_index=True) if val else pd.DataFrame()
        elif isinstance(val, pd.DataFrame):
            out[key] = val
        else:
            out[key] = pd.DataFrame()
    return out


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2A — Macaque: format age table
# Dataset-specific: reads Garin_macaque.xlsx
# ══════════════════════════════════════════════════════════════════════════════

def format_macaque_age(excel_path: str) -> pd.DataFrame:
    """
    Reads Garin_macaque.xlsx and returns a standardised age table:

        subject | session | age_years

    Logic:
      - 'subject' = animal ID (matches sub-<ID> in BIDS)
      - 'session' = chronological session number (1, 2, 3...)
      - 'age_years' = (session_date - birthdate) / 365.25
    """
    ani = pd.read_excel(excel_path, sheet_name='animalinfo')
    mri = pd.read_excel(excel_path, sheet_name='MRIsessions')
    try:
        mri2 = pd.read_excel(excel_path, sheet_name='MRIsessions12')
        mri  = pd.concat([mri, mri2], ignore_index=True)
    except Exception:
        pass

    # Merge animal info (birthdate) with session dates
    meta = pd.merge(ani, mri, on='ID', how='inner')

    # ── Find birthdate column ────────────────────────────────────────────────
    dob_col = next((c for c in meta.columns
                    if c.upper() in ('DOB', 'DATE_NAISSANCE', 'BIRTHDATE', 'NAISSANCE')),
                   None)
    if dob_col is None:
        raise KeyError(f"No birthdate column in animalinfo. Columns: {list(ani.columns)}")
    meta['_dob'] = pd.to_datetime(meta[dob_col], errors='coerce')

    # ── Find session date column ─────────────────────────────────────────────
    dos_col = next((c for c in meta.columns
                    if c.upper() in ('DOS', 'DATE', 'SESSION_DATE', 'SCAN_DATE')),
                   None)
    if dos_col is None:
        raise KeyError(f"No session date column. Columns: {list(mri.columns)}")
    meta['_dos'] = pd.to_datetime(meta[dos_col], errors='coerce')

    meta = meta.dropna(subset=['_dob', '_dos'])
    meta['age_days']  = (meta['_dos'] - meta['_dob']).dt.days.astype(float)
    meta['age_years'] = meta['age_days'] / 365.25

    # Subject = ID (strip sub- prefix issues)
    meta['subject'] = meta['ID'].astype(str)

    # Session = chronological rank per animal
    meta = meta.sort_values(['subject', '_dos'])
    meta['session'] = (meta.groupby('subject').cumcount() + 1).astype(str)

    result = meta[['subject', 'session', 'age_years']].dropna().reset_index(drop=True)
    print(f"  [Macaque age] {len(result)} sessions  "
          f"{result['subject'].nunique()} subjects  "
          f"age range [{result['age_years'].min():.1f}, "
          f"{result['age_years'].max():.1f}] years")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2B — Human: format age table from participants.tsv
# Dataset-specific: wide format with AgeMRI_W1, AgeMRI_W2 ...
# ══════════════════════════════════════════════════════════════════════════════

def format_human_age(bids_dirs: list,
                     age_prefix: str = 'AgeMRI_W') -> pd.DataFrame:
    """
    Reads participants.tsv (wide format) and returns standardised age table:

        subject | session | age

    Wide columns like AgeMRI_W1, AgeMRI_W2 → session 1, 2, ...
    participant_id is normalised to strip 'sub-' prefix.
    """
    frames = []
    for bd in bids_dirs:
        tsv = opj(bd, 'participants.tsv')
        if not os.path.exists(tsv):
            warnings.warn(f"participants.tsv not found: {tsv}")
            continue
        df = pd.read_csv(tsv, sep='\t')

        # Normalise subject column
        if 'participant_id' in df.columns:
            df['subject'] = (df['participant_id']
                             .str.replace('sub-', '', regex=False)
                             .astype(str))
        elif 'subject' not in df.columns:
            warnings.warn(f"No subject/participant_id in {tsv}")
            continue

        # Find wave columns: e.g. AgeMRI_W1, AgeMRI_W2, ...
        wave_cols = sorted([c for c in df.columns if c.startswith(age_prefix)])
        if not wave_cols:
            # Fallback: single 'age' column → all sessions share same age
            if 'age' in df.columns:
                long = df[['subject', 'age']].copy()
                long = long.rename(columns={'age': 'age'})
                long['session'] = '1'
                frames.append(long[['subject', 'session', 'age']])
            else:
                warnings.warn(f"No age columns found in {tsv}")
            continue

        # Melt wide → long
        id_cols  = ['subject']
        age_wide = df[id_cols + wave_cols].copy()
        long     = age_wide.melt(id_vars='subject',
                                  value_vars=wave_cols,
                                  var_name='wave_col',
                                  value_name='age')
        # Extract wave number → session label
        long['session'] = (long['wave_col']
                           .str.replace(age_prefix, '', regex=False)
                           .astype(str))
        long = long.drop(columns=['wave_col'])
        # Replace 'n/a' strings and drop missing
        long['age'] = pd.to_numeric(long['age'], errors='coerce')
        long = long.dropna(subset=['age'])
        frames.append(long[['subject', 'session', 'age']])

    if not frames:
        return pd.DataFrame(columns=['subject', 'session', 'age'])

    result = pd.concat(frames, ignore_index=True).drop_duplicates()
    print(f"  [Human age] {len(result)} sessions  "
          f"{result['subject'].nunique()} subjects  "
          f"age range [{result['age'].min():.1f}, {result['age'].max():.1f}] years")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    os.makedirs(OUT_DIR, exist_ok=True)

    all_stats = []

    # ── Atlas paths ───────────────────────────────────────────────────────────
    atlas_paths, _ = _build_atlas_paths(['Macaque'])

    # ══════════════════════════════════════════════════════════════════════════
    # MACAQUE
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 64)
    print("  MACAQUE")
    print("=" * 64)

    mac_morph = load_morphometry(
        {'Macaque': MACAQUE_BIDS_DIRS}, atlas_paths)

    mac_age = format_macaque_age(MACAQUE_EXCEL)

    for mod, metric_col, metric_label in [
        ('surface',   'surface_area_mm2', 'Surface area (mm²)'),
        ('thickness', 'thickness_mm',     'Cortical thickness (mm)'),
        ('volume',    'volume_mm3',        'Volume (mm³)'),
    ]:
        df_mod = mac_morph.get(mod, pd.DataFrame())
        if df_mod is None or df_mod.empty:
            print(f"  [Macaque/{mod}] no data")
            continue

        stats = run_longitudinal(
            label        = f'Macaque_{mod}',
            df_morph     = df_mod,
            df_age       = mac_age,        # subject | session | age_years
            metric_col   = metric_col,
            metric_label = metric_label,
            age_col      = 'age_years',
            age_label    = 'Age (years)',
            out_dir      = opj(OUT_DIR, 'Macaque', mod),
            atlas_level  = ATLAS_LEVEL,
            hemisphere   = 'left',
            regions      = REGIONS,
            min_sessions = MIN_SESSIONS,
        )
        all_stats.append(stats)

    # ══════════════════════════════════════════════════════════════════════════
    # HUMAN
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 64)
    print("  HUMAN")
    print("=" * 64)
    os.makedirs(OUT_DIR, exist_ok=True)

    all_stats = []

    # ── Atlas paths ───────────────────────────────────────────────────────────
    atlas_paths, _ = _build_atlas_paths(['Human'])

    hum_morph = load_morphometry(
        {'Human': HUMAN_BIDS_DIRS}, atlas_paths)

    hum_age = format_human_age(HUMAN_BIDS_DIRS, age_prefix=HUMAN_AGE_PREFIX)

    for mod, metric_col, metric_label in [
        ('surface',   'surface_area_mm2', 'Surface area (mm²)'),
        ('thickness', 'thickness_mm',     'Cortical thickness (mm)'),
        ('volume',    'volume_mm3',        'Volume (mm³)'),
    ]:
        df_mod = hum_morph.get(mod, pd.DataFrame())
        if df_mod is None or df_mod.empty:
            print(f"  [Human/{mod}] no data")
            continue

        stats = run_longitudinal(
            label        = f'Human_{mod}',
            df_morph     = df_mod,
            df_age       = hum_age,         # subject | session | age
            metric_col   = metric_col,
            metric_label = metric_label,
            age_col      = 'age',
            age_label    = 'Age (years)',
            out_dir      = opj(OUT_DIR, 'Human', mod),
            atlas_level  = ATLAS_LEVEL,
            hemisphere   = 'left',
            regions      = REGIONS,
            min_sessions = MIN_SESSIONS,
        )
        all_stats.append(stats)

    # ── Combined summary ──────────────────────────────────────────────────────
    if all_stats:
        summary = pd.concat([s for s in all_stats if not s.empty],
                            ignore_index=True)
        out_path = opj(OUT_DIR, 'longitudinal_summary.xlsx')
        summary.to_excel(out_path, index=False)
        print(f"\n  Summary → {out_path}")

    print("\n  Done.")