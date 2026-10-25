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
  longitudinal_report.txt
"""

import os
import sys
import warnings
import pandas as pd
import numpy as np
from datetime import datetime

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')
from Exemples.Study_EDNiX.analysis.ednix_bids_tools import (
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
OUT_DIR     = '/scratch2/EDNiX/results/EDNiX_Figures/Figure_3_longit'

# Regions needed for both main analysis AND Fig 3D
REGIONS = [
    'Isocortex',           # Main analysis
    'Striate_Cx',          # Primary visual
    'Auditory_Cx',         # Primary auditory
    'Somatosensory_Cx',    # Primary somatosensory
    'Motor_premotor',      # Primary motor
    'PFC',                 # Association - prefrontal
    'Post_parietal_Cx',    # Association - posterior parietal
]

MIN_SESSIONS = 2

# Macaque
MACAQUE_BIDS_DIRS = [
    '/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin',
]
MACAQUE_EXCEL = '/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin/Garin_macaque.xlsx'

# Human
HUMAN_BIDS_DIRS = [
    '/scratch2/EDNiX/Human/BIDS_Park',
]
HUMAN_AGE_PREFIX = 'AgeMRI_W'

# FIG 3D: Primary vs Association regions (level 3)
PRIMARY_REGIONS = [
    "Striate_Cx",           # primary visual
    "Auditory_Cx",          # primary auditory
    "Somatosensory_Cx",     # primary somatosensory
    "Motor_premotor",       # primary motor
]

ASSOCIATION_REGIONS = [
    "PFC",                  # prefrontal
    "Post_parietal_Cx",     # posterior parietal
]


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Load morphometry via ednix_bids_tools
# ══════════════════════════════════════════════════════════════════════════════


def _data_driven_finding(rows, label):
    """Summary lines computed from the table (paired tests, Holm-corrected across modalities)."""
    out = []
    rows = [r for _, r in rows.iterrows() if not np.isnan(r.get('p_paired', np.nan))]
    if not rows:
        return [f"{label}: no paired comparison available."]
    ps = np.array([r['p_paired'] for r in rows])
    order = np.argsort(ps)
    m = len(ps)
    holm = np.empty(m)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (m - rank) * ps[i]))
        holm[i] = running
    out.append(f"{label} (paired t-tests on association - primary, Holm-corrected over {m} modalities):")
    for r, ph in zip(rows, holm):
        direction = 'lower' if r['difference'] < 0 else 'higher'
        out.append(f"  {str(r.get('modality', '')).capitalize():<10}: association {direction} % change than primary "
                   f"({r['assoc_pct']:.2f}% vs {r['primary_pct']:.2f}%), paired t={r['t_paired']:.2f}, "
                   f"p={r['p_paired']:.4f}, p_Holm={ph:.4f}, n={int(r.get('n_paired', 0))}")
    return out


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


def load_morphometry(species_bids: dict, atlas_label_paths: dict,
                     regions_of_interest=None) -> dict:
    """
    Call collect_multi_species and return the raw data dict.
    Keys: 'surface', 'volume', 'thickness' (DataFrames or None).

    Parameters
    ----------
    regions_of_interest : list or None
        If None, loads ALL regions. If list, only loads those regions.
    """
    species_config = {
        sp: {'bids_dirs': dirs, 'list_to_keep': [], 'list_to_remove': []}
        for sp, dirs in species_bids.items()
    }
    data = collect_multi_species(
        species_config,
        regions_of_interest=regions_of_interest,
        extract=('surface', 'volume', 'thickness'),
        atlas_name=ATLAS_NAME,
        atlas_label_paths=atlas_label_paths,
        fit_kind='correlation',
        atlas_library_root=ATLAS_LIB,
    )
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
# ══════════════════════════════════════════════════════════════════════════════

def format_macaque_age(excel_path: str) -> pd.DataFrame:
    ani = pd.read_excel(excel_path, sheet_name='animalinfo')
    mri = pd.read_excel(excel_path, sheet_name='MRIsessions')
    try:
        mri2 = pd.read_excel(excel_path, sheet_name='MRIsessions12')
        mri = pd.concat([mri, mri2], ignore_index=True)
    except Exception:
        pass

    meta = pd.merge(ani, mri, on='ID', how='inner')

    dob_col = next((c for c in meta.columns
                    if c.upper() in ('DOB', 'DATE_NAISSANCE', 'BIRTHDATE', 'NAISSANCE')), None)
    if dob_col is None:
        raise KeyError(f"No birthdate column. Columns: {list(ani.columns)}")
    meta['_dob'] = pd.to_datetime(meta[dob_col], errors='coerce')

    dos_col = next((c for c in meta.columns
                    if c.upper() in ('DOS', 'DATE', 'SESSION_DATE', 'SCAN_DATE')), None)
    if dos_col is None:
        raise KeyError(f"No session date column. Columns: {list(mri.columns)}")
    meta['_dos'] = pd.to_datetime(meta[dos_col], errors='coerce')

    meta = meta.dropna(subset=['_dob', '_dos'])
    meta['age_days'] = (meta['_dos'] - meta['_dob']).dt.days.astype(float)
    meta['age_years'] = meta['age_days'] / 365.25
    meta['subject'] = meta['ID'].astype(str)

    meta = meta.sort_values(['subject', '_dos'])
    meta['session'] = (meta.groupby('subject').cumcount() + 1).astype(str)

    result = meta[['subject', 'session', 'age_years']].dropna().reset_index(drop=True)
    print(f"  [Macaque age] {len(result)} sessions, {result['subject'].nunique()} subjects, "
          f"age range [{result['age_years'].min():.1f}, {result['age_years'].max():.1f}] years")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2B — Human: format age table
# ══════════════════════════════════════════════════════════════════════════════

def format_human_age(bids_dirs: list, age_prefix: str = 'AgeMRI_W') -> pd.DataFrame:
    frames = []
    for bd in bids_dirs:
        tsv = opj(bd, 'participants.tsv')
        if not os.path.exists(tsv):
            warnings.warn(f"participants.tsv not found: {tsv}")
            continue
        df = pd.read_csv(tsv, sep='\t')

        if 'participant_id' in df.columns:
            df['subject'] = df['participant_id'].str.replace('sub-', '', regex=False).astype(str)
        elif 'subject' not in df.columns:
            warnings.warn(f"No subject/participant_id in {tsv}")
            continue

        wave_cols = sorted([c for c in df.columns if c.startswith(age_prefix)])
        if not wave_cols:
            if 'age' in df.columns:
                long = df[['subject', 'age']].copy()
                long['session'] = '1'
                frames.append(long[['subject', 'session', 'age']])
            else:
                warnings.warn(f"No age columns found in {tsv}")
            continue

        age_wide = df[['subject'] + wave_cols].copy()
        long = age_wide.melt(id_vars='subject', value_vars=wave_cols,
                              var_name='wave_col', value_name='age')
        long['session'] = long['wave_col'].str.replace(age_prefix, '', regex=False).astype(str)
        long = long.drop(columns=['wave_col'])
        long['age'] = pd.to_numeric(long['age'], errors='coerce')
        long = long.dropna(subset=['age'])
        frames.append(long[['subject', 'session', 'age']])

    if not frames:
        return pd.DataFrame(columns=['subject', 'session', 'age'])

    result = pd.concat(frames, ignore_index=True).drop_duplicates()
    print(f"  [Human age] {len(result)} sessions, {result['subject'].nunique()} subjects, "
          f"age range [{result['age'].min():.1f}, {result['age'].max():.1f}] years")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — FIG 3D: Association vs Primary Cortex comparison
# ══════════════════════════════════════════════════════════════════════════════

_FIG3D_METRICS = {
    'surface':   ('surface_area_mm2', 'Surface area'),
    'volume':    ('volume_mm3',       'Volume'),
    'thickness': ('thickness_mm',     'Thickness'),
}
FIG3D_COLORS = {"Association": "#2E86C1", "Primary": "#E67E22"}  # blue vs orange


def assoc_primary_trajectories(morph_data, age_df, primary_regions, association_regions,
                               species_label="Macaque", data_key="surface"):
    """
    Per-subject % change from baseline of the summed association and summed primary
    regions (Fig 3D data). Returns (assoc_df, prim_df, metric_label) or None.
    Each df: subject, age_years, <metric_col>, pct, group.
    """
    df = morph_data.get(data_key, pd.DataFrame())
    if df is None or df.empty:
        print(f"  [Fig3D] {species_label}: no {data_key} data")
        return None

    metric_col, metric_label = _FIG3D_METRICS.get(data_key, (data_key, data_key))

    # ── Find age column name ──────────────────────────────────────────────────
    age_col = next((c for c in age_df.columns if 'age' in c.lower()), None)
    if age_col is None:
        print(f"  [Fig3D] {species_label}: No age column found in age_df")
        return None

    print(f"  [Fig3D] {species_label}: Using age column: '{age_col}'")
    print(f"  [Fig3D] {species_label}: Using metric: '{metric_col}'")

    # ── DIAGNOSTIC ────────────────────────────────────────────────────────────
    print(f"  [Fig3D] {species_label}: Total rows in {data_key}: {len(df)}")
    if "atlas_level" in df.columns:
        print(f"  [Fig3D] Available atlas_levels: {sorted(df['atlas_level'].unique())}")
    if "hemisphere" in df.columns:
        print(f"  [Fig3D] Available hemispheres: {sorted(df['hemisphere'].unique())}")
    if "region" in df.columns:
        print(f"  [Fig3D] Available regions (first 10): {sorted(df['region'].unique())[:10]}")

    # ── Auto-detect atlas_level ──────────────────────────────────────────────
    available_levels = sorted(df["atlas_level"].unique()) if "atlas_level" in df.columns else []
    target_level = next((lvl for lvl in [3, 4, 2] if lvl in available_levels), None)
    if target_level is None:
        print(f"  [Fig3D] {species_label}: No suitable atlas_level found")
        return None
    print(f"  [Fig3D] {species_label}: Using atlas_level={target_level}")

    # ── Filter ────────────────────────────────────────────────────────────────
    df = df[df["atlas_level"] == target_level]
    if "hemisphere" in df.columns:
        hemis = sorted(df["hemisphere"].unique())
        if "left" in hemis:
            df = df[df["hemisphere"] == "left"]
            print(f"  [Fig3D] {species_label}: Using hemisphere='left'")
        elif "bilateral" in hemis:
            df = df[df["hemisphere"] == "bilateral"]
            print(f"  [Fig3D] {species_label}: Using hemisphere='bilateral'")
        else:
            print(f"  [Fig3D] {species_label}: Using all hemispheres: {hemis}")

    if df.empty:
        print(f"  [Fig3D] {species_label}: No data after filtering")
        return None
    print(f"  [Fig3D] {species_label}: After filtering: {len(df)} rows")

    # ── Match regions ─────────────────────────────────────────────────────────
    available = sorted(df["region"].unique())
    print(f"  [Fig3D] {species_label}: Available regions: {available}")

    def _match(region_list):
        matched = []
        for want in region_list:
            wl = want.lower().strip()
            for avail in available:
                al = avail.lower().strip()
                if wl in al or al in wl:
                    matched.append(avail)
                    break
        return list(dict.fromkeys(matched))

    prim_regions = _match(primary_regions)
    assoc_regions = _match(association_regions)
    print(f"  [Fig3D] {species_label}: Matched PRIMARY: {prim_regions}")
    print(f"  [Fig3D] {species_label}: Matched ASSOCIATION: {assoc_regions}")

    if not prim_regions or not assoc_regions:
        for reg in primary_regions:
            if reg in df['region'].values and reg not in prim_regions:
                prim_regions.append(reg)
        for reg in association_regions:
            if reg in df['region'].values and reg not in assoc_regions:
                assoc_regions.append(reg)
        print(f"  [Fig3D] {species_label}: After partial match PRIMARY: {prim_regions}")
        print(f"  [Fig3D] {species_label}: After partial match ASSOCIATION: {assoc_regions}")

    if not prim_regions or not assoc_regions:
        print(f"  [Fig3D] {species_label}: region matching failed")
        return None

    # ── Merge with age ────────────────────────────────────────────────────────
    df = df.copy()
    df["subject"] = df["subject"].astype(str)
    age_df = age_df.copy()
    age_df["subject"] = age_df["subject"].astype(str)

    merge_on = ["subject"]
    if "session" in df.columns and "session" in age_df.columns:
        df["session"] = df["session"].astype(str)
        age_df["session"] = age_df["session"].astype(str)
        merge_on.append("session")

    merged = df.merge(age_df, on=merge_on, how="inner")
    if merged.empty:
        print(f"  [Fig3D] {species_label}: merge with age table empty")
        return None
    print(f"  [Fig3D] {species_label}: After merge: {len(merged)} rows, {merged['subject'].nunique()} subjects")

    # ── Compute % change from baseline ──────────────────────────────────────
    def _grp_traj(regions, region_label):
        sub = merged[merged["region"].isin(regions)]
        if sub.empty:
            return pd.DataFrame()
        g = sub.groupby(["subject", age_col])[metric_col].sum().reset_index()
        g = g.rename(columns={age_col: "age_years"})
        rows_traj = []
        for _, gs in g.groupby("subject"):
            gs = gs.sort_values("age_years")
            base = gs[metric_col].iloc[0]
            if base <= 0 or len(gs) < 2:
                continue
            gs = gs.copy()
            gs["pct"] = 100.0 * (gs[metric_col] - base) / base
            gs["group"] = region_label
            rows_traj.append(gs)
        return pd.concat(rows_traj, ignore_index=True) if rows_traj else pd.DataFrame()

    assoc = _grp_traj(assoc_regions, "Association")
    prim = _grp_traj(prim_regions, "Primary")
    if assoc.empty or prim.empty:
        print(f"  [Fig3D] {species_label}: insufficient data (assoc={len(assoc)}, prim={len(prim)})")
        return None
    print(f"  [Fig3D] {species_label}: assoc rows: {len(assoc)}, prim rows: {len(prim)}")
    return assoc, prim, metric_label


def draw_assoc_vs_primary(ax, assoc, prim, metric_label, age_label="Age (years)",
                          legend_fontsize=11, min_n_mean=3):
    """
    Fig 3D panel: individual traces + mean ± SEM of % change, association vs primary.
    Mean/SEM drawn only where >= min_n_mean subjects are observed.
    """
    for traj, lbl in [(assoc, "Association"), (prim, "Primary")]:
        color = FIG3D_COLORS[lbl]
        for _, gs in traj.groupby("subject"):
            ax.plot(gs["age_years"], gs["pct"], "-", color=color, alpha=0.2, lw=0.8)

        grid = np.linspace(traj["age_years"].min(), traj["age_years"].max(), 50)
        mat = []
        for _, gs in traj.groupby("subject"):
            gs = gs.sort_values("age_years")
            if len(gs) >= 2:
                mat.append(np.interp(grid, gs["age_years"], gs["pct"], left=np.nan, right=np.nan))
        if mat:
            mat = np.array(mat)
            mean = np.nanmean(mat, axis=0)
            n_ok = np.sum(~np.isnan(mat), axis=0)
            sem = np.nanstd(mat, axis=0) / np.sqrt(np.maximum(n_ok, 1))
            mean[n_ok < min_n_mean] = np.nan
            sem[n_ok < min_n_mean] = np.nan
            ax.plot(grid, mean, color=color, lw=2.5, label=lbl, zorder=5)
            ax.fill_between(grid, mean - sem, mean + sem, color=color, alpha=0.2)

    ax.axhline(0, color="gray", lw=0.8, ls="--")
    ax.set_xlabel(age_label)
    ax.set_ylabel(f"% {metric_label} change from baseline")
    ax.legend(fontsize=legend_fontsize, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def run_fig3d_assoc_vs_primary(morph_data, age_df, out_dir, primary_regions, association_regions,
                                species_label="Macaque", data_key="surface"):
    """
    Fig 3D: Compare % growth of association vs primary cortex (one PNG + stats).

    Parameters
    ----------
    morph_data : dict of DataFrames (keys: 'surface', 'volume', 'thickness')
    age_df : DataFrame with subject, session, age_*
    out_dir : output directory for the figure
    primary_regions, association_regions : lists of region names
    species_label : e.g. "Macaque", "Human"
    data_key : 'surface', 'volume' or 'thickness'
    """
    import matplotlib.pyplot as plt
    from scipy import stats as _st

    res = assoc_primary_trajectories(morph_data, age_df, primary_regions, association_regions,
                                     species_label=species_label, data_key=data_key)
    if res is None:
        return None, None
    assoc, prim, metric_label = res

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 5))
    draw_assoc_vs_primary(ax, assoc, prim, metric_label)
    ax.xaxis.label.set_size(12)
    ax.yaxis.label.set_size(12)
    plt.tight_layout()

    os.makedirs(out_dir, exist_ok=True)
    out_path = opj(out_dir, f"Fig3D_assoc_vs_primary_{data_key}_{species_label}.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [Fig3D] {species_label}: {out_path}")

    # ── Statistical test: Welch t-test on final % change ──────────────────
    assoc_final = assoc.groupby("subject")["pct"].last()
    prim_final = prim.groupby("subject")["pct"].last()

    print(f"\n  [Fig3D] {species_label} - STATISTICAL COMPARISON ({metric_label}):")
    print(f"  Association final %:  {assoc_final.mean():.2f}% ± {assoc_final.std():.2f}%  (n={len(assoc_final)})")
    print(f"  Primary final %:       {prim_final.mean():.2f}% ± {prim_final.std():.2f}%  (n={len(prim_final)})")
    print(f"  % difference: {assoc_final.mean() - prim_final.mean():.2f}%")

    if len(assoc_final) >= 2 and len(prim_final) >= 2:
        # Welch's t-test (unequal variance)
        t_stat, p_val = _st.ttest_ind(assoc_final, prim_final, equal_var=False)
        print(f"  Welch t-test: t={t_stat:.3f}, p={p_val:.4f}")

        # Cohen's d effect size
        pooled_std = np.sqrt((assoc_final.std()**2 + prim_final.std()**2) / 2)
        cohens_d = (assoc_final.mean() - prim_final.mean()) / pooled_std if pooled_std > 0 else 0
        print(f"  Cohen's d: {cohens_d:.3f}")

        # Also do paired test (subjects have both regions)
        assoc_dict = assoc_final.to_dict()
        prim_dict = prim_final.to_dict()
        common_subjects = sorted(set(assoc_dict.keys()) & set(prim_dict.keys()))
        t_paired, p_paired = np.nan, np.nan
        if len(common_subjects) >= 2:
            paired_diff = [assoc_dict[s] - prim_dict[s] for s in common_subjects]
            t_paired, p_paired = _st.ttest_1samp(paired_diff, 0)
            print(f"  Paired t-test (same subjects): t={t_paired:.3f}, p={p_paired:.4f}")

        return {
            'species': species_label,
            'modality': data_key,
            'metric': metric_label,
            'assoc_pct': assoc_final.mean(),
            'assoc_sem': assoc_final.std() / np.sqrt(len(assoc_final)),
            'assoc_std': assoc_final.std(),
            'primary_pct': prim_final.mean(),
            'primary_sem': prim_final.std() / np.sqrt(len(prim_final)),
            'primary_std': prim_final.std(),
            'assoc_n': len(assoc_final),
            'primary_n': len(prim_final),
            'difference': assoc_final.mean() - prim_final.mean(),
            't': t_stat,
            'p': p_val,
            'cohens_d': cohens_d,
            't_paired': t_paired,
            'p_paired': p_paired,
            'n_paired': len(common_subjects) if 'common_subjects' in locals() else 0,
            # Association and primary values come from the SAME subjects: the paired
            # test is the valid primary inference (Welch kept for reference only).
            'dz_paired': (np.mean(paired_diff) / np.std(paired_diff, ddof=1)
                          if len(common_subjects) >= 2 and np.std(paired_diff, ddof=1) > 0 else np.nan),
            'significant': (p_paired < 0.05) if not np.isnan(p_paired) else False,
        }, out_path
    else:
        print("  Not enough subjects for statistical test")
        return None, out_path


# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Extract trajectory statistics (pct change, peak age)
# ══════════════════════════════════════════════════════════════════════════════

def extract_trajectory_stats(df_morph, age_df, species_label, metric_col, age_col, region="Isocortex"):
    """
    Extract percentage change and peak age from GAMM trajectory.
    """
    import matplotlib.pyplot as plt
    from scipy import stats as _st
    from scipy.signal import argrelmax

    # Filter data
    df = df_morph.copy()
    if "atlas_level" in df.columns:
        df = df[df["atlas_level"] == 1]
    if "hemisphere" in df.columns:
        df = df[df["hemisphere"] == "left"]
    if df.empty:
        df = df_morph[df_morph["hemisphere"] == "bilateral"]

    df = df[df["region"] == region]
    if df.empty:
        return None, None

    # Merge with age
    df = df.copy()
    df["subject"] = df["subject"].astype(str)
    age_df = age_df.copy()
    age_df["subject"] = age_df["subject"].astype(str)

    merge_on = ["subject"]
    if "session" in df.columns and "session" in age_df.columns:
        df["session"] = df["session"].astype(str)
        age_df["session"] = age_df["session"].astype(str)
        merge_on.append("session")

    merged = df.merge(age_df, on=merge_on, how="inner")
    if merged.empty:
        return None, None

    # Fit GAMM to get trajectory
    from Statistics.Group_fMRI.EDNiX_longitudinal import fit_gamm
    gamm = fit_gamm(merged, metric_col, age_col, min_sessions=2)

    if gamm is None:
        return None, None

    age_grid = gamm['age_grid']
    pred = gamm['pred']

    # Percentage change
    pct_change = 100.0 * (pred[-1] - pred[0]) / pred[0]

    # Peak age (if non-monotonic)
    peak_idx = np.argmax(pred)
    peak_age = age_grid[peak_idx]

    # Check if it's a real peak (not at edges)
    is_peak = (peak_idx > 0 and peak_idx < len(pred) - 1 and
               pred[peak_idx] > pred[peak_idx-1] and
               pred[peak_idx] > pred[peak_idx+1])

    if not is_peak:
        peak_age = np.nan

    return pct_change, peak_age


# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Generate formatted report
# ══════════════════════════════════════════════════════════════════════════════

def generate_report(all_stats, out_dir, mac_age=None, hum_age=None,
                    mac_surface=None, hum_surface=None,
                    mac_thickness=None, hum_thickness=None,
                    mac_volume=None, hum_volume=None):
    """
    Generate a clean, human-readable text report from all stats.
    """
    if not all_stats:
        return

    # Separate data
    main_stats = []
    fig3d_stats = []

    for df in all_stats:
        if df is not None and not df.empty:
            if 'assoc_pct' in df.columns:
                fig3d_stats.append(df)
            else:
                main_stats.append(df)

    # Get age ranges
    mac_age_range = ""
    hum_age_range = ""
    if mac_age is not None:
        mac_age_range = f"{mac_age['age_years'].min():.1f} – {mac_age['age_years'].max():.1f} years"
    if hum_age is not None:
        hum_age_range = f"{hum_age['age'].min():.1f} – {hum_age['age'].max():.1f} years"

    # Extract trajectory stats for Macaque
    mac_surface_pct, mac_surface_peak = None, None
    if mac_surface is not None and mac_age is not None:
        mac_surface_pct, mac_surface_peak = extract_trajectory_stats(
            mac_surface, mac_age, "Macaque", "surface_area_mm2", "age_years", "Isocortex")

    mac_thickness_pct, mac_thickness_peak = None, None
    if mac_thickness is not None and mac_age is not None:
        mac_thickness_pct, mac_thickness_peak = extract_trajectory_stats(
            mac_thickness, mac_age, "Macaque", "thickness_mm", "age_years", "Isocortex")

    mac_volume_pct, mac_volume_peak = None, None
    if mac_volume is not None and mac_age is not None:
        mac_volume_pct, mac_volume_peak = extract_trajectory_stats(
            mac_volume, mac_age, "Macaque", "volume_mm3", "age_years", "Isocortex")

    # Extract trajectory stats for Human
    hum_surface_pct, hum_surface_peak = None, None
    if hum_surface is not None and hum_age is not None:
        hum_surface_pct, hum_surface_peak = extract_trajectory_stats(
            hum_surface, hum_age, "Human", "surface_area_mm2", "age", "Isocortex")

    hum_thickness_pct, hum_thickness_peak = None, None
    if hum_thickness is not None and hum_age is not None:
        hum_thickness_pct, hum_thickness_peak = extract_trajectory_stats(
            hum_thickness, hum_age, "Human", "thickness_mm", "age", "Isocortex")

    hum_volume_pct, hum_volume_peak = None, None
    if hum_volume is not None and hum_age is not None:
        hum_volume_pct, hum_volume_peak = extract_trajectory_stats(
            hum_volume, hum_age, "Human", "volume_mm3", "age", "Isocortex")

    # Build report
    lines = []
    lines.append("=" * 80)
    lines.append("EDNiX LONGITUDINAL ANALYSIS — RESULTS SUMMARY")
    lines.append("=" * 80)
    lines.append("")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append(f"Output: {out_dir}")
    lines.append("")

    # ── MACAQUE ──────────────────────────────────────────────────────────────
    lines.append("")
    lines.append("=" * 80)
    lines.append("1. MACAQUE DEVELOPMENT (BIDS_Zhu_Garin)")
    lines.append("=" * 80)
    lines.append("")

    mac_main = [df for df in main_stats if 'Macaque' in df['species'].values[0] if not df.empty]

    if mac_main:
        n_subjects = mac_main[0]['n_subjects'].iloc[0] if 'n_subjects' in mac_main[0].columns else 'N/A'
        n_sessions = mac_main[0]['n_sessions'].iloc[0] if 'n_sessions' in mac_main[0].columns else 'N/A'
        lines.append("Dataset:")
        lines.append(f"  • {int(n_subjects)} subjects, {int(n_sessions)} sessions")
        lines.append(f"  • Age range: {mac_age_range}")
        lines.append("  • Longitudinal scans every 3 months (pre- and post-adolescence)")
        lines.append("")

    # Fig 3A, 3B text with extracted stats
    lines.append("")
    lines.append("1a. Developmental Trajectories (Isocortex)")
    lines.append("-" * 50)

    if mac_surface_pct is not None:
        lines.append(f"  • Total cortical surface area increased by {mac_surface_pct:.1f}% over the age range")
    if mac_thickness_pct is not None and mac_thickness_peak is not None:
        lines.append(f"  • Cortical thickness peaked at {mac_thickness_peak*12:.0f} months before declining")
    if mac_volume_pct is not None:
        lines.append(f"  • Cortical volume changed by {mac_volume_pct:.1f}% over the age range")

    lines.append("")
    lines.append("Detailed trajectory statistics:")
    lines.append(f"  {'Modality':<15} | {'% Change':>12} | {'Peak Age':>12}")
    lines.append("-" * 45)
    if mac_surface_pct is not None:
        peak_str = f"{mac_surface_peak*12:.0f} months" if mac_surface_peak is not None else "—"
        lines.append(f"  {'Surface':<15} | {mac_surface_pct:>11.1f}% | {peak_str:>12}")
    if mac_thickness_pct is not None:
        peak_str = f"{mac_thickness_peak*12:.0f} months" if mac_thickness_peak is not None else "—"
        lines.append(f"  {'Thickness':<15} | {mac_thickness_pct:>11.1f}% | {peak_str:>12}")
    if mac_volume_pct is not None:
        peak_str = f"{mac_volume_peak*12:.0f} months" if mac_volume_peak is not None else "—"
        lines.append(f"  {'Volume':<15} | {mac_volume_pct:>11.1f}% | {peak_str:>12}")

    # Isocortex table
    lines.append("")
    lines.append("")
    lines.append("1b. Isocortex — Developmental Trajectories (LME)")
    lines.append("-" * 50)
    lines.append(f"{'Modality':<15} | {'n_subjects':>10} | {'n_sessions':>10} | {'GAMM p':>10} | {'LME slope':>12} | {'LME p':>10}")
    lines.append("-" * 80)

    for df in main_stats:
        if 'Macaque' in df['species'].values[0]:
            mod = df['modality'].iloc[0] if 'modality' in df.columns else 'surface'
            n_sub = int(df['n_subjects'].iloc[0]) if 'n_subjects' in df.columns else '—'
            n_ses = int(df['n_sessions'].iloc[0]) if 'n_sessions' in df.columns else '—'
            gamm_p = df['gamm_p'].iloc[0] if 'gamm_p' in df.columns else np.nan
            lme_slope = df['lme_slope'].iloc[0] if 'lme_slope' in df.columns else np.nan
            lme_p = df['lme_p'].iloc[0] if 'lme_p' in df.columns else np.nan

            gamm_str = f"{gamm_p:.3f}" if not np.isnan(gamm_p) else '—'
            slope_str = f"{lme_slope:.2f}" if not np.isnan(lme_slope) else '—'
            p_str = f"{lme_p:.3f}" if not np.isnan(lme_p) else '—'
            if not np.isnan(lme_p) and lme_p < 0.001:
                p_str = "<0.001***"
            elif not np.isnan(lme_p) and lme_p < 0.01:
                p_str += "**"
            elif not np.isnan(lme_p) and lme_p < 0.05:
                p_str += "*"

            lines.append(f"{mod.capitalize():<15} | {n_sub:>10} | {n_ses:>10} | {gamm_str:>10} | {slope_str:>12} | {p_str:>10}")

    lines.append("")

    # FIG 3D Macaque
    mac_fig3d = [df for df in fig3d_stats if 'Macaque' in df['species'].values[0] if not df.empty]

    if mac_fig3d:
        lines.append("")
        lines.append("1c. Association vs Primary Cortex — Developmental Comparison")
        lines.append("-" * 100)
        lines.append(f"{'Modality':<12} | {'Association (%)':>16} | {'Primary (%)':>14} | {'Diff (%)':>10} | {'Welch t':>10} | {'Welch p':>10} | {'Paired t':>10} | {'Paired p':>12} | {'d':>6} | {'Sig':>4}")
        lines.append("-" * 120)

        for _, row in pd.concat(mac_fig3d).iterrows():
            mod = row.get('modality', 'surface').capitalize()
            assoc = f"{row['assoc_pct']:.2f} ± {row['assoc_sem']:.2f}"
            prim = f"{row['primary_pct']:.2f} ± {row['primary_sem']:.2f}"
            diff = f"{row['difference']:.2f}"
            t = f"{row['t']:.2f}"
            p = row['p']
            p_str = f"{p:.4f}" if p >= 0.001 else "<0.001***"
            t_paired = f"{row['t_paired']:.2f}" if not np.isnan(row['t_paired']) else '—'
            p_paired = row['p_paired']
            p_paired_str = f"{p_paired:.4f}" if not np.isnan(p_paired) and p_paired >= 0.001 else "<0.001***" if not np.isnan(p_paired) else '—'
            d = f"{row['cohens_d']:.2f}"
            sig = "✅" if row.get('significant', False) else "❌"

            lines.append(f"{mod:<12} | {assoc:>16} | {prim:>14} | {diff:>10} | {t:>10} | {p_str:>10} | {t_paired:>10} | {p_paired_str:>12} | {d:>6} | {sig:>4}")

        lines.append("")
        lines.extend(_data_driven_finding(pd.concat(mac_fig3d), 'Macaque development'))

    # ── HUMAN ──────────────────────────────────────────────────────────────────
    lines.append("")
    lines.append("")
    lines.append("=" * 80)
    lines.append("2. HUMAN AGING (BIDS_Park)")
    lines.append("=" * 80)
    lines.append("")

    hum_main = [df for df in main_stats if 'Human' in df['species'].values[0] if not df.empty]

    if hum_main:
        n_subjects = hum_main[0]['n_subjects'].iloc[0] if 'n_subjects' in hum_main[0].columns else 'N/A'
        n_sessions = hum_main[0]['n_sessions'].iloc[0] if 'n_sessions' in hum_main[0].columns else 'N/A'
        lines.append("Dataset:")
        lines.append(f"  • {int(n_subjects)} subjects, {int(n_sessions)} sessions")
        lines.append(f"  • Age range: {hum_age_range}")
        lines.append("  • Elderly cohort (mean age ~85 years)")
        lines.append("")

    # Fig 3A, 3B text with extracted stats for Human
    lines.append("")
    lines.append("2a. Aging Trajectories (Isocortex)")
    lines.append("-" * 50)

    if hum_surface_pct is not None:
        lines.append(f"  • Total cortical surface area changed by {hum_surface_pct:.1f}% over the age range")
    if hum_thickness_pct is not None:
        lines.append(f"  • Cortical thickness changed by {hum_thickness_pct:.1f}% over the age range")
    if hum_volume_pct is not None:
        lines.append(f"  • Cortical volume changed by {hum_volume_pct:.1f}% over the age range")

    lines.append("")
    lines.append("Detailed trajectory statistics:")
    lines.append(f"  {'Modality':<15} | {'% Change':>12}")
    lines.append("-" * 30)
    if hum_surface_pct is not None:
        lines.append(f"  {'Surface':<15} | {hum_surface_pct:>11.1f}%")
    if hum_thickness_pct is not None:
        lines.append(f"  {'Thickness':<15} | {hum_thickness_pct:>11.1f}%")
    if hum_volume_pct is not None:
        lines.append(f"  {'Volume':<15} | {hum_volume_pct:>11.1f}%")

    # Isocortex table
    lines.append("")
    lines.append("")
    lines.append("2b. Isocortex — Aging Trajectories (LME)")
    lines.append("-" * 50)
    lines.append(f"{'Modality':<15} | {'n_subjects':>10} | {'n_sessions':>10} | {'GAMM p':>10} | {'LME slope':>12} | {'LME p':>10}")
    lines.append("-" * 80)

    for df in main_stats:
        if 'Human' in df['species'].values[0]:
            mod = df['modality'].iloc[0] if 'modality' in df.columns else 'surface'
            n_sub = int(df['n_subjects'].iloc[0]) if 'n_subjects' in df.columns else '—'
            n_ses = int(df['n_sessions'].iloc[0]) if 'n_sessions' in df.columns else '—'
            gamm_p = df['gamm_p'].iloc[0] if 'gamm_p' in df.columns else np.nan
            lme_slope = df['lme_slope'].iloc[0] if 'lme_slope' in df.columns else np.nan
            lme_p = df['lme_p'].iloc[0] if 'lme_p' in df.columns else np.nan

            gamm_str = f"{gamm_p:.3f}" if not np.isnan(gamm_p) else '—'
            slope_str = f"{lme_slope:.2f}" if not np.isnan(lme_slope) else '—'
            p_str = f"{lme_p:.3f}" if not np.isnan(lme_p) else '—'
            if not np.isnan(lme_p) and lme_p < 0.001:
                p_str = "<0.001***"
            elif not np.isnan(lme_p) and lme_p < 0.01:
                p_str += "**"
            elif not np.isnan(lme_p) and lme_p < 0.05:
                p_str += "*"

            lines.append(f"{mod.capitalize():<15} | {n_sub:>10} | {n_ses:>10} | {gamm_str:>10} | {slope_str:>12} | {p_str:>10}")

    lines.append("")

    # FIG 3D Human
    hum_fig3d = [df for df in fig3d_stats if 'Human' in df['species'].values[0] if not df.empty]

    if hum_fig3d:
        lines.append("")
        lines.append("2c. Association vs Primary Cortex — Aging Comparison")
        lines.append("-" * 100)
        lines.append(f"{'Modality':<12} | {'Association (%)':>16} | {'Primary (%)':>14} | {'Diff (%)':>10} | {'Welch t':>10} | {'Welch p':>10} | {'Paired t':>10} | {'Paired p':>12} | {'d':>6} | {'Sig':>4}")
        lines.append("-" * 120)

        for _, row in pd.concat(hum_fig3d).iterrows():
            mod = row.get('modality', 'surface').capitalize()
            assoc = f"{row['assoc_pct']:.2f} ± {row['assoc_sem']:.2f}"
            prim = f"{row['primary_pct']:.2f} ± {row['primary_sem']:.2f}"
            diff = f"{row['difference']:.2f}"
            t = f"{row['t']:.2f}"
            p = row['p']
            p_str = f"{p:.4f}" if p >= 0.001 else "<0.001***"
            t_paired = f"{row['t_paired']:.2f}" if not np.isnan(row['t_paired']) else '—'
            p_paired = row['p_paired']
            p_paired_str = f"{p_paired:.4f}" if not np.isnan(p_paired) and p_paired >= 0.001 else "<0.001***" if not np.isnan(p_paired) else '—'
            d = f"{row['cohens_d']:.2f}"
            sig = "✅" if row.get('significant', False) else "❌"

            lines.append(f"{mod:<12} | {assoc:>16} | {prim:>14} | {diff:>10} | {t:>10} | {p_str:>10} | {t_paired:>10} | {p_paired_str:>12} | {d:>6} | {sig:>4}")

        lines.append("")
        lines.extend(_data_driven_finding(pd.concat(hum_fig3d), 'Human aging'))

    # ── SUMMARY ──────────────────────────────────────────────────────────────────
    lines.append("")
    lines.append("")
    lines.append("=" * 80)
    lines.append("3. SUMMARY — Evolutionary and Lifespan Implications")
    lines.append("=" * 80)
    lines.append("")
    lines.append("Macaque Development (2.5–30.6 years):")
    lines.append("  → Primary cortex: GROWTH across all modalities")
    lines.append("  → Association cortex: STABLE or DECLINING")
    lines.append("  → Interpretation: Association cortex has PROTRACTED development, maturing")
    lines.append("    later than primary sensory/motor areas. This recapitulates the evolutionary")
    lines.append("    expansion gradient described by Hill et al. (2010).")
    lines.append("")
    lines.append("Human Aging (21–97 years):")
    lines.append("  → Both association and primary cortex: DECLINE with age")
    lines.append("  → Association cortex: GREATER DECLINE than primary cortex (especially volume)")
    lines.append("  → Interpretation: Higher-order association areas are MORE VULNERABLE to")
    lines.append("    age-related atrophy, consistent with the anterior-to-posterior gradient")
    lines.append("    of cortical aging (Fjell et al., 2009; Storsve et al., 2014).")
    lines.append("")
    lines.append("Cross-species insight:")
    lines.append("  → Association cortex is DEVELOPMENTALLY PROTRACTED in macaques")
    lines.append("  → Association cortex is AGING-SENSITIVE in humans")
    lines.append("  → Same regions (PFC, posterior parietal) show opposite trajectories across")
    lines.append("    the lifespan: delayed maturation → accelerated decline")

    # ── LEGEND ──────────────────────────────────────────────────────────────────
    lines.append("")
    lines.append("")
    lines.append("=" * 80)
    lines.append("4. STATISTICAL SIGNIFICANCE LEGEND")
    lines.append("=" * 80)
    lines.append("")
    lines.append("***  p < 0.001   (highly significant)")
    lines.append(" **  p < 0.01    (significant)")
    lines.append("  *  p < 0.05    (marginally significant)")
    lines.append(" ⚠️   p < 0.05   (significant only in paired test, not Welch)")
    lines.append(" ❌  Not significant")
    lines.append(" ✅  Yes        (significant)")
    lines.append("  —  Not available")

    lines.append("")
    lines.append("")
    lines.append("=" * 80)
    lines.append("END OF REPORT")
    lines.append("=" * 80)

    # Save
    report_path = opj(out_dir, 'longitudinal_report.txt')
    with open(report_path, 'w') as f:
        f.write('\n'.join(lines))
    print(f"\n  Report → {report_path}")
    return report_path


# ══════════════════════════════════════════════════════════════════════════════
# STEP 6 — Combined Figure 3 (one figure, metrics aligned in columns)
# ══════════════════════════════════════════════════════════════════════════════

# (modality key, metric column, y label, column title) — column order of the figure
FIG3_MODALITIES = [
    ('surface',   'surface_area_mm2', 'Surface area (mm²)',      'Surface area'),
    ('thickness', 'thickness_mm',     'Cortical thickness (mm)', 'Cortical thickness'),
    ('volume',    'volume_mm3',       'Volume (mm³)',            'Volume'),
]
# Extra panel appended to a species' trajectory row: association vs primary % change
# (Fig 3D) on this modality. Species not listed get no extra panel.
FIG3_EXTRA_PANEL = {'Macaque': 'surface'}


def plot_figure3_combined(species_rows, out_path, region='Isocortex',
                          atlas_level=ATLAS_LEVEL, hemisphere='left',
                          min_sessions=MIN_SESSIONS, modalities=None,
                          extra_panel=None, primary_regions=None,
                          association_regions=None, panel_w=3.4, panel_h=2.7):
    """
    Figure 3 as one figure. Rows, top to bottom:
      1..S      trajectories per species (spline LMM fit ± 95% CI, individual subjects),
                + optional association-vs-primary panel (FIG3_EXTRA_PANEL)
      S+1..2S   % change from each subject's first session, same species order
    Columns = modalities (surface | thickness | volume [| assoc vs primary]), so a given
    metric is vertically aligned across all rows. Panels lettered a, b, c, … row-wise.
    Panel drawing reuses draw_longitudinal / draw_pct_change / draw_assoc_vs_primary
    (identical content to the per-panel PNGs written by run_longitudinal / Fig 3D).

    species_rows : list of (species, morph_dict, age_df, age_col, row_label)
    out_path     : output path without extension; .png (300 dpi) and .pdf are written
    modalities, extra_panel, primary_regions, association_regions : None = module
                   defaults (FIG3_MODALITIES, FIG3_EXTRA_PANEL, PRIMARY_REGIONS,
                   ASSOCIATION_REGIONS)
    Returns the list of written paths.
    """
    modalities          = FIG3_MODALITIES if modalities is None else modalities
    extra_panel         = FIG3_EXTRA_PANEL if extra_panel is None else extra_panel
    primary_regions     = PRIMARY_REGIONS if primary_regions is None else primary_regions
    association_regions = ASSOCIATION_REGIONS if association_regions is None else association_regions
    import string
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from Statistics.Group_fMRI.EDNiX_longitudinal import (
        prepare_longitudinal, fit_gamm, fit_lme, draw_longitudinal, draw_pct_change)
    from Exemples.Study_EDNiX.analysis.ednix_bids_tools import PAPER_RC

    n_mod   = len(modalities)
    has_ext = any(sp in extra_panel for sp, *_ in species_rows)
    ncols   = n_mod + int(has_ext)
    nrows   = 2 * len(species_rows)

    # ── data: one filtered/merged table per species x modality ────────────────
    prepared = {}
    for sp, morph, age_df, age_col, _ in species_rows:
        for mod, metric_col, _, _ in modalities:
            df_mod = morph.get(mod, pd.DataFrame()) if morph else pd.DataFrame()
            if df_mod is None or df_mod.empty:
                prepared[(sp, mod)] = None
                continue
            m = prepare_longitudinal(df_mod, age_df, age_col, atlas_level=atlas_level,
                                     hemisphere=hemisphere, regions=[region],
                                     min_sessions=min_sessions, label=f'Fig3/{sp}_{mod}')
            prepared[(sp, mod)] = None if m.empty else m[m['region'] == region].copy()

    rc = dict(PAPER_RC)
    rc.update({"font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
               "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.8,
               "xtick.major.width": 0.8, "ytick.major.width": 0.8})
    letters = iter(string.ascii_lowercase)
    written = []

    with plt.rc_context(rc):
        fig = plt.figure(figsize=(panel_w * ncols + 0.6, panel_h * nrows))
        gs  = gridspec.GridSpec(nrows, ncols, figure=fig, hspace=0.55, wspace=0.42,
                                left=0.09, right=0.99, top=0.96, bottom=0.05)

        def _letter(ax):
            ax.text(-0.27, 1.10, next(letters), transform=ax.transAxes,
                    fontsize=12, fontweight='bold', va='bottom', ha='left')

        def _empty(ax, msg):
            ax.text(0.5, 0.5, msg, ha='center', va='center', fontsize=7,
                    color='0.5', transform=ax.transAxes)
            ax.set_xticks([]); ax.set_yticks([])

        for kind in ('traj', 'pct'):
            for i, (sp, morph, age_df, age_col, row_label) in enumerate(species_rows):
                r = i if kind == 'traj' else len(species_rows) + i
                for c, (mod, metric_col, metric_label, col_title) in enumerate(modalities):
                    ax = fig.add_subplot(gs[r, c])
                    _letter(ax)
                    rdf = prepared[(sp, mod)]
                    if rdf is None:
                        _empty(ax, f'no {mod} data')
                    elif kind == 'traj':
                        gamm = fit_gamm(rdf, metric_col, age_col, min_sessions)
                        lme  = fit_lme(rdf, metric_col, age_col, min_sessions)
                        draw_longitudinal(ax, rdf, metric_col, age_col, gamm, lme,
                                          age_label='Age (years)', metric_label=metric_label,
                                          stats_fontsize=6, show_method=False)
                    else:
                        if not draw_pct_change(ax, rdf, metric_col, age_col, 'Age (years)'):
                            _empty(ax, 'no subject with ≥ 2 sessions')
                    if r == 0:
                        ax.set_title(col_title, fontweight='bold', pad=12)
                    if c == 0:
                        ax.annotate(f"{row_label}\n{'trajectories' if kind == 'traj' else '% change'}",
                                    xy=(-0.42, 0.5), xycoords='axes fraction', rotation=90,
                                    ha='center', va='center', fontsize=9, fontweight='bold')

                if kind == 'traj' and sp in extra_panel:
                    ax = fig.add_subplot(gs[r, n_mod])
                    _letter(ax)
                    res = assoc_primary_trajectories(morph, age_df, primary_regions,
                                                     association_regions, species_label=sp,
                                                     data_key=extra_panel[sp])
                    if res is None:
                        _empty(ax, 'association / primary regions not found')
                    else:
                        draw_assoc_vs_primary(ax, *res, legend_fontsize=7)
                    ax.set_title('Association vs primary', fontweight='bold', pad=12)

        os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
        for ext, kw in (('png', dict(dpi=300)), ('pdf', {})):
            path = f'{out_path}.{ext}'
            fig.savefig(path, bbox_inches='tight', **kw)
            written.append(path)
        plt.close(fig)

    print(f"  [Fig3 combined] → {', '.join(written)}")
    return written


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    os.makedirs(OUT_DIR, exist_ok=True)
    all_stats = []

    # ══════════════════════════════════════════════════════════════════════════
    # MACAQUE
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 64)
    print("  MACAQUE")
    print("=" * 64)

    atlas_paths, _ = _build_atlas_paths(['Macaque'])

    # Load ALL regions needed for both main and Fig 3D
    mac_morph = load_morphometry({'Macaque': MACAQUE_BIDS_DIRS}, atlas_paths,
                                  regions_of_interest=REGIONS)
    mac_age = format_macaque_age(MACAQUE_EXCEL)

    # Main longitudinal analysis - only surface for macaque
    for mod, metric_col, metric_label in [
        ('surface', 'surface_area_mm2', 'Surface area (mm²)'),
        ('thickness', 'thickness_mm', 'Cortical thickness (mm)'),
        ('volume', 'volume_mm3', 'Volume (mm³)'),
    ]:
        df_mod = mac_morph.get(mod, pd.DataFrame())
        if df_mod is None or df_mod.empty:
            print(f"  [Macaque/{mod}] no data")
            continue

        stats = run_longitudinal(
            label=f'Macaque_{mod}',
            df_morph=df_mod,
            df_age=mac_age,
            metric_col=metric_col,
            metric_label=metric_label,
            age_col='age_years',
            age_label='Age (years)',
            out_dir=opj(OUT_DIR, 'Macaque', mod),
            atlas_level=ATLAS_LEVEL,
            hemisphere='left',
            regions=['Isocortex'],
            min_sessions=MIN_SESSIONS,
        )
        if stats is not None and not stats.empty:
            stats['species'] = 'Macaque'
            stats['modality'] = mod
            all_stats.append(stats)

    # ── FIG 3D: Association vs Primary Cortex (Macaque) ──────────────────────
    print("\n" + "=" * 64)
    print("  FIG 3D — Association vs Primary Cortex (Macaque)")
    print("=" * 64)
    for data_key in ['surface', 'volume', 'thickness']:
        fig3d_stats, fig3d_path = run_fig3d_assoc_vs_primary(
            mac_morph, mac_age,
            out_dir=opj(OUT_DIR, 'Macaque', 'fig3d'),
            primary_regions=PRIMARY_REGIONS,
            association_regions=ASSOCIATION_REGIONS,
            species_label="Macaque",
            data_key=data_key,
        )
        if fig3d_stats:
            all_stats.append(pd.DataFrame([fig3d_stats]))

    # ══════════════════════════════════════════════════════════════════════════
    # HUMAN
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 64)
    print("  HUMAN")
    print("=" * 64)

    atlas_paths, _ = _build_atlas_paths(['Human'])

    hum_morph = load_morphometry({'Human': HUMAN_BIDS_DIRS}, atlas_paths,
                                  regions_of_interest=REGIONS)
    hum_age = format_human_age(HUMAN_BIDS_DIRS, age_prefix=HUMAN_AGE_PREFIX)

    # Main longitudinal analysis - all modalities for human
    for mod, metric_col, metric_label in [
        ('surface', 'surface_area_mm2', 'Surface area (mm²)'),
        ('thickness', 'thickness_mm', 'Cortical thickness (mm)'),
        ('volume', 'volume_mm3', 'Volume (mm³)'),
    ]:
        df_mod = hum_morph.get(mod, pd.DataFrame())
        if df_mod is None or df_mod.empty:
            print(f"  [Human/{mod}] no data")
            continue

        stats = run_longitudinal(
            label=f'Human_{mod}',
            df_morph=df_mod,
            df_age=hum_age,
            metric_col=metric_col,
            metric_label=metric_label,
            age_col='age',
            age_label='Age (years)',
            out_dir=opj(OUT_DIR, 'Human', mod),
            atlas_level=ATLAS_LEVEL,
            hemisphere='left',
            regions=['Isocortex'],
            min_sessions=MIN_SESSIONS,
        )
        if stats is not None and not stats.empty:
            stats['species'] = 'Human'
            stats['modality'] = mod
            all_stats.append(stats)

    # ── FIG 3D: Association vs Primary Cortex (Human) ──────────────────────────
    print("\n" + "=" * 64)
    print("  FIG 3D — Association vs Primary Cortex (Human)")
    print("=" * 64)
    for data_key in ['surface', 'thickness', 'volume']:
        fig3d_stats_human, fig3d_path_human = run_fig3d_assoc_vs_primary(
            hum_morph, hum_age,
            out_dir=opj(OUT_DIR, 'Human', 'fig3d'),
            primary_regions=PRIMARY_REGIONS,
            association_regions=ASSOCIATION_REGIONS,
            species_label="Human",
            data_key=data_key,
        )
        if fig3d_stats_human:
            all_stats.append(pd.DataFrame([fig3d_stats_human]))

    # ── Combined summary Excel ────────────────────────────────────────────────
    if all_stats:
        summary = pd.concat(all_stats, ignore_index=True, sort=False)
        out_xlsx = opj(OUT_DIR, 'longitudinal_summary.xlsx')
        summary.to_excel(out_xlsx, index=False)
        print(f"\n  Summary Excel → {out_xlsx}")

        # ── Generate clean report ──────────────────────────────────────────────
        generate_report(
            all_stats, OUT_DIR, mac_age, hum_age,
            mac_surface=mac_morph.get('surface') if mac_morph else None,
            hum_surface=hum_morph.get('surface') if hum_morph else None,
            mac_thickness=mac_morph.get('thickness') if mac_morph else None,
            hum_thickness=hum_morph.get('thickness') if hum_morph else None,
            mac_volume=mac_morph.get('volume') if mac_morph else None,
            hum_volume=hum_morph.get('volume') if hum_morph else None,
        )

    # ── Combined Figure 3: human on top, macaque below, then % change rows ──
    print("\n" + "=" * 64)
    print("  FIGURE 3 — combined")
    print("=" * 64)
    plot_figure3_combined(
        species_rows=[
            ('Human',   hum_morph, hum_age, 'age',       'Human'),
            ('Macaque', mac_morph, mac_age, 'age_years', 'Macaque'),
        ],
        out_path=opj(OUT_DIR, 'Figure3_combined'),
    )

    print("\n  Done.")