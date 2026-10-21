"""
run_SBA_by_category.py
──────────────────────
Seed-based analysis (SBA) rendering script for EDNiX species,
run separately for each classification category (Specific, Unspecific, No, Spurious).
Threshold is based on the total number of BOLD runs (files), not subjects.
"""

import os
import sys
import glob
import argparse
import random
import shutil
from pathlib import Path
import pandas as pd
import numpy as np

# ── EDNiX environment ─────────────────────────────────────────────────────────
MAIN_PATH = '/home/cgarin/PycharmProjects/EDNiX/'
ATLAS_LIB = os.path.join(MAIN_PATH, 'Atlases_library')
sys.path.insert(0, MAIN_PATH)
sys.path.insert(0, os.path.join(MAIN_PATH, 'code'))

# ── Imports ───────────────────────────────────────────────────────────────────
from Tools import Load_subject_with_BIDS, load_bids
from Exemples.Study_EDNiX.analysis.ednix_bids_tools import (
    extract_bold_paths, parse_label_file, get_atlas_label_path,
)
from Tools import niftotoWBsurface
from Statistics.Group_fMRI._Group_anal_3dTtest import _3dttest_EDNiX
from Statistics.Group_fMRI import _Group_anal_3dLMEr_SBA as _lmer_sba
from Statistics.Group_fMRI import _Group_anal_3dLMEr_Mirror as _lmer_mirror

# ── Singularity / wb_command ─────────────────────────────────────────────────
SING_WB = (
    'singularity exec --cleanenv '
    f'--bind {ATLAS_LIB}/,/scratch2/,/scratch/ '
    f'/home/cgarin/PycharmProjects/EDNiX/Singularity_library/Singularity/workbench_2.1.0.sif'
)

NORMALIZE_PERCENTILE = 95

# ── Output directory ─────────────────────────────────────────────────────────
OUTPUT_BASE = '/scratch2/EDNiX/results/EDNiX_Figures/Figure_5_SBA'

# ── Dataset definitions ───────────────────────────────────────────────────────
SPECIES_BIDS_DICT = {
    'Rat':        '/scratch2/EDNiX/Rat/BIDS_Grandjean',
    'Mouse':      '/scratch2/EDNiX/Mouse/BIDS_Grandjean2',
    'Dog':        '/scratch2/EDNiX/Dog/BIDS_Boch_K9',
    'Marmoset':   '/scratch2/EDNiX/Marmoset/BIDS_Tian',
    'Mouselemur': '/scratch2/EDNiX/Mouselemur/BIDS_Garin',
}
SPECIES_MULTI_BIDS = {
    'Macaque': [
        '/scratch2/EDNiX/Macaque/BIDS_Tremblay',
        '/scratch2/EDNiX/Macaque/BIDS_BenHamed',
        '/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin',
    ],
    'Human': [
        '/scratch2/EDNiX/Human/BIDS_Merida',
        '/scratch2/EDNiX/Human/BIDS_Park',
        '/scratch2/EDNiX/Human/BIDS_Castrillon',
    ],
    'Marmoset': [
    '/scratch2/EDNiX/Marmoset/BIDS_Tian',
    '/scratch2/EDNiX/Marmoset/BIDS_Hori',
    ],
}

# ── Per-species configuration overrides ──────────────────────────────────────
SPECIES_OVERRIDES = {
    'Mouse': {
        'species_fragment': 'atlas/mammals/rodents/Mouse',
        'anat_suffix':      'T2w',
        'alpha':            0.1,
        'method_mask_func':  'onlyprovidedmask',
        'cut_coords':       [-6, -5, -4, -2, -1, 1, 3, 4, 5, 6],
        'projection':       'ribbon_dil',
    },
    'Rat': {
        'species_fragment': 'atlas/mammals/rodents/Rat',
        'anat_suffix':      'T2w',
        'alpha':            0.1,
        'method_mask_func':  'onlyprovidedmask',
        'cut_coords':       [-6, -5, -4, -2, -1, 1, 3, 4, 5, 6],
        'projection':       'ribbon_dil',
    },
    'Dog': {
        'species_fragment': 'atlas/mammals/carnivorans/Dog',
        'alpha':            0.1,
        'cut_coords':       [-6, -5, -4, -2, -1, 1, 3, 4, 5, 6],
        'method_mask_func':  'onlyprovidedmask',
        'projection':       'ribbon_dil',
        'delete_surf_file':  False,
    },
    'Marmoset': {
        'species_fragment': 'atlas/mammals/primates/Marmoset',
        'alpha':            0.1,
        'anat_suffix':      'T2w',
        'method_mask_func':  'onlyprovidedmask',
        'cut_coords':       [-6, -5, -4, -2, -1, 1, 3, 4, 5, 6],
        'projection':       'ribbon_dil',
    },
    'Mouselemur': {
        'species_fragment': 'atlas/mammals/primates/Mouselemur',
        'alpha':            0.1,
        'anat_suffix':      'T2w',
        'method_mask_func':  'onlyprovidedmask',
        'cut_coords':       [-6, -5, -4, -2, -1, 1, 3, 4, 5, 6],
        'projection':       'ribbon_dil',
    },
    'Macaque': {
        'species_fragment':  'atlas/mammals/primates/Macaque',
        'method_mask_func':  'onlyprovidedmask',
        'projection':       'ribbon_dil',
        'alpha': 0.1,
        'delete_surf_file':  False,
    },
    'Human': {
        'species_fragment':  'atlas/mammals/primates/Human',
        'method_mask_func':  'onlyprovidedmask',
        'projection':       'ribbon_dil',
        'alpha': 0.1,
        'delete_surf_file':  False,
    },
}

# ── Regions of interest ──────────────────────────────────────────────────────
# ── Regions of interest ──────────────────────────────────────────────────────
REGIONS_OF_INTEREST = [
    'BA_26_29_30__retrosplenial',
    'Posterior_medial_cortex__PMC',
]

CATEGORIES = ["Specific", "Unspecific", "No", "Spurious", "All"]

# =============================================================================
# CLASSIFICATION FILTERING
# =============================================================================

def load_classification_subjects(excel_path, species=None, bids_dir=None,
                                 categories=None, min_n=4):
    """Load classification Excel and extract subjects for each category."""
    if not os.path.exists(excel_path):
        print(f"  [ERROR] Classification file not found: {excel_path}")
        return None, None

    # Build mapping from BIDS directory name to full path
    bids_name_to_path = {}
    for sp, d in SPECIES_BIDS_DICT.items():
        bids_name_to_path[Path(d).name] = d
    for sp, dirs in SPECIES_MULTI_BIDS.items():
        for d in dirs:
            bids_name_to_path[Path(d).name] = d

    try:
        df_subject = pd.read_excel(excel_path, sheet_name='Subject_Summary')
    except:
        df_subject = pd.read_excel(excel_path, sheet_name='Run_Classification')
        df_subject = df_subject.groupby(['species', 'bids_dir', 'subject']).agg({
            'fp_category': lambda x: list(x)
        }).reset_index()

        def get_primary_category(cat_list):
            counts = {c: 0 for c in CATEGORIES if c != 'All'}
            for c in cat_list:
                if c in counts:
                    counts[c] += 1
            return max(counts, key=counts.get) if counts else 'No'

        df_subject['primary_category'] = df_subject['fp_category'].apply(get_primary_category)

    if species is not None:
        df_subject = df_subject[df_subject['species'] == species]
    if bids_dir is not None:
        df_subject = df_subject[df_subject['bids_dir'] == bids_dir]

    if df_subject.empty:
        return None, None

    cat_col = None
    for col in ['primary_category', 'fp_category']:
        if col in df_subject.columns:
            cat_col = col
            break

    if cat_col is None:
        print(f"  [ERROR] No category column found in {excel_path}")
        return None, None

    category_subjects = {}
    categories_to_use = categories if categories is not None else CATEGORIES
    real_categories = [c for c in CATEGORIES if c != 'All']

    # If 'All' was requested, we still need every real category loaded to build it,
    # even if the caller only asked for 'All' itself.
    need_all = 'All' in categories_to_use
    categories_to_load = real_categories if need_all else [c for c in categories_to_use if c != 'All']

    for cat in categories_to_load:
        cat_df = df_subject[df_subject[cat_col] == cat]
        if cat_df.empty:
            category_subjects[cat] = []
            continue

        subjects = []
        for _, row in cat_df.iterrows():
            sp = str(row['species'])
            bd_name = str(row['bids_dir'])
            subj = str(row['subject'])
            ses = str(row.get('session', '1'))

            if bd_name in bids_name_to_path:
                bd_full = bids_name_to_path[bd_name]
            else:
                bd_full = bd_name
                print(f"  [WARN] No full path found for BIDS directory: {bd_name}")

            subjects.append((sp, bd_full, subj, ses))

        category_subjects[cat] = subjects

    # Now handle the 'All' category - combine all subjects from all real categories
    if need_all:
        all_subjects = []
        for cat in real_categories:
            all_subjects.extend(category_subjects.get(cat, []))

        # Remove duplicates
        seen = set()
        unique_subjects = []
        for subject in all_subjects:
            key = (subject[1], subject[2], subject[3])  # (bids_dir, subject, session)
            if key not in seen:
                seen.add(key)
                unique_subjects.append(subject)

        category_subjects['All'] = unique_subjects
        print(f"  [INFO] Created 'All' category with {len(unique_subjects)} subjects")

    # Trim back down to only what the caller actually asked for
    category_subjects = {c: category_subjects.get(c, []) for c in categories_to_use}

    return category_subjects

    # Now handle the 'All' category - combine all subjects from all categories
    if 'All' in categories_to_use:
        all_subjects = []
        for cat in CATEGORIES:
            if cat != 'All' and cat in category_subjects:
                all_subjects.extend(category_subjects[cat])

        # Remove duplicates
        seen = set()
        unique_subjects = []
        for subject in all_subjects:
            key = (subject[1], subject[2], subject[3])  # (bids_dir, subject, session)
            if key not in seen:
                seen.add(key)
                unique_subjects.append(subject)

        category_subjects['All'] = unique_subjects
        print(f"  [INFO] Created 'All' category with {len(unique_subjects)} subjects")

    return category_subjects


def subsample_categories(category_subjects, target_n=None, random_seed=42):
    """Subsample categories to have the same number of subjects. Skip 'All' category."""
    random.seed(random_seed)
    np.random.seed(random_seed)

    n_per_category = {}
    for cat, subjects in category_subjects.items():
        n_per_category[cat] = len(subjects)

    # Only consider categories with >= 4 subjects, EXCLUDING 'All'
    valid_cats = {cat: n for cat, n in n_per_category.items() if n >= 4 and cat != 'All'}

    if not valid_cats:
        return {}, 0

    if target_n is None:
        target_n = min(valid_cats.values())
        print(f"  [subsample] Target N = {target_n} (smallest category with >=4)")

    subsampled = {}
    for cat, subjects in category_subjects.items():
        if cat == 'All':
            subsampled[cat] = []
            print(f"  [subsample] {cat}: SKIPPED (not used for homogenization)")
            continue

        n = len(subjects)
        if n >= 4:
            if n <= target_n:
                subsampled[cat] = subjects
                print(f"  [subsample] {cat}: {n} → {n} (kept all)")
            else:
                sampled = random.sample(subjects, target_n)
                subsampled[cat] = sampled
                print(f"  [subsample] {cat}: {n} → {target_n} (sampled)")
        else:
            subsampled[cat] = []
            print(f"  [subsample] {cat}: {n} → 0 (skipped, n<4)")

    return subsampled, target_n

# =============================================================================
# SBA PROCESSING
# =============================================================================

def build_species_cfg(species: str) -> dict:
    defaults = dict(
        species_fragment  = f'atlas/mammals/{species.lower()}/{species}',
        anat_suffix       = 'T1w',
        alpha             = 0.05,
        lower_cutoff      = 0.90,
        upper_cutoff      = 0.95,
        cut_coords        = 5,
        method_mask_func  = None,
        projection        = None,
        delete_surf_file  = None,
    )
    defaults.update(SPECIES_OVERRIDES.get(species, {}))
    return defaults


def find_template_low(bids_dirs: list[str], species: str) -> str | None:
    pattern = '**/templates/EDNiX/preprocessing/BASE_SS_fMRI.nii.gz'
    for bids_dir in bids_dirs:
        matches = sorted(glob.glob(os.path.join(bids_dir, pattern), recursive=True))
        if matches:
            return matches[0]
    return None


def collect_bold_and_means(bids_dirs: list[str], list_to_keep=None) -> tuple[dict, list[str]]:
    """Collect BOLD paths and mean images for subjects in list_to_keep."""
    keep_set = None
    if list_to_keep is not None:
        keep_set = set()
        for item in list_to_keep:
            if len(item) >= 2:
                keep_set.add((item[0], item[1]))
            elif len(item) == 1:
                keep_set.add((item[0], '1'))

    all_bold = {'bold_path': [], 'subject': [], 'session': [], 'run': []}
    mean_imgs = []
    seen = set()

    for bids_dir in bids_dirs:
        print(f"  [DEBUG] Searching for BOLD files in: {bids_dir}")

        if not os.path.exists(bids_dir):
            print(f"  [ERROR] BIDS directory does not exist: {bids_dir}")
            continue

        for sub, ses in keep_set:
            # Try multiple session formats
            ses_formats = [ses]
            if ses.isdigit():
                ses_formats.extend([f"{int(ses):02d}", f"{int(ses):03d}"])
            ses_formats = list(set(ses_formats))

            # Find all existing sessions for this subject
            subj_pattern = os.path.join(bids_dir, f'sub-{sub}', 'ses-*')
            existing_sessions = glob.glob(subj_pattern)
            if existing_sessions:
                for ses_path in existing_sessions:
                    ses_name = os.path.basename(ses_path).replace('ses-', '')
                    if ses_name not in ses_formats:
                        ses_formats.append(ses_name)
                print(f"  [DEBUG] Found existing sessions for {sub}: {ses_formats}")

            found_mean = False
            found_bold = False
            best_ses_fmt = None

            for ses_fmt in ses_formats:
                ses_dir = os.path.join(bids_dir, f'sub-{sub}', f'ses-{ses_fmt}')
                if not os.path.exists(ses_dir):
                    continue

                # Search for mean images
                postpro_dirs = [
                    os.path.join(ses_dir, 'func', 'templates', 'EDNiX', 'postprocessed_rs'),
                    os.path.join(ses_dir, 'func', 'templates', 'EDNiX', 'preprocessed_rs'),
                ]

                mean_patterns = [
                    '*fMRI_Mean_Image_SS.nii.gz',
                    '*desc-fMRI_Mean_Image_SS.nii.gz',
                    'all_runs_space-template-func_desc-fMRI_Mean_Image_SS.nii.gz',
                ]

                for postpro_dir in postpro_dirs:
                    for pattern in mean_patterns:
                        glob_pattern = os.path.join(postpro_dir, pattern)
                        matches = glob.glob(glob_pattern)
                        if matches:
                            key = (bids_dir, sub, ses_fmt)
                            if key not in seen:
                                seen.add(key)
                                mean_imgs.append(matches[0])
                                print(f"  [mean] Found mean image for {sub}/{ses_fmt}: {os.path.basename(matches[0])}")
                                found_mean = True
                                best_ses_fmt = ses_fmt
                                break
                    if found_mean:
                        break
                if found_mean:
                    break

            if not found_mean:
                print(f"  [WARN] No mean image found for {sub}/{ses} in {bids_dir}")
                continue

            ses_dir = os.path.join(bids_dir, f'sub-{sub}', f'ses-{best_ses_fmt}')
            bold_search_dirs = [
                os.path.join(ses_dir, 'func', 'templates', 'EDNiX', 'postprocessed_rs'),
                os.path.join(ses_dir, 'func', 'templates', 'EDNiX', 'preprocessed_rs'),
                os.path.join(ses_dir, 'func'),
            ]

            bold_patterns = [
                '*_task-rest_*_desc-fMRI_residual.nii.gz',
                '*_task-rest_*_bold_space-template_desc-fMRI_residual.nii.gz',
                '*_desc-fMRI_residual.nii.gz',
                '*_desc-preproc_bold.nii.gz',
                '*_bold.nii.gz',
                '*bold*.nii.gz',
                '*residual.nii.gz',
            ]

            for search_dir in bold_search_dirs:
                if not os.path.exists(search_dir):
                    continue
                for pattern in bold_patterns:
                    glob_pattern = os.path.join(search_dir, pattern)
                    matches = glob.glob(glob_pattern)
                    if matches:
                        for bold_path in matches:
                            run = '01'
                            filename = os.path.basename(bold_path)
                            if 'run-RL-' in filename or 'run-RL_' in filename:
                                if 'run-RL-' in filename:
                                    run_parts = filename.split('run-RL-')
                                    if len(run_parts) > 1:
                                        run = 'RL-' + run_parts[1].split('_')[0]
                                elif 'run-RL_' in filename:
                                    run_parts = filename.split('run-RL_')
                                    if len(run_parts) > 1:
                                        run = 'RL-' + run_parts[1].split('_')[0]
                            elif 'run-LR-' in filename or 'run-LR_' in filename:
                                if 'run-LR-' in filename:
                                    run_parts = filename.split('run-LR-')
                                    if len(run_parts) > 1:
                                        run = 'LR-' + run_parts[1].split('_')[0]
                                elif 'run-LR_' in filename:
                                    run_parts = filename.split('run-LR_')
                                    if len(run_parts) > 1:
                                        run = 'LR-' + run_parts[1].split('_')[0]
                            elif 'run-' in filename:
                                run = filename.split('run-')[1].split('_')[0]
                            elif 'run-' in bold_path:
                                parts = bold_path.split(os.sep)
                                for part in parts:
                                    if 'run-' in part:
                                        run = part.split('run-')[1].split('_')[0]
                                        break

                            if bold_path not in all_bold['bold_path']:
                                all_bold['bold_path'].append(bold_path)
                                all_bold['subject'].append(sub)
                                all_bold['session'].append(best_ses_fmt)
                                all_bold['run'].append(run)
                                print(f"  [bold] Found BOLD file for {sub}/{best_ses_fmt}/run-{run}: {os.path.basename(bold_path)}")
                                found_bold = True
                        if found_bold:
                            break
                if found_bold:
                    break

            if not found_bold:
                func_dir = os.path.join(ses_dir, 'func')
                if os.path.exists(func_dir):
                    for root, dirs, files in os.walk(func_dir):
                        for filename in files:
                            if not filename.endswith('.nii.gz'):
                                continue
                            if 'mean' in filename.lower() or 'mask' in filename.lower():
                                continue
                            if 'bold' in filename.lower() or 'BOLD' in filename or 'residual' in filename:
                                nii_file = os.path.join(root, filename)
                                run = '01'
                                if 'run-RL-' in filename:
                                    run_parts = filename.split('run-RL-')
                                    if len(run_parts) > 1:
                                        run = 'RL-' + run_parts[1].split('_')[0]
                                elif 'run-LR-' in filename:
                                    run_parts = filename.split('run-LR-')
                                    if len(run_parts) > 1:
                                        run = 'LR-' + run_parts[1].split('_')[0]
                                elif 'run-' in filename:
                                    run = filename.split('run-')[1].split('_')[0]

                                if nii_file not in all_bold['bold_path']:
                                    all_bold['bold_path'].append(nii_file)
                                    all_bold['subject'].append(sub)
                                    all_bold['session'].append(best_ses_fmt)
                                    all_bold['run'].append(run)
                                    print(f"  [bold] Found BOLD file in func dir for {sub}/{best_ses_fmt}: {filename}")
                                    found_bold = True
                                    break
                        if found_bold:
                            break

            if not found_bold:
                print(f"  [WARN] No BOLD file found for {sub}/{best_ses_fmt} in {bids_dir}")

    print(f"  [DEBUG] Found {len(all_bold['bold_path'])} BOLD files total")
    print(f"  [DEBUG] Found {len(mean_imgs)} mean images total")

    return all_bold, mean_imgs


def process_species_category(species: str, bids_dirs: list[str],
                             category: str, list_to_keep: list,
                             bids_label: str,
                             output_suffix: str = "",
                             skip_analyses: list = None,
                             results_base_dir: str = None,
                             min_runs: int = 15,
                             exclude_lateralized: bool = True) -> bool:   # NEW) -> bool:
    """
    Run SBA for one species and one category.
    Returns True if analysis was attempted, False if skipped due to insufficient runs.
    """
    if skip_analyses is None:
        skip_analyses = []

    if not list_to_keep:
        print(f"\n  [SKIP] {species} - {category}: No subjects to process")
        return False

    # ─── Group subjects by BIDS directory ──────────────────────────────────
    subjects_by_bids = {}
    for bids_dir in bids_dirs:
        subjects_by_bids[bids_dir] = []

    for subject_info in list_to_keep:
        bids_dir_from_excel = subject_info[1]
        matched = False
        for bids_dir in bids_dirs:
            if Path(bids_dir_from_excel).name in bids_dir or Path(bids_dir).name in bids_dir_from_excel:
                subjects_by_bids[bids_dir].append(subject_info)
                matched = True
                break
        if not matched:
            for bids_dir in bids_dirs:
                subj = subject_info[2]
                subj_dir = os.path.join(bids_dir, f'sub-{subj}')
                if os.path.exists(subj_dir):
                    subjects_by_bids[bids_dir].append(subject_info)
                    matched = True
                    break
        if not matched:
            print(f"  [WARN] Could not match subject {subject_info[2]} to any BIDS directory")

    subjects_by_bids = {bids_dir: subjects for bids_dir, subjects in subjects_by_bids.items() if subjects}
    if not subjects_by_bids:
        print(f"  [SKIP] {species} - {category}: No subjects matched to any BIDS directory")
        return False

    print(f"\n  Subjects by BIDS directory:")
    for bids_dir, subjects in subjects_by_bids.items():
        print(f"    {Path(bids_dir).name}: {len(subjects)} subjects")

    # ─── Process each BIDS directory separately ────────────────────────────
    for bids_dir, subjects in subjects_by_bids.items():
        subject_sessions = [(subj, ses) for _, _, subj, ses in subjects]
        if not subject_sessions:
            continue

        # ─── Set up output directory ──────────────────────────────────────
        if results_base_dir is None:
            results_base_dir = OUTPUT_BASE

        bids_name = Path(bids_dir).name
        output_dir = os.path.join(
            results_base_dir,
            species,
            bids_name,
            output_suffix if output_suffix else "Full",
            category
        )
        os.makedirs(output_dir, exist_ok=True)

        print(f"\n{'#' * 70}")
        print(f"  SPECIES: {species}  |  CATEGORY: {category}  ({len(subject_sessions)} subjects)")
        print(f"  BIDS: {bids_name}")
        print(f"  OUTPUT: {output_dir}")
        print(f"{'#' * 70}")

        cfg = build_species_cfg(species)

        # ── Derived atlas paths ────────────────────────────────────────────
        atlas_ednix = os.path.join(ATLAS_LIB, cfg['species_fragment'], 'EDNiX')
        base_ss = os.path.join(atlas_ednix, 'volumes',
                               f'{species}_space-acpc_desc-SS_T1w.nii.gz')
        mask_func = os.path.join(atlas_ednix, 'volumes', 'masks',
                                 f'{species}_desc-Cerebral_Gray_mask.nii.gz')
        scene_template = os.path.join(atlas_ednix, 'surfaces', 'Native_resol',
                                      f'{species}_thickness.scene')
        template_low = find_template_low([bids_dir], species)

        # ── Collect BOLD paths and mean images ────────────────────────────
        bold_paths, mean_imgs = collect_bold_and_means([bids_dir], subject_sessions)
        n_runs = len(bold_paths['bold_path'])
        print(f"  {n_runs} BOLD files found.")
        print(f"  {len(mean_imgs)} mean images collected.")

        # ─── Check minimum runs ────────────────────────────────────────────
        if n_runs < min_runs:
            print(f"  [WARN] Only {n_runs} BOLD runs found, below minimum {min_runs}")
            print(f"  [WARN] Skipping SBA for {bids_name} - {category}")
            return False

        if not mean_imgs:
            print(f"  [WARN] No mean images found for {bids_name} - {category}")
            print(f"  [WARN] Skipping SBA for this BIDS directory")
            return False

        # ── Atlas label file ───────────────────────────────────────────────
        label_path = get_atlas_label_path(ATLAS_LIB, cfg['species_fragment'], 'EDNIxCSC')
        label_df = parse_label_file(label_path)
        print(f"  {len(label_df)} regions in label file.")

        # ── Surface rendering config ───────────────────────────────────────
        surface_cfg = dict(
            species_name=species,
            atlas_lib_path=ATLAS_LIB,
            species_fragment=cfg['species_fragment'],
            sing_wb=SING_WB,
            surface_type='inflated',
            palette_name='RY-BC-BL',
            pos_max=90,
            neg_max=10,
            colorbar_mode='abs_pct',
            scene_template=scene_template,
            scene_label='New Scene 1',
            display_zero=False,  # ← cache les voxels à valeur nulle
        )
        if cfg.get('projection') is not None:
            surface_cfg['projection'] = cfg['projection']
        if cfg.get('delete_surf_file') is not None:
            surface_cfg['delete_surf_file'] = cfg['delete_surf_file']

        # ── Shared analysis kwargs ────────────────────────────────────────
        shared = dict(
            templatehigh=base_ss,
            templatelow=template_low,
            oversample_map=False,
            mask_func=mask_func,
            label_df=label_df,
            MAIN_PATH=MAIN_PATH,
            bold_paths=bold_paths,
            mean_imgs=mean_imgs,
            regions_of_interest=REGIONS_OF_INTEREST,
            alpha=cfg['alpha'],
            lower_cutoff=cfg['lower_cutoff'],
            upper_cutoff=cfg['upper_cutoff'],
            cut_coords=cfg['cut_coords'],
            normalize_target=NORMALIZE_PERCENTILE,
        )
        if cfg.get('method_mask_func') is not None:
            shared['method_mask_func'] = cfg['method_mask_func']

        tag = Path(bids_dir).name
        category_tag = f"{category}{output_suffix}"

        # ── ANALYSIS 1: 3dTtest++ ─────────────────────────────────────────
        if 'analysis1' not in skip_analyses:
            print(f"\n  {'=' * 66}")
            print(f"  ANALYSIS 1  —  3dTtest++ seed-based  [{tag} - {category_tag}]")
            print(f"  {'=' * 66}")
            try:
                shared_with_bids = dict(shared, bids_dir=bids_dir)
                _3dttest_EDNiX(**shared_with_bids)

                source_dir = os.path.join(bids_dir, 'Results', 'Grp_SBA_3dTTEST')
                target_dir = os.path.join(output_dir, f'Grp_SBA_3dTTEST_{category_tag}')
                if os.path.exists(source_dir):
                    if os.path.exists(target_dir):
                        shutil.rmtree(target_dir)
                    shutil.copytree(source_dir, target_dir)
                    stat_suffixes = ['ttest-stat_fisher_zmap', 'ttest-stat_fisher_cmap']
                    if NORMALIZE_PERCENTILE is not None:
                        stat_suffixes.append('correlation_norm')
                    niftotoWBsurface.render_stat_maps(
                        results_dir=output_dir,
                        results_subdir=f'Grp_SBA_3dTTEST_{category_tag}',
                        stat_suffixes=stat_suffixes,
                        regions_of_interest=REGIONS_OF_INTEREST,
                        surface_cfg=surface_cfg,
                    )
            except Exception as e:
                print(f"  [WARN] Analysis 1 failed: {e}")
                import traceback
                traceback.print_exc()

        # ── ANALYSIS 2: 3dLMEr SBA (no mirroring) ────────────────────────
        if 'analysis2' not in skip_analyses:
            print(f"\n  {'=' * 66}")
            print(f"  ANALYSIS 2  —  3dLMEr SBA (no mirroring)  [{tag} - {category_tag}]")
            print(f"  {'=' * 66}")
            try:
                shared_with_bids = dict(shared, bids_dir=bids_dir)
                _lmer_sba._3dLMEr_EDNiX(
                    **shared_with_bids,
                    model="'(1|Subj)'",
                    gltCode="-gltCode groupeffect 'Subj : 1'",
                    stat_img=['groupeffect'],
                )

                source_dir = os.path.join(bids_dir, 'Results', 'Grp_SBA_3dLME_network')
                target_dir = os.path.join(output_dir, f'Grp_SBA_3dLME_network_{category_tag}')
                if os.path.exists(source_dir):
                    if os.path.exists(target_dir):
                        shutil.rmtree(target_dir)
                    shutil.copytree(source_dir, target_dir)
                    niftotoWBsurface.render_stat_maps(
                        results_dir=output_dir,
                        results_subdir=f'Grp_SBA_3dLME_network_{category_tag}',
                        stat_suffixes=['groupeffect_zstat', 'groupeffect_coef'],
                        regions_of_interest=REGIONS_OF_INTEREST,
                        surface_cfg=surface_cfg,
                    )
            except Exception as e:
                print(f"  [WARN] Analysis 2 failed: {e}")
                import traceback
                traceback.print_exc()

        # ── ANALYSIS 3: 3dLMEr Mirror (hemisphere) ──────────────────────
        if 'analysis3' not in skip_analyses:
            print(f"\n  {'=' * 66}")
            print(f"  ANALYSIS 3  —  3dLMEr Mirror  [{tag} - {category_tag}]")
            print(f"  {'=' * 66}")
            try:
                shared_with_bids = dict(shared, bids_dir=bids_dir)
                _lmer_mirror._3dLMEr_EDNiX(
                    **shared_with_bids,
                    model='(1|Subj)*Hemisphere',
                    glt_spec=[
                        (1, 'Overall_Effect', 'Hemisphere : 0.5*L +0.5*R'),
                        (2, 'Hemisphere_Effect', 'Hemisphere : 1*L -1*R'),
                    ],
                    contrast_names=['Overall_Effect', 'Hemisphere_Effect'],
                    midline_x=0.0,
                    visualize='threshold',
                    percent=5,
                )

                source_dir = os.path.join(bids_dir, 'Results', 'Grp_SBA_3dLME_network')
                target_dir = os.path.join(output_dir, f'Grp_SBA_3dLME_network_{category_tag}')
                if os.path.exists(source_dir):
                    if os.path.exists(target_dir):
                        shutil.rmtree(target_dir)
                    shutil.copytree(source_dir, target_dir)
                    stat_suffixes = ['Overall_Effect_zstat', 'Overall_Effect_coef']
                    if NORMALIZE_PERCENTILE is not None:
                        stat_suffixes += ['Overall_Effect_zstat_norm', 'Overall_Effect_coef_norm']
                    niftotoWBsurface.render_stat_maps(
                        results_dir=output_dir,
                        results_subdir=f'Grp_SBA_3dLME_network_{category_tag}',
                        stat_suffixes=stat_suffixes,
                        regions_of_interest=['BA_26_29_30__retrosplenial'],
                        surface_cfg=surface_cfg,
                    )
            except Exception as e:
                print(f"  [WARN] Analysis 3 failed: {e}")
                import traceback
                traceback.print_exc()

    return True

def main() -> None:
    parser = argparse.ArgumentParser(
        description='Run SBA rendering for EDNiX species — "All" category only.'
    )
    parser.add_argument(
        '--species', nargs='+', default=None,
        help='Subset of species to process (default: all).',
    )
    parser.add_argument(
        '--classification', type=str,
        default='/scratch2/EDNiX/results/EDNiX_Figures/Figure_4_rsfMRI_QC/lvl2/all/classification_results.xlsx',
        help='Path to classification_results.xlsx',
    )
    parser.add_argument(
        '--min-n', type=int, default=4,
        help='Minimum number of BOLD runs required to process (default: 4).',
    )
    parser.add_argument(
        '--skip-analyses', nargs='+', default=None,
        help='Skip specific analyses (e.g., analysis2 analysis3 for 3dttest++ only).',
    )
    parser.add_argument(
        '--output', type=str, default=OUTPUT_BASE,
        help='Base output directory for results.',
    )
    parser.add_argument(
        '--exclude-lateralized', action='store_true', default=True,
        help='Only run bilateral seeds (drop any region prefixed L_/R_) when a ROI substring matches multiple label rows.',
    )
    args, _ = parser.parse_known_args()

    # Build full species→[bids_dirs] map
    all_species: dict[str, list[str]] = {}
    for sp, d in SPECIES_BIDS_DICT.items():
        all_species[sp] = [d]
    for sp, dirs in SPECIES_MULTI_BIDS.items():
        all_species[sp] = dirs

    if args.species:
        unknown = set(args.species) - set(all_species)
        if unknown:
            print(f"[ERROR] Unknown species: {unknown}")
            sys.exit(1)
        to_run = {sp: all_species[sp] for sp in args.species}
    else:
        to_run = all_species

    print(f"\n{'=' * 70}")
    print(f"  EDNiX SBA — 'All' category only (threshold based on BOLD runs)")
    print(f"  Species: {', '.join(to_run)}")
    print(f"  Minimum runs: {args.min_n}")
    print(f"  Output: {args.output}")
    print(f"{'=' * 70}\n")

    # ── Load classification — only ask for 'All' ─────────────────────────
    # load_classification_subjects still needs to internally load every
    # sub-category to BUILD the 'All' union (that's its own logic), but we
    # only ever iterate/process the 'All' result below — Specific/
    # Unspecific/No/Spurious are never processed individually.
    #
    # ⚠ On a genuine failure (missing file, no matching rows), this
    # function returns the TUPLE (None, None), not a dict — .get() on
    # that raises AttributeError, so the type is checked before calling it.
    all_category_subjects = {}
    for species, bids_dirs in to_run.items():
        print(f"\n{'=' * 70}")
        print(f"  Loading classification for {species}")
        print(f"{'=' * 70}")

        category_subjects = load_classification_subjects(
            args.classification, species=species, bids_dir=None, categories=['All'],
        )
        if not isinstance(category_subjects, dict) or not category_subjects.get('All'):
            print(f"  [WARN] No 'All'-category subjects for {species}")
            continue

        all_category_subjects[species] = category_subjects
        print(f"  {species}: {len(category_subjects['All'])} subjects in 'All'")

    # ── Process each species — 'All' category only ────────────────────────
    print(f"\n{'=' * 70}")
    print(f"  Starting SBA processing (min runs = {args.min_n})")
    print(f"{'=' * 70}\n")

    all_processed, all_skipped = [], []
    for species, bids_dirs in to_run.items():
        if species not in all_category_subjects:
            continue
        subjects = all_category_subjects[species]['All']
        if not subjects:
            all_skipped.append((species, 0, "no subjects"))
            continue
        try:
            success = process_species_category(
                species, bids_dirs, 'All', subjects,
                bids_label=Path(bids_dirs[0]).name,
                output_suffix="Full",
                skip_analyses=args.skip_analyses if args.skip_analyses else ['analysis3', 'analysis2'],
                results_base_dir=args.output,
                min_runs=args.min_n,
                exclude_lateralized=args.exclude_lateralized,
            )
            if success:
                all_processed.append((species, len(subjects), "full"))
            else:
                all_skipped.append((species, len(subjects), "insufficient runs"))
        except Exception as exc:
            print(f"\n[ERROR] {species} - All failed: {exc}")
            import traceback
            traceback.print_exc()
            all_skipped.append((species, len(subjects), "error"))

    # ── Final summary ─────────────────────────────────────────────────────
    print(f"\n{'=' * 70}")
    print(f"  SBA COMPLETE")
    print(f"{'=' * 70}")
    print(f"\n  Processed:")
    for species, n_subj, mode in all_processed:
        print(f"    ✅ {species}: {n_subj} subjects ({mode})")
    if all_skipped:
        print(f"\n  Skipped:")
        for species, n_subj, reason in all_skipped:
            print(f"    ❌ {species}: {n_subj} subjects ({reason})")
    print(f"\n  Results saved to: {args.output}")
    print(f"\n{'#' * 70}\n  Done.\n{'#' * 70}\n")


if __name__ == '__main__':
    main()


'''
# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description='Run SBA rendering for EDNiX species by classification category (threshold based on number of BOLD runs).'
    )
    parser.add_argument(
        '--species', nargs='+', default=None,
        help='Subset of species to process (default: all).',
    )
    parser.add_argument(
        '--classification', type=str,
        default='/scratch2/EDNiX/results/EDNiX_Figures/Figure_4_rsfMRI_QC/lvl2/all/classification_results.xlsx',
        help='Path to classification_results.xlsx',
    )
    parser.add_argument(
        '--categories', nargs='+', default=None,
        help='Categories to process (default: all).',
    )
    parser.add_argument(
        '--min-n', type=int, default=4,
        help='Minimum number of BOLD runs required to process a category (default: 4).',
    )
    parser.add_argument(
        '--homogenize', action='store_true',
        help='Subsample all categories to the same N (smallest category with >=4 subjects).',
    )
    parser.add_argument(
        '--skip-analyses', nargs='+', default=None,
        help='Skip specific analyses (e.g., analysis1 analysis2).',
    )
    parser.add_argument(
        '--output', type=str, default=OUTPUT_BASE,
        help='Base output directory for results.',
    )
    parser.add_argument(
        '--exclude-lateralized', action='store_true', default=True,
        help='Only run bilateral seeds (drop any region prefixed L_/R_) when a ROI substring matches multiple label rows.',
    )
    args, _ = parser.parse_known_args()

    # Build full species→[bids_dirs] map
    all_species: dict[str, list[str]] = {}
    for sp, d in SPECIES_BIDS_DICT.items():
        all_species[sp] = [d]
    for sp, dirs in SPECIES_MULTI_BIDS.items():
        all_species[sp] = dirs

    if args.species:
        unknown = set(args.species) - set(all_species)
        if unknown:
            print(f"[ERROR] Unknown species: {unknown}")
            sys.exit(1)
        to_run = {sp: all_species[sp] for sp in args.species}
    else:
        to_run = all_species

    categories_to_process = args.categories if args.categories is not None else CATEGORIES

    print(f"\n{'=' * 70}")
    print(f"  EDNiX SBA by Category (threshold based on BOLD runs)")
    print(f"  Species: {', '.join(to_run)}")
    print(f"  Categories: {categories_to_process}")
    print(f"  Minimum runs: {args.min_n}")
    print(f"  Homogenize: {args.homogenize}")
    print(f"  Output: {args.output}")
    print(f"{'=' * 70}\n")

    # ── Load classification ───────────────────────────────────────────────────
    all_category_subjects = {}

    for species, bids_dirs in to_run.items():
        print(f"\n{'=' * 70}")
        print(f"  Loading classification for {species}")
        print(f"{'=' * 70}")

        category_subjects = load_classification_subjects(
            args.classification,
            species=species,
            bids_dir=None,
            categories=categories_to_process,
        )

        if category_subjects is None:
            print(f"  [WARN] No classification data for {species}")
            continue

        all_category_subjects[species] = category_subjects

        print(f"\n  Summary for {species}:")
        total_subj = sum(len(subj) for subj in category_subjects.values() if subj is not None)
        print(f"    Total subjects (all categories): {total_subj}")
        for cat in categories_to_process:
            n = len(category_subjects.get(cat, []))
            print(f"    {cat}: {n} subjects")

    # ── Process each species and category ────────────────────────────────────
    print(f"\n{'=' * 70}")
    print(f"  Starting SBA processing (min runs = {args.min_n})")
    print(f"{'=' * 70}\n")

    all_processed = []
    all_skipped = []

    # First pass: full dataset (no homogenization), includes 'All'
    if not args.homogenize:
        for species, bids_dirs in to_run.items():
            if species not in all_category_subjects:
                continue

            category_subjects = all_category_subjects[species]
            categories_to_process_with_all = list(set(categories_to_process + ['All']))

            for cat in categories_to_process_with_all:
                subjects = category_subjects.get(cat, [])
                if not subjects:
                    all_skipped.append((species, cat, 0, "no subjects"))
                    continue

                # We'll let process_species_category check the run count
                try:
                    success = process_species_category(
                        species, bids_dirs, cat, subjects,
                        bids_label=Path(bids_dirs[0]).name,
                        output_suffix="Full",
                        skip_analyses=args.skip_analyses if args.skip_analyses else ['analysis1', 'analysis2'],
                        results_base_dir=args.output,
                        min_runs=args.min_n,
                        exclude_lateralized=args.exclude_lateralized,   # NEW — was missing
                    )
                    if success:
                        all_processed.append((species, cat, len(subjects), "full"))
                    else:
                        all_skipped.append((species, cat, len(subjects), "insufficient runs"))
                except Exception as exc:
                    print(f"\n[ERROR] {species} - {cat} failed: {exc}")
                    import traceback
                    traceback.print_exc()
                    all_skipped.append((species, cat, len(subjects), "error"))

    # Second pass: homogenized (EXCLUDING 'All')
    if args.homogenize:
        print(f"\n{'=' * 70}")
        print(f"  Starting HOMOGENIZED SBA processing")
        print(f"  (Subsampling subjects, min runs still applies)")
        print(f"{'=' * 70}\n")

        for species, bids_dirs in to_run.items():
            if species not in all_category_subjects:
                continue

            category_subjects = all_category_subjects[species]
            categories_for_homogenization = [cat for cat in categories_to_process if cat != 'All']
            category_subjects_filtered = {cat: category_subjects[cat] for cat in categories_for_homogenization}

            subsampled, target_n = subsample_categories(
                category_subjects_filtered,
                target_n=None,
                random_seed=42
            )

            if target_n < 1:
                print(f"\n  [SKIP] {species}: no categories with >=4 subjects")
                continue

            for cat in categories_for_homogenization:
                subjects = subsampled.get(cat, [])
                if not subjects:
                    all_skipped.append((species, cat, 0, "no subjects after subsampling"))
                    continue

                try:
                    success = process_species_category(
                        species, bids_dirs, cat, subjects,
                        bids_label=Path(bids_dirs[0]).name,
                        output_suffix="Homogenized",
                        skip_analyses=args.skip_analyses if args.skip_analyses else ['analysis2'],
                        results_base_dir=args.output,
                        min_runs=args.min_n,
                    )
                    if success:
                        all_processed.append((species, cat, len(subjects), "homogenized"))
                    else:
                        all_skipped.append((species, cat, len(subjects), "insufficient runs (homogenized)"))
                except Exception as exc:
                    print(f"\n[ERROR] {species} - {cat} failed (homogenized): {exc}")
                    import traceback
                    traceback.print_exc()
                    all_skipped.append((species, cat, len(subjects), "error"))

    # ── Final summary ─────────────────────────────────────────────────────────
    print(f"\n{'=' * 70}")
    print(f"  SBA by Category COMPLETE")
    print(f"{'=' * 70}")

    print(f"\n  Processed runs:")
    for species, cat, n_subj, mode in all_processed:
        print(f"    ✅ {species} - {cat}: {n_subj} subjects ({mode})")

    if all_skipped:
        print(f"\n  Skipped runs:")
        for species, cat, n_subj, reason in all_skipped:
            print(f"    ❌ {species} - {cat}: {n_subj} subjects ({reason})")

    print(f"\n  Results saved to: {args.output}")
    print(f"\n{'#' * 70}")
    print("  Done.")
    print(f"{'#' * 70}\n")


if __name__ == '__main__':
    main()
'''