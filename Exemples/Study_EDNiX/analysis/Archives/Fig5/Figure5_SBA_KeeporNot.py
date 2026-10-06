"""
launch_SBA_MultiConfig.py
────────────────────────────
Runs the DMN region-profile QC pipeline across MULTIPLE BIDS directories,
each with its own config (species, label path, min_correlation threshold,
seeds) — then builds cross-BIDS comparison figures on top:

  1. Per-BIDS: kept vs outlier average profile (impact of the selection).
  2. Inter-species: overlay each BIDS's KEPT-only average profile on the
     common EDNiX lvl3 regions, to see how close the profiles are.
  3. Inter-species: distance of every BIDS's profile to a reference
     (e.g. Human), on the same common regions.

Edit BIDS_CONFIGS below, then run this file directly or call launch()
from a console.
"""

import os
import re
import glob
from pathlib import Path

import pandas as pd
import nibabel as nib

from Exemples.Study_EDNiX.analysis.Fig5 import SBA_Region_Profile_Index as sri

opj = os.path.join
ope = os.path.exists

# ═══════════════════════════════════════════════════════════════════════════
# CONFIG — one entry per BIDS directory. Fill in real values before running.
# ═══════════════════════════════════════════════════════════════════════════

# REFERENCE MAP CONFIG — optional per-BIDS keys (each entry in
# BIDS_CONFIGS below can set these; if unset, the group map is used,
# same as before):
#   reference_subject: None (default) uses the GROUP map. Set to a
#     subject id (e.g. "01") to use THAT ONE individual subject's raw
#     correlation map as the reference instead — e.g. to check every
#     other run against one hand-picked "gold standard" individual
#     rather than the group average.
#   reference_session / reference_run: only relevant when
#     reference_subject is set — pin down WHICH session/run of that
#     subject to use. None (default) auto-picks the first one found
#     (printed either way, so which exact file was used is always
#     visible).

BIDS_CONFIGS = [
    # Single-species datasets (one BIDS directory per species)
        # Single-species datasets (one BIDS directory per species)
    dict(
        label="Dog",
        species="Dog",
        bids_dir="/scratch2/EDNiX/Dog/BIDS_Boch_K9",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "carnivorans/Dog/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        # Example (commented out): use ONE individual subject's raw
        # correlation map as the reference instead of the group map —
        # see REFERENCE MAP CONFIG below for the full explanation.
        reference_subject=None,
        reference_session=None,   # None = auto-pick first session found
        reference_run=None,       # None = auto-pick first run found
    ),
    dict(
        label="Macaque_BenHamed",
        species="Macaque",
        bids_dir="/scratch2/EDNiX/Macaque/BIDS_BenHamed",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "primates/Macaque/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject=None,
        reference_session=None,   # None = auto-pick first session found
        reference_run=None,       # None = auto-pick first run found
    ),
    dict(
        label="Macaque_ZhuGarin",
        species="Macaque",
        bids_dir="/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "primates/Macaque/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject=None,
        reference_session=None,   # None = auto-pick first session found
        reference_run=None,       # None = auto-pick first run found
    ),
    dict(
        label="Rat",
        species="Rat",
        bids_dir="/scratch2/EDNiX/Rat/BIDS_Grandjean",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "rodents/Rat/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject=None,
        reference_session=None,   # None = auto-pick first session found
        reference_run=None,       # None = auto-pick first run found
    ),
    dict(
        label="Mouse",
        species="Mouse",
        bids_dir="/scratch2/EDNiX/Mouse/BIDS_Grandjean2",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "rodents/Mouse/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject=None,
        reference_session=None,   # None = auto-pick first session found
        reference_run=None,       # None = auto-pick first run found
    ),
    dict(
        label="Marmoset_Tian",
        species="Marmoset",
        bids_dir="/scratch2/EDNiX/Marmoset/BIDS_Tian",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "primates/Marmoset/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject=None,
        reference_session=None,   # None = auto-pick first session found
        reference_run=None,       # None = auto-pick first run found
    ),
    dict(
        label="Marmoset_Hori",
        species="Marmoset",
        bids_dir="/scratch2/EDNiX/Marmoset/BIDS_Hori",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "primates/Marmoset/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject=None,
        reference_session=None,   # None = auto-pick first session found
        reference_run=None,       # None = auto-pick first run found
    ),
    dict(
        label="Mouselemur",
        species="Mouselemur",
        bids_dir="/scratch2/EDNiX/Mouselemur/BIDS_Garin",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "primates/Mouselemur/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject=None,
        reference_session=None,   # None = auto-pick first session found
        reference_run=None,       # None = auto-pick first run found
    ),

    dict(
        label="Human_Castrillon",
        species="Human",
        bids_dir="/scratch2/EDNiX/Human/BIDS_Castrillon",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "primates/Human/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject=None,
        reference_session=None,   # None = auto-pick first session found
        reference_run=None,       # None = auto-pick first run found
    ),
    dict(
        label="Human_Park",
        species="Human",
        bids_dir="/scratch2/EDNiX/Human/BIDS_Park",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "primates/Human/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject=None,
        reference_session=None,   # None = auto-pick first session found
        reference_run=None,       # None = auto-pick first run found
    ),
]

'''
BIDS_CONFIGS = [
    # Single-species datasets (one BIDS directory per species)
        # Single-species datasets (one BIDS directory per species)
    dict(
        label="Dog",
        species="Dog",
        bids_dir="/scratch2/EDNiX/Dog/BIDS_Boch_K9",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "carnivorans/Dog/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        # Example (commented out): use ONE individual subject's raw
        # correlation map as the reference instead of the group map —
        # see REFERENCE MAP CONFIG below for the full explanation.
        reference_subject="22",
        reference_session="1",  # None = auto-pick first session found
        reference_run=None,  # None = auto-pick first run found
    ),
    dict(
        label="Macaque_BenHamed",
        species="Macaque",
        bids_dir="/scratch2/EDNiX/Macaque/BIDS_BenHamed",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "primates/Macaque/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject="Samourai",
        reference_session="1",   # None = auto-pick first session found
        reference_run="2",       # None = auto-pick first run found
    ),
    dict(
        label="Macaque_ZhuGarin",
        species="Macaque",
        bids_dir="/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "primates/Macaque/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject="Oliver",
        reference_session="11",   # None = auto-pick first session found
        reference_run="01",       # None = auto-pick first run found
    ),
    dict(
        label="Rat",
        species="Rat",
        bids_dir="/scratch2/EDNiX/Rat/BIDS_Grandjean",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "rodents/Rat/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject="300101",
        reference_session="1",   # None = auto-pick first session found
        reference_run=None,       # None = auto-pick first run found
    ),
    dict(
        label="Mouse",
        species="Mouse",
        bids_dir="/scratch2/EDNiX/Mouse/BIDS_Grandjean2",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "rodents/Mouse/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject="jgrAesMEDISOc2NT",
        reference_session="2",   # None = auto-pick first session found
        reference_run="2",       # None = auto-pick first run found
    ),
    dict(
        label="Marmoset_Tian",
        species="Marmoset",
        bids_dir="/scratch2/EDNiX/Marmoset/BIDS_Tian",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "primates/Marmoset/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject="NIHm27",
        reference_session="02",  # None = auto-pick first session found
        reference_run="3",  # None = auto-pick first run found
    ),
    dict(
        label="Marmoset_Hori",
        species="Marmoset",
        bids_dir="/scratch2/EDNiX/Marmoset/BIDS_Hori",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "primates/Marmoset/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject="03",
        reference_session="01",   # None = auto-pick first session found
        reference_run="06",       # None = auto-pick first run found
    ),
    dict(
        label="Mouselemur",
        species="Mouselemur",
        bids_dir="/scratch2/EDNiX/Mouselemur/BIDS_Garin",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "primates/Mouselemur/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject="283CCA",
        reference_session="02",   # None = auto-pick first session found
        reference_run=None,       # None = auto-pick first run found
    ),

    dict(
        label="Human_Castrillon",
        species="Human",
        bids_dir="/scratch2/EDNiX/Human/BIDS_Castrillon",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "primates/Human/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject="007",
        reference_session="open",  # None = auto-pick first session found
        reference_run=None,  # None = auto-pick first run found
    ),
    dict(
        label="Human_Park",
        species="Human",
        bids_dir="/scratch2/EDNiX/Human/BIDS_Park",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "primates/Human/label_code/EDNIxCSCLR_ITK.txt"),
        percentile=25,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject="1600",
        reference_session="2",   # None = auto-pick first session found
        reference_run="2",       # None = auto-pick first run found
    ),
]
'''
'''
# Multi‑directory species (Macaque)
dict(
    label="Macaque_Tremblay",
    species="Macaque",
    bids_dir="/scratch2/EDNiX/Macaque/BIDS_Tremblay",
    label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                "primates/Macaque/label_code/EDNIxCSCLR_ITK.txt"),
    min_correlation=0.6,
    seeds=["BA_26_29_30__retrosplenial"],
),

'''

'''
BIDS_CONFIGS = [
    # Single-species datasets (one BIDS directory per species)
        # Single-species datasets (one BIDS directory per species)
    dict(
        label="Dog",
        species="Dog",
        bids_dir="/scratch2/EDNiX/Dog/BIDS_Boch_K9",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "carnivorans/Dog/label_code/EDNIxCSCLR_ITK.txt"),
        min_correlation=0.5,
        seeds=["BA_26_29_30__retrosplenial"],
        # Example (commented out): use ONE individual subject's raw
        # correlation map as the reference instead of the group map —
        # see REFERENCE MAP CONFIG below for the full explanation.
        reference_subject="22",
        reference_session="1",  # None = auto-pick first session found
        reference_run=None,  # None = auto-pick first run found
    ),
    dict(
        label="Macaque_BenHamed",
        species="Macaque",
        bids_dir="/scratch2/EDNiX/Macaque/BIDS_BenHamed",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "primates/Macaque/label_code/EDNIxCSCLR_ITK.txt"),
        min_correlation=0.6,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject="Samourai",
        reference_session="1",   # None = auto-pick first session found
        reference_run="2",       # None = auto-pick first run found
    ),
    dict(
        label="Macaque_ZhuGarin",
        species="Macaque",
        bids_dir="/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "primates/Macaque/label_code/EDNIxCSCLR_ITK.txt"),
        min_correlation=0.6,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject="Oliver",
        reference_session="11",   # None = auto-pick first session found
        reference_run="01",       # None = auto-pick first run found
    ),
    dict(
        label="Rat",
        species="Rat",
        bids_dir="/scratch2/EDNiX/Rat/BIDS_Grandjean",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "rodents/Rat/label_code/EDNIxCSCLR_ITK.txt"),
        min_correlation=0.8,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject="300101",
        reference_session="1",   # None = auto-pick first session found
        reference_run=None,       # None = auto-pick first run found
    ),
    dict(
        label="Mouse",
        species="Mouse",
        bids_dir="/scratch2/EDNiX/Mouse/BIDS_Grandjean2",
        label_path=("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
                    "rodents/Mouse/label_code/EDNIxCSCLR_ITK.txt"),
        min_correlation=0.7,
        seeds=["BA_26_29_30__retrosplenial"],
        reference_subject="jgrAesMEDISOc2NT",
        reference_session="2",   # None = auto-pick first session found
        reference_run="2",       # None = auto-pick first run found
    ),

]
'''


CLASSIFICATION_PATH = ("/scratch2/EDNiX/results/EDNiX_Figures/Figure_4_rsfMRI_QC/"
                       "lvl2/all/classification_results.xlsx")
CATEGORY = "All"
VALUE_COL = "median_r"
METRIC = "correlation"

# ⚠ Same reasoning as run_SBA_Region_Profile.py — None auto-detects the
# atlas level PER REGION. Set an int to force one level everywhere.
ATLAS_LEVEL = 3

# Interspecies comparison figures use this level's regions for the common-
# region intersection — matches "EDNiX lvl3" from the request. Only
# matters for readability of the comparison; profile extraction itself
# still respects ATLAS_LEVEL above.
INTERSPECIES_ATLAS_LEVEL_LABEL = "lvl2"

# Which BIDS_CONFIGS entry to use as the reference for the
# distance-to-reference figure. Accepts EITHER an exact "label"
# (e.g. "Human_Park") OR a "species" name (e.g. "Human") — if given a
# species name, the FIRST BIDS_CONFIGS entry of that species (in
# BIDS_CONFIGS order) with kept data is used automatically; with
# multiple BIDS for that species (e.g. Human_Park + Human_Castrillon),
# the others are still shown in the comparison figure, just not used as
# the reference point. Must resolve to something with kept data, or that
# figure is skipped with a clear message.
REFERENCE_LABEL = "Human"

RESULTS_BASE_DIR = "/scratch2/EDNiX/results/EDNiX_Figures/Figure_5_SBA"
OUTPUT_SUFFIX = "Full"
OUTPUT_DIR = "/scratch2/EDNiX/results/EDNiX_Figures/Figure_5_SBA"

# Which stat map to build each BIDS's reference curve from:
#   "ttest-stat_fisher_cmap" — raw (non-normalized) correlation (default)
#   "correlation_norm"       — normalized r, clipped to [-1, 1]
#   "ttest-stat_fisher_zmap" — Fisher z-stat (NOT bounded — if you switch
#                              to this, every BIDS's min_correlation
#                              threshold no longer means the same thing,
#                              since z-stats aren't on a [-1,1] scale)
# Applies to every BIDS_CONFIGS entry uniformly — set per-BIDS by editing
# _group_map_path's call below if you ever need it to differ by BIDS.
STAT_SUFFIX = "ttest-stat_fisher_cmap"

# Whether _3dttest_EDNiX should drop lateralized (L_.../R_...) regions
# and keep only the merged/bilateral one — set to True to study
# BA_26_29_30__retrosplenial as ONE bilateral region across every BIDS,
# without also processing L_.../R_... separately. Confirmed (via
# report_lateralization()) that every BIDS's label file has BOTH a
# merged entry AND lateralized ones for this seed — so True is safe
# everywhere here and doesn't drop any BIDS's results to zero. My
# earlier per-BIDS caution (a "some BIDS only have lateralized entries"
# warning) was based on stale evidence from before an unrelated
# atlas-file-name bugfix — re-run report_lateralization() any time you
# add a new BIDS/seed to confirm this still holds before assuming True
# is safe for it too.
# This is the DEFAULT used for any BIDS_CONFIGS entry that doesn't set
# its own "exclude_lateralized" key — applies to RERUN_GROUP_ANALYSIS /
# rerun_group_analysis_for_all only (the group map RESOLUTION for the
# main profiling step already auto-detects lateralized-vs-merged via
# find_group_map_path).
EXCLUDE_LATERALIZED = True

# Copy each run's SBA preview jpg into <bids_out_dir>/images/kept|outliers/,
# for visual side-by-side comparison. Off by default (file-copy side
# effect, not just computation).
COPY_QC_IMAGES = True

# Re-run the group 3dttest++ separately on the balanced kept vs outlier
# runs (sri.rerun_balanced_by_run) for EACH BIDS/seed, to statistically
# compare the two groups' resulting group maps. Off by default — this
# launches real AFNI jobs per BIDS, not just computation on what's
# already on disk.
RERUN_GROUP_ANALYSIS = True

# Use sri.classify_runs_iterative_robust for EVERY BIDS_CONFIGS entry,
# instead of whichever classification_method each entry already sets
# (or defaults to). False (default) leaves every entry's own
# classification_method untouched -- set True to switch every BIDS to
# the iterative-robust classification for this launch() call, without
# editing BIDS_CONFIGS itself. Each entry still needs its own
# "percentile" or "min_correlation" set (the per-iteration threshold) --
# this flag only overrides WHICH method is used, not the threshold.
USE_ITERATIVE_ROBUST = False

# Optional label identifying THIS run's configuration (e.g. "percentile25",
# "iterative_robust", "group_ref_v2") -- when set, every output this
# script produces (Excel/plots/images under OUTPUT_DIR, AND the AFNI
# group-rerun results under RESULTS_BASE_DIR) is nested under this tag
# instead of overwriting whatever the LAST run with a different
# configuration wrote to the same paths. None (default, unchanged
# behavior) keeps the original flat layout -- set this whenever you
# want to compare two configurations (e.g. different percentile/
# min_correlation, or classification_method) side by side afterward,
# rather than losing the earlier one.
RUN_TAG = "percentile25_mean"

# Defaults for classification_method="progressive_shrink" (see
# run_one_bids_config's docstring) -- used by any BIDS_CONFIGS entry
# that doesn't set its own "keep_fraction_schedule" / "min_n_kept".
# DEFAULT_KEEP_FRACTION_SCHEDULE: fractions tried in order, typically
# loosest to strictest -- at each step the reference is rebuilt from
# the CURRENT kept set (same self-purifying mechanism as
# classify_runs_iterative_robust), then reclassified at the NEXT
# fraction. DEFAULT_MIN_N_KEPT: an ABSOLUTE COUNT floor -- a step that
# would keep fewer runs than this is refused, and the PREVIOUS step's
# result is kept instead. This protects small BIDS from shrinking all
# the way down to a near-individual reference (which would recreate the
# "arbitrary individual reference" problem this whole approach avoids),
# while letting large BIDS shrink much further in real terms.
DEFAULT_KEEP_FRACTION_SCHEDULE = [0.75, 0.5, 0.35, 0.25, 0.15]
DEFAULT_MIN_N_KEPT = 5

# ═══════════════════════════════════════════════════════════════════════════


def _group_map_path(species, bids_dir, seed):
    """
    Resolves via sri.find_group_map_path — does NOT assume whether the
    group results folder uses a lateralized name (L_BA_.../R_BA_...) or
    a merged bilateral name (BA_..., no prefix). Confirmed real layout
    varies by BIDS: Dog's Grp_SBA_3dTTEST_AllFull/ contains the
    unprefixed "BA_26_29_30__retrosplenial" folder, while other BIDS use
    separate L_.../R_... folders — this works with `seed` given in
    either form. stat_suffix comes from module-level STAT_SUFFIX (read
    at call time — see COPY_QC_IMAGES/RERUN_GROUP_ANALYSIS's reasoning
    elsewhere in this file for why globals are read this way).
    """
    bids_name = Path(bids_dir).name
    return sri.find_group_map_path(RESULTS_BASE_DIR, species, bids_name, seed,
                                   stat_suffix=STAT_SUFFIX,
                                   category=CATEGORY, output_suffix=OUTPUT_SUFFIX)


def _run_token_matches(actual_token, requested_run):
    """
    Flexible run-token match for the reference-run filter below: exact
    match first, then match on trailing digits only after stripping an
    "RL-"/"LR-" acquisition-direction prefix if present on EITHER side.
    Confirmed real gap this fixes: Marmoset's run tokens are always
    "RL-<N>"/"LR-<N>" (SBA_Region_Profile_Index._extract_run_id never
    returns a bare digit for a multi-run session there), so a
    BIDS_CONFIGS entry with reference_run="3" never matched the real
    "RL-3" token under the old exact-string comparison — the reference
    lookup silently fell through to "no matching raw correlation map
    found" for every session tried. Never matches the "NA" sentinel
    (no run token in the filename at all).
    """
    if actual_token == "NA":
        return False
    requested_str = str(requested_run)
    if actual_token == requested_str:
        return True
    actual_digits = re.sub(r"^(RL|LR)-", "", actual_token)
    requested_digits = re.sub(r"^(RL|LR)-", "", requested_str)
    return actual_digits == requested_digits


def _grids_match(map_path, atlas_path):
    """
    True if map_path and atlas_path share the same spatial voxel grid
    shape. Confirmed real failure mode this guards against: different
    subjects within the SAME BIDS can have their EDNIxCSCLR atlas
    volume on a slightly different grid (confirmed real case: Human_Park
    sub-1225's atlas was 53x63x52 vs a reference map's 53x63x53, a
    1-slice difference) — silently substituting one subject's atlas for
    another's map (the pick_representative_atlas_volume fallback used
    when a specific subject's own atlas can't be found) can pair
    mismatched grids, which previously only surfaced much later as an
    opaque ValueError deep inside compute_region_profile, far from the
    actual cause. Checking here lets _resolve_reference_map refuse a
    mismatched fallback immediately with a clear, specific message.
    """
    try:
        map_shape = nib.load(map_path).shape[:3]
        atlas_shape = nib.load(atlas_path).shape[:3]
        return tuple(map_shape) == tuple(atlas_shape)
    except Exception as e:
        print(f"  [WARN] could not check grid match ({map_path} vs {atlas_path}): {e}")
        return True  # don't block on a check failure itself; let downstream code raise its own error


def _find_matching_atlas_volume(bids_dir, map_path, preferred_atlas_path):
    """
    Return an atlas volume path whose spatial grid matches map_path's —
    preferring preferred_atlas_path if it ALREADY matches (fast path, no
    behavior change for the common case where every subject in a BIDS
    shares one grid, which is most of them). If it doesn't match, scans
    every subject under bids_dir for ANY atlas volume on a matching
    grid, since the confirmed real failure (Human_Park: sub-1225's atlas
    one slice narrower than the group map) affects only SOME subjects,
    not necessarily all — using a different subject's matching atlas
    lets this BIDS still succeed instead of failing outright, exactly
    the same way pick_representative_atlas_volume itself already picks
    an arbitrary subject's atlas as a stand-in for the whole BIDS (this
    just also checks that the stand-in is actually usable).

    Used for BOTH _resolve_reference_map branches (group-average
    reference AND an individual reference_subject's own atlas lookup
    failing) — the group-average branch previously skipped this check
    entirely (returned default_atlas_vol_path unconditionally), which is
    the confirmed real gap this fixes: the individual-reference branch
    already refused a mismatched fallback, but the group-average branch
    crashed deep inside compute_region_profile instead, since it never
    checked at all.

    Returns None (with a printed reason) if nothing under bids_dir
    matches map_path's grid.
    """
    if preferred_atlas_path is not None and _grids_match(map_path, preferred_atlas_path):
        return preferred_atlas_path

    print(f"  [INFO] preferred atlas volume ({preferred_atlas_path}) does not match "
         f"{Path(map_path).name}'s grid -- scanning other subjects under {bids_dir} "
         f"for one that does")
    sub_dirs = sorted(glob.glob(opj(bids_dir, "sub-*")))
    for sub_dir in sub_dirs:
        subject = Path(sub_dir).name.replace("sub-", "")
        ses_dirs = sorted(glob.glob(opj(sub_dir, "ses-*")))
        sessions = [Path(d).name.replace("ses-", "") for d in ses_dirs] or ["1"]
        for ses in sessions:
            candidate = sri.find_atlas_volume_path(bids_dir, subject, ses)
            if candidate is not None and _grids_match(map_path, candidate):
                print(f"  [INFO] found a matching atlas volume from sub-{subject}/ses-{ses}: "
                     f"{candidate}")
                return candidate

    print(f"  [WARN] no atlas volume anywhere under {bids_dir} matches "
         f"{Path(map_path).name}'s grid — giving up on this reference.")
    return None


def _resolve_reference_map(cfg, seed, bids_dir, default_atlas_vol_path):
    """
    Resolves the reference map path + the atlas volume to parcellate it
    with, for ONE BIDS config. Two modes, selected by
    cfg.get("reference_subject"):

      - None (default): the GROUP map, via _group_map_path.
        default_atlas_vol_path (pick_representative_atlas_volume's
        result) is used IF its grid matches the group map's — otherwise
        _find_matching_atlas_volume scans every subject for one that
        does, rather than pairing a mismatched grid and crashing later
        (confirmed real gap: this branch previously used
        default_atlas_vol_path unconditionally, with no grid check at
        all).

      - a subject id: that ONE INDIVIDUAL subject's raw (non-Fisher-z)
        correlation map for this seed is used as the reference instead
        of the group map — e.g. to check every other run against one
        specific, hand-picked "gold standard" individual rather than
        the group average. cfg["reference_session"] /
        cfg["reference_run"] optionally pin down WHICH session/run of
        that subject to use (auto-picks the first one found, printed
        either way, if not given). Uses THAT subject's own atlas volume
        (find_atlas_volume_path) if its grid matches, else
        _find_matching_atlas_volume scans for a substitute the same way
        as the group-average branch. reference_run is matched via
        _run_token_matches (see above) rather than an exact string
        comparison, so a bare digit still matches an "RL-"/"LR-"-tagged
        token.

    Returns (map_path, atlas_vol_path) — map_path is None if nothing
    could be resolved (with a printed reason, including the run tokens
    that WERE found for that subject/session if the run filter is what
    failed — rather than a bare "nothing found" that doesn't say why).
    """
    reference_subject = cfg.get("reference_subject")
    if reference_subject is None:
        group_map = _group_map_path(cfg["species"], bids_dir, seed)
        if group_map is None:
            return None, default_atlas_vol_path
        atlas_vol_path = _find_matching_atlas_volume(bids_dir, group_map, default_atlas_vol_path)
        if atlas_vol_path is None:
            # No atlas anywhere in this BIDS matches the group map's grid --
            # there is nothing usable to parcellate with, so this is a full
            # resolution failure (None map), not "map found, atlas missing".
            return None, None
        return group_map, atlas_vol_path

    reference_session = cfg.get("reference_session")
    reference_run = cfg.get("reference_run")

    sessions = [reference_session] if reference_session is not None else None
    if sessions is None:
        ses_dirs = sorted(glob.glob(opj(bids_dir, f"sub-{reference_subject}", "ses-*")))
        sessions = [Path(d).name.replace("ses-", "") for d in ses_dirs] or ["1"]

    for ses in sessions:
        candidates = sri.find_raw_corr_maps_for_run(bids_dir, reference_subject, ses, seed)
        if not candidates:
            continue
        if reference_run is not None:
            run_tokens = {c: sri._extract_run_id(Path(c).name) for c in candidates}
            matched = [c for c in candidates if _run_token_matches(run_tokens[c], reference_run)]
            if not matched:
                print(f"  [INFO] reference_subject={reference_subject!r} ses-{ses}: "
                     f"run={reference_run!r} not found -- available run token(s) here: "
                     f"{sorted(set(run_tokens.values()))}")
                continue
            candidates = matched
        chosen = sorted(candidates)[0]
        preferred_atlas = sri.find_atlas_volume_path(bids_dir, reference_subject, ses)
        print(f"  [INFO] reference_subject={reference_subject!r}: using "
             f"{Path(chosen).name} (ses-{ses}) as the reference map "
             f"instead of the group map")
        if preferred_atlas is None:
            print(f"  [WARN] no atlas volume found for reference subject "
                 f"{reference_subject}/ses-{ses} — will look for a substitute "
                 f"with a matching grid")

        atlas_vol_path = _find_matching_atlas_volume(bids_dir, chosen, preferred_atlas)
        if atlas_vol_path is None:
            print(f"  [WARN] reference_subject={reference_subject!r} ses-{ses}: no atlas "
                 f"volume with a matching grid found anywhere in this BIDS -- trying "
                 f"the next session/run candidate, if any.")
            continue
        return chosen, atlas_vol_path

    print(f"  [WARN] reference_subject={reference_subject!r} "
         f"(session={reference_session!r}, run={reference_run!r}): no matching "
         f"raw correlation map found for seed={seed!r} — check the subject id, "
         f"session, and run.")
    return None, default_atlas_vol_path


def run_one_bids_config(cfg, seed, output_dir=None):
    """
    Full per-BIDS pipeline for ONE seed: reference curve, run profiles,
    classification, plot, Excel.

    Classification method — controlled by cfg.get("classification_method"),
    per-BIDS (each BIDS_CONFIGS entry can use a different one):
      - "absolute" (default, unchanged behavior): sri.classify_runs_by_absolute_correlation,
        using this cfg's own "min_correlation" — a fixed correlation-to-
        reference cutoff, meaning the same thing regardless of how many
        runs this BIDS has. Needs a single fixed reference (group map or
        an individual reference_subject) resolved via
        _resolve_reference_map — skipped if that fails.
      - "percentile": sri.classify_runs_by_percentile, using this cfg's
        own "percentile" (0-100, e.g. 50 keeps the best half of THIS
        BIDS's own runs) — always keeps the same FRACTION of runs for
        this BIDS specifically, regardless of the absolute correlation
        values reached. Also needs a single fixed reference, same as
        "absolute".
      - "iterative_robust": sri.classify_runs_iterative_robust — needs
        EXACTLY ONE of "min_correlation" or "percentile" set (same
        meaning as the corresponding method above, applied at EVERY
        iteration). Unlike the other two methods, this one needs NO
        single fixed reference at all: it starts from the median profile
        across every run and iteratively purifies it using only the
        currently-kept runs — no group map to resolve, no individual
        subject to justify. _resolve_reference_map / build_reference_curve
        are skipped entirely for this method, so a BIDS whose group map
        has a grid-mismatched atlas (see _find_matching_atlas_volume) can
        still be classified this way even if "absolute"/"percentile"
        would have to skip it.
    "percentile" is interpreted as keep_fraction = percentile / 100 in
    every method that uses it.

    Returns dict(label_df=, atlas_vol_path=, ref_profile=, profiles=,
    df_class=) or None if profiles (or, for "absolute"/"percentile", the
    reference map) couldn't be resolved.
    """
    output_dir = output_dir if output_dir is not None else OUTPUT_DIR
    if RUN_TAG:
        output_dir = opj(output_dir, RUN_TAG)
    label, species, bids_dir = cfg["label"], cfg["species"], cfg["bids_dir"]
    print(f"\n{'=' * 70}\n  {label} ({species} / {Path(bids_dir).name}) — {seed}\n{'=' * 70}")

    label_df = sri.parse_itk_snap_label_file(cfg["label_path"])
    rep_atlas_vol_path = sri.pick_representative_atlas_volume(bids_dir)
    if rep_atlas_vol_path is None:
        print(f"  [WARN] no atlas volume found under {bids_dir} — skipping {label}")
        return None
    atlas_vol_path = rep_atlas_vol_path  # used below for the per-run profiling step

    classification_method = cfg.get("classification_method")
    if classification_method is None:
        # Auto-infer: a BIDS_CONFIGS entry with "percentile" set but no
        # explicit "classification_method" almost certainly meant
        # percentile mode (confirmed real symptom this fixes: an entry
        # with percentile set but classification_method left unset fell
        # through to the "absolute" default, which then KeyError'd on
        # min_correlation being absent -- a confusing failure for an
        # easy-to-make config mistake). Only auto-infers when
        # min_correlation is ALSO absent, so an entry that genuinely
        # has both keys set (e.g. left over from switching modes back
        # and forth) still defaults to the explicit, unambiguous
        # "absolute" behavior rather than guessing.
        if "percentile" in cfg and "min_correlation" not in cfg:
            classification_method = "percentile"
            print(f"  [INFO] {label}: 'classification_method' not set, but "
                 f"'percentile' is present and 'min_correlation' is not -- "
                 f"inferring classification_method='percentile'")
        else:
            classification_method = "absolute"

    # profiles are needed by EVERY classification method (including as
    # iterative_robust's/progressive_shrink's own reference source) --
    # compute once, before resolving any single fixed reference (which
    # neither of those two methods needs at all).
    profiles, long_df = sri.compute_all_run_profiles(
        species, bids_dir, label_df, CLASSIFICATION_PATH,
        category=CATEGORY, seeds=[seed], atlas_level=ATLAS_LEVEL)
    if not profiles:
        print(f"  [WARN] no run profiles computed for {label}/{seed} — skipping")
        return None

    bids_out_dir = opj(output_dir, label, seed)
    os.makedirs(bids_out_dir, exist_ok=True)

    if classification_method == "iterative_robust":
        percentile = cfg.get("percentile")
        min_correlation = cfg.get("min_correlation")
        if (percentile is None) == (min_correlation is None):
            raise ValueError(
                f"{label}: classification_method='iterative_robust' needs "
                f"EXACTLY ONE of 'min_correlation' or 'percentile' set in "
                f"this BIDS_CONFIGS entry (the per-iteration threshold)."
            )

        def _save_iteration(i, df_class_i, ref_profile_i):
            # Cheap-only per-iteration snapshot: Excel + the profile plot
            # (matplotlib, no AFNI/surface) -- so you can see HOW the
            # classification evolved and its impact on the profile curve,
            # iteration by iteration, without paying for a full group
            # rerun at every step. The final converged result still gets
            # the FULL finalization (Excel/plots/jpg copy/AFNI+surface)
            # below, exactly like "absolute"/"percentile" -- this is
            # purely additional, intermediate bookkeeping.
            iter_dir = opj(bids_out_dir, "iterations", f"iter_{i}")
            os.makedirs(iter_dir, exist_ok=True)
            excel_path_i = opj(iter_dir, f"{label}_{seed}_iter{i}.xlsx")
            with pd.ExcelWriter(excel_path_i, engine="openpyxl") as writer:
                df_class_i.to_excel(writer, sheet_name="classification", index=False)
                ref_profile_i.to_excel(writer, sheet_name="reference_profile", index=False)

            outlier_map_i = {
                (row.subject, row.session, row.run, row.seed): bool(row.is_outlier)
                for row in df_class_i.itertuples()
            }
            sri.plot_region_profile_curve(
                ref_profile_i, profiles, value_col=VALUE_COL, is_outlier_map=outlier_map_i,
                title=f"{label} — {seed} — iteration {i} ({VALUE_COL})",
                output_path=opj(iter_dir, f"{label}_{seed}_iter{i}_profile.png"),
            )
            try:
                sri.plot_kept_vs_outlier_average(
                    profiles, df_class_i, value_col=VALUE_COL,
                    title=f"{label} — {seed} — iteration {i} — kept vs outlier average",
                    output_path=opj(iter_dir, f"{label}_{seed}_iter{i}_kept_vs_outlier.png"),
                )
            except ValueError as e:
                print(f"  [iter {i}] kept-vs-outlier plot skipped: {e}")

            n_kept_i = int((~df_class_i["is_outlier"]).sum())
            print(f"  [iter {i}] saved -> {iter_dir} ({n_kept_i}/{len(df_class_i)} kept)")

        if percentile is not None:
            result = sri.classify_runs_iterative_robust(
                profiles, value_col=VALUE_COL, metric=METRIC,
                keep_fraction=percentile / 100, use_median=True,
                on_iteration=_save_iteration)
        else:
            result = sri.classify_runs_iterative_robust(
                profiles, value_col=VALUE_COL,
                min_correlation=min_correlation, use_median=True,
                on_iteration=_save_iteration)
        df_class = result["df_class"]
        ref_profile = result["ref_profile"]  # the CONVERGED reference, not a fixed map
        print(f"  [INFO] iterative_robust: converged={result['converged']} in "
             f"{result['n_iterations']} iteration(s)")

    elif classification_method == "progressive_shrink":
        keep_fraction_schedule = cfg.get("keep_fraction_schedule", DEFAULT_KEEP_FRACTION_SCHEDULE)
        min_n_kept = cfg.get("min_n_kept", DEFAULT_MIN_N_KEPT)

        def _save_step(step_i, frac, df_class_i, ref_profile_i):
            # Same cheap-only bookkeeping as iterative_robust's
            # _save_iteration above, one folder per SCHEDULE STEP
            # (not per within-step convergence pass, since
            # classify_runs_progressive_shrink itself only calls this
            # once per fraction actually applied). Folder/file names use
            # the PERCENTAGE (e.g. "iter_75") rather than the step index,
            # since that's what actually varies from step to step here.
            pct = int(round(frac * 100))
            step_dir = opj(bids_out_dir, "iterations", f"iter_{pct}")
            os.makedirs(step_dir, exist_ok=True)
            excel_path_i = opj(step_dir, f"{label}_{seed}_iter{pct}.xlsx")
            with pd.ExcelWriter(excel_path_i, engine="openpyxl") as writer:
                df_class_i.to_excel(writer, sheet_name="classification", index=False)
                ref_profile_i.to_excel(writer, sheet_name="reference_profile", index=False)

            outlier_map_i = {
                (row.subject, row.session, row.run, row.seed): bool(row.is_outlier)
                for row in df_class_i.itertuples()
            }
            sri.plot_region_profile_curve(
                ref_profile_i, profiles, value_col=VALUE_COL, is_outlier_map=outlier_map_i,
                title=f"{label} — {seed} — iter_{pct} (keep_fraction={frac}, {VALUE_COL})",
                output_path=opj(step_dir, f"{label}_{seed}_iter{pct}_profile.png"),
            )
            try:
                sri.plot_kept_vs_outlier_average(
                    profiles, df_class_i, value_col=VALUE_COL,
                    title=f"{label} — {seed} — iter_{pct} (keep_fraction={frac}) "
                         f"— kept vs outlier average",
                    output_path=opj(step_dir, f"{label}_{seed}_iter{pct}_kept_vs_outlier.png"),
                )
            except ValueError as e:
                print(f"  [iter_{pct}] kept-vs-outlier plot skipped: {e}")

            n_kept_i = int((~df_class_i["is_outlier"]).sum())
            print(f"  [iter_{pct}] saved -> {step_dir} ({n_kept_i}/{len(df_class_i)} kept)")

        result = sri.classify_runs_progressive_shrink(
            profiles, keep_fraction_schedule=keep_fraction_schedule,
            min_n_kept=min_n_kept, value_col=VALUE_COL, metric=METRIC,
            use_median=True, on_step=_save_step)
        df_class = result["df_class"]
        ref_profile = result["ref_profile"]
        print(f"  [INFO] progressive_shrink: stopped ({result['stopped_reason']}) after "
             f"{result['n_steps']} step(s) -- final keep_fraction={result['keep_fraction']}")

        # Deviation of EVERY run (not just the final kept set) to the
        # "last average" over the final surviving group --
        # final_kept_reference is ONE MORE average/median pass over
        # whoever was still kept at the very last successful step (see
        # classify_runs_progressive_shrink's own docstring).
        sri.plot_deviation_to_reference(
            profiles, result["final_kept_reference"], df_class=df_class,
            value_col=VALUE_COL, metric=METRIC,
            title=f"{label} — {seed} — deviation to final reference "
                 f"(keep_fraction={result['keep_fraction']})",
            output_path=opj(bids_out_dir, f"{label}_{seed}_deviation_to_final_reference.png"),
        )

    else:
        ref_map_path, ref_atlas_vol_path = _resolve_reference_map(
            cfg, seed, bids_dir, rep_atlas_vol_path)
        if not ref_map_path or not os.path.exists(ref_map_path):
            print(f"  [WARN] reference map not found, skipping: {ref_map_path}")
            return None
        if not ref_atlas_vol_path or not os.path.exists(ref_atlas_vol_path):
            # Defensive backstop: _resolve_reference_map should already refuse
            # to return a None/missing atlas (see its own grid-matching logic),
            # but a missing atlas paired with nib.load() crashes with a
            # confusing TypeError far from the actual cause -- catch it here
            # too, with a message that says what's actually missing.
            print(f"  [WARN] no usable atlas volume resolved for the reference map "
                 f"({ref_map_path}), skipping: {ref_atlas_vol_path}")
            return None

        ref_profile = sri.build_reference_curve(ref_map_path, ref_atlas_vol_path, label_df,
                                                atlas_level=ATLAS_LEVEL)

        if classification_method == "percentile":
            percentile = cfg.get("percentile")
            if percentile is None:
                raise ValueError(
                    f"{label}: classification_method='percentile' but no "
                    f"'percentile' key set in this BIDS_CONFIGS entry."
                )
            print(f"  [INFO] classification_method='percentile': keeping the best "
                 f"{percentile:.0f}% of {label}'s own runs (keep_fraction="
                 f"{percentile / 100:.2f})")
            df_class = sri.classify_runs_by_percentile(
                profiles, ref_profile, value_col=VALUE_COL, metric=METRIC,
                keep_fraction=percentile / 100)
        elif classification_method == "absolute":
            min_correlation = cfg.get("min_correlation")
            if min_correlation is None:
                raise ValueError(
                    f"{label}: classification_method='absolute' (the default) but "
                    f"no 'min_correlation' key set in this BIDS_CONFIGS entry -- "
                    f"either add 'min_correlation', or set "
                    f"classification_method='percentile' with a 'percentile' key, "
                    f"or classification_method='iterative_robust'/'progressive_shrink', "
                    f"if you meant one of those instead."
                )
            df_class = sri.classify_runs_by_absolute_correlation(
                profiles, ref_profile, value_col=VALUE_COL,
                min_correlation=min_correlation)
        else:
            raise ValueError(
                f"{label}: unknown classification_method={classification_method!r} "
                f"-- must be 'absolute', 'percentile', 'iterative_robust', "
                f"or 'progressive_shrink'."
            )

    if classification_method == "percentile":
        classif_desc = f"percentile={cfg.get('percentile')}"
    elif classification_method == "iterative_robust":
        classif_desc = (f"iterative_robust(percentile={cfg.get('percentile')}, "
                        f"min_correlation={cfg.get('min_correlation')})")
    elif classification_method == "progressive_shrink":
        classif_desc = f"progressive_shrink(final_keep_fraction={result['keep_fraction']})"
    else:
        classif_desc = f"min_correlation={cfg.get('min_correlation')}"

    outlier_map = {
        (row.subject, row.session, row.run, row.seed): bool(row.is_outlier)
        for row in df_class.itertuples()
    }
    sri.plot_region_profile_curve(
        ref_profile, profiles, value_col=VALUE_COL, is_outlier_map=outlier_map,
        title=f"{label} — {seed} ({VALUE_COL}, {classif_desc})",
        output_path=opj(bids_out_dir, f"{label}_{seed}_profile.png"),
    )

    try:
        sri.plot_kept_vs_outlier_average(
            profiles, df_class, value_col=VALUE_COL,
            title=f"{label} — {seed} — kept vs outlier average",
            output_path=opj(bids_out_dir, f"{label}_{seed}_kept_vs_outlier.png"),
        )
    except ValueError as e:
        print(f"  [WARN] kept-vs-outlier plot skipped: {e}")

    excel_path = opj(bids_out_dir, f"{label}_{seed}_region_profile.xlsx")
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        df_class.to_excel(writer, sheet_name="classification", index=False)
        ref_profile.to_excel(writer, sheet_name="reference_profile", index=False)
    print(f"  [saved] {excel_path}")

    n_out = int(df_class["is_outlier"].sum())
    print(f"  {label}/{seed}: {n_out}/{len(df_class)} runs flagged as outliers ({classif_desc})")

    if COPY_QC_IMAGES:
        images_dir = opj(bids_out_dir, "images")
        sri.copy_run_qc_images(df_class, bids_dir, images_dir)

    if RERUN_GROUP_ANALYSIS:
        exclude_lateralized = cfg.get("exclude_lateralized", EXCLUDE_LATERALIZED)
        tag_suffix = f"_{RUN_TAG}" if RUN_TAG else ""
        try:
            sri.rerun_balanced_by_run(
                df_class, species, bids_dir, seed=seed,
                output_suffix_kept=f"{seed}_KeptRunBalanced{tag_suffix}",
                output_suffix_outlier=f"{seed}_OutlierRunBalanced{tag_suffix}",
                exclude_lateralized=exclude_lateralized,
            )
        except Exception as e:
            print(f"  [WARN] group rerun failed for {label}/{seed}: {e}")

    return dict(label_df=label_df, atlas_vol_path=atlas_vol_path,
               ref_profile=ref_profile, profiles=profiles, df_class=df_class)


def launch(bids_configs=None, seeds_override=None):
    """
    Run every BIDS config (each with its own min_correlation), then build
    the cross-BIDS comparison figures. Returns
    {(label, seed): result_dict} plus writes the comparison figures to
    OUTPUT_DIR/_interspecies_comparison/.

    If USE_ITERATIVE_ROBUST is True, every cfg's classification_method
    is overridden to "iterative_robust" for this call (each entry's own
    "percentile"/"min_correlation" is left as-is and used as the
    per-iteration threshold).
    """
    bids_configs = bids_configs if bids_configs is not None else BIDS_CONFIGS

    if USE_ITERATIVE_ROBUST:
        for cfg in bids_configs:
            cfg["classification_method"] = "iterative_robust"
        print("  [INFO] USE_ITERATIVE_ROBUST=True -- every BIDS_CONFIGS entry's "
             "classification_method overridden to 'iterative_robust' for this launch()")

    all_results = {}
    for cfg in bids_configs:
        seeds = seeds_override if seeds_override is not None else cfg["seeds"]
        for seed in seeds:
            result = run_one_bids_config(cfg, seed)
            if result is not None:
                all_results[(cfg["label"], seed)] = result

    if not all_results:
        print("\n[WARN] no BIDS config produced results — nothing to compare")
        return all_results

    _build_comparison_figures(all_results, bids_configs)
    return all_results


def _resolve_reference_label(kept_avg_profiles, bids_configs, reference_label):
    """
    reference_label may be an exact BIDS_CONFIGS "label" (e.g.
    "Human_Park") or a "species" name (e.g. "Human") — if it doesn't
    match a label directly, fall back to the FIRST BIDS_CONFIGS entry (in
    list order) whose species equals reference_label and that has kept
    data. Returns the resolved label (a key of kept_avg_profiles), or
    None if nothing matches.
    """
    if reference_label in kept_avg_profiles:
        return reference_label
    candidates = [cfg["label"] for cfg in bids_configs
                 if cfg["species"] == reference_label and cfg["label"] in kept_avg_profiles]
    if candidates:
        extra = f" (other {reference_label} BIDS also in the comparison: {candidates[1:]})" \
                if len(candidates) > 1 else ""
        print(f"  [INFO] REFERENCE_LABEL={reference_label!r} matched by species -> "
             f"using {candidates[0]!r}{extra}")
        return candidates[0]
    return None


def _build_comparison_figures(all_results, bids_configs):
    """
    Builds, per seed (comparisons only make sense within the SAME seed):
      - inter-species profile overlay (kept-only average, common regions)
      - distance-to-REFERENCE_LABEL bar chart
    """
    comp_dir = opj(OUTPUT_DIR, "_interspecies_comparison")
    os.makedirs(comp_dir, exist_ok=True)

    seeds_present = sorted({seed for (_, seed) in all_results})
    for seed in seeds_present:
        print(f"\n{'#' * 70}\n  Inter-species comparison — {seed}\n{'#' * 70}")

        label_dfs, kept_avg_profiles = {}, {}
        for cfg in bids_configs:
            key = (cfg["label"], seed)
            if key not in all_results:
                continue
            r = all_results[key]
            label_dfs[cfg["label"]] = r["label_df"]

            kept_keys = {
                (str(row.subject), str(row.session), str(row.run), str(row.seed))
                for row in r["df_class"].itertuples() if not bool(row.is_outlier)
            }
            profiles_str_keys = {tuple(str(x) for x in k): v for k, v in r["profiles"].items()}
            avg = sri.compute_group_average_profile(profiles_str_keys, kept_keys)
            if avg is None:
                print(f"  [WARN] {cfg['label']}: no kept runs for {seed} — excluded "
                      f"from comparison")
                continue
            kept_avg_profiles[cfg["label"]] = avg

        if len(kept_avg_profiles) < 2:
            print(f"  [WARN] fewer than 2 BIDS with kept data for {seed} — "
                 f"skipping comparison figures")
            continue

        common_bases = sri.find_common_base_regions(label_dfs)
        print(f"  {len(common_bases)} common base region(s) across "
             f"{list(kept_avg_profiles)} ({INTERSPECIES_ATLAS_LEVEL_LABEL})")

        sri.plot_interspecies_profile_comparison(
            kept_avg_profiles, common_bases, value_col=VALUE_COL,
            title=f"Inter-species retrosplenial profile — {seed} "
                 f"({INTERSPECIES_ATLAS_LEVEL_LABEL}, kept-only average)",
            output_path=opj(comp_dir, f"interspecies_profile_{seed}.png"),
        )

        resolved_ref = _resolve_reference_label(kept_avg_profiles, bids_configs, REFERENCE_LABEL)
        if resolved_ref is None:
            print(f"  [WARN] REFERENCE_LABEL={REFERENCE_LABEL!r} has no kept data "
                 f"for {seed} (no matching label or species) — distance-to-reference "
                 f"figure skipped")
            continue

        sri.plot_distance_to_reference_species(
            kept_avg_profiles, common_bases, reference_label=resolved_ref,
            value_col=VALUE_COL, metric=METRIC,
            title=f"Distance to {resolved_ref} — {seed} "
                 f"({INTERSPECIES_ATLAS_LEVEL_LABEL}, kept-only average)",
            output_path=opj(comp_dir, f"distance_to_{resolved_ref}_{seed}.png"),
        )


def rerun_group_analysis_for_all(all_results, min_n_per_group=1):
    """
    Loop the balanced kept-vs-outlier group rerun (sri.rerun_balanced_by_run)
    across every (label, seed) already in `all_results` (e.g. from a prior
    launch() call) — lets you trigger the (heavy, AFNI-launching) group
    reruns as a SEPARATE pass, after inspecting the classification first,
    instead of only via the RERUN_GROUP_ANALYSIS flag during launch().

    Each call is restricted to its OWN seed (never the module's full
    REGIONS_OF_INTEREST list), and uses that BIDS_CONFIGS entry's own
    "exclude_lateralized" override if it has one, else the module-level
    EXCLUDE_LATERALIZED default — same reasoning as run_one_bids_config's
    inline RERUN_GROUP_ANALYSIS path.

    min_n_per_group: skip a (label, seed) if either its kept or outlier
    group has fewer runs than this (rerun_balanced_by_run itself already
    requires both non-empty, but very small groups make for a weak group
    map — this lets you set a stricter floor).

    Returns {(label, seed): (success_kept, success_outlier, kept_runs, outlier_runs)}
    for whichever (label, seed) were actually attempted.
    """
    cfg_by_label = {cfg["label"]: cfg for cfg in BIDS_CONFIGS}
    rerun_results = {}

    for (label, seed), r in all_results.items():
        cfg = cfg_by_label.get(label)
        if cfg is None:
            print(f"  [WARN] {label}/{seed}: no matching BIDS_CONFIGS entry — skipped")
            continue

        df_class = r["df_class"]
        n_kept = int((~df_class["is_outlier"]).sum())
        n_out = int(df_class["is_outlier"].sum())
        if n_kept < min_n_per_group or n_out < min_n_per_group:
            print(f"  [SKIP] {label}/{seed}: kept={n_kept}, outlier={n_out} "
                 f"(min_n_per_group={min_n_per_group})")
            continue

        exclude_lateralized = cfg.get("exclude_lateralized", EXCLUDE_LATERALIZED)
        tag_suffix = f"_{RUN_TAG}" if RUN_TAG else ""
        print(f"\n{'#' * 70}\n  Group rerun — {label} ({cfg['species']}) — {seed} "
             f"(exclude_lateralized={exclude_lateralized})\n{'#' * 70}")
        try:
            rerun_results[(label, seed)] = sri.rerun_balanced_by_run(
                df_class, cfg["species"], cfg["bids_dir"], seed=seed,
                output_suffix_kept=f"{seed}_KeptRunBalanced{tag_suffix}",
                output_suffix_outlier=f"{seed}_OutlierRunBalanced{tag_suffix}",
                exclude_lateralized=exclude_lateralized,
            )
        except Exception as e:
            print(f"  [WARN] group rerun failed for {label}/{seed}: {e}")

    return rerun_results


def report_lateralization(bids_configs=None, seed=None):
    """
    For each BIDS_CONFIGS entry, load the SAME "EDNIxCSC" label file
    rerun_balanced_by_run itself uses for seed matching, and report
    whether `seed` matches a merged/bilateral entry, lateralized (L_/R_)
    entries, or both — WITHOUT launching any AFNI job. Use this to
    decide exclude_lateralized per-BIDS from actual evidence instead of
    guessing/forcing one value everywhere: a BIDS with ONLY a merged
    entry means exclude_lateralized is irrelevant there (nothing to
    exclude either way); a BIDS with ONLY lateralized entries would
    produce NOTHING if exclude_lateralized=True is forced on it (every
    match gets filtered out); a BIDS with BOTH needs a real choice
    depending on whether you want one merged profile or two
    hemisphere-specific ones.

    seed: which seed to check — defaults to each BIDS_CONFIGS entry's
    own cfg["seeds"][0] if not given.

    Returns a DataFrame: label, species, seed, n_merged, n_lateralized,
    recommendation. Also printed.
    """
    from Statistics.Group_fMRI._Group_anal_3dTtest import _format_seed

    bids_configs = bids_configs if bids_configs is not None else BIDS_CONFIGS
    rows = []
    for cfg in bids_configs:
        this_seed = seed if seed is not None else cfg["seeds"][0]
        species = cfg["species"]
        try:
            species_cfg = sri.fig5.build_species_cfg(species)
            label_path = sri.fig5.get_atlas_label_path(
                sri.fig5.ATLAS_LIB, species_cfg["species_fragment"], "EDNIxCSC")
            label_df = sri._load_label_df_for_seed_matching(label_path, [this_seed])
        except Exception as e:
            rows.append(dict(label=cfg["label"], species=species, seed=this_seed,
                             n_merged=None, n_lateralized=None,
                             recommendation=f"ERROR: {e}"))
            continue

        matched = label_df[label_df["region_name"].apply(
            lambda b: any(r in _format_seed(b) for r in [this_seed]))]
        lateralized = matched[matched["region_name"].apply(
            lambda b: _format_seed(b).startswith(("L_", "R_")))]
        merged = matched[~matched.index.isin(lateralized.index)]
        n_merged, n_lat = len(merged), len(lateralized)

        if n_merged and not n_lat:
            rec = "exclude_lateralized doesn't matter (nothing lateralized to exclude)"
        elif n_lat and not n_merged:
            rec = "exclude_lateralized=True would produce NOTHING here -- use False"
        elif n_merged and n_lat:
            rec = "both exist -- True=merged only, False=merged + both hemispheres separately"
        else:
            rec = "seed not found at all in this label file"
        rows.append(dict(label=cfg["label"], species=species, seed=this_seed,
                         n_merged=n_merged, n_lateralized=n_lat, recommendation=rec))

    df = pd.DataFrame(rows)
    print(df.to_string(index=False))
    return df


def load_all_results_from_disk(bids_configs=None, output_dir=None):
    """
    Reload {(label, seed): {"df_class": ...}} from the classification
    Excel files a PRIOR launch() call already wrote to disk
    (<output_dir>/<label>/<seed>/<label>_<seed>_region_profile.xlsx,
    "classification" sheet) — instead of rerunning the entire profiling
    pipeline (SBA profile extraction, classification, jpg copying, and
    any AFNI reruns if RERUN_GROUP_ANALYSIS=True) just to rebuild
    all_results for a downstream report like
    kept_outlier_cross_species_report.run_cross_species_report.

    Only df_class is reconstructed this way — it's the only piece that
    report (and rerun_group_analysis_for_all) actually reads from
    all_results. label_df / atlas_vol_path / ref_profile / profiles are
    NOT reloaded (launch() doesn't save those to a reloadable format);
    if you need those too, either rerun launch() or read the same
    Excel's "reference_profile" sheet yourself.

    Returns {(label, seed): {"df_class": df_class}} for every BIDS/seed
    whose Excel file exists — a missing or unreadable file is printed
    and skipped rather than raising, so a partial/incomplete launch()
    run doesn't block reloading what IS there.

    If RUN_TAG was set when the classification was produced, pass the
    SAME output_dir=opj(OUTPUT_DIR, that_tag) here (or just set RUN_TAG
    again before calling this with no output_dir) — otherwise this looks
    in the untagged OUTPUT_DIR and won't find the tagged run's files.
    """
    bids_configs = bids_configs if bids_configs is not None else BIDS_CONFIGS
    output_dir = output_dir if output_dir is not None else OUTPUT_DIR
    if output_dir == OUTPUT_DIR and RUN_TAG:
        output_dir = opj(output_dir, RUN_TAG)

    all_results = {}
    for cfg in bids_configs:
        label = cfg["label"]
        for seed in cfg["seeds"]:
            excel_path = opj(output_dir, label, seed, f"{label}_{seed}_region_profile.xlsx")
            if not ope(excel_path):
                print(f"  [WARN] {label}/{seed}: no classification file found at "
                     f"{excel_path} -- skipped")
                continue
            try:
                # keep_default_na=False: WITHOUT this, pandas treats the
                # literal string "NA" as a missing-value marker and
                # silently turns it into a real float NaN on reload —
                # but "NA" is OUR OWN sentinel string (from
                # _extract_run_id) meaning "no run token found in this
                # filename", not a missing value. Confirmed real
                # symptom: every "NA"-run row came back as "run=nan",
                # breaking find_lvl3_matrix_path's run-token matching
                # entirely for those rows.
                df_class = pd.read_excel(excel_path, sheet_name="classification",
                                         keep_default_na=False, na_values=[])
                for col in ("subject", "session", "run", "seed"):
                    if col in df_class.columns:
                        df_class[col] = df_class[col].astype(str)
            except Exception as e:
                print(f"  [WARN] {label}/{seed}: could not read {excel_path}: {e}")
                continue
            all_results[(label, seed)] = dict(df_class=df_class)
            print(f"  {label}/{seed}: loaded {len(df_class)} run(s) from {excel_path}")
    return all_results


if __name__ == "__main__":
    launch()