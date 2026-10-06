"""
run_SBA_Region_Profile.py
────────────────────────────
Clean, complete driver for the DMN (retrosplenial-seed) region-profile
quality index — ties together everything in SBA_Region_Profile_Index.py
into one runnable pipeline: parse the atlas labels, build the reference
curve from the existing group map, profile every individual run, classify
outliers, plot, and export to Excel.

One seed (L_ or R_ retrosplenial) is processed at a time — they use
different group maps and can have a different number of runs, so keeping
them separate avoids silently mixing two different reference curves into
one distance number. Loop over both if you want both.

Edit the CONFIG block below, then either run this file directly or import
and call run_one_seed() / run_all() from a console.
"""

import os
from pathlib import Path

import pandas as pd

from Exemples.Study_EDNiX.analysis.Fig5 import SBA_Region_Profile_Index as sri

opj = os.path.join

# ═══════════════════════════════════════════════════════════════════════════
# CONFIG — edit this block
# ═══════════════════════════════════════════════════════════════════════════

SPECIES = "Macaque"
BIDS_NAME = "BIDS_Zhu_Garin"
BIDS_DIR = "/scratch2/EDNiX/" + SPECIES + "/" + BIDS_NAME
# Which seeds to process (each gets its own reference curve + classification).
SEEDS = ["BA_26_29_30__retrosplenial"]

LABEL_PATH = ("/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/"
             "carnivorans/Dog/label_code/EDNIxCSCLR_ITK.txt")
CLASSIFICATION_PATH = ("/scratch2/EDNiX/results/EDNiX_Figures/Figure_4_rsfMRI_QC/"
                       "lvl2/all/classification_results.xlsx")



# ⚠ The atlas is 4D / hierarchical (4 levels). None = auto-detect the
# right level PER REGION (recommended — see SBA_Region_Profile_Index.py's
# diagnose_unmatched_regions if coverage looks low for your atlas). Set an
# int (1-4) to force every region to that one level instead.
ATLAS_LEVEL = None

VALUE_COL = "median_r"          # which per-region stat drives the distance/plot
METRIC = "correlation"          # distance metric: correlation | cosine | euclidean
OUTLIER_THRESHOLD = 0.0001        # MAD modified-z cutoff, one-sided ("upper" = worse)
CATEGORY = "All"

# Plot appearance — SMOOTH_WINDOW=1 shows the raw, unsmoothed reference
# curve (no rolling-average flattening); higher values smooth more.
# CI_PCT=None disables the confidence band entirely.
SMOOTH_WINDOW = 5
CI_PCT = 95

# Copy each run's SBA preview jpg into <OUTPUT_DIR>/<species>_<seed>_images/
# kept|outliers/, for visual side-by-side comparison. Off by default since
# it's a file-copy side effect, not just computation.
COPY_QC_IMAGES = False

# Group-map path template — matches the confirmed real layout:
# .../Figure_5_SBA/<species>/<bids_name>/Full/All/Grp_SBA_3dTTEST_AllFull/<seed>/<seed>_correlation_norm.nii.gz
RESULTS_BASE_DIR = "/scratch2/EDNiX/results/EDNiX_Figures/Figure_5_SBA"
OUTPUT_SUFFIX = "Full"

OUTPUT_DIR = "/scratch2/EDNiX/results/EDNiX_Figures/Figure_5_SBA_RegionProfile"

# ═══════════════════════════════════════════════════════════════════════════


def _group_map_path(species, bids_dir, seed):
    bids_name = Path(bids_dir).name
    return opj(RESULTS_BASE_DIR, species, bids_name, OUTPUT_SUFFIX, CATEGORY,
              "Grp_SBA_3dTTEST_AllFull", seed, f"{seed}_correlation_norm.nii.gz")


def run_one_seed(seed, label_df, atlas_vol_path, output_dir=None):
    """Full pipeline for ONE seed: reference curve, run profiles,
    classification, plot, Excel. Returns (ref_profile, profiles, df_class).
    output_dir=None reads the CURRENT value of the module-level OUTPUT_DIR
    at call time (not at function-definition time), so overriding
    run_mod.OUTPUT_DIR after import actually takes effect."""
    output_dir = output_dir if output_dir is not None else OUTPUT_DIR
    print(f"\n{'=' * 70}\n  {SPECIES} / {BIDS_DIR} — {seed}\n{'=' * 70}")

    group_map_path = _group_map_path(SPECIES, BIDS_DIR, seed)
    if not os.path.exists(group_map_path):
        print(f"  [WARN] group map not found, skipping this seed: {group_map_path}")
        return None, {}, pd.DataFrame()

    ref_profile = sri.build_reference_curve(
        group_map_path, atlas_vol_path, label_df, atlas_level=ATLAS_LEVEL)

    profiles, long_df = sri.compute_all_run_profiles(
        SPECIES, BIDS_DIR, label_df, CLASSIFICATION_PATH,
        category=CATEGORY, seeds=[seed], atlas_level=ATLAS_LEVEL)

    if not profiles:
        print(f"  [WARN] no run profiles computed for {seed} — nothing to classify/plot")
        return ref_profile, {}, pd.DataFrame()

    df_class = sri.classify_runs_by_profile_distance(
        profiles, ref_profile, value_col=VALUE_COL, metric=METRIC,
        outlier_threshold=OUTLIER_THRESHOLD)

    os.makedirs(output_dir, exist_ok=True)

    outlier_map = {
        (row.subject, row.session, row.run, row.seed): bool(row.is_outlier)
        for row in df_class.itertuples()
    }
    plot_path = opj(output_dir, f"{SPECIES}_{seed}_profile.png")
    sri.plot_region_profile_curve(
        ref_profile, profiles, value_col=VALUE_COL, is_outlier_map=outlier_map,
        title=f"{SPECIES} / {Path(BIDS_DIR).name} — {seed} ({VALUE_COL})",
        output_path=plot_path, smooth_window=SMOOTH_WINDOW, ci_pct=CI_PCT,
    )

    excel_path = opj(output_dir, f"{SPECIES}_{seed}_region_profile.xlsx")
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        df_class.to_excel(writer, sheet_name="classification", index=False)
        ref_profile.to_excel(writer, sheet_name="reference_profile", index=False)
        long_df.to_excel(writer, sheet_name="all_run_profiles", index=False)
    print(f"  [saved] {excel_path}")

    # Plain CSV of the reference curve too — quicker to open/grep/diff
    # than digging into one sheet of an .xlsx when you just want the
    # numbers.
    ref_csv_path = opj(output_dir, f"{SPECIES}_{seed}_reference_profile.csv")
    ref_profile.to_csv(ref_csv_path, index=False)
    print(f"  [saved] {ref_csv_path}")

    n_out = int(df_class["is_outlier"].sum())
    print(f"  {seed}: {n_out}/{len(df_class)} runs flagged as outliers")

    if COPY_QC_IMAGES:
        images_dir = opj(output_dir, f"{SPECIES}_{seed}_images")
        sri.copy_run_qc_images(df_class, BIDS_DIR, images_dir)

    return ref_profile, profiles, df_class


def run_all(seeds=None):
    """Parse the atlas once, find the representative atlas volume once,
    then run_one_seed() for every seed. Returns {seed: (ref_profile,
    profiles, df_class)}. seeds=None reads the CURRENT value of the
    module-level SEEDS at call time — same reasoning as run_one_seed's
    output_dir."""
    seeds = seeds if seeds is not None else SEEDS
    label_df = sri.parse_itk_snap_label_file(LABEL_PATH)
    atlas_vol_path = sri.pick_representative_atlas_volume(BIDS_DIR)
    if atlas_vol_path is None:
        raise FileNotFoundError(f"No atlas volume found under {BIDS_DIR} — "
                                f"check pick_representative_atlas_volume()'s glob pattern")
    print(f"  Using atlas volume: {atlas_vol_path}")

    results = {}
    for seed in seeds:
        results[seed] = run_one_seed(seed, label_df, atlas_vol_path)
    return results


if __name__ == "__main__":
    results = run_all()
    ref_profile, profiles, df_class = results[str(SEEDS[0])]

df_class = sri.classify_runs_by_absolute_correlation(profiles, ref_profile, value_col=VALUE_COL, min_correlation=0.4)


sri.copy_run_qc_images(df_class, BIDS_DIR, "/scratch2/EDNiX/results/EDNiX_Figures/Figure_5_SBA_RegionProfile/" + BIDS_NAME+ "_" + str(SEEDS[0]) + "_" + "df_class")
sk, so, kept_runs, outlier_runs = sri.rerun_balanced_by_run(
    df_class, SPECIES, BIDS_DIR,
)