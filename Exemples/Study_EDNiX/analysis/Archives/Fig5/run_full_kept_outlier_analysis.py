"""
run_full_kept_outlier_analysis.py
────────────────────────────────────
Consolidated, top-to-bottom version of everything built across this
conversation for the Kept/Outlier paper analysis. Each STEP prints what
it's doing and SKIPS gracefully (with a clear message) if its
prerequisites aren't met, instead of crashing the rest of the script —
run it whole, or copy individual STEPs into a console.

Toggle USE_FRESH_LAUNCH below depending on whether you need `profiles`/
`ref_profile` (only STEP 7's region boxplots need them): loading from
disk (load_all_results_from_disk) is fast but only gives df_class;
relaunching (launch()) is slower but gives everything, including
profiles/ref_profile.
"""

import os
from pathlib import Path

import Exemples.Study_EDNiX.analysis.Archives.Fig5.Figure5_SBA_KeeporNot as launcher
from Exemples.Study_EDNiX.analysis.Archives.Fig5.outlier_cross_species_report import (
    run_cross_species_report, plot_reference_correlation_distributions,
)
from Exemples.Study_EDNiX.analysis.Fig5.external_qc_validation import (
    load_qc_from_excel, compare_kept_outlier_on_external_metrics,
    plot_external_qc_comparison, compare_reference_strategies,
    plot_reference_strategy_comparison,
)

# ═══════════════════════════════════════════════════════════════════════════
# CONFIG
# ═══════════════════════════════════════════════════════════════════════════

OUT_DIR = "/scratch2/EDNiX/results/EDNiX_Figures"
QC_RAW_XLSX = "/scratch2/EDNiX/results/EDNiX_Figures/Suppl_Figure_4-5-6_QC/qc_raw_per_individual.xlsx"

# False (default): load_all_results_from_disk() -- fast, df_class only.
# True: launcher.launch() -- slower (re-runs the profiling pipeline),
# but needed for STEP 7 (region boxplots, which need `profiles`/`ref_profile`).
USE_FRESH_LAUNCH = False

# Level/atlas for the matrix-based steps (fingerprint, cross-species report).
LEVEL = 2
ATLAS_NAME = "EDNIxCSCLR"
FDR_ALPHA = 0.05
N_PERM = 5000  # 500 for a quick check, 5000-10000 for the paper

REGION_PATTERNS = [("Somatosensory", "Somatosensory"), ("Prefrontal", "Prefrontal")]

opj = os.path.join
os.makedirs(OUT_DIR, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════════
# STEP 1 — Load classification results
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("  STEP 1 — Loading classification results")
print("=" * 70)

if USE_FRESH_LAUNCH:
    all_results = launcher.launch()
else:
    all_results = launcher.load_all_results_from_disk()

bids_dirs_by_label = {cfg["label"]: cfg["bids_dir"] for cfg in launcher.BIDS_CONFIGS}
print(f"  {len(all_results)} (label, seed) result(s) loaded")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 2 — Cross-species report (fig04: FC matrices, fig05/06: fingerprint,
#          fig07: combined connectivity+QC ranking)
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("  STEP 2 — Cross-species report")
print("=" * 70)

out = run_cross_species_report(
    all_results, bids_dirs_by_label,
    output_dir=opj(OUT_DIR, "KeptOutlier_Report"),
    level=LEVEL, atlas_name=ATLAS_NAME, fdr_alpha=FDR_ALPHA, n_top=10,
)

'''
# ═══════════════════════════════════════════════════════════════════════════
# STEP 3 — Permutation test on the fingerprint (fig08)
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("  STEP 3 — Permutation test")
print("=" * 70)

if out.get("fingerprint") is None:
    print("  [SKIP] no fingerprint available from STEP 2 -- need >=2 BIDS with kept+outlier data")
else:
    perm = permutation_test_fingerprint(
        out["matrices_by_bids"], out["common_rois"], n_perm=N_PERM, fdr_alpha=FDR_ALPHA,
    )
    plot_permutation_null_distribution(perm, output_path=opj(OUT_DIR, "08_permutation_test.png"))
'''

# ═══════════════════════════════════════════════════════════════════════════
# STEP 4 — Reference-correlation distributions per BIDS (fig09)
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("  STEP 4 — Reference-correlation distributions")
print("=" * 70)

thresholds = {cfg["label"]: cfg.get("min_correlation") for cfg in launcher.BIDS_CONFIGS}
plot_reference_correlation_distributions(
    all_results, min_correlation_by_label=thresholds,
    output_path=opj(OUT_DIR, "09_reference_correlation_distributions.png"),
)


# ═══════════════════════════════════════════════════════════════════════════
# STEP 5 — External QC validation (fig10): Kept vs Outlier on motion/tSNR/DVARS,
#          metrics NEVER used to build the classification
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("  STEP 5 — External QC validation")
print("=" * 70)

if not os.path.exists(QC_RAW_XLSX):
    print(f"  [SKIP] {QC_RAW_XLSX} not found -- run Suppl_Figure_4-5-6_QC.py "
         f"with its df_qc.to_excel(...) save added first")
    df_qc_by_label = None
else:
    df_qc = load_qc_from_excel(QC_RAW_XLSX)
    print("  bids_dir values found in df_qc:", sorted(df_qc["bids_dir"].astype(str).unique()))

    df_qc_by_label = {
        label: df_qc[df_qc["bids_dir"].astype(str).str.contains(Path(bids_dir).name)]
        for label, bids_dir in bids_dirs_by_label.items()
    }
    empty_labels = [lbl for lbl, d in df_qc_by_label.items() if d.empty]
    if empty_labels:
        print(f"  [WARN] no external QC rows matched for: {empty_labels} -- "
             f"check the bids_dir values printed above against bids_dirs_by_label")

    comp = compare_kept_outlier_on_external_metrics(all_results, df_qc_by_label)
    if not comp.empty:
        plot_external_qc_comparison(comp, output_path=opj(OUT_DIR, "10_external_qc_validation.png"))
    else:
        print("  [SKIP] no comparison rows produced -- nothing to plot")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 6 — Reference-strategy comparison (fig11): individual reference vs
#          group-average reference, validated against the SAME external
#          QC metrics as STEP 5. Reruns the profiling pipeline with
#          reference_subject=None for every BIDS -- heavier than the
#          steps above, and off by default.
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("  STEP 6 — Reference-strategy comparison")
print("=" * 70)

RUN_REFERENCE_STRATEGY_COMPARISON = True  # set True to actually run this (see note above)

if not RUN_REFERENCE_STRATEGY_COMPARISON:
    print("  [SKIP] RUN_REFERENCE_STRATEGY_COMPARISON is False -- set True to run "
         "(reruns the full profiling pipeline with reference_subject=None for every BIDS)")
elif df_qc_by_label is None:
    print("  [SKIP] STEP 5's external QC data wasn't available")
else:
    group_avg_output_dir = launcher.OUTPUT_DIR + "_GroupAvgRef"

    # avoid heavy/destructive side effects for this comparison-only rerun
    launcher.COPY_QC_IMAGES = False
    launcher.RERUN_GROUP_ANALYSIS = False

    all_results_individual = all_results
    all_results_group_avg = {}
    for cfg in launcher.BIDS_CONFIGS:
        cfg_group_avg = dict(cfg)
        cfg_group_avg["reference_subject"] = None  # force the group-average reference
        for seed in cfg["seeds"]:
            result = launcher.run_one_bids_config(cfg_group_avg, seed, output_dir=group_avg_output_dir)
            if result is not None:
                all_results_group_avg[(cfg["label"], seed)] = result

    combined = compare_reference_strategies(
        all_results_individual, all_results_group_avg, df_qc_by_label)
    if not combined.empty:
        plot_reference_strategy_comparison(
            combined, output_path=opj(OUT_DIR, "11_reference_strategy_comparison.png"))
    else:
        print("  [SKIP] no comparison rows produced -- nothing to plot")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 7 — Region boxplots (Somatosensory / Prefrontal), Kept vs Outlier,
#          with the reference marked. Needs `profiles`/`ref_profile`,
#          which are only available when USE_FRESH_LAUNCH=True (STEP 1).
# ═══════════════════════════════════════════════════════════════════════════

from Exemples.Study_EDNiX.analysis.Archives.Fig5.extract_region_median_r import extract_region_median_r_per_run, plot_median_r_kept_vs_outlier, find_region_names_by_species
import Exemples.Study_EDNiX.analysis.Fig5.SBA_Region_Profile_Index as sri

matches = find_region_names_by_species(sri, launcher.BIDS_CONFIGS, "PFC")

for cfg in launcher.BIDS_CONFIGS:
    label, bids_dir = cfg["label"], cfg["bids_dir"]
    label_df = sri.parse_itk_snap_label_file(cfg["label_path"])
    for seed in cfg["seeds"]:
        df_class = all_results[(label, seed)]["df_class"]
        long_df = extract_region_median_r_per_run(
            sri, bids_dir, label_df, df_class,
            region_patterns=[("L_Somatosensory_cortex", "R_Somatosensory_cortex"), ("L_Prefrontal_cortex_(PFC)", "R_Prefrontal_cortex_(PFC)")],
        )
        plot_median_r_kept_vs_outlier(
            long_df, title=f"{label} — {seed}",
            output_path=f"/scratch2/EDNiX/results/EDNiX_Figures/{label}_{seed}_somato_pfc_median_r.png",
        )