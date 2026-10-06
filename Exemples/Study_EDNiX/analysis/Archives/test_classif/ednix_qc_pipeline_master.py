"""
EDNiX — Master QC Orchestrator
==================================

Drop this at the end of your launcher (after your existing
run_pipeline_multilevel(...) call), or run it standalone. It uses the SAME
species_bids_dict / species_multi_bids / ALL_BIDS config pattern as your
production script, and for every (species, bids_dir):

  1. runs the dual-threshold MATRIX classification (ednix_batch_classify_lvl2)
  2. runs the biological SBA QC (ednix_sba_biological_qc)
  3. merges both on (species, bids_dir, subject, session)
  4. writes one master xlsx under OUT_DIR, plus a per-species summary sheet

The lvl2-matrix manifest is discovered ONCE up front (not per bids_dir) and
then filtered per dataset — avoids re-walking all of /scratch2/EDNiX for
every entry in ALL_BIDS.

*** SET THESE THREE PATHS FOR YOUR ENVIRONMENT before running: ***
  BIDS_ROOT, LUT_PATH, OUT_DIR
"""

import os
import warnings

import pandas as pd

from Exemples.Study_EDNiX.analysis.Archives.test_classif import ednix_scoring_core as core
from Exemples.Study_EDNiX.analysis.Archives.test_classif import ednix_batch_classify_lvl2 as batch_matrix
from Exemples.Study_EDNiX.analysis.Archives.test_classif import ednix_sba_biological_qc as bio

opj = os.path.join

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION — mirrors your launcher's species_bids_dict / species_multi_bids
# ═══════════════════════════════════════════════════════════════════════════════

BIDS_ROOT = "/scratch2/EDNiX"
LUT_PATH = "/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/primates/Macaque/label_code/EDNIxCSCLR_label.txt"
OUT_DIR = "/scratch2/EDNiX/results/EDNiX_Figures/Figure_4_rsfMRI_QC"

species_bids_dict = {
    'Rat':        '/scratch2/EDNiX/Rat/BIDS_Grandjean',
    'Mouse':      '/scratch2/EDNiX/Mouse/BIDS_Grandjean2',
    'Dog':        '/scratch2/EDNiX/Dog/BIDS_Boch_K9',
    'Marmoset':   '/scratch2/EDNiX/Marmoset/BIDS_Tian',
    'Mouselemur': '/scratch2/EDNiX/Mouselemur/BIDS_Garin',
}
species_multi_bids = {
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
}

ALL_BIDS = (
    [(sp, bd) for sp, bd in species_bids_dict.items()] +
    [(sp, bd) for sp, bds in species_multi_bids.items() for bd in bds]
)

THRESH_INTRA = core.DEFAULT_THRESH_INTRA
THRESH_INTER = core.DEFAULT_THRESH_INTER
COMMON_ROIS = core.DEFAULT_COMMON_ROIS

# ═══════════════════════════════════════════════════════════════════════════════
# PER-DATASET STEPS
# ═══════════════════════════════════════════════════════════════════════════════

def matrix_classification_from_manifest(full_manifest, species, bids_dir_name, verbose=False):
    """Filters the already-discovered manifest to one dataset, scores + classifies."""
    manifest = full_manifest[(full_manifest["species"] == species) &
                              (full_manifest["bids_dir"] == bids_dir_name)]
    if manifest.empty:
        warnings.warn(f"[matrix] no lvl2 matrices found for {species}/{bids_dir_name}")
        return pd.DataFrame()
    df_scored = batch_matrix.score_all_runs(manifest, common_rois=COMMON_ROIS, verbose=verbose)
    if df_scored.empty:
        return df_scored
    cats = batch_matrix.classify_scored_runs(df_scored, thresh_intra=THRESH_INTRA,
                                              thresh_inter=THRESH_INTER, verbose=verbose)
    df_scored["fp_category"] = cats.values
    return df_scored


def biological_qc_for_dataset(species, bids_dir_name, verbose=False):
    try:
        return bio.batch_compute_sba_qc(BIDS_ROOT, LUT_PATH, species=species,
                                         bids_dir=bids_dir_name, verbose=verbose)
    except Exception as e:
        warnings.warn(f"[bio] failed for {species}/{bids_dir_name}: {e}")
        return pd.DataFrame()


def merge_matrix_and_bio(df_matrix, df_bio):
    if df_matrix.empty:
        return df_bio
    if df_bio.empty:
        return df_matrix

    print("Matrix columns:", df_matrix.columns.tolist())
    print("Bio columns:", df_bio.columns.tolist())
    print("Matrix run_idx sample:", df_matrix["run_idx"].head(5).tolist())
    print("Bio run_idx sample:", df_bio["run_idx"].head(5).tolist())

    key_cols = ["species", "bids_dir", "subject", "session"]
    if "run_idx" in df_matrix.columns and "run_idx" in df_bio.columns:
        key_cols.append("run_idx")
        print("Merging on run_idx")
    elif "run" in df_matrix.columns and "run" in df_bio.columns:
        key_cols.append("run")
        print("Merging on run (fallback)")
    else:
        print("Merging only on subject/session (no run identifier)")

    merged = df_matrix.merge(df_bio, on=key_cols, how="left", suffixes=("", "_bio"))
    return merged

# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def run_full_qc_pipeline(all_bids=ALL_BIDS, bids_root=BIDS_ROOT, out_dir=OUT_DIR, verbose=True):
    os.makedirs(out_dir, exist_ok=True)

    if verbose:
        print("Discovering all lvl2 matrices under", bids_root, "...")
    full_manifest = batch_matrix.discover_all_lvl2_runs(bids_root, verbose=verbose)

    all_results = []
    for species, bids_dir_path in all_bids:
        bids_dir_name = os.path.basename(bids_dir_path.rstrip('/'))
        if verbose:
            print(f"\n{'=' * 60}\n  {species} — {bids_dir_name}\n{'=' * 60}")

        df_matrix = matrix_classification_from_manifest(full_manifest, species, bids_dir_name, verbose=False)
        df_bio = biological_qc_for_dataset(species, bids_dir_name, verbose=False)
        merged = merge_matrix_and_bio(df_matrix, df_bio)

        if merged.empty:
            if verbose:
                print(f"  [skip] no results for {species}/{bids_dir_name}")
            continue

        all_results.append(merged)
        if verbose:
            n_specific = (merged["fp_category"] == "Specific").sum() if "fp_category" in merged.columns else "n/a"
            print(f"  matrix rows: {len(df_matrix):4d}   bio rows: {len(df_bio):4d}   "
                  f"merged: {len(merged):4d}   Specific: {n_specific}")

    if not all_results:
        raise RuntimeError("No results produced for any species/bids_dir — check the warnings above.")

    df_all = pd.concat(all_results, ignore_index=True)

    xlsx_path = opj(out_dir, "combined_qc_matrix_and_biological.xlsx")
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        df_all.to_excel(writer, sheet_name="Combined_Per_Run", index=False)

        summary_rows = []
        for sp in sorted(df_all["species"].unique()):
            sub = df_all[df_all["species"] == sp]
            row = {"species": sp, "n_runs": len(sub)}
            if "fp_category" in sub.columns:
                row["frac_specific"] = round(float((sub["fp_category"] == "Specific").mean()), 3)
            for col in ["median_r_somato", "contiguity_dice", "anticorr_distance", "frac_negative_somato"]:
                if col in sub.columns:
                    row[f"{col}_mean"] = round(float(sub[col].mean()), 3)
            summary_rows.append(row)
        pd.DataFrame(summary_rows).to_excel(writer, sheet_name="Species_Summary", index=False)

        thresholds_note = pd.DataFrame({
            "Parameter": ["thresh_intra (T1)", "thresh_inter (T2)", "common_rois"],
            "Value": [THRESH_INTRA, THRESH_INTER, ", ".join(COMMON_ROIS)],
        })
        thresholds_note.to_excel(writer, sheet_name="Thresholds_Used", index=False)

    if verbose:
        print(f"\n[done] wrote {xlsx_path}  ({len(df_all)} total rows across {len(all_results)} datasets)")

    return df_all


if __name__ == "__main__":
    run_full_qc_pipeline()
