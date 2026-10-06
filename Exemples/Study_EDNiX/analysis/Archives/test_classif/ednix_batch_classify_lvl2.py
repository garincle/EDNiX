"""
EDNiX — Batch lvl2 Classification (scan -> score -> classify -> xlsx)
========================================================================

Two-phase design, deliberately split so changing a threshold or the
common-ROI list doesn't force you to re-walk /scratch2/EDNiX/ and reload
every matrix from disk:

  PHASE A (slow, I/O bound):  discover_all_lvl2_runs()  + score_all_runs()
      Walks the whole tree once, loads every lvl2 correlation matrix,
      computes species_intra / primary_inter / primary_intra (and the
      descriptive species_inter_all) for every run, and writes a
      "raw scores" csv you keep around.

  PHASE B (fast, no I/O):     classify_scored_runs() + export_to_xlsx()
      Re-applies classify_one() to the already-computed scores. If you
      only change thresh_intra / thresh_inter, this is instant. If you
      change common_rois you DO need to re-run Phase A (species_intra
      depends on which ROIs are pooled).

Usage
-----
Full scan + classify + export (first run, or after changing common_rois):
    python ednix_batch_classify_lvl2.py --scan --output classification_lvl2.xlsx

Re-classify only (after changing thresholds, reusing existing raw scores):
    python ednix_batch_classify_lvl2.py --reclassify raw_scores_lvl2.csv \
        --thresh-intra 0.12 --thresh-inter 0.20 --output classification_lvl2_v2.xlsx

As a library:
    import ednix_batch_classify_lvl2 as batch
    manifest = batch.discover_all_lvl2_runs("/scratch2/EDNiX/")
    df_scored = batch.score_all_runs(manifest, common_rois=batch.core.DEFAULT_COMMON_ROIS)
    df_scored.to_csv("raw_scores_lvl2.csv", index=False)   # keep this around
    cats = batch.classify_scored_runs(df_scored)           # fast, re-runnable
    batch.export_to_xlsx(df_scored, cats, "classification_lvl2.xlsx")
"""

import os
import re
import glob
import argparse
import warnings

import numpy as np
import pandas as pd

from Exemples.Study_EDNiX.analysis.Archives.test_classif import ednix_scoring_core as core

# =============================================================================
# PHASE A — DISCOVER + SCORE
# =============================================================================

DEFAULT_BIDS_ROOT = "/scratch2/EDNiX"
DEFAULT_ATLAS_NAME = "EDNIxCSC"
DEFAULT_LEVEL = 2
DEFAULT_LR = True
DEFAULT_FIT_KIND = "correlation"

# Path shape (mirrors _find_all_run_paths in the FC threshold explorer):
#   {bids_root}/{species}/{bids_dir}/sub-{subject}/[ses-{session}/]func/acpc-func/
#       Stats/Correl_matrix/{atlas_name}{LR}/{fit_kind}/
#       {atlas_name}{LR}_{level}_run_{n}[_flattened]_{fit_kind}_matrix.csv
_FNAME_RE_TEMPLATE = (
    r"^{atlas}{lr}_{level}_run_(\d+)(?:_flattened)?_{fit}_matrix\.csv$"
)


def discover_all_lvl2_runs(bids_root=DEFAULT_BIDS_ROOT,
                            atlas_name=DEFAULT_ATLAS_NAME,
                            atlas_level=DEFAULT_LEVEL,
                            use_lr=DEFAULT_LR,
                            fit_kind=DEFAULT_FIT_KIND,
                            verbose=True):
    """
    Walk bids_root and return a manifest DataFrame of every lvl2
    correlation-matrix CSV found, with species / bids_dir / subject /
    session / run parsed from the path.

    species and bids_dir are taken from the two path components right
    under bids_root: {bids_root}/{species}/{bids_dir}/...
    """
    lr = "LR" if use_lr else ""
    fname_re = re.compile(_FNAME_RE_TEMPLATE.format(
        atlas=re.escape(atlas_name), lr=lr, level=atlas_level, fit=re.escape(fit_kind)))
    excl = ("check_fit", "pval", "tstat", "flattened")

    pattern = os.path.join(bids_root, "*", "*", "sub-*", "**",
                            "Correl_matrix", atlas_name + lr, fit_kind, "*.csv")
    all_paths = glob.glob(pattern, recursive=True)

    rows = []
    for p in sorted(set(all_paths)):
        fn = os.path.basename(p)
        if any(x in fn for x in excl):
            continue
        m = fname_re.match(fn)
        if not m:
            continue
        run_idx = int(m.group(1))

        rel = os.path.relpath(p, bids_root)
        parts = rel.split(os.sep)
        # parts[0]=species, parts[1]=bids_dir, parts[2]='sub-XXX', ...
        if len(parts) < 3:
            continue
        species = parts[0]
        bids_dir = parts[1]
        subject = parts[2].replace("sub-", "", 1)
        ses_match = re.search(r"ses-([^/\\]+)", rel)
        session = ses_match.group(1) if ses_match else ""

        rows.append({
            "species": species,
            "bids_dir": bids_dir,
            "subject": subject,
            "session": session,
            "run_idx": run_idx,
            "run": fn,
            "corr_matrix_path": p,
        })

    manifest = pd.DataFrame(rows)
    if manifest.empty:
        raise ValueError(
            f"No lvl{atlas_level} matrices found under {bids_root} matching "
            f"pattern {atlas_name}{lr}_{atlas_level}_run_N_{fit_kind}_matrix.csv — "
            f"check bids_root / atlas_name / atlas_level / fit_kind."
        )
    manifest = manifest.drop_duplicates(subset=["species", "bids_dir", "subject", "session", "run"])

    if verbose:
        print(f"[discover] found {len(manifest)} lvl{atlas_level} runs "
              f"across {manifest['species'].nunique()} species, "
              f"{manifest.groupby(['species','bids_dir']).ngroups} bids dirs")
        for sp in sorted(manifest["species"].unique()):
            print(f"    {sp}: {(manifest['species'] == sp).sum()} runs")

    return manifest


def score_all_runs(manifest, common_rois=None,
                    n_top_homotopic=core.DEFAULT_N_TOP_HOMOTOPIC,
                    use_inter_all=True, verbose=True):
    """
    Phase A scoring: load every matrix in `manifest` and compute the raw
    dual scores (species_intra, primary_inter, primary_intra, ...).
    Does NOT classify — that's Phase B (classify_scored_runs), so the
    same raw-scores table can be reclassified cheaply under different
    thresholds.
    """
    if common_rois is None:
        common_rois = core.DEFAULT_COMMON_ROIS

    out_rows = []
    n_failed = 0
    for _, row in manifest.iterrows():
        path = row["corr_matrix_path"]
        try:
            rois, mat = core.load_matrix(path)
            scores = core.compute_dual_scores_from_matrix(
                mat, rois, common_rois=common_rois,
                n_top_homotopic=n_top_homotopic, use_inter_all=use_inter_all)
        except Exception as e:
            warnings.warn(f"[score] failed on {path}: {e}")
            n_failed += 1
            continue

        out_row = dict(row)
        out_row.update(scores)
        out_rows.append(out_row)

    if verbose:
        print(f"[score] scored {len(out_rows)}/{len(manifest)} runs "
              f"({n_failed} failed to load)")

    df_scored = pd.DataFrame(out_rows)
    return df_scored

# =============================================================================
# PHASE B — CLASSIFY (fast, no disk I/O) + EXPORT
# =============================================================================

def classify_scored_runs(df_scored,
                          thresh_intra=core.DEFAULT_THRESH_INTRA,
                          thresh_inter=core.DEFAULT_THRESH_INTER,
                          verbose=True):
    """Re-applies classify_one() to an already-scored DataFrame. Cheap —
    safe to call repeatedly while you tune thresholds."""
    cats = pd.Series("No", index=df_scored.index, dtype=object)
    for idx, row in df_scored.iterrows():
        cats[idx] = core.classify_one(
            row["species_intra"], row["primary_inter"], row["primary_intra"],
            thresh_intra=thresh_intra, thresh_inter=thresh_inter)

    if verbose:
        print(f"\n[classify] T1(species_intra) > {thresh_intra}, "
              f"T2(primary_inter) < {thresh_inter}")
        if "species" in df_scored.columns:
            for sp in sorted(df_scored["species"].unique()):
                mask = df_scored["species"] == sp
                c = cats[mask].value_counts()
                print(f"    {sp:14s} n={mask.sum():3d}   " +
                      "  ".join(f"{cat}={c.get(cat, 0)}" for cat in core.CATS))
    return cats


def export_to_xlsx(df_scored, cats, output_path,
                    thresh_intra=core.DEFAULT_THRESH_INTRA,
                    thresh_inter=core.DEFAULT_THRESH_INTER,
                    common_rois=None, verbose=True):
    """
    Writes a classification workbook readable by downstream tools
    (e.g. the Step-3 noise-robustness script): one row per run with its
    scores + fp_category, plus per-species and per-subject summaries.
    """
    df = df_scored.copy()
    df["fp_category"] = cats.values if hasattr(cats, "values") else cats

    run_cols = ["species", "bids_dir", "subject", "session", "run", "run_idx",
                "corr_matrix_path", "fp_category",
                "species_intra", "primary_inter", "primary_intra",
                "species_inter_all", "n_homotopic_pairs_all", "n_homotopic_pairs_primary"]
    run_cols = [c for c in run_cols if c in df.columns]
    df_runs = df[run_cols].copy()

    subject_cols = ["species", "bids_dir", "subject"]
    df_subject = df.groupby(subject_cols).agg(
        species_intra=("species_intra", "mean"),
        primary_inter=("primary_inter", "mean"),
        total_runs=("fp_category", "count"),
    ).reset_index()
    for cat in core.CATS:
        df_subject[f"n_{cat}"] = df.groupby(subject_cols)["fp_category"].apply(
            lambda s: (s == cat).sum()).values
    df_subject["frac_specific"] = df_subject["n_Specific"] / df_subject["total_runs"]

    species_summary = df.groupby("species").agg(
        n_subjects=("subject", "nunique"),
        total_runs=("fp_category", "count"),
        species_intra_mean=("species_intra", "mean"),
        primary_inter_mean=("primary_inter", "mean"),
    ).reset_index()
    for cat in core.CATS:
        species_summary[f"n_{cat}"] = df.groupby("species")["fp_category"].apply(
            lambda s: (s == cat).sum()).values
    species_summary["frac_specific"] = species_summary["n_Specific"] / species_summary["total_runs"]

    thresholds_note = pd.DataFrame({
        "Parameter": ["thresh_intra (T1)", "thresh_inter (T2)", "common_rois"],
        "Value": [thresh_intra, thresh_inter,
                  ", ".join(common_rois) if common_rois else ", ".join(core.DEFAULT_COMMON_ROIS)],
    })

    os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        df_runs.to_excel(writer, sheet_name="Run_Classification", index=False)
        df_subject.to_excel(writer, sheet_name="Subject_Summary", index=False)
        species_summary.to_excel(writer, sheet_name="Species_Summary", index=False)
        thresholds_note.to_excel(writer, sheet_name="Thresholds_Used", index=False)

    if verbose:
        print(f"\n[export] wrote {output_path}")
        print(f"    {len(df_runs)} runs, {df_subject.shape[0]} subjects, "
              f"{species_summary.shape[0]} species")

    return df_runs, df_subject, species_summary

# =============================================================================
# CONVENIENCE: full pipeline in one call
# =============================================================================

def run_full_scan_and_classify(bids_root=DEFAULT_BIDS_ROOT,
                                common_rois=None,
                                thresh_intra=core.DEFAULT_THRESH_INTRA,
                                thresh_inter=core.DEFAULT_THRESH_INTER,
                                raw_scores_csv=None,
                                output_xlsx="classification_lvl2.xlsx",
                                verbose=True):
    """Phase A + Phase B in one call. Pass raw_scores_csv to also save the
    raw (unclassified) scores table for fast re-classification later."""
    manifest = discover_all_lvl2_runs(bids_root, verbose=verbose)
    df_scored = score_all_runs(manifest, common_rois=common_rois, verbose=verbose)
    if raw_scores_csv:
        df_scored.to_csv(raw_scores_csv, index=False)
        if verbose:
            print(f"[save] raw scores -> {raw_scores_csv} (reuse for fast reclassification)")
    cats = classify_scored_runs(df_scored, thresh_intra=thresh_intra,
                                 thresh_inter=thresh_inter, verbose=verbose)
    return export_to_xlsx(df_scored, cats, output_xlsx,
                           thresh_intra=thresh_intra, thresh_inter=thresh_inter,
                           common_rois=common_rois, verbose=verbose)

# =============================================================================
# CLI
# =============================================================================

def main():
    p = argparse.ArgumentParser(description="EDNiX lvl2 batch classifier (dual-threshold).")
    p.add_argument("--bids-root", default=DEFAULT_BIDS_ROOT)
    p.add_argument("--scan", action="store_true",
                    help="Full rescan of bids-root + score + classify + export.")
    p.add_argument("--reclassify", metavar="RAW_SCORES_CSV",
                    help="Reuse an existing raw-scores CSV (from a prior --scan) and just "
                         "reclassify with new thresholds — no disk rescan.")
    p.add_argument("--thresh-intra", type=float, default=core.DEFAULT_THRESH_INTRA)
    p.add_argument("--thresh-inter", type=float, default=core.DEFAULT_THRESH_INTER)
    p.add_argument("--raw-scores-out", default="raw_scores_lvl2.csv",
                    help="Where to save raw scores when using --scan.")
    p.add_argument("--output", default="classification_lvl2.xlsx")
    args = p.parse_args()

    if args.reclassify:
        df_scored = pd.read_csv(args.reclassify)
        cats = classify_scored_runs(df_scored, thresh_intra=args.thresh_intra,
                                     thresh_inter=args.thresh_inter)
        export_to_xlsx(df_scored, cats, args.output,
                        thresh_intra=args.thresh_intra, thresh_inter=args.thresh_inter)
    elif args.scan:
        run_full_scan_and_classify(
            bids_root=args.bids_root,
            thresh_intra=args.thresh_intra, thresh_inter=args.thresh_inter,
            raw_scores_csv=args.raw_scores_out, output_xlsx=args.output)
    else:
        p.error("Specify either --scan (full rescan) or --reclassify RAW_SCORES_CSV (fast).")


if __name__ == "__main__":
    main()
