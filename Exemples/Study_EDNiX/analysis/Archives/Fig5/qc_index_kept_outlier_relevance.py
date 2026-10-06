"""
qc_index_kept_outlier_relevance.py
────────────────────────────────────
Given the current pipeline's Kept/Outlier classification (df_class, from
SBA_Region_Profile_Index.classify_runs_by_absolute_correlation or
classify_runs_by_profile_distance) and its per-run region profiles
(profiles, from compute_all_run_profiles), rank a set of candidate QC
variables by how well each one separates Kept from Outlier — same
approach as the earlier ednix_qc_variable_selection.py: Cohen's d
(direction + magnitude) and AUC (Mann-Whitney U — separation strength
regardless of direction).

Candidate variables, computed PER RUN from its own profile_df (already-
available data — no new file I/O needed beyond what SBA_Region_Profile_Index
already produced):
  - global_<stat>: each DEFAULT_REGION_STATS column (median_r, std_r,
    frac_positive, frac_negative, iqr_r if present), aggregated as the
    MEDIAN across every region in that run's profile
  - seed_region_median_r: median_r restricted to the seed's OWN region
    row (self-correlation sanity check — should be high; a low value
    here is a strong signal something's wrong with that specific run).
    Only added if seed_region_hint is given.
  - n_regions_valid: how many regions had ANY voxels at all in that run
    (a low count can flag a badly-masked/cropped/misaligned run)

⚠ Two biological QC metrics from earlier work are NOT recreated here —
spatial contiguity (Dice overlap of the run's significant cluster with
the seed's own anatomical mask) and anticorrelation distance (median r
in the seed mask minus median r in a reference mask, e.g. somatosensory
cortex). Both need each run's THRESHOLDED significance map (the earlier
convention was a `higher_threshold.nii.gz` sibling of the raw
correlation map) and a specific reference mask — neither is part of
what `profiles` already holds. Tell me the current confirmed path for
that thresholded map and the reference mask you want, and I'll add both
as extra candidate columns in build_candidate_qc_table().
"""

import numpy as np
import pandas as pd
from scipy import stats


# =============================================================================
# 1) Build the candidate QC table from already-computed profiles
# =============================================================================

def build_candidate_qc_table(df_class, profiles, seed_region_hint=None):
    """
    One row per run, with candidate QC variables pulled from its
    already-computed profile_df (compute_all_run_profiles's `profiles`
    dict, keyed by (subject, session, run, seed)).

    seed_region_hint: optional case-insensitive substring identifying
    the seed's OWN region row within each profile (e.g.
    "retrosplenial") — if given, adds a seed_region_median_r column.
    Skipped silently if None, since not every use case needs it.
    """
    rows = []
    n_no_profile = 0
    for row in df_class.itertuples():
        key = (str(row.subject), str(row.session), str(row.run), str(row.seed))
        prof = profiles.get(key)
        if prof is None:
            # profiles dict keys may not be str-cast the same way
            # df_class's columns are (int vs str, etc.) -- try the raw
            # (non-str-cast) tuple too before giving up on this run.
            key_raw = (row.subject, row.session, row.run, row.seed)
            prof = profiles.get(key_raw)
        if prof is None:
            n_no_profile += 1
            continue

        out = dict(subject=row.subject, session=row.session, run=row.run,
                  seed=row.seed, is_outlier=bool(row.is_outlier))

        for stat_col in ("median_r", "std_r", "frac_positive", "frac_negative", "iqr_r"):
            if stat_col in prof.columns:
                out[f"global_{stat_col}"] = float(prof[stat_col].median(skipna=True))

        if "n_voxels_valid" in prof.columns:
            out["n_regions_valid"] = int((prof["n_voxels_valid"] > 0).sum())

        if seed_region_hint is not None and "region_name" in prof.columns:
            seed_rows = prof[prof["region_name"].str.contains(seed_region_hint, case=False, na=False)]
            if not seed_rows.empty and "median_r" in seed_rows.columns:
                out["seed_region_median_r"] = float(seed_rows["median_r"].median())

        rows.append(out)

    if n_no_profile:
        print(f"  [WARN] {n_no_profile} run(s) in df_class had no matching entry in "
             f"profiles — skipped (check that df_class and profiles came from the "
             f"same compute_all_run_profiles() call)")

    df = pd.DataFrame(rows)
    if df.empty:
        print("  [WARN] candidate QC table is empty -- nothing to rank")
        return df
    print(f"  Built candidate QC table: {len(df)} run(s), "
         f"{int(df['is_outlier'].sum())} outlier / {int((~df['is_outlier']).sum())} kept")
    return df


# =============================================================================
# 2) Rank candidate variables by separation strength
# =============================================================================

def _cohens_d(a, b):
    """Cohen's d, pooled SD. Positive = group a has the higher mean."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    na, nb = len(a), len(b)
    pooled_sd = np.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    if pooled_sd == 0:
        return np.nan
    return float((np.mean(a) - np.mean(b)) / pooled_sd)


def rank_qc_variables(qc_table, label_col="is_outlier",
                      exclude_cols=("subject", "session", "run", "seed")):
    """
    For every numeric column in qc_table (except exclude_cols/label_col),
    rank how well it separates the two groups defined by label_col
    (kept = False, outlier = True) using:
      - Cohen's d (kept vs outlier means, pooled SD) — sign shows
        direction: POSITIVE = kept has the higher value
      - AUC (Mann-Whitney U / (n_kept*n_outlier)) — 0.5 = no
        separation; towards 0 or 1 = strong separation, regardless of
        direction (use the sign of cohens_d to see which way)

    Returns a DataFrame sorted by |Cohen's d| descending: variable,
    cohens_d, auc, p_value, n_kept, n_outlier, n_missing.
    """
    if qc_table.empty:
        return pd.DataFrame(columns=["variable", "cohens_d", "auc", "p_value",
                                     "n_kept", "n_outlier", "n_missing"])

    candidate_cols = [c for c in qc_table.columns
                      if c not in exclude_cols and c != label_col
                      and pd.api.types.is_numeric_dtype(qc_table[c])]

    kept_mask = ~qc_table[label_col].astype(bool)
    outlier_mask = qc_table[label_col].astype(bool)

    rows = []
    for col in candidate_cols:
        kept_vals = qc_table.loc[kept_mask, col].dropna().values
        outlier_vals = qc_table.loc[outlier_mask, col].dropna().values
        n_missing = int(qc_table[col].isna().sum())

        d = _cohens_d(kept_vals, outlier_vals)

        if len(kept_vals) >= 1 and len(outlier_vals) >= 1:
            try:
                u_stat, p = stats.mannwhitneyu(kept_vals, outlier_vals, alternative="two-sided")
                auc = float(u_stat / (len(kept_vals) * len(outlier_vals)))
            except ValueError:
                auc, p = np.nan, np.nan
        else:
            auc, p = np.nan, np.nan

        rows.append(dict(variable=col, cohens_d=d, auc=auc, p_value=p,
                         n_kept=len(kept_vals), n_outlier=len(outlier_vals), n_missing=n_missing))

    ranked = pd.DataFrame(rows)
    ranked["_abs_cohens_d"] = ranked["cohens_d"].abs()
    ranked = ranked.sort_values("_abs_cohens_d", ascending=False).drop(columns="_abs_cohens_d")
    return ranked.reset_index(drop=True)


# =============================================================================
# 3) Plot
# =============================================================================

def plot_qc_variable_ranking(ranked, top_n=15, output_path=None, dpi=300):
    """Horizontal bar chart of |Cohen's d| for the top_n candidate variables."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if ranked.empty:
        print("  [WARN] nothing to plot -- ranked table is empty")
        return None

    top = ranked.head(top_n).iloc[::-1]  # reverse so the strongest is on top
    fig, ax = plt.subplots(figsize=(8, max(4, 0.4 * len(top))))
    colors = ["#2ca02c" if d > 0 else "#d62728" for d in top["cohens_d"]]
    ax.barh(top["variable"], top["cohens_d"].abs(), color=colors)
    ax.set_xlabel("|Cohen's d|  (green = kept higher, red = outlier higher)")
    ax.set_title(f"Top {len(top)} QC variables by separation strength (Kept vs Outlier)")
    for i, (d, auc) in enumerate(zip(top["cohens_d"], top["auc"])):
        ax.text(abs(d) + 0.02, i, f"AUC={auc:.2f}", va="center", fontsize=8)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi)
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


# =============================================================================
# 4) Orchestration
# =============================================================================

def run_qc_variable_ranking(df_class, profiles, seed_region_hint=None,
                            output_dir=None, label=""):
    """
    One-call pipeline: build the candidate QC table, rank it, plot the
    top variables, and (if output_dir given) save the ranking as Excel
    and the plot as PNG. Returns (qc_table, ranked_df).
    """
    qc_table = build_candidate_qc_table(df_class, profiles, seed_region_hint=seed_region_hint)
    ranked = rank_qc_variables(qc_table)

    if output_dir and not ranked.empty:
        import os
        os.makedirs(output_dir, exist_ok=True)
        excel_path = os.path.join(output_dir, f"{label}_qc_variable_ranking.xlsx" if label
                                  else "qc_variable_ranking.xlsx")
        with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
            ranked.to_excel(writer, sheet_name="ranking", index=False)
            qc_table.to_excel(writer, sheet_name="raw_values", index=False)
        print(f"  [saved] {excel_path}")

        plot_qc_variable_ranking(
            ranked, output_path=os.path.join(
                output_dir, f"{label}_qc_variable_ranking.png" if label else "qc_variable_ranking.png"))

    return qc_table, ranked
