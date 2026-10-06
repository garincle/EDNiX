"""
compare_reference_strategies.py
──────────────────────────────────
Empirically decides WHICH reference-construction strategy works best,
per BIDS, judged against EXTERNAL QC metrics (motion/tSNR/DVARS) that
never enter the reference construction or the classification -- so the
comparison is non-circular, unlike judging a reference by how well runs
correlate to it (which the reference itself defines).

Strategies compared (all from SBA_Region_Profile_Index):
  - "uniform_mean":   compute_group_average_profile over every run
  - "uniform_median": compute_group_median_profile over every run
  - "medoid":         compute_medoid_reference (single most central run)
  - "weighted":       compute_weighted_reference (centrality-weighted)
  - "kcore":          compute_kcore_reference (medoid + k-1 neighbours)

For each strategy: build the reference, classify runs against it at a
fixed threshold, then measure how well that Kept/Outlier split separates
on each external QC metric. The strategy with the strongest separation
on metrics it never saw is the best-supported one.

⚠ INTERPRETING THE SIGN: Cohen's d here is Kept-minus-Outlier. For
func_TSNR, HIGHER is better, so a POSITIVE d is the "good" direction.
For func_mean_fd / func_mean_dvars, LOWER is better, so a NEGATIVE d is
the good direction. rank_strategies() below accounts for this via
METRIC_GOOD_DIRECTION rather than ranking on |d| alone, which would
reward a strategy that separates strongly in the WRONG direction.

⚠ WHAT THIS CANNOT SETTLE: external QC metrics measure data quality
(motion, temporal SNR), not whether the retained network is
biologically correct. A strategy that converges on a shared STRUCTURED
artifact (anesthesia depth, respiratory/cardiac signature, coil or
sequence effects -- all common in animal fMRI) could still score well
here, because such runs can be both mutually consistent AND
low-motion. Treat a good score as necessary, not sufficient, and
corroborate with something independent of run-to-run consistency
(anatomical priors, cross-species convergence).
"""

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as _stats

opj = os.path.join

# +1: higher is better for Kept (expect positive Cohen's d)
# -1: lower is better for Kept  (expect negative Cohen's d)
METRIC_GOOD_DIRECTION = {
    "func_TSNR": +1,
    "func_mean_fd": -1,
    "func_mean_dvars": -1,
    "func_avg_snr_gray": +1,
    "func_cnr": +1,
}

DEFAULT_STRATEGIES = ("uniform_mean", "uniform_median", "medoid", "weighted", "kcore")


def build_reference(sri, profiles, strategy, value_col="median_r", metric="correlation",
                    weight_power=3.0, k=5, verbose=False):
    """
    Build ONE reference profile from `profiles` using the named strategy.
    Returns (ref_profile, detail) where detail carries whatever is
    strategy-specific and worth reporting (the medoid key, the weights,
    the core keys) or None.
    """
    all_keys = set(profiles.keys())
    if strategy == "uniform_mean":
        return sri.compute_group_average_profile(profiles, all_keys), None
    if strategy == "uniform_median":
        return sri.compute_group_median_profile(profiles, all_keys), None
    if strategy == "medoid":
        ref, medoid_key, _, _ = sri.compute_medoid_reference(
            profiles, value_col=value_col, metric=metric, verbose=verbose)
        return ref, dict(medoid_key=medoid_key)
    if strategy == "weighted":
        ref, weights, _, _ = sri.compute_weighted_reference(
            profiles, value_col=value_col, metric=metric,
            weight_power=weight_power, verbose=verbose)
        return ref, dict(weights=weights)
    if strategy == "kcore":
        ref, core_keys, _, _ = sri.compute_kcore_reference(
            profiles, k=k, value_col=value_col, metric=metric, verbose=verbose)
        return ref, dict(core_keys=core_keys)
    raise ValueError(f"Unknown strategy {strategy!r} -- must be one of {DEFAULT_STRATEGIES}")


def _cohens_d(a, b):
    """Cohen's d, pooled SD. Positive = group a (Kept) has the higher mean."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    na, nb = len(a), len(b)
    pooled = np.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    if pooled == 0:
        return np.nan
    return float((np.mean(a) - np.mean(b)) / pooled)


def evaluate_strategy_against_qc(df_class, df_qc_this_bids, metrics, merge_fn):
    """
    Given ONE strategy's classification, measure Kept-vs-Outlier
    separation on each external QC metric.

    merge_fn: the merge_qc_with_classification function from
    external_qc_validation (passed in rather than imported, so this
    module works even if that one isn't importable in a given session).

    Returns a list of dicts (one per metric) with cohens_d, auc,
    p_value, n_kept, n_outlier, join_level.
    """
    merged, level = merge_fn(df_qc_this_bids, df_class, "strategy-eval")
    if merged is None:
        return []

    kept_mask = ~merged["is_outlier"].astype(bool)
    outlier_mask = merged["is_outlier"].astype(bool)

    rows = []
    for metric in metrics:
        if metric not in merged.columns:
            continue
        kept_vals = pd.to_numeric(merged.loc[kept_mask, metric], errors="coerce").dropna().values
        outlier_vals = pd.to_numeric(merged.loc[outlier_mask, metric], errors="coerce").dropna().values
        if len(kept_vals) < 2 or len(outlier_vals) < 2:
            continue
        d = _cohens_d(kept_vals, outlier_vals)
        try:
            u, p = _stats.mannwhitneyu(kept_vals, outlier_vals, alternative="two-sided")
            auc = float(u / (len(kept_vals) * len(outlier_vals)))
        except ValueError:
            auc, p = np.nan, np.nan
        rows.append(dict(metric=metric, cohens_d=d, auc=auc, p_value=p,
                        n_kept=len(kept_vals), n_outlier=len(outlier_vals),
                        join_level=level))
    return rows


def compare_strategies_for_bids(sri, profiles, df_qc_this_bids, merge_fn,
                                keep_fraction=0.25, metrics=None,
                                strategies=DEFAULT_STRATEGIES,
                                value_col="median_r", metric="correlation",
                                weight_power=3.0, k=5, label="", verbose=True):
    """
    Run every strategy on ONE BIDS's profiles and evaluate each against
    that BIDS's external QC. Classification uses a FIXED keep_fraction
    for every strategy, so the only thing differing between them is the
    REFERENCE -- not how many runs get kept.

    Returns a long DataFrame: strategy, metric, cohens_d, auc, p_value,
    n_kept, n_outlier, join_level, label.
    """
    metrics = metrics or list(METRIC_GOOD_DIRECTION)
    rows = []
    for strategy in strategies:
        try:
            ref_profile, detail = build_reference(
                sri, profiles, strategy, value_col=value_col, metric=metric,
                weight_power=weight_power, k=k, verbose=verbose)
        except Exception as e:
            print(f"  [WARN] {label}/{strategy}: reference construction failed: {e}")
            continue
        if ref_profile is None:
            print(f"  [WARN] {label}/{strategy}: reference is None -- skipped")
            continue

        df_class = sri.classify_runs_by_percentile(
            profiles, ref_profile, value_col=value_col, metric=metric,
            keep_fraction=keep_fraction)

        for r in evaluate_strategy_against_qc(df_class, df_qc_this_bids, metrics, merge_fn):
            r.update(strategy=strategy, label=label)
            rows.append(r)

    return pd.DataFrame(rows)


def rank_strategies(comparison_df, metric_good_direction=None):
    """
    Rank strategies by how strongly, and in the RIGHT DIRECTION, they
    separate Kept from Outlier on the external QC metrics.

    Scoring: each (strategy, metric) contributes cohens_d * expected_sign
    -- so separating in the good direction scores positive, and
    separating strongly in the WRONG direction scores NEGATIVE rather
    than being rewarded (which ranking on |d| alone would wrongly do).
    Metrics with no known direction are skipped rather than guessed.

    Returns a DataFrame: strategy, mean_directional_d, n_metrics,
    n_significant (p<0.05 AND in the good direction), sorted best-first.
    """
    metric_good_direction = metric_good_direction or METRIC_GOOD_DIRECTION
    if comparison_df.empty:
        return pd.DataFrame(columns=["strategy", "mean_directional_d", "n_metrics", "n_significant"])

    df = comparison_df.copy()
    df["expected_sign"] = df["metric"].map(metric_good_direction)
    unknown = df["expected_sign"].isna()
    if unknown.any():
        skipped = sorted(df.loc[unknown, "metric"].unique())
        print(f"  [INFO] no known good-direction for {skipped} -- excluded from ranking "
             f"(add them to METRIC_GOOD_DIRECTION to include them)")
        df = df[~unknown]
    if df.empty:
        return pd.DataFrame(columns=["strategy", "mean_directional_d", "n_metrics", "n_significant"])

    df["directional_d"] = df["cohens_d"] * df["expected_sign"]
    df["good_and_significant"] = (df["p_value"] < 0.05) & (df["directional_d"] > 0)

    ranked = (df.groupby("strategy")
             .agg(mean_directional_d=("directional_d", "mean"),
                  n_metrics=("directional_d", "count"),
                  n_significant=("good_and_significant", "sum"))
             .reset_index()
             .sort_values("mean_directional_d", ascending=False))
    return ranked


def plot_strategy_comparison(comparison_df, metric_good_direction=None, output_path=None):
    """
    One panel per QC metric: grouped bars of Cohen's d per strategy,
    per BIDS. A line at 0, and the panel title states which direction is
    the GOOD one for that metric, so a reader can't misread a
    strongly-separating-but-wrong-direction bar as a success.
    """
    metric_good_direction = metric_good_direction or METRIC_GOOD_DIRECTION
    if comparison_df.empty:
        print("  [plot] nothing to plot")
        return None

    metrics = sorted(comparison_df["metric"].unique())
    strategies = [s for s in DEFAULT_STRATEGIES if s in set(comparison_df["strategy"])]
    labels = sorted(comparison_df["label"].unique())

    fig, axes = plt.subplots(len(metrics), 1,
                             figsize=(max(9, len(labels) * 1.6), 4 * len(metrics)),
                             squeeze=False)
    axes = axes[:, 0]
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, max(1, len(strategies))))

    for ax, met in zip(axes, metrics):
        sub = comparison_df[comparison_df["metric"] == met]
        x = np.arange(len(labels))
        width = 0.8 / max(1, len(strategies))
        for si, (strategy, color) in enumerate(zip(strategies, colors)):
            vals = []
            for lab in labels:
                cell = sub[(sub["label"] == lab) & (sub["strategy"] == strategy)]["cohens_d"]
                vals.append(cell.mean() if len(cell) else np.nan)
            ax.bar(x + si * width - 0.4 + width / 2, vals, width, label=strategy, color=color)
        ax.axhline(0, color="#888888", lw=0.9)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("Cohen's d (Kept - Outlier)", fontsize=8)
        direction = metric_good_direction.get(met)
        dir_txt = ("higher d = better" if direction == +1
                   else "LOWER (negative) d = better" if direction == -1
                   else "direction unknown")
        ax.set_title(f"{met}  ({dir_txt})", fontsize=10, fontweight="bold")
        ax.legend(fontsize=7, ncol=max(1, len(strategies)))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def run_full_strategy_comparison(sri, launcher, df_qc_by_label, merge_fn,
                                 keep_fraction=0.25, metrics=None,
                                 strategies=DEFAULT_STRATEGIES,
                                 weight_power=3.0, k=5, output_dir=None,
                                 bids_configs=None):
    """
    Whole-pipeline driver: for every BIDS_CONFIGS entry, compute its run
    profiles once, compare every reference strategy on it, then pool and
    rank across BIDS.

    sri: the SBA_Region_Profile_Index module.
    launcher: the Figure5_SBA_KeeporNot module (for its CONFIG constants).
    df_qc_by_label: {label: df_qc slice for that BIDS} -- same structure
    external_qc_validation.compare_kept_outlier_on_external_metrics takes.
    merge_fn: external_qc_validation.merge_qc_with_classification.

    Returns (comparison_df, ranked_df).
    """
    bids_configs = bids_configs if bids_configs is not None else launcher.BIDS_CONFIGS
    frames = []

    for cfg in bids_configs:
        label = cfg["label"]
        seed = cfg["seeds"][0]
        print(f"\n{'=' * 70}\n  {label} / {seed}\n{'=' * 70}")

        df_qc_this = df_qc_by_label.get(label)
        if df_qc_this is None or df_qc_this.empty:
            print(f"  [SKIP] no external QC rows for {label}")
            continue

        try:
            label_df = sri.parse_itk_snap_label_file(cfg["label_path"])
            profiles, _ = sri.compute_all_run_profiles(
                cfg["species"], cfg["bids_dir"], label_df, launcher.CLASSIFICATION_PATH,
                category=launcher.CATEGORY, seeds=[seed], atlas_level=launcher.ATLAS_LEVEL)
            if not profiles:
                print(f"  [SKIP] no run profiles computed")
                continue

            df_bids = compare_strategies_for_bids(
                sri, profiles, df_qc_this, merge_fn,
                keep_fraction=keep_fraction, metrics=metrics, strategies=strategies,
                value_col=launcher.VALUE_COL, metric=launcher.METRIC,
                weight_power=weight_power, k=k, label=label, verbose=True)
            if not df_bids.empty:
                frames.append(df_bids)
        except Exception as e:
            print(f"  [ERROR] {label}: {e}")

    if not frames:
        print("\n[WARN] no BIDS produced comparable results")
        return pd.DataFrame(), pd.DataFrame()

    comparison_df = pd.concat(frames, ignore_index=True)
    ranked_df = rank_strategies(comparison_df)

    print(f"\n\n{'=' * 70}\n  STRATEGY RANKING (across all BIDS)\n{'=' * 70}")
    print(ranked_df.to_string(index=False))
    print("\n(mean_directional_d: positive = separates in the GOOD direction on average;\n"
          " n_significant: (strategy, BIDS, metric) cells with p<0.05 AND good direction)")

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        excel_path = opj(output_dir, "reference_strategy_comparison.xlsx")
        with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
            ranked_df.to_excel(writer, sheet_name="ranking", index=False)
            comparison_df.to_excel(writer, sheet_name="per_bids_per_metric", index=False)
        print(f"\n  [saved] {excel_path}")
        plot_strategy_comparison(
            comparison_df, output_path=opj(output_dir, "reference_strategy_comparison.png"))

    return comparison_df, ranked_df
