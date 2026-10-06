"""
external_qc_validation.py
────────────────────────────
Validates the Kept/Outlier classification against EXTERNAL, independent
QC metrics (motion, tSNR, DVARS, ...) — the metrics computed by EDNiX's
own QC pipeline (collect_multi_species(..., extract=('qc',)), same
source as Suppl_Figure_4-5-6_QC.py's df_qc), which are NEVER used
anywhere in building df_class (that comes from SBA correlation-to-
reference). This breaks the circularity a same-metric comparison would
have: "is Kept better than Outlier" here means "on a metric that had
NOTHING to do with how Kept/Outlier was decided."

⚠ You collect df_qc YOURSELF first (same as Suppl_Figure_4-5-6_QC.py's
Step 1-2, via collect_multi_species(..., extract=('qc',))) and pass it
in — this module never calls collect_multi_species itself, since its
required species_config/atlas_label_paths setup is specific to your
own script and every kwarg it needs hasn't been independently verified
here.

⚠ UNCERTAIN GRANULARITY: whether df_qc has one row per RUN, per SESSION,
or per SUBJECT depends on how the underlying QC pipeline computed these
metrics (motion/tSNR are usually per-run in fMRI, but
Suppl_Figure_4-5-6_QC.py's own session-averaging step for df_qc_plot
suggests the RAW df_qc might already only reach session granularity for
some metrics — its own groupby explicitly doesn't even check for a
"run" column). merge_qc_with_classification() below tries progressively
coarser (subject, session, run) -> (subject, session) -> (subject)
joins and reports which one actually matched real rows, rather than
assuming — check that printed join_level before trusting the results:
a "subject"-only join means every run of that subject shares the SAME
external QC value, which is a much weaker test than a true per-run join.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as _stats

DEFAULT_EXTERNAL_METRICS = ["func_mean_fd", "func_TSNR", "func_mean_dvars"]


def load_qc_from_excel(path, sheet_name=0):
    """
    Load the raw per-individual QC spreadsheet saved by
    Suppl_Figure_4-5-6_QC.py (df_qc.to_excel(path, index=False), added
    right after df_qc is collected — BEFORE any session averaging, so
    this reload has the same per-run/session granularity the original
    collect_multi_species(..., extract=('qc',)) call produced).

    Reads subject/session/run with dtype=str FORCED AT READ TIME (not
    cast afterward) — confirmed real gap: a plain read followed by
    .astype(str) is NOT enough, because pandas' Excel engine already
    numeric-infers a column like subject="01" DURING parsing (openpyxl
    doesn't know it should stay a string), so the leading zero is gone
    before any post-hoc cast could ever see it — the same "01" -> 1.0
    loss confirmed elsewhere in this pipeline, but this time surviving
    even the .astype(str) fix that was sufficient for CSV/other reads.
    keep_default_na=False additionally stops "NA" from being read as a
    real NaN (same reasoning as SBA_Region_Profile_Index's
    load_all_results_from_disk).
    """
    header_df = pd.read_excel(path, sheet_name=sheet_name, nrows=0)
    id_cols = [c for c in ("subject", "session", "run") if c in header_df.columns]
    dtype_map = {c: str for c in id_cols}

    df_qc = pd.read_excel(path, sheet_name=sheet_name, dtype=dtype_map,
                          keep_default_na=False, na_values=[])
    print(f"  Loaded {len(df_qc)} row(s), {len(df_qc.columns)} column(s) from {path}")
    return df_qc


def _cohens_d(a, b):
    """Cohen's d, pooled SD. Positive = group a (Kept) has the higher mean."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    na, nb = len(a), len(b)
    pooled_sd = np.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    if pooled_sd == 0:
        return np.nan
    return float((np.mean(a) - np.mean(b)) / pooled_sd)


def merge_qc_with_classification(df_qc, df_class, bids_label):
    """
    Merge one BIDS's df_class (subject/session/run/is_outlier) with the
    EXTERNAL df_qc (subject/session[/run], func_mean_fd/func_TSNR/...),
    trying progressively coarser join keys and reporting which one
    actually produced matches — since df_qc's real granularity isn't
    confirmed (see module docstring).

    df_qc is expected to already be filtered/restricted to this ONE
    BIDS's rows (e.g. df_qc[df_qc["bids_dir"].str.contains(bids_name)])
    — this function does not do that filtering itself, since the exact
    bids_dir string format in df_qc (full path? basename?) isn't
    confirmed either.

    Returns (merged_df, join_level) where join_level is one of
    "subject+session+run", "subject+session", "subject", or (None, None)
    if nothing matched at any granularity.
    """
    df_class = df_class.copy()
    for col in ("subject", "session", "run"):
        if col in df_class.columns:
            df_class[col] = df_class[col].astype(str)

    df_qc = df_qc.copy()
    for col in ("subject", "session", "run"):
        if col in df_qc.columns:
            df_qc[col] = df_qc[col].astype(str)

    candidate_keys = []
    if "run" in df_qc.columns and "run" in df_class.columns:
        candidate_keys.append(["subject", "session", "run"])
    if "session" in df_qc.columns and "session" in df_class.columns:
        candidate_keys.append(["subject", "session"])
    candidate_keys.append(["subject"])

    for keys in candidate_keys:
        if not all(k in df_class.columns and k in df_qc.columns for k in keys):
            continue
        merged = df_class.merge(df_qc, on=keys, how="inner", suffixes=("", "_qc"))
        if not merged.empty:
            level = "+".join(keys)
            print(f"  [{bids_label}] matched {len(merged)}/{len(df_class)} df_class rows "
                 f"to external QC at join level: {level}")
            return merged, level

    print(f"  [WARN] {bids_label}: no join key (tried "
         f"{['+'.join(k) for k in candidate_keys]}) matched anything between "
         f"df_class and df_qc — check subject/session/run value formats on both sides "
         f"(e.g. zero-padding, or df_qc's subject ids not matching df_class's at all).")
    return None, None


def compare_kept_outlier_on_external_metrics(all_results, df_qc_by_label,
                                              metrics=None, verbose=True):
    """
    For every BIDS in all_results, merge its df_class with the matching
    slice of df_qc (df_qc_by_label[label]), then compare Kept vs Outlier
    on each external metric with Mann-Whitney U + Cohen's d — the
    non-circular validation the classification needs: does Outlier show
    worse motion/tSNR/DVARS on a metric NEVER used to build the
    classification?

    df_qc_by_label: {label: df_qc_slice_for_this_bids} — YOU split the
    full df_qc by BIDS yourself (this function doesn't assume how
    "bids_dir" is spelled in df_qc), e.g.:
        df_qc_by_label = {label: df_qc[df_qc["bids_dir"].astype(str).str.contains(bids_name)]
                          for label, bids_name in ...}

    Returns a long DataFrame: one row per (label, metric) with n_kept,
    n_outlier, kept_mean, outlier_mean, cohens_d, auc, p_value,
    join_level — ready for a table or the ranking figure below.
    """
    metrics = metrics or DEFAULT_EXTERNAL_METRICS
    rows = []
    for (label, seed), result in all_results.items():
        df_qc_this = df_qc_by_label.get(label)
        if df_qc_this is None or df_qc_this.empty:
            if verbose:
                print(f"  [SKIP] {label}: no external QC data given")
            continue
        merged, level = merge_qc_with_classification(df_qc_this, result["df_class"], label)
        if merged is None:
            continue

        kept_mask = ~merged["is_outlier"].astype(bool)
        outlier_mask = merged["is_outlier"].astype(bool)
        for metric in metrics:
            if metric not in merged.columns:
                continue
            kept_vals = pd.to_numeric(merged.loc[kept_mask, metric], errors="coerce").dropna().values
            outlier_vals = pd.to_numeric(merged.loc[outlier_mask, metric], errors="coerce").dropna().values
            if len(kept_vals) < 2 or len(outlier_vals) < 2:
                continue
            d = _cohens_d(kept_vals, outlier_vals)
            try:
                u_stat, p = _stats.mannwhitneyu(kept_vals, outlier_vals, alternative="two-sided")
                auc = float(u_stat / (len(kept_vals) * len(outlier_vals)))
            except ValueError:
                auc, p = np.nan, np.nan
            rows.append(dict(label=label, metric=metric, join_level=level,
                            n_kept=len(kept_vals), n_outlier=len(outlier_vals),
                            kept_mean=float(np.mean(kept_vals)),
                            outlier_mean=float(np.mean(outlier_vals)),
                            cohens_d=d, auc=auc, p_value=p))

    result_df = pd.DataFrame(rows)
    if verbose and not result_df.empty:
        print()
        print(result_df.sort_values("p_value").to_string(index=False))
    return result_df


def plot_external_qc_comparison(comparison_df, metrics=None, output_path=None):
    """
    One panel per external metric: for every BIDS, a bar showing
    Cohen's d (Kept vs Outlier). Sign shows direction (positive = Kept
    has the higher value — remember to flip interpretation for metrics
    where LOWER is better, e.g. func_mean_fd/func_mean_dvars: for those,
    a NEGATIVE bar here — Outlier higher — is the "good" direction).
    Bars annotated with significance stars from p_value.
    """
    metrics = metrics or sorted(comparison_df["metric"].unique())
    n = len(metrics)
    if n == 0:
        print("  [plot] nothing to plot")
        return None
    fig, axes = plt.subplots(1, n, figsize=(4.5 * n, 4.5), squeeze=False)
    axes = axes[0]

    for ax, metric in zip(axes, metrics):
        sub = comparison_df[comparison_df["metric"] == metric].sort_values("label")
        if sub.empty:
            ax.axis("off")
            continue
        colors = ["#0072B2" if d > 0 else "#D55E00" for d in sub["cohens_d"]]
        ax.bar(sub["label"], sub["cohens_d"], color=colors, edgecolor="white")
        for i, (d, p) in enumerate(zip(sub["cohens_d"], sub["p_value"])):
            if not np.isfinite(d) or not np.isfinite(p):
                continue
            sig = "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else ""
            ax.text(i, d + (0.05 if d >= 0 else -0.05), sig, ha="center",
                   va="bottom" if d >= 0 else "top", fontsize=10, fontweight="bold")
        ax.axhline(0, color="#888888", lw=0.8)
        ax.set_title(metric, fontsize=10, fontweight="bold")
        ax.set_ylabel("Cohen's d (Kept − Outlier)", fontsize=8)
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def compare_reference_strategies(all_results_individual_ref, all_results_group_ref,
                                 df_qc_by_label, metrics=None, verbose=True):
    """
    THE reference-strategy comparison discussed in this conversation:
    for each BIDS, run compare_kept_outlier_on_external_metrics TWICE —
    once for a classification built with an individual reference
    subject, once for the SAME BIDS classified against the group-average
    reference — and put the two |Cohen's d| side by side per metric.
    Whichever strategy produces the LARGER |Cohen's d| (stronger
    separation on a metric NEVER used to classify) is the
    better-justified reference choice for that BIDS — this is the
    non-circular comparison the "individual vs group-average reference"
    question needs.

    all_results_individual_ref / all_results_group_ref: two separate
    all_results-shaped dicts for the SAME BIDS_CONFIGS, one built with
    each BIDS's cfg["reference_subject"] set as usual, one with it set
    to None (group average) — e.g. run launch() twice with a temporarily
    modified BIDS_CONFIGS (reference_subject=None for the second run),
    saving each to a separate OUTPUT_DIR, then load_all_results_from_disk()
    against each directory.

    Returns a long DataFrame: one row per (label, metric, strategy) with
    the same columns as compare_kept_outlier_on_external_metrics, plus
    a "strategy" column ("individual" / "group_average").
    """
    df_ind = compare_kept_outlier_on_external_metrics(
        all_results_individual_ref, df_qc_by_label, metrics=metrics, verbose=False)
    df_ind["strategy"] = "individual"
    df_grp = compare_kept_outlier_on_external_metrics(
        all_results_group_ref, df_qc_by_label, metrics=metrics, verbose=False)
    df_grp["strategy"] = "group_average"

    combined = pd.concat([df_ind, df_grp], ignore_index=True)
    if verbose and not combined.empty:
        summary = combined.pivot_table(index=["label", "metric"], columns="strategy",
                                       values="cohens_d")
        print(summary.to_string())
    return combined


def plot_reference_strategy_comparison(combined_df, metrics=None, output_path=None):
    """
    Grouped bar chart: for each (BIDS, metric), |Cohen's d| for the
    individual-reference classification next to the group-average
    classification's own |Cohen's d| — the taller bar in each pair is
    the better-separated (hence better-justified) reference choice for
    that BIDS/metric.
    """
    metrics = metrics or sorted(combined_df["metric"].unique())
    n = len(metrics)
    if n == 0:
        print("  [plot] nothing to plot")
        return None
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4.5), squeeze=False)
    axes = axes[0]

    for ax, metric in zip(axes, metrics):
        sub = combined_df[combined_df["metric"] == metric]
        labels = sorted(sub["label"].unique())
        x = np.arange(len(labels))
        width = 0.35
        for offset, strategy, color in [(-width / 2, "individual", "#0072B2"),
                                        (width / 2, "group_average", "#D55E00")]:
            vals = [sub.loc[(sub["label"] == l) & (sub["strategy"] == strategy), "cohens_d"]
                   .abs().mean() for l in labels]
            ax.bar(x + offset, vals, width, label=strategy, color=color, edgecolor="white")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("|Cohen's d| (Kept vs Outlier)", fontsize=8)
        ax.set_title(metric, fontsize=10, fontweight="bold")
        ax.legend(fontsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig