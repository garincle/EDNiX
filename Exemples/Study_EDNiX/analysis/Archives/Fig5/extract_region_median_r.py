"""
extract_region_median_r.py
─────────────────────────────
Lightweight, per-run extraction of median_r for SPECIFIC named regions
(e.g. Somatosensory, Prefrontal), Kept vs Outlier — WITHOUT needing a
full compute_all_run_profiles() relaunch (which recomputes every
region's profile for every run in the BIDS). Each run's own raw SBA
correlation map + atlas volume are located directly and
compute_region_profile is called with label_df already restricted to
just the matched region(s), so this stays cheap even for a full BIDS.

⚠ Step 1 first: region names in the atlas label file are NOT guessed —
find_region_names() looks them up for real and prints them, so you
confirm the actual names ("Somatosensory", "SS", "S1", whatever this
species' EDNIxCSCLR label file really calls it) before extracting
anything.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# =============================================================================
# 1) Find the REAL region names — no guessing
# =============================================================================

def find_region_names(label_df, pattern):
    """
    Look up region_name values in label_df (from
    sri.parse_itk_snap_label_file) that CASE-INSENSITIVELY contain
    `pattern` — the real names, not a guess. Prints and returns the
    sorted unique matches.
    """
    matched = sorted(label_df[label_df["region_name"].str.contains(
        pattern, case=False, na=False)]["region_name"].unique().tolist())
    print(f"  '{pattern}' matched {len(matched)} region(s): {matched}")
    return matched


def find_region_names_by_species(sri, bids_configs, pattern):
    """
    Same lookup as find_region_names, but run for EVERY BIDS_CONFIGS
    entry's own label_path — since different species' EDNIxCSCLR label
    files can name the same anatomical concept differently (or not at
    all). Returns {label: [matched region names]}.
    """
    out = {}
    for cfg in bids_configs:
        label_df = sri.parse_itk_snap_label_file(cfg["label_path"])
        print(f"\n{cfg['label']}:")
        out[cfg["label"]] = find_region_names(label_df, pattern)
    return out


# =============================================================================
# 2) Per-run median_r extraction, restricted to the matched region(s)
# =============================================================================

def _find_run_map_path(sri, bids_dir, subject, session, run, seed):
    """
    Locate ONE run's raw SBA correlation map — same building blocks
    _resolve_reference_map already uses (find_raw_corr_maps_for_run +
    _extract_run_id-based token matching), just for an arbitrary
    (subject, session, run) instead of a configured reference subject.
    Returns None if nothing matches this exact run token.
    """
    candidates = sri.find_raw_corr_maps_for_run(bids_dir, subject, session, seed)
    if not candidates:
        return None
    run_str = str(run)
    exact = [c for c in candidates if sri._extract_run_id(Path(c).name) == run_str]
    if exact:
        return sorted(exact)[0]
    # same RL-/LR- flexible fallback confirmed necessary elsewhere in this pipeline
    def _digits(tok):
        return re.sub(r"^(RL|LR)-", "", tok)
    flexible = [c for c in candidates
               if _digits(sri._extract_run_id(Path(c).name)) == _digits(run_str)]
    return sorted(flexible)[0] if flexible else None


def extract_region_median_r_per_run(sri, bids_dir, label_df, df_class, region_patterns,
                                    atlas_level=None, verbose=True):
    """
    For each run in df_class, locate its OWN raw SBA correlation map +
    atlas volume directly, and compute median_r restricted to just the
    region(s) matching each pattern in region_patterns — without a full
    compute_all_run_profiles() pass.

    region_patterns: list of (label, substring) tuples, e.g.
    [("Somatosensory", "Somatosensory"), ("Prefrontal", "Prefrontal")] —
    substring matched case-insensitively. Run find_region_names first to
    confirm these substrings actually match something real for this
    species' label file.

    Returns a long DataFrame: subject, session, run, seed, is_outlier,
    pattern, median_r (averaged across every region matched by that
    pattern, for that run), n_regions_matched. Runs where the map/atlas
    couldn't be found, or where compute_region_profile produced no
    usable value, are skipped (counted and reported at the end, not
    per-row, to avoid flooding the console for a full BIDS).
    """
    # Restrict label_df once per pattern -- computing on a 2-row label_df
    # instead of the full ~276 is what keeps this cheap.
    matched_by_pattern = {}
    for pat_label, pat_substr in region_patterns:
        sub_df = label_df[label_df["region_name"].str.contains(pat_substr, case=False, na=False)]
        matched_by_pattern[pat_label] = sub_df
        if sub_df.empty and verbose:
            print(f"  [WARN] pattern '{pat_label}' ('{pat_substr}') matched NO region "
                 f"in this label_df -- run find_region_names first to check the real names")

    rows = []
    n_no_map, n_no_atlas, n_profile_fail = 0, 0, 0

    for row in df_class.itertuples():
        subject, session, run, seed = (str(row.subject), str(row.session),
                                       str(row.run), str(row.seed))
        map_path = _find_run_map_path(sri, bids_dir, subject, session, run, seed)
        if map_path is None:
            n_no_map += 1
            continue
        atlas_vol_path = sri.find_atlas_volume_path(bids_dir, subject, session)
        if atlas_vol_path is None:
            n_no_atlas += 1
            continue

        for pat_label, sub_label_df in matched_by_pattern.items():
            if sub_label_df.empty:
                continue
            try:
                prof = sri.compute_region_profile(map_path, atlas_vol_path, sub_label_df,
                                                  atlas_level=atlas_level)
            except Exception as e:
                n_profile_fail += 1
                if verbose:
                    print(f"  [WARN] {subject}/{session}/{run} — {pat_label}: "
                         f"compute_region_profile failed: {e}")
                continue
            vals = prof["median_r"].dropna() if "median_r" in prof.columns else pd.Series(dtype=float)
            if vals.empty:
                continue
            rows.append(dict(subject=subject, session=session, run=run, seed=seed,
                            is_outlier=bool(row.is_outlier), pattern=pat_label,
                            median_r=float(vals.mean()), n_regions_matched=len(vals)))

    if verbose:
        print(f"  Extracted {len(rows)} (run, pattern) value(s) from {len(df_class)} run(s) "
             f"-- {n_no_map} map(s) not found, {n_no_atlas} atlas(es) not found, "
             f"{n_profile_fail} profile computation failure(s)")

    return pd.DataFrame(rows)


# =============================================================================
# 3) Plot: Kept vs Outlier boxplot per pattern
# =============================================================================

def plot_median_r_kept_vs_outlier(long_df, title=None, output_path=None, dpi=300, jitter=0.08,
                                  random_seed=0):
    """
    Boxplot(s) of median_r from extract_region_median_r_per_run's long
    DataFrame, one panel per pattern, Kept vs Outlier, with individual
    run points overlaid.
    """
    if long_df.empty:
        print("  [plot] nothing to plot -- long_df is empty")
        return None

    patterns = list(long_df["pattern"].unique())
    rng = np.random.default_rng(random_seed)
    fig, axes = plt.subplots(1, len(patterns), figsize=(4.5 * len(patterns), 5), squeeze=False)
    axes = axes[0]

    for ax, pattern in zip(axes, patterns):
        sub = long_df[long_df["pattern"] == pattern]
        kept_vals = sub.loc[~sub["is_outlier"], "median_r"].dropna().tolist()
        outlier_vals = sub.loc[sub["is_outlier"], "median_r"].dropna().tolist()

        data = [kept_vals, outlier_vals]
        positions = [1, 2]
        bp = ax.boxplot(data, positions=positions, widths=0.5, showfliers=False, patch_artist=True)
        for patch, color in zip(bp["boxes"], ["#009E73", "#D55E00"]):
            patch.set_facecolor(color)
            patch.set_alpha(0.25)
            patch.set_edgecolor(color)
        for median in bp["medians"]:
            median.set_color("#333333")

        for pos, vals, color in zip(positions, data, ["#009E73", "#D55E00"]):
            if not vals:
                continue
            x_jitter = pos + rng.uniform(-jitter, jitter, size=len(vals))
            ax.scatter(x_jitter, vals, color=color, alpha=0.6, s=22,
                      edgecolors="white", lw=0.4, zorder=3)

        ax.set_xticks(positions)
        ax.set_xticklabels([f"Kept\n(n={len(kept_vals)})", f"Outlier\n(n={len(outlier_vals)})"],
                           fontsize=8)
        ax.set_ylabel("median_r", fontsize=9)
        ax.set_title(pattern, fontsize=9, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle(title or "median_r by region, Kept vs Outlier", fontsize=11)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig
