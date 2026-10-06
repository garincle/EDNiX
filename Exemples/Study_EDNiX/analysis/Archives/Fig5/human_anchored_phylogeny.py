"""
human_anchored_phylogeny.py
──────────────────────────────
Tests the hypothesis: if runs in each non-human species are selected for
maximum resemblance to a HUMAN reference profile, does the RESIDUAL
distance to human still track evolutionary divergence time?

The reasoning: even when explicitly picking the most human-like runs
available, a species-specific FLOOR should remain -- a mouse brain
cannot resemble a human beyond some biological limit no matter which
run is chosen. If that residual floor correlates with divergence time,
that is informative about conserved network organization.

Pipeline:
  1. Build a HUMAN reference profile from the human BIDS (two methods
     compared side by side: "kcore" and "uniform_median").
  2. Find the lvl2 regions common to EVERY species' atlas.
  3. Compute each run's profile in every BIDS, restricted to those
     common regions.
  4. Per BIDS, select the N runs whose profiles correlate best with the
     human reference.
  5. Run the group 3dttest (via rerun_balanced_by_run) on those selected
     runs, producing a group map per BIDS.
  6. Measure each species' group-map distance to the human group map, in
     the same common regions.
  7. Plot that distance against divergence time from human, with a
     correlation test.

⚠ CIRCULARITY, STATED PLAINLY: runs are SELECTED for human-likeness and
then MEASURED for human-likeness. Every species will therefore look
more human-like than it would unselected, and the absolute distances are
not interpretable on their own. What the phylogenetic test asks is
narrower: whether the RESIDUAL differences that survive this selection
still rank species by divergence time. A positive result is suggestive,
not confirmatory, and needs one of the held-out controls (select on one
seed/region set, test on another) before it can carry a strong claim.

⚠ SAMPLE-SIZE BIAS, NOT CORRECTED HERE: BIDS differ enormously in run
count (e.g. ~195 for Marmoset vs ~25 for Dog). With more runs to choose
from, a species gets more chances to contain a human-like run BY
CHANCE, so it can appear closer to human for reasons unrelated to
biology. subsample_n below (off by default) equalizes the pool size
before selection; use it before drawing any phylogenetic conclusion.

⚠ HUMAN IS EXCLUDED from the phylogenetic correlation: it sits at
distance ~0 from itself by construction, so including it would inflate
the correlation artificially. That leaves few points (~6), enough for a
trend, not for a well-powered test.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as _stats

opj = os.path.join

# Divergence time from human, in millions of years (TimeTree-style
# estimates).
# ⚠ NOT INDEPENDENTLY VERIFIED -- these vary by source, calibration and
# tree version. Check them against TimeTree (timetree.org) for your own
# species set before publishing anything that rests on them.
# ⚠ Note the resolution limit this implies: rodents (~87) and carnivores
# (~94) diverge at nearly the same depth, so this design mostly
# separates primates from non-primates rather than resolving a fine
# gradient.
DIVERGENCE_TIME_MYA = {
    "Human": 0.0,
    "Macaque": 29.0,
    "Marmoset": 43.0,
    "Mouselemur": 74.0,
    "Rat": 87.0,
    "Mouse": 87.0,
    "Dog": 94.0,
}

SPECIES_COLORS = {
    "Human": "#D35400", "Macaque": "#F5B041", "Marmoset": "#F4D03F",
    "Mouselemur": "#58D68D", "Rat": "#2E86C1", "Mouse": "#1A5276",
    "Dog": "#8E44AD",
}


# =============================================================================
# 1) Human reference
# =============================================================================

def build_human_reference(sri, launcher, human_label, method="kcore", k=5,
                          weight_power=3.0, verbose=True):
    """
    Build the human reference profile from ONE human BIDS, using either
    "kcore" (medoid + k-1 nearest neighbours) or "uniform_median" (plain
    median over every human run).

    Returns (ref_profile, profiles, detail) -- `profiles` is that BIDS's
    full per-run profile dict (reused later so it isn't recomputed).
    """
    cfg = next(c for c in launcher.BIDS_CONFIGS if c["label"] == human_label)
    seed = cfg["seeds"][0]
    label_df = sri.parse_itk_snap_label_file(cfg["label_path"])
    profiles, _ = sri.compute_all_run_profiles(
        cfg["species"], cfg["bids_dir"], label_df, launcher.CLASSIFICATION_PATH,
        category=launcher.CATEGORY, seeds=[seed], atlas_level=launcher.ATLAS_LEVEL)
    if not profiles:
        raise ValueError(f"{human_label}: no run profiles computed -- cannot build a reference.")

    if method == "kcore":
        ref, core_keys, _, _ = sri.compute_kcore_reference(
            profiles, k=k, value_col=launcher.VALUE_COL, metric=launcher.METRIC,
            verbose=verbose)
        detail = dict(core_keys=core_keys)
    elif method == "uniform_median":
        ref = sri.compute_group_median_profile(profiles, set(profiles.keys()))
        detail = dict(n_runs=len(profiles))
    else:
        raise ValueError(f"Unknown method {method!r} -- use 'kcore' or 'uniform_median'.")

    if verbose:
        print(f"  Human reference ({method}) built from {human_label} "
             f"({len(profiles)} run(s) available)")
    return ref, profiles, detail


def build_reference_for_profiles(sri, profiles, method="kcore", k=5, value_col="median_r",
                                 metric="correlation", verbose=False):
    """
    Build a reference profile from an ALREADY-COMPUTED profiles dict,
    using the SAME method ("kcore" or "uniform_median") for whichever
    BIDS it's called on -- the human reference AND every species' own
    reference are built this identical way, so they are directly
    comparable (same construction, same method, just different input
    runs).
    """
    if method == "kcore":
        ref, core_keys, _, _ = sri.compute_kcore_reference(
            profiles, k=k, value_col=value_col, metric=metric, verbose=verbose)
        return ref, dict(core_keys=core_keys)
    if method == "uniform_median":
        ref = sri.compute_group_median_profile(profiles, set(profiles.keys()))
        return ref, dict(n_runs=len(profiles))
    raise ValueError(f"Unknown method {method!r} -- use 'kcore' or 'uniform_median'.")


def build_run_correlation_table(sri, profiles_by_label, own_ref_by_label, human_ref,
                                common_base_names, selections=None,
                                value_col="median_r", metric="correlation"):
    """
    For EVERY run of EVERY BIDS, compute its Pearson r to (a) its OWN
    species' reference and (b) the human reference, both restricted to
    the fixed common-region list. This is the actual metric selection is
    based on (b) alongside the complementary "how typical for its own
    species" question (a) -- put side by side so both can be inspected
    directly rather than only read off a plot.

    selections: optional {label: [key, ...]} (e.g.
    results[method]["selections"]) -- adds an is_selected column marking
    which runs the human-anchored selection actually kept.

    Returns a long DataFrame: label, species, subject, session, run, seed,
    r_to_species_ref, r_to_human_ref, is_selected.
    """
    selections = selections or {}
    rows = []
    for label, profiles in profiles_by_label.items():
        own_ref = own_ref_by_label.get(label)
        selected_str = set(selections.get(label, []))
        for key, prof in profiles.items():
            str_key = tuple(str(x) for x in key)
            r_species = r_human = np.nan
            if own_ref is not None:
                try:
                    _, r_species, _ = sri.profile_distance(
                        prof, own_ref, value_col=value_col, metric=metric)
                except Exception:
                    pass
            try:
                _, r_human, _ = sri.profile_distance(
                    prof, human_ref, value_col=value_col, metric=metric)
            except Exception:
                pass
            rows.append(dict(label=label, species=None, subject=key[0], session=key[1],
                            run=key[2], seed=key[3], r_to_species_ref=r_species,
                            r_to_human_ref=r_human,
                            is_selected=str_key in selected_str))
    return pd.DataFrame(rows)


# =============================================================================
# 2) Common regions + restriction
# =============================================================================

def get_common_lvl_regions(sri, launcher, bids_configs=None, verbose=True):
    """
    Base region names shared by EVERY BIDS_CONFIGS entry's atlas label
    file (via sri.find_common_base_regions).

    ⚠ Name-intersection based -- see find_common_base_regions' own
    docstring: if your Atlas_nomenclature module exposes a cross-species
    homology/alias table, that would be more correct than bare name
    matching, since homologous regions can be named differently across
    species atlases.
    """
    bids_configs = bids_configs if bids_configs is not None else launcher.BIDS_CONFIGS
    label_dfs = {}
    for cfg in bids_configs:
        try:
            label_dfs[cfg["label"]] = sri.parse_itk_snap_label_file(cfg["label_path"])
        except Exception as e:
            print(f"  [WARN] {cfg['label']}: could not parse label file: {e}")
    common = sri.find_common_base_regions(label_dfs)
    if verbose:
        print(f"  {len(common)} common base region(s) across {len(label_dfs)} BIDS")
    return common


def restrict_profile(profile, common_base_names):
    """Keep only rows whose base_name is in common_base_names."""
    if profile is None or "base_name" not in profile.columns:
        return profile
    mask = profile["base_name"].str.lower().str.strip().isin(common_base_names)
    return profile[mask].copy()


# =============================================================================
# 3) Select the runs closest to the human reference
# =============================================================================

def select_runs_closest_to_reference(sri, profiles, ref_profile, n_select=5,
                                     common_base_names=None, value_col="median_r",
                                     metric="correlation", subsample_n=None,
                                     random_seed=42, verbose=True):
    """
    Rank every run in `profiles` by correlation to `ref_profile`
    (restricted to common_base_names if given) and return the top
    n_select.

    subsample_n: if given, randomly draw this many runs from `profiles`
    BEFORE ranking -- the sample-size control described in the module
    docstring. Without it, a BIDS with 195 runs has far more chances to
    contain an unusually human-like run than one with 25, so the
    resulting "closeness to human" partly reflects run count rather than
    biology. None (default) uses every run, which is fine for a first
    exploratory look but not for the phylogenetic claim.

    Returns (selected_keys, ranking_df) where ranking_df has every run's
    pearson_r to the reference, sorted best-first.
    """
    ref = restrict_profile(ref_profile, common_base_names) if common_base_names else ref_profile

    keys = list(profiles.keys())
    if subsample_n is not None and subsample_n < len(keys):
        rng = np.random.default_rng(random_seed)
        idx = rng.choice(len(keys), size=subsample_n, replace=False)
        keys = [keys[i] for i in sorted(idx)]
        if verbose:
            print(f"  [subsample] drew {subsample_n} of {len(profiles)} run(s) before ranking")

    rows = []
    for key in keys:
        prof = restrict_profile(profiles[key], common_base_names) if common_base_names else profiles[key]
        try:
            _, r, n = sri.profile_distance(prof, ref, value_col=value_col, metric=metric)
        except Exception:
            continue
        rows.append(dict(subject=key[0], session=key[1], run=key[2], seed=key[3],
                        pearson_r=r, n_regions=n, _key=key))

    ranking_df = pd.DataFrame(rows).sort_values("pearson_r", ascending=False).reset_index(drop=True)
    if ranking_df.empty:
        return [], ranking_df

    selected_keys = ranking_df.head(n_select)["_key"].tolist()
    if verbose:
        top_r = ranking_df.head(n_select)["pearson_r"]
        print(f"  Selected {len(selected_keys)} run(s), r to human reference: "
             f"{top_r.min():.4f} to {top_r.max():.4f}")
    return selected_keys, ranking_df


def selected_keys_to_df_class(profiles, selected_keys):
    """
    Build a df_class (the shape rerun_balanced_by_run expects) marking
    selected_keys as kept and everything else as outlier.
    """
    selected_str = {tuple(str(x) for x in k) for k in selected_keys}
    rows = []
    for key in profiles:
        str_key = tuple(str(x) for x in key)
        rows.append(dict(subject=key[0], session=key[1], run=key[2], seed=key[3],
                        is_outlier=str_key not in selected_str))
    return pd.DataFrame(rows)


# =============================================================================
# 4) Group maps + distance to human
# =============================================================================

def run_group_analysis_for_selection(sri, cfg, seed, df_class, output_suffix,
                                     render_surface=False, verbose=True):
    """
    Run the group 3dttest on one BIDS's selected runs via
    rerun_balanced_by_run. render_surface defaults to False here since
    the phylogeny test only needs the volumetric group map; set True for
    the final figures.

    Returns True on success, False if the rerun raised.
    """
    try:
        sri.rerun_balanced_by_run(
            df_class, cfg["species"], cfg["bids_dir"], seed=seed,
            output_suffix_kept=output_suffix,
            output_suffix_outlier=f"{output_suffix}_NotSelected",
            exclude_lateralized=cfg.get("exclude_lateralized", True),
            render_surface=render_surface,
        )
        return True
    except Exception as e:
        print(f"  [WARN] group rerun failed for {cfg['label']}: {e}")
        return False


def group_map_profile(sri, launcher, cfg, seed, output_suffix, common_base_names,
                      verbose=True):
    """
    Locate the group map produced by run_group_analysis_for_selection and
    turn it into a region profile restricted to common_base_names.

    Returns the profile DataFrame, or None (with a printed reason) if the
    map or a grid-matching atlas volume couldn't be resolved.
    """
    bids_name = Path(cfg["bids_dir"]).name
    map_path = sri.find_group_map_path(
        launcher.RESULTS_BASE_DIR, cfg["species"], bids_name, seed,
        stat_suffix=launcher.STAT_SUFFIX, category=launcher.CATEGORY,
        output_suffix=output_suffix)
    if not map_path or not os.path.exists(map_path):
        print(f"  [WARN] {cfg['label']}: group map not found at {map_path}")
        return None

    rep_atlas = sri.pick_representative_atlas_volume(cfg["bids_dir"])
    atlas_vol_path = launcher._find_matching_atlas_volume(cfg["bids_dir"], map_path, rep_atlas)
    if atlas_vol_path is None:
        print(f"  [WARN] {cfg['label']}: no atlas volume matching the group map's grid")
        return None

    label_df = sri.parse_itk_snap_label_file(cfg["label_path"])
    prof = sri.compute_region_profile(map_path, atlas_vol_path, label_df,
                                      atlas_level=launcher.ATLAS_LEVEL)
    return restrict_profile(prof, common_base_names)


# =============================================================================
# 5) Phylogenetic test + figure
# =============================================================================

def test_phylogenetic_signal(distance_df, divergence_times=None, exclude_human=True,
                             verbose=True):
    """
    Correlate each species' distance-to-human against its divergence time
    from human.

    exclude_human: True (default) drops the human row -- it is at
    distance ~0 from itself BY CONSTRUCTION, and including a point that
    is fixed at (0, 0) inflates the correlation without evidence.

    Returns a dict with pearson_r/p and spearman_rho/p, plus the
    DataFrame actually used.
    """
    divergence_times = divergence_times or DIVERGENCE_TIME_MYA
    df = distance_df.copy()
    df["divergence_mya"] = df["species"].map(divergence_times)

    missing = df[df["divergence_mya"].isna()]["species"].unique()
    if len(missing):
        print(f"  [WARN] no divergence time for {sorted(missing)} -- excluded "
             f"(add them to DIVERGENCE_TIME_MYA)")
    df = df[df["divergence_mya"].notna()]

    if exclude_human:
        df = df[df["species"] != "Human"]

    df = df[df["distance_to_human"].notna()]
    if len(df) < 3:
        print(f"  [WARN] only {len(df)} species with both a distance and a divergence "
             f"time -- too few for a correlation")
        return dict(n=len(df), pearson_r=np.nan, pearson_p=np.nan,
                   spearman_rho=np.nan, spearman_p=np.nan, data=df)

    pr, pp = _stats.pearsonr(df["divergence_mya"], df["distance_to_human"])
    sr, sp = _stats.spearmanr(df["divergence_mya"], df["distance_to_human"])

    if verbose:
        print(f"\n  Phylogenetic signal (n={len(df)} species, human "
             f"{'excluded' if exclude_human else 'included'}):")
        print(f"    Pearson  r={pr:.4f}, p={pp:.4f}")
        print(f"    Spearman rho={sr:.4f}, p={sp:.4f}")
        if len(df) < 6:
            print(f"    [NOTE] with n={len(df)}, this is a trend at best -- "
                 f"a single species shifts it substantially.")

    return dict(n=len(df), pearson_r=pr, pearson_p=pp,
               spearman_rho=sr, spearman_p=sp, data=df)


def plot_distance_vs_divergence(distance_df, test_result=None, divergence_times=None,
                                title=None, output_path=None, dpi=300):
    """
    Scatter of distance-to-human against divergence time, one point per
    species, colored by species, with a fitted line and the correlation
    stats annotated. Human (distance 0 by construction) is drawn as an
    open marker to make clear it anchors the axis rather than
    contributing evidence.
    """
    divergence_times = divergence_times or DIVERGENCE_TIME_MYA
    df = distance_df.copy()
    df["divergence_mya"] = df["species"].map(divergence_times)
    df = df[df["divergence_mya"].notna() & df["distance_to_human"].notna()]
    if df.empty:
        print("  [plot] nothing to plot")
        return None

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for _, row in df.iterrows():
        is_human = row["species"] == "Human"
        ax.scatter(row["divergence_mya"], row["distance_to_human"],
                  s=140, color=SPECIES_COLORS.get(row["species"], "#888888"),
                  edgecolors="black", lw=1.0, zorder=3,
                  facecolors="none" if is_human else None,
                  label=row["species"])
        ax.annotate(row["species"], (row["divergence_mya"], row["distance_to_human"]),
                   textcoords="offset points", xytext=(7, 5), fontsize=8)

    fit_df = df[df["species"] != "Human"]
    if len(fit_df) >= 3:
        slope, intercept, r, p, se = _stats.linregress(
            fit_df["divergence_mya"], fit_df["distance_to_human"])
        xs = np.linspace(fit_df["divergence_mya"].min(), fit_df["divergence_mya"].max(), 100)
        ax.plot(xs, intercept + slope * xs, "k--", lw=1.4, alpha=0.7, zorder=2)

    if test_result:
        ax.text(0.03, 0.97,
               f"Pearson r={test_result['pearson_r']:.3f}, p={test_result['pearson_p']:.3f}\n"
               f"Spearman rho={test_result['spearman_rho']:.3f}, p={test_result['spearman_p']:.3f}\n"
               f"n={test_result['n']} species (human excluded from fit)",
               transform=ax.transAxes, va="top", fontsize=8,
               bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))

    ax.set_xlabel("Divergence time from human (Mya)", fontsize=10)
    ax.set_ylabel("Distance of group map to human group map", fontsize=10)
    ax.set_title(title or "Residual distance to human vs phylogenetic divergence",
                fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


# =============================================================================
# 6) Driver
# =============================================================================

def run_human_anchored_analysis(sri, launcher, human_label="Human_Park",
                                reference_methods=("kcore", "uniform_median"),
                                n_select=5, k_human=5, subsample_n=None,
                                run_group_analysis=True, render_surface=False, copy_qc_images=True,
                                output_dir=None, bids_configs=None,
                                value_col=None, metric=None):
    """
    Full pipeline, run once per reference method so the two can be
    compared side by side.

    run_group_analysis: True runs the real AFNI 3dttest per BIDS (slow).
    False skips it and reuses whatever group maps already exist at the
    expected paths -- useful to re-plot without recomputing.

    Returns {method: dict(distance_df=, test_result=, rankings=,
    selections=)}.
    """
    bids_configs = bids_configs if bids_configs is not None else launcher.BIDS_CONFIGS
    value_col = value_col or launcher.VALUE_COL
    metric = metric or launcher.METRIC

    # PASS 1: compute every BIDS's profiles ONCE, up front. Needed so the
    # common-region set can be derived from the profiles that actually
    # exist (see get_common_regions_from_profiles) rather than from the
    # label files, which list every atlas level and so vastly overstate
    # the overlap.
    print(f"\n{'=' * 70}\n  PASS 1: computing profiles for every BIDS\n{'=' * 70}")
    profiles_by_label = {}
    for cfg in bids_configs:
        label, seed = cfg["label"], cfg["seeds"][0]
        try:
            label_df = sri.parse_itk_snap_label_file(cfg["label_path"])
            profs, _ = sri.compute_all_run_profiles(
                cfg["species"], cfg["bids_dir"], label_df, launcher.CLASSIFICATION_PATH,
                category=launcher.CATEGORY, seeds=[seed], atlas_level=launcher.ATLAS_LEVEL)
            if profs:
                profiles_by_label[label] = profs
            else:
                print(f"  [SKIP] {label}: no run profiles computed")
        except Exception as e:
            print(f"  [ERROR] {label}: {e}")

    if not profiles_by_label:
        raise ValueError("No BIDS produced run profiles -- nothing to compare.")

    print(f"\n{'=' * 70}\n  Common regions across all BIDS' PROFILES\n{'=' * 70}")
    common = get_common_regions_from_profiles(profiles_by_label, value_col=value_col)
    if not common:
        raise ValueError("No common base regions across BIDS profiles -- cannot compare.")

    results = {}
    for method in reference_methods:
        print(f"\n{'#' * 70}\n  REFERENCE METHOD: {method}\n{'#' * 70}")

        human_profiles = profiles_by_label.get(human_label)
        if human_profiles is None:
            raise ValueError(f"{human_label} produced no profiles -- cannot anchor on it.")
        human_ref, _ = build_reference_for_profiles(
            sri, human_profiles, method=method, k=k_human,
            value_col=value_col, metric=metric, verbose=True)
        print(f"  Human reference ({method}) from {human_label} "
             f"({len(human_profiles)} run(s))")

        # Every BIDS's OWN reference (same method, its own runs) -- "what
        # a typical run from this species looks like", independent of the
        # human-anchored selection. Used for figure 3 (alongside
        # human_ref) and for the run-correlation table below.
        own_ref_by_label = {}
        for lbl, profs in profiles_by_label.items():
            try:
                own_ref_by_label[lbl], _ = build_reference_for_profiles(
                    sri, profs, method=method, k=k_human,
                    value_col=value_col, metric=metric, verbose=False)
            except Exception as e:
                print(f"  [WARN] {lbl}: could not build its own reference: {e}")

        rankings, selections, distance_rows = {}, {}, []
        bids_ref_profiles = {}

        for cfg in bids_configs:
            label = cfg["label"]
            seed = cfg["seeds"][0]
            print(f"\n{'=' * 70}\n  {label} / {seed}  [{method}]\n{'=' * 70}")

            try:
                profiles = profiles_by_label.get(label)
                if not profiles:
                    print(f"  [SKIP] no run profiles computed")
                    continue

                selected_keys, ranking_df = select_runs_closest_to_reference(
                    sri, profiles, human_ref, n_select=n_select,
                    common_base_names=common, value_col=value_col, metric=metric,
                    subsample_n=subsample_n)
                if not selected_keys:
                    print(f"  [SKIP] no runs could be ranked against the human reference")
                    continue

                rankings[label] = ranking_df.drop(columns="_key")
                selections[label] = [tuple(str(x) for x in kk) for kk in selected_keys]

                output_suffix = f"{seed}_HumanAnchored_{method}"
                df_class_sel = selected_keys_to_df_class(profiles, selected_keys)
                if run_group_analysis:
                    df_class = df_class_sel
                    ok = run_group_analysis_for_selection(
                        sri, cfg, seed, df_class, output_suffix,
                        render_surface=render_surface)
                    if not ok:
                        continue

                group_prof = group_map_profile(sri, launcher, cfg, seed, output_suffix, common)
                if group_prof is None:
                    continue

                if label == human_label:
                    human_group_prof = group_prof
                    distance_rows.append(dict(label=label, species=cfg["species"],
                                             distance_to_human=0.0, pearson_r=1.0,
                                             n_selected=len(selected_keys)))
                else:
                    distance_rows.append(dict(label=label, species=cfg["species"],
                                             _group_prof=group_prof,
                                             n_selected=len(selected_keys)))

                bids_ref_profiles[label] = group_prof

                # Copy each run's SBA preview jpg into selected/not_selected,
                # same as the main pipeline's COPY_QC_IMAGES step -- so the
                # human-anchored selection can be eyeballed run by run.
                if copy_qc_images and output_dir:
                    images_dir = opj(output_dir, method, "images", label)
                    try:
                        sri.copy_run_qc_images(df_class_sel, cfg["bids_dir"], images_dir)
                    except Exception as e:
                        print(f"  [WARN] {label}: jpg copy failed: {e}")

                # FIGURE 3 (per BIDS): species reference AND human
                # reference (distinct colors), plus every individual run,
                # in the same common-region space.
                if output_dir and label in own_ref_by_label:
                    per_bids_dir = opj(output_dir, method, "per_bids")
                    os.makedirs(per_bids_dir, exist_ok=True)
                    plot_bids_ref_with_individuals(
                        profiles, own_ref_by_label[label], human_ref, common,
                        selected_keys=selected_keys, value_col=value_col,
                        title=f"{label} - species vs human reference ({method})",
                        output_path=opj(per_bids_dir, f"{label}_ref_and_individuals.png"))
            except Exception as e:
                print(f"  [ERROR] {label}: {e}")

        # Resolve distances against the human group map (computed above)
        human_rows = [r for r in distance_rows if r.get("distance_to_human") == 0.0]
        if not human_rows:
            print(f"\n[WARN] {method}: no human group map -- cannot compute distances")
            results[method] = dict(distance_df=pd.DataFrame(), test_result=None,
                                  rankings=rankings, selections=selections)
            continue

        final_rows = []
        for r in distance_rows:
            if r.get("distance_to_human") == 0.0:
                final_rows.append({kk: vv for kk, vv in r.items() if not kk.startswith("_")})
                continue
            try:
                dist, pr, n = sri.profile_distance(
                    r["_group_prof"], human_group_prof, value_col=value_col, metric=metric)
                final_rows.append(dict(label=r["label"], species=r["species"],
                                      distance_to_human=dist, pearson_r=pr,
                                      n_selected=r["n_selected"]))
            except Exception as e:
                print(f"  [WARN] {r['label']}: distance to human failed: {e}")

        # Narrow the plotting axis to regions present in EVERY group map,
        # not just in the run profiles. Confirmed real symptom otherwise:
        # figure 2 drew gaps for several species, because a region can have
        # voxels in the individual run profiles (which is what
        # get_common_regions_from_profiles intersects) yet be absent from
        # the thresholded 3dttest GROUP map that figure 2 actually plots.
        # The region set is FIXED at step 1 and used unchanged everywhere --
        # no post-hoc narrowing. If a group map is nevertheless missing a
        # region, that is a real data problem to surface, not something to
        # silently filter away, so report exactly which ones and where.
        for lab, gp in [(human_label, human_group_prof)] + sorted(bids_ref_profiles.items()):
            vals = _ordered_common_values(gp, common, value_col)
            missing = [n for n, v in zip(common, vals) if not np.isfinite(v)]
            if missing:
                print(f"  [WARN] {lab}: group map has no value for {len(missing)}/{len(common)} "
                     f"step-1 common region(s): {missing}")

        distance_df = pd.DataFrame(final_rows)

        # STEP 5 SANITY CHECK: the OTHER human BIDS should fit the human
        # reference best of all. It is the same species, so if some animal
        # BIDS beats it, the pipeline is not measuring human-likeness --
        # it is measuring something else (profile flatness, amplitude, or
        # a shared artifact), and the phylogenetic reading below is not
        # trustworthy. Reported rather than assumed.
        if not distance_df.empty and "distance_to_human" in distance_df.columns:
            ranked_fit = distance_df[distance_df["label"] != human_label].sort_values(
                "distance_to_human")
            if not ranked_fit.empty:
                best = ranked_fit.iloc[0]
                other_humans = ranked_fit[ranked_fit["species"] == "Human"]
                print(f"\n  [step-5 check] closest to the human reference: "
                     f"{best['label']} (d={best['distance_to_human']:.4f})")
                if not other_humans.empty:
                    oh = other_humans.iloc[0]
                    rank = int(ranked_fit.reset_index(drop=True).index[
                        ranked_fit.reset_index(drop=True)["label"] == oh["label"]][0]) + 1
                    print(f"  [step-5 check] other human BIDS ({oh['label']}) ranks "
                         f"{rank}/{len(ranked_fit)} (d={oh['distance_to_human']:.4f})")
                    if rank != 1:
                        print(f"  [WARN] a NON-human BIDS fits the human reference better "
                             f"than the other human dataset -- the selection is likely "
                             f"not capturing human-likeness. Treat the phylogenetic "
                             f"result below as uninterpretable until this is resolved.")
                else:
                    print(f"  [INFO] no second human BIDS available -- step-5 check skipped")
        test_result = test_phylogenetic_signal(distance_df) if not distance_df.empty else None

        if output_dir:
            method_dir = opj(output_dir, method)
            os.makedirs(method_dir, exist_ok=True)

        # Always build the correlation table, whether or not output_dir is
        # given -- it's cheap (pure profile correlations, no AFNI) and
        # genuinely useful returned in memory even without saving to disk.
        species_lut = {c["label"]: c["species"] for c in bids_configs}
        corr_table = build_run_correlation_table(
            sri, profiles_by_label, own_ref_by_label, human_ref, common,
            selections=selections, value_col=value_col, metric=metric)
        corr_table["species"] = corr_table["label"].map(species_lut)

        if output_dir:
            # FIGURE 1: the human reference over the common regions.
            plot_common_regions_with_human_ref(
                restrict_profile(human_ref, common), common, value_col=value_col,
                title=f"Human reference ({human_label}, {method}) over "
                      f"{len(common)} common regions",
                output_path=opj(method_dir, f"01_human_reference_common_regions_{method}.png"))

            # FIGURE 2: every BIDS reference overlaid on the human one.
            if bids_ref_profiles:
                plot_all_bids_refs_vs_human(
                    bids_ref_profiles, human_group_prof, common,
                    value_col=value_col,
                    title=f"Each BIDS reference vs human reference ({method})",
                    output_path=opj(method_dir, f"02_all_refs_vs_human_{method}.png"))
            excel_path = opj(method_dir, f"human_anchored_{method}.xlsx")
            with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
                distance_df.to_excel(writer, sheet_name="distance_to_human", index=False)
                corr_table.to_excel(writer, sheet_name="run_correlations", index=False)
                for label, rdf in rankings.items():
                    rdf.head(50).to_excel(writer, sheet_name=f"rank_{label[:24]}", index=False)
            print(f"  [saved] {excel_path}")

            plot_distance_vs_divergence(
                distance_df, test_result=test_result,
                title=f"Residual distance to human vs divergence time ({method})",
                output_path=opj(method_dir, f"distance_vs_divergence_{method}.png"))

        results[method] = dict(distance_df=distance_df, test_result=test_result,
                              rankings=rankings, selections=selections,
                              correlation_table=corr_table,
                              own_ref_by_label=own_ref_by_label, human_ref=human_ref)

    return results


# =============================================================================
# 7) Profile figures (common regions / refs overlay / per-BIDS individuals)
# =============================================================================

def _drop_all_nan_regions(profiles_list, common_base_names, value_col="median_r"):
    """
    Keep only regions that have a real value in AT LEAST ONE of the given
    profiles. Confirmed real symptom without this: the figures drew all
    138 candidate region ticks while only a handful carried data, making
    the plots unreadable and implying coverage that does not exist.
    """
    keep = []
    for i, name in enumerate(common_base_names):
        vals = [_ordered_common_values(p, [name], value_col)[0] for p in profiles_list
                if p is not None]
        if any(np.isfinite(v) for v in vals):
            keep.append(name)
    return keep


def _sort_regions_by_value(reference_profile, common_base_names, value_col="median_r"):
    """
    Order regions by their value in `reference_profile`, ascending
    (lowest median_r first) -- the same convention as the per-BIDS
    profile figures, where a monotonically rising reference curve makes
    each run's deviation from it immediately readable. Regions with no
    value in the reference keep their original relative order at the end.
    """
    vals = _ordered_common_values(reference_profile, common_base_names, value_col)
    finite = [(v, n) for v, n in zip(vals, common_base_names) if np.isfinite(v)]
    missing = [n for v, n in zip(vals, common_base_names) if not np.isfinite(v)]
    return [n for _, n in sorted(finite, key=lambda t: t[0])] + missing


def _keep_regions_present_in_all(profiles_list, common_base_names, value_col="median_r"):
    """
    Keep only regions that have a real value in EVERY given profile --
    the strict counterpart to _drop_all_nan_regions (which keeps a region
    if ANY profile has it). Needed for the cross-species overlay figure:
    a region missing from even one species' group map leaves a visible
    GAP in that species' line, so the shared axis must be the strict
    intersection, not the union.
    """
    keep = []
    for name in common_base_names:
        vals = [_ordered_common_values(p, [name], value_col)[0]
                for p in profiles_list if p is not None]
        if vals and all(np.isfinite(v) for v in vals):
            keep.append(name)
    return keep


def _ordered_common_values(profile, common_base_names, value_col="median_r"):
    """
    Values of `profile` for each name in common_base_names, in THAT
    order (NaN where the profile has no such region), so every species
    plots on the same x-axis positions and is directly comparable.
    """
    if profile is None or "base_name" not in profile.columns:
        return np.full(len(common_base_names), np.nan)
    lut = (profile.assign(_b=profile["base_name"].str.lower().str.strip())
                  .groupby("_b")[value_col].mean())
    return np.array([lut.get(b, np.nan) for b in common_base_names], dtype=float)


def plot_common_regions_with_human_ref(human_ref, common_base_names,
                                       value_col="median_r", title=None,
                                       output_path=None, dpi=300,
                                       sort_by_reference=True):
    """
    FIGURE 1 -- the human reference profile across the regions common to
    every species: what the anchor actually looks like, region by region,
    on the exact axis every other species will be compared on.
    """
    if sort_by_reference:
        common_base_names = _sort_regions_by_value(human_ref, common_base_names, value_col)
    vals = _ordered_common_values(human_ref, common_base_names, value_col)
    x = np.arange(len(common_base_names))

    fig, ax = plt.subplots(figsize=(max(10, len(common_base_names) * 0.22), 5))
    ax.plot(x, vals, color="#D35400", lw=1.8, marker="o", ms=3, label="Human reference")
    ax.axhline(0, color="#999999", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(common_base_names, rotation=90, fontsize=4)
    ax.set_ylabel(value_col, fontsize=10)
    ax.set_xlabel(f"{len(common_base_names)} common regions"
                 f"{' (ordered by reference value)' if sort_by_reference else ''}", fontsize=9)
    ax.set_title(title or "Human reference profile over cross-species common regions",
                fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_all_bids_refs_vs_human(bids_ref_profiles, human_ref, common_base_names,
                                value_col="median_r", title=None,
                                output_path=None, dpi=300,
                                sort_by_reference=True):
    """
    FIGURE 2 -- every BIDS's own selected-group reference overlaid on the
    human reference, all on the common-region axis.

    bids_ref_profiles: {label: profile} -- typically each BIDS's group
    map profile from group_map_profile(), i.e. the reference built from
    its selected runs.

    Human is drawn thick/black so the others read as deviations from it.
    """
    if sort_by_reference:
        # ordered by the HUMAN reference, so every species is plotted on
        # one shared ordering and deviations from human read directly
        common_base_names = _sort_regions_by_value(human_ref, common_base_names, value_col)
    x = np.arange(len(common_base_names))
    fig, ax = plt.subplots(figsize=(max(11, len(common_base_names) * 0.24), 6))

    ax.plot(x, _ordered_common_values(human_ref, common_base_names, value_col),
           color="black", lw=2.6, label="Human reference", zorder=5)

    for label, prof in sorted(bids_ref_profiles.items()):
        species = next((s for s in SPECIES_COLORS if label.startswith(s)), None)
        ax.plot(x, _ordered_common_values(prof, common_base_names, value_col),
               lw=1.2, alpha=0.85, label=label,
               color=SPECIES_COLORS.get(species, None))

    ax.axhline(0, color="#999999", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(common_base_names, rotation=90, fontsize=4)
    ax.set_ylabel(value_col, fontsize=10)
    ax.set_title(title or "Each BIDS reference vs the human reference (common regions)",
                fontsize=11, fontweight="bold")
    ax.legend(fontsize=7, ncol=2)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_bids_ref_with_individuals(profiles, species_ref, human_ref, common_base_names,
                                   selected_keys=None, value_col="median_r",
                                   title=None, output_path=None, dpi=300,
                                   sort_by_reference=True):
    """
    FIGURE 3 (per BIDS) -- BOTH the species' own reference AND the human
    reference, plus every individual run, in the same common-region space.

    Confirmed real problem this fixes: the previous version plotted only
    the POST-HOC group map (built from the 3dttest++ on the selected
    runs) as a single black "Reference" line -- which is NOT the human
    reference the selection was actually based on, is DIFFERENT for every
    species (so it looked wrong that the black line "should be the same
    for everyone" -- it should, but wasn't, because it wasn't the human
    reference at all), and can have GAPS where the thresholded group map
    lost coverage in a region the run-level profiles still had. Using two
    PROFILE-based references (species_ref, human_ref -- both built the
    same way as the common-region list itself, so complete by
    construction) fixes both: distinct colors for what is actually being
    compared, and no coverage gaps.

    species_ref: this BIDS's OWN reference (e.g. compute_kcore_reference
    over ALL its own runs) -- "what does a typical run from this species
    look like".
    human_ref: the fixed human reference every species was ranked against
    for selection -- "what were runs actually selected for resembling".

    selected_keys: if given, runs in this set are drawn highlighted
    (green) and the rest faint grey, showing WHICH runs the human-anchored
    selection actually kept relative to the overall spread.
    """
    if sort_by_reference:
        common_base_names = _sort_regions_by_value(human_ref, common_base_names, value_col)
    x = np.arange(len(common_base_names))
    selected_str = {tuple(str(v) for v in k) for k in (selected_keys or [])}

    fig, ax = plt.subplots(figsize=(max(11, len(common_base_names) * 0.24), 6))

    n_sel = n_other = 0
    for key, prof in profiles.items():
        is_sel = tuple(str(v) for v in key) in selected_str
        ax.plot(x, _ordered_common_values(prof, common_base_names, value_col),
               lw=0.8 if is_sel else 0.5,
               color="#009E73" if is_sel else "#BBBBBB",
               alpha=0.85 if is_sel else 0.35,
               zorder=3 if is_sel else 1)
        n_sel += is_sel
        n_other += (not is_sel)

    ax.plot(x, _ordered_common_values(species_ref, common_base_names, value_col),
           color="#1A5276", lw=2.6, zorder=5, label="Species reference", ls="-")
    ax.plot(x, _ordered_common_values(human_ref, common_base_names, value_col),
           color="#D35400", lw=2.6, zorder=6, label="Human reference", ls="--")

    handles = [
        plt.Line2D([], [], color="#1A5276", lw=2.6, label="Species reference"),
        plt.Line2D([], [], color="#D35400", lw=2.6, ls="--", label="Human reference"),
        plt.Line2D([], [], color="#009E73", lw=1.2, label=f"Selected runs (n={n_sel})"),
        plt.Line2D([], [], color="#BBBBBB", lw=1.2, label=f"Other runs (n={n_other})"),
    ]
    ax.legend(handles=handles, fontsize=8)

    ax.axhline(0, color="#999999", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(common_base_names, rotation=90, fontsize=4)
    ax.set_ylabel(value_col, fontsize=10)
    ax.set_title(title or "Species vs human reference, with individual runs (common regions)",
                fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def get_common_regions_from_profiles(profiles_by_label, value_col="median_r",
                                     require_nonnan=True, verbose=True):
    """
    Base regions common to every BIDS's ACTUALLY COMPUTED profiles --
    replaces get_common_lvl_regions for anything that compares profile
    values.

    Confirmed real problem this fixes: get_common_lvl_regions intersects
    the LABEL FILES, which list every region at every atlas level (276
    entries), while the profiles are computed at ONE level (ATLAS_LEVEL).
    The label-file intersection therefore returned ~138 names of which
    only ~12 existed in the profiles, so plots showed 138 x-axis ticks
    with a dozen points, and -- much worse -- run selection ranked runs
    on a dozen regions instead of the full overlap.

    require_nonnan: True (default) additionally drops regions that are
    NaN (no voxels) in some BIDS, since those contribute nothing to a
    correlation and just pad the axis.

    profiles_by_label: {label: {key: profile_df}} -- every BIDS's full
    per-run profile dict.
    """
    per_label_sets = []
    for label, profiles in profiles_by_label.items():
        names = set()
        for prof in profiles.values():
            if prof is None or "base_name" not in prof.columns:
                continue
            sub = prof
            if require_nonnan and value_col in prof.columns:
                sub = prof[prof[value_col].notna()]
            names |= set(sub["base_name"].str.lower().str.strip())
        if names:
            per_label_sets.append(names)
            if verbose:
                print(f"    {label}: {len(names)} region(s) with data")

    if not per_label_sets:
        return []
    common = sorted(set.intersection(*per_label_sets))
    if verbose:
        print(f"  {len(common)} region(s) common to all {len(per_label_sets)} BIDS' profiles")
    return common
