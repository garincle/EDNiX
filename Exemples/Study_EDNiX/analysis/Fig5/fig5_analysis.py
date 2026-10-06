"""
fig5_analysis.py
================
All PURE-ANALYSIS functions for Figure 5, consolidated into one module.

Consolidated from what used to be spread over three files:
  - consensus_anchored_selection.py   (consensus building + selection)
  - human_anchored_phylogeny.py       (region sets, references, figures)
  - Figure5_SBA_KeeporNot.py          (used only as a constants holder --
                                       those constants now live in
                                       fig5_pipeline.CONFIG)

WHAT IS AND IS NOT HERE
-----------------------
Everything in this file is pure numpy/pandas/matplotlib: profiles in,
numbers and figures out. Nothing here touches AFNI, nibabel, atlas files
or the filesystem layout of a BIDS tree.

That work stays in SBA_Region_Profile_Index (imported as `sri` by the
caller), which remains the one external dependency. Exactly five of its
functions are used:
    parse_itk_snap_label_file   read an ITK-SNAP atlas label file
    compute_all_run_profiles    per-run regional profiles from SBA maps
    rerun_balanced_by_run       run-level 3dttest++ (+ optional surfaces)
    find_group_map_path         locate a group map on disk
    copy_run_qc_images          copy SBA preview jpgs
Those wrap AFNI, nibabel and the atlas library; inlining them here would
mean duplicating that machinery, so they are deliberately left in place.

Call graph is one-directional: fig5_pipeline -> fig5_analysis -> sri.
No module here calls back into the pipeline.
"""

import os
import math

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as _stats

opj = os.path.join



# === Constants ===========================================================

# Same 8-color, colorblind-safe (Okabe-Ito) palette as Supp_Figure_4-6_QC.py's
# PALETTE, so every figure in this module uses the same color code the QC
# figures already use -- kept as a literal copy (not imported) since
# fig5_analysis.py has no dependency on ednix_bids_tools.py.
PALETTE = [
    "#0072B2", "#E69F00", "#009E73", "#CC79A7",
    "#56B4E9", "#D55E00", "#F0E442", "#000000",
]


def _palette_color_map(names):
    """{name: color}, names sorted then cycled through PALETTE -- same
    assignment rule Supp_Figure_4-6_QC.py uses for its GLOBAL_BIDS_COL
    (sorted(set(...)), PALETTE[i % len(PALETTE)]). With more than 8
    distinct names, colors repeat; harmless for reading one panel, but
    two BIDS can end up sharing a color once there are >8 -- true here
    with 10 BIDS, see the two BIDS falling on Rat and Mouse's colors."""
    return {n: PALETTE[i % len(PALETTE)] for i, n in enumerate(sorted(set(names)))}


# Species-level color map (F0/F1/F4/multilevel-F6, where a point is one
# species, not one BIDS): same PALETTE, keyed by species instead of BIDS,
# so these panels also follow the QC figures' color code even though they
# have no BIDS-level color -- swap SPECIES_COLORS with your own dict if
# you want species colors independent of this palette.
SPECIES_COLORS = _palette_color_map(
    ["Human", "Macaque", "Marmoset", "Mouselemur", "Rat", "Mouse", "Dog"])


# Divergence time (Mya) between human and each species = age of their most
# recent common ancestor. PROVENANCE: entered by hand as approximate TimeTree
# median estimates (Kumar et al. 2022, TimeTree 5, timetree.org); they were
# NOT downloaded or queried by this code. Before publication, look each pair
# up on timetree.org, replace the values, and report the access date and the
# TimeTree version in the methods. fig5_phylo.SPECIES_CLADES must stay
# consistent with this dict (checked automatically at run time).
DIVERGENCE_TIME_MYA = {
    "Human": 0.0,
    "Macaque": 29.0,
    "Marmoset": 43.0,
    "Mouselemur": 74.0,
    "Rat": 87.0,
    "Mouse": 87.0,
    "Dog": 94.0,
}



# === Low-level helpers ===================================================

VALID_NORMALIZATIONS = (True, False, None, "zscore", "scale", "none")


def _normalize_values(vals, mode):
    """
    Per-profile normalization applied before profiles are combined into
    a consensus or compared. `mode` is cfg["ZSCORE"]:
      True / "zscore" -- subtract the profile's own mean, divide by its
          SD. Equalizes amplitude across datasets but RE-CENTERS each
          profile on zero, so roughly half of every profile becomes
          "negative" regardless of the real sign of the connectivity.
      "scale" -- divide by the profile's root-mean-square, NO centering.
          Equalizes amplitude across datasets (so a high-amplitude
          dataset such as Marmoset_Tian cannot dictate the consensus at
          every region) while keeping the true sign: a positive region
          stays positive, a real anti-correlation stays negative, and 0
          stays 0.
      False / None / "none" -- raw values. Sign preserved, but amplitude
          NOT equalized: a region-wise median across species with very
          different amplitudes returns, region by region, whichever
          species happens to sit in the middle there -- a patchwork
          whose SHAPE mixes datasets, not a shared pattern.
    Pearson correlation is invariant to all three, so run RANKING under
    SIMILARITY="correlation" is unaffected; what changes is the
    consensus each run is ranked against.
    """
    if mode is True or mode == "zscore":
        return _zscore_values(vals)
    if mode == "scale":
        v = np.asarray(vals, dtype=float)
        f = np.isfinite(v)
        if f.sum() < 2:
            return v
        rms = float(np.sqrt(np.mean(v[f] ** 2)))
        if rms == 0:
            return np.full_like(v, np.nan, dtype=float)
        out = np.full_like(v, np.nan, dtype=float)
        out[f] = v[f] / rms
        return out
    if mode is False or mode is None or mode == "none":
        return np.asarray(vals, dtype=float)
    raise ValueError(f"Unknown normalization {mode!r} -- use one of {VALID_NORMALIZATIONS}")


def _norm_label(mode):
    """Human-readable axis/legend prefix for a normalization mode."""
    if mode is True or mode == "zscore":
        return "z-scored "
    if mode == "scale":
        return "RMS-scaled "
    return ""


def _zscore_values(vals):
    """Z-score a value vector, ignoring NaN. Returns all-NaN unchanged if
    there is no spread to normalize (constant profile)."""
    vals = np.asarray(vals, dtype=float)
    finite = np.isfinite(vals)
    if finite.sum() < 2:
        return vals
    mu, sd = np.nanmean(vals[finite]), np.nanstd(vals[finite])
    if sd == 0:
        return np.full_like(vals, np.nan, dtype=float)
    out = np.full_like(vals, np.nan, dtype=float)
    out[finite] = (vals[finite] - mu) / sd
    return out


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


def restrict_profile(profile, common_base_names):
    """Keep only rows whose base_name is in common_base_names."""
    if profile is None or "base_name" not in profile.columns:
        return profile
    mask = profile["base_name"].str.lower().str.strip().isin(common_base_names)
    return profile[mask].copy()



# === Region sets and references ==========================================

def get_common_regions_from_profiles(profiles_by_label, value_col="median_r",
                                     require_nonnan=False, verbose=True):
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

    require_nonnan: False (default, changed from the original True) --
    a region is kept as soon as it is STRUCTURALLY present (its
    base_name appears somewhere in the profile, i.e. the atlas actually
    defines it for that BIDS), regardless of whether any individual
    run's value is NaN there. Every consumer downstream that reads off
    `common` (_ordered_common_values, compute_region_invariance,
    profile_distance, build_consensus_pattern, ...) already fills a
    missing/NaN region with np.nan and skips it with an isfinite/
    nanmedian mask at the point of use, region by region -- so an
    occasional NaN run never breaks anything downstream, it is simply
    ignored exactly where it occurs. Confirmed real problem this fixes:
    with the old True default, a SINGLE NaN in a SINGLE run, for a
    SINGLE species, anywhere, dropped that region from `common` for
    EVERY species -- at atlas_level=3 with UNIT="run" (hundreds of
    individual run profiles per dataset) this left only 18/276 regions,
    most of the anatomy silently discarded from every figure and from
    the region-invariance computation that decides which regions the
    run selection scores on.

    Pass True to restore the old, strictly conservative all-or-nothing
    behavior (drop a region globally the instant ANY run anywhere is
    NaN for it) if you ever need to double-check a result against the
    stricter criterion.

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


def get_common_regions_from_group_maps(group_profiles, value_col="median_r",
                                       verbose=True):
    """
    Common regions taken from the ALL-INDIVIDUALS group maps themselves:
    a region is kept when it has a real (non-NaN) value in EVERY
    dataset's group-map profile.

    This single criterion enforces both conditions that matter, because
    a group-map profile is computed by compute_region_profile at ONE
    atlas_level:
      - the region must exist at that atlas level in that dataset's
        atlas (otherwise no label_id is painted there -> NaN), so the
        result is automatically the intersection ACROSS ATLASES at the
        requested level;
      - the region must actually be sampled by the group z-map/c-map
        (voxels present and finite -> non-NaN), so regions the group map
        simply does not cover are dropped.

    Why this replaces the per-RUN criterion: the old
    get_common_regions_from_profiles intersected the per-run profiles,
    so ONE NaN run, for ONE species, anywhere among several hundred
    runs, removed that region from every species (confirmed: 18/276
    regions survived at atlas_level=3 with UNIT="run"). Evaluated on the
    group map instead, a region only disappears when it is genuinely
    absent or unsampled for that whole dataset -- which is the condition
    that actually matters, since the group map is what every figure
    downstream is drawn from.

    Note there is deliberately no separate enumeration of "all regions
    defined at level L" from the label files: parse_itk_snap_label_file
    returns no level column, because level membership is a property of
    the 4D atlas VOLUME (one 3D label slice per level), knowable only by
    reading it -- which is exactly what computing the group-map profile
    at atlas_level already does.

    group_profiles: {label: (species, profile_df)} -- the output of
    fig5_units.profiles_from_group_maps.

    Returns the sorted list of common base names.
    """
    per_label_sets = []
    for label, entry in sorted(group_profiles.items()):
        prof = entry[1] if isinstance(entry, (tuple, list)) else entry
        if prof is None or "base_name" not in prof.columns:
            continue
        sub = prof[prof[value_col].notna()] if value_col in prof.columns else prof
        names = set(sub["base_name"].str.lower().str.strip())
        if names:
            per_label_sets.append(names)
            if verbose:
                print(f"    {label}: {len(names)} region(s) present in its group map")

    if not per_label_sets:
        print("  [WARN] no group-map profile could be read -- no common regions")
        return []
    common = sorted(set.intersection(*per_label_sets))
    if verbose:
        print(f"  {len(common)} region(s) common to all {len(per_label_sets)} "
              f"group map(s)")
    return common


def build_reference_for_profiles(sri, profiles, method="kcore", k=5, value_col="median_r",
                                 metric="correlation", verbose=False):
    """
    Build a reference profile from an ALREADY-COMPUTED profiles dict,
    using the SAME method ("kcore" or "uniform_median") for whichever
    BIDS it's called on -- every dataset's reference (all runs, kept
    only, not-kept only, any selection round) is built this identical way,
    so they are directly comparable (same construction, same method, just
    different input runs).
    """
    if method == "kcore":
        ref, core_keys, _, _ = sri.compute_kcore_reference(
            profiles, k=k, value_col=value_col, metric=metric, verbose=verbose)
        return ref, dict(core_keys=core_keys)
    if method == "uniform_median":
        ref = sri.compute_group_median_profile(profiles, set(profiles.keys()))
        return ref, dict(n_runs=len(profiles))
    raise ValueError(f"Unknown method {method!r} -- use 'kcore' or 'uniform_median'.")


def species_profile(sri, profiles, source="all_runs", method="kcore", k=5,
                    value_col="median_r", metric="correlation", verbose=False):
    """
    One species' contribution to the consensus.

    source="all_runs": median across EVERY run of that species -- no
    selection at any point, so nothing about this profile reflects a
    prior choice. Preferred when the consensus will later be used to
    select runs.
    source="reference": that species' own kcore/median reference (see
    SBA_Region_Profile_Index.compute_kcore_reference) -- sharper, but the
    medoid/kcore step is itself a selection, so a consensus built this
    way is one step less independent of the runs it will later judge.
    """
    if source == "all_runs":
        return sri.compute_group_median_profile(profiles, set(profiles.keys()))
    if source == "reference":
        if method == "kcore":
            ref, _, _, _ = sri.compute_kcore_reference(
                profiles, k=k, value_col=value_col, metric=metric, verbose=verbose)
            return ref
        if method == "uniform_median":
            return sri.compute_group_median_profile(profiles, set(profiles.keys()))
        raise ValueError(f"Unknown method {method!r}")
    raise ValueError(f"Unknown source {source!r} -- use 'all_runs' or 'reference'")



# === Consensus and run selection =========================================

def diagnose_consensus_sign(species_profiles, common_base_names, value_col="median_r",
                            species_of_label=None, verbose=True):
    """
    Ground truth for "the consensus is never negative -- is that a bug?":
    for every region, the per-species (or per-BIDS, if species_of_label is
    None) contributor values that build_consensus_pattern's median is
    actually computed from, plus how many are positive/negative/NaN.

    A region where the median stays (barely) positive despite several
    negative contributors is expected, not a bug -- median only needs a
    majority on one side, and individual species/BIDS dipping negative
    from noise doesn't change that as long as most don't. This function
    lets you SEE which situation you're in per region instead of
    inferring it from the aggregate range alone.

    Returns a DataFrame: region, median (raw, no z-score -- matches what
    F4's Cross-species consensus line actually plots), n_pos, n_neg,
    n_nan, n_total, one column per species/BIDS with its raw value.
    """
    contributors = {}
    for label, prof in species_profiles.items():
        key = (species_of_label or {}).get(label, label)
        vals = _ordered_common_values(prof, common_base_names, value_col)
        contributors.setdefault(key, []).append(vals)
    with np.errstate(invalid="ignore"):
        by_key = {k: np.nanmedian(np.vstack(v), axis=0) for k, v in contributors.items()}

    rows = []
    for j, region in enumerate(common_base_names):
        vals = {k: by_key[k][j] for k in by_key}
        finite = [v for v in vals.values() if np.isfinite(v)]
        row = dict(region=region,
                   median=float(np.median(finite)) if finite else np.nan,
                   n_pos=sum(v > 0 for v in finite), n_neg=sum(v < 0 for v in finite),
                   n_nan=sum(not np.isfinite(v) for v in vals.values()),
                   n_total=len(vals))
        row.update({f"val_{k}": vals[k] for k in sorted(vals)})
        rows.append(row)
    df = pd.DataFrame(rows).sort_values("median")
    if verbose:
        close_calls = df[(df["median"] > 0) & (df["n_neg"] >= df["n_pos"])]
        print(f"  [diagnose_consensus_sign] {len(df)} region(s), "
              f"{(df['median'] < 0).sum()} with negative median, "
              f"{(df['median'] > 0).sum()} positive")
        if len(close_calls):
            print(f"    {len(close_calls)} region(s) positive median despite "
                 f">= as many negative as positive contributors (median barely "
                 f"tipping positive, not a clean majority):")
            print(close_calls[["region", "median", "n_pos", "n_neg", "n_nan"]]
                  .to_string(index=False))
    return df


def build_consensus_pattern(species_profiles, common_base_names, value_col="median_r",
                            exclude_label=None, zscore=True, threshold_abs=None,
                            species_of_label=None, verbose=True):
    """
    Combine per-BIDS profiles into ONE consensus pattern over
    common_base_names.

    species_of_label: {label: species}. When given, BIDS belonging to the
    SAME species are first collapsed (median) into one line, and the
    consensus is then the median ACROSS SPECIES. Without this, a species
    contributed by several datasets (e.g. 2 Macaque BIDS, 2 Marmoset, 2
    Human) counts two or three times as much as one represented by a
    single dataset -- the consensus then describes "the species we
    happen to have most datasets for", not what is conserved across
    evolution. Strongly recommended; None keeps the old per-BIDS
    behaviour only for backward comparison.

    exclude_label: the BIDS to leave OUT (leave-one-out). When
    species_of_label is given, the WHOLE species that BIDS belongs to is
    excluded, not just that one dataset -- otherwise a Macaque BIDS
    would still be scored against a consensus containing the other
    Macaque BIDS, which is exactly the leak leave-one-out exists to
    prevent.
    zscore: z-score each contributor before combining (default True) so
    amplitude differences don't let one dominate.
    threshold_abs: if given, consensus regions whose |value| is below
    this are set to NaN. None (default) keeps every region.

    Returns (consensus_values, used_labels).
    """
    exclude_species = None
    if species_of_label and exclude_label is not None:
        exclude_species = species_of_label.get(exclude_label)

    # 1) per-BIDS vectors, z-scored if requested
    by_label = {}
    for label, prof in sorted(species_profiles.items()):
        if exclude_species is not None:
            if species_of_label.get(label) == exclude_species:
                continue
        elif exclude_label is not None and label == exclude_label:
            continue
        vals = _ordered_common_values(prof, common_base_names, value_col)
        if not np.isfinite(vals).any():
            continue
        by_label[label] = _normalize_values(vals, zscore)

    if not by_label:
        return np.full(len(common_base_names), np.nan), []

    # 2) collapse BIDS of the same species, so each species weighs once
    if species_of_label:
        grouped = {}
        for label, vec in by_label.items():
            grouped.setdefault(species_of_label.get(label, label), []).append(vec)
        with np.errstate(invalid="ignore"):
            contributors = [np.nanmedian(np.vstack(v), axis=0) for v in grouped.values()]
        contributor_names = sorted(grouped)
    else:
        contributors = list(by_label.values())
        contributor_names = sorted(by_label)

    arr = np.vstack(contributors)
    with np.errstate(invalid="ignore"):
        consensus = np.nanmedian(arr, axis=0)

    if threshold_abs is not None:
        n_before = int(np.isfinite(consensus).sum())
        consensus = np.where(np.abs(consensus) >= threshold_abs, consensus, np.nan)
        n_after = int(np.isfinite(consensus).sum())
        if verbose:
            print(f"    threshold_abs={threshold_abs}: kept {n_after}/{n_before} region(s)")

    if verbose:
        if exclude_species is not None:
            excl = f" (excluding species {exclude_species})"
        elif exclude_label is not None:
            excl = f" (excluding {exclude_label})"
        else:
            excl = ""
        unit = "species" if species_of_label else "BIDS"
        print(f"    consensus from {len(contributor_names)} {unit}{excl}: {contributor_names}")
    return consensus, sorted(by_label)


def compute_loo_species_distances(profiles_by_label, common_base_names, value_col="median_r",
                                  zscore=True, verbose=True):
    """
    PHYLOGENY-SIDE ONLY (feeds F6 and PGLS's distance_to_consensus). Has
    NO leave_one_out parameter -- LOO is unconditional here, and this is
    NOT controlled by cfg["LEAVE_ONE_OUT"] (that flag only affects
    selection, see build_selection_target above). The two are
    deliberately independent: toggling selection's LOO never changes
    what F6/PGLS see, and vice versa. Do not disable this one to "match"
    a selection-side experiment -- it is what makes the phylogenetic
    signal claim non-circular in the first place; there is no
    engineering-tradeoff argument for turning it off the way there is
    for selection.

    For every label (BIDS) in profiles_by_label ({label: (species,
    profile_df)}), correlate its profile against a leave-one-SPECIES-out
    consensus (build_consensus_pattern with species_of_label, which
    excludes every BIDS of that label's own species, not just the label
    itself -- see build_consensus_pattern's docstring for why a
    same-species BIDS left in would still leak). This is what makes a
    species' own distance non-circular regardless of how many BIDS it
    has, so nothing needs to be excluded from the fit afterward.

    Returns a DataFrame, one row per BIDS: label, species,
    distance_to_human (1 - r), pearson_r, n_regions, n_species_in_ref.
    ("distance_to_human" kept as the column name for compatibility with
    plot_distance_vs_divergence/test_phylogenetic_signal, even though
    the reference is a leave-one-species-out cross-species consensus,
    not literally human.)
    """
    dist_rows = []
    species_of_label = {lab: sp for lab, (sp, p) in profiles_by_label.items()}
    profiles_only = {lab: p for lab, (sp, p) in profiles_by_label.items()}
    for lab, (sp, p) in sorted(profiles_by_label.items()):
        v = _ordered_common_values(p, common_base_names, value_col)
        ref_loo, used = build_consensus_pattern(
            profiles_only, common_base_names, value_col=value_col,
            exclude_label=lab, zscore=zscore, threshold_abs=None,
            species_of_label=species_of_label, verbose=False)
        both = np.isfinite(v) & np.isfinite(ref_loo)
        if both.sum() < 3:
            if verbose:
                print(f"  [SKIP] {lab}: only {both.sum()} finite common region(s) "
                      f"vs its leave-one-species-out reference, need >= 3")
            continue
        rr = float(np.corrcoef(v[both], ref_loo[both])[0, 1])
        dist_rows.append(dict(label=lab, species=sp, distance_to_human=1 - rr,
                              pearson_r=rr, n_regions=int(both.sum()),
                              n_species_in_ref=len(used)))
    return pd.DataFrame(dist_rows)


def compute_run_distance_to_reference(runs_by_label, ref_profile_by_label, common_base_names,
                                      value_col="median_r", kept_keys_by_label=None, verbose=True):
    """
    For every label with a reference profile (ref_profile_by_label --
    typically each BIDS' "kept" group-map profile), correlate EVERY run
    in runs_by_label[label] against it individually -- kept, not-kept,
    and anything outside the balanced 3dttest++ alike, since the point
    is to see the full spread each BIDS actually produced, not just the
    runs that made the cut.

    kept_keys_by_label: optional {label: set-of-run-keys} (e.g. the
    step-2 `selections` dict) -- when given, an `is_kept` bool column
    marks which rows were in the selected/kept subset, so a plot can
    show where they sit within the full spread.

    Returns a long DataFrame: label, species, subject, session, run,
    seed, pearson_r, distance (1-r), n_regions, is_kept.
    """
    rows = []
    for label, (species, ref_profile) in ref_profile_by_label.items():
        profs = runs_by_label.get(label)
        if not profs:
            continue
        kept_keys = {tuple(str(x) for x in k) for k in
                    (kept_keys_by_label or {}).get(label, [])}
        ref_v = _ordered_common_values(ref_profile, common_base_names, value_col)
        for run_key, run_profile in profs.items():
            v = _ordered_common_values(run_profile, common_base_names, value_col)
            both = np.isfinite(v) & np.isfinite(ref_v)
            if both.sum() < 3:
                continue
            rr = float(np.corrcoef(v[both], ref_v[both])[0, 1])
            rows.append(dict(label=label, species=species,
                             subject=run_key[0], session=run_key[1],
                             run=run_key[2], seed=run_key[3],
                             pearson_r=rr, distance=1 - rr,
                             n_regions=int(both.sum()),
                             is_kept=tuple(str(x) for x in run_key) in kept_keys))
    result = pd.DataFrame(rows)
    if verbose:
        print(f"  compute_run_distance_to_reference: {len(result)} run(s) scored "
              f"across {result['label'].nunique() if not result.empty else 0} BIDS")
    return result


def _qc_norm_id(x):
    """Session/run id normalisation for the QC join: '01' == '1', and
    'NA' / '' / NaN (pandas turns the string 'NA' into NaN when reading an
    Excel file, even with dtype=str) all become 'NA'."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "NA"
    t = str(x).strip()
    if t in ("", "nan", "NaN", "None", "NA", "<NA>"):
        return "NA"
    return t.lstrip("0") or "0"


def correlate_qc_with_distance(dist_df, df_qc_by_label, qc_metrics=None,
                               id_cols=("subject", "session"), verbose=True):
    """
    Within each BIDS, Pearson r between each QC metric and each unit's
    distance to its dataset reference (dist_df: one row per run).

    Granularity is decided PER METRIC, from the QC table itself:
      - metric constant within every (subject, session) -> session-level
        metric (e.g. tSNR/FD computed once per session and repeated on
        each run's row): QC and distance are both reduced to ONE row per
        session (distance = mean over that session's runs);
      - metric varying within sessions AND run ids matching >= 90% of the
        distance rows -> run-level join (one row per run);
      - otherwise -> session-level (QC averaged per session).

    Confirmed bug this replaces: distance was averaged per session, then
    joined on (subject, session) to a QC table with one row PER RUN, so
    every session was duplicated by its number of runs (real log:
    "matched 40/4", "169/57", "195/46", "52/31"). n, hence every p-value
    and every meta-analytic weight (n - 3), was inflated by pseudo-
    replication. A merge can no longer produce more rows than its left
    side (asserted).

    Returns a long DataFrame: label, metric, level ("run"/"session"), n,
    pearson_r, pearson_p.
    """
    if dist_df is None or dist_df.empty or not df_qc_by_label:
        if verbose:
            print("  [SKIP] correlate_qc_with_distance: no distance data or no QC table")
        return pd.DataFrame()

    no_qc = sorted(set(dist_df["label"]) - set(df_qc_by_label))
    if no_qc and verbose:
        print(f"  [WARN] no QC table for {no_qc} -- excluded from the QC analysis")

    rows = []
    for label, df_qc in df_qc_by_label.items():
        d = dist_df[dist_df["label"] == label].copy()
        if d.empty or df_qc is None or df_qc.empty or "subject" not in df_qc.columns:
            continue
        q = df_qc.copy()
        for df_ in (d, q):
            df_["_sub"] = df_["subject"].astype(str).str.strip()
            df_["_ses"] = (df_["session"].map(_qc_norm_id) if "session" in df_.columns else "NA")
            df_["_run"] = (df_["run"].map(_qc_norm_id) if "run" in df_.columns else "NA")
        d_ses = (d.groupby(["_sub", "_ses"], as_index=False)
                  .agg(distance=("distance", "mean"), n_runs=("run", "nunique")))

        run_join = False
        if "run" in q.columns and (q["_run"] != "NA").any():
            kd = set(zip(d["_sub"], d["_ses"], d["_run"]))
            kq = set(zip(q["_sub"], q["_ses"], q["_run"]))
            run_join = len(kd & kq) >= 0.9 * len(kd)

        metrics = qc_metrics or [c for c in q.columns if
                                 c.startswith(("func_", "anat_", "net_", "homotopic_",
                                               "cross_", "sp_", "contam_"))
                                 and pd.api.types.is_numeric_dtype(q[c])]
        n_by_level = {}
        for metric in metrics:
            if metric not in q.columns:
                continue
            qv = q[["_sub", "_ses", "_run", metric]].copy()
            qv[metric] = pd.to_numeric(qv[metric], errors="coerce")
            per_ses_unique = qv.dropna(subset=[metric]).groupby(["_sub", "_ses"])[metric].nunique()
            session_level = per_ses_unique.empty or bool((per_ses_unique <= 1).all())
            if run_join and not session_level:
                qm = qv.groupby(["_sub", "_ses", "_run"], as_index=False)[metric].mean()
                m = d.merge(qm, on=["_sub", "_ses", "_run"], how="inner")
                level, n_left = "run", len(d)
            else:
                qm = qv.groupby(["_sub", "_ses"], as_index=False)[metric].mean()
                m = d_ses.merge(qm, on=["_sub", "_ses"], how="inner")
                level, n_left = "session", len(d_ses)
            assert len(m) <= n_left, f"{label}/{metric}: QC join duplicated rows"
            n_by_level.setdefault(level, (len(m), n_left))
            x = pd.to_numeric(m["distance"], errors="coerce")
            y = m[metric]
            ok = x.notna() & y.notna()
            if ok.sum() < 4 or y[ok].nunique() < 2:
                continue
            r, p = _stats.pearsonr(x[ok], y[ok])
            rows.append(dict(label=label, metric=metric, level=level, n=int(ok.sum()),
                             pearson_r=float(r), pearson_p=float(p)))
        if verbose:
            desc = ", ".join(f"{lvl}-level: {a}/{b} matched" for lvl, (a, b) in n_by_level.items())
            print(f"  [{label}] QC join -- {desc or 'no metric'}")

    result = pd.DataFrame(rows)
    if verbose and not result.empty:
        print(result.sort_values("pearson_p").head(30).to_string(index=False))
    return result


def meta_analyze_qc_correlations(per_bids_corr_df, min_bids=3, alpha=0.05, verbose=True):
    """
    Combines correlate_qc_with_distance's per-BIDS rows into ONE number
    per metric, WITHOUT ever pooling raw rows across BIDS -- a standard
    fixed-effect correlation meta-analysis: Fisher-z each BIDS's r
    (z_i = arctanh(r_i)), weight by its precision (n_i - 3), and test
    whether the weighted-average z is consistently non-zero across BIDS.
    This is what "significant in inter-species" should mean here: not a
    big pooled r (which the between-BIDS confound can produce on its
    own, see correlate_qc_with_distance's docstring), but a QC-distance
    relationship that shows up with the same sign across independent
    datasets.

    min_bids: a metric needs at least this many BIDS with a computed r
    to even attempt the meta-test (3 is the same floor used elsewhere in
    this pipeline for phylogenetic signal -- below that, one dataset can
    flip the result and "significant" is not a meaningful claim).

    Returns a DataFrame: metric, n_bids, meta_r, meta_p, n_same_sign
    (out of n_bids), sorted by meta_p. Only rows with n_bids >= min_bids
    are returned at all.
    """
    if per_bids_corr_df is None or per_bids_corr_df.empty:
        if verbose:
            print("  [SKIP] meta_analyze_qc_correlations: nothing to combine")
        return pd.DataFrame()

    rows = []
    for metric, g in per_bids_corr_df.groupby("metric"):
        g = g.dropna(subset=["pearson_r", "n"])
        g = g[g["n"] > 3]
        if len(g) < min_bids:
            continue
        z = np.arctanh(g["pearson_r"].clip(-0.999, 0.999).values)
        w = (g["n"].values - 3).astype(float)
        z_meta = float(np.sum(w * z) / np.sum(w))
        se_meta = float(1.0 / np.sqrt(np.sum(w)))
        z_stat = z_meta / se_meta
        p_meta = float(2 * _stats.norm.sf(abs(z_stat)))
        r_meta = float(np.tanh(z_meta))
        n_same_sign = int((np.sign(g["pearson_r"].values) == np.sign(r_meta)).sum())
        rows.append(dict(metric=metric, n_bids=len(g), meta_r=r_meta, meta_p=p_meta,
                         n_same_sign=n_same_sign))

    result = pd.DataFrame(rows)
    if result.empty:
        if verbose:
            print(f"  [meta] no metric had >= {min_bids} BIDS with a computed correlation")
        return result
    result = result.sort_values("meta_p")
    # Benjamini-Hochberg across metrics: ~40 QC metrics are tested, so at
    # alpha=0.05 about 2 would pass by chance with no correction.
    p = result["meta_p"].values
    order = np.argsort(p)
    q = np.empty(len(p))
    running = 1.0
    for rank, idx in reversed(list(enumerate(order, start=1))):
        running = min(running, p[idx] * len(p) / rank)
        q[idx] = running
    result["meta_q_fdr"] = q
    n_sig = int((result["meta_q_fdr"] < alpha).sum())
    if verbose:
        print(f"  [meta] {n_sig}/{len(result)} metric(s) significant at q_FDR<{alpha} "
              f"(BH across metrics) across >= {min_bids} BIDS "
              f"({int((result['meta_p'] < alpha).sum())} at uncorrected p<{alpha})")
        print(result.head(20).to_string(index=False))
    return result


def compute_region_invariance(species_profiles, common_base_names, value_col="median_r",
                              species_of_label=None, zscore=True, verbose=True):
    """
    How INVARIANT is each region across species?

    Motivation (a real flaw in selecting on the plain consensus): the
    consensus is a median across species, so a region where species
    genuinely disagree still contributes to it -- at whatever value the
    median happens to land on. Runs then get selected for matching an
    average of disagreement, which is noise, not conserved biology. Only
    regions where species actually AGREE carry information about what is
    conserved.

    Invariance is measured as 1 - (dispersion across species), where
    dispersion is the median absolute deviation of the per-species values
    around their median, rescaled to [0, 1] across regions. A region
    where every species reports the same z-scored value scores near 1; a
    region where they scatter widely scores near 0.

    Returns (invariance, per_species_matrix, species_names):
    invariance is aligned to common_base_names.
    """
    by_species = {}
    for label, prof in species_profiles.items():
        vals = _ordered_common_values(prof, common_base_names, value_col)
        if not np.isfinite(vals).any():
            continue
        vec = _normalize_values(vals, zscore)
        sp = (species_of_label or {}).get(label, label)
        by_species.setdefault(sp, []).append(vec)

    if not by_species:
        return np.full(len(common_base_names), np.nan), None, []

    with np.errstate(invalid="ignore"):
        per_species = {sp: np.nanmedian(np.vstack(v), axis=0) for sp, v in by_species.items()}
    names = sorted(per_species)
    mat = np.vstack([per_species[s] for s in names])

    with np.errstate(invalid="ignore"):
        med = np.nanmedian(mat, axis=0)
        mad = np.nanmedian(np.abs(mat - med), axis=0)

    finite = np.isfinite(mad)
    invariance = np.full(len(common_base_names), np.nan)
    if finite.sum():
        lo, hi = np.nanmin(mad[finite]), np.nanmax(mad[finite])
        spread = (mad - lo) / (hi - lo) if hi > lo else np.zeros_like(mad)
        invariance[finite] = 1.0 - spread[finite]

    if verbose and finite.sum():
        order = np.argsort(-np.nan_to_num(invariance, nan=-1))
        top = [common_base_names[i] for i in order[:5]]
        bot = [common_base_names[i] for i in order[-5:]]
        print(f"    invariance across {len(names)} species: "
             f"most invariant {top}, least {bot}")
    return invariance, mat, names


def select_invariant_regions(invariance, common_base_names, mode="top_fraction",
                             value=0.5, verbose=True):
    """
    Turn the invariance scores into the region subset selection runs on.

    mode="top_fraction": keep the most invariant `value` fraction
        (0.5 = the top half). Scale-free, so it adapts to however many
        common regions exist -- preferred, since an absolute cutoff on a
        rescaled score has no fixed meaning.
    mode="threshold": keep regions with invariance >= value.
    mode="all": keep everything (the previous behaviour, kept so the
        effect of restricting can be measured rather than assumed).

    Returns a boolean mask aligned to common_base_names.
    """
    inv = np.asarray(invariance, dtype=float)
    finite = np.isfinite(inv)
    if mode == "all" or not finite.any():
        mask = finite.copy() if finite.any() else np.ones(len(inv), bool)
    elif mode == "top_fraction":
        frac = float(value)
        if frac >= 1.0:
            # Explicit "keep everything": every region with a finite
            # score, no ranking applied. Must not go through the max(3, .)
            # floor below, which would be a no-op here but would also hide
            # the intent.
            mask = finite.copy()
        elif frac <= 0.0:
            # Explicit "keep nothing". The floor of 3 exists to stop a
            # careless fraction from leaving too few regions to correlate
            # on -- it must not override an explicit 0, otherwise the knob
            # silently refuses its own endpoint (confirmed by test: 0.0
            # used to return every region).
            mask = np.zeros(len(inv), bool)
            print("    [WARN] value=0 -- no region is used for the run "
                  "selection; scoring runs needs at least 3 regions, so "
                  "STEP 2 will not be able to rank anything.")
        else:
            k = max(3, int(round(finite.sum() * frac)))
            order = np.argsort(-np.where(finite, inv, -np.inf))
            mask = np.zeros(len(inv), bool)
            mask[order[:k]] = True
    elif mode == "threshold":
        mask = finite & (inv >= float(value))
        if mask.sum() < 3:
            print(f"    [WARN] threshold {value} leaves only {mask.sum()} region(s) "
                 f"-- falling back to the top 3 most invariant")
            order = np.argsort(-np.where(finite, inv, -np.inf))
            mask = np.zeros(len(inv), bool); mask[order[:3]] = True
    else:
        raise ValueError(f"Unknown mode {mode!r}")

    if verbose:
        kept = [n for n, m in zip(common_base_names, mask) if m]
        print(f"    selecting on {mask.sum()}/{len(mask)} region(s) "
             f"(mode={mode}, value={value}): {kept}")
    return mask


def plot_region_invariance(invariance, common_base_names, mask=None, title=None,
                           output_path=None, dpi=300):
    """Bar chart of per-region invariance, regions used for selection
    highlighted. Makes explicit which regions the selection actually
    rests on."""
    inv = np.asarray(invariance, dtype=float)
    order = np.argsort(-np.nan_to_num(inv, nan=-1))
    names = [common_base_names[i] for i in order]
    vals = inv[order]
    used = (np.asarray(mask, bool)[order] if mask is not None
            else np.ones(len(vals), bool))

    fig, ax = plt.subplots(figsize=(max(8, len(names) * 0.3), 4.5))
    ax.bar(range(len(vals)), vals,
          color=["#009E73" if u else "#CCCCCC" for u in used], edgecolor="white")
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=90, fontsize=6)
    ax.set_ylabel("invariance across species", fontsize=10)
    ax.set_title(title or "Region invariance (green = used for selection)",
                fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def rank_runs_against_consensus(sri, profiles, consensus_values, common_base_names,
                                value_col="median_r", zscore=True, similarity="correlation"):
    """
    Score every run's profile against the consensus vector.

    similarity controls WHAT is being matched:
      "correlation" (default) -- Pearson r: matches the SHAPE of the
          profile and is completely blind to amplitude. Note that
          Pearson is invariant to any affine rescaling, so the zscore
          argument has NO effect on this metric (verified: r is
          identical to 6 decimals with and without z-scoring, and after
          multiplying one side by 100). zscore only really matters for
          building the consensus itself, where the median across species
          IS scale-sensitive.
      "cosine" -- cosine similarity WITHOUT mean-centering: unlike
          Pearson it is sensitive to sign/offset, so a profile that is
          uniformly positive scores differently from one centered on
          zero. Still amplitude-invariant (it normalizes vector length).
      "neg_euclidean" -- negative Euclidean distance on the (optionally
          z-scored) values: the ONLY option here that is sensitive to
          AMPLITUDE, so a run with the right shape but a much weaker
          signal scores worse than one matching in both. Use it when
          "weak but correctly shaped" should NOT count as equivalent to
          "strong and correctly shaped". With zscore=True this
          amplitude-sensitivity is removed again per run, so pair it
          with zscore=False if amplitude is what you want to capture.

    Returns a DataFrame sorted best-first: subject, session, run, seed,
    pearson_r (the score under the chosen metric -- named this way for
    backward compatibility with callers), n_regions, _key.
    """
    rows = []
    for key, prof in profiles.items():
        vals = _ordered_common_values(prof, common_base_names, value_col)
        vals = _normalize_values(vals, zscore)
        vals = np.asarray(vals, dtype=float)
        both = np.isfinite(vals) & np.isfinite(consensus_values)
        if both.sum() < 3:
            continue
        v, c = vals[both], consensus_values[both]

        if similarity == "correlation":
            score = float(np.corrcoef(v, c)[0, 1])
        elif similarity == "cosine":
            denom = np.linalg.norm(v) * np.linalg.norm(c)
            score = float(np.dot(v, c) / denom) if denom > 0 else np.nan
        elif similarity == "neg_euclidean":
            # negated so that, like the other two, HIGHER is always better
            score = -float(np.sqrt(np.mean((v - c) ** 2)))
        else:
            raise ValueError(f"Unknown similarity {similarity!r} -- use 'correlation', "
                             f"'cosine', or 'neg_euclidean'")

        rows.append(dict(subject=key[0], session=key[1], run=key[2], seed=key[3],
                        pearson_r=score, n_regions=int(both.sum()), _key=key))
    if not rows:
        return pd.DataFrame(columns=["subject", "session", "run", "seed",
                                     "pearson_r", "n_regions", "_key"])
    return pd.DataFrame(rows).sort_values("pearson_r", ascending=False).reset_index(drop=True)


def build_selection_target(species_profiles, common_base_names, label, leave_one_out,
                           zscore, threshold_abs, species_of_label, value_col, verbose):
    """
    SELECTION-SIDE ONLY. Builds the consensus a single BIDS (`label`)
    is scored against for run selection (round 1 and round >= 3).

    Controlled entirely by cfg["LEAVE_ONE_OUT"] (selection's own flag,
    default True). Deliberately kept separate from compute_loo_species_distances
    below (the phylogeny-side function, used by F6/PGLS): the two answer
    different questions and have different stakes if disabled --
    disabling LOO here is a selection/engineering choice with a small,
    symmetric effect (every species gets the same ~1-vote self-nudge on
    its own median, real but bounded, since with n=7 species dropping
    one from a median IS a non-trivial change, not just 1/7 dilution --
    it matters most for species that are themselves outliers, exactly
    where an independent check is most valuable). Disabling LOO on the
    phylogeny side would remove exactly the guard that makes F6/PGLS's
    "phylogenetic signal" claim non-circular -- these are not
    interchangeable levers, hence two functions, not one flag.

    Returns (consensus_values, used_labels) -- same contract as
    build_consensus_pattern.
    """
    return build_consensus_pattern(
        species_profiles, common_base_names, value_col=value_col,
        exclude_label=label if leave_one_out else None,
        zscore=zscore, threshold_abs=threshold_abs,
        species_of_label=species_of_label, verbose=verbose)


def _select_against_species_profiles(sri, species_profiles, profiles_by_label,
                                     common_base_names, n_select, leave_one_out,
                                     zscore, threshold_abs, similarity, species_of_label,
                                     region_mask, value_col, metric, verbose):
    """
    Shared core of run_consensus_selection and run_iterative_selection:
    given a {label: profile} dict to build the LOO cross-species
    consensus FROM, rank and select the top n_select of EVERY run in
    profiles_by_label[label] (the full, un-truncated pool) against it.
    Factored out so round >=2 of the iterative loop scores against a
    consensus built from a DIFFERENT species_profiles dict (kept-only)
    using the exact same consensus/ranking code as round 1, rather than
    a re-derived copy that could silently drift out of sync.
    """
    selections, rankings = {}, {}
    consensus_by_label, used_labels_by_label = {}, {}
    for label, profiles in profiles_by_label.items():
        if verbose:
            print(f"\n  {label}:")
        consensus, used = build_selection_target(
            species_profiles, common_base_names, label, leave_one_out,
            zscore, threshold_abs, species_of_label, value_col, verbose)
        if region_mask is not None:
            consensus = np.where(region_mask, consensus, np.nan)
        if not np.isfinite(consensus).any():
            print(f"  [SKIP] {label}: consensus has no usable region")
            continue

        ranking = rank_runs_against_consensus(
            sri, profiles, consensus, common_base_names,
            value_col=value_col, zscore=zscore, similarity=similarity)
        if ranking.empty:
            print(f"  [SKIP] {label}: no run could be ranked")
            continue

        selected = ranking.head(n_select)["_key"].tolist()
        selections[label] = selected
        rankings[label] = ranking.drop(columns="_key")
        consensus_by_label[label] = consensus
        used_labels_by_label[label] = used
        if verbose:
            top = ranking.head(n_select)["pearson_r"]
            print(f"    selected {len(selected)} run(s), r = {top.min():.4f} to {top.max():.4f}")
    return selections, rankings, consensus_by_label, used_labels_by_label


def run_consensus_selection(sri, cfg, profiles_by_label, common_base_names,
                            n_select=5, leave_one_out=True, consensus_source="all_runs",
                            method="kcore", k=5, zscore=True, threshold_abs=None,
                            similarity="correlation", species_of_label=None,
                            region_mask=None,
                            value_col=None, metric=None, output_dir=None, verbose=True):
    """
    Select runs in every BIDS against a cross-species consensus.

    profiles_by_label: {label: {key: profile_df}} -- compute these once
    (e.g. the PASS 1 loop in human_anchored_phylogeny's own driver) and
    pass them in, so nothing is recomputed here.
    common_base_names: the FIXED cross-species region list (step 1).

    Returns dict(selections={label: [key,...]}, rankings={label: df},
    consensus_by_label={label: consensus vector actually used for it},
    species_profiles={label: profile}, used_labels_by_label={label: [...]}).
    """
    # cfg is the pipeline CONFIG dict; the launcher module used to be
    # passed here purely to read these two values off it.
    value_col = value_col or cfg["VALUE_COL"]
    metric = metric or cfg["METRIC"]

    print(f"\n{'=' * 70}")
    print(f"  Consensus selection (leave_one_out={leave_one_out}, "
          f"source={consensus_source}, zscore={zscore}, threshold_abs={threshold_abs}, "
          f"similarity={similarity})")
    print(f"{'=' * 70}")

    species_profiles = {}
    for label, profiles in profiles_by_label.items():
        try:
            species_profiles[label] = species_profile(
                sri, profiles, source=consensus_source, method=method, k=k,
                value_col=value_col, metric=metric)
        except Exception as e:
            print(f"  [WARN] {label}: could not build its species profile: {e}")

    # region_mask restricts BOTH the consensus and the run scoring to the
    # regions that are actually invariant across species. Without it, runs
    # are ranked partly on regions where species genuinely disagree, i.e.
    # on cross-species noise rather than conserved signal.
    if region_mask is not None:
        region_mask = np.asarray(region_mask, dtype=bool)

    selections, rankings, consensus_by_label, used_labels_by_label = (
        _select_against_species_profiles(
            sri, species_profiles, profiles_by_label, common_base_names,
            n_select, leave_one_out, zscore, threshold_abs, similarity,
            species_of_label, region_mask, value_col, metric, verbose))

    if output_dir and rankings:
        os.makedirs(output_dir, exist_ok=True)
        tag = (f"loo{int(leave_one_out)}_{consensus_source}_{similarity}"
               f"_z{int(zscore)}_thr{threshold_abs if threshold_abs is not None else 'none'}")
        excel_path = opj(output_dir, f"consensus_selection_{tag}.xlsx")
        with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
            summary = pd.DataFrame([
                dict(label=lab, n_runs=len(profiles_by_label[lab]),
                     n_selected=len(selections.get(lab, [])),
                     consensus_from=", ".join(used_labels_by_label.get(lab, [])),
                     n_consensus_regions=int(np.isfinite(consensus_by_label[lab]).sum()))
                for lab in rankings
            ])
            summary.to_excel(writer, sheet_name="summary", index=False)
            for lab, rdf in rankings.items():
                rdf.head(100).to_excel(writer, sheet_name=f"rank_{lab[:24]}", index=False)
        print(f"\n  [saved] {excel_path}")

        plot_consensus_pattern(
            consensus_by_label, common_base_names,
            title=f"Consensus pattern per species ({tag})",
            output_path=opj(output_dir, f"consensus_pattern_{tag}.png"))

    return dict(selections=selections, rankings=rankings,
               consensus_by_label=consensus_by_label,
               species_profiles=species_profiles,
               used_labels_by_label=used_labels_by_label)


def run_round2_selection(sri, cfg, profiles_by_label, selections, common_base_names,
                         n_select, region_mask=None, verbose=True):
    """
    Second selection round, run AFTER run_consensus_selection.

    Round 1 anchors every species' selection to the cross-species
    consensus, which can favor runs that resemble OTHER species over
    ones that are simply typical examples of THIS species. Round 2
    re-injects the species' own signal: build a reference from round
    1's kept subset (kcore over those runs only), then re-rank the FULL
    candidate pool -- every run that passed QC, not just round 1's kept
    8 -- against THIS species-anchored target instead of the
    cross-species one.

    A run round 1 rejected for looking insufficiently cross-species-
    typical can be picked back up here if it matches this species' own
    signal well; a run round 1 kept can drop out if it turns out to be
    atypical FOR THE SPECIES despite resembling the consensus.

    On circularity: round 2's target is built from round 1's own kept
    subset of the SAME pool being re-ranked, so it is not fully
    independent of that pool -- but it is considerably less circular
    than ranking from scratch against nothing, since round 1's kept
    subset was itself gated by an external, cross-species anchor. Treat
    round 2 as a deliberately species-anchored ALTERNATIVE view, not a
    strictly "more correct" one -- this pipeline reports both rather
    than replacing round 1's outputs.

    Returns dict(selections={label: [key,...]},
    round1_kept_ref={label: profile_df} -- round 2's scoring target,
    rankings={label: df}).
    """
    selections_r2, round1_kept_ref, rankings_r2, criterion_by_label = {}, {}, {}, {}
    for label, profs in profiles_by_label.items():
        kept_keys_r1 = {tuple(str(x) for x in k) for k in selections.get(label, [])}
        kept_profs_r1 = {k: v for k, v in profs.items()
                         if tuple(str(x) for x in k) in kept_keys_r1}
        if len(kept_profs_r1) < 2:
            if verbose:
                print(f"  [{label}] round 2 skipped: only {len(kept_profs_r1)} "
                      f"round-1 kept run(s), need >= 2 to build a reference")
            continue
        try:
            kept_ref_r1, _ = build_reference_for_profiles(
                sri, kept_profs_r1, method=cfg["REFERENCE_METHOD"],
                k=min(cfg["K_REF"], len(kept_profs_r1)),
                value_col=cfg["VALUE_COL"], metric=cfg["METRIC"])
        except Exception as e:
            print(f"  [WARN] {label}: round-1 kept reference failed: {e}")
            continue
        round1_kept_ref[label] = kept_ref_r1

        # Target on the same footing as round 1's: normalised over ALL common
        # regions with the same mode as the run vectors, THEN masked. The
        # previous version hardcoded similarity="correlation" and left the
        # target raw, so SIMILARITY="cosine"/"neg_euclidean" was silently
        # ignored in round 2 (and a raw target is not comparable to
        # normalised run vectors for the euclidean metric).
        consensus_values = _normalize_values(
            _ordered_common_values(kept_ref_r1, common_base_names, cfg["VALUE_COL"]),
            cfg["ZSCORE"])
        if region_mask is not None:
            consensus_values = np.where(np.asarray(region_mask, dtype=bool),
                                        consensus_values, np.nan)
        criterion_by_label[label] = consensus_values
        ranking = rank_runs_against_consensus(
            sri, profs, consensus_values, common_base_names,
            value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"], similarity=cfg["SIMILARITY"])
        if ranking.empty:
            print(f"  [SKIP] {label}: round 2 -- no run could be ranked")
            continue
        n = min(n_select, len(ranking))
        top = ranking.head(n)
        selections_r2[label] = top["_key"].tolist()
        rankings_r2[label] = ranking.drop(columns="_key")
        if verbose:
            overlap = len({tuple(str(x) for x in k) for k in top["_key"]} & kept_keys_r1)
            print(f"    round 2: selected {n} run(s), r = {top['pearson_r'].min():.4f} "
                  f"to {top['pearson_r'].max():.4f}, {overlap}/{n} overlap with round-1 kept")
    return dict(selections=selections_r2, round1_kept_ref=round1_kept_ref,
               rankings=rankings_r2, criterion_by_label=criterion_by_label)


def run_iterative_selection(sri, cfg, profiles_by_label, common_base_names, n_select,
                            start_selections, start_round, n_rounds,
                            species_of_label=None, region_mask=None, verbose=True):
    """
    Continues the selection chain from `start_round` (its kept runs are
    `start_selections`) up to round `n_rounds` inclusive. Each round r:
      (a) rebuild each BIDS's representative profile from round r-1's
          kept runs (build_reference_for_profiles, same
          REFERENCE_METHOD/K_REF/METRIC as everywhere else);
      (b) rebuild the cross-species consensus from these kept-only
          profiles (build_selection_target: LOO follows
          cfg["LEAVE_ONE_OUT"], exactly like round 1);
      (c) re-rank the FULL original pool (a run rejected earlier can come
          back) against that consensus; keep the top n_select.

    Chain the pipeline builds: round 1 = cross-species consensus from all
    runs; round 2 (RUN_ROUND2_SELECTION) = species-anchored, each BIDS
    scored against its own round-1 kept mean; rounds 3..N (this
    function) = cross-species consensus rebuilt from the previous
    round's kept means. When round 2 is off, start_round=1 and this
    function produces rounds 2..N with the same rule.

    Confirmed bug this replaces: the previous version always restarted
    from round 1's selection even when round 2 had run, so a run with
    both options produced TWO different "round 2"s (species-anchored in
    03b maps, cross-species in F1e) and the loop never used kept2 at all.

    ON CIRCULARITY: each round moves one step further from the only
    external anchor (round 1's all-runs consensus). The loop converges on
    a self-consistent pattern; it cannot tell a shared biological signal
    from a shared artefact among the kept runs. overlap_trajectory
    (plotted in F1e) shows whether it actually converged.

    Returns dict(rounds={r: {"selections", "species_profiles",
    "consensus_by_label"}}, overlap_trajectory=DataFrame(round, label,
    n_select, n_overlap_with_previous)). rounds contains only the rounds
    computed here (start_round + 1 .. n_rounds).
    """
    value_col, metric = cfg["VALUE_COL"], cfg["METRIC"]
    rounds, overlap_rows = {}, []
    prev_selections = start_selections or {}

    for r in range(start_round + 1, n_rounds + 1):
        if verbose:
            print(f"\n  --- round {r} (cross-species consensus from round {r - 1} kept) ---")
        species_profiles_r = {}
        for label, profs in profiles_by_label.items():
            kept_keys = {tuple(str(x) for x in k) for k in prev_selections.get(label, [])}
            kept_profs = {k: v for k, v in profs.items()
                         if tuple(str(x) for x in k) in kept_keys}
            if len(kept_profs) < 2:
                if verbose:
                    print(f"  [{label}] round {r}: only {len(kept_profs)} kept run(s) "
                          f"from round {r - 1} -- excluded from this round's consensus")
                continue
            try:
                species_profiles_r[label], _ = build_reference_for_profiles(
                    sri, kept_profs, method=cfg["REFERENCE_METHOD"],
                    k=min(cfg["K_REF"], len(kept_profs)),
                    value_col=value_col, metric=metric)
            except Exception as e:
                print(f"  [WARN] {label}: round {r} reference failed: {e}")

        selections_r, _, consensus_r, _ = _select_against_species_profiles(
            sri, species_profiles_r, profiles_by_label, common_base_names,
            n_select, cfg["LEAVE_ONE_OUT"], cfg["ZSCORE"], cfg["THRESHOLD_ABS"],
            cfg["SIMILARITY"], species_of_label, region_mask, value_col, metric,
            verbose=False)

        for label, keys in selections_r.items():
            cur = {tuple(str(x) for x in k) for k in keys}
            prev = {tuple(str(x) for x in k) for k in prev_selections.get(label, [])}
            overlap_rows.append(dict(round=r, label=label, n_select=len(cur),
                                     n_overlap_with_previous=len(cur & prev)))
            if verbose:
                print(f"  [{label}] round {r}: {len(cur & prev)}/{len(cur)} "
                      f"unchanged from round {r - 1}")
        missing = sorted(set(prev_selections) - set(selections_r))
        if missing:
            print(f"  [WARN] round {r}: no selection for {missing} -- they keep "
                  f"their round {r - 1} selection downstream")
            for label in missing:
                selections_r[label] = prev_selections[label]

        rounds[r] = dict(selections=selections_r, species_profiles=species_profiles_r,
                         consensus_by_label=consensus_r)
        prev_selections = selections_r

    return dict(rounds=rounds, overlap_trajectory=pd.DataFrame(overlap_rows))


def plot_iterative_selection_convergence(overlap_df, output_path=None, dpi=300):
    """One line per BIDS: fraction of the n_select kept runs unchanged
    from the previous round, across rounds. Flat at 1.0 = converged;
    still moving at the last round = n_rounds was too small to trust,
    or this BIDS is oscillating (also worth knowing)."""
    if overlap_df is None or overlap_df.empty:
        print("  [plot] nothing to plot")
        return None
    df = overlap_df.copy()
    df["frac_overlap"] = df["n_overlap_with_previous"] / df["n_select"]
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = _palette_color_map(df["label"].unique())
    for label, g in df.groupby("label"):
        g = g.sort_values("round")
        ax.plot(g["round"], g["frac_overlap"], marker="o", color=colors[label], label=label)
    ax.set_xlabel("round", fontsize=10)
    ax.set_ylabel("fraction of kept runs unchanged\nfrom previous round", fontsize=9)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xticks(sorted(df["round"].unique()))
    ax.set_title("Iterative selection: convergence by round", fontsize=11, fontweight="bold")
    ax.legend(fontsize=7, ncol=2)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def compare_reference_methods_on_round2(sri, cfg, profiles_by_label, selections,
                                        common_base_names, n_select, region_mask=None,
                                        verbose=True):
    """
    Runs round-2 selection TWICE against the SAME round-1 kept pool --
    once with REFERENCE_METHOD="kcore", once with "uniform_median" --
    and reports how much the selected 8 runs actually change.

    This is the one place in the pipeline where REFERENCE_METHOD already
    has a direct effect on which runs get selected: round 2's target
    (round1_kept_ref) is built with build_reference_for_profiles(method=...).
    Round 1's own selection does NOT depend on REFERENCE_METHOD under the
    default CONSENSUS_SOURCE="all_runs" (species_profile only reads
    `method` when source="reference") -- so this comparison is scoped to
    round 2, not round 1, unless CONSENSUS_SOURCE is changed too.

    Pure Python, no AFNI calls, cheap to run both.

    Returns a DataFrame: label, n_select, n_overlap (out of n_select),
    kcore_only / uniform_only (run keys picked by one method but not
    the other).
    """
    results = {}
    for method in ("kcore", "uniform_median"):
        cfg_m = dict(cfg)
        cfg_m["REFERENCE_METHOD"] = method
        r = run_round2_selection(sri, cfg_m, profiles_by_label, selections,
                                 common_base_names, n_select=n_select,
                                 region_mask=region_mask, verbose=False)
        results[method] = r["selections"]
        if verbose:
            print(f"  [{method}] round-2 selection done for {len(r['selections'])} BIDS")

    rows = []
    for label in profiles_by_label:
        kcore_keys = {tuple(str(x) for x in k) for k in results["kcore"].get(label, [])}
        unif_keys = {tuple(str(x) for x in k)
                    for k in results["uniform_median"].get(label, [])}
        if not kcore_keys and not unif_keys:
            continue
        overlap = len(kcore_keys & unif_keys)
        n_sel = max(len(kcore_keys), len(unif_keys))
        rows.append(dict(label=label, n_select=n_sel, n_overlap=overlap,
                         kcore_only=sorted(kcore_keys - unif_keys),
                         uniform_only=sorted(unif_keys - kcore_keys)))
    result = pd.DataFrame(rows)
    if verbose and not result.empty:
        print(result[["label", "n_select", "n_overlap"]].to_string(index=False))
    return result


def plot_reference_method_overlap(overlap_df, output_path=None, dpi=300):
    """
    Bar chart: for each BIDS, how many of the round-2 selected runs
    agree between REFERENCE_METHOD="kcore" and "uniform_median" (out of
    n_select). A full-height bar means the method choice doesn't matter
    for that BIDS; a short bar means it changes which runs get kept.
    """
    if overlap_df is None or overlap_df.empty:
        print("  [plot] nothing to plot")
        return None
    df = overlap_df.sort_values("n_overlap")
    fig, ax = plt.subplots(figsize=(max(6, 0.7 * len(df)), 5))
    colors = ["#D55E00" if r.n_overlap < r.n_select else "#0072B2" for r in df.itertuples()]
    ax.bar(df["label"], df["n_overlap"], color=colors, edgecolor="black")
    for i, row in enumerate(df.itertuples()):
        ax.text(i, row.n_overlap + 0.1, f"{row.n_overlap}/{row.n_select}",
               ha="center", va="bottom", fontsize=8)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["label"], rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("Runs agreeing between kcore and uniform_median\n(round-2 selection)",
                 fontsize=9)
    ax.set_title("Effect of REFERENCE_METHOD on round-2 run selection",
                fontsize=11, fontweight="bold")
    ax.set_ylim(0, df["n_select"].max() + 1)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig



# === Statistics ==========================================================

def _label_rng(random_seed, label):
    """One independent, reproducible RNG per dataset: seeded by the global
    seed AND the dataset label, so a dataset's random draws do not depend on
    its position in BIDS_CONFIGS (a shared RNG consumed label after label
    made results change when the config order changed)."""
    import zlib
    return np.random.default_rng([int(random_seed), zlib.crc32(str(label).encode())])


def balanced_kept_vs_notkept(profiles_by_label, selections, common_base_names,
                             consensus_by_label=None, value_col="median_r",
                             zscore=True, n_draws=1000, random_seed=42,
                             species_of_label=None, verbose=True):
    """
    Compare selected ("kept") runs against non-selected ones, with the
    two groups matched in SIZE.

    Why this is needed with a fixed n_select: keeping 5 of 195 Marmoset
    runs but 5 of 25 Dog runs leaves wildly unequal comparison groups
    (5 vs 190 against 5 vs 20). An unmatched test is then driven largely
    by the size imbalance -- the larger group's median is estimated far
    more precisely, and effect sizes are not comparable across BIDS.
    Here the non-selected group is randomly subsampled to exactly the
    size of the selected group, the comparison is repeated n_draws
    times, and the reported statistic is the distribution across draws.
    That removes the imbalance without discarding data arbitrarily (no
    single subsample decides the answer).

    The compared quantity is each run's correlation to the consensus
    that BIDS was scored against (consensus_by_label). If that isn't
    given, it falls back to correlation with that BIDS's own median
    profile -- but note that fallback is circular for judging a
    selection that was itself made on consensus similarity, so pass
    consensus_by_label whenever you have it.

    Returns a DataFrame, one row per BIDS: n_kept, n_notkept_total,
    median_kept, median_notkept (mean over draws), delta, auc, p_value
    (median across draws), p_frac_sig (fraction of draws reaching
    p<0.05).
    """
    from scipy import stats as _st
    rows = []

    for label, profiles in profiles_by_label.items():
        rng = _label_rng(random_seed, label)
        selected = {tuple(str(x) for x in k) for k in selections.get(label, [])}
        if not selected:
            continue

        ref_vec = None
        if consensus_by_label is not None:
            ref_vec = consensus_by_label.get(label)
        if ref_vec is None:
            all_vals = [_ordered_common_values(p, common_base_names, value_col)
                        for p in profiles.values()]
            with np.errstate(invalid="ignore"):
                ref_vec = np.nanmedian(np.vstack(all_vals), axis=0)
            ref_vec = _normalize_values(ref_vec, zscore)

        kept_scores, notkept_scores = [], []
        for key, prof in profiles.items():
            vals = _ordered_common_values(prof, common_base_names, value_col)
            vals = _normalize_values(vals, zscore)
            vals = np.asarray(vals, dtype=float)
            both = np.isfinite(vals) & np.isfinite(ref_vec)
            if both.sum() < 3:
                continue
            r = float(np.corrcoef(vals[both], ref_vec[both])[0, 1])
            (kept_scores if tuple(str(x) for x in key) in selected
             else notkept_scores).append(r)

        n_k, n_n = len(kept_scores), len(notkept_scores)
        if n_k < 2 or n_n < n_k:
            if verbose:
                print(f"  [SKIP] {label}: kept={n_k}, not-kept={n_n} "
                     f"(need >= 2 kept and at least as many not-kept)")
            continue

        deltas, aucs, pvals = [], [], []
        kept_arr = np.asarray(kept_scores)
        notkept_arr = np.asarray(notkept_scores)
        for _ in range(n_draws):
            draw = rng.choice(notkept_arr, size=n_k, replace=False)
            deltas.append(np.median(kept_arr) - np.median(draw))
            try:
                u, p = _st.mannwhitneyu(kept_arr, draw, alternative="two-sided")
                aucs.append(u / (n_k * n_k))
                pvals.append(p)
            except ValueError:
                continue

        if not pvals:
            continue
        row = dict(label=label,
                  species=(species_of_label or {}).get(label),
                  n_kept=n_k, n_notkept_total=n_n,
                  median_kept=float(np.median(kept_arr)),
                  median_notkept=float(np.median(notkept_arr)),
                  delta_median=float(np.mean(deltas)),
                  auc=float(np.mean(aucs)),
                  n_draws=len(pvals),
                  note=("descriptive only: kept runs were CHOSEN on this very score, "
                        "so any test of kept vs not-kept on it is circular"))
        # No p-value is reported: the kept runs are, by construction, the
        # top-n on the score compared here, so a Mann-Whitney p is guaranteed
        # to be small and tests nothing (circular analysis, Kriegeskorte et
        # al. 2009). delta/AUC remain useful as a description of how strongly
        # the selection separated the two groups.
        rows.append(row)
        if verbose:
            print(f"  {label}: kept={n_k} vs {n_k} drawn from {n_n} | "
                 f"delta={row['delta_median']:+.4f}, AUC={row['auc']:.3f} "
                 f"(descriptive; no p-value -- circular by construction)")

    return pd.DataFrame(rows)


def validate_selection_against_qc(profiles_by_label, selections, df_qc_by_label,
                                  merge_fn, metrics=None, n_draws=500,
                                  random_seed=42, verbose=True):
    """
    The non-circular check: do the SELECTED runs differ from the rest on
    external QC metrics (motion, tSNR, DVARS) that never entered the
    consensus or the selection?

    Groups are size-matched by subsampling the non-selected side to the
    selected size, repeated n_draws times -- with 5 selected out of 195
    vs 5 out of 25, an unmatched test would be driven by group size
    rather than by any real difference (confirmed: identical effects
    produced p-values five orders of magnitude apart when unmatched).

    merge_fn: external_qc_validation.merge_qc_with_classification.

    Returns a DataFrame: label, metric, n_kept, delta (kept - notkept
    mean over draws), auc, p_value (median over draws), p_frac_sig.
    """
    from scipy import stats as _st
    metrics = metrics or ["func_TSNR", "func_mean_fd", "func_mean_dvars"]
    rows = []

    for label, profiles in profiles_by_label.items():
        rng = _label_rng(random_seed, label)
        df_qc = df_qc_by_label.get(label)
        if df_qc is None or df_qc.empty:
            continue
        selected = {tuple(str(x) for x in k) for k in selections.get(label, [])}
        if not selected:
            continue

        df_class = pd.DataFrame([
            dict(subject=k[0], session=k[1], run=k[2], seed=k[3],
                 is_outlier=tuple(str(x) for x in k) not in selected)
            for k in profiles
        ])
        merged, level = merge_fn(df_qc, df_class, label)
        if merged is None:
            continue

        kept_mask = ~merged["is_outlier"].astype(bool)
        for metric in metrics:
            if metric not in merged.columns:
                continue
            kept = pd.to_numeric(merged.loc[kept_mask, metric], errors="coerce").dropna().values
            rest = pd.to_numeric(merged.loc[~kept_mask, metric], errors="coerce").dropna().values
            if len(kept) < 2 or len(rest) < len(kept):
                continue
            deltas, aucs, pvals = [], [], []
            for _ in range(n_draws):
                draw = rng.choice(rest, size=len(kept), replace=False)
                deltas.append(np.mean(kept) - np.mean(draw))
                try:
                    u, p = _st.mannwhitneyu(kept, draw, alternative="two-sided")
                    aucs.append(u / (len(kept) ** 2)); pvals.append(p)
                except ValueError:
                    continue
            if not pvals:
                continue
            rows.append(dict(label=label, metric=metric, join_level=level,
                            n_kept=len(kept), n_rest=len(rest),
                            delta=float(np.mean(deltas)), auc=float(np.mean(aucs)),
                            p_value=float(np.median(pvals)),
                            p_frac_sig=float(np.mean(np.asarray(pvals) < 0.05))))
            if verbose:
                print(f"  {label} {metric}: delta={rows[-1]['delta']:+.4f}, "
                     f"AUC={rows[-1]['auc']:.3f}, median p={rows[-1]['p_value']:.4f}")
    return pd.DataFrame(rows)


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


def permutation_test_phylogeny(distance_df, divergence_times, n_perm=10000,
                               random_seed=42, verbose=True):
    """
    Is the distance-vs-divergence correlation stronger than chance?

    With only ~6-8 species, a Pearson p-value rests on normality
    assumptions that cannot be checked at that n, and one species moves
    it a lot. Shuffling the divergence times across species and
    recomputing gives an exact, assumption-free null for exactly this
    sample size. p is computed with the (1 + count) / (n_perm + 1)
    correction (Phipson & Smyth 2010): the observed arrangement is
    itself one member of the null, so a p of exactly 0 is never
    reportable.

    Confirmed real bug this fixes: this function used to hard-drop
    divergence_mya == 0 (Human) regardless of what exclude_human was
    set to in test_phylogenetic_signal(), so F6 (parametric) and F7
    (this function) silently ran on different samples (n=7 vs n=6) and
    could report contradictory r/p with no indication why. The
    filter's original justification -- "distance 0 from itself carries
    no information" -- was written when Human's distance actually was
    circularly ~0 by construction; since the leave-one-species-out fix,
    Human's distance is a genuine, non-trivial measurement (often one
    of the LARGEST in the dataset), so dropping it is no longer
    self-evidently correct. This function now uses distance_df exactly
    as given -- pass it test_phylogenetic_signal()'s own tr["data"]
    (not a fresh, unfiltered distance_df) so the two figures can never
    diverge on which species were included again.

    Returns dict(observed_r, p_perm, null_r, n_species) or None if fewer
    than 4 species have both a distance and a divergence time.
    """
    from scipy import stats as _st
    rng = np.random.default_rng(random_seed)

    df = distance_df.copy()
    if "divergence_mya" not in df.columns:
        df["divergence_mya"] = df["species"].map(divergence_times)
    df = df[df["divergence_mya"].notna() & df["distance_to_human"].notna()]
    if len(df) < 4:
        if verbose:
            print(f"  [SKIP] permutation test: only {len(df)} species with both values")
        return None

    x = df["divergence_mya"].values
    y = df["distance_to_human"].values
    observed_r = float(_st.pearsonr(x, y)[0])

    null_r = np.empty(n_perm)
    for i in range(n_perm):
        null_r[i] = _st.pearsonr(rng.permutation(x), y)[0]

    p = (1 + int(np.sum(np.abs(null_r) >= abs(observed_r)))) / (n_perm + 1)
    if verbose:
        print(f"  Phylogeny permutation test (n={len(df)} species, {n_perm} shuffles):")
        print(f"    observed r={observed_r:+.4f}, permutation p={p:.4f} (two-sided)")
        if len(df) < 6:
            print(f"    [NOTE] with n={len(df)}, a single species shifts this substantially.")
    return dict(observed_r=observed_r, p_perm=float(p), null_r=null_r, n_species=len(df))



# === Figures =============================================================

def plot_consensus_pattern(consensus_by_label, common_base_names, title=None,
                           output_path=None, dpi=300):
    """
    Plot each species' consensus vector. With leave_one_out=True these
    differ slightly per species (each excludes itself), so overlaying
    them shows how much any single species moves the consensus -- a
    large spread means the consensus is fragile and driven by a few
    datasets, which matters for interpreting anything built on it.
    """
    if not consensus_by_label:
        print("  [plot] nothing to plot")
        return None

    ref_vec = np.nanmedian(np.vstack(list(consensus_by_label.values())), axis=0)
    order = np.argsort(np.where(np.isfinite(ref_vec), ref_vec, np.inf))
    names_sorted = [common_base_names[i] for i in order]
    x = np.arange(len(names_sorted))

    fig, ax = plt.subplots(figsize=(max(10, len(names_sorted) * 0.26), 5.5))
    for label, vec in sorted(consensus_by_label.items()):
        species = next((s for s in SPECIES_COLORS if label.startswith(s)), None)
        ax.plot(x, np.asarray(vec)[order], lw=1.2, alpha=0.8,
               color=SPECIES_COLORS.get(species), label=label)
    ax.plot(x, ref_vec[order], color="black", lw=2.4, zorder=5,
           label="median across variants")

    ax.axhline(0, color="#999999", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(names_sorted, rotation=90, fontsize=4)
    ax.set_ylabel("consensus (z-scored median_r)", fontsize=9)
    ax.set_title(title or "Cross-species consensus pattern", fontsize=11, fontweight="bold")
    ax.legend(fontsize=6, ncol=2)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_consensus_by_bids(species_profiles, common_base_names, species_of_label,
                           value_col="median_r", zscore=True, title=None,
                           output_path=None, dpi=300):
    """
    Companion to plot_consensus_by_species: same black consensus line
    (species-collapsed, so a 2-BIDS species still counts once -- NOT
    recomputed here, taken as-is from the species-level figure so both
    plots are directly comparable), but one colored line per BIDS
    instead of per species (10 lines here vs 7 there) -- lets you see,
    e.g., whether Macaque_BenHamed and Macaque_ZhuGarin actually agree
    with each other before they get collapsed into one "Macaque" line.
    Colors follow the same per-BIDS PALETTE as F6/F8/F9 (_palette_color_map),
    not SPECIES_COLORS.
    """
    by_species, per_bids = {}, {}
    for label, prof in species_profiles.items():
        vals = _ordered_common_values(prof, common_base_names, value_col)
        if not np.isfinite(vals).any():
            continue
        vec = _normalize_values(vals, zscore)
        per_bids[label] = vec
        by_species.setdefault(species_of_label.get(label, label), []).append(vec)

    if not per_bids:
        print("  [plot] nothing to plot")
        return None

    with np.errstate(invalid="ignore"):
        per_species = {sp: np.nanmedian(np.vstack(v), axis=0) for sp, v in by_species.items()}
        consensus = np.nanmedian(np.vstack(list(per_species.values())), axis=0)

    order = np.argsort(np.where(np.isfinite(consensus), consensus, np.inf))
    names_sorted = [common_base_names[i] for i in order]
    x = np.arange(len(names_sorted))
    bids_colors = _palette_color_map(per_bids.keys())

    fig, ax = plt.subplots(figsize=(max(10, len(names_sorted) * 0.3), 6))
    for label, vec in sorted(per_bids.items()):
        ax.plot(x, vec[order], lw=1.6, alpha=0.9, color=bids_colors[label], label=label)
    ax.plot(x, consensus[order], color="black", lw=3.0, zorder=5,
           label=f"consensus ({len(per_species)} species)")

    ax.axhline(0, color="#999999", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(names_sorted, rotation=90, fontsize=6)
    ax.set_ylabel(f"{_norm_label(zscore)}{value_col}", fontsize=10)
    ax.set_title(title or "Consensus pattern, one line per BIDS", fontsize=12, fontweight="bold")
    ax.legend(fontsize=7, ncol=2)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_consensus_by_species(species_profiles, common_base_names, species_of_label,
                              value_col="median_r", zscore=True, title=None,
                              output_path=None, dpi=300):
    """
    FIGURE: the consensus with ONE LINE PER SPECIES (BIDS of the same
    species collapsed by median first), plus the overall consensus in
    black.

    Why per species and not per BIDS: with 2 Macaque, 2 Marmoset and 2
    Human datasets against 1 Dog and 1 Mouse, a per-BIDS plot (and a
    per-BIDS consensus) silently weights some species two or three times
    more than others. Collapsing within species first makes every
    species count once, which is what a claim about conserved
    organisation requires.
    """
    by_species = {}
    for label, prof in species_profiles.items():
        vals = _ordered_common_values(prof, common_base_names, value_col)
        if not np.isfinite(vals).any():
            continue
        vec = _normalize_values(vals, zscore)
        by_species.setdefault(species_of_label.get(label, label), []).append(vec)

    if not by_species:
        print("  [plot] nothing to plot")
        return None

    with np.errstate(invalid="ignore"):
        per_species = {sp: np.nanmedian(np.vstack(v), axis=0) for sp, v in by_species.items()}
        consensus = np.nanmedian(np.vstack(list(per_species.values())), axis=0)

    order = np.argsort(np.where(np.isfinite(consensus), consensus, np.inf))
    names_sorted = [common_base_names[i] for i in order]
    x = np.arange(len(names_sorted))

    fig, ax = plt.subplots(figsize=(max(10, len(names_sorted) * 0.3), 6))
    for sp, vec in sorted(per_species.items()):
        n_bids = len(by_species[sp])
        ax.plot(x, vec[order], lw=1.8, alpha=0.9, color=SPECIES_COLORS.get(sp),
               label=f"{sp} ({n_bids} BIDS)")
    ax.plot(x, consensus[order], color="black", lw=3.0, zorder=5,
           label=f"consensus ({len(per_species)} species)")

    ax.axhline(0, color="#999999", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(names_sorted, rotation=90, fontsize=6)
    ax.set_ylabel(f"{_norm_label(zscore)}{value_col}", fontsize=10)
    ax.set_title(title or "Consensus pattern, one line per species",
                fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, ncol=2)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_common_regions_with_human_ref(human_ref, common_base_names,
                                       value_col="median_r", title=None,
                                       output_path=None, dpi=300,
                                       sort_by_reference=True):
    """
    DEPRECATED -- not called by the current pipeline (selection is no longer
    anchored to a human reference); kept only for old scripts.

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
    DEPRECATED -- not called by the current pipeline (selection is no longer
    anchored to a human reference); kept only for old scripts.

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


def plot_bids_ref_with_individuals(profiles, species_ref, consensus_ref, common_base_names,
                                   selected_keys=None, kept_ref=None, notkept_ref=None,
                                   grand_mean_ref=None, kept_ref_round2=None,
                                   actual_selection_target=None,
                                   final_round_label="Kept final round",
                                   target_label="Round-1 selection target (right axis)",
                                   value_col="median_r",
                                   title=None, output_path=None, dpi=300,
                                   sort_by_reference=True):
    """
    FIGURE 3 (per BIDS) -- the species' own reference, the cross-species
    consensus, and every individual run, in the same common-region space.

    Confirmed real problem this fixes: the previous version plotted only
    the POST-HOC group map (built from the 3dttest++ on the selected
    runs) as a single black "Reference" line -- which is NOT what the
    selection was actually based on, is DIFFERENT for every species (so
    it looked wrong that the black line "should be the same for
    everyone" -- it should, but wasn't, because it wasn't the shared
    reference at all), and can have GAPS where the thresholded group map
    lost coverage in a region the run-level profiles still had. Using
    PROFILE-based references (built the same way as the common-region
    list itself, so complete by construction) fixes both.

    Also confirmed real problem: the previous version of this function
    called the second reference "human_ref" and labeled it "Human
    reference" in the legend, on EVERY panel including non-human species
    -- leftover naming from an earlier version of the pipeline where
    selection really was anchored to a literal human reference. It is
    now the generic cross-species leave-one-out consensus (identical
    curve on every panel, human or not), so it is renamed consensus_ref
    / "Cross-species consensus" throughout. This is exactly why on the
    two Human panels the "Human reference" line didn't look like an
    average of the two Human datasets -- it wasn't one; it's the
    consensus across all 7 species, so it has no reason to sit between
    Castrillon's and Park's own curves at any given region.

    species_ref: this BIDS's OWN reference (kcore over ALL its own
    runs, kept and not-kept alike) -- "what does a typical run from this
    dataset look like".
    consensus_ref: the shared cross-species consensus every species was
    ranked against for selection (same curve on every panel).
    kept_ref / notkept_ref: optional -- kcore reference built from just
    the kept, or just the not-kept, runs of THIS BIDS (round 1). Lets
    you see whether the selection actually shifted the average, or
    whether kept and not-kept look about the same (weak selection).
    kept_ref_round2: optional -- kcore reference built from round 2's
    kept runs (species-anchored re-selection, see run_round2_selection):
    round 1's kept-only average re-ranked against itself rather than
    the cross-species consensus. Shown alongside kept_ref so you can see
    how much re-anchoring to the species' own signal moved the average.
    grand_mean_ref: optional -- unweighted median of every BIDS's own
    species_ref (each of the 10 BIDS counted once, unlike consensus_ref
    which collapses same-species BIDS together first).
    selected_keys: if given, runs in this set are drawn highlighted
    (green) and the rest faint grey, showing WHICH runs the selection
    actually kept relative to the overall spread.
    actual_selection_target: optional -- the EXACT vector selection
    actually scored this BIDS's runs against (sel["consensus_by_label"][label]
    from run_consensus_selection): leave-one-species-out AND z-scored,
    unlike consensus_ref above (all species included, always raw scale,
    identical curve on every panel -- confirmed real gap: consensus_ref
    was never literally what any panel's runs were judged against).
    Because z-scoring happens BEFORE the per-species median that builds
    the consensus, not after, this can differ from consensus_ref in
    SHAPE, not just scale -- rescaling consensus_ref would not
    reproduce it. Plotted on a separate right-hand axis (different
    units: z-scored) rather than forced onto the left axis's raw scale.
    Must be given in the ORIGINAL common_base_names order (the order
    this function received, before any sort_by_reference re-sort).
    """
    original_names = list(common_base_names)
    if sort_by_reference:
        common_base_names = _sort_regions_by_value(consensus_ref, common_base_names, value_col)
    region_reorder_idx = [original_names.index(r) for r in common_base_names]
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
           color="#1A5276", lw=2.6, zorder=5, label="Species reference (all runs)", ls="-")
    ax.plot(x, _ordered_common_values(consensus_ref, common_base_names, value_col),
           color="#D35400", lw=2.6, zorder=6, label="Cross-species consensus", ls="--")

    handles = [
        plt.Line2D([], [], color="#1A5276", lw=2.6, label="Species reference (all runs)"),
        plt.Line2D([], [], color="#D35400", lw=2.6, ls="--", label="Cross-species consensus"),
        plt.Line2D([], [], color="#009E73", lw=1.2, label=f"Selected runs (n={n_sel})"),
        plt.Line2D([], [], color="#BBBBBB", lw=1.2, label=f"Other runs (n={n_other})"),
    ]

    if kept_ref is not None:
        ax.plot(x, _ordered_common_values(kept_ref, common_base_names, value_col),
               color="#009E73", lw=2.0, zorder=7, ls=":")
        handles.append(plt.Line2D([], [], color="#009E73", lw=2.0, ls=":",
                                  label="Kept-only reference"))
    if notkept_ref is not None:
        ax.plot(x, _ordered_common_values(notkept_ref, common_base_names, value_col),
               color="#888888", lw=2.0, zorder=7, ls=":")
        handles.append(plt.Line2D([], [], color="#888888", lw=2.0, ls=":",
                                  label="Not-kept-only reference"))
    if kept_ref_round2 is not None:
        ax.plot(x, _ordered_common_values(kept_ref_round2, common_base_names, value_col),
               color="#CC79A7", lw=2.2, zorder=8, ls="-")
        handles.append(plt.Line2D([], [], color="#CC79A7", lw=2.2,
                                  label=final_round_label))
    if grand_mean_ref is not None:
        ax.plot(x, _ordered_common_values(grand_mean_ref, common_base_names, value_col),
               color="#000000", lw=1.8, zorder=4, ls="-.")
        handles.append(plt.Line2D([], [], color="#000000", lw=1.8, ls="-.",
                                  label="Grand mean (all 10 BIDS, unweighted)"))

    ax2 = None
    if actual_selection_target is not None:
        ax2 = ax.twinx()
        tvals = np.asarray(actual_selection_target, dtype=float)
        ax2.plot(x, tvals[region_reorder_idx], color="#8E44AD", lw=2.2, zorder=9, ls="-")
        ax2.set_ylabel(target_label, fontsize=9, color="#8E44AD")
        ax2.tick_params(axis="y", colors="#8E44AD")
        ax2.spines["top"].set_visible(False)
        handles.append(plt.Line2D([], [], color="#8E44AD", lw=2.2, label=target_label))

    ax.legend(handles=handles, fontsize=8, loc="upper left")

    ax.axhline(0, color="#999999", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(common_base_names, rotation=90, fontsize=4)
    ax.set_ylabel(value_col, fontsize=10)
    ax.set_title(title or "Species vs consensus reference, with individual runs (common regions)",
                fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_score_distributions(rankings, selections, title=None, output_path=None,
                             dpi=300, ncols=3):
    """
    Per BIDS: histogram of every run's consensus score, with the actual
    selection cutoff drawn. Shows whether the selected runs sit in a
    distinct upper mode (a real subgroup) or are just the right tail of
    one continuous distribution (an arbitrary cut) -- which matters a lot
    for how the selection can be described in a paper.
    """
    labels = sorted(rankings)
    if not labels:
        print("  [plot] nothing to plot")
        return None
    ncols = min(ncols, len(labels))
    nrows = int(np.ceil(len(labels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4.2, nrows * 3.2), squeeze=False)

    for i, label in enumerate(labels):
        ax = axes[i // ncols][i % ncols]
        vals = pd.to_numeric(rankings[label]["pearson_r"], errors="coerce").dropna()
        if vals.empty:
            ax.axis("off"); continue
        n_sel = len(selections.get(label, []))
        ax.hist(vals, bins=30, color="#4472c4", alpha=0.75, edgecolor="white")
        if n_sel and n_sel <= len(vals):
            cutoff = float(np.sort(vals)[::-1][n_sel - 1])
            ax.axvline(cutoff, color="#D55E00", lw=2.2)
            ax.text(cutoff, ax.get_ylim()[1] * 0.85, f" cutoff\n n={n_sel}",
                   color="#D55E00", fontsize=7, va="top")
        ax.set_title(f"{label} (n={len(vals)})", fontsize=9, fontweight="bold")
        ax.set_xlabel("score vs consensus", fontsize=7)
        ax.set_ylabel("runs", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    for j in range(len(labels), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.suptitle(title or "Consensus score distribution per BIDS", fontsize=11, fontweight="bold")
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def _within_across_species_pairs(profiles_by_label, common_base_names, species_of_label,
                                 value_col="median_r", zscore=True, max_pairs=20000,
                                 random_seed=0):
    """
    Pairwise run-to-run profile correlations, with, for every pair, the two
    datasets and species it comes from. Returns a dict of arrays:
    r, da, db (dataset index), same_dataset, same_species, plus
    labels (dataset names) and species_of (species per dataset index).
    Pairs are subsampled above max_pairs.
    """
    rng = np.random.default_rng(random_seed)
    labels = sorted(profiles_by_label)
    vecs, ds = [], []
    for i, label in enumerate(labels):
        for key, prof in profiles_by_label[label].items():
            vals = _ordered_common_values(prof, common_base_names, value_col)
            if not np.isfinite(vals).any():
                continue
            vecs.append(_normalize_values(vals, zscore))
            ds.append(i)
    out = dict(labels=labels, species_of=np.array([species_of_label.get(l, l) for l in labels]),
               r=np.array([]), da=np.array([], int), db=np.array([], int))
    n = len(vecs)
    if n < 4:
        return out
    arr, ds = np.vstack(vecs), np.array(ds)
    ii, jj = np.triu_indices(n, k=1)
    if len(ii) > max_pairs:
        sel = rng.choice(len(ii), size=max_pairs, replace=False)
        ii, jj = ii[sel], jj[sel]
    r, keep = np.empty(len(ii)), np.zeros(len(ii), bool)
    for t, (i, j) in enumerate(zip(ii, jj)):
        both = np.isfinite(arr[i]) & np.isfinite(arr[j])
        if both.sum() >= 3:
            r[t] = np.corrcoef(arr[i][both], arr[j][both])[0, 1]
            keep[t] = True
    out.update(r=r[keep], da=ds[ii[keep]], db=ds[jj[keep]])
    return out


def _within_across_stats(pairs, n_perm=2000, random_seed=0):
    """
    Three pair classes: same dataset / same species but another dataset /
    different species. Descriptive AUC (same species vs different species)
    plus a permutation p that respects the dependence between pairs: the
    species labels are shuffled ACROSS DATASETS (each dataset keeps all its
    runs), the pair correlations stay fixed, and the statistic
    median(same species) - median(different species) is recomputed.
    Pairs are not independent (each run enters hundreds of pairs), so a
    Mann-Whitney p on pairs is anti-conservative by orders of magnitude --
    the reason it was replaced. With few datasets per species the
    permutation null is coarse and the test is weak: read the three
    medians, not only p.
    """
    r, da, db = pairs["r"], pairs["da"], pairs["db"]
    sp = pairs["species_of"]
    if len(r) < 4:
        return None
    same_ds = da == db
    same_sp = sp[da] == sp[db]
    within, across = r[same_sp], r[~same_sp]
    if len(within) < 2 or len(across) < 2:
        return None
    u, _ = _stats.mannwhitneyu(within, across, alternative="two-sided")
    auc = u / (len(within) * len(across))
    obs = np.median(within) - np.median(across)
    rng = np.random.default_rng(random_seed)
    null = []
    for _ in range(n_perm):
        spp = rng.permutation(sp)
        w = spp[da] == spp[db]
        if w.all() or not w.any():
            continue
        null.append(np.median(r[w]) - np.median(r[~w]))
    null = np.asarray(null)
    p_perm = (1 + int(np.sum(null >= obs - 1e-12))) / (len(null) + 1) if len(null) else np.nan
    cross_ds_same_sp = r[same_sp & ~same_ds]
    return dict(n_same_dataset=int(same_ds.sum()),
                n_same_species_other_dataset=int(len(cross_ds_same_sp)),
                n_different_species=int(len(across)),
                median_same_dataset=float(np.median(r[same_ds])) if same_ds.any() else np.nan,
                median_same_species_other_dataset=(float(np.median(cross_ds_same_sp))
                                                   if len(cross_ds_same_sp) else np.nan),
                median_different_species=float(np.median(across)),
                median_within=float(np.median(within)), median_across=float(np.median(across)),
                auc=float(auc), p_perm_dataset_level=float(p_perm), n_perm=int(len(null)))


def _draw_within_across_panel(ax, pairs, title=None):
    """Three histograms (same dataset / same species other dataset /
    different species) + descriptive AUC and dataset-level permutation p."""
    st = _within_across_stats(pairs)
    if st is None:
        ax.axis("off")
        return None
    r, da, db, sp = pairs["r"], pairs["da"], pairs["db"], pairs["species_of"]
    same_ds, same_sp = da == db, sp[da] == sp[db]
    bins = np.linspace(-1, 1, 60)
    groups = [(r[same_ds], "#009E73", "same dataset", st["median_same_dataset"]),
              (r[same_sp & ~same_ds], "#0072B2", "same species, other dataset",
               st["median_same_species_other_dataset"]),
              (r[~same_sp], "#D55E00", "different species", st["median_different_species"])]
    for vals, col, lab, med in groups:
        if len(vals):
            ax.hist(vals, bins=bins, alpha=0.55, color=col, density=True,
                    label=f"{lab} (n={len(vals)}, median={med:.3f})")
            ax.axvline(med, color=col, lw=2, ls="--")
    ax.set_xlabel("run-to-run profile correlation", fontsize=10)
    ax.set_ylabel("density", fontsize=10)
    ax.set_title((title + "\n" if title else "") +
                 f"same vs different species: AUC={st['auc']:.3f}, "
                 f"dataset-level permutation p={st['p_perm_dataset_level']:.3g}",
                 fontsize=10, fontweight="bold")
    ax.legend(fontsize=7)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    return st


def plot_intra_vs_inter_species_distance(profiles_by_label, common_base_names,
                                         species_of_label, value_col="median_r",
                                         zscore=True, title=None, output_path=None,
                                         dpi=300, max_pairs=20000, random_seed=0):
    """
    FIGURE F2: run-to-run profile correlation for pairs from the same
    dataset, from the same species but another dataset, and from different
    species. The middle class is the one that isolates species from
    dataset/protocol and animal identity; if it is not clearly above the
    different-species class, species is not the dominant structure.
    """
    pairs = _within_across_species_pairs(profiles_by_label, common_base_names,
                                         species_of_label, value_col=value_col, zscore=zscore,
                                         max_pairs=max_pairs, random_seed=random_seed)
    fig, ax = plt.subplots(figsize=(9, 5))
    st = _draw_within_across_panel(ax, pairs, title=title or "Within- vs across-species similarity")
    if st is None:
        plt.close(fig)
        print("  [plot] not enough pairs")
        return None, None
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    print(f"    medians: same dataset={st['median_same_dataset']:.4f}, same species other "
          f"dataset={st['median_same_species_other_dataset']:.4f}, different species="
          f"{st['median_different_species']:.4f}; AUC={st['auc']:.3f}, "
          f"dataset-level permutation p={st['p_perm_dataset_level']:.3g}")
    return fig, st


def plot_within_vs_across_all_vs_selected(all_runs_by_label, selected_by_label,
                                          common_base_names, species_of_label,
                                          value_col="median_r", zscore=True, title=None,
                                          output_path=None, dpi=300, max_pairs=20000,
                                          random_seed=0, selected_name="Selected runs"):
    """F2 computed on every run (left) and on the selected runs (right)."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 5.5))
    stats = {}
    for ax, (name, pbl) in zip(axes, [("All runs", all_runs_by_label),
                                      (selected_name, selected_by_label)]):
        pairs = _within_across_species_pairs(pbl, common_base_names, species_of_label,
                                             value_col=value_col, zscore=zscore,
                                             max_pairs=max_pairs, random_seed=random_seed)
        n_runs = sum(len(p) for p in pbl.values())
        stats[name] = _draw_within_across_panel(ax, pairs, title=f"{name} (n={n_runs} runs)")
    fig.suptitle(title or "Within- vs across-species similarity: all runs vs selected",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    for name, st in stats.items():
        if st:
            print(f"    [{name}] same species other dataset={st['median_same_species_other_dataset']:.4f}"
                  f" vs different species={st['median_different_species']:.4f}, "
                  f"AUC={st['auc']:.3f}, permutation p={st['p_perm_dataset_level']:.3g}")
    return fig, stats

def _draw_distance_divergence_panel(ax, distance_df, test_result=None, divergence_times=None,
                                    exclude_human=False, title=None, bids_colors=None):
    """
    Shared drawing routine behind plot_distance_vs_divergence and
    plot_distance_vs_divergence_multilevel, factored out so both draw
    identically instead of risking the two styles drifting apart.

    exclude_human MUST match whatever exclude_human was passed to
    test_phylogenetic_signal() to build test_result -- otherwise the
    drawn regression line and the printed r/p come from two different
    fits. (Confirmed real bug in the previous version of this function:
    fit_df = df[species != "Human"] was hardcoded regardless of what
    test_result was actually computed from, so once the caller started
    passing exclude_human=False the annotation still said "human
    excluded from fit" while the line kept excluding it anyway.)
    """
    divergence_times = divergence_times or DIVERGENCE_TIME_MYA
    df = distance_df.copy()
    df["divergence_mya"] = df["species"].map(divergence_times)
    df = df[df["divergence_mya"].notna() & df["distance_to_human"].notna()]
    if df.empty:
        ax.axis("off")
        return

    for _, row in df.iterrows():
        is_human = row["species"] == "Human"
        point_color = (bids_colors.get(row["label"], "#888888")
                       if bids_colors is not None and "label" in df.columns
                       else SPECIES_COLORS.get(row["species"], "#888888"))
        ax.scatter(row["divergence_mya"], row["distance_to_human"],
                  s=140, color=point_color,
                  edgecolors="black", lw=1.0, zorder=3,
                  facecolors="none" if (is_human and exclude_human) else None,
                  label=row["species"])
        # Multiple BIDS of the same species (Human, Macaque, Marmoset) would
        # otherwise print the exact same label text on top of each other --
        # append the BIDS name whenever this df has more than one row for
        # that species (a species-level df, one row per species, never
        # triggers this: label is then identical to species anyway).
        same_species_n = (df["species"] == row["species"]).sum()
        label_txt = (f"{row['species']} ({row['label']})"
                    if "label" in df.columns and same_species_n > 1
                    else row["species"])
        ax.annotate(label_txt, (row["divergence_mya"], row["distance_to_human"]),
                   textcoords="offset points", xytext=(7, 5), fontsize=8)

    fit_df = df[df["species"] != "Human"] if exclude_human else df
    if len(fit_df) >= 3:
        slope, intercept, r, p, se = _stats.linregress(
            fit_df["divergence_mya"], fit_df["distance_to_human"])
        xs = np.linspace(fit_df["divergence_mya"].min(), fit_df["divergence_mya"].max(), 100)
        ax.plot(xs, intercept + slope * xs, "k--", lw=1.4, alpha=0.7, zorder=2)

    if test_result:
        human_note = "human excluded from fit" if exclude_human else "human included in fit"
        ax.text(0.03, 0.97,
               f"Pearson r={test_result['pearson_r']:.3f}, p={test_result['pearson_p']:.3f}\n"
               f"Spearman rho={test_result['spearman_rho']:.3f}, p={test_result['spearman_p']:.3f}\n"
               f"n={test_result['n']} species ({human_note})",
               transform=ax.transAxes, va="top", fontsize=8,
               bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))

    ax.set_xlabel("Divergence time from human (Mya)", fontsize=9)
    ax.set_ylabel("Distance to leave-one-species-out consensus", fontsize=9)
    if title:
        ax.set_title(title, fontsize=10, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_distance_vs_divergence(distance_df, test_result=None, divergence_times=None,
                                exclude_human=False, title=None, output_path=None, dpi=300):
    """
    Scatter of distance-to-consensus against divergence time, one point
    per species, colored by species, with a fitted line and the
    correlation stats annotated.

    exclude_human: pass the SAME value used for test_phylogenetic_signal
    when building test_result (see _draw_distance_divergence_panel).
    Default False matches compute_loo_species_distances' leave-one-
    SPECIES-out reference, which makes Human's own distance non-
    circular -- there is no longer a reason to exclude it by default.
    """
    dt = divergence_times or DIVERGENCE_TIME_MYA
    check = distance_df.copy()
    check["divergence_mya"] = check["species"].map(dt)
    if check["divergence_mya"].notna().sum() == 0 or distance_df.empty:
        print("  [plot] nothing to plot")
        return None

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    bids_colors = (_palette_color_map(distance_df["label"])
                  if "label" in distance_df.columns else None)
    _draw_distance_divergence_panel(
        ax, distance_df, test_result, divergence_times, exclude_human=exclude_human,
        title=title or "Distance to leave-one-species-out consensus vs divergence",
        bids_colors=bids_colors)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_distance_vs_divergence_multilevel(results_by_level, divergence_times=None,
                                           output_path=None, dpi=300):
    """
    One panel per atlas level (results_by_level: {level: {"ddf_species":
    df, "test_result": dict, "n_common_regions": int}}), same drawing
    routine as plot_distance_vs_divergence, to compare how segmentation
    granularity changes the phylogenetic slope/fit.
    """
    levels = sorted(results_by_level)
    n = len(levels)
    if n == 0:
        print("  [plot] nothing to plot")
        return None
    n_cols = min(2, n)
    n_rows = math.ceil(n / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7 * n_cols, 5 * n_rows), squeeze=False)
    axes_flat = axes.ravel()
    for ax, level in zip(axes_flat, levels):
        res = results_by_level[level]
        _draw_distance_divergence_panel(
            ax, res["ddf_species"], res.get("test_result"), divergence_times,
            exclude_human=False,
            title=(f"Atlas level {level} ({res.get('n_common_regions', '?')} common regions)"
                   + (f"\nMantel r={res['mantel']['mantel_pearson']['r']:+.2f}, "
                      f"exact p(>)={res['mantel']['mantel_pearson']['p_greater']:.3f}"
                      if res.get("mantel") else "")))
    for ax in axes_flat[n:]:
        ax.set_visible(False)
    fig.suptitle("Impact of segmentation level on distance-to-consensus vs divergence",
                fontsize=12, fontweight="bold")
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_run_distance_distribution(dist_df, value_col="pearson_r", output_path=None,
                                   title=None, dpi=300):
    """
    One violin+jittered-strip per BIDS, showing every run's correlation
    (or distance) to its BIDS' reference map. Self-contained rather than
    calling EDNiX_figures' private _violin_strip_quartiles (a different
    package), but same visual spirit: grey violin backdrop, colored
    jittered points, quartile ticks.
    """
    if dist_df is None or dist_df.empty:
        print("  [plot] nothing to plot")
        return None
    labels = sorted(dist_df["label"].unique())
    bids_colors = _palette_color_map(labels)
    fig, ax = plt.subplots(figsize=(max(6, 1.1 * len(labels)), 5))
    rng = np.random.default_rng(42)
    for xi, lab in enumerate(labels):
        sub = dist_df.loc[dist_df["label"] == lab, [value_col] + (
            ["is_kept"] if "is_kept" in dist_df.columns else [])].dropna(subset=[value_col])
        vals = sub[value_col].values
        if len(vals) == 0:
            continue
        color = bids_colors[lab]
        if len(vals) >= 2:
            parts = ax.violinplot(vals, positions=[xi], widths=0.7, showextrema=False)
            for pc in parts["bodies"]:
                pc.set_facecolor("#dddddd")
                pc.set_alpha(0.6)
        jitter = rng.uniform(-0.15, 0.15, len(vals))
        ax.scatter(xi + jitter, vals, s=22, color=color, edgecolors="black",
                  lw=0.3, alpha=0.85, zorder=3)
        if "is_kept" in sub.columns:
            kept_mask = sub["is_kept"].fillna(False).values.astype(bool)
            if kept_mask.any():
                ax.scatter((xi + jitter)[kept_mask], vals[kept_mask], marker="x", s=70,
                          color="black", linewidths=1.6, zorder=5)
        for q in (25, 50, 75):
            ax.hlines(np.percentile(vals, q), xi - 0.25, xi + 0.25,
                     color="black", lw=1.4 if q == 50 else 0.8, zorder=4)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("Pearson r to BIDS reference" if value_col == "pearson_r"
                 else "Distance (1-r) to BIDS reference", fontsize=10)
    ax.set_title(title or "Per-run distance to the kept-group reference map",
                fontsize=11, fontweight="bold")
    if "is_kept" in dist_df.columns and dist_df["is_kept"].any():
        ax.legend(handles=[plt.Line2D([], [], marker="x", color="black", lw=0,
                                      markersize=9, markeredgewidth=1.6,
                                      label="Selected (kept) run")],
                  fontsize=8, loc="best")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_run_distance_std_by_bids(dist_df, value_col="pearson_r", output_path=None, dpi=300):
    """
    Bar chart: SD of each BIDS' per-run correlation to its reference --
    a compact single-number view of within-BIDS run-to-run consistency,
    companion to plot_run_distance_distribution's full spread.
    """
    if dist_df is None or dist_df.empty:
        print("  [plot] nothing to plot")
        return None
    summary = (dist_df.groupby(["label", "species"], as_index=False)[value_col]
              .agg(std="std", n="count")).sort_values("std", ascending=False)
    bids_colors = _palette_color_map(summary["label"])
    fig, ax = plt.subplots(figsize=(max(6, 1.1 * len(summary)), 4.5))
    colors = [bids_colors[lab] for lab in summary["label"]]
    ax.bar(summary["label"], summary["std"], color=colors, edgecolor="black", lw=0.5)
    for i, (n, s) in enumerate(zip(summary["n"], summary["std"])):
        ax.text(i, s, f"n={n}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(len(summary)))
    ax.set_xticklabels(summary["label"], rotation=35, ha="right", fontsize=8)
    ax.set_ylabel(f"SD of per-run {value_col} to reference", fontsize=10)
    ax.set_title("Run-to-run consistency by BIDS", fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_qc_distance_correlation(meta_df, output_path=None, dpi=300, alpha=0.05, max_metrics=40):
    """
    Bar chart of the meta-analytic r (meta_analyze_qc_correlations) between
    each QC metric and distance to reference, for EVERY metric tested (top
    max_metrics by |r|). Metrics significant after Benjamini-Hochberg
    correction (meta_q_fdr < alpha) are coloured (blue r>0, orange r<0);
    the others are grey. Previously only significant metrics were drawn,
    so a run where nothing passed produced no figure at all -- which reads
    as a crash and hides a genuine null result.
    """
    if meta_df is None or meta_df.empty:
        print("  [plot] nothing to plot")
        return None
    pcol = "meta_q_fdr" if "meta_q_fdr" in meta_df.columns else "meta_p"
    sub = meta_df.reindex(meta_df["meta_r"].abs().sort_values(ascending=False).index)
    sub = sub.head(max_metrics).reset_index(drop=True)
    sig = sub[pcol] < alpha
    colors = [("#0072B2" if r > 0 else "#D55E00") if g else "#BBBBBB"
              for r, g in zip(sub["meta_r"], sig)]
    fig, ax = plt.subplots(figsize=(max(6, 0.45 * len(sub)), 5))
    ax.bar(range(len(sub)), sub["meta_r"], color=colors, edgecolor="white")
    for i, (r, n_bids, q) in enumerate(zip(sub["meta_r"], sub["n_bids"], sub[pcol])):
        ax.text(i, r + (0.02 if r >= 0 else -0.02), f"n={n_bids}\nq={q:.2g}", ha="center",
                va="bottom" if r >= 0 else "top", fontsize=6)
    ax.axhline(0, color="#888888", lw=0.8)
    ax.set_xticks(range(len(sub)))
    ax.set_xticklabels(sub["metric"], rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("Meta-analytic r (Fisher-z weighted across BIDS)", fontsize=10)
    lab = "BH-FDR q" if pcol == "meta_q_fdr" else "meta p"
    ax.set_title(f"QC metrics vs distance to reference -- {int(sig.sum())}/{len(meta_df)} "
                 f"significant across BIDS (coloured: {lab}<{alpha})",
                 fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_phylogeny_permutation(perm_result, title=None, output_path=None, dpi=300):
    """Null distribution of the distance-vs-divergence correlation, with
    the observed value marked."""
    if perm_result is None:
        return None
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(perm_result["null_r"], bins=50, color="#888888", alpha=0.75, edgecolor="white")
    ax.axvline(perm_result["observed_r"], color="#D55E00", lw=2.5,
              label=f"observed r={perm_result['observed_r']:+.3f} "
                    f"(p={perm_result['p_perm']:.4f})")
    ax.set_xlabel("Pearson r (distance vs divergence time)", fontsize=10)
    ax.set_ylabel("permutations", fontsize=10)
    ax.set_title(title or f"Phylogeny permutation test "
                          f"(n={perm_result['n_species']} species)",
                fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig