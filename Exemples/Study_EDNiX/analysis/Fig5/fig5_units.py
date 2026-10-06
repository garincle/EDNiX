"""
fig5_units.py
=============
THE UNIT OF ANALYSIS layer for the Figure 5 pipeline.

WHY THIS EXISTS
---------------
Everything downstream of `compute_all_run_profiles` used to operate on
RUNS. That silently weights datasets by how many runs they happen to
contain, at three different places:

1. Within a dataset, the "species profile" was a median over every RUN
   (`compute_group_median_profile(profiles, all keys)`). A subject with
   20 runs then contributes 20 vectors to that median while a subject
   with 1 run contributes 1 -- so the dataset's profile is really
   "whatever the most-scanned animals look like".

2. The selection kept the N best RUNS per dataset. Nothing stopped all
   N from coming from the same animal, in which case the "group map"
   for that species is one animal measured N times.

3. `3dttest++` received one input per RUN, so its degrees of freedom
   counted runs, not animals, and runs of one subject (not independent)
   were treated as independent samples.

Note what was ALREADY correct and is left untouched:
`build_consensus_pattern` collapses datasets of the same species by
median before taking the median across species, so a species with 2 or 3
BIDS already counted once. This module fixes the level BELOW that.

WHAT THIS MODULE DOES
---------------------
`aggregate_profiles_to_unit()` collapses the raw per-run profile dict
into one profile per UNIT (subject by default, optionally session), by
median across that unit's runs -- reusing
`sri.compute_group_median_profile`, the same function the rest of the
pipeline already uses to combine profiles, so no new averaging
convention is introduced.

Every analysis function in fig5_analysis.py works on a generic
{(a, b, c, seed): profile_df} dict and never interprets those first
three fields, so aggregating here changes the unit of the WHOLE
downstream analysis (common regions, species profiles, invariance,
scoring, selection, balanced stats) without touching a single line of
fig5_analysis.py.

The aggregated keys keep the 4-tuple shape expected everywhere:
    unit="subject"  ->  (subject, "ALL", "ALL", seed)
    unit="session"  ->  (subject, session, "ALL", seed)
    unit="run"      ->  unchanged (the old behaviour)

`expand_units_to_runs()` goes back the other way: AFNI still needs the
real per-run BOLD files, so a selection made at subject level is
expanded to every run of the selected subjects before `df_class` is
built. The averaging over those runs then happens inside AFNI (3dMean),
which is what "moyenne calculee avec AFNI" means here.
"""

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

VALID_UNITS = ("run", "session", "subject")


# ═══════════════════════════════════════════════════════════════════════════
# Keys
# ═══════════════════════════════════════════════════════════════════════════

def unit_key(run_key, unit="subject"):
    """Map a raw run key (subject, session, run, seed) to its unit key.

    The placeholder is the literal string "ALL" rather than None or "" so
    that every consumer that does str() on these fields (the selection
    spreadsheet, the QC merge, the figure labels) shows something
    readable instead of "None", and so a unit key can never be confused
    with a real run whose id happens to be missing -- that case already
    uses the sentinel "NA" (compute_all_run_profiles' fallback).
    """
    subject, session, run, seed = (str(x) for x in run_key)
    if unit == "run":
        return (subject, session, run, seed)
    if unit == "session":
        return (subject, session, "ALL", seed)
    if unit == "subject":
        return (subject, "ALL", "ALL", seed)
    raise ValueError(f"Unknown unit {unit!r} -- use one of {VALID_UNITS}")


def aggregate_profiles_to_unit(sri, profiles, unit="subject", verbose=True,
                               label=""):
    """Collapse {run_key: profile_df} into {unit_key: profile_df}.

    Median across the unit's runs, via sri.compute_group_median_profile
    (the pipeline's existing combiner). Returns (unit_profiles, members)
    where members maps each unit key to the list of run keys behind it --
    needed later to expand a selection back to real BOLD files.
    """
    if unit not in VALID_UNITS:
        raise ValueError(f"Unknown unit {unit!r} -- use one of {VALID_UNITS}")

    members = {}
    for run_key in profiles:
        members.setdefault(unit_key(run_key, unit), []).append(run_key)

    if unit == "run":
        return dict(profiles), members

    unit_profiles = {}
    for ukey, run_keys in members.items():
        prof = sri.compute_group_median_profile(profiles, set(run_keys))
        if prof is None or prof.empty:
            if verbose:
                print(f"    [WARN] {label} {ukey[:2]}: median over "
                      f"{len(run_keys)} run(s) produced nothing -- skipped")
            continue
        unit_profiles[ukey] = prof

    if verbose:
        counts = np.array([len(v) for v in members.values()])
        print(f"  {label}: {len(profiles)} run(s) -> {len(unit_profiles)} "
              f"{unit}(s)  (runs per {unit}: min={counts.min()}, "
              f"median={int(np.median(counts))}, max={counts.max()})")
    return unit_profiles, members


def expand_units_to_runs(unit_keys, members):
    """Every run key belonging to the given unit keys.

    Used to build the df_class that AFNI consumes: the selection is made
    at subject level, but 3dttest still needs each subject's actual runs
    so it can average them itself.
    """
    out = []
    for ukey in unit_keys:
        out.extend(members.get(tuple(str(x) for x in ukey), []))
    return out


def check_unit_counts(profiles_by_label, n_select, unit="subject"):
    """Warn when a dataset has fewer units than the selection asks for.

    With unit="run" a dataset with 171 runs could always supply N_SELECT;
    with unit="subject" that same dataset may only have 8 animals, and a
    dataset with 4 animals simply cannot give 5. Silently returning 4 is
    exactly the kind of unequal weighting this change exists to remove,
    so it is reported explicitly and the effective N is returned.
    """
    rows = []
    for label, profs in sorted(profiles_by_label.items()):
        n = len(profs)
        rows.append(dict(label=label, n_units=n, unit=unit,
                         n_requested=n_select, n_available=min(n, n_select),
                         short=n < n_select))
    df = pd.DataFrame(rows)
    if df.empty:
        return df, n_select
    short = df[df["short"]]
    print(f"\n  Units available per dataset (unit={unit}, N_SELECT={n_select}):")
    print(df[["label", "n_units", "n_available"]].to_string(index=False))
    if not short.empty:
        print(f"\n  [WARN] {len(short)} dataset(s) have fewer than {n_select} "
              f"{unit}(s): {', '.join(short['label'])}")
        print(f"  [WARN] those datasets will contribute fewer {unit}s than the "
              f"others, which reintroduces exactly the unequal weighting this "
              f"unit-level analysis removes. Either lower N_SELECT to "
              f"{int(df['n_units'].min())} (the smallest dataset's count, so "
              f"every dataset contributes the same number) or drop those "
              f"datasets.")
    return df, n_select


# ═══════════════════════════════════════════════════════════════════════════
# STEP 1 figure -- species profiles read back off the AFNI group maps
# ═══════════════════════════════════════════════════════════════════════════

def find_matching_atlas_volume(sri, cb, run_keys, map_path, verbose=True,
                               label=""):
    """
    Pick an atlas volume whose spatial grid actually matches the group
    map's grid, instead of blindly using whichever candidate resolves
    first.

    Confirmed real bug this fixes: the previous version took the FIRST
    run key it could resolve an atlas path for (or the dataset's generic
    "representative" volume), with no check that this particular file
    was on the SAME grid as the group map. compute_region_profile
    itself already refuses a mismatched pair ("Grid mismatch: map ...
    vs atlas volume spatial shape ..."), but by then the damage is
    done -- the whole species silently drops out of the figure. Case
    that surfaced this: Rat, UNIT="run" -- the first candidate happened
    to be sub-300101/ses-1's NATIVE-space atlas (48, 82, 48) while the
    group map (built by 3dMean/3dttest++ on the template grid) was
    (48, 48, 30).

    This tries every run's atlas path (deduplicated), then the dataset's
    representative volume, and returns the first one whose shape (first
    3 spatial dims -- a 4th dim, if present, is the multi-level axis,
    not a spatial one) matches the map's. Only a cheap header read
    (nib.load(...).shape) per candidate, not a full profile computation.

    Returns the matching atlas path, or None if no candidate matches (in
    which case the caller's compute_region_profile call is left to
    raise its own explicit "Grid mismatch" error on whatever fallback
    it is given, rather than silently skipping the species here).
    """
    try:
        map_shape = tuple(sri.nib.load(map_path).shape[:3])
    except Exception as e:
        if verbose:
            print(f"    [WARN] {label}: could not read the group map's own "
                  f"shape from {map_path}: {e}")
        return None

    candidates = []
    for run_key in run_keys:
        atlas = sri.find_atlas_volume_path(cb["bids_dir"], run_key[0], run_key[1])
        if atlas and atlas not in candidates:
            candidates.append(atlas)
    rep = sri.pick_representative_atlas_volume(cb["bids_dir"])
    if rep and rep not in candidates:
        candidates.append(rep)

    for atlas in candidates:
        try:
            atlas_shape = tuple(sri.nib.load(atlas).shape[:3])
        except Exception:
            continue
        if atlas_shape == map_shape:
            if verbose and atlas != candidates[0]:
                print(f"    [{label}] grid match found on candidate "
                      f"{candidates.index(atlas) + 1}/{len(candidates)} "
                      f"({os.path.basename(atlas)})")
            return atlas

    if verbose:
        print(f"    [WARN] {label}: none of {len(candidates)} atlas candidate(s) "
              f"match the group map's grid {map_shape} -- tried: "
              f"{[os.path.basename(a) for a in candidates]}")
    return None


def profiles_from_group_maps(sri, map_index, bids_configs, common_regions,
                             atlas_level, members_by_label=None, verbose=True):
    """Read each step-1 AFNI group map back into a regional profile.

    map_index: DataFrame with label, species, map_path, found -- as
    written by step 1. The map itself is the AFNI subject-level mean
    (3dMean within subject, then 3dttest++ across subjects), so the
    profile computed from it is "the species' mean pattern according to
    AFNI", not a re-derivation from the per-run profiles.

    Returns {label: (species, profile_df)}.
    """
    cfg_by_label = {c["label"]: c for c in bids_configs}
    out = {}
    for r in map_index.itertuples():
        if not getattr(r, "found", False) or not r.map_path:
            continue
        cb = cfg_by_label.get(r.label)
        if cb is None:
            continue
        try:
            # members_by_label[label] is {unit_key: [run_key, ...]} --
            # flatten to the real (subject, session, run, seed) run keys
            # behind this label, whatever UNIT was used to build it, so
            # find_matching_atlas_volume always has real run/session ids
            # to resolve paths from (not a unit_key, which under
            # UNIT="subject"/"session" carries "ALL" placeholders that
            # never resolve to a real atlas file).
            run_keys = []
            if members_by_label and r.label in members_by_label:
                for rks in members_by_label[r.label].values():
                    run_keys.extend(rks)
            atlas = find_matching_atlas_volume(
                sri, cb, run_keys, r.map_path, verbose=verbose, label=r.label)
            atlas = atlas or sri.pick_representative_atlas_volume(cb["bids_dir"])
            label_df = sri.parse_itk_snap_label_file(cb["label_path"])
            prof = sri.compute_region_profile(r.map_path, atlas, label_df,
                                              atlas_level=atlas_level)
            out[r.label] = (r.species, prof)
            if verbose:
                print(f"    {r.label}: profile from {os.path.basename(r.map_path)}")
        except Exception as e:
            print(f"    [WARN] {r.label}: could not profile the group map: {e}")
    return out

def reject_bad_runs(sri, fa, profs, method, k, value_col, metric,
                     mad_threshold, label, verbose=True):
    """
    QC au niveau RUN, appliqué AVANT aggregate_profiles_to_unit.

    Chaque run est évalué indépendamment contre une référence robuste du
    dataset (kcore/médoïde sur TOUS les runs bruts, method=method), avec
    un cutoff MAD (classify_runs_by_profile_distance). Un run rejeté ici
    ne peut plus être caché dans la médiane d'un sujet "kept" en aval.

    Retourne (profs_propres, df_qc_runs) -- df_qc_runs a une ligne par
    run avec sa distance/is_outlier, pour export/inspection.
    """
    n_runs_in = len(profs)
    n_subj_in = len({k_[0] for k_ in profs})
    if verbose:
        print(f"  [{label}] run-level QC: {n_runs_in} run(s) across "
              f"{n_subj_in} subject(s), method={method}, k={k}, "
              f"mad_threshold={mad_threshold}")

    if n_runs_in < 4:
        if verbose:
            print(f"  [{label}] [SKIP-QC] only {n_runs_in} run(s), too few "
                  f"for a robust MAD cutoff -- all runs kept unchanged")
        return dict(profs), pd.DataFrame()

    ref_runs, ref_info = fa.build_reference_for_profiles(
        sri, profs, method=method, k=k, value_col=value_col, metric=metric)
    if verbose:
        print(f"  [{label}] reference built: {ref_info}")

    df_qc = sri.classify_runs_by_profile_distance(
        profs, ref_runs, value_col=value_col, metric=metric,
        outlier_threshold=mad_threshold)
    if df_qc.empty:
        if verbose:
            print(f"  [{label}] [WARN] classify_runs_by_profile_distance "
                  f"returned empty -- all runs kept unchanged")
        return dict(profs), df_qc

    bad_keys = {(str(r.subject), str(r.session), str(r.run), str(r.seed))
                for r in df_qc.itertuples() if r.is_outlier}
    clean = {k_: v for k_, v in profs.items()
             if tuple(str(x) for x in k_) not in bad_keys}

    n_subj_out = len({k_[0] for k_ in clean})
    if verbose:
        print(f"  [{label}] run-level QC result: {len(bad_keys)}/{n_runs_in} "
              f"run(s) rejected, {len(clean)} kept "
              f"({n_subj_in} -> {n_subj_out} subject(s) with >=1 valid run)")
        if bad_keys:
            per_subj_rejected = {}
            for k_ in bad_keys:
                per_subj_rejected.setdefault(k_[0], 0)
                per_subj_rejected[k_[0]] += 1
            print(f"  [{label}] rejected runs per subject: {per_subj_rejected}")
        if n_subj_out < n_subj_in:
            print(f"  [{label}] [WARN] {n_subj_in - n_subj_out} subject(s) "
                  f"lost ALL their runs to run-level QC -- they will be "
                  f"absent from the unit-level profiles entirely")

    return clean, df_qc

def plot_species_profiles_from_afni(afni_profiles, common_base_names,
                                    region_mask=None, invariance=None,
                                    value_col="median_r", zscore=True,
                                    title=None, output_path=None, dpi=300,
                                    analysis_fa=None, group_by="species"):
    """FIGURE: one line per species, straight off the step-1 AFNI group
    maps, with the regions that define the cross-species profile shaded.

    What makes this different from F1 (plot_consensus_by_species): F1 is
    built from the per-unit profiles computed in Python, this one is
    built from the volumetric group maps AFNI actually produced on ALL
    units, with no selection anywhere in the chain. If the two disagree
    in shape, the discrepancy is between the profile arithmetic and the
    AFNI group statistics, and that is worth knowing before any of it
    goes in a figure.

    region_mask: boolean over common_base_names -- the regions the run
    selection was restricted to (the invariant core). These are shaded
    green and their names are printed in bold, so the figure answers
    "which regions are actually being used to define the cross-species
    profile" directly.
    """
    if analysis_fa is None:
        import fig5_analysis as analysis_fa
    fa = analysis_fa

    # group_by="species" (default): the BIDS of one species are collapsed
    # to one line -- the published F0.
    # group_by="label": one line PER BIDS, coloured by species, so the two
    # datasets of a species (different lab, scanner and protocol) can be
    # read against each other. That is the reproducibility check: if
    # Human_Castrillon and Human_Park track each other more closely than
    # either tracks another species, the cross-species differences are not
    # protocol artefacts.
    if group_by not in ("species", "label"):
        raise ValueError(f"group_by must be 'species' or 'label', got {group_by!r}")

    by_species = {}
    species_of_group = {}
    for label, (species, prof) in afni_profiles.items():
        vals = fa._ordered_common_values(prof, common_base_names, value_col)
        if not np.isfinite(vals).any():
            continue
        vec = fa._normalize_values(vals, zscore)
        gkey = species if group_by == "species" else label
        by_species.setdefault(gkey, []).append(vec)
        species_of_group[gkey] = species

    if not by_species:
        print("  [plot] no AFNI group-map profile could be read -- nothing to plot")
        return None

    with np.errstate(invalid="ignore"):
        per_species = {sp: np.nanmedian(np.vstack(v), axis=0)
                       for sp, v in by_species.items()}
        # The consensus is ALWAYS the median across SPECIES, never across
        # groups -- otherwise, with group_by="label", a species carried by
        # two BIDS would weigh twice as much in the black reference line.
        by_true_species = {}
        for gkey, vec in per_species.items():
            by_true_species.setdefault(species_of_group[gkey], []).append(vec)
        consensus = np.nanmedian(
            np.vstack([np.nanmedian(np.vstack(v), axis=0)
                       for v in by_true_species.values()]), axis=0)
        n_species_in_consensus = len(by_true_species)

    order = np.argsort(np.where(np.isfinite(consensus), consensus, np.inf))
    names_sorted = [common_base_names[i] for i in order]
    x = np.arange(len(names_sorted))

    mask_sorted = None
    if region_mask is not None:
        mask_sorted = np.asarray(region_mask, dtype=bool)[order]

    n_panels = 2 if invariance is not None else 1
    fig, axes = plt.subplots(
        n_panels, 1, sharex=True,
        figsize=(max(10, len(names_sorted) * 0.32), 6.5 if n_panels == 1 else 8.5),
        gridspec_kw=dict(height_ratios=[3, 1]) if n_panels == 2 else None)
    axes = np.atleast_1d(axes)
    ax = axes[0]

    # Shade the selected (invariant) regions BEFORE the lines, so the
    # shading reads as background rather than as another data series.
    if mask_sorted is not None:
        for i, keep in enumerate(mask_sorted):
            if keep:
                ax.axvspan(i - 0.5, i + 0.5, color="#2E8B57", alpha=0.10, lw=0)

    # With group_by="label", sibling datasets of the same species share the
    # species colour and are told apart by line style, so proximity within
    # a colour reads as reproducibility at a glance.
    styles = ["-", "--", ":", "-."]
    seen_per_species = {}
    for gkey, vec in sorted(per_species.items()):
        sp = species_of_group[gkey]
        if group_by == "species":
            label_txt = f"{sp} ({len(by_species[gkey])} BIDS)"
            style = "-"
        else:
            i = seen_per_species.get(sp, 0)
            seen_per_species[sp] = i + 1
            style = styles[i % len(styles)]
            label_txt = gkey
        ax.plot(x, vec[order], lw=1.8, alpha=0.9, linestyle=style,
                color=fa.SPECIES_COLORS.get(sp), label=label_txt)
    ax.plot(x, consensus[order], color="black", lw=3.0, zorder=5,
            label=f"consensus ({n_species_in_consensus} species)")

    ax.axhline(0, color="#999999", lw=0.8)
    ax.set_ylabel(f"{fa._norm_label(zscore)}{value_col}\n(AFNI group map)",
                  fontsize=10)
    ax.set_title(title or "Species profiles from the all-units AFNI group maps",
                 fontsize=12, fontweight="bold")
    handles, labels = ax.get_legend_handles_labels()
    if mask_sorted is not None:
        handles.append(plt.Rectangle((0, 0), 1, 1, color="#2E8B57", alpha=0.25))
        labels.append(f"used to define the cross-species profile "
                      f"({int(mask_sorted.sum())}/{len(mask_sorted)} regions)")
    ax.legend(handles, labels, fontsize=8, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if invariance is not None:
        ax2 = axes[1]
        inv_sorted = np.asarray(invariance, dtype=float)[order]
        colors = ["#2E8B57" if (mask_sorted is None or mask_sorted[i]) else "#BBBBBB"
                  for i in range(len(inv_sorted))]
        ax2.bar(x, inv_sorted, color=colors, edgecolor="white", linewidth=0.3)
        ax2.set_ylabel("invariance\nacross species", fontsize=9)
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)

    ax_last = axes[-1]
    ax_last.set_xticks(x)
    ax_last.set_xticklabels(names_sorted, rotation=90, fontsize=6)
    if mask_sorted is not None:
        for tick, keep in zip(ax_last.get_xticklabels(), mask_sorted):
            if keep:
                tick.set_fontweight("bold")
                tick.set_color("#1E5631")

    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig
