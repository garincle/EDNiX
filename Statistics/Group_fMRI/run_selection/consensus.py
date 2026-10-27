"""
Consensus run selection -- profile helpers, cross-species consensus, iterative
selection, region invariance, kept-vs-not-kept balance, QC validation.

Moved VERBATIM from Exemples/Study_EDNiX/analysis/Fig5/fig5_analysis.py (2026-10-27):
same code, same defaults, same results. fig5_analysis imports these names back, so
every caller keeps working unchanged.
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


# Same 8-color, colorblind-safe (Okabe-Ito) palette as Supp_Figure_4-6_QC.py's
# PALETTE, so every figure in this module uses the same color code the QC
# figures already use -- kept as a literal copy (not imported) since
# fig5_analysis.py has no dependency on ednix_bids_tools.py.
PALETTE = [
    "#0072B2", "#E69F00", "#009E73", "#CC79A7",
    "#56B4E9", "#D55E00", "#F0E442", "#000000",
]


def _palette_color_map(names):
    """{name: colour} in the colour code of the article (analysis/paper_style.py):
    a species name gets its species colour (= Figure_2_anat.SPECIES_COLORS), a dataset
    label ('Macaque_BenHamed', 'Rat', ...) the colour of its species, so that a
    species has one colour in every figure. The two datasets of a species share the
    colour; plots that draw them side by side tell them apart by line style or
    marker (dataset_style). Unknown names fall back on PALETTE."""
    names = sorted(set(names))
    out = {}
    for i, n in enumerate(names):
        c = article_color(n)
        out[n] = c if c is not None else PALETTE[i % len(PALETTE)]
    return out


def article_color(name):
    """Article colour of a species or dataset label, or None."""
    from Exemples.Study_EDNiX.analysis import paper_style as _ps
    if name in _ps.SPECIES_COLORS:
        return _ps.SPECIES_COLORS[name]
    sp = _ps.dataset_info(name)[0] or next(
        (s for s in _ps.SPECIES_COLORS if str(name).startswith(s)), None)
    return _ps.SPECIES_COLORS.get(sp)


# Species-level colour map: the article's species colours.
SPECIES_COLORS = _palette_color_map(
    ["Human", "Macaque", "Marmoset", "Mouselemur", "Rat", "Mouse", "Dog"])


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


def _channel_block(name):
    """0 for a plain region or the 'high' channel of a two-channel profile
    ('x|pos', 'x|q80'), 1 for the 'low' channel ('x|neg', 'x|q20')."""
    ch = str(name).rpartition("|")[2]
    if ch == "neg":
        return 1
    if len(ch) > 1 and ch[0] == "q" and ch[1:].isdigit():
        return 1 if int(ch[1:]) < 50 else 0
    return 0


def _value_order(values, names):
    """Indices sorting `names` by increasing value (NaN last); with a
    two-channel profile ('x|pos'/'x|neg' or 'x|q80'/'x|q20'), all 'high'
    channels first then all 'low' ones, each block sorted -- so a profile
    curve never jumps between the two channels."""
    key = np.where(np.isfinite(values), values, np.inf)
    block = np.array([_channel_block(n) for n in names])
    return np.lexsort((key, block))


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
    order = _value_order(ref_vec, common_base_names)
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
