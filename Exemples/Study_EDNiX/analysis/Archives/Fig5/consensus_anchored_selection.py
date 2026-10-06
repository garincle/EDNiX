"""
consensus_anchored_selection.py
----------------------------------
Selects runs against a CROSS-SPECIES CONSENSUS pattern instead of a
single anchor species (the human-anchored approach in
human_anchored_phylogeny). Two orthogonal options control how defensible
the result is against the "you found what you were looking for"
objection:

  leave_one_out : when True, the consensus each species is scored
      against is rebuilt EXCLUDING that species. That species therefore
      never contributed to the criterion that selects its own runs, so
      any resulting cross-species structure is an out-of-sample
      prediction rather than something guaranteed by construction. This
      is the setting that answers the circularity objection; False is
      only for seeing how much the answer depends on it.

  consensus_source : "all_runs" (default) builds each species'
      contribution from EVERY run it has -- no selection anywhere, so the
      consensus is contaminated only by unsystematic noise, never by a
      prior choice. "reference" instead uses each species' own
      medoid/kcore/median reference, which is cleaner per species but
      re-introduces a selection step upstream of the consensus.

Why z-scoring matters here: species differ a lot in overall correlation
amplitude. Averaging raw profiles lets a high-amplitude species dominate
the consensus regardless of shape. Each species' profile is therefore
z-scored across regions before combining, so the consensus reflects the
common SHAPE, which is what "conserved across evolution" should mean.

On thresholding near-zero regions (threshold_abs, default None = keep
everything): a region sitting near zero IS information -- "the
retrosplenial is NOT coupled to motor cortex in any species" is a real
constraint, and a run showing strong coupling there is genuinely
atypical. Dropping those regions also costs statistical power you cannot
spare (~18-25 common regions to begin with) and adds one more arbitrary
parameter to justify. The option exists so you can test sensitivity, not
because thresholding is recommended.
"""

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from Exemples.Study_EDNiX.analysis.Archives.Fig5.human_anchored_phylogeny import (
    SPECIES_COLORS, _ordered_common_values,
)

opj = os.path.join


# =============================================================================
# 1) Per-species contribution to the consensus
# =============================================================================

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


# =============================================================================
# 2) The consensus itself
# =============================================================================

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
        by_label[label] = _zscore_values(vals) if zscore else np.asarray(vals, dtype=float)

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
        if zscore:
            vals = _zscore_values(vals)
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


# =============================================================================
# 3) Driver
# =============================================================================

def run_consensus_selection(sri, launcher, profiles_by_label, common_base_names,
                            n_select=5, leave_one_out=True, consensus_source="all_runs",
                            method="kcore", k=5, zscore=True, threshold_abs=None,
                            similarity="correlation", species_of_label=None,
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
    value_col = value_col or launcher.VALUE_COL
    metric = metric or launcher.METRIC

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

    selections, rankings = {}, {}
    consensus_by_label, used_labels_by_label = {}, {}

    for label, profiles in profiles_by_label.items():
        if verbose:
            print(f"\n  {label}:")
        consensus, used = build_consensus_pattern(
            species_profiles, common_base_names, value_col=value_col,
            exclude_label=label if leave_one_out else None,
            zscore=zscore, threshold_abs=threshold_abs,
            species_of_label=species_of_label, verbose=verbose)
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
        vec = _zscore_values(vals) if zscore else np.asarray(vals, dtype=float)
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
    ax.set_ylabel(f"{'z-scored ' if zscore else ''}{value_col}", fontsize=10)
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


def plot_intra_vs_inter_species_distance(profiles_by_label, common_base_names,
                                         species_of_label, value_col="median_r",
                                         zscore=True, title=None, output_path=None,
                                         dpi=300, max_pairs=20000, random_seed=0):
    """
    FIGURE: distribution of run-to-run correlation for pairs WITHIN the
    same species vs pairs ACROSS species.

    This is the sanity check the whole cross-species premise rests on: if
    within-species pairs are not clearly more similar than across-species
    pairs, then "species" is not the dominant structure in these profiles
    and any phylogenetic reading is on thin ice. Reported as two
    histograms plus a Mann-Whitney test, so the answer is quantitative
    rather than visual.

    max_pairs: pairs are subsampled above this count (all-pairs grows
    quadratically and gets slow past a few hundred runs); the subsample
    is drawn separately for the within and across sets so neither is
    starved.
    """
    rng = np.random.default_rng(random_seed)

    vecs, labels_of_vec = [], []
    for label, profiles in profiles_by_label.items():
        sp = species_of_label.get(label, label)
        for key, prof in profiles.items():
            vals = _ordered_common_values(prof, common_base_names, value_col)
            if not np.isfinite(vals).any():
                continue
            vecs.append(_zscore_values(vals) if zscore else np.asarray(vals, dtype=float))
            labels_of_vec.append(sp)

    n = len(vecs)
    if n < 4:
        print("  [plot] too few runs for a within/across comparison")
        return None, None

    arr = np.vstack(vecs)
    sp_arr = np.array(labels_of_vec)

    idx_i, idx_j = np.triu_indices(n, k=1)
    if len(idx_i) > max_pairs:
        sel = rng.choice(len(idx_i), size=max_pairs, replace=False)
        idx_i, idx_j = idx_i[sel], idx_j[sel]

    within, across = [], []
    for i, j in zip(idx_i, idx_j):
        a, b = arr[i], arr[j]
        both = np.isfinite(a) & np.isfinite(b)
        if both.sum() < 3:
            continue
        r = float(np.corrcoef(a[both], b[both])[0, 1])
        (within if sp_arr[i] == sp_arr[j] else across).append(r)

    if len(within) < 2 or len(across) < 2:
        print("  [plot] not enough within- or across-species pairs")
        return None, None

    from scipy import stats as _st
    u, p = _st.mannwhitneyu(within, across, alternative="two-sided")
    auc = u / (len(within) * len(across))

    fig, ax = plt.subplots(figsize=(8, 5))
    bins = np.linspace(-1, 1, 60)
    ax.hist(within, bins=bins, alpha=0.65, color="#009E73",
           label=f"within species (n={len(within)}, median={np.median(within):.3f})")
    ax.hist(across, bins=bins, alpha=0.65, color="#D55E00",
           label=f"across species (n={len(across)}, median={np.median(across):.3f})")
    ax.axvline(np.median(within), color="#009E73", lw=2, ls="--")
    ax.axvline(np.median(across), color="#D55E00", lw=2, ls="--")
    ax.set_xlabel("run-to-run profile correlation", fontsize=10)
    ax.set_ylabel("pairs", fontsize=10)
    ax.set_title(title or f"Within- vs across-species similarity "
                          f"(AUC={auc:.3f}, p={p:.2e})", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")

    stats = dict(n_within=len(within), n_across=len(across),
                median_within=float(np.median(within)),
                median_across=float(np.median(across)),
                auc=float(auc), p_value=float(p))
    print(f"    within median={stats['median_within']:.4f}, "
         f"across median={stats['median_across']:.4f}, AUC={auc:.3f}, p={p:.2e}")
    return fig, stats


# =============================================================================
# 4) Balanced kept vs not-kept statistics
# =============================================================================

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
    rng = np.random.default_rng(random_seed)
    rows = []

    for label, profiles in profiles_by_label.items():
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
            if zscore:
                ref_vec = _zscore_values(ref_vec)

        kept_scores, notkept_scores = [], []
        for key, prof in profiles.items():
            vals = _ordered_common_values(prof, common_base_names, value_col)
            if zscore:
                vals = _zscore_values(vals)
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
                  p_value=float(np.median(pvals)),
                  p_frac_sig=float(np.mean(np.asarray(pvals) < 0.05)),
                  n_draws=len(pvals))
        rows.append(row)
        if verbose:
            print(f"  {label}: kept={n_k} vs {n_k} drawn from {n_n} | "
                 f"delta={row['delta_median']:+.4f}, AUC={row['auc']:.3f}, "
                 f"median p={row['p_value']:.4f} "
                 f"({row['p_frac_sig']*100:.0f}% of draws p<0.05)")

    return pd.DataFrame(rows)


# =============================================================================
# 5) Score distributions, external QC validation, phylogeny permutation
# =============================================================================

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
    rng = np.random.default_rng(random_seed)
    metrics = metrics or ["func_TSNR", "func_mean_fd", "func_mean_dvars"]
    rows = []

    for label, profiles in profiles_by_label.items():
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

    Returns dict(observed_r, p_perm, null_r, n_species) or None if fewer
    than 4 species have both a distance and a divergence time.
    """
    from scipy import stats as _st
    rng = np.random.default_rng(random_seed)

    df = distance_df.copy()
    df["divergence_mya"] = df["species"].map(divergence_times)
    df = df[df["divergence_mya"].notna() & df["distance_to_human"].notna()]
    df = df[df["divergence_mya"] > 0]  # a species at distance 0 from itself carries no information
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
