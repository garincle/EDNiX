"""
run_consensus_full_pipeline.py
---------------------------------
End-to-end driver for the consensus-anchored approach. Everything lands
in ONE output tree, and every heavy step is optional so a fast
exploratory pass costs nothing but profile maths.

Output layout under output_dir/:
    00_log.txt                        full console log of the run
    01_common_regions.xlsx            the fixed step-1 region list
    02_consensus_selection_<tag>.xlsx per-run rankings + selection summary
    03_run_correlations.xlsx          every run: r to consensus, selected or not
    04_distance_to_consensus.xlsx     per-species group-map distances
    figures/
        01_consensus_pattern.png      the consensus each species was scored against
        02_all_species_vs_consensus.png
        03_per_species/<label>.png    species ref + consensus + all runs
        04_distance_vs_divergence.png phylogenetic test
    jpgs/<label>/{selected,not_selected}/   SBA previews per run
    (group maps and surfaces are written by AFNI under RESULTS_BASE_DIR,
     not here -- their paths are recorded in 04_distance_to_consensus.xlsx)

The pipeline never silently skips: anything that fails for a BIDS is
printed with the reason and that BIDS is excluded from later steps,
rather than being quietly absent from the figures.
"""

import os
import sys

import numpy as np
import pandas as pd

from Exemples.Study_EDNiX.analysis.Archives.Fig5.human_anchored_phylogeny import (
    get_common_regions_from_profiles, selected_keys_to_df_class,
    run_group_analysis_for_selection, group_map_profile,
    test_phylogenetic_signal, plot_distance_vs_divergence,
    plot_all_bids_refs_vs_human, plot_bids_ref_with_individuals,
    build_reference_for_profiles, )
from Exemples.Study_EDNiX.analysis.Archives.Fig5.consensus_anchored_selection import (
    run_consensus_selection, plot_consensus_pattern,
    plot_consensus_by_species, plot_intra_vs_inter_species_distance,
    balanced_kept_vs_notkept, _zscore_values,
    plot_score_distributions, validate_selection_against_qc,
    permutation_test_phylogeny, plot_phylogeny_permutation,
)

opj = os.path.join


class _Tee:
    """Mirror stdout to a log file so 00_log.txt holds exactly what the
    console showed -- the debug prints are the main way to tell whether a
    step did what it claimed, so they belong in the saved output too."""

    def __init__(self, path):
        self.file = open(path, "w", encoding="utf-8")
        self.stdout = sys.stdout

    def write(self, data):
        self.stdout.write(data)
        self.file.write(data)

    def flush(self):
        self.stdout.flush()
        self.file.flush()

    def close(self):
        self.file.close()


def copy_selection_jpgs(sri, cfg, profiles, selected_keys, output_dir, verbose=True):
    """
    Copy each run's SBA preview jpg into selected/ or not_selected/,
    reusing SBA_Region_Profile_Index.copy_run_qc_images (which already
    knows the real on-disk jpg convention) by handing it a df_class where
    "kept" means "selected by the consensus".
    """
    df_class = selected_keys_to_df_class(profiles, selected_keys)
    try:
        sri.copy_run_qc_images(df_class, cfg["bids_dir"], output_dir)
    except Exception as e:
        print(f"    [WARN] jpg copy failed: {e}")


def run_full_consensus_pipeline(
        sri, launcher, output_dir,
        bids_configs=None,
        n_select=5,
        leave_one_out=True,
        consensus_source="all_runs",
        similarity="correlation",
        zscore=True,
        threshold_abs=None,
        reference_method="kcore",
        k_ref=5,
        run_group_analysis=True,
        render_surface=False,
        copy_jpgs=True,
        make_figures=True,
        divergence_times=None,
        df_qc_by_label=None,
        n_perm_phylogeny=10000,
):
    """
    Full pipeline. Heavy steps are individually switchable:

    run_group_analysis: real AFNI 3dttest++ per BIDS on the selected runs.
        False skips it AND everything that depends on group maps (the
        distance/phylogeny step), leaving selection + profile figures.
    render_surface: wb_command surface rendering on top of each group
        map. Only meaningful with run_group_analysis=True; off by default
        since it is the slowest step and only needed for final figures.
    copy_jpgs / make_figures: cheap, on by default.

    Returns a dict with every intermediate: common_regions, profiles,
    selection (the run_consensus_selection result), distance_df,
    test_result, group_suffixes.
    """
    bids_configs = bids_configs if bids_configs is not None else launcher.BIDS_CONFIGS
    os.makedirs(output_dir, exist_ok=True)
    fig_dir = opj(output_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    tee = _Tee(opj(output_dir, "00_log.txt"))
    sys.stdout = tee
    try:
        return _run(sri, launcher, output_dir, fig_dir, bids_configs, n_select,
                    leave_one_out, consensus_source, similarity, zscore, threshold_abs,
                    reference_method, k_ref, run_group_analysis, render_surface,
                    copy_jpgs, make_figures, divergence_times, df_qc_by_label,
                    n_perm_phylogeny)
    finally:
        sys.stdout = tee.stdout
        tee.close()
        print(f"\n[log] {opj(output_dir, '00_log.txt')}")


def _run(sri, launcher, output_dir, fig_dir, bids_configs, n_select,
         leave_one_out, consensus_source, similarity, zscore, threshold_abs,
         reference_method, k_ref, run_group_analysis, render_surface,
         copy_jpgs, make_figures, divergence_times, df_qc_by_label,
         n_perm_phylogeny):

    tag = (f"loo{int(leave_one_out)}_{consensus_source}_{similarity}"
           f"_z{int(zscore)}_thr{threshold_abs if threshold_abs is not None else 'none'}")
    print(f"\n{'#' * 70}\n  CONSENSUS PIPELINE -- {tag}\n{'#' * 70}")

    # ---- STEP 1: profiles for every BIDS -------------------------------
    print(f"\n{'=' * 70}\n  STEP 1: per-run profiles\n{'=' * 70}")
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
                print(f"  {label}: {len(profs)} run profile(s)")
            else:
                print(f"  [SKIP] {label}: no run profiles computed")
        except Exception as e:
            print(f"  [ERROR] {label}: {e}")
    if not profiles_by_label:
        raise ValueError("No BIDS produced run profiles -- nothing to do.")

    # ---- STEP 2: fixed common-region list -------------------------------
    print(f"\n{'=' * 70}\n  STEP 2: common regions (fixed from here on)\n{'=' * 70}")
    common = get_common_regions_from_profiles(profiles_by_label,
                                              value_col=launcher.VALUE_COL)
    if not common:
        raise ValueError("No common regions across BIDS profiles.")
    pd.DataFrame({"common_region": common}).to_excel(
        opj(output_dir, "01_common_regions.xlsx"), index=False)

    species_of_label = {c["label"]: c["species"] for c in bids_configs}

    # ---- STEP 3: consensus + selection ----------------------------------
    print(f"\n{'=' * 70}\n  STEP 3: consensus + run selection\n{'=' * 70}")
    sel = run_consensus_selection(
        sri, launcher, profiles_by_label, common,
        n_select=n_select, leave_one_out=leave_one_out,
        consensus_source=consensus_source, similarity=similarity,
        method=reference_method, k=k_ref, zscore=zscore,
        threshold_abs=threshold_abs, species_of_label=species_of_label,
        output_dir=None, verbose=True)
    selections = sel["selections"]
    if not selections:
        raise ValueError("Consensus selection produced nothing.")

    # per-run correlation table, with selected flag
    rows = []
    species_lut = {c["label"]: c["species"] for c in bids_configs}
    for label, ranking in sel["rankings"].items():
        chosen = {tuple(str(x) for x in k) for k in selections.get(label, [])}
        for r in ranking.itertuples():
            key = (str(r.subject), str(r.session), str(r.run), str(r.seed))
            rows.append(dict(label=label, species=species_lut.get(label),
                            subject=r.subject, session=r.session, run=r.run,
                            seed=r.seed, score=r.pearson_r, n_regions=r.n_regions,
                            is_selected=key in chosen))
    corr_table = pd.DataFrame(rows)
    with pd.ExcelWriter(opj(output_dir, f"02_consensus_selection_{tag}.xlsx"),
                        engine="openpyxl") as w:
        pd.DataFrame([
            dict(label=lab, n_runs=len(profiles_by_label[lab]),
                 n_selected=len(selections.get(lab, [])),
                 consensus_from=", ".join(sel["used_labels_by_label"].get(lab, [])),
                 n_consensus_regions=int(np.isfinite(sel["consensus_by_label"][lab]).sum()))
            for lab in sel["rankings"]
        ]).to_excel(w, sheet_name="summary", index=False)
        for lab, rdf in sel["rankings"].items():
            rdf.head(200).to_excel(w, sheet_name=f"rank_{lab[:24]}", index=False)
    corr_table.to_excel(opj(output_dir, "03_run_correlations.xlsx"), index=False)
    print(f"  [saved] 02_consensus_selection_{tag}.xlsx, 03_run_correlations.xlsx")

    # ---- STEP 3b: balanced kept vs not-kept -----------------------------
    print(f"\n{'=' * 70}\n  STEP 3b: kept vs not-kept (size-matched)\n{'=' * 70}")
    balanced_df = balanced_kept_vs_notkept(
        profiles_by_label, selections, common,
        consensus_by_label=sel["consensus_by_label"],
        value_col=launcher.VALUE_COL, zscore=zscore,
        species_of_label=species_of_label)
    if not balanced_df.empty:
        balanced_df.to_excel(opj(output_dir, "05_kept_vs_notkept_balanced.xlsx"), index=False)
        print(f"  [saved] 05_kept_vs_notkept_balanced.xlsx")

    # ---- STEP 3c: external QC validation (non-circular) -----------------
    if df_qc_by_label:
        print(f"\n{'=' * 70}\n  STEP 3c: external QC validation (size-matched)\n{'=' * 70}")
        try:
            from Exemples.Study_EDNiX.analysis.Fig5.external_qc_validation import merge_qc_with_classification
            qc_df = validate_selection_against_qc(
                profiles_by_label, selections, df_qc_by_label,
                merge_qc_with_classification)
            if not qc_df.empty:
                qc_df.to_excel(opj(output_dir, "06_external_qc_validation.xlsx"), index=False)
                print(f"  [saved] 06_external_qc_validation.xlsx")
        except Exception as e:
            print(f"  [WARN] external QC validation failed: {e}")
    else:
        print(f"\n  [SKIP] STEP 3c: no df_qc_by_label given "
             f"(pass it to get the non-circular QC validation)")

    # ---- STEP 4: jpgs ---------------------------------------------------
    if copy_jpgs:
        print(f"\n{'=' * 70}\n  STEP 4: SBA preview jpgs\n{'=' * 70}")
        cfg_by_label = {c["label"]: c for c in bids_configs}
        for label, keys in selections.items():
            print(f"  {label}:")
            copy_selection_jpgs(sri, cfg_by_label[label], profiles_by_label[label],
                                keys, opj(output_dir, "jpgs", label))

    # ---- STEP 5: group 3dttest (+ optional surfaces) --------------------
    group_suffixes, group_profiles = {}, {}
    cfg_by_label = {c["label"]: c for c in bids_configs}
    if run_group_analysis:
        print(f"\n{'=' * 70}\n  STEP 5: group 3dttest"
              f"{' + surface rendering' if render_surface else ''}\n{'=' * 70}")
        failed_ttest, failed_readback = [], []
        for label, keys in selections.items():
            cfg = cfg_by_label[label]
            seed = cfg["seeds"][0]
            suffix = f"{seed}_Consensus_{tag}"
            df_class = selected_keys_to_df_class(profiles_by_label[label], keys)
            ok = run_group_analysis_for_selection(
                sri, cfg, seed, df_class, suffix, render_surface=render_surface)
            if not ok:
                failed_ttest.append(label)
                continue
            prof = group_map_profile(sri, launcher, cfg, seed, suffix, common)
            if prof is None:
                failed_readback.append(label)
                continue
            group_suffixes[label] = suffix
            group_profiles[label] = prof

        # Explicit accounting: without this, a BIDS whose 3dttest failed or
        # whose group map could not be read back just vanished from the
        # later steps, and the phylogeny figure silently never appeared
        # with no indication of why.
        print(f"\n  group maps usable: {len(group_profiles)}/{len(selections)}")
        if failed_ttest:
            print(f"  [WARN] 3dttest failed for: {failed_ttest}")
        if failed_readback:
            print(f"  [WARN] group map written but could NOT be read back "
                 f"(missing file, or no atlas matching its grid) for: {failed_readback}")
    else:
        print(f"\n  [SKIP] STEP 5 (run_group_analysis=False) -- no group maps, "
             f"so STEP 6's distance/phylogeny is skipped too")

    # ---- STEP 6: distance to the consensus-anchored human/group map -----
    distance_df, test_result = pd.DataFrame(), None
    if run_group_analysis and not group_profiles:
        print(f"\n  [WARN] STEP 6 skipped: run_group_analysis was True but NO group map "
             f"could be used (see the warnings just above) -- this is why "
             f"04_distance_vs_divergence.png is not produced.")
    if group_profiles:
        print(f"\n{'=' * 70}\n  STEP 6: distance between species group maps\n{'=' * 70}")
        # Reference = the consensus itself, expressed as the median of
        # every species' group map. Unlike anchoring on one species, no
        # single dataset defines the target here.
        stack = np.vstack([
            np.asarray([
                (p.assign(_b=p["base_name"].str.lower().str.strip())
                  .groupby("_b")[launcher.VALUE_COL].mean().get(b, np.nan))
                for b in common], dtype=float)
            for p in group_profiles.values()])
        with np.errstate(invalid="ignore"):
            ref_vec = np.nanmedian(stack, axis=0)

        rows = []
        for label, prof in group_profiles.items():
            vals = np.asarray([
                (prof.assign(_b=prof["base_name"].str.lower().str.strip())
                     .groupby("_b")[launcher.VALUE_COL].mean().get(b, np.nan))
                for b in common], dtype=float)
            both = np.isfinite(vals) & np.isfinite(ref_vec)
            if both.sum() < 3:
                print(f"  [WARN] {label}: too few overlapping regions")
                continue
            r = float(np.corrcoef(vals[both], ref_vec[both])[0, 1])
            rows.append(dict(label=label, species=species_lut.get(label),
                            distance_to_human=1 - r, pearson_r=r,
                            n_regions=int(both.sum()),
                            group_map_suffix=group_suffixes[label]))
            print(f"  {label}: r to consensus = {r:.4f}")
        distance_df = pd.DataFrame(rows)
        distance_df.to_excel(opj(output_dir, "04_distance_to_consensus.xlsx"), index=False)

        if len(distance_df) >= 3:
            test_result = test_phylogenetic_signal(
                distance_df, divergence_times=divergence_times, exclude_human=False)

    # ---- STEP 7: figures ------------------------------------------------
    if make_figures:
        print(f"\n{'=' * 70}\n  STEP 7: figures\n{'=' * 70}")
        plot_consensus_pattern(
            sel["consensus_by_label"], common,
            title=f"Consensus pattern, per BIDS variant ({tag})",
            output_path=opj(fig_dir, "01_consensus_pattern.png"))

        # one line per SPECIES (BIDS of a species collapsed first), so
        # species with several datasets don't count several times
        plot_consensus_by_species(
            sel["species_profiles"], common, species_of_label,
            value_col=launcher.VALUE_COL, zscore=zscore,
            title=f"Consensus by species ({tag})",
            output_path=opj(fig_dir, "01b_consensus_by_species.png"))

        # is "species" actually the dominant structure in these profiles?
        plot_intra_vs_inter_species_distance(
            profiles_by_label, common, species_of_label,
            value_col=launcher.VALUE_COL, zscore=zscore,
            title=f"Within- vs across-species run similarity ({tag})",
            output_path=opj(fig_dir, "05_within_vs_across_species.png"))

        own_refs = {}
        for label, profs in profiles_by_label.items():
            try:
                own_refs[label], _ = build_reference_for_profiles(
                    sri, profs, method=reference_method, k=k_ref,
                    value_col=launcher.VALUE_COL, metric=launcher.METRIC)
            except Exception as e:
                print(f"  [WARN] {label}: own reference failed: {e}")

        # Confirmed real problem being fixed here: the consensus is
        # z-scored (roughly -1.5..+1.5) while raw species references sit
        # around -0.1..+0.4, so plotting them on one axis made the
        # consensus look like a huge outlier and flattened every species
        # into a band near zero. Z-score the species references the same
        # way so the two are actually comparable -- the figure is about
        # SHAPE agreement, which is exactly what z-scoring isolates.
        consensus_as_profile = None
        if own_refs:
            if zscore:
                for lab, ref in list(own_refs.items()):
                    vals = np.asarray([
                        (ref.assign(_b=ref["base_name"].str.lower().str.strip())
                            .groupby("_b")[launcher.VALUE_COL].mean().get(b, np.nan))
                        for b in common], dtype=float)
                    lut_z = dict(zip(common, _zscore_values(vals)))
                    ref = ref.copy()
                    ref[launcher.VALUE_COL] = (
                        ref["base_name"].str.lower().str.strip().map(lut_z))
                    own_refs[lab] = ref

            any_label = next(iter(own_refs))
            consensus_as_profile = own_refs[any_label].copy()
            med = np.nanmedian(np.vstack(list(sel["consensus_by_label"].values())), axis=0)
            lut = dict(zip(common, med))
            consensus_as_profile[launcher.VALUE_COL] = (
                consensus_as_profile["base_name"].str.lower().str.strip().map(lut))
            plot_all_bids_refs_vs_human(
                own_refs, consensus_as_profile, common,
                value_col=launcher.VALUE_COL,
                title=f"Each species reference vs consensus ({tag})",
                output_path=opj(fig_dir, "02_all_species_vs_consensus.png"))

        per_dir = opj(fig_dir, "03_per_species")
        os.makedirs(per_dir, exist_ok=True)
        for label, profs in profiles_by_label.items():
            if label not in own_refs or consensus_as_profile is None:
                continue
            plot_bids_ref_with_individuals(
                profs, own_refs[label], consensus_as_profile, common,
                selected_keys=selections.get(label), value_col=launcher.VALUE_COL,
                title=f"{label} - species ref vs consensus ({tag})",
                output_path=opj(per_dir, f"{label}.png"))

        plot_score_distributions(
            sel["rankings"], selections,
            title=f"Consensus score distribution per BIDS ({tag})",
            output_path=opj(fig_dir, "06_score_distributions.png"))

        if not distance_df.empty:
            from human_anchored_phylogeny import DIVERGENCE_TIME_MYA
            perm = permutation_test_phylogeny(
                distance_df, divergence_times or DIVERGENCE_TIME_MYA,
                n_perm=n_perm_phylogeny)
            if perm is not None:
                plot_phylogeny_permutation(
                    perm, title=f"Phylogeny permutation test ({tag})",
                    output_path=opj(fig_dir, "07_phylogeny_permutation.png"))

        if distance_df.empty:
            print(f"  [SKIP] 04_distance_vs_divergence.png -- distance_df is empty "
                 f"(no usable group maps; see STEP 5/6 warnings)")
        elif test_result is None:
            print(f"  [SKIP] 04_distance_vs_divergence.png -- only {len(distance_df)} "
                 f"species with a distance, need >= 3 for the correlation")
        if not distance_df.empty and test_result:
            plot_distance_vs_divergence(
                distance_df, test_result=test_result, divergence_times=divergence_times,
                title=f"Distance to consensus vs divergence time ({tag})",
                output_path=opj(fig_dir, "04_distance_vs_divergence.png"))

    print(f"\n{'#' * 70}\n  DONE -- everything under {output_dir}\n{'#' * 70}")
    return dict(common_regions=common, profiles_by_label=profiles_by_label,
               selection=sel, correlation_table=corr_table,
               balanced_kept_vs_notkept=balanced_df,
               distance_df=distance_df, test_result=test_result,
               group_suffixes=group_suffixes)
