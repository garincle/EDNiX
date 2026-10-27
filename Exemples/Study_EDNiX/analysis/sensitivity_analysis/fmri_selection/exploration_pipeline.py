"""
fig5_pipeline.py
================
ONE entry point for the whole Figure 5 analysis. Replaces the scattered
run_*.py drivers: everything is controlled from the CONFIG block below,
and every output lands in one predictable tree under OUTPUT_ROOT.

WHY THIS EXISTS
---------------
Four concrete problems in the previous organisation, all fixed here:

1. THE UNIT OF ANALYSIS IS NO LONGER THE RUN (`UNIT` in CONFIG, default
   "subject"). Everything -- the per-dataset profile, the cross-species
   consensus, the selection, and the 3dttest++ group maps -- now treats
   one ANIMAL as one sample. Before, a subject scanned 20 times counted
   20 times while a subject scanned once counted once, so a dataset's
   profile described its most-scanned animals, N_SELECT "best runs"
   could all come from a single animal, and 3dttest++ spent its degrees
   of freedom on non-independent scans. Runs are still the raw material:
   they are averaged within each subject, in Python for the profiles
   (median) and by AFNI 3dMean for the group maps.

2. STEP 1 ACTUALLY PRODUCES A MAP. It used to call
   `rerun_balanced_by_run` -- a kept-vs-outlier function -- with one
   arbitrary run flagged as a fake outlier just to satisfy it. That
   function then BALANCED the two groups, so min(N-1, 1) = 1: the
   "all runs" map was computed on a single run, AFNI refused it
   ("[SKIP] ...: only 1 run(s)"), and every species got an empty output
   folder while summary.xlsx still claimed 24/39/80/170 runs. Step 1 now
   passes balance=False and keeps every unit, and the `_unused` folder is
   gone.

3. Two different folder-naming conventions existed at once --
   SBA_group_maps.py (former Figure_5_SBA.py) wrote `Grp_SBA_3dTTEST_All<suffix>` (category glued
   on) while rerun_balanced_by_run wrote `Grp_SBA_3dTTEST_<suffix>`.
   Lookups guessed wrong and silently reported "group map not found at
   None". Here, ONE convention is used everywhere and the resolved path
   is printed and recorded for every map.

4. Outputs were spread across several roots with opaque names. Here the
   tree is fixed (see OUTPUT LAYOUT), every figure gets a .txt legend,
   and a methods.txt is written with the exact parameters used.

PART 1 (selection) NOW LIVES IN Statistics/Group_fMRI/run_selection/ (moved verbatim
2026-10-27: selection.py, consensus.py, units.py, profiles.py, group_maps.py ...) and is
imported back here; run(..., part=...) works as before.

DEPENDS ON
----------
- fig5_analysis.py   (unchanged: it works on generic profile dicts)
- fig5_units.py      (new: the run -> subject aggregation, and the
                      step-1 figure read off the AFNI group maps)
- SBA_Region_Profile_Index.py, patched by patch_SBA_unit_level.py
  (adds `unit=` and `balance=` to rerun_balanced_by_run). Running this
  pipeline against an unpatched module raises a TypeError immediately,
  by design -- silently falling back to run-level would undo point 1.

OUTPUT LAYOUT (under OUTPUT_ROOT/)
----------------------------------
    methods.txt                     parameters + method text for the paper
    00_log.txt                      full console log
    01_all_units/                   group maps from ALL units, per BIDS
        summary.xlsx                what actually entered each map
        map_index.xlsx              every map's resolved path
        F0_species_profiles_afni.png   species profiles read off those maps,
        F0_species_profiles_afni.txt   with the cross-species core shaded
    02_selection/                   the consensus selection itself
        selection.xlsx              per-unit scores, kept flag
        common_regions.xlsx
        region_invariance.xlsx
        unit_counts.xlsx            units available per dataset
        balanced_stats.xlsx         kept vs not-kept, size-matched (descriptive, no p)
        external_qc.xlsx            non-circular QC validation (optional)
    03_group_maps/                  unit-level 3dttest for kept AND not-kept
        map_index.xlsx              every map's resolved path
    04_jpg/                         SBA previews
    05_surfaces/                    every surface .png, renamed unambiguously
        <label>__<group>__<stat>.png
    06_figures/                     all figures + one .txt legend each
"""

import os
import sys
import glob
import shutil
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

import Exemples.Study_EDNiX.analysis.fmri_fingerprints.profile_stats as fa   # pure numpy/pandas/matplotlib
from Statistics.Group_fMRI.run_selection import units as fu      # unit aggregation + step-1 figure
import Exemples.Study_EDNiX.analysis.fmri_fingerprints.phylo as fp      # Mantel / Blomberg K / PGLS
import Exemples.Study_EDNiX.analysis.fmri_fingerprints.fingerprints as ff  # F6ab, fingerprints, F15-F18
import Exemples.Study_EDNiX.analysis.fmri_fingerprints.paper_panels as fpp         # panels of Figures 5 and 6
import Exemples.Study_EDNiX.analysis.fmri_fingerprints.supp_panels as fsp    # print-size panels, Figs 4-5 / Supp 7-11
from Statistics.Group_fMRI.run_selection.selection import (  # part 1, moved verbatim 2026-10-27
    CONFIG,
    FIGURE_KEYS,
    LEGACY_CONFIG_KEYS,
    PARTS,
    STATE_FILE,
    _DERIVED_KEYS,
    _Tee,
    _archive_previous_outputs,
    _check_map_inputs,
    _close_log,
    _copy_selection_jpgs,
    _criterion_in_space,
    _delete_old_map_files,
    _expected_balanced_subsample,
    _kept_pair_reusable,
    _map_input_runs,
    _map_suffix,
    _merge_figure_config,
    _norm_id,
    _norm_sentence,
    _open_log,
    _resolve_map,
    _run_region_indices,
    _run_selection,
    _seed_regions_to_exclude,
    _state_path,
    _tree,
    _validate_config,
    _vec_to_profile,
    _write_methods,
    load_state,
    save_state,
    step1_all_units_maps,
    step3_group_maps,
    step3b_final_group_maps,
    step5_collect_surfaces)
from Exemples.Study_EDNiX.analysis.fmri_fingerprints.context import (  # moved verbatim 2026-10-27
    _DIVERGENCE_DEFAULT,
    _kept_references,
    _load_time_tree)

# ═══════════════════════════════════════════════════════════════════════════
# CONFIG -- everything controllable lives here
# ═══════════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════════
# Small helpers
# ═══════════════════════════════════════════════════════════════════════════


def _write_legend(path, title, what, how, reading):
    """Write the .txt legend that accompanies each figure."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{title}\n{'=' * len(title)}\n\n")
        f.write(f"WHAT IS PLOTTED\n---------------\n{what}\n\n")
        f.write(f"HOW IT WAS COMPUTED\n-------------------\n{how}\n\n")
        f.write(f"HOW TO READ IT\n--------------\n{reading}\n")


def _write_profile_reliability_legend(png, rel, tests, cfg):
    lines = []
    for r in rel.sort_values(["label", "method"]).itertuples():
        lines.append(f"{r.label} / {r.method}: discriminability {r.discriminability:.3f} "
                     f"(p vs chance {r.p_discriminability:.3f}); r inter-run "
                     f"{r.r_inter_run:.3f}, same animal {r.r_same_animal:.3f}, different "
                     f"animals {r.r_inter_animal:.3f} (gap {r.gap_animal:+.3f}); "
                     f"{r.n_runs} runs, {r.n_animals} animals, {r.n_values} values")
    for k, d in tests.items():
        for m, t in d.items():
            lines.append(f"PAIRED {k}: {m} - median_r = {t['mean_diff']:+.4f}, higher in "
                         f"{t['n_better']}/{t['n']} datasets, p = {t['p']:.3f} (exact sign-flip, "
                         f"two-sided)")
    _write_legend(
        os.path.splitext(png)[0] + ".txt",
        "F22 -- Which construction of the profile is the most reliable?",
        "Every run of every dataset with several runs per animal, BEFORE any selection. For "
        "each profile method: left, the discriminability (Bridgeford et al. 2021 PLoS Comput "
        "Biol 17:e1009279) = probability that a run is closer (1 - r between profiles) to "
        "another run of the same animal than to a run of another animal, 0.5 = chance; right, "
        "median r between runs of the same animal minus median r between runs of different "
        "animals. One line per dataset, black = mean. Methods: median_r (median of the "
        f"region's voxels); fraction_posneg (% of the region's voxels in the top / bottom "
        f"{100 * cfg['TAIL_Q']:g} % of the map, percentiles over the common regions); "
        f"percentiles ({100 - 100 * cfg['TAIL_Q']:g}th / {100 * cfg['TAIL_Q']:g}th "
        f"percentile of the region's voxels). All three are computed from the same voxels "
        f"of the same runs and normalised as in the analysis (ZSCORE={cfg['ZSCORE']!r}).",
        "p of a dataset: animal labels shuffled across its runs (500 permutations, "
        "one-sided). Method comparison: difference with median_r per dataset, exact sign-flip "
        "over datasets (2^n), two-sided. The criterion was fixed before looking at the "
        "cross-species results and does not involve them.\n" + "\n".join(lines),
        "Discriminability is rank-based, hence comparable between profiles of different "
        "length (16 values for median_r, 32 for the two-channel methods); the r values "
        "themselves are not. A method that is more discriminable measures the animal more "
        "reliably from one run. For Mouse and Marmoset_Hori, runs of one animal span arousal "
        "states, which lowers every method alike.")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 1 -- group map from ALL units (no selection)
# ═══════════════════════════════════════════════════════════════════════════


def step1b_species_profiles_figure(sri, cfg, paths, afni_profiles,
                                   common, region_mask, invariance):
    """The two figures drawn off the all-individuals group maps: F0a with
    one line per BIDS, F0 with one line per species, both with the
    regions that define the cross-species profile shaded.

    Reading the profile back OFF the volumetric group map (rather than
    re-deriving it from the per-unit profiles, which figure F1 already
    does) is the point: it is the only check that the AFNI group
    statistics and the Python profile arithmetic describe the same
    pattern. If F0 and F1 disagree in shape, one of the two is wrong and
    it is worth knowing before either goes in a paper.

    afni_profiles is passed in already computed (the common regions are
    now derived from these same profiles, so they must exist before this
    is called) -- this function no longer reads the maps itself.
    """
    if not afni_profiles:
        print("  [SKIP] F0: no group map could be profiled")
        return None

    # F0a -- one line per BIDS: the reproducibility check. Two datasets of
    # one species come from different labs, scanners and protocols, so
    # their agreement (or lack of it) bounds how much of any cross-species
    # difference can be attributed to acquisition rather than biology.
    png_bids = os.path.join(paths["all_units"], "F0a_bids_profiles_afni.png")
    fu.plot_species_profiles_from_afni(
        afni_profiles, common, region_mask=region_mask, invariance=invariance,
        value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"], group_by="label",
        title=(f"Per-BIDS profiles from the all-{cfg['UNIT']}s AFNI group maps "
               f"(lvl {cfg['ATLAS_LEVEL']}, seed {cfg['SEED']})"),
        output_path=png_bids, analysis_fa=fa)

    png = os.path.join(paths["all_units"], "F0_species_profiles_afni.png")
    fu.plot_species_profiles_from_afni(
        afni_profiles, common, region_mask=region_mask, invariance=invariance,
        value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"], group_by="species",
        title=(f"Species profiles from the all-{cfg['UNIT']}s AFNI group maps "
               f"(lvl {cfg['ATLAS_LEVEL']}, seed {cfg['SEED']})"),
        output_path=png, analysis_fa=fa)

    n_sel = int(np.asarray(region_mask, dtype=bool).sum()) if region_mask is not None else len(common)
    _write_legend(
        os.path.join(paths["all_units"], "F0_species_profiles_afni.txt"),
        "F0/F0a -- Species and per-BIDS profiles from the all-units group maps",
        f"F0a draws one line per BIDS (species colour, line style per dataset) so "
        f"the datasets of one species can be read against each other; F0 collapses "
        f"them into one line per species. "
        f"Top panel: one line per species, over the {len(common)} regions common "
        f"to every dataset, read directly off the step-1 AFNI group maps (every "
        f"{cfg['UNIT']}, no selection). Black = the median across species. The "
        f"green shading and the bold region names mark the "
        f"{n_sel}/{len(common)} regions that define the cross-species profile, "
        f"i.e. the ones the selection was restricted to. Bottom panel: each "
        f"region's invariance across species, green where it was retained.",
        f"Common regions are taken from these group maps themselves: present at "
        f"atlas level {cfg['ATLAS_LEVEL']} in every dataset's atlas AND sampled by "
        f"every group map. "
        f"Each dataset's group map is a 3dttest++ over its {cfg['UNIT']}-level "
        f"Fisher-z maps (runs averaged within {cfg['UNIT']} by AFNI 3dMean). "
        f"The map was then summarised into a regional profile with the same "
        f"{cfg['VALUE_COL']}-per-region routine used for the individual "
        f"profiles, datasets of one species collapsed by median, and each "
        f"profile {fa._norm_label(cfg['ZSCORE']).strip() or 'left unscaled'} before "
        f"combining. Regions are ordered by the consensus value. Invariance and "
        f"the retained set come from the profile-level analysis "
        f"(mode={cfg['INVARIANCE_MODE']}, value={cfg['INVARIANCE_VALUE']}).",
        "Read the green band as the answer to 'which regions is the "
        "cross-species claim actually made on'. Inside it, species lines should "
        "track the black line closely -- that is what 'invariant' means, and it "
        "is the region set the selection scored runs against. Outside it, "
        "divergence between species is expected and is precisely why those "
        "regions were excluded. Compare the shape of this figure with F1: F1 is "
        "computed from the Python profiles, F0 from the AFNI volumes, and they "
        "should agree. If they do not, the discrepancy is in the pipeline, not "
        "in the biology.")
    return png


# ═══════════════════════════════════════════════════════════════════════════
# STEP 3 -- group maps for kept AND not-kept, with resolved-path index
# ═══════════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════════
# STEP 3b -- kept / not-kept maps for the FINAL selection round (2..N).
# Own suffix (Kept<N>/NotKept<N>) and own folder, so round 1's maps and the
# all-units maps are never touched.
# ═══════════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════════
# STEP 5 -- collect surface PNGs under unambiguous names
# ═══════════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════════
# STEP 7 -- F6 at several atlas levels, from the SAME already-built
# all-units maps (no AFNI recompute)
# ═══════════════════════════════════════════════════════════════════════════

def step6_multilevel_f6(sri, fa, fu, cfg, paths, bids_configs, all_units_index,
                        members_by_label, levels=(1, 2, 3, 4)):
    """
    Re-summarizes the SAME step-1 all-units cmaps at several atlas
    levels (compute_region_profile just re-groups the same volume by a
    different label hierarchy -- no 3dttest++/ClustSim involved), then
    redoes the leave-one-species-out F6 distance/fit at each level, to
    see whether the phylogenetic signal survives coarser or finer
    segmentation. Returns {level: {"ddf_bids", "ddf_species",
    "test_result", "n_common_regions"}}.
    """
    cfg_by_label = {c["label"]: c for c in bids_configs}
    results = {}
    for level in levels:
        print(f"\n  --- F6 @ atlas_level={level} ---")
        profiles_lvl = {}
        # one map per BIDS: the index may carry several groups (e.g. kept /
        # not-kept) -- a later row would silently overwrite an earlier one
        idx = all_units_index
        if "group" in idx.columns and idx["group"].nunique() > 1:
            keep_group = "kept_final" if "kept_final" in set(idx["group"]) else idx["group"].iloc[0]
            idx = idx[idx["group"] == keep_group]
        for r in idx.itertuples():
            if not r.found or not r.map_path:
                continue
            cb = cfg_by_label.get(r.label)
            if cb is None:
                continue
            try:
                run_keys = []
                if r.label in members_by_label:
                    for rks in members_by_label[r.label].values():
                        run_keys.extend(rks)
                atlas = fu.find_matching_atlas_volume(
                    sri, cb, run_keys, r.map_path, verbose=False,
                    label=f"F6lvl{level}/{r.label}")
                atlas = atlas or sri.pick_representative_atlas_volume(cb["bids_dir"])
                ldf = sri.parse_itk_snap_label_file(cb["label_path"])
                prof = sri.compute_region_profile(r.map_path, atlas, ldf, atlas_level=level)
                profiles_lvl[r.label] = (r.species, prof)
            except Exception as e:
                print(f"    [WARN] {r.label} @ lvl{level}: {e}")

        if not profiles_lvl:
            print(f"    [SKIP] lvl{level}: no profile could be built")
            continue

        common_lvl = fa.get_common_regions_from_profiles(
            {lab: {"map": p} for lab, (sp, p) in profiles_lvl.items()},
            value_col=cfg["VALUE_COL"],
            # ALWAYS True here, unlike the single-level pipeline's
            # COMMON_REGIONS_REQUIRE_NONNAN=False default: at a forced
            # atlas_level, label_df still lists all 276 ITK-SNAP rows
            # regardless of whether a region has any voxels at that
            # specific level (compute_region_profile tries every region
            # at the forced level and returns NaN, not a dropped row, for
            # ones with none there) -- so require_nonnan=False here would
            # count every region as "common" at every level structurally,
            # even ones with no real data at that level, exactly what
            # produced the identical 138-region count across levels 1-4.
            require_nonnan=True, verbose=False)
        print(f"    {len(common_lvl)} common region(s) at level {level}")
        if len(common_lvl) < 3:
            print(f"    [SKIP] lvl{level}: fewer than 3 common regions")
            continue

        ddf_bids_lvl = fa.compute_loo_species_distances(
            profiles_lvl, common_lvl, value_col=cfg["VALUE_COL"],
            zscore=cfg["ZSCORE"], verbose=False)
        if len(ddf_bids_lvl) < 3:
            print(f"    [SKIP] lvl{level}: only {len(ddf_bids_lvl)} usable BIDS")
            continue

        ddf_species_lvl = (ddf_bids_lvl.groupby("species", as_index=False)
                           .agg(distance_to_human=("distance_to_human", "mean"),
                                pearson_r=("pearson_r", "mean")))
        tr_lvl = fa.test_phylogenetic_signal(ddf_species_lvl, exclude_human=False)
        mres_lvl = (fp.run_mantel(profiles_lvl, common_lvl, value_col=cfg["VALUE_COL"],
                                  zscore=cfg["ZSCORE"], label=f" lvl{level}")
                    if cfg["RUN_PHYLO_COMPARATIVE"] else None)
        results[level] = dict(ddf_bids=ddf_bids_lvl, ddf_species=ddf_species_lvl,
                              test_result=tr_lvl, n_common_regions=len(common_lvl),
                              mantel=mres_lvl)
        print(f"    lvl{level}: r={tr_lvl['pearson_r']:.3f}, p={tr_lvl['pearson_p']:.4f}, "
              f"n={tr_lvl['n']} species")

    if not results:
        print("  [SKIP] multilevel F6: no level produced usable results")
        return results

    out_xlsx = os.path.join(paths["group_maps"], "distance_to_consensus_multilevel.xlsx")
    with pd.ExcelWriter(out_xlsx) as xw:
        for level, res in results.items():
            res["ddf_bids"].to_excel(xw, sheet_name=f"lvl{level}_per_bids", index=False)
            res["ddf_species"].to_excel(xw, sheet_name=f"lvl{level}_per_species", index=False)
        mantel_rows = [dict(level=lv, n_regions=r["n_common_regions"],
                            **{f"{meth}_{k}": r["mantel"][f"mantel_{meth}"][k]
                               for meth in ("pearson", "spearman")
                               for k in ("r", "p_greater", "p_two_sided")})
                       for lv, r in results.items() if r.get("mantel")]
        if mantel_rows:
            pd.DataFrame(mantel_rows).to_excel(xw, sheet_name="mantel_by_level", index=False)
    print(f"  [saved] {os.path.basename(out_xlsx)}")

    fa.plot_distance_vs_divergence_multilevel(
        results, output_path=os.path.join(paths["figures"],
                                          "F6_multilevel_segmentation_impact.png"))
    return results


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════════
# The two parts: SELECTION (part 1) and FIGURES (part 2)
# ═══════════════════════════════════════════════════════════════════════════
# Part 1 = everything that depends on the data and on the selection settings:
# run profiles, run QC, AFNI all-units maps, common regions, regional indices of
# every run, F22 table, selection rounds, AFNI kept / not-kept maps, jpgs,
# surface renders, profiles of every group map. It ends by writing
# <OUTPUT_ROOT>/07_state/selection_state.pkl.
# Part 2 = every statistic and figure, from that file only (no AFNI, no NIfTI
# read except STEP 7, which is off by default): it can be rerun as often as
# needed, with other figure settings, in minutes.
# run(..., part="both") does both in one go (the default).


def run(sri, bids_configs, config=None, df_qc_by_label=None, part="both"):
    """Run the Figure 5 pipeline.

    sri          : SBA_Region_Profile_Index module -- the only external
                   dependency (AFNI / nibabel / atlas work). Must be the
                   version patched by patch_SBA_unit_level.py.
    bids_configs : list of per-dataset dicts (label, species, bids_dir,
                   label_path, seeds, optional exclude_lateralized, state).
    config       : dict overriding CONFIG entries
    df_qc_by_label : {label: df_qc slice} to enable the external QC check
    part         : "both" (default) = selection then figures;
                   "selection" = part 1 only (ends with 07_state/selection_state.pkl);
                   "figures" = part 2 only, from that file (only FIGURE_KEYS may
                   differ from the selection run; bids_configs and df_qc_by_label
                   are then taken from the file).
    """
    if part not in PARTS:
        raise ValueError(f"part must be one of {PARTS}, got {part!r}")
    if part == "figures":
        root = (config or {}).get("OUTPUT_ROOT", CONFIG["OUTPUT_ROOT"])
        state = load_state(root)
        cfg = _merge_figure_config(state["cfg"], config)
        paths = _tree(root)
        archived = _archive_previous_outputs(paths, figures_only=True)
        tee, old_cwd = _open_log(paths, "00_log_figures.txt")
        if archived:
            print(f"  Previous figures moved to {archived}")
        print(f"  part 2 (figures) from {_state_path(root)} (saved {state.get('saved')})")
        _load_time_tree(cfg, state["bids_configs"])
        if hasattr(sri, "set_profile_method"):
            sri.set_profile_method(cfg["PROFILE_METHOD"], cfg["TAIL_Q"], keep_sketch=False)
            if cfg["PROFILE_METHOD"] != "median_r":
                sri.set_tail_reference(state["anat_common"])
        try:
            return _run_figures(sri, dict(state, cfg=cfg, paths=paths))
        finally:
            if hasattr(sri, "set_profile_method"):
                sri.set_profile_method("median_r")
            _close_log(tee, old_cwd, paths, "00_log_figures.txt")

    cfg = dict(CONFIG)
    if config:
        cfg.update(config)
    config_warnings = _validate_config(cfg, config)
    if cfg["UNIT"] not in fu.VALID_UNITS:
        raise ValueError(f"UNIT must be one of {fu.VALID_UNITS}, got {cfg['UNIT']!r}")

    paths = _tree(cfg["OUTPUT_ROOT"])
    if cfg["MAP_TAG"] is None:
        import hashlib
        cfg["MAP_TAG"] = "r" + hashlib.md5(
            os.path.abspath(cfg["OUTPUT_ROOT"]).encode()).hexdigest()[:6]
    if cfg["RESULTS_BASE_DIR"] is None:
        cfg["RESULTS_BASE_DIR"] = paths["group_maps"]
    # AFNI/singularity fail with "stale file handle" if the process cwd is
    # deleted under them (seen in the real log): start from a directory this
    # run owns. All pipeline paths are absolute, so this changes nothing else.
    archived = _archive_previous_outputs(paths)
    tee, old_cwd = _open_log(paths, "00_log.txt")
    if archived:
        print(f"  Previous run's figures/tables/log moved to {archived} "
              f"(AFNI maps and *.json balance records left in place for reuse)")
    print(f"  part = {part!r}")
    print(f"  AFNI statmaps stored under RESULTS_BASE_DIR={cfg['RESULTS_BASE_DIR']!r}, "
          f"suffix tag MAP_TAG={cfg['MAP_TAG']!r}, "
          f"{'THRESHOLDED' if cfg['STAT_MAP_USE_THRESHOLDED'] else 'unthresholded'} "
          f"group maps")
    for w in config_warnings:
        print(f"  [CONFIG WARNING] {w}")
    _load_time_tree(cfg, bids_configs)
    # The profile method is a module-level setting of sri: set it for this run
    # and put it back afterwards, so the next run of the same process (control,
    # matrix) starts clean. The voxel values of every region are kept until the
    # regional indices of every run (6d) and F22 are computed.
    if hasattr(sri, "set_profile_method"):
        sri.set_profile_method(cfg["PROFILE_METHOD"], cfg["TAIL_Q"],
                               keep_sketch=cfg.get("PROFILE_KIND", "seed") == "seed")
    print(f"  Profile method: {cfg['PROFILE_METHOD']}"
          + (f" (q = {cfg['TAIL_Q']:g}, value column {cfg['VALUE_COL']})"
             if cfg["PROFILE_METHOD"] != "median_r" else ""))
    try:
        state = _run_selection(sri, cfg, paths, bids_configs, df_qc_by_label)
        save_state(state)
        if part == "selection":
            _write_methods(cfg, paths, sorted(state["profiles_by_label"]),
                           state["n_regions"], state["unit_counts"],
                           final_round=state["final_round"])
            print(f"\n{'#' * 70}\n  PART 1 DONE -- {paths['root']}\n  figures: "
                  f"run(..., part='figures')\n{'#' * 70}")
            return state
        return _run_figures(sri, state)
    finally:
        if hasattr(sri, "set_profile_method"):
            sri.set_profile_method("median_r")
        _close_log(tee, old_cwd, paths, "00_log.txt")


def _write_methods_article(st, cfg, paths, phylo_source=None):
    """Methods_article.txt: the Methods of Figures 5-6 written as for a
    manuscript, with the numbers of this run (datasets, runs, animals,
    regions, selection, tests), followed by the references cited. Rewritten
    at the end of every part 2."""
    import scipy
    import matplotlib
    unit = cfg["UNIT"]
    labels = sorted(st["profiles_by_label"])
    species_of = st["species_of_label"]
    cfg_by = {c["label"]: c for c in st["bids_configs"]}
    rows = []
    for lab in labels:
        res = ff.StateResolver(cfg_by.get(lab, {}))
        keys = list(st["runs_by_label"].get(lab, {}))
        animals = {res.animal(tuple(str(x) for x in k)) for k in keys}
        try:
            state = ff.dataset_state(res, keys)
        except Exception:
            state = "unknown"
        rows.append((lab, species_of.get(lab), len(keys), len(animals), state))
    n_species = len({r[1] for r in rows})
    n_runs = sum(r[2] for r in rows)
    n_anim = sum(r[3] for r in rows)
    ds_txt = "; ".join(f"{ff.pretty_label(lab)} ({sp}, {nr} runs, {na} animals, {state})"
                       for lab, sp, nr, na, state in rows)
    method = cfg["PROFILE_METHOD"]
    q = float(cfg["TAIL_Q"])
    lo, hi = int(round(100 * q)), int(round(100 * (1 - q)))
    anat = st.get("anat_common") or st["common"]
    seed_found = (st.get("seed_regions") or cfg.get("_seed_regions_found")
                  or cfg.get("_seed_regions_excluded") or [])
    out_ana = bool(cfg.get("_seed_regions_excluded"))
    out_sel = bool(cfg.get("_seed_regions_excluded_selection",
                           cfg.get("_seed_regions_excluded")))
    if seed_found and out_ana and out_sel:
        seed_txt = ("; the region containing the seed (" + ", ".join(seed_found) + "), found "
                    "from the seed voxels in each atlas, was removed from every profile because "
                    "its value reflects the autocorrelation of the seed rather than coupling")
    elif seed_found and out_ana:
        seed_txt = ("; the region containing the seed (" + ", ".join(seed_found) + "), found "
                    "from the seed voxels in each atlas, was kept in the profiles scored by the "
                    "selection and removed from every subsequent analysis, because its value "
                    "reflects the autocorrelation of the seed rather than coupling")
    elif seed_found and out_sel:
        seed_txt = ("; the region containing the seed (" + ", ".join(seed_found) + ") was "
                    "removed from the profiles scored by the selection only")
    else:
        seed_txt = ""
    # the selection score: regions and mask of the selection space
    ss = st.get("selection_space") or dict(common=st["common"], anat_common=anat,
                                           region_mask=st.get("region_mask"))
    n_val = len(ss["common"])
    n_inv = (int(np.asarray(ss["region_mask"], bool).sum())
             if ss.get("region_mask") is not None else n_val)
    sel_anat = ss.get("anat_common") or anat
    inv_txt = (f"{n_inv} of the {n_val} profile values ({len(sel_anat)} regions x 2 channels)"
               if n_val != len(sel_anat) else f"{n_inv} of the {n_val} regions")
    n_kept = {lab: len(v) for lab, v in st["selections_final"].items()}
    if method == "percentiles":
        prof = (f"two values per region: the {hi}th and the {lo}th percentiles of the "
                f"region's voxel values (its most and least coupled {lo} %), left and right "
                f"hemispheres pooled by voxel count. Each of the two channels was standardised "
                f"across the common regions (mean 0, s.d. 1), so that the Pearson correlation "
                f"between two profiles is the mean of the two channel correlations and a global "
                f"offset between channels cannot create similarity")
    elif method == "fraction_posneg":
        prof = (f"two values per region: the percentage of the region's voxels in the top "
                f"{lo} % and in the bottom {lo} % of the map, the map percentiles being taken "
                f"over the voxels of the common regions (both hemispheres); rank-based, hence "
                f"invariant to any monotonic transform of the map. Each channel was standardised "
                f"across the common regions")
    else:
        prof = "the median of the region's voxel values"
    norm = {True: "z-scored across regions", "zscore": "z-scored across regions",
            "scale": "divided by their root mean square across regions (sign preserved)",
            False: "left on their raw scale", None: "left on their raw scale",
            "none": "left on their raw scale"}.get(cfg["ZSCORE"], str(cfg["ZSCORE"]))
    rounds = st["final_round"]
    chain = (f"Round 1 ranked the {unit}s of each dataset by the Pearson correlation of their "
             f"profile with the leave-one-species-out consensus (median across the other "
             f"species of the species profiles; datasets of one species first collapsed by "
             f"median, so that each species counts once)"
             + ("; round 2 re-ranked each dataset against the median of its own round-1 "
                "selection" if cfg["RUN_ROUND2_SELECTION"] else "")
             + ((f"; round 3 rebuilt" if rounds == 3 else f"; rounds 3 to {rounds} each rebuilt")
                + f" the cross-species consensus from the previous round's selected {unit}s and "
                  f"re-ranked the whole pool" if rounds >= 3 else "")
             + f". The {cfg['N_SELECT']} best-ranked {unit}s of each dataset were kept at every "
             f"round (final selection: round {rounds}).")
    versions = (f"Python {sys.version.split()[0]}, NumPy {np.__version__}, pandas "
                f"{pd.__version__}, SciPy {scipy.__version__}, Matplotlib "
                f"{matplotlib.__version__}")
    g6 = " and ".join((cfg.get("FIG6D") or {}).get("group", ("Human", "Macaque")))
    regs6 = ", ".join(ff.short_region(r) for r in (cfg.get("FIG6D") or {}).get(
        "regions", ()))
    pairs5b = ("pairs of runs from the same species but two different datasets (available for "
               "the species with two datasets) and pairs from different species; pairs from the "
               "same dataset, which also share protocol, scanner and usually the animal, were "
               "not compared"
               if cfg.get("WITHIN_ACROSS_PAIRS", "other_dataset") == "other_dataset" else
               "intra-species pairs (same species, same or another dataset) and inter-species "
               "pairs")
    test5b = ("one-sided Wilcoxon rank-sum test on the pairs, not corrected for multiple "
              "comparisons (single primary hypothesis); because every run enters many pairs, "
              "the pairs are not independent and this test is liberal, so a permutation test "
              "of the species labels across datasets (each dataset keeping its runs; 2,000 "
              "permutations; statistic, difference between the intra- and inter-species "
              "medians), which respects that dependence, is reported with it"
              if cfg.get("WITHIN_ACROSS_TEST", "wilcoxon") == "wilcoxon" else
              "permutation test of the species labels across datasets (each dataset keeping its "
              "runs; 2,000 permutations; statistic, difference between the intra- and "
              "inter-species medians), one-sided; a Wilcoxon rank-sum test on the pairs is "
              "reported for comparison")
    txt = f"""METHODS -- Cross-species retrosplenial connectivity profiles (Figures 5 and 6)
(written by fig5_pipeline from the run in {paths['root']}; numbers are those of this run)

Datasets
{len(labels)} resting-state fMRI datasets from {n_species} species ({n_runs} runs, {n_anim} animals) were analysed: {ds_txt}. Every run was preprocessed and registered to its species template with EDNiX, and a seed-based correlation map (Pearson r) of the retrosplenial seed ({cfg['SEED']}) was computed in template space.

Regional profiles
Each run's map was summarised over the regions of the EDNiX atlas at level {cfg['ATLAS_LEVEL']} into a regional profile: {prof}. Only regions present in the atlas of every dataset and sampled by every dataset's group map were kept ({len(anat)} common regions in the analyses){seed_txt}. {'Runs were screened beforehand by their distance to the most central run of their dataset (median absolute deviation threshold ' + str(cfg['RUN_QC_MAD_THRESHOLD']) + ').' if cfg['RUN_QC_ENABLE'] else 'No run was excluded before the selection.'} Profiles were {norm} before any comparison. The unit of the selection and of the group maps was the {unit}.

Consensus-based selection
{chain} The score was computed on the {inv_txt} that were the most invariant across species (top fraction {cfg['INVARIANCE_VALUE']:g}, ranked by the median absolute deviation of the species values); every common region was used everywhere else. Leave-one-species-out ensured that no dataset contributed to the criterion that selected its own runs. Selected per dataset: {', '.join(f'{ff.pretty_label(k)} {v}' for k, v in sorted(n_kept.items()))}.

Group maps and surfaces
For each dataset, group maps were computed with AFNI 3dttest++ (Cox 1996) on the Fisher-z transformed run maps of (i) every {unit}, (ii) the {unit}s of the final selection and (iii) the same number of {unit}s drawn at random among those not selected. They were rendered on each species' inflated cortical surface with Connectome Workbench (wb_command; Marcus et al. 2011), with one symmetric colour range per dataset shared by its three maps (± the {cfg.get('SURFACE_COMMON_PCT', 90):g}th percentile of |z| over their non-zero voxels).

Reliability of the profile constructions
The three constructions of the profile (median; upper and lower percentiles; tail fractions) were compared on every run before the selection, within the datasets holding several runs per animal, by their discriminability (Bridgeford et al. 2021): the probability that a run is closer (1 − r) to another run of the same animal than to a run of a different animal (0.5, chance). Its significance against chance was assessed by permuting the animal labels across runs (500 permutations), and the methods were compared with an exact paired sign-flip test across datasets. The criterion was fixed before, and independently of, the cross-species analyses.

Statistics
Figure 5b. Pearson correlations between the profiles of every pair of {unit}s were split into {pairs5b}, among every {unit} and among the selected {unit}s; the two distributions were compared with a {test5b}.
Figure 5c. The agreement between species profiles was the mean Pearson r over all pairs of species, tested against {10000:,} permutations of the region labels within each species (one-sided); this permutation ignores the spatial autocorrelation between neighbouring regions and is therefore liberal.
Figure 5d. Whether the selection reorders the datasets was tested by the Spearman correlation between the dataset-to-consensus correlations before and after selection ({10000:,} permutations, two-sided).
Figure 6a. Pearson correlations between the species-mean similarity to the leave-one-species-out consensus (or to Human) and divergence time (time-calibrated tree{': ' + os.path.basename(cfg['TREE_PATH']) if cfg.get('TREE_PATH') else ''}) or log10 brain mass (Herculano-Houzel et al. 2015; Jardim-Messeder et al. 2017) were tested exactly over all assignments of the x values to the species (two-sided), with and without Human. Profiles: {phylo_source or 'group maps'}.
Figure 6c. Within the datasets holding both arousal states in the same animals, the per-run correlation with the selection criterion was compared between anaesthetised and awake runs by the mean over animals of the per-animal difference, tested with an exact sign-flip of the animals (2^n assignments).
Figure 6d. For each of six regions ({regs6}), the regional value of every selected run (one panel per index: upper and lower percentile, median, % of voxels in the top and bottom tail of the map, % of voxels with r > 0 and r < 0) was averaged over the runs of each animal, then over the animals of each dataset, then over the datasets of each species, so that neither an animal with more runs nor a dataset with more animals weighs more. The species is the unit of inference -- runs of one species share its brain, protocol and scanner, and testing runs would be pseudoreplication (Hurlbert 1984; Lazic 2010; Aarts et al. 2014). The difference between the mean of {g6} and the mean of the other species was tested exactly over every assignment of 2 of the 7 species to the group (21 assignments; two-sided; smallest possible P = 1/21 = 0.048), one test per region, not corrected for the number of regions. Robustness: the same test with one dataset per species, for every combination of the datasets of the species that have two. Region size: the group difference of each region was compared with the linear fit, across the other common regions, of the group difference against the group difference in log10 relative region size.
Significance marks on the figures: *** P < 0.001, ** P < 0.01, * P < 0.05, n.s. not significant; exact P values, n and statistics are given in the .txt file of each panel. Permutation P values are (1 + number of null statistics at least as extreme) / (1 + number of permutations) except for exact enumerations.

Software
{versions}; AFNI (3dttest++, 3dMean); Connectome Workbench (wb_command); nibabel.

References
Aarts E, Verhage M, Veenvliet JV, Dolan CV, van der Sluis S (2014) A solution to dependency: using multilevel analysis to accommodate nested data. Nat Neurosci 17:491-496.
Bridgeford EW, et al. (2021) Eliminating accidental deviations to minimize generalization error and maximize replicability: applications in connectomics and genomics. PLoS Comput Biol 17:e1009279.
Cox RW (1996) AFNI: software for analysis and visualization of functional magnetic resonance neuroimages. Comput Biomed Res 29:162-173.
Herculano-Houzel S, Catania K, Manger PR, Kaas JH (2015) Mammalian brains are made of these: a dataset of the numbers and densities of neuronal and nonneuronal cells in the brain of glires, primates, scandentia, eulipotyphlans, afrotherians and artiodactyls, and their relationship with body mass. Brain Behav Evol 86:145-163.
Hurlbert SH (1984) Pseudoreplication and the design of ecological field experiments. Ecol Monogr 54:187-211.
Jardim-Messeder D, et al. (2017) Dogs have the most neurons, though not the largest brain: trade-off between body mass and number of neurons in the cerebral cortex of large carnivoran species. Front Neuroanat 11:118.
Lazic SE (2010) The problem of pseudoreplication in neuroscientific studies: is it affecting your analysis? BMC Neurosci 11:5.
Mantel N (1967) The detection of disease clustering and a generalized regression approach. Cancer Res 27:209-220.
Marcus DS, et al. (2011) Informatics and data mining tools and strategies for the Human Connectome Project. Front Neuroinform 5:4.
Nili H, et al. (2014) A toolbox for representational similarity analysis. PLoS Comput Biol 10:e1003553.
"""
    path = os.path.join(paths["root"], "Methods_article.txt")
    with open(path, "w", encoding="utf-8") as f:
        f.write(txt)
    print(f"  [saved] Methods_article.txt")
    return path


def _f4_per_dataset(sri, cfg, fd, sel, selections, final_round, profiles_by_label, common,
                    species_profs, species_of_label, kept_ref_by_label,
                    kept_ref_final_by_label, space_note=""):
    """F4: one figure per dataset -- its own reference, the consensus it was
    scored against, every unit (kept / not kept), the kept references. Drawn in
    the region space of the SELECTION (it shows the selection's target).
    Returns own_refs ({label: profile}); empty when no reference could be built."""
    own_refs = {}
    for label, profs in profiles_by_label.items():
        try:
            own_refs[label], _ = fa.build_reference_for_profiles(
                sri, profs, method=cfg["REFERENCE_METHOD"], k=cfg["K_REF"],
                value_col=cfg["VALUE_COL"], metric=cfg["METRIC"])
        except Exception as e:
            print(f"  [WARN] {label}: reference failed: {e}")
    if not own_refs:
        return own_refs
    unit = cfg["UNIT"]
    # cons_prof must stay on the RAW value_col scale to be comparable to
    # own_refs and the individual units. The consensus in
    # sel["consensus_by_label"] is z-scored and leave-one-out (correct for the
    # selection SCORING) -- rebuild a separate display-only consensus here: all
    # species, no LOO, no z-score.
    display_consensus, _ = fa.build_consensus_pattern(
        species_profs, common, value_col=cfg["VALUE_COL"],
        exclude_label=None, zscore=False, threshold_abs=cfg["THRESHOLD_ABS"],
        species_of_label=species_of_label, verbose=False)
    finite = np.isfinite(display_consensus)
    print(f"  [F4] display consensus (raw scale): "
          f"{int(finite.sum())}/{len(display_consensus)} region(s) finite, "
          f"range=[{np.nanmin(display_consensus):.3f}, "
          f"{np.nanmax(display_consensus):.3f}]")
    cons_prof = own_refs[next(iter(own_refs))].copy()
    cons_prof[cfg["VALUE_COL"]] = (
        cons_prof["base_name"].str.lower().str.strip()
        .map(dict(zip(common, display_consensus))))
    # Unweighted grand mean: every BIDS's own reference counted once, unlike
    # display_consensus which collapses same-species BIDS together first.
    grand_mean_vals = np.nanmedian(
        np.vstack([fa._ordered_common_values(p, common, cfg["VALUE_COL"])
                   for p in own_refs.values()]), axis=0)
    grand_mean_prof = own_refs[next(iter(own_refs))].copy()
    grand_mean_prof[cfg["VALUE_COL"]] = (
        grand_mean_prof["base_name"].str.lower().str.strip()
        .map(dict(zip(common, grand_mean_vals))))

    per = os.path.join(fd, "F4_per_dataset")
    os.makedirs(per, exist_ok=True)
    for label, profs in profiles_by_label.items():
        if label not in own_refs:
            continue
        own_vals = own_refs[label][cfg["VALUE_COL"]].values
        print(f"  [F4] {label}: own_ref range=[{np.nanmin(own_vals):.3f}, "
              f"{np.nanmax(own_vals):.3f}]  (should be same order of "
              f"magnitude as the consensus range above)")
        kept_keys = {tuple(str(x) for x in k) for k in selections.get(label, [])}
        notkept_profs = {k: v for k, v in profs.items()
                         if tuple(str(x) for x in k) not in kept_keys}
        notkept_ref = None
        try:
            if len(notkept_profs) >= 2:
                notkept_ref, _ = fa.build_reference_for_profiles(
                    sri, notkept_profs, method=cfg["REFERENCE_METHOD"],
                    k=min(cfg["K_REF"], len(notkept_profs)),
                    value_col=cfg["VALUE_COL"], metric=cfg["METRIC"])
        except Exception as e:
            print(f"  [WARN] {label}: not-kept-only reference failed: {e}")
        fa.plot_bids_ref_with_individuals(
            profs, own_refs[label], cons_prof, common,
            selected_keys=selections.get(label),
            kept_ref=kept_ref_by_label.get(label), notkept_ref=notkept_ref,
            grand_mean_ref=grand_mean_prof,
            kept_ref_round2=kept_ref_final_by_label.get(label),
            actual_selection_target=sel["consensus_by_label"].get(label),
            final_round_label=(f"Kept final round (round {final_round})"
                               if final_round >= 2 else "Kept final round"),
            target_label=(f"Round-1 selection target ("
                          f"{fa._norm_label(cfg['ZSCORE']).strip() or 'raw'}, "
                          f"{'LOO' if cfg['LEAVE_ONE_OUT'] else 'no LOO'}, "
                          f"right axis)"),
            value_col=cfg["VALUE_COL"],
            title=f"{label} (lvl {cfg['ATLAS_LEVEL']}, {unit} level)",
            output_path=os.path.join(per, f"{label}.png"))
    _write_legend(
        os.path.join(fd, "F4_per_dataset.txt"),
        "F4 -- Per-dataset reference, consensus and individual units",
        f"One figure per dataset: its own reference (blue), the consensus it "
        f"was scored against (orange dashed), and every individual {unit} "
        f"(green = kept, grey = not kept).{space_note}",
        f"Dataset reference built with {cfg['REFERENCE_METHOD']} over ALL its "
        f"{unit}s, independent of the selection. Regions ordered by consensus "
        f"value.",
        f"Kept {unit}s should track the orange line more closely than grey "
        f"ones. If they instead look flatter, the metric is rewarding low "
        f"amplitude rather than matching shape -- try "
        f"SIMILARITY='neg_euclidean' with ZSCORE=False.")
    return own_refs


def _paper_out(paths, name):
    """06_figures/Paper/<name>: print-size panels laid out by analysis/Figure_4_*.py and
    analysis/Supplementary/Supp_Fig_07-11*.py."""
    d = os.path.join(paths["figures"], "Paper")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, name)


def _paper(fn, *a, **k):
    """A print-size panel failing must not stop the pipeline."""
    try:
        return fn(*a, **k)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"  [paper panel] [WARN] {getattr(fn, '__name__', fn)} failed: {e}")
        return None


def _run_figures(sri, st):
    """Part 2: every statistic and figure, from the part-1 state `st` only.

    Two region spaces when the seed region is kept for the selection and
    removed from the analyses (state['selection_space'] set): F3, F4 and F5
    describe the selection and are drawn in its regions; every other output
    uses the analysis regions (the state's usual keys)."""
    cfg, paths, bids_configs = st["cfg"], st["paths"], st["bids_configs"]
    unit = cfg["UNIT"]
    species_of_label = st["species_of_label"]
    profiles_by_label, members_by_label = st["profiles_by_label"], st["members_by_label"]
    runs_by_label, unit_counts = st["runs_by_label"], st["unit_counts"]
    df_qc_by_label = st["df_qc_by_label"]
    all_units_index, afni_profiles = st["all_units_index"], st["afni_profiles"]
    common, n_regions = st["common"], st["n_regions"]
    species_profs, invariance, region_mask = (st["species_profs"], st["invariance"],
                                              st["region_mask"])
    sel, selections = st["sel"], st["selections"]
    all_round_selections, final_round = st["all_round_selections"], st["final_round"]
    selections_final, final_criterion = st["selections_final"], st["final_criterion"]
    map_index, map_index_final = st["map_index"], st["map_index_final"]
    group_map_profiles = st["group_map_profiles"]
    # ---- region space of the selection vs of the analyses -------------------
    split = st.get("selection_space") is not None
    ss = st["selection_space"] if split else dict(
        common=common, profiles_by_label=profiles_by_label, species_profs=species_profs,
        invariance=invariance, region_mask=region_mask, final_criterion=final_criterion)
    # species profiles of every analysis (F1, F1c, fingerprints, 5c): those of the
    # analysis space (with one space, the selection's own, as before)
    ana_species_profiles = species_profs if split else sel["species_profiles"]
    seed_txt = ", ".join(st.get("seed_regions") or cfg.get("_seed_regions_found") or [])
    sel_note = (f" Drawn in the regions of the SELECTION ({len(ss['common'])} values): the "
                f"region containing the seed ({seed_txt}) is "
                f"{'included' if not cfg.get('SEED_REGION_EXCLUDE_SELECTION') else 'excluded'}"
                f", as it was for the selection; the analyses use {len(common)} values."
                if split else "")
    # selections are unit keys; every run-level step needs the RUN keys behind
    # them (with UNIT="subject" a unit key (subject, 'ALL', 'ALL', seed) never
    # equals a run key, so kept flags of runs were all False before)
    kept_runs_final = fpp.selected_run_keys(selections_final, members_by_label)
    kept_runs_r1 = fpp.selected_run_keys(selections, members_by_label)
    print(f"\n{'#' * 70}\n  PART 2: statistics and figures\n{'#' * 70}")
    if split:
        print(f"  region spaces: selection {len(ss['common'])} values (seed region "
              f"{seed_txt} {'excluded' if cfg.get('SEED_REGION_EXCLUDE_SELECTION') else 'kept'})"
              f", analyses {len(common)} values (seed region "
              f"{'excluded' if cfg.get('SEED_REGION_EXCLUDE_ANALYSES') else 'kept'}); F3, F4, "
              f"F5 in the selection's regions, everything else in the analysis regions")

    # ---- figures computed in part 1, drawn here ------------------------------
    if cfg["MAKE_FIGURES"]:
        if afni_profiles:
            print(f"\n{'=' * 70}\n  F0 figures (per BIDS, then per species)\n{'=' * 70}")
            step1b_species_profiles_figure(sri, cfg, paths, afni_profiles,
                                           common, region_mask, invariance)
        if st["overlap_df"] is not None and not st["overlap_df"].empty:
            fa.plot_reference_method_overlap(
                st["overlap_df"], output_path=os.path.join(
                    paths["figures"], "F1d_reference_method_effect.png"))
        if st["conv_df"] is not None and not st["conv_df"].empty:
            fa.plot_iterative_selection_convergence(
                st["conv_df"], output_path=os.path.join(
                    paths["figures"], "F1e_selection_convergence.png"))
            if cfg.get("RUN_PAPER_FIGURES", True):
                _paper(fsp.f1e_convergence, st["conv_df"], _paper_out(paths, "F1e_convergence"))
        rel = st.get("reliability")
        if rel and not rel["rel"].empty:
            png = os.path.join(paths["figures"], "F22_profile_reliability_by_method.png")
            ff.plot_profile_reliability(rel["rel"], species_of_label, png,
                                        current_method=cfg["PROFILE_METHOD"])
            _write_profile_reliability_legend(png, rel["rel"], rel["tests"], cfg)

    # ---- phylogeny source (F6, F7, F11-F13, STEP 7) -----------------------
    # Resolved ONCE here so every phylogenetic output uses the same maps.
    phylo_profiles = afni_profiles
    phylo_index = all_units_index
    phylo_source = f"all-{unit}s group maps (no selection)"
    # Every dataset that has runs must be represented; compare against that,
    # not against afni_profiles (which can itself be incomplete).
    expected = set(profiles_by_label)
    missing_all_units = sorted(expected - set(afni_profiles))
    if missing_all_units:
        print(f"  [PHYLO] [WARN] no all-{unit}s profile for {missing_all_units}")
    if cfg["PHYLOGENY_USE_FINAL_ROUND"] and final_round >= 2:
        kf = (map_index_final[(map_index_final.group == "kept_final") & map_index_final.found]
              if not map_index_final.empty else pd.DataFrame(columns=["label"]))
        missing = sorted(expected - set(kf["label"]))
        if missing:
            print(f"  [PHYLO] [WARN] PHYLOGENY_USE_FINAL_ROUND=True but no round-"
                  f"{final_round} kept map for {missing} -- keeping all-{unit}s maps "
                  f"(mixing sources across species would bias F6)")
        else:
            cand = {lab: v for lab, v in group_map_profiles.get("kept_final", {}).items()
                    if lab in expected}
            not_profiled = sorted(expected - set(cand))
            if not_profiled:
                print(f"  [PHYLO] [WARN] could not profile round-{final_round} kept maps "
                      f"for {not_profiled} -- keeping all-{unit}s maps")
            else:
                phylo_profiles, phylo_index = cand, kf
                phylo_source = f"round-{final_round} KEPT group maps (selection applied)"
    phylo_kind = "afni_kept_final" if phylo_profiles is not afni_profiles else "afni_all"
    if not phylo_profiles:
        # No AFNI group map at all (RUN_ALL_UNITS_MAPS=False, e.g. the matrix
        # analysis of fig5_matrices): fall back to the median over every unit
        # of each dataset -- still no selection, only a different estimator.
        phylo_profiles = {
            lab: (species_of_label.get(lab), sri.compute_group_median_profile(profs, set(profs)))
            for lab, profs in profiles_by_label.items()}
        phylo_source = f"median of all {unit} profiles (no AFNI group map, no selection)"
        phylo_kind = "python_all"
    if set(phylo_profiles) != expected:
        print(f"  [PHYLO] [WARN] phylogenetic analyses cover {len(phylo_profiles)}/"
              f"{len(expected)} datasets -- missing {sorted(expected - set(phylo_profiles))}")
    print(f"\n  [PHYLO] phylogenetic analyses use: {phylo_source}")

    # ---- STEP 6: figures + legends --------------------------------------
    fingerprints = {}
    kept_ref_by_label, kept_ref_final_by_label = {}, {}
    if cfg["MAKE_FIGURES"]:
        print(f"\n{'=' * 70}\n  STEP 6: figures\n{'=' * 70}")
        fd = paths["figures"]

        fa.plot_consensus_by_species(
            ana_species_profiles, common, species_of_label,
            value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"],
            title=f"Consensus by species -- all runs (lvl {cfg['ATLAS_LEVEL']}, {unit} level)",
            output_path=os.path.join(fd, "F1_consensus_by_species.png"))
        _write_legend(
            os.path.join(fd, "F1_consensus_by_species.txt"),
            "F1 -- Consensus pattern by species (all runs)",
            f"One line per species, over the {len(common)} regions common to every "
            f"dataset. Black = the consensus (median across species). Built from "
            f"EVERY run (kept and not-kept alike) -- source=all_runs, no selection.",
            f"Built from the {unit}-level profiles (runs averaged within {unit} by "
            f"median). Datasets of the same species were collapsed by median "
            f"first, so a species with several datasets counts once. Profiles were "
            f"{fa._norm_label(cfg['ZSCORE']).strip() or 'left unscaled'} before combining.",
            "Regions are ordered by the consensus value. Species lines tracking the "
            "black line indicates a shared profile shape; a species departing from "
            "it is atypical for that region. F0 shows the same thing computed from "
            "the AFNI group maps instead -- the two should agree. See F1a/F1b for "
            "the same plot built from round-1 / final-round KEPT runs only instead of "
            "every run.")

        # Kept-round-1 and kept-final per-BIDS references (analysis space), reused
        # by F1a/F1b and the fingerprints; F4 has its own, in the selection space.
        kept_ref_by_label, kept_ref_final_by_label = _kept_references(
            sri, cfg, profiles_by_label, selections, selections_final, final_round)

        if kept_ref_by_label:
            fa.plot_consensus_by_species(
                kept_ref_by_label, common, species_of_label,
                value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"],
                title=f"Consensus by species -- round-1 kept only "
                      f"(lvl {cfg['ATLAS_LEVEL']}, {unit} level)",
                output_path=os.path.join(fd, "F1a_consensus_by_species_kept1.png"))
            _write_legend(
                os.path.join(fd, "F1a_consensus_by_species_kept1.txt"),
                "F1a -- Consensus pattern by species (round-1 kept only)",
                "Same as F1, but each BIDS's line is now its round-1 KEPT "
                "reference (cross-species-anchored selection) instead of the "
                "all-runs median.",
                f"{cfg['REFERENCE_METHOD']} reference over each BIDS's round-1 "
                f"kept runs only (k={cfg['K_REF']}, metric={cfg['METRIC']}).",
                "Compare against F1: if a species' line moves a lot from F1 to "
                "F1a, round-1 selection changed what that species looks like "
                "here -- either cleaning it up (closer to consensus) or biasing "
                "it toward the consensus by construction, which F6's use of the "
                "all-units map (not this kept map) by default is meant to avoid.")

        if kept_ref_final_by_label:
            fa.plot_consensus_by_species(
                kept_ref_final_by_label, common, species_of_label,
                value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"],
                title=f"Consensus by species -- final round (round {final_round}) kept only "
                      f"(lvl {cfg['ATLAS_LEVEL']}, {unit} level)",
                output_path=os.path.join(fd, "F1b_consensus_by_species_kept_final.png"))
            _write_legend(
                os.path.join(fd, "F1b_consensus_by_species_kept_final.txt"),
                f"F1b -- Consensus pattern by species (final round = round {final_round})",
                f"Same as F1a, but each BIDS's line is its KEPT reference from the LAST "
                f"selection round (round {final_round}; chain: round 1 cross-species"
                f"{', round 2 species-anchored' if 2 in all_round_selections else ''}"
                f"{f', rounds 3-{final_round} cross-species from previous kept' if final_round >= 3 else ''}).",
                f"{cfg['REFERENCE_METHOD']} reference over each BIDS's round-{final_round} "
                f"kept runs only.",
                "Compare against F1a: if F1b and F1a are nearly identical, the refinement "
                "rounds barely changed that species. F1e shows round-to-round stability.")

        fa.plot_consensus_by_bids(
            ana_species_profiles, common, species_of_label,
            value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"],
            title=f"Consensus by BIDS -- all runs (lvl {cfg['ATLAS_LEVEL']}, {unit} level)",
            output_path=os.path.join(fd, "F1c_consensus_by_bids.png"))
        if cfg.get("RUN_PAPER_FIGURES", True):
            _paper(fsp.f1c_consensus_by_dataset, ana_species_profiles, common,
                   species_of_label, cfg["VALUE_COL"], cfg["ZSCORE"],
                   _paper_out(paths, "F1c_consensus_by_dataset"))
        _write_legend(
            os.path.join(fd, "F1c_consensus_by_bids.txt"),
            "F1c -- Consensus pattern by BIDS (all runs)",
            f"Companion to F1: one line per BIDS (10 here) instead of one per "
            f"species (7 in F1) -- same black consensus line, unchanged.",
            f"Same {unit}-level profiles and z-scoring as F1, just not collapsed "
            f"by species before plotting.",
            "Lets you see whether the 2-BIDS species (Human, Macaque, Marmoset) "
            "actually agree with each other before F1 collapses them into one "
            "line -- two BIDS of the same species tracking each other closely "
            "but diverging from the consensus points at a real species effect; "
            "two BIDS of the same species diverging from EACH OTHER points at a "
            "dataset/protocol effect instead.")

        pairs_def = cfg.get("WITHIN_ACROSS_PAIRS", "other_dataset")
        fa.plot_intra_vs_inter_species_distance(
            profiles_by_label, common, species_of_label,
            value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"],
            title=f"Within- vs across-species similarity (lvl {cfg['ATLAS_LEVEL']})",
            output_path=os.path.join(fd, "F2_within_vs_across_species.png"),
            intra=pairs_def)
        compared = ("'same species, other dataset' vs 'different species' (same-dataset "
                    "pairs shown, not compared)" if pairs_def == "other_dataset" else
                    "every same-species pair (same dataset included) vs 'different species'")
        _write_legend(
            os.path.join(fd, "F2_within_vs_across_species.txt"),
            "F2 -- Within- vs across-species similarity",
            f"Distribution of {unit}-to-{unit} profile correlation for three kinds of "
            f"pairs: same dataset; same species but another dataset; different species.",
            f"All {unit} pairs (subsampled above 20000 pairs). Classes compared "
            f"(WITHIN_ACROSS_PAIRS={pairs_def!r}): {compared}. AUC is descriptive. p comes "
            f"from permuting species labels ACROSS DATASETS (each dataset keeps its {unit}s; "
            f"pair correlations fixed): pairs share {unit}s and are not independent, so a "
            f"Mann-Whitney p on pairs would be anti-conservative by orders of magnitude."
            + (f" With UNIT='run', same-dataset pairs include repeated scans of the same "
               f"animal and the same protocol, so they overestimate species similarity."
               if unit == "run" else ""),
            "The informative comparison is 'same species, other dataset' vs 'different "
            "species': it separates species from dataset/protocol and animal identity. "
            "If it is not clearly above the different-species class, species is not the "
            "dominant structure and phylogenetic interpretation is unsupported. Only "
            "species with several datasets contribute to that middle class, so the "
            "permutation test has few distinct configurations and limited power.")

        selected_by_label = {}
        for label, profs in profiles_by_label.items():
            kept_keys = {tuple(str(x) for x in k) for k in selections_final.get(label, [])}
            selected_by_label[label] = {k: v for k, v in profs.items()
                                        if tuple(str(x) for x in k) in kept_keys}
        fa.plot_within_vs_across_all_vs_selected(
            profiles_by_label, selected_by_label, common, species_of_label,
            value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"],
            title=f"Within- vs across-species similarity: all vs final-round selection "
                  f"(lvl {cfg['ATLAS_LEVEL']})",
            selected_name=f"Selected, round {final_round}",
            output_path=os.path.join(fd, "F2b_within_vs_across_all_vs_selected.png"),
            intra=pairs_def)
        _write_legend(
            os.path.join(fd, "F2b_within_vs_across_all_vs_selected.txt"),
            "F2b -- Within- vs across-species similarity, all runs vs selected",
            f"Same analysis as F2: left, every {unit}; right, the {unit}s kept by the "
            f"final selection round (round {final_round}). Paper panel 5b draws the two "
            f"compared classes of this figure with the same pairs, bins and statistics.",
            f"Same pair classes, AUC and dataset-level permutation test as F2 (compared: "
            f"{compared}), computed independently on each pool.",
            "The left panel is the unbiased estimate. Selection moves the right panel "
            "in both directions: rounds scored against a cross-species consensus pull "
            "runs toward a shared pattern and should REDUCE species separation, while a "
            "species-anchored round (round 2) pulls each dataset toward its own mean "
            "and can INFLATE it. A separation seen only on the right is therefore not "
            "evidence of species structure.")

        fa.plot_score_distributions(
            sel["rankings"], selections,
            title=f"Consensus score per dataset (lvl {cfg['ATLAS_LEVEL']}, {unit} level)",
            output_path=os.path.join(fd, "F3_score_distributions.png"))
        if cfg.get("RUN_PAPER_FIGURES", True):
            _paper(fsp.f3_score_distributions, sel["rankings"], selections,
                   _paper_out(paths, "F3_score_distributions"))
        _write_legend(
            os.path.join(fd, "F3_score_distributions.txt"),
            "F3 -- Selection score distributions",
            f"Per dataset, the distribution of every {unit}'s similarity to the "
            f"consensus, with the selection cutoff marked in orange.{sel_note}",
            f"Similarity metric: {cfg['SIMILARITY']}. Cutoff is the score of the "
            f"{cfg['N_SELECT']}th best {unit}.",
            f"If the kept {unit}s form a separate upper mode, the selection captures "
            f"a real subgroup. If the distribution is smooth across the cutoff, the "
            f"selection is a cut through a continuum and should be described as such. "
            f"Note the x-axis is now much shorter than at run level -- each dataset "
            f"has far fewer {unit}s than runs, so a single {unit} moves the cutoff "
            f"more.")

        # F4 shows the selection's own target: drawn in the selection's regions
        if split:
            f4_kept, f4_kept_final = _kept_references(
                sri, cfg, ss["profiles_by_label"], selections, selections_final, final_round)
        else:
            f4_kept, f4_kept_final = kept_ref_by_label, kept_ref_final_by_label
        own_refs = _f4_per_dataset(
            sri, cfg, fd, sel, selections, final_round, ss["profiles_by_label"], ss["common"],
            ss["species_profs"], species_of_label, f4_kept, f4_kept_final,
            space_note=sel_note)
        if own_refs:
            sign_df = fa.diagnose_consensus_sign(
                species_profs, common, value_col=cfg["VALUE_COL"],
                species_of_label=species_of_label)
            sign_df.to_excel(os.path.join(
                paths["group_maps"], "consensus_sign_diagnostic.xlsx"), index=False)

        fa.plot_region_invariance(
            ss["invariance"], ss["common"], mask=ss["region_mask"],
            title=f"Region invariance across species (lvl {cfg['ATLAS_LEVEL']})",
            output_path=os.path.join(fd, "F5_region_invariance.png"))
        if cfg.get("RUN_PAPER_FIGURES", True):
            _paper(fsp.f5_region_invariance, ss["invariance"], ss["common"], ss["region_mask"],
                   _paper_out(paths, "F5_region_invariance"))
        _write_legend(
            os.path.join(fd, "F5_region_invariance.txt"),
            "F5 -- Cross-species region invariance",
            "Per region, how much species AGREE on its value. Green bars are the "
            f"regions the selection was restricted to.{sel_note}",
            f"Invariance = 1 - rescaled median absolute deviation of the per-species "
            f"values, computed on the {unit}-level profiles. Selection used "
            f"mode={cfg['INVARIANCE_MODE']}, value={cfg['INVARIANCE_VALUE']}.",
            "Grey regions vary across species, so matching them would mean matching "
            "cross-species noise. If very few regions are green, the conserved core "
            "is small and any cross-species claim rests on little.")

        # --- phylogenetic regression (F6/F7) --------------------------------
        if phylo_profiles:
            print(f"  [F6] using {phylo_source} ({len(phylo_profiles)} dataset(s)), "
                  f"leave-one-SPECIES-out reference per point")
            ddf_bids = fa.compute_loo_species_distances(
                phylo_profiles, common, value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"])
            if not ddf_bids.empty:
                ddf_bids["source"] = phylo_source
            if len(ddf_bids) >= 3:
                # One row per SPECIES for the fit: two BIDS of one species
                # are not independent observations.
                ddf = (ddf_bids.groupby("species", as_index=False)
                       .agg(distance_to_human=("distance_to_human", "mean"),
                            pearson_r=("pearson_r", "mean"),
                            n_bids=("label", "count")))
                ddf["source"] = "all_units"

                # Intra- vs inter-species contrast
                intra = (ddf_bids.groupby("species")["distance_to_human"]
                         .agg(lambda s: float(np.abs(s.max() - s.min())) if len(s) > 1 else np.nan)
                         .dropna())
                d = ddf["distance_to_human"].values
                inter = float(np.nanmean(np.abs(d[:, None] - d[None, :])[np.triu_indices(len(d), 1)])) \
                    if len(d) > 1 else np.nan
                if len(intra):
                    print(f"  [F6] intra-species spread (|BIDS1-BIDS2|): "
                          f"mean={intra.mean():.4f} over {len(intra)} species "
                          f"{dict(intra.round(4))}")
                    print(f"  [F6] inter-species mean pairwise distance: {inter:.4f}")
                    print(f"  [F6] ratio inter/intra = {inter / intra.mean():.2f} "
                          f"(>1 means datasets of a species agree more with each "
                          f"other than species do with each other)")
                pd.DataFrame(dict(metric=["intra_species_spread_mean",
                                          "inter_species_pairwise_mean",
                                          "ratio_inter_over_intra"],
                                  value=[float(intra.mean()) if len(intra) else np.nan,
                                         inter,
                                         float(inter / intra.mean()) if len(intra) else np.nan])
                             ).to_excel(os.path.join(paths["group_maps"],
                                                     "intra_vs_inter_species.xlsx"), index=False)

                with pd.ExcelWriter(os.path.join(paths["group_maps"],
                                                 "distance_to_consensus.xlsx")) as xw:
                    # "distance_to_human" is the historical column name: it is the
                    # distance to the leave-one-species-out consensus.
                    for df_ in (ddf_bids, ddf):
                        df_["divergence_mya_used"] = df_["species"].map(fa.DIVERGENCE_TIME_MYA)
                        df_["divergence_source"] = cfg["TREE_PATH"] or "hand-entered (approximate)"
                    ddf_bids.to_excel(xw, sheet_name="per_bids", index=False)
                    ddf.to_excel(xw, sheet_name="per_species", index=False)
                print(f"  [F6] distance_to_consensus.xlsx: {len(ddf_bids)} BIDS row(s), "
                      f"{len(ddf)} species row(s) -- fit on species rows only")
                tr = fa.test_phylogenetic_signal(ddf, exclude_human=False)
                fa.plot_distance_vs_divergence(
                    ddf_bids, test_result=tr, exclude_human=False,
                    title=f"Distance to consensus vs divergence (lvl {cfg['ATLAS_LEVEL']})",
                    output_path=os.path.join(fd, "F6_distance_vs_divergence.png"))
                perm = fa.permutation_test_phylogeny(
                    tr["data"], fa.DIVERGENCE_TIME_MYA, n_perm=cfg["N_PERM_PHYLOGENY"])
                if perm is not None:
                    fa.plot_phylogeny_permutation(
                        perm, title=f"Phylogeny permutation (lvl {cfg['ATLAS_LEVEL']})",
                        output_path=os.path.join(fd, "F7_phylogeny_permutation.png"))
                _write_legend(
                    os.path.join(fd, "F6_distance_vs_divergence.txt"),
                    "F6/F7 -- Distance to consensus vs phylogenetic divergence",
                    f"Distance (1 - r) of each dataset's group map ({phylo_source}) to "
                    f"the leave-one-species-out median across species, against "
                    f"divergence time from human. "
                    "One point per BIDS, so the two datasets of a species appear "
                    "separately; the regression is fitted on the per-species means "
                    "only.",
                    (f"Group maps from {unit}-level 3dttest on EVERY {unit}, with no "
                     f"selection -- deliberately not the kept maps, whose selection "
                     f"pressure varies with pool size across datasets and would bias "
                     f"the slope. " if phylo_kind == "afni_all" else
                     f"Median over EVERY {unit} profile of each dataset, no selection "
                     f"(no AFNI group map in this run). " if phylo_kind == "python_all" else
                     f"SENSITIVITY ANALYSIS: round-{final_round} KEPT maps "
                     f"(PHYLOGENY_USE_FINAL_ROUND=True). Selection pressure varies with "
                     f"pool size, which is confounded with clade (winner's curse); "
                     f"{'LEAVE_ONE_OUT=False also made the selection circular. ' if not cfg['LEAVE_ONE_OUT'] else ''}"
                     f"Compare with the all-units run before interpreting. ")
                    + f"Distance computed on the {len(common)} common "
                    f"regions. Fit on species means (per-BIDS rows are not "
                    f"independent). F7 shuffles divergence times for an exact null.",
                    "A positive slope means more distantly related species differ more "
                    "from the consensus. Check intra_vs_inter_species.xlsx first: if "
                    "the two BIDS of a species sit closer to each other than species "
                    "do to one another (ratio > 1), the cross-species differences are "
                    "not protocol artefacts. With ~7 species one point moves the fit "
                    "a lot -- trust F7's permutation p over the parametric one.")
            else:
                print(f"  [SKIP] phylogeny figures: only {len(ddf_bids)} usable "
                      f"group-map profile(s), need >= 3")

        # ---- F6ab, fingerprints, F15-F18 (fig5_fingerprints) ------------
        if cfg["RUN_FINGERPRINTS"] and phylo_profiles:
            print(f"\n{'=' * 70}\n  Fingerprints, two-reference F6, proximity, "
                  f"before/after, arousal state\n{'=' * 70}")
            # Before/after needs the SAME estimator on both sides: AFNI
            # all-units vs AFNI kept-final maps when every dataset has both,
            # otherwise the median of all vs kept unit profiles.
            kept_final_keys = {lab: {tuple(str(x) for x in k)
                                     for k in selections_final.get(lab, [])}
                               for lab in profiles_by_label}
            py_all = {lab: (species_of_label.get(lab),
                            sri.compute_group_median_profile(profs, set(profs)))
                      for lab, profs in profiles_by_label.items()}
            py_sel = {}
            for lab, profs in profiles_by_label.items():
                keys = {k for k in profs if tuple(str(x) for x in k) in kept_final_keys[lab]}
                if keys:
                    py_sel[lab] = (species_of_label.get(lab),
                                   sri.compute_group_median_profile(profs, keys))
            ba_all, ba_sel = py_all, py_sel
            ba_source = (f"median of {unit} profiles: all vs round-{final_round} kept")
            if (afni_profiles and set(afni_profiles) == expected and final_round >= 2
                    and not map_index_final.empty):
                kf = map_index_final[(map_index_final.group == "kept_final")
                                     & map_index_final.found]
                if set(kf["label"]) >= expected:
                    afni_sel = {lab: v for lab, v in
                                group_map_profiles.get("kept_final", {}).items()
                                if lab in expected}
                    if set(afni_sel) == expected:
                        ba_all, ba_sel = afni_profiles, afni_sel
                        ba_source = (f"AFNI group maps: all {unit}s vs round-{final_round} "
                                     f"kept")
            fingerprints = ff.make_fingerprint_figures(
                group_all=phylo_profiles, group_all_source=phylo_source,
                ba_all=ba_all, ba_sel=ba_sel, ba_source=ba_source,
                f1_all=ana_species_profiles,
                f1_sel=kept_ref_final_by_label or kept_ref_by_label,
                common=common, region_mask=region_mask, species_of_label=species_of_label,
                bids_configs=bids_configs, runs_by_label=runs_by_label,
                selections_final=kept_runs_final, final_round=final_round, cfg=cfg,
                fig_dir=fd, table_dir=paths["group_maps"],
                divergence=dict(fa.DIVERGENCE_TIME_MYA),
                ref_species=cfg["REFERENCE_SPECIES"], n_perm=cfg["N_PERM_PHYLOGENY"])

            if cfg["RUN_PAPER_FIGURES"]:
                print(f"\n{'=' * 70}\n  Paper panels: Figures 5 and 6\n{'=' * 70}")
                try:
                    fingerprints["paper"] = fpp.make_paper_figures(
                        sri=sri, fu=fu, cfg=cfg, paths=paths, bids_configs=bids_configs,
                        species_of_label=species_of_label, common=common,
                        region_mask=region_mask, final_round=final_round,
                        profiles_by_label=profiles_by_label, runs_by_label=runs_by_label,
                        members_by_label=members_by_label, selections_final=selections_final,
                        rankings=sel["rankings"], selections_r1=selections,
                        f1_all=ana_species_profiles,
                        f1_sel=kept_ref_final_by_label or kept_ref_by_label,
                        group_all=phylo_profiles, group_all_source=phylo_source,
                        fingerprints=fingerprints, all_units_index=all_units_index,
                        map_index_final=map_index_final, phylo_index=phylo_index,
                        df_qc_by_label=df_qc_by_label,
                        divergence=dict(fa.DIVERGENCE_TIME_MYA),
                        surface_pics=st["surface_pics"], surface_ranges=st["surface_ranges"],
                        run_index=st["run_index"], anat_common=st["anat_common"])
                except Exception as e:
                    import traceback
                    traceback.print_exc()
                    print(f"  [WARN] paper panels failed: {e}")

        # ---- F11-F13: comparative phylogenetic tests on the full tree ----
        if cfg["RUN_PHYLO_COMPARATIVE"] and phylo_profiles:
            print(f"\n{'=' * 70}\n  Phylogenetic comparative tests (Mantel / K / PGLS)\n{'=' * 70}")
            fp.check_tree_against_divergence_times(fa.DIVERGENCE_TIME_MYA)
            mres = fp.run_mantel(phylo_profiles, common, value_col=cfg["VALUE_COL"],
                                 zscore=cfg["ZSCORE"])
            if mres is not None:
                with pd.ExcelWriter(os.path.join(paths["group_maps"],
                                                 "phylo_mantel.xlsx")) as xw:
                    pd.DataFrame(mres["D_func"], index=mres["species"],
                                 columns=mres["species"]).to_excel(xw, sheet_name="D_functional")
                    pd.DataFrame(mres["D_phylo"], index=mres["species"],
                                 columns=mres["species"]).to_excel(xw, sheet_name="D_mrca_Mya")
                    mres["bids_pairs"].to_excel(xw, sheet_name="bids_pairs", index=False)
                    pd.DataFrame([{k: v for k, v in m.items() if k != "null"}
                                  for m in (mres["mantel_pearson"], mres["mantel_spearman"])]
                                 ).to_excel(xw, sheet_name="mantel", index=False)
                fp.plot_mantel(mres, title=f"Pairwise functional vs phylogenetic distance "
                                           f"(lvl {cfg['ATLAS_LEVEL']})",
                               output_path=os.path.join(fd, "F11_mantel.png"))
                _write_legend(
                    os.path.join(fd, "F11_mantel.txt"),
                    "F11 -- Mantel test: functional vs phylogenetic distance",
                    "Top: species x species functional distance (mean 1 - r over "
                    f"every cross-species BIDS pair of {phylo_source} profiles) "
                    "and MRCA age from the time-calibrated tree. Bottom-left: every "
                    "species pair (blue, tested) and every BIDS pair (grey); "
                    "within-species BIDS pairs in orange at 0 Mya (reported, not "
                    "tested). Bottom-right: exact permutation null.",
                    f"Mantel correlation between upper triangles; null = all "
                    f"{mres['mantel_pearson']['n_perm']} permutations of species "
                    f"labels (exact). Pearson and Spearman. {len(common)} common "
                    f"regions, profiles z-scored within dataset.",
                    "Primary phylogenetic test: uses the whole tree, not only the "
                    "distance to human, so Rat/Mouse are no longer pseudo-replicates. "
                    "p(>) is the pre-specified directional test (more distant -> "
                    "more different). Within-species BIDS distances as large as "
                    "cross-species ones mean dataset effects rival species effects. "
                    "Mantel tests can be liberal on phylogenetic data (Harmon & "
                    "Glor 2010) -- read with F12/F13.")

            kdf = fp.regionwise_phylo_signal(phylo_profiles, common,
                                             value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"])
            if not kdf.empty:
                kdf.to_excel(os.path.join(paths["group_maps"], "phylo_blomberg_K.xlsx"),
                             index=False)
                fp.plot_regionwise_k(kdf, output_path=os.path.join(fd, "F12_blomberg_K.png"))
                _write_legend(
                    os.path.join(fd, "F12_blomberg_K.txt"),
                    "F12 -- Phylogenetic signal per region (Blomberg's K)",
                    "Per common region, Blomberg's K of the species-level value "
                    "(median over that species' BIDS of the z-scored group-map "
                    "profile). Green = q_FDR < 0.05.",
                    "K from the Brownian covariance of the time-calibrated tree; "
                    "exact permutation p (all tip shuffles); Benjamini-Hochberg "
                    "across regions.",
                    "K = 1: as much signal as Brownian motion; K << 1: less. With 7 "
                    "species, power is low -- a non-significant region is not "
                    "evidence of absence. Useful to see WHICH regions drive F11.")

            ddf_pg = fa.compute_loo_species_distances(
                phylo_profiles, common, value_col=cfg["VALUE_COL"],
                zscore=cfg["ZSCORE"], verbose=False)
            if len(ddf_pg) >= 4:
                sp_df = (ddf_pg.groupby("species", as_index=False)["distance_to_human"]
                         .mean().rename(columns={"distance_to_human": "distance_to_consensus"}))
                n_runs = {sp: sum(len(runs_by_label[l]) for l in runs_by_label
                                  if species_of_label.get(l) == sp)
                          for sp in sp_df["species"]}
                sp_df["log_n_runs"] = sp_df["species"].map(
                    lambda s: np.log(n_runs[s]) if n_runs.get(s) else np.nan)
                covs = ["log_n_runs"]
                qc_cov = fp.species_covariates_from_qc(
                    df_qc_by_label, species_of_label, list(cfg["PGLS_QC_COVARIATES"]))
                if not qc_cov.empty:
                    sp_df = sp_df.merge(qc_cov, on="species", how="left")
                    covs += [c for c in cfg["PGLS_QC_COVARIATES"] if c in sp_df.columns]
                for cname, values in (cfg["SPECIES_COVARIATES"] or {}).items():
                    sp_df[cname] = sp_df["species"].map(values)
                    covs.append(cname)
                pg = fp.pgls_battery(sp_df, "distance_to_consensus", covs)
                if not pg.empty:
                    with pd.ExcelWriter(os.path.join(paths["group_maps"],
                                                     "phylo_pgls.xlsx")) as xw:
                        pg.to_excel(xw, sheet_name="pgls", index=False)
                        sp_df.to_excel(xw, sheet_name="species_data", index=False)
                    print(f"  [saved] phylo_pgls.xlsx")

        print(f"  [saved] figures + legends -> 06_figures/")

    # ---- STEP 7: F6 at several atlas levels ------------------------------
    multilevel_results = {}
    if cfg["RUN_MULTILEVEL_F6"]:
        print(f"\n{'=' * 70}\n  STEP 7: F6 at levels {cfg['MULTILEVEL_F6_LEVELS']}\n{'=' * 70}")
        print(f"  [STEP 7] source: {phylo_source}")
        multilevel_results = step6_multilevel_f6(
            sri, fa, fu, cfg, paths, bids_configs, phylo_index, members_by_label,
            levels=cfg["MULTILEVEL_F6_LEVELS"])

    # ---- STEP 8: per-run distance to the kept reference -------------------
    dist_df = pd.DataFrame()
    corr_df = pd.DataFrame()
    if cfg["RUN_DISTANCE_TO_REF"]:
        # Dispersion is measured around the dataset's reference from the
        # FINAL selection round (its kept map), whatever
        # PHYLOGENY_USE_FINAL_ROUND says.
        use_final = (final_round >= 2 and not map_index_final.empty
                     and cfg["DISTANCE_REF_GROUP"] == "kept")
        ref_group = "kept_final" if use_final else cfg["DISTANCE_REF_GROUP"]
        ref_selections = kept_runs_final if use_final else kept_runs_r1
        ref_name = (f"round-{final_round} kept" if use_final
                    else f"round-1 {cfg['DISTANCE_REF_GROUP']}")
        print(f"\n{'=' * 70}\n  STEP 8: per-run distance to each BIDS' "
              f"{ref_name} group-map reference\n{'=' * 70}")
        # profiles of these maps were read in part 1 (07_state)
        ref_profiles = dict(group_map_profiles.get(ref_group, {}))
        if not ref_profiles:
            print(f"  [STEP 8] [WARN] no '{ref_group}' group-map profile in the saved state "
                  f"(groups: {sorted(group_map_profiles)})")
        dist_df = fa.compute_run_distance_to_reference(
            runs_by_label, ref_profiles, common, value_col=cfg["VALUE_COL"],
            kept_keys_by_label=ref_selections)
        if not dist_df.empty:
            dist_df.to_excel(os.path.join(paths["group_maps"], "run_distance_to_reference.xlsx"),
                            index=False)
            fa.plot_run_distance_distribution(
                dist_df, output_path=os.path.join(
                    paths["figures"], "F8_run_distance_to_reference.png"),
                title=f"Per-run correlation to each BIDS' {ref_name} group-map reference "
                      f"(x = round-{final_round if use_final else 1} kept)")
            fa.plot_run_distance_std_by_bids(
                dist_df, output_path=os.path.join(
                    paths["figures"], "F9_run_distance_std_by_bids.png"))
            if cfg.get("RUN_PAPER_FIGURES", True):
                _paper(fsp.f9_run_distance_sd, dist_df, _paper_out(paths, "F9_run_distance_sd"))

            # ---- STEP 9: correlate against external QC metrics ------------
            if df_qc_by_label:
                print(f"\n{'=' * 70}\n  STEP 9: QC metrics vs distance to reference\n{'=' * 70}")
                corr_df = fa.correlate_qc_with_distance(
                    dist_df, df_qc_by_label, qc_metrics=cfg["QC_METRICS_FOR_CORRELATION"])
                if not corr_df.empty:
                    corr_df.to_excel(os.path.join(
                        paths["group_maps"], "qc_vs_distance_correlation_per_bids.xlsx"),
                        index=False)
                    meta_df = fa.meta_analyze_qc_correlations(
                        corr_df, min_bids=cfg["QC_META_MIN_BIDS"], alpha=cfg["QC_META_ALPHA"])
                    if not meta_df.empty:
                        meta_df.to_excel(os.path.join(
                            paths["group_maps"], "qc_vs_distance_meta_analysis.xlsx"),
                            index=False)
                        fa.plot_qc_distance_correlation(
                            meta_df, alpha=cfg["QC_META_ALPHA"], output_path=os.path.join(
                                paths["figures"], "F10_qc_vs_distance_correlation.png"))
                        if cfg.get("RUN_PAPER_FIGURES", True):
                            _paper(fsp.f10_qc_vs_distance, meta_df,
                                   _paper_out(paths, "F10_qc_vs_distance"),
                                   alpha=cfg["QC_META_ALPHA"])
            else:
                print("  [SKIP] STEP 9: no df_qc_by_label given to run() -- "
                      "pass QC data to correlate it against the per-run distances")

    if cfg["RUN_DISTANCE_TO_REF"]:
        # ---- F8b: distance to the FINAL round's selection criterion ----
        # The profile the final round scored each BIDS's runs against (round 1:
        # LOO/non-LOO consensus; round 2: the BIDS's own round-1 kept mean;
        # round >= 3: the consensus rebuilt from the previous round's kept
        # means). With two region spaces it is the same criterion rebuilt in the
        # analysis regions (_criterion_in_space, part 1).
        crit_note = None
        if split:
            crit_note = (f"The criterion is rebuilt in the analysed regions ({len(common)} "
                         f"values; seed region {seed_txt} removed), with the same construction "
                         f"and the same kept runs at every round as the selection, which "
                         f"scored {len(ss['common'])} values (seed region included); the r "
                         f"is therefore close to, not identical with, the selection score.")
        ref_profiles_criterion = {
            label: (species_of_label.get(label), ref)
            for label, ref in final_criterion.items() if ref is not None}
        dist_df_b = fa.compute_run_distance_to_reference(
            runs_by_label, ref_profiles_criterion, common,
            value_col=cfg["VALUE_COL"], kept_keys_by_label=kept_runs_final)
        if not dist_df_b.empty:
            dist_df_b.to_excel(os.path.join(
                paths["group_maps"], "run_distance_to_selection_criterion.xlsx"),
                index=False)
            # F8b with the arousal state of every run (circles awake,
            # triangles anaesthetised; Mouse and Marmoset_Hori split and
            # tested within animal); the same figure is paper panel 6c.
            try:
                fpp.fig6c_f8b(paths, dist_df_b, bids_configs, final_round,
                              n_perm=cfg["N_PERM_PHYLOGENY"],
                              paper=cfg.get("RUN_PAPER_FIGURES", True),
                              criterion_note=crit_note)
            except Exception as e:
                import traceback
                traceback.print_exc()
                print(f"  [F8b] [WARN] state version failed ({e}) -- plain F8b drawn")
                fa.plot_run_distance_distribution(
                    dist_df_b, output_path=os.path.join(
                        paths["figures"], "F8b_run_distance_to_selection_criterion.png"),
                    title=f"Per-run correlation to the round-{final_round} SELECTION "
                          f"CRITERION (what the ranking used, not the AFNI map)")

    if cfg["MAKE_FIGURES"] and cfg["RUN_PAPER_FIGURES"]:
        fpp.organise_figures(paths["figures"])
    _write_methods(cfg, paths, sorted(profiles_by_label), n_regions, unit_counts,
                   final_round=final_round, phylo_source=phylo_source,
                   n_species=len({species_of_label.get(l) for l in phylo_profiles}))
    try:
        _write_methods_article(st, cfg, paths, phylo_source=phylo_source)
    except Exception as e:
        print(f"  [WARN] Methods_article.txt failed: {e}")
    print(f"\n{'#' * 70}\n  DONE -- {paths['root']}\n{'#' * 70}")
    return dict(config=cfg, paths=paths, common_regions=common,
                multilevel_f6=multilevel_results,
                run_distance_to_ref=dist_df, qc_vs_distance=corr_df,
                profiles_by_label=profiles_by_label,
                members_by_label=members_by_label,
                runs_by_label=runs_by_label,
                unit_counts=unit_counts,
                selection=sel, invariance=invariance, region_mask=region_mask,
                all_units_index=all_units_index, map_index=map_index,
                all_round_selections=all_round_selections, final_round=final_round,
                selections_final=selections_final, phylo_profiles=phylo_profiles,
                phylo_source=phylo_source, fingerprints=fingerprints,
                run_index=st["run_index"])
