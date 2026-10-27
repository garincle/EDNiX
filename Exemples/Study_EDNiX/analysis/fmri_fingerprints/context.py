"""
Inputs of every fMRI fingerprint statistic, from the saved selection (07_state).

load_state(root)  the state of Statistics/Group_fMRI/run_selection + the statistics
                  settings (settings.FIGURE_SETTINGS), time tree loaded, profile method
                  of the selection restored -- what exploration_pipeline.run(part="figures")
                  does before its first figure
prepare(st)       the objects the statistics share (phylogeny source, kept references,
                  selected units, region spaces), built with the code of the first part of
                  exploration_pipeline._run_figures

_DIVERGENCE_DEFAULT, _load_time_tree and _kept_references were moved VERBATIM from
sensitivity_analysis/fmri_selection/exploration_pipeline.py (2026-10-27), which imports
them back.
"""
import os

import numpy as np
import pandas as pd

import Exemples.Study_EDNiX.analysis.fmri_fingerprints.profile_stats as fa
import Exemples.Study_EDNiX.analysis.fmri_fingerprints.phylo as fp
import Exemples.Study_EDNiX.analysis.fmri_fingerprints.fingerprints as ff
import Exemples.Study_EDNiX.analysis.fmri_fingerprints.paper_panels as fpp
from Exemples.Study_EDNiX.analysis.fmri_fingerprints.settings import FIGURE_SETTINGS
from Statistics.Group_fMRI.run_selection import profiles as sri
from Statistics.Group_fMRI.run_selection.selection import load_state as _load_saved
from Statistics.Group_fMRI.run_selection.selection import _merge_figure_config


_DIVERGENCE_DEFAULT = dict(fa.DIVERGENCE_TIME_MYA)   # hand-entered values, at import


def _load_time_tree(cfg, bids_configs):
    """
    Load the time-calibrated Newick tree (cfg["TREE_PATH"]) and make it the
    single source of divergence times: phylo uses it for every pairwise
    age (Mantel, K, PGLS), and fa.DIVERGENCE_TIME_MYA (F6/F7, multi-level)
    is overwritten with the MRCA age of each species with Human.
    """
    # reset module state first: a tree loaded by an earlier run() in the same
    # Python session must not leak into this one
    fp._TREE_PAIR_AGES, fp._TREE_SOURCE = None, None
    fa.DIVERGENCE_TIME_MYA.clear()
    fa.DIVERGENCE_TIME_MYA.update(_DIVERGENCE_DEFAULT)
    if not cfg["TREE_PATH"]:
        print("  [tree] TREE_PATH=None -- using the hand-entered approximate ages "
              "(fa.DIVERGENCE_TIME_MYA / fp.SPECIES_CLADES); not for publication")
        return
    if not os.path.isfile(cfg["TREE_PATH"]):
        # "è" can be stored as one code point (NFC) or two (NFD) depending on
        # how the file was copied: look for the same name under the other form.
        import unicodedata
        d, name = os.path.split(cfg["TREE_PATH"])
        if os.path.isdir(d or "."):
            want = unicodedata.normalize("NFC", name)
            hit = [f for f in os.listdir(d or ".") if unicodedata.normalize("NFC", f) == want]
            if hit:
                print(f"  [tree] {cfg['TREE_PATH']!r} matched on disk as {hit[0]!r} "
                      f"(different Unicode normalisation of the same name)")
                cfg["TREE_PATH"] = os.path.join(d, hit[0])
    if not os.path.isfile(cfg["TREE_PATH"]):
        raise FileNotFoundError(
            f"TREE_PATH {cfg['TREE_PATH']!r} not found. Point TREE_PATH to the "
            f"time-calibrated Newick file, or set TREE_PATH=None to fall back to "
            f"approximate hand-entered ages.")
    species = sorted({c["species"] for c in bids_configs})
    unmapped = [sp for sp in species if sp not in cfg["SPECIES_TO_TREE_TIP"]]
    if unmapped:
        raise KeyError(f"SPECIES_TO_TREE_TIP has no tree tip for {unmapped}")
    mapping = {sp: cfg["SPECIES_TO_TREE_TIP"][sp] for sp in species}
    fp.use_newick_tree(cfg["TREE_PATH"], mapping)
    old = dict(fa.DIVERGENCE_TIME_MYA)
    for sp in species:
        fa.DIVERGENCE_TIME_MYA[sp] = 0.0 if sp == "Human" else fp.mrca_age("Human", sp)
    for sp in species:
        print(f"  [tree] {sp:11s} ({mapping[sp]}): divergence from Human = "
              f"{fa.DIVERGENCE_TIME_MYA[sp]:7.3f} Mya (hand-entered value was "
              f"{old.get(sp, float('nan'))})")


def _kept_references(sri, cfg, profiles_by_label, selections, selections_final, final_round):
    """Per dataset, the reference (REFERENCE_METHOD / K_REF) of its round-1
    kept units and of its final-round kept units -- F1a/F1b, F4, fingerprints."""
    kept_ref_by_label, kept_ref_final_by_label = {}, {}
    for label, profs in profiles_by_label.items():
        kept_keys_r1 = {tuple(str(x) for x in k) for k in selections.get(label, [])}
        kept_profs_r1 = {k: v for k, v in profs.items()
                         if tuple(str(x) for x in k) in kept_keys_r1}
        try:
            if len(kept_profs_r1) >= 2:
                kept_ref_by_label[label], _ = fa.build_reference_for_profiles(
                    sri, kept_profs_r1, method=cfg["REFERENCE_METHOD"],
                    k=min(cfg["K_REF"], len(kept_profs_r1)),
                    value_col=cfg["VALUE_COL"], metric=cfg["METRIC"])
        except Exception as e:
            print(f"  [WARN] {label}: round-1 kept reference failed: {e}")
        if final_round >= 2:
            kept_keys_r2 = {tuple(str(x) for x in k)
                            for k in selections_final.get(label, [])}
            kept_profs_r2 = {k: v for k, v in profs.items()
                             if tuple(str(x) for x in k) in kept_keys_r2}
            try:
                if len(kept_profs_r2) >= 2:
                    kept_ref_final_by_label[label], _ = fa.build_reference_for_profiles(
                        sri, kept_profs_r2, method=cfg["REFERENCE_METHOD"],
                        k=min(cfg["K_REF"], len(kept_profs_r2)),
                        value_col=cfg["VALUE_COL"], metric=cfg["METRIC"])
            except Exception as e:
                print(f"  [WARN] {label}: round-{final_round} kept reference failed: {e}")
    return kept_ref_by_label, kept_ref_final_by_label


# ═══════════════════════════════════════════════════════════════════════════
# Entry points of the paper statistics
# ═══════════════════════════════════════════════════════════════════════════
def load_state(root, settings=FIGURE_SETTINGS):
    """Saved selection of `root` + statistics settings, ready for the statistics:
    the same steps as exploration_pipeline.run(part="figures") before _run_figures
    (selection settings checked, time tree loaded, profile method of the selection
    set on the profiles module). No output folder is created."""
    state = _load_saved(root)
    cfg = _merge_figure_config(state["cfg"], settings)
    _load_time_tree(cfg, state["bids_configs"])
    if hasattr(sri, "set_profile_method"):
        sri.set_profile_method(cfg["PROFILE_METHOD"], cfg["TAIL_Q"], keep_sketch=False)
        if cfg["PROFILE_METHOD"] != "median_r":
            sri.set_tail_reference(state["anat_common"])
    return dict(state, cfg=cfg)


def prepare(st):
    """Objects shared by the statistics, built exactly as at the start of
    exploration_pipeline._run_figures (same code, same order):
      ss                     selection-space objects (F3 / F5 are drawn in the regions
                             the selection scored)
      ana_species_profiles   per-dataset profiles of every unit (F1c, region order)
      kept_runs_final / _r1  run keys behind the kept units
      phylo_profiles         group profiles of Figure 4 (PHYLOGENY_USE_FINAL_ROUND:
                             kept maps of the final round when every dataset has one)
      kept_ref_by_label / kept_ref_final_by_label   references of the kept units
      selected_units         {label: {unit key: profile}} kept by the final round (4b)
    """
    cfg, bids_configs = st["cfg"], st["bids_configs"]
    unit = cfg["UNIT"]
    species_of_label = st["species_of_label"]
    profiles_by_label, members_by_label = st["profiles_by_label"], st["members_by_label"]
    all_units_index, afni_profiles = st["all_units_index"], st["afni_profiles"]
    common = st["common"]
    species_profs, invariance, region_mask = (st["species_profs"], st["invariance"],
                                              st["region_mask"])
    sel, selections = st["sel"], st["selections"]
    final_round = st["final_round"]
    selections_final, final_criterion = st["selections_final"], st["final_criterion"]
    map_index_final = st["map_index_final"]
    group_map_profiles = st["group_map_profiles"]
    split = st.get("selection_space") is not None
    ss = st["selection_space"] if split else dict(
        common=common, profiles_by_label=profiles_by_label, species_profs=species_profs,
        invariance=invariance, region_mask=region_mask, final_criterion=final_criterion)
    ana_species_profiles = species_profs if split else sel["species_profiles"]
    kept_runs_final = fpp.selected_run_keys(selections_final, members_by_label)
    kept_runs_r1 = fpp.selected_run_keys(selections, members_by_label)

    # ---- phylogeny source (same code as _run_figures) ------------------------
    phylo_profiles = afni_profiles
    phylo_index = all_units_index
    phylo_source = f"all-{unit}s group maps (no selection)"
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
    if not phylo_profiles:
        phylo_profiles = {
            lab: (species_of_label.get(lab), sri.compute_group_median_profile(profs, set(profs)))
            for lab, profs in profiles_by_label.items()}
        phylo_source = f"median of all {unit} profiles (no AFNI group map, no selection)"
    print(f"  [PHYLO] Figure 4 group profiles: {phylo_source}")

    kept_ref_by_label, kept_ref_final_by_label = _kept_references(
        sri, cfg, profiles_by_label, selections, selections_final, final_round)

    # units kept by the final round (paper_panels.make_paper_figures, panel 5b)
    selected_units = {}
    for lab, profs in profiles_by_label.items():
        kept = {fpp._key(x) for x in selections_final.get(lab, [])}
        selected_units[lab] = {k: v for k, v in profs.items() if fpp._key(k) in kept}

    return dict(st, ss=ss, split=split, ana_species_profiles=ana_species_profiles,
                kept_runs_final=kept_runs_final, kept_runs_r1=kept_runs_r1,
                phylo_profiles=phylo_profiles, phylo_index=phylo_index,
                phylo_source=phylo_source, kept_ref_by_label=kept_ref_by_label,
                kept_ref_final_by_label=kept_ref_final_by_label,
                selected_units=selected_units, divergence=dict(fa.DIVERGENCE_TIME_MYA))

