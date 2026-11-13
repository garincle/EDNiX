"""
paper_stats.py -- every statistic of the EDNiX paper, in figure order.

One function per figure. Each one computes the numbers of its figure, writes them to
<PAPER_DIR>/<Figure>/<Figure>_stats.xlsx (one sheet per table) and _stats.txt (result
lines), and returns them with what paper_figures.py needs to draw the figure.

  figure1   atlas hierarchy                 no statistic
  figure2   cortical morphometry, level 2   per species median, IQR, n subjects
  figure3   longitudinal morphometry        spline LMM, linear LME, association vs primary
  figure4   selection of the fMRI runs      a  r of every dataset with the LOSO consensus
                                            b, d, e selection vs 10 random runs per dataset:
                                               noise (pair classes), species identity (d),
                                               phylogenetic gradient (e) (former 5b histograms
                                               and per-dataset convergence: tables only)
                                            c  dataset x dataset distance + sibling test
  figure5   seed maps after selection       runs kept / not kept per dataset, colour ranges
  figure6   multimodal profile (PET, FC,    a-b bratio of 2 human + 1 macaque FDG datasets, each human
            morphometry), Human / Macaque      dataset vs macaque (permutation, BH-FDR); c agreement of the
                                               human datasets (r, CCC, noise ceiling, pair classes);
                                               d z-scored profiles, cross-species / cross-modal r
  supp1-3   pipeline schematics (anat MRI,  no statistic
            fMRI, FDG-PET with the CT path)
  supp4     acquisition parameters          the <= 4 parameters reported for every fMRI dataset
  supp5     anatomical QC                   per dataset summary of every panel
  supp6     functional QC                   per dataset summary of every panel
  (supp7    FreeSurfer vs EDNiX, human      validation/supp_freesurfer_comparison.py, not here)
  supp8     longitudinal vs cross-sectional same scans processed with and without the longitudinal
            (Fig. 3 datasets)                  option: residual SD, slope SE and bias (subject bootstrap)
  supp9     human association vs primary    paired t-tests (BIDS_Park)
  supp10    behaviour of the selection      convergence, consensus by dataset, invariance
  supp11    selection vs variability levels animal / session gaps vs random subsets;
                                            b raw species gap, selected vs random
  supp12    proximity of the datasets       sibling rank of every dataset
  supp13    scores and arousal state        selection scores; r with the criterion by state;
                                            units kept per state
  supp14    run dispersion and QC           SD of r with the kept map; QC meta-analysis
  supp15    FDG profile similarity          pairs of subjects by class, label-permutation test
  stable1   acquisition parameters          per dataset: n, field, scanner, TR, TE, flip, voxel,
                                            volumes, duration, arousal state; implausible TRs flagged

fMRI statistics (figure4, figure5, supp10-supp14) read only the saved run selection
(Statistics/Group_fMRI/run_selection, <FMRI_SELECTION_ROOT>/07_state), with the settings
of fmri_fingerprints/settings.py: the functions and the settings are those of the former
Fig5 pipeline (sensitivity_analysis/fmri_selection/exploration_pipeline.py), so every
number is the one that pipeline wrote (checked on synthetic data: validation/).

    python paper_stats.py                  # every figure
    python paper_stats.py figure4 supp12   # some figures
Everything (statistics and figures): run_paper.py.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_data as pdat

opj = os.path.join

SECTIONS = ("figure1", "figure2", "figure3", "figure4", "figure5", "figure6", "supp1", "supp2",
            "supp3", "supp4", "supp5", "supp6", "supp8", "supp9", "supp10", "supp11", "supp12",
            "supp13", "supp14", "supp15", "stable1")
FOLDER = dict(figure1="Figure_1", figure2="Figure_2", figure3="Figure_3",
              figure4="Figure_4", figure5="Figure_5", figure6="Figure_6", supp1="Supp_Fig_01",
              supp2="Supp_Fig_02", supp3="Supp_Fig_03", supp4="Supp_Fig_04",
              supp5="Supp_Fig_05",
              supp6="Supp_Fig_06", supp9="Supp_Fig_09", supp10="Supp_Fig_10",
              supp13="Supp_Fig_13", supp12="Supp_Fig_12", supp11="Supp_Fig_11",
              supp14="Supp_Fig_14", supp15="Supp_Fig_15", supp8="Supp_Fig_08",
              stable1="Supp_Table_01")


def out_base(section, paper_dir=None):
    """<PAPER_DIR>/<Figure>/<Figure>: path of a figure without extension."""
    name = FOLDER[section]
    return opj(paper_dir or pdat.PAPER_DIR, name, name)


# ═══════════════════════════════════════════════════════════════════════════
# Inputs, loaded once and shared by the sections
# ═══════════════════════════════════════════════════════════════════════════
class Inputs:
    """Lazy loader: each input is read the first time a section asks for it.
    fmri_root: folder of the run selection (07_state inside)."""

    def __init__(self, fmri_root=None, refresh=False):
        self.fmri_root = fmri_root or pdat.FMRI_SELECTION_ROOT
        self.refresh = refresh
        self._cache = {}

    def _get(self, key, fn):
        if key not in self._cache:
            self._cache[key] = fn()
        return self._cache[key]

    def morphometry(self):
        from Exemples.Study_EDNiX.analysis.library import morphometry_figure as mf

        def load():
            regions = pdat.CORTICAL_REGIONS[mf.LEVEL]
            morpho = pdat.load_morphometry(level=mf.LEVEL, refresh=self.refresh)
            morpho = {m: pdat.select_regions(df, regions) if not df.empty else df
                      for m, df in morpho.items()}
            return morpho, regions
        return self._get("morphometry", load)

    def figure3_rows(self):
        from Exemples.Study_EDNiX.analysis.library import longitudinal as lg
        return self._get("figure3_rows", lg.load_figure3_data)

    def longitudinal_ablation(self):
        from Exemples.Study_EDNiX.analysis.library import longitudinal_ablation as la
        return self._get("longitudinal_ablation", la.load_ablation_data)

    def human_park(self):
        from Exemples.Study_EDNiX.analysis.library import human_assoc_primary as hap
        return self._get("human_park", hap.load_supp6_data)

    def qc(self, kind):
        from Exemples.Study_EDNiX.analysis.library import qc_figures as qf
        if kind == "anat":
            return self._get("qc_anat", lambda: qf.add_subcortical_snr(
                pdat.load_qc(pdat.ANAT_DATASETS, 'qc_anat', refresh=self.refresh)))
        return self._get("qc_func", lambda: pdat.load_qc(pdat.FUNC_DATASETS, 'qc_func',
                                                         refresh=self.refresh))

    def acquisition(self):
        from Exemples.Study_EDNiX.analysis.library import acquisition_table as acq
        return self._get("acquisition", lambda: acq.load(refresh=self.refresh))

    def multimodal(self):
        from Exemples.Study_EDNiX.analysis.library import multimodal_figure as mmf
        return self._get("multimodal", lambda: mmf.figure6_stats(fmri=self.fmri))

    def fmri(self):
        """Saved selection + shared objects (fmri_fingerprints.context) + the group
        fingerprints used by Figure 4 and Supp. 10 (fmri_group)."""
        def load():
            from Exemples.Study_EDNiX.analysis.fmri_fingerprints import context
            c = context.prepare(context.load_state(self.fmri_root))
            c.update(fmri_group(c))
            return c
        return self._get("fmri", load)


# ═══════════════════════════════════════════════════════════════════════════
# Shared fMRI objects (same code as fingerprints.make_fingerprint_figures)
# ═══════════════════════════════════════════════════════════════════════════
def fmri_group(c):
    """Group fingerprints of every dataset and what Figure 4 and Supp. 10 share:
      vecs_all   {label: scaled group profile} of the Figure 4 group maps
                 (context.prepare: round-final KEPT maps, PHYLOGENY_USE_FINAL_ROUND)
      dist_all   per dataset r / 1 - r with the LOSO consensus and with Human
      order      region order of the radars (lobes, Figure 1)
      D, sibling dataset x dataset 1 - r, sibling test
      resolvers, dataset_state   arousal state of every run / dataset"""
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import fingerprints as ff
    cfg = c["cfg"]
    value_col, norm = cfg["VALUE_COL"], cfg["ZSCORE"]
    sol, common, div = c["species_of_label"], c["common"], c["divergence"]
    ref = cfg["REFERENCE_SPECIES"]
    cfg_by_label = {b["label"]: b for b in c["bids_configs"]}
    resolvers = {lab: ff.StateResolver(cfg_by_label.get(lab, {})) for lab in c["runs_by_label"]}
    ds_state = {lab: ff.dataset_state(resolvers[lab], list(c["runs_by_label"][lab]))
                for lab in c["runs_by_label"]}
    vecs_all = ff.group_vectors(c["phylo_profiles"], common, value_col, norm)
    dist_all = ff.dataset_distances(vecs_all, sol, div, ref)
    matrix_mode = cfg.get("PROFILE_KIND", "seed") == "matrix" or ff.is_edge_profile(common)
    vecs_f1 = ff.group_vectors(c["ana_species_profiles"], common, value_col, norm)
    if matrix_mode or not vecs_f1:
        order = np.arange(len(common))
    else:
        order = ff._region_order(ff.consensus(vecs_f1, sol), common)
    D, sibling = None, None
    if len(vecs_all) >= 3:
        labels = ff.order_labels(vecs_all, sol, div)
        D = ff.distance_matrix(vecs_all, labels)
        sibling = ff.sibling_test(D, sol, n_perm=cfg["N_PERM_PHYLOGENY"])
    return dict(vecs_all=vecs_all, dist_all=dist_all, order=order, matrix_mode=matrix_mode,
                D=D, sibling=sibling, resolvers=resolvers, dataset_state=ds_state)


def _sibling_tables(sibling):
    if not sibling:
        return {}
    return dict(sibling_ranks=sibling["per_dataset"],
                sibling_test=pd.DataFrame([{k: v for k, v in sibling.items()
                                            if k != "per_dataset"}]))


def _sibling_lines(sibling):
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import fingerprints as ff
    if not sibling:
        return ["sibling test: not computed (fewer than 3 datasets or no species with "
                "two datasets)"]
    lines = [f"sibling test: mean 1 - r same species {sibling['mean_same']:.3f} "
             f"({sibling['n_same']} pairs) vs different species {sibling['mean_cross']:.3f} "
             f"({sibling['n_cross']} pairs); difference {sibling['statistic']:+.3f}, "
             f"one-sided p = {ff._fmt_p(sibling['p_one_sided'])} ({sibling['n_perm']} "
             f"permutations of the species labels across datasets)"]
    for r in sibling["per_dataset"].itertuples():
        if r.sibling is not None:
            lines.append(f"  {r.label}: sibling {r.sibling} at 1 - r = {r.d_sibling:.3f}, "
                         f"rank {r.sibling_rank:.0f} of {r.n_others}; nearest "
                         f"{r.nearest_neighbour} ({r.d_nearest:.3f})")
    return lines


# ═══════════════════════════════════════════════════════════════════════════
# Figures
# ═══════════════════════════════════════════════════════════════════════════
def figure1(inp):
    """Atlas hierarchy: drawing only."""
    return dict(tables={}, lines=["Figure 1 has no statistic (atlas hierarchy)."])


def figure2(inp):
    """Per species, region and modality: n subjects, median, quartiles (the values of
    the violins; one value per subject, sessions averaged, bilateral)."""
    from Exemples.Study_EDNiX.analysis.library import morphometry_figure as mf
    morpho, regions = inp.morphometry()
    table = mf.figure2_table(morpho, regions)
    lines = [f"atlas level {mf.LEVEL}, {len(regions)} regions, bilateral values, one value per "
             f"subject (sessions averaged); y axis log10"]
    for mod, g in table.groupby("modality", sort=False) if not table.empty else ():
        n = g.groupby("species")["n"].max()
        lines.append(f"{mod}: subjects per species " + ", ".join(f"{k} {v}" for k, v in n.items()))
    return dict(tables=dict(values_by_species=table), lines=lines,
                draw=dict(morpho=morpho, regions=regions))


def figure3(inp):
    """Trajectories (spline LMM + linear LME) and association vs primary (paired t)."""
    from Exemples.Study_EDNiX.analysis.library import longitudinal as lg
    st = lg.figure3_stats(inp.figure3_rows())
    return dict(tables=dict(panels=st["table"]), lines=st["report"], text=st["text"],
                draw=dict(stats=st))


def figure4(inp):
    """a  r of each dataset's group fingerprint with the leave-one-species-out consensus (the
          dashed profile of panel a; r with the Human datasets kept as a second column),
          Figure 4 group maps
       5b run-to-run profile correlation: same species / other dataset vs different
          species, every run and the runs kept by the final round (tables and lines only)
       b, d, e selection vs 10 random runs per dataset (selection_validation.py): b noise
          removed (pair classes), d phylogenetic gradient (Mantel, raw profiles), e species
          identity kept (residual profiles); per-dataset convergence kept in the tables
       c  dataset x dataset similarity r (drawn) / distance 1 - r and sibling test (F16a)"""
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import profile_stats as fa
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import paper_panels as fpp
    c = inp.fmri()
    cfg = c["cfg"]
    value_col, norm = cfg["VALUE_COL"], cfg["ZSCORE"]
    sol, common = c["species_of_label"], c["common"]
    intra = cfg.get("WITHIN_ACROSS_PAIRS", "other_dataset")
    test = cfg.get("WITHIN_ACROSS_TEST", "wilcoxon")
    lines = [f"group maps: {c['phylo_source']}; {len(common)} regional values "
             f"({cfg['VALUE_COL']}, {fa._norm_label(norm).strip() or 'unscaled'})"]
    # a
    dist = c["dist_all"]
    lines.append("a  r with the leave-one-species-out consensus (median of the other species), "
                 "panel order")
    for r in dist.sort_values("r_consensus", ascending=False).itertuples():
        lines.append(f"  {r.label}: r = {r.r_consensus:.3f} (r with the Human datasets = "
                     f"{r.r_reference:.3f})")
    # b (same pairs and statistics as fa.plot_within_vs_across_all_vs_selected)
    pools = (("All runs", c["profiles_by_label"]), ("Selected", c["selected_units"]))
    pairs, stats = {}, {}
    for name, pbl in pools:
        pairs[name] = fa._within_across_species_pairs(pbl, common, sol, value_col=value_col,
                                                      zscore=norm, max_pairs=20000,
                                                      random_seed=0)
        stats[name] = fa._within_across_stats(pairs[name], intra=intra)
    named = {("all runs" if k == "All runs" else "selected"): v for k, v in stats.items()}
    lines.append(f"b  pairs compared: {intra}; test on the figure: {test}")
    lines += ["  " + ln for ln in fpp.within_across_result_lines(named, test=test,
                                                                 intra=intra).splitlines()]
    n_units = {k: sum(len(p) for p in pool.values()) for k, pool in pools}
    lines.append(f"  units: all {n_units['All runs']}, selected {n_units['Selected']}")
    # c
    lines.append("c  dataset x dataset similarity: r of the group fingerprints (matrix of the "
                 "figure = 1 - distance below)")
    lines += ["  " + ln for ln in _sibling_lines(c["sibling"])]
    # b, d, e (Clément 2026-10-30): what the selection buys vs 10 random runs per dataset
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import selection_validation as sv
    val = sv.compute(c)
    lines.append("b, d, e  selection vs random choice of the same number of runs "
                 "(fmri_fingerprints/selection_validation.py); the histograms above are kept "
                 "in the tables, no longer drawn")
    lines += ["  " + ln for ln in sv.result_lines(val)]
    tables = dict(a_distances=dist,
                  heldout_random_draws=pd.DataFrame(val["heldout"]["rand"],
                                                      columns=["between_datasets",
                                                               "with_other_species", "G"]),
                  heldout_selected=pd.DataFrame([dict(zip(
                      ["between_datasets", "with_other_species", "G"], val["heldout"]["sel"]),
                      **{f"p_{k}": v for k, v in zip(
                          ["between_datasets", "with_other_species", "G"],
                          val["heldout"]["p"])})]),
                  c_heldout_per_dataset=val["heldout"]["per_dataset"],
                  b_pair_levels=val["levels"],
                  d_species_signal=pd.DataFrame(dict(
                      pool=["all runs", "selected"] + ["random"] * len(val["species"]["rand"]),
                      gap=[val["species"]["gap_all"], val["species"]["gap_sel"],
                           *val["species"]["rand"]],
                      p_label_permutation=[val["species"]["p_all"], val["species"]["p_sel"],
                                           *val["species"]["rand_p"]])),
                  e_phylo_heldout=pd.DataFrame(dict(
                      pool=["selected"] + ["random"] * len(val["phylo_heldout"]["random"]),
                      rho=[val["phylo_heldout"]["selected"], *val["phylo_heldout"]["random"]],
                      p=[val["phylo_heldout"]["p"]] + [np.nan]
                      * len(val["phylo_heldout"]["random"]))),
                  e_phylo_delta=pd.DataFrame([dict(profiles=k, **v_) for k, v_ in
                                              val["species"]["mantel_delta"].items()]),
                  e_phylo_mantel=pd.DataFrame(
                      [dict(profiles=k, pool=pool, rho=rho, p=p_)
                       for k, m in val["species"]["mantel"].items()
                       for pool, (rho, p_) in (("all runs", m["all"]),
                                               ("selected", m["selected"]))]
                      + [dict(profiles=k, pool="random", rho=r_, p=p_)
                         for k, m in val["species"]["mantel"].items()
                         for r_, p_ in zip(m["random"], m["random_p"])]),
                  b_pair_class_p=pd.DataFrame(dict(
                      pair_class=["same dataset", "same species, other dataset",
                                  "different species"],
                      p_selected_vs_random=val["heldout"]["p_classes"])),
                  b_within_across=pd.DataFrame([dict(pool=k, **v) for k, v in named.items()
                                                if v]),
                  c_similarity_matrix=1 - c["D"] if c["D"] is not None else pd.DataFrame(),
                  c_distance_matrix=c["D"] if c["D"] is not None else pd.DataFrame(),
                  **{f"c_{k}": v for k, v in _sibling_tables(c["sibling"]).items()})
    return dict(tables=tables, lines=lines,
                draw=dict(pairs=pairs, stats=stats, test=test, intra=intra, validation=val))


def figure5(inp):
    """Seed maps: per dataset, units kept by the final round, maps drawn (all / kept /
    same number not kept) and the colour range shared by the three maps."""
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import paper_panels as fpp
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import fingerprints as ff
    c = inp.fmri()
    maps = fpp.surface_maps_for_5a(c["all_units_index"], c["map_index_final"])
    ranges = c.get("surface_ranges") or {}
    rows = []
    for lab in ff.order_labels(list(c["profiles_by_label"]), c["species_of_label"],
                               c["divergence"]):
        rows.append(dict(label=lab, species=c["species_of_label"].get(lab),
                         n_units=len(c["profiles_by_label"][lab]),
                         n_kept_final=len(c["selections_final"].get(lab, [])),
                         map_all=maps.get(lab, {}).get("all"),
                         map_kept_final=maps.get(lab, {}).get("kept_final"),
                         map_not_kept_final=maps.get(lab, {}).get("not_kept_final"),
                         colour_range=ranges.get(lab, np.nan)))
    t = pd.DataFrame(rows)
    lines = [f"final round {c['final_round']}; one symmetric colour range per dataset "
             f"(+/- P{c['cfg'].get('SURFACE_COMMON_PCT', 90.0):g} of |z|, shared by its maps)"]
    lines += [f"  {r.label}: {r.n_kept_final} kept of {r.n_units} {c['cfg']['UNIT']}s; range "
              f"+/-{r.colour_range:.3g}" for r in t.itertuples()]
    return dict(tables=dict(datasets=t), lines=lines,
                draw=dict(maps=maps, pics=c.get("surface_pics") or {}, ranges=ranges))


def figure6(inp):
    """Multimodal profile: a-b bratio of the FDG datasets (two human, one macaque) and difference with the macaque,
    c reproducibility of the human datasets, d multimodal fingerprints (library/multimodal_figure.py,
    library/pet_datasets.py)."""
    from Exemples.Study_EDNiX.analysis.library import multimodal_figure as mmf
    st = inp.multimodal()
    pet = st["pet"]
    w, sb = pet["w"], pet["stats"]
    lines = [f"level {mmf.LEVEL}, {len(pet['regions'])} regions (a-c), {len(st['regions_c'])} regions (d); "
             f"test {mmf.TEST} ({mmf.N_PERM} permutations), BH-FDR within dataset; reference {mmf.REFERENCE_DATASET}"]
    lines.append("B  n subjects: " + ", ".join(f"{d} {w.loc[d].shape[0]}" for d in pet["labels"]
                                               if d in w.index.get_level_values(0)))
    lines += [f"  {r.dataset} | {r.region}: {r.dataset} - {r.reference} = {r.difference:+.3f} "
              f"(g = {r.hedges_g:+.2f}), q = {r.q_fdr:.3g}" for r in sb.itertuples()]
    g = pet["agree"].iloc[0]
    lines.append(f"C  agreement {g.dataset_a} vs {g.dataset_b} ({int(g.n_regions)} regions, n {int(g.n_a)} / "
                 f"{int(g.n_b)}): Pearson r = {g.pearson_r:.3f} [{g.pearson_r_ci_low:.3f}, {g.pearson_r_ci_high:.3f}], "
                 f"Spearman = {g.spearman_r:.3f}, Lin CCC = {g.lin_ccc:.3f} [{g.lin_ccc_ci_low:.3f}, "
                 f"{g.lin_ccc_ci_high:.3f}], mean |diff| = {g.mean_abs_difference:.3f}; split-half "
                 f"{g.reliability_a:.3f} / {g.reliability_b:.3f}, disattenuated r = {g.pearson_r_disattenuated:.3f}")
    bt = pet["between"]
    lines.append(f"C  {bt.dataset.iat[0]} - {bt.reference.iat[0]}: {int((bt.q_fdr < mmf.ALPHA).sum())} of {len(bt)} "
                 f"regions q < {mmf.ALPHA}")
    lines += [f"  {r.region}: {r.difference:+.3f} (g = {r.hedges_g:+.2f}), q = {r.q_fdr:.3g}"
              for r in bt.itertuples() if r.q_fdr < mmf.ALPHA]
    pr = pet["pairs"]
    if len(pr):
        lines += [f"C  pairs {c}: n = {len(v)}, median r = {v.median():.3f} [IQR {v.quantile(.25):.3f}-"
                  f"{v.quantile(.75):.3f}]" for c, v in pr.groupby("pair_class", sort=False)["r"]]
    for t in pet["ptest"].itertuples():
        lines.append(f"C  {t.test}: {t.statistic:+.4f} (null {t.null_mean:+.4f} +- {t.null_sd:.4f}), "
                     f"P = {t.p_permutation:.3g}, {t.n_subjects} subjects")
    lines += [f"D  {r.species} {r.modality_a}-{r.modality_b}: Spearman r = {r.spearman_r:.2f}"
              for r in st["cmc"].itertuples()]
    summ = st.get("summ")
    if summ is not None and len(summ):
        lines.append("D  summary (Spearman over regions; P permutation, CI bootstrap, Holm over all rows)")
        lines += [f"  {r.analysis} | {r.species} | {r.modality} (n={r.n_regions}): r = {r.r:+.2f} "
                  f"[{r.ci_low:+.2f}, {r.ci_high:+.2f}], P = {r.p:.3g}, Holm P = {r.p_holm:.3g}"
                  for r in summ.itertuples()]
    return dict(tables=mmf.stats_tables(st), lines=lines, draw=dict(stats=st))


def _schematic(n, what):
    return dict(tables={}, lines=[f"Supplementary Figure {n} has no statistic (schematic of the {what})."])


def supp1(inp):
    """Anatomical MRI module: drawing only."""
    return _schematic(1, "anatomical MRI module")


def supp2(inp):
    """fMRI module: drawing only."""
    return _schematic(2, "fMRI module")


def supp3(inp):
    """FDG-PET module with the CT path: drawing only."""
    return _schematic(3, "FDG-PET module with the CT path")


def _qc_section(df, metrics, unit):
    from Exemples.Study_EDNiX.analysis.library import qc_figures as qf
    summ = qf.qc_summary(df, metrics)
    lines = [f"one point per {unit}; n per dataset and median [IQR] of every panel"]
    for (col, label), g in summ.groupby(["metric", "label"], sort=False):
        lines.append(f"{label} ({col}): " + "; ".join(
            f"{r.bids_dir} n={r.n} {r.median:.3g} [{r.q25:.3g}-{r.q75:.3g}]"
            for r in g.itertuples()))
    missing = [m[0] for m in metrics if m is not None and m[0] not in df.columns]
    if missing:
        lines.append(f"no column for: {missing} (panel empty: re-run the EDNiX QC step)")
    return summ, lines


def supp4(inp):
    """Acquisition parameters of the fMRI datasets: at most 4, only those reported for
    every dataset (library/acquisition_table.ACQ_CANDIDATES gives the order of interest)."""
    from Exemples.Study_EDNiX.analysis.library import acquisition_table as acq
    f = acq.acq_runs(inp.acquisition())
    chosen, dropped = acq.acq_metrics(f)
    summ = acq.summary_per_metric(f, chosen)
    return dict(tables=dict(summary=summ), lines=acq.acq_lines(f, chosen, dropped),
                draw=dict(runs=f, chosen=chosen))


def supp5(inp):
    from Exemples.Study_EDNiX.analysis.library import qc_figures as qf
    df = inp.qc("anat")
    metrics = qf.anat_metrics(df)
    summ, lines = _qc_section(df, metrics, "(subject, session)")
    cols = ['species', 'bids_dir', 'subject', 'session'] + [m for m, _ in metrics]
    return dict(tables=dict(summary=summ, values=df[[c for c in cols if c in df.columns]]),
                lines=lines, draw=dict(df=df))


def supp6(inp):
    from Exemples.Study_EDNiX.analysis.library import qc_figures as qf
    df = inp.qc("func")
    metrics = qf.func_metrics(df)
    summ, lines = _qc_section(df, metrics, "BOLD run")
    cols = ['species', 'bids_dir', 'subject', 'session', 'run'] + [
        m[0] for m in metrics if m is not None]
    return dict(tables=dict(summary=summ, values=df[[c for c in cols if c in df.columns]]),
                lines=lines, draw=dict(df=df))


def supp9(inp):
    from Exemples.Study_EDNiX.analysis.library import human_assoc_primary as hap
    morph, age = inp.human_park()
    st = hap.supp6_stats(morph, age)
    return dict(tables=dict(panels=st["table"]), lines=st["lines"], text=hap.stats_text(st),
                draw=dict(stats=st))


def supp10(inp):
    """a  convergence of the iterative selection (F1e)
       b  regional profile of every dataset (all runs) and the consensus (F1c)
       c  invariance of every region; regions used by the score (F5)"""
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import supp_panels as fsp
    c = inp.fmri()
    cfg, ss = c["cfg"], c["ss"]
    conv = c.get("conv_df")
    conv = conv if conv is not None else pd.DataFrame()
    lines = []
    if not conv.empty:
        last = conv[conv["round"] == conv["round"].max()]
        lines.append(f"a  final round {int(conv['round'].max())}: kept runs unchanged from the "
                     f"previous round " + ", ".join(
                         f"{r.label} {r.n_overlap_with_previous}/{r.n_select}"
                         for r in last.itertuples()))
    per, cons = fsp.f1c_values(c["ana_species_profiles"], c["common"], c["species_of_label"],
                               cfg["VALUE_COL"], cfg["ZSCORE"])
    f1c = pd.DataFrame(per, index=c["common"]) if per else pd.DataFrame()
    if per:
        f1c["consensus"] = cons
        lines.append(f"b  {len(per)} datasets, {len(c['common'])} regional values (all runs)")
    inv = pd.DataFrame(dict(region=ss["common"], invariance=np.asarray(ss["invariance"], float),
                            used_by_selection=(np.asarray(ss["region_mask"], bool)
                                               if ss["region_mask"] is not None else True)))
    lines.append(f"c  {int(inv['used_by_selection'].sum())} of {len(inv)} regional values used "
                 f"by the selection score ({cfg['INVARIANCE_MODE']} = {cfg['INVARIANCE_VALUE']})")
    return dict(tables=dict(a_convergence=conv, b_consensus_by_dataset=f1c.reset_index(
        names="region") if per else f1c, c_region_invariance=inv), lines=lines, draw={})


def supp13(inp):
    """a  selection score of every unit, per dataset (F3)
       b  r of every run with its dataset's final-round selection criterion, by arousal
          state; within-animal sign-flip test (F8b / 6c)
       c  units kept in each state (datasets holding both states)"""
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import profile_stats as fa
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import paper_panels as fpp
    c = inp.fmri()
    cfg, sol = c["cfg"], c["species_of_label"]
    rk = []
    for lab, df in c["sel"]["rankings"].items():
        kept = {fpp._key(k) for k in c["selections"].get(lab, [])}
        d = df.copy()
        d.insert(0, "label", lab)
        d["kept_round1"] = [fpp._key(k) in kept for k in
                            zip(d["subject"], d["session"], d["run"], d["seed"])]
        rk.append(d)
    rk = pd.concat(rk, ignore_index=True) if rk else pd.DataFrame()
    lines = ["a  round-1 score (Pearson r with the leave-one-species-out consensus); cutoff = "
             f"score of the {cfg['N_SELECT']}th best unit"]
    for lab, df in c["sel"]["rankings"].items():
        v = pd.to_numeric(df["pearson_r"], errors="coerce").dropna()
        lines.append(f"  {lab}: {len(v)} units, median r = {v.median():.3f}, kept "
                     f"{len(c['selections'].get(lab, []))}")
    ref_crit = {lab: (sol.get(lab), ref) for lab, ref in c["final_criterion"].items()
                if ref is not None}
    dist_b = fa.compute_run_distance_to_reference(c["runs_by_label"], ref_crit, c["common"],
                                                  value_col=cfg["VALUE_COL"],
                                                  kept_keys_by_label=c["kept_runs_final"])
    d, tests, pooled, slines = (fpp.f8b_state_stats(dist_b, c["bids_configs"],
                                                    n_perm=cfg["N_PERM_PHYLOGENY"])
                                if not dist_b.empty else (pd.DataFrame(), {}, None, []))
    lines.append(f"b  r with the round-{c['final_round']} selection criterion, by state")
    lines += ["  " + ln for ln in slines]
    state_t = pd.DataFrame([dict(label=k, **{kk: vv for kk, vv in v.items()
                                             if kk != "subject_diffs"})
                            for k, v in tests.items()])
    kept_t, klines = _kept_by_state(d, list(tests))
    lines += klines
    return dict(tables=dict(a_selection_scores=rk, b_run_r_by_state=d, b_state_tests=state_t,
                            c_kept_by_state=kept_t),
                lines=lines, draw=dict(d=d, tests=tests))


def _kept_by_state(d, labels):
    """Number of units (run x seed) kept by the selection in each arousal state, for the
    datasets holding both states (Clément, V4 comment: 'XX awake and XX anaesthetised runs
    kept'). Fisher exact test (two-sided) of kept vs not kept x awake vs anaesthetised per
    dataset; it treats units as independent, so it is descriptive only (runs are nested in
    animals; the within-animal test of section b is the inferential one)."""
    from scipy.stats import fisher_exact
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import fingerprints as ff
    if d is None or d.empty or "is_kept" not in d or "state_class" not in d:
        return pd.DataFrame(), ["c  kept units by state: not computed (no state or kept column)"]
    rows, lines = [], ["c  units kept by the selection, by arousal state (datasets with both "
                       "states; Fisher exact test, descriptive)"]
    for lab in labels:
        g = d[d["label"] == lab]
        cnt = {}
        for st in (ff.AWAKE, ff.ANESTH):
            s_ = g[g["state_class"] == st]
            cnt[st] = (int(s_["is_kept"].sum()), int(len(s_)))
        (ka, na), (kn, nn) = cnt[ff.AWAKE], cnt[ff.ANESTH]
        p = (fisher_exact([[ka, na - ka], [kn, nn - kn]])[1]
             if min(na, nn) > 0 else np.nan)
        rows.append(dict(label=lab, kept_awake=ka, n_awake=na, kept_anaesthetised=kn,
                         n_anaesthetised=nn, p_fisher=p))
        lines.append(f"  {lab}: kept {ka}/{na} awake, {kn}/{nn} anaesthetised; "
                     f"Fisher p = {ff._fmt_p(p) if p == p else 'n/a'}")
    tot = pd.DataFrame(rows)
    if not tot.empty:
        lines.append(f"  total: kept {int(tot['kept_awake'].sum())} awake and "
                     f"{int(tot['kept_anaesthetised'].sum())} anaesthetised units")
    return tot, lines


def supp12(inp):
    """Proximity of every dataset to its sibling and to the other species (F16b)."""
    c = inp.fmri()
    return dict(tables=dict(distances=c["dist_all"],
                            distance_matrix=c["D"] if c["D"] is not None else pd.DataFrame(),
                            **_sibling_tables(c["sibling"])),
                lines=_sibling_lines(c["sibling"]), draw={})


def supp11(inp):
    """a  inter-run / inter-session / inter-animal similarity, selected runs vs random
          subsets of the same size (F20)
       b  species gap on the raw group profiles (common profile kept), selected vs random"""
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import fingerprints as ff
    c = inp.fmri()
    cfg = c["cfg"]
    if c["matrix_mode"]:
        return dict(tables={}, lines=["matrix profiles: not computed"], draw={})
    vt, vdraws = ff.variability_levels(c["runs_by_label"], c["kept_runs_final"],
                                       c["resolvers"], c["common"], cfg["VALUE_COL"],
                                       cfg["ZSCORE"], n_draw=min(cfg["N_PERM_PHYLOGENY"], 2000))
    pooled = ({k: ff.pooled_gap_test(vdraws, k) for k in ("gap_animal", "gap_session")}
              if not vt.empty else {})
    lines = ff.variability_result_lines(vt, pooled) if not vt.empty else ["no dataset with "
                                                                          "several runs per animal"]
    pooled_t = pd.DataFrame([dict(gap=k, **t) for k, t in pooled.items() if t])
    # b (Clément, V4 comment): species gap on the raw group profiles, selected vs random
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import selection_validation as sv
    sp = sv.compute(c)["species"]
    lines = ["a  " + lines[0]] + lines[1:] + sv.raw_gap_lines(sp)
    return dict(tables=dict(variability_levels=vt, pooled=pooled_t,
                            b_species_gap_raw=sv.raw_gap_table(sp)), lines=lines,
                draw=dict(vt=vt, vdraws=vdraws, species=sp))


def supp14(inp):
    """a  SD across runs of r with the dataset's kept group map (F9)
       b  external QC vs that r: per-dataset correlations, meta-analysis (F10)"""
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import profile_stats as fa
    c = inp.fmri()
    cfg = c["cfg"]
    mif = c["map_index_final"]
    use_final = (c["final_round"] >= 2 and not mif.empty and cfg["DISTANCE_REF_GROUP"] == "kept")
    ref_group = "kept_final" if use_final else cfg["DISTANCE_REF_GROUP"]
    ref_sel = c["kept_runs_final"] if use_final else c["kept_runs_r1"]
    ref_profiles = dict(c["group_map_profiles"].get(ref_group, {}))
    lines = [f"reference = {ref_group} group map of each dataset"]
    if not ref_profiles:
        lines.append(f"no '{ref_group}' group-map profile in the saved selection")
        return dict(tables={}, lines=lines, draw={})
    dist = fa.compute_run_distance_to_reference(c["runs_by_label"], ref_profiles, c["common"],
                                                value_col=cfg["VALUE_COL"],
                                                kept_keys_by_label=ref_sel)
    tables, draw = dict(a_run_r_to_kept_map=dist), dict(dist=dist)
    if not dist.empty:
        sd = dist.groupby("label")["pearson_r"].agg(["count", "mean", "std"]).reset_index()
        tables["a_sd_by_dataset"] = sd
        lines += [f"a  {r.label}: {r.count} runs, mean r = {r.mean:.3f}, SD = {r.std:.3f}"
                  for r in sd.itertuples()]
        if c["df_qc_by_label"]:
            corr = fa.correlate_qc_with_distance(dist, c["df_qc_by_label"],
                                                 qc_metrics=cfg["QC_METRICS_FOR_CORRELATION"])
            tables["b_qc_correlation_per_dataset"] = corr
            if not corr.empty:
                meta = fa.meta_analyze_qc_correlations(corr, min_bids=cfg["QC_META_MIN_BIDS"],
                                                       alpha=cfg["QC_META_ALPHA"])
                tables["b_qc_meta_analysis"] = meta
                draw["meta"] = meta
                if not meta.empty:
                    lines.append("b  meta-analysis across datasets (Fisher z, random effects)")
                    lines += ["  " + ", ".join(f"{k}={v:.3g}" if isinstance(v, float) else
                                               f"{k}={v}" for k, v in r.items())
                              for r in meta.to_dict("records")]
        else:
            lines.append("b  no external QC table in the saved selection")
    return dict(tables=tables, lines=lines, draw=draw)


def supp15(inp):
    """Similarity of the regional FDG profiles between pairs of subjects (former right half of Figure 6c):
    classes same subject / other session, within a dataset, between the human datasets, human vs macaque."""
    pet = inp.multimodal()["pet"]
    pr, lines = pet["pairs"], []
    if len(pr):
        lines += [f"{c}: n = {len(v)} pairs, median r = {v.median():.3f} [IQR {v.quantile(.25):.3f}-"
                  f"{v.quantile(.75):.3f}]" for c, v in pr.groupby("pair_class", sort=False)["r"]]
    for t in pet["ptest"].itertuples():
        lines.append(f"{t.test}: {t.statistic:+.4f} (null {t.null_mean:+.4f} +- {t.null_sd:.4f}), "
                     f"P = {t.p_permutation:.3g}, {t.n_subjects} subjects, {t.n_perm} permutations")
    return dict(tables=dict(pairs=pr, test=pet["ptest"], sessions_used=pet["sessions"]), lines=lines,
                draw=dict(pet=pet))


def supp8(inp):
    """Longitudinal vs cross-sectional processing of the Fig. 3 datasets
    (library/longitudinal_ablation.py): same scans, same LMM as the Fig. 3 inset."""
    from Exemples.Study_EDNiX.analysis.library import longitudinal_ablation as la
    st = la.ablation_stats(inp.longitudinal_ablation())
    return dict(tables=dict(summary=st["summary"], assoc_primary=st["assoc"],
                            paired_values=st["paired"]),
                lines=st["lines"], draw=dict(stats=st))


def stable1(inp):
    """Supplementary Table 1: acquisition parameters of every dataset (BIDS sidecars and
    NIfTI headers), fMRI then anatomy; implausible sidecar TRs flagged, not corrected."""
    from Exemples.Study_EDNiX.analysis.library import acquisition_table as acq
    df = inp.acquisition()
    s = acq.summary(df)
    flagged = df[df["TR_check"].fillna("") != ""] if "TR_check" in df else df.iloc[:0]
    return dict(tables=dict(summary=s, per_image=df, TR_flagged=flagged),
                lines=acq.lines(s, df), draw=dict(summary=s))


# ═══════════════════════════════════════════════════════════════════════════
def write(section, res, paper_dir=None):
    """<Figure>_stats.xlsx (one sheet per table) and <Figure>_stats.txt."""
    base = out_base(section, paper_dir)
    os.makedirs(os.path.dirname(base), exist_ok=True)
    tables = {k: v for k, v in res.get("tables", {}).items()
              if isinstance(v, pd.DataFrame) and not v.empty}
    written = []
    if tables:
        with pd.ExcelWriter(base + "_stats.xlsx") as xw:
            for name, t in tables.items():
                t.to_excel(xw, sheet_name=name[:31], index=name.endswith("matrix"))
        written.append(base + "_stats.xlsx")
    text = res.get("text")
    if text is None:
        text = f"{FOLDER[section]} -- statistics (paper_stats.{section})\n\n" + \
               "\n".join(res.get("lines", [])) + "\n"
    with open(base + "_stats.txt", "w", encoding="utf-8") as fh:
        fh.write(text)
    written.append(base + "_stats.txt")
    return written


def run(sections=SECTIONS, inputs=None, paper_dir=None, write_files=True):
    """Statistics of `sections`; {section: result}. A section failing (missing input)
    is reported and the others run."""
    inp = inputs or Inputs()
    out = {}
    for s in sections:
        if s not in SECTIONS:
            raise KeyError(f"unknown section {s!r}; one of {SECTIONS}")
        print(f"\n{'=' * 70}\n  statistics: {FOLDER[s]}\n{'=' * 70}")
        try:
            out[s] = globals()[s](inp)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"  [{s}] [ERROR] {e}")
            continue
        if write_files:
            for p in write(s, out[s], paper_dir):
                print(f"  [saved] {p}")
    return out


if __name__ == "__main__":
    run(tuple(sys.argv[1:]) or SECTIONS)
