"""
paper_figures.py -- every figure of the EDNiX paper, in figure order, drawn from the
results of paper_stats.py (nothing is computed here).

One function per figure, same names as in paper_stats. Multi-panel figures: the
panels are drawn at print size into <PAPER_DIR>/<Figure>/panels/ and laid out, with
their letters, as <PAPER_DIR>/<Figure>/<Figure>.pdf|png (paper_style.compose_images).
Nature widths: 89 / 120 / 183 mm; text 5-7 pt; PDF with TrueType fonts, PNG 600 dpi.

  figure1  atlas hierarchy (sunburst)                     library/atlas_sunburst.py
  figure2  cortical morphometry, level 2                  library/morphometry_figure.py
  figure3  longitudinal morphometry                       library/longitudinal.py
  figure4  a fingerprint over Human (F15b), b run-to-run r (5b), c distance matrix (F16a)
  figure5  seed maps, selected vs excluded runs (5a)
  figure6  multimodal profile, Human / Macaque (former Figure XX) library/multimodal_figure.py
  supp4/5  anatomical / functional QC                     library/qc_figures.py
  supp6    human association vs primary                   library/human_assoc_primary.py
  supp7    a convergence (F1e), b consensus by dataset (F1c), c region invariance (F5)
  supp8    a selection scores (F3), b r with the criterion by arousal state (6c)
  supp9    proximity rows (F16b)
  supp10   selection vs variability levels (F20)
  supp11   a SD of run r with the kept map (F9), b QC vs that r (F10)
fMRI panels: fmri_fingerprints/ (fingerprints.py, profile_stats.py, paper_panels.py,
supp_panels.py), the drawing functions of the former Fig5 pipeline.

Run with the statistics: run_paper.py.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_style as ps
from Exemples.Study_EDNiX.analysis import paper_stats as pst

opj = os.path.join


def _panel(section, name, paper_dir=None):
    d = opj(os.path.dirname(pst.out_base(section, paper_dir)), "panels")
    os.makedirs(d, exist_ok=True)
    return opj(d, name)


def _png(path):
    """<path>.png when the panel was written, else None (empty cell)."""
    p = path + ".png"
    return p if os.path.isfile(p) else None


def _try(name, fn, *a, **k):
    """A panel failing leaves an empty cell; the error is printed."""
    try:
        return fn(*a, **k)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"  [{name}] [WARN] panel failed: {e}")
        return None


# ═══════════════════════════════════════════════════════════════════════════
def figure1(res, inp, paper_dir=None):
    from Exemples.Study_EDNiX.analysis.library import atlas_sunburst as asb
    base = pst.out_base("figure1", paper_dir)
    asb.write_sunburst_figures(out_dir=opj(os.path.dirname(base), "parts"), paper_out=base)
    return base


def figure2(res, inp, paper_dir=None):
    from Exemples.Study_EDNiX.analysis.library import morphometry_figure as mf
    d = res["draw"]
    return mf.plot_figure2(d["morpho"], d["regions"], pst.out_base("figure2", paper_dir))


def figure3(res, inp, paper_dir=None):
    from Exemples.Study_EDNiX.analysis.library import longitudinal as lg
    base = pst.out_base("figure3", paper_dir)
    return lg.draw_figure3(res["draw"]["stats"], os.path.dirname(base),
                           name=os.path.basename(base))


def figure4(res, inp, paper_dir=None):
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import fingerprints as ff
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import profile_stats as fa
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import supp_panels as fsp
    c, d = inp.fmri(), res["draw"]
    cfg, sol = c["cfg"], c["species_of_label"]
    pa, pb, pc = (_panel("figure4", n, paper_dir) for n in
                  ("a_F15b_fingerprint_vs_human", "b_5b_intra_vs_inter_species",
                   "c_F16a_proximity_matrix"))
    if not c["matrix_mode"] and len(c["dist_all"]) >= 3:
        _try("4a", ff.plot_distance_fingerprint_grid, c["vecs_all"], c["dist_all"], sol,
             "reference", c["common"], c["order"], c["region_mask"], pa,
             cfg["REFERENCE_SPECIES"], c["divergence"], c["phylo_source"], paper=True,
             print_width_mm=183, label_first=True)
    _try("4b", fa.paper_within_across_figure, d["pairs"], d["stats"],
         ["All runs", "Selected"], d["test"], d["intra"], pb)
    if c["D"] is not None:
        _try("4c", fsp.f16a_proximity_matrix, c["D"], sol, pc)
    base = pst.out_base("figure4", paper_dir)
    ps.compose_images([_png(pa), _png(pb), _png(pc)], base, layout=[[0], [1, 2]])
    return base


def figure5(res, inp, paper_dir=None):
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import fingerprints as ff
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import paper_panels as fpp
    c, d = inp.fmri(), res["draw"]
    cfg, sol = c["cfg"], c["species_of_label"]
    labels = ff.order_labels([b["label"] for b in c["bids_configs"] if b["label"] in sol],
                             sol, c["divergence"])
    p = _panel("figure5", "5a_surface_maps_selected_vs_excluded", paper_dir)
    _try("5", fpp.fig5a_surfaces, p, cfg, d["maps"], d["pics"], d["ranges"], labels, sol,
         tuple(cfg.get("SURFACE_PNG_KEEP") or ("",)), supp_out=None)
    base = pst.out_base("figure5", paper_dir)
    ps.compose_images([_png(p)], base, letters=" ")
    return base


def figure6(res, inp, paper_dir=None):
    from Exemples.Study_EDNiX.analysis.library import multimodal_figure as mmf
    return mmf.draw_figure6(res["draw"]["stats"], pst.out_base("figure6", paper_dir))


def supp4(res, inp, paper_dir=None):
    from Exemples.Study_EDNiX.analysis.library import qc_figures as qf
    return qf.plot_supp4(res["draw"]["df"], pst.out_base("supp4", paper_dir))


def supp5(res, inp, paper_dir=None):
    from Exemples.Study_EDNiX.analysis.library import qc_figures as qf
    return qf.plot_supp5(res["draw"]["df"], pst.out_base("supp5", paper_dir))


def supp6(res, inp, paper_dir=None):
    from Exemples.Study_EDNiX.analysis.library import human_assoc_primary as hap
    return hap.draw_supp6(res["draw"]["stats"], pst.out_base("supp6", paper_dir))


def supp7(res, inp, paper_dir=None):
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import supp_panels as fsp
    c = inp.fmri()
    cfg, ss = c["cfg"], c["ss"]
    pa, pb, pc = (_panel("supp7", n, paper_dir) for n in
                  ("a_F1e_convergence", "b_F1c_consensus_by_dataset", "c_F5_region_invariance"))
    if c.get("conv_df") is not None and not c["conv_df"].empty:
        _try("S7a", fsp.f1e_convergence, c["conv_df"], pa)
    _try("S7b", fsp.f1c_consensus_by_dataset, c["ana_species_profiles"], c["common"],
         c["species_of_label"], cfg["VALUE_COL"], cfg["ZSCORE"], pb)
    _try("S7c", fsp.f5_region_invariance, ss["invariance"], ss["common"], ss["region_mask"], pc)
    base = pst.out_base("supp7", paper_dir)
    ps.compose_images([_png(pa), _png(pb), _png(pc)], base, layout=[[0, 1], [2]])
    return base


def supp8(res, inp, paper_dir=None):
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import fingerprints as ff
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import supp_panels as fsp
    c, d = inp.fmri(), res["draw"]
    pa, pb = (_panel("supp8", n, paper_dir) for n in
              ("a_F3_score_distributions", "b_6c_arousal_state_F8b"))
    _try("S8a", fsp.f3_score_distributions, c["sel"]["rankings"], c["selections"], pa)
    if not d["d"].empty:
        _try("S8b", ff.plot_f8b_by_state, d["d"], d["tests"], pb + ".png", paper=True)
    base = pst.out_base("supp8", paper_dir)
    ps.compose_images([_png(pa), _png(pb)], base, layout=[[0], [1]])
    return base


def supp9(res, inp, paper_dir=None):
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import supp_panels as fsp
    c = inp.fmri()
    p = _panel("supp9", "F16b_proximity_rows", paper_dir)
    if c["D"] is not None:
        _try("S9", fsp.f16b_proximity_rows, c["D"], c["dist_all"], c["species_of_label"],
             c["sibling"], p)
    base = pst.out_base("supp9", paper_dir)
    ps.compose_images([_png(p)], base, letters=" ")
    return base


def supp10(res, inp, paper_dir=None):
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import supp_panels as fsp
    c, d = inp.fmri(), res["draw"]
    p = _panel("supp10", "F20_variability_levels", paper_dir)
    if d and not d["vt"].empty:
        _try("S10", fsp.f20_variability_levels, d["vt"], d["vdraws"], c["species_of_label"], p)
    base = pst.out_base("supp10", paper_dir)
    ps.compose_images([_png(p)], base, letters=" ")
    return base


def supp11(res, inp, paper_dir=None):
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import supp_panels as fsp
    c, d = inp.fmri(), res["draw"]
    pa, pb = (_panel("supp11", n, paper_dir) for n in
              ("a_F9_run_distance_sd", "b_F10_qc_vs_distance"))
    if d.get("dist") is not None and not d["dist"].empty:
        _try("S11a", fsp.f9_run_distance_sd, d["dist"], pa)
    if d.get("meta") is not None and not d["meta"].empty:
        _try("S11b", fsp.f10_qc_vs_distance, d["meta"], pb, alpha=c["cfg"]["QC_META_ALPHA"])
    base = pst.out_base("supp11", paper_dir)
    ps.compose_images([_png(pa), _png(pb)], base, layout=[[0, 1]])
    return base


def run(results, inputs, paper_dir=None):
    """Draw every figure of `results` (output of paper_stats.run)."""
    out = {}
    for s in pst.SECTIONS:
        if s not in results:
            continue
        print(f"\n{'=' * 70}\n  figure: {pst.FOLDER[s]}\n{'=' * 70}")
        out[s] = _try(s, globals()[s], results[s], inputs, paper_dir)
    return out
