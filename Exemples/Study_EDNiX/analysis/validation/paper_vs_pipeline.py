"""
Are the numbers and panels of paper_stats / paper_figures those of the former Fig5
pipeline? Synthetic check, no AFNI.

    python paper_vs_pipeline.py <old_repo> <synthetic_root> <work_dir>

  <old_repo>        checkout with the former pipeline (Exemples/Study_EDNiX/analysis/Fig5,
                    e.g. git worktree add /tmp/ednix_old 0707c94)
  <synthetic_root>  output of Statistics/Group_fMRI/run_selection/validation/
                    synthetic_pipeline.py (its 07_state is read)
  <work_dir>        scratch folder (emptied)

1. The synthetic selection has no AFNI group map. Group maps are simulated as the median
   profile of the units of each group (all, kept / not kept at round 1 and at the final
   round), every surface render as a small random PNG, and two datasets get both arousal
   states in the same animals (run 1 awake, run 2 anaesthetised), so that the paper settings
   (fmri_fingerprints.settings.FIGURE_SETTINGS: kept maps of the final round for
   Figure 4, distance to the kept map for Supp. 12) run as on real data.
2. The former pipeline (part 2, "figures", from that state) runs in <old_repo>.
3. paper_stats + paper_figures run in this checkout on the same state.
4. Every table, result line and print-size panel of the paper is compared.
Exit code 0 = identical.
"""
import os
import sys
import pickle
import shutil
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, REPO)
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd

opj = os.path.join
GROUPS = ("all", "kept", "not_kept", "kept_final", "not_kept_final")
MIXED_STATE = ("Mouse", "Marmoset_Hori")


# ═══════════════════════════════════════════════════════════════════════════
def augment(src_root, dst_root):
    """Synthetic state + simulated group maps, written to dst_root/07_state."""
    from Statistics.Group_fMRI.run_selection import profiles as sri
    from Statistics.Group_fMRI.run_selection.selection import load_state, save_state
    st = load_state(src_root)
    cfg = dict(st["cfg"], OUTPUT_ROOT=dst_root)
    sri.set_profile_method(cfg["PROFILE_METHOD"], cfg["TAIL_Q"], keep_sketch=False)
    if cfg["PROFILE_METHOD"] != "median_r":
        sri.set_tail_reference(st["anat_common"])

    def key(k):
        return tuple(str(x) for x in k)

    sol = st["species_of_label"]
    gmp = {g: {} for g in GROUPS}
    rows_all, rows_r1, rows_fin = [], [], []
    pics, ranges = {}, {}
    pic_dir = opj(dst_root, "05b_surfaces_final")
    os.makedirs(pic_dir, exist_ok=True)
    rng = np.random.default_rng(3)
    import matplotlib.pyplot as plt
    for lab, profs in sorted(st["profiles_by_label"].items()):
        k1 = {key(k) for k in st["selections"].get(lab, [])}
        kf = {key(k) for k in st["selections_final"].get(lab, [])}
        pools = dict(all=set(profs),
                     kept={k for k in profs if key(k) in k1},
                     not_kept={k for k in profs if key(k) not in k1},
                     kept_final={k for k in profs if key(k) in kf},
                     not_kept_final={k for k in profs if key(k) not in kf})
        for g, keys in pools.items():
            if keys:
                gmp[g][lab] = (sol.get(lab), sri.compute_group_median_profile(profs, keys))
        for g, rows in (("all", rows_all), ("kept", rows_r1), ("not_kept", rows_r1),
                        ("kept_final", rows_fin), ("not_kept_final", rows_fin)):
            rows.append(dict(label=lab, species=sol.get(lab), group=g,
                             round=st["final_round"] if "final" in g else 1,
                             output_suffix=f"synthetic_{g}",
                             map_path=f"/synthetic/{lab}_{g}.nii.gz", found=True))
        for g in ("all", "kept_final", "not_kept_final"):
            p = opj(pic_dir, f"{lab}_{g}.png")
            plt.imsave(p, rng.random((90, 160, 3)))
            pics[(lab, g)] = p
        ranges[lab] = float(rng.uniform(0.5, 2))
    # two datasets with both arousal states in the same animals (as Mouse and
    # Marmoset_Hori): run 1 awake, run 2 anaesthetised, read from BIDS task names
    bids_configs = [dict(b) for b in st["bids_configs"]]
    for b in bids_configs:
        if b["label"] in MIXED_STATE:
            b["bids_dir"] = opj(dst_root, "_bids", b["label"])
            b["state"] = {"from": "task", "map": {"restawake": "awake",
                                                  "restanest": "anesthetized"}}
            for k in st["runs_by_label"][b["label"]]:
                subj, ses, run = (str(x) for x in k[:3])
                d = opj(b["bids_dir"], f"sub-{subj}", f"ses-{ses}", "func")
                os.makedirs(d, exist_ok=True)
                task = "restawake" if run == "1" else "restanest"
                open(opj(d, f"sub-{subj}_ses-{ses}_task-{task}_run-{run}_bold.nii.gz"),
                     "w").close()
    st = dict(st, cfg=cfg, afni_profiles=dict(gmp["all"]), group_map_profiles=gmp,
              bids_configs=bids_configs,
              all_units_index=pd.DataFrame(rows_all), map_index=pd.DataFrame(rows_r1),
              map_index_final=pd.DataFrame(rows_fin), surface_pics=pics,
              surface_ranges=ranges)
    save_state(st)
    return st


OLD_RUNNER = r'''
import os, sys
sys.path.insert(0, sys.argv[1])
os.environ["MPLBACKEND"] = "Agg"
import json
from Exemples.Study_EDNiX.analysis.Fig5 import fig5_pipeline as pipe
from Exemples.Study_EDNiX.analysis.Fig5 import SBA_Region_Profile_Index as sri
settings = json.loads(sys.argv[3])
pipe.run(sri, None, config=dict(OUTPUT_ROOT=sys.argv[2], **settings), part="figures")
'''


def settings_for_validation():
    """The paper settings, with the brain-mass table off (not on this machine; not used
    by the paper statistics) and the full-tree comparative tests off (no tree)."""
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints.settings import FIGURE_SETTINGS
    return dict(FIGURE_SETTINGS, BRAIN_MASS_XLSX=None, RUN_PHYLO_COMPARATIVE=False)


# ═══════════════════════════════════════════════════════════════════════════
def _frame_equal(a, b, name, rtol=1e-12):
    a, b = a.reset_index(drop=True), b.reset_index(drop=True)
    if list(a.columns) != list(b.columns):
        return f"{name}: columns differ {list(a.columns)} vs {list(b.columns)}"
    if a.shape != b.shape:
        return f"{name}: shape {a.shape} vs {b.shape}"
    for c in a.columns:
        x, y = a[c], b[c]
        xn, yn = pd.to_numeric(x, errors="coerce"), pd.to_numeric(y, errors="coerce")
        if xn.notna().any() or yn.notna().any():
            if not np.allclose(xn.to_numpy(float), yn.to_numpy(float), rtol=rtol, atol=0,
                               equal_nan=True):
                return f"{name}: column {c} differs"
            if (xn.isna() != yn.isna()).any():
                return f"{name}: column {c} differs (missing values)"
        if (xn.isna() & yn.isna()).any():
            xs = x[xn.isna()].astype(object).where(x[xn.isna()].notna(), "<NA>").astype(str)
            ys = y[yn.isna()].astype(object).where(y[yn.isna()].notna(), "<NA>").astype(str)
            if list(xs) != list(ys):
                return f"{name}: column {c} differs (text)"
    return None


def _xl(path, sheet=0, **kw):
    return pd.read_excel(path, sheet_name=sheet, **kw)


def _png_equal(a, b):
    import matplotlib.image as mpimg
    if not (os.path.isfile(a) and os.path.isfile(b)):
        return False, "missing " + (a if not os.path.isfile(a) else b)
    x, y = mpimg.imread(a), mpimg.imread(b)
    if x.shape != y.shape:
        return False, f"size {x.shape} vs {y.shape}"
    d = float(np.abs(x.astype(float) - y.astype(float)).max())
    return d == 0, f"max pixel difference {d:.3g}"


def compare(old_root, new_paper):
    """[(check, ok, detail)]"""
    gm, fig = opj(old_root, "03_group_maps"), opj(old_root, "06_figures")
    out = []

    def table(check, new_xlsx, new_sheet, old_xlsx, old_sheet=0, index_col=None):
        try:
            n = _xl(new_xlsx, new_sheet, index_col=index_col)
            o = _xl(old_xlsx, old_sheet, index_col=index_col)
            if index_col is not None:
                n, o = n.reset_index(), o.reset_index()
            err = _frame_equal(n, o, check)
            out.append((check, err is None, err or f"{len(n)} rows identical"))
        except Exception as e:
            out.append((check, False, f"{type(e).__name__}: {e}"))

    def lines(check, new_txt, old_txt, marker, after=None):
        """result lines of the new text starting with `marker` (after the line starting
        with `after`, when given) must all be in the former text"""
        try:
            txt = [ln.strip() for ln in open(new_txt, encoding="utf-8")]
            if after is not None:
                txt = txt[next(i for i, ln in enumerate(txt) if ln.startswith(after)) + 1:]
            new = [ln for ln in txt if ln.startswith(marker)]
            old = open(old_txt, encoding="utf-8").read()
            miss = [ln for ln in new if ln not in old]
            out.append((check, bool(new) and not miss,
                        f"{len(new)} result lines found in the former text" if new and not miss
                        else f"not found: {miss[:2]}" if miss else "no line"))
        except Exception as e:
            out.append((check, False, f"{type(e).__name__}: {e}"))

    def panel(check, new_png, old_png):
        ok, det = _png_equal(new_png, old_png)
        out.append((check, ok, det))

    P = lambda f: opj(new_paper, f, f + "_stats.xlsx")
    T = lambda f: opj(new_paper, f, f + "_stats.txt")
    N = lambda f, n: opj(new_paper, f, "panels", n + ".png")
    fp = opj(gm, "fingerprint_tables.xlsx")
    # Figure 4
    table("Fig 4a distances", P("Figure_4"), "a_distances", fp, "distances_all")
    lines("Fig 4b statistics", T("Figure_4"), opj(fig, "Fig5", "5b_intra_vs_inter_species.txt"),
          ("all runs:", "selected:"))
    table("Fig 4c distance matrix", P("Figure_4"), "c_distance_matrix", fp,
          "dataset_distance_matrix", index_col=0)
    table("Fig 4c sibling ranks", P("Figure_4"), "c_sibling_ranks", fp, "sibling_ranks")
    table("Fig 4c sibling test", P("Figure_4"), "c_sibling_test", fp, "sibling_test")
    panel("Fig 4a panel", N("Figure_4", "a_F15b_fingerprint_vs_human"),
          opj(fig, "Paper", "F15b_fingerprint_vs_human.png"))
    panel("Fig 4b panel", N("Figure_4", "b_5b_intra_vs_inter_species"),
          opj(fig, "Fig5", "5b_intra_vs_inter_species.png"))
    panel("Fig 4c panel", N("Figure_4", "c_F16a_proximity_matrix"),
          opj(fig, "Paper", "F16a_proximity_matrix.png"))
    # Figure 5
    panel("Fig 5 panel", N("Figure_5", "5a_surface_maps_selected_vs_excluded"),
          opj(fig, "Fig5", "5a_surface_maps_selected_vs_excluded.png"))
    # Supp 8
    for n in ("F1e_convergence", "F1c_consensus_by_dataset", "F5_region_invariance"):
        panel(f"Supp 8 {n}", N("Supp_Fig_08", {"F1e_convergence": "a_", "F5_region_invariance":
                                               "c_"}.get(n, "b_") + n),
              opj(fig, "Paper", n + ".png"))
    # Supp 9
    table("Supp 9b run r by state", P("Supp_Fig_09"), "b_run_r_by_state",
          opj(gm, "run_r_to_selection_criterion_by_state.xlsx"))
    lines("Supp 9b statistics", T("Supp_Fig_09"), opj(fig, "Fig6", "6c_arousal_state_F8b.txt"),
          ("Mouse", "Marmoset", "Rat", "Dog", "Macaque", "Human", "Mouselemur", "POOLED"),
          after="b  ")
    panel("Supp 9a panel", N("Supp_Fig_09", "a_F3_score_distributions"),
          opj(fig, "Paper", "F3_score_distributions.png"))
    panel("Supp 9b panel", N("Supp_Fig_09", "b_6c_arousal_state_F8b"),
          opj(fig, "Fig6", "6c_arousal_state_F8b.png"))
    # Supp 10
    table("Supp 10 sibling ranks", P("Supp_Fig_10"), "sibling_ranks", fp, "sibling_ranks")
    panel("Supp 10 panel", N("Supp_Fig_10", "F16b_proximity_rows"),
          opj(fig, "Paper", "F16b_proximity_rows.png"))
    # Supp 11
    table("Supp 11 variability levels", P("Supp_Fig_11"), "variability_levels", fp,
          "variability_levels")
    lines("Supp 11 statistics", T("Supp_Fig_11"),
          opj(fig, "supplementary", "F20_selection_vs_variability_levels.txt"),
          ("Dog", "Macaque", "Mouse", "Rat", "Marmoset", "Human", "Mouselemur", "POOLED"))
    panel("Supp 11 panel", N("Supp_Fig_11", "F20_variability_levels"),
          opj(fig, "Paper", "F20_variability_levels.png"))
    # Supp 12
    table("Supp 12a run r to kept map", P("Supp_Fig_12"), "a_run_r_to_kept_map",
          opj(gm, "run_distance_to_reference.xlsx"))
    table("Supp 12b QC correlations", P("Supp_Fig_12"), "b_qc_correlation_per_dataset",
          opj(gm, "qc_vs_distance_correlation_per_bids.xlsx"))
    table("Supp 12b QC meta-analysis", P("Supp_Fig_12"), "b_qc_meta_analysis",
          opj(gm, "qc_vs_distance_meta_analysis.xlsx"))
    panel("Supp 12a panel", N("Supp_Fig_12", "a_F9_run_distance_sd"),
          opj(fig, "Paper", "F9_run_distance_sd.png"))
    panel("Supp 12b panel", N("Supp_Fig_12", "b_F10_qc_vs_distance"),
          opj(fig, "Paper", "F10_qc_vs_distance.png"))
    return out


def main(old_repo, synthetic_root, work):
    if os.path.isdir(work):
        shutil.rmtree(work)
    root_old, root_new = opj(work, "pipeline_old"), opj(work, "paper_new")
    augment(synthetic_root, root_old)
    shutil.copytree(root_old, root_new)
    import json
    settings = settings_for_validation()
    runner = opj(work, "old_runner.py")
    with open(runner, "w") as f:
        f.write(OLD_RUNNER)
    print("  former pipeline (part 2) ...")
    r = subprocess.run([sys.executable, runner, old_repo, root_old, json.dumps(settings)],
                       capture_output=True, text=True)
    with open(opj(work, "old_pipeline.log"), "w") as f:
        f.write(r.stdout + r.stderr)
    if r.returncode:
        print(r.stderr[-3000:])
        raise SystemExit("former pipeline failed")
    # the paper on the same state (part-2 settings of the validation)
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import context
    from Exemples.Study_EDNiX.analysis import paper_stats as pst
    from Exemples.Study_EDNiX.analysis import paper_figures as pfg
    paper = opj(work, "Paper")
    inp = pst.Inputs(fmri_root=root_new)
    inp._cache["fmri"] = None
    c = context.prepare(context.load_state(root_new, settings=settings))
    c.update(pst.fmri_group(c))
    inp._cache["fmri"] = c
    sections = ("figure4", "figure5", "supp8", "supp9", "supp10", "supp11", "supp12")
    res = pst.run(sections, inp, paper_dir=paper)
    pfg.run(res, inp, paper_dir=paper)
    checks = compare(root_old, paper)
    print(f"\n{'=' * 70}")
    for name, ok, det in checks:
        print(f"  [{'SAME' if ok else 'DIFF'}] {name}: {det}")
    n_ok = sum(ok for _, ok, _ in checks)
    print(f"\n  {n_ok}/{len(checks)} identical")
    return 0 if n_ok == len(checks) else 1


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:4]))
