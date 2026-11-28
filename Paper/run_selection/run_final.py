"""Night run of the fMRI part of the paper for one INVARIANCE_VALUE, with the new Figure 4.

    python Paper/run_selection/run_final.py 0.8          # select -> paper -> fig4
    python Paper/run_selection/run_final.py 0.5 fig4     # one step only
    steps: select | paper | fig4 | all (default)   options: --profiles all|selected

  select  run selection (main + common-target control) into NEW roots
              <OUTPUT_ROOT of config_ednix>_final<XX>[_control_common_target]
          (XX = 100 x INVARIANCE_VALUE; the reference OUTPUT_ROOT is not touched)
  paper   fMRI sections of the paper on that root -> <FIGURES_ROOT>/Paper_final<XX>/
  fig4    new Figure 4, harmonisation figure and supplements (Paper/analysis/fig4_maps)
              -> <FIGURES_ROOT>/Paper_final<XX>/Figure4_maps/

Only the two configuration dicts of config_ednix are changed (INVARIANCE_VALUE and
OUTPUT_ROOT), before any module reads them; nothing else differs from launcher.py.
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from Paper.run_selection import config_ednix as cfg  # noqa: E402

PAPER_SECTIONS = ("figure4", "figure5", "figure6", "supp10", "supp11", "supp12", "supp13",
                  "supp14", "supp15", "supp16")
BASE_ROOT = cfg.SELECTION_CONFIG["OUTPUT_ROOT"]


def configure(inv):
    """Point the selection and the paper at the root of this invariance value."""
    if not 0 < inv <= 1:
        raise SystemExit(f"INVARIANCE_VALUE must be in (0, 1], got {inv}")
    tag = f"final{int(round(100 * inv))}"
    root = f"{BASE_ROOT}_{tag}"
    cfg.SELECTION_CONFIG.update(INVARIANCE_VALUE=inv, OUTPUT_ROOT=root)
    cfg.CONTROL_SELECTION_CONFIG.update(INVARIANCE_VALUE=inv,
                                        OUTPUT_ROOT=root + "_control_common_target")
    from Paper.analysis import paper_data as pdat
    pdat.FMRI_SELECTION_ROOT = root
    paper_dir = os.path.join(pdat.FIGURES_ROOT, f"Paper_{tag}")
    pdat.PAPER_DIR = paper_dir
    return tag, root, paper_dir


def select(root):
    from Statistics.Group_fMRI.run_selection import profiles as sri
    from Statistics.Group_fMRI.run_selection.selection import run_selection
    qc = cfg.load_qc_by_label()
    run_selection(sri, cfg.BIDS_CONFIGS, config=cfg.SELECTION_CONFIG, df_qc_by_label=qc)
    regions = cfg.read_common_regions(root)
    if regions is None:
        print("  [CONTROL] no common_regions.xlsx from the main run -- control skipped")
        return
    run_selection(sri, cfg.BIDS_CONFIGS,
                  config=dict(cfg.CONTROL_SELECTION_CONFIG, COMMON_REGIONS_OVERRIDE=regions),
                  df_qc_by_label=qc)


def paper(root, paper_dir):
    from Paper.analysis import run_paper
    run_paper.main(PAPER_SECTIONS, fmri_root=root, paper_dir=paper_dir)


def fig4(root, paper_dir, pool):
    from Paper.analysis import paper_stats as pst
    from Paper.analysis.fig4_maps import run as f4run
    c = pst.Inputs(fmri_root=root).fmri()
    f4run.run(c, os.path.join(paper_dir, "Figure4_maps"), pool=pool,
              title=os.path.basename(root))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("invariance", type=float)
    ap.add_argument("step", nargs="?", default="all", choices=("select", "paper", "fig4", "all"))
    ap.add_argument("--profiles", choices=("all", "selected"), default="all",
                    help="Fig. 4b/c from the all-runs maps (default) or the selected runs")
    a = ap.parse_args(argv)
    tag, root, paper_dir = configure(a.invariance)
    print(f"INVARIANCE_VALUE = {a.invariance}; selection {root}; paper {paper_dir}; "
          f"step {a.step}", flush=True)
    steps = (("select", lambda: select(root)), ("paper", lambda: paper(root, paper_dir)),
             ("fig4", lambda: fig4(root, paper_dir, a.profiles)))
    for name, fn in steps:
        if a.step not in (name, "all"):
            continue
        if name != "select" and not os.path.isdir(os.path.join(root, "07_state")):
            raise SystemExit(f"{root}: no 07_state -- run the 'select' step first")
        t0 = time.time()
        print(f"\n=== {name} ===", flush=True)
        fn()
        print(f"=== {name} done ({(time.time() - t0) / 3600:.1f} h)", flush=True)


if __name__ == "__main__":
    main()
