"""New Figure 4, harmonisation figure and supplements from one finished selection.

    python -m Paper.analysis.fig4_maps.run <OUTPUT_ROOT of the selection> <out_dir>
        [--profiles all|selected]   Fig. 4b/c from the all-runs maps (default, as 4a) or
                                    from the selected (kept) runs
        [--draws 100]               random draws of the 3D-map size control
        [--perm 100]                random partitions of the noise-bin null
        [--amp-draws 1000]          null draws of the linear (amplitude, matrix) analyses
        [--bins 5] [--no-matrices] [--no-surfaces]

Writes into <out_dir>: Figure4.pdf/png (+ panels/), Figure4_harmonisation.pdf/png,
Supp_Figure4_noise_bins, Supp_Figure4_bins_profiles, Supp_Figure4_surfaces,
Figure4_stats.xlsx (one sheet per analysis) and Figure4_stats.txt (every number printed).
"""
import argparse
import os
import sys
import time

import pandas as pd

from Paper.analysis.fig4_maps import analyses as fa4
from Paper.analysis.fig4_maps import figures as ff4


def run(c, out_dir, pool="all", n_draws=100, n_perm=100, amp_draws=1000, n_bins=5,
        do_matrices=True, do_surfaces=True, title=""):
    os.makedirs(out_dir, exist_ok=True)
    log = fa4.Log()
    log(f"=== Figure 4 map analyses {title}\n  4b/c profiles from the {pool}-runs maps; "
        f"{n_draws} draws, {n_perm} partitions, {amp_draws} linear null draws, B = {n_bins}")
    t0 = time.time()
    ctx = fa4.Context(c, n_bins=n_bins, log=log)
    ph = fa4.phylogeny(ctx, n_draws)
    pr = fa4.profiles(ctx, ph, pool)
    nb = fa4.noise_bins(ctx, ph, n_perm)
    am = fa4.amplitude(ctx, amp_draws)
    mx = None
    if do_matrices:
        try:
            mx = fa4.matrices(ctx, amp_draws)
        except Exception as e:                     # matrices are optional (step 10 outputs)
            import traceback
            traceback.print_exc()
            log(f"  [WARN] correlation matrices skipped: {e}")
    sheets = {"phylogeny_3D": ph["table"], "noise_bins": nb["table"],
              "amplitude_seed": am["table"], "amplitude_per_dataset": am["amp"].reset_index()
              .rename(columns={"index": "dataset"}),
              "fig4b_dataset_r": pr["dist"]}
    if mx is not None:
        sheets.update(matrices=mx["table"],
                      matrix_amp_per_dataset=mx["amp"].reset_index().rename(columns={"index": "dataset"}))
    xlsx = os.path.join(out_dir, "Figure4_stats.xlsx")
    with pd.ExcelWriter(xlsx) as xw:
        for name, df in sheets.items():
            df.to_excel(xw, sheet_name=name[:31], index=False)
    log(f"\n[saved] {xlsx}")
    ff4.figure4(ctx, ph, pr, out_dir)
    ff4._try(ff4.harmonisation, ctx, am, mx, out_dir)
    ff4._try(ff4.supp_noise_bins, ctx, ph, nb, out_dir)
    ff4._try(ff4.supp_bins_profiles, ctx, ph, nb, out_dir)
    if do_surfaces:
        ff4._try(ff4.supp_surfaces, ctx, ph, nb, out_dir)
    log(f"\n  done in {(time.time() - t0) / 60:.0f} min -> {out_dir}")
    with open(os.path.join(out_dir, "Figure4_stats.txt"), "w") as fh:
        fh.write(log.text())
    return dict(ctx=ctx, phylogeny=ph, profiles=pr, noise_bins=nb, amplitude=am, matrices=mx)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("root", help="OUTPUT_ROOT of a finished selection (with 07_state)")
    ap.add_argument("out_dir")
    ap.add_argument("--profiles", choices=("all", "selected"), default="all")
    ap.add_argument("--draws", type=int, default=100)
    ap.add_argument("--perm", type=int, default=100)
    ap.add_argument("--amp-draws", type=int, default=1000)
    ap.add_argument("--bins", type=int, default=5)
    ap.add_argument("--no-matrices", action="store_true")
    ap.add_argument("--no-surfaces", action="store_true")
    a = ap.parse_args(argv)
    if not os.path.isdir(os.path.join(a.root, "07_state")):
        raise SystemExit(f"{a.root}: no 07_state -- run the selection first")
    from Paper.analysis import paper_stats as pst
    c = pst.Inputs(fmri_root=a.root).fmri()
    run(c, a.out_dir, pool=a.profiles, n_draws=a.draws, n_perm=a.perm, amp_draws=a.amp_draws,
        n_bins=a.bins, do_matrices=not a.no_matrices, do_surfaces=not a.no_surfaces,
        title=os.path.basename(a.root.rstrip("/")))


if __name__ == "__main__":
    sys.exit(main())
