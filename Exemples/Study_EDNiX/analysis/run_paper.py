"""
run_paper.py -- statistics then figures of the EDNiX paper.

    python run_paper.py                    # every figure
    python run_paper.py figure4 supp10     # some figures
    python run_paper.py --stats-only       # numbers only, no drawing
    python run_paper.py --refresh          # re-read the EDNiX outputs (no cache)

Before: the run selection of Figures 4-5 / Supp. 8-12 must exist
(python Statistics/Group_fMRI/run_selection/launcher.py, writes its 07_state).
Output: <PAPER_DIR>/<Figure>/<Figure>.pdf|png, <Figure>_stats.xlsx|txt, panels/.
"""
import sys

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_stats as pst
from Exemples.Study_EDNiX.analysis import paper_figures as pfg


def main(sections=pst.SECTIONS, stats_only=False, refresh=False, fmri_root=None,
         paper_dir=None):
    inputs = pst.Inputs(fmri_root=fmri_root, refresh=refresh)
    results = pst.run(sections, inputs, paper_dir=paper_dir)
    figures = {} if stats_only else pfg.run(results, inputs, paper_dir=paper_dir)
    failed = [s for s in sections if s not in results or (not stats_only and not figures.get(s))]
    print(f"\n{'#' * 70}\n  paper: {len(results)}/{len(sections)} statistics"
          + ("" if stats_only else f", {sum(1 for v in figures.values() if v)} figures")
          + (f"; FAILED: {failed}" if failed else "") + f"\n{'#' * 70}")
    return results, figures


if __name__ == "__main__":
    args = sys.argv[1:]
    flags = {a for a in args if a.startswith("--")}
    unknown = flags - {"--stats-only", "--refresh"}
    if unknown:
        raise SystemExit(f"unknown option(s) {sorted(unknown)}")
    main(tuple(a for a in args if not a.startswith("--")) or pst.SECTIONS,
         stats_only="--stats-only" in flags, refresh="--refresh" in flags)
