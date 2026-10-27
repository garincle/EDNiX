"""Compare every table (.xlsx, .csv) and text (.txt) under two output roots."""
import os, sys, glob
import pandas as pd
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..')))
from Statistics.Group_fMRI.run_selection.compare_states import _same

def main():
    """Run from the command line only (the test suite imports every module)."""
    a, b = sys.argv[1], sys.argv[2]
    skip = ("00_log", "_previous_runs", "selection_state.txt", "_inputs")
    def files(root):
        out = set()
        for ext in ("xlsx", "csv", "txt"):
            for p in glob.glob(os.path.join(root, "**", f"*.{ext}"), recursive=True):
                r = os.path.relpath(p, root)
                if not any(s in r for s in skip):
                    out.add(r)
        return out
    fa_, fb_ = files(a), files(b)
    print(f"{len(fa_)} files OLD, {len(fb_)} NEW")
    for r in sorted(fa_ ^ fb_):
        print("  only in", "OLD" if r in fa_ else "NEW", r)
    nd = 0
    for r in sorted(fa_ & fb_):
        pa, pb = os.path.join(a, r), os.path.join(b, r)
        if r.endswith(".txt"):
            ok = open(pa, errors="ignore").read() == open(pb, errors="ignore").read()
        elif r.endswith(".csv"):
            ok = _same(pd.read_csv(pa), pd.read_csv(pb))
        else:
            xa, xb = pd.read_excel(pa, sheet_name=None), pd.read_excel(pb, sheet_name=None)
            ok = set(xa) == set(xb) and all(_same(xa[s], xb[s]) for s in xa)
        if not ok:
            nd += 1
            print("  DIFFERS", r)
    print(f"{len(fa_ & fb_) - nd}/{len(fa_ & fb_)} common files identical")
    sys.exit(1 if nd or (fa_ ^ fb_) else 0)


if __name__ == "__main__":
    main()
