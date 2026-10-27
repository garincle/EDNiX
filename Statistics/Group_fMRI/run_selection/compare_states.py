"""
Is a new selection identical to an earlier one?

    python compare_states.py                      # newest archived state vs current state
    python compare_states.py OLD NEW              # two selection_state.pkl, or two OUTPUT_ROOTs

Without arguments: OUTPUT_ROOT of config_ednix.SELECTION_CONFIG; OLD = the most recent
<OUTPUT_ROOT>/_previous_runs/<date>/07_state/selection_state.pkl (every new run moves the
previous state there), NEW = <OUTPUT_ROOT>/07_state/selection_state.pkl.

Every entry of the two states is compared (tables cell by cell, arrays element-wise,
exact equality; floats with rtol=1e-12 to absorb the last-bit noise of a BLAS change),
except 'saved' / 'versions' (time stamps) and 'paths' (folders). In 'cfg' the figure
settings (selection.FIGURE_KEYS) are reported but not counted: they do not act on the
selection. Exit code 0 = same selection.
"""
import os
import sys
import glob
import pickle

import numpy as np
import pandas as pd

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

IGNORED = {"saved", "versions", "paths"}
RTOL = 1e-12


def _state_file(p):
    if os.path.isdir(p):
        p = os.path.join(p, "07_state", "selection_state.pkl")
    return p


def latest_archived(root):
    hits = sorted(glob.glob(os.path.join(root, "_previous_runs", "*", "07_state",
                                         "selection_state.pkl")))
    return hits[-1] if hits else None


def _same(a, b):
    if isinstance(a, pd.DataFrame) or isinstance(b, pd.DataFrame):
        if not (isinstance(a, pd.DataFrame) and isinstance(b, pd.DataFrame)):
            return False
        if list(a.columns) != list(b.columns) or len(a) != len(b):
            return False
        try:
            pd.testing.assert_frame_equal(a, b, check_exact=False, rtol=RTOL, atol=0)
            return True
        except AssertionError:
            return False
    if isinstance(a, pd.Series) or isinstance(b, pd.Series):
        try:
            pd.testing.assert_series_equal(a, b, check_exact=False, rtol=RTOL, atol=0)
            return True
        except (AssertionError, TypeError):
            return False
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        a, b = np.asarray(a), np.asarray(b)
        if a.shape != b.shape:
            return False
        if a.dtype.kind in "fc" or b.dtype.kind in "fc":
            return bool(np.allclose(a, b, rtol=RTOL, atol=0, equal_nan=True))
        return bool(np.array_equal(a, b))
    if isinstance(a, float) and isinstance(b, float):
        return (np.isnan(a) and np.isnan(b)) or np.isclose(a, b, rtol=RTOL, atol=0)
    if isinstance(a, dict) and isinstance(b, dict):
        return set(a) == set(b) and all(_same(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    try:
        return bool(a == b)
    except Exception:
        return repr(a) == repr(b)


def _diff(a, b, path, out):
    """Paths of the differing leaves (recurses into dicts only)."""
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            p = f"{path}[{k!r}]"
            if k not in a or k not in b:
                out.append(f"{p}: only in {'NEW' if k not in a else 'OLD'}")
            else:
                _diff(a[k], b[k], p, out)
    elif not _same(a, b):
        out.append(f"{path}: differs")
    return out


def compare(old_path, new_path):
    from Statistics.Group_fMRI.run_selection.selection import FIGURE_KEYS
    with open(old_path, "rb") as f:
        old = pickle.load(f)
    with open(new_path, "rb") as f:
        new = pickle.load(f)
    print(f"OLD {old_path}  (saved {old.get('saved')})\nNEW {new_path}  (saved {new.get('saved')})")
    diffs = []
    for k in sorted(set(old) | set(new)):
        if k in IGNORED:
            continue
        if k == "cfg":
            for c in sorted(set(old["cfg"]) | set(new["cfg"])):
                if not _same(old["cfg"].get(c), new["cfg"].get(c)):
                    line = f"cfg[{c!r}]: {old['cfg'].get(c)!r} -> {new['cfg'].get(c)!r}"
                    if c in FIGURE_KEYS:
                        print(f"  (figure setting, not counted) {line}")
                    else:
                        diffs.append(line)
            continue
        if k not in old or k not in new:
            diffs.append(f"[{k!r}]: only in {'NEW' if k not in old else 'OLD'}")
            continue
        _diff(old[k], new[k], f"[{k!r}]", diffs)
    kept_old = {l: len(v) for l, v in old.get("selections_final", {}).items()}
    kept_new = {l: len(v) for l, v in new.get("selections_final", {}).items()}
    print(f"kept per dataset, OLD: {kept_old}\nkept per dataset, NEW: {kept_new}")
    if diffs:
        print(f"\n{len(diffs)} DIFFERENCE(S):\n  " + "\n  ".join(diffs[:200]))
    else:
        print("\nSAME SELECTION: every entry of the two states is identical.")
    return diffs


if __name__ == "__main__":
    if len(sys.argv) == 3:
        old_p, new_p = _state_file(sys.argv[1]), _state_file(sys.argv[2])
    else:
        from Statistics.Group_fMRI.run_selection import config_ednix as cfg
        root = cfg.SELECTION_CONFIG["OUTPUT_ROOT"]
        old_p, new_p = latest_archived(root), _state_file(root)
        if old_p is None:
            sys.exit(f"no archived state under {root}/_previous_runs")
    sys.exit(1 if compare(old_p, new_p) else 0)
