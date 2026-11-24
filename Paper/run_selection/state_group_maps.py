"""Group seed maps of each arousal state, inside the datasets that hold both states
(Supplementary figure 'awake vs anaesthetised', paper section supp18).

    python Paper/run_selection/state_group_maps.py                 # kept runs of the final round
    python Paper/run_selection/state_group_maps.py --pool all      # every run of each state
    python Paper/run_selection/state_group_maps.py --plan           # list the runs, run nothing

For every dataset whose runs are both awake and anaesthetised (state rule of
config_ednix.BIDS_CONFIGS; every state other than "awake" counts as anaesthetised, as in
Supp. 15), the runs of the pool are split by state and each half goes through the SAME
AFNI 3dttest++ as the kept / not-kept maps of Figure 5 (profiles.rerun_balanced_by_run,
balance=False: every run of the state, no subsampling; the map drawn is the mean Fisher z,
sub-brick 0, so its scale does not grow with n). The two maps of a dataset are then
rendered on its inflated surface with ONE symmetric colour range (render_common_scale,
+/- P<pct> of |z| over the two maps) so that the states are directly comparable.

Needs the saved selection (<OUTPUT_ROOT>/07_state, Paper/run_selection/launcher.py), AFNI
and Workbench (Singularity images of group_maps). Cost: one 3dttest++ per state and
dataset (2 datasets x 2 states); 3dClustSim only for n > 14 runs.

Output: <OUTPUT_ROOT>/08_state_maps/
    state_map_index.xlsx    label, state, n_runs, n_animals, map_path, found, png, range
    state_runs.xlsx         every run of the pool with its state and animal
    state_maps.json         same index, read by paper_stats.supp18
    surfaces/<label>__<state>/   wb_command renders
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from Paper.run_selection import config_ednix  # noqa: E402  (configures group_maps)

STATES = ("awake", "anesthetized")
SUFFIX_CORE = {"awake": "StateAwake", "anesthetized": "StateAnesth"}
OUT_FOLDER = "08_state_maps"


def output_dir(root):
    return os.path.join(root, OUT_FOLDER)


def pool_runs(st, pool="kept"):
    """One row per run of the datasets holding both states: label, species, subject,
    session, run, seed, animal, state (raw rule value), state_class (awake /
    anesthetized), kept (final round). pool='kept': kept runs only."""
    from Paper.analysis.fmri_fingerprints import fingerprints as ff
    from Paper.analysis.fmri_fingerprints import paper_panels as fpp
    cfg_by = {c["label"]: c for c in st["bids_configs"]}
    sol = st["species_of_label"]
    kept_runs = fpp.selected_run_keys(st["selections_final"], st["members_by_label"])
    rows = []
    for lab, runs in st["runs_by_label"].items():
        res = ff.StateResolver(cfg_by.get(lab, {}))
        if ff.dataset_state(res, list(runs)) != "mixed":
            continue
        for k in runs:
            key = tuple(str(x) for x in k)
            state = res.state(key)
            rows.append(dict(label=lab, species=sol.get(lab), subject=key[0], session=key[1],
                             run=key[2], seed=key[3] if len(key) > 3 else "",
                             animal=res.animal(key), state=state,
                             state_class=ff.state_class(state),
                             kept=key in kept_runs.get(lab, set())))
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df = df[df["state_class"].isin(STATES)]
    return df[df["kept"]].reset_index(drop=True) if pool == "kept" else df.reset_index(drop=True)


def _voxel_r(path_a, path_b):
    """Pearson r between two group maps over the voxels non-zero and finite in both."""
    import nibabel as nib
    a = np.asarray(nib.load(path_a).get_fdata()).ravel()
    b = np.asarray(nib.load(path_b).get_fdata()).ravel()
    if a.shape != b.shape:
        return np.nan, 0
    m = np.isfinite(a) & np.isfinite(b) & (a != 0) & (b != 0)
    if m.sum() < 10:
        return np.nan, int(m.sum())
    return float(np.corrcoef(a[m], b[m])[0, 1]), int(m.sum())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", default=config_ednix.SELECTION_CONFIG["OUTPUT_ROOT"],
                    help="selection OUTPUT_ROOT (07_state inside)")
    ap.add_argument("--pool", choices=("kept", "all"), default="kept",
                    help="runs kept by the final round (default) or every run")
    ap.add_argument("--min-runs", type=int, default=3,
                    help="no map for a state with fewer runs (AFNI needs >= 2; default 3)")
    ap.add_argument("--pct", type=float, default=90.0, help="colour range = +/- P<pct> of |z|")
    ap.add_argument("--no-render", action="store_true", help="AFNI maps only, no surface render")
    ap.add_argument("--plan", action="store_true", help="print the runs of each map and stop")
    a = ap.parse_args(argv)

    from Statistics.Group_fMRI.run_selection import profiles as sri
    from Statistics.Group_fMRI.run_selection import surfaces
    from Statistics.Group_fMRI.run_selection.selection import (load_state, _map_suffix,
                                                               _resolve_map)
    st = load_state(a.root)
    cfg = st["cfg"]
    out = output_dir(a.root)
    os.makedirs(out, exist_ok=True)
    df = pool_runs(st, a.pool)
    if df.empty:
        print("  no dataset holds both states -- nothing to do")
        return 1
    df.to_excel(os.path.join(out, "state_runs.xlsx"), index=False)
    cfg_by = {c["label"]: c for c in st["bids_configs"]}
    print(f"  pool = {a.pool} (final round {st['final_round']}), datasets "
          f"{sorted(df['label'].unique())}")
    for (lab, sc), g in df.groupby(["label", "state_class"]):
        print(f"    {lab:16s} {sc:13s} {len(g):3d} runs, {g['animal'].nunique():3d} animals "
              f"({', '.join(f'{k} {v}' for k, v in g['state'].value_counts().items())})")
    if a.plan:
        return 0

    rows, maps = [], {}
    for lab, g in df.groupby("label"):
        cb = cfg_by[lab]
        for state in STATES:
            gs = g[g["state_class"] == state]
            core = f"{SUFFIX_CORE[state]}{a.pool.capitalize()}R{st['final_round']}"
            suffix = _map_suffix(cfg, core)
            row = dict(label=lab, species=cb["species"], state=state, pool=a.pool,
                       n_runs=len(gs), n_animals=int(gs["animal"].nunique()),
                       n_sessions=int(gs[["subject", "session"]].drop_duplicates().shape[0]),
                       suffix=suffix, map_path="", found=False)
            if len(gs) < a.min_runs:
                print(f"  [SKIP] {lab}/{state}: {len(gs)} run(s) < --min-runs {a.min_runs}")
                rows.append(row)
                continue
            df_class = gs[["subject", "session", "run", "seed"]].assign(is_outlier=False)
            print(f"\n  {lab} / {state}: 3dttest++ on {len(gs)} run(s) ({suffix})")
            try:
                sri.rerun_balanced_by_run(
                    df_class, cb["species"], cb["bids_dir"], seed=cfg["SEED"],
                    results_base_dir=cfg["RESULTS_BASE_DIR"], output_suffix_kept=suffix,
                    output_suffix_outlier=suffix + "_unused", balance=False, unit=cfg["UNIT"],
                    exclude_lateralized=cb.get("exclude_lateralized", True),
                    render_surface=False)
            except Exception as e:
                print(f"  [ERROR] {lab}/{state}: 3dttest++ failed ({e})")
                rows.append(row)
                continue
            path, found = _resolve_map(sri, cfg, cb, suffix)
            row.update(map_path=path or "", found=bool(found))
            if found:
                maps.setdefault(lab, {})[state] = path
            rows.append(row)

    pics, ranges = {}, {}
    if maps and not a.no_render:
        pics, ranges = surfaces.render_common_scale(
            sri, os.path.join(out, "surfaces"), maps, st["species_of_label"], pct=a.pct)
    idx = pd.DataFrame(rows)
    idx["png"] = [pics.get((r.label, r.state), "") for r in idx.itertuples()]
    idx["colour_range"] = [ranges.get(r.label, np.nan) for r in idx.itertuples()]
    vr = {}
    for lab, m in maps.items():
        if all(s in m for s in STATES):
            vr[lab] = _voxel_r(surfaces._oversampled(m[STATES[0]]),
                               surfaces._oversampled(m[STATES[1]]))
    idx["voxel_r_between_states"] = [vr.get(r.label, (np.nan, 0))[0] for r in idx.itertuples()]
    idx["voxel_n"] = [vr.get(r.label, (np.nan, 0))[1] for r in idx.itertuples()]
    idx.to_excel(os.path.join(out, "state_map_index.xlsx"), index=False)
    with open(os.path.join(out, "state_maps.json"), "w") as f:
        json.dump(dict(pool=a.pool, final_round=st["final_round"], pct=a.pct,
                       rows=json.loads(idx.to_json(orient="records"))), f, indent=2)
    print(f"\n  [saved] {out}/state_map_index.xlsx, state_maps.json")
    print(idx[["label", "state", "n_runs", "n_animals", "found", "colour_range",
               "voxel_r_between_states"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
