"""
Supplementary figure: awake vs anaesthetised, inside the datasets that hold both states
(paper section supp16; group maps computed by Paper/run_selection/state_group_maps.py).

One row per dataset (Mouse, Marmoset_Hori):
  a  retrosplenial seed map (mean Fisher z) of the runs of each state, on the dataset's
     inflated surface, ONE symmetric colour range for the two states;
     n runs / n animals under each map
  b  fingerprint (regional profile, Figure 4 radar layout): median run profile of each
     state + the leave-one-species-out (LOSO) consensus
  c  how far apart the states are, compared with how far apart the animals are:
     Pearson r between animal-mean profiles for three kinds of pairs (other animal / same
     state; same animal / other state; other animal / other state), and r between the two
     state profiles of b (diamond) -- and, when the maps exist, the voxelwise r between the
     two maps of a.

Descriptive by design: state is nested in sessions (and, in Mouse, in the subject label)
and the pools are small (the final-round selection keeps 10 runs per dataset), so no
awake - anaesthetised test is drawn here; the within-animal tests are Supp. 13.

Pure functions (no file reading except the PNG renders), so the figure can be drawn on
synthetic inputs (Paper/tests/test_state_maps_figure.py, preview script).
"""
import itertools
import json
import os

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

import Paper.analysis.fmri_fingerprints.fingerprints as ff
import Paper.analysis.fmri_fingerprints.profile_stats as fa
from Paper.analysis import paper_style as ps

AWAKE, ANESTH = ff.AWAKE, ff.ANESTH
STATES = (AWAKE, ANESTH)
STATE_NAME = {AWAKE: "awake", ANESTH: "anaesthetised"}
PAIR_KINDS = ("other animal, same state", "same animal, other state",
              "other animal, other state")
PAIR_COLORS = {PAIR_KINDS[0]: "#888888", PAIR_KINDS[1]: "#7B3294",
               PAIR_KINDS[2]: "#C2A5CF"}
# approximation of Workbench's RY-BC-BL palette, for the colour bar only
RYBCBL = LinearSegmentedColormap.from_list(
    "RYBCBL", [(0.0, "#00FFFF"), (0.25, "#0000FF"), (0.5, "#FFFFFF"),
               (0.75, "#FF0000"), (1.0, "#FFFF00")])


# ═══════════════════════════════════════════════════════════════════════════
# Inputs
# ═══════════════════════════════════════════════════════════════════════════
def runs_from_context(c, pool="kept"):
    """(runs table, {label: {run key: scaled vector}}) for the datasets of `c`
    (paper_stats.Inputs.fmri()) that hold both states; pool 'kept' = runs behind the
    final-round kept units (c['kept_runs_final']), 'all' = every run."""
    cfg = c["cfg"]
    value_col, norm = cfg["VALUE_COL"], cfg["ZSCORE"]
    rows, vecs = [], {}
    for lab, runs in c["runs_by_label"].items():
        res = c["resolvers"][lab]
        if c["dataset_state"].get(lab) != "mixed":
            continue
        kept = {tuple(str(x) for x in k) for k in c["kept_runs_final"].get(lab, set())}
        for k, prof in runs.items():
            key = tuple(str(x) for x in k)
            state = res.state(key)
            sc = ff.state_class(state)
            if sc not in STATES or (pool == "kept" and key not in kept):
                continue
            vecs.setdefault(lab, {})[key] = ff._run_vector(prof, c["common"], value_col, norm)
            rows.append(dict(label=lab, species=c["species_of_label"].get(lab),
                             subject=key[0], session=key[1], run=key[2], key=key,
                             animal=res.animal(key), state=state, state_class=sc,
                             kept=key in kept))
    return pd.DataFrame(rows), vecs


def read_map_index(root):
    """Rows of <root>/08_state_maps/state_maps.json (state_group_maps.py), or None."""
    p = os.path.join(root, "08_state_maps", "state_maps.json")
    if not os.path.isfile(p):
        return None
    with open(p) as f:
        return json.load(f)


# ═══════════════════════════════════════════════════════════════════════════
# Statistics (descriptive)
# ═══════════════════════════════════════════════════════════════════════════
def _median(vs):
    if not vs:
        return None
    with np.errstate(invalid="ignore"):
        return np.nanmedian(np.vstack(vs), axis=0)


def compute(runs, run_vecs, group_vecs, species_of_label, map_rows=None):
    """runs: table of runs_from_context; run_vecs: {label: {key: vector}};
    group_vecs: {label: group vector} of every dataset (Figure 4), for the LOSO
    consensus. Returns dict(summary, pairs, profiles, lines)."""
    maps = {(r["label"], r["state"]): r for r in (map_rows or [])}
    summary, pairs, profiles = [], [], {}
    for lab, g in runs.groupby("label", sort=False):
        sp = g["species"].iloc[0]
        prof = {s: _median([run_vecs[lab][k] for k in g.loc[g["state_class"] == s, "key"]])
                for s in STATES}
        prof["consensus"] = ff.consensus(group_vecs, species_of_label, exclude_species=sp)
        profiles[lab] = prof
        # animal-mean profile per state
        am = {}
        for (a, s), ga in g.groupby(["animal", "state_class"]):
            am[(a, s)] = _median([run_vecs[lab][k] for k in ga["key"]]) if len(ga) > 1 \
                else run_vecs[lab][ga["key"].iloc[0]]
        for (k1, v1), (k2, v2) in itertools.combinations(sorted(am.items()), 2):
            (a1, s1), (a2, s2) = k1, k2
            if a1 == a2 and s1 == s2:
                continue
            kind = (PAIR_KINDS[0] if s1 == s2 else
                    PAIR_KINDS[1] if a1 == a2 else PAIR_KINDS[2])
            pairs.append(dict(label=lab, kind=kind, animal_1=a1, state_1=s1, animal_2=a2,
                              state_2=s2, r=ff._r(v1, v2)[0]))
        r_states, n_reg = ff._r(prof[AWAKE], prof[ANESTH])
        row = dict(label=lab, species=sp, r_states=r_states, n_regions=n_reg)
        for s in STATES:
            gs = g[g["state_class"] == s]
            row.update({f"n_runs_{s}": len(gs), f"n_animals_{s}": int(gs["animal"].nunique()),
                        f"states_{s}": ", ".join(f"{k} {v}" for k, v in
                                                 gs["state"].value_counts().items()),
                        f"r_consensus_{s}": ff._r(prof[s], prof["consensus"])[0]})
            m = maps.get((lab, s), {})
            row[f"map_{s}"] = m.get("map_path", "")
            row[f"png_{s}"] = m.get("png", "")
            row["colour_range"] = m.get("colour_range", row.get("colour_range", np.nan))
            row["voxel_r_between_states"] = m.get("voxel_r_between_states",
                                                  row.get("voxel_r_between_states", np.nan))
        both = set(g.loc[g["state_class"] == AWAKE, "animal"]) & \
            set(g.loc[g["state_class"] == ANESTH, "animal"])
        row["n_animals_both_states"] = len(both)
        summary.append(row)
    summary = pd.DataFrame(summary)
    pairs = pd.DataFrame(pairs, columns=["label", "kind", "animal_1", "state_1", "animal_2",
                                         "state_2", "r"])
    lines = ["Awake vs anaesthetised inside the datasets holding both states (descriptive; "
             "every state other than awake = anaesthetised)"]
    for r in summary.itertuples():
        lines.append(f"  {r.label}: awake {r.n_runs_awake} runs / {r.n_animals_awake} animals "
                     f"({r.states_awake}); anaesthetised {r.n_runs_anesthetized} runs / "
                     f"{r.n_animals_anesthetized} animals ({r.states_anesthetized}); "
                     f"{r.n_animals_both_states} animals in both states")
        lines.append(f"    r(awake, anaesthetised) median profiles = {r.r_states:.3f} "
                     f"({r.n_regions} values); r with LOSO consensus: awake "
                     f"{r.r_consensus_awake:.3f}, anaesthetised {r.r_consensus_anesthetized:.3f}"
                     + (f"; voxelwise r between the two maps = {r.voxel_r_between_states:.3f}"
                        if np.isfinite(r.voxel_r_between_states) else
                        "; voxelwise r: maps not computed (state_group_maps.py)"))
        p = pairs[pairs["label"] == r.label]
        for kind in PAIR_KINDS:
            v = p.loc[p["kind"] == kind, "r"].dropna()
            lines.append(f"    {kind}: {len(v)} pairs of animal profiles, median r = "
                         + (f"{v.median():.3f} [IQR {v.quantile(.25):.3f}, {v.quantile(.75):.3f}]"
                            if len(v) else "n/a"))
    return dict(summary=summary, pairs=pairs, profiles=profiles, lines=lines)


# ═══════════════════════════════════════════════════════════════════════════
# Drawing
# ═══════════════════════════════════════════════════════════════════════════
def _load_pngs(row):
    from Paper.analysis.fmri_fingerprints import paper_panels as fpp
    raw = {}
    for s in STATES:
        p = row.get(f"png_{s}") or ""
        raw[s] = plt.imread(p) if p and os.path.isfile(p) else None
    if not any(im is not None for im in raw.values()):
        return {}
    return fpp._compress_gaps(fpp._aligned_crops(raw))


def draw(stats, common, order, region_mask, out_base, pct=90.0):
    """Supp figure at print size (183 mm wide). Per dataset: a strip of the two state
    maps (a) above the fingerprint (b) and the profile similarity (c).
    out_base: path without extension."""
    summ, pairs, profs = stats["summary"], stats["pairs"], stats["profiles"]
    labels = list(summ["label"])
    if not labels:
        print("  [supp16] no dataset holding both states -- skipped")
        return None
    o = np.arange(len(common)) if order is None else np.asarray(order)
    mask = (np.ones(len(common), bool) if region_mask is None
            else np.asarray(region_mask, bool))
    with ff._print_size():
        W = ps.WIDTH_2COL
        maps_h, low_h, name_h, gap_h, foot_h = 34, 74, 6, 6, 14      # mm
        H = len(labels) * (name_h + maps_h + low_h + gap_h) + foot_h
        fig = plt.figure(figsize=(W, H * ps.MM))
        fy = lambda mm: 1 - mm / H                                   # mm from top -> fig y
        y = 0.0
        for i, lab in enumerate(labels):
            row = summ[summ["label"] == lab].iloc[0].to_dict()
            sp = row["species"]
            fig.text(0.01, fy(y + 1), ff.pretty_label(lab), fontsize=7, fontweight="bold",
                     color=ps.text_color(sp), ha="left", va="top")
            y += name_h
            # a: maps
            imgs = _load_pngs(row)
            rng_ = row.get("colour_range", np.nan)
            ax_maps = []
            for j, s in enumerate(STATES):
                x0 = 0.035 + j * 0.485
                ax = fig.add_axes([x0, fy(y + maps_h - 9), 0.455, (maps_h - 13) / H])
                ax_maps.append(ax)
                im = imgs.get(s)
                if im is not None:
                    ax.imshow(im, interpolation="lanczos", aspect="equal")
                else:
                    ax.text(0.5, 0.5, "map not computed (state_group_maps.py)", ha="center",
                            va="center", fontsize=5.5, color="#999999", transform=ax.transAxes)
                ax.set_axis_off()
                ax.set_title(f"{STATE_NAME[s]}: {row[f'n_runs_{s}']} runs, "
                             f"{row[f'n_animals_{s}']} animals", fontsize=6,
                             color=ff._state_color(s), pad=2)
            if np.isfinite(rng_):
                cax = fig.add_axes([0.38, fy(y + maps_h - 4.5), 0.24, 1.6 / H])
                cb = fig.colorbar(plt.cm.ScalarMappable(Normalize(-rng_, rng_), RYBCBL),
                                  cax=cax, orientation="horizontal")
                cb.set_ticks([-rng_, 0, rng_])
                cb.set_ticklabels([f"{-rng_:.2g}", "0", f"{rng_:.2g}"])
                cb.ax.tick_params(labelsize=5, length=1.5, pad=1)
                cb.outline.set_linewidth(0.4)
                cb.set_label("mean Fisher z (same range for both states)", fontsize=5,
                             labelpad=1)
            y += maps_h
            # b: radar (square, labels need ~14 mm around)
            side = low_h - 30
            ax_radar = fig.add_axes([0.12, fy(y + 15 + side), side / (W / ps.MM), side / H],
                                    projection="polar")
            p = profs[lab]
            arrays = [p[s][o] for s in STATES if p.get(s) is not None] + (
                [p["consensus"][o]] if p.get("consensus") is not None else [])
            theta = ff._spider_axes(ax_radar, [ff.short_region(common[k]) for k in o],
                                    ff._rlim(arrays), [not mask[k] for k in o],
                                    regions=[common[k] for k in o])
            for s in STATES:
                if p.get(s) is not None:
                    ff._poly(ax_radar, theta, p[s][o], color=ff._state_color(s), lw=0.9,
                             fill=True)
            if p.get("consensus") is not None:
                ff._poly(ax_radar, theta, p["consensus"][o], color="black", lw=0.7, ls="--")
            # c: similarity
            ax = fig.add_axes([0.63, fy(y + low_h - 16), 0.34, (low_h - 26) / H])
            pp = pairs[pairs["label"] == lab]
            rs = np.random.default_rng(i)
            for k, kind in enumerate(PAIR_KINDS):
                v = pp.loc[pp["kind"] == kind, "r"].dropna().values
                if len(v):
                    ax.scatter(k + rs.uniform(-0.15, 0.15, len(v)), v, s=4, lw=0,
                               color=PAIR_COLORS[kind], alpha=0.85, zorder=3)
                    ax.plot([k - 0.28, k + 0.28], [np.median(v)] * 2, color="black", lw=0.8,
                            zorder=4)
            if np.isfinite(row["r_states"]):
                ax.scatter([3.0], [row["r_states"]], marker="D", s=14, color="#7B3294",
                           edgecolor="black", lw=0.4, zorder=5)
            vr = row.get("voxel_r_between_states", np.nan)
            if vr is not None and np.isfinite(vr):
                ax.scatter([3.0], [vr], marker="s", s=12, facecolor="white",
                           edgecolor="#7B3294", lw=0.8, zorder=5)
            npair = [int((pp["kind"] == k).sum()) for k in PAIR_KINDS]
            ax.set_xticks(range(4))
            ax.set_xticklabels([f"other animal,\nsame state\n(n={npair[0]})",
                                f"same animal,\nother state\n(n={npair[1]})",
                                f"other animal,\nother state\n(n={npair[2]})",
                                "state\ngroups"], fontsize=5)
            ax.set_xlim(-0.6, 3.5)
            lo = np.nanmin(np.r_[pp["r"].values, row["r_states"], 0.0])
            ax.set_ylim(min(-0.1, lo - 0.05), 1.02)
            ax.axhline(0, color="#BBBBBB", lw=0.4, zorder=1)
            ax.set_ylabel("Pearson r between profiles", fontsize=5.5)
            if i == 0:
                for axl, letter, dx in ((ax_maps[0], "a", 0.03), (ax_radar, "b", 0.10),
                                        (ax, "c", 0.06)):
                    bb = axl.get_position()
                    fig.text(bb.x0 - dx, bb.y1 + (0.004 if letter != "b" else 0.03), letter,
                             fontsize=ps.PANEL_LABEL_SIZE, fontweight="bold", ha="left",
                             va="bottom")
            y += low_h + gap_h
        handles = [Line2D([], [], color=ff._state_color(s), lw=1.2, label=STATE_NAME[s])
                   for s in STATES]
        handles += [Line2D([], [], color="black", lw=0.8, ls="--",
                           label="leave-one-species-out consensus")]
        handles += [Line2D([], [], marker="o", ls="none", ms=2.5, color=PAIR_COLORS[k],
                           label=f"pairs: {k}") for k in PAIR_KINDS]
        handles += [Line2D([], [], marker="D", ls="none", ms=3, color="#7B3294", mec="black",
                           mew=0.4, label="state profiles of b"),
                    Line2D([], [], marker="s", ls="none", ms=3, mfc="white", mec="#7B3294",
                           mew=0.8, label="voxelwise, maps of a")]
        fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=5, frameon=False,
                   bbox_to_anchor=(0.5, 0.0))
        ps.save(fig, out_base)
        plt.close(fig)
    _write_legend(stats, out_base, pct)
    return out_base


def _write_legend(stats, out_base, pct):
    s = stats["summary"]
    n_txt = "; ".join(
        f"{ff.pretty_label(r.label)}: awake {r.n_runs_awake} runs from {r.n_animals_awake} "
        f"animals, anaesthetised {r.n_runs_anesthetized} runs from "
        f"{r.n_animals_anesthetized} animals ({r.states_anesthetized})"
        for r in s.itertuples())
    txt = (
        "Supplementary Figure | Seed connectivity of the selected runs in the awake and the "
        "anaesthetised state, in the two datasets that contain both. "
        "a, Retrosplenial seed maps (group mean Fisher z, AFNI 3dttest++ on the runs kept by "
        "the final selection round, split by arousal state) on each dataset's inflated "
        "cortical surface; the two maps of a dataset share one symmetric colour range "
        f"(± the {pct:g}th percentile of |z| over both maps). "
        "b, Regional fingerprints: median scaled profile of the runs of each state, with the "
        "leave-one-species-out consensus (dashed; the round-1 selection target). "
        "c, Pearson r between animal-mean profiles, for pairs of different animals in the "
        "same state, of the same animal in the two states and of different animals in "
        "different states (dots = pairs, bars = medians); diamond, r between the two "
        "state profiles of b; square, voxelwise r between the two maps of a. "
        f"{n_txt}. Arousal state is nested in sessions and acquisition protocol, and every "
        "anaesthetic regimen is pooled as 'anaesthetised'; the panels are descriptive and "
        "no awake versus anaesthetised test is made here (within-animal tests: "
        "Supplementary Fig. 13).")
    with open(out_base + "_legend.txt", "w", encoding="utf-8") as f:
        f.write(txt + "\n")
