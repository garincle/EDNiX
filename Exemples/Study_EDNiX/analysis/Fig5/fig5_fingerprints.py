"""
fig5_fingerprints.py
====================
Interpretive figures for Figure 5, computed from profiles the pipeline already
holds. Pure numpy/pandas/matplotlib; the only file access is the listing of
each BIDS func/ folder, used to read the arousal state of a run from its task
name (Marmoset_Hori).

    F6ab   distance of every dataset to TWO references vs divergence time
           (a) leave-one-species-out (LOSO) consensus -- magnitude of departure
           (b) the Human datasets                     -- direction of departure
    F1_fp  F1  redrawn as a connectivity fingerprint, all runs | selected
    F1c_fp F1c redrawn as a connectivity fingerprint, all runs | selected
    F15    each dataset's fingerprint over its reference, ordered by distance
    F16    dataset x dataset distances + sibling test (species vs dataset)
    F17    group means before vs after selection
    F18    arousal state, contrasted INSIDE the datasets that contain both
           awake and anaesthetised scans

Connectivity fingerprint = the vector of a seed's connectivity with a set of
target regions, drawn as a polar plot (Passingham, Stephan & Kotter 2002 Nat
Rev Neurosci 3:606; Mars, Passingham & Jbabdi 2018 Trends Cogn Sci 22:1026).
With PROFILE_KIND="matrix" (edge profiles from correlation matrices) the polar
plot is replaced by a region x region matrix.

Arousal-state rules live in the launcher's BIDS_CONFIGS, key "state":
    "awake" | "anesthetized" | "isoflurane" | ...   whole dataset
    {"from": "task", "map": {"restawake": "awake", "restanest": "anesthetized"}}
    {"from": "subject", "map": [["MEDISO", "medetomidine+isoflurane"],
                                ["ISO", "isoflurane"], ["AW", "awake"]]}
    None                                             unknown
Every state other than "awake" counts as anaesthetised in the binary contrast.
"""

import os
import re
import glob
import itertools
import contextlib

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats as st

import Exemples.Study_EDNiX.analysis.Fig5.fig5_analysis as fa

opj = os.path.join
opb = os.path.basename

REGION_SHORT = {
    "pallidum": "Pallidum",
    "striatum": "Striatum",
    "thalamus": "Thalamus",
    "amygdala": "Amygdala",
    "hippocampal_formation": "Hippocampus",
    "olfactory_cortex": "Olfactory",
    "periarchicortex": "Periarchicortex",
    "insula_and_others_in_lateral_sulcus": "Insula",
    "somatosensory_cortex": "Somatosensory",
    "motor_and_premotor": "Motor",
    "auditory_cortex_(superior_temporal)": "Auditory",
    "middle_temporal_inferior_temporal_temporal_pole_(mipt)": "MT/IT/TP",
    "prefrontal_cortex_(pfc)": "PFC",
    "posterior_parietal_cortex": "Post. parietal",
    "visual_pre_and_extra_striate_cortex": "Extrastriate",
    "visual_striate_cortex": "V1",
}
NODE_ORDER = list(REGION_SHORT)          # subcortical -> allocortex -> isocortex
EDGE_LINK = " <> "
LINESTYLES = ["-", (0, (5, 2)), (0, (1, 1.2)), "-."]
AWAKE, ANESTH = "awake", "anesthetized"


CHANNEL_SIGN = {"pos": "+", "neg": "\u2212"}   # fraction_posneg channels


def pretty_label(label):
    """Dataset name for display: 'Macaque_BenHamed' -> 'Macaque BenHamed'."""
    return str(label).replace("_", " ")


def short_region(name):
    """Short label of a region; a two-channel profile entry reads 'X +' /
    'X \u2212' (fraction_posneg) or 'X P80' / 'X P20' (percentiles)."""
    key = str(name).lower().strip()
    base, sep, ch = key.rpartition("|")
    if sep and (ch in CHANNEL_SIGN or re.fullmatch(r"q\d{1,2}", ch)):
        tag = CHANNEL_SIGN.get(ch, "P" + ch[1:] if ch.startswith("q") else ch)
        return f"{REGION_SHORT.get(base, base.replace('_', ' ')[:18])} {tag}"
    return REGION_SHORT.get(key, key.replace("_", " ")[:18])


def _legend(path, title, what, how, reading):
    """Same layout as fig5_pipeline._write_legend."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{title}\n{'=' * len(title)}\n\n")
        f.write(f"WHAT IS PLOTTED\n---------------\n{what}\n\n")
        f.write(f"HOW IT WAS COMPUTED\n-------------------\n{how}\n\n")
        f.write(f"HOW TO READ IT\n--------------\n{reading}\n")


def _save(fig, path, dpi=300, pdf=False):
    """pdf=True (paper panels): also writes <path>.pdf, text kept as text."""
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    if pdf:
        with plt.rc_context({"pdf.fonttype": 42}):
            fig.savefig(os.path.splitext(path)[0] + ".pdf", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  [saved] {path}" + (" (+ .pdf)" if pdf else ""))


# ── print size of the paper panels ────────────────────────────────────────────
# Inside `with _print_size():` the figure is drawn at its printed size
# (fa.paper_style: 5-7 pt text, thin lines) and the shared helpers below pick
# their print sizes through _P(screen value, print value).
_PRINT = {"on": False}


@contextlib.contextmanager
def _print_size():
    old = _PRINT["on"]
    _PRINT["on"] = True
    try:
        with fa.paper_style():
            yield
    finally:
        _PRINT["on"] = old


def _P(screen, printed):
    return printed if _PRINT["on"] else screen


def _color(sp):
    return fa.SPECIES_COLORS.get(sp, "#555555")


def _tcolor(sp):
    """Text version of the species colour: the palette's yellow (Rat) is
    unreadable as text on white, so light colours are darkened."""
    from matplotlib.colors import to_rgb
    rgb = np.array(to_rgb(_color(sp)))
    lum = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
    return tuple(rgb * 0.62) if lum > 0.7 else tuple(rgb)


# ═══════════════════════════════════════════════════════════════════════════
# Arousal state
# ═══════════════════════════════════════════════════════════════════════════

def run_token(filename):
    """Run id from a file name -- same rule as
    SBA_Region_Profile_Index._extract_run_id (kept in sync by a unit test)."""
    m = re.search(r"run-(RL|LR)-(\d+)", filename)
    if m:
        return f"{m.group(1)}-{m.group(2)}"
    m = re.search(r"run-(\d+)", filename)
    return m.group(1) if m else "NA"


def task_of(filename):
    m = re.search(r"task-([A-Za-z0-9]+)", filename)
    return m.group(1) if m else None


def state_class(state):
    """'awake' / 'anesthetized' / None (unknown, ambiguous, mixed)."""
    if state is None or str(state).lower() in ("", "unknown", "ambiguous", "mixed"):
        return None
    return AWAKE if str(state).lower() == AWAKE else ANESTH


class StateResolver:
    """State of one unit key (subject, session, run, seed) from the rule in a
    BIDS_CONFIGS entry. With rule {"from": "task"} the raw BIDS func/ folder of
    that subject/session is listed once and cached."""

    def __init__(self, bids_config):
        self.rule = bids_config.get("state")
        self.bids_dir = bids_config.get("bids_dir")
        self._tasks = {}
        self.animal_re = (self.rule.get("animal") if isinstance(self.rule, dict) else None)

    def animal(self, key):
        """Animal identity of a unit. Default: the BIDS subject. With
        rule["animal"] = regex, the regex is removed from the subject label
        (case-insensitive result), e.g. Mouse "jgrAesISOc11L" -> "c11l": the
        state is encoded in the subject label but the animal is the same."""
        subj = str(key[0])
        if not self.animal_re:
            return subj
        return (re.sub(self.animal_re, "", subj) or subj).lower()

    def _task_pairs(self, subj, ses):
        key = (subj, ses)
        if key not in self._tasks:
            files = []
            for d in (opj(self.bids_dir, f"sub-{subj}", f"ses-{ses}", "func"),
                      opj(self.bids_dir, f"sub-{subj}", "func")):
                if os.path.isdir(d):
                    files = sorted(glob.glob(opj(d, "*_bold.nii*")))
                    if files:
                        break
            self._tasks[key] = [(run_token(opb(f)), task_of(opb(f))) for f in files]
        return self._tasks[key]

    def state(self, key):
        subj, ses, run = str(key[0]), str(key[1]), str(key[2])
        rule = self.rule
        if rule is None:
            return "unknown"
        if isinstance(rule, str):
            return rule
        src, mapping = rule.get("from"), rule.get("map")
        if src == "subject":
            items = mapping.items() if isinstance(mapping, dict) else mapping
            for pattern, state in items:          # first match wins: MEDISO before ISO
                if pattern in subj:
                    return state
            return "unknown"
        if src == "task":
            tasks = {t for r, t in self._task_pairs(subj, ses)
                     if t and (run == "ALL" or r == run)}
            states = {mapping.get(t, "unknown") for t in tasks}
            if not states:
                return "unknown"
            if len(states) == 1:
                return states.pop()
            return "mixed" if run == "ALL" else "ambiguous"
        raise ValueError(f"state rule {rule!r}: 'from' must be 'task' or 'subject'")


def dataset_state(resolver, keys):
    """'awake' | 'anesthetized' | 'mixed' | 'unknown' for a whole dataset."""
    classes = {state_class(resolver.state(k)) for k in keys} - {None}
    if not classes:
        return "unknown"
    return classes.pop() if len(classes) == 1 else "mixed"


STATE_STYLE = {
    AWAKE: dict(fillstyle="full"),
    ANESTH: dict(fillstyle="none"),
    "mixed": dict(fillstyle="left"),
    "unknown": dict(fillstyle="none"),
}


def _state_marker(ax, x, y, color, state, size=8, zorder=4):
    kw = STATE_STYLE.get(state, STATE_STYLE["unknown"])
    ax.plot([x], [y], marker="D" if state == "unknown" else "o", ms=_P(size, 0.5 * size),
            color=color, mec=color, mew=_P(1.6, 0.7), linestyle="none", zorder=zorder,
            markerfacecoloralt="white", **kw)


def _state_legend_handles():
    out = []
    for st_, lab in ((AWAKE, "awake"), (ANESTH, "anaesthetised"),
                     ("mixed", "both states in the dataset"), ("unknown", "state unknown")):
        kw = STATE_STYLE[st_]
        out.append(Line2D([], [], marker="D" if st_ == "unknown" else "o", ms=_P(7, 3.6),
                          color="#444444", mec="#444444", mew=_P(1.4, 0.7), linestyle="none",
                          markerfacecoloralt="white", label=lab, **kw))
    return out


# ═══════════════════════════════════════════════════════════════════════════
# Vectors, references, distances
# ═══════════════════════════════════════════════════════════════════════════

def group_vectors(profiles, common, value_col, norm):
    """{label: normalised vector over `common`} from {label: df} or
    {label: (species, df)}; datasets with < 3 finite regions are dropped."""
    out = {}
    for lab, v in profiles.items():
        prof = v[1] if isinstance(v, tuple) else v
        if prof is None:
            continue
        vec = fa._normalize_values(fa._ordered_common_values(prof, common, value_col), norm)
        vec = np.asarray(vec, dtype=float)
        if np.isfinite(vec).sum() >= 3:
            out[lab] = vec
    return out


def species_medians(vecs, species_of_label, exclude_species=None):
    grouped = {}
    for lab, v in sorted(vecs.items()):
        sp = species_of_label.get(lab, lab)
        if sp != exclude_species:
            grouped.setdefault(sp, []).append(v)
    with np.errstate(invalid="ignore"):
        return {sp: np.nanmedian(np.vstack(v), axis=0) for sp, v in grouped.items()}


def consensus(vecs, species_of_label, exclude_species=None):
    """Median across species of the per-species medians -- identical to
    fa.build_consensus_pattern(..., species_of_label=...) on normalised input."""
    spm = species_medians(vecs, species_of_label, exclude_species)
    if not spm:
        return None
    with np.errstate(invalid="ignore"):
        return np.nanmedian(np.vstack(list(spm.values())), axis=0)


def reference_vector(vecs, species_of_label, ref_species="Human", exclude_label=None):
    vs = [v for lab, v in sorted(vecs.items())
          if species_of_label.get(lab) == ref_species and lab != exclude_label]
    if not vs:
        return None
    with np.errstate(invalid="ignore"):
        return np.nanmedian(np.vstack(vs), axis=0)


def _r(a, b):
    if a is None or b is None:
        return np.nan, 0
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3 or np.std(a[m]) == 0 or np.std(b[m]) == 0:
        return np.nan, int(m.sum())
    return float(np.corrcoef(a[m], b[m])[0, 1]), int(m.sum())


def order_labels(labels, species_of_label, divergence):
    return sorted(labels, key=lambda l: (divergence.get(species_of_label.get(l), 1e9),
                                         species_of_label.get(l, ""), l))


def dataset_distances(vecs, species_of_label, divergence, ref_species="Human"):
    """One row per dataset: 1 - r to its LOSO consensus (F6a, identical to
    fa.compute_loo_species_distances) and 1 - r to the reference species
    (F6b; for a dataset OF the reference species: the other dataset(s) of that
    species -- the within-species floor)."""
    rows = []
    for lab in order_labels(vecs, species_of_label, divergence):
        sp = species_of_label.get(lab, lab)
        v = vecs[lab]
        cons = consensus(vecs, species_of_label, exclude_species=sp)
        r_c, n = _r(v, cons)
        ref = reference_vector(vecs, species_of_label, ref_species, exclude_label=lab)
        r_h, _ = _r(v, ref)
        rows.append(dict(label=lab, species=sp, divergence_mya=divergence.get(sp, np.nan),
                         d_consensus=1 - r_c, r_consensus=r_c,
                         d_reference=1 - r_h, r_reference=r_h,
                         reference_is_own_species=sp == ref_species,
                         n_regions=n))
    return pd.DataFrame(rows)


def species_table(dist_df, ref_species="Human"):
    g = dist_df.groupby("species", as_index=False).agg(
        divergence_mya=("divergence_mya", "first"),
        d_consensus=("d_consensus", "mean"),
        d_reference=("d_reference", "mean"),
        r_consensus=("r_consensus", "mean"),
        r_reference=("r_reference", "mean"),
        n_datasets=("label", "count"))
    return g.sort_values("divergence_mya").reset_index(drop=True)


def perm_corr(x, y, max_exact=8, n_mc=20000, seed=0):
    """Pearson r with a permutation p over ALL n! orderings of x when n <= 8
    (exact), Monte Carlo otherwise. p_two = P(|r*| >= |r|), p_pos = P(r* >= r).
    The identity permutation is part of the null, so p >= 1/n!."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = len(x)
    if n < 3 or np.std(x) == 0 or np.std(y) == 0:
        return dict(r=np.nan, p_two=np.nan, p_pos=np.nan, n=n, n_perm=0, exact=False)
    r_obs = float(np.corrcoef(x, y)[0, 1])
    if n <= max_exact:
        idx = np.array(list(itertools.permutations(range(n))))
        exact = True
    else:
        rng = np.random.default_rng(seed)
        idx = np.array([rng.permutation(n) for _ in range(n_mc)])
        exact = False
    X = x[idx]
    Xc = X - X.mean(axis=1, keepdims=True)
    yc = y - y.mean()
    rs = (Xc @ yc) / (np.sqrt((Xc ** 2).sum(axis=1)) * np.sqrt((yc ** 2).sum()))
    eps = 1e-12
    return dict(r=r_obs, p_two=float(np.mean(np.abs(rs) >= abs(r_obs) - eps)),
                p_pos=float(np.mean(rs >= r_obs - eps)), n=n, n_perm=len(rs), exact=exact)


def _fmt_p(p):
    if not np.isfinite(p):
        return "n/a"
    return f"{p:.3f}" if p >= 0.001 else f"{p:.1e}"


# ═══════════════════════════════════════════════════════════════════════════
# Label placement (no overlap) for scatter panels
# ═══════════════════════════════════════════════════════════════════════════

def _label_points(ax, xs, ys, texts, colors, dx, min_gap, xlim_right=None):
    """Labels to the right of each point; points within 1.5*dx in x share a
    column whose labels are pushed apart vertically by at least min_gap.
    At print size the labels are only recorded here and placed without overlap
    once the layout is final (_place_labels, just before saving)."""
    if _PRINT["on"]:
        ax._fig5_labels = getattr(ax, "_fig5_labels", []) + [
            (float(x), float(y), pretty_label(t), c) for x, y, t, c in zip(xs, ys, texts, colors)]
        return
    pts = sorted(zip(xs, ys, texts, colors), key=lambda t: (t[0], t[1]))
    cols = []
    for p in pts:
        if cols and abs(p[0] - cols[-1][0][0]) <= 1.5 * dx:
            cols[-1].append(p)
        else:
            cols.append([p])
    for col in cols:
        col = sorted(col, key=lambda t: t[1])
        ly = [c[1] for c in col]
        for i in range(1, len(ly)):
            ly[i] = max(ly[i], ly[i - 1] + min_gap)
        shift = (np.mean([c[1] for c in col]) - np.mean(ly))
        ly = [v + shift for v in ly]
        for (x, y, t, c), yl in zip(col, ly):
            xt = max(c_[0] for c_ in col) + dx
            ax.annotate(pretty_label(t), xy=(x, y), xytext=(xt, yl), fontsize=_P(7, 5),
                        color=c, va="center",
                        arrowprops=dict(arrowstyle="-", color="#BBBBBB", lw=_P(0.6, 0.4),
                                        shrinkA=0, shrinkB=_P(3, 1.5)))


def _place_labels(fig, fs=5):
    """Point labels recorded by _label_points at print size, placed once the
    layout is final (after tight_layout, before saving)."""
    for ax in fig.axes:
        items = getattr(ax, "_fig5_labels", None)
        if items:
            fa.repel_labels(ax, items, fs=fs, obstacles=getattr(ax, "_fig5_obstacles", ()))
            ax._fig5_labels = []


# ═══════════════════════════════════════════════════════════════════════════
# F6ab -- two references
# ═══════════════════════════════════════════════════════════════════════════

def _fit_line(ax, sp_df, ycol, color, ls, label, xcol="divergence_mya"):
    x, y = sp_df[xcol].values, sp_df[ycol].values
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return
    b1, b0 = np.polyfit(x[m], y[m], 1)
    xx = np.array([x[m].min(), x[m].max()])
    ax.plot(xx, b0 + b1 * xx, color=color, ls=ls, lw=_P(1.4, 0.8), zorder=1, label=label)


def _f6_panel(ax, dist_df, ycol, states, title, ylab, ref_species, floor=None,
              xcol="divergence_mya", xlab="Divergence time from Human (Mya)", paper=False):
    """paper=True: no panel title, no text on the grey floor (both go to the
    .txt legend); everything else identical."""
    sp_df = species_table(dist_df, ref_species)
    if xcol not in sp_df:
        sp_df[xcol] = sp_df["species"].map(dist_df.groupby("species")[xcol].first())
    is_ref = ycol in ("d_reference", "r_reference")
    for r in dist_df.itertuples():
        if is_ref and r.species == ref_species:
            continue
        _state_marker(ax, getattr(r, xcol), getattr(r, ycol), _color(r.species),
                      states.get(r.label, "unknown"))
    for r in sp_df.itertuples():
        if is_ref and r.species == ref_species:
            continue
        ax.scatter([getattr(r, xcol)], [getattr(r, ycol)], s=_P(230, 55), facecolor="none",
                   edgecolor=_color(r.species), lw=_P(1.0, 0.6), alpha=0.6, zorder=2)
        if _PRINT["on"]:                       # labels must not cover the rings
            ax._fig5_obstacles = getattr(ax, "_fig5_obstacles", []) + [
                (float(getattr(r, xcol)), float(getattr(r, ycol)), 4.0)]
    if floor is not None and np.isfinite(floor).any():
        if ycol.startswith("r_"):          # r between the two reference datasets: ceiling
            lo, hi = float(np.nanmin(floor)), 1.0
        else:                              # 1 - r: floor
            lo, hi = 0.0, float(np.nanmax(floor))
        ax.axhspan(lo, hi, color="#BBBBBB", alpha=0.35, lw=0, zorder=0)
        if not paper:
            ax.text(0.01, lo, f" within-{ref_species}: each {ref_species} dataset vs the "
                              f"other one ({', '.join(f'{v:.2f}' for v in floor)})",
                    transform=ax.get_yaxis_transform(), ha="left", va="top",
                    fontsize=6.5, color="#555555")
    sub = dist_df[dist_df["species"] != ref_species] if is_ref else dist_df
    yv = sub[ycol].values
    yall = np.r_[yv, floor] if floor is not None else yv
    if ycol.startswith("r_"):
        yr = (np.nanmax(yall) - np.nanmin(yall)) or 0.1
    else:
        yr = (np.nanmax(yall) - min(0, np.nanmin(yall))) if len(yall) else 1.0
    xv = dist_df[xcol].values.astype(float)
    x_lo, x_hi = (np.nanmin(xv), np.nanmax(xv)) if np.isfinite(xv).any() else (0.0, 1.0)
    xr = (x_hi - x_lo) or 1.0
    _label_points(ax, sub[xcol].values, yv, list(sub["label"]),
                  [_tcolor(s) for s in sub["species"]], dx=0.037 * xr, min_gap=0.045 * yr)
    if title and not paper:
        ax.set_title(title, fontsize=10, fontweight="bold", loc="left")
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    if xcol == "divergence_mya":
        ax.set_xlim(-6, max(x_hi, 1) * 1.33)
    else:
        ax.set_xlim(x_lo - 0.06 * xr, x_hi + 0.33 * xr)
    if ycol.startswith("r_"):
        top = 1.0 if floor is not None and np.isfinite(floor).any() else np.nanmax(yall)
        ax.set_ylim(np.nanmin(yall) - 0.12 * yr, top + 0.06 * yr)
    else:
        ax.set_ylim(min(0, np.nanmin(yall)) - 0.03 * yr, np.nanmax(yall) + 0.12 * yr)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return sp_df


def f6_tests(dist_df, ref_species="Human"):
    """Divergence tests; brain-mass tests too when dist_df has log_brain_g."""
    sp_df = species_table(dist_df, ref_species)
    nonref = sp_df[sp_df["species"] != ref_species]
    out = dict(
        consensus_with_ref=perm_corr(sp_df["divergence_mya"], sp_df["r_consensus"]),
        consensus_without_ref=perm_corr(nonref["divergence_mya"], nonref["r_consensus"]),
        reference=perm_corr(nonref["divergence_mya"], nonref["r_reference"]))
    if "log_brain_g" in dist_df and dist_df["log_brain_g"].notna().sum() >= 4:
        sp_df["log_brain_g"] = sp_df["species"].map(
            dist_df.groupby("species")["log_brain_g"].first())
        nonref = sp_df[sp_df["species"] != ref_species]
        out["brain_mass_with_ref"] = perm_corr(sp_df["log_brain_g"], sp_df["r_consensus"])
        out["brain_mass_without_ref"] = perm_corr(nonref["log_brain_g"], nonref["r_consensus"])
    return out, sp_df


def _np_header(t, unit="species"):
    """'n = 7 species, p = 0.015' (supplementary panels)."""
    if t is None or not np.isfinite(t.get("r", np.nan)):
        return "n/a"
    return f"n = {t['n']} {unit}, p = {_fmt_p(t['p_two'])}"


def _r_stars(t):
    """'r = +0.62 *' -- the statistics written on a paper panel (n, test and
    exact p in the .txt)."""
    if t is None or not np.isfinite(t.get("r", np.nan)):
        return "n/a"
    return f"r = {t['r']:+.2f} {fa.stars(t['p_two'])}"


def _test_line(name, t):
    if not np.isfinite(t["r"]):
        return f"{name}: n/a"
    kind = f"exact, {t['n_perm']} perm." if t["exact"] else f"{t['n_perm']} MC perm."
    return (f"{name}: r = {t['r']:+.2f}, p = {_fmt_p(t['p_two'])} two-sided "
            f"(n = {t['n']} species, {kind})")


def plot_f6_two_references(dist_df, states, out_png, source, ref_species="Human",
                           level=None, paper=False, brain_mass=None):
    """F6ab. paper=True (6a): same panels at print size (183 mm), without
    suptitle, panel titles, statistics boxes or floor text (all in the .txt);
    top right of each panel: r and its stars (exact permutation p). With
    brain_mass ({species: g}) a third panel plots the distance to the LOSO
    consensus against log10 brain mass."""
    d = dist_df.copy()
    has_bm = False
    if brain_mass:
        d["log_brain_g"] = np.log10(d["species"].map(brain_mass).astype(float))
        has_bm = d.groupby("species")["log_brain_g"].first().notna().sum() >= 4
        if not has_bm:
            d = d.drop(columns=["log_brain_g"])
    tests, sp_df = f6_tests(d, ref_species)
    n_pan = 3 if (paper and has_bm) else 2
    ctx = _print_size() if paper else contextlib.nullcontext()
    G = f"without {ref_species}"

    def _stats_text(ax, lines):
        ax.text(1.0, 1.0, "\n".join(lines), transform=ax.transAxes, ha="right", va="bottom",
                fontsize=6, linespacing=1.15)

    with ctx:
        fig, axes = plt.subplots(1, n_pan, figsize=((fa.PAPER_WIDTH["double"], 64 * fa.MM)
                                                    if paper else (7.25 * n_pan, 6.2)))
        ax = axes[0]
        _f6_panel(ax, d, "r_consensus", states,
                  "a  Similarity to the leave-one-species-out consensus",
                  "r with the LOSO consensus" if paper else
                  "r  (dataset vs median of the OTHER species)", ref_species, paper=paper)
        _fit_line(ax, sp_df, "r_consensus", "black", "-", "fit, all species")
        _fit_line(ax, sp_df[sp_df["species"] != ref_species], "r_consensus", "#777777", "--",
                  f"fit, {G}")
        if paper:
            _stats_text(ax, [f"all: {_r_stars(tests['consensus_with_ref'])}",
                             f"{G}: {_r_stars(tests['consensus_without_ref'])}"])
        else:
            ax.text(0.01, 0.99, _test_line(f"with {ref_species}", tests["consensus_with_ref"])
                    + "\n" + _test_line(G, tests["consensus_without_ref"]),
                    transform=ax.transAxes, va="top", fontsize=7.2,
                    bbox=dict(fc="white", ec="#CCCCCC", lw=0.6))
        ax.legend(loc="best", fontsize=_P(7, 5.5), frameon=False)

        ax = axes[1]
        floor = d.loc[d["species"] == ref_species, "r_reference"].values
        _f6_panel(ax, d, "r_reference", states,
                  f"b  Similarity to {ref_species}",
                  f"r with {ref_species}" if paper else
                  f"r  (dataset vs median of the {ref_species} datasets)", ref_species,
                  floor=floor, paper=paper)
        _fit_line(ax, sp_df[sp_df["species"] != ref_species], "r_reference", "black", "-",
                  f"fit, species other than {ref_species}")
        if paper:
            _stats_text(ax, [_r_stars(tests["reference"])])
        else:
            ax.text(0.01, 0.99, _test_line("fit", tests["reference"]), transform=ax.transAxes,
                    va="top", fontsize=7.2, bbox=dict(fc="white", ec="#CCCCCC", lw=0.6))
        ax.legend(loc="lower right", fontsize=_P(7, 5.5), frameon=False)

        if n_pan == 3:
            ax = axes[2]
            _f6_panel(ax, d, "r_consensus", states, None,
                      "r with the LOSO consensus", ref_species,
                      xcol="log_brain_g", xlab="Brain mass (log$_{10}$ g)", paper=True)
            _fit_line(ax, sp_df, "r_consensus", "black", "-", "fit, all species",
                      xcol="log_brain_g")
            _fit_line(ax, sp_df[sp_df["species"] != ref_species], "r_consensus", "#777777",
                      "--", f"fit, {G}", xcol="log_brain_g")
            _stats_text(ax, [f"all: {_r_stars(tests['brain_mass_with_ref'])}",
                             f"{G}: {_r_stars(tests['brain_mass_without_ref'])}"])
            ax.legend(loc="lower right", fontsize=_P(7, 5.5), frameon=False)

        handles = _state_legend_handles() + [
            Line2D([], [], marker="o", ms=_P(13, 6.5), mfc="none", mec="#888888", ls="none",
                   mew=_P(1.0, 0.6), label="species mean (fitted point)")]
        fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=_P(7.5, 6),
                   frameon=False, bbox_to_anchor=(0.5, -0.02))
        if not paper:
            fig.suptitle(f"Retrosplenial profile: distance to two references vs divergence"
                         f"{f' (atlas level {level})' if level else ''} -- source: {source}",
                         fontsize=11, fontweight="bold", y=1.0)
            fig.tight_layout(rect=(0, 0.05, 1, 0.97))
            _save(fig, out_png)
        else:
            fig.tight_layout(rect=(0, 0.07, 1, 1), w_pad=1.5)
            _place_labels(fig)
            fa.save_paper(fig, out_png)
    return tests, sp_df


# ═══════════════════════════════════════════════════════════════════════════
# Fingerprints -- polar (seed) or matrix (edges)
# ═══════════════════════════════════════════════════════════════════════════

def _spider_axes(ax, labels, rlim, excluded, show_labels=True):
    """Polar axes with one spoke per entry, labelled with the full entry name:
    with a two-channel profile every spoke keeps its channel ('Insula P80',
    'Insula P20' -- or '+' / '−' with fraction_posneg), as in V12, so the
    two channels of every region are visible. Grey + '*' = not used by the
    selection score. show_labels=False: no spoke label at all (6b)."""
    n = len(labels)
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks(theta)
    ax.set_ylim(*rlim)
    texts = [f"{l}*" if ex else str(l) for l, ex in zip(labels, excluded)]
    if not show_labels:
        ax.set_xticklabels([])
    elif _PRINT["on"]:
        # print size: labels written along their spoke (radial), so that the
        # labels of neighbouring spokes never overlap, even with 32 spokes
        ax.set_xticklabels([])
        lab_r = rlim[1] + 0.05 * (rlim[1] - rlim[0])
        for th, txt, ex in zip(theta, texts, excluded):
            phi = (90.0 - np.degrees(th)) % 360.0
            phi = phi - 360.0 if phi > 180.0 else phi
            if -90.0 <= phi <= 90.0:
                rot, ha = phi, "left"
            else:
                rot, ha = (phi - 180.0 if phi > 0 else phi + 180.0), "right"
            ax.text(th, lab_r, txt, rotation=rot, rotation_mode="anchor", ha=ha, va="center",
                    fontsize=5, color="#999999" if ex else "#222222", clip_on=False)
    else:
        ax.set_xticklabels(texts, fontsize=7)
        for tl, ex in zip(ax.get_xticklabels(), excluded):
            tl.set_color("#999999" if ex else "#222222")
    ring = np.linspace(0, 2 * np.pi, 200)
    ax.plot(ring, np.zeros_like(ring), color="#666666", lw=_P(0.9, 0.5), ls=":", zorder=1)
    ticks = [t for t in (-1.0, 0.0, 1.0, 2.0) if rlim[0] < t < rlim[1]]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{t:g}" for t in ticks], fontsize=_P(6, 4.5), color="#777777")
    ax.set_rlabel_position(180.0 / max(n, 1))          # between the first two spokes
    ax.grid(color="#DDDDDD", lw=_P(0.6, 0.35))
    ax.spines["polar"].set_color("#BBBBBB")
    ax.spines["polar"].set_linewidth(_P(0.8, 0.5))
    return theta


def _texts_top(ax):
    """Top of the texts of a polar axes (spoke labels), in axes fraction."""
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    tops = [t.get_window_extent(r).y1 for t in ax.texts if t.get_text()]
    if not tops:
        return 1.0
    return float(ax.transAxes.inverted().transform((0, max(tops)))[1])


def _texts_bottom(ax):
    """Bottom of the texts of a polar axes (spoke labels), in figure fraction."""
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    lows = [t.get_window_extent(r).y0 for t in ax.texts if t.get_text()]
    if not lows:
        return float(ax.get_position().y0)
    return float(fig.transFigure.inverted().transform((0, min(lows)))[1])


def _close(theta, v):
    """Close a polar curve. The last point is at theta[0] + 2 pi, not theta[0]:
    going back from 2 pi - step to 0 makes fill_between sweep the whole circle
    (the light wedges seen on some fingerprints)."""
    return np.r_[theta, theta[0] + 2 * np.pi], np.r_[v, v[0]]


def _poly(ax, theta, v, fill=False, **kw):
    t, r = _close(theta, v)
    line, = ax.plot(t, r, **kw)
    if fill:
        ax.fill(t, np.where(np.isfinite(r), r, 0.0), color=line.get_color(), alpha=0.10, lw=0)
    return line


def _rlim(arrays):
    vals = np.concatenate([a[np.isfinite(a)] for a in arrays if a is not None])
    lo, hi = float(np.min(vals)), float(np.max(vals))
    pad = 0.08 * (hi - lo if hi > lo else 1.0)
    return min(lo, 0.0) - pad, hi + pad


def _region_order(cons_all, common):
    """Regions by increasing consensus value; with a two-channel profile, all
    'high' channels first then all 'low' ones, each block sorted."""
    key = np.where(np.isfinite(cons_all), cons_all, np.inf)
    block = np.array([fa._channel_block(c) for c in common])
    return np.lexsort((key, block))


def is_edge_profile(common):
    return bool(common) and all(" |" in c for c in common)


def parse_edge(name):
    body, _, cls = name.rpartition(" |")
    if cls == "homotopic" or EDGE_LINK not in body:
        return body, body, cls
    a, b = body.split(EDGE_LINK, 1)
    return a, b, cls


def edge_nodes(common):
    nodes = set()
    for c in common:
        a, b, _ = parse_edge(c)
        nodes |= {a, b}
    known = [n for n in NODE_ORDER if n in nodes]
    return known + sorted(nodes - set(known))


def edges_to_matrix(vec, common, nodes):
    """Upper triangle = intra-hemispheric edges (L-L and R-R averaged),
    lower triangle = inter-hemispheric heterotopic edges, diagonal = homotopic
    edge (L_x - R_x). Bilateral (non-lateralised) nodes go to the upper half."""
    idx = {n: i for i, n in enumerate(nodes)}
    M = np.full((len(nodes), len(nodes)), np.nan)
    for v, c in zip(vec, common):
        a, b, cls = parse_edge(c)
        i, j = idx[a], idx[b]
        if cls == "homotopic":
            M[i, i] = v
        elif cls == "inter":
            M[max(i, j), min(i, j)] = v
        else:
            M[min(i, j), max(i, j)] = v
    return M


def _entries(vecs, species_of_label, by, divergence):
    """[(name, vector, color, linestyle)] -- per species median or per BIDS."""
    out = []
    if by == "species":
        spm = species_medians(vecs, species_of_label)
        n_by = pd.Series([species_of_label.get(l, l) for l in vecs]).value_counts()
        for sp in sorted(spm, key=lambda s: (divergence.get(s, 1e9), s)):
            out.append((f"{sp} ({int(n_by.get(sp, 0))} BIDS)", spm[sp], _color(sp), "-"))
    else:
        seen = {}
        for lab in order_labels(vecs, species_of_label, divergence):
            sp = species_of_label.get(lab, lab)
            k = seen.get(sp, 0)
            seen[sp] = k + 1
            out.append((lab, vecs[lab], _color(sp), LINESTYLES[k % len(LINESTYLES)]))
    return out


def plot_fingerprints(all_profiles, sel_profiles, common, species_of_label, value_col, norm,
                      by, region_mask, out_png, divergence, sel_name, title, paper=False,
                      notes=None):
    """F1/F1c as polar fingerprints: left = every run (the F1/F1c input),
    right = final-round selection (the F1b input). Same region order on both
    sides (ascending all-runs consensus, as on F1).
    paper=True (5c): same data at print size (183 mm, 5-7 pt), panel headers
    'all runs' / 'selected' followed by notes[k] (effect and stars), legend
    without dataset counts; .png (600 dpi) + .pdf."""
    vecs_all = group_vectors(all_profiles, common, value_col, norm)
    vecs_sel = group_vectors(sel_profiles or {}, common, value_col, norm)
    if not vecs_all:
        print("  [F1_fp] nothing to plot")
        return None
    cons_all = consensus(vecs_all, species_of_label)
    order = _region_order(cons_all, common)
    labels = [short_region(common[i]) for i in order]
    mask = np.ones(len(common), bool) if region_mask is None else np.asarray(region_mask, bool)
    excluded = [not mask[i] for i in order]
    panels = ([("all runs" if paper else "All runs", vecs_all)]
              + ([("selected" if paper else sel_name, vecs_sel)] if vecs_sel else []))
    allv = [v[order] for _, vv in panels for v in vv.values()]
    rlim = _rlim(allv + [cons_all[order]])
    ctx = _print_size() if paper else contextlib.nullcontext()
    with ctx:
        fig = plt.figure(figsize=((fa.PAPER_WIDTH["double"], 98 * fa.MM) if paper
                                  else (6.6 * len(panels), 7.4)))
        handles = None
        axes_hdr = []
        for k, (name, vecs) in enumerate(panels):
            ax = fig.add_subplot(1, len(panels), k + 1, projection="polar")
            theta = _spider_axes(ax, labels, rlim, excluded)
            hs = []
            for lab, vec, col, ls in _entries(vecs, species_of_label, by, divergence):
                hs.append(_poly(ax, theta, vec[order], color=col, ls=ls, lw=_P(1.5, 0.8),
                                alpha=0.95, label=lab.split(" (")[0] if paper else lab))
            cons = consensus(vecs, species_of_label)
            hs.append(_poly(ax, theta, cons[order], color="black", lw=_P(2.8, 1.5), zorder=5,
                            label="consensus" if paper else
                            f"consensus ({len(species_medians(vecs, species_of_label))} "
                            f"species)"))
            if paper:
                note = (notes or {}).get(k)
                axes_hdr.append((ax, name, note))
            else:
                ax.set_title(name, fontsize=11, fontweight="bold", pad=22)
            handles = hs
        if not paper:
            fig.legend(handles=handles, loc="lower center", ncol=min(8, len(handles)),
                       fontsize=_P(8, 6), frameon=False, bbox_to_anchor=(0.5, 0.0))
            fig.suptitle(title, fontsize=12, fontweight="bold", y=0.99)
            fig.text(0.5, 0.935, f"{fa._norm_label(norm) or 'raw '}{value_col}; dotted ring = 0; "
                     f"regions ordered by the all-runs consensus (clockwise from the top), as "
                     f"on F1" + ("; * = region not used for the selection (invariance filter)"
                                 if any(excluded) else ""),
                     ha="center", fontsize=8, color="#555555")
            fig.subplots_adjust(top=0.84, bottom=0.13, left=0.06, right=0.94, wspace=0.42)
            _save(fig, out_png)
        else:
            fig.subplots_adjust(top=0.84, bottom=0.12, left=0.08, right=0.92, wspace=0.55)
            # header and note above the spoke labels (texts are not in the axes box)
            top = max(_texts_top(ax) for ax, _, _ in axes_hdr) if axes_hdr else 1.0
            for ax, name, note in axes_hdr:
                y = top + 0.015
                if note:
                    ax.text(0.5, y, note, transform=ax.transAxes, ha="center", va="bottom",
                            fontsize=6)
                    y += 0.055
                ax.text(0.5, y, name, transform=ax.transAxes, ha="center", va="bottom",
                        fontsize=7, fontweight="bold")
            # legend under the lowest spoke label (labels carry their channel, P80 / P20,
            # and reach below the axes box)
            low = min(_texts_bottom(ax) for ax, _, _ in axes_hdr) if axes_hdr else 0.06
            fig.legend(handles=handles, loc="upper center", ncol=min(8, len(handles)),
                       fontsize=6, frameon=False, bbox_to_anchor=(0.5, low - 0.01))
            fa.save_paper(fig, out_png)
    return dict(order=[common[i] for i in order])


def _matrix_grid(n_panels, ncols, panel=3.3, top_in=1.3):
    nrows = int(np.ceil(n_panels / ncols))
    h = panel * nrows + top_in
    fig = plt.figure(figsize=(panel * ncols + 1.6, h))
    gs = fig.add_gridspec(nrows, ncols, hspace=0.32, wspace=0.08, top=1 - top_in / h,
                          bottom=0.9 / h, left=1.0 / (panel * ncols + 1.6), right=0.9)
    return fig, gs, nrows


def _matrix_axis(ax, nodes, show_x, show_y):
    ax.set_xticks(range(len(nodes)))
    ax.set_yticks(range(len(nodes)))
    ax.set_xticklabels([short_region(n) for n in nodes] if show_x else [], rotation=90,
                       fontsize=5.5)
    ax.set_yticklabels([short_region(n) for n in nodes] if show_y else [], fontsize=5.5)
    ax.plot([-0.5, len(nodes) - 0.5], [-0.5, len(nodes) - 0.5], color="#888888", lw=0.5)


def plot_matrix_fingerprints(all_profiles, sel_profiles, common, species_of_label, value_col,
                             norm, by, out_png, divergence, sel_name, title):
    """Matrix analogue of plot_fingerprints: one region x region matrix per
    entry, one figure per pool (all runs / selected), one colour scale."""
    nodes = edge_nodes(common)
    written = []
    for tag, profiles in (("all", all_profiles), ("selected", sel_profiles)):
        vecs = group_vectors(profiles or {}, common, value_col, norm)
        if not vecs:
            continue
        ent = _entries(vecs, species_of_label, by, divergence)
        ent.append((f"consensus ({len(species_medians(vecs, species_of_label))} species)",
                    consensus(vecs, species_of_label), "black", "-"))
        mats = [edges_to_matrix(v, common, nodes) for _, v, _, _ in ent]
        vmax = np.nanpercentile(np.abs(np.concatenate([m[np.isfinite(m)] for m in mats])), 99)
        ncols = 4
        fig, gs, nrows = _matrix_grid(len(ent), ncols)
        im = None
        for k, ((name, _, col, _), M) in enumerate(zip(ent, mats)):
            ax = fig.add_subplot(gs[k // ncols, k % ncols])
            im = ax.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax, interpolation="nearest")
            ax.set_title(name, fontsize=9, fontweight="bold",
                         color="black" if col == "black" else _tcolor(
                             species_of_label.get(name, name.split(" (")[0])))
            last_row = k // ncols == nrows - 1 or k + ncols >= len(ent)
            _matrix_axis(ax, nodes, show_x=last_row, show_y=k % ncols == 0)
        cax = fig.add_axes([0.92, 0.3, 0.012, 0.4])
        fig.colorbar(im, cax=cax, label=f"{fa._norm_label(norm)}r")
        fig.suptitle(f"{title} -- {'All runs' if tag == 'all' else sel_name}\n"
                     f"upper triangle: intra-hemispheric edges (L-L and R-R averaged); lower "
                     f"triangle: inter-hemispheric (L-R); diagonal: homotopic (L_x-R_x)",
                     fontsize=10, fontweight="bold")
        path = out_png.replace(".png", f"_{tag}.png")
        _save(fig, path)
        written.append(path)
    return written


def _reference_for(vecs, species_of_label, which, lab, sp, ref_species):
    if which == "consensus":
        return consensus(vecs, species_of_label, exclude_species=sp), "LOSO consensus"
    ref = reference_vector(vecs, species_of_label, ref_species, exclude_label=lab)
    return ref, (f"other {ref_species} dataset" if sp == ref_species
                 else f"{ref_species} datasets")


def plot_distance_fingerprint_grid(vecs, dist_df, species_of_label, which, common, order,
                                   region_mask, out_png, ref_species, divergence, source,
                                   paper=False):
    """F15: one panel per dataset over its reference, closest first. Both
    vectors are z-scored, so the mean squared gap between them equals
    2 x (1 - r) exactly: the area between the curves (polar) or the energy of
    the difference matrix (matrix) IS the distance plotted on F6.
    paper=True (6b): the V12 layout -- same radar size, line widths and titles
    as the supplementary F15a, five panels per row -- without the spoke labels
    (same region order as Figure 5c, which has them), so the gaps that held
    them are closed, and without the suptitle; .png (600 dpi) + .pdf."""
    col = "r_consensus" if which == "consensus" else "r_reference"
    rows = dist_df.dropna(subset=[col]).sort_values(col, ascending=False)
    if rows.empty:
        return None
    what = ("leave-one-species-out consensus (median of the other species)"
            if which == "consensus" else f"{ref_species} reference")
    ncols = 5
    if is_edge_profile(common):
        nodes = edge_nodes(common)
        diffs = []
        for r in rows.itertuples():
            ref, ref_name = _reference_for(vecs, species_of_label, which, r.label, r.species,
                                           ref_species)
            diffs.append((r, ref_name, edges_to_matrix(
                fa._zscore_values(vecs[r.label]) - fa._zscore_values(ref), common, nodes)))
        lim = np.nanpercentile(np.abs(np.concatenate(
            [D[np.isfinite(D)] for _, _, D in diffs])), 99)
        fig, gs, nrows = _matrix_grid(len(diffs), ncols, panel=3.4, top_in=1.5)
        im = None
        for k, (r, ref_name, D) in enumerate(diffs):
            ax = fig.add_subplot(gs[k // ncols, k % ncols])
            im = ax.imshow(D, cmap="PuOr_r", vmin=-lim, vmax=lim, interpolation="nearest")
            _matrix_axis(ax, nodes, show_x=k + ncols >= len(diffs), show_y=k % ncols == 0)
            ax.set_title(f"{pretty_label(r.label)}\nr = {getattr(r, col):.2f}"
                         + ("" if paper else f"  (vs {ref_name})"),
                         fontsize=8.5, fontweight="bold", color=_tcolor(r.species))
        cax = fig.add_axes([0.92, 0.3, 0.012, 0.4])
        fig.colorbar(im, cax=cax, label="z(dataset) - z(reference)")
        if not paper:
            fig.suptitle(f"Each dataset minus its {what}, most similar first -- {source}\n"
                         f"z-scored edge profiles, one colour scale; mean of the squared "
                         f"differences = 2 x (1 - r). Triangles as on F1_fingerprint.",
                         fontsize=10.5, fontweight="bold")
        _save(fig, out_png, pdf=paper)
        return rows

    labels = [short_region(common[i]) for i in order]
    mask = np.ones(len(common), bool) if region_mask is None else np.asarray(region_mask, bool)
    excluded = [not mask[i] for i in order]
    nrows = int(np.ceil(len(rows) / ncols))
    # ONE radial range for every panel: on a polar plot the position of the
    # centre changes the apparent shape, so a per-panel range would draw the
    # same reference (e.g. the Human median behind both macaques) differently
    pairs = []
    for r in rows.itertuples():
        ref, ref_name = _reference_for(vecs, species_of_label, which, r.label, r.species,
                                       ref_species)
        pairs.append((r, ref_name, fa._zscore_values(vecs[r.label]), fa._zscore_values(ref)))
    allv = np.concatenate([np.r_[v[np.isfinite(v)], ref[np.isfinite(ref)]]
                           for _, _, v, ref in pairs])
    rlim = (min(allv.min(), 0) - 0.3, allv.max() + 0.3)
    if paper:
        # same radar size as V12 (~2.9 in); without spoke labels the gaps that
        # held them are closed
        top_in = 0.45
        h = 3.45 * nrows + top_in + 0.05
        fig = plt.figure(figsize=(3.3 * ncols, h))
        gs = fig.add_gridspec(nrows, ncols, hspace=0.22, wspace=0.12, top=1 - top_in / h,
                              bottom=0.01, left=0.01, right=0.99)
    else:
        top_in = 1.25
        h = 4.5 * nrows + top_in + 0.05
        fig = plt.figure(figsize=(4.4 * ncols, h))
        gs = fig.add_gridspec(nrows, ncols, hspace=0.38, wspace=0.55, top=1 - top_in / h,
                              bottom=0.02, left=0.03, right=0.98)
    for k, (r, ref_name, v, ref) in enumerate(pairs):
        sp = r.species
        ax = fig.add_subplot(gs[k // ncols, k % ncols], projection="polar")
        theta = _spider_axes(ax, labels, rlim, excluded, show_labels=not paper)
        for tl in ax.get_xticklabels():
            tl.set_fontsize(5.5)
        t, a = _close(theta, v[order])
        _, b = _close(theta, ref[order])
        # shaded only where the dataset exceeds its reference
        ax.fill_between(t, a, b, where=a > b, interpolate=True, color=_color(sp), alpha=0.25,
                        lw=0)
        ax.plot(t, b, color="black", ls="--", lw=1.4)
        ax.plot(t, a, color=_color(sp), lw=1.8)
        ax.set_title(f"{pretty_label(r.label)}\nr = {getattr(r, col):.2f}"
                     + ("" if paper else f"  (vs {ref_name})"),
                     fontsize=8.5, fontweight="bold", color=_tcolor(sp),
                     pad=6 if paper else 14)
    if not paper:
        fig.suptitle(
            f"Each dataset (coloured) over its {what} (black dashed), most similar first -- "
            f"{source}\nBoth curves z-scored: mean squared gap between them = 2 x (1 - r). "
            "Pearson r ignores a common offset or gain -- only the shape across regions counts",
            fontsize=10.5, fontweight="bold")
        _save(fig, out_png)
    else:
        fa.save_paper(fig, out_png)
    return rows


# ═══════════════════════════════════════════════════════════════════════════
# F16 -- dataset x dataset proximity
# ═══════════════════════════════════════════════════════════════════════════

def distance_matrix(vecs, labels):
    n = len(labels)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = 1 - _r(vecs[labels[i]], vecs[labels[j]])[0]
    return pd.DataFrame(D, index=labels, columns=labels)


def sibling_test(D, species_of_label, n_perm=10000, seed=0):
    """Do two datasets of the same species resemble each other more than
    datasets of different species? Statistic = mean d(same-species pairs) -
    mean d(cross-species pairs); null = species labels permuted across
    datasets (group sizes preserved). One-sided p (more similar)."""
    labels = list(D.index)
    sp = np.array([species_of_label.get(l, l) for l in labels])
    iu = np.triu_indices(len(labels), 1)
    d = D.values[iu]

    def stat(s):
        same = s[iu[0]] == s[iu[1]]
        if same.sum() == 0 or (~same).sum() == 0:
            return np.nan
        return d[same].mean() - d[~same].mean()

    obs = stat(sp)
    if not np.isfinite(obs):
        return None
    rng = np.random.default_rng(seed)
    null = np.array([stat(rng.permutation(sp)) for _ in range(n_perm)])
    p = (1 + int(np.sum(null <= obs + 1e-12))) / (n_perm + 1)
    same = sp[iu[0]] == sp[iu[1]]
    rows = []
    for i, lab in enumerate(labels):
        others = [(D.iloc[i, j], labels[j]) for j in range(len(labels)) if j != i]
        others.sort()
        sib = [(dd, l) for dd, l in others if species_of_label.get(l) == sp[i]]
        rows.append(dict(label=lab, species=sp[i], nearest_neighbour=others[0][1],
                         d_nearest=others[0][0],
                         sibling=sib[0][1] if sib else None,
                         d_sibling=sib[0][0] if sib else np.nan,
                         sibling_rank=(1 + [l for _, l in others].index(sib[0][1]))
                         if sib else np.nan,
                         n_others=len(others)))
    return dict(statistic=float(obs), mean_same=float(d[same].mean()),
                mean_cross=float(d[~same].mean()), n_same=int(same.sum()),
                n_cross=int((~same).sum()), p_one_sided=float(p), n_perm=n_perm,
                per_dataset=pd.DataFrame(rows))


def plot_dataset_proximity(D, dist_df, species_of_label, test, out_png, source):
    labels = list(D.index)
    n = len(labels)
    fig = plt.figure(figsize=(17, 0.55 * n + 3.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.25], wspace=0.35)
    ax = fig.add_subplot(gs[0])
    Dm = D.values.astype(float).copy()
    np.fill_diagonal(Dm, np.nan)
    cmap = plt.get_cmap("viridis_r").with_extremes(bad="#EEEEEE")
    im = ax.imshow(Dm, cmap=cmap, vmin=0, vmax=np.nanmax(Dm))
    for i in range(n):
        for j in range(n):
            if i != j:
                ax.text(j, i, f"{D.values[i, j]:.2f}", ha="center", va="center", fontsize=6.5,
                        color="white" if D.values[i, j] > 0.55 * np.nanmax(D.values) else "black")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels([pretty_label(l) for l in labels], rotation=60, ha="right", fontsize=7.5)
    ax.set_yticklabels([pretty_label(l) for l in labels], fontsize=7.5)
    for tl in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        tl.set_color(_tcolor(species_of_label.get(tl.get_text())))
    sp = [species_of_label.get(l) for l in labels]
    for i in range(1, n):
        if sp[i] != sp[i - 1]:
            ax.axhline(i - 0.5, color="white", lw=1.5)
            ax.axvline(i - 0.5, color="white", lw=1.5)
    fig.colorbar(im, ax=ax, shrink=0.75, label="1 - r between dataset profiles")
    ax.set_title("a  Distance between every pair of datasets", fontsize=10,
                 fontweight="bold", loc="left")

    ax = fig.add_subplot(gs[1])
    per = test["per_dataset"].set_index("label") if test else None
    cons = dist_df.set_index("label")["d_consensus"]
    for i, lab in enumerate(labels):
        y = n - 1 - i
        for j, other in enumerate(labels):
            if other == lab:
                continue
            same = species_of_label.get(other) == species_of_label.get(lab)
            ax.scatter(D.loc[lab, other], y, s=170 if same else 34,
                       marker="*" if same else "o", color=_color(species_of_label.get(other)),
                       edgecolor="black" if same else "none", lw=0.7, zorder=3 if same else 2,
                       alpha=1 if same else 0.8)
        if lab in cons.index and np.isfinite(cons[lab]):
            ax.scatter(cons[lab], y, marker="|", s=260, color="black", lw=2, zorder=4)
        if per is not None and lab in per.index:
            nn = per.loc[lab, "nearest_neighbour"]
            ax.scatter(D.loc[lab, nn], y, s=170, facecolor="none", edgecolor="#222222",
                       lw=1.0, zorder=5)
            sr = per.loc[lab, "sibling_rank"]
            txt = f"NN: {nn}" + (f" | sibling rank {int(sr)}/{int(per.loc[lab, 'n_others'])}"
                                 if np.isfinite(sr) else "")
            ax.text(1.01, y, txt, transform=ax.get_yaxis_transform(), fontsize=6.8,
                    va="center", color="#444444")
    ax.set_yticks(range(n))
    ax.set_yticklabels([pretty_label(l) for l in labels[::-1]], fontsize=8)
    for tl in ax.get_yticklabels():
        tl.set_color(_tcolor(species_of_label.get(tl.get_text())))
    ax.set_xlabel("1 - r to the dataset on that row")
    ax.set_xlim(0, np.nanmax(D.values) * 1.05)
    ax.grid(axis="x", color="#EEEEEE")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    t = (f"same-species pairs: mean 1 - r = {test['mean_same']:.2f} (n = {test['n_same']}) vs "
         f"cross-species {test['mean_cross']:.2f} (n = {test['n_cross']}); permutation "
         f"p = {_fmt_p(test['p_one_sided'])} (one-sided, {test['n_perm']} shuffles of species "
         f"labels across datasets)") if test else "no species with two datasets"
    ax.set_title("b  Proximity of each dataset to every other one\n" + t, fontsize=9,
                 fontweight="bold", loc="left")
    sps = sorted(set(sp), key=lambda s: list(sp).index(s))
    h = [Line2D([], [], marker="o", color=_color(s), ls="none", ms=6, label=s) for s in sps]
    h += [Line2D([], [], marker="*", color="#888888", mec="black", ls="none", ms=12,
                 label="same species, other dataset"),
          Line2D([], [], marker="|", color="black", ls="none", ms=14, mew=2,
                 label="LOSO consensus (F6a)"),
          Line2D([], [], marker="o", mfc="none", mec="#222222", ls="none", ms=11,
                 label="nearest neighbour")]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=6,
              fontsize=7, frameon=False)
    fig.suptitle(f"Dataset proximity -- {source}", fontsize=11, fontweight="bold")
    _save(fig, out_png)


# ═══════════════════════════════════════════════════════════════════════════
# F17 -- before vs after selection
# ═══════════════════════════════════════════════════════════════════════════

def _rank_perm_p(a, b, n_perm=10000, seed=0):
    """Spearman rho between two columns, p from permuting one of them
    (two-sided; exact enough at n = 10, unlike the t approximation)."""
    rho = st.spearmanr(a, b)[0]
    rng = np.random.default_rng(seed)
    null = np.array([st.spearmanr(a, rng.permutation(b))[0] for _ in range(n_perm)])
    return float(rho), (1 + int(np.sum(np.abs(null) >= abs(rho) - 1e-12))) / (n_perm + 1)


def _before_after_panel(ax, sub, col, states, final_round):
    for r in sub.itertuples():
        a, b = getattr(r, f"{col}_all"), getattr(r, f"{col}_sel")
        c = _color(r.species)
        ax.plot([0, 1], [a, b], color=c, lw=_P(1.3, 0.7), alpha=0.9)
        _state_marker(ax, 0, a, c, states.get(r.label, "unknown"), size=7)
        _state_marker(ax, 1, b, c, states.get(r.label, "unknown"), size=7)
    yv = np.r_[sub[f"{col}_sel"].values]
    both = np.r_[yv, sub[f"{col}_all"].values]
    yr = np.nanmax(both) - np.nanmin(both)
    _label_points(ax, np.ones(len(sub)), yv, list(sub["label"]),
                  [_tcolor(s) for s in sub["species"]], dx=0.08, min_gap=0.05 * (yr or 1))
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["all runs", f"selected\n(round {final_round})"])
    ax.set_xlim(-0.25, 1.75)
    ax.set_ylabel("r")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_before_after(dist_all, dist_sel, states, out_png, source, final_round,
                      ref_species="Human", paper=False):
    """F17. paper=True (5d): panel a alone (distance to the LOSO consensus, all
    runs -> selected) at print size (65 x 70 mm), no title; top right: Spearman
    rho between the two columns and its stars (p from 10000 permutations,
    n and exact p in the .txt); .png (600 dpi) + .pdf."""
    m = dist_all.merge(dist_sel, on=["label", "species", "divergence_mya"],
                       suffixes=("_all", "_sel"))
    if paper:
        rho, p = _rank_perm_p(m["r_consensus_all"].values, m["r_consensus_sel"].values)
        with _print_size():
            fig, ax = plt.subplots(figsize=(65 * fa.MM, 70 * fa.MM))
            _before_after_panel(ax, m, "r_consensus", states, final_round)
            ax.set_xticklabels(["all runs", "selected"])
            ax.set_ylabel("r with the leave-one-species-out consensus")
            ax.text(0.0, 1.01, f"\u03c1 = {rho:.2f} {fa.stars(p)}", transform=ax.transAxes,
                    ha="left", va="bottom", fontsize=6)
            fig.tight_layout()
            _place_labels(fig)
            fa.save_paper(fig, out_png)
        return dict(r_consensus=dict(spearman_rho=rho, p_perm=p, n=len(m),
                                     median_change=float((m["r_consensus_sel"]
                                                          - m["r_consensus_all"]).median()))), m
    fig, axes = plt.subplots(1, 3, figsize=(17, 6.0),
                             gridspec_kw=dict(width_ratios=[1, 1, 1.35]))
    res = {}
    for ax, col, ttl in ((axes[0], "r_consensus", "a  Similarity to the LOSO consensus"),
                         (axes[1], "r_reference", f"b  Similarity to {ref_species}")):
        sub = m if col == "r_consensus" else m[m["species"] != ref_species]
        _before_after_panel(ax, sub, col, states, final_round)
        rho, p = (st.spearmanr(sub[f"{col}_all"], sub[f"{col}_sel"])
                  if len(sub) >= 3 else (np.nan, np.nan))
        delta = (sub[f"{col}_sel"] - sub[f"{col}_all"]).median()
        res[col] = dict(spearman_rho=float(rho), spearman_p=float(p),
                        median_change=float(delta), n=len(sub))
        ax.set_title(f"{ttl}\nrank across datasets: rho = {rho:.2f}; median change "
                     f"{delta:+.3f}", fontsize=9, fontweight="bold", loc="left")
    ax = axes[2]
    out_t = {}
    for tag, df, ls, filled in (("all", dist_all, "--", False), ("selected", dist_sel, "-", True)):
        t, sp_df = f6_tests(df, ref_species)
        out_t[tag] = t
        for r in sp_df.itertuples():
            ax.scatter(r.divergence_mya, r.r_consensus, s=60, color=_color(r.species),
                       facecolor=_color(r.species) if filled else "white", lw=1.4, zorder=3)
        _fit_line(ax, sp_df, "r_consensus", "black" if filled else "#888888", ls,
                  f"{tag}: r = {t['consensus_with_ref']['r']:+.2f}, p = "
                  f"{_fmt_p(t['consensus_with_ref']['p_two'])}; without {ref_species}: r = "
                  f"{t['consensus_without_ref']['r']:+.2f}, p = "
                  f"{_fmt_p(t['consensus_without_ref']['p_two'])}")
    sp_sel = species_table(dist_sel, ref_species)
    sp_all = species_table(dist_all, ref_species).set_index("species")
    for r in sp_sel.itertuples():
        if r.species in sp_all.index:
            ax.annotate("", xy=(r.divergence_mya, r.r_consensus),
                        xytext=(r.divergence_mya, sp_all.loc[r.species, "r_consensus"]),
                        arrowprops=dict(arrowstyle="->", color=_color(r.species), lw=0.9,
                                        alpha=0.7, shrinkA=4, shrinkB=4))
    for r in sp_sel.itertuples():
        ax.annotate(r.species, (r.divergence_mya, r.r_consensus), xytext=(5, 4),
                    textcoords="offset points", fontsize=7, color=_tcolor(r.species))
    ax.set_xlabel("Divergence time from Human (Mya)")
    ax.set_ylabel("r with the LOSO consensus (species mean)")
    ax.set_title("c  F6a gradient, all runs (open) vs selected (filled)\n"
                 "exact permutation p over species", fontsize=9, fontweight="bold", loc="left")
    ax.legend(fontsize=7, loc="upper left", frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.suptitle(f"Group means before vs after selection -- {source}", fontsize=11,
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _save(fig, out_png)
    res["f6"] = out_t
    return res, m


# ═══════════════════════════════════════════════════════════════════════════
# F18 -- arousal state
# ═══════════════════════════════════════════════════════════════════════════

STATE_COLORS = {AWAKE: "#0072B2", "isoflurane": "#D55E00",
                "medetomidine+isoflurane": "#CC79A7", "anesthetized": "#D55E00",
                ANESTH: "#D55E00"}


def _state_color(state):
    return STATE_COLORS.get(state, STATE_COLORS.get(state_class(state), "#777777"))


def _run_vector(prof, common, value_col, norm):
    return np.asarray(fa._normalize_values(fa._ordered_common_values(prof, common, value_col),
                                           norm), float)


def run_level_similarity(runs_by_label, group_vecs, species_of_label, common, value_col,
                         norm, labels, kept_by_label, animal_of):
    """Per run of `labels`:
    r_consensus = r with the leave-one-species-out consensus built from the
                  GROUP profiles of the OTHER species (round-1 selection target,
                  independent of the dataset);
    r_selected  = r with the dataset's own selected reference: median of the
                  final-round KEPT runs of the same dataset, leaving out every
                  run of the same animal (so no run is compared with a reference
                  containing itself or its own animal).
    Returns (table, {label: {key: vector}})."""
    rows, vecs = [], {}
    for lab in labels:
        sp = species_of_label.get(lab)
        cons = consensus(group_vecs, species_of_label, exclude_species=sp)
        profs = runs_by_label.get(lab, {})
        v_lab = {tuple(str(x) for x in k): _run_vector(p, common, value_col, norm)
                 for k, p in profs.items()}
        vecs[lab] = v_lab
        kept = [k for k in (tuple(str(x) for x in kk) for kk in kept_by_label.get(lab, []))
                if k in v_lab]
        kept_animal = {k: animal_of(lab, k) for k in kept}
        for key, v in v_lab.items():
            a = animal_of(lab, key)
            others = [v_lab[k] for k in kept if kept_animal[k] != a]
            with np.errstate(invalid="ignore"):
                ref = np.nanmedian(np.vstack(others), axis=0) if others else None
            rows.append(dict(label=lab, species=sp, subject=str(key[0]), session=str(key[1]),
                             run=str(key[2]), key=key, kept=key in set(kept),
                             r_consensus=_r(v, cons)[0], r_selected=_r(v, ref)[0]))
    return pd.DataFrame(rows), vecs


def state_contrast(df, measure, n_perm=10000, seed=0):
    """Effect = anaesthetised - awake, in units of `measure` (r). Two designs:
    within-subject (>= 2 animals scanned in both states): the animal is the
    unit; p from sign-flipping each animal's difference (exact when 2^n <=
    n_perm). A run-level p (state labels shuffled among each animal's runs) is
    reported too: it treats runs of an animal as exchangeable, i.e. ignores the
    session structure, and is therefore liberal.
    between-subject: one mean per animal; p from permuting state across
    animals (exact when the number of assignments <= n_perm)."""
    from math import comb
    d = df[df["state_class"].isin([AWAKE, ANESTH]) & np.isfinite(df[measure])]
    if d.empty or d["state_class"].nunique() < 2:
        return None
    rng = np.random.default_rng(seed)
    per_subj = d.groupby(["subject", "state_class"])[measure].mean().unstack()
    both = per_subj.dropna()
    if len(both) >= 2:
        diffs = (both[ANESTH] - both[AWAKE]).values
        n = len(diffs)
        if 2 ** n <= n_perm:
            signs = np.array(list(itertools.product([1, -1], repeat=n)))
            exact = True
        else:
            signs = rng.choice([1, -1], size=(n_perm, n))
            exact = False
        obs = float(diffs.mean())
        null = (signs * diffs).mean(axis=1)
        p_subj = float(np.mean(np.abs(null) >= abs(obs) - 1e-12))
        if not exact:
            p_subj = (1 + int(np.sum(np.abs(null) >= abs(obs) - 1e-12))) / (n_perm + 1)
        sub = d[d["subject"].isin(both.index)]
        groups = [(g[measure].values, (g["state_class"].values == ANESTH))
                  for _, g in sub.groupby("subject")]

        def stat(gs):
            return np.mean([v[a].mean() - v[~a].mean() for v, a in gs])
        nulls = np.array([stat([(v, rng.permutation(a)) for v, a in groups])
                          for _ in range(n_perm)])
        p_run = (1 + int(np.sum(np.abs(nulls) >= abs(obs) - 1e-12))) / (n_perm + 1)
        return dict(design="within-subject", measure=measure, effect=obs,
                    n_subjects=n, n_runs=int(len(sub)), p_subject=float(p_subj),
                    p_subject_method=("exact sign-flip" if exact else "MC sign-flip")
                    + f", {n} animals",
                    p_run=float(p_run), subject_diffs=diffs)
    one_state = d.groupby("subject")["state_class"].nunique() == 1
    subj = (d[d["subject"].isin(one_state[one_state].index)]
            .groupby("subject").agg(v=(measure, "mean"), s=("state_class", "first")))
    if subj["s"].nunique() < 2:
        return None
    v, s = subj["v"].values, subj["s"].values == ANESTH
    obs = float(v[s].mean() - v[~s].mean())
    n, k = len(v), int(s.sum())
    if comb(n, k) <= n_perm:
        masks = []
        for idx in itertools.combinations(range(n), k):
            m = np.zeros(n, bool)
            m[list(idx)] = True
            masks.append(m)
        masks = np.array(masks)
        method = f"exact, {len(masks)} assignments"
    else:
        masks = np.array([rng.permutation(s) for _ in range(n_perm)])
        method = f"{n_perm} MC permutations"
    null = np.array([v[m].mean() - v[~m].mean() for m in masks])
    p = float(np.mean(np.abs(null) >= abs(obs) - 1e-12))
    if method.endswith("permutations"):
        p = (1 + int(np.sum(np.abs(null) >= abs(obs) - 1e-12))) / (n_perm + 1)
    return dict(design="between-subject", measure=measure, effect=obs,
                n_subjects=n, n_awake=int(n - k), n_anesth=k,
                n_runs=int(len(d)), p_subject=float(p),
                p_subject_method=f"state permuted across animals, {method}", p_run=np.nan)


def kept_share_test(df, n_perm=10000, seed=0):
    """Is being kept by the selection associated with the state, within animal?
    T = number of kept runs that are anaesthetised. Null: state labels permuted
    among the runs of each animal (each animal keeps its number of runs per
    state; the kept set is fixed), so animal identity cannot create the effect.
    Two-sided: |T* - E[T*]| >= |T - E[T*]|. Animals seen in one state only add
    nothing to the null variance. Works for one dataset or several pooled
    (animal ids must then be unique across datasets)."""
    d = df[df["state_class"].isin([AWAKE, ANESTH])]
    if d.empty or d["state_class"].nunique() < 2:
        return None
    an = (d["state_class"] == ANESTH).values
    kept = d["kept"].values.astype(bool)
    groups = [np.where(d["subject"].values == g)[0] for g in d["subject"].unique()]
    groups = [g for g in groups if 0 < an[g].sum() < len(g)]
    if not groups:
        return None
    rng = np.random.default_rng(seed)
    T = int(np.sum(kept & an))
    null = np.empty(n_perm)
    for i in range(n_perm):
        a = an.copy()
        for g in groups:
            a[g] = rng.permutation(a[g])
        null[i] = np.sum(kept & a)
    E = float(null.mean())
    p = (1 + int(np.sum(np.abs(null - E) >= abs(T - E) - 1e-9))) / (n_perm + 1)
    return dict(n_runs=int(len(d)), n_kept=int(kept.sum()), kept_anesth=T, expected=E,
                frac_kept_awake=float(kept[~an].mean()), frac_kept_anesth=float(kept[an].mean()),
                n_animals_both_states=len(groups), p=float(p), n_perm=n_perm)


def state_group_similarity(df_runs, vecs, group_vecs, species_of_label, value_col=None):
    """r of each state's MEDIAN run profile with the LOSO consensus (the same
    kind of number as F6a, within one dataset)."""
    rows = []
    for lab, g in df_runs.groupby("label"):
        sp = species_of_label.get(lab)
        cons = consensus(group_vecs, species_of_label, exclude_species=sp)
        for state, gs in g.groupby("state"):
            with np.errstate(invalid="ignore"):
                med = np.nanmedian(np.vstack([vecs[lab][k] for k in gs["key"]]), axis=0)
            rows.append(dict(label=lab, state=state, n_runs=len(gs),
                             r_consensus=_r(med, cons)[0]))
    return pd.DataFrame(rows)


def state_fingerprints(df_runs, vecs, group_vecs, species_of_label):
    """{label: dict(states={state: median vector}, kept=median of kept runs,
    consensus=LOSO consensus)} -- what the fingerprint column of F18 draws."""
    out = {}
    for lab, g in df_runs.groupby("label"):
        sp = species_of_label.get(lab)
        with np.errstate(invalid="ignore"):
            st_ = {s: np.nanmedian(np.vstack([vecs[lab][k] for k in gs["key"]]), axis=0)
                   for s, gs in g.groupby("state")}
            kk = list(g.loc[g["kept"], "key"])
            kept = np.nanmedian(np.vstack([vecs[lab][k] for k in kk]), axis=0) if kk else None
        out[lab] = dict(states=st_, kept=kept,
                        consensus=consensus(group_vecs, species_of_label, exclude_species=sp))
    return out


def state_region_effects(df_runs, vecs, common, n_perm=10000, seed=0):
    """Per dataset and region: mean over animals of (anaesthetised - awake) of
    the normalised regional value; p = sign-flip over animals (exact when
    possible); Benjamini-Hochberg across regions within a dataset."""
    rows = []
    rng = np.random.default_rng(seed)
    for lab, g in df_runs.groupby("label"):
        g = g[g["state_class"].isin([AWAKE, ANESTH])]
        M = np.vstack([vecs[lab][k] for k in g["key"]])
        per = {}
        for (subj, sc), idx in g.groupby(["subject", "state_class"]).groups.items():
            pos = [g.index.get_loc(i) for i in idx]
            with np.errstate(invalid="ignore"):
                per[(subj, sc)] = np.nanmean(M[pos], axis=0)
        animals = sorted({s for s, c in per if (s, AWAKE) in per and (s, ANESTH) in per})
        if len(animals) < 2:
            continue
        D = np.vstack([per[(a, ANESTH)] - per[(a, AWAKE)] for a in animals])
        n = len(animals)
        signs = (np.array(list(itertools.product([1, -1], repeat=n))) if 2 ** n <= n_perm
                 else rng.choice([1, -1], size=(n_perm, n)))
        ps = []
        for j, reg in enumerate(common):
            dj = D[:, j]
            if not np.isfinite(dj).all():
                ps.append(np.nan)
                continue
            null = (signs * dj).mean(axis=1)
            ps.append(float(np.mean(np.abs(null) >= abs(dj.mean()) - 1e-12)))
            rows.append(dict(label=lab, region=reg, effect=float(dj.mean()), n_animals=n,
                             p=ps[-1]))
        sub = [r for r in rows if r["label"] == lab]
        pv = np.array([r["p"] for r in sub])
        order = np.argsort(pv)
        q = np.empty(len(pv))
        prev = 1.0
        for rank, i in reversed(list(enumerate(order, start=1))):
            prev = min(prev, pv[i] * len(pv) / rank)
            q[i] = prev
        for r, qq in zip(sub, q):
            r["q_bh"] = float(qq)
    return pd.DataFrame(rows)


def selection_by_state(df_runs, kept_by_label=None, n_perm=10000):
    """Share of each state's runs kept by the final round, per dataset; the p
    attached (attrs['p']) is the within-animal permutation test of
    kept_share_test; Fisher's p (runs treated as independent) is kept in
    attrs['fisher_p'] for reference."""
    out = {}
    for lab, g in df_runs.groupby("label"):
        if kept_by_label is not None and "kept" not in g:
            kept = {tuple(str(x) for x in k) for k in kept_by_label.get(lab, [])}
            g = g.assign(kept=g["key"].isin(kept))
        sc = (g.groupby("state").agg(n_runs=("kept", "size"), n_kept=("kept", "sum"))
              .reset_index())
        sc["frac_kept"] = sc["n_kept"] / sc["n_runs"]
        b = g[g["state_class"].isin([AWAKE, ANESTH])]
        fp = np.nan
        if b["state_class"].nunique() == 2:
            tab = pd.crosstab(b["state_class"], b["kept"]).reindex(
                index=[AWAKE, ANESTH], columns=[True, False], fill_value=0)
            fp = float(st.fisher_exact(tab.values)[1])
        t = kept_share_test(g, n_perm=n_perm)
        sc.attrs["fisher_p"] = fp
        sc.attrs["p"] = t["p"] if t else np.nan
        sc.attrs["test"] = t
        out[lab] = sc
    return out


def plot_state_effect(df_runs, tests, grp, sel_counts, out_png, fps=None, common=None,
                      order=None, region_mask=None, paper=False):
    """F18 (paper: 6c). One row per dataset holding both states.
    1  r of each run with the LOSO consensus (dots = runs, joined markers =
       animals, star = median profile of the state);
    2  r of each run with the dataset's own selected reference (kept runs of
       the OTHER animals);
    3  share of each state's runs kept by the final selection;
    4  fingerprint: median profile of each state, the selected reference
       (dashed) and the LOSO consensus (black).
    paper=True: no suptitle or statistics boxes; dataset name heads the row,
    'n, p' heads each panel; also saved as .pdf."""
    labs = sorted(df_runs["label"].unique())
    ncol = 4 if fps else 3
    fig = plt.figure(figsize=(5.0 * ncol, 4.8 * len(labs) + 0.8))
    gs = fig.add_gridspec(len(labs), ncol, wspace=0.32, hspace=0.42,
                          width_ratios=[1, 1, 0.9, 1.1][:ncol])
    rng = np.random.default_rng(1)
    for i, lab in enumerate(labs):
        d = df_runs[df_runs["label"] == lab]
        states = sorted(d["state"].unique(), key=lambda s: (s != AWAKE, s))
        xpos = {s: k for k, s in enumerate(states)}
        col = _color(d["species"].iloc[0])
        for j, (measure, name) in enumerate((("r_consensus", "r with the LOSO consensus"),
                                             ("r_selected",
                                              "r with the selected reference\n"
                                              "(kept runs of the other animals)"))):
            ax = fig.add_subplot(gs[i, j])
            x = d["state"].map(xpos).values + rng.uniform(-0.12, 0.12, len(d))
            ax.scatter(x, d[measure], s=10, color=col, alpha=0.35, lw=0)
            sm = d.groupby(["subject", "state"])[measure].mean().reset_index()
            for subj, g in sm.groupby("subject"):
                g = g.sort_values("state", key=lambda s: s.map(xpos))
                ax.plot(g["state"].map(xpos), g[measure], color=col, lw=0.9, alpha=0.7,
                        marker="o", ms=4.5)
            if measure == "r_consensus" and grp is not None and not grp.empty:
                for r in grp[grp["label"] == lab].itertuples():
                    ax.scatter(xpos[r.state] + 0.3, r.r_consensus, marker="*", s=160,
                               color="black", zorder=5)
            t = tests.get((lab, measure))
            if t and paper:
                ax.set_title(f"n = {t['n_subjects']} animals, p = {_fmt_p(t['p_subject'])}",
                             loc="right", fontsize=10)
            elif t:
                txt = (f"anaesth - awake = {t['effect']:+.3f}\n"
                       f"p = {_fmt_p(t['p_subject'])} ({t['p_subject_method']})")
                if np.isfinite(t.get("p_run", np.nan)):
                    txt += f"\nrun-level p = {_fmt_p(t['p_run'])} (liberal)"
                ax.text(0.02, 0.02, txt, transform=ax.transAxes, va="bottom", fontsize=7,
                        bbox=dict(fc="white", ec="#CCCCCC", lw=0.6))
            ax.set_xticks(range(len(states)))
            ax.set_xticklabels([f"{s}\n(n={int((d['state'] == s).sum())} runs)"
                                for s in states], fontsize=8)
            ax.set_xlim(-0.5, len(states) - 0.2)
            ax.set_ylabel(name)
            if j == 0:
                ax.set_title(pretty_label(lab), fontsize=10, fontweight="bold", loc="left",
                             color=_tcolor(d["species"].iloc[0]))
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        ax = fig.add_subplot(gs[i, 2])
        sc = sel_counts.get(lab)
        if sc is not None and not sc.empty:
            p = sc.attrs.get("p", np.nan)
            sc = sc.set_index("state").reindex(states).reset_index()
            ax.bar(range(len(sc)), 100 * sc["frac_kept"], color=[_state_color(s)
                                                                 for s in sc["state"]],
                   alpha=0.8)
            for k, r in enumerate(sc.itertuples()):
                ax.text(k, 100 * r.frac_kept + 0.5, f"{int(r.n_kept)}/{int(r.n_runs)}",
                        ha="center", va="bottom", fontsize=8)
            ax.set_xticks(range(len(sc)))
            ax.set_xticklabels(sc["state"], fontsize=8)
            ax.set_title(f"n = {int(sc['n_runs'].sum())} runs, p = {_fmt_p(p)}"
                         + ("" if paper else "\n(state permuted within animal)"),
                         loc="right", fontsize=10 if paper else 9)
            ax.set_ylabel("% of runs kept")
            ax.set_ylim(0, max(100 * sc["frac_kept"].max() * 1.3, 10))
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        else:
            ax.axis("off")
        if fps and lab in fps:
            f = fps[lab]
            ax = fig.add_subplot(gs[i, 3], projection="polar")
            o = np.arange(len(common)) if order is None else np.asarray(order)
            labels = [short_region(common[k]) for k in o]
            mask = (np.ones(len(common), bool) if region_mask is None
                    else np.asarray(region_mask, bool))
            arrays = [v[o] for v in f["states"].values()] + [
                a[o] for a in (f["kept"], f["consensus"]) if a is not None]
            theta = _spider_axes(ax, labels, _rlim(arrays), [not mask[k] for k in o])
            for tl in ax.get_xticklabels():
                tl.set_fontsize(6)
            hs = []
            for s in states:
                if s in f["states"]:
                    hs.append(_poly(ax, theta, f["states"][s][o], color=_state_color(s), lw=1.6,
                                    label=s))
            if f["kept"] is not None:
                hs.append(_poly(ax, theta, f["kept"][o], color="#555555", lw=1.4, ls="--",
                                label="selected reference"))
            if f["consensus"] is not None:
                hs.append(_poly(ax, theta, f["consensus"][o], color="black", lw=2.2,
                                label="LOSO consensus"))
            ax.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=2,
                      fontsize=7, frameon=False)
    handles = [Line2D([], [], marker="o", color="#666666", ls="-", ms=4.5,
                      label="animal mean (lines join the same animal)"),
               Line2D([], [], marker="o", color="#666666", alpha=0.4, ls="none", ms=4,
                      label="single run"),
               Line2D([], [], marker="*", color="black", ls="none", ms=12,
                      label="median profile of the state (comparable to Fig. 6a)")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8, frameon=False)
    if not paper:
        fig.suptitle("Arousal state inside the datasets that contain both states",
                     fontsize=11, fontweight="bold")
    fig.subplots_adjust(left=0.05, right=0.98, top=0.93 if not paper else 0.96,
                        bottom=0.12 if len(labs) < 3 else 0.08)
    _save(fig, out_png, pdf=paper)


# ═══════════════════════════════════════════════════════════════════════════
# F8b by arousal state (paper 6c)
# ═══════════════════════════════════════════════════════════════════════════

def annotate_states(dist_df, bids_configs):
    """Adds state, state_class and animal to a per-run table (label, subject,
    session, run, seed) using each dataset's BIDS_CONFIGS 'state' rule."""
    cfg_by = {c["label"]: c for c in bids_configs}
    res = {lab: StateResolver(cfg_by.get(lab, {})) for lab in dist_df["label"].unique()}
    d = dist_df.copy()
    keys = [tuple(str(x) for x in k) for k in
            zip(d["subject"], d["session"], d["run"], d.get("seed", d["run"]))]
    d["state"] = [res[l].state(k) for l, k in zip(d["label"], keys)]
    d["state_class"] = d["state"].map(state_class)
    d["animal"] = [res[l].animal(k) for l, k in zip(d["label"], keys)]
    return d


def state_r_tests(d, value_col="pearson_r", n_perm=10000):
    """Per dataset holding both states: within-animal contrast of value_col
    (anaesthetised - awake). Uses state_contrast (sign-flip over animals; run-
    level permutation reported as liberal). Returns {label: test}."""
    out = {}
    for lab, g in d.groupby("label"):
        if (g["state_class"] == AWAKE).any() and (g["state_class"] == ANESTH).any():
            g = g.copy()
            g["subject"] = g["animal"].values      # the unit of the test is the animal
            t = state_contrast(g, value_col, n_perm=n_perm)
            if t:
                out[lab] = t
    return out


def plot_f8b_by_state(d, tests, out_png, value_col="pearson_r", title=None, paper=False):
    """F8b with the arousal state of every run: one column per dataset (violin
    + jittered runs, crosses = kept, quartile ticks), as F8b; anaesthetised
    runs drawn as triangles, awake runs as circles. Datasets holding both
    states get two sub-columns (awake left, anaesthetised right) with their
    own violin, joined by a bracket with the stars of the within-animal test
    (n, effect and exact p in the .txt).
    paper=True (6c): print size (183 mm), no title, no n or p on the figure
    (they are in the .txt), brackets and stars only."""
    labels = sorted(d["label"].unique())
    # species colours, as on every other panel (two datasets of a species share it;
    # the x labels name the dataset)
    if "species" in d and d["species"].notna().all():
        sp_of = d.groupby("label")["species"].first()
        colors = {lab: _color(sp_of[lab]) for lab in labels}
    else:
        colors = fa._palette_color_map(labels)
    ctx = fa.paper_style() if paper else plt.rc_context({})
    with ctx:
        if paper:
            fig, ax = plt.subplots(figsize=(fa.PAPER_WIDTH["double"], 62 * fa.MM))
        else:
            fig, ax = plt.subplots(figsize=(max(6.0, 1.25 * len(labels)), 5.2))
        rng = np.random.default_rng(42)
        vals_all = d[value_col].values.astype(float) if len(d) else np.array([1.0])
        ymax = float(np.nanmax(vals_all))
        ymin = float(np.nanmin(vals_all))
        span = (ymax - ymin) or 0.1
        ms_o, ms_t, ms_x = (6, 8, 16) if paper else (24, 30, 70)
        lw_q = (0.8, 0.5) if paper else (1.4, 0.8)
        for xi, lab in enumerate(labels):
            g = d[d["label"] == lab].dropna(subset=[value_col])
            if g.empty:
                continue
            mixed = lab in tests
            groups = ([(AWAKE, -0.2, 0.34), (ANESTH, 0.2, 0.34)] if mixed else
                      [(None, 0.0, 0.7)])
            for sc_, dx, w in groups:
                sub = g if sc_ is None else g[g["state_class"] == sc_]
                if sub.empty:
                    continue
                vals = sub[value_col].values
                if len(vals) >= 2:
                    parts = ax.violinplot(vals, positions=[xi + dx], widths=w,
                                          showextrema=False)
                    for pc in parts["bodies"]:
                        pc.set_facecolor("#dddddd")
                        pc.set_edgecolor("none")
                        pc.set_alpha(0.6)
                jit = rng.uniform(-w / 4.5, w / 4.5, len(vals))
                an = (sub["state_class"] == ANESTH).values
                for mask, mk in ((~an, "o"), (an, "^")):
                    if mask.any():
                        ax.scatter((xi + dx + jit)[mask], vals[mask],
                                   s=ms_o if mk == "o" else ms_t, marker=mk, color=colors[lab],
                                   edgecolors="black", lw=0.2 if paper else 0.3, alpha=0.85,
                                   zorder=3)
                kept = (sub["is_kept"].fillna(False).values.astype(bool) if "is_kept" in sub
                        else np.zeros(len(sub), bool))
                if kept.any():
                    ax.scatter((xi + dx + jit)[kept], vals[kept], marker="x", s=ms_x,
                               color="black", linewidths=0.7 if paper else 1.6, zorder=5)
                for q in (25, 50, 75):
                    ax.hlines(np.percentile(vals, q), xi + dx - w / 3, xi + dx + w / 3,
                              color="black", lw=lw_q[0] if q == 50 else lw_q[1], zorder=4)
            if mixed:
                t = tests[lab]
                top = float(np.nanmax(g[value_col].values))
                y = top + 0.05 * span
                fa.bracket(ax, xi - 0.2, xi + 0.2, y, t["p_subject"], dy=0.03 * span,
                           fs=6 if paper else 9, lw=0.5 if paper else 0.9)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels([pretty_label(l) for l in labels], rotation=35, ha="right",
                           fontsize=None if paper else 8)
        ax.set_ylabel(("r with the selection criterion" if paper else
                       "Pearson r to BIDS reference") if value_col == "pearson_r"
                      else value_col, fontsize=None if paper else 10)
        ax.set_ylim(top=ymax + (0.16 if paper else 0.2) * (span if paper else abs(ymax) + 0.1))
        ax.set_xlim(-0.6, len(labels) - 0.4)
        if not paper:
            ax.set_title(title or "Per-run correlation to the selection criterion, by arousal "
                         "state", fontsize=11, fontweight="bold")
        ms = 3.5 if paper else 6
        handles = [Line2D([], [], marker="o", color="#777777", mec="black", mew=0.3, lw=0,
                          ms=ms, label="awake"),
                   Line2D([], [], marker="^", color="#777777", mec="black", mew=0.3, lw=0,
                          ms=ms, label="anaesthetised"),
                   Line2D([], [], marker="x", color="black", lw=0, ms=ms + 1,
                          mew=0.8 if paper else 1.6, label="selected run")]
        if paper:
            ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.0),
                      frameon=False, ncol=3, borderaxespad=0.3)
        else:
            ax.legend(handles=handles, fontsize=8, loc="lower left", frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        if paper:
            fa.save_paper(fig, out_png)
        else:
            _save(fig, out_png)


# ═══════════════════════════════════════════════════════════════════════════
# F20 -- does the selection smooth inter-run, inter-session or inter-animal
# variability?
# ═══════════════════════════════════════════════════════════════════════════

LEVELS = ("inter-run", "inter-session", "inter-animal")


def _level_stats(R, C, idx):
    """Median r per level over the pairs inside idx, and the two gaps."""
    idx = np.asarray(sorted(idx))
    if len(idx) < 3:
        return None
    sub_R = R[np.ix_(idx, idx)]
    sub_C = C[np.ix_(idx, idx)]
    iu = np.triu_indices(len(idx), 1)
    r, c = sub_R[iu], sub_C[iu]
    ok = np.isfinite(r)
    r, c = r[ok], c[ok]
    out = {}
    for k, name in enumerate(LEVELS):
        v = r[c == k]
        out[name] = float(np.median(v)) if len(v) else np.nan
        out[f"n_{name}"] = int(len(v))
    same = r[c <= 1]
    out["same_animal"] = float(np.median(same)) if len(same) else np.nan
    out["gap_animal"] = out["same_animal"] - out["inter-animal"]
    out["gap_session"] = out["inter-run"] - out["inter-session"]
    return out


def variability_levels(runs_by_label, selections_final, resolvers, common, value_col, norm,
                       n_draw=2000, seed=0):
    """Per dataset with several runs per animal: pairwise r between run
    profiles, classed as inter-run (same animal, same session), inter-session
    (same animal, other session) and inter-animal. Statistics on every run, on
    the runs kept by the final round, and on n_draw random subsets of the same
    size (the null of 'a selection that ignores animals and sessions').
    gap_animal = median r same animal - median r different animals: if the
    selection erases the identity of the animal, it falls below the null.
    Returns (table, {label: dict(null=..., obs_sel=..., obs_all=...)})."""
    rng = np.random.default_rng(seed)
    rows, draws = [], {}
    for lab, profs in runs_by_label.items():
        keys = sorted(profs, key=lambda k: tuple(str(x) for x in k))
        if len(keys) < 6:
            continue
        V = np.vstack([_run_vector(profs[k], common, value_col, norm) for k in keys])
        ok = np.isfinite(V).all(axis=0)
        if ok.sum() < 3:
            continue
        V = V[:, ok]
        R = np.corrcoef(V)
        res = resolvers.get(lab)
        skeys = [tuple(str(x) for x in k) for k in keys]
        animal = np.array([res.animal(k) if res else k[0] for k in skeys])
        session = np.array([f"{k[0]}|{k[1]}" for k in skeys])
        C = np.where(animal[:, None] != animal[None, :], 2,
                     np.where(session[:, None] != session[None, :], 1, 0))
        all_idx = np.arange(len(keys))
        obs_all = _level_stats(R, C, all_idx)
        if obs_all is None or obs_all["n_inter-run"] + obs_all["n_inter-session"] < 3 \
                or len(set(animal)) < 2:
            continue                      # one run per animal: nothing to test
        kept = {tuple(str(x) for x in k) for k in selections_final.get(lab, [])}
        sel_idx = [i for i, k in enumerate(skeys) if k in kept]
        obs_sel = _level_stats(R, C, sel_idx)
        if obs_sel is None:
            continue
        null = [_level_stats(R, C, rng.choice(len(keys), len(sel_idx), replace=False))
                for _ in range(n_draw)]
        null = pd.DataFrame([n for n in null if n is not None])
        row = dict(label=lab, n_runs=len(keys), n_kept=len(sel_idx),
                   n_animals=len(set(animal)), n_sessions=len(set(session)),
                   n_animals_kept=len(set(animal[sel_idx])))
        for k in LEVELS + ("same_animal", "gap_animal", "gap_session"):
            row[f"{k}_all"] = obs_all[k]
            row[f"{k}_selected"] = obs_sel[k]
            nv = null[k].dropna().values
            row[f"{k}_null_median"] = float(np.median(nv)) if len(nv) else np.nan
            if k.startswith("gap") and len(nv) and np.isfinite(obs_sel[k]):
                c = np.median(nv)
                row[f"p_{k}"] = (1 + int(np.sum(np.abs(nv - c) >= abs(obs_sel[k] - c)
                                                 - 1e-12))) / (len(nv) + 1)
        for k in LEVELS:
            row[f"n_pairs_{k}_selected"] = obs_sel[f"n_{k}"]
        rows.append(row)
        draws[lab] = dict(null=null, obs_sel=obs_sel, obs_all=obs_all)
    return pd.DataFrame(rows), draws


def pooled_gap_test(draws, key="gap_animal"):
    """Mean over datasets of the z-score of the selected gap within its own
    null; null of that mean = the same mean over the draws (draw i of every
    dataset taken together)."""
    labs = [l for l, d in draws.items() if np.isfinite(d["obs_sel"].get(key, np.nan))
            and d["null"][key].notna().sum() > 10]
    if len(labs) < 2:
        return None
    n = min(len(draws[l]["null"]) for l in labs)
    Z, zobs = [], []
    for l in labs:
        nv = draws[l]["null"][key].values[:n]
        mu, sd = np.nanmean(nv), np.nanstd(nv) or np.nan
        zobs.append((draws[l]["obs_sel"][key] - mu) / sd)
        Z.append((nv - mu) / sd)
    Z = np.nanmean(np.vstack(Z), axis=0)
    obs = float(np.nanmean(zobs))
    p = (1 + int(np.sum(np.abs(Z) >= abs(obs) - 1e-12))) / (len(Z) + 1)
    return dict(key=key, mean_z=obs, p=float(p), n_datasets=len(labs), n_draws=n)


def plot_variability_levels(table, draws, species_of_label, out_png, pooled=None):
    labs = list(table["label"])
    if not labs:
        return None
    ncol = min(5, len(labs))
    nrow = int(np.ceil(len(labs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 3.6 * nrow + 0.8),
                             squeeze=False, sharey=True)
    x = np.arange(len(LEVELS))
    for k, lab in enumerate(labs):
        ax = axes[k // ncol, k % ncol]
        d = draws[lab]
        col = _color(species_of_label.get(lab))
        lo = [np.nanpercentile(d["null"][lv], 2.5) for lv in LEVELS]
        hi = [np.nanpercentile(d["null"][lv], 97.5) for lv in LEVELS]
        md = [np.nanmedian(d["null"][lv]) for lv in LEVELS]
        ax.vlines(x, lo, hi, color="#BBBBBB", lw=6, alpha=0.8, zorder=1)
        ax.scatter(x, md, color="#888888", marker="_", s=120, zorder=2)
        ax.plot(x, [d["obs_all"][lv] for lv in LEVELS], color="#555555", marker="s", ms=5,
                lw=1.0, ls="--", zorder=3)
        ax.plot(x, [d["obs_sel"][lv] for lv in LEVELS], color=col, marker="o", ms=7, lw=1.8,
                zorder=4)
        row = table[table["label"] == lab].iloc[0]
        ax.set_title(pretty_label(lab), loc="left", fontsize=9, fontweight="bold",
                     color=_tcolor(species_of_label.get(lab)))
        ax.set_title(f"n = {int(row['n_kept'])} runs, p = {_fmt_p(row.get('p_gap_animal', np.nan))}",
                     loc="right", fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels(LEVELS, fontsize=8, rotation=20)
        if k % ncol == 0:
            ax.set_ylabel("median run-to-run profile r")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    for k in range(len(labs), nrow * ncol):
        axes[k // ncol, k % ncol].axis("off")
    handles = [Line2D([], [], color="#555555", marker="s", ls="--", label="all runs"),
               Line2D([], [], color="#444444", marker="o", lw=1.8,
                      label="selected (final round)"),
               Line2D([], [], color="#BBBBBB", lw=6,
                      label="random subsets of the same size (95%)")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8, frameon=False)
    sup = "Does the selection smooth inter-run, inter-session or inter-animal variability?"
    if pooled:
        sup += (f"  -- pooled over {pooled['n_datasets']} datasets: mean z of the animal gap "
                f"= {pooled['mean_z']:+.2f}, p = {_fmt_p(pooled['p'])}")
    fig.suptitle(sup, fontsize=10, fontweight="bold")
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    _save(fig, out_png)


# ═══════════════════════════════════════════════════════════════════════════
# F22 -- reliability of the profile, by profile method
# ═══════════════════════════════════════════════════════════════════════════

METHOD_LABELS = {"median_r": "median r", "fraction_posneg": "fraction pos / neg",
                 "percentiles": "percentiles"}


def discriminability(D, animal):
    """Bridgeford et al. 2021 (PLoS Comput Biol 17:e1009279): probability that
    a run is closer to another run of the SAME animal than to a run of a
    different animal. D = distance matrix (runs x runs), animal = label per
    run. Mean over every ordered same-animal pair (i, j) of the fraction of
    runs k of other animals with D[i, k] > D[i, j] (ties count 1/2). 0.5 =
    chance, 1 = every run identifies its animal."""
    animal = np.asarray(animal)
    num, den = 0.0, 0
    for i in range(len(animal)):
        same = np.flatnonzero((animal == animal[i]) & (np.arange(len(animal)) != i))
        other = np.flatnonzero(animal != animal[i])
        if not same.size or not other.size:
            continue
        dij = D[i, same][:, None]
        dik = D[i, other][None, :]
        num += float(np.sum((dik > dij).mean(axis=1) + 0.5 * (dik == dij).mean(axis=1)))
        den += same.size
    return num / den if den else np.nan


def _method_matrix(sri, profs, keys, common, method, norm, q):
    """Runs x values for one profile method, computed from the raw profiles
    (which carry the voxel sketch)."""
    if method == "median_r":
        names, vcol = list(common), "median_r"
        rows = [fa._ordered_common_values(profs[k], names, vcol) for k in keys]
    else:
        names, vcol = sri.channel_names(common, method, q), sri.profile_value_col(method, q)
        rows = [fa._ordered_common_values(sri.to_channel_profile(profs[k], method,
                                                                 reference=common, q=q),
                                          names, vcol) for k in keys]
    return np.vstack([fa._normalize_values(np.asarray(r, float), norm) for r in rows])


def profile_reliability(sri, runs_by_label, common, resolvers, norm, q=0.2,
                        methods=("median_r", "fraction_posneg", "percentiles"),
                        n_perm=500, seed=0):
    """Per dataset with several runs per animal and per profile method, on
    EVERY run (before any selection): median r inter-run (same session),
    same animal, different animals, the animal gap, and the discriminability
    (with its permutation p against chance: animal labels shuffled across
    runs). Returns a long DataFrame (one row per dataset x method)."""
    rng = np.random.default_rng(seed)
    rows = []
    for lab, profs in runs_by_label.items():
        keys = sorted(profs, key=lambda k: tuple(str(x) for x in k))
        if len(keys) < 6:
            continue
        skeys = [tuple(str(x) for x in k) for k in keys]
        res = resolvers.get(lab)
        animal = np.array([res.animal(k) if res else k[0] for k in skeys])
        session = np.array([f"{k[0]}|{k[1]}" for k in skeys])
        if len(set(animal)) < 2 or len(set(animal)) == len(animal):
            continue                                  # no repeated runs per animal
        C = np.where(animal[:, None] != animal[None, :], 2,
                     np.where(session[:, None] != session[None, :], 1, 0))
        iu = np.triu_indices(len(keys), 1)
        perms = [rng.permutation(animal) for _ in range(n_perm)]
        for method in methods:
            try:
                V = _method_matrix(sri, profs, keys, common, method, norm, q)
            except ValueError as e:
                print(f"  [F22] {lab} / {method}: {e}")
                continue
            ok = np.isfinite(V).all(axis=0)
            if ok.sum() < 3:
                continue
            R = np.corrcoef(V[:, ok])
            r, c = R[iu], C[iu]
            med = {k: float(np.median(r[c == v])) if np.any(c == v) else np.nan
                   for k, v in (("inter_run", 0), ("inter_session", 1), ("inter_animal", 2))}
            same = r[c <= 1]
            dist = 1.0 - R
            d_obs = discriminability(dist, animal)
            null = np.array([discriminability(dist, a) for a in perms])
            rows.append(dict(
                label=lab, method=method, n_runs=len(keys), n_animals=len(set(animal)),
                n_values=int(ok.sum()), r_inter_run=med["inter_run"],
                r_inter_session=med["inter_session"], r_inter_animal=med["inter_animal"],
                r_same_animal=float(np.median(same)) if same.size else np.nan,
                gap_animal=(float(np.median(same)) - med["inter_animal"]) if same.size
                else np.nan,
                discriminability=d_obs,
                p_discriminability=(1 + int(np.sum(null >= d_obs - 1e-12))) / (n_perm + 1)))
    return pd.DataFrame(rows)


def compare_methods(table, measure="discriminability", reference="median_r"):
    """Paired comparison of each method with `reference` over datasets: mean
    difference and exact two-sided sign-flip p (2^n assignments)."""
    out = {}
    if table.empty:
        return out
    w = table.pivot(index="label", columns="method", values=measure)
    if reference not in w:
        return out
    for m in w.columns:
        if m == reference:
            continue
        d = (w[m] - w[reference]).dropna().values
        if len(d) < 2:
            continue
        signs = np.array(list(itertools.product((1, -1), repeat=len(d))))
        null = np.abs((signs * d).mean(axis=1))
        out[m] = dict(mean_diff=float(d.mean()), n=len(d),
                      p=float(np.mean(null >= abs(d.mean()) - 1e-12)),
                      n_better=int(np.sum(d > 0)))
    return out


def plot_profile_reliability(table, species_of_label, out_png, current_method=None):
    """F22: discriminability and animal gap per dataset, one column per method."""
    if table.empty:
        return None
    methods = [m for m in METHOD_LABELS if m in set(table["method"])]
    x = np.arange(len(methods))
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    for ax, measure, ylab, ref_line in (
            (axes[0], "discriminability", "discriminability (0.5 = chance)", 0.5),
            (axes[1], "gap_animal", "median r same animal - different animals", 0.0)):
        w = table.pivot(index="label", columns="method", values=measure).reindex(columns=methods)
        for lab, row in w.iterrows():
            col = _color(species_of_label.get(lab, lab))
            ax.plot(x, row.values, "-o", color=col, lw=1.4, ms=6, alpha=0.9)
            ax.text(x[-1] + 0.08, row.values[-1], pretty_label(lab), fontsize=7, color=col,
                    va="center")
        ax.plot(x, w.mean(axis=0).values, "-s", color="black", lw=2.6, ms=8, zorder=5,
                label="mean over datasets")
        ax.axhline(ref_line, color="#999999", lw=0.8, ls=":")
        tests = compare_methods(table, measure)
        ax.set_xticks(x)
        ax.set_xticklabels([METHOD_LABELS[m] + ("\n(this run)" if m == current_method else "")
                            for m in methods])
        ax.set_xlim(-0.3, len(methods) - 0.2)
        ax.set_ylabel(ylab)
        ax.set_title("  ".join(f"{METHOD_LABELS[m]} vs median r: {t['n_better']}/{t['n']}, "
                               f"p = {_fmt_p(t['p'])}" for m, t in tests.items()),
                     loc="right", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.legend(handles=[Line2D([], [], color="black", marker="s", lw=2.6, ms=8,
                               label="mean over datasets"),
                        Line2D([], [], color="#999999", lw=0.8, ls=":",
                               label="chance (left) / no animal effect (right)")],
               loc="lower center", ncol=2, frameon=False, fontsize=8)
    fig.suptitle("F22 -- How reliable is the profile? Every run, before selection; datasets "
                 "with several runs per animal", fontsize=10, fontweight="bold")
    fig.tight_layout(rect=(0, 0.06, 1, 0.93))
    _save(fig, out_png)
    return True


# ═══════════════════════════════════════════════════════════════════════════
# Entry point called by fig5_pipeline
# ═══════════════════════════════════════════════════════════════════════════

def make_fingerprint_figures(*, group_all, group_all_source, ba_all, ba_sel, ba_source,
                             f1_all, f1_sel, common, region_mask, species_of_label,
                             bids_configs, runs_by_label, selections_final, final_round, cfg,
                             fig_dir, table_dir, divergence, ref_species="Human",
                             n_perm=10000):
    """
    group_all        : {label: (species, profile_df)} -- what F6 uses (all units,
                       AFNI group maps when they exist).
    ba_all / ba_sel  : {label: (species, profile_df)} -- group mean of every unit /
                       of the final-round selection, SAME estimator on both sides
                       (F17); ba_source names it.
    f1_all / f1_sel  : {label: profile_df} -- exactly what F1 / F1b plot.
    runs_by_label         : {label: {run_key: profile_df}} (state test).
    selections_final      : {label: [unit keys]} -- final-round kept units.
    divergence            : {species: Mya from Human} (fa.DIVERGENCE_TIME_MYA at call time).
    """
    value_col, norm = cfg["VALUE_COL"], cfg["ZSCORE"]
    level = cfg.get("ATLAS_LEVEL")
    matrix_mode = cfg.get("PROFILE_KIND", "seed") == "matrix" or is_edge_profile(common)
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(table_dir, exist_ok=True)
    cfg_by_label = {c["label"]: c for c in bids_configs}
    out = {}

    # ---- arousal state of every run --------------------------------------
    resolvers = {lab: StateResolver(cfg_by_label.get(lab, {})) for lab in runs_by_label}
    run_states = {lab: {tuple(str(x) for x in k): resolvers[lab].state(k) for k in runs}
                  for lab, runs in runs_by_label.items()}
    ds_state = {lab: dataset_state(resolvers[lab], list(runs_by_label[lab]))
                for lab in runs_by_label}
    st_rows = []
    for lab, d in run_states.items():
        vc = pd.Series(list(d.values())).value_counts()
        st_rows.append(dict(label=lab, species=species_of_label.get(lab),
                            dataset_state=ds_state[lab],
                            rule=repr(cfg_by_label.get(lab, {}).get("state")),
                            **{f"n_runs_{k}": int(v) for k, v in vc.items()}))
    states_df = pd.DataFrame(st_rows)
    print("  [states] " + "; ".join(f"{r.label}={r.dataset_state}"
                                    for r in states_df.itertuples()))
    unknown = states_df.loc[states_df["dataset_state"] == "unknown", "label"].tolist()
    if unknown:
        print(f"  [states] [WARN] no arousal state for {unknown} -- add 'state' to their "
              f"BIDS_CONFIGS entry")

    # ---- F6ab -------------------------------------------------------------
    vecs_all = group_vectors(group_all, common, value_col, norm)
    dist_all = dataset_distances(vecs_all, species_of_label, divergence, ref_species)
    if len(dist_all) >= 3:
        tests, sp_df = plot_f6_two_references(
            dist_all, ds_state, opj(fig_dir, "F6ab_distance_two_references.png"),
            group_all_source, ref_species, level)
        out["f6_tests"] = tests
        t = tests
        _legend(opj(fig_dir, "F6ab_distance_two_references.txt"),
                "F6ab -- Similarity to two references vs divergence time",
                f"(a) r between each dataset's group profile ({group_all_source}) and "
                f"the leave-one-species-out (LOSO) consensus: median, across the OTHER "
                f"species, of the per-species median profile -- the quantity of the original "
                f"F6 (shown as r instead of 1 - r; same p, opposite sign of the slope). "
                f"(b) r between each dataset and the median of the {ref_species} "
                f"datasets; the {ref_species} datasets are compared with each other and "
                f"define the grey band (their r). Marker fill = arousal state (filled awake, open "
                f"anaesthetised, half both, diamond unknown); large rings = species means, "
                f"the points that are fitted.",
                f"{len(common)} common regions, profiles {fa._norm_label(norm) or 'unscaled'}"
                f"(Pearson r is insensitive to it). Fits and p on species means (two datasets "
                f"of a species are not independent). p = exact permutation over all n! "
                f"assignments of divergence times to species (n = 7: 5040, n = 6: 720). "
                f"Divergence times: {cfg.get('TREE_PATH') or 'hand-entered'}.\n"
                f"Correlation of r with divergence. (a) with {ref_species}: r = "
                f"{t['consensus_with_ref']['r']:+.3f}, "
                f"p = {_fmt_p(t['consensus_with_ref']['p_two'])}; without: r = "
                f"{t['consensus_without_ref']['r']:+.3f}, p = "
                f"{_fmt_p(t['consensus_without_ref']['p_two'])}. (b) r = "
                f"{t['reference']['r']:+.3f}, p = {_fmt_p(t['reference']['p_two'])}.",
                "(a) measures HOW MUCH a dataset shares the pattern of the other species, "
                "not in which direction it departs: two species can have the same r while "
                "departing in opposite directions. The consensus is a median across the "
                "SAMPLED species: it sits wherever most of them agree, so it is not an "
                "ancestral state and it depends on which species were sampled. A positive "
                f"slope means the species closest to {ref_species} share the LEAST with the "
                "others (the consensus resembles the distant species); a negative slope the "
                "reverse. Report the fit without "
                f"{ref_species} next to the full fit: {ref_species} sits alone at 0 Mya and "
                "has the largest leverage. (b) asks the directional "
                f"question 'how {ref_species}-like is each species'; it is dominated by "
                f"dataset noise when the r between the two {ref_species} datasets is low. "
                "Divergence time "
                "is collinear with brain size, arousal state (most awake datasets are "
                "primates) and effective resolution in voxels per region: neither panel "
                "can separate them with 7 species, and the Dog (large brain, distant from "
                "Human) is the only species that pulls them apart. The within-species state "
                "contrast (F18) is the only direct estimate of the state effect."
                + (" SOURCE WARNING: these are the KEPT maps of the final selection round "
                   "(PHYLOGENY_USE_FINAL_ROUND=True); the unselected values are in F17 and "
                   "fingerprint_tables.xlsx/before_after." if "KEPT" in group_all_source
                   else ""))
        out["f6_species"] = sp_df

    # ---- F1 / F1c as fingerprints -----------------------------------------
    sel_name = f"Selected (final round = {final_round})"
    if matrix_mode:
        for by, name in (("species", "F1_fingerprint_species"), ("bids", "F1c_fingerprint_bids")):
            plot_matrix_fingerprints(
                f1_all, f1_sel, common, species_of_label, value_col, norm, by,
                opj(fig_dir, f"{name}.png"), divergence, sel_name,
                f"{'F1' if by == 'species' else 'F1c'} as connectivity matrices, one per "
                f"{'species' if by == 'species' else 'BIDS'}")
        order = np.arange(len(common))
    else:
        r1 = plot_fingerprints(f1_all, f1_sel, common, species_of_label, value_col, norm,
                               "species", region_mask, opj(fig_dir, "F1_fingerprint_species.png"),
                               divergence, sel_name,
                               f"F1 as a connectivity fingerprint: one line per species "
                               f"(atlas level {level})")
        plot_fingerprints(f1_all, f1_sel, common, species_of_label, value_col, norm,
                          "bids", region_mask, opj(fig_dir, "F1c_fingerprint_bids.png"),
                          divergence, sel_name,
                          f"F1c as a connectivity fingerprint: one line per BIDS "
                          f"(atlas level {level})")
        order = np.array([common.index(c) for c in r1["order"]]) if r1 else np.arange(len(common))
    for name, by in (("F1_fingerprint_species", "species"), ("F1c_fingerprint_bids", "BIDS")):
        _legend(opj(fig_dir, f"{name}.txt"),
                f"{name} -- {'F1' if by == 'species' else 'F1c'} drawn as connectivity "
                f"fingerprints",
                f"Left: one {by} line per {'species (BIDS of a species collapsed by median)' if by == 'species' else 'dataset (colour = species, line style = 1st/2nd dataset of that species)'}"
                f", built from every run -- the exact input of {'F1' if by == 'species' else 'F1c'}. "
                f"Right: the same after the final selection round (round {final_round}) -- the "
                f"input of F1b. Black: consensus (median across species)."
                + (" Matrices: upper triangle intra-hemispheric edges, lower triangle "
                   "inter-hemispheric, diagonal homotopic." if matrix_mode else ""),
                f"Profiles {fa._norm_label(norm) or 'unscaled'}{value_col}; regions ordered "
                f"clockwise by the all-runs consensus (as on F1), same order on both sides.",
                "A fingerprint shows the whole pattern at once (Passingham et al. 2002; Mars "
                "et al. 2018): species sharing the consensus shape overlap the black line; a "
                "region where one species bulges out or caves in is where it departs. Compare "
                "left vs right: selection is scored against a cross-species consensus, so the "
                "right side is expected to be tighter around the black line by construction.")

    # ---- F15: distance fingerprints ---------------------------------------
    if len(dist_all) >= 3:
        for which, tag in (("consensus", "F15a_fingerprint_vs_consensus"),
                           ("reference", "F15b_fingerprint_vs_reference")):
            plot_distance_fingerprint_grid(vecs_all, dist_all, species_of_label, which, common,
                                           order, region_mask, opj(fig_dir, f"{tag}.png"),
                                           ref_species, divergence, group_all_source)
        _legend(opj(fig_dir, "F15_fingerprint_vs_reference.txt"),
                "F15a/b -- Each dataset's fingerprint over its reference",
                f"One panel per dataset ({group_all_source}), closest first. F15a: reference "
                f"= LOSO consensus (F6a). F15b: reference = {ref_species} datasets (F6b); a "
                f"{ref_species} dataset is drawn over the other {ref_species} dataset. Shaded "
                f"area = where the dataset exceeds its reference (regions where it is below are "
                f"left unshaded).",
                "Both curves z-scored over the common regions. For z-scored vectors the mean "
                "squared difference equals 2(1 - r), so the shaded gap is what lowers r on "
                "F6, region by region.",
                "Panels are ordered from the most to the least similar (r in each title, as "
                "on F6). The shape of "
                "the shaded area tells WHICH regions make a dataset distant: a uniform thin "
                "band means a small global departure, one or two large lobes mean the "
                "distance is carried by a few regions -- in that case check those regions "
                "on F5 (invariance) and in the atlas (size, partial volume) before "
                "interpreting.")

    # ---- F16: proximity ----------------------------------------------------
    if len(vecs_all) >= 3:
        labels = order_labels(vecs_all, species_of_label, divergence)
        D = distance_matrix(vecs_all, labels)
        test = sibling_test(D, species_of_label, n_perm=n_perm)
        plot_dataset_proximity(D, dist_all, species_of_label, test,
                               opj(fig_dir, "F16_dataset_proximity.png"), group_all_source)
        out["distance_matrix"], out["sibling_test"] = D, test
        _legend(opj(fig_dir, "F16_dataset_proximity.txt"),
                "F16 -- Proximity of each dataset to its sibling and to the other species",
                "(a) 1 - r between the group profiles of every pair of datasets, ordered by "
                "divergence from Human. (b) one row per dataset: each dot is its distance to "
                "another dataset (colour = that dataset's species), star = the other dataset "
                "of the same species, black bar = its distance to the LOSO consensus (F6a), "
                "orange ring = nearest neighbour; the sibling's rank among the 9 others is "
                "written on the right.",
                f"Group profiles: {group_all_source}. Test: mean 1 - r over same-species "
                f"pairs minus mean over cross-species pairs; null = species labels permuted "
                f"across the datasets ({n_perm} shuffles, group sizes kept); one-sided.",
                "If species identity dominates dataset/protocol, stars sit at the far left "
                "and the sibling is the nearest neighbour. Only three species have two "
                "datasets, so the test has few distinct configurations and low power; the "
                "per-dataset ranks are the informative part. A sibling farther away than "
                "other species means the protocol (state, field, preprocessing, coverage) "
                "moves the profile more than the species does -- the two Macaque datasets "
                "also differ in arousal state, so their distance is an upper bound on a "
                "pure site effect.")

    # ---- F17: before / after ---------------------------------------------
    va = group_vectors(ba_all or {}, common, value_col, norm)
    vs = group_vectors(ba_sel or {}, common, value_col, norm)
    shared = sorted(set(va) & set(vs))
    if len(shared) >= 3:
        dist_a = dataset_distances({l: va[l] for l in shared}, species_of_label, divergence,
                                   ref_species)
        dist_s = dataset_distances({l: vs[l] for l in shared}, species_of_label, divergence,
                                   ref_species)
        res, merged = plot_before_after(dist_a, dist_s, ds_state,
                                        opj(fig_dir, "F17_before_after_selection.png"),
                                        ba_source, final_round, ref_species)
        out["before_after"], out["before_after_table"] = res, merged
        out["ba_dist_all"], out["ba_dist_sel"], out["ba_source"] = dist_a, dist_s, ba_source
        _legend(opj(fig_dir, "F17_before_after_selection.txt"),
                "F17 -- Group means before vs after selection",
                f"(a, b) per dataset, r of the group mean of ALL runs (left) and of "
                f"the runs kept by the final round {final_round} (right) to the LOSO "
                f"consensus and to {ref_species}; each reference is rebuilt from the same "
                f"pool (all vs selected) of the other datasets. (c) F6a gradient computed on "
                f"each pool.",
                f"Source: {ba_source}. Spearman rho across datasets tells whether the "
                f"ordering of datasets survives selection.",
                "A rise of r with the consensus is EXPECTED: runs were chosen for "
                "resembling a cross-species consensus. That rise is therefore not a result. "
                "What is informative: (1) whether the ranking of datasets is kept (rho close "
                "to 1); (2) whether the gradient of (c) survives. Selection toward a shared "
                "consensus compresses between-species differences, so a gradient that "
                "persists after selection is conservative with respect to the selection -- "
                "but the kept pools differ in size (top 10 of 25 vs of 170 runs), and the "
                "expected similarity of a top-N subset grows with pool size (order "
                "statistics), so (c) is a sensitivity analysis, never the primary test.")

    # ---- F18: arousal state ----------------------------------------------
    mixed = [lab for lab, s in ds_state.items() if s == "mixed"]
    if mixed and len(vecs_all) >= 3:
        def animal_of(lab, key):
            return resolvers[lab].animal(key)
        df_runs, run_vecs = run_level_similarity(runs_by_label, vecs_all, species_of_label,
                                                 common, value_col, norm, mixed,
                                                 selections_final, animal_of)
        df_runs["state"] = [run_states[l].get(k, "unknown")
                            for l, k in zip(df_runs["label"], df_runs["key"])]
        # the contrast is run on ANIMALS: the BIDS subject unless the state
        # rule says how to read the animal out of it (Mouse)
        df_runs["bids_subject"] = df_runs["subject"]
        df_runs["subject"] = [animal_of(l, k) for l, k in zip(df_runs["label"], df_runs["key"])]
        df_runs["state_class"] = df_runs["state"].map(state_class)
        tests = {}
        for lab in mixed:
            for measure in ("r_consensus", "r_selected"):
                t = state_contrast(df_runs[df_runs["label"] == lab], measure, n_perm=n_perm)
                if t:
                    tests[(lab, measure)] = t
        grp = state_group_similarity(df_runs, run_vecs, vecs_all, species_of_label)
        sel_counts = selection_by_state(df_runs, n_perm=n_perm)
        fps = None if matrix_mode else state_fingerprints(df_runs, run_vecs, vecs_all,
                                                          species_of_label)
        region_fx = state_region_effects(df_runs, run_vecs, common, n_perm=n_perm)
        plot_state_effect(df_runs, tests, grp, sel_counts,
                          opj(fig_dir, "F18_arousal_state.png"), fps=fps, common=common,
                          order=order, region_mask=region_mask)
        out["state_runs"], out["state_tests"], out["state_group"] = df_runs, tests, grp
        out["state_selection"], out["state_fingerprints"] = sel_counts, fps
        out["state_region_effects"] = region_fx
        f6_range = (np.nanmax(dist_all["r_consensus"]) - np.nanmin(dist_all["r_consensus"])
                    if len(dist_all) else np.nan)
        _legend(opj(fig_dir, "F18_arousal_state.txt"),
                "F18 -- Arousal state, contrasted within the datasets that contain both",
                f"Datasets with awake and anaesthetised scans: {mixed}. Column 1: r of each "
                f"run's regional profile with the leave-one-species-out consensus (group "
                f"profiles of the OTHER species = the round-1 selection target); star = r of "
                f"the state's median profile. Column 2: r with the dataset's own selected "
                f"reference = median of the final-round KEPT runs of the same dataset, "
                f"leaving out every run of the run's own animal. Column 3: share of each "
                f"state's runs kept by the final selection. Column 4: fingerprint of each "
                f"state (median run profile), the selected reference (dashed) and the "
                f"consensus (black). Dots = runs; joined markers = one animal.",
                "Columns 1-2, within-subject design (animals scanned in both states): effect "
                "= mean over animals of (mean anaesthetised - mean awake), p from "
                "sign-flipping animals (exact when possible); the run-level p shuffles state "
                "among an animal's runs and ignores sessions -- liberal. Column 3: T = number "
                "of kept runs that are anaesthetised; null = state labels permuted among the "
                "runs of each animal (kept set fixed), two-sided; Fisher's p (runs treated as "
                "independent) is in fingerprint_tables.xlsx only. Region by region "
                "(state_region_effects sheet): sign-flip over animals, Benjamini-Hochberg "
                "across regions.",
                f"This is the only contrast where state varies while species, site, coil and "
                f"pipeline are held fixed. Compare the effect with the spread of F6a across "
                f"datasets ({f6_range:.3f} in r here): an effect of the same order means "
                f"arousal state alone could produce differences as large as those attributed "
                f"to phylogeny, because state is collinear with divergence across datasets. "
                f"Column 2 is partly circular in one direction: the selected reference is "
                f"built from the runs the selection preferred, so a state that is rarely kept "
                f"is compared with a reference made mostly of the other state. Column 1 has "
                f"no such link. Isoflurane alters cortical functional connectivity in a "
                f"dose-dependent way in macaques (Hutchison et al. 2014 Hum Brain Mapp "
                f"35:5754) and differs between awake and anaesthetised marmosets (Hori et al. "
                f"2020 Cereb Cortex 30:5943); the anaesthetic changes mouse connectivity "
                f"patterns (Grandjean et al. 2014 NeuroImage 102:838).")

    # ---- F20: selection vs variability levels -----------------------------
    if not matrix_mode:
        vt, vdraws = variability_levels(runs_by_label, selections_final, resolvers, common,
                                        value_col, norm, n_draw=min(n_perm, 2000))
        if not vt.empty:
            pooled_v = {k: pooled_gap_test(vdraws, k) for k in ("gap_animal", "gap_session")}
            plot_variability_levels(vt, vdraws, species_of_label,
                                    opj(fig_dir, "F20_selection_vs_variability_levels.png"),
                                    pooled=pooled_v["gap_animal"])
            out["variability_levels"], out["variability_pooled"] = vt, pooled_v
            res_lines = []
            for r in vt.itertuples():
                res_lines.append(
                    f"{r.label}: {r.n_kept} kept of {r.n_runs} runs ({r.n_animals_kept} of "
                    f"{r.n_animals} animals). Animal gap (median r same animal - different "
                    f"animals): all runs {r.gap_animal_all:+.3f}, selected "
                    f"{r.gap_animal_selected:+.3f}, random subsets {r.gap_animal_null_median:+.3f}, "
                    f"p = {_fmt_p(getattr(r, 'p_gap_animal', np.nan))}. Session gap (inter-run - "
                    f"inter-session): selected {r.gap_session_selected:+.3f}, random "
                    f"{r.gap_session_null_median:+.3f}, p = "
                    f"{_fmt_p(getattr(r, 'p_gap_session', np.nan))}.")
            for k, t in pooled_v.items():
                if t:
                    res_lines.append(f"POOLED {k}: mean z = {t['mean_z']:+.2f}, p = "
                                     f"{_fmt_p(t['p'])} ({t['n_datasets']} datasets)")
            _legend(opj(fig_dir, "F20_selection_vs_variability_levels.txt"),
                    "F20 -- Does the selection smooth inter-run, inter-session or inter-animal "
                    "variability?",
                    "Datasets with several runs per animal only. For every pair of runs, r "
                    "between their regional profiles, classed as inter-run (same animal, same "
                    "session), inter-session (same animal, other session; for Mouse a "
                    "different session is also a different arousal state) and inter-animal. "
                    "Dashed grey: median per class over every run; coloured: over the runs kept "
                    "by the final round; grey bars: 95% of 2000 random subsets of the same "
                    "number of runs drawn from the dataset.",
                    "Animal gap = median r same animal (inter-run and inter-session pairs) - "
                    "median r different animals; session gap = inter-run - inter-session. p: "
                    "the selected gap against the random subsets (two-sided). Any selection "
                    "toward a common target raises every class; the random subsets of the same "
                    "size are the right comparison for the GAP, not for the level. Pooled: "
                    "mean over datasets of the z of the selected gap within its own null.\n"
                    + "\n".join(res_lines),
                    "A selected gap below the random subsets = the selection makes runs of "
                    "different animals as alike as runs of the same animal: it erases "
                    "individual differences (or keeps few runs per animal). A gap above = it "
                    "keeps the identity of the animal. Few kept runs per animal give few "
                    "same-animal pairs (n in fingerprint_tables.xlsx, variability_levels): "
                    "read the gap with its n.")

    # ---- tables -----------------------------------------------------------
    with pd.ExcelWriter(opj(table_dir, "fingerprint_tables.xlsx")) as xw:
        states_df.to_excel(xw, sheet_name="arousal_state", index=False)
        dist_all.to_excel(xw, sheet_name="distances_all", index=False)
        if "f6_species" in out:
            out["f6_species"].to_excel(xw, sheet_name="F6_species_means", index=False)
            pd.DataFrame([dict(test=k, **{kk: vv for kk, vv in v.items()})
                          for k, v in out["f6_tests"].items()]).to_excel(
                xw, sheet_name="F6_tests", index=False)
        if "distance_matrix" in out:
            out["distance_matrix"].to_excel(xw, sheet_name="dataset_distance_matrix")
            if out["sibling_test"]:
                out["sibling_test"]["per_dataset"].to_excel(xw, sheet_name="sibling_ranks",
                                                            index=False)
                pd.DataFrame([{k: v for k, v in out["sibling_test"].items()
                               if k != "per_dataset"}]).to_excel(
                    xw, sheet_name="sibling_test", index=False)
        if "before_after_table" in out:
            out["before_after_table"].to_excel(xw, sheet_name="before_after", index=False)
        if "state_runs" in out:
            out["state_runs"].drop(columns=["key"]).to_excel(xw, sheet_name="state_runs",
                                                            index=False)
            pd.DataFrame([dict(label=k[0], **{kk: vv for kk, vv in v.items()
                                              if kk != "subject_diffs"})
                          for k, v in out["state_tests"].items()]).to_excel(
                xw, sheet_name="state_tests", index=False)
            out["state_group"].to_excel(xw, sheet_name="state_group_profiles", index=False)
            pd.concat([sc.assign(label=lab, p_within_animal=sc.attrs.get("p"),
                                 fisher_p=sc.attrs.get("fisher_p"))
                       for lab, sc in out["state_selection"].items()]).to_excel(
                xw, sheet_name="selection_by_state", index=False)
            if out.get("state_region_effects") is not None and \
                    not out["state_region_effects"].empty:
                out["state_region_effects"].to_excel(xw, sheet_name="state_region_effects",
                                                     index=False)
        if out.get("variability_levels") is not None:
            out["variability_levels"].to_excel(xw, sheet_name="variability_levels", index=False)
    print(f"  [saved] {opj(table_dir, 'fingerprint_tables.xlsx')}")
    out["dist_all"], out["dataset_state"] = dist_all, ds_state
    out["vecs_all"], out["region_order"] = vecs_all, order
    out["group_all_source"] = group_all_source
    out["common"], out["region_mask"] = common, region_mask
    return out
