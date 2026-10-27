"""
profile_stats.py
================
All PURE-ANALYSIS functions for Figure 5, consolidated into one module.

Consolidated from what used to be spread over three files:
  - consensus_anchored_selection.py   (consensus building + selection)
  - human_anchored_phylogeny.py       (region sets, references, figures)
  - Figure5_SBA_KeeporNot.py          (used only as a constants holder --
                                       those constants now live in
                                       exploration_pipeline.CONFIG)

WHAT IS AND IS NOT HERE
-----------------------
Everything in this file is pure numpy/pandas/matplotlib: profiles in,
numbers and figures out. Nothing here touches AFNI, nibabel, atlas files
or the filesystem layout of a BIDS tree.

That work stays in SBA_Region_Profile_Index (imported as `sri` by the
caller), which remains the one external dependency. Exactly five of its
functions are used:
    parse_itk_snap_label_file   read an ITK-SNAP atlas label file
    compute_all_run_profiles    per-run regional profiles from SBA maps
    rerun_balanced_by_run       run-level 3dttest++ (+ optional surfaces)
    find_group_map_path         locate a group map on disk
    copy_run_qc_images          copy SBA preview jpgs
Those wrap AFNI, nibabel and the atlas library; inlining them here would
mean duplicating that machinery, so they are deliberately left in place.

Call graph is one-directional: exploration_pipeline -> profile_stats -> sri.
No module here calls back into the pipeline.
"""

import os
import math

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as _stats
from Statistics.Group_fMRI.run_selection.consensus import (  # moved verbatim, 2026-10-27
    PALETTE,
    SPECIES_COLORS,
    VALID_NORMALIZATIONS,
    _channel_block,
    _label_rng,
    _norm_label,
    _normalize_values,
    _ordered_common_values,
    _palette_color_map,
    _select_against_species_profiles,
    _value_order,
    _zscore_values,
    article_color,
    balanced_kept_vs_notkept,
    build_consensus_pattern,
    build_reference_for_profiles,
    build_selection_target,
    compare_reference_methods_on_round2,
    compute_region_invariance,
    get_common_regions_from_group_maps,
    get_common_regions_from_profiles,
    opj,
    plot_consensus_pattern,
    rank_runs_against_consensus,
    run_consensus_selection,
    run_iterative_selection,
    run_round2_selection,
    select_invariant_regions,
    species_profile,
    validate_selection_against_qc)


# === Constants ===========================================================


def dataset_style(name):
    """(marker, line style) of a dataset inside its species -- same rank rule as
    paper_style.dataset_marker."""
    from Exemples.Study_EDNiX.analysis import paper_style as _ps
    m = _ps.dataset_marker(name)
    ls = {"o": "-", "s": (0, (4, 1.5)), "^": (0, (1, 1)), "D": "-.", "v": (0, (3, 1, 1, 1))}[m]
    return m, ls


# === Publication style of the paper panels (Figures 5, 6, S5a, S6d) ======
# Drawn at their final printed size: widths of a Nature page (89 mm one
# column, 120 mm one and a half, 183 mm two columns), text 5-7 pt (Nature:
# 5-7 pt at final size), lines 0.5-1 pt, fonts embedded as TrueType in the
# PDF (editable text in Illustrator / Inkscape), PNG at 600 dpi.
MM = 1.0 / 25.4
PAPER_WIDTH = dict(single=89 * MM, onehalf=120 * MM, double=183 * MM)
INTRA_COLOR, INTER_COLOR = "#0072B2", "#D55E00"   # Okabe-Ito blue / vermillion


def paper_font():
    """First available of Arial, Helvetica, Liberation Sans (metric clone of
    Arial), Nimbus Sans, DejaVu Sans -- asking matplotlib for a missing family
    prints a warning at every text."""
    from matplotlib import font_manager
    have = {f.name for f in font_manager.fontManager.ttflist}
    for name in ("Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"):
        if name in have:
            return name
    return "sans-serif"


def paper_rc():
    return {
        "font.family": [paper_font()], "font.size": 7, "axes.labelsize": 7,
        "axes.titlesize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
        "legend.fontsize": 6, "legend.title_fontsize": 6, "legend.frameon": False,
        "legend.handlelength": 1.4, "legend.handletextpad": 0.4, "legend.columnspacing": 1.0,
        "legend.borderaxespad": 0.2, "axes.linewidth": 0.5, "axes.labelpad": 2.0,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2.5,
        "ytick.major.size": 2.5, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
        "xtick.minor.width": 0.4, "ytick.minor.width": 0.4, "lines.linewidth": 0.8,
        "lines.markersize": 3.5, "patch.linewidth": 0.5, "hatch.linewidth": 0.5,
        "axes.spines.top": False, "axes.spines.right": False, "axes.titlepad": 3.0,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.dpi": 600, "figure.dpi": 150, "savefig.pad_inches": 0.02,
        "mathtext.default": "regular"}


def paper_style():
    """with paper_style(): ...  -- the rcParams of every paper panel."""
    return plt.rc_context(paper_rc())


def stars(p):
    """Significance marks of the paper panels: *** p < 0.001, ** p < 0.01,
    * p < 0.05, n.s. otherwise ('' when there is no test). The exact p, the
    test and n are in the .txt of the panel."""
    if p is None or not np.isfinite(p):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."


def bracket(ax, x1, x2, y, p=None, text=None, dy=None, fs=6, lw=0.5, color="black",
            clip=False):
    """Significance bracket between x1 and x2 (data units) at height y, with
    the stars of p (or `text`) centred above. dy = arm length (data units,
    default 2 % of the y range)."""
    lab = stars(p) if text is None else text
    if not lab:
        return
    if dy is None:
        lo, hi = ax.get_ylim()
        dy = 0.02 * (hi - lo)
    ax.plot([x1, x1, x2, x2], [y, y + dy, y + dy, y], color=color, lw=lw,
            clip_on=clip, solid_capstyle="butt", zorder=6)
    ax.text((x1 + x2) / 2, y + dy, lab, ha="center", va="bottom", color=color,
            fontsize=fs if lab != "n.s." else fs - 0.5, clip_on=clip, zorder=6)


def repel_labels(ax, items, fs=5, side="right", pad_pt=2.5, leader=True, obstacles=()):
    """Labels of points without overlap, placed AFTER the final layout (call
    just before saving). items = [(x, y, text, color)] in data units. Each
    label goes to the right of its point (then to the left), shifted
    vertically by steps of one line until its box overlaps neither another
    label nor a point (nor an `obstacles` marker: (x, y, radius in pt)) and
    stays inside the axes; a thin leader joins it to its point when it moved."""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    px = fig.dpi / 72.0
    axbb = ax.get_window_extent(renderer)
    pts = [(ax.transData.transform((x, y)), t, c, x, y) for x, y, t, c in items
           if np.isfinite(x) and np.isfinite(y)]
    pr = 2.0 * px                                  # half-size of a point's box
    extra = []
    for ox, oy, rad_pt in obstacles:               # other markers (x, y, radius in pt)
        if np.isfinite(ox) and np.isfinite(oy):
            X, Y = ax.transData.transform((ox, oy))
            extra.append((X - rad_pt * px, Y - rad_pt * px, X + rad_pt * px, Y + rad_pt * px))
    obstacles = [(p[0] - pr, p[1] - pr, p[0] + pr, p[1] + pr) for p, *_ in pts] + extra
    placed = []

    def free(b, X, Y):
        if b[0] < axbb.x0 or b[2] > axbb.x1 or b[1] < axbb.y0 or b[3] > axbb.y1:
            return False
        for o in placed + obstacles:
            own = o[0] <= X <= o[2] and o[1] <= Y <= o[3]     # the label's own marker / ring
            if not own and b[0] < o[2] and o[0] < b[2] and b[1] < o[3] and o[1] < b[3]:
                return False
        return True

    order = sorted(range(len(pts)), key=lambda i: (pts[i][0][0], pts[i][0][1]))
    for i in order:
        (X, Y), txt, col, x, y = pts[i]
        tmp = ax.text(0, 0, txt, fontsize=fs)
        bb = tmp.get_window_extent(renderer)
        tmp.remove()
        w, h = bb.width, bb.height
        choice = None
        sides = (1, -1) if side == "right" else (-1, 1)
        for k in (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6):
            for sgn in sides:
                x0 = X + pad_pt * px if sgn > 0 else X - pad_pt * px - w
                yc = Y + k * 1.05 * h
                b = (x0, yc - h / 2, x0 + w, yc + h / 2)
                if free(b, X, Y):
                    choice = (b, sgn, k)
                    break
            if choice:
                break
        if choice is None:                         # nowhere free: right of the point
            x0 = X + pad_pt * px
            choice = ((x0, Y - h / 2, x0 + w, Y + h / 2), 1, 0)
        b, sgn, k = choice
        placed.append(b)
        inv = ax.transData.inverted()
        xt, yt = inv.transform((b[0] if sgn > 0 else b[2], (b[1] + b[3]) / 2))
        ax.annotate(txt, (x, y), xytext=(xt, yt), textcoords="data", fontsize=fs, color=col,
                    va="center", ha="left" if sgn > 0 else "right",
                    arrowprops=(dict(arrowstyle="-", lw=0.3, color="#AAAAAA", shrinkA=0,
                                     shrinkB=1.2) if (leader and k != 0) else None))


def save_paper(fig, path, svg=False):
    """Paper panel: .png (600 dpi) and .pdf (TrueType fonts), optional .svg;
    `path` with or without the .png extension."""
    base = path[:-4] if path.endswith(".png") else path
    exts = ("png", "pdf") + (("svg",) if svg else ())
    with plt.rc_context({"pdf.fonttype": 42, "svg.fonttype": "none",
                         "image.composite_image": False}):
        for ext in exts:
            fig.savefig(f"{base}.{ext}", dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print(f"  [saved] {base}." + "/.".join(exts))


# Divergence time (Mya) between human and each species = age of their most
# recent common ancestor. PROVENANCE: entered by hand as approximate TimeTree
# median estimates (Kumar et al. 2022, TimeTree 5, timetree.org); they were
# NOT downloaded or queried by this code. Before publication, look each pair
# up on timetree.org, replace the values, and report the access date and the
# TimeTree version in the methods. phylo.SPECIES_CLADES must stay
# consistent with this dict (checked automatically at run time).
DIVERGENCE_TIME_MYA = {
    "Human": 0.0,
    "Macaque": 29.0,
    "Marmoset": 43.0,
    "Mouselemur": 74.0,
    "Rat": 87.0,
    "Mouse": 87.0,
    "Dog": 94.0,
}


# === Low-level helpers ===================================================


def _sort_regions_by_value(reference_profile, common_base_names, value_col="median_r"):
    """
    Order regions by their value in `reference_profile`, ascending
    (lowest median_r first) -- the same convention as the per-BIDS
    profile figures, where a monotonically rising reference curve makes
    each run's deviation from it immediately readable. Regions with no
    value in the reference keep their original relative order at the end.
    """
    vals = _ordered_common_values(reference_profile, common_base_names, value_col)
    return [common_base_names[i] for i in _value_order(vals, common_base_names)]


def _drop_all_nan_regions(profiles_list, common_base_names, value_col="median_r"):
    """
    Keep only regions that have a real value in AT LEAST ONE of the given
    profiles. Confirmed real symptom without this: the figures drew all
    138 candidate region ticks while only a handful carried data, making
    the plots unreadable and implying coverage that does not exist.
    """
    keep = []
    for i, name in enumerate(common_base_names):
        vals = [_ordered_common_values(p, [name], value_col)[0] for p in profiles_list
                if p is not None]
        if any(np.isfinite(v) for v in vals):
            keep.append(name)
    return keep


def _keep_regions_present_in_all(profiles_list, common_base_names, value_col="median_r"):
    """
    Keep only regions that have a real value in EVERY given profile --
    the strict counterpart to _drop_all_nan_regions (which keeps a region
    if ANY profile has it). Needed for the cross-species overlay figure:
    a region missing from even one species' group map leaves a visible
    GAP in that species' line, so the shared axis must be the strict
    intersection, not the union.
    """
    keep = []
    for name in common_base_names:
        vals = [_ordered_common_values(p, [name], value_col)[0]
                for p in profiles_list if p is not None]
        if vals and all(np.isfinite(v) for v in vals):
            keep.append(name)
    return keep


def restrict_profile(profile, common_base_names):
    """Keep only rows whose base_name is in common_base_names."""
    if profile is None or "base_name" not in profile.columns:
        return profile
    mask = profile["base_name"].str.lower().str.strip().isin(common_base_names)
    return profile[mask].copy()


# === Region sets and references ==========================================


# === Consensus and run selection =========================================

def diagnose_consensus_sign(species_profiles, common_base_names, value_col="median_r",
                            species_of_label=None, verbose=True):
    """
    Ground truth for "the consensus is never negative -- is that a bug?":
    for every region, the per-species (or per-BIDS, if species_of_label is
    None) contributor values that build_consensus_pattern's median is
    actually computed from, plus how many are positive/negative/NaN.

    A region where the median stays (barely) positive despite several
    negative contributors is expected, not a bug -- median only needs a
    majority on one side, and individual species/BIDS dipping negative
    from noise doesn't change that as long as most don't. This function
    lets you SEE which situation you're in per region instead of
    inferring it from the aggregate range alone.

    Returns a DataFrame: region, median (raw, no z-score -- matches what
    F4's Cross-species consensus line actually plots), n_pos, n_neg,
    n_nan, n_total, one column per species/BIDS with its raw value.
    """
    contributors = {}
    for label, prof in species_profiles.items():
        key = (species_of_label or {}).get(label, label)
        vals = _ordered_common_values(prof, common_base_names, value_col)
        contributors.setdefault(key, []).append(vals)
    with np.errstate(invalid="ignore"):
        by_key = {k: np.nanmedian(np.vstack(v), axis=0) for k, v in contributors.items()}

    rows = []
    for j, region in enumerate(common_base_names):
        vals = {k: by_key[k][j] for k in by_key}
        finite = [v for v in vals.values() if np.isfinite(v)]
        row = dict(region=region,
                   median=float(np.median(finite)) if finite else np.nan,
                   n_pos=sum(v > 0 for v in finite), n_neg=sum(v < 0 for v in finite),
                   n_nan=sum(not np.isfinite(v) for v in vals.values()),
                   n_total=len(vals))
        row.update({f"val_{k}": vals[k] for k in sorted(vals)})
        rows.append(row)
    df = pd.DataFrame(rows).sort_values("median")
    if verbose:
        close_calls = df[(df["median"] > 0) & (df["n_neg"] >= df["n_pos"])]
        print(f"  [diagnose_consensus_sign] {len(df)} region(s), "
              f"{(df['median'] < 0).sum()} with negative median, "
              f"{(df['median'] > 0).sum()} positive")
        if len(close_calls):
            print(f"    {len(close_calls)} region(s) positive median despite "
                 f">= as many negative as positive contributors (median barely "
                 f"tipping positive, not a clean majority):")
            print(close_calls[["region", "median", "n_pos", "n_neg", "n_nan"]]
                  .to_string(index=False))
    return df


def compute_loo_species_distances(profiles_by_label, common_base_names, value_col="median_r",
                                  zscore=True, verbose=True):
    """
    PHYLOGENY-SIDE ONLY (feeds F6 and PGLS's distance_to_consensus). Has
    NO leave_one_out parameter -- LOO is unconditional here, and this is
    NOT controlled by cfg["LEAVE_ONE_OUT"] (that flag only affects
    selection, see build_selection_target above). The two are
    deliberately independent: toggling selection's LOO never changes
    what F6/PGLS see, and vice versa. Do not disable this one to "match"
    a selection-side experiment -- it is what makes the phylogenetic
    signal claim non-circular in the first place; there is no
    engineering-tradeoff argument for turning it off the way there is
    for selection.

    For every label (BIDS) in profiles_by_label ({label: (species,
    profile_df)}), correlate its profile against a leave-one-SPECIES-out
    consensus (build_consensus_pattern with species_of_label, which
    excludes every BIDS of that label's own species, not just the label
    itself -- see build_consensus_pattern's docstring for why a
    same-species BIDS left in would still leak). This is what makes a
    species' own distance non-circular regardless of how many BIDS it
    has, so nothing needs to be excluded from the fit afterward.

    Returns a DataFrame, one row per BIDS: label, species,
    distance_to_human (1 - r), pearson_r, n_regions, n_species_in_ref.
    ("distance_to_human" kept as the column name for compatibility with
    plot_distance_vs_divergence/test_phylogenetic_signal, even though
    the reference is a leave-one-species-out cross-species consensus,
    not literally human.)
    """
    dist_rows = []
    species_of_label = {lab: sp for lab, (sp, p) in profiles_by_label.items()}
    profiles_only = {lab: p for lab, (sp, p) in profiles_by_label.items()}
    for lab, (sp, p) in sorted(profiles_by_label.items()):
        v = _ordered_common_values(p, common_base_names, value_col)
        ref_loo, used = build_consensus_pattern(
            profiles_only, common_base_names, value_col=value_col,
            exclude_label=lab, zscore=zscore, threshold_abs=None,
            species_of_label=species_of_label, verbose=False)
        both = np.isfinite(v) & np.isfinite(ref_loo)
        if both.sum() < 3:
            if verbose:
                print(f"  [SKIP] {lab}: only {both.sum()} finite common region(s) "
                      f"vs its leave-one-species-out reference, need >= 3")
            continue
        rr = float(np.corrcoef(v[both], ref_loo[both])[0, 1])
        # `used` lists the DATASETS behind the reference (the column used to
        # report that count under the name n_species_in_ref)
        dist_rows.append(dict(label=lab, species=sp, distance_to_human=1 - rr,
                              pearson_r=rr, n_regions=int(both.sum()),
                              n_species_in_ref=len({species_of_label[u] for u in used}),
                              n_datasets_in_ref=len(used)))
    return pd.DataFrame(dist_rows)


def compute_run_distance_to_reference(runs_by_label, ref_profile_by_label, common_base_names,
                                      value_col="median_r", kept_keys_by_label=None, verbose=True):
    """
    For every label with a reference profile (ref_profile_by_label --
    typically each BIDS' "kept" group-map profile), correlate EVERY run
    in runs_by_label[label] against it individually -- kept, not-kept,
    and anything outside the balanced 3dttest++ alike, since the point
    is to see the full spread each BIDS actually produced, not just the
    runs that made the cut.

    kept_keys_by_label: optional {label: set-of-run-keys} (e.g. the
    step-2 `selections` dict) -- when given, an `is_kept` bool column
    marks which rows were in the selected/kept subset, so a plot can
    show where they sit within the full spread.

    Returns a long DataFrame: label, species, subject, session, run,
    seed, pearson_r, distance (1-r), n_regions, is_kept.
    """
    rows = []
    for label, (species, ref_profile) in ref_profile_by_label.items():
        profs = runs_by_label.get(label)
        if not profs:
            continue
        kept_keys = {tuple(str(x) for x in k) for k in
                    (kept_keys_by_label or {}).get(label, [])}
        ref_v = _ordered_common_values(ref_profile, common_base_names, value_col)
        for run_key, run_profile in profs.items():
            v = _ordered_common_values(run_profile, common_base_names, value_col)
            both = np.isfinite(v) & np.isfinite(ref_v)
            if both.sum() < 3:
                continue
            rr = float(np.corrcoef(v[both], ref_v[both])[0, 1])
            rows.append(dict(label=label, species=species,
                             subject=run_key[0], session=run_key[1],
                             run=run_key[2], seed=run_key[3],
                             pearson_r=rr, distance=1 - rr,
                             n_regions=int(both.sum()),
                             is_kept=tuple(str(x) for x in run_key) in kept_keys))
    result = pd.DataFrame(rows)
    if verbose:
        print(f"  compute_run_distance_to_reference: {len(result)} run(s) scored "
              f"across {result['label'].nunique() if not result.empty else 0} BIDS")
    return result


def _qc_norm_id(x):
    """Session/run id normalisation for the QC join: '01' == '1', and
    'NA' / '' / NaN (pandas turns the string 'NA' into NaN when reading an
    Excel file, even with dtype=str) all become 'NA'."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "NA"
    t = str(x).strip()
    if t in ("", "nan", "NaN", "None", "NA", "<NA>"):
        return "NA"
    return t.lstrip("0") or "0"


def correlate_qc_with_distance(dist_df, df_qc_by_label, qc_metrics=None,
                               id_cols=("subject", "session"), verbose=True):
    """
    Within each BIDS, Pearson r between each QC metric and each unit's
    distance to its dataset reference (dist_df: one row per run).

    Granularity is decided PER METRIC, from the QC table itself:
      - metric constant within every (subject, session) -> session-level
        metric (e.g. tSNR/FD computed once per session and repeated on
        each run's row): QC and distance are both reduced to ONE row per
        session (distance = mean over that session's runs);
      - metric varying within sessions AND run ids matching >= 90% of the
        distance rows -> run-level join (one row per run);
      - otherwise -> session-level (QC averaged per session).

    Confirmed bug this replaces: distance was averaged per session, then
    joined on (subject, session) to a QC table with one row PER RUN, so
    every session was duplicated by its number of runs (real log:
    "matched 40/4", "169/57", "195/46", "52/31"). n, hence every p-value
    and every meta-analytic weight (n - 3), was inflated by pseudo-
    replication. A merge can no longer produce more rows than its left
    side (asserted).

    Returns a long DataFrame: label, metric, level ("run"/"session"), n,
    pearson_r, pearson_p.
    """
    if dist_df is None or dist_df.empty or not df_qc_by_label:
        if verbose:
            print("  [SKIP] correlate_qc_with_distance: no distance data or no QC table")
        return pd.DataFrame()

    no_qc = sorted(set(dist_df["label"]) - set(df_qc_by_label))
    if no_qc and verbose:
        print(f"  [WARN] no QC table for {no_qc} -- excluded from the QC analysis")

    rows = []
    for label, df_qc in df_qc_by_label.items():
        d = dist_df[dist_df["label"] == label].copy()
        if d.empty or df_qc is None or df_qc.empty or "subject" not in df_qc.columns:
            continue
        q = df_qc.copy()
        for df_ in (d, q):
            df_["_sub"] = df_["subject"].astype(str).str.strip()
            df_["_ses"] = (df_["session"].map(_qc_norm_id) if "session" in df_.columns else "NA")
            df_["_run"] = (df_["run"].map(_qc_norm_id) if "run" in df_.columns else "NA")
        d_ses = (d.groupby(["_sub", "_ses"], as_index=False)
                  .agg(distance=("distance", "mean"), n_runs=("run", "nunique")))

        run_join = False
        if "run" in q.columns and (q["_run"] != "NA").any():
            kd = set(zip(d["_sub"], d["_ses"], d["_run"]))
            kq = set(zip(q["_sub"], q["_ses"], q["_run"]))
            run_join = len(kd & kq) >= 0.9 * len(kd)

        metrics = qc_metrics or [c for c in q.columns if
                                 c.startswith(("func_", "anat_", "net_", "homotopic_",
                                               "cross_", "sp_", "contam_"))
                                 and pd.api.types.is_numeric_dtype(q[c])]
        n_by_level = {}
        for metric in metrics:
            if metric not in q.columns:
                continue
            qv = q[["_sub", "_ses", "_run", metric]].copy()
            qv[metric] = pd.to_numeric(qv[metric], errors="coerce")
            per_ses_unique = qv.dropna(subset=[metric]).groupby(["_sub", "_ses"])[metric].nunique()
            session_level = per_ses_unique.empty or bool((per_ses_unique <= 1).all())
            if run_join and not session_level:
                qm = qv.groupby(["_sub", "_ses", "_run"], as_index=False)[metric].mean()
                m = d.merge(qm, on=["_sub", "_ses", "_run"], how="inner")
                level, n_left = "run", len(d)
            else:
                qm = qv.groupby(["_sub", "_ses"], as_index=False)[metric].mean()
                m = d_ses.merge(qm, on=["_sub", "_ses"], how="inner")
                level, n_left = "session", len(d_ses)
            assert len(m) <= n_left, f"{label}/{metric}: QC join duplicated rows"
            n_by_level.setdefault(level, (len(m), n_left))
            x = pd.to_numeric(m["distance"], errors="coerce")
            y = m[metric]
            ok = x.notna() & y.notna()
            if ok.sum() < 4 or y[ok].nunique() < 2:
                continue
            r, p = _stats.pearsonr(x[ok], y[ok])
            rows.append(dict(label=label, metric=metric, level=level, n=int(ok.sum()),
                             pearson_r=float(r), pearson_p=float(p)))
        if verbose:
            desc = ", ".join(f"{lvl}-level: {a}/{b} matched" for lvl, (a, b) in n_by_level.items())
            print(f"  [{label}] QC join -- {desc or 'no metric'}")

    result = pd.DataFrame(rows)
    if verbose and not result.empty:
        print(result.sort_values("pearson_p").head(30).to_string(index=False))
    return result


def meta_analyze_qc_correlations(per_bids_corr_df, min_bids=3, alpha=0.05, verbose=True):
    """
    Combines correlate_qc_with_distance's per-BIDS rows into ONE number
    per metric, WITHOUT ever pooling raw rows across BIDS -- a standard
    fixed-effect correlation meta-analysis: Fisher-z each BIDS's r
    (z_i = arctanh(r_i)), weight by its precision (n_i - 3), and test
    whether the weighted-average z is consistently non-zero across BIDS.
    This is what "significant in inter-species" should mean here: not a
    big pooled r (which the between-BIDS confound can produce on its
    own, see correlate_qc_with_distance's docstring), but a QC-distance
    relationship that shows up with the same sign across independent
    datasets.

    min_bids: a metric needs at least this many BIDS with a computed r
    to even attempt the meta-test (3 is the same floor used elsewhere in
    this pipeline for phylogenetic signal -- below that, one dataset can
    flip the result and "significant" is not a meaningful claim).

    Returns a DataFrame: metric, n_bids, meta_r, meta_p, n_same_sign
    (out of n_bids), sorted by meta_p. Only rows with n_bids >= min_bids
    are returned at all.
    """
    if per_bids_corr_df is None or per_bids_corr_df.empty:
        if verbose:
            print("  [SKIP] meta_analyze_qc_correlations: nothing to combine")
        return pd.DataFrame()

    rows = []
    for metric, g in per_bids_corr_df.groupby("metric"):
        g = g.dropna(subset=["pearson_r", "n"])
        g = g[g["n"] > 3]
        if len(g) < min_bids:
            continue
        z = np.arctanh(g["pearson_r"].clip(-0.999, 0.999).values)
        w = (g["n"].values - 3).astype(float)
        z_meta = float(np.sum(w * z) / np.sum(w))
        se_meta = float(1.0 / np.sqrt(np.sum(w)))
        z_stat = z_meta / se_meta
        p_meta = float(2 * _stats.norm.sf(abs(z_stat)))
        r_meta = float(np.tanh(z_meta))
        n_same_sign = int((np.sign(g["pearson_r"].values) == np.sign(r_meta)).sum())
        rows.append(dict(metric=metric, n_bids=len(g), meta_r=r_meta, meta_p=p_meta,
                         n_same_sign=n_same_sign))

    result = pd.DataFrame(rows)
    if result.empty:
        if verbose:
            print(f"  [meta] no metric had >= {min_bids} BIDS with a computed correlation")
        return result
    result = result.sort_values("meta_p")
    # Benjamini-Hochberg across metrics: ~40 QC metrics are tested, so at
    # alpha=0.05 about 2 would pass by chance with no correction.
    p = result["meta_p"].values
    order = np.argsort(p)
    q = np.empty(len(p))
    running = 1.0
    for rank, idx in reversed(list(enumerate(order, start=1))):
        running = min(running, p[idx] * len(p) / rank)
        q[idx] = running
    result["meta_q_fdr"] = q
    n_sig = int((result["meta_q_fdr"] < alpha).sum())
    if verbose:
        print(f"  [meta] {n_sig}/{len(result)} metric(s) significant at q_FDR<{alpha} "
              f"(BH across metrics) across >= {min_bids} BIDS "
              f"({int((result['meta_p'] < alpha).sum())} at uncorrected p<{alpha})")
        print(result.head(20).to_string(index=False))
    return result


def plot_region_invariance(invariance, common_base_names, mask=None, title=None,
                           output_path=None, dpi=300):
    """Bar chart of per-region invariance, regions used for selection
    highlighted. Makes explicit which regions the selection actually
    rests on."""
    inv = np.asarray(invariance, dtype=float)
    order = np.argsort(-np.nan_to_num(inv, nan=-1))
    names = [common_base_names[i] for i in order]
    vals = inv[order]
    used = (np.asarray(mask, bool)[order] if mask is not None
            else np.ones(len(vals), bool))

    fig, ax = plt.subplots(figsize=(max(8, len(names) * 0.3), 4.5))
    ax.bar(range(len(vals)), vals,
          color=["#009E73" if u else "#CCCCCC" for u in used], edgecolor="white")
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=90, fontsize=6)
    ax.set_ylabel("invariance across species", fontsize=10)
    ax.set_title(title or "Region invariance (green = used for selection)",
                fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_iterative_selection_convergence(overlap_df, output_path=None, dpi=300):
    """One line per BIDS: fraction of the n_select kept runs unchanged
    from the previous round, across rounds. Flat at 1.0 = converged;
    still moving at the last round = n_rounds was too small to trust,
    or this BIDS is oscillating (also worth knowing)."""
    if overlap_df is None or overlap_df.empty:
        print("  [plot] nothing to plot")
        return None
    df = overlap_df.copy()
    df["frac_overlap"] = df["n_overlap_with_previous"] / df["n_select"]
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = _palette_color_map(df["label"].unique())
    for label, g in df.groupby("label"):
        g = g.sort_values("round")
        ax.plot(g["round"], g["frac_overlap"], marker="o", color=colors[label], label=label)
    ax.set_xlabel("round", fontsize=10)
    ax.set_ylabel("fraction of kept runs unchanged\nfrom previous round", fontsize=9)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xticks(sorted(df["round"].unique()))
    ax.set_title("Iterative selection: convergence by round", fontsize=11, fontweight="bold")
    ax.legend(fontsize=7, ncol=2)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_reference_method_overlap(overlap_df, output_path=None, dpi=300):
    """
    Bar chart: for each BIDS, how many of the round-2 selected runs
    agree between REFERENCE_METHOD="kcore" and "uniform_median" (out of
    n_select). A full-height bar means the method choice doesn't matter
    for that BIDS; a short bar means it changes which runs get kept.
    """
    if overlap_df is None or overlap_df.empty:
        print("  [plot] nothing to plot")
        return None
    df = overlap_df.sort_values("n_overlap")
    fig, ax = plt.subplots(figsize=(max(6, 0.7 * len(df)), 5))
    colors = ["#D55E00" if r.n_overlap < r.n_select else "#0072B2" for r in df.itertuples()]
    ax.bar(df["label"], df["n_overlap"], color=colors, edgecolor="black")
    for i, row in enumerate(df.itertuples()):
        ax.text(i, row.n_overlap + 0.1, f"{row.n_overlap}/{row.n_select}",
               ha="center", va="bottom", fontsize=8)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["label"], rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("Runs agreeing between kcore and uniform_median\n(round-2 selection)",
                 fontsize=9)
    ax.set_title("Effect of REFERENCE_METHOD on round-2 run selection",
                fontsize=11, fontweight="bold")
    ax.set_ylim(0, df["n_select"].max() + 1)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


# === Statistics ==========================================================


def test_phylogenetic_signal(distance_df, divergence_times=None, exclude_human=True,
                             verbose=True):
    """
    Correlate each species' distance-to-human against its divergence time
    from human.

    exclude_human: True (default) drops the human row -- it is at
    distance ~0 from itself BY CONSTRUCTION, and including a point that
    is fixed at (0, 0) inflates the correlation without evidence.

    Returns a dict with pearson_r/p and spearman_rho/p, plus the
    DataFrame actually used.
    """
    divergence_times = divergence_times or DIVERGENCE_TIME_MYA
    df = distance_df.copy()
    df["divergence_mya"] = df["species"].map(divergence_times)

    missing = df[df["divergence_mya"].isna()]["species"].unique()
    if len(missing):
        print(f"  [WARN] no divergence time for {sorted(missing)} -- excluded "
             f"(add them to DIVERGENCE_TIME_MYA)")
    df = df[df["divergence_mya"].notna()]

    if exclude_human:
        df = df[df["species"] != "Human"]

    df = df[df["distance_to_human"].notna()]
    if len(df) < 3:
        print(f"  [WARN] only {len(df)} species with both a distance and a divergence "
             f"time -- too few for a correlation")
        return dict(n=len(df), pearson_r=np.nan, pearson_p=np.nan,
                   spearman_rho=np.nan, spearman_p=np.nan, data=df)

    pr, pp = _stats.pearsonr(df["divergence_mya"], df["distance_to_human"])
    sr, sp = _stats.spearmanr(df["divergence_mya"], df["distance_to_human"])

    if verbose:
        print(f"\n  Phylogenetic signal (n={len(df)} species, human "
             f"{'excluded' if exclude_human else 'included'}):")
        print(f"    Pearson  r={pr:.4f}, p={pp:.4f}")
        print(f"    Spearman rho={sr:.4f}, p={sp:.4f}")
        if len(df) < 6:
            print(f"    [NOTE] with n={len(df)}, this is a trend at best -- "
                 f"a single species shifts it substantially.")

    return dict(n=len(df), pearson_r=pr, pearson_p=pp,
               spearman_rho=sr, spearman_p=sp, data=df)


def permutation_test_phylogeny(distance_df, divergence_times, n_perm=10000,
                               random_seed=42, verbose=True):
    """
    Is the distance-vs-divergence correlation stronger than chance?

    With only ~6-8 species, a Pearson p-value rests on normality
    assumptions that cannot be checked at that n, and one species moves
    it a lot. Shuffling the divergence times across species and
    recomputing gives an exact, assumption-free null for exactly this
    sample size. p is computed with the (1 + count) / (n_perm + 1)
    correction (Phipson & Smyth 2010): the observed arrangement is
    itself one member of the null, so a p of exactly 0 is never
    reportable.

    Confirmed real bug this fixes: this function used to hard-drop
    divergence_mya == 0 (Human) regardless of what exclude_human was
    set to in test_phylogenetic_signal(), so F6 (parametric) and F7
    (this function) silently ran on different samples (n=7 vs n=6) and
    could report contradictory r/p with no indication why. The
    filter's original justification -- "distance 0 from itself carries
    no information" -- was written when Human's distance actually was
    circularly ~0 by construction; since the leave-one-species-out fix,
    Human's distance is a genuine, non-trivial measurement (often one
    of the LARGEST in the dataset), so dropping it is no longer
    self-evidently correct. This function now uses distance_df exactly
    as given -- pass it test_phylogenetic_signal()'s own tr["data"]
    (not a fresh, unfiltered distance_df) so the two figures can never
    diverge on which species were included again.

    Returns dict(observed_r, p_perm, null_r, n_species) or None if fewer
    than 4 species have both a distance and a divergence time.
    """
    from scipy import stats as _st
    rng = np.random.default_rng(random_seed)

    df = distance_df.copy()
    if "divergence_mya" not in df.columns:
        df["divergence_mya"] = df["species"].map(divergence_times)
    df = df[df["divergence_mya"].notna() & df["distance_to_human"].notna()]
    if len(df) < 4:
        if verbose:
            print(f"  [SKIP] permutation test: only {len(df)} species with both values")
        return None

    x = df["divergence_mya"].values
    y = df["distance_to_human"].values
    observed_r = float(_st.pearsonr(x, y)[0])

    null_r = np.empty(n_perm)
    for i in range(n_perm):
        null_r[i] = _st.pearsonr(rng.permutation(x), y)[0]

    p = (1 + int(np.sum(np.abs(null_r) >= abs(observed_r)))) / (n_perm + 1)
    if verbose:
        print(f"  Phylogeny permutation test (n={len(df)} species, {n_perm} shuffles):")
        print(f"    observed r={observed_r:+.4f}, permutation p={p:.4f} (two-sided)")
        if len(df) < 6:
            print(f"    [NOTE] with n={len(df)}, a single species shifts this substantially.")
    return dict(observed_r=observed_r, p_perm=float(p), null_r=null_r, n_species=len(df))


# === Figures =============================================================


def plot_consensus_by_bids(species_profiles, common_base_names, species_of_label,
                           value_col="median_r", zscore=True, title=None,
                           output_path=None, dpi=300):
    """
    Companion to plot_consensus_by_species: same black consensus line
    (species-collapsed, so a 2-BIDS species still counts once -- NOT
    recomputed here, taken as-is from the species-level figure so both
    plots are directly comparable), but one colored line per BIDS
    instead of per species (10 lines here vs 7 there) -- lets you see,
    e.g., whether Macaque_BenHamed and Macaque_ZhuGarin actually agree
    with each other before they get collapsed into one "Macaque" line.
    Colors follow the same per-BIDS PALETTE as F6/F8/F9 (_palette_color_map),
    not SPECIES_COLORS.
    """
    by_species, per_bids = {}, {}
    for label, prof in species_profiles.items():
        vals = _ordered_common_values(prof, common_base_names, value_col)
        if not np.isfinite(vals).any():
            continue
        vec = _normalize_values(vals, zscore)
        per_bids[label] = vec
        by_species.setdefault(species_of_label.get(label, label), []).append(vec)

    if not per_bids:
        print("  [plot] nothing to plot")
        return None

    with np.errstate(invalid="ignore"):
        per_species = {sp: np.nanmedian(np.vstack(v), axis=0) for sp, v in by_species.items()}
        consensus = np.nanmedian(np.vstack(list(per_species.values())), axis=0)

    order = _value_order(consensus, common_base_names)
    names_sorted = [common_base_names[i] for i in order]
    x = np.arange(len(names_sorted))
    bids_colors = _palette_color_map(per_bids.keys())

    fig, ax = plt.subplots(figsize=(max(10, len(names_sorted) * 0.3), 6))
    for label, vec in sorted(per_bids.items()):
        ax.plot(x, vec[order], lw=1.6, alpha=0.9, color=bids_colors[label], label=label)
    ax.plot(x, consensus[order], color="black", lw=3.0, zorder=5,
           label=f"consensus ({len(per_species)} species)")

    ax.axhline(0, color="#999999", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(names_sorted, rotation=90, fontsize=6)
    ax.set_ylabel(f"{_norm_label(zscore)}{value_col}", fontsize=10)
    ax.set_title(title or "Consensus pattern, one line per BIDS", fontsize=12, fontweight="bold")
    ax.legend(fontsize=7, ncol=2)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_consensus_by_species(species_profiles, common_base_names, species_of_label,
                              value_col="median_r", zscore=True, title=None,
                              output_path=None, dpi=300):
    """
    FIGURE: the consensus with ONE LINE PER SPECIES (BIDS of the same
    species collapsed by median first), plus the overall consensus in
    black.

    Why per species and not per BIDS: with 2 Macaque, 2 Marmoset and 2
    Human datasets against 1 Dog and 1 Mouse, a per-BIDS plot (and a
    per-BIDS consensus) silently weights some species two or three times
    more than others. Collapsing within species first makes every
    species count once, which is what a claim about conserved
    organisation requires.
    """
    by_species = {}
    for label, prof in species_profiles.items():
        vals = _ordered_common_values(prof, common_base_names, value_col)
        if not np.isfinite(vals).any():
            continue
        vec = _normalize_values(vals, zscore)
        by_species.setdefault(species_of_label.get(label, label), []).append(vec)

    if not by_species:
        print("  [plot] nothing to plot")
        return None

    with np.errstate(invalid="ignore"):
        per_species = {sp: np.nanmedian(np.vstack(v), axis=0) for sp, v in by_species.items()}
        consensus = np.nanmedian(np.vstack(list(per_species.values())), axis=0)

    order = _value_order(consensus, common_base_names)
    names_sorted = [common_base_names[i] for i in order]
    x = np.arange(len(names_sorted))

    fig, ax = plt.subplots(figsize=(max(10, len(names_sorted) * 0.3), 6))
    for sp, vec in sorted(per_species.items()):
        n_bids = len(by_species[sp])
        ax.plot(x, vec[order], lw=1.8, alpha=0.9, color=SPECIES_COLORS.get(sp),
               label=f"{sp} ({n_bids} BIDS)")
    ax.plot(x, consensus[order], color="black", lw=3.0, zorder=5,
           label=f"consensus ({len(per_species)} species)")

    ax.axhline(0, color="#999999", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(names_sorted, rotation=90, fontsize=6)
    ax.set_ylabel(f"{_norm_label(zscore)}{value_col}", fontsize=10)
    ax.set_title(title or "Consensus pattern, one line per species",
                fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, ncol=2)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_common_regions_with_human_ref(human_ref, common_base_names,
                                       value_col="median_r", title=None,
                                       output_path=None, dpi=300,
                                       sort_by_reference=True):
    """
    DEPRECATED -- not called by the current pipeline (selection is no longer
    anchored to a human reference); kept only for old scripts.

    FIGURE 1 -- the human reference profile across the regions common to
    every species: what the anchor actually looks like, region by region,
    on the exact axis every other species will be compared on.
    """
    if sort_by_reference:
        common_base_names = _sort_regions_by_value(human_ref, common_base_names, value_col)
    vals = _ordered_common_values(human_ref, common_base_names, value_col)
    x = np.arange(len(common_base_names))

    fig, ax = plt.subplots(figsize=(max(10, len(common_base_names) * 0.22), 5))
    ax.plot(x, vals, color="#D95F02", lw=1.8, marker="o", ms=3, label="Human reference")
    ax.axhline(0, color="#999999", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(common_base_names, rotation=90, fontsize=4)
    ax.set_ylabel(value_col, fontsize=10)
    ax.set_xlabel(f"{len(common_base_names)} common regions"
                 f"{' (ordered by reference value)' if sort_by_reference else ''}", fontsize=9)
    ax.set_title(title or "Human reference profile over cross-species common regions",
                fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_all_bids_refs_vs_human(bids_ref_profiles, human_ref, common_base_names,
                                value_col="median_r", title=None,
                                output_path=None, dpi=300,
                                sort_by_reference=True):
    """
    DEPRECATED -- not called by the current pipeline (selection is no longer
    anchored to a human reference); kept only for old scripts.

    FIGURE 2 -- every BIDS's own selected-group reference overlaid on the
    human reference, all on the common-region axis.

    bids_ref_profiles: {label: profile} -- typically each BIDS's group
    map profile from group_map_profile(), i.e. the reference built from
    its selected runs.

    Human is drawn thick/black so the others read as deviations from it.
    """
    if sort_by_reference:
        # ordered by the HUMAN reference, so every species is plotted on
        # one shared ordering and deviations from human read directly
        common_base_names = _sort_regions_by_value(human_ref, common_base_names, value_col)
    x = np.arange(len(common_base_names))
    fig, ax = plt.subplots(figsize=(max(11, len(common_base_names) * 0.24), 6))

    ax.plot(x, _ordered_common_values(human_ref, common_base_names, value_col),
           color="black", lw=2.6, label="Human reference", zorder=5)

    for label, prof in sorted(bids_ref_profiles.items()):
        species = next((s for s in SPECIES_COLORS if label.startswith(s)), None)
        ax.plot(x, _ordered_common_values(prof, common_base_names, value_col),
               lw=1.2, alpha=0.85, label=label,
               color=SPECIES_COLORS.get(species, None))

    ax.axhline(0, color="#999999", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(common_base_names, rotation=90, fontsize=4)
    ax.set_ylabel(value_col, fontsize=10)
    ax.set_title(title or "Each BIDS reference vs the human reference (common regions)",
                fontsize=11, fontweight="bold")
    ax.legend(fontsize=7, ncol=2)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_bids_ref_with_individuals(profiles, species_ref, consensus_ref, common_base_names,
                                   selected_keys=None, kept_ref=None, notkept_ref=None,
                                   grand_mean_ref=None, kept_ref_round2=None,
                                   actual_selection_target=None,
                                   final_round_label="Kept final round",
                                   target_label="Round-1 selection target (right axis)",
                                   value_col="median_r",
                                   title=None, output_path=None, dpi=300,
                                   sort_by_reference=True):
    """
    FIGURE 3 (per BIDS) -- the species' own reference, the cross-species
    consensus, and every individual run, in the same common-region space.

    Confirmed real problem this fixes: the previous version plotted only
    the POST-HOC group map (built from the 3dttest++ on the selected
    runs) as a single black "Reference" line -- which is NOT what the
    selection was actually based on, is DIFFERENT for every species (so
    it looked wrong that the black line "should be the same for
    everyone" -- it should, but wasn't, because it wasn't the shared
    reference at all), and can have GAPS where the thresholded group map
    lost coverage in a region the run-level profiles still had. Using
    PROFILE-based references (built the same way as the common-region
    list itself, so complete by construction) fixes both.

    Also confirmed real problem: the previous version of this function
    called the second reference "human_ref" and labeled it "Human
    reference" in the legend, on EVERY panel including non-human species
    -- leftover naming from an earlier version of the pipeline where
    selection really was anchored to a literal human reference. It is
    now the generic cross-species leave-one-out consensus (identical
    curve on every panel, human or not), so it is renamed consensus_ref
    / "Cross-species consensus" throughout. This is exactly why on the
    two Human panels the "Human reference" line didn't look like an
    average of the two Human datasets -- it wasn't one; it's the
    consensus across all 7 species, so it has no reason to sit between
    Castrillon's and Park's own curves at any given region.

    species_ref: this BIDS's OWN reference (kcore over ALL its own
    runs, kept and not-kept alike) -- "what does a typical run from this
    dataset look like".
    consensus_ref: the shared cross-species consensus every species was
    ranked against for selection (same curve on every panel).
    kept_ref / notkept_ref: optional -- kcore reference built from just
    the kept, or just the not-kept, runs of THIS BIDS (round 1). Lets
    you see whether the selection actually shifted the average, or
    whether kept and not-kept look about the same (weak selection).
    kept_ref_round2: optional -- kcore reference built from round 2's
    kept runs (species-anchored re-selection, see run_round2_selection):
    round 1's kept-only average re-ranked against itself rather than
    the cross-species consensus. Shown alongside kept_ref so you can see
    how much re-anchoring to the species' own signal moved the average.
    grand_mean_ref: optional -- unweighted median of every BIDS's own
    species_ref (each of the 10 BIDS counted once, unlike consensus_ref
    which collapses same-species BIDS together first).
    selected_keys: if given, runs in this set are drawn highlighted
    (green) and the rest faint grey, showing WHICH runs the selection
    actually kept relative to the overall spread.
    actual_selection_target: optional -- the EXACT vector selection
    actually scored this BIDS's runs against (sel["consensus_by_label"][label]
    from run_consensus_selection): leave-one-species-out AND z-scored,
    unlike consensus_ref above (all species included, always raw scale,
    identical curve on every panel -- confirmed real gap: consensus_ref
    was never literally what any panel's runs were judged against).
    Because z-scoring happens BEFORE the per-species median that builds
    the consensus, not after, this can differ from consensus_ref in
    SHAPE, not just scale -- rescaling consensus_ref would not
    reproduce it. Plotted on a separate right-hand axis (different
    units: z-scored) rather than forced onto the left axis's raw scale.
    Must be given in the ORIGINAL common_base_names order (the order
    this function received, before any sort_by_reference re-sort).
    """
    original_names = list(common_base_names)
    if sort_by_reference:
        common_base_names = _sort_regions_by_value(consensus_ref, common_base_names, value_col)
    region_reorder_idx = [original_names.index(r) for r in common_base_names]
    x = np.arange(len(common_base_names))
    selected_str = {tuple(str(v) for v in k) for k in (selected_keys or [])}

    fig, ax = plt.subplots(figsize=(max(11, len(common_base_names) * 0.24), 6))

    n_sel = n_other = 0
    for key, prof in profiles.items():
        is_sel = tuple(str(v) for v in key) in selected_str
        ax.plot(x, _ordered_common_values(prof, common_base_names, value_col),
               lw=0.8 if is_sel else 0.5,
               color="#009E73" if is_sel else "#BBBBBB",
               alpha=0.85 if is_sel else 0.35,
               zorder=3 if is_sel else 1)
        n_sel += is_sel
        n_other += (not is_sel)

    ax.plot(x, _ordered_common_values(species_ref, common_base_names, value_col),
           color="#1A5276", lw=2.6, zorder=5, label="Species reference (all runs)", ls="-")
    ax.plot(x, _ordered_common_values(consensus_ref, common_base_names, value_col),
           color="#D35400", lw=2.6, zorder=6, label="Cross-species consensus", ls="--")

    handles = [
        plt.Line2D([], [], color="#1A5276", lw=2.6, label="Species reference (all runs)"),
        plt.Line2D([], [], color="#D35400", lw=2.6, ls="--", label="Cross-species consensus"),
        plt.Line2D([], [], color="#009E73", lw=1.2, label=f"Selected runs (n={n_sel})"),
        plt.Line2D([], [], color="#BBBBBB", lw=1.2, label=f"Other runs (n={n_other})"),
    ]

    if kept_ref is not None:
        ax.plot(x, _ordered_common_values(kept_ref, common_base_names, value_col),
               color="#009E73", lw=2.0, zorder=7, ls=":")
        handles.append(plt.Line2D([], [], color="#009E73", lw=2.0, ls=":",
                                  label="Kept-only reference"))
    if notkept_ref is not None:
        ax.plot(x, _ordered_common_values(notkept_ref, common_base_names, value_col),
               color="#888888", lw=2.0, zorder=7, ls=":")
        handles.append(plt.Line2D([], [], color="#888888", lw=2.0, ls=":",
                                  label="Not-kept-only reference"))
    if kept_ref_round2 is not None:
        ax.plot(x, _ordered_common_values(kept_ref_round2, common_base_names, value_col),
               color="#CC79A7", lw=2.2, zorder=8, ls="-")
        handles.append(plt.Line2D([], [], color="#CC79A7", lw=2.2,
                                  label=final_round_label))
    if grand_mean_ref is not None:
        ax.plot(x, _ordered_common_values(grand_mean_ref, common_base_names, value_col),
               color="#000000", lw=1.8, zorder=4, ls="-.")
        handles.append(plt.Line2D([], [], color="#000000", lw=1.8, ls="-.",
                                  label="Grand mean (all 10 BIDS, unweighted)"))

    ax2 = None
    if actual_selection_target is not None:
        ax2 = ax.twinx()
        tvals = np.asarray(actual_selection_target, dtype=float)
        ax2.plot(x, tvals[region_reorder_idx], color="#8E44AD", lw=2.2, zorder=9, ls="-")
        ax2.set_ylabel(target_label, fontsize=9, color="#8E44AD")
        ax2.tick_params(axis="y", colors="#8E44AD")
        ax2.spines["top"].set_visible(False)
        handles.append(plt.Line2D([], [], color="#8E44AD", lw=2.2, label=target_label))

    ax.legend(handles=handles, fontsize=8, loc="upper left")

    ax.axhline(0, color="#999999", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(common_base_names, rotation=90, fontsize=4)
    ax.set_ylabel(value_col, fontsize=10)
    ax.set_title(title or "Species vs consensus reference, with individual runs (common regions)",
                fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_score_distributions(rankings, selections, title=None, output_path=None,
                             dpi=300, ncols=3):
    """
    Per BIDS: histogram of every run's consensus score, with the actual
    selection cutoff drawn. Shows whether the selected runs sit in a
    distinct upper mode (a real subgroup) or are just the right tail of
    one continuous distribution (an arbitrary cut) -- which matters a lot
    for how the selection can be described in a paper.
    """
    labels = sorted(rankings)
    if not labels:
        print("  [plot] nothing to plot")
        return None
    ncols = min(ncols, len(labels))
    nrows = int(np.ceil(len(labels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4.2, nrows * 3.2), squeeze=False)

    for i, label in enumerate(labels):
        ax = axes[i // ncols][i % ncols]
        vals = pd.to_numeric(rankings[label]["pearson_r"], errors="coerce").dropna()
        if vals.empty:
            ax.axis("off"); continue
        n_sel = len(selections.get(label, []))
        ax.hist(vals, bins=30, color="#4472c4", alpha=0.75, edgecolor="white")
        if n_sel and n_sel <= len(vals):
            cutoff = float(np.sort(vals)[::-1][n_sel - 1])
            ax.axvline(cutoff, color="#D55E00", lw=2.2)
            ax.text(cutoff, ax.get_ylim()[1] * 0.85, f" cutoff\n n={n_sel}",
                   color="#D55E00", fontsize=7, va="top")
        ax.set_title(f"{str(label).replace('_', ' ')} (n={len(vals)})", fontsize=9,
                     fontweight="bold")
        ax.set_xlabel("score vs consensus", fontsize=7)
        ax.set_ylabel("runs", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    for j in range(len(labels), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.suptitle(title or "Consensus score distribution per BIDS", fontsize=11, fontweight="bold")
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def _within_across_species_pairs(profiles_by_label, common_base_names, species_of_label,
                                 value_col="median_r", zscore=True, max_pairs=20000,
                                 random_seed=0):
    """
    Pairwise run-to-run profile correlations, with, for every pair, the two
    datasets and species it comes from. Returns a dict of arrays:
    r, da, db (dataset index), same_dataset, same_species, plus
    labels (dataset names) and species_of (species per dataset index).
    Pairs are subsampled above max_pairs.
    """
    rng = np.random.default_rng(random_seed)
    labels = sorted(profiles_by_label)
    vecs, ds = [], []
    for i, label in enumerate(labels):
        for key, prof in profiles_by_label[label].items():
            vals = _ordered_common_values(prof, common_base_names, value_col)
            if not np.isfinite(vals).any():
                continue
            vecs.append(_normalize_values(vals, zscore))
            ds.append(i)
    out = dict(labels=labels, species_of=np.array([species_of_label.get(l, l) for l in labels]),
               r=np.array([]), da=np.array([], int), db=np.array([], int))
    n = len(vecs)
    if n < 4:
        return out
    arr, ds = np.vstack(vecs), np.array(ds)
    ii, jj = np.triu_indices(n, k=1)
    if len(ii) > max_pairs:
        sel = rng.choice(len(ii), size=max_pairs, replace=False)
        ii, jj = ii[sel], jj[sel]
    r, keep = np.empty(len(ii)), np.zeros(len(ii), bool)
    for t, (i, j) in enumerate(zip(ii, jj)):
        both = np.isfinite(arr[i]) & np.isfinite(arr[j])
        if both.sum() >= 3:
            r[t] = np.corrcoef(arr[i][both], arr[j][both])[0, 1]
            keep[t] = True
    out.update(r=r[keep], da=ds[ii[keep]], db=ds[jj[keep]])
    return out


INTRA_PAIRS = ("other_dataset", "all")


def _within_across_stats(pairs, n_perm=2000, random_seed=0, intra="other_dataset"):
    """
    Three pair classes: same dataset / same species but another dataset /
    different species. The two classes compared (AUC, Wilcoxon, permutation)
    are set by `intra`:
      "other_dataset" (default): intra = same species, ANOTHER dataset;
          inter = different species; same-dataset pairs are left out of the
          test (they also share protocol, scanner and usually the animal, so
          they cannot separate species from dataset);
      "all": intra = every same-species pair, same dataset included (previous
          5b) -- dominated by the same-dataset pairs, which are the majority.
    Descriptive AUC (intra vs inter) plus a permutation p that respects the
    dependence between pairs: the species labels are shuffled ACROSS DATASETS
    (each dataset keeps all its runs), the pair correlations stay fixed, and
    median(intra) - median(inter) is recomputed with the same definition of
    intra. Pairs are not independent (each run enters hundreds of pairs), so a
    Mann-Whitney p on pairs is anti-conservative by orders of magnitude. With
    few datasets per species the permutation null is coarse and the test is
    weak: read the three medians, not only p.
    """
    if intra not in INTRA_PAIRS:
        raise ValueError(f"intra must be one of {INTRA_PAIRS}, got {intra!r}")
    r, da, db = pairs["r"], pairs["da"], pairs["db"]
    sp = pairs["species_of"]
    if len(r) < 4:
        return None
    same_ds = da == db
    same_sp = sp[da] == sp[db]
    intra_mask = (same_sp & ~same_ds) if intra == "other_dataset" else same_sp
    within, across = r[intra_mask], r[~same_sp]
    if len(within) < 2 or len(across) < 2:
        return None
    u, _ = _stats.mannwhitneyu(within, across, alternative="two-sided")
    auc = u / (len(within) * len(across))
    # Wilcoxon rank-sum (Mann-Whitney U), one-sided intra > inter, on the pairs:
    # requested as the test of the primary hypothesis. Pairs sharing a run are
    # not independent, so this p is optimistic; p_perm_dataset_level respects
    # the dependence.
    _, p_wilcoxon = _stats.mannwhitneyu(within, across, alternative="greater")
    obs = np.median(within) - np.median(across)
    rng = np.random.default_rng(random_seed)
    null = []
    for _ in range(n_perm):
        spp = rng.permutation(sp)
        w = spp[da] == spp[db]
        wi = (w & ~same_ds) if intra == "other_dataset" else w
        if not wi.any() or w.all():
            continue
        null.append(np.median(r[wi]) - np.median(r[~w]))
    null = np.asarray(null)
    p_perm = (1 + int(np.sum(null >= obs - 1e-12))) / (len(null) + 1) if len(null) else np.nan
    cross_ds_same_sp = r[same_sp & ~same_ds]
    return dict(intra_definition=intra,
                n_same_dataset=int(same_ds.sum()),
                n_same_species_other_dataset=int(len(cross_ds_same_sp)),
                n_different_species=int(len(across)),
                median_same_dataset=float(np.median(r[same_ds])) if same_ds.any() else np.nan,
                median_same_species_other_dataset=(float(np.median(cross_ds_same_sp))
                                                   if len(cross_ds_same_sp) else np.nan),
                median_different_species=float(np.median(across)),
                median_within=float(np.median(within)), median_across=float(np.median(across)),
                auc=float(auc), p_perm_dataset_level=float(p_perm), n_perm=int(len(null)),
                p_wilcoxon=float(p_wilcoxon), n_within=int(len(within)),
                n_across=int(len(across)))


_INTRA_TXT = {"other_dataset": "same species other dataset vs different species",
              "all": "intra (same species, any dataset) vs inter species"}
WITHIN_ACROSS_BINS = np.linspace(-1, 1, 41)     # 5b and F2/F2b: same bins


def _draw_within_across_panel(ax, pairs, title=None, paper=False, header=None,
                              n_runs=None, legend=True, intra="other_dataset"):
    """Three histograms (same dataset / same species other dataset /
    different species) + descriptive AUC and dataset-level permutation p
    (classes compared: `intra`, see _within_across_stats).
    paper=True: same drawing, but no bold title and no statistics in the
    legend; `header` (e.g. 'all runs') at the top left and 'n = .. runs,
    p = ..' at the top right -- the rest goes to the figure legend (.txt)."""
    st = _within_across_stats(pairs, intra=intra)
    if st is None:
        ax.axis("off")
        return None
    r, da, db, sp = pairs["r"], pairs["da"], pairs["db"], pairs["species_of"]
    same_ds, same_sp = da == db, sp[da] == sp[db]
    bins = WITHIN_ACROSS_BINS
    groups = [(r[same_ds], "#009E73", "same dataset", st["median_same_dataset"]),
              (r[same_sp & ~same_ds], INTRA_COLOR, "same species, other dataset",
               st["median_same_species_other_dataset"]),
              (r[~same_sp], INTER_COLOR, "different species", st["median_different_species"])]
    for vals, col, lab, med in groups:
        if len(vals):
            ax.hist(vals, bins=bins, alpha=0.55, color=col, density=True,
                    label=lab if paper else f"{lab} (n={len(vals)}, median={med:.3f})")
            ax.axvline(med, color=col, lw=2, ls="--")
    ax.set_xlabel("run-to-run profile correlation", fontsize=10)
    ax.set_ylabel("density", fontsize=10)
    if paper:
        if header:
            ax.set_title(header, loc="left", fontsize=10)
        p = st["p_perm_dataset_level"]
        p_txt = f"{p:.3f}" if p >= 0.001 else f"{p:.1e}"
        ax.set_title((f"n = {n_runs} runs, " if n_runs is not None else "") + f"p = {p_txt}",
                     loc="right", fontsize=10)
    else:
        ax.set_title((title + "\n" if title else "") +
                     f"{_INTRA_TXT[intra]}: AUC={st['auc']:.3f}\nWilcoxon p="
                     f"{st['p_wilcoxon']:.3g} (pairs), dataset-level permutation p="
                     f"{st['p_perm_dataset_level']:.3g}",
                     fontsize=10, fontweight="bold")
    if legend:
        ax.legend(fontsize=8 if paper else 7, loc="upper left", frameon=not paper)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    return st


def plot_intra_vs_inter_species_distance(profiles_by_label, common_base_names,
                                         species_of_label, value_col="median_r",
                                         zscore=True, title=None, output_path=None,
                                         dpi=300, max_pairs=20000, random_seed=0,
                                         intra="other_dataset"):
    """
    FIGURE F2: run-to-run profile correlation for pairs from the same
    dataset, from the same species but another dataset, and from different
    species. The middle class is the one that isolates species from
    dataset/protocol and animal identity; if it is not clearly above the
    different-species class, species is not the dominant structure.
    """
    pairs = _within_across_species_pairs(profiles_by_label, common_base_names,
                                         species_of_label, value_col=value_col, zscore=zscore,
                                         max_pairs=max_pairs, random_seed=random_seed)
    fig, ax = plt.subplots(figsize=(9, 5))
    st = _draw_within_across_panel(ax, pairs, title=title or "Within- vs across-species similarity",
                                   intra=intra)
    if st is None:
        plt.close(fig)
        print("  [plot] not enough pairs")
        return None, None
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    print(f"    medians: same dataset={st['median_same_dataset']:.4f}, same species other "
          f"dataset={st['median_same_species_other_dataset']:.4f}, different species="
          f"{st['median_different_species']:.4f}; {_INTRA_TXT[intra]}: AUC={st['auc']:.3f}, "
          f"dataset-level permutation p={st['p_perm_dataset_level']:.3g}")
    return fig, st


def _draw_intra_inter_paper(ax, pairs, st, test, header, legend, intra="other_dataset"):
    """Paper 5b: two distributions, dashed = medians, bracket with the stars of
    `test` ('wilcoxon' or 'permutation') between the medians. intra =
    'other_dataset': same species, another dataset (blue) vs different species
    (orange), same-dataset pairs not drawn; 'all': every same-species pair."""
    r, da, db, sp = pairs["r"], pairs["da"], pairs["db"], pairs["species_of"]
    same_sp = sp[da] == sp[db]
    intra_mask = (same_sp & (da != db)) if intra == "other_dataset" else same_sp
    intra_lab = ("same species, other dataset" if intra == "other_dataset"
                 else "intra-species pairs")
    inter_lab = "different species" if intra == "other_dataset" else "inter-species pairs"
    bins = WITHIN_ACROSS_BINS
    top = 0.0
    for vals, col, lab in ((r[~same_sp], INTER_COLOR, inter_lab),
                           (r[intra_mask], INTRA_COLOR, intra_lab)):
        h, _, _ = ax.hist(vals, bins=bins, density=True, histtype="stepfilled", color=col,
                          alpha=0.28, lw=0)
        ax.hist(vals, bins=bins, density=True, histtype="step", color=col, lw=0.8, label=lab)
        top = max(top, float(np.max(h)) if len(h) else 0.0)
    m_in, m_out = st["median_within"], st["median_across"]
    for m, col in ((m_out, INTER_COLOR), (m_in, INTRA_COLOR)):
        ax.plot([m, m], [0, top * 1.04], color=col, lw=0.8, ls=(0, (3, 2)))
    ax.set_ylim(0, top * 1.32)
    p = st["p_wilcoxon"] if test == "wilcoxon" else st["p_perm_dataset_level"]
    bracket(ax, min(m_in, m_out), max(m_in, m_out), top * 1.08, p, dy=top * 0.05)
    ax.set_xlim(-1, 1)
    ax.set_xlabel("Run-to-run profile correlation (r)")
    ax.set_ylabel("Density")
    ax.set_title(header, loc="left", fontweight="bold")
    if legend:
        ax.legend(loc="upper left", frameon=False)


def paper_within_across_figure(pairs_by, stats, names, test, intra, output_path):
    """5b / Figure 4b at print size (120 mm): one panel per pool in `names` ('all runs',
    'selected'), the two compared classes and the stars of `test`; pairs and statistics
    from _within_across_species_pairs / _within_across_stats."""
    with paper_style():
        fig, axes = plt.subplots(1, 2, figsize=(PAPER_WIDTH["onehalf"], 46 * MM))
        for k, (ax, name) in enumerate(zip(axes, names)):
            if stats[name] is None:
                ax.axis("off")
                continue
            _draw_intra_inter_paper(ax, pairs_by[name], stats[name], test,
                                    "all runs" if k == 0 else "selected", legend=k == 0,
                                    intra=intra)
        fig.tight_layout(w_pad=1.5)
        if output_path:
            save_paper(fig, output_path)
        else:
            plt.close(fig)


def plot_within_vs_across_all_vs_selected(all_runs_by_label, selected_by_label,
                                          common_base_names, species_of_label,
                                          value_col="median_r", zscore=True, title=None,
                                          output_path=None, dpi=300, max_pairs=20000,
                                          random_seed=0, selected_name="Selected runs",
                                          paper=False, test="wilcoxon",
                                          intra="other_dataset"):
    """F2 computed on every run (left) and on the selected runs (right).
    paper=True (5b): same pairs, same statistics (same `intra` definition),
    drawn at print size (120 mm) with the two compared classes only and the
    stars of `test` ('wilcoxon': rank-sum on the pairs, one-sided; or
    'permutation': species permuted across datasets); also saved as .pdf."""
    pools = [("All runs", all_runs_by_label), (selected_name, selected_by_label)]
    stats, pairs_by = {}, {}
    for name, pbl in pools:
        pairs_by[name] = _within_across_species_pairs(
            pbl, common_base_names, species_of_label, value_col=value_col, zscore=zscore,
            max_pairs=max_pairs, random_seed=random_seed)
        stats[name] = _within_across_stats(pairs_by[name], intra=intra)
    if paper:
        paper_within_across_figure(pairs_by, stats, [name for name, _ in pools], test, intra,
                                   output_path)
        return None, stats
    fig, axes = plt.subplots(1, 2, figsize=(16, 5.5))
    for ax, (name, pbl) in zip(axes, pools):
        n_runs = sum(len(p) for p in pbl.values())
        _draw_within_across_panel(ax, pairs_by[name], title=f"{name} (n={n_runs} runs)",
                                  intra=intra)
    fig.suptitle(title or "Within- vs across-species similarity: all runs vs selected",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    for name, st in stats.items():
        if st:
            print(f"    [{name}] {_INTRA_TXT[intra]}: {st['median_within']:.4f} vs "
                  f"{st['median_across']:.4f} (same dataset {st['median_same_dataset']:.4f}), "
                  f"AUC={st['auc']:.3f}, Wilcoxon p={st['p_wilcoxon']:.3g}, permutation p="
                  f"{st['p_perm_dataset_level']:.3g}")
    return fig, stats


def _draw_distance_divergence_panel(ax, distance_df, test_result=None, divergence_times=None,
                                    exclude_human=False, title=None, bids_colors=None):
    """
    Shared drawing routine behind plot_distance_vs_divergence and
    plot_distance_vs_divergence_multilevel, factored out so both draw
    identically instead of risking the two styles drifting apart.

    exclude_human MUST match whatever exclude_human was passed to
    test_phylogenetic_signal() to build test_result -- otherwise the
    drawn regression line and the printed r/p come from two different
    fits. (Confirmed real bug in the previous version of this function:
    fit_df = df[species != "Human"] was hardcoded regardless of what
    test_result was actually computed from, so once the caller started
    passing exclude_human=False the annotation still said "human
    excluded from fit" while the line kept excluding it anyway.)
    """
    divergence_times = divergence_times or DIVERGENCE_TIME_MYA
    df = distance_df.copy()
    df["divergence_mya"] = df["species"].map(divergence_times)
    df = df[df["divergence_mya"].notna() & df["distance_to_human"].notna()]
    if df.empty:
        ax.axis("off")
        return

    for _, row in df.iterrows():
        is_human = row["species"] == "Human"
        point_color = (bids_colors.get(row["label"], "#888888")
                       if bids_colors is not None and "label" in df.columns
                       else SPECIES_COLORS.get(row["species"], "#888888"))
        ax.scatter(row["divergence_mya"], row["distance_to_human"],
                  s=140, color=point_color,
                  edgecolors="black", lw=1.0, zorder=3,
                  facecolors="none" if (is_human and exclude_human) else None,
                  label=row["species"])
        # Multiple BIDS of the same species (Human, Macaque, Marmoset) would
        # otherwise print the exact same label text on top of each other --
        # append the BIDS name whenever this df has more than one row for
        # that species (a species-level df, one row per species, never
        # triggers this: label is then identical to species anyway).
        same_species_n = (df["species"] == row["species"]).sum()
        label_txt = (f"{row['species']} ({row['label']})"
                    if "label" in df.columns and same_species_n > 1
                    else row["species"])
        ax.annotate(label_txt, (row["divergence_mya"], row["distance_to_human"]),
                   textcoords="offset points", xytext=(7, 5), fontsize=8)

    fit_df = df[df["species"] != "Human"] if exclude_human else df
    if len(fit_df) >= 3:
        slope, intercept, r, p, se = _stats.linregress(
            fit_df["divergence_mya"], fit_df["distance_to_human"])
        xs = np.linspace(fit_df["divergence_mya"].min(), fit_df["divergence_mya"].max(), 100)
        ax.plot(xs, intercept + slope * xs, "k--", lw=1.4, alpha=0.7, zorder=2)

    if test_result:
        human_note = "human excluded from fit" if exclude_human else "human included in fit"
        ax.text(0.03, 0.97,
               f"Pearson r={test_result['pearson_r']:.3f}, p={test_result['pearson_p']:.3f}\n"
               f"Spearman rho={test_result['spearman_rho']:.3f}, p={test_result['spearman_p']:.3f}\n"
               f"n={test_result['n']} species ({human_note})",
               transform=ax.transAxes, va="top", fontsize=8,
               bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))

    ax.set_xlabel("Divergence time from human (Mya)", fontsize=9)
    ax.set_ylabel("Distance to leave-one-species-out consensus", fontsize=9)
    if title:
        ax.set_title(title, fontsize=10, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_distance_vs_divergence(distance_df, test_result=None, divergence_times=None,
                                exclude_human=False, title=None, output_path=None, dpi=300):
    """
    Scatter of distance-to-consensus against divergence time, one point
    per species, colored by species, with a fitted line and the
    correlation stats annotated.

    exclude_human: pass the SAME value used for test_phylogenetic_signal
    when building test_result (see _draw_distance_divergence_panel).
    Default False matches compute_loo_species_distances' leave-one-
    SPECIES-out reference, which makes Human's own distance non-
    circular -- there is no longer a reason to exclude it by default.
    """
    dt = divergence_times or DIVERGENCE_TIME_MYA
    check = distance_df.copy()
    check["divergence_mya"] = check["species"].map(dt)
    if check["divergence_mya"].notna().sum() == 0 or distance_df.empty:
        print("  [plot] nothing to plot")
        return None

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    bids_colors = (_palette_color_map(distance_df["label"])
                  if "label" in distance_df.columns else None)
    _draw_distance_divergence_panel(
        ax, distance_df, test_result, divergence_times, exclude_human=exclude_human,
        title=title or "Distance to leave-one-species-out consensus vs divergence",
        bids_colors=bids_colors)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_distance_vs_divergence_multilevel(results_by_level, divergence_times=None,
                                           output_path=None, dpi=300):
    """
    One panel per atlas level (results_by_level: {level: {"ddf_species":
    df, "test_result": dict, "n_common_regions": int}}), same drawing
    routine as plot_distance_vs_divergence, to compare how segmentation
    granularity changes the phylogenetic slope/fit.
    """
    levels = sorted(results_by_level)
    n = len(levels)
    if n == 0:
        print("  [plot] nothing to plot")
        return None
    n_cols = min(2, n)
    n_rows = math.ceil(n / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7 * n_cols, 5 * n_rows), squeeze=False)
    axes_flat = axes.ravel()
    for ax, level in zip(axes_flat, levels):
        res = results_by_level[level]
        _draw_distance_divergence_panel(
            ax, res["ddf_species"], res.get("test_result"), divergence_times,
            exclude_human=False,
            title=(f"Atlas level {level} ({res.get('n_common_regions', '?')} common regions)"
                   + (f"\nMantel r={res['mantel']['mantel_pearson']['r']:+.2f}, "
                      f"exact p(>)={res['mantel']['mantel_pearson']['p_greater']:.3f}"
                      if res.get("mantel") else "")))
    for ax in axes_flat[n:]:
        ax.set_visible(False)
    fig.suptitle("Impact of segmentation level on distance-to-consensus vs divergence",
                fontsize=12, fontweight="bold")
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_run_distance_distribution(dist_df, value_col="pearson_r", output_path=None,
                                   title=None, dpi=300):
    """
    One violin+jittered-strip per BIDS, showing every run's correlation
    (or distance) to its BIDS' reference map. Self-contained rather than
    calling EDNiX_figures' private _violin_strip_quartiles (a different
    package), but same visual spirit: grey violin backdrop, colored
    jittered points, quartile ticks.
    """
    if dist_df is None or dist_df.empty:
        print("  [plot] nothing to plot")
        return None
    labels = sorted(dist_df["label"].unique())
    bids_colors = _palette_color_map(labels)
    fig, ax = plt.subplots(figsize=(max(6, 1.1 * len(labels)), 5))
    rng = np.random.default_rng(42)
    for xi, lab in enumerate(labels):
        sub = dist_df.loc[dist_df["label"] == lab, [value_col] + (
            ["is_kept"] if "is_kept" in dist_df.columns else [])].dropna(subset=[value_col])
        vals = sub[value_col].values
        if len(vals) == 0:
            continue
        color = bids_colors[lab]
        if len(vals) >= 2:
            parts = ax.violinplot(vals, positions=[xi], widths=0.7, showextrema=False)
            for pc in parts["bodies"]:
                pc.set_facecolor("#dddddd")
                pc.set_alpha(0.6)
        jitter = rng.uniform(-0.15, 0.15, len(vals))
        ax.scatter(xi + jitter, vals, s=22, color=color, edgecolors="black",
                  lw=0.3, alpha=0.85, zorder=3)
        if "is_kept" in sub.columns:
            kept_mask = sub["is_kept"].fillna(False).values.astype(bool)
            if kept_mask.any():
                ax.scatter((xi + jitter)[kept_mask], vals[kept_mask], marker="x", s=70,
                          color="black", linewidths=1.6, zorder=5)
        for q in (25, 50, 75):
            ax.hlines(np.percentile(vals, q), xi - 0.25, xi + 0.25,
                     color="black", lw=1.4 if q == 50 else 0.8, zorder=4)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([str(l).replace("_", " ") for l in labels], rotation=35, ha="right",
                       fontsize=8)
    ax.set_ylabel("Pearson r to BIDS reference" if value_col == "pearson_r"
                 else "Distance (1-r) to BIDS reference", fontsize=10)
    ax.set_title(title or "Per-run distance to the kept-group reference map",
                fontsize=11, fontweight="bold")
    if "is_kept" in dist_df.columns and dist_df["is_kept"].any():
        ax.legend(handles=[plt.Line2D([], [], marker="x", color="black", lw=0,
                                      markersize=9, markeredgewidth=1.6,
                                      label="Selected (kept) run")],
                  fontsize=8, loc="best")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_run_distance_std_by_bids(dist_df, value_col="pearson_r", output_path=None, dpi=300):
    """
    Bar chart: SD of each BIDS' per-run correlation to its reference --
    a compact single-number view of within-BIDS run-to-run consistency,
    companion to plot_run_distance_distribution's full spread.
    """
    if dist_df is None or dist_df.empty:
        print("  [plot] nothing to plot")
        return None
    summary = (dist_df.groupby(["label", "species"], as_index=False)[value_col]
              .agg(std="std", n="count")).sort_values("std", ascending=False)
    bids_colors = _palette_color_map(summary["label"])
    fig, ax = plt.subplots(figsize=(max(6, 1.1 * len(summary)), 4.5))
    colors = [bids_colors[lab] for lab in summary["label"]]
    ax.bar(summary["label"], summary["std"], color=colors, edgecolor="black", lw=0.5)
    for i, (n, s) in enumerate(zip(summary["n"], summary["std"])):
        ax.text(i, s, f"n={n}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(len(summary)))
    ax.set_xticklabels(summary["label"], rotation=35, ha="right", fontsize=8)
    ax.set_ylabel(f"SD of per-run {value_col} to reference", fontsize=10)
    ax.set_title("Run-to-run consistency by BIDS", fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_qc_distance_correlation(meta_df, output_path=None, dpi=300, alpha=0.05, max_metrics=40):
    """
    Bar chart of the meta-analytic r (meta_analyze_qc_correlations) between
    each QC metric and distance to reference, for EVERY metric tested (top
    max_metrics by |r|). Metrics significant after Benjamini-Hochberg
    correction (meta_q_fdr < alpha) are coloured (blue r>0, orange r<0);
    the others are grey. Previously only significant metrics were drawn,
    so a run where nothing passed produced no figure at all -- which reads
    as a crash and hides a genuine null result.
    """
    if meta_df is None or meta_df.empty:
        print("  [plot] nothing to plot")
        return None
    pcol = "meta_q_fdr" if "meta_q_fdr" in meta_df.columns else "meta_p"
    sub = meta_df.reindex(meta_df["meta_r"].abs().sort_values(ascending=False).index)
    sub = sub.head(max_metrics).reset_index(drop=True)
    sig = sub[pcol] < alpha
    colors = [("#0072B2" if r > 0 else "#D55E00") if g else "#BBBBBB"
              for r, g in zip(sub["meta_r"], sig)]
    fig, ax = plt.subplots(figsize=(max(6, 0.45 * len(sub)), 5))
    ax.bar(range(len(sub)), sub["meta_r"], color=colors, edgecolor="white")
    for i, (r, n_bids, q) in enumerate(zip(sub["meta_r"], sub["n_bids"], sub[pcol])):
        ax.text(i, r + (0.02 if r >= 0 else -0.02), f"n={n_bids}\nq={q:.2g}", ha="center",
                va="bottom" if r >= 0 else "top", fontsize=6)
    ax.axhline(0, color="#888888", lw=0.8)
    ax.set_xticks(range(len(sub)))
    ax.set_xticklabels(sub["metric"], rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("Meta-analytic r (Fisher-z weighted across BIDS)", fontsize=10)
    lab = "BH-FDR q" if pcol == "meta_q_fdr" else "meta p"
    ax.set_title(f"QC metrics vs distance to reference -- {int(sig.sum())}/{len(meta_df)} "
                 f"significant across BIDS (coloured: {lab}<{alpha})",
                 fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_phylogeny_permutation(perm_result, title=None, output_path=None, dpi=300):
    """Null distribution of the distance-vs-divergence correlation, with
    the observed value marked."""
    if perm_result is None:
        return None
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(perm_result["null_r"], bins=50, color="#888888", alpha=0.75, edgecolor="white")
    ax.axvline(perm_result["observed_r"], color="#D55E00", lw=2.5,
              label=f"observed r={perm_result['observed_r']:+.3f} "
                    f"(p={perm_result['p_perm']:.4f})")
    ax.set_xlabel("Pearson r (distance vs divergence time)", fontsize=10)
    ax.set_ylabel("permutations", fontsize=10)
    ax.set_title(title or f"Phylogeny permutation test "
                          f"(n={perm_result['n_species']} species)",
                fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig