"""
Supplementary Figure 10 -- region-by-region trajectories of the Figure 3 composites.

Does one region drive the effect of a composite (Fig. 3 a-f, Supp. Fig. 9), or does the
effect follow the category (association vs primary)? Same subjects, sessions, hemisphere
and models as Figure 3; every region of both categories, both species, three measures.

  a-c / d-f  Human / Macaque: spline-LMM trajectory of every region, in % of the value
             fitted at the youngest age (log scale: a region's size does not weigh in);
             association regions in blues, primary regions in oranges
  (archive, not in the figure since 2026-11-24; numbers kept in the stats)
  g-i / j-l  Human / Macaque: annual change (%/year, 95% CI) of a linear mixed model on
             log values, in two blocks: top, every region alone (circles); bottom, every
             composite (squares) and every composite without one region (open diamonds,
             "Primary without Striate (V1)", leave-one-region-out, LORO)

Annual change = 100 * (exp(b) - 1), b = slope per year of log(value) ~ age + (1|subject)
(ML), CI from b +/- 1.96 SE. Single-region P values: BH-FDR over the regions of one
species x measure. LORO only where >= 2 regions remain: association has 2 regions, so
"association without one region" is the other region alone, already in the top block.
An estimate that cannot be computed carries its reason (region absent, values missing
or <= 0, < 2 sessions, model failed) in the stats.

  region_stats()   every number (paper_stats.supp10)
  draw_regions()   the drawing (paper_figures.supp10)
"""
import os

import numpy as np
import pandas as pd

from Paper.analysis.library import longitudinal as LG

LABELS = {"PFC": "Prefrontal", "Post_parietal_Cx": "Post. parietal",
          "Striate_Cx": "Striate (V1)", "Auditory_Cx": "Auditory",
          "Somatosensory_Cx": "Somatosensory", "Motor_premotor": "Motor / premotor"}
CATEGORY_OF = {**{r: "Association" for r in LG.ASSOCIATION_REGIONS},
               **{r: "Primary" for r in LG.PRIMARY_REGIONS}}
# shades of the Fig. 3 g-i colours, one per region
_SHADES = {"Association": ["#1B4F72", "#5DADE2"],
           "Primary": ["#A04000", "#E67E22", "#F5B041", "#F8C471"]}
COLORS = {r: _SHADES[c][i] for c in ("Association", "Primary")
          for i, r in enumerate([x for x in CATEGORY_OF if CATEGORY_OF[x] == c])}


def _bh(p):
    """Benjamini-Hochberg adjusted P (NaN kept)."""
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    m = ok.sum()
    if not m:
        return out
    q = p[ok]
    o = np.argsort(q)
    adj = q[o] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    res = np.empty(m)
    res[o] = np.clip(adj, 0, 1)
    out[ok] = res
    return out


def _long_rows(morph, mod, name, regions, age_df, age_col, hemisphere, min_sessions):
    """Sessions of one region / composite merged with age, subjects with >= min_sessions;
    log of the measure in 'log_y'. (DataFrame, composite info, reason when empty)."""
    from Statistics.Group_fMRI.EDNiX_longitudinal import prepare_longitudinal
    col = LG._FIG3D_METRICS[mod][0]
    rows, info = LG.composite_region_rows(morph, mod, name, regions=regions,
                                          hemisphere=hemisphere)
    if rows.empty:
        why = ("region(s) absent: " + ", ".join(info["missing"]) if info["missing"]
               else "no finite value (all NaN or a region missing in every session)")
        return pd.DataFrame(), info, why
    n_bad = int((rows[col] <= 0).sum())
    rows = rows[rows[col] > 0]
    if rows.empty:
        return pd.DataFrame(), info, f"all {n_bad} values <= 0"
    m = prepare_longitudinal(rows, age_df, age_col, atlas_level=info["level"],
                             hemisphere=hemisphere, regions=[name],
                             min_sessions=min_sessions, label=f"SuppRegions/{name}_{mod}")
    if m.empty:
        return m, info, f"< {min_sessions} sessions with age" + (
            f" ({n_bad} values <= 0 dropped)" if n_bad else "")
    m = m.copy()
    m["log_y"] = np.log(m[col])
    return m, info, (f"{n_bad} values <= 0 dropped" if n_bad else "")


def _annual(m, age_col, min_sessions, why=""):
    import warnings
    from Statistics.Group_fMRI.EDNiX_longitudinal import fit_lme
    base = dict(pct_per_year=np.nan, ci_low=np.nan, ci_high=np.nan, p=np.nan,
                n_subjects=m["subject"].nunique() if not m.empty else 0, n_sessions=len(m),
                note=why)
    if m.empty:
        return base
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        lm = fit_lme(m, "log_y", age_col, min_sessions)
    if not lm:
        wsd = m.groupby("subject")["log_y"].std().median()
        msg = next((str(x.message) for x in w if "fit_lme" in str(x.message)), "no estimate")
        base["note"] = (f"model failed ({msg}); median within-subject SD of log value = "
                        f"{wsd:.2g}" + (f"; {why}" if why else ""))
        return base
    b, se = lm["slope_per_unit"], lm["se_per_unit"]
    return dict(pct_per_year=100 * np.expm1(b), ci_low=100 * np.expm1(b - 1.96 * se),
                ci_high=100 * np.expm1(b + 1.96 * se), p=lm["p_slope"],
                n_subjects=m["subject"].nunique(), n_sessions=len(m), note=why)


def region_stats(species_rows, modalities=None, hemisphere="left",
                 min_sessions=LG.MIN_SESSIONS, letters="abcdefghijkl"):
    """Every number of Supp. Fig. 10. species_rows = longitudinal.load_figure3_data().
    Returns dict(curves, table (one row per estimate), lines, modalities, species, labels)."""
    from Statistics.Group_fMRI.EDNiX_longitudinal import fit_gamm
    modalities = LG.FIG3_MODALITIES if modalities is None else modalities
    regions = LG.ASSOCIATION_REGIONS + LG.PRIMARY_REGIONS
    curves, rows, lines = {}, [], []
    for sp, morph, age_df, age_col, label in species_rows:
        for mod, _, _, title in modalities:
            block = []
            for reg in regions:                                   # single regions
                m, info, why = _long_rows(morph, mod, reg, [reg], age_df, age_col,
                                          hemisphere, min_sessions)
                g = (fit_gamm(m, "log_y", age_col, min_sessions, spline_df=LG.FIG3_SPLINE_DF)
                     if not m.empty else None)
                if g:
                    curves[(sp, mod, reg)] = dict(age=g["age_grid"],
                                                  pct=100 * np.expm1(g["pred"] - g["pred"][0]))
                block.append(dict(species=sp, modality=mod, kind="region",
                                  category=CATEGORY_OF[reg], estimate=reg, left_out="",
                                  label=LABELS.get(reg, reg),
                                  gamm_p_age=g["p_age"] if g else np.nan,
                                  **_annual(m, age_col, min_sessions, why)))
            q = _bh([r["p"] for r in block])
            for r, qq in zip(block, q):
                r["p_fdr"] = qq
            for cat, regs in (("Association", LG.ASSOCIATION_REGIONS),
                              ("Primary", LG.PRIMARY_REGIONS)):   # composites + LORO
                outs = [None] + (list(regs) if len(regs) > 2 else [])
                for out in outs:
                    keep = [r for r in regs if r != out]
                    m, _, why = _long_rows(morph, mod, cat, keep, age_df, age_col,
                                           hemisphere, min_sessions)
                    block.append(dict(species=sp, modality=mod,
                                      kind="composite" if out is None else "loro",
                                      category=cat,
                                      estimate=cat if out is None else f"{cat} without {out}",
                                      label=(f"{cat} (all regions)" if out is None else
                                             f"{cat} without {LABELS.get(out, out)}"),
                                      left_out=out or "", regions=", ".join(keep),
                                      **_annual(m, age_col, min_sessions, why)))
            rows.extend(block)
            lines.extend(_block_lines(sp, title, block))
    table = pd.DataFrame(rows)
    return dict(curves=curves, table=table, lines=lines, modalities=modalities,
                species=[r[0] for r in species_rows],
                labels={r[0]: r[4] for r in species_rows}, letters=letters)


def _block_lines(sp, title, block):
    reg = [r for r in block if r["kind"] == "region" and np.isfinite(r["pct_per_year"])]
    out = [f"{sp} {title}:"]
    for r in block:
        if not np.isfinite(r["pct_per_year"]):
            out.append(f"  {r['label']}: NO ESTIMATE ({r['note'] or 'unknown'})")
            continue
        extra = f", P_FDR={r['p_fdr']:.3g}" if r["kind"] == "region" else ""
        out.append(f"  {r['label']}: {r['pct_per_year']:+.3f} %/year "
                   f"[{r['ci_low']:+.3f}, {r['ci_high']:+.3f}], P={r['p']:.3g}{extra}; "
                   f"n={r['n_subjects']} subjects, {r['n_sessions']} sessions"
                   + (f" [{r['note']}]" if r["note"] else ""))
    missing = [r["label"] for r in block if r["kind"] == "region"
               and not np.isfinite(r["pct_per_year"])]
    a = [r["pct_per_year"] for r in reg if r["category"] == "Association"]
    p = [r["pct_per_year"] for r in reg if r["category"] == "Primary"]
    if a and p:
        if max(a) < min(p):
            sep = "every association region has a lower annual change than every primary region"
        elif min(a) > max(p):
            sep = "every association region has a higher annual change than every primary region"
        else:
            sep = "association and primary single-region slopes OVERLAP"
        out.append(f"  ranges: association [{min(a):+.3f}, {max(a):+.3f}], primary "
                   f"[{min(p):+.3f}, {max(p):+.3f}] %/year: {sep}"
                   + (f" (WITHOUT {', '.join(missing)}: no estimate)" if missing else ""))
    for cat in ("Association", "Primary"):
        regs = [r for r in reg if r["category"] == cat]
        allr = next((r for r in block if r["kind"] == "composite" and r["category"] == cat), None)
        if allr is None or not np.isfinite(allr["pct_per_year"]) or not regs:
            continue
        sign = np.sign(allr["pct_per_year"])
        against = [r["label"] for r in regs if np.sign(r["pct_per_year"]) != sign]
        out.append(f"  {cat}: composite {allr['pct_per_year']:+.3f} %/year; regions changing in "
                   f"the opposite direction: {', '.join(against) or 'none'}")
        lo = [r for r in block if r["kind"] == "loro" and r["category"] == cat
              and np.isfinite(r["pct_per_year"])]
        if lo:
            flips = [LABELS.get(r["left_out"], r["left_out"]) for r in lo
                     if np.sign(r["pct_per_year"]) != sign]
            lost = [LABELS.get(r["left_out"], r["left_out"]) for r in lo
                    if allr["p"] < 0.05 <= r["p"]]
            out.append(f"    leave-one-region-out: range [{min(r['pct_per_year'] for r in lo):+.3f}, "
                       f"{max(r['pct_per_year'] for r in lo):+.3f}] %/year; sign change without: "
                       f"{', '.join(flips) or 'none'}; P >= 0.05 without: {', '.join(lost) or 'none'}")
    return out


def stats_text(st):
    return ("Supplementary Figure 10 -- region-by-region trajectories of the Figure 3 "
            "composites (left hemisphere)\n"
            "Annual change = 100*(exp(b)-1), b = slope/year of log(value) ~ age + (1|subject), ML; "
            "95% CI from b +/- 1.96 SE; single regions: BH-FDR over the 6 regions of one "
            "species x measure. Leave-one-region-out only where >= 2 regions remain (primary); "
            "association has 2 regions: without one = the other region alone.\n\n"
            + "\n".join(st["lines"]) + "\n")


def draw_regions(st, out_path, panel_w=3.4, panel_h=2.7, archive_annual=True):
    """Supp. Fig. 10 = a-f, the region curves (rows Human, Macaque; columns = Figure 3
    measures). The annual-change panels (former g-l, 2026-11-24: removed from the figure
    for now) are drawn alone in <dir>/archive/<name>_annual_change_ARCHIVE.png/.pdf when
    archive_annual. Returns the paths."""
    paths = _draw(st, out_path, panel_w, panel_h, curves=True, annual=False)
    if archive_annual:
        d, name = os.path.split(out_path)
        paths += _draw(st, os.path.join(d, "archive", f"{name}_annual_change_ARCHIVE"),
                       panel_w, panel_h, curves=False, annual=True)
    return paths


def _draw(st, out_path, panel_w, panel_h, curves=True, annual=True):
    """Region curves (one row per species) and / or annual-change panels (one row per
    species), columns = Figure 3 measures."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.lines import Line2D
    from Tools.ednix_bids_tools import PAPER_RC

    mods, species, tab = st["modalities"], st["species"], st["table"]
    nrows, ncols = (int(curves) + int(annual)) * len(species), len(mods)
    r0 = len(species) if curves else 0                            # first annual-change row
    regions = LG.ASSOCIATION_REGIONS + LG.PRIMARY_REGIONS
    rc = dict(PAPER_RC); rc.update(LG._FIG_RC)
    with plt.rc_context(rc):
        fig = plt.figure(figsize=(panel_w * ncols + 1.0,
                                  panel_h * nrows * 0.95 + (0.5 if annual else 0)))
        gs = gridspec.GridSpec(nrows, ncols, figure=fig, hspace=0.45,
                               wspace=0.75 if annual else 0.42,
                               left=0.15 if annual else 0.11, right=0.99,
                               top=0.96 if nrows > 2 else 0.92,
                               bottom=(0.07 if nrows > 2 else 0.14) if annual else 0.1)
        letter, empty, row_label = LG._panel_tools(fig)
        for i, sp in enumerate(species if curves else []):        # trajectories
            for c, (mod, _, _, title) in enumerate(mods):
                ax = fig.add_subplot(gs[i, c])
                letter(ax)
                drawn = False
                for reg in regions:
                    cv = st["curves"].get((sp, mod, reg))
                    if cv is None:
                        continue
                    ax.plot(cv["age"], cv["pct"], color=COLORS[reg], lw=1.4,
                            ls="-" if CATEGORY_OF[reg] == "Association" else "--")
                    drawn = True
                if not drawn:
                    empty(ax, "no data")
                    continue
                ax.axhline(0, color="0.5", lw=0.6, ls=":")
                ax.set_xlabel("Age (years)")
                ax.set_ylabel(f"% {title.lower()} change")
                ax.spines[["top", "right"]].set_visible(False)
                if i == 0:
                    ax.set_title(title, fontweight="bold", pad=12)
                if c == 0:
                    row_label(ax, f"{st['labels'][sp]}\neach region")
                if i == 0 and c == ncols - 1:
                    h = [Line2D([], [], color=COLORS[r],
                                ls="-" if CATEGORY_OF[r] == "Association" else "--",
                                label=f"{LABELS.get(r, r)} ("
                                      f"{'assoc.' if CATEGORY_OF[r] == 'Association' else 'primary'})")
                         for r in regions]
                    ax.legend(handles=h, fontsize=5.5, frameon=False, loc="best")
        for i, sp in enumerate(species if annual else []):        # annual change
            for c, (mod, _, _, title) in enumerate(mods):
                ax = fig.add_subplot(gs[r0 + i, c])
                letter(ax)
                b = tab[(tab["species"] == sp) & (tab["modality"] == mod)] if len(tab) else tab
                if b.empty:
                    empty(ax, "no data")
                    continue
                top = [b[(b.kind == "region") & (b.estimate == r)] for r in regions]
                top = pd.concat([t for t in top if len(t)])
                bottom = pd.concat([b[(b.kind != "region") & (b.category == cat)]
                                    for cat in ("Association", "Primary")])
                n_top, n_bot = len(top), len(bottom)
                ys_top = np.arange(n_top)[::-1] + n_bot + 1.5     # gap between the blocks
                ys_bot = np.arange(n_bot)[::-1]
                ticks, labels = [], []
                for ys, blk in ((ys_top, top), (ys_bot, bottom)):
                    for yy, (_, r) in zip(ys, blk.iterrows()):
                        ticks.append(yy)
                        labels.append(r["label"])
                        if not np.isfinite(r["pct_per_year"]):
                            ax.text(0.02, yy, "no estimate", fontsize=5, color="0.5",
                                    va="center", transform=ax.get_yaxis_transform())
                            continue
                        col = (COLORS[r["estimate"]] if r["kind"] == "region"
                               else LG.FIG3D_COLORS[r["category"]])
                        mk = {"region": "o", "composite": "s", "loro": "D"}[r["kind"]]
                        face = "white" if r["kind"] == "loro" else col
                        ax.plot([r["ci_low"], r["ci_high"]], [yy, yy], color=col, lw=1)
                        ax.plot(r["pct_per_year"], yy, mk, color=col, mfc=face, ms=3.5, mew=0.9)
                ax.axhline(n_bot + 0.75, color="0.7", lw=0.6)
                ax.set_yticks(ticks)
                ax.set_yticklabels(labels, fontsize=5.5)
                ax.set_ylim(-0.8, n_bot + n_top + 1.0)
                ax.axvline(0, color="0.5", lw=0.6, ls=":")
                ax.set_xlabel("Annual change (%/year, 95% CI)")
                if i == 0 and not curves:
                    ax.set_title(title, fontweight="bold", pad=12)
                ax.spines[["top", "right"]].set_visible(False)
                if c == 0:
                    ax.annotate(f"{st['labels'][sp]}\nannual change", xy=(-0.95, 0.5),
                                xycoords="axes fraction", rotation=90, ha="center",
                                va="center", fontsize=9, fontweight="bold")
        h = [] if not annual else [
             Line2D([], [], marker="o", color="0.3", ls="", ms=4, label="one region"),
             Line2D([], [], marker="s", color="0.3", ls="", ms=4,
                    label="composite, all regions of the category"),
             Line2D([], [], marker="D", color="0.3", mfc="white", ls="", ms=4,
                    label="composite without one region (leave-one-region-out)")]
        if h:
            fig.legend(handles=h, loc="lower center", ncol=3, fontsize=7, frameon=False,
                       bbox_to_anchor=(0.55, 0.0))
        paths = LG._save(fig, out_path)
        plt.close(fig)
    return paths


def plot_regions(species_rows, out_path, **kw):
    """region_stats + draw_regions + <out_path>_stats.txt (stand-alone entry point)."""
    st = region_stats(species_rows, **kw)
    paths = draw_regions(st, out_path)
    with open(out_path + "_stats.txt", "w") as fh:
        fh.write(stats_text(st))
    return paths + [out_path + "_stats.txt"]
