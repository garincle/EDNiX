"""
Supplementary Figure 10 -- region-by-region trajectories of the Figure 3 composites.

Does one region drive the effect of a composite (Fig. 3 a-f, Supp. Fig. 9), or does the
effect follow the category (association vs primary)? Same subjects, sessions, hemisphere
and models as Figure 3; every region of both categories, both species, three measures.

  a-c / d-f  Human / Macaque: spline-LMM trajectory of every region, in % of the value
             fitted at the youngest age (log scale: a region's size does not weigh in);
             association regions in blues, primary regions in oranges
  g-i / j-l  Human / Macaque: annual change (%/year, 95% CI) of a linear mixed model on
             log values, for every single region, every composite, and every composite
             with one region left out (leave-one-region-out, LORO)

Annual change = 100 * (exp(b) - 1), b = slope per year of log(value) ~ age + (1|subject)
(ML), CI from b +/- 1.96 SE. Single-region P values: BH-FDR over the regions of one
species x measure. With two association regions, leaving one out leaves a single region:
those LORO rows equal the single-region rows of the other region (said in the stats).

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
    log of the measure in 'log_y'. (DataFrame, composite info)."""
    from Statistics.Group_fMRI.EDNiX_longitudinal import prepare_longitudinal
    rows, info = LG.composite_region_rows(morph, mod, name, regions=regions,
                                          hemisphere=hemisphere)
    if rows.empty:
        return pd.DataFrame(), info
    m = prepare_longitudinal(rows, age_df, age_col, atlas_level=info["level"],
                             hemisphere=hemisphere, regions=[name],
                             min_sessions=min_sessions, label=f"SuppRegions/{name}_{mod}")
    if m.empty:
        return m, info
    col = LG._FIG3D_METRICS[mod][0]
    m = m[m[col] > 0].copy()
    m["log_y"] = np.log(m[col])
    return m, info


def _annual(m, age_col, min_sessions):
    from Statistics.Group_fMRI.EDNiX_longitudinal import fit_lme
    lm = fit_lme(m, "log_y", age_col, min_sessions) if not m.empty else {}
    if not lm:
        return dict(pct_per_year=np.nan, ci_low=np.nan, ci_high=np.nan, p=np.nan,
                    n_subjects=m["subject"].nunique() if not m.empty else 0, n_sessions=len(m))
    b, se = lm["slope_per_unit"], lm["se_per_unit"]
    return dict(pct_per_year=100 * np.expm1(b), ci_low=100 * np.expm1(b - 1.96 * se),
                ci_high=100 * np.expm1(b + 1.96 * se), p=lm["p_slope"],
                n_subjects=m["subject"].nunique(), n_sessions=len(m))


def region_stats(species_rows, modalities=None, hemisphere="left",
                 min_sessions=LG.MIN_SESSIONS, letters="abcdefghijkl"):
    """Every number of Supp. Fig. 10. species_rows = longitudinal.load_figure3_data().
    Returns dict(curves, table (one row per estimate), lines, nrows/ncols, species)."""
    from Statistics.Group_fMRI.EDNiX_longitudinal import fit_gamm
    modalities = LG.FIG3_MODALITIES if modalities is None else modalities
    regions = LG.ASSOCIATION_REGIONS + LG.PRIMARY_REGIONS
    curves, rows, lines = {}, [], []
    for sp, morph, age_df, age_col, label in species_rows:
        for mod, _, _, title in modalities:
            block = []
            for reg in regions:                                   # single regions
                m, info = _long_rows(morph, mod, reg, [reg], age_df, age_col, hemisphere,
                                     min_sessions)
                g = fit_gamm(m, "log_y", age_col, min_sessions) if not m.empty else None
                if g:
                    curves[(sp, mod, reg)] = dict(age=g["age_grid"],
                                                  pct=100 * np.expm1(g["pred"] - g["pred"][0]))
                block.append(dict(species=sp, modality=mod, kind="region", category=CATEGORY_OF[reg],
                                  estimate=reg, left_out="", found=bool(info["matched"]),
                                  gamm_p_age=g["p_age"] if g else np.nan,
                                  **_annual(m, age_col, min_sessions)))
            q = _bh([r["p"] for r in block])
            for r, qq in zip(block, q):
                r["p_fdr"] = qq
            for cat, regs in (("Association", LG.ASSOCIATION_REGIONS),
                              ("Primary", LG.PRIMARY_REGIONS)):   # composites + LORO
                for out in [None] + list(regs):
                    keep = [r for r in regs if r != out]
                    m, _ = _long_rows(morph, mod, cat, keep, age_df, age_col, hemisphere,
                                      min_sessions)
                    block.append(dict(species=sp, modality=mod,
                                      kind="composite" if out is None else "loro",
                                      category=cat, estimate=cat if out is None else f"{cat} - {out}",
                                      left_out=out or "", regions=", ".join(keep),
                                      single_region_left=len(keep) == 1,
                                      **_annual(m, age_col, min_sessions)))
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
            out.append(f"  {r['estimate']}: no data")
            continue
        extra = (f", P_FDR={r['p_fdr']:.3g}" if r["kind"] == "region" else
                 " (single region left)" if r.get("single_region_left") else "")
        out.append(f"  {r['estimate']}: {r['pct_per_year']:+.3f} %/year "
                   f"[{r['ci_low']:+.3f}, {r['ci_high']:+.3f}], P={r['p']:.3g}{extra}; "
                   f"n={r['n_subjects']} subjects, {r['n_sessions']} sessions")
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
                   f"[{min(p):+.3f}, {max(p):+.3f}] %/year: {sep}")
    for cat in ("Association", "Primary"):
        allr = next((r for r in block if r["kind"] == "composite" and r["category"] == cat), None)
        lo = [r for r in block if r["kind"] == "loro" and r["category"] == cat
              and np.isfinite(r["pct_per_year"])]
        if allr is None or not lo or not np.isfinite(allr["pct_per_year"]):
            continue
        sign = np.sign(allr["pct_per_year"])
        flips = [r["left_out"] for r in lo if np.sign(r["pct_per_year"]) != sign]
        lost = [r["left_out"] for r in lo if allr["p"] < 0.05 <= r["p"]]
        out.append(f"  LORO {cat}: range [{min(r['pct_per_year'] for r in lo):+.3f}, "
                   f"{max(r['pct_per_year'] for r in lo):+.3f}] vs all {allr['pct_per_year']:+.3f}"
                   f" %/year; sign change when leaving out: {', '.join(flips) or 'none'}; "
                   f"P >= 0.05 when leaving out: {', '.join(lost) or 'none'}")
    return out


def stats_text(st):
    return ("Supplementary Figure 10 -- region-by-region trajectories of the Figure 3 "
            "composites (left hemisphere)\n"
            "Annual change = 100*(exp(b)-1), b = slope/year of log(value) ~ age + (1|subject), ML; "
            "95% CI from b +/- 1.96 SE; single regions: BH-FDR over the 6 regions of one "
            "species x measure. LORO = composite without one region (association has 2 "
            "regions: leaving one out = the other single region).\n\n"
            + "\n".join(st["lines"]) + "\n")


def draw_regions(st, out_path, panel_w=3.4, panel_h=2.7):
    """4 x 3 figure (rows: Human curves, Macaque curves, Human annual change, Macaque
    annual change; columns = Figure 3 measures). Writes <out_path>.png/.pdf."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.lines import Line2D
    from Tools.ednix_bids_tools import PAPER_RC

    mods, species, tab = st["modalities"], st["species"], st["table"]
    nrows, ncols = 2 * len(species), len(mods)
    rc = dict(PAPER_RC); rc.update(LG._FIG_RC)
    with plt.rc_context(rc):
        fig = plt.figure(figsize=(panel_w * ncols + 0.6, panel_h * nrows * 0.95))
        gs = gridspec.GridSpec(nrows, ncols, figure=fig, hspace=0.45, wspace=0.55,
                               left=0.13, right=0.99, top=0.96, bottom=0.05)
        letter, empty, row_label = LG._panel_tools(fig)
        regions = LG.ASSOCIATION_REGIONS + LG.PRIMARY_REGIONS
        for i, sp in enumerate(species):                           # trajectories
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
                    row_label(ax, f"{st['labels'][sp]}\nper region")
                if i == 0 and c == ncols - 1:
                    h = [Line2D([], [], color=COLORS[r],
                                ls="-" if CATEGORY_OF[r] == "Association" else "--",
                                label=LABELS.get(r, r)) for r in regions]
                    ax.legend(handles=h, fontsize=5.5, frameon=False, loc="best")
        for i, sp in enumerate(species):                           # annual change
            for c, (mod, _, _, title) in enumerate(mods):
                ax = fig.add_subplot(gs[len(species) + i, c])
                letter(ax)
                b = tab[(tab["species"] == sp) & (tab["modality"] == mod)] if len(tab) else tab
                b = b[np.isfinite(b["pct_per_year"])] if len(b) else b
                if b.empty:
                    empty(ax, "no data")
                    continue
                order = ([r for r in regions if r in set(b.loc[b.kind == "region", "estimate"])] +
                         ["Association"] + [f"Association - {r}" for r in LG.ASSOCIATION_REGIONS] +
                         ["Primary"] + [f"Primary - {r}" for r in LG.PRIMARY_REGIONS])
                b = b.set_index("estimate").reindex([o for o in order if o in set(b["estimate"])])
                y = np.arange(len(b))[::-1]
                for yy, (est, r) in zip(y, b.iterrows()):
                    cat = r["category"]
                    col = (COLORS[est] if r["kind"] == "region"
                           else LG.FIG3D_COLORS[cat])
                    mk = {"region": "o", "composite": "s", "loro": "D"}[r["kind"]]
                    face = col if r["kind"] != "loro" else "white"
                    ax.plot([r["ci_low"], r["ci_high"]], [yy, yy], color=col, lw=1)
                    ax.plot(r["pct_per_year"], yy, mk, color=col, mfc=face, ms=3.5, mew=0.9)
                labels = [LABELS.get(e, e) if b.loc[e, "kind"] == "region" else
                          (e if b.loc[e, "kind"] == "composite" else
                           f"  - {LABELS.get(b.loc[e, 'left_out'], b.loc[e, 'left_out'])}")
                          for e in b.index]
                ax.set_yticks(y)
                ax.set_yticklabels(labels, fontsize=5.5)
                ax.axvline(0, color="0.5", lw=0.6, ls=":")
                ax.set_xlabel("Annual change (%/year, 95% CI)")
                ax.spines[["top", "right"]].set_visible(False)
                if c == 0:
                    ax.annotate(f"{st['labels'][sp]}\nannual change", xy=(-0.62, 0.5),
                                xycoords="axes fraction", rotation=90, ha="center",
                                va="center", fontsize=9, fontweight="bold")
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
