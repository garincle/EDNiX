"""
Figure 4 of the paper, drawn at print size (183 mm wide). Clément, 2026-10-27:

  a  the common profile: left, one large radar of the cross-species consensus (median
     across species of the species medians, every dataset; grey band = interquartile range
     across species) with the region names on its spokes, no title; right, one small radar per dataset, its group fingerprint (species colour)
     over the consensus it was compared with (black dashed: leave-one-species-out
     consensus, the median of the OTHER species, so that a dataset is never compared
     with a profile it contributed to), ranked by that r. Same spokes, same radial range
     everywhere.
  b  run-to-run profile correlation, same species / other dataset vs different species,
     all runs and selected runs (as before; the legend is above the panels, it used to
     sit on the n.s. bracket)
  c  dataset x dataset similarity, r (no longer 1 - r) of the group fingerprints, dataset
     names above the matrix, so that its lower edge is level with the x axis of b.

2026-10-30: b is replaced by b-d (fig4_validation_and_similarity, selection_validation.py:
held-out test, per dataset, pair classes) and the matrix becomes e; fig4bc is kept for the
former layout.
b and c are one figure (fig4bc) so that their axes are aligned exactly. The numbers are
those of paper_stats.figure4 (same vectors, same r).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import Exemples.Study_EDNiX.analysis.fmri_fingerprints.profile_stats as fa
import Exemples.Study_EDNiX.analysis.fmri_fingerprints.fingerprints as ff
from Exemples.Study_EDNiX.analysis.fmri_fingerprints import supp_panels as fsp
from Exemples.Study_EDNiX.analysis import paper_style as ps

MM = 1.0 / 25.4
LETTER_PT = ps.PANEL_LABEL_SIZE


def _radar_curve(ax, theta, v, order, **kw):
    t, a = ff._close(theta, v[order])
    return ax.plot(t, a, **kw)


def _common_with_band(vecs, species_of_label):
    """Common profile (median across species of the species medians), z-scored, and the
    interquartile range across species of the species medians, in the same units (the
    affine z-transform of the common profile applied to the quartiles)."""
    spm = ff.species_medians(vecs, species_of_label)
    mat = np.vstack(list(spm.values()))
    with np.errstate(invalid="ignore"):
        c = np.nanmedian(mat, axis=0)
        q25, q75 = np.nanpercentile(mat, [25, 75], axis=0)
    ok = np.isfinite(c)
    mu, sd = np.nanmean(c[ok]), np.nanstd(c[ok])
    sd = sd if sd > 0 else 1.0
    return (c - mu) / sd, (q25 - mu) / sd, (q75 - mu) / sd, len(spm)


def _grid(n, width_mm, big_mm, gap_mm=3.0, aspect=1.12):
    """Columns of the small-radar grid whose height is closest to the large radar's."""
    best = None
    for nc in (3, 4, 5, 6):
        cell = (width_mm - big_mm - gap_mm) / nc
        nr = int(np.ceil(n / nc))
        d = abs(nr * cell * aspect - big_mm)
        if best is None or d < best[0]:
            best = (d, nc, nr, cell)
    return best[1:]


def fig4a_common_profile(vecs, dist_df, species_of_label, common, order, region_mask, out,
                         width_mm=183, big_mm=80):
    """Panel a. vecs: {label: group fingerprint}; dist_df: ff.dataset_distances (column
    r_consensus = r with the leave-one-species-out consensus).
    Left: common profile (black) with the interquartile range across species (grey band),
    region names on the spokes. Right: one small radar per dataset, ranked by r; the
    number of columns is the one whose grid height matches the large radar; the legend
    fills the empty cells of the last row (or goes under the grid)."""
    rows = dist_df.dropna(subset=["r_consensus"]).sort_values("r_consensus", ascending=False)
    if rows.empty:
        return None
    order = np.asarray(order)
    labels = [ff.short_region(common[i]) for i in order]
    regions = [common[i] for i in order]
    mask = np.ones(len(common), bool) if region_mask is None else np.asarray(region_mask, bool)
    excluded = [not mask[i] for i in order]
    glob_c, lo, hi, n_sp = _common_with_band(vecs, species_of_label)
    pairs = []
    for r in rows.itertuples():
        ref, _ = ff._reference_for(vecs, species_of_label, "consensus", r.label, r.species, None)
        pairs.append((r, fa._zscore_values(vecs[r.label]), fa._zscore_values(ref)))
    allv = np.concatenate([x[np.isfinite(x)] for x in (glob_c, lo, hi)] +
                          [np.r_[v[np.isfinite(v)], c[np.isfinite(c)]] for _, v, c in pairs])
    rlim = (min(allv.min(), 0) - 0.3, allv.max() + 0.3)

    n = len(pairs)
    ncols, nrows, cell = _grid(n, width_mm, big_mm)
    step_y = cell * 1.12
    H_mm = max(big_mm, nrows * step_y)
    W, H = width_mm * MM, H_mm * MM
    with ff._print_size():
        fig = plt.figure(figsize=(W, H))
        # large radar, square, centred vertically, with room for the spoke names
        lab_mm = 15.0
        r_mm = big_mm - 2 * lab_mm
        y0 = (H_mm - big_mm) / 2 + lab_mm
        bax = fig.add_axes([lab_mm / width_mm, y0 / H_mm, r_mm / width_mm, r_mm / H_mm],
                           projection="polar")
        theta = ff._spider_axes(bax, labels, rlim, excluded, show_labels=True, regions=regions)
        t, a = ff._close(theta, lo[order])
        _, b = ff._close(theta, hi[order])
        bax.fill_between(t, a, b, color="0.55", alpha=0.35, lw=0, zorder=2)
        _radar_curve(bax, theta, glob_c, order, color="black", lw=1.0, zorder=3)
        bax.tick_params(axis="y", labelsize=4)
        # small radars, best r first
        x0 = big_mm + 3.0
        top = (H_mm - nrows * step_y) / 2
        for k, (r, v, c) in enumerate(pairs):
            i, j = divmod(k, ncols)
            side = cell * 0.80
            ax = fig.add_axes([(x0 + j * cell + (cell - side) / 2) / width_mm,
                               1 - (top + (i + 1) * step_y - 0.8) / H_mm,
                               side / width_mm, side / H_mm], projection="polar")
            th = ff._spider_axes(ax, labels, rlim, excluded, show_labels=False, regions=regions)
            ax.set_yticklabels([])
            t, a = ff._close(th, v[order])
            _, b = ff._close(th, c[order])
            col = ff._color(r.species)
            ax.fill_between(t, a, b, where=a > b, interpolate=True, color=col, alpha=0.25, lw=0)
            ax.plot(t, b, color="black", ls=(0, (2.5, 1.5)), lw=0.6)
            ax.plot(t, a, color=col, lw=0.8)
            ax.set_title(f"{ff.pretty_label(r.label)}\nr = {r.r_consensus:.2f}", fontsize=5,
                         color=ff._tcolor(r.species), pad=1.5)
        handles = [Line2D([], [], color="black", lw=1.0, label="common profile"),
                   Patch(facecolor="0.55", alpha=0.35, lw=0,
                         label=f"IQR across species (n = {n_sp})"),
                   Line2D([], [], color="black", lw=0.6, ls=(0, (2.5, 1.5)),
                          label="leave-one-species-out consensus"),
                   Line2D([], [], color="0.4", lw=0.8, label="dataset (species colour)")]
        free = nrows * ncols - n
        if free >= 1:
            # in the empty cells of the last row
            j0 = n - (nrows - 1) * ncols
            fig.legend(handles=handles, loc="center left", ncol=1, fontsize=5, frameon=False,
                       bbox_to_anchor=((x0 + j0 * cell + 1.0) / width_mm,
                                       1 - (top + (nrows - 0.5) * step_y) / H_mm))
        else:
            fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=5, frameon=False,
                       bbox_to_anchor=((x0 + ncols * cell / 2) / width_mm, 0.0))
        fa.save_paper(fig, out)
    return out


def _matrix_r(ax, cax, D, species_of_label, annotate=True):
    """Dataset x dataset r = 1 - D, names above the matrix and on the left."""
    labels = fsp._labels_sorted(list(D.index), species_of_label)
    M = 1.0 - D.loc[labels, labels].to_numpy(float)
    np.fill_diagonal(M, np.nan)
    n = len(labels)
    lo, hi = np.nanmin(M), np.nanmax(M)
    cmap = plt.get_cmap("viridis").with_extremes(bad="white")
    im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=lo, vmax=hi, interpolation="none")
    if annotate:
        for i in range(n):
            for j in range(n):
                if i != j and np.isfinite(M[i, j]):
                    ax.text(j, i, f"{M[i, j]:.2f}".replace("0.", ".", 1), ha="center",
                            va="center", fontsize=3.6,
                            color="white" if M[i, j] < lo + 0.45 * (hi - lo) else "black")
    sp = [fsp._sp(l, species_of_label) for l in labels]
    for i in range(1, n):
        if sp[i] != sp[i - 1]:
            ax.axhline(i - 0.5, color="white", lw=0.8)
            ax.axvline(i - 0.5, color="white", lw=0.8)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.xaxis.tick_top()
    ax.set_xticklabels([fsp._lab(l) for l in labels], rotation=55, ha="left",
                       rotation_mode="anchor", fontsize=4.8)
    ax.set_yticklabels([fsp._lab(l) for l in labels], fontsize=4.8)
    for t, s in zip(list(ax.get_xticklabels()) + list(ax.get_yticklabels()), sp + sp):
        t.set_color(fsp._tcol(s))
    ax.tick_params(length=0, pad=1)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    cb = ax.figure.colorbar(im, cax=cax)
    cb.set_label("r", fontsize=5.5)
    cb.ax.tick_params(labelsize=5, length=1.5)
    cb.outline.set_linewidth(0.3)


def fig4bc(pairs_by, stats, names, test, intra, D, species_of_label, out, width_mm=183,
           height_mm=68):
    """Panels b (two histograms) and c (similarity matrix) in one figure, the lower edge of
    the matrix level with the x axes of b. Letter 'c' drawn here; 'b' is added by the
    layout (paper_style.compose_images)."""
    with fa.paper_style():
        fig = plt.figure(figsize=(width_mm * MM, height_mm * MM))
        bottom, top_b = 13.0, height_mm - 9.0                      # mm
        hb = top_b - bottom
        # matrix: square, same lower edge as b, its height = what is left under the names
        names_mm = 17.0
        # height 62 -> 68 mm (Clément 2026-10-29): matrix side 35 -> 41 mm
        side = min(hb + 9.0 - names_mm + 3.0, 46.0)
        mx = width_mm - side - 9.0
        axes_b = []
        left_b, gap, right_b = 10.0, 9.0, mx - 30.0
        wb = (right_b - left_b - gap) / 2
        for k in range(2):
            axes_b.append(fig.add_axes([(left_b + k * (wb + gap)) / width_mm, bottom / height_mm,
                                        wb / width_mm, hb / height_mm]))
        for k, (ax, name) in enumerate(zip(axes_b, names)):
            if stats.get(name) is None:
                ax.axis("off")
                continue
            fa._draw_intra_inter_paper(ax, pairs_by[name], stats[name], test,
                                       "all runs" if k == 0 else "selected", legend=False,
                                       intra=intra)
            if k:
                ax.set_ylabel("")
        h, l = axes_b[0].get_legend_handles_labels()
        fig.legend(h, l, loc="upper left", ncol=2, frameon=False, fontsize=5.5,
                   bbox_to_anchor=(left_b / width_mm, 1.0))
        if D is not None:
            ax_c = fig.add_axes([mx / width_mm, bottom / height_mm, side / width_mm,
                                 side / height_mm])
            cax = fig.add_axes([(mx + side + 1.5) / width_mm, bottom / height_mm,
                                1.2 / width_mm, side / height_mm])
            _matrix_r(ax_c, cax, D, species_of_label)
            fig.canvas.draw()
            r = fig.canvas.get_renderer()
            ylab = ax_c.get_yticklabels()
            x_c = min(t.get_window_extent(r).x0 for t in ylab) / fig.bbox.width if ylab else mx / width_mm
            # top of the figure = top of the b+c image once cropped: level with letter b
            fig.text(x_c - 1.5 / width_mm, 1.0, "c",
                     fontsize=LETTER_PT, fontweight="bold", va="top", ha="right")
        fa.save_paper(fig, out)
    return out


def fig4_validation_and_similarity(val, D, species_of_label, out, width_mm=183, height_mm=68):
    """Panels b-d (what the selection buys, selection_validation.py) and e (similarity
    matrix, former c) in one figure; lower edges aligned. Clément, 2026-10-30: the
    histograms of the former b are replaced by b (held-out test), c (per dataset) and d
    (pair classes). Letter 'b' is added by the layout; c, d, e are drawn here."""
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import selection_validation as sv
    fs = 5.5
    with fa.paper_style():
        fig = plt.figure(figsize=(width_mm * MM, height_mm * MM))
        bottom, top = 13.0, height_mm - 9.0
        hb = top - bottom
        W, H = width_mm, height_mm
        # height 68 mm and matrix side 41 mm as fig4bc (Clément 2026-10-29)
        ax_sch = fig.add_axes([9.0 / W, (top + 1.5) / H, 30.0 / W, 3.5 / H])
        ax_b = fig.add_axes([9.0 / W, bottom / H, 30.0 / W, (hb - 2.0) / H])
        ax_c = fig.add_axes([51.0 / W, bottom / H, 17.0 / W, (hb - 2.0) / H])
        ax_d = fig.add_axes([80.0 / W, bottom / H, 32.0 / W, (hb - 2.0) / H])
        h = val["heldout"]
        sv.draw_b_schematic(ax_sch, n=h["n_values"], fs=fs)
        sv.draw_b(ax_b, h, fs=fs)
        sv.draw_c(ax_c, h, ps.species_color, ps.dataset_marker, fs=fs)
        sv.draw_d(ax_d, val["levels"], fs=fs)
        fig.canvas.draw()
        rend = fig.canvas.get_renderer()
        for ax, letter in ((ax_c, "c"), (ax_d, "d")):
            x0 = ax.get_tightbbox(rend).x0 / fig.bbox.width
            fig.text(x0 - 0.5 / W, 1.0, letter, fontsize=LETTER_PT, fontweight="bold",
                     va="top", ha="right")
        if D is not None:
            names_mm = 17.0
            side = min(hb + 9.0 - names_mm + 3.0, 46.0)
            mx = W - side - 9.0
            ax_e = fig.add_axes([mx / W, bottom / H, side / W, side / H])
            cax = fig.add_axes([(mx + side + 1.5) / W, bottom / H, 1.2 / W, side / H])
            _matrix_r(ax_e, cax, D, species_of_label)
            fig.canvas.draw()
            ylab = ax_e.get_yticklabels()
            x_e = (min(t.get_window_extent(rend).x0 for t in ylab) / fig.bbox.width
                   if ylab else mx / W)
            fig.text(x_e - 1.5 / W, 1.0, "e", fontsize=LETTER_PT, fontweight="bold",
                     va="top", ha="right")
        fa.save_paper(fig, out)
    return out
