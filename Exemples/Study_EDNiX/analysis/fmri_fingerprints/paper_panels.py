"""
fig5_paper.py
=============
Panels of the paper's Figures 5 and 6, one file per panel (PNG 300 dpi + PDF).

Each panel is drawn by the SAME function as its supplementary figure, with
the same size, fonts, line widths and resolution, called with paper=True:
no suptitle, no panel title, no statistics box -- only the data, axis
labels, legends, n and p. Same function also means same numbers: the p on a
paper panel is the p of its supplementary twin. Everything else (question,
test, full statistics, reading) is in the .txt written next to each panel.
At the end of a run, every other figure of 06_figures/ is moved to
06_figures/supplementary/.

Figure 5 -- what the selection does
    5a  surface maps, all runs | selected, one colour range per dataset
    5b  run-to-run similarity, all runs | selected            (= F2b)
    5e  selection x species, dataset pairs, main vs control run (launcher, last step)
    5c  connectivity fingerprints, all runs | selected        (= F1_fingerprint)
    5d  distance to the LOSO consensus, all runs -> selected   (= F17 a)

Figure 6 -- what the comparison shows
    6a  distance to the LOSO consensus vs divergence | distance to Human vs
        divergence | distance to the LOSO consensus vs brain mass  (= F6ab + c)
    6b  each dataset's fingerprint over its LOSO consensus      (= F15a)
    6c  arousal state within the datasets holding both states   (= F18)
    6d  fragmentation: posterior parietal vs somatosensory coverage and peak,
        vs brain mass (no supplementary twin; same style)

Brain mass: read from the Herculano-Houzel workbook (title rows, 'mean±SD'
strings and thousand separators handled; values printed), or CONFIG
BRAIN_MASS_G.
"""

import os
import math
import re
import glob
import shutil
import itertools

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats as st

import Exemples.Study_EDNiX.analysis.fmri_fingerprints.profile_stats as fa
import Exemples.Study_EDNiX.analysis.fmri_fingerprints.fingerprints as ff
from Statistics.Group_fMRI.run_selection.surfaces import (  # moved verbatim, 2026-10-27
    _oversampled,
    map_paths_by_dataset,
    opj,
    render_common_scale,
    surface_maps_for_5a)

opb = os.path.basename

DPI = 300                     # same as every supplementary figure
MAIN_DIRS = ("Fig5", "Fig6", "Paper", "supplementary")
BINOMIAL_ALIASES = {"Dog": ["canis lupus familiaris", "canis familiaris", "canis lupus"]}


# ═══════════════════════════════════════════════════════════════════════════
# helpers
# ═══════════════════════════════════════════════════════════════════════════

def _save(fig, path_noext, svg=False):
    """png + pdf (text kept as text); svg=True adds an .svg in which every image
    and every text is a separate object, to be rearranged in Inkscape."""
    exts = ("png", "pdf") + (("svg",) if svg else ())
    with plt.rc_context({"pdf.fonttype": 42, "svg.fonttype": "none",
                         "image.composite_image": False}):
        for ext in exts:
            fig.savefig(f"{path_noext}.{ext}", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  [saved] {path_noext}." + "/.".join(exts))


def _p(p):
    return f"p = {ff._fmt_p(p)}" if p is not None and np.isfinite(p) else "p = n/a"


def _clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _legend_txt(path_noext, question, what, test, result, reading, legend=None):
    """The .txt of a panel: LEGEND (the figure legend, manuscript style; also
    gathered in 06_figures/Figure_legends.txt), then question, what is plotted,
    test, full result and how to read it."""
    with open(f"{path_noext}.txt", "w", encoding="utf-8") as f:
        if legend:
            f.write(f"LEGEND\n------\n{legend}\n\n")
        f.write(f"QUESTION\n--------\n{question}\n\nWHAT IS PLOTTED\n---------------\n{what}\n\n"
                f"TEST\n----\n{test}\n\nRESULT\n------\n{result}\n\n"
                f"HOW TO READ IT\n--------------\n{reading}\n")


def _key(k):
    return tuple(str(x) for x in k)


def _run_number(run):
    d = re.findall(r"\d+", str(run))
    return int(d[-1]) if d else 1


# ═══════════════════════════════════════════════════════════════════════════
# Brain mass (Herculano-Houzel workbook)
# ═══════════════════════════════════════════════════════════════════════════

_SP_COL = re.compile(r"species|binomial|scientific|taxon", re.I)
_BM_EXCLUDE = re.compile(r"cort|cereb|neuron|densit|ratio|log|%|body|other|rest|olfact|bulb",
                         re.I)


def _is_brain_mass_col(c):
    c = str(c)
    return bool(re.search(r"brain", c, re.I)
                and re.search(r"mass|weight|\(g\)|_g\b|\bg\b", c, re.I)
                and not _BM_EXCLUDE.search(c))


def _norm_name(s):
    """Binomial key: lower case, letters only, common name in brackets removed
    ('Canis familiaris (Domestic dog)' -> 'canisfamiliaris')."""
    return re.sub(r"[^a-z]", "", re.sub(r"\(.*?\)", "", str(s)).lower())


def _parse_number(v):
    """'1,508.910±299.140' -> 1508.91; '4,600' -> 4600; '0,42' -> 0.42."""
    if isinstance(v, (int, float, np.integer, np.floating)):
        return float(v) if np.isfinite(v) else np.nan
    s = str(v).strip()
    if not s or s.lower() in ("nan", "-", "n/a", "na"):
        return np.nan
    s = re.split(r"±|\+/-|\+-", s)[0].replace(" ", "").replace("\xa0", "").replace(" ", "")
    if "," in s and "." in s:
        s = s.replace(",", "")
    elif re.fullmatch(r"-?\d{1,3}(,\d{3})+", s):
        s = s.replace(",", "")
    else:
        s = s.replace(",", ".")
    m = re.match(r"-?\d+(\.\d+)?([eE]-?\d+)?", s)
    return float(m.group(0)) if m else np.nan


def _brain_mass_tables(path, max_header_row=8):
    """Every sheet holding a species column and a whole-brain mass column,
    whatever row the header is on (the HH workbook has a title row first).
    Returns [(sheet, species column name, mass column name, DataFrame(name, g))]."""
    out = []
    for sheet, raw in pd.read_excel(path, sheet_name=None, header=None).items():
        for i in range(min(max_header_row, len(raw))):
            cells = [str(c) for c in raw.iloc[i].tolist()]
            sp = [j for j, c in enumerate(cells) if _SP_COL.search(c)]
            bm = [j for j, c in enumerate(cells) if _is_brain_mass_col(c)]
            if sp and bm:
                col = cells[bm[0]]
                scale = 1e-3 if re.search(r"\bmg\b|\(mg\)", col, re.I) else (
                    1e3 if re.search(r"\bkg\b|\(kg\)", col, re.I) else 1.0)
                body = raw.iloc[i + 1:]
                t = pd.DataFrame(dict(name=body.iloc[:, sp[0]].astype(str),
                                      g=body.iloc[:, bm[0]].map(_parse_number) * scale))
                out.append((sheet, cells[sp[0]], col, t.dropna(subset=["g"])))
                break
    return out


def load_brain_mass(path, species, species_to_tip, override=None):
    """{species: whole-brain mass in g}. override (CONFIG BRAIN_MASS_G) wins.
    Otherwise every sheet with a species column and a whole-brain mass column
    is read (header row found automatically); several rows of one species ->
    median. Prints, for every species, the value and the sheet it came from."""
    if override:
        print(f"  [brain mass] CONFIG BRAIN_MASS_G used: {override}")
        return {sp: float(override[sp]) for sp in species if sp in override}
    if not path or not os.path.isfile(path):
        print(f"  [brain mass] [WARN] {path!r} not found -- brain-mass panels skipped")
        return {}
    tables = _brain_mass_tables(path)
    if not tables:
        listing = {s: [str(c) for c in d.iloc[:3].values.ravel()]
                   for s, d in pd.read_excel(path, sheet_name=None, header=None).items()}
        print(f"  [brain mass] [WARN] no sheet with a species column and a whole-brain mass "
              f"column in {path}; first cells: {listing} -- set BRAIN_MASS_G in CONFIG")
        return {}
    out, src = {}, {}
    for sp in species:
        tip = species_to_tip.get(sp, sp).replace("_", " ")
        keys = {_norm_name(c) for c in [tip] + BINOMIAL_ALIASES.get(sp, [])}
        gs = " ".join(tip.split()[:2])
        vals, where = [], []
        for sheet, _, col, t in tables:
            names = t["name"].map(_norm_name)
            hit = names.isin(keys)
            if not hit.any():
                hit = names.str.startswith(_norm_name(gs)) & (names.str.len() > 0)
            if hit.any():
                vals += list(t.loc[hit, "g"])
                where.append(f"{sheet!r}/{col!r}")
        if vals:
            out[sp] = float(np.median(vals))
            src[sp] = ", ".join(where)
    print("  [brain mass] " + "; ".join(f"{sp} = {m:.4g} g ({src[sp]})" for sp, m in out.items()))
    missing = sorted(set(species) - set(out))
    if missing:
        print(f"  [brain mass] [WARN] no row for {missing} -- add them to BRAIN_MASS_G")
    if "Human" in out and not 800 < out["Human"] < 2000:
        print(f"  [brain mass] [WARN] Human = {out['Human']:.4g} g: check the unit")
    return out


# ═══════════════════════════════════════════════════════════════════════════
# Figure 5
# ═══════════════════════════════════════════════════════════════════════════

def _surface_png(map_path, keep):
    if not map_path or not isinstance(map_path, str):
        return None
    pngs = sorted(glob.glob(opj(os.path.dirname(map_path), "*.png")))
    for k in keep:
        hits = [p for p in pngs if k in opb(p)]
        if hits:
            return hits[0]
    return None


def _content_bbox(img, tol=0.03):
    """(r0, r1, c0, c1) of everything that differs from the background colour
    (top-left pixel), or None for a uniform image."""
    a = img[..., :3] if img.ndim == 3 else img
    diff = np.abs(a - a[0, 0]).max(axis=-1) if a.ndim == 3 else np.abs(a - a[0, 0])
    rows = np.where(diff.max(axis=1) > tol)[0]
    cols = np.where(diff.max(axis=0) > tol)[0]
    if rows.size and cols.size:
        return rows[0], rows[-1] + 1, cols[0], cols[-1] + 1
    return None


def _autocrop(img, tol=0.03):
    """Remove the uniform border (background colour of the top-left pixel)."""
    b = _content_bbox(img, tol)
    return img[b[0]:b[1], b[2]:b[3]] if b else img


def _valid_png(path, min_bytes=2000):
    """A render that exists, is not a stub and is not a uniform image."""
    if not path or not os.path.isfile(path) or os.path.getsize(path) < min_bytes:
        return False
    try:
        return _content_bbox(plt.imread(path)) is not None
    except Exception:
        return False


def _aligned_crops(images):
    """Crop the images of ONE dataset with the SAME box (union of their content
    boxes). Renders of one dataset come from the same scene and camera, so the
    brains sit on the same pixels: a common box keeps them aligned, whereas
    cropping each image on its own shifts them. Images of different sizes are
    first padded (background colour, centred) to the largest size."""
    keys = [k for k, im in images.items() if im is not None]
    if not keys:
        return {}
    H = max(images[k].shape[0] for k in keys)
    W = max(images[k].shape[1] for k in keys)
    padded = {}
    for k in keys:
        im = images[k]
        if im.shape[:2] != (H, W):
            canvas = np.empty((H, W) + im.shape[2:], dtype=im.dtype)
            canvas[...] = im[0, 0]
            r0, c0 = (H - im.shape[0]) // 2, (W - im.shape[1]) // 2
            canvas[r0:r0 + im.shape[0], c0:c0 + im.shape[1]] = im
            im = canvas
        padded[k] = im
    boxes = [b for b in (_content_bbox(im) for im in padded.values()) if b]
    if not boxes:
        return padded
    r0, r1 = min(b[0] for b in boxes), max(b[1] for b in boxes)
    c0, c1 = min(b[2] for b in boxes), max(b[3] for b in boxes)
    return {k: im[r0:r1, c0:c1] for k, im in padded.items()}


def _compress_gaps(images, max_gap_frac=0.04, tol=0.03):
    """Shrink the blank bands BETWEEN the views of a render (lateral, medial,
    dorsal) to at most max_gap_frac of the image height. The bands are found
    on the union of the images of ONE dataset and the same rows / columns are
    removed from all of them, so they stay aligned. Scenes whose views are far
    apart (Dog, Macaque) otherwise end up with small brains once fitted in a
    cell."""
    keys = [k for k, im in images.items() if im is not None]
    if not keys:
        return images
    ref = images[keys[0]]
    if any(images[k].shape[:2] != ref.shape[:2] for k in keys):
        return images
    blank = np.ones(ref.shape[:2], bool)
    for k in keys:
        a = images[k][..., :3] if images[k].ndim == 3 else images[k]
        diff = (np.abs(a - a[0, 0]).max(axis=-1) if a.ndim == 3 else np.abs(a - a[0, 0]))
        blank &= diff <= tol
    gap = max(2, int(round(max_gap_frac * ref.shape[0])))

    def keep_index(is_blank):
        keep = np.ones(is_blank.size, bool)
        i = 0
        while i < is_blank.size:
            if is_blank[i]:
                j = i
                while j < is_blank.size and is_blank[j]:
                    j += 1
                if j - i > gap:
                    keep[i + gap // 2: j - (gap - gap // 2)] = False
                i = j
            else:
                i += 1
        return keep

    rows = keep_index(blank.all(axis=1))
    cols = keep_index(blank.all(axis=0))
    return {k: (im[rows][:, cols] if im is not None else None) for k, im in images.items()}


pretty_label = ff.pretty_label


def _surface_montage(out, labels, species_of_label, pics, groups, group_names,
                     per_row=2, mirror=False, name_color=None, zoom=None):
    """`per_row` datasets per row; for each dataset its groups side by side,
    all cropped with one common box so the brains are aligned across groups,
    and the blank space between the views shrunk (same rows / columns for
    every group of a dataset). Layout in inches (no overlap): one header line
    with the group names (bold) at the top, then for every row a line with the
    dataset name and the images; a thin vertical line separates the groups of
    a dataset.
    mirror=True: the left dataset of a row shows its groups in reverse order, so
    with groups (kept, excluded) a row reads  excluded | kept || kept | excluded.
    name_color: one colour for every dataset name (None = species colour).
    zoom: {species: relative size} -- the images of a species are drawn at
    zoom / max(zoom) of their cell, centred (1 = fills the cell).
    Saved as png, pdf and svg (each image a separate object for Inkscape)."""
    imgs = {}
    for lab in labels:
        raw = {g: (plt.imread(pics[(lab, g)]) if (lab, g) in pics and
                   os.path.isfile(pics[(lab, g)]) else None) for g in groups}
        for g, im in _compress_gaps(_aligned_crops(raw)).items():
            if im is not None:
                imgs[(lab, g)] = im
    if not imgs:
        return False
    zoom = {k: float(v) for k, v in (zoom or {}).items() if float(v) > 0}
    zmax = max(zoom.values()) if zoom else 1.0
    aspect = float(np.median([im.shape[1] / im.shape[0] for im in imgs.values()]))
    # print size: two Nature columns (183 mm), text 6-7 pt
    n_pairs, W = max(1, int(per_row)), fa.PAPER_WIDTH["double"]
    gap_pair, gap_img = 0.20, 0.07
    head_h, name_h, row_gap = 0.17, 0.15, 0.07
    nrows = int(np.ceil(len(labels) / n_pairs))
    ng = len(groups)
    cell_w = (W - gap_pair * (n_pairs - 1) - gap_img * (ng - 1) * n_pairs) / (ng * n_pairs)
    cell_h = cell_w / aspect
    H = head_h + nrows * (name_h + cell_h + row_gap)

    def x_of(pair, j):
        return pair * (ng * cell_w + gap_img * (ng - 1) + gap_pair) + j * (cell_w + gap_img)

    def order_of(pair):
        idx = list(range(ng))
        return idx[::-1] if (mirror and pair % 2 == 0 and n_pairs > 1) else idx

    fig = plt.figure(figsize=(W, H))
    for pair in range(min(n_pairs, len(labels))):
        for j, k in enumerate(order_of(pair)):
            fig.text((x_of(pair, j) + cell_w / 2) / W, 1 - 0.5 * head_h / H, group_names[k],
                     ha="center", va="center", fontsize=7, fontweight="bold", color="#222222",
                     family=fa.paper_font())
        for j in range(ng - 1):
            xs = (x_of(pair, j) + cell_w + gap_img / 2) / W
            fig.add_artist(Line2D([xs, xs], [row_gap / H, 1 - head_h / H], color="#BBBBBB",
                                  lw=0.5, transform=fig.transFigure))
    for i, lab in enumerate(labels):
        row, pair = divmod(i, n_pairs)
        y_top = H - head_h - row * (name_h + cell_h + row_gap)
        fig.text(x_of(pair, 0) / W, (y_top - 0.5 * name_h) / H, pretty_label(lab),
                 ha="left", va="center", fontsize=6.5, fontweight="bold",
                 color=name_color or ff._tcolor(species_of_label.get(lab)),
                 family=fa.paper_font())
        z = zoom.get(species_of_label.get(lab), zmax) / zmax
        for j, g in enumerate([groups[k] for k in order_of(pair)]):
            w_, h_ = cell_w * z, cell_h * z
            x0 = x_of(pair, j) + (cell_w - w_) / 2
            y0 = y_top - name_h - cell_h + (cell_h - h_) / 2
            ax = fig.add_axes([x0 / W, y0 / H, w_ / W, h_ / H])
            ax.axis("off")
            im = imgs.get((lab, g))
            if im is not None:
                # same box for every group of a dataset: 'equal' keeps pixels
                # square, so the brains stay aligned
                ax.imshow(im, interpolation="lanczos", aspect="equal")
            else:
                ax.text(0.5, 0.5, "no map", ha="center", va="center", fontsize=6,
                        color="#999999", transform=ax.transAxes)
    fa.save_paper(fig, out, svg=True)
    return True


def fig5a_surfaces(out, cfg, maps, pics, ranges, labels, species_of_label, keep,
                   supp_out=None):
    """Surface maps of the group seed maps, rendered in part 1 (render_common_scale,
    one colour range per dataset shared by its three maps: all runs, kept by the
    final round, same number not kept), drawn twice:
      5a  (out)      : selected | excluded, two datasets per row, mirrored
                       (excluded | selected || selected | excluded);
      S5a (supp_out) : all runs | selected | excluded, two datasets per row.
    A map without a common-scale render falls back to its own per-map render.
    Brains aligned within a dataset (one crop box), dataset names in black."""
    if not maps:
        print("  [5a] no group map -- skipped")
        return None
    pics = dict(pics or {})
    per_map, used = [], []
    for lab, gs in maps.items():
        for g, p in gs.items():
            if (lab, g) in pics and not _valid_png(pics[(lab, g)]):
                print(f"  [5a] [WARN] {lab}/{g}: common-scale render missing or empty "
                      f"({pics[(lab, g)]}) -- per-map render used instead")
                pics.pop((lab, g))
            if (lab, g) not in pics:
                fallback = _surface_png(p, keep)
                if fallback and _valid_png(fallback):
                    pics[(lab, g)] = fallback
                    per_map.append(f"{lab}/{g}")
            if (lab, g) in pics:
                used.append(f"{lab}/{g}: {pics[(lab, g)]}")
    per_row = cfg.get("SURFACE_DATASETS_PER_ROW", 2)
    pct = cfg.get("SURFACE_COMMON_PCT", 90.0)
    scale = ("None (illustration). Colour range: one fixed symmetric range per dataset, "
             "shared by its three maps (all runs, kept, not kept; the same images are used "
             f"in 5a and S5a) = +/- the {pct:g}th percentile of |Fisher z| over the non-zero "
             "voxels of those maps (palette RY-BC-BL).")
    res = ("Ranges: " + ", ".join(f"{k} +/-{v:.3g}" for k, v in (ranges or {}).items())
           + (f". PER-MAP SCALE (common-scale render failed or disabled) for: "
              f"{', '.join(per_map)}" if per_map else "")
           + "\nImages used:\n" + "\n".join(used))
    read = ("Within a dataset the colours are directly comparable between its maps; across "
            "datasets they are not (each dataset has its own range).")
    surf = ("Group seed map (mean Fisher z) of every dataset on its own inflated surface "
            "(lateral, medial and dorsal views as rendered by the EDNiX scene). ")
    done = {}
    zoom = cfg.get("SURFACE_ZOOM") or {}
    legend5a = ("5a | Retrosplenial seed maps of the selected runs and of the runs the "
                "selection excluded. Group mean Fisher z of each dataset on its inflated "
                "cortical surface (lateral, medial and dorsal views); centre of each row, runs "
                "kept by the final selection round; outer columns, the same number of runs "
                "drawn from those not kept. One symmetric colour range per dataset, shared by "
                f"its maps (± the {pct:g}th percentile of |z|); colours are not comparable "
                "across datasets.")
    ok = _surface_montage(out, labels, species_of_label, pics, ("kept_final", "not_kept_final"),
                          ("selected", "excluded"), per_row=per_row, mirror=True,
                          name_color="black", zoom=zoom)
    if ok:
        _legend_txt(out, "(illustration) What do the retrosplenial seed maps of the selected "
                         "runs look like next to the runs the selection left out?",
                    surf + "Runs kept by the final selection round ('selected', centre of each "
                    "row) and the same number of runs drawn from those the selection left out "
                    "('excluded', outer columns). "
                    f"{per_row} datasets per row, mirrored (excluded | selected || selected "
                    "| excluded), a thin line between the two groups. The two images of a "
                    "dataset share one crop box and the blank space between the views is "
                    "shrunk identically in both, so the brains are aligned; image size per "
                    f"species: SURFACE_ZOOM = {zoom or 'all 1'}. Printed width 183 mm. An .svg "
                    "is written next to the png/pdf: each image and each text is a separate "
                    "object (Inkscape). The version with all runs is supplementary S5a.",
                    scale, res, read, legend=legend5a)
        done["5a"] = out
    else:
        print("  [5a] no kept / not-kept surface png available -- skipped")
    if supp_out:
        ok = _surface_montage(supp_out, labels, species_of_label, pics,
                              ("all", "kept_final", "not_kept_final"),
                              ("all runs", "selected", "excluded"), per_row=per_row,
                              mirror=False, name_color="black", zoom=zoom)
        if ok:
            _legend_txt(supp_out, "(illustration) What do the retrosplenial seed maps look "
                                  "like before and after the selection, and what did the "
                                  "selection leave out?",
                        surf + "Every run; the runs kept by the final selection round; the "
                        "same number of runs drawn from those not kept. The three images of a "
                        "dataset share one crop box, so the brains are aligned. An .svg is "
                        "written next to the png/pdf (Inkscape).", scale, res, read,
                        legend=legend5a.replace("5a |", "S5a |", 1).replace(
                            "centre of each row, runs kept by the final selection round; outer "
                            "columns, the same number of runs drawn from those not kept",
                            "left to right, every run, the runs kept by the final selection "
                            "round, the same number of runs drawn from those not kept"))
            done["S5a"] = supp_out
        else:
            print("  [S5a] no surface png available -- skipped")
    if not done:
        return None
    return dict(ranges=ranges, per_map=per_map, written=done)


def _value_txt(value_col, norm):
    """Plain-language name of a profile column for the legends."""
    m = re.fullmatch(r"q(\d+)q(\d+)_r", str(value_col))
    t = re.fullmatch(r"pct_tail(\d+)", str(value_col))
    if m:
        txt = (f"upper (P{m.group(2)}) and lower (P{m.group(1)}) percentiles of each region's "
               f"voxels, each channel standardised across regions")
    elif t:
        txt = (f"% of each region's voxels in the top (+) and bottom (\u2212) {t.group(1)} % of "
               f"the map, each channel standardised across regions")
    else:
        txt = "median r of each region's voxels"
    n = fa._norm_label(norm).strip()
    return txt + (f", {n.replace('RMS-scaled', 'divided by its root mean square')}"
                  if n else "")


def fig5b_within_across(out, runs_all, runs_sel, common, sol, norm, value_col,
                        test="wilcoxon", unit="run", intra="other_dataset"):
    """= F2b (fa.plot_within_vs_across_all_vs_selected) in paper mode: same
    pairs and statistics, the two compared classes, the stars of `test` on the
    figure; both tests in the .txt. intra='other_dataset' (default): same
    species in ANOTHER dataset vs different species, same-dataset pairs left
    out; 'all': every same-species pair (same dataset included) vs different
    species."""
    _, stats = fa.plot_within_vs_across_all_vs_selected(
        runs_all, runs_sel, common, sol, value_col=value_col, zscore=norm,
        output_path=out + ".png", paper=True, test=test, intra=intra)
    stats = {("all runs" if k == "All runs" else "selected"): v for k, v in stats.items()}
    shown = "Wilcoxon rank-sum" if test == "wilcoxon" else "dataset-level permutation"
    cross = intra == "other_dataset"
    in_name = "same species, other dataset" if cross else "intra-species"
    out_name = "different species" if cross else "inter-species"
    def _res_line(k, v):
        if cross:
            pairs_txt = (f"{in_name} median r = {v['median_within']:.3f} (n = {v['n_within']} "
                         f"pairs), {out_name} median r = {v['median_across']:.3f} (n = "
                         f"{v['n_across']} pairs); not compared: same dataset median r = "
                         f"{v['median_same_dataset']:.3f} (n = {v['n_same_dataset']} pairs)")
        else:
            pairs_txt = (f"{in_name} median r = {v['median_within']:.3f} (n = {v['n_within']} "
                         f"pairs: same dataset {v['median_same_dataset']:.3f}, n = "
                         f"{v['n_same_dataset']}; same species other dataset "
                         f"{v['median_same_species_other_dataset']:.3f}, n = "
                         f"{v['n_same_species_other_dataset']}), {out_name} median r = "
                         f"{v['median_across']:.3f} (n = {v['n_across']} pairs)")
        p_fig = v["p_wilcoxon"] if test == "wilcoxon" else v["p_perm_dataset_level"]
        return (f"{k}: {pairs_txt}; AUC = {v['auc']:.3f}; Wilcoxon rank-sum one-sided p = "
                f"{ff._fmt_p(v['p_wilcoxon'])}; dataset-level permutation p = "
                f"{ff._fmt_p(v['p_perm_dataset_level'])} ({v['n_perm']} permutations); on the "
                f"figure ({shown}): {fa.stars(p_fig)}")

    res = "\n".join(_res_line(k, v) for k, v in stats.items() if v)
    n_units = {k: sum(len(p) for p in pool.values()) for k, pool in
               (("all runs", runs_all), ("selected", runs_sel))}
    if cross:
        classes = (f"pairs of {unit}s from the same species but two different datasets (blue; "
                   f"available for the species with two datasets) and pairs from different "
                   f"species (orange); pairs from the same dataset, which also share protocol, "
                   f"scanner and usually the animal, are not shown")
    else:
        classes = ("intra-species pairs (blue; same species, same or another dataset) and "
                   "inter-species pairs (orange)")
    legend = ("5b | Run-to-run similarity of the regional profiles, within and across species. "
              f"Distribution of the Pearson correlation between the profiles of two {unit}s for "
              f"{classes}; dashed lines, medians. Left, every {unit} (n = "
              f"{n_units['all runs']}); right, the {unit}s kept by the final selection round "
              f"(n = {n_units['selected']}). Brackets: one-sided "
              + ("Wilcoxon rank-sum test on the pairs, not corrected for multiple comparisons"
                 if test == "wilcoxon" else
                 "permutation test, species labels permuted across datasets")
              + "; *** P < 0.001, ** P < 0.01, * P < 0.05, n.s. not significant.")
    reading = (
        "Only 'same species, other dataset' pairs separate species from dataset: two runs "
        "of one dataset share the protocol, the scanner and usually the animal, so their "
        "similarity says nothing about species. "
        + ("This is why same-dataset pairs are left out here (WITHIN_ACROSS_PAIRS="
           "'other_dataset'); they are shown, with the two compared classes, on F2b. "
           if cross else
           "Here they are pooled with the other intra-species pairs (WITHIN_ACROSS_PAIRS="
           "'all') and, being the majority, they dominate the intra-species distribution and "
           "the Wilcoxon test; the permutation test is barely affected because same-dataset "
           "pairs stay intra-species under every permutation. ")
        + "After selection, two datasets of one species are scored against the same "
          "leave-one-species-out consensus, so part of their gain in similarity is "
          "mechanical; the control run with one common target (LEAVE_ONE_OUT=False, "
          "*_control_common_target) removes that asymmetry.")
    _legend_txt(out, "Are the profiles of two runs of the same species more similar than those "
                     "of two runs of different species (primary hypothesis)?",
                "Distribution of the correlation between the regional profiles of two "
                f"{unit}s: {classes}; dashed = medians. Left: every unit; right: units kept by "
                "the final selection round. Same pairs, bins and statistics as F2b. Printed "
                "width 120 mm.",
                "Primary test, on the figure: " + (
                    "Wilcoxon rank-sum (Mann-Whitney U), one-sided (intra > inter), on the "
                    "pairs, not corrected for multiple comparisons (one hypothesis). The pairs "
                    "share units (each unit enters many pairs), so they are not independent "
                    "and this p is optimistic (Mantel 1967 Cancer Res 27:209; Nili et al. 2014 "
                    "PLoS Comput Biol 10:e1003553). Reported with it: the dataset-level "
                    "permutation test (statistic = median(intra) - median(inter), same pair "
                    "classes; species labels permuted across DATASETS, each dataset keeps its "
                    "units, pair correlations fixed; 2000 permutations, one-sided), which "
                    "respects that dependence." if test == "wilcoxon" else
                    "dataset-level permutation (statistic = median(intra) - median(inter), same "
                    "pair classes; species labels permuted across DATASETS, each dataset keeps "
                    "its units; 2000 permutations, one-sided). Reported with it: the Wilcoxon "
                    "rank-sum on the pairs, optimistic because pairs share units."),
                res, reading, legend=legend)
    return stats


def _region_perm_test(spm, n_perm=10000, seed=0):
    names = sorted(spm)
    M = np.vstack([spm[s] for s in names])
    ok = np.isfinite(M).all(axis=0)
    M = M[:, ok]
    iu = np.triu_indices(len(names), 1)
    obs = np.corrcoef(M)[iu].mean()
    rng = np.random.default_rng(seed)
    null = np.array([np.corrcoef(np.vstack([rng.permutation(m) for m in M]))[iu].mean()
                     for _ in range(n_perm)])
    return dict(r=float(obs), p=(1 + int(np.sum(null >= obs))) / (n_perm + 1),
                n_species=len(names), n_regions=int(ok.sum()))


def fig5c_fingerprint(out, f1_all, f1_sel, common, sol, value_col, norm, region_mask,
                      divergence, final_round):
    """= F1_fingerprint_species (ff.plot_fingerprints) in paper mode, with
    'n = .. species, p = ..' under each panel header."""
    tests = {}
    for k, prof in enumerate((f1_all, f1_sel or {})):
        vv = ff.group_vectors(prof, common, value_col, norm)
        if vv:
            tests[k] = _region_perm_test(ff.species_medians(vv, sol))
    if not tests:
        return None
    notes = {k: f"mean r between species = {t['r']:.2f} {fa.stars(t['p'])}"
             for k, t in tests.items()}
    ff.plot_fingerprints(f1_all, f1_sel, common, sol, value_col, norm, "species", region_mask,
                         out + ".png", divergence, f"Selected (final round = {final_round})",
                         "", paper=True, notes=notes)
    names = {0: "all runs", 1: "selected"}
    _legend_txt(out, "Is the retrosplenial profile shared across species?",
                "Connectivity fingerprint (Passingham et al. 2002 Nat Rev Neurosci 3:606; Mars "
                "et al. 2018 eLife 7:e35237): one line per species (median of its datasets), "
                "black = consensus (median across species). Left: every run (input of F1); "
                f"right: final-round selection (round {final_round}, input of F1b). "
                f"{fa._norm_label(norm)}{value_col}; regions grouped by lobe as in Figure 1 "
                "(counter-clockwise from the top; rim = lobe colour); dotted ring = 0; grey label with * = region not used "
                "by the selection score (invariance filter, 02_selection/region_invariance.xlsx).",
                "Mean Pearson r between all pairs of species profiles, against 10000 "
                "permutations of the region labels within each species (one-sided).",
                "; ".join(f"{names[k]}: mean r between species = {v['r']:.3f}, p = "
                          f"{ff._fmt_p(v['p'])}, n = {v['n_species']} species, "
                          f"{v['n_regions']} regions" for k, v in tests.items()),
                "A shared shape means the species lines follow the black line. The region-label "
                "permutation ignores spatial autocorrelation between neighbouring regions and is "
                "therefore liberal: read the size of r, not only p. The right panel is tighter by "
                "construction (runs were selected for resembling a cross-species consensus); "
                "the regions marked * were not used for that selection and are the fair "
                "comparison between the two panels.",
                legend=("5c | Connectivity fingerprints of the retrosplenial seed. Polar plot of "
                        "the regional profile of each species (median of its datasets; "
                        f"{_value_txt(value_col, norm)}) over the common "
                        "regions, grouped by lobe as in Figure 1 (counter-clockwise from the top; "
                        "rim, lobe colour); "
                        "black, consensus (median across species); dotted ring, 0; * region not "
                        "used by the selection score. Left, every run; right, the runs kept by "
                        f"the final selection round (round {final_round}). Above each panel, "
                        "mean Pearson r between all pairs of species profiles and its "
                        "significance against 10,000 permutations of the region labels within "
                        "each species (one-sided; *** P < 0.001, ** P < 0.01, * P < 0.05, n.s. "
                        f"not significant); n = {next(iter(tests.values()))['n_species']} "
                        "species."))
    return tests


def fig5d_rank(out, fingerprints, states, final_round, ref="Human"):
    """= F17 panel a (ff.plot_before_after) in paper mode."""
    da, ds = fingerprints.get("ba_dist_all"), fingerprints.get("ba_dist_sel")
    if da is None or ds is None or len(da) < 3:
        return None
    res, m = ff.plot_before_after(da, ds, states, out + ".png",
                                  fingerprints.get("ba_source", ""), final_round, ref,
                                  paper=True)
    r = res["r_consensus"]
    _legend_txt(out, "Does the selection reorder the datasets?",
                "Each line = one dataset (colour = species; marker fill = arousal state: filled "
                "awake, open anaesthetised, half both): r between its group profile and the "
                "leave-one-species-out consensus, every run (left) vs runs kept by the final "
                f"round {final_round} (right). Each consensus is rebuilt from the same pool (all "
                f"/ selected) of the other species. Source: {fingerprints.get('ba_source', '')}.",
                "Spearman rho between the two columns across datasets; p from 10000 "
                "permutations of one column (two-sided).",
                f"rho = {r['spearman_rho']:.3f}, p = {ff._fmt_p(r['p_perm'])}, n = {r['n']} "
                f"datasets; median change of r {r['median_change']:+.3f}.",
                "A higher r after selection is expected by construction and is not a result. A "
                "rank close to 1 means the selection did not reorder the datasets: the "
                "between-dataset differences of Figure 6 are not created by the selection.",
                legend=("5d | The selection does not reorder the datasets. Pearson r between "
                        "the group profile of each dataset (line; colour, species; marker fill, "
                        "arousal state: filled awake, open anaesthetised, half both, diamond "
                        "unknown) and the leave-one-species-out consensus, every run (left) and "
                        f"runs kept by the final selection round (right; round {final_round}). "
                        "Top right, Spearman rho between the two columns (two-sided, 10,000 "
                        f"permutations; *** P < 0.001, ** P < 0.01, * P < 0.05); n = {r['n']} "
                        "datasets."))
    return r


# ═══════════════════════════════════════════════════════════════════════════
# Figure 6
# ═══════════════════════════════════════════════════════════════════════════

def fig6a_distances(out, dist_df, states, brain_mass, source, ref="Human"):
    """= F6ab (ff.plot_f6_two_references) in paper mode + brain-mass panel."""
    tests, sp = ff.plot_f6_two_references(dist_df, states, out + ".png", source, ref,
                                          paper=True, brain_mass=brain_mass)
    names = dict(consensus_with_ref="r with the LOSO consensus vs divergence, all species",
                 consensus_without_ref=f"same without {ref}",
                 reference=f"r with {ref} vs divergence (species other than {ref})",
                 brain_mass_with_ref="r with the LOSO consensus vs log10 brain mass, all species",
                 brain_mass_without_ref=f"same without {ref}")
    res = "\n".join(f"{names.get(k, k)}: correlation = {v['r']:+.3f}, p = "
                    f"{ff._fmt_p(v['p_two'])} (n = {v['n']} species, "
                    f"{'exact, ' + str(v['n_perm']) + ' permutations' if v['exact'] else 'MC'})"
                    for k, v in tests.items())
    band = dist_df.loc[dist_df["species"] == ref, "r_reference"].values
    if brain_mass:
        res += "\nbrain mass used (g): " + ", ".join(f"{k} {v:.4g}" for k, v in
                                                     sorted(brain_mass.items()))
    _legend_txt(out, "Left: do species share the pattern of the other species according to "
                     "their divergence from Human? Middle: are species closer to Human in time "
                     "also more similar to Human? Right: does the shared part follow brain size?",
                "Small markers = datasets (colour = species; fill = arousal state: filled awake, "
                "open anaesthetised, half both, diamond unknown), labelled with the dataset name; "
                "large rings = species means, the points that are fitted. Left and right: r "
                "between a dataset's group profile and the median of the OTHER species "
                "(leave-one-species-out consensus). Middle: r with the median of the "
                f"{ref} datasets; grey band = r between the two {ref} datasets "
                f"({', '.join(f'{v:.2f}' for v in band)}), the ceiling for any other species. "
                f"Solid line: fit on all species; dashed: without {ref}. Brain mass: "
                "Herculano-Houzel et al. 2015 Brain Behav Evol 86:145 (Table 5) and "
                "Jardim-Messeder et al. 2017 Front Neuroanat 11:118 (dog), whole brain without "
                f"olfactory bulbs. Source of the profiles: {source}.",
                "Pearson correlation on species means; p exact over all n! assignments of the x "
                "values to species (7! = 5040, 6! = 720), two-sided. n on the panels = species. "
                "Same p as with 1 - r (linear transform); the sign of the slope is reversed.",
                res,
                "The consensus is a median across the sampled species, not an ancestral state. "
                f"A gradient carried by {ref} alone (significant with, not without {ref}) says "
                f"that {ref} is the outlier, not that the shared part changes with time. "
                "Divergence from Human and brain mass are collinear in this sample; the Dog "
                "(large brain, distant from Human) is the species that separates them. Arousal "
                "state is also collinear with both (Figure 6c measures its effect directly).",
                legend=("6a | Similarity of each dataset to the other species and to Human. "
                        "Small markers, datasets (colour, species; fill, arousal state: filled "
                        "awake, open anaesthetised, half both, diamond unknown); rings, species "
                        "means, the fitted points. Left, Pearson r between a dataset's group "
                        "profile and the leave-one-species-out consensus against divergence "
                        f"time from {ref}; middle, r with the median of the {ref} datasets "
                        f"(grey band, r between the two {ref} datasets: the ceiling); right, r "
                        "with the leave-one-species-out consensus against log10 brain mass "
                        "(Herculano-Houzel et al. 2015; Jardim-Messeder et al. 2017). Solid "
                        f"line, fit over all species; dashed, without {ref}. Top right, Pearson "
                        "r over species means, exact two-sided permutation over all "
                        "assignments of the x values to species (*** P < 0.001, ** P < 0.01, * "
                        "P < 0.05, n.s. not significant); n = "
                        f"{tests['consensus_with_ref']['n']} species."))
    return tests


def fig6b_fingerprint_grid(out, fingerprints, sol, common, region_mask, divergence, ref="Human"):
    """= F15a (ff.plot_distance_fingerprint_grid, which='consensus') in paper mode."""
    vecs, dist_df = fingerprints.get("vecs_all"), fingerprints.get("dist_all")
    order = fingerprints.get("region_order")
    if not vecs or dist_df is None or len(dist_df) < 3:
        return None
    ff.plot_distance_fingerprint_grid(vecs, dist_df, sol, "consensus", common, order,
                                      region_mask, out + ".png", ref, divergence,
                                      fingerprints.get("group_all_source", ""), paper=True)
    labels = [ff.short_region(common[i]) for i in order]
    _legend_txt(out, "Which regions carry each dataset's departure from the shared profile?",
                "One fingerprint per dataset (colour) over its leave-one-species-out consensus "
                "(black dashed), most similar first; r under each name = its value on Figure 6a. "
                "Both curves z-scored over the common regions, so the mean squared gap between "
                "them equals 2 x (1 - r). Shading marks the regions where the dataset EXCEEDS "
                "its reference; where it falls below, the gap is left unshaded. "
                "Region order (counter-clockwise from the top, by lobe) as in Figure 5c: "
                + ", ".join(labels) + "; * = not used by the selection score.",
                "Descriptive (no test).",
                "See 03_group_maps/fingerprint_tables.xlsx (distances_all) for the values.",
                "Shaded lobes = regions more coupled with the seed than in the other species; "
                "a dashed line outside the coloured one = regions less coupled. Pearson r "
                "ignores a common offset or gain.",
                legend=("6b | Where each dataset departs from the shared profile. Fingerprint of "
                        "each dataset (colour) over its leave-one-species-out consensus (black "
                        "dashed), both z-scored over the common regions (so the mean squared gap "
                        "between them is 2 x (1 - r)), ordered from the most to the least "
                        "similar; r, Pearson r of Figure 6a; shading, regions where the dataset "
                        "exceeds its reference. Regions in the order of Figure 5c, where they are "
                        "labelled (counter-clockwise from the top; rim, lobe colour)."))
    return True


# ---- 6c arousal state -------------------------------------------------------

def fig6c_f8b(paths, dist_df_b, bids_configs, final_round, n_perm=10000, paper=True,
              criterion_note=None):
    """F8b by arousal state (supplementary F8b) and its paper version (6c).
    Every run of every dataset: Pearson r with the final round's selection
    criterion of its dataset (all runs, kept included). Mouse and
    Marmoset_Hori: awake vs anaesthetised, tested within animal.
    criterion_note: a sentence added to the .txt and to the legend when the
    criterion was rebuilt in other regions than the selection used (seed
    region kept for the selection, removed from the analyses)."""
    d = ff.annotate_states(dist_df_b, bids_configs)
    tests = ff.state_r_tests(d, "pearson_r", n_perm=n_perm)
    ff.plot_f8b_by_state(d, tests, opj(paths["figures"],
                                       "F8b_run_distance_to_selection_criterion.png"),
                         title=f"Per-run correlation to the round-{final_round} selection "
                               f"criterion, by arousal state")
    out6 = opj(paths["figures"], "Fig6", "6c_arousal_state_F8b")
    if paper:
        os.makedirs(opj(paths["figures"], "Fig6"), exist_ok=True)
        ff.plot_f8b_by_state(d, tests, out6 + ".png", paper=True)
    pooled = None
    mixed = d[d["label"].isin(list(tests))].copy()
    if not mixed.empty:
        mixed["animal"] = mixed["label"].astype(str) + ":" + mixed["animal"].astype(str)
        mixed["subject"] = mixed["animal"]
        pooled = ff.state_contrast(mixed, "pearson_r", n_perm=n_perm)
    lines = []
    for lab, t in tests.items():
        g = d[d["label"] == lab]
        med = g.groupby("state_class")["pearson_r"].median()
        lines.append(f"{lab}: median r awake {med.get(ff.AWAKE, np.nan):.3f}, anaesthetised "
                     f"{med.get(ff.ANESTH, np.nan):.3f}; mean over animals of (anaesthetised - "
                     f"awake) = {t['effect']:+.4f}, p = {ff._fmt_p(t['p_subject'])} "
                     f"({t['p_subject_method']}); run-level p = {ff._fmt_p(t['p_run'])} "
                     f"(liberal)")
        for agent in sorted(set(g.loc[g["state_class"] == ff.ANESTH, "state"])):
            sub = g[(g["state_class"] == ff.AWAKE) | (g["state"] == agent)].copy()
            sub["subject"] = sub["animal"]
            ta = ff.state_contrast(sub, "pearson_r", n_perm=n_perm)
            if ta and g["state"].nunique() > 2:
                lines.append(f"   {agent} vs awake: {ta['effect']:+.4f}, p = "
                             f"{ff._fmt_p(ta['p_subject'])} ({ta['n_subjects']} animals)")
    if pooled:
        lines.append(f"POOLED ({', '.join(tests)}): {pooled['effect']:+.4f}, p = "
                     f"{ff._fmt_p(pooled['p_subject'])} ({pooled['p_subject_method']})")
    d.to_excel(opj(paths["group_maps"], "run_r_to_selection_criterion_by_state.xlsx"),
               index=False)
    txt = dict(question="Are the runs of anaesthetised animals more (or less) similar to "
                        "their dataset's selection criterion than the awake runs of the "
                        "same animals?",
               what="Every run of every dataset: Pearson r between its regional profile and "
                    f"the profile the final round ({final_round}) ranked that dataset's runs "
                    "against (as F8b), all runs, kept ones included (crosses). Circles = awake "
                    "runs, triangles = anaesthetised runs (all anaesthetics). Datasets holding "
                    "both states in the same animals (Mouse, Marmoset_Hori) are split into two "
                    "sub-columns, awake left and anaesthetised right; violin and quartile ticks "
                    "per sub-column; bracket = within-animal test (stars), the only statistics "
                    "on the figure." + (" " + criterion_note if criterion_note else ""),
               test="Within animal: effect = mean over animals of (mean r anaesthetised - mean "
                    "r awake); p = exact sign-flip of the per-animal differences (2^n "
                    "assignments; with 4 animals the smallest possible p is 0.125). The run-level "
                    "p (state shuffled among the runs of each animal) treats runs as independent "
                    "and is liberal. Pooled: the same sign-flip over the animals of both "
                    "datasets.",
               reading="The criterion is what decides selection, so a lower r under "
                       "anaesthesia means anaesthetised runs are less likely to be kept. The "
                       "criterion of the final round is built from the runs kept at the "
                       "previous round: if those were mostly awake, the comparison favours "
                       "awake runs by construction. In Marmoset_Hori the state is also the "
                       "session (anaesthetised = ses-01).")
    n_an = {lab: t["n_subjects"] for lab, t in tests.items()}
    legend = ("6c | Similarity of every run to its dataset's selection criterion, by arousal "
              "state. Pearson r between each run's regional profile and the profile the final "
              f"selection round (round {final_round}) ranked that dataset's runs against; "
              "circles, awake runs; triangles, anaesthetised runs; crosses, selected runs; "
              "violins and ticks, distribution and quartiles. Datasets holding both states in "
              "the same animals are split into awake (left) and anaesthetised (right) "
              "sub-columns; brackets, exact sign-flip test of the per-animal difference "
              "(anaesthetised − awake; " + "; ".join(
                  f"{pretty_label(k)}, n = {v} animals" for k, v in n_an.items())
              + "); * P < 0.05, ** P < 0.01, *** P < 0.001, n.s. not significant."
              + (" " + criterion_note if criterion_note else ""))
    for path in [opj(paths["figures"], "F8b_run_distance_to_selection_criterion")] + (
            [out6] if paper else []):
        _legend_txt(path, txt["question"], txt["what"], txt["test"], "\n".join(lines),
                    txt["reading"], legend=legend.replace("6c |", "F8b |", 1)
                    if "F8b" in path else legend)
    return dict(per_dataset=tests, pooled=pooled)


def anesthesia_table(state_runs, rankings=None):
    """Run-level table for the pooled model of 6c, built from the F18 runs:
    anesth (0/1), animal_uid (dataset:animal), run_order (rank of the run number
    within subject and session), the two r of F18 and the round-1 score."""
    if state_runs is None or state_runs.empty:
        return pd.DataFrame()
    df = state_runs[state_runs["state_class"].isin([ff.AWAKE, ff.ANESTH])].copy()
    df["anesth"] = (df["state_class"] == ff.ANESTH).astype(int)
    df["animal_uid"] = df["label"].astype(str) + ":" + df["subject"].astype(str)
    df["run_order"] = df.groupby(["label", "bids_subject", "session"])["run"].transform(
        lambda s: s.map(_run_number).rank(method="first"))
    sc = []
    for lab in df["label"].unique():
        rk = (rankings or {}).get(lab)
        if rk is None or rk.empty:
            continue
        sc.append(pd.DataFrame(dict(
            label=lab, key=[(str(a), str(b), str(c), str(d)) for a, b, c, d in
                            zip(rk["subject"], rk["session"], rk["run"], rk["seed"])],
            score=rk["pearson_r"].values)))
    if sc:
        df = df.merge(pd.concat(sc), on=["label", "key"], how="left")
    return df


def _hier(df, value, extra=()):
    """Hierarchical means: run -> animal (mean of its runs) -> dataset (mean
    of its animals) -> species (mean of its datasets), so an animal with many
    runs, or a dataset with many animals, does not weigh more. df needs the
    columns species, label, animal, `value` and `extra`. Returns (per dataset
    DataFrame [species, label, *extra, value], per species Series indexed by
    [species, *extra])."""
    keys = list(extra)
    a = df.groupby(["species", "label", "animal"] + keys)[value].mean().reset_index()
    d = a.groupby(["species", "label"] + keys)[value].mean().reset_index()
    s = d.groupby(["species"] + keys)[value].mean()
    return d, s


def _anatomical_regions(common):
    """Anatomical names of `common` (two-channel entries 'x|pos' -> 'x')."""
    out = []
    for c in common:
        base, sep, ch = str(c).rpartition("|")
        name = base if sep and re.fullmatch(r"pos|neg|q\d{1,2}", ch) else str(c)
        if name not in out:
            out.append(name)
    return out


def group_contrast_test(values, species, group, alternative="two-sided"):
    """Exact test of mean(group) - mean(others) over species: null = every
    assignment of len(group) species to the group (C(n, k); 21 for 2 of 7,
    hence p >= 1/21 = 0.048 whatever the alternative). alternative:
    'two-sided' (|null| >= |obs|), 'greater' or 'less'. values: {species: value}."""
    sp = [s_ for s_ in species if np.isfinite(values.get(s_, np.nan))]
    k = sum(s_ in group for s_ in sp)
    if k == 0 or k == len(sp):
        return dict(diff=np.nan, p=np.nan, n_assign=0, n=len(sp), k=k, alternative=alternative)
    v = np.array([values[s_] for s_ in sp], float)

    def stat(idx):
        m = np.zeros(len(sp), bool)
        m[list(idx)] = True
        return v[m].mean() - v[~m].mean()

    obs = stat([i for i, s_ in enumerate(sp) if s_ in group])
    null = np.array([stat(c) for c in itertools.combinations(range(len(sp)), k)])
    eps = 1e-12 * max(1.0, float(np.max(np.abs(null))))
    if alternative == "greater":
        p = np.mean(null >= obs - eps)
    elif alternative == "less":
        p = np.mean(null <= obs + eps)
    else:
        p = np.mean(np.abs(null) >= abs(obs) - eps)
    return dict(diff=float(obs), p=float(p), n_assign=len(null), n=len(sp), k=k,
                alternative=alternative)


def _one_dataset_tests(d, value, species, group, alternative="two-sided"):
    """group_contrast_test with ONE dataset per species, for every
    combination (species with two datasets: Human, Macaque, Marmoset ->
    2 x 2 x 2 = 8). d: per dataset [species, label, value]. Asks whether the
    result rests on a single dataset."""
    choices = []
    for s_ in species:
        v = [float(x) for x in d.loc[d["species"] == s_, value] if np.isfinite(x)]
        choices.append(v or [np.nan])
    res = [group_contrast_test(dict(zip(species, c)), species, group, alternative)
           for c in itertools.product(*choices)]
    diffs = [r["diff"] for r in res if np.isfinite(r["diff"])]
    ps = [r["p"] for r in res if np.isfinite(r["p"])]
    if not diffs:
        return dict(n=0, diff_min=np.nan, diff_max=np.nan, p_min=np.nan, p_max=np.nan)
    return dict(n=len(res), diff_min=float(min(diffs)), diff_max=float(max(diffs)),
                p_min=float(min(ps)), p_max=float(max(ps)))


# Indices of the regional tables (sri.region_indices), with their unit, the
# label of the y axis and the sentence of the legend. `q` = TAIL_Q.
def region_index_specs(q):
    """Every regional index of the run table (sri.region_indices) that 6d can
    show: column -> file tag, axis label, legend text (+ chance level)."""
    lo, hi = int(round(100 * q)), int(round(100 * (1 - q)))
    return {
        "upper_percentile_r": dict(
            tag=f"P{hi}", file=f"upper_P{hi}", method="percentiles",
            ylabel=f"P{hi} of the region's voxels (r)",
            what=f"the {hi}th percentile of the region's voxel values in the run's seed map "
                 f"(r; its most coupled {lo} %)"),
        "lower_percentile_r": dict(
            tag=f"P{lo}", file=f"lower_P{lo}", method="percentiles",
            ylabel=f"P{lo} of the region's voxels (r)",
            what=f"the {lo}th percentile of the region's voxel values in the run's seed map "
                 f"(r; its least coupled {lo} %)"),
        "median_r": dict(
            tag="median r", file="median_r", method="median_r",
            ylabel="median r of the region's voxels",
            what="the median of the region's voxel values in the run's seed map (r)"),
        "top_tail_pct": dict(
            tag=f"top {lo} %", file=f"top{lo}_fraction", method="fraction_posneg",
            ylabel=f"% of voxels in the map's top {lo} %",
            what=f"the % of the region's voxels in the top {lo} % of the run's seed map "
                 f"(map percentiles over the voxels of the common regions; {lo} % without "
                 f"regional structure)", chance=float(lo)),
        "bottom_tail_pct": dict(
            tag=f"bottom {lo} %", file=f"bottom{lo}_fraction", method="fraction_posneg",
            ylabel=f"% of voxels in the map's bottom {lo} %",
            what=f"the % of the region's voxels in the bottom {lo} % of the run's seed map "
                 f"(map percentiles over the voxels of the common regions; {lo} % without "
                 f"regional structure)", chance=float(lo)),
        "positive_pct": dict(
            tag="r > 0", file="positive_fraction", method=None,
            ylabel="% of voxels with r > 0",
            what="the % of the region's voxels with a positive correlation with the seed "
                 "(r > 0) in the run's seed map; depends on the global level of correlation "
                 "of the map (no global signal regression)"),
        "negative_pct": dict(
            tag="r < 0", file="negative_fraction", method=None,
            ylabel="% of voxels with r < 0",
            what="the % of the region's voxels anticorrelated with the seed (r < 0) in the "
                 "run's seed map; depends on the global level of correlation of the map (no "
                 "global signal regression)"),
    }


def selected_run_keys(selections_final, members_by_label):
    """{label: set of run keys} behind the final-round selected units (with
    UNIT='subject' every run of each selected animal)."""
    out = {}
    for lab, sel in (selections_final or {}).items():
        mem = {_key(k): v for k, v in (members_by_label or {}).get(lab, {}).items()}
        out[lab] = {_key(r) for u in sel for r in mem.get(_key(u), [u])}
    return out


def region_index_table(run_index, keys_by_label, animal_of):
    """The rows of run_index (one per run x region) of the selected runs, with
    the animal of each run."""
    if run_index is None or run_index.empty:
        return pd.DataFrame()
    t = run_index.copy()
    key = list(zip(t["subject"].astype(str), t["session"].astype(str), t["run"].astype(str),
                   t["seed"].astype(str)))
    keep = np.array([k in keys_by_label.get(lab, ()) for lab, k in zip(t["label"], key)])
    t = t[keep].copy()
    if t.empty:
        return t
    kk = [k for k, ok in zip(key, keep) if ok]
    t["animal"] = [str(animal_of(lab, k)) for lab, k in zip(t["label"], kk)]
    t["run_id"] = ["_".join(k) for k in kk]
    return t.reset_index(drop=True)


def fig6d_region_values(out, tab, species_of_label, divergence, index_col, spec, regions,
                        group=("Human", "Macaque"), all_regions=None, source="",
                        panel="6d"):
    """6d / S6d: per region, the value of `index_col` in each SELECTED run (its
    own seed map); species value = mean over the runs of each animal, then over
    the animals of each dataset, then over the datasets of the species. One
    exact two-sided test per region, Human + Macaque vs the other species (the
    species is the unit: 21 assignments of 2 among 7, p >= 0.048), not
    corrected for the number of regions. Right panel: across every common
    region, the group difference of the index against the group difference in
    log10 relative size (region voxels / voxels of all common regions); dashed
    = fit through the regions not shown on the left, band = +/-2 SD: a tested
    region far from the band is not explained by its size."""
    if tab is None or tab.empty or index_col not in tab:
        print(f"  [{panel}] no {index_col} in the run table -- skipped")
        return None
    t = tab.dropna(subset=[index_col]).copy()
    t = t[t["species"].notna()]
    tot = t.groupby(["label", "run_id"])["n_voxels"].transform("sum")
    t["log_rel_size"] = np.log10(t["n_voxels"] / tot)
    ds_val, sp_val = _hier(t.rename(columns={index_col: "v"}), "v", ["region"])
    _, sp_size = _hier(t, "log_rel_size", ["region"])
    val, size = sp_val.unstack(), sp_size.unstack()
    species = sorted(val.index, key=lambda s_: divergence.get(s_, np.inf))
    val, size = val.loc[species], size.loc[species]
    group = [g for g in group if g in species]
    show = [r for r in regions if r in val.columns and val[r].notna().sum() >= 3]
    if not show or len(group) in (0, len(species)):
        print(f"  [{panel}] regions {regions} / group {group} not testable -- skipped")
        return None
    in_g = np.array([s_ in group for s_ in species])
    tests = {}
    for reg in show:
        tests[reg] = dict(
            group=group_contrast_test(val[reg].to_dict(), species, group),
            one_dataset=_one_dataset_tests(ds_val[ds_val["region"] == reg], "v", species,
                                           group),
            size=group_contrast_test(size[reg].to_dict(), species, group))
    present = [r for r in val.columns if val[r].notna().sum() >= 3]
    d_val = {r: float(val.loc[in_g, r].mean() - val.loc[~in_g, r].mean()) for r in present}
    d_size = {r: float(size.loc[in_g, r].mean() - size.loc[~in_g, r].mean()) for r in present}
    others = [r for r in present if r not in show]
    size_fit = None
    if len(others) >= 4:
        x_o = np.array([d_size[r] for r in others])
        y_o = np.array([d_val[r] for r in others])
        if np.ptp(x_o) > 0:
            b1, b0 = np.polyfit(x_o, y_o, 1)
            res_o = y_o - (b0 + b1 * x_o)
            sd = float(np.std(res_o, ddof=2))
            size_fit = dict(slope=float(b1), intercept=float(b0), sd=sd, n=len(others),
                            r=float(np.corrcoef(x_o, y_o)[0, 1]))
            for reg in show:
                resid = d_val[reg] - (b0 + b1 * d_size[reg])
                size_fit[reg] = dict(residual=float(resid),
                                     z=float(resid / sd) if sd else np.nan,
                                     n_regions_beyond=int(np.sum(np.abs(res_o) >= abs(resid))))
    n_runs = int(t[["label", "run_id"]].drop_duplicates().shape[0])
    n_animals = int(t[["label", "animal"]].drop_duplicates().shape[0])

    # ---- figure (print size) ------------------------------------------------
    order = [s_ for s_ in species if s_ in group] + [s_ for s_ in species if s_ not in group]
    k, n_o = len(group), len(species) - len(group)
    xo = list(np.linspace(-0.36, -0.24, k) if k > 1 else [-0.30])
    xo += list(np.linspace(-0.10, 0.36, n_o) if n_o > 1 else [0.13])
    xoff = dict(zip(order, xo))
    with fa.paper_style():
        fig, axes = plt.subplots(1, 2, figsize=(fa.PAPER_WIDTH["double"], 64 * fa.MM),
                                 gridspec_kw=dict(width_ratios=[1.75, 1]))
        ax = axes[0]
        rng = np.random.default_rng(0)
        ymin = float(np.nanmin(t.loc[t["region"].isin(show), index_col]))
        ymax = float(np.nanmax(t.loc[t["region"].isin(show), index_col]))
        span = (ymax - ymin) or 1.0
        for i, reg in enumerate(show):
            sub = t[t["region"] == reg]
            if i > 0:
                ax.axvline(i - 0.5, color="#E6E6E6", lw=0.5, zorder=0)
            for s_ in order:
                yv = val.loc[s_, reg]
                if not np.isfinite(yv):
                    continue
                x0 = i + xoff[s_]
                rv = sub.loc[sub["species"] == s_, index_col].values
                ax.scatter(x0 + rng.uniform(-0.022, 0.022, len(rv)), rv, s=2.2, lw=0,
                           alpha=0.45, color=ff._color(s_), zorder=2, rasterized=True)
                if s_ in group:
                    ax.scatter(x0, yv, s=16, zorder=4, lw=0.4, edgecolor="black",
                               facecolor=ff._color(s_))
                else:
                    ax.scatter(x0, yv, s=16, zorder=4, lw=0.9, edgecolor=ff._color(s_),
                               facecolor="white")
            means = []
            for members in (group, [s_ for s_ in order if s_ not in group]):
                xs = [i + xoff[s_] for s_ in members if np.isfinite(val.loc[s_, reg])]
                m = float(val.loc[members, reg].mean())
                ax.hlines(m, min(xs) - 0.04, max(xs) + 0.04, color="black", lw=1.0, zorder=5)
                means.append(np.mean(xs))
            top = float(np.nanmax(sub[index_col].values))
            fa.bracket(ax, means[0], means[1], top + 0.04 * span,
                       tests[reg]["group"]["p"], dy=0.025 * span)
        if spec.get("chance") is not None:
            ax.axhline(spec["chance"], color="#AAAAAA", lw=0.5, ls=(0, (2, 2)), zorder=0)
        if "percentile" in index_col or index_col == "median_r":
            ax.axhline(0, color="#AAAAAA", lw=0.5, ls=(0, (2, 2)), zorder=0)
        ax.set_xticks(range(len(show)))
        ax.set_xticklabels([ff.short_region(r) for r in show])
        ax.set_xlim(-0.5, len(show) - 0.5)
        ax.set_ylim(ymin - 0.04 * span, ymax + 0.16 * span)
        ax.set_ylabel(spec["ylabel"])
        ax.tick_params(axis="x", length=0)

        ax = axes[1]
        for r in present:
            tested = r in show
            ax.scatter(d_size[r], d_val[r], s=14 if tested else 7, lw=0, zorder=3,
                       color="black" if tested else "#B0B0B0")
        if size_fit:
            xx = np.linspace(min(d_size.values()), max(d_size.values()), 50)
            yy = size_fit["intercept"] + size_fit["slope"] * xx
            ax.fill_between(xx, yy - 2 * size_fit["sd"], yy + 2 * size_fit["sd"],
                            color="#EBEBEB", lw=0, zorder=0)
            ax.plot(xx, yy, color="#888888", lw=0.7, ls=(0, (3, 2)), zorder=1)
        ax.axhline(0, color="#CCCCCC", lw=0.5, zorder=0)
        ax.axvline(0, color="#CCCCCC", lw=0.5, zorder=0)
        G = " + ".join(group)
        ax.set_xlabel(f"Δ log$_{{10}}$ relative region size ({G} − others)")
        ax.set_ylabel(f"Δ {spec['tag']} ({G} − others)")
        h = [Line2D([], [], color="#888888", marker="o", mec="black", mew=0.4, ms=4, lw=0,
                    label=f"{G} (species mean)"),
             Line2D([], [], color="#888888", marker="o", mfc="white", mew=0.9, ms=4, lw=0,
                    label="other species (species mean)"),
             Line2D([], [], color="#888888", marker="o", ms=1.8, lw=0, alpha=0.6,
                    label="selected run")]
        h += [Line2D([], [], color=ff._color(s_), marker="o", ms=4, lw=0, label=s_)
              for s_ in species]
        h += [Line2D([], [], color="#888888", ls=(0, (3, 2)), lw=0.7,
                     label="size fit through the other regions ±2 SD")]
        fig.legend(handles=h, loc="lower center", ncol=6, frameon=False,
                   bbox_to_anchor=(0.5, -0.005))
        fig.tight_layout(rect=(0, 0.13, 1, 1), w_pad=2.0)
        # labels of the tested regions only (grey = the other regions, values in the
        # .xlsx), pushed apart vertically when they would overlap
        fa.repel_labels(axes[1], [(d_size[r], d_val[r], ff.short_region(r), "black")
                                  for r in show], fs=5)
        fa.save_paper(fig, out)

    # ---- tables and legend ----------------------------------------------------
    G = " + ".join(group)
    lines, trows = [], []
    for reg, tt in tests.items():
        g = tt["group"]
        od = tt["one_dataset"]
        line = (f"{ff.short_region(reg)}: {G} - others = {g['diff']:+.4f}, p = "
                f"{ff._fmt_p(g['p'])} ({fa.stars(g['p'])}; exact two-sided, "
                f"{g['n_assign']} assignments)")
        if od["n"] > 1:
            pr = (f"p = {ff._fmt_p(od['p_max'])}" if od["p_min"] == od["p_max"] else
                  f"p {ff._fmt_p(od['p_min'])} to {ff._fmt_p(od['p_max'])}")
            line += (f"; one dataset per species ({od['n']} combinations): {od['diff_min']:+.4f} "
                     f"to {od['diff_max']:+.4f}, {pr}")
        sz = tt["size"]
        line += (f"; log10 relative size {G} - others = {sz['diff']:+.3f} "
                 f"(p = {ff._fmt_p(sz['p'])})")
        if size_fit and reg in size_fit:
            line += (f"; residual from the size fit {size_fit[reg]['residual']:+.4f} "
                     f"({size_fit[reg]['z']:+.1f} SD, {size_fit[reg]['n_regions_beyond']} of "
                     f"{size_fit['n']} other regions at least as far)")
        lines.append(line)
        for name, v in tt.items():
            trows.append(dict(region=reg, test=name, **v))
    if size_fit:
        lines.append(f"Size fit through the {size_fit['n']} regions not shown on the left: "
                     f"slope {size_fit['slope']:+.4f} per log10 unit, r = {size_fit['r']:+.2f}, "
                     f"residual SD {size_fit['sd']:.4f}")
    per_sp = (t.assign(_r=t["label"] + "/" + t["run_id"], _a=t["label"] + "/" + t["animal"])
              .groupby("species").agg(runs=("_r", "nunique"), animals=("_a", "nunique"),
                                      datasets=("label", "nunique")))
    lines.append(f"n = {len(species)} species ({G} vs {n_o} others); {n_runs} selected runs, "
                 f"{n_animals} animals: " + "; ".join(
                     f"{s_} {per_sp.loc[s_, 'runs']} runs / {per_sp.loc[s_, 'animals']} "
                     f"animals / {per_sp.loc[s_, 'datasets']} dataset(s)" for s_ in species))
    lines.append("species values: " + "; ".join(
        f"{s_}: " + ", ".join(f"{ff.short_region(r)} {val.loc[s_, r]:.4g}" for r in show)
        for s_ in species))
    with pd.ExcelWriter(out + ".xlsx") as xw:
        t.to_excel(xw, sheet_name="runs", index=False)
        ds_val.pivot_table(index=["species", "label"], columns="region", values="v"
                           ).reset_index().to_excel(xw, sheet_name="datasets", index=False)
        val.reset_index().to_excel(xw, sheet_name="species", index=False)
        pd.DataFrame(trows).to_excel(xw, sheet_name="tests", index=False)
    na = math.comb(len(species), len(group))
    regs_txt = ", ".join(ff.short_region(r) for r in show)
    legend = (
        f"{panel} | {spec['tag']} of six regions, {G} versus the other species. Left, per "
        f"region: small dots, the {n_runs} selected runs ({spec['what']}); large markers, "
        f"species values (filled, {G}; open, other species), mean over the runs of each animal, "
        f"then over the animals of each dataset, then over the datasets of the species; "
        f"horizontal bars, means of the two groups; brackets, exact two-sided test of the "
        f"group difference with the species as unit ({na} assignments of {len(group)} species "
        f"among {len(species)}, smallest possible P = 1/{na} = {1 / na:.3f}; not corrected for "
        f"the {len(show)} regions; * P < 0.05, n.s. not significant). Right, every common "
        f"region (the region containing the seed excluded): group difference of the same index "
        f"against the group difference in log10 relative size of the region; dashed line and "
        f"grey band, linear fit through the regions not shown on the left ±2 s.d. "
        f"n = {len(species)} species, {n_animals} animals, {n_runs} runs.")
    _legend_txt(
        out, f"Is the {spec['tag']} of {regs_txt} different in {' and '.join(group)} from the "
             f"other species, beyond what the size of the region predicts?",
        f"{spec['what'][0].upper() + spec['what'][1:]}, in every selected run ({source}); "
        f"species value = run -> animal -> dataset -> species means (an animal with more runs, "
        f"or a dataset with more animals, does not weigh more). Left: small dots = runs, large "
        f"markers = species (filled = {G}, open = others), bars = group means, bracket = test. "
        f"Right: size control (see TEST).",
        f"One test per region, exact two-sided, over the {na} ways of choosing {len(group)} of "
        f"the {len(species)} species (p >= 1/{na} = {1 / na:.3f}); statistic = mean of {G} minus "
        f"mean of the others; NOT corrected for the {len(show)} regions (each region is its own "
        f"question). The species is the unit of inference: runs of one species share its "
        f"brain, protocol and scanner, so a test over runs would count {n_runs} runs as "
        f"{n_runs} independent observations of a species difference (pseudoreplication; "
        f"Hurlbert 1984 Ecol Monogr 54:187; Lazic 2010 BMC Neurosci 11:5; Aarts et al. 2014 Nat "
        f"Neurosci 17:491). Robustness: the same test with one dataset per species, every "
        f"combination of the datasets of the species that have two. Size: the group "
        f"difference of a region is compared with the linear fit, across the regions not "
        f"tested, of the group difference against the group difference in log10 relative "
        f"size (region voxels / voxels of all common regions, per run, same hierarchy).",
        "\n".join(lines),
        "Hypothesis: in anthropoid primates with a large association cortex, posterior parietal "
        "(and prefrontal) cortex contains territory anti-correlated with the default-mode "
        "network -- the dorsal attention and frontoparietal systems (Fox et al. 2005 PNAS "
        "102:9673; Vincent et al. 2007 Nature 447:83; Mantini et al. 2011 J Neurosci 31:12954) "
        "--, whereas in rodents posterior parietal cortex is coupled with retrosplenial cortex "
        "within the default-mode-like network (Lu et al. 2012 PNAS 109:3979; Whitesell et al. "
        "2021 Neuron 109:545); the dorsal premotor cortex holds the frontal eye field, frontal "
        "node of the dorsal attention network (Corbetta & Shulman 2002 Nat Rev Neurosci 3:201). "
        "Somatosensory cortex, V1 and striatum are the comparison regions. A single-run map "
        "is noisier than a group map: noise moves every percentile of a region away from its "
        "median and pulls tail fractions towards their chance value, more in noisier datasets. "
        "The group contrast confounds brain size, gyrification, arousal state and acquisition; "
        "Human and Macaque are the only catarrhines sampled, so a difference is that of the two "
        "catarrhines, not of primates (Marmoset and Mouselemur are in the other group).",
        legend=legend)
    return dict(tests=tests, size_fit=size_fit, n_runs=n_runs, n_animals=n_animals,
                species_values=val)


def fig6d_all(d6, sup, run_index, selections_final, members_by_label, bids_configs,
              species_of_label, divergence, cfg):
    """6d for EVERY regional index of the run table, all in 06_figures/Fig6
    (6d_<index>.png/.pdf/.txt/.xlsx): upper and lower percentile (P80, P20 with
    TAIL_Q = 0.2), median r, % in the top / bottom tail of the map, % of voxels
    with r > 0 and r < 0. Same runs, regions, group and test in every panel;
    the index of PROFILE_METHOD comes first. `sup` is kept for the call
    signature (S6d is no longer written: every index is a 6d)."""
    if run_index is None or run_index.empty:
        print("  [6d] no regional index table (07_state) -- skipped")
        return None
    q = float(run_index["tail_q"].iloc[0]) if "tail_q" in run_index else cfg["TAIL_Q"]
    specs = region_index_specs(q)
    fig6d = cfg.get("FIG6D") or {}
    regions = [str(r).lower() for r in fig6d.get(
        "regions", ("posterior_parietal_cortex", "prefrontal_cortex_(pfc)",
                    "somatosensory_cortex", "visual_striate_cortex", "striatum",
                    "motor_and_premotor"))]
    group = tuple(fig6d.get("group", ("Human", "Macaque")))
    resolvers = {c["label"]: ff.StateResolver(c) for c in bids_configs}

    def animal_of(lab, key):
        r = resolvers.get(lab)
        return r.animal(key) if r is not None else str(key[0])

    tab = region_index_table(run_index, selected_run_keys(selections_final, members_by_label),
                             animal_of)
    if tab.empty:
        print("  [6d] no selected run in the regional index table -- skipped")
        return None
    # shares of voxels > 0 / < 0 in % (the run table stores them as fractions)
    for frac, pct in (("frac_positive", "positive_pct"), ("frac_negative", "negative_pct")):
        if frac in tab:
            tab[pct] = 100.0 * tab[frac].astype(float)
    absent = [r for r in regions if r not in set(tab["region"])]
    if absent:
        print(f"  [6d] [WARN] FIG6D regions not among the common regions: {absent}")
    method = cfg["PROFILE_METHOD"]
    order = ([c for c, s_ in specs.items() if s_["method"] == method]
             + [c for c, s_ in specs.items() if s_["method"] != method])
    src_txt = (f"runs of the {cfg['N_SELECT']} {cfg.get('UNIT', 'run')}s kept per dataset by "
               f"the final selection round")
    out = {}
    for col in order:
        s_ = specs[col]
        if col not in tab or not tab[col].notna().any():
            print(f"  [6d] {col}: not in the run table -- skipped")
            continue
        out[col] = fig6d_region_values(opj(d6, f"6d_{s_['file']}"), tab, species_of_label,
                                       divergence, col, s_, regions, group=group,
                                       source=src_txt, panel="6d")
    print(f"  [6d] {len(out)} panel(s) in Fig6: "
          + ", ".join(f"6d_{specs[c]['file']}" for c in out))
    return out


# ═══════════════════════════════════════════════════════════════════════════
# Diagnostics and folder organisation
# ═══════════════════════════════════════════════════════════════════════════

def global_level_table(profiles_by_label, sel_r1, sel_final, common, value_col):
    """Mean of the raw regional profile (all common regions) of each run: the
    'global level' that Pearson-based selection ignores. Kept vs not kept."""
    rows = []
    for lab, profs in profiles_by_label.items():
        lv = {_key(k): float(np.nanmean(fa._ordered_common_values(p, common, value_col)))
              for k, p in profs.items()}
        s = pd.Series(lv)
        row = dict(label=lab, global_level_all_runs=float(s.median()))
        for tag, sel in (("round1", sel_r1), ("final", sel_final)):
            kept = {_key(k) for k in sel.get(lab, [])}
            m = s.index.isin(list(kept))
            if m.any() and (~m).any():
                u = st.mannwhitneyu(s[m], s[~m], alternative="two-sided")[0]
                row.update({f"global_level_kept_{tag}": float(s[m].median()),
                            f"global_level_not_kept_{tag}": float(s[~m].median()),
                            f"auc_kept_vs_not_{tag}": float(u / (m.sum() * (~m).sum()))})
        rows.append(row)
    return pd.DataFrame(rows)


def dataset_diagnostics(path, D, runs_by_label, df_qc_by_label, global_levels=None):
    rows = []
    if D is not None and not D.empty:
        md = D.where(~np.eye(len(D), dtype=bool)).mean(axis=1)
        mad = 1.4826 * (md - md.median()).abs().median()
        z = (md - md.median()) / mad if mad else md * np.nan
        for lab in D.index:
            rows.append(dict(label=lab, mean_distance_to_other_datasets=md[lab],
                             robust_z=z[lab]))
    diag = pd.DataFrame(rows)
    cov_rows = []
    for lab, profs in runs_by_label.items():
        frames = [p[["base_name", "n_voxels", "n_voxels_valid"]] for p in profs.values()
                  if {"n_voxels", "n_voxels_valid"} <= set(p.columns)]
        if frames:
            f = pd.concat(frames)
            f = f[f["n_voxels"] > 0]
            if len(f):
                cov_rows.append(dict(label=lab, median_fraction_of_region_with_signal=float(
                    (f["n_voxels_valid"] / f["n_voxels"]).median())))
    if cov_rows:
        diag = diag.merge(pd.DataFrame(cov_rows), on="label", how="outer") if not diag.empty \
            else pd.DataFrame(cov_rows)
    qrows = []
    for lab, dq in (df_qc_by_label or {}).items():
        if dq is None or dq.empty:
            continue
        r = dict(label=lab)
        for c in ("func_TSNR", "func_mean_fd", "func_mean_dvars"):
            if c in dq:
                r[c] = float(pd.to_numeric(dq[c], errors="coerce").median())
        qrows.append(r)
    if qrows:
        diag = diag.merge(pd.DataFrame(qrows), on="label", how="outer") if not diag.empty \
            else pd.DataFrame(qrows)
    if global_levels is not None and not global_levels.empty:
        diag = diag.merge(global_levels, on="label", how="outer") if not diag.empty \
            else global_levels
    if not diag.empty:
        diag.to_excel(path, index=False)
        print(f"  [saved] {path}")
    return diag


def collect_legends(fig_dir):
    """06_figures/Figure_legends.txt: the LEGEND section of every panel .txt of
    Fig5, Fig6 and supplementary (S5a, S6d), in panel order."""
    out = []
    for sub, title in (("Fig5", "Figure 5"), ("Fig6", "Figure 6"),
                       ("supplementary", "Supplementary panels")):
        d = opj(fig_dir, sub)
        if not os.path.isdir(d):
            continue
        block = []
        for f in sorted(os.listdir(d)):
            if not f.endswith(".txt"):
                continue
            txt = open(opj(d, f), encoding="utf-8").read()
            if txt.startswith("LEGEND\n"):
                body = txt.split("\n\n", 1)[0].split("\n", 2)[-1].strip()
                block.append(body)
        if block:
            out.append(f"{title}\n{'=' * len(title)}\n\n" + "\n\n".join(sorted(block)))
    if out:
        path = opj(fig_dir, "Figure_legends.txt")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("FIGURE LEGENDS (exact P, n and test of each panel: its own .txt)\n\n"
                     + "\n\n\n".join(out) + "\n")
        print(f"  [saved] {path}")


def organise_figures(fig_dir):
    sup = opj(fig_dir, "supplementary")
    os.makedirs(sup, exist_ok=True)
    n = 0
    for name in os.listdir(fig_dir):
        if name in MAIN_DIRS or name == "Figure_legends.txt":
            continue
        src = opj(fig_dir, name)
        dst = opj(sup, name)
        if os.path.isdir(dst):
            shutil.rmtree(dst)
        elif os.path.exists(dst):
            os.remove(dst)
        shutil.move(src, dst)
        n += 1
    print(f"  [figures] {n} item(s) moved to 06_figures/supplementary/; paper panels in "
          f"06_figures/Fig5 and 06_figures/Fig6")
    collect_legends(fig_dir)


# ═══════════════════════════════════════════════════════════════════════════
# 5e -- does the selection make the species stand out, beyond its construction?
# ═══════════════════════════════════════════════════════════════════════════
# Unit of the test: the PAIR OF DATASETS. Two runs of a dataset pair share the two
# datasets (protocol, scanner, animals, selection), so the hundreds of run pairs
# of one dataset pair are one observation, not hundreds. With 10 datasets there
# are 45 dataset pairs, 3 of which are "same species" (Human, Macaque, Marmoset).
# Statistic: gap = mean over same-species dataset pairs of their median run-to-run
# r - the same mean over different-species dataset pairs; gain = gap(selected) -
# gap(all runs) (the selection x species interaction). Null: every assignment of
# the species to the datasets that keeps the observed species structure (two
# datasets for 3 species, one for 4: 3150 assignments), enumerated exactly.

_5E_BEGIN = "=== SELECTION x SPECIES (Figure 5e) ==="
_5E_END = "=== end of Figure 5e ==="


def _n_partitions(n, sizes):
    """Number of partitions of n items into blocks of the given sizes."""
    num = math.factorial(n)
    for s in sizes:
        num //= math.factorial(s)
    for s in set(sizes):
        num //= math.factorial(sizes.count(s))
    return num


def _partitions(items, sizes):
    """Every partition of `items` into blocks whose sizes are the multiset
    `sizes`, each exactly once (blocks ordered by their first item)."""
    items = list(items)
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for s in sorted(set(sizes)):
        rem = list(sizes)
        rem.remove(s)
        for comp in itertools.combinations(rest, s - 1):
            left = [x for x in rest if x not in comp]
            for tail in _partitions(left, rem):
                yield [(first,) + comp] + tail


def _pairwise_r(X):
    """Pearson r between every two rows of X over their common finite columns
    (pairwise complete); NaN when fewer than 3 common values."""
    F = np.isfinite(X).astype(float)
    Z = np.where(F > 0, X, 0.0)
    n = F @ F.T
    sx = Z @ F.T
    sy = sx.T
    sxx = (Z ** 2) @ F.T
    syy = sxx.T
    sxy = Z @ Z.T
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = sxy - sx * sy / n
        vx = sxx - sx ** 2 / n
        vy = syy - sy ** 2 / n
        R = cov / np.sqrt(vx * vy)
    R[n < 3] = np.nan
    return R


def dataset_pair_medians(profiles_by_label, labels, common, value_col, norm, keep=None):
    """M[i, j] = median Pearson r between the profiles of the units of dataset
    labels[i] and those of labels[j] (i != j: every cross pair; i == j: every
    pair of two different units of the dataset). keep: {label: set of unit
    keys} to restrict to (None = every unit). Returns (M, n_units)."""
    vecs, ds = [], []
    for i, lab in enumerate(labels):
        profs = profiles_by_label.get(lab, {})
        keys = [k for k in profs if keep is None or _key(k) in keep.get(lab, set())]
        for k in keys:
            v = fa._ordered_common_values(profs[k], common, value_col)
            if np.isfinite(v).sum() >= 3:
                vecs.append(fa._normalize_values(v, norm))
                ds.append(i)
    n = len(labels)
    M = np.full((n, n), np.nan)
    n_units = np.bincount(np.asarray(ds, int), minlength=n) if ds else np.zeros(n, int)
    if not vecs:
        return M, n_units
    R = _pairwise_r(np.vstack(vecs))
    ds = np.asarray(ds)
    for i in range(n):
        ii = np.flatnonzero(ds == i)
        for j in range(i, n):
            jj = np.flatnonzero(ds == j)
            if i == j:
                block = R[np.ix_(ii, ii)][np.triu_indices(len(ii), 1)]
            else:
                block = R[np.ix_(ii, jj)].ravel()
            block = block[np.isfinite(block)]
            if block.size:
                M[i, j] = M[j, i] = float(np.median(block))
    return M, n_units


def _same_species(species, n):
    return np.array([[species[i] == species[j] for j in range(n)] for i in range(n)])


def _gap(M, same):
    """mean over same-species dataset pairs - mean over different-species pairs
    (off-diagonal cells only)."""
    iu = np.triu_indices(M.shape[0], 1)
    v, s = M[iu], same[iu]
    ok = np.isfinite(v)
    if not (s & ok).any() or not (~s & ok).any():
        return np.nan
    return float(np.mean(v[s & ok]) - np.mean(v[~s & ok]))


def species_gain_test(M_all, M_sel, species, max_exact=200000, n_mc=20000, seed=0):
    """gap(all), gap(selected), gain = gap(selected) - gap(all), and their null
    over every assignment of the species to the datasets that keeps the observed
    structure (exact enumeration up to max_exact assignments, Monte Carlo
    otherwise). One-sided p (>=), the observed assignment included."""
    n = len(species)
    obs_same = _same_species(species, n)
    obs = dict(gap_all=_gap(M_all, obs_same), gap_sel=_gap(M_sel, obs_same))
    obs["gain"] = obs["gap_sel"] - obs["gap_all"]
    sizes = sorted(pd.Series(species).value_counts().tolist(), reverse=True)
    n_part = _n_partitions(n, sizes)
    null = {k: [] for k in obs}

    def _add(same):
        ga, gs = _gap(M_all, same), _gap(M_sel, same)
        null["gap_all"].append(ga)
        null["gap_sel"].append(gs)
        null["gain"].append(gs - ga)

    if n_part <= max_exact:
        exact = True
        for part in _partitions(range(n), sizes):
            sp = np.empty(n, int)
            for b, block in enumerate(part):
                sp[list(block)] = b
            _add(_same_species(sp, n))
    else:
        exact = False
        rng = np.random.default_rng(seed)
        for _ in range(n_mc):
            _add(_same_species(rng.permutation(np.asarray(species)), n))
    null = {k: np.asarray(v, float) for k, v in null.items()}
    out = dict(obs, exact=exact, n_null=int(len(null["gain"])), n_assignments=int(n_part),
               null=null)
    for k in ("gap_all", "gap_sel", "gain"):
        nl = null[k][np.isfinite(null[k])]
        ge = int(np.sum(nl >= obs[k] - 1e-12))
        out[f"p_{k}"] = (ge / len(nl)) if exact else (1 + ge) / (len(nl) + 1)
        out[f"min_p_{k}"] = 1.0 / len(nl) if exact else 1.0 / (len(nl) + 1)
    return out


def _load_state_file(root):
    import pickle
    path = opj(root, "07_state", "selection_state.pkl")
    if not os.path.isfile(path):
        return None
    with open(path, "rb") as f:
        return pickle.load(f)


def _species_gain_run(st, labels):
    """Matrices and test of one run (its state)."""
    cfg = st["cfg"]
    sol = st["species_of_label"]
    keep = selected_run_keys(st["selections_final"], st["members_by_label"])
    keep_units = {lab: {_key(k) for k in v} for lab, v in st["selections_final"].items()}
    profs = st["profiles_by_label"]
    # unit-level profiles: the selected UNITS (with UNIT='run' units are runs)
    M_all, n_all = dataset_pair_medians(profs, labels, st["common"], cfg["VALUE_COL"],
                                        cfg["ZSCORE"])
    M_sel, n_sel = dataset_pair_medians(profs, labels, st["common"], cfg["VALUE_COL"],
                                        cfg["ZSCORE"], keep=keep_units)
    species = [sol[l] for l in labels]
    res = species_gain_test(M_all, M_sel, species)
    res.update(M_all=M_all, M_sel=M_sel, n_all=n_all, n_sel=n_sel, species=species,
               n_runs_kept=sum(len(v) for v in keep.values()),
               leave_one_out=bool(cfg.get("LEAVE_ONE_OUT", True)),
               final_round=st.get("final_round"))
    return res


def _slopegraph(ax, res, labels, header, ylim):
    n = len(labels)
    same = _same_species(res["species"], n)
    iu = np.triu_indices(n, 1)
    a, b, s = res["M_all"][iu], res["M_sel"][iu], same[iu]
    for x0, x1, intra in zip(a, b, s):
        if not intra:
            ax.plot([0, 1], [x0, x1], color=fa.INTER_COLOR, lw=0.4, alpha=0.35, zorder=1)
    ends = []
    for (i, j), x0, x1, intra in zip(zip(*iu), a, b, s):
        if intra:
            ax.plot([0, 1], [x0, x1], color=fa.INTRA_COLOR, lw=0.9, zorder=3,
                    marker="o", ms=2.2, mew=0)
            ends.append([x1, res["species"][i]])
    # species names at the right end, pushed apart so that they never overlap
    ends.sort(key=lambda e: e[0])
    step = 0.045 * (ylim[1] - ylim[0])
    ys = [e[0] for e in ends]
    for k in range(1, len(ys)):
        ys[k] = max(ys[k], ys[k - 1] + step)
    over = ys[-1] - (ylim[1] - step / 2) if ys else 0
    if over > 0:
        ys = [y - over for y in ys]
    for (y_end, name), y_txt in zip(ends, ys):
        ax.text(1.05, y_txt, name, color=fa.INTRA_COLOR, fontsize=5, va="center", ha="left")
    m = {k: (float(np.nanmean(v[s])), float(np.nanmean(v[~s]))) for k, v in (("a", a), ("b", b))}
    for x, (mi, mo) in ((0, m["a"]), (1, m["b"])):
        ax.plot([x], [mo], marker="s", ms=4.2, color=fa.INTER_COLOR, mec="white", mew=0.4,
                zorder=5)
        ax.plot([x], [mi], marker="o", ms=4.6, color=fa.INTRA_COLOR, mec="white", mew=0.4,
                zorder=5)
        xg = x - 0.09 if x == 0 else x + 0.46
        ax.annotate("", xy=(xg, mi), xytext=(xg, mo),
                    arrowprops=dict(arrowstyle="<->", lw=0.6, color="#333333",
                                    shrinkA=0, shrinkB=0))
        ax.text(xg + (-0.03 if x == 0 else 0.03), (mi + mo) / 2, "gap",
                ha="right" if x == 0 else "left", va="center", fontsize=5, color="#333333")
    ax.plot([0, 1], [m["a"][1], m["b"][1]], color=fa.INTER_COLOR, lw=1.6, zorder=4)
    ax.plot([0, 1], [m["a"][0], m["b"][0]], color=fa.INTRA_COLOR, lw=1.6, zorder=4)
    ax.set_xlim(-0.28, 1.72)
    ax.set_ylim(*ylim)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["all runs", "selected"])
    ax.set_ylabel("Median r between the runs\nof two datasets")
    ax.set_title(header, loc="left", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _matrix_panel(ax, M, labels, species, vmin, vmax, header, show_y):
    n = len(labels)
    D = M.copy()
    np.fill_diagonal(D, np.nan)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#E6E6E6")
    im = ax.imshow(D, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
    same = _same_species(species, n)
    for i in range(n):
        for j in range(n):
            if i != j and same[i, j]:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                           ec="white", lw=1.4, zorder=3))
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                           ec=fa.INTRA_COLOR, lw=0.8, zorder=4))
    names = [ff.pretty_label(l) for l in labels]
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(names, rotation=90, fontsize=5)
    ax.set_yticklabels(names if show_y else [], fontsize=5)
    ax.tick_params(length=0)
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
    ax.set_title(header, loc="left", fontweight="bold")
    return im


def _null_panel(ax, res, color, header, xlim, show_x):
    nl = res["null"]["gain"]
    nl = nl[np.isfinite(nl)]
    bins = np.linspace(xlim[0], xlim[1], 41)
    ax.hist(np.clip(nl, xlim[0], xlim[1]), bins=bins, color="#BDBDBD", lw=0)
    top = ax.get_ylim()[1]
    ax.axvline(res["gain"], color=color, lw=1.2)
    ax.axvline(0, color="#777777", lw=0.5, ls=(0, (2, 2)))
    ax.text(res["gain"], top * 0.97, " " + (fa.stars(res["p_gain"]) or "n.s."), color=color,
            ha="left", va="top", fontsize=6, fontweight="bold")
    ax.set_xlim(*xlim)
    ax.set_yticks([])
    ax.set_title(header, loc="left", fontweight="bold")
    for k in ("top", "right", "left"):
        ax.spines[k].set_visible(False)
    if show_x:
        ax.set_xlabel("Gain = gap(selected) − gap(all runs)")
    else:
        ax.set_xticklabels([])


def fig5e_selection_species_gain(main_root, control_root=None):
    """Figure 5e (06_figures/Fig5/5e_selection_species_gain.png/.pdf/.txt/.xlsx)
    from the saved part-1 states of the main run and, when it exists, of the
    control run (one common target, LEAVE_ONE_OUT=False). Also appends the
    result to Fig5/5b_intra_vs_inter_species.txt and refreshes
    Figure_legends.txt. Returns the results dict."""
    st = _load_state_file(main_root)
    if st is None:
        print(f"  [5e] no 07_state/selection_state.pkl under {main_root} -- skipped")
        return None
    ctl = _load_state_file(control_root) if control_root else None
    sol = st["species_of_label"]
    labels = ff.order_labels([l for l in sol if l in st["profiles_by_label"]], sol,
                             dict(fa.DIVERGENCE_TIME_MYA))
    print(f"\n{'=' * 70}\n  Figure 5e: selection x species (dataset pairs)\n{'=' * 70}")
    res = {"main": _species_gain_run(st, labels)}
    if ctl is not None:
        if ctl.get("cfg", {}).get("LEAVE_ONE_OUT", True):
            print("  [5e] [WARN] the control state was made with LEAVE_ONE_OUT=True -- not a "
                  "common-target control, ignored")
        elif set(ctl["profiles_by_label"]) != set(labels) or ctl["common"] != st["common"]:
            print("  [5e] [WARN] the control run has other datasets or regions than the main "
                  "run -- ignored (rerun RUN_CONTROL)")
        else:
            res["control"] = _species_gain_run(ctl, labels)
    else:
        print(f"  [5e] no control state under {control_root} -- main run only "
              f"(RUN_CONTROL=True draws the control)")
    for k, r in res.items():
        print(f"  [5e] {k}: gap all runs {r['gap_all']:+.4f} (p = {r['p_gap_all']:.4f}), "
              f"selected {r['gap_sel']:+.4f} (p = {r['p_gap_sel']:.4f}), gain {r['gain']:+.4f} "
              f"(p = {r['p_gain']:.4f}; {'exact, ' if r['exact'] else ''}{r['n_null']} "
              f"assignments)")

    fig_dir = opj(st["cfg"]["OUTPUT_ROOT"], "06_figures")
    d5 = opj(fig_dir, "Fig5")
    os.makedirs(d5, exist_ok=True)
    out = opj(d5, "5e_selection_species_gain")
    main, cont = res["main"], res.get("control")
    mats = [main["M_all"], main["M_sel"]] + ([cont["M_sel"]] if cont else [])
    offd = np.concatenate([m[~np.eye(len(labels), dtype=bool)] for m in mats])
    offd = offd[np.isfinite(offd)]
    vmin, vmax = float(np.percentile(offd, 1)), float(np.percentile(offd, 99))
    ylo = min(float(np.nanmin(m[np.triu_indices(len(labels), 1)])) for m in mats)
    yhi = max(float(np.nanmax(m[np.triu_indices(len(labels), 1)])) for m in mats)
    pad = 0.06 * (yhi - ylo)
    ylim = (ylo - pad, yhi + pad)
    gains = np.concatenate([r["null"]["gain"] for r in res.values()]
                           + [np.array([r["gain"] for r in res.values()])])
    gains = gains[np.isfinite(gains)]
    lim = max(abs(float(np.percentile(gains, 0.5))), abs(float(np.percentile(gains, 99.5))),
              *(abs(r["gain"]) * 1.15 for r in res.values()))
    xlim = (-lim, lim)
    with fa.paper_style():
        fig = plt.figure(figsize=(fa.PAPER_WIDTH["double"], 132 * fa.MM))
        outer = fig.add_gridspec(2, 1, height_ratios=[1.05, 1.0], hspace=0.62)
        top = outer[0].subgridspec(1, 5, width_ratios=[1, 1, 0.05, 0.42, 0.9], wspace=0.12)
        ax_a1 = fig.add_subplot(top[0, 0])
        ax_a2 = fig.add_subplot(top[0, 1])
        im = _matrix_panel(ax_a1, main["M_all"], labels, main["species"], vmin, vmax,
                           "a  all runs", True)
        _matrix_panel(ax_a2, main["M_sel"], labels, main["species"], vmin, vmax,
                      "selected", False)
        cax = fig.add_subplot(top[0, 2])
        cb = fig.colorbar(im, cax=cax)
        cb.set_label("Median r between the runs\nof the two datasets")
        cb.outline.set_linewidth(0.4)
        ax_key = fig.add_subplot(top[0, 4])
        ax_key.axis("off")
        key = ("Each cell: one pair of datasets.\n"
               "Blue frames: the 3 pairs of datasets\n"
               "of the SAME species (the only ones\n"
               "that can separate species from\n"
               "dataset). Grey diagonal: runs of\n"
               "one dataset (not used).\n\n"
               "gap = mean of the blue cells −\n"
               "mean of the other cells.\n"
               "gain = gap(selected) − gap(all runs):\n"
               "does the selection make the\n"
               "species stand out?")
        ax_key.text(0.0, 0.98, key, transform=ax_key.transAxes, ha="left", va="top",
                    fontsize=5.5, linespacing=1.25, color="#333333")
        bot = outer[1].subgridspec(2, 3, width_ratios=[1, 1, 1.15], wspace=0.55, hspace=0.75)
        ax_b = fig.add_subplot(bot[:, 0])
        _slopegraph(ax_b, main, labels, "b  main selection\n(target: the other species)",
                    ylim)
        ax_c = fig.add_subplot(bot[:, 1])
        if cont:
            _slopegraph(ax_c, cont, labels, "c  control selection\n(one common target)", ylim)
            ax_c.set_ylabel("")
        else:
            ax_c.axis("off")
            ax_c.text(0.5, 0.5, "c  control run not available\n(RUN_CONTROL=True)",
                      ha="center", va="center", fontsize=6, transform=ax_c.transAxes)
        ax_d1 = fig.add_subplot(bot[0, 2])
        _null_panel(ax_d1, main, fa.INTRA_COLOR, "d  main", xlim, show_x=cont is None)
        if cont:
            ax_d2 = fig.add_subplot(bot[1, 2])
            _null_panel(ax_d2, cont, "#7B3294", "control", xlim, show_x=True)
        handles = [Line2D([], [], color=fa.INTRA_COLOR, marker="o", ms=4, lw=1.6,
                          label="same species, two datasets (mean)"),
                   Line2D([], [], color=fa.INTER_COLOR, marker="s", ms=4, lw=1.6,
                          label="different species (mean)"),
                   plt.Rectangle((0, 0), 1, 1, color="#BDBDBD",
                                 label="gain under every assignment of the species "
                                       "to the datasets")]
        fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
                   bbox_to_anchor=(0.5, -0.01), fontsize=5.5)
        fig.subplots_adjust(left=0.10, right=0.98, top=0.95, bottom=0.10)
        fa.save_paper(fig, out)

    # ---- tables -------------------------------------------------------------------
    names = [ff.pretty_label(l) for l in labels]
    with pd.ExcelWriter(out + ".xlsx") as xw:
        rows = []
        for k, r in res.items():
            rows.append(dict(run=k, leave_one_out=r["leave_one_out"],
                             final_round=r["final_round"], n_datasets=len(labels),
                             n_same_species_dataset_pairs=int(
                                 _same_species(r["species"], len(labels))[
                                     np.triu_indices(len(labels), 1)].sum()),
                             gap_all_runs=r["gap_all"], p_gap_all_runs=r["p_gap_all"],
                             gap_selected=r["gap_sel"], p_gap_selected=r["p_gap_sel"],
                             gain=r["gain"], p_gain=r["p_gain"], min_possible_p=r["min_p_gain"],
                             exact=r["exact"], n_assignments=r["n_null"]))
        pd.DataFrame(rows).to_excel(xw, sheet_name="tests", index=False)
        for k, r in res.items():
            for tag, M in (("all", r["M_all"]), ("selected", r["M_sel"])):
                if k == "control" and tag == "all":
                    continue
                pd.DataFrame(M, index=names, columns=names).to_excel(
                    xw, sheet_name=f"{k}_{tag}_median_r")
    print(f"  [saved] {out}.xlsx")

    # ---- legend and .txt -------------------------------------------------------------
    def _line(k, r):
        return (f"{k}: gap all runs = {r['gap_all']:+.3f} (p = {ff._fmt_p(r['p_gap_all'])}), "
                f"gap selected = {r['gap_sel']:+.3f} (p = {ff._fmt_p(r['p_gap_sel'])}), "
                f"gain = {r['gain']:+.3f}, p = {ff._fmt_p(r['p_gain'])} (one-sided; "
                f"{'exact over' if r['exact'] else 'Monte Carlo,'} {r['n_null']} assignments; "
                f"smallest possible p = {ff._fmt_p(r['min_p_gain'])}); "
                f"{r['n_runs_kept']} selected runs")
    result = "\n".join(_line(k, r) for k, r in res.items())
    n_pairs = len(labels) * (len(labels) - 1) // 2
    n_same = int(_same_species(main["species"], len(labels))[
        np.triu_indices(len(labels), 1)].sum())
    stars_main = fa.stars(main["p_gain"]) or "n.s."
    legend = (
        "5e | Does the selection make the species stand out? (a) Median Pearson r between the "
        f"regional profiles of the runs of every pair of the {len(labels)} datasets, all runs "
        f"(left) and runs kept by the final selection round (right); blue frames, the {n_same} "
        f"pairs of datasets of the same species, the only ones that separate species from "
        f"dataset; grey diagonal, pairs of runs of one dataset (not used). (b) Each line, one of "
        f"the {n_pairs} pairs of datasets, from all runs to selected runs (blue, same species; "
        f"orange, different species); large markers, means; gap, mean same-species minus mean "
        f"different-species. (c) Same for the control selection, in which every dataset is "
        f"ranked against one common target, so that two datasets of a species no longer share a "
        f"target of their own. (d) Gain = gap(selected) − gap(all runs) (vertical line) and its "
        f"distribution over every assignment of the species to the datasets that keeps the "
        f"observed species structure (grey; exact enumeration); one-sided; "
        f"*** P < 0.001, ** P < 0.01, * P < 0.05, n.s. not significant.")
    question = ("Does the selection make two datasets of the same species more similar than two "
                "datasets of different species -- beyond what it does to every pair of datasets, "
                "and beyond what its construction alone produces?")
    what = ("(a) For every pair of datasets, the median r between the profiles of their runs "
            "(dataset x dataset matrix), all runs and selected runs; the same-species pairs are "
            "framed. (b, c) One line per dataset pair, all runs -> selected, main and control "
            "selections. (d) Observed gain and its null distribution, main and control.")
    test = ("Unit: the pair of datasets (the run pairs of one dataset pair share the two "
            "datasets' protocol, scanner, animals and selection and are not independent; a "
            "test on run pairs, like the Wilcoxon of 5b, counts them as thousands of "
            "observations). Statistic: gap = mean over same-species dataset pairs of their "
            "median run-to-run r - mean over different-species dataset pairs; gain = "
            "gap(selected) - gap(all runs), i.e. the selection x species interaction. Null: "
            "every assignment of the species to the datasets that keeps the observed structure "
            "(how many datasets each species has), enumerated exactly; the same assignment is "
            "applied to both pools. One-sided (gain > 0, gap > 0). Each dataset pair weighs the "
            "same whatever its number of runs.")
    reading = (
        "The selection raises the similarity of every pair of datasets (all lines go up: every "
        "dataset is pulled toward a shared cross-species pattern), so the gap after selection "
        "alone says little; the gain is the quantity of interest. In the main selection two "
        "datasets of one species are ranked against the SAME leave-one-species-out target "
        "(rounds 1 and 3 onward), which brings them closer by construction, with no biology: "
        "part of the main gain can be mechanical. The control removes that shared target. A "
        "gain in the main run that persists in the control is not produced by the construction "
        "of the selection; a gain that vanishes in the control was. With 3 same-species dataset "
        "pairs the test has little power, and the smallest possible p is set by the number of "
        "assignments: read the gap of each line, not only p.")
    _legend_txt(out, question, what, test, result, reading, legend=legend)

    # ---- appended to the .txt of 5b ------------------------------------------------------
    p5b = opj(d5, "5b_intra_vs_inter_species.txt")
    if os.path.isfile(p5b):
        txt = open(p5b, encoding="utf-8").read()
        if _5E_BEGIN in txt:
            txt = txt[:txt.index(_5E_BEGIN)].rstrip() + "\n"
        txt += (f"\n{_5E_BEGIN}\nDoes the selection increase the separation of the species "
                f"beyond chance and beyond its construction? Dataset-pair level (see "
                f"5e_selection_species_gain.txt):\n{result}\n{_5E_END}\n")
        with open(p5b, "w", encoding="utf-8") as f:
            f.write(txt)
        print(f"  [5e] result added to {p5b}")
    collect_legends(fig_dir)
    print(f"  [5e] main gain {main['gain']:+.3f} {stars_main}"
          + (f"; control gain {cont['gain']:+.3f} {fa.stars(cont['p_gain']) or 'n.s.'}"
             if cont else ""))
    return res


# ═══════════════════════════════════════════════════════════════════════════
# Entry point
# ═══════════════════════════════════════════════════════════════════════════

def make_paper_figures(*, sri, fu, cfg, paths, bids_configs, species_of_label, common,
                       region_mask, final_round, profiles_by_label, runs_by_label,
                       members_by_label, selections_final, rankings, f1_all, f1_sel,
                       selections_r1=None,
                       group_all, group_all_source, fingerprints, all_units_index,
                       map_index_final, phylo_index, df_qc_by_label, divergence,
                       surface_pics=None, surface_ranges=None, run_index=None,
                       anat_common=None):
    """Every paper panel (06_figures/Fig5, Fig6) and its supplementary twins
    S5a, S6d. surface_pics / surface_ranges: the 5a renders made in part 1;
    run_index: the regional indices of every run (part 1)."""
    value_col, norm = cfg["VALUE_COL"], cfg["ZSCORE"]
    ref = cfg["REFERENCE_SPECIES"]
    seed_kind = cfg.get("PROFILE_KIND", "seed") == "seed"
    d5, d6 = opj(paths["figures"], "Fig5"), opj(paths["figures"], "Fig6")
    sup = opj(paths["figures"], "supplementary")
    for d in (d5, d6, sup):
        os.makedirs(d, exist_ok=True)
    out = {}
    states = fingerprints.get("dataset_state", {})
    dist_all = fingerprints.get("dist_all", pd.DataFrame())
    labels = ff.order_labels([c["label"] for c in bids_configs if c["label"] in species_of_label],
                             species_of_label, divergence)
    species = sorted({species_of_label[l] for l in labels})
    brain_mass = load_brain_mass(cfg.get("BRAIN_MASS_XLSX"), species,
                                 cfg["SPECIES_TO_TREE_TIP"], cfg.get("BRAIN_MASS_G"))
    out["brain_mass"] = brain_mass

    def _run(name, fn, *a, **k):
        """One panel failing must not cost the others."""
        try:
            out[name] = fn(*a, **k)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"  [{name}] [WARN] failed: {e}")

    # ---- Figure 5 ----------------------------------------------------------
    if seed_kind:
        _run("5a", fig5a_surfaces, opj(d5, "5a_surface_maps_selected_vs_excluded"), cfg,
             surface_maps_for_5a(all_units_index, map_index_final), surface_pics or {},
             surface_ranges or {}, labels, species_of_label,
             tuple(cfg.get("SURFACE_PNG_KEEP") or ("",)),
             supp_out=opj(sup, "S5a_surface_maps_all_selected_excluded"))
    selected_units = {}
    for lab, profs in profiles_by_label.items():
        kept = {_key(x) for x in selections_final.get(lab, [])}
        selected_units[lab] = {k: v for k, v in profs.items() if _key(k) in kept}
    _run("5b", fig5b_within_across, opj(d5, "5b_intra_vs_inter_species"), profiles_by_label,
         selected_units, common, species_of_label, norm, value_col,
         test=cfg.get("WITHIN_ACROSS_TEST", "wilcoxon"), unit=cfg.get("UNIT", "run"),
         intra=cfg.get("WITHIN_ACROSS_PAIRS", "other_dataset"))
    if seed_kind:
        _run("5c", fig5c_fingerprint, opj(d5, "5c_fingerprint_species"), f1_all, f1_sel, common,
             species_of_label, value_col, norm, region_mask, divergence, final_round)
    _run("5d", fig5d_rank, opj(d5, "5d_rank_before_after_selection"), fingerprints, states,
         final_round, ref)

    # ---- Figure 6 ----------------------------------------------------------
    if len(dist_all) >= 3:
        _run("6a", fig6a_distances, opj(d6, "6a_distance_divergence_brain_mass"), dist_all,
             states, brain_mass, group_all_source, ref)
        if seed_kind:
            _run("6b", fig6b_fingerprint_grid, opj(d6, "6b_fingerprint_vs_consensus"),
                 fingerprints, species_of_label, common, region_mask, divergence, ref)
        try:
            an = anesthesia_table(fingerprints.get("state_runs"), rankings)
            if not an.empty:
                an.drop(columns=["key"]).to_excel(opj(paths["group_maps"],
                                                      "arousal_state_runs.xlsx"), index=False)
        except Exception as e:
            print(f"  [6c] [WARN] run-level table failed: {e}")
        # 6c is drawn at STEP 8 (fig6c_f8b), where the F8b table exists; the F18
        # analysis stays in supplementary.
    if seed_kind:
        _run("6d", fig6d_all, d6, sup, run_index, selections_final, members_by_label,
             bids_configs, species_of_label, divergence, cfg)
    gl = global_level_table(profiles_by_label, selections_r1 or {}, selections_final,
                            common, value_col)
    dataset_diagnostics(opj(paths["group_maps"], "dataset_diagnostics.xlsx"),
                        fingerprints.get("distance_matrix"), runs_by_label, df_qc_by_label,
                        global_levels=gl)
    return out
