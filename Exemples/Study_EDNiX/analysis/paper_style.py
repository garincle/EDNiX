"""
paper_style.py -- one visual vocabulary for every figure of the EDNiX paper.

Every paper script (Figure_*.py, Supplementary/Supp_Fig_*.py) takes its colours,
species order, fonts, sizes and panel letters from here, so that a species, a
dataset or a lobe looks the same in every figure.

  species  -> colour            SPECIES_COLORS (= Figure_2_anat.SPECIES_COLORS, unchanged)
  dataset  -> colour            DATASET_COLORS, one fixed colour per BIDS dataset for the
                                whole paper; a species with one dataset keeps its species
                                colour, the other datasets of a species get colours far
                                from it, so that datasets side by side are told apart.
                                Every point is a circle.
  lobe     -> colour            Atlas_nomenclature.LOBE_COLORS (= Figure 1 sunburst)
  regions  -> order             Atlas_nomenclature.atlas_order (lobe order of Figure 1)
  species  -> order             PHYLO_ORDER (rodents -> primates -> carnivores)

Sizes follow Nature's figure guide: 89 mm (1 column) or 183 mm (2 columns) wide,
5-7 pt text, bold lower-case panel letters at 8 pt, vector PDF (TrueType fonts) and
a 600 dpi PNG for review.
"""
import os
import re

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from Exemples.Study_EDNiX.atlases import Atlas_nomenclature as AN

MM = 1 / 25.4
WIDTH_1COL = 89 * MM
WIDTH_15COL = 120 * MM
WIDTH_2COL = 183 * MM
MAX_HEIGHT = 247 * MM

# ── species ──────────────────────────────────────────────────────────────────
PHYLO_ORDER = ["Mouse", "Rat", "Mouselemur", "Bat", "Marmoset",
               "Macaque", "Human", "Dog", "Cat", "Pig"]

SPECIES_COLORS = {
    "Mouse": "#1A5276", "Rat": "#2E86C1", "Bat": "#00A896",
    "Mouselemur": "#58D68D", "Marmoset": "#F4D03F", "Macaque": "#F5B041",
    "Human": "#D35400", "Dog": "#8E44AD", "Cat": "#C0392B", "Pig": "#6E2C00",
}

# ── datasets: (species, BIDS folder, label used by the Fig 5 pipeline) ───────
# Order inside a species fixes its marker in every figure.
DATASETS = [
    ("Mouse", "BIDS_Grandjean2", "Mouse"),
    ("Rat", "BIDS_Grandjean", "Rat"),
    ("Mouselemur", "BIDS_Garin", "Mouselemur"),
    ("Marmoset", "BIDS_Tian", "Marmoset_Tian"),
    ("Marmoset", "BIDS_Hori", "Marmoset_Hori"),
    ("Macaque", "BIDS_BenHamed", "Macaque_BenHamed"),
    ("Macaque", "BIDS_Zhu_Garin", "Macaque_ZhuGarin"),
    ("Macaque", "BIDS_Tremblay", "Macaque_Tremblay"),
    ("Human", "BIDS_Castrillon", "Human_Castrillon"),
    ("Human", "BIDS_Babayan", "Human_Babayan"),
    ("Human", "BIDS_Park", "Human_Park"),
    ("Human", "BIDS_Merida", "Human_Merida"),
    ("Dog", "BIDS_Boch_K9", "Dog"),
]
# One colour per dataset, the same in every figure. The first dataset of a species keeps
# the species colour; the others are chosen far from it and from each other (Okabe-Ito /
# Tol colours), because they are drawn side by side at the same x.
DATASET_COLORS = {
    "BIDS_Grandjean2": "#1A5276",   # Mouse (species colour)
    "BIDS_Grandjean": "#2E86C1",    # Rat (species colour)
    "BIDS_Garin": "#58D68D",        # Mouse lemur (species colour)
    "BIDS_Tian": "#F4D03F",         # Marmoset (species colour)
    "BIDS_Hori": "#882255",         # Marmoset, wine
    "BIDS_BenHamed": "#F5B041",     # Macaque (species colour)
    "BIDS_Zhu_Garin": "#44AA99",    # Macaque, teal
    "BIDS_Tremblay": "#332288",     # Macaque, indigo
    "BIDS_Castrillon": "#D35400",   # Human (species colour)
    "BIDS_Babayan": "#117733",      # Human, dark green
    "BIDS_Park": "#CC79A7",         # Human, pink
    "BIDS_Merida": "#000000",       # Human, black
    "BIDS_Boch_K9": "#8E44AD",      # Dog (species colour)
}
DATASET_MARKERS = ["o"]

# ── lobes ────────────────────────────────────────────────────────────────────
LOBE_COLORS = {k: tuple(c / 255 for c in v) for k, v in AN.LOBE_COLORS.items()}

# ── typography ───────────────────────────────────────────────────────────────
PAPER_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 6, "axes.labelsize": 6, "axes.titlesize": 6,
    "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
    "xtick.minor.size": 1.2, "ytick.minor.size": 1.2,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelpad": 1.5, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
    "legend.frameon": False, "legend.handletextpad": 0.3, "legend.columnspacing": 0.8,
    "legend.borderaxespad": 0.2, "legend.handlelength": 1.0,
    "lines.linewidth": 0.8, "lines.markersize": 2.5,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 150, "savefig.dpi": 600,
    "figure.facecolor": "white", "axes.facecolor": "white",
}
PANEL_LABEL_SIZE = 8


def _key(text):
    return re.sub(r"[^a-z0-9]", "", str(text).lower().replace("bids_", ""))


_DATASET_BY_KEY = {}
for _sp, _bids, _lab in DATASETS:
    for _k in (_bids, _lab, _bids.replace("BIDS_", ""), f"{_sp}_{_bids.replace('BIDS_', '')}"):
        _DATASET_BY_KEY[_key(_k)] = (_sp, _bids, _lab)


def dataset_info(name, species=None):
    """(species, BIDS folder, pipeline label) of a dataset given as a path, a BIDS
    folder name ('BIDS_Zhu_Garin') or a Fig 5 label ('Macaque_ZhuGarin'). Unknown
    names return (species, name, name)."""
    base = os.path.basename(str(name).rstrip("/"))
    hit = _DATASET_BY_KEY.get(_key(base))
    if hit is None and species is not None:
        hit = _DATASET_BY_KEY.get(_key(f"{species}_{base}"))
    return hit if hit is not None else (species, base, base)


def dataset_short(name, species=None):
    """Name printed on figures: 'Grandjean2', 'Zhu_Garin', ..."""
    return dataset_info(name, species)[1].replace("BIDS_", "")


def species_color(species):
    return SPECIES_COLORS.get(species, "#7F8C8D")


def dataset_color(name, species=None):
    """Fixed colour of a dataset (DATASET_COLORS); unknown dataset -> species colour."""
    sp, bids, _ = dataset_info(name, species)
    return DATASET_COLORS.get(bids, species_color(sp or species))


def dataset_marker(name, species=None):
    """Every dataset is drawn with circles (datasets differ by colour)."""
    return "o"


def dataset_offsets(names, spread=0.26):
    """x offset of each dataset inside one species column: datasets side by side,
    centred on the species tick (one dataset -> 0)."""
    names = list(names)
    if len(names) <= 1:
        return {n: 0.0 for n in names}
    xs = np.linspace(-spread, spread, len(names))
    return dict(zip(names, xs))


def phylo_sort(species):
    species = list(dict.fromkeys(species))
    known = [s for s in PHYLO_ORDER if s in species]
    return known + sorted(s for s in species if s not in PHYLO_ORDER)


def dataset_sort(names, species_of=None):
    """Datasets grouped by species (PHYLO_ORDER), in DATASETS order inside a species."""
    order = {b: i for i, (_, b, _) in enumerate(DATASETS)}

    def k(n):
        sp, b, _ = dataset_info(n, (species_of or {}).get(n))
        return (PHYLO_ORDER.index(sp) if sp in PHYLO_ORDER else 99, order.get(b, 99), str(n))
    return sorted(dict.fromkeys(names), key=k)


def region_order(names):
    return AN.atlas_order(names)


def lobe_color(region):
    lobe = AN.lobe_of(region)
    return LOBE_COLORS.get(lobe, (0.6, 0.6, 0.6)) if lobe else (0.6, 0.6, 0.6)


def region_label(name):
    """Readable region name: 'Somatosensory_Cx' -> 'Somatosensory Cx'."""
    return str(name).replace("_lobe", "").replace("_", " ")


# ── drawing helpers ──────────────────────────────────────────────────────────
def panel_label(ax, letter, x=-0.02, y=1.0, fig=None):
    """Bold lower-case panel letter at the top-left corner of `ax` (axes fraction,
    shifted left by the y tick labels)."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=PANEL_LABEL_SIZE,
            fontweight="bold", va="bottom", ha="right")


def row_label(fig, ax, letter, dx_mm=1.0):
    """Panel letter in the left margin of the figure, level with the top of `ax`."""
    bb = ax.get_position()
    fig.text(dx_mm * MM / fig.get_figwidth(), bb.y1, letter, fontsize=PANEL_LABEL_SIZE,
             fontweight="bold", va="top", ha="left")


def lobe_strip(ax, region, height=0.035):
    """Thin bar in the lobe colour of Figure 1 above the axes of `region`."""
    ax.add_patch(matplotlib.patches.Rectangle(
        (0, 1.0), 1, height, transform=ax.transAxes, clip_on=False,
        facecolor=lobe_color(region), edgecolor="none"))


VIOLIN_GREY = "#A6A6A6"


def violin(ax, values, x, color, width=0.7, alpha=0.25):
    """Half-transparent KDE violin with median (thick) and quartiles (thin)."""
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return
    if v.size < 3 or np.std(v) == 0:
        ax.hlines(np.median(v), x - width / 4, x + width / 4, colors=color, lw=1.0, zorder=3)
        return
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(v)
    y = np.linspace(v.min(), v.max(), 200)
    d = kde(y)
    d = d / d.max() * width / 2
    ax.fill_betweenx(y, x - d, x + d, color=color, alpha=alpha, lw=0, zorder=1)
    for q, lw in zip(np.percentile(v, [25, 50, 75]), (0.5, 1.0, 0.5)):
        dq = float(kde(q)[0]) / kde(y).max() * width / 2
        ax.hlines(q, x - dq, x + dq, colors=color, lw=lw, zorder=3)


def dots(ax, values, x, color, marker="o", spread=0.12, size=3.0, rng=None, edge=True):
    """Jittered points (one per subject / run)."""
    rng = rng if rng is not None else np.random.default_rng(42)
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return
    ax.scatter(x + rng.uniform(-spread, spread, v.size), v, s=size, marker=marker,
               facecolor=color, edgecolor="black" if edge else "none", linewidths=0.2,
               zorder=4, alpha=0.9)


def species_dataset_legend(datasets, species_of=None, ncol=None, marker_size=3.5):
    """Handles + labels: one entry per dataset (dataset colour, circle)."""
    handles, labels = [], []
    for d in dataset_sort(datasets, species_of):
        sp = dataset_info(d, (species_of or {}).get(d))[0]
        handles.append(Line2D([], [], ls="", marker="o",
                              markerfacecolor=dataset_color(d, sp), markeredgecolor="black",
                              markeredgewidth=0.2, markersize=marker_size))
        labels.append(f"{sp} {dataset_short(d, sp)}")
    return handles, labels


def legend_below(fig, axes, handles, labels, ncol=7, gap_mm=2.0, x=0.5, **kw):
    """Figure legend placed gap_mm under the lowest tick label of `axes` (whatever the
    label length), centred on x (figure fraction). save() crops to the tight bbox."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    y0 = min(ax.get_tightbbox(r).y0 for ax in np.ravel(axes) if ax.get_visible())
    y = fig.transFigure.inverted().transform((0, y0))[1] - gap_mm * MM / fig.get_figheight()
    return fig.legend(handles, labels, loc="upper center", ncol=min(len(handles), ncol),
                      bbox_to_anchor=(x, y), handlelength=0.8, frameon=False, **kw)


def species_legend(species):
    return ([Patch(facecolor=species_color(s), edgecolor="none") for s in phylo_sort(species)],
            phylo_sort(species))


def save(fig, path_no_ext, formats=("pdf", "png")):
    """Write <path>.pdf (vector, TrueType) and <path>.png (600 dpi); returns the paths."""
    os.makedirs(os.path.dirname(os.path.abspath(path_no_ext)), exist_ok=True)
    out = []
    for f in formats:
        p = f"{path_no_ext}.{f}"
        fig.savefig(p, dpi=600 if f == "png" else None, bbox_inches="tight", pad_inches=0.01)
        out.append(p)
    print("  [figure] " + ", ".join(out))
    return out


def compose_images(panels, out_path_no_ext, width=WIDTH_2COL, layout=None,
                   wspace_mm=2.0, hspace_mm=2.0, letters=None, crop=True, gutter_mm=3.0):
    """Lay out existing PNG panels (drawn at print size) as one lettered figure.

    panels : list of PNG paths, in reading order (None or a missing file = empty cell).
    layout : list of rows, each a list of panel indices, e.g. [[0], [1, 2]]; default:
             one panel per row.
    In a row, panels share one height and their widths follow their aspect ratios, so
    that the row fills `width` exactly (no empty band, nothing stretched). Panels
    drawn at print size for that place keep a scale close to 1 (fonts stay 5-7 pt).
    Letters a, b, c... follow the order of `panels`. Writes .pdf + .png."""
    import matplotlib.image as mpimg
    imgs = []
    for p in panels:
        if p is None or not os.path.isfile(p):
            if p is not None:
                print(f"  [compose] [WARN] missing panel: {p}")
            imgs.append(None)
            continue
        im = mpimg.imread(p)
        imgs.append(_autocrop(im) if crop else im)
    layout = layout or [[i] for i in range(len(imgs))]
    letters = letters or "abcdefghijklmnopqrstuvwxyz"
    ws, hs, gut = wspace_mm * MM, hspace_mm * MM, gutter_mm * MM
    rows = []
    for row in layout:
        asp = [(imgs[i].shape[1] / imgs[i].shape[0]) if imgs[i] is not None else 1.0 for i in row]
        h = (width - ws * (len(row) - 1) - gut * len(row)) / sum(asp)
        rows.append((row, asp, h))
    height = sum(h for _, _, h in rows) + hs * (len(rows) - 1)
    with plt.rc_context(PAPER_RC):
        fig = plt.figure(figsize=(width, height))
        y = height
        for row, asp, h in rows:
            y -= h
            x = 0.0
            for i, a in zip(row, asp):
                w = a * h
                ax = fig.add_axes([(x + gut) / width, y / height, w / width, h / height])
                ax.axis("off")
                if imgs[i] is not None:
                    ax.imshow(imgs[i], interpolation="lanczos")
                    fig.text(x / width, (y + h) / height, letters[i], fontsize=PANEL_LABEL_SIZE,
                             fontweight="bold", va="top", ha="left")
                x += gut + w + ws
            y -= hs
        paths = save(fig, out_path_no_ext)
        plt.close(fig)
    return paths


def _autocrop(im, tol=0.02):
    """Remove uniform white margins of an RGB(A) image array."""
    rgb = im[..., :3].astype(float)
    if rgb.max() > 1.5:
        rgb = rgb / 255.0
    nonwhite = (1.0 - rgb.min(axis=2)) > tol
    if im.shape[-1] == 4:
        nonwhite &= im[..., 3] > 0.02 * (255 if im.max() > 1.5 else 1)
    ys, xs = np.where(nonwhite)
    if ys.size == 0:
        return im
    return im[ys.min():ys.max() + 1, xs.min():xs.max() + 1]


def plot_by_dataset(df, metrics, out_path, ncols=None, letters=None,
                    height_mm=45, value_unit="dataset"):
    """One panel per metric; x = datasets grouped by species (PHYLO_ORDER); grey violin
    = distribution of the dataset, points in the dataset colour (DATASET_COLORS).
    Values beyond the axis limits are drawn on the axis limit as open circles, and their
    number is written in the corner, so that one aberrant session does not flatten the
    panel (axis limits = quartiles -/+ 4 IQR, within the data range).
    metrics: list of (column, axis label); a column absent from df gives an empty panel
    that says so (the letter is kept, so the numbering of the legend stays right); None
    leaves the cell empty (no letter), e.g. to start a new row."""
    import math
    metrics = list(metrics)
    if not metrics:
        print("  [plot_by_dataset] no metric -- nothing drawn")
        return None
    sp_of = dict(zip(df['bids_dir'], df['species']))
    dsets = dataset_sort(list(sp_of), sp_of)
    n = len(metrics)
    ncols = ncols or n
    nrows = math.ceil(n / ncols)
    letters = letters or "abcdefghijklmnop"
    rng = np.random.default_rng(42)
    with plt.rc_context(PAPER_RC):
        fig, axes = plt.subplots(nrows, ncols, squeeze=False,
                                 figsize=(WIDTH_2COL, (height_mm * nrows + 14) * MM))
        fig.subplots_adjust(left=0.06, right=0.995, top=0.97,
                            bottom=0.05, wspace=0.42, hspace=0.45)
        n_letter = 0
        for k, item in enumerate(metrics):
            ax = axes.flat[k]
            if item is None:
                ax.set_visible(False)
                continue
            col, label = item
            panel_label(ax, letters[n_letter], x=-0.16, y=1.0)
            n_letter += 1
            ax.set_ylabel(label)
            v_all = pd_to_float(df[col]) if col in df.columns else np.array([np.nan])
            if not np.isfinite(v_all).any():
                ax.text(0.5, 0.5, "not computed", ha="center", va="center",
                        transform=ax.transAxes, color="0.5")
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            q1, q3 = np.nanpercentile(v_all, (25, 75))
            iqr = q3 - q1
            lo = max(np.nanmin(v_all), q1 - 4 * iqr)      # Tukey-like fence, wide (4 IQR)
            hi = min(np.nanmax(v_all), q3 + 4 * iqr)
            pad = 0.05 * (hi - lo if hi > lo else abs(hi) + 1)
            lo_l, hi_l = lo - pad, hi + pad
            n_clip = 0
            for x, d in enumerate(dsets):
                sp = sp_of[d]
                c = dataset_color(d, sp)
                v = pd_to_float(df.loc[df['bids_dir'] == d, col])
                v = v[np.isfinite(v)]
                if v.size == 0:
                    continue
                inside = (v >= lo_l) & (v <= hi_l)
                n_clip += int((~inside).sum())
                violin(ax, v[inside], x, VIOLIN_GREY, width=0.75, alpha=0.45)
                dots(ax, v[inside], x, c, spread=0.15, size=2.5, rng=rng)
                if (~inside).any():
                    ax.scatter(np.full((~inside).sum(), x), np.clip(v[~inside], lo_l, hi_l),
                               s=5, marker="o", facecolor="none", edgecolor=c,
                               linewidths=0.5, zorder=5, clip_on=False)
            ax.set_ylim(lo_l, hi_l)
            ax.set_xlim(-0.6, len(dsets) - 0.4)
            ax.set_xticks(range(len(dsets)))
            ax.set_xticklabels([dataset_short(d, sp_of[d]) for d in dsets], rotation=55,
                               ha="right", rotation_mode="anchor", fontsize=5)
            for t, d in zip(ax.get_xticklabels(), dsets):
                t.set_color(dataset_color(d, sp_of[d]))
            if n_clip:
                ax.text(0.98, 0.99, f"{n_clip} off-scale", transform=ax.transAxes, ha="right",
                        va="top", fontsize=4.5, color="0.35",
                        bbox=dict(facecolor="white", edgecolor="none", pad=0.3, alpha=0.8), zorder=6)
        for ax in list(axes.flat)[n:]:
            ax.set_visible(False)
        h, l = species_dataset_legend(dsets, sp_of)
        legend_below(fig, axes, h, l)
        paths = save(fig, out_path)
        plt.close(fig)
    return paths


def pd_to_float(series):
    import pandas as _pd
    return _pd.to_numeric(series, errors="coerce").to_numpy(float)
