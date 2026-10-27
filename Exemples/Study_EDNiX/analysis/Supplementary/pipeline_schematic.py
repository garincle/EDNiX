"""
pipeline_schematic.py -- drawing engine of the EDNiX pipeline schematics (Supplementary Figs 1-3).

A schematic is drawn on one canvas in millimetres (183 mm wide, Nature double column):

  strip   a row of image frames, one per processing stage, laid on coloured bands that name the
          coordinate space the image lives in (scanner, native AC-PC, functional, PET, CT, template).
          Arrows between frames say what links them: an ordinary processing step (grey), a
          registration that is ESTIMATED (thick vermilion, transform written above), or transforms
          that are only APPLIED (dashed purple).
  cards   the steps of the module as numbered cards (step number = Skip_step of the launcher),
          grouped by processing block, coloured by the kind of operation.
  graph   small node-link diagrams of registration chains.

Images come from a processed EDNiX subject when the files are found (Thumb.from_files), else
from a synthetic brain (SyntheticBrain, built from the MNI ICBM152 2009 template shipped with
nilearn, or from an analytic phantom when nilearn is missing). Synthetic frames carry the word
'synthetic' in their corner, so that an illustration is never mistaken for data.
"""
import os
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle
from scipy import ndimage as ndi

from Exemples.Study_EDNiX.analysis import paper_style as ps

MM = ps.MM
W = 183.0                      # figure width (mm), Nature double column

# ── colours ──────────────────────────────────────────────────────────────────
# coordinate spaces: light fill of the band, darker colour of its label
SPACES = {
    'scanner':  dict(fill='#F0F0F0', ink='#555555', name='Scanner space'),
    'native':   dict(fill='#E5EEF8', ink='#1F4E79', name='Native AC-PC space'),
    'func':     dict(fill='#E3F3EF', ink='#00695C', name='Functional space'),
    'funcacpc': dict(fill='#E3EEF3', ink='#1B5E6B', name='Functional AC-PC space'),
    'pet':      dict(fill='#FCEEDC', ink='#9A5300', name='PET space'),
    'ct':       dict(fill='#EFE8F4', ink='#5E3A87', name='CT space'),
    'study':    dict(fill='#F5F0E1', ink='#7A6A2E', name='Study template'),
    'template': dict(fill='#ECE7F5', ink='#45358A', name='Reference template'),
}
# kinds of operation (Okabe & Ito colour-blind safe palette)
OPS = {
    'prepro':   dict(color='#0072B2', name='Intensity / geometry'),
    'mask':     dict(color='#009E73', name='Brain extraction / masks'),
    'reg':      dict(color='#D55E00', name='Registration (estimated)'),
    'apply':    dict(color='#CC79A7', name='Transforms applied'),
    'check':    dict(color='#E69F00', name='Visual check (optional)'),
    'analysis': dict(color='#56B4E9', name='Signal / analysis'),
}
FLOW = '#7F7F7F'
TEXT = '#222222'
SUBTLE = '#666666'

FS_TITLE = 6.5      # block headers
FS_TEXT = 5.6       # card titles / text (detail lines FS_TEXT - 0.5)
FS_SMALL = 5.0      # frame titles, captions, transform labels (Nature minimum: 5 pt)


def label_lut(n=40, seed=3):
    """Qualitative colours for atlas labels (tab20 repeated, shuffled with a fixed seed)."""
    base = list(plt.get_cmap('tab20').colors) + list(plt.get_cmap('tab20b').colors)
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(base))
    cols = [base[i] for i in idx]
    return [cols[i % len(cols)] for i in range(n)]


# ── canvas in millimetres ────────────────────────────────────────────────────
class Canvas:
    """Figure of W x H mm with a full-size axes in mm, y downwards (0 = top)."""

    def __init__(self, height_mm, width_mm=W):
        self.w, self.h = width_mm, height_mm
        self._rc = plt.rc_context(ps.PAPER_RC)
        self._rc.__enter__()
        self.fig = plt.figure(figsize=(width_mm * MM, height_mm * MM))
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set_xlim(0, width_mm)
        self.ax.set_ylim(height_mm, 0)
        self.ax.axis('off')

    # primitives -----------------------------------------------------------------
    def rect(self, x, y, w, h, fc='none', ec='none', lw=0.5, r=1.2, ls='-', z=1, alpha=1.0):
        p = FancyBboxPatch((x, y), w, h, boxstyle=f'round,pad=0,rounding_size={r}', fc=fc, ec=ec, lw=lw,
                           ls=ls, zorder=z, alpha=alpha, mutation_aspect=1)
        self.ax.add_patch(p)
        return p

    def text(self, x, y, s, size=FS_TEXT, weight='normal', ha='left', va='top', color=TEXT, z=10, **kw):
        return self.ax.text(x, y, s, fontsize=size, fontweight=weight, ha=ha, va=va, color=color, zorder=z,
                            linespacing=1.15, **kw)

    def arrow(self, x0, y0, x1, y1, color=FLOW, lw=0.6, ls='-', head=3.0, rad=0.0, z=6, shrink=0.0):
        a = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=f'-|>,head_length={head * 0.35},head_width={head * 0.2}',
                            connectionstyle=f'arc3,rad={rad}', color=color, lw=lw, ls=ls, zorder=z,
                            shrinkA=shrink, shrinkB=shrink, mutation_scale=1, capstyle='round')
        self.ax.add_patch(a)
        return a

    def badge(self, x, y, s, color, r=1.55, size=4.6, z=12):
        """Step number in a filled circle centred on (x, y); a pill for longer labels ('itk_1', '10-15'),
        which then starts at x - r. Returns the right edge (mm)."""
        s = str(s)
        if len(s) <= 2:
            self.ax.add_patch(Circle((x, y), r, fc=color, ec='white', lw=0.4, zorder=z))
            self.text(x, y + 0.05, s, size=size, weight='bold', ha='center', va='center', color='white', z=z + 1)
            return x + r
        w = 0.62 * size * 0.3528 * len(s) + 1.2
        self.ax.add_patch(FancyBboxPatch((x - r, y - r), w, 2 * r, boxstyle=f'round,pad=0,rounding_size={r}',
                                         fc=color, ec='white', lw=0.4, zorder=z))
        self.text(x - r + w / 2, y + 0.05, s, size=size - 0.4, weight='bold', ha='center', va='center',
                  color='white', z=z + 1)
        return x - r + w

    def panel(self, letter, x, y, title=None):
        self.text(x, y, letter, size=ps.PANEL_LABEL_SIZE, weight='bold', va='top')
        if title:
            self.text(x + 3.6, y + 0.45, title, size=FS_TITLE + 0.5, weight='bold', va='top')

    def axes_mm(self, x, y, w, h):
        """Matplotlib axes covering the mm box (x, y, w, h) of the canvas."""
        return self.fig.add_axes([x / self.w, 1 - (y + h) / self.h, w / self.w, h / self.h])

    def save(self, path_no_ext):
        out = ps.save(self.fig, path_no_ext)
        plt.close(self.fig)
        self._rc.__exit__(None, None, None)
        return out

    # composite elements ---------------------------------------------------------
    def space_band(self, x, y, w, h, space, label=True, dashed=False, label_text=None):
        s = SPACES[space]
        self.rect(x, y, w, h, fc=s['fill'], ec=s['ink'] if dashed else 'none', lw=0.4,
                  ls=(0, (2, 1.5)) if dashed else '-', r=1.6, z=0)
        if label:
            self.text(x + w / 2, y + 1.0, label_text or s['name'], size=FS_SMALL + 0.6, weight='bold',
                      ha='center', va='top', color=s['ink'], z=2)

    def transition(self, x0, x1, y, kind='flow', label=None, label_dy=-1.2, rad=0.0, sub=None):
        """Arrow between two frames. kind: 'flow' (processing), 'reg' (registration estimated),
        'apply' (transforms applied)."""
        if kind == 'reg':
            self.arrow(x0, y, x1, y, color=OPS['reg']['color'], lw=1.1, head=4.2, rad=rad)
        elif kind == 'apply':
            self.arrow(x0, y, x1, y, color=OPS['apply']['color'], lw=0.8, ls=(0, (2.2, 1.2)), head=3.6, rad=rad)
        else:
            self.arrow(x0, y, x1, y, color=FLOW, lw=0.6, head=3.0, rad=rad)
        xm = (x0 + x1) / 2
        if label:
            col = OPS['reg']['color'] if kind == 'reg' else OPS['apply']['color'] if kind == 'apply' else SUBTLE
            self.text(xm, y + label_dy, label, size=FS_SMALL, ha='center', va='bottom', color=col,
                      weight='bold' if kind == 'reg' else 'normal')
        if sub:
            self.text(xm, y + 1.2, sub, size=FS_SMALL, ha='center', va='top', color=SUBTLE)

    def card(self, x, y, w, step, title, lines, op, h=None, dashed=False, tool=None):
        """Step card: coloured left edge (kind of operation), step badge, bold title, detail lines,
        optional tool line in grey italics. Returns the card height (mm)."""
        n = len(lines) + (1 if tool else 0)
        LS = 2.3
        h = h or 5.0 + LS * n
        col = OPS[op]['color']
        self.rect(x, y, w, h, fc='white', ec='#BFBFBF' if not dashed else col, lw=0.35,
                  ls=(0, (2, 1.2)) if dashed else '-', r=0.9, z=3)
        self.ax.add_patch(Rectangle((x, y + 0.25), 0.9, h - 0.5, fc=col, ec='none', zorder=4))
        tx = x + 2.0
        if step is not None and step != '':
            tx = self.badge(x + 3.4, y + 2.35, step, col, r=1.5, size=4.8) + 1.0
        self.text(tx, y + 0.95, title, size=FS_TEXT + 0.3, weight='bold', z=5)
        yy = y + 3.9
        for ln in lines:
            self.text(x + 2.0, yy, ln, size=FS_TEXT - 0.5, color=TEXT, z=5)
            yy += LS
        if tool:
            self.text(x + 2.0, yy, tool, size=FS_TEXT - 0.6, color=SUBTLE, style='italic', z=5)
        return h

    def elbow(self, x, y_from, y_to, dx, color, lw=0.7, ls=(0, (2.2, 1.2)), head=3.2):
        """Arrow leaving (x, y_from) to the right, running vertically at x + dx, entering (x, y_to)."""
        xr = x + dx
        self.ax.plot([x, xr, xr], [y_from, y_from, y_to], color=color, lw=lw, ls=ls, zorder=6,
                     solid_capstyle='round', dash_capstyle='round')
        self.arrow(xr, y_to, x + 0.3, y_to, color=color, lw=lw, head=head)
        return xr

    def trim(self, used):
        """Shrink the canvas height to `used` mm, keeping every element at its mm position."""
        H = self.h
        boxes = [(ax, ax.get_position(original=True).bounds) for ax in self.fig.axes]
        self.fig.set_size_inches(self.w * MM, used * MM)
        for ax, (x0, y0, w, h) in boxes:
            if ax is self.ax:
                ax.set_position([0, 0, 1, 1])
                ax.set_ylim(used, 0)
                continue
            top_mm, h_mm = (1 - y0 - h) * H, h * H
            ax.set_position([x0, 1 - (top_mm + h_mm) / used, w, h_mm / used])
        self.h = used

    def block_header(self, x, y, w, title, sub=None, color=TEXT):
        self.text(x, y, title, size=FS_TITLE, weight='bold', color=color)
        if sub:
            self.text(x, y + 2.9, sub, size=FS_SMALL, color=SUBTLE, style='italic')
        self.ax.plot([x, x + w], [y + (5.4 if sub else 2.9)] * 2, color='#BFBFBF', lw=0.4, zorder=2)

    def ops_legend(self, x, y, keys=None):
        """One row: coloured ticks + names of the kinds of operation, then the arrow conventions."""
        keys = keys or list(OPS)
        xx = x
        for k in keys:
            self.ax.add_patch(Rectangle((xx, y - 0.9), 0.9, 1.8, fc=OPS[k]['color'], ec='none', zorder=4))
            self.text(xx + 1.6, y, OPS[k]['name'], size=FS_SMALL, va='center')
            xx += 3.4 + 1.02 * len(OPS[k]['name'])
        return xx

    def arrow_legend(self, x, y):
        items = [('flow', 'processing'), ('reg', 'registration estimated'), ('apply', 'transforms applied')]
        xx = x
        for kind, name in items:
            self.transition(xx, xx + 6, y, kind=kind)
            self.text(xx + 7.0, y, name, size=FS_SMALL, va='center')
            xx += 9.5 + 1.02 * len(name)
        return xx


# ── frames (thumbnails) ──────────────────────────────────────────────────────
class Thumb:
    """One image frame: 2-D image + optional contours and label map, or a custom plot function."""

    def __init__(self, img=None, cmap='gray', vmin=None, vmax=None, contours=(), labels=None, label_alpha=0.75,
                 aspect=1.0, plot=None, synthetic=False, bg='black', overlay=None,
                 overlay_cmap='hot', overlay_alpha=0.85, overlay_thr=None, lines=()):
        self.img, self.cmap, self.vmin, self.vmax = img, cmap, vmin, vmax
        self.contours = list(contours)          # [(mask2d, colour, lw, linestyle)]
        self.labels, self.label_alpha = labels, label_alpha
        self.aspect = aspect                     # pixel height / pixel width (mm)
        self.plot = plot                         # callable(ax) for non-image frames
        self.synthetic = synthetic
        self.bg = bg
        self.overlay, self.overlay_cmap, self.overlay_alpha, self.overlay_thr = overlay, overlay_cmap, overlay_alpha, overlay_thr
        self.lines = list(lines)                 # [((x0, y0), (x1, y1), colour, label)] in pixels

    def draw(self, ax):
        if self.plot is not None:
            self.plot(ax)
            if self.synthetic:
                _synthetic_tag(ax, dark=False)
            return
        ax.set_facecolor(self.bg)
        img = self.img
        vmin = self.vmin if self.vmin is not None else np.nanpercentile(img, 1)
        vmax = self.vmax if self.vmax is not None else np.nanpercentile(img, 99.5)
        ax.imshow(img, cmap=self.cmap, vmin=vmin, vmax=vmax, interpolation='bilinear', aspect=self.aspect,
                  origin='upper')
        if self.overlay is not None:
            ov = np.ma.masked_where(~np.isfinite(self.overlay) | (self.overlay < (self.overlay_thr or 0)),
                                    self.overlay)
            ax.imshow(ov, cmap=self.overlay_cmap, alpha=self.overlay_alpha, interpolation='bilinear',
                      aspect=self.aspect, vmin=self.overlay_thr, vmax=np.nanmax(self.overlay))
        if self.labels is not None:
            lab = self.labels.astype(int)
            vals = [v for v in np.unique(lab) if v > 0]
            if vals:
                lut = label_lut(max(vals) + 1)
                rgba = np.zeros(lab.shape + (4,))
                for v in vals:
                    rgba[lab == v, :3] = lut[v][:3] if len(lut[v]) >= 3 else lut[v]
                    rgba[lab == v, 3] = self.label_alpha
                ax.imshow(rgba, interpolation='nearest', aspect=self.aspect)
                # thin black borders between labels
                edges = (ndi.grey_dilation(lab, size=2) != ndi.grey_erosion(lab, size=2)) & (lab > 0)
                e = np.zeros(lab.shape + (4,))
                e[edges, 3] = 0.35
                ax.imshow(e, interpolation='nearest', aspect=self.aspect)
        for c in self.contours:
            m, col, lw = c[0], c[1], c[2]
            ls = c[3] if len(c) > 3 else '-'
            if m is not None and np.any(m > 0.5) and np.any(m < 0.5):
                ax.contour(m.astype(float), levels=[0.5], colors=[col], linewidths=lw, linestyles=[ls])
        for (x0, y0), (x1, y1), col, lab in self.lines:
            ax.plot([x0, x1], [y0, y1], color=col, lw=0.5, ls=(0, (2.5, 1.5)))
            if lab:
                ax.text(x0, y0 + 4, lab, color=col, fontsize=4.2, ha='left', va='top')
        ax.set_xlim(-0.5, img.shape[1] - 0.5)
        ax.set_ylim(img.shape[0] - 0.5, -0.5)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        if self.synthetic:
            _synthetic_tag(ax)


def _synthetic_tag(ax, dark=True):
    ax.text(0.97, 0.03, 'synthetic', transform=ax.transAxes, ha='right', va='bottom', fontsize=3.4,
            color='#BBBBBB' if dark else '#999999', style='italic', zorder=20)


def draw_frame(cv, x, y, w, h, thumb, title=None, step=None, op=None, caption=None, frame_color=None):
    """Frame at (x, y) of size w x h mm: title above (with step badge), image, caption below."""
    if title:
        tx = x
        if step is not None:
            tx = cv.badge(x + 1.5, y - 2.1, step, OPS[op]['color'] if op else FLOW, r=1.45, size=4.6) + 0.9
        cv.text(tx, y - 2.1, title, size=FS_SMALL + 0.4, weight='bold', va='center')
    ax = cv.axes_mm(x, y, w, h)
    thumb.draw(ax)
    if frame_color:
        cv.rect(x - 0.25, y - 0.25, w + 0.5, h + 0.5, ec=frame_color, lw=0.6, r=0.4, z=8)
    if caption:
        cv.text(x + w / 2, y + h + 0.9, caption, size=FS_SMALL, ha='center', va='top', color=SUBTLE)
    return ax


def fit_crop(img, aspect_target, support=None, pad=0.06):
    """Crop a 2-D image (and return the slices) around `support` (bool mask of the content) so that
    its height / width ratio equals aspect_target, with `pad` relative margin."""
    if support is None:
        support = np.isfinite(img) & (img > np.nanpercentile(img, 5))
    ys, xs = np.where(support)
    if ys.size == 0:
        return slice(None), slice(None)
    y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
    h, w = y1 - y0, x1 - x0
    h, w = h * (1 + 2 * pad), w * (1 + 2 * pad)
    cy, cx = (y0 + y1) / 2, (x0 + x1) / 2
    if h / w > aspect_target:
        w = h / aspect_target
    else:
        h = w * aspect_target
    sy = slice(int(round(cy - h / 2)), int(round(cy + h / 2)))
    sx = slice(int(round(cx - w / 2)), int(round(cx + w / 2)))
    return sy, sx


def crop_pad(a, sy, sx, fill=0):
    """a[sy, sx] with out-of-range rows / columns filled with `fill`."""
    H, Wd = a.shape[:2]
    out = np.full((sy.stop - sy.start, sx.stop - sx.start) + a.shape[2:], fill, dtype=float)
    ys0, ys1 = max(sy.start, 0), min(sy.stop, H)
    xs0, xs1 = max(sx.start, 0), min(sx.stop, Wd)
    if ys1 > ys0 and xs1 > xs0:
        out[ys0 - sy.start:ys1 - sy.start, xs0 - sx.start:xs1 - sx.start] = a[ys0:ys1, xs0:xs1]
    return out


# ── strip of frames on space bands ───────────────────────────────────────────
def draw_strip(cv, frames, x0, y0, width, frame_w, frame_h, band_pad=1.6, title_h=4.6, label_band=True):
    """frames: list of dicts with keys
         thumb, title, step, op, space, caption, link (arrow INTO this frame: dict(kind, label, sub)),
         band_label (optional text of the space band, on the first frame of a band).
    Frames are spread evenly over `width`. Consecutive frames in the same space share one band.
    Returns the list of frame boxes (x, y, w, h)."""
    n = len(frames)
    gap = (width - n * frame_w) / (n - 1) if n > 1 else 0
    xs = [x0 + i * (frame_w + gap) for i in range(n)]
    band_top = y0
    fy = y0 + title_h + 3.2                                 # frame top (space label + title above)
    band_h = (fy - band_top) + frame_h + band_pad + 3.2     # caption room below
    # bands
    i = 0
    while i < n:
        j = i
        while j + 1 < n and frames[j + 1]['space'] == frames[i]['space']:
            j += 1
        bx = xs[i] - band_pad
        bw = xs[j] + frame_w - xs[i] + 2 * band_pad
        cv.space_band(bx, band_top, bw, band_h, frames[i]['space'], label=label_band,
                      dashed=frames[i].get('dashed', False), label_text=frames[i].get('band_label'))
        i = j + 1
    boxes = []
    for k, (x, f) in enumerate(zip(xs, frames)):
        draw_frame(cv, x, fy, frame_w, frame_h, f['thumb'], f.get('title'), f.get('step'), f.get('op'),
                   f.get('caption'), f.get('frame_color'))
        boxes.append((x, fy, frame_w, frame_h))
        link = f.get('link')
        if k > 0 and link:
            xa, xb = xs[k - 1] + frame_w + 0.35, x - 0.35
            cv.transition(xa, xb, fy + frame_h / 2, kind=link.get('kind', 'flow'), label=link.get('label'),
                          sub=link.get('sub'), label_dy=link.get('label_dy', -1.0))
    return boxes


# ── data: real files ─────────────────────────────────────────────────────────
def find_session(bids_dir, patterns, prefer=None):
    """First sub-*/ses-* (or sub-*) folder of `bids_dir` where every glob of `patterns` ({key: relative
    glob}) matches a file. Returns ({key: path}, folder) or (None, None). `prefer`: subject id tried first."""
    if not bids_dir or not os.path.isdir(bids_dir):
        return None, None
    subs = sorted(glob.glob(os.path.join(bids_dir, 'sub-*')))
    if prefer:
        subs = sorted(subs, key=lambda s: os.path.basename(s) != prefer)
    for sub in subs:
        sess = sorted(glob.glob(os.path.join(sub, 'ses-*'))) or [sub]
        for ses in sess:
            found = {}
            for k, pat in patterns.items():
                hits = sorted(glob.glob(os.path.join(ses, pat)))
                if not hits:
                    break
                found[k] = hits[0]
            if len(found) == len(patterns):
                return found, ses
    return None, None


def missing_report(bids_dir, patterns):
    """For the console: patterns that match nothing anywhere in the dataset."""
    out = []
    for k, pat in patterns.items():
        if not glob.glob(os.path.join(bids_dir, 'sub-*', 'ses-*', pat)) and \
                not glob.glob(os.path.join(bids_dir, 'sub-*', pat)):
            out.append(f'{k}: {pat}')
    return out


class Volume:
    """NIfTI image in RAS+ orientation (nibabel.as_closest_canonical)."""

    def __init__(self, path, frame=None):
        import nibabel as nib
        img = nib.as_closest_canonical(nib.load(path))
        data = np.asanyarray(img.dataobj)
        if data.ndim == 4:
            if frame == 'mean':
                data = np.asarray(data[..., :min(data.shape[3], 20)], dtype=np.float32).mean(-1)
            else:
                data = data[..., frame or 0]
        self.data = np.asarray(data, dtype=np.float32)
        self.affine = img.affine
        self.zooms = img.header.get_zooms()[:3]
        self.path = path
        self.img = nib.Nifti1Image(self.data, self.affine)

    def resample_like(self, other, order=0):
        """This image on the grid of `other` (labels / masks: order 0)."""
        from nibabel.processing import resample_from_to
        r = resample_from_to(self.img, (other.data.shape, other.affine), order=order)
        v = Volume.__new__(Volume)
        v.data = np.asarray(r.dataobj, dtype=np.float32)
        v.affine, v.zooms, v.path, v.img = other.affine, other.zooms, self.path, r
        return v

    def sag_index(self, support=None, offset_mm=0.0):
        """Sagittal index at the centre of mass of `support` (or of the image) + offset (mm, + = right)."""
        s = support if support is not None else (self.data > np.percentile(self.data[self.data > 0], 30)
                                                 if np.any(self.data > 0) else np.ones_like(self.data, bool))
        cx = ndi.center_of_mass(s.astype(float))[0]
        return int(np.clip(round(cx + offset_mm / self.zooms[0]), 0, self.data.shape[0] - 1))

    def sagittal(self, i):
        """2-D sagittal slice: rows = superior -> inferior, columns = posterior -> anterior."""
        return np.rot90(self.data[i, :, :])

    @property
    def sag_aspect(self):
        return float(self.zooms[2] / self.zooms[1])


# ── data: synthetic brain ────────────────────────────────────────────────────
class SyntheticBrain:
    """2-D parasagittal illustration images, 1 mm pixels, derived from the MNI ICBM152 2009c template and
    its grey / white matter probability maps (nilearn package data), or from an analytic phantom.
    Everything is generated in 2-D in the sagittal plane: scanner obliquity is an in-plane pitch,
    EPI distortion an in-plane shift along the phase-encoding (antero-posterior) axis."""

    def __init__(self, x_offset_mm=-14.0, seed=7):
        self.rng = np.random.default_rng(seed)
        try:
            from nilearn import datasets
            t1 = np.asarray(datasets.load_mni152_template(resolution=1).dataobj, dtype=np.float32)
            gm = np.asarray(datasets.load_mni152_gm_template(resolution=1).dataobj, dtype=np.float32)
            wm = np.asarray(datasets.load_mni152_wm_template(resolution=1).dataobj, dtype=np.float32)
            i = int(round(t1.shape[0] / 2 + x_offset_mm))
            self.t1 = np.rot90(t1[i])
            self.gm, self.wm = np.rot90(gm[i]), np.rot90(wm[i])
            self.source = 'MNI ICBM152 2009c (nilearn)'
        except Exception:
            self.t1, self.gm, self.wm = _phantom()
            self.source = 'analytic phantom'
        self.t1 = self.t1 / np.percentile(self.t1[self.t1 > 0], 99)
        self.brain = ndi.binary_fill_holes(ndi.binary_closing(self.t1 > 0.08, iterations=2))
        self.brain = _largest(self.brain)
        self.csf = np.clip(self.brain * 1.0 - self.gm - self.wm, 0, 1)
        pad = 40
        for k in ('t1', 'gm', 'wm', 'csf'):
            setattr(self, k, np.pad(getattr(self, k), pad))
        self.brain = np.pad(self.brain, pad)
        self.shape = self.t1.shape
        # smooth envelope of the brain (the skull does not follow the sulci); skull 0-5 mm, scalp and fat 7-7+t mm (thicker over the vertex and occiput)
        yy, xx = np.mgrid[0:self.shape[0], 0:self.shape[1]]
        env = ndi.binary_closing(self.brain, structure=_disk(14), iterations=1)
        env = ndi.binary_fill_holes(env)
        env = ndi.gaussian_filter(env.astype(float), 3) > 0.5
        d_env = ndi.distance_transform_edt(~env)
        cy, cx = ndi.center_of_mass(self.brain)
        ang = np.arctan2(-(yy - cy), xx - cx)
        thick = 6.0 + 2.5 * np.sin(ang) + 1.5 * np.cos(2 * ang)
        self.d_out = d_env
        self.csf_rim = env & ~self.brain
        self.skull = (d_env > 0) & (d_env <= 5)
        self.scalp = (d_env > 5) & (d_env <= 5 + thick)
        self.head = env | (d_env <= 5 + thick)
        self.neck_soft = np.zeros(self.shape, bool)
        self._lobes = None

    # helpers ----------------------------------------------------------------
    def smooth(self, a, fwhm_mm):
        return ndi.gaussian_filter(a, fwhm_mm / 2.3548)

    def noise(self, a, sd):
        return a + self.rng.normal(0, sd, a.shape) * (a > -1e9)

    def rotate(self, a, deg, shift=(0, 0), order=1):
        r = ndi.rotate(a, deg, reshape=False, order=order, mode='constant', cval=float(np.min(a)) if order else 0)
        if shift != (0, 0):
            r = ndi.shift(r, shift, order=order, mode='constant', cval=float(np.min(a)) if order else 0)
        return r

    def acpc_line(self, deg=0.0, shift=(0, 0), half=0.42):
        """End points (x, y) of an illustrative AC-PC plane: horizontal line through the centre of the brain,
        then the same rotation / shift as SyntheticBrain.rotate."""
        cy, cx = ndi.center_of_mass(self.brain)
        cy += 0.06 * np.ptp(np.where(self.brain)[0])
        w = half * np.ptp(np.where(self.brain)[1])
        m = np.zeros(self.shape)
        m[int(cy) - 1:int(cy) + 2, int(cx - w):int(cx + w)] = 1
        if deg or shift != (0, 0):
            m = self.rotate(m, deg, shift)
        ys, xs = np.where(m > 0.5)
        i0, i1 = np.argmin(xs), np.argmax(xs)
        return (xs[i0], ys[i0]), (xs[i1], ys[i1])

    def bias(self, strength=0.45):
        yy, xx = np.mgrid[0:self.shape[0], 0:self.shape[1]]
        yy, xx = yy / self.shape[0], xx / self.shape[1]
        return 1 + strength * (0.9 * (xx - 0.5) - 0.7 * (yy - 0.5) + 0.6 * (xx - 0.5) * (yy - 0.5))

    # images -----------------------------------------------------------------
    def t1_head(self):
        """T1-weighted head: brain (template), dark CSF rim and skull, bright scalp with fine texture."""
        tex = self.smooth(self.rng.normal(0, 1, self.shape), 2.5)
        h = self.t1 * self.brain
        h = h + self.csf_rim * 0.10 + self.skull * (0.05 + 0.03 * tex) + self.scalp * (0.66 + 0.12 * tex)
        h = np.where(self.neck_soft, 0.40 + 0.08 * tex, h)
        return np.clip(self.smooth(h, 0.8), 0, None)

    def lobes(self):
        """Illustrative cortical parcellation of the grey matter: angular sectors around the brain centre
        (frontal, parietal, occipital, temporal...) + cerebellum. Labels 1-9."""
        if self._lobes is not None:
            return self._lobes
        cy, cx = ndi.center_of_mass(self.brain)
        yy, xx = np.mgrid[0:self.shape[0], 0:self.shape[1]]
        ang = np.degrees(np.arctan2(-(yy - cy), xx - cx))        # 0 = anterior, 90 = superior
        gm = (self.gm > 0.35) & self.brain
        lab = np.zeros(self.shape, int)
        sectors = [(-30, 20, 1), (20, 60, 2), (60, 95, 3), (95, 130, 4), (130, 165, 5), (165, 200, 6)]
        for a0, a1, v in sectors:
            sel = ((ang >= a0) & (ang < a1)) | ((ang + 360 >= a0) & (ang + 360 < a1))
            lab[gm & sel] = v
        lab[gm & (ang < -30) & (ang > -100)] = 7                  # temporal
        lab[gm & (ang <= -100) & (ang > -170)] = 8                 # cerebellum / occipito-temporal
        lab[(self.wm > 0.5) & self.brain & (yy > cy + 0.25 * (np.ptp(np.where(self.brain)[0])))] = 0
        self._lobes = lab
        return lab

    def epi(self, distort=0.0, dropout=0.0):
        """T2*-weighted EPI (CSF bright, grey > white), 3 mm pixels, optional susceptibility distortion
        (shift along the antero-posterior phase-encoding axis, mm) and signal loss near the orbitofrontal
        cortex."""
        e = (0.62 * self.gm + 0.42 * self.wm + 0.95 * self.csf) * self.brain + 0.10 * self.csf_rim
        e = self.smooth(e, 2.5)
        if distort or dropout:
            yy, xx = np.indices(self.shape, dtype=float)
            ys, xs = np.where(self.brain)
            # air-tissue interface: inferior frontal (ventral, anterior) region
            y_c, x_c = np.percentile(ys, 78), np.percentile(xs, 88)
            field = np.exp(-(((yy - y_c) / 16) ** 2 + ((xx - x_c) / 22) ** 2))
            e = e * (1 - dropout * field)
            e = ndi.map_coordinates(e, [yy, xx - distort * field], order=1)
        return block_mean(e, 3)

    def fdg(self, fwhm=6.5, extracerebral=0.14):
        """FDG-like uptake (GM 1, WM 0.35, CSF 0.05; Phelps et al. 1979), extracerebral uptake in the scalp,
        smoothed to the PET resolution, 2 mm pixels."""
        a = (self.gm * 1.0 + self.wm * 0.35 + self.csf * 0.05) * self.brain + extracerebral * self.scalp
        return block_mean(self.smooth(a, fwhm), 2)

    def synthetic_target(self, fwhm=6.5):
        """coreg_target='synthetic' of the PET pipeline: GM 1, WM 0.35, CSF 0.05 inside the brain mask, Gaussian
        of the PET FWHM (no extracerebral signal)."""
        a = np.where(self.brain, 0.05, 0.0)
        a = np.where(self.wm > 0.5, 0.35, a)
        a = np.where(self.gm > 0.5, 1.0, a)
        return self.smooth(a, fwhm)

    def ct(self):
        """CT in Hounsfield units: air -1000, soft tissue 40, skull 1000-1600 (diploe darker), brain 25-40, CSF 5."""
        tex = self.smooth(self.rng.normal(0, 1, self.shape), 2.0)
        h = np.full(self.shape, -1000.0)
        h[self.scalp] = 40 + 15 * tex[self.scalp]
        sk = self.skull
        dip = sk & (np.abs(self.d_out - 2.5) < 0.8)
        h[sk] = 1450 + 120 * tex[sk]
        h[dip] = 700
        h[self.csf_rim] = 8
        b = self.brain
        h[b] = 5 + 32 * self.gm[b] + 24 * self.wm[b] + 3 * tex[b]
        h[self.neck_soft] = 35 + 10 * tex[self.neck_soft]
        return h


def block_mean(a, k):
    """Downsample a 2-D image by block averaging (k x k pixels)."""
    H, Wd = (a.shape[0] // k) * k, (a.shape[1] // k) * k
    return a[:H, :Wd].reshape(H // k, k, Wd // k, k).mean(axis=(1, 3))


def upsample(a, k, order=0):
    return ndi.zoom(a, k, order=order)


def _disk(r):
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    return (yy ** 2 + xx ** 2) <= r ** 2


def _largest(m):
    lab, n = ndi.label(m)
    if n <= 1:
        return m
    sizes = ndi.sum(m, lab, range(1, n + 1))
    return lab == (1 + int(np.argmax(sizes)))


def _phantom(n=230):
    """Analytic sagittal 'brain' (ellipse + cerebellum + folded ribbon), when nilearn is not installed."""
    yy, xx = np.indices((n, n), dtype=float)
    cy, cx = n * 0.45, n * 0.5
    r = np.sqrt(((yy - cy) / 70) ** 2 + ((xx - cx) / 88) ** 2)
    cereb = np.sqrt(((yy - cy - 62) / 26) ** 2 + ((xx - cx + 52) / 34) ** 2) < 1
    brain = (r < 1) | cereb
    th = np.arctan2(yy - cy, xx - cx)
    ribbon = (r > 0.80 + 0.06 * np.sin(14 * th)) & brain
    gm = ribbon.astype(float) + cereb * 0.6
    wm = (brain & ~ribbon & ~cereb).astype(float)
    vent = np.sqrt(((yy - cy + 2) / 10) ** 2 + ((xx - cx) / 30) ** 2) < 1
    wm[vent] = 0
    t1 = 0.55 * gm + 0.85 * wm + 0.15 * vent
    return t1 * brain, gm * brain, wm * brain


# ── small helpers for plotted frames ─────────────────────────────────────────
def style_plot_axes(ax, bg='white'):
    ax.set_facecolor(bg)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_linewidth(0.35)
    ax.tick_params(length=1.2, width=0.35, labelsize=3.6, pad=0.8)


def edges_of(a, level):
    """Binary mask of a > level (for iso-intensity contours of a registration QC)."""
    return (a > level).astype(float)
