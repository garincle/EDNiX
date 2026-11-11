"""
Figure 1 a: atlas levels (Lvl2, Lvl3, Lvl4) and cortical thickness on the template surface of each species,
rendered from the Workbench scenes stored in the atlas library (atlases/install_figure_scenes.py), then laid out
with the sunbursts b (cortical) and c (non-cortical) and the abbreviation table (library/atlas_sunburst.py).

    render_species(library, out_dir, variant)     wb_command -show-scene of every <Species>_Figure1[_white].scene;
                                                  falls back to the reference PNG of the library if rendering fails
    draw_panel_a(images, out_no_ext, background)  5 x 2 species grid, group brackets (Primate, Ungulate, Rodents,
                                                  Carnivoran), binomial names in italics
    figure1(parts_dir, out_no_ext)                full Figure 1 (a | b over c, abbreviations below), white background;
                                                  also writes panel a alone on white and on black

Run alone (renders + assembles, sunbursts included):
    python -m Exemples.Study_EDNiX.analysis.library.figure1_surfaces
"""
import os
import shutil
import subprocess
from os.path import join as opj

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

from Exemples.Study_EDNiX.analysis import paper_data as pdat
from Exemples.Study_EDNiX.analysis import paper_style as ps

# (library folder, binomial name), in the reading order of panel a (row-major, 2 columns)
SPECIES = [('Human', 'Homo sapiens'), ('Chimpanzee', 'Pan troglodytes'),
           ('Macaque', 'Macaca mulatta'), ('Marmoset', 'Callithrix jacchus'),
           ('Mouselemur', 'Microcebus murinus'), ('Pig', 'Sus scrofa'),
           ('Rat', 'Rattus norvegicus'), ('Mouse', 'Mus musculus'),
           ('Dog', 'Canis familiaris'), ('Cat', 'Felis catus')]
NCOL = 2
# group label -> grid cells (row, col) it spans
GROUPS = [('Primate', [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0)]), ('Ungulate', [(2, 1)]),
          ('Rodents', [(3, 0), (3, 1)]), ('Carnivoran', [(4, 0), (4, 1)])]
# labels of the 2 x 2 views of every species image (Workbench tiles), drawn on the first species only
VIEW_LABELS = ['Lvl2', 'Lvl3', 'Lvl4', 'Thickness']
SCENE = {'white': 'Figure1_white', 'black': 'Figure1'}
SCENE_NUMBER = 1              # scene of the .scene file to render (name or 1-based number)
RENDER_SIZE = None            # None: window size saved in the scene (same framing as the hand export); or (w, h) px
OUT_DIR = opj(pdat.PAPER_DIR, 'Figure_1', 'surfaces')


def species_scene_dirs(library=pdat.ATLAS_LIB):
    """{species: <library>/atlas/.../<Species>/EDNiX/scenes}."""
    out = {}
    for dirpath, dirnames, _ in os.walk(opj(library, 'atlas')):
        if os.path.basename(dirpath) == 'freesurfer':
            dirnames[:] = []
        elif os.path.basename(dirpath) == 'scenes' and os.path.basename(os.path.dirname(dirpath)) == 'EDNiX':
            out[os.path.basename(os.path.dirname(os.path.dirname(dirpath)))] = dirpath
            dirnames[:] = []
    return out


def workbench(library=pdat.ATLAS_LIB, bind=pdat.PAPER_DIR):
    """wb_command prefix (EDNiX Workbench container, headless: no vglrun)."""
    from Tools import Load_EDNiX_requirement
    main_path = os.path.dirname(os.path.normpath(library)) + '/'
    sing_wb = Load_EDNiX_requirement.load_requirement(main_path, library, bind, 'yes')[4]
    return sing_wb.replace('vglrun ', '')


def render_species(library=pdat.ATLAS_LIB, out_dir=OUT_DIR, variant='white', sing_wb=None):
    """{species: png} rendered from <Species>_<SCENE[variant]>.scene; reference PNG when rendering fails."""
    os.makedirs(out_dir, exist_ok=True)
    dirs = species_scene_dirs(library)
    if sing_wb is None:
        sing_wb = workbench(library, os.path.dirname(os.path.normpath(out_dir)))
    out = {}
    for sp, _ in SPECIES:
        d = dirs.get(sp)
        if d is None:
            print(f'  [figure1a] {sp}: no scenes/ folder in the library')
            continue
        scene = opj(d, f'{sp}_{SCENE[variant]}.scene')
        ref = opj(d, f'{sp}_{SCENE[variant]}_reference.png')
        png = opj(out_dir, f'{sp}_{SCENE[variant]}.png')
        if os.path.isfile(scene):
            w, h = RENDER_SIZE or (1600, 1200)
            cmd = (f'{sing_wb}wb_command -show-scene "{scene}" "{SCENE_NUMBER}" "{png}" {w} {h}'
                   + ('' if RENDER_SIZE else ' -use-window-size'))
            r = subprocess.run(cmd, shell=True, cwd=d, capture_output=True, text=True)
            if r.returncode == 0 and os.path.isfile(png):
                out[sp] = png
                continue
            print(f'  [figure1a] {sp}: show-scene failed ({r.returncode}): {(r.stderr or r.stdout)[-400:]}')
        if os.path.isfile(ref):
            shutil.copy2(ref, png)
            out[sp] = png
            print(f'  [figure1a] {sp}: reference PNG used ({ref})')
        else:
            print(f'  [figure1a] {sp}: no scene and no reference PNG')
    return out


def autocrop(im, tol=0.03, margin=4):
    """Crop the uniform background (colour of the top-left pixel) of an RGB(A) image."""
    rgb = im[..., :3].astype(float)
    if rgb.max() > 1.5:
        rgb = rgb / 255.0
    keep = np.abs(rgb - rgb[0, 0]).max(axis=2) > tol
    if im.shape[-1] == 4:
        keep &= im[..., 3] > 0.02 * (255 if im.max() > 1.5 else 1)
    ys, xs = np.where(keep)
    if ys.size == 0:
        return im
    y0, y1 = max(ys.min() - margin, 0), min(ys.max() + margin + 1, im.shape[0])
    x0, x1 = max(xs.min() - margin, 0), min(xs.max() + margin + 1, im.shape[1])
    return im[y0:y1, x0:x1]


def _draw_grid(fig, rect, images, background='white', title_pt=7, group_pt=8, view_pt=5):
    """Species grid in fig coordinates rect = (x, y, w, h); returns nothing."""
    fg = 'black' if background == 'white' else 'white'
    x0, y0, w, h = rect
    nrow = int(np.ceil(len(SPECIES) / NCOL))
    gutter = 0.06 * w                       # band left of EVERY column for the group labels
    colw, rowh = (w - NCOL * gutter) / NCOL, h / nrow
    title_h = 0.16 * rowh
    for k, (sp, latin) in enumerate(SPECIES):
        r, c = divmod(k, NCOL)
        cx, cy = x0 + (c + 1) * gutter + c * colw, y0 + h - (r + 1) * rowh
        fig.text(cx + 0.02 * colw, cy + rowh - 0.02 * rowh, latin, fontsize=title_pt, style='italic',
                 color=fg, va='top', ha='left')
        ax = fig.add_axes([cx + 0.06 * colw, cy, 0.92 * colw, rowh - title_h])
        ax.axis('off')
        ax.set_facecolor(background)
        if sp in images and os.path.isfile(images[sp]):
            ax.imshow(autocrop(mpimg.imread(images[sp])), interpolation='lanczos')
            if k == 0 and VIEW_LABELS:
                for i, lab in enumerate(VIEW_LABELS):      # 2 x 2 Workbench tiles
                    ax.text(0.5 * (i % 2) + 0.48, 0.5 - 0.5 * (i // 2) + 0.02, lab, transform=ax.transAxes,
                            fontsize=view_pt, color=fg, ha='right', va='bottom')
        else:
            ax.text(0.5, 0.5, f'{sp}: missing', transform=ax.transAxes, ha='center', va='center',
                    fontsize=6, color='0.5')
    for name, cells in GROUPS:
        rows = [r for r, _ in cells]
        col = min(c for _, c in cells)
        top, bot = y0 + h - min(rows) * rowh - title_h, y0 + h - (max(rows) + 1) * rowh
        lx = x0 + col * (gutter + colw) + 0.4 * gutter          # inside the gutter of its own column
        fig.add_artist(plt.Line2D([lx + 0.4 * gutter] * 2, [bot + 0.03 * rowh, top], transform=fig.transFigure,
                                  color=fg, lw=0.8))
        fig.text(lx, 0.5 * (top + bot), name, rotation=90, fontsize=group_pt, fontweight='bold', color=fg,
                 ha='center', va='center')


def draw_panel_a(images, out_no_ext, background='white', width=ps.WIDTH_2COL * 0.58, height=None):
    """Panel a alone (.pdf + .png 600 dpi)."""
    height = height or width * 1.2
    with plt.rc_context(ps.PAPER_RC):
        fig = plt.figure(figsize=(width, height), facecolor=background)
        _draw_grid(fig, (0.0, 0.0, 1.0, 1.0), images, background)
        os.makedirs(os.path.dirname(out_no_ext), exist_ok=True)
        for ext in ('pdf', 'png'):
            fig.savefig(f'{out_no_ext}.{ext}', dpi=600 if ext == 'png' else None, facecolor=background,
                        bbox_inches='tight', pad_inches=0.02)
        plt.close(fig)
    print(f'  [figure1a] {out_no_ext}.pdf|png')
    return out_no_ext


def figure1(parts_dir, out_no_ext, library=pdat.ATLAS_LIB, sing_wb=None):
    """Full Figure 1, white background: a (species grid) | b cortical over c non-cortical sunburst, abbreviations
    below. parts_dir holds ednix_sunburst_{cortical,subcortical,legend}.png (atlas_sunburst.write_sunburst_figures)."""
    surf_dir = opj(os.path.dirname(out_no_ext), 'surfaces')
    images = {}
    for variant in ('white', 'black'):
        images[variant] = render_species(library, surf_dir, variant, sing_wb)
        draw_panel_a(images[variant], opj(surf_dir, f'Figure_1a_{variant}'), variant)
    width = ps.WIDTH_2COL
    height = width * 842 / 941                       # proportions of the current Figure 1
    leg_h, a_w = 0.16, 0.60
    with plt.rc_context(ps.PAPER_RC):
        fig = plt.figure(figsize=(width, height), facecolor='white')
        _draw_grid(fig, (0.0, leg_h + 0.01, a_w, 1.0 - leg_h - 0.01), images['white'], 'white')
        for letter, name, rect in (('b', 'cortical', (a_w, 0.585, 1 - a_w, 0.415)),
                                   ('c', 'subcortical', (a_w, leg_h + 0.01, 1 - a_w, 0.415)),
                                   (None, 'legend', (0.0, 0.0, 1.0, leg_h))):
            ax = fig.add_axes(rect)
            ax.axis('off')
            p = opj(parts_dir, f'ednix_sunburst_{name}.png')
            if os.path.isfile(p):
                ax.imshow(ps._autocrop(mpimg.imread(p)), interpolation='lanczos')
            if letter:
                fig.text(rect[0], rect[1] + rect[3], letter, fontsize=ps.PANEL_LABEL_SIZE, fontweight='bold',
                         va='top', ha='left')
        fig.text(0.0, 1.0, 'a', fontsize=ps.PANEL_LABEL_SIZE, fontweight='bold', va='top', ha='left')
        paths = ps.save(fig, out_no_ext)
        plt.close(fig)
    return paths


if __name__ == '__main__':
    from Exemples.Study_EDNiX.analysis.library import atlas_sunburst as asb
    base = opj(pdat.PAPER_DIR, 'Figure_1', 'Figure_1')
    parts = opj(pdat.PAPER_DIR, 'Figure_1', 'parts')
    asb.write_sunburst_figures(out_dir=parts, paper_out=opj(parts, 'Figure_1bc_sunburst'))
    figure1(parts, base)
