"""
Figure 6 -- Multimodal, multi-species imaging approaches to build an intersectional brain function profile.

 A  surface maps of the group-mean brain ratio of FDG uptake (bratio), human left, macaque right, same colour
    scale (Statistics/Group_PET: pet_group.run_group_mean + pet_surface.render_pet_surfaces, ribbon projection);
 B  cross-species comparison of glucose metabolism in the regions of LEVEL (bilateral), in bratio (bars
    drawn from 0, the values of the maps of A): bars and statistics of Plot_interspecies_regional_profiles.py;
 C  per-region multimodal profile (morphometry, FC, metabolism) on radar axes, one radar per species: each
    modality z-scored across the regions within each subject, mean +- SEM across subjects.
Regions (B, C) grouped by lobe in the order of the atlas sunburst of Figure 1 (Atlas_nomenclature.atlas_order),
radars counter-clockwise from the top as the sunburst, lobe colours of Figure 1 on the rim and under the bars.

Data, units and level conventions: multimodal_profiles.py. Colours of the species: those of the article figures
(paper_style.py); colours of the modalities: Okabe-Ito.
  figure6_stats()        every table and test of B and C (paper_stats.figure6)
  draw_figure6()         panel A renders + the figure + its legend (paper_figures.figure6)
Former multimodal_PET/Figure_XX_multimodal.py (Figure XX, now Figure 6, Clément 2026-10-27):
main() split into these two functions, computations unchanged.
"""
import os
import sys
import glob
import textwrap

opj = os.path.join
MAIN_PATH = '/home/cgarin/PycharmProjects/EDNiX/'
sys.path.insert(1, MAIN_PATH)

import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Exemples.Study_EDNiX.analysis.library import multimodal_profiles as MP
from Exemples.Study_EDNiX.analysis.ednix_bids_tools import PAPER_RC
from Plotting import Plot_interspecies_regional_profiles as PIRP

# ═══════════════════════════════════════════════════════════════════════════════════════════════════ CONFIG
ATLAS_LIB = opj(MAIN_PATH, 'Atlases_library')
ATLAS = 'EDNIxCSCLR'
LEVEL = 3                                   # hierarchy level (1-4, 1-based) of the Figure 4 connectivity fingerprints
                                            # (run_selection ATLAS_LEVEL=3) = volume 2 of the 4D atlases, the FC
                                            # matrices <atlas>_2_run_* of Figure_2_anat (CORR_ATLAS_LEVEL=2, 0-based)
CORTICAL_ONLY = False                       # False: every brain region of LEVEL (cortical + subcortical: striatum, thalamus...)
REGIONS = None                              # None: every region of LEVEL (see CORTICAL_ONLY) present in all modalities and
                                            # species (fingerprint names = label base names); or a list of names
SPECIES = ['Human', 'Macaque']              # first = reference of the statistics of panel B
DESC = ('bratio', 'bratio-mean')            # static human run / dynamic macaque run (window SUV_window)
OUT_DIR = '/scratch2/EDNiX/results/EDNiX_Figures/Figure_XX_multimodal'

SOURCES = {
    'morphometry': {
        'Human': ['/scratch2/EDNiX/Human/BIDS_Merida',
                  '/scratch2/EDNiX/Human/BIDS_Castrillon', '/scratch2/EDNiX/Human/BIDS_Babayan'],
        'Macaque': ['/scratch2/EDNiX/Macaque/BIDS_Tremblay', '/scratch2/EDNiX/Macaque/BIDS_BenHamed',
                    '/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin'],
    },
    'FC': {
        'Human': ['/scratch2/EDNiX/Human/BIDS_Castrillon',
                  '/scratch2/EDNiX/Human/BIDS_Babayan'],
        'Macaque': ['/scratch2/EDNiX/Macaque/BIDS_BenHamed', '/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin',],
    },
    'metabolism': {
        'Human': ['/scratch2/EDNiX/Human/BIDS_Merida'],
        'Macaque': ['/scratch2/EDNiX/Macaque/BIDS_Tremblay'],
    },
}
FIT_KIND = 'correlation'
PANEL_C = 'summary'                         # 'summary': heatmap of the profiles + cross-species conservation and superior-inferior
                                            # position (forest plot); 'radar': per-species radars + consensus
FC_SOURCE = 'sba'                           # connectivity of panel c: 'sba' = retrosplenial seed Q80/Q20 of the runs kept by the
                                            # Figure 4 selection (10 per dataset); 'matrix' = strength from the correlation matrices
SBA_INDEX = 'tails'                         # 'tails': positive / negative tail of the seed map in r units ('pos' / 'neg');
                                            # 'z': Q80 / Q20 standardised across regions (nearly identical profiles)
SBA_WHICH = 'selected'                      # 'selected' (kept runs) | 'all' (every run) for FC_SOURCE = 'sba'
FC_STRENGTH = 'positive'                    # 'positive' | 'absolute' | 'signed' (MP.STRENGTHS); the others go to the sensitivity sheet
MIN_VOXELS = 10                             # PET voxels per hemisphere label, every subject

# panel A
RUN_SURFACES = True                         # False: reuse the PNGs of a previous run
SURFACE_TEMPLATE = 'EDNiX'                  # `reference` of the PET launchers (pet/templates/EDNiX)
SPECIES_FRAGMENT = {'Human': 'atlas/mammals/primates/Human', 'Macaque': 'atlas/mammals/primates/Macaque'}
BRATIO_RANGE = None                         # (low, high) shared scale; None = 1st-99th percentile of both maps

# statistics of panel B
TEST, N_PERM, ALPHA = 'permutation', 10000, 0.05

# colours of the article figures
from Exemples.Study_EDNiX.analysis.paper_style import SPECIES_COLORS  # noqa: E402

MM = 1 / 25.4
FIG_W = 183 * MM                            # Nature double column


# ═════════════════════════════════════════════════════════════════════════════════════════════════ panel A
def surface_pngs():
    """Group mean bratio of each species in the template space, rendered on its surfaces with one shared scale.
    Returns ({species: png}, (low, high))."""
    from Statistics.Group_PET import pet_group, pet_surface
    from Tools import Load_EDNiX_requirement
    means = {}
    for sp in SPECIES:
        bids = MP._as_list(SOURCES['metabolism'][sp])[0]
        atlas_ednix = opj(ATLAS_LIB, SPECIES_FRAGMENT[sp], 'EDNiX')
        if RUN_SURFACES:
            means[sp] = pet_group.run_group_mean(
                bids, SURFACE_TEMPLATE, 'bratio', rel_threshold=0,
                tissue_mask=opj(atlas_ednix, 'volumes', 'masks', sp + '_desc-Whole_Brain_Gray_mask.nii.gz'),
                background=opj(atlas_ednix, 'volumes', sp + '_space-acpc_desc-SS_T1w.nii.gz'))
        else:
            out = opj(bids, 'Results', 'Grp_PET_mean_desc-bratio')
            means[sp] = {'output_dir': out, 'maps': {'mean_all': opj(out, 'group-all_mean.nii.gz')}}
    rng = BRATIO_RANGE
    if rng is None:
        vals = []
        for m in means.values():
            d = np.asarray(nib.load(m['maps']['mean_all']).dataobj, float)
            vals.append(d[np.isfinite(d) & (d != 0)])
        v = np.concatenate(vals)
        rng = (float(np.percentile(v, 1)), float(np.percentile(v, 99)))
    pngs = {}
    for sp, m in means.items():
        out = opj(m['output_dir'], 'surfaces_figureXX')
        if RUN_SURFACES:
            bids = MP._as_list(SOURCES['metabolism'][sp])[0]
            sing_wb = Load_EDNiX_requirement.load_requirement(MAIN_PATH, ATLAS_LIB, bids, 'yes')[4]
            atlas_ednix = opj(ATLAS_LIB, SPECIES_FRAGMENT[sp], 'EDNiX')
            done = pet_surface.render_pet_surfaces(
                m, sp, ATLAS_LIB, SPECIES_FRAGMENT[sp], sing_wb,
                scene_template=opj(atlas_ednix, 'surfaces', 'Native_resol', sp + '_thickness.scene'),
                maps=['mean_all'], projection='ribbon', surface_type='inflated', output_dir=out,
                ranges={'mean_all': rng})
            pngs[sp] = done['mean_all'].replace('.dscalar.nii', '.png') if 'mean_all' in done else None
        else:
            found = sorted(glob.glob(opj(out, '*.png')))
            pngs[sp] = found[0] if found else None
    return pngs, rng


# ═════════════════════════════════════════════════════════════════════════════════════════════════ panel B
def select_regions(tables):
    """REGIONS, or every region of LEVEL found in the tables (cortical only if CORTICAL_ONLY) (Atlas_nomenclature order)."""
    if REGIONS:
        return list(REGIONS)
    names = pd.concat([t['region'] for t in tables.values()]).unique()
    return MP.cortical_regions(names, LEVEL) if CORTICAL_ONLY else MP.brain_regions(names, LEVEL)


def metabolism_stats(met, regions):
    """Subjects x regions bratio (regions present in every subject of both species) and the species comparison of
    Plot_interspecies_regional_profiles."""
    regions, cov = MP.complete_regions(met, regions, SPECIES, 'metabolism')
    if len(regions) < 3:
        raise ValueError(f'metabolism: only {len(regions)} region(s) with a value in every subject; see the coverage '
                         f'table (MIN_VOXELS = {MIN_VOXELS}, LEVEL = {LEVEL})')
    w = MP.subject_matrix(met, regions).reindex(SPECIES, level=0)
    stats = PIRP.compare_species(w, SPECIES[0], TEST, N_PERM)
    return w, stats, cov


# ═════════════════════════════════════════════════════════════════════════════════════════════════ figure
def _axes_mm(fig, x, top, w, h, **kw):
    """Axes at x mm from the left and `top` mm from the top of the figure, w x h mm."""
    W, H = fig.get_size_inches() / MM
    return fig.add_axes([x / W, 1 - (top + h) / H, w / W, h / H], **kw)


def _text_mm(fig, texts, fontsize):
    """Width (mm) of the widest of `texts` at `fontsize` (multi-line texts allowed)."""
    r = fig.canvas.get_renderer()
    w = 0.0
    for t in texts:
        a = fig.text(0, 0, t, fontsize=fontsize)
        w = max(w, a.get_window_extent(r).width / fig.dpi / MM)
        a.remove()
    return w


def _panel_c_summary(fig, W, L, top, h, prof_c, regions_c, summ, si):
    """Left: heatmap (regions x [volume, + coupling, - coupling, FDG, superior-inferior position] per species,
    z of the mean profile); right: forest plot of the Spearman r (95 % CI over regions; filled = Holm P < 0.05):
    cross-species conservation of each modality, and relation of each modality to the superior-inferior position."""
    from matplotlib.colors import TwoSlopeNorm
    species = SPECIES[:2]
    mods = [m for m in MP.SUMMARY_MODALITIES if all(m in prof_c[sp] for sp in species)]
    short = {'morphometry': 'volume', 'pos': '+ coupling', 'neg': '- coupling', 'metabolism': 'FDG'}
    cols, heads = [], []
    for k, sp in enumerate(species):
        for m in mods:
            z = prof_c[sp][m].reindex(columns=regions_c).mean().values
            cols.append(z)
            heads.append((sp, short[m]))
        if si is not None:
            v = si[si.species == sp].set_index('region')['si'].reindex(regions_c).values * 2
            cols.append(v)
            heads.append((sp, 'S-I'))
        if k == 0:
            cols.append(np.full(len(regions_c), np.nan))
            heads.append(('', ''))
    mat = np.vstack(cols).T
    lab_w = 40.0
    ax = _axes_mm(fig, L + lab_w, top + 13, 70.0, h - 13)
    ax.imshow(mat, aspect='auto', cmap='RdBu_r', norm=TwoSlopeNorm(0, -2, 2), interpolation='nearest')
    ax.set_yticks(range(len(regions_c)))
    ax.set_yticklabels([MP.radar_name(r, 60).replace('\n', ' ') for r in regions_c], fontsize=5)
    ax.set_xticks(range(len(heads)))
    ax.set_xticklabels([h_[1] for h_ in heads], rotation=60, ha='left', rotation_mode='anchor', fontsize=5)
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for sp_ in ('top', 'right', 'left', 'bottom'):
        ax.spines[sp_].set_visible(False)
    n1 = len(mods) + (1 if si is not None else 0)
    for k, sp in enumerate(species):
        x0 = k * (n1 + 1)
        fig.text((L + lab_w + 70.0 * (x0 + n1 / 2) / len(heads)) / W,
                 1 - (top - 0.5) / (fig.get_figheight() / MM), sp, ha='center', va='top', fontsize=7,
                 fontweight='bold', color=SPECIES_COLORS.get(sp, 'black'))
    # forest plot
    fx = L + lab_w + 70.0 + 28.0
    axf = _axes_mm(fig, fx, top + 13, W - fx - 4, h - 13)
    rows = summ.reset_index(drop=True)
    y = np.arange(len(rows))[::-1].astype(float)
    # gaps between analyses / species
    gap = np.cumsum([0] + [0.8 if (rows.loc[i, 'analysis'], rows.loc[i, 'species']) != (rows.loc[i - 1, 'analysis'], rows.loc[i - 1, 'species']) else 0.0 for i in range(1, len(rows))])
    y = y - gap
    for yi, r in zip(y, rows.itertuples()):
        col = MP.MODALITY_COLORS[r.modality]
        axf.plot([r.ci_low, r.ci_high], [yi, yi], color=col, lw=1.0)
        axf.plot(r.r, yi, 'o', ms=3.2, color=col, mfc=col if r.p_holm < 0.05 else 'white', mew=0.8)
    axf.axvline(0, color='0.3', lw=0.5)
    axf.set_xlim(-1.05, 1.05)
    axf.set_yticks(y)
    axf.set_yticklabels([f'{short[r.modality]}' for r in rows.itertuples()], fontsize=5)
    axf.set_xlabel('Spearman r [95% CI over regions]\nfilled: Holm P < 0.05', fontsize=5.5)
    axf.tick_params(labelsize=5)
    for sp_ in ('top', 'right'):
        axf.spines[sp_].set_visible(False)
    last = None
    for yi, r in zip(y, rows.itertuples()):
        key = (r.analysis, r.species)
        if key != last:
            axf.text(1.07, yi, ('conservation: ' if r.analysis == 'conservation' else 'vs S-I position: ') + r.species,
                     fontsize=5, va='center', ha='left', clip_on=False, color='0.2')
            last = key
    fig.text((fx - 6) / W, 1 - (top - 1) / (fig.get_figheight() / MM), 'd', fontsize=8, fontweight='bold', va='top')


def build_figure(pngs, rng, prof_b, stats_b, prof_c, regions_c, summ=None, si=None):
    """A: the two species side by side (human left, macaque right), image height set by the aspect of the
    renders so that no space is left empty; B: bars of the bratio drawn from 0; C: the two radars with the
    modality legend between them, their size computed from the width of the region names so that names and
    legend never overlap. Positions in mm (Nature double column, FIG_W wide); the height follows the content."""
    with plt.rc_context(PAPER_RC):
        imgs = {sp: (MP.autocrop(plt.imread(pngs[sp])) if pngs.get(sp) and os.path.isfile(pngs[sp]) else None)
                for sp in SPECIES}
        aspect = max([im.shape[1] / im.shape[0] for im in imgs.values() if im is not None] or [3.0])
        W = FIG_W / MM
        L, R, GAP = 6.0, 1.5, 3.0                                     # margins (mm)
        cw = (W - L - R - GAP * (len(SPECIES) - 1)) / len(SPECIES)    # width of one surface image
        ha = cw / aspect
        fs_lab = 5.0
        names_c = [MP.radar_name(r, 13) for r in regions_c]
        fig = plt.figure(figsize=(FIG_W, 1.0))
        lab_w = _text_mm(fig, names_c, fs_lab) + 2.5                  # region name + tick pad, each side
        pt = 25.4 / 72                                                # mm per point
        lab_h = max(n.count('\n') + 1 for n in names_c) * fs_lab * 1.2 * pt + 2.5   # name above / below
        title_h = 2 * (fs_lab + 2) * 1.25 * pt + 2                   # species name and n of the radar titles
        leg_labels = [MP.MODALITY_LABEL[m].replace(' (', '\n(') for m in MP.MODALITIES]
        leg_w = _text_mm(fig, leg_labels + ['z across regions', '(within subject)'], 6) + 8
        d = (W - L - R) / 2 - 2 * lab_w - 2                            # diameter of the three radars (two columns)
        names_b = [MP.short_name(r, 30) for r in prof_b.columns]
        b_lab = _text_mm(fig, [n.split('\n')[0] for n in names_b], 5.5) * np.sin(np.radians(40)) + 4
        plt.close(fig)

        top_a = 6.0
        top_b = top_a + 5 + ha + 9                                   # species name, image, letter + legend
        hb = 48.0
        top_c = top_b + hb + b_lab + 2
        rad_top = top_c + title_h + lab_h                             # top of the circle
        row_h = d + 2 * lab_h + title_h + 4                               # one stacked radar (species left)
        H = rad_top + d + lab_h + row_h + 2
        if PANEL_C == 'summary' and summ is not None:
            hh = 3.0 * len(regions_c) + 22.0
            H = top_c + hh + 4
        else:
            row_h = d + 2 * lab_h + title_h + 4                               # one stacked radar (species left)
            H = rad_top + d + lab_h + row_h + 2
        fig = plt.figure(figsize=(FIG_W, H * MM))

        # A: human left, macaque right
        axa = []
        for k, sp in enumerate(SPECIES):
            ax = _axes_mm(fig, L + k * (cw + GAP), top_a + 5, cw, ha)
            ax.axis('off')
            if imgs[sp] is not None:
                ax.imshow(imgs[sp], aspect='auto')
            else:
                ax.text(0.5, 0.5, 'surface not rendered', ha='center', va='center', fontsize=6, color='0.4',
                        transform=ax.transAxes)
            n = prof_b.loc[sp].shape[0] if sp in prof_b.index.get_level_values(0) else 0
            ax.text(0.5, 1.02, f'{sp} (n = {n})', transform=ax.transAxes, ha='center', va='bottom', fontsize=8,
                    fontweight='bold', color=SPECIES_COLORS.get(sp, 'black'))
            axa.append(ax)
        # (no caption under A: the measure and the shared colour range are in the legend, Clément 2026-10-27)

        # B: bratio, bars from 0 (the values of the surface maps), dashed line = whole-brain mean
        axb = _axes_mm(fig, L + 6, top_b, W - L - 6 - R, hb)
        vals = prof_b.values[np.isfinite(prof_b.values)]
        PIRP.draw_bars(axb, prof_b, stats_b, SPECIES[0], 'bratio', alpha=ALPHA, colors=SPECIES_COLORS,
                       fontsize=6.5, legend=False, bar_from=0, ylim=(0, float(np.max(vals)) * 1.12),
                       bar_alpha=1.0)                    # full species colour, as in the other figures
        axb.set_ylabel('bratio (ratio to whole-brain mean)', fontsize=6.5)
        axb.set_xticklabels(names_b, rotation=40, ha='right', rotation_mode='anchor', fontsize=5.5)
        axb.legend(fontsize=6.5, frameon=False, ncol=len(SPECIES), loc='lower right', bbox_to_anchor=(1, 1.0))

        # C: summary (heatmap + forest plot)
        if PANEL_C == 'summary' and summ is not None:
            _panel_c_summary(fig, W, L, top_c + 6, hh - 6, prof_c, regions_c, summ, si)
        else:
            # C: left, one radar per species (human top, macaque bottom); right, the consensus across species
            rmax = max(np.nanmax(np.abs(m.mean().values)) for mods in prof_c.values() for m in mods.values())
            rmax = float(np.ceil(rmax * 1.2 * 2) / 2)
            cons = MP.consensus_profiles(prof_c, SPECIES[:2], regions_c)
            rmax_c = float(np.ceil(max(np.nanmax(np.abs(np.r_[v[1], v[2]])) for v in cons.values()) * 1.2 * 2) / 2)
            xl = L + lab_w
            xr = L + (W - L - R) / 2 + lab_w
            mods_present = [m for m in MP.MODALITIES if any(m in prof_c[sp] for sp in SPECIES[:2])]
    
            def _yt(ax, r):
                ax.set_yticks([-r / 2, 0, r / 2])
                ax.set_yticklabels([f'{-r / 2:g}', '0', f'{r / 2:g}'], fontsize=5, color='0.4')
                ax.set_rlabel_position(180 / len(regions_c))
    
            for k, sp in enumerate(SPECIES[:2]):
                top = rad_top + k * row_h
                ax = _axes_mm(fig, xl, top, d, d, projection='polar')
                ns = ' / '.join(str(len(prof_c[sp][m])) for m in mods_present if m in prof_c[sp])
                MP.draw_radar(ax, prof_c[sp], regions_c, f'{sp}\nn = {ns}', SPECIES_COLORS.get(sp, 'black'),
                              fontsize=fs_lab, rlim=(-rmax, rmax), width=13)
                _yt(ax, rmax)
            top = rad_top + (row_h - 0) / 2
            ax = _axes_mm(fig, xr, top, d, d, projection='polar')
            MP.draw_radar_consensus(ax, cons, regions_c, fontsize=fs_lab, rlim=(-rmax_c, rmax_c), width=13)
            _yt(ax, rmax_c)
            ax.set_title('Consensus across species\n(mean; shading, range of the species means)', fontsize=fs_lab + 2,
                         fontweight='bold', pad=4 + 1.25 * fs_lab * 3)
            axl = _axes_mm(fig, xr - lab_w, rad_top + (row_h - 0) / 2 + d + lab_h + 3, d + 2 * lab_w, 18)
            axl.axis('off')
            handles = [plt.Line2D([], [], color=MP.MODALITY_COLORS[m], lw=1.6, ls=MP.MODALITY_LINESTYLE[m])
                       for m in mods_present]
            axl.legend(handles, [MP.MODALITY_LABEL[m] for m in mods_present], loc='center', ncol=2, fontsize=6,
                       frameon=False, handlelength=2.4, title='z across regions (within run / subject)',
                       title_fontsize=6, columnspacing=1.5)
        # bold lower-case letters at 8 pt, as every figure of the paper (Nature style)
        for lab, top in zip('abc', [top_a, top_b - 9, top_c]):
            fig.text(1 / W, 1 - top / H, lab, fontsize=8, fontweight='bold', va='top')
        return fig


def legend_text(prof_b, stats_b, prof_c, regions_c, rng, cmc):
    n = {sp: prof_b.loc[sp].shape[0] for sp in SPECIES}
    nc = {sp: {m: len(v) for m, v in prof_c[sp].items()} for sp in SPECIES}
    sig = stats_b[stats_b.q_fdr < ALPHA]
    txt = (
        f'Figure 6 | Multimodal, multi-species imaging approaches to build an intersectional brain function '
        f'profile. a, Group-mean [18F]FDG brain ratio (bratio: regional activity / mean activity of the brain, '
        f'proportional scaling) projected on the inflated cortical surface of each species (ribbon-constrained '
        f'mapping), same colour scale for both species ({rng[0]:.2f}-{rng[1]:.2f}); '
        + ', '.join(f'{sp} n = {n[sp]}' for sp in SPECIES) + '. '
        f'b, Cross-species comparison of glucose metabolism (bratio) in the cortical and subcortical regions of the EDNiX atlas, '
        f'level {LEVEL}, bilateral (valid-voxel-weighted mean of the left and right labels, native PET space), in '
        f'bratio as in a. Bars, mean across subjects +- SEM; dashed line, whole-brain mean (1); dots, one subject (sessions averaged); stars, '
        f'{TEST} test of {", ".join(SPECIES[1:])} vs {SPECIES[0]}, Benjamini-Hochberg q < {ALPHA} (*), '
        f'< {ALPHA / 5:g} (**), < {ALPHA / 50:g} (***) across {prof_b.shape[1]} regions '
        f'({len(sig)} significant). '
        f'c, Per-region multimodal profile (level {LEVEL}, bilateral, {len(regions_c)} regions) of each species: heatmap of '
        f'the mean profile of each modality (z across regions within subject or run, then averaged; red, high; blue, low), '
        f'in atlas order: volume, % of total brain volume (ventricles and CSF excluded); + coupling and - coupling, '
        + ('positive coupling (regional 80th percentile of the voxel values of the retrosplenial seed correlation map, '
           'negative values set to 0) and negative coupling (depth of the regional 20th percentile when below 0, |r|), in the '
           + ('runs kept by the final round of the run selection of Figure 4 (10 per dataset)' if SBA_WHICH == 'selected' else 'runs')
           if FC_SOURCE == 'sba' and SBA_INDEX == 'tails' else 'connectivity') +
        f'; FDG, bratio; S-I, superior-inferior position of the centroid of the region (-1 inferior, +1 superior, two times '
        f'scaled for display). d, Spearman r with its 95% CI (bootstrap over regions) between the human and macaque mean '
        f'profiles of each modality (conservation) and between each mean profile and the superior-inferior position, for '
        f'each species; filled, Holm P < 0.05 over all the tests (permutation of the regions). These are descriptive '
        f'region-level statistics: regions are few, coarse and spatially dependent, and the modalities come from different '
        f'datasets ('
        + '; '.join(f'{sp}: ' + ', '.join(f'{m} n = {k}' for m, k in nc[sp].items()) for sp in SPECIES)
        + '). Region names: '
        f'those of the connectivity fingerprints (EDNiX atlas labels, level {LEVEL}); in b and c, regions are '
        f'grouped by lobe in the order of the atlas of Figure 1 (radars: counter-clockwise from the top; colour '
        f'strip and rim, lobe colour of Figure 1).'
    )
    if len(cmc):
        txt += (' Cross-modal Spearman r between mean profiles (descriptive): '
                + '; '.join(f'{r.species} {r.modality_a}-{r.modality_b} {r.spearman_r:.2f}' for r in cmc.itertuples())
                + '.')
    return txt


def _fc_variants(tables, regions, fmri):
    """Cross-modal Spearman r of the connectivity variants (sensitivity): matrix strength (positive / absolute /
    signed) and, when the fMRI selection is available, retrosplenial Q80/Q20 of the selected and of all runs."""
    sens = []

    def run(label, fc):
        t2 = {k: v for k, v in tables.items() if k not in ('FC', 'Q80', 'Q20', 'pos', 'neg')}
        t2.update(fc)
        t2 = MP.harmonize_region_names(t2)
        p2, _ = MP.multimodal_profiles(t2, regions, SPECIES)
        sens.append(MP.cross_modal_correlations(p2).assign(connectivity=label))

    for kind in MP.STRENGTHS:
        run(f'matrix {kind} strength', {'FC': MP.fc_table(SOURCES['FC'], ATLAS, LEVEL, FIT_KIND, kind)})
    if fmri is not None:
        try:
            c = fmri() if callable(fmri) else fmri
            variants = {f'SBA {idx}, {which} runs': MP.sba_tables(c, SPECIES, which, idx)
                        for idx in ('tails', 'z') for which in ('selected', 'all')}
        except Exception as e:                      # noqa: BLE001
            if FC_SOURCE == 'sba':
                raise
            print(f'  [sensitivity] SBA variants skipped: {type(e).__name__}: {e}')
            variants = {}
        for label, fc in variants.items():
            run(label, fc)
    return pd.concat(sens, ignore_index=True)


def figure6_stats(fmri=None):
    """Tables of every modality, regions, panel B (bratio per subject and region, species test, region
    coverage) and panel C (z-scored multimodal profiles, cross-modal Spearman r). `fmri` = Inputs().fmri (or its
    result): needed for FC_SOURCE = 'sba' (the runs kept by the Figure 4 selection)."""
    cfg = dict(SOURCES=SOURCES, ATLAS=ATLAS, LEVEL=LEVEL, ATLAS_LIB=ATLAS_LIB, DESC=DESC, FIT_KIND=FIT_KIND,
               FC_STRENGTH=FC_STRENGTH,
               MIN_VOXELS=MIN_VOXELS)
    tables = MP.collect_modalities(cfg)
    if FC_SOURCE == 'sba':
        if fmri is None:
            raise ValueError("FC_SOURCE = 'sba' needs the fMRI selection: figure6_stats(fmri=Inputs().fmri)")
        tables.pop('FC', None)
        tables.update(MP.sba_tables(fmri() if callable(fmri) else fmri, SPECIES, SBA_WHICH, SBA_INDEX))
    tables = MP.harmonize_region_names(tables)
    regions = select_regions(tables)
    prof_b, stats_b, cov_b = metabolism_stats(tables['metabolism'], regions)
    prof_c, regions_c = MP.multimodal_profiles(tables, regions, SPECIES)
    cmc = MP.cross_modal_correlations(prof_c)
    sens = _fc_variants(tables, regions, fmri)
    try:
        si = MP.si_table(SOURCES['morphometry'], ATLAS, LEVEL, ATLAS_LIB)
    except Exception as e:                          # noqa: BLE001
        print(f'  [S-I] superior-inferior position not computed: {type(e).__name__}: {e}')
        si = None
    summ = MP.summary_stats(prof_c, regions_c, SPECIES, si)
    return dict(tables=tables, regions=regions, prof_b=prof_b, stats_b=stats_b, cov_b=cov_b,
                prof_c=prof_c, regions_c=regions_c, cmc=cmc, sens=sens, si=si, summ=summ)


def stats_tables(st):
    """{sheet: DataFrame}: every table behind the panels (former Figure_XX_multimodal_data.xlsx)."""
    out = {f'{k}_subjects': v for k, v in st['tables'].items()}
    out.update(B_bratio=st['prof_b'].reset_index(), B_stats=st['stats_b'], B_region_coverage=st['cov_b'])
    for sp in SPECIES:
        for m, v in st['prof_c'][sp].items():
            out[f'C_{sp}_{m}'[:31]] = v.reset_index()
    out['C_cross_modal_r'] = st['cmc']
    out['C_connectivity_sensitivity'] = st['sens']
    out['C_summary_stats'] = st['summ']
    if st.get('si') is not None:
        out['C_superior_inferior_position'] = st['si']
    return out


def draw_figure6(st, pre):
    """Panel A (group-mean surfaces, rendered here), the figure (<pre>.png/.svg/.pdf) and its legend
    (<pre>_legend.txt, which needs the colour range of A)."""
    os.makedirs(os.path.dirname(pre), exist_ok=True)
    pngs, rng = surface_pngs()
    fig = build_figure(pngs, rng, st['prof_b'], st['stats_b'], st['prof_c'], st['regions_c'], st.get('summ'), st.get('si'))
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(f'{pre}.{ext}', dpi=300 if ext == 'png' else None, bbox_inches='tight')
    plt.close(fig)
    with open(pre + '_legend.txt', 'w') as f:
        f.write(textwrap.fill(legend_text(st['prof_b'], st['stats_b'], st['prof_c'], st['regions_c'], rng,
                                          st['cmc']), 110) + '\n')
    print(f'  [FIG] {pre}.png')
    return pre


def main():
    """Stand-alone run into OUT_DIR (the paper goes through run_paper.py)."""
    os.makedirs(OUT_DIR, exist_ok=True)
    st = figure6_stats()
    pre = opj(OUT_DIR, 'Figure_6_multimodal')
    draw_figure6(st, pre)
    with pd.ExcelWriter(pre + '_data.xlsx') as xw:
        for k, v in stats_tables(st).items():
            v.to_excel(xw, sheet_name=k, index=False)


if __name__ == '__main__':
    main()
