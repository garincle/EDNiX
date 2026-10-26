"""
Figure XX -- Multimodal, multi-species imaging approaches to build an intersectional brain function profile.

 A  surface maps of the group-mean brain ratio of FDG uptake (bratio), human and macaque, same colour scale
    (Statistics/Group_PET: pet_group.run_group_mean + pet_surface.render_pet_surfaces, ribbon projection);
 B  cross-species comparison of glucose metabolism in the cortical regions of LEVEL (bilateral): bars and
    statistics of Plotting/Plot_interspecies_regional_profiles.py (same rendering as F2_bars);
 C  per-region multimodal profile (morphometry, FC, metabolism) on radar axes, one radar per species: each
    modality z-scored across the regions within each subject, mean +- SEM across subjects.

Data, units and level conventions: multimodal_profiles.py. Colours of the species: those of the article figures
(Figure_2_anat.py, Supp_Figure_4-6_QC.py); colours of the modalities: Okabe-Ito.
Outputs in OUT_DIR: Figure_XX_multimodal.png/.svg/.pdf, Figure_XX_multimodal.txt (legend),
Figure_XX_multimodal_data.xlsx (every table behind the panels).
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
import matplotlib.gridspec as gridspec

from Exemples.Study_EDNiX.atlases import Atlas_nomenclature as AN
from Exemples.Study_EDNiX.analysis import multimodal_profiles as MP
from Exemples.Study_EDNiX.analysis.ednix_bids_tools import PAPER_RC
from Plotting import Plot_interspecies_regional_profiles as PIRP

# ═══════════════════════════════════════════════════════════════════════════════════════════════════ CONFIG
ATLAS_LIB = opj(MAIN_PATH, 'Atlases_library')
ATLAS = 'EDNIxCSCLR'
LEVEL = 3                                   # hierarchy level (1-4, 1-based) of the Figure 5 connectivity fingerprints
                                            # (Launcher_Fig5 ATLAS_LEVEL=3) = volume 2 of the 4D atlases, the FC
                                            # matrices <atlas>_2_run_* of Figure_2_anat (CORR_ATLAS_LEVEL=2, 0-based)
REGIONS = None                              # None: every cortical region of LEVEL present in all modalities and
                                            # species (fingerprint names = label base names); or a list of names
SPECIES = ['Human', 'Macaque']              # first = reference of the statistics of panel B
DESC = ('bratio', 'bratio-mean')            # static human run / dynamic macaque run (window SUV_window)
OUT_DIR = '/scratch2/EDNiX/results/EDNiX_Figures/Figure_XX_multimodal'

SOURCES = {
    'morphometry': {
        'Human': ['/scratch2/EDNiX/Human/BIDS_Merida', '/scratch2/EDNiX/Human/BIDS_Park',
                  '/scratch2/EDNiX/Human/BIDS_Castrillon', '/scratch2/EDNiX/Human/BIDS_Babayan'],
        'Macaque': ['/scratch2/EDNiX/Macaque/BIDS_Tremblay', '/scratch2/EDNiX/Macaque/BIDS_BenHamed',
                    '/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin'],
    },
    'FC': {
        'Human': ['/scratch2/EDNiX/Human/BIDS_Park', '/scratch2/EDNiX/Human/BIDS_Castrillon',
                  '/scratch2/EDNiX/Human/BIDS_Babayan'],
        'Macaque': ['/scratch2/EDNiX/Macaque/BIDS_BenHamed', '/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin',
                    '/scratch2/EDNiX/Macaque/BIDS_Tremblay'],
    },
    'metabolism': {
        'Human': ['/scratch2/EDNiX/Human/BIDS_Merida'],
        'Macaque': ['/scratch2/EDNiX/Macaque/BIDS_Tremblay'],
    },
}
FIT_KIND = 'correlation'
MIN_VOXELS = 10                             # PET voxels per hemisphere label, every subject

# panel A
RUN_SURFACES = True                         # False: reuse the PNGs of a previous run
SURFACE_TEMPLATE = 'EDNiX'                  # `reference` of the PET launchers (pet/templates/EDNiX)
SPECIES_FRAGMENT = {'Human': 'atlas/mammals/primates/Human', 'Macaque': 'atlas/mammals/primates/Macaque'}
BRATIO_RANGE = None                         # (low, high) shared scale; None = 1st-99th percentile of both maps

# statistics of panel B
TEST, N_PERM, ALPHA = 'permutation', 10000, 0.05

# colours of the article figures (Figure_2_anat.SPECIES_COLORS)
SPECIES_COLORS = {"Mouse": "#1A5276", "Rat": "#2E86C1", "Bat": "#00A896", "Mouselemur": "#58D68D",
                  "Marmoset": "#F4D03F", "Macaque": "#F5B041", "Human": "#D35400", "Dog": "#8E44AD",
                  "Cat": "#C0392B", "Pig": "#6E2C00"}

MM = 1 / 25.4
FIG_W = 183 * MM                            # Nature double column
FIG_H = 245 * MM                            # full page (Nature maximum 247 mm)


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
                bids, SURFACE_TEMPLATE, 'bratio', rel_threshold=0.8,
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
    """REGIONS, or every cortical region of LEVEL found in the tables (Atlas_nomenclature order)."""
    if REGIONS:
        return list(REGIONS)
    names = pd.concat([t['region'] for t in tables.values()]).unique()
    return MP.cortical_regions(names, LEVEL)


def metabolism_stats(met, regions):
    """Subjects x regions bratio (regions present in every subject of both species) and the species comparison of
    Plot_interspecies_regional_profiles."""
    w = MP.subject_matrix(met, regions)
    regions = [r for r in regions if r in w.columns and w[r].notna().all()]
    w = w[regions].reindex(SPECIES, level=0)
    stats = PIRP.compare_species(w, SPECIES[0], TEST, N_PERM)
    return w, stats


# ═════════════════════════════════════════════════════════════════════════════════════════════════ figure
def build_figure(pngs, rng, prof_b, stats_b, prof_c, regions_c):
    """A: one full-width row per species (surfaces as large as the page allows); B: bars; C: the two radars with
    the modality legend between them."""
    with plt.rc_context(PAPER_RC):
        n_sp = len(SPECIES)
        fig = plt.figure(figsize=(FIG_W, FIG_H))
        # rows: surfaces of each species, bars, room for the bar labels, radars
        gs = gridspec.GridSpec(n_sp + 3, 1, figure=fig, height_ratios=[0.8] * n_sp + [1.0, 0.5, 1.35], hspace=0.12,
                               left=0.07, right=0.99, top=0.98, bottom=0.03)
        # A: species name on the left, surface image over the whole width
        axa = []
        for k, sp in enumerate(SPECIES):
            ax = fig.add_subplot(gs[k, 0])
            n = prof_b.loc[sp].shape[0] if sp in prof_b.index.get_level_values(0) else 0
            MP.draw_surface(ax, pngs.get(sp), '')
            ax.text(-0.01, 0.5, f'{sp}\nn = {n}', transform=ax.transAxes, rotation=90, ha='right', va='center',
                    fontsize=8, fontweight='bold', color=SPECIES_COLORS.get(sp, 'black'))
            axa.append(ax)
        axa[-1].text(0.0, -0.02, f'bratio {rng[0]:.2f} - {rng[1]:.2f}, same scale for both species',
                     transform=axa[-1].transAxes, ha='left', va='top', fontsize=6.5)
        # B
        axb = fig.add_subplot(gs[n_sp, 0])
        PIRP.draw_bars(axb, prof_b, stats_b, SPECIES[0], 'ratio to whole-brain mean', alpha=ALPHA,
                       colors=SPECIES_COLORS, fontsize=6.5, legend=False)
        axb.set_xticklabels([MP.short_name(r, 30) for r in prof_b.columns], rotation=40, ha='right',
                            rotation_mode='anchor', fontsize=5.5)
        axb.legend(fontsize=6.5, frameon=False, ncol=len(SPECIES), loc='lower right', bbox_to_anchor=(1, 1.0))
        # C: radar | legend | radar
        # side columns leave room for the region names, which extend outside the polar axes
        gsc = gridspec.GridSpecFromSubplotSpec(1, 5, subplot_spec=gs[n_sp + 2, 0],
                                               width_ratios=[0.22, 1, 0.95, 1, 0.22], wspace=0.0)
        rmax = max(np.nanmax(np.abs(m.mean().values) + np.nan_to_num(m.sem().values))
                   for mods in prof_c.values() for m in mods.values())
        rmax = float(np.ceil(rmax * 2) / 2)
        axc = []
        for k, sp in enumerate(SPECIES[:2]):
            ax = fig.add_subplot(gsc[0, 1 + 2 * k], projection='polar')
            ns = ' / '.join(str(len(prof_c[sp][m])) for m in MP.MODALITIES if m in prof_c[sp])
            MP.draw_radar(ax, prof_c[sp], regions_c, f'{sp}\nn = {ns}', SPECIES_COLORS.get(sp, 'black'),
                          fontsize=5.5, rlim=(-rmax, rmax))
            ax.set_xticklabels([MP.short_name(r, 15) for r in regions_c], fontsize=5)
            ax.set_yticks([-rmax / 2, 0, rmax / 2])
            ax.set_yticklabels([f'{-rmax / 2:g}', '0', f'{rmax / 2:g}'], fontsize=5, color='0.4')
            ax.set_rlabel_position(90 / len(regions_c))
            axc.append(ax)
        axl = fig.add_subplot(gsc[0, 2])
        axl.axis('off')
        handles = [plt.Line2D([], [], color=MP.MODALITY_COLORS[m], lw=1.8) for m in MP.MODALITIES]
        axl.legend(handles, [MP.MODALITY_LABEL[m].replace(' (', '\n(') for m in MP.MODALITIES], loc='center',
                   fontsize=6, frameon=False, title='z across regions\n(within subject)', title_fontsize=6,
                   handlelength=1.5, labelspacing=1.2)
        for lab, ax in zip('ABC', [axa[0], axb, axc[0]]):
            bb = ax.get_position()
            fig.text(0.02, bb.y1 + 0.012, lab, fontsize=12, fontweight='bold', va='bottom')
        return fig


def legend_text(prof_b, stats_b, prof_c, regions_c, rng, cmc):
    n = {sp: prof_b.loc[sp].shape[0] for sp in SPECIES}
    nc = {sp: {m: len(v) for m, v in prof_c[sp].items()} for sp in SPECIES}
    sig = stats_b[stats_b.q_fdr < ALPHA]
    txt = (
        f'Figure XX | Multimodal, multi-species imaging approaches to build an intersectional brain function '
        f'profile. (A) Group-mean [18F]FDG brain ratio (bratio: regional activity / mean activity of the brain, '
        f'proportional scaling) projected on the inflated cortical surface of each species (ribbon-constrained '
        f'mapping), same colour scale for both species ({rng[0]:.2f}-{rng[1]:.2f}); '
        + ', '.join(f'{sp} n = {n[sp]}' for sp in SPECIES) + '. '
        f'(B) Cross-species comparison of glucose metabolism (bratio) in the cortical regions of the EDNiX atlas, '
        f'level {LEVEL}, bilateral (valid-voxel-weighted mean of the left and right labels, native PET space). Bars, mean '
        f'across subjects +- SEM drawn from the whole-brain mean (1); dots, one subject (sessions averaged); stars, '
        f'{TEST} test of {", ".join(SPECIES[1:])} vs {SPECIES[0]}, Benjamini-Hochberg q < {ALPHA} (*), '
        f'< {ALPHA / 5:g} (**), < {ALPHA / 50:g} (***) across {prof_b.shape[1]} regions '
        f'({len(sig)} significant). '
        f'(C) Per-region multimodal profile of each species on radar axes (level {LEVEL}, bilateral, {len(regions_c)} '
        f'cortical regions): morphometry, region volume in % of total brain volume (ventricles and CSF excluded); '
        f'FC, average connectivity (mean Fisher z of the resting-state correlations of the region with every other '
        f'region); metabolism, bratio. Each modality is z-scored across regions within each subject, so radii '
        f'compare the shape of the profiles, not absolute levels; line, mean; band, SEM across subjects ('
        + '; '.join(f'{sp}: ' + ', '.join(f'{m} n = {k}' for m, k in nc[sp].items()) for sp in SPECIES)
        + '). Modalities come from different datasets of each species (multimodal_profiles.py). Region names: '
        f'those of the connectivity fingerprints (EDNiX atlas labels, level {LEVEL}).'
    )
    if len(cmc):
        txt += (' Cross-modal Spearman r between mean profiles (descriptive): '
                + '; '.join(f'{r.species} {r.modality_a}-{r.modality_b} {r.spearman_r:.2f}' for r in cmc.itertuples())
                + '.')
    return txt


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    cfg = dict(SOURCES=SOURCES, ATLAS=ATLAS, LEVEL=LEVEL, ATLAS_LIB=ATLAS_LIB, DESC=DESC, FIT_KIND=FIT_KIND,
               MIN_VOXELS=MIN_VOXELS)
    tables = MP.collect_modalities(cfg)
    regions = select_regions(tables)
    prof_b, stats_b = metabolism_stats(tables['metabolism'], regions)
    prof_c, regions_c = MP.multimodal_profiles(tables, regions, SPECIES)
    cmc = MP.cross_modal_correlations(prof_c)
    pngs, rng = surface_pngs()
    fig = build_figure(pngs, rng, prof_b, stats_b, prof_c, regions_c)
    pre = opj(OUT_DIR, 'Figure_XX_multimodal')
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(f'{pre}.{ext}', dpi=300 if ext == 'png' else None, bbox_inches='tight')
    plt.close(fig)
    with open(pre + '.txt', 'w') as f:
        f.write(textwrap.fill(legend_text(prof_b, stats_b, prof_c, regions_c, rng, cmc), 110) + '\n')
    with pd.ExcelWriter(pre + '_data.xlsx') as xw:
        for k, v in tables.items():
            v.to_excel(xw, sheet_name=f'{k}_subjects', index=False)
        prof_b.reset_index().to_excel(xw, sheet_name='B_bratio', index=False)
        stats_b.to_excel(xw, sheet_name='B_stats', index=False)
        for sp in SPECIES:
            for m, v in prof_c[sp].items():
                v.reset_index().to_excel(xw, sheet_name=f'C_{sp}_{m}'[:31], index=False)
        cmc.to_excel(xw, sheet_name='C_cross_modal_r', index=False)
    print(f'  [FIG] {pre}.png')


if __name__ == '__main__':
    main()
