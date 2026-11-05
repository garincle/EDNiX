"""
Figure 6 -- Multimodal, multi-species imaging approaches to build an intersectional brain function profile.

 a  top: group-mean brain ratio of FDG uptake (bratio) of each FDG dataset (two human datasets: Merida = CERMEP-IDB-MRXFDG,
    static frame; Castrillon = simultaneous PET/MR, late frames 20-45 min; one macaque dataset: Tremblay), projected on
    the inflated surface of its species, one shared colour scale (Statistics/Group_PET: pet_group.run_group_mean +
    pet_surface.render_pet_surfaces, ribbon projection). Bottom: regional difference with the macaque (mean bratio of the
    human dataset - mean bratio of the macaque, per homologous EDNiX region) painted on the human atlas and surface;
    symmetric colour scale shared by the two human datasets.
 b  regional bratio of the three datasets (bars from 0 = the values of the maps of a), every subject as a dot; stars =
    difference of each human dataset with the macaque (permutation test, BH-FDR across regions within a dataset).
 c  reproducibility of the two human datasets: left, regional group means of one dataset against the other (identity
    line; Pearson r, Lin's concordance, disattenuated r with the split-half reliabilities); right, Pearson r between the
    regional profiles of pairs of subjects, by class (same subject / other session = eyes open vs closed in Castrillon;
    between subjects within a dataset; between the two human datasets; human vs macaque), label-permutation test of
    within - between human datasets (library/pet_datasets.py).
 d  multimodal fingerprint of each species (radar, PANEL_D = 'radar'): volume, positive coupling to the retrosplenial
    seed and FDG, each z-scored across the regions within each subject, mean +- SEM; cross-species Spearman r of each
    modality (conservation) in the legend. PANEL_D = 'summary': heatmap + forest plot of the former panels c-d.
Regions grouped by lobe in the order of the atlas sunburst of Figure 1 (Atlas_nomenclature.atlas_order), radars
counter-clockwise from the top as the sunburst, lobe colours of Figure 1 on the rim and under the bars.

Data, units and level conventions: multimodal_profiles.py, pet_datasets.py. Colours: species colour of the article
(paper_style.py); the datasets of a species differ by marker (paper_style.dataset_marker) and, for the bars, by hatching.
  figure6_stats()        every table and test (paper_stats.figure6)
  draw_figure6()         panel a renders + the figure + its legend (paper_figures.figure6)
Clément 2026-11-05: Castrillon added (a, b), difference with the macaque, reproducibility of the human datasets (c),
multimodal fingerprint (d).
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
from Exemples.Study_EDNiX.analysis.library import pet_datasets as PD
from Exemples.Study_EDNiX.analysis.ednix_bids_tools import PAPER_RC
from Exemples.Study_EDNiX.analysis import paper_style as PS
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
SPECIES = ['Human', 'Macaque']              # order of the radars of d (conservation = SPECIES[0] vs SPECIES[1])
DESC = ('bratio', 'bratio-mean')            # static run / dynamic run (late-frame window SUV_window of the launcher)
OUT_DIR = '/scratch2/EDNiX/results/EDNiX_Figures/Figure_XX_multimodal'

# FDG datasets of panels a-c: (label, species, BIDS). Label = paper_style.DATASETS label (colour + marker).
PET_DATASETS = [
    ('Human_Merida', 'Human', '/scratch2/EDNiX/Human/BIDS_Merida'),
    ('Human_Castrillon', 'Human', '/scratch2/EDNiX/Human/BIDS_Castrillon'),
    ('Macaque_Tremblay', 'Macaque', '/scratch2/EDNiX/Macaque/BIDS_Tremblay'),
]
REFERENCE_DATASET = 'Macaque_Tremblay'      # difference maps / stars of a-b: human dataset - this dataset
HUMAN_DATASETS = ('Human_Merida', 'Human_Castrillon')   # the pair of panel c
SESSION_PREFER = {'Human_Castrillon': 'open'}   # panel c, between-subject classes: ONE session per subject (eyes open,
                                            # as most resting PET); the other session gives the same-subject class

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
    'metabolism': {                         # d: the two human FDG datasets pooled (first = map of fig6_energetic_cost)
        'Human': ['/scratch2/EDNiX/Human/BIDS_Merida', '/scratch2/EDNiX/Human/BIDS_Castrillon'],
        'Macaque': ['/scratch2/EDNiX/Macaque/BIDS_Tremblay'],
    },
}
FIT_KIND = 'correlation'
PANEL_D = 'radar'                           # 'radar': one multimodal fingerprint per species; 'summary': heatmap of the
                                            # profiles + conservation / superior-inferior forest plot (former c-d)
D_MODALITIES = ('morphometry', 'pos', 'metabolism')   # radar lines of d ('neg' ~ mirror of 'pos', r ~ -0.8: left out)
FC_SOURCE = 'sba'                           # connectivity of panel d: 'sba' = retrosplenial seed Q80/Q20 of the runs kept by the
                                            # Figure 4 selection (10 per dataset); 'matrix' = strength from the correlation matrices
SBA_INDEX = 'tails'                         # 'tails': positive / negative tail of the seed map in r units ('pos' / 'neg');
                                            # 'z': Q80 / Q20 standardised across regions (nearly identical profiles)
SBA_WHICH = 'selected'                      # 'selected' (kept runs) | 'all' (every run) for FC_SOURCE = 'sba'
FC_STRENGTH = 'positive'                    # 'positive' | 'absolute' | 'signed' (MP.STRENGTHS); the others go to the sensitivity sheet
MIN_VOXELS = 10                             # PET voxels per hemisphere label, every subject

# panel a
RUN_SURFACES = True                         # False: reuse the PNGs of a previous run
SURFACE_TEMPLATE = 'EDNiX'                  # `reference` of the PET launchers (pet/templates/EDNiX)
SPECIES_FRAGMENT = {'Human': 'atlas/mammals/primates/Human', 'Macaque': 'atlas/mammals/primates/Macaque'}
BRATIO_RANGE = None                         # (low, high) shared scale; None = 1st-99th percentile of all the maps
DIFF_RANGE = None                           # +/- max of the difference maps; None = max |difference| over the regions
SURF_DIR = 'surfaces_figure6'               # sub-folder of each group-mean result
DIFF_DIR = 'Grp_PET_regional_difference_vs_{ref}'   # in <bids>/Results/ of each human dataset

# statistics of panels a-c
TEST, N_PERM, ALPHA = 'permutation', 10000, 0.05

# colours of the article figures
from Exemples.Study_EDNiX.analysis.paper_style import SPECIES_COLORS  # noqa: E402

MM = 1 / 25.4
FIG_W = 183 * MM                            # Nature double column


def ds_name(label):
    """'Human_Castrillon' -> 'Human, Castrillon' (figure text)."""
    sp = next((s for l_, s, _ in PET_DATASETS if l_ == label), None)
    return f'{sp}, {PS.dataset_short(label, sp)}' if sp else label


def ds_species(label):
    return next((s for l_, s, _ in PET_DATASETS if l_ == label), label.split('_')[0])


def ds_style(label):
    """(colour, marker, hatch): species colour, marker of paper_style, hatch for the second dataset of a species."""
    sp = ds_species(label)
    same = [l_ for l_, s, _ in PET_DATASETS if s == sp]
    return (SPECIES_COLORS.get(sp, 'black'), PS.dataset_marker(label, sp),
            '//////' if same.index(label) % 2 else None)


# ═════════════════════════════════════════════════════════════════════════════════════════════════ panel a
def _group_mean(label, sp, bids):
    from Statistics.Group_PET import pet_group
    atlas_ednix = opj(ATLAS_LIB, SPECIES_FRAGMENT[sp], 'EDNiX')
    if RUN_SURFACES:
        return pet_group.run_group_mean(
            bids, SURFACE_TEMPLATE, 'bratio', rel_threshold=0,
            tissue_mask=opj(atlas_ednix, 'volumes', 'masks', sp + '_desc-Whole_Brain_Gray_mask.nii.gz'),
            background=opj(atlas_ednix, 'volumes', sp + '_space-acpc_desc-SS_T1w.nii.gz'))
    found = sorted(glob.glob(opj(bids, 'Results', 'Grp_PET_mean_desc-bratio*', 'group-all_mean.nii.gz')))
    if not found:
        raise FileNotFoundError(f'{label}: no Results/Grp_PET_mean_desc-bratio*/group-all_mean.nii.gz '
                                f'(RUN_SURFACES = True the first time)')
    out = os.path.dirname(found[0])
    return {'output_dir': out, 'maps': {'mean_all': found[0]}}


def _render(result, key, sp, rng, out):
    """PNG of map `key` of a pet_group-like result on the surfaces of sp (or the PNG of a previous run)."""
    if not RUN_SURFACES:
        found = sorted(glob.glob(opj(out, '*.png')))
        return found[0] if found else None
    from Statistics.Group_PET import pet_surface
    from Tools import Load_EDNiX_requirement
    bids = result['bids']
    sing_wb = Load_EDNiX_requirement.load_requirement(MAIN_PATH, ATLAS_LIB, bids, 'yes')[4]
    atlas_ednix = opj(ATLAS_LIB, SPECIES_FRAGMENT[sp], 'EDNiX')
    done = pet_surface.render_pet_surfaces(
        result, sp, ATLAS_LIB, SPECIES_FRAGMENT[sp], sing_wb,
        scene_template=opj(atlas_ednix, 'surfaces', 'Native_resol', sp + '_thickness.scene'),
        maps=[key], projection='ribbon', surface_type='inflated', output_dir=out, ranges={key: rng})
    return done[key].replace('.dscalar.nii', '.png') if key in done else None


def difference_volume(mean_result, diff, sp, out_dir):
    """Regional difference (Series region -> human - macaque bratio) painted into the template atlas of the human
    dataset, on the grid and inside the mask (non-zero voxels) of its group-mean map. Returns the NIfTI path."""
    from Statistics.Group_PET import pet_group
    from Tools.atlas_labels import parse_label_file, get_atlas_label_path, find_species_path
    ref = nib.load(mean_result['maps']['mean_all'])
    m = np.asarray(ref.dataobj, float)
    images = pet_group.collect_images(mean_result['bids'], SURFACE_TEMPLATE, 'bratio')
    atlas_file = pet_group.template_files(images)['atlases'].get(ATLAS)
    if not atlas_file:
        raise FileNotFoundError(f'no {ATLAS} atlas next to the template-space PET images of {mean_result["bids"]}')
    label_df = parse_label_file(get_atlas_label_path(ATLAS_LIB, find_species_path(ATLAS_LIB, sp), ATLAS))
    img = PD.paint_regions(diff, atlas_file, LEVEL - 1, ref, label_df, mask=np.isfinite(m) & (m != 0))
    os.makedirs(out_dir, exist_ok=True)
    path = opj(out_dir, 'regional_difference.nii.gz')
    nib.save(img, path)
    return path


def surface_pngs(diffs=None):
    """Group mean bratio of every dataset of PET_DATASETS in its template space, rendered on the surfaces of its
    species with one shared scale, and (diffs = {human dataset: Series region -> difference with the macaque}) the
    regional difference maps on the human surface with one symmetric scale.
    Returns ({dataset: png}, (low, high), {dataset: png}, (-d, d))."""
    means = {}
    for label, sp, bids in PET_DATASETS:
        means[label] = dict(_group_mean(label, sp, bids), bids=bids)
    rng = BRATIO_RANGE
    if rng is None:
        vals = []
        for m in means.values():
            d = np.asarray(nib.load(m['maps']['mean_all']).dataobj, float)
            vals.append(d[np.isfinite(d) & (d != 0)])
        v = np.concatenate(vals)
        rng = (float(np.percentile(v, 1)), float(np.percentile(v, 99)))
    pngs = {label: _render(means[label], 'mean_all', sp, rng, opj(means[label]['output_dir'], SURF_DIR))
            for label, sp, _ in PET_DATASETS}
    dpngs, drng = {}, None
    if diffs:
        dmax = DIFF_RANGE or float(max(np.nanmax(np.abs(s.values)) for s in diffs.values()))
        drng = (-dmax, dmax)
        for label, sp, bids in PET_DATASETS:
            if label not in diffs:
                continue
            out = opj(bids, 'Results', DIFF_DIR.format(ref=REFERENCE_DATASET))
            if RUN_SURFACES:
                path = difference_volume(means[label], diffs[label], sp, out)
                res = {'output_dir': out, 'bids': bids, 'maps': {'regional_difference': path}}
            else:
                res = {'output_dir': out, 'bids': bids, 'maps': {}}
            dpngs[label] = _render(res, 'regional_difference', sp, drng, opj(out, SURF_DIR))
    return pngs, rng, dpngs, drng


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


# ═════════════════════════════════════════════════════════════════════════════════════ panels b-c (datasets)
def pet_dataset_stats(regions=None):
    """Every table of panels a-c from the ROI tables of PET_DATASETS:
      w        subjects (sessions averaged) x regions, index (dataset, subject)
      stats    each human dataset vs REFERENCE_DATASET, per region (difference, Hedges g, P, BH q)
      between  Castrillon vs Merida, per region (systematic difference between the human datasets)
      agree    agreement of the regional group means of the two human datasets
      pairs / ptest   pairwise subject-profile similarity and the label-permutation test (one session per subject)
      sessions the session used for every subject in the pair classes."""
    df = PD.dataset_table(PET_DATASETS, ATLAS, LEVEL, DESC, MIN_VOXELS, sessions=False)
    dfs = PD.dataset_table(PET_DATASETS, ATLAS, LEVEL, DESC, MIN_VOXELS, sessions=True)
    labels = [l_ for l_, _, _ in PET_DATASETS]
    if regions is None:
        regions = select_regions({'metabolism': df})
    regions, cov = PD.complete_regions(dfs, regions, 'metabolism (datasets)')
    if len(regions) < 3:
        raise ValueError(f'metabolism: only {len(regions)} region(s) with a value in every session of every dataset')
    w = PD.dataset_matrix(df, regions, labels)
    stats = PD.compare_to_reference(w, REFERENCE_DATASET, TEST, N_PERM)
    a, b = HUMAN_DATASETS
    between = PD.compare_to_reference(w.loc[[a, b]], a, TEST, N_PERM)
    agree = PD.agreement(w, a, b)
    ws = PD.session_matrix(dfs, regions, labels)
    one, chosen = PD.one_session(ws, SESSION_PREFER)
    multi = ws[[ds in SESSION_PREFER for ds in ws.index.get_level_values(0)]]
    pairs, ptest = PD.pair_similarity(one, multi, humans=HUMAN_DATASETS, reference=REFERENCE_DATASET, n_perm=N_PERM)
    return dict(w=w, stats=stats, between=between, agree=agree, pairs=pairs, ptest=ptest, sessions=chosen,
                coverage=cov, regions=regions, labels=labels)


def regional_differences(pet):
    """{human dataset: Series region -> mean difference with the reference dataset} (difference maps of a)."""
    s = pet['stats']
    return {d: g.set_index('region')['difference'] for d, g in s.groupby('dataset', sort=False)}


PAIR_ORDER = ['same subject, other session', 'within {a}', 'within {b}', 'between human datasets', 'human vs macaque']


def _draw_bars(ax, w, stats, regions, fontsize=5.5):
    """Grouped bars (mean +- SEM, from 0), one per dataset of PET_DATASETS; dots = subjects (marker of the dataset);
    stars above the bar of each human dataset = its difference with REFERENCE_DATASET (BH q)."""
    labels = [l_ for l_ in (x for x, _, _ in PET_DATASETS) if l_ in w.index.get_level_values(0)]
    nd = len(labels)
    bw = 0.8 / nd
    x = np.arange(len(regions))
    st = stats.set_index(['dataset', 'region'])
    rng = np.random.default_rng(0)
    top = float(np.nanmax(w.values)) * 1.16
    for i, lab in enumerate(labels):
        c, mk, hatch = ds_style(lab)
        d = w.loc[lab, regions]
        xs = x - 0.4 + bw * (i + 0.5)
        ax.bar(xs, d.mean().values, bw * 0.9, yerr=d.sem().values, color='white' if hatch else c, edgecolor=c,
               hatch=hatch, lw=0.6, error_kw=dict(lw=0.6, capsize=1.0, ecolor='0.15'),
               label=f'{ds_name(lab)} (n = {d.shape[0]})')
        for j, reg in enumerate(regions):
            v = d[reg].values
            ax.scatter(xs[j] + rng.uniform(-bw * 0.22, bw * 0.22, len(v)), v, s=3.5, marker=mk, facecolor=c,
                       edgecolors='white', linewidths=0.25, zorder=3)
            if lab != REFERENCE_DATASET and (lab, reg) in st.index:
                s_ = PIRP._stars(st.loc[(lab, reg), 'q_fdr'], ALPHA)
                if s_:
                    ax.text(xs[j], top * 0.985, s_, rotation=90, ha='center', va='top', fontsize=fontsize - 0.5,
                            color=PS.text_color(ds_species(lab)), linespacing=0.6)
    ax.axhline(1.0, color='k', lw=0.5, ls='--')
    ax.set_xticks(x)
    ax.set_xlim(-0.6, len(regions) - 0.4)
    ax.set_ylim(0, top)
    PIRP.draw_lobe_strip(ax, regions)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='y', labelsize=fontsize)


def _draw_agreement(ax, w, agree, regions, fontsize=5.5):
    """Regional group means of the two human datasets (x: first, y: second), SEM bars, lobe colours, identity line."""
    a, b = HUMAN_DATASETS
    ma, mb = w.loc[a, regions].mean(), w.loc[b, regions].mean()
    sa, sb = w.loc[a, regions].sem(), w.loc[b, regions].sem()
    cols = [c if c is not None else (0.5, 0.5, 0.5) for c in PIRP.lobe_colors(regions)]
    lo = float(min(ma.min() - sa.max(), mb.min() - sb.max())) - 0.03
    hi = float(max(ma.max() + sa.max(), mb.max() + sb.max())) + 0.03
    ax.plot([lo, hi], [lo, hi], color='0.5', lw=0.5, ls='--', zorder=0)
    ax.errorbar(ma, mb, xerr=sa, yerr=sb, fmt='none', ecolor='0.6', elinewidth=0.4, capsize=0, zorder=1)
    ax.scatter(ma, mb, s=10, c=cols, edgecolors='0.2', linewidths=0.3, zorder=2)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect('equal')
    ax.set_xlabel(f'{ds_name(a)} bratio', fontsize=fontsize + 0.5)
    ax.set_ylabel(f'{ds_name(b)} bratio', fontsize=fontsize + 0.5)
    ax.tick_params(labelsize=fontsize)
    ax.spines[['top', 'right']].set_visible(False)
    g = agree.iloc[0]
    txt = (f"r = {g.pearson_r:.2f} [{g.pearson_r_ci_low:.2f}, {g.pearson_r_ci_high:.2f}]\n"
           f"CCC = {g.lin_ccc:.2f} [{g.lin_ccc_ci_low:.2f}, {g.lin_ccc_ci_high:.2f}]\n"
           f"r / ceiling = {g.pearson_r_disattenuated:.2f}\n{int(g.n_regions)} regions")
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, ha='left', va='top', fontsize=fontsize - 0.3, linespacing=1.25)


def _draw_pairs(ax, pairs, ptest, fontsize=5.5):
    """Pearson r between subject profiles by pair class (box: median, IQR; dots: pairs)."""
    a, b = HUMAN_DATASETS
    order = [c.format(a=a, b=b) for c in PAIR_ORDER]
    order = [c for c in order if (pairs.pair_class == c).any()]
    names = {f'within {a}': f'within\n{PS.dataset_short(a)}', f'within {b}': f'within\n{PS.dataset_short(b)}',
             'same subject, other session': f'same subject\n{PS.dataset_short(b)}\nopen/closed',
             'between human datasets': f'{PS.dataset_short(a)}\nvs {PS.dataset_short(b)}',
             'human vs macaque': 'human\nvs macaque'}
    colour = {f'within {a}': SPECIES_COLORS['Human'], f'within {b}': SPECIES_COLORS['Human'],
              'same subject, other session': SPECIES_COLORS['Human'], 'between human datasets': SPECIES_COLORS['Human'],
              'human vs macaque': SPECIES_COLORS['Macaque']}
    rng = np.random.default_rng(1)
    data = [pairs.loc[pairs.pair_class == c, 'r'].values for c in order]
    bp = ax.boxplot(data, positions=range(len(order)), widths=0.55, showfliers=False, patch_artist=True,
                    medianprops=dict(color='black', lw=0.8), boxprops=dict(lw=0.5), whiskerprops=dict(lw=0.5),
                    capprops=dict(lw=0.5))
    for patch, c in zip(bp['boxes'], order):
        patch.set_facecolor('white')
        patch.set_edgecolor('0.2')
    for k, (c, v) in enumerate(zip(order, data)):
        n = len(v)
        sub = v if n <= 400 else rng.choice(v, 400, replace=False)          # dots: at most 400 pairs per class
        ax.scatter(k + rng.uniform(-0.2, 0.2, len(sub)), sub, s=1.5, color=colour[c], alpha=0.45, lw=0, zorder=0)
        ax.text(k, 1.005, f'{n}', ha='center', va='bottom', fontsize=fontsize - 1, color='0.35',
                transform=ax.get_xaxis_transform())
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([names.get(c, c) for c in order], fontsize=fontsize - 0.3)
    ax.set_ylabel('profile similarity between\nsubjects (Pearson r)', fontsize=fontsize + 0.5)
    ax.tick_params(axis='y', labelsize=fontsize)
    ax.spines[['top', 'right']].set_visible(False)
    if len(ptest):
        t = ptest.iloc[0]
        ax.text(0.99, 0.03, f'within - between human datasets: {t.statistic:+.3f}, P = {t.p_permutation:.2g}',
                transform=ax.transAxes, ha='right', va='bottom', fontsize=fontsize - 0.5, color='0.2')


def _draw_radars(fig, W, L, R, top, d, lab_w, lab_h, title_h, prof_c, regions_c, summ, fs_lab):
    """d: one multimodal radar per species (D_MODALITIES), side by side, modality legend + conservation r below."""
    mods_present = [m for m in D_MODALITIES if all(m in prof_c[sp] for sp in SPECIES[:2])]
    rmax = max(np.nanmax(np.abs(prof_c[sp][m].reindex(columns=regions_c).mean().values))
               for sp in SPECIES[:2] for m in mods_present)
    rmax = float(np.ceil(rmax * 1.2 * 2) / 2)
    half = (W - L - R) / 2
    for k, sp in enumerate(SPECIES[:2]):
        ax = _axes_mm(fig, L + k * half + (half - d) / 2, top + title_h + lab_h, d, d, projection='polar')
        ns = ' / '.join(str(len(prof_c[sp][m])) for m in mods_present)
        MP.draw_radar(ax, {m: prof_c[sp][m] for m in mods_present}, regions_c, f'{sp}\nn = {ns}',
                      SPECIES_COLORS.get(sp, 'black'), fontsize=fs_lab, rlim=(-rmax, rmax), width=13)
        ax.set_yticks([-rmax / 2, 0, rmax / 2])
        ax.set_yticklabels([f'{-rmax / 2:g}', '0', f'{rmax / 2:g}'], fontsize=fs_lab - 0.5, color='0.4')
        ax.set_rlabel_position(180 / len(regions_c))
    cons = {}
    if summ is not None and len(summ):
        c = summ[summ.analysis == 'conservation'].set_index('modality')
        cons = {m: c.loc[m] for m in mods_present if m in c.index}
    short = {'morphometry': 'volume', 'pos': 'positive coupling', 'neg': 'negative coupling', 'metabolism': 'FDG (bratio)',
             'FC': 'connectivity strength'}
    labels_ = []
    for m in mods_present:
        s_ = short.get(m, m)
        if m in cons:
            r = cons[m]
            s_ += f'   {SPECIES[0]}-{SPECIES[1]} r = {r.r:+.2f} [{r.ci_low:+.2f}, {r.ci_high:+.2f}]' + \
                  (' *' if r.p_holm < 0.05 else '')
        labels_.append(s_)
    axl = _axes_mm(fig, L, top + title_h + 2 * lab_h + d + 1, W - L - R, 14)
    axl.axis('off')
    handles = [plt.Line2D([], [], color=MP.MODALITY_COLORS[m], lw=1.4, ls=MP.MODALITY_LINESTYLE[m])
               for m in mods_present]
    axl.legend(handles, labels_, loc='center', ncol=1, fontsize=5.5, frameon=False, handlelength=2.6, title='z across regions within subject (mean +- SEM); r: cross-species Spearman '
                                        'r of the mean profiles [95% CI], * Holm P < 0.05', title_fontsize=5.5)


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


def build_figure(pngs, rng, pet, prof_c, regions_c, summ=None, si=None, dpngs=None, drng=None):
    """a: surfaces of the datasets (top) and the regional differences with the macaque (bottom, under their human
    dataset); b: bars of the three datasets; c: agreement of the human datasets (scatter) and pairwise profile
    similarity (boxes); d: multimodal radars (or the summary). Positions in mm (Nature double column, FIG_W wide);
    the height follows the content."""
    dpngs = dpngs or {}
    w, regions = pet['w'], pet['regions']
    labels = [l_ for l_, _, _ in PET_DATASETS]
    with plt.rc_context(PAPER_RC):
        imgs = {k: (MP.autocrop(plt.imread(p)) if p and os.path.isfile(p) else None)
                for k, p in list(pngs.items()) + [(('diff', k), v) for k, v in dpngs.items()]}
        aspect = max([im.shape[1] / im.shape[0] for im in imgs.values() if im is not None] or [3.0])
        W = FIG_W / MM
        L, R, GAP = 6.0, 1.5, 3.0
        cw = (W - L - R - GAP * (len(labels) - 1)) / len(labels)
        ha = cw / aspect
        fs_lab = 5.0
        names_c = [MP.radar_name(r, 13) for r in regions_c]
        fig = plt.figure(figsize=(FIG_W, 1.0))
        lab_w = _text_mm(fig, names_c, fs_lab) + 2.5
        pt = 25.4 / 72
        lab_h = max(n.count('\n') + 1 for n in names_c) * fs_lab * 1.2 * pt + 2.5
        title_h = 2 * (fs_lab + 2) * 1.25 * pt + 2
        names_b = [MP.short_name(r, 30) for r in regions]
        b_lab = _text_mm(fig, [n.split('\n')[0] for n in names_b], 5.5) * np.sin(np.radians(40)) + 4
        plt.close(fig)

        top_a = 6.0
        row2 = top_a + 4 + ha + 3
        has_diff = any(imgs.get(('diff', k)) is not None for k in labels) or bool(dpngs)
        top_b = (row2 + 4 + ha + 9) if has_diff else (top_a + 4 + ha + 9)
        hb = 40.0
        top_c = top_b + hb + b_lab + 4
        hc = 44.0
        top_d = top_c + hc + 7
        half = (W - L - R) / 2
        d = float(min(52.0, half - 2 * lab_w - 2))
        if PANEL_D == 'summary' and summ is not None:
            hh = 3.0 * len(regions_c) + 22.0
            H = top_d + hh + 4
        else:
            H = top_d + title_h + 2 * lab_h + d + 17
        fig = plt.figure(figsize=(FIG_W, H * MM))

        # a: one surface per dataset, then the difference maps under their human dataset
        for k, lab in enumerate(labels):
            sp = ds_species(lab)
            ax = _axes_mm(fig, L + k * (cw + GAP), top_a + 4, cw, ha)
            ax.axis('off')
            if imgs.get(lab) is not None:
                ax.imshow(imgs[lab], aspect='auto')
            else:
                ax.text(0.5, 0.5, 'surface not rendered', ha='center', va='center', fontsize=6, color='0.4',
                        transform=ax.transAxes)
            n = w.loc[lab].shape[0] if lab in w.index.get_level_values(0) else 0
            ax.text(0.5, 1.02, f'{ds_name(lab)} (n = {n})', transform=ax.transAxes, ha='center', va='bottom',
                    fontsize=7, fontweight='bold', color=PS.text_color(sp))
            if has_diff and lab != REFERENCE_DATASET:
                ax2 = _axes_mm(fig, L + k * (cw + GAP), row2 + 4, cw, ha)
                ax2.axis('off')
                im = imgs.get(('diff', lab))
                if im is not None:
                    ax2.imshow(im, aspect='auto')
                else:
                    ax2.text(0.5, 0.5, 'difference map not rendered', ha='center', va='center', fontsize=6,
                             color='0.4', transform=ax2.transAxes)
                ax2.text(0.5, 1.02, f'{ds_name(lab)} - {ds_name(REFERENCE_DATASET).split(", ")[0]}',
                         transform=ax2.transAxes, ha='center', va='bottom', fontsize=6.5, color='0.15')
        if has_diff:
            k = labels.index(REFERENCE_DATASET) if REFERENCE_DATASET in labels else len(labels) - 1
            axk = _axes_mm(fig, L + k * (cw + GAP), row2 + 4, cw, ha)
            axk.axis('off')
            axk.text(0.5, 0.5, 'bratio difference, human - macaque,\nper homologous region (EDNiX level '
                     f'{LEVEL})\npainted on the human atlas\nscale +/- {drng[1]:.2f}' if drng else '',
                     ha='center', va='center', fontsize=5.5, color='0.3', transform=axk.transAxes)

        # b: bratio of the three datasets, stars = human dataset vs macaque
        axb = _axes_mm(fig, L + 6, top_b, W - L - 6 - R, hb)
        _draw_bars(axb, w, pet['stats'], regions)
        axb.set_ylabel('bratio (ratio to whole-brain mean)', fontsize=6)
        axb.set_xticklabels(names_b, rotation=40, ha='right', rotation_mode='anchor', fontsize=5.5)
        axb.legend(fontsize=6, frameon=False, ncol=len(labels), loc='lower right', bbox_to_anchor=(1, 1.0))

        # c: reproducibility of the human datasets
        axc1 = _axes_mm(fig, L + 9, top_c, hc - 4, hc - 4)
        _draw_agreement(axc1, w, pet['agree'], regions)
        xc2 = L + 9 + hc + 18
        axc2 = _axes_mm(fig, xc2, top_c, W - xc2 - R - 2, hc - 10)
        _draw_pairs(axc2, pet['pairs'], pet['ptest'])

        # d: multimodal fingerprint
        if PANEL_D == 'summary' and summ is not None:
            _panel_c_summary(fig, W, L, top_d + 6, hh - 6, prof_c, regions_c, summ, si)
        else:
            _draw_radars(fig, W, L, R, top_d, d, lab_w, lab_h, title_h, prof_c, regions_c, summ, fs_lab)
        for lab, top in zip('abcd', [top_a, top_b - 7, top_c - 2, top_d]):
            fig.text(1 / W, 1 - top / H, lab, fontsize=8, fontweight='bold', va='top')
        return fig


def legend_text(pet, prof_c, regions_c, rng, cmc, summ=None, drng=None):
    w, st, g = pet['w'], pet['stats'], pet['agree'].iloc[0]
    labels = [l_ for l_ in pet['labels'] if l_ in w.index.get_level_values(0)]
    n = {l_: w.loc[l_].shape[0] for l_ in labels}
    nc = {sp: {m: len(v) for m, v in prof_c[sp].items()} for sp in SPECIES}
    a, b = HUMAN_DATASETS
    nsig = {d_: int((gg.q_fdr < ALPHA).sum()) for d_, gg in st.groupby('dataset', sort=False)}
    t = pet['ptest'].iloc[0] if len(pet['ptest']) else None
    txt = (
        f'Figure 6 | Multimodal, multi-species imaging approaches to build an intersectional brain function profile. '
        f'a, Top, group-mean [18F]FDG brain ratio (bratio: regional activity / mean activity of the brain, proportional '
        f'scaling) of each dataset, projected on the inflated cortical surface of its species (ribbon-constrained '
        f'mapping), one colour scale ({rng[0]:.2f}-{rng[1]:.2f}); '
        + ', '.join(f'{ds_name(l_)} n = {n[l_]}' for l_ in labels) + '. '
        f'{ds_name(a)}: static 10-min frame; {ds_name(b)}: simultaneous PET/MR, mean of the 20-45 min frames, eyes-open '
        f'and eyes-closed sessions averaged; macaque: anaesthetised. '
        + (f'Bottom, regional difference of each human dataset with the macaque (mean bratio, human - macaque, per '
           f'homologous EDNiX region, level {LEVEL}) painted on the human atlas and surface (scale +/- {drng[1]:.2f}). '
           if drng else '') +
        f'b, Regional bratio (level {LEVEL}, bilateral: valid-voxel-weighted mean of the left and right labels, native '
        f'PET space; {len(pet["regions"])} regions with a value in every session). Bars, mean +- SEM (hatched: '
        f'{ds_name(b)}); dots, subjects; dashed line, whole-brain mean (1). Stars, {TEST} test of each human dataset vs '
        f'the macaque, Benjamini-Hochberg q < {ALPHA} (*), < {ALPHA / 5:g} (**), < {ALPHA / 50:g} (***) across regions '
        f'within a dataset (' + ', '.join(f'{ds_name(d_)} {k} significant' for d_, k in nsig.items()) + '); with '
        f'n = {n.get(REFERENCE_DATASET, 0)} macaques the test is coarse and the difference descriptive. '
        f'c, Reproducibility of the two human datasets. Left, regional group means of {ds_name(a)} vs {ds_name(b)} '
        f'(dot colour, lobe of Figure 1; bars, SEM; dashed, identity): Pearson r = {g.pearson_r:.2f} '
        f'[{g.pearson_r_ci_low:.2f}, {g.pearson_r_ci_high:.2f}], Lin concordance CCC = {g.lin_ccc:.2f} '
        f'[{g.lin_ccc_ci_low:.2f}, {g.lin_ccc_ci_high:.2f}] (95% CI, bootstrap of the subjects); split-half '
        f'reliabilities {g.reliability_a:.2f} and {g.reliability_b:.2f} (Spearman-Brown), disattenuated r = '
        f'{g.pearson_r_disattenuated:.2f}. Right, Pearson r between the regional bratio profiles of pairs of subjects '
        f'(one session per subject: {ds_name(b)} eyes open), by class: the same {ds_name(b)} subject in its two sessions '
        f'(eyes open vs closed), two subjects of the same dataset, one subject of each human dataset, a human and a '
        f'macaque; box, median and interquartile range; number of pairs above. '
        + (f'Mean within-dataset minus between-dataset r = {t.statistic:+.3f}, P = {t.p_permutation:.2g} '
           f'(permutation of the dataset labels of the {int(t.n_subjects)} human subjects, {int(t.n_perm)} permutations). '
           if t is not None else '') +
        f'd, Multimodal fingerprint of each species (level {LEVEL}, {len(regions_c)} regions): volume (% of brain '
        f'volume), positive coupling to the retrosplenial seed (regional 80th percentile of the seed correlation map, '
        + ('runs kept by the final round of the run selection of Figure 4' if SBA_WHICH == 'selected' else 'all runs') +
        f') and FDG (bratio, the two human datasets pooled), each z-scored across regions within subject, mean +- SEM; '
        f'lobe colours on the rim. Cross-species Spearman r of the mean profiles with 95% CI (bootstrap over regions), '
        f'P by permutation of the regions, Holm-corrected; descriptive: regions are few, coarse and spatially dependent, '
        f'and the modalities come from different subjects ('
        + '; '.join(f'{sp}: ' + ', '.join(f'{m} n = {k}' for m, k in nc[sp].items() if m in D_MODALITIES)
                    for sp in SPECIES) + '). '
        f'Regions grouped by lobe in the order of the atlas of Figure 1.'
    )
    if summ is not None and len(summ):
        c = summ[(summ.analysis == 'conservation') & summ.modality.isin(D_MODALITIES)]
        txt += (' Conservation: ' + '; '.join(f'{r.modality} r = {r.r:+.2f} [{r.ci_low:+.2f}, {r.ci_high:+.2f}], '
                                              f'Holm P = {r.p_holm:.2g}' for r in c.itertuples()) + '.')
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
    """Tables of panels a-c (datasets, pet_dataset_stats) and d (z-scored multimodal profiles per species, cross-modal
    and cross-species Spearman r). `fmri` = Inputs().fmri (or its result): needed for FC_SOURCE = 'sba' (the runs kept
    by the Figure 4 selection)."""
    cfg = dict(SOURCES=SOURCES, ATLAS=ATLAS, LEVEL=LEVEL, ATLAS_LIB=ATLAS_LIB, DESC=DESC, FIT_KIND=FIT_KIND,
               FC_STRENGTH=FC_STRENGTH, MIN_VOXELS=MIN_VOXELS)
    tables = MP.collect_modalities(cfg)
    if FC_SOURCE == 'sba':
        if fmri is None:
            raise ValueError("FC_SOURCE = 'sba' needs the fMRI selection: figure6_stats(fmri=Inputs().fmri)")
        tables.pop('FC', None)
        tables.update(MP.sba_tables(fmri() if callable(fmri) else fmri, SPECIES, SBA_WHICH, SBA_INDEX))
    tables = MP.harmonize_region_names(tables)
    regions = select_regions(tables)
    pet = pet_dataset_stats(None if REGIONS is None else list(REGIONS))
    prof_c, regions_c = MP.multimodal_profiles(tables, regions, SPECIES)
    cmc = MP.cross_modal_correlations(prof_c)
    sens = _fc_variants(tables, regions, fmri)
    try:
        si = MP.si_table(SOURCES['morphometry'], ATLAS, LEVEL, ATLAS_LIB)
    except Exception as e:                          # noqa: BLE001
        print(f'  [S-I] superior-inferior position not computed: {type(e).__name__}: {e}')
        si = None
    summ = MP.summary_stats(prof_c, regions_c, SPECIES, si)
    return dict(tables=tables, regions=regions, pet=pet, prof_c=prof_c, regions_c=regions_c, cmc=cmc, sens=sens,
                si=si, summ=summ)


def stats_tables(st):
    """{sheet: DataFrame}: every table behind the panels."""
    out = {f'{k}_subjects': v for k, v in st['tables'].items()}
    pet = st['pet']
    out.update(B_bratio=pet['w'].reset_index(), B_vs_macaque=pet['stats'], B_region_coverage=pet['coverage'],
               C_human_datasets_regions=pet['between'], C_agreement=pet['agree'], C_pair_similarity=pet['pairs'],
               C_pair_test=pet['ptest'], C_sessions_used=pet['sessions'])
    for sp in SPECIES:
        for m, v in st['prof_c'][sp].items():
            out[f'D_{sp}_{m}'[:31]] = v.reset_index()
    out['D_cross_modal_r'] = st['cmc']
    out['D_connectivity_sensitivity'] = st['sens']
    out['D_summary_stats'] = st['summ']
    if st.get('si') is not None:
        out['D_superior_inferior_position'] = st['si']
    return out


def draw_figure6(st, pre):
    """Panel a (group-mean and difference surfaces, rendered here), the figure (<pre>.png/.svg/.pdf) and its legend
    (<pre>_legend.txt, which needs the colour ranges of a)."""
    os.makedirs(os.path.dirname(pre), exist_ok=True)
    pngs, rng, dpngs, drng = surface_pngs(regional_differences(st['pet']))
    fig = build_figure(pngs, rng, st['pet'], st['prof_c'], st['regions_c'], st.get('summ'), st.get('si'),
                       dpngs, drng)
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(f'{pre}.{ext}', dpi=300 if ext == 'png' else None, bbox_inches='tight')
    plt.close(fig)
    with open(pre + '_legend.txt', 'w') as f:
        f.write(textwrap.fill(legend_text(st['pet'], st['prof_c'], st['regions_c'], rng, st['cmc'], st.get('summ'),
                                          drng), 110) + '\n')
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
