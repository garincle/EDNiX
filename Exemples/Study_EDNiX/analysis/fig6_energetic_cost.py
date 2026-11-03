"""
fig6_energetic_cost.py -- energetic cost of signalling across species (Castrillon et al., Sci. Adv. 2023).

    python fig6_energetic_cost.py [--refresh]

  1. Coupling FDG ~ degree of functional connectivity, voxel-wise, in each species (cortex only):
       degree = sum of the positive correlations (p < 0.001) with every other cortical voxel, per run (z-scored),
       averaged over the N_RUNS runs kept by the Figure 4 selection (taken in turn from the datasets of the species),
       compared with the group-mean bratio map of the PET (same template space). Spearman r, P from phase-randomised
       surrogate maps of the degree map (spatial autocorrelation), residual map = "energetic cost".
  2. Human - macaque difference in FDG vs relative expansion of the region (log2 volume ratio of our morphometry).

Output: <PAPER_DIR>/Figure_6/Energetic_cost/ energetic_cost_stats.xlsx / .txt, coupling_<species>.png,
residual maps <species>_energetic_cost.nii.gz (grid of the fMRI template).
"""
import glob
import os
import re
import sys

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Exemples.Study_EDNiX.analysis import paper_data as pdat
from Exemples.Study_EDNiX.analysis import paper_stats as pst
from Exemples.Study_EDNiX.analysis.library import energetic_cost as EC
from Exemples.Study_EDNiX.analysis.library import multimodal_figure as mmf
from Exemples.Study_EDNiX.analysis.library import multimodal_profiles as MP

opj = os.path.join

SPECIES = ('Human', 'Macaque')
N_RUNS = 10                                  # runs per species (the best = kept by the final selection round)
FWHM_MM = {'Human': 6.0, 'Macaque': 3.0}     # 6 mm in Castrillon et al.; macaque scaled by the brain-size ratio (~2)
N_SURR = 1000
P_DEGREE = 0.001
PET_MAP = opj('Results', 'Grp_PET_mean_desc-bratio', 'group-all_mean.nii.gz')
OUT = opj(pdat.PAPER_DIR, 'Figure_6', 'Energetic_cost')


# --------------------------------------------------------------------------------------------- runs
def select_runs(c, species, n=N_RUNS):
    """[(bids_dir, subject, session, run or None)] : the runs kept by the final selection round, taken in turn from
    the datasets of `species` (sorted) until n runs are found."""
    per = []
    for cfg in c['bids_configs']:
        if cfg['species'] != species:
            continue
        keys = sorted(c['kept_runs_final'].get(cfg['label'], []))
        per.append([(cfg['bids_dir'], cfg['label']) + tuple(k) for k in keys])
    out = []
    while len(out) < n and any(per):
        for lst in per:
            if lst and len(out) < n:
                out.append(lst.pop(0))
    return out


def find_bold(bids_dir, sub, ses, run=None):
    base = opj(bids_dir, f'sub-{sub}', f'ses-{ses}')
    if not os.path.isdir(base):
        base = opj(bids_dir, f'sub-{sub}')
    pats = [opj(base, 'func', 'templates', 'EDNiX', 'postprocessed_rs', '*desc-fMRI_residual.nii.gz'),
            opj(base, 'func', '*_space-template_desc-fMRI_residual.nii.gz')]
    files = sorted(f for p in pats for f in glob.glob(p))
    if run is not None:
        sel = [f for f in files if f'run-{run}' in os.path.basename(f)]
        files = sel or files
    return files[0] if files else None


def cortex_mask(bids_dir, sub, ses, label_path, level=mmf.LEVEL):
    """Boolean cortical mask (grid of the fMRI template) and {label id: region} of the level, from the subject's
    EDNIxCSCLR dseg in template space."""
    ses_dir = opj(bids_dir, f'sub-{sub}', f'ses-{ses}')
    if not os.path.isdir(ses_dir):
        ses_dir = opj(bids_dir, f'sub-{sub}')
    dseg = opj(ses_dir, 'func', 'templates', 'EDNiX', 'labels', f'{sub}_seg-{mmf.ATLAS}_dseg.nii.gz')
    if not os.path.exists(dseg):
        return None, None, None
    lab = np.asarray(nib.load(dseg).dataobj)
    lab = np.rint(lab[..., level - 1] if lab.ndim == 4 else lab).astype(np.int64)
    names = {}                                   # ITK-SNAP label description: ID R G B A VIS MSH "NAME"
    with open(label_path, encoding='utf-8-sig') as f:
        for line in f:
            m = re.match(r'^\s*(\d+)\s+\d+\s+\d+\s+\d+\s+[\d.]+\s+\d+\s+\d+\s+"([^"]*)"\s*$', line)
            if m and int(m.group(1)) > 0:
                names[int(m.group(1))] = MP.region_key(m.group(2))
    cort = set(MP.cortical_regions(sorted(set(names.values())), level))
    ids = [i for i, n in names.items() if n in cort]
    return np.isin(lab, ids), lab, names


# ------------------------------------------------------------------------------------------ analysis 1
def coupling_for_species(c, sp):
    runs = select_runs(c, sp)
    if not runs:
        raise RuntimeError(f'{sp}: no selected run')
    print(f'  [{sp}] {len(runs)} runs: ' + ', '.join(f'{r[1]}:{r[2]}/{r[3]}' for r in runs))
    maps, aff, mask, lab, names, used = [], None, None, None, None, []
    for bids_dir, label, sub, ses, *run in runs:
        bold = find_bold(bids_dir, sub, ses, run[0] if run else None)
        if bold is None:
            print(f'    [skip] no residual BOLD for {label} sub-{sub} ses-{ses}')
            continue
        cfg = next(x for x in c['bids_configs'] if x['label'] == label)
        m, lb, nm = cortex_mask(bids_dir, sub, ses, cfg['label_path'])
        if m is None:
            print(f'    [skip] no template-space dseg for {label} sub-{sub} ses-{ses}')
            continue
        dm, aff = EC.degree_map(bold, mask=m, p=P_DEGREE)
        maps.append(dm)
        used.append(dict(species=sp, dataset=label, subject=sub, session=ses, bold=bold))
        mask, lab, names = (m, lb, nm) if mask is None else (mask | m, lab, names)
        print(f'    [degree] {os.path.basename(bold)}')
    if len(maps) < 3:
        raise RuntimeError(f'{sp}: only {len(maps)} usable run(s)')
    deg = EC.group_degree_map(maps)
    mask = mask & np.isfinite(deg)
    pet_bids = mmf.SOURCES['metabolism'][sp]
    pet_bids = pet_bids[0] if isinstance(pet_bids, (list, tuple)) else pet_bids
    pet_img = nib.load(opj(pet_bids, PET_MAP))
    pet = EC.resample_to_grid(np.asarray(pet_img.dataobj, np.float32).squeeze(), pet_img.affine, deg.shape, aff)
    overlap = float(np.isfinite(pet[mask]).mean())
    print(f'  [{sp}] cortical voxels {int(mask.sum())}, with PET {overlap:.0%}')
    if overlap < 0.5:
        print(f'  [{sp}] [WARN] less than half of the cortical voxels have a PET value: check that the PET group map '
              f'and the fMRI template are in the same space')
    fw = FWHM_MM[sp]
    deg_s, pet_s = EC.smooth_masked(np.where(mask, deg, np.nan), fw, aff), EC.smooth_masked(np.where(mask, pet, np.nan), fw, aff)
    stt, res = EC.coupling_analysis(deg_s, pet_s, mask, n_surr=N_SURR)
    stt.update(species=sp, n_runs=len(maps), fwhm_mm=fw, pet_overlap=overlap)
    reg = EC.regional_means(res, lab, names).rename('residual_bratio').rename_axis('region').reset_index()
    reg['species'] = sp
    nib.save(nib.Nifti1Image(res, aff), opj(OUT, f'{sp}_energetic_cost.nii.gz'))
    fig, ax = plt.subplots(figsize=(3.2, 3.0))
    m = mask & np.isfinite(deg_s) & np.isfinite(pet_s)
    idx = np.random.default_rng(0).choice(np.flatnonzero(m.ravel()), min(5000, int(m.sum())), replace=False)
    ax.scatter(deg_s.ravel()[idx], pet_s.ravel()[idx], s=2, alpha=.3, color='0.3', rasterized=True)
    ax.set_xlabel('degree of connectivity (z)')
    ax.set_ylabel('bratio')
    ax.set_title(f"{sp}: r = {stt['r']:+.2f}, P_surr = {stt['p_surrogate']:.3g}", fontsize=8)
    fig.tight_layout()
    fig.savefig(opj(OUT, f'coupling_{sp}.png'), dpi=200)
    plt.close(fig)
    return stt, reg, pd.DataFrame(used)


# ------------------------------------------------------------------------------------------ analysis 2
def expansion(inp):
    st = inp.multimodal()
    t = st['tables']
    tab = EC.expansion_table(t['morphometry'], t['metabolism'], SPECIES)
    cort = MP.cortical_regions(list(tab.index), mmf.LEVEL)
    res = EC.expansion_analysis(tab, cortical=cort)
    return tab.reset_index(), res


def main(refresh=False):
    os.makedirs(OUT, exist_ok=True)
    inp = pst.Inputs(refresh=refresh)
    lines, sheets = [], {}
    tab, res = expansion(inp)
    sheets['expansion_points'], sheets['expansion_stats'] = tab, res
    lines.append('2. species difference in FDG (log2 human/macaque bratio) vs relative expansion (log2 volume ratio)')
    lines += [f'  {r.regions} (n={r.n_regions}): r = {r.r:+.2f} [{r.ci_low:+.2f}, {r.ci_high:+.2f}], P = {r.p:.3g}, '
              f'Holm P = {r.p_holm:.3g}' for r in res.itertuples()]
    c = inp.fmri()
    stats, regs, used = [], [], []
    for sp in SPECIES:
        try:
            s, r, u = coupling_for_species(c, sp)
        except Exception as e:                                   # noqa: BLE001
            import traceback
            traceback.print_exc()
            print(f'  [{sp}] [ERROR] {type(e).__name__}: {e}')
            continue
        stats.append(s), regs.append(r), used.append(u)
    if stats:
        cs = pd.DataFrame(stats)
        sheets['coupling_stats'] = cs
        sheets['residual_by_region'] = pd.concat(regs, ignore_index=True)
        sheets['runs_used'] = pd.concat(used, ignore_index=True)
        lines.insert(0, '1. FDG (bratio) ~ positive degree of connectivity, cortical voxels (Spearman; P from phase-'
                        'randomised surrogates of the degree map)')
        lines[1:1] = [f'  {r.species} ({r.n_runs} runs, {r.n_voxels} voxels, FWHM {r.fwhm_mm:g} mm): r = {r.r:+.2f}, '
                      f'P_surr = {r.p_surrogate:.3g}, slope = {r.slope:+.3f}' for r in cs.itertuples()]
    with pd.ExcelWriter(opj(OUT, 'energetic_cost_stats.xlsx')) as w:
        for k, v in sheets.items():
            v.to_excel(w, sheet_name=k[:31], index=False)
    with open(opj(OUT, 'energetic_cost_stats.txt'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    print(f'  [saved] {OUT}')


if __name__ == '__main__':
    main(refresh='--refresh' in sys.argv)
