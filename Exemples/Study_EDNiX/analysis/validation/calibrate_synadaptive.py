"""
Calibration of the SyNadaptive thresholds on independent labels (BIDS_Park, FreeSurfer recon-all benchmark subjects).

Question
--------
SyNadaptive (modalities/anat/registration_adaptive.py) refines the standard SyNCC (candidate A) with a compliant
masked SyN (candidate B) when A fails intensity-based checks, and keeps B when it raises the k-means tissue Dice.
The tissue Dice is computed from the same intensities that SyN optimises, so a more compliant warp can raise it
without any anatomical gain. This script measures, on labels that the registration never sees, whether
  1. B is better than A at all (paired, per subject),
  2. the intensity gain of B (delta tissue Dice) predicts its anatomical gain (delta label Dice),
  3. the checks on A (tissue Dice, mask Dice) identify the subjects in which B helps,
and which thresholds of the gate maximise the label overlap of the adaptive policy.

Design
------
- Subjects: the BIDS_Park scans with a completed recon-all of the FreeSurfer benchmark (same scan as EDNiX).
- For each subject, A and B are BOTH computed (not gated), exactly as in EDNiX step 6: fixed = EDNiX human template
  (skull-stripped), moving = the EDNiX native skull-stripped T1 (space-acpc), translation first.
- Independent labels: the recon-all aseg.mgz of the subject, brought to the EDNiX native space by a rigid ANTs
  registration of recon-all brainmask.mgz on the EDNiX skull-stripped T1 (same scan, 6 dof; quality logged).
- Evaluation in native space (where EDNiX uses the atlases, step 8): the template seg-4FS (FreeSurfer codes) is
  warped to native space with the inverse transforms of A and of B (nearest neighbour), and compared with the
  recon-all aseg: Dice per structure and hemisphere (cortex, cerebral WM, lateral ventricles, thalamus, caudate,
  putamen, pallidum, hippocampus, amygdala, cerebellar cortex and WM, brainstem), Klein et al. (2009) style.
- Outputs (OUT_DIR): per-subject table (calibration_subjects.csv), per-structure table (calibration_structures.csv),
  statistics and recommended thresholds (calibration_summary.json), figure (calibration.png / .pdf), and a QC
  figure per subject (QC/<subject>_QC.png) written as soon as that subject is finished. Registrations are cached: re-running only redoes the analysis.

n = 10 (Park, ages 87-97): the analysis is descriptive (paired Wilcoxon, Spearman, bootstrap CI); the thresholds it
recommends must be read with that sample size in mind.

Run on the machine that has the data:
    python Exemples/Study_EDNiX/analysis/validation/calibrate_synadaptive.py            # all subjects
    python .../calibrate_synadaptive.py --analyse-only                                     # tables + figure only
"""
import os
import sys
import glob
import json
import time
import argparse
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

MAIN_PATH = '/home/cgarin/PycharmProjects/EDNiX/'
sys.path.insert(0, MAIN_PATH)

opj = os.path.join
opi = os.path.isfile

# ----------------------------------------------------------------------------------------------------------- settings
BIDS_DIR = '/scratch2/EDNiX/Human/BIDS_Park'
FS_SUBJECTS_DIR = '/scratch2/EDNiX/FreeSurfer_benchmark/BIDS_Park'      # recon-all: sub-<S>_ses-<T>
LIBRARY = opj(os.path.dirname(MAIN_PATH), 'Atlases_library')            # same rule as the pipeline
SPECIES, REFERENCE = 'Human', 'EDNiX'
OUT_DIR = '/scratch2/EDNiX/SyNadaptive_calibration/BIDS_Park'
N_JOBS = 5               # subjects in parallel
THREADS_PER_JOB = 8      # ITK threads per subject
FS_ALIGN_MIN_CORR = 0.80  # below: the recon-all -> EDNiX alignment is flagged and the subject excluded

# FreeSurfer aseg codes (left, right); the template seg-4FS uses the same codes.
STRUCTURES = {'cortex': (3, 42), 'cerebral_WM': (2, 41), 'lateral_ventricle': (4, 43), 'thalamus': (10, 49),
              'caudate': (11, 50), 'putamen': (12, 51), 'pallidum': (13, 52), 'hippocampus': (17, 53),
              'amygdala': (18, 54), 'cerebellum_cortex': (8, 47), 'cerebellum_WM': (7, 46), 'brainstem': (16,)}


# ------------------------------------------------------------------------------------------------------------ inputs
def template_files():
    """EDNiX human template (skull-stripped T1w) and its seg-4FS, found as atlases/templatefeat finds them."""
    hits = sorted(glob.glob(opj(LIBRARY, 'atlas', '**', SPECIES, REFERENCE, 'volumes',
                                f'{SPECIES}_space-acpc_desc-SS_T1w.nii.gz'), recursive=True),
                  key=lambda p: (p.count(os.sep), p))
    hits = [h for h in hits if f'{os.sep}freesurfer{os.sep}' not in h]
    if not hits:
        sys.exit(f'No {SPECIES}/{REFERENCE} template under {LIBRARY}/atlas')
    vol = os.path.dirname(hits[0])
    seg = opj(vol, 'labels', f'{SPECIES}_seg-4FS_dseg.nii.gz')
    if not opi(seg):
        sys.exit(f'Template segmentation missing: {seg}')
    return hits[0], seg


def subjects():
    """[(fsid, sub, ses, ednix_T1_SS, fs_dir)] for every recon-all subject with an EDNiX skull-stripped T1."""
    out = []
    for fs_dir in sorted(glob.glob(opj(FS_SUBJECTS_DIR, 'sub-*_ses-*'))):
        fsid = os.path.basename(fs_dir)
        if not opi(opj(fs_dir, 'scripts', 'recon-all.done')) or opi(opj(fs_dir, 'scripts', 'recon-all.error')):
            continue
        sub, ses = fsid[4:].split('_ses-')
        t1 = sorted(glob.glob(opj(BIDS_DIR, f'sub-{sub}', f'ses-{ses}', 'anat', 'native', 'volumes',
                                  '*_space-acpc_desc-SS_*T1w.nii.gz')))
        t1 = [t for t in t1 if 'desc-SS-step' not in t]
        if t1:
            out.append((fsid, sub, ses, t1[0], fs_dir))
        else:
            print(f'WARNING: {fsid}: no EDNiX space-acpc_desc-SS T1w, skipped')
    return out


# ------------------------------------------------------------------------------------------------------------ labels
def dice(a, b):
    s = a.sum() + b.sum()
    return float(2 * np.logical_and(a, b).sum() / s) if s else np.nan


def label_dice(ref, test):
    """Dice per structure (mean of the two hemispheres) between two label arrays with FreeSurfer codes."""
    out = {}
    for name, codes in STRUCTURES.items():
        vals = [dice(ref == c, test == c) for c in codes if (ref == c).any() and (test == c).any()]
        out[name] = float(np.mean(vals)) if vals else np.nan
    return out


def mgz_to_ants(path, tmp_nii):
    import nibabel as nib
    import ants
    img = nib.load(path)
    nib.save(nib.Nifti1Image(np.asanyarray(img.dataobj), img.affine), tmp_nii)
    return ants.image_read(tmp_nii)


# ---------------------------------------------------------------------------------------------------------- per subject
def run_subject(args):
    fsid, sub, ses, t1_path, fs_dir, tpl_path, seg_path = args
    for k in ('ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[k] = str(THREADS_PER_JOB)
    import ants
    from modalities.anat import registration_adaptive as RA

    d = opj(OUT_DIR, fsid)
    os.makedirs(d, exist_ok=True)
    res_json = opj(d, 'result.json')
    if opi(res_json):
        with open(res_json) as f:
            return json.load(f)
    row = dict(fsid=fsid, subject=sub, session=ses, ednix_T1=t1_path)

    moving = ants.image_read(t1_path)               # subject (EDNiX native, skull-stripped)
    fixed = ants.image_read(tpl_path)               # template
    tpl_seg = ants.image_read(seg_path)

    # recon-all aseg -> EDNiX native space (same scan: rigid)
    fs_brain = mgz_to_ants(opj(fs_dir, 'mri', 'brainmask.mgz'), opj(d, 'fs_brainmask.nii.gz'))
    fs_aseg = mgz_to_ants(opj(fs_dir, 'mri', 'aseg.mgz'), opj(d, 'fs_aseg.nii.gz'))
    rig = ants.registration(fixed=moving, moving=fs_brain, type_of_transform='Rigid', aff_metric='mattes',
                            outprefix=opj(d, 'fs_to_ednix_'))
    fs_brain_n = ants.apply_transforms(moving, fs_brain, rig['fwdtransforms'], interpolator='linear')
    aseg_n = ants.apply_transforms(moving, fs_aseg, rig['fwdtransforms'], interpolator='nearestNeighbor')
    ants.image_write(aseg_n, opj(d, 'fs_aseg_in_ednix_native.nii.gz'))
    m = (moving.numpy() > 0) & (fs_brain_n.numpy() > 0)
    row['fs_align_corr'] = float(np.corrcoef(moving.numpy()[m], fs_brain_n.numpy()[m])[0, 1]) if m.sum() > 100 else np.nan

    # A: SyNCC from a translation, exactly as norm2template.norm (step 6, no masks)
    t0 = time.time()
    tr = ants.registration(fixed=fixed, moving=moving, type_of_transform='Translation', aff_metric='mattes',
                           outprefix=opj(d, 'shift_'))
    preA = opj(d, 'A_')
    if not opi(preA + '1Warp.nii.gz'):
        ants.registration(fixed=fixed, moving=moving, type_of_transform='SyNCC',
                          initial_transform=tr['fwdtransforms'], outprefix=preA)
    resA = RA._result(preA)
    row['A_minutes'] = (time.time() - t0) / 60

    # B: the SyNadaptive refinement, from A's affine
    t0 = time.time()
    preB = opj(d, 'B_')
    if not opi(preB + '1Warp.nii.gz'):
        RA.refine(fixed, moving, resA['fwdtransforms'][1], preB)
    resB = RA._result(preB)
    row['B_minutes'] = (time.time() - t0) / 60

    ref = aseg_n.numpy().astype(int)
    for tag, res in (('A', resA), ('B', resB)):
        # label-free checks used by SyNadaptive (template space)
        for k, v in RA.evaluate(fixed, moving, res['fwdtransforms'], None, None).items():
            row[f'{tag}_{k}'] = v
        # independent check: template labels -> native, against recon-all aseg
        seg_n = ants.apply_transforms(moving, tpl_seg, res['invtransforms'], whichtoinvert=[True, False],
                                      interpolator='nearestNeighbor')
        ants.image_write(seg_n, opj(d, f'{tag}_template_seg_in_native.nii.gz'))
        ld = label_dice(ref, seg_n.numpy().astype(int))
        for k, v in ld.items():
            row[f'{tag}_label_{k}'] = v
        row[f'{tag}_label_mean'] = float(np.nanmean(list(ld.values())))
    with open(res_json, 'w') as f:
        json.dump(row, f, indent=1)
    return row



# --------------------------------------------------------------------------------------- per-subject QC figure
import nibabel as nib  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

# FreeSurfer codes merged into a few classes for display
CLASSES = {'cortex': ((3, 42), '#D62728'), 'WM': ((2, 41), '#1F77B4'), 'ventricles': ((4, 43, 5, 44), '#FFD700'),
           'subcortical': ((10, 49, 11, 50, 12, 51, 13, 52), '#2CA02C'), 'hippo/amyg': ((17, 53, 18, 54), '#9467BD'),
           'cerebellum': ((7, 46, 8, 47), '#FF7F0E')}


def load(p):  # noqa: E302
    return np.asanyarray(nib.load(p).dataobj)


def classes(seg):
    out = np.zeros(seg.shape, np.int16)
    for i, (codes, _) in enumerate(CLASSES.values(), 1):
        out[np.isin(seg, codes)] = i
    return out


def slices(vol, c):
    return [np.rot90(vol[c[0], :, :]), np.rot90(vol[:, c[1], :]), np.rot90(vol[:, :, c[2]])]


def subject_qc(d):
    sub = os.path.basename(d)
    with open(opj(d, 'result.json')) as f:
        r = json.load(f)
    t1 = load(r['ednix_T1']).astype(float)
    ref = classes(load(opj(d, 'fs_aseg_in_ednix_native.nii.gz')))
    segs = {t: classes(load(opj(d, f'{t}_template_seg_in_native.nii.gz'))) for t in ('A', 'B')}
    nz = np.argwhere(t1 > 0)
    c = ((nz.min(0) + nz.max(0)) // 2) if len(nz) else np.array(t1.shape) // 2
    vmax = np.percentile(t1[t1 > 0], 99) if (t1 > 0).any() else 1
    colors = [col for _, col in CLASSES.values()]
    fig, ax = plt.subplots(3, 3, figsize=(12, 12), facecolor='k')
    rows = [('recon-all aseg', ref, None), ('A (SyNCC)', segs['A'], ref), ('B (refine)', segs['B'], ref)]
    for i, (name, seg, filled) in enumerate(rows):
        for j, (img, s) in enumerate(zip(slices(t1, c), slices(seg, c))):
            a = ax[i, j]
            a.imshow(img, cmap='gray', vmin=0, vmax=vmax)
            if filled is not None:
                f = slices(filled, c)[j].astype(float)
                f[f == 0] = np.nan
                a.imshow(f, cmap=matplotlib.colors.ListedColormap(colors), vmin=0.5, vmax=len(colors) + 0.5,
                         alpha=0.25, interpolation='nearest')
            for k, col in enumerate(colors, 1):
                if (s == k).any():
                    a.contour(s == k, levels=[0.5], colors=[col], linewidths=0.6)
            a.set_axis_off()
        tag = name.split()[0]
        txt = (f"fs_align r={r.get('fs_align_corr', np.nan):.3f}" if i == 0 else
               f"label Dice {r.get(f'{tag}_label_mean', np.nan):.3f} | tissue {r.get(f'{tag}_tissue_dice', np.nan):.3f}"
               f" | mask {r.get(f'{tag}_mask_dice', np.nan):.3f} | folding {r.get(f'{tag}_folding', np.nan):.1e}")
        ax[i, 0].set_title(f'{name}: {txt}', color='w', loc='left', fontsize=9)
    handles = [matplotlib.patches.Patch(color=col, label=n) for n, (_, col) in CLASSES.items()]
    fig.legend(handles=handles, loc='lower center', ncol=6, fontsize=8, facecolor='k', labelcolor='w')
    fig.suptitle(f'{sub}: contours = labels, fill = recon-all aseg', color='w')
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    os.makedirs(opj(OUT_DIR, 'QC'), exist_ok=True)
    out = opj(OUT_DIR, 'QC', f'{sub}_QC.png')
    fig.savefig(out, dpi=110, facecolor='k')
    plt.close(fig)
    return r, out



# ------------------------------------------------------------------------------------------------------------ analysis
def boot_ci(x, n=10000, seed=0):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return [np.nan, np.nan]
    rng = np.random.default_rng(seed)
    means = rng.choice(x, (n, len(x))).mean(1)
    return [float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))]


def policy_label(df, t_tissue, t_mask, min_gain, max_folding):
    """Label Dice of SyNadaptive for given thresholds (B used only if A fails and B passes the gain test)."""
    fail = (df['A_tissue_dice'] < t_tissue) | (df['A_mask_dice'] < t_mask) | (df['A_folding'] > max_folding)
    take_b = fail & ((df['B_tissue_dice'] - df['A_tissue_dice']) > min_gain) & (df['B_folding'] <= max_folding)
    return np.where(take_b, df['B_label_mean'], df['A_label_mean']), take_b


def analyse(rows):
    import pandas as pd
    from scipy import stats
    df = pd.DataFrame(rows).sort_values('fsid')
    df.to_csv(opj(OUT_DIR, 'calibration_subjects.csv'), index=False)
    bad = df['fs_align_corr'] < FS_ALIGN_MIN_CORR
    if bad.any():
        print(f'WARNING: excluded (recon-all -> EDNiX alignment r < {FS_ALIGN_MIN_CORR}): '
              + ', '.join(df.loc[bad, 'fsid']))
    df = df[~bad].reset_index(drop=True)
    n = len(df)
    df['d_label'] = df['B_label_mean'] - df['A_label_mean']
    df['d_tissue'] = df['B_tissue_dice'] - df['A_tissue_dice']

    # per structure
    srows = []
    for s in STRUCTURES:
        a, b = df[f'A_label_{s}'], df[f'B_label_{s}']
        ok = a.notna() & b.notna()
        p = stats.wilcoxon(b[ok], a[ok]).pvalue if ok.sum() >= 5 and (b[ok] != a[ok]).any() else np.nan
        srows.append(dict(structure=s, n=int(ok.sum()), A_mean=a[ok].mean(), B_mean=b[ok].mean(),
                          delta_mean=(b - a)[ok].mean(), delta_CI95=boot_ci((b - a)[ok]),
                          n_B_better=int(((b - a)[ok] > 0).sum()), wilcoxon_p=p))
    sdf = pd.DataFrame(srows)
    sdf.to_csv(opj(OUT_DIR, 'calibration_structures.csv'), index=False)

    def spear(x, y):
        ok = df[x].notna() & df[y].notna()
        if ok.sum() < 4:
            return dict(rho=np.nan, p=np.nan, n=int(ok.sum()))
        r = stats.spearmanr(df.loc[ok, x], df.loc[ok, y])
        return dict(rho=float(r.correlation), p=float(r.pvalue), n=int(ok.sum()))

    s = dict(n_subjects=n,
             A_label_mean=float(df['A_label_mean'].mean()), B_label_mean=float(df['B_label_mean'].mean()),
             delta_label_mean=float(df['d_label'].mean()), delta_label_CI95=boot_ci(df['d_label']),
             n_B_better=int((df['d_label'] > 0).sum()),
             wilcoxon_B_vs_A_p=float(stats.wilcoxon(df['B_label_mean'], df['A_label_mean']).pvalue)
             if n >= 5 and (df['d_label'] != 0).any() else np.nan,
             delta_tissue_mean=float(df['d_tissue'].mean()),
             B_folding_max=float(df['B_folding'].max()), A_folding_max=float(df['A_folding'].max()),
             minutes_A_median=float(df['A_minutes'].median()), minutes_B_median=float(df['B_minutes'].median()),
             spearman_dtissue_vs_dlabel=spear('d_tissue', 'd_label'),
             spearman_Atissue_vs_dlabel=spear('A_tissue_dice', 'd_label'),
             spearman_Amask_vs_dlabel=spear('A_mask_dice', 'd_label'))

    # threshold sweep of the gate (tissue threshold x min_gain; mask threshold and folding kept at the defaults)
    from modalities.anat import registration_adaptive as RA
    th0 = RA.THRESHOLDS
    grid_t = np.unique(np.r_[0.0, np.round(df['A_tissue_dice'].values + 1e-4, 4), 1.01])
    sweep = []
    for t in grid_t:
        for g in (0.0, 0.005, 0.01, 0.02):
            lab, take = policy_label(df, t, th0['mask_dice'], g, th0['max_folding'])
            sweep.append(dict(tissue_threshold=float(t), min_gain=g, n_B=int(take.sum()),
                              label_mean=float(np.mean(lab))))
    sw = pd.DataFrame(sweep)
    sw.to_csv(opj(OUT_DIR, 'calibration_threshold_sweep.csv'), index=False)
    best = sw.sort_values(['label_mean', 'n_B'], ascending=[False, True]).iloc[0].to_dict()
    cur, _ = policy_label(df, th0['tissue_dice'], th0['mask_dice'], th0['min_gain'], th0['max_folding'])
    s['policy'] = dict(always_A=s['A_label_mean'], always_B=s['B_label_mean'],
                       current_thresholds=float(np.mean(cur)), current_thresholds_values=th0,
                       best_in_sweep=best, oracle=float(np.maximum(df['A_label_mean'], df['B_label_mean']).mean()))
    s['recommendation'] = recommend(s)
    with open(opj(OUT_DIR, 'calibration_summary.json'), 'w') as f:
        json.dump(s, f, indent=1, default=float)
    figure(df, sdf, sw, s)
    print(json.dumps({k: s[k] for k in ('n_subjects', 'A_label_mean', 'B_label_mean', 'delta_label_mean',
                                        'delta_label_CI95', 'n_B_better', 'wilcoxon_B_vs_A_p',
                                        'spearman_dtissue_vs_dlabel', 'recommendation')}, indent=1, default=float))
    return s


def recommend(s):
    """Plain-language reading of the results (thresholds are only suggested when the data support a gate)."""
    ci = s['delta_label_CI95']
    r = s['spearman_dtissue_vs_dlabel']
    if np.isfinite(ci[1]) and ci[1] < 0:
        return 'B lowers the independent label overlap: keep SyNCC (do not use SyNadaptive).'
    if np.isfinite(ci[0]) and ci[0] > 0 and s['B_folding_max'] <= 1e-4:
        return ('B raises the label overlap in the group without folding: consider B for every subject '
                '(gate not needed); check per structure.')
    if not (np.isfinite(r['rho']) and r['rho'] > 0.5 and r['p'] < 0.05):
        return ('No group gain and the tissue-Dice gain does not predict the label gain: the intensity gate is not '
                'validated. Keep SyNCC as default; SyNadaptive should not be claimed as an improvement.')
    b = s['policy']['best_in_sweep']
    return (f"The tissue-Dice gain predicts the label gain (rho={r['rho']:.2f}). Suggested gate: tissue_dice < "
            f"{b['tissue_threshold']:.3f}, min_gain = {b['min_gain']} (B used in {b['n_B']} subject(s)); "
            f"confirm on a second dataset before using it in the paper.")


def figure(df, sdf, sw, s):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 4, figsize=(17, 4.2))
    # a: per-structure paired Dice
    x = np.arange(len(sdf))
    for i, st in enumerate(sdf['structure']):
        a, b = df[f'A_label_{st}'], df[f'B_label_{st}']
        for ya, yb in zip(a, b):
            ax[0].plot([i - 0.15, i + 0.15], [ya, yb], color='0.75', lw=0.6)
        ax[0].scatter(np.full(len(a), i - 0.15), a, s=10, color='#1F77B4', zorder=3)
        ax[0].scatter(np.full(len(b), i + 0.15), b, s=10, color='#D62728', zorder=3)
    ax[0].set_xticks(x, sdf['structure'], rotation=60, ha='right', fontsize=7)
    ax[0].set_ylabel('Dice vs recon-all aseg')
    ax[0].set_title('a  A (blue) vs B (red), per structure', loc='left', fontsize=9)
    # b: intensity gain vs label gain
    ax[1].axhline(0, color='0.6', lw=0.6)
    ax[1].axvline(0, color='0.6', lw=0.6)
    ax[1].scatter(df['d_tissue'], df['d_label'], color='k', s=18)
    r = s['spearman_dtissue_vs_dlabel']
    ax[1].set_xlabel(r'$\Delta$ tissue Dice (B$-$A, intensity)')
    ax[1].set_ylabel(r'$\Delta$ label Dice (B$-$A, independent)')
    ax[1].set_title(f"b  rho = {r['rho']:.2f}, p = {r['p']:.2g}", loc='left', fontsize=9)
    # c: A check vs label gain
    ax[2].axhline(0, color='0.6', lw=0.6)
    ax[2].scatter(df['A_tissue_dice'], df['d_label'], color='k', s=18)
    r = s['spearman_Atissue_vs_dlabel']
    ax[2].set_xlabel('tissue Dice of A (gate variable)')
    ax[2].set_ylabel(r'$\Delta$ label Dice (B$-$A)')
    ax[2].set_title(f"c  rho = {r['rho']:.2f}, p = {r['p']:.2g}", loc='left', fontsize=9)
    # d: policy sweep
    for g, sub in sw.groupby('min_gain'):
        ax[3].plot(sub['tissue_threshold'], sub['label_mean'], marker='.', label=f'min_gain {g}')
    ax[3].axhline(s['A_label_mean'], color='#1F77B4', ls='--', lw=0.8, label='always A (SyNCC)')
    ax[3].axhline(s['B_label_mean'], color='#D62728', ls='--', lw=0.8, label='always B')
    ax[3].set_xlabel('tissue Dice threshold of the gate')
    ax[3].set_ylabel('mean label Dice of the policy')
    ax[3].legend(fontsize=6)
    ax[3].set_title('d  SyNadaptive policy', loc='left', fontsize=9)
    fig.suptitle(f"SyNadaptive calibration, BIDS_Park (n = {s['n_subjects']})", fontsize=10)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(opj(OUT_DIR, f'calibration.{ext}'), dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------------------- main
def main(analyse_only=False):
    os.makedirs(OUT_DIR, exist_ok=True)
    tpl, seg = template_files()
    subs = subjects()
    print(f'template: {tpl}\n{len(subs)} subject(s) with recon-all and EDNiX outputs', flush=True)
    rows = []
    if analyse_only:
        for f in sorted(glob.glob(opj(OUT_DIR, '*', 'result.json'))):
            with open(f) as fh:
                rows.append(json.load(fh))
            subject_qc(os.path.dirname(f))
    else:
        with ProcessPoolExecutor(max_workers=N_JOBS) as ex:
            futs = {ex.submit(run_subject, (*sb, tpl, seg)): sb[0] for sb in subs}
            for f in as_completed(futs):
                try:
                    r = f.result()
                    rows.append(r)
                    print(f"{futs[f]}: done | label Dice A {r['A_label_mean']:.3f} B {r['B_label_mean']:.3f} | "
                          f"tissue A {r['A_tissue_dice']:.3f} B {r['B_tissue_dice']:.3f} | "
                          f"fs_align r {r['fs_align_corr']:.3f}", flush=True)
                    try:
                        print(f'   QC: {subject_qc(opj(OUT_DIR, futs[f]))[1]}', flush=True)
                    except Exception as e:
                        print(f'   QC figure failed: {e}', flush=True)
                except Exception:
                    print(f'{futs[f]}: FAILED\n{traceback.format_exc()}', flush=True)
    if len(rows) < 3:
        sys.exit(f'{len(rows)} subject(s) available: not enough for the analysis.')
    analyse(rows)
    print(f'Outputs: {OUT_DIR}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--analyse-only', action='store_true', help='only redo tables and figure from cached results')
    main(ap.parse_args().analyse_only)
