"""
Voxel-wise group analysis of static PET images ([18F]FDG SUV, activity, PVC images...) in a template space.

Input: the images exported by step 10 of the PET pipeline,
    <BIDS>/sub-<ID>/ses-<ses>/pet/templates/<template>/<base>_space-<template>_desc-<desc>_pet.nii.gz
(not skull-stripped; the analysis mask is built here), with the template atlases and brain mask written on
the same grid (templates/<template>/labels, masks).

Steps (run_group_analysis):
 1. images and subject table (participants.tsv / sessions.tsv columns: groups, covariates);
 2. common grid (images exported at the PET voxel size of each scanner are resampled, linearly, to the
    grid of the first image), coverage (voxels finite and > 0 in every image);
 3. intensity normalisation of every image, before smoothing (`normalization`):
      'global'       proportional scaling: image / its mean in the brain mask (Fox et al. 1988; SPM);
      reference      SUVR: image / its mean in a reference region, atlas labels or a mask, with the
                     `reference_signal` syntax of the PET pipeline (atlases/getlabels.py), e.g.
                     {'type': 'atlas', 'atlas': 'EDNIxCSCLR', 'level': 0, 'labels': 'default', 'name': 'Cer'}
                     (cerebellum) or {'type': 'mask', 'file': '<mask in the template space>', 'name': 'pons'};
      'none'         SUV (or activity) as they are.
    The normalising value of every image is written (inputs.tsv), and compared between groups: when it
    differs, the SUVR differences are partly differences of the reference (Borghammer et al. 2008;
    Yakushev et al. 2009): check this test before interpreting a two-group result;
 4. analysis mask = template brain mask (and `tissue_mask`, e.g. grey matter) x coverage x relative
    threshold: voxels below `rel_threshold` x the mean of the image in the brain, in any image, are
    excluded (SPM PET default 0.8): it removes the CSF and the background, whose spill-over is not a
    metabolic signal;
 5. smoothing restricted to the mask (`fwhm` mm; normalised convolution, Knutsson & Westin 1993): the
    signal outside the brain (CSF, scalp, background) does not enter the cortex;
 6. GLM (Statistics/Group_PET/_glm.py): 'one-sample', 'two-sample', 'regression' (covariate effect),
    'paired' (difference between two sessions of the same subjects), with covariates;
 7. inference: voxel-wise threshold `voxel_p` (two-sided, 0.001 by default as the fMRI group analyses)
    and cluster-level family-wise error `alpha`, by permutation (`cluster_method='permutation'`,
    Freedman-Lane, Winkler et al. 2014; Eklund et al. 2016) or with AFNI 3dClustSim on the ACF of the
    residuals (`'clustsim'`, as the fMRI 3dLMEr / 3dttest++ analyses). The permutation also gives
    voxel-level FWE p-values (maximum |t|);
 8. maps (NIfTI + PNG), cluster table, summary.json; surface projection: pet_surface.py.

Outputs: <BIDS>/Results/Grp_PET_<design>_desc-<desc>_norm-<normalization>[_<suffix>]/ (output_dir).
"""
import os
import re
import glob
import json

import numpy as np
import pandas as pd
import nibabel as nib
from scipy import ndimage, stats

from Statistics.Group_PET import _glm

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
opi = os.path.isfile

DESC_PRIORITY = ('SUV-mean', 'SUV', 'activity-mean', 'activity')
CLUSTER_METHODS = ('permutation', 'clustsim')


def _log(msg):
    print(msg, flush=True)


# ======================================================================================================= inputs
def _norm_id(x, prefix):
    x = str(x)
    return x[len(prefix):] if x.startswith(prefix) else x


def _entities(name):
    ent = {}
    for part in opb(name).split('_'):
        if '-' in part:
            k, v = part.split('-', 1)
            ent.setdefault(k, v)
    return ent


def _path_entities(path):
    """Entities of a BIDS file: those of its name, with sub / ses taken from the folders when the name does not
    carry them (a run named sub-X_trc-FDG_pet.nii.gz in sub-X/ses-Y/pet/: session Y)."""
    ent = _entities(path)
    for part in os.path.normpath(path).split(os.sep)[:-1]:
        for k in ('sub', 'ses'):
            if part.startswith(k + '-') and not ent.get(k):
                ent[k] = part[len(k) + 1:]
    return ent


def desc_alternatives(desc):
    """The desc asked, then its static / dynamic counterpart: step 7 writes desc-<x> for a static run and
    desc-<x>-mean (mean of the frames of SUV_window) for a dynamic run; both are the same quantity."""
    return [desc, desc[:-len('-mean')]] if desc.endswith('-mean') else [desc, desc + '-mean']


def pick_desc(df, desc):
    """For every (subject, session) of df (columns subject, session, desc), the row(s) of `desc`, or of its
    static / dynamic counterpart when `desc` is absent. Returns the selected rows and a log message (or '')."""
    alts = desc_alternatives(desc)
    parts, used = [], {}
    for _, g in df.groupby(['subject', 'session'], sort=False):
        for d in alts:
            sel = g[g.desc == d]
            if len(sel):
                parts.append(sel)
                used[d] = used.get(d, 0) + 1
                break
    out = pd.concat(parts) if parts else df.iloc[:0]
    msg = ''
    if any(d != desc for d in used):
        msg = (f'INFO: desc-{desc} resolved to ' + ', '.join(f'desc-{d} ({n} session(s))' for d, n in used.items())
               + ' (desc-<x>: static run; desc-<x>-mean: mean of the frames of SUV_window of a dynamic run)')
    return out.sort_index(), msg


def collect_images(bids_dir, template, desc=None, list_to_keep=None, list_to_remove=None):
    """
    Exported template-space images of every session. Returns a DataFrame: subject, session, run, base,
    image, data_path (BIDS session folder) and the `desc` used.
    desc: one desc of the PET pipeline ('SUV', 'SUV-mean', 'SUV-PVC', 'bratio'...), or None: the first of
    DESC_PRIORITY available for every session (the same desc for all the subjects). 'SUV' and 'SUV-mean'
    (static / dynamic run, see desc_alternatives) are interchangeable: each session uses the one it has.
    list_to_keep / list_to_remove: (ID, session) tuples, as in the launchers.
    """
    pats = [opj(bids_dir, 'sub-*', 'ses-*', 'pet', 'templates', template, f'*_space-{template}_desc-*_pet.nii.gz'),
            opj(bids_dir, 'sub-*', 'pet', 'templates', template, f'*_space-{template}_desc-*_pet.nii.gz')]
    rows = []
    for f in sorted(set(sum((glob.glob(p) for p in pats), []))):
        ent = _path_entities(f)
        d = ent.get('desc', '')
        if d.endswith('-SS') or 'grid' in d:
            continue
        base = '_'.join(p for p in opb(f).split('_') if p.split('-')[0] not in ('space', 'desc', 'pet.nii.gz'))
        rows.append(dict(subject=ent.get('sub'), session=ent.get('ses', ''), run=ent.get('run', ''), desc=d,
                         base=base, image=f, data_path=opd(opd(opd(opd(f))))))
    df = pd.DataFrame(rows, columns=['subject', 'session', 'run', 'desc', 'base', 'image', 'data_path'])
    if df.empty:
        raise FileNotFoundError(f'no exported PET image in {bids_dir}/sub-*/[ses-*/]pet/templates/{template}: run '
                                f'step 10 of the PET pipeline (template name = folder name under templates/)')
    df = _filter(df, list_to_keep, list_to_remove)
    keys = df[['subject', 'session']].drop_duplicates()
    if desc is None:
        for d in DESC_PRIORITY:
            have = df[df.desc == d][['subject', 'session']].drop_duplicates()
            if len(have) == len(keys):
                desc = d
                break
        if desc is None:
            raise ValueError(f'none of {DESC_PRIORITY} is available for every session: give desc '
                             f'(available: {sorted(df.desc.unique())})')
    out, msg = pick_desc(df, desc)
    if msg:
        _log(msg)
    out = out.reset_index(drop=True)
    missing = keys.merge(out[['subject', 'session']].drop_duplicates(), how='left', indicator=True)
    missing = missing[missing._merge == 'left_only']
    if out.empty:
        raise FileNotFoundError(f'no image with desc-{desc} or desc-{desc_alternatives(desc)[1]} '
                                f'(available: {sorted(df.desc.unique())})')
    if len(missing):
        _log(f'WARNING: desc-{desc} not found for {len(missing)} session(s): '
             + ', '.join(f'sub-{s}' + (f'_ses-{e}' if e else '') for s, e in zip(missing.subject, missing.session)))
    out.attrs['desc'] = desc
    return out


def _filter(df, list_to_keep, list_to_remove):
    def key(i, s):
        return _norm_id(i, 'sub-'), _norm_id(s, 'ses-')
    k = df.apply(lambda r: key(r.subject, r.session), axis=1)
    if list_to_keep:
        keep = {key(*t) if isinstance(t, (tuple, list)) else (key(t, '')[0], None) for t in list_to_keep}
        ok = [(i, s) in keep or (i, None) in keep for i, s in k]
        df = df[ok]
    if list_to_remove:
        rem = {key(*t) if isinstance(t, (tuple, list)) else (key(t, '')[0], None) for t in list_to_remove}
        ok = [not ((i, s) in rem or (i, None) in rem) for i, s in k[df.index]]
        df = df[ok]
    return df.reset_index(drop=True)


SUBJECT_ID_COLUMNS = ('participant_id', 'subject_id', 'subject', 'sub', 'participant', 'id', 'ID', 'Subject')
SESSION_ID_COLUMNS = ('session_id', 'session', 'ses', 'Session')


def _id_column(df, candidates, path):
    """Column of `df` holding the IDs: the first of `candidates` present (case-insensitive), else a first
    column whose values look like BIDS IDs (sub-xx / ses-xx)."""
    lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    first = df.columns[0]
    prefix = 'ses-' if candidates is SESSION_ID_COLUMNS else 'sub-'
    if df[first].astype(str).str.startswith(prefix).all():
        return first
    raise KeyError(f'{path}: no ID column (expected one of {list(candidates)}, or a first column of {prefix}<label> '
                   f'values); columns found: {list(df.columns)}')


def subject_table(images, bids_dir, columns=None):
    """
    images (collect_images) + the columns of sub-<ID>/sub-<ID>_sessions.tsv (by session) and of
    participants.tsv (by subject); a session value wins over a participant value. columns: the columns
    needed (KeyError if one is found nowhere). No column needed (one-sample or paired on sessions without
    covariates, groups given by group_assignment): the tables are not read.
    The ID column is participant_id / session_id (BIDS), or one of SUBJECT_ID_COLUMNS / SESSION_ID_COLUMNS.
    """
    t = images.copy()
    if not columns:
        return t
    part = opj(bids_dir, 'participants.tsv')
    if opi(part):
        p = pd.read_csv(part, sep='\t', dtype=str)
        idc = _id_column(p, SUBJECT_ID_COLUMNS, part)
        p['subject'] = p[idc].map(lambda x: _norm_id(x, 'sub-'))
        p = p.drop(columns=[idc] if idc != 'subject' else [])
        p = p.apply(lambda c: pd.to_numeric(c, errors='coerce') if c.name != 'subject' and _numeric(c) else c)
        t = t.merge(p, on='subject', how='left', suffixes=('', '_participants'))
    ses_rows = []
    for sub in t.subject.unique():
        f = opj(bids_dir, f'sub-{sub}', f'sub-{sub}_sessions.tsv')
        if opi(f):
            s = pd.read_csv(f, sep='\t', dtype=str)
            idc = _id_column(s, SESSION_ID_COLUMNS, f)
            s['session'] = s[idc].map(lambda x: _norm_id(x, 'ses-'))
            s = s.drop(columns=[idc] if idc != 'session' else [])
            s = s.apply(lambda c: pd.to_numeric(c, errors='coerce') if c.name != 'session' and _numeric(c) else c)
            s['subject'] = sub
            ses_rows.append(s)
    if ses_rows:
        s = pd.concat(ses_rows, ignore_index=True)
        t = t.merge(s, on=['subject', 'session'], how='left', suffixes=('', '_ses'))
        for c in [c for c in t.columns if c.endswith('_ses')]:
            base = c[:-4]
            t[base] = t[c].combine_first(t[base])
            t = t.drop(columns=[c])
    missing = [c for c in columns if c not in t.columns]
    if missing:
        avail = [c for c in t.columns if c not in ('subject', 'session', 'run', 'desc', 'base', 'image', 'data_path')]
        raise KeyError(f'column(s) {missing} found neither in participants.tsv nor in the sessions.tsv files '
                       f'(available: {avail})' + ('' if opi(part) else f'; no {part}'))
    return t


def _numeric(col):
    """Text column of a TSV read as str whose non-empty values are all numbers (n/a, empty: missing)."""
    v = col.dropna()
    v = v[~v.isin(['n/a', 'NA', 'NaN', ''])]
    return len(v) > 0 and pd.to_numeric(v, errors='coerce').notna().all()


def assign_groups(table, group_assignment, group_col='group'):
    """Group column from {group: [subject IDs]} (no participants.tsv needed). Subjects in no list: NaN
    (left out of a two-sample design); a subject in two lists or an ID with no image: error."""
    owner = {}
    for g, ids in group_assignment.items():
        for i in ids:
            i = _norm_id(i, 'sub-')
            if i in owner:
                raise ValueError(f'subject {i} is in groups {owner[i]!r} and {g!r}')
            owner[i] = str(g)
    absent = sorted(set(owner) - set(table.subject))
    if absent:
        raise ValueError(f'group_assignment: no image for subject(s) {absent} (subjects: {sorted(table.subject.unique())})')
    t = table.copy()
    t[group_col] = t.subject.map(owner)
    return t


# ======================================================================================================= images
def load_on_grid(files, ref_img=None):
    """4D array (X, Y, Z, n) of 3D images on the grid of ref_img (default: the first image). Images on
    another grid are resampled linearly (images exported at the PET resolution of different scanners)."""
    ref = nib.load(ref_img or files[0])
    out = np.zeros(ref.shape[:3] + (len(files),), np.float32)
    n_rs = 0
    for k, f in enumerate(files):
        im = nib.load(f)
        if im.ndim == 4 and im.shape[3] == 1:
            im = nib.Nifti1Image(np.asarray(im.dataobj)[..., 0], im.affine, im.header)
        if im.ndim != 3:
            raise ValueError(f'{f} is not a 3D image')
        if im.shape[:3] != ref.shape[:3] or not np.allclose(im.affine, ref.affine, atol=1e-3):
            im = resample_to_ref(im, ref, order=1)
            n_rs += 1
        out[..., k] = np.asarray(im.dataobj, np.float32)
    if n_rs:
        _log(f'INFO: {n_rs} image(s) resampled (linear) to the grid of {opb(ref_img or files[0])}')
    return out, ref


def resample_to_ref(img, ref, order=1):
    """img resampled to the grid of ref (3D), order 0 = nearest neighbour, 1 = trilinear; 0 outside the field
    of view. nibabel.processing (stable since nibabel 2.1) rather than nilearn.image.resample_to_img, whose
    keywords changed between versions (force_resample / copy_header exist only from nilearn 0.11)."""
    from nibabel.processing import resample_from_to
    return resample_from_to(img, (ref.shape[:3], ref.affine), order=order, mode='constant', cval=0.0)


def mask_on_grid(mask_file, ref, interp='nearest'):
    m = nib.load(mask_file)
    if m.ndim == 4:
        m = nib.Nifti1Image(np.asarray(m.dataobj)[..., 0], m.affine, m.header)
    if m.shape[:3] != ref.shape[:3] or not np.allclose(m.affine, ref.affine, atol=1e-3):
        if interp == 'nearest':
            m = resample_to_ref(m, ref, order=0)
            return np.asarray(m.dataobj) > 0
        m = resample_to_ref(m, ref, order=1)
        return np.asarray(m.dataobj) >= 0.5
    return np.asarray(m.dataobj) > 0


def labels_on_grid(atlas_file, level, ref):
    """Integer label volume of one level of an atlas, resampled (nearest) to the grid of ref."""
    a = nib.load(atlas_file)
    data = np.asarray(a.dataobj)
    if data.ndim == 4:
        data = data[..., level]
    elif level not in (0, None):
        raise ValueError(f'{atlas_file} has one level; level {level} requested')
    a3 = nib.Nifti1Image(np.rint(data).astype(np.int32), a.affine)
    if a3.shape != ref.shape[:3] or not np.allclose(a3.affine, ref.affine, atol=1e-3):
        a3 = resample_to_ref(a3, ref, order=0)
    return np.rint(np.asarray(a3.dataobj)).astype(np.int64)


def template_files(images):
    """Brain mask and atlases written by the PET export next to the first image (templates/<template>/)."""
    d = opd(images.image.iloc[0])
    masks = sorted(glob.glob(opj(d, 'masks', '*.nii.gz')))
    brain = [m for m in masks if re.search(r'(_mask|brain_mask|desc-brain)', opb(m)) and
             not re.search(r'(Gray|White|Vent|gray|white|csf|erod)', opb(m))]
    atlases = {re.sub(r'_dseg$', '', pre.split('_seg-')[-1]): f for f in
               sorted(glob.glob(opj(d, 'labels', '*_seg-*_dseg.nii.gz'))) for pre in [opb(f)[:-7]]}
    return dict(brain=brain[0] if brain else (masks[0] if masks else None), masks=masks, atlases=atlases)


# ======================================================================================================= normalisation
def reference_mask(normalization, ref, files, label_tables=None):
    """Boolean volume of the reference region of `normalization` (dict) on the grid of ref.
    files: template_files() (atlases of the template space)."""
    from atlases import getlabels
    r = dict(normalization)
    if r.get('type') == 'mask':
        f = r.get('file') or (files['brain'] if r.get('name') == 'brain' else None)
        if not f or not opi(f):
            raise FileNotFoundError(f"reference mask {r}: give 'file' (a mask in the template space)")
        return mask_on_grid(f, ref, interp='linear')
    if r.get('type') != 'atlas':
        raise ValueError(f"normalization={normalization!r}: expected 'global', 'none' or a reference region "
                         f"{{'type': 'atlas', ...}} / {{'type': 'mask', 'file': ...}}")
    atlas = r.get('atlas')
    if atlas not in files['atlases']:
        raise FileNotFoundError(f"atlas {atlas!r} not exported in the template space (available: "
                                f"{sorted(files['atlases'])}); check atlas_followers / selected_atlases")
    level = int(r.get('level', 0))
    lab = labels_on_grid(files['atlases'][atlas], level, ref)
    present = sorted(int(x) for x in np.unique(lab) if x > 0)
    table = (label_tables or {}).get(atlas)
    if isinstance(r.get('labels'), str) and r['labels'] != 'default' and table is None:
        raise ValueError(f"labels given as a name pattern ({r['labels']!r}): give label_tables={{'{atlas}': "
                         f"<label table>}} or the label numbers")
    res = getlabels.resolve(r, {'name': atlas, 'level': level}, present, table, None)
    m = np.isin(lab, res['labels'])
    if not m.any():
        raise ValueError(f'reference region {res}: no voxel in the template space')
    return m


def normalise(data, method, brain, ref_mask=None):
    """Divides every image by its mean in the brain (global) or in the reference region. Returns
    (normalised data, values) — values: the normalising value of each image (NaN for 'none')."""
    n = data.shape[-1]
    if method == 'none':
        return data, np.full(n, np.nan)
    region = brain if method == 'global' else ref_mask & brain
    if region.sum() == 0:
        raise ValueError('the normalisation region has no voxel inside the brain and the coverage of all images')
    values = np.array([float(np.mean(data[..., k][region])) for k in range(n)])
    if np.any(~np.isfinite(values)) or np.any(values <= 0):
        raise ValueError(f'non-positive normalising value(s): {values}')
    return data / values[None, None, None, :], values


# ======================================================================================================= mask, smoothing
def analysis_mask(data, brain, rel_threshold=0.8, tissue=None):
    """brain x tissue x coverage (finite and > 0 in every image) x relative threshold (value >= rel_threshold x
    the mean of the image in brain x coverage, in every image; None: not applied). Returns (mask, report)."""
    cover = np.all(np.isfinite(data) & (data > 0), axis=-1)
    m = brain & cover
    rep = dict(brain=int(brain.sum()), covered=int(m.sum()))
    if tissue is not None:
        m &= tissue
        rep['tissue'] = int(m.sum())
    if rel_threshold:
        base = brain & cover
        for k in range(data.shape[-1]):
            v = data[..., k]
            m &= v >= rel_threshold * float(v[base].mean())
        rep['relative_threshold'] = int(m.sum())
    if m.sum() == 0:
        raise ValueError(f'empty analysis mask ({rep}): check the brain mask, the coverage and rel_threshold')
    return m, rep


def smooth_in_mask(data, mask, fwhm_mm, zooms):
    """Gaussian smoothing restricted to the mask (normalised convolution): sum of the weights inside the mask
    only, so that values outside the brain do not enter it. fwhm_mm: scalar or (x, y, z) mm; 0/None: none."""
    if not fwhm_mm:
        return np.where(mask[..., None], data, 0).astype(np.float32)
    fwhm = np.broadcast_to(np.asarray(fwhm_mm, float), (3,))
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2))) / np.asarray(zooms[:3], float)
    w = ndimage.gaussian_filter(mask.astype(np.float64), sigma, mode='constant')
    out = np.zeros(data.shape, np.float32)
    for k in range(data.shape[-1]):
        num = ndimage.gaussian_filter(np.where(mask, data[..., k], 0).astype(np.float64), sigma, mode='constant')
        with np.errstate(invalid='ignore', divide='ignore'):
            out[..., k] = np.where(mask, num / w, 0)
    return out


# ======================================================================================================= clusters
def cluster_table(z3d, t3d, supra_thr, ref, mask, null_max=None, labels=None, label_names=None):
    """Clusters of |t| >= supra_thr (positive and negative separately, face connectivity): size, volume,
    peak z and its world coordinates, region of the peak, cluster FWE p (permutation)."""
    vox = float(np.prod(ref.header.get_zooms()[:3]))
    rows = []
    for sign, sup in ((1, (t3d >= supra_thr) & mask), (-1, (t3d <= -supra_thr) & mask)):
        lab, n = ndimage.label(sup)
        for c in range(1, n + 1):
            vox_c = lab == c
            size = int(vox_c.sum())
            zc = np.where(vox_c, z3d * sign, -np.inf)
            pk = np.unravel_index(int(np.argmax(zc)), z3d.shape)
            xyz = nib.affines.apply_affine(ref.affine, pk)
            row = dict(sign='+' if sign > 0 else '-', n_voxels=size, volume_mm3=round(size * vox, 3),
                       peak_z=round(float(z3d[pk]), 4), peak_x=round(float(xyz[0]), 2),
                       peak_y=round(float(xyz[1]), 2), peak_z_mm=round(float(xyz[2]), 2))
            if null_max is not None:
                row['p_FWE_cluster'] = float(_glm.fwe_p([size], null_max)[0])
            if labels is not None:
                l = int(labels[pk])
                row['peak_label'] = l
                row['peak_region'] = (label_names or {}).get(l, '' if l == 0 else f'label_{l}')
                inside = labels[vox_c]
                vals, cnt = np.unique(inside[inside > 0], return_counts=True)
                row['regions'] = '; '.join(f"{(label_names or {}).get(int(v), f'label_{int(v)}')} ({c_})"
                                           for v, c_ in sorted(zip(vals, cnt), key=lambda x: -x[1])[:5])
            rows.append((row, vox_c))
    rows.sort(key=lambda r: -r[0]['n_voxels'])
    table = pd.DataFrame([r for r, _ in rows])
    return table, [m for _, m in rows]


# ======================================================================================================= plots
def plot_map(path, bg, title, out_png, threshold=None, cut_coords=7, display_mode='z', symmetric=True,
             cmap='cold_hot'):
    """Slices of a map over the template (nilearn); skipped when the map is empty."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        from nilearn import plotting
        data = np.asarray(nib.load(path).dataobj)
        if not np.any(data[np.isfinite(data)]):
            _log(f'  [plot] {opb(path)} is empty: no figure')
            return None
        vmax = float(np.nanpercentile(np.abs(data[data != 0]), 99.5)) or None
        disp = plotting.plot_stat_map(path, bg_img=bg, display_mode=display_mode, cut_coords=cut_coords,
                                      threshold=threshold, title=title, colorbar=True, vmax=vmax,
                                      symmetric_cbar=symmetric, cmap=cmap)
        disp.savefig(out_png, dpi=200)
        disp.close()
        return out_png
    except Exception as e:  # figures must not stop the analysis
        _log(f'  [plot] WARNING: {opb(path)}: {e}')
        return None


# ======================================================================================================= main
def _needed_columns(design, group_col, groups, covariates, paired_col, group_assignment):
    """(group_col, groups, columns to read from participants.tsv / sessions.tsv)."""
    if group_assignment:
        if design != 'two-sample':
            raise ValueError("group_assignment is for design='two-sample'")
        group_col = group_col or 'group'
        groups = list(groups or group_assignment)
        if len(groups) != 2:
            raise ValueError(f'group_assignment: two groups expected, got {groups}')
    needed = list(covariates or []) + ([group_col] if design == 'two-sample' and not group_assignment else []) + \
        ([paired_col] if design == 'paired' and paired_col not in ('session', 'subject') else [])
    return group_col, groups, needed


def run_group_analysis(bids_dir, template, design, desc=None, group_col=None, groups=None, covariates=None,
                       test_covariate=None, paired_col='session', paired_levels=None, normalization='global',
                       label_tables=None, brain_mask=None, tissue_mask=None, rel_threshold=0.8, fwhm=None,
                       aggregate='subject', voxel_p=0.001, alpha=0.05, cluster_method='permutation',
                       n_perm=5000, seed=0, sing_afni='', list_to_keep=None, list_to_remove=None,
                       output_dir=None, output_suffix=None, background=None, peak_atlas=None, peak_level=0,
                       group_assignment=None):
    """
    Voxel-wise group analysis of PET images in a template space (see the module docstring).

    Parameters
    ----------
    bids_dir : str
        BIDS folder processed by the PET pipeline (step 10 run).
    template : str
        Template space: folder name under pet/templates/ (the `reference` of the PET launcher, e.g. 'EDNiX',
        or the study template). Surfaces (pet_surface.py) need the reference template of the atlas library.
    design : str
        'one-sample', 'two-sample', 'regression' or 'paired'.
    desc : str or None
        Image of the PET pipeline ('SUV', 'SUV-mean', 'SUV-PVC', 'activity'...); None: first of
        DESC_PRIORITY available for every session.
    group_col, groups : str, [str, str]
        Two-sample: column of participants.tsv / sessions.tsv and the two groups (contrast groups[0] - groups[1]).
    group_assignment : dict or None
        Two-sample without participants.tsv: {group: [subject IDs]} (group_col defaults to 'group', groups to
        the two keys). Subjects in no list are left out.
    covariates : list of str or None
        Columns of participants.tsv / sessions.tsv (numeric: centred; text: categorical). None / []: no
        covariate; when no column is needed the tables are not read.
    test_covariate : str
        Regression: the covariate tested (default: the first one).
    paired_col, paired_levels : str, [a, b]
        Paired: column that distinguishes the two conditions (default the session) and their values;
        the difference a - b of every subject having both is tested (sign flipping).
    normalization : 'global' | 'none' | dict
        Intensity normalisation (module docstring). Default 'global' (proportional scaling).
    label_tables : dict
        {atlas: label table (DataFrame label_id, region_name)} for names in the reference region and the
        cluster table (Tools.atlas_labels.parse_label_file).
    brain_mask, tissue_mask : str
        Brain mask (default: the mask exported in templates/<template>/masks) and an optional tissue mask
        (e.g. the grey-matter mask of the template) intersected with it.
    rel_threshold : float or None
        Relative threshold of the analysis mask (0.8 = SPM PET default; None: off).
    fwhm : float, (fx, fy, fz) or None
        Smoothing FWHM in mm, inside the mask. None: twice the voxel size (no smoothing: 0). The final
        resolution is about sqrt(PET resolution^2 + fwhm^2); scale it with the brain (8-12 mm in humans
        for FDG, Signorini et al. 1999; a few mm in macaques).
    aggregate : 'subject' | 'session'
        Images of the same subject averaged (after normalisation) into one observation ('subject',
        default: the observations of a GLM must be independent) or kept per session ('session': only
        with one session per subject). Runs of one session are always averaged.
    voxel_p, alpha : float
        Cluster-forming two-sided voxel threshold and family-wise error rate.
    cluster_method : 'permutation' | 'clustsim'
        Permutation (default; no AFNI needed) or AFNI 3dClustSim with the ACF of the residuals
        (sing_afni: prefix of the AFNI container, '' for a native AFNI).
    n_perm, seed : int
        Number of permutations / sign flips (exhaustive when there are fewer) and random seed.
    list_to_keep, list_to_remove : list of (ID, session)
        Selection, as in the launchers.
    output_dir, output_suffix : str
        Output folder (default <bids_dir>/Results/Grp_PET_<design>_desc-<desc>_norm-<norm>[_<suffix>]).
    background : str
        Background image of the figures (e.g. the skull-stripped template; default: the group mean image).
    peak_atlas, peak_level : str, int
        Atlas of the template space used to name the clusters (default: first atlas exported).

    Returns
    -------
    dict: output_dir, contrast, maps (paths), clusters (DataFrame), summary (dict).
    """
    if cluster_method not in CLUSTER_METHODS:
        raise ValueError(f'cluster_method={cluster_method!r}: expected one of {CLUSTER_METHODS}')
    if aggregate not in ('subject', 'session'):
        raise ValueError(f"aggregate={aggregate!r}: expected 'subject' or 'session'")
    covariates = list(covariates or [])

    # ------------------------------------------------------------------------------------------- 1 inputs
    images = collect_images(bids_dir, template, desc, list_to_keep, list_to_remove)
    desc = images.attrs.get('desc') or images.desc.iloc[0]
    group_col, groups, needed = _needed_columns(design, group_col, groups, covariates, paired_col, group_assignment)
    table = subject_table(images, bids_dir, needed)
    if group_assignment:
        table = assign_groups(table, group_assignment, group_col)
    files = template_files(images)
    norm_tag = normalization if isinstance(normalization, str) else (normalization.get('name') or 'ref')
    out = output_dir or opj(bids_dir, 'Results', f'Grp_PET_{design}_desc-{desc}_norm-{norm_tag}'
                            + (f'_{output_suffix}' if output_suffix else ''))
    os.makedirs(out, exist_ok=True)
    _log(f'PET group analysis: {len(images)} image(s), desc-{desc}, template {template}, design {design}\n'
         f'  output: {out}')

    # ------------------------------------------------------------------------------------------- 2 grid
    data, ref = load_on_grid(list(table.image))
    zooms = ref.header.get_zooms()[:3]
    brain_file = brain_mask or files['brain']
    if not brain_file:
        raise FileNotFoundError(f'no brain mask in {opd(images.image.iloc[0])}/masks: give brain_mask')
    _log(f'  brain mask: {brain_file}')
    brain = mask_on_grid(brain_file, ref, interp='linear')
    cover = np.all(np.isfinite(data) & (data > 0), axis=-1)
    brain_cov = brain & cover

    # ------------------------------------------------------------------------------------------- 3 normalisation
    ref_mask = None
    if isinstance(normalization, dict):
        ref_mask = reference_mask(normalization, ref, files, label_tables)
        method = 'reference'
    elif normalization in ('global', 'none'):
        method = normalization
    else:
        raise ValueError(f"normalization={normalization!r}: expected 'global', 'none' or a reference dict")
    data, values = normalise(np.nan_to_num(data, nan=0.0), method, brain_cov, ref_mask)
    table['normalising_value'] = values
    if ref_mask is not None:
        table['reference_voxels'] = int((ref_mask & brain_cov).sum())

    # ------------------------------------------------------------------------------------------- 4 mask
    tissue = mask_on_grid(tissue_mask, ref, interp='linear') if tissue_mask else None
    mask, mrep = analysis_mask(data, brain, rel_threshold, tissue)
    _log(f'  analysis mask (voxels): {mrep}')

    # ------------------------------------------------------------------------------------------- 5 smoothing
    if fwhm is None:
        fwhm = 2 * float(min(zooms))
    sm = smooth_in_mask(data, mask, fwhm, zooms)
    del data

    # ------------------------------------------------------------------------------------------- observations
    table['_col'] = np.arange(len(table))
    obs, Y = _observations(table, sm, mask, design, aggregate, paired_col, paired_levels)
    obs.drop(columns=['_col'], errors='ignore').to_csv(opj(out, 'inputs.tsv'), sep='\t', index=False)
    _write(sm[..., list(table._col)] * 1, ref, opj(out, 'inputs_4D.nii.gz'))
    _write(mask.astype(np.uint8), ref, opj(out, 'mask_group.nii.gz'))

    # ------------------------------------------------------------------------------------------- 6 GLM
    design_glm = 'one-sample' if design == 'paired' else design
    X, c, names, rows = _glm.design_matrix(obs, design_glm, group_col, groups, covariates, test_covariate)
    Y = Y[rows]
    obs = obs[rows].reset_index(drop=True)
    cname = _glm.contrast_name(design, c, names, groups)
    glm = _glm.GLM(X, c)
    cope, t, resid = glm.fit(Y)
    z = _glm.t_to_z(t, glm.df)
    _log(f'  GLM: {X.shape[0]} observations, columns {names}, contrast {cname}, df = {glm.df}')

    maps = {}
    for k, v in (('cope', cope), ('tstat', t), ('zstat', z)):
        maps[k] = _write(_vol(v, mask), ref, opj(out, f'{cname}_{k}.nii.gz'))
    maps.update(_descriptive(obs, Y, mask, ref, out, design, group_col, groups))

    # ------------------------------------------------------------------------------------------- 7 inference
    summary = dict(design=design, contrast=cname, columns=names, n=int(X.shape[0]), df=int(glm.df), desc=desc,
                   template=template, normalization=normalization, fwhm_mm=fwhm, rel_threshold=rel_threshold,
                   mask_voxels=int(mask.sum()), mask_report=mrep, voxel_p=voxel_p, alpha=alpha,
                   cluster_method=cluster_method, aggregate=aggregate)
    thr = _glm.t_threshold(voxel_p, glm.df)
    null_max = None
    if cluster_method == 'permutation':
        sign_flip = design_glm == 'one-sample'
        perm = _glm.permutation_test(glm, Y, mask, voxel_p, n_perm, sign_flip, seed,
                                     progress=lambda k, m: _log(f'  permutation {k}/{m}'))
        null_max = perm['max_cluster']
        k_min = _glm.critical_value(null_max, alpha, integer=True)
        p_vox = _glm.fwe_p(np.abs(t), perm['max_t'])
        maps['logp_fwe_voxel'] = _write(_vol(-np.log10(np.maximum(p_vox, 1.0 / perm['n_perm'])), mask), ref,
                                        opj(out, f'{cname}_logp-FWE-voxel.nii.gz'))
        summary.update(n_perm=perm['n_perm'], exhaustive=perm['exhaustive'], cluster_k_min=k_min,
                       t_FWE_voxel=_glm.critical_value(perm['max_t'], alpha),
                       min_p_attainable=1.0 / perm['n_perm'],
                       method_reference='Freedman & Lane 1983; Winkler et al. 2014 NeuroImage 92:381')
        np.savetxt(opj(out, f'{cname}_null_max_cluster.txt'), null_max, fmt='%d')
        if 1.0 / perm['n_perm'] > alpha:
            _log(f'WARNING: {perm["n_perm"]} permutations only: the smallest p is {1.0 / perm["n_perm"]:.3f} > '
                 f'alpha = {alpha}: nothing can be significant (too few subjects)')
    else:
        from Statistics.Group_fMRI._Group_anal_3dLMEr_SBA import _clustsim_acf, _extract_cluster_threshold
        rpath = _write(np.stack([_vol(r, mask) for r in resid], -1), ref, opj(out, f'{cname}_residuals.nii.gz'))
        table_cs, acf = _clustsim_acf(sing_afni, opj(out, 'mask_group.nii.gz'), rpath, opj(out, 'ClustSim'))
        # positive and negative clusters are formed separately here: AFNI 'bi-sided' table
        table_cs = table_cs.replace('NN1_2sided', 'NN1_bisided')
        k_min = _extract_cluster_threshold(sing_afni, table_cs, voxel_p, alpha)
        summary.update(cluster_k_min=k_min, acf=acf, method_reference='Cox et al. 2017 Brain Connect 7:152')
    summary['t_cluster_forming'] = thr

    # ------------------------------------------------------------------------------------------- 8 outputs
    t3d, z3d = _vol(t, mask), _vol(z, mask)
    labels, names_map = None, None
    atl = peak_atlas or (sorted(files['atlases'])[0] if files['atlases'] else None)
    if atl and atl in files['atlases']:
        labels = labels_on_grid(files['atlases'][atl], peak_level, ref)
        lt = (label_tables or {}).get(atl)
        names_map = dict(zip(lt['label_id'], lt['region_name'])) if lt is not None else None
    clusters, cmasks = cluster_table(z3d, t3d, thr, ref, mask, null_max, labels, names_map)
    surv = np.zeros(mask.shape, bool)
    if len(clusters):
        clusters['significant'] = clusters['n_voxels'] >= k_min
        for keep, cm in zip(clusters['significant'], cmasks):
            if keep:
                surv |= cm
    clusters.to_csv(opj(out, f'{cname}_clusters.tsv'), sep='\t', index=False)
    maps['zstat_thresholded'] = _write(np.where(surv, z3d, 0), ref, opj(out, f'{cname}_zstat_thresholded.nii.gz'))
    maps['cope_thresholded'] = _write(np.where(surv, _vol(cope, mask), 0), ref,
                                      opj(out, f'{cname}_cope_thresholded.nii.gz'))
    summary.update(n_clusters=int(len(clusters)), n_significant=int(clusters['significant'].sum()) if len(clusters)
                   else 0)
    if design == 'two-sample' and method != 'none':
        summary['normalising_value_test'] = _reference_test(obs, group_col, groups)
    with open(opj(out, f'{cname}_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2, default=lambda o: o if not isinstance(o, np.generic) else o.item())
    _log(f'  voxel threshold |t| >= {thr:.3f} (p < {voxel_p} two-sided), clusters >= {k_min} voxels '
         f'(FWE {alpha}, {cluster_method}): {summary["n_significant"]} of {summary["n_clusters"]} cluster(s)')

    bg = background or maps.get('mean_all')
    plot_map(maps['zstat'], bg, f'{cname} z (unthresholded)', maps['zstat'].replace('.nii.gz', '.png'))
    plot_map(maps['zstat_thresholded'], bg, f'{cname} z (p<{voxel_p}, k>={k_min:g}, FWE {alpha}, {cluster_method})',
             maps['zstat_thresholded'].replace('.nii.gz', '.png'), threshold=1e-6)
    return dict(output_dir=out, contrast=cname, maps=maps, clusters=clusters, summary=summary)


def _observations(table, sm, mask, design, aggregate, paired_col, paired_levels):
    """One row per observation and Y (n x voxels in the mask)."""
    flat = sm[mask]                                              # voxels x images
    # runs of a session first
    ses = table.groupby(['subject', 'session'], sort=True, as_index=False)
    rows, cols = [], []
    for _, g in ses:
        r = g.iloc[0].drop(labels=['image', 'run', 'base']).to_dict()
        r['n_images'] = len(g)
        r['normalising_value'] = float(np.nanmean(g['normalising_value']))
        rows.append(r)
        cols.append(flat[:, list(g._col)].mean(1))
    obs = pd.DataFrame(rows)
    Y = np.stack(cols, 0)
    if design == 'paired':
        if not paired_levels or len(paired_levels) != 2:
            raise ValueError("design='paired' needs paired_levels=[condition a, condition b] (values of paired_col)")
        a, b = (str(x) for x in paired_levels)
        key = obs[paired_col].astype(str)
        ia = {s: i for i, s in enumerate(obs.subject) if key.iloc[i] == a}
        ib = {s: i for i, s in enumerate(obs.subject) if key.iloc[i] == b}
        subs = sorted(set(ia) & set(ib))
        if len(subs) < 3:
            raise ValueError(f'{len(subs)} subject(s) with both {paired_col}={a} and {b}')
        Y = np.stack([Y[ia[s]] - Y[ib[s]] for s in subs], 0)
        obs = obs.loc[[ia[s] for s in subs]].reset_index(drop=True)
        obs['pair'] = f'{a}-{b}'
        return obs, Y
    if aggregate == 'subject' and obs.subject.duplicated().any():
        g = obs.groupby('subject', sort=True)
        new_rows, new_Y = [], []
        for _, gi in g:
            r = gi.iloc[0].to_dict()
            r['session'] = ','.join(str(x) for x in gi.session)
            r['n_images'] = int(gi.n_images.sum())
            for c in gi.columns:
                if c not in ('subject', 'session', 'n_images', 'desc', 'data_path') and \
                        pd.api.types.is_numeric_dtype(gi[c]):
                    r[c] = float(gi[c].mean())
            new_rows.append(r)
            new_Y.append(Y[gi.index].mean(0))
        obs, Y = pd.DataFrame(new_rows), np.stack(new_Y, 0)
        _log(f'  sessions of the same subject averaged: {len(obs)} subjects')
    elif aggregate == 'session' and obs.subject.duplicated().any():
        raise ValueError("aggregate='session' with several sessions per subject: the observations are not "
                         "independent; use aggregate='subject' or design='paired'")
    return obs, Y


def _descriptive(obs, Y, mask, ref, out, design, group_col, groups):
    maps = {}
    sets = [('all', np.ones(len(obs), bool))]
    if design == 'two-sample':
        sets += [(str(g), (obs[group_col].astype(str) == str(g)).to_numpy()) for g in groups]
    for name, sel in sets:
        if sel.sum() == 0:
            continue
        m = Y[sel].mean(0)
        maps[f'mean_{name}'] = _write(_vol(m, mask), ref, opj(out, f'group-{name}_mean.nii.gz'))
        if sel.sum() > 1:
            sd = Y[sel].std(0, ddof=1)
            maps[f'sd_{name}'] = _write(_vol(sd, mask), ref, opj(out, f'group-{name}_sd.nii.gz'))
    if design == 'two-sample':
        a, b = (Y[(obs[group_col].astype(str) == str(g)).to_numpy()].mean(0) for g in groups)
        with np.errstate(divide='ignore', invalid='ignore'):
            pc = np.where(b != 0, 100 * (a - b) / b, 0)
        maps['percent_difference'] = _write(_vol(pc, mask), ref,
                                            opj(out, f'{groups[0]}-vs-{groups[1]}_percent-difference.nii.gz'))
    return maps


def _reference_test(obs, group_col, groups):
    """Welch t-test of the normalising value between the two groups (bias of the normalisation)."""
    va = obs.loc[obs[group_col].astype(str) == str(groups[0]), 'normalising_value'].to_numpy(float)
    vb = obs.loc[obs[group_col].astype(str) == str(groups[1]), 'normalising_value'].to_numpy(float)
    if len(va) < 2 or len(vb) < 2 or not np.all(np.isfinite(np.r_[va, vb])):
        return None
    tt = stats.ttest_ind(va, vb, equal_var=False)
    res = dict(mean_a=float(va.mean()), mean_b=float(vb.mean()), ratio=float(va.mean() / vb.mean()),
               t=float(tt.statistic), p=float(tt.pvalue))
    if tt.pvalue < 0.05:
        _log(f'WARNING: the normalising value differs between {groups[0]} and {groups[1]} (ratio '
             f'{res["ratio"]:.3f}, p = {tt.pvalue:.3g}): the normalised differences are biased by '
             f'{100 * (1 / res["ratio"] - 1):+.1f} % (Borghammer et al. 2008); choose another reference region')
    return res


def _vol(v, mask):
    out = np.zeros(mask.shape, np.float32)
    out[mask] = v
    return out


def _write(arr, ref, path):
    img = nib.Nifti1Image(np.asarray(arr, np.float32 if arr.dtype != np.uint8 else np.uint8), ref.affine)
    img.header.set_xyzt_units('mm')
    nib.save(img, path)
    return path
