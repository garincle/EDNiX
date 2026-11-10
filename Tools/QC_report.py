"""
EDNiX quality-control reports, three levels, written in <BIDS>/QC/ (the existing QC/<check>/ folders are left as
they are; the reports only read them):

  1. one recap per subject / session      QC/subjects/sub-<ID>[_ses-<S>].html
       processing status, errors / warnings of the logs, metrics (anat, surfaces, every fMRI run, FC specificity,
       PET motion) and every QC image of that subject (QC/<check>/ images, anat/QC_anat, func/QC, pet/QC, ct/QC)
  2. one page gathering all the subjects  QC/QC_all_subjects.html
       one row per subject / session with the key metrics, outliers (robust z > 3.5, median / MAD) highlighted,
       link to each recap, then one thumbnail grid per QC check
  3. the dataset (BIDS) level              QC/QC_dataset.html + QC/EDNiX_BIDS_summary.xlsx
       subjects, sessions, images per modality, acquisition parameters read from the raw BIDS files (TR from the
       sidecar and from the NIfTI header, flagged when they disagree, volumes, duration, voxel size, field...),
       processing status and all the QC metrics, one sheet each.

The surface QC (surface_qc) is the figure of the surface tests (white / pial contours on orig.mgz on 4 planes,
thickness histogram in real mm, distance of the surfaces to the aseg boundaries). It only reads the FreeSurfer
folder, so it works with any surface reconstruction parameters. The anatomical pipeline calls it at step 16.

    from Tools import QC_report
    QC_report.build(bids_dir)                       # the three levels for one dataset
    QC_report.study_summary([bids1, bids2], out)    # one Excel file across datasets

Everything is read-only on the data; a report that cannot be built is reported, never fatal.
"""
import os
import re
import glob
import json
import html
import shutil
import datetime
import traceback

import numpy as np

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
opi = os.path.isfile
ope = os.path.exists

IMG_EXT = ('.png', '.jpg', '.jpeg', '.gif', '.svg')
HEMIS = ('lh', 'rh')
CSF_CODES = [4, 5, 43, 44, 14, 15, 24, 72, 30, 62, 31, 63]          # preFS.CSF_edit
Z_THRESH = 3.5
# metrics shown in the all-subjects table (all the others are in the recaps and the Excel file)
KEY_METRICS = {
    'func': ['mean_fd', 'censor_fraction', 'TSNR', 'mean_dvars', 'gcor', 'avg_snr_gray', 'cc', 'nmi_brain'],
    'surface': ['thick_median_mm', 'white_within_1vox', 'pial_within_1vox', 'n_vertices_white'],
}
ERROR_RE = re.compile(r'\bERROR\b|Traceback \(most recent call last\)|exited with code [1-9]|^\s*(?:\*\* )?error:|'
                      r'^\s*\w*(?:Error|Exception)\b(?::|$)|\bFATAL\b|Segmentation fault', re.MULTILINE)
WARN_RE = re.compile(r'\bWARNING\b|^\s*(?:\*\* )?warning:', re.MULTILINE)

CSS = """
body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:24px;color:#1d1d1f;background:#fff}
h1{font-size:22px;margin:0 0 4px} h2{font-size:17px;margin:28px 0 8px;border-bottom:1px solid #ddd;padding-bottom:4px}
h3{font-size:14px;margin:18px 0 6px} .muted{color:#6e6e73;font-size:12px}
table{border-collapse:collapse;font-size:12px;margin:6px 0} td,th{border:1px solid #ddd;padding:3px 7px;text-align:left}
th{background:#f5f5f7;position:sticky;top:0} td.num{text-align:right;font-variant-numeric:tabular-nums}
td.out{background:#ffd7d5;font-weight:600} td.ok{background:#e3f6e8} td.no{background:#fbeaea} td.na{color:#aaa}
.grid{display:flex;flex-wrap:wrap;gap:8px} .card{width:260px;font-size:11px;word-break:break-all}
.card img{width:260px;border:1px solid #ddd;background:#000} .big img{width:520px}
.badge{display:inline-block;padding:1px 6px;border-radius:8px;font-size:11px;margin-right:4px}
.err{background:#ffd7d5} .warn{background:#fff1c2} .fine{background:#e3f6e8}
a{color:#0a58ca;text-decoration:none} a:hover{text-decoration:underline}
"""


# =====================================================================================================================
# helpers
# =====================================================================================================================
def _esc(x):
    return html.escape(str(x))


def _rel(path, start):
    return os.path.relpath(path, start).replace(os.sep, '/')


def _load_json(p):
    try:
        with open(p) as f:
            return json.load(f)
    except Exception:
        return None


def _flatten(d, prefix=''):
    """Scalar values of a (nested) dict; lists are summarised by their mean when numeric."""
    out = {}
    if not isinstance(d, dict):
        return out
    for k, v in d.items():
        key = f'{prefix}{k}'
        if isinstance(v, dict):
            out.update(_flatten(v, key + '.'))
        elif isinstance(v, (int, float, str, bool)) or v is None:
            out[key] = v
        elif isinstance(v, (list, tuple)) and v and all(isinstance(x, (int, float)) for x in v):
            out[key + '_mean'] = float(np.mean(v))
    return out


def _fmt(v):
    if isinstance(v, float):
        if not np.isfinite(v):
            return 'nan'
        return f'{v:.4g}'
    return '' if v is None else str(v)


def _robust_z(values):
    v = np.asarray([x if isinstance(x, (int, float)) and not isinstance(x, bool) else np.nan for x in values],
                   dtype=float)
    ok = np.isfinite(v)
    z = np.full(len(v), np.nan)
    if ok.sum() < 5:
        return z
    med = np.median(v[ok])
    mad = np.median(np.abs(v[ok] - med)) * 1.4826
    if mad > 0:
        z[ok] = (v[ok] - med) / mad
    return z


def _page(title, body, subtitle=''):
    now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
    return (f'<!doctype html><html><head><meta charset="utf-8"><title>{_esc(title)}</title>'
            f'<meta name="viewport" content="width=device-width,initial-scale=1"><style>{CSS}</style></head><body>'
            f'<h1>{_esc(title)}</h1><div class="muted">{_esc(subtitle)} · generated {now} by Tools/QC_report.py</div>'
            f'{body}</body></html>')


def _table(rows, columns, outliers=None, links=None):
    """rows: list of dicts. outliers: {(i, col)}; links: {i: href} applied to the first column."""
    outliers = outliers or set()
    links = links or {}
    h = ['<table><tr>' + ''.join(f'<th>{_esc(c)}</th>' for c in columns) + '</tr>']
    for i, r in enumerate(rows):
        cells = []
        for j, c in enumerate(columns):
            v = r.get(c, '')
            cls = 'num' if isinstance(v, (int, float)) and not isinstance(v, bool) else ''
            if (i, c) in outliers:
                cls += ' out'
            if isinstance(v, bool):
                cls, v = ('ok', 'yes') if v else ('no', 'no')
            txt = _esc(_fmt(v))
            if j == 0 and i in links:
                txt = f'<a href="{_esc(links[i])}">{txt}</a>'
            cells.append(f'<td class="{cls.strip()}">{txt}</td>')
        h.append('<tr>' + ''.join(cells) + '</tr>')
    h.append('</table>')
    return ''.join(h)


def _images_html(paths, start, big=False):
    if not paths:
        return '<p class="muted">none</p>'
    cards = []
    for p in paths:
        r = _rel(p, start)
        cards.append(f'<div class="card{" big" if big else ""}"><a href="{_esc(r)}" target="_blank">'
                     f'<img loading="lazy" src="{_esc(r)}"></a><br>{_esc(opb(p))}</div>')
    return '<div class="grid">' + ''.join(cards) + '</div>'


# =====================================================================================================================
# dataset layout
# =====================================================================================================================
def sessions(bids_dir):
    """[(ID, Session, data_path)] from the BIDS tree (Session '' when the dataset has no ses- level)."""
    out = []
    for sub in sorted(glob.glob(opj(bids_dir, 'sub-*'))):
        if not os.path.isdir(sub):
            continue
        ID = opb(sub)[4:]
        ses = sorted(d for d in glob.glob(opj(sub, 'ses-*')) if os.path.isdir(d))
        if ses:
            out += [(ID, opb(s)[4:], s) for s in ses]
        else:
            out.append((ID, '', sub))
    return out


def _key(ID, Session):
    return f'sub-{ID}' + (f'_ses-{Session}' if Session != '' else '')


def _qc_check_images(bids_dir):
    """{check folder: [image paths]} of <BIDS>/QC/<check>/ (the report folders written here are excluded)."""
    root = opj(bids_dir, 'QC')
    out = {}
    if not os.path.isdir(root):
        return out
    for d in sorted(os.listdir(root)):
        p = opj(root, d)
        if not os.path.isdir(p) or d in ('subjects', 'dataset'):
            continue
        imgs = sorted(f for f in glob.glob(opj(p, '**', '*'), recursive=True) if f.lower().endswith(IMG_EXT))
        if imgs:
            out[d] = imgs
    return out


def assign_images(images, subjects_sessions):
    """
    Attribute the files of QC/<check>/ to subjects / sessions. The QC writers name their files differently
    ('<ID>_<ses>_...', 'sub-<ID>_ses-<S>_task-...', '<ID>meanIMG...'), so the subject is the LONGEST subject ID that
    starts the name (after an optional 'sub-') or appears as 'sub-<ID>'; the session is the one written in the name
    ('ses-<S>' or '<ID>_<S>_'), else the image goes to every session of the subject.
    Returns {(ID, Session): [paths]}.
    """
    by_id = {}
    for ID, S, _ in subjects_sessions:
        by_id.setdefault(str(ID), []).append(str(S))
    ids = sorted(by_id, key=len, reverse=True)
    out = {}
    for p in images:
        name = opb(p)
        base = name[4:] if name.startswith('sub-') else name
        ID = next((i for i in ids if base.startswith(i) or f'sub-{i}_' in name or f'sub-{i}.' in name), None)
        if ID is None:
            continue
        sess = by_id[ID]
        rest = base[len(ID):] if base.startswith(ID) else name.split(f'sub-{ID}', 1)[1]
        hit = [S for S in sess if S != '' and (f'ses-{S}_' in rest or f'ses-{S}.' in rest
                                                or rest.startswith(f'_{S}_') or rest.startswith(f'_{S}.'))]
        for S in (hit or sess):
            out.setdefault((ID, S), []).append(p)
    return out


# =====================================================================================================================
# surface QC (from the surface variant tests: test_surface_variants_library.py)
# =====================================================================================================================
def _fs_scale(prepro_dir, ID):
    """Real mm per FreeSurfer mm: 1 / s, s = header scale of smallbrain.get (first diagonal term of the .dat file)."""
    info = _load_json(opj(prepro_dir, ID + '_rescale.json')) or {}
    datfile = info.get('datfile', '')
    if not datfile or not opi(datfile):
        return 1.0
    with open(datfile) as fh:
        lines = [l.strip() for l in fh if l.strip()]
    try:
        s = float(lines[4].split()[0])
        return 1.0 / s if s > 0 else 1.0
    except (IndexError, ValueError):
        return 1.0


def _surf_path(FS_dir, sp, h, kind):
    s = opj(FS_dir, sp, 'surf')
    names = (h + '.pial', h + '.pial.T1') if kind == 'pial' else (h + '.' + kind,)
    return next((opj(s, n) for n in names if opi(opj(s, n))), None)


def _thickness_real(FS_dir, sp, k):
    import nibabel as nib
    vals = []
    for h in HEMIS:
        f = opj(FS_dir, sp, 'surf', h + '.thickness')
        if not opi(f):
            continue
        t = nib.freesurfer.read_morph_data(f)
        lab = opj(FS_dir, sp, 'label', h + '.cortex.label')
        if opi(lab):
            t = t[nib.freesurfer.read_label(lab)]
        vals.append(t[t > 0])
    return np.concatenate(vals) * k if vals else np.array([])


def _plane_segments(v, f, axis, c):
    """Segments where the triangle mesh (v in voxel coords) crosses the plane v[:, axis] == c."""
    d = v[:, axis] - c
    s = d[f] > 0
    f = f[(s.sum(1) == 1) | (s.sum(1) == 2)]
    if not len(f):
        return np.zeros((0, 2, 2))
    pts = []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        a_, b_ = f[:, i], f[:, j]
        m = (d[a_] > 0) != (d[b_] > 0)
        t = np.where(m, d[a_] / np.where(m, d[a_] - d[b_], 1), np.nan)
        pts.append(v[a_] + t[:, None] * (v[b_] - v[a_]))
    P = np.stack(pts, 1)
    ok = ~np.isnan(P[:, :, 0])
    segs = [P[r][ok[r]][:2] for r in range(len(P)) if ok[r].sum() == 2]
    if not segs:
        return np.zeros((0, 2, 2))
    keep = [k for k in range(3) if k != axis]
    return np.stack(segs)[:, :, keep]


def aseg_fidelity(FS_dir, sp, k=1.0):
    """
    Distance of the white / pial vertices to the boundaries they should follow, in real mm: white -> boundary of
    filled.mgz > 0; pial -> outer boundary of the aseg brain (aseg > 0 without CSF codes, holes filled).
    Returns {'white'|'pial': (median |d| in real mm, fraction of vertices with |d| <= 1 voxel)}.
    """
    import nibabel as nib
    from scipy import ndimage as ndi
    mri = opj(FS_dir, sp, 'mri')
    out = {}
    if not (opi(opj(mri, 'filled.mgz')) and opi(opj(mri, 'aseg.mgz'))):
        return out
    fi = nib.load(opj(mri, 'filled.mgz'))
    vox = float(np.mean(fi.header.get_zooms()[:3]))
    inv = np.linalg.inv(fi.header.get_vox2ras_tkr())
    aseg = np.asarray(nib.load(opj(mri, 'aseg.mgz')).dataobj).astype(np.int32)
    masks = {'white': np.asarray(fi.dataobj) > 0,
             'pial': ndi.binary_fill_holes((aseg > 0) & ~np.isin(aseg, CSF_CODES))}
    for kind, m in masks.items():
        sd = (ndi.distance_transform_edt(~m) - ndi.distance_transform_edt(m)).astype(np.float32)
        d = []
        for h in HEMIS:
            p = _surf_path(FS_dir, sp, h, kind)
            if p:
                v = nib.freesurfer.read_geometry(p)[0]
                ijk = (inv @ np.c_[v, np.ones(len(v))].T)[:3]
                d.append(np.abs(ndi.map_coordinates(sd, ijk, order=1, mode='nearest')))
        if d:
            d = np.concatenate(d)
            out[kind] = (float(np.median(d)) * vox * k, float(np.mean(d <= 1.0)))
    return out


def surface_qc(FS_dir, sp, out_png, real_mm_per_fs_mm=1.0, title=''):
    """
    Surface QC of one FreeSurfer subject: white (blue) / pial (red) on orig.mgz at coronal 30/50/70 % and axial 55 %
    of the brain extent, thickness histogram in real mm, aseg fidelity. Writes out_png and out_png .json (metrics).
    Returns the metrics dict, {} when there is no surface.
    """
    import nibabel as nib
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    orig_f = next((opj(FS_dir, sp, 'mri', n) for n in ('orig.mgz', 'brain.mgz') if opi(opj(FS_dir, sp, 'mri', n))),
                  None)
    if orig_f is None or _surf_path(FS_dir, sp, 'lh', 'white') is None:
        return {}
    k = float(real_mm_per_fs_mm)
    img = nib.load(orig_f)
    vol = np.asarray(img.dataobj).astype(np.float32)
    inv = np.linalg.inv(img.header.get_vox2ras_tkr())
    nz = np.argwhere(vol > 0)
    lo, hi = nz.min(0), nz.max(0)
    meshes, n_white = {}, 0
    for kind in ('white', 'pial'):
        for h in HEMIS:
            p = _surf_path(FS_dir, sp, h, kind)
            if p:
                vv, ff = nib.freesurfer.read_geometry(p)
                if kind == 'white':
                    n_white += len(vv)
                vv = (inv @ np.c_[vv, np.ones(len(vv))].T).T[:, :3]
                meshes[(kind, h)] = (vv, ff)

    t = _thickness_real(FS_dir, sp, k)
    try:
        fid = aseg_fidelity(FS_dir, sp, k)
    except Exception:
        fid = {}
    metrics = dict(real_mm_per_fs_mm=k, n_vertices_white=int(n_white),
                   has_pial=any(kk == 'pial' for kk, _ in meshes),
                   has_sphere_reg=all(opi(opj(FS_dir, sp, 'surf', h + '.sphere.reg')) for h in HEMIS))
    for name, q in (('thick_median_mm', 50), ('thick_p05_mm', 5), ('thick_p95_mm', 95)):
        metrics[name] = round(float(np.percentile(t, q)), 3) if len(t) else None
    for kind in ('white', 'pial'):
        md, fr = fid.get(kind, (None, None))
        metrics[kind + '_dist_med_mm'] = round(md, 3) if md is not None else None
        metrics[kind + '_within_1vox'] = round(fr, 3) if fr is not None else None

    planes = [('coronal 30 %', 2, 0.30), ('coronal 50 %', 2, 0.50), ('coronal 70 %', 2, 0.70), ('axial 55 %', 1, 0.55)]
    fig, axes = plt.subplots(1, len(planes) + 1, figsize=(3.2 * (len(planes) + 1), 3.6))
    for j, (ptitle, axis, frac) in enumerate(planes):
        ax = axes[j]
        c = int(round(lo[axis] + frac * (hi[axis] - lo[axis])))
        ax.imshow(np.take(vol, c, axis=axis).T, cmap='gray', origin='upper', interpolation='nearest')
        for (kind, h), (vv, ff) in meshes.items():
            seg = _plane_segments(vv, ff, axis, c)
            if len(seg):
                ax.add_collection(LineCollection(seg, colors='tab:blue' if kind == 'white' else 'tab:red',
                                                 linewidths=0.6))
        keep = [kk for kk in range(3) if kk != axis]
        ax.set_xlim(lo[keep[0]] - 3, hi[keep[0]] + 3)
        ax.set_ylim(hi[keep[1]] + 3, lo[keep[1]] - 3)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(ptitle, fontsize=9)
    ax = axes[-1]
    if len(t):
        ax.hist(t, bins=60, color='tab:gray')
        ax.axvline(np.median(t), color='k', lw=1)
        ax.set_title(f'median {np.median(t):.2f} mm (real)', fontsize=9)
    else:
        ax.text(0.5, 0.5, 'no thickness', ha='center', va='center', transform=ax.transAxes)
    ax.set_xlabel('thickness (real mm)', fontsize=8)
    ax.tick_params(labelsize=7)
    fid_txt = '  '.join(f"{kk} ≤1 vox {metrics[kk + '_within_1vox']:.0%}" for kk in ('white', 'pial')
                        if metrics[kk + '_within_1vox'] is not None)
    fig.suptitle(f'{title or sp}: white (blue) / pial (red) on {opb(orig_f)}   {fid_txt}', fontsize=10)
    fig.tight_layout()
    os.makedirs(opd(out_png) or '.', exist_ok=True)
    fig.savefig(out_png, dpi=110)
    plt.close(fig)
    with open(os.path.splitext(out_png)[0] + '.json', 'w') as fh:
        json.dump(metrics, fh, indent=1)
    return metrics


def surface_qc_subject(data_path, ID, Session='', bids_dir=''):
    """
    Surface QC of one subject of an EDNiX dataset: anat/QC_anat/<ID>_surface_QC.png (+ .json), copied to
    <BIDS>/QC/surfaces/. Returns the metrics ({} without surfaces).
    """
    path_anat = opj(data_path, 'anat')
    FS_dir = opj(path_anat, 'freesurfer')
    k = _fs_scale(opj(path_anat, 'preprocessing'), ID)
    out = opj(path_anat, 'QC_anat', f'{ID}_surface_QC.png')
    m = surface_qc(FS_dir, ID, out, k, title=_key(ID, Session))
    if m and bids_dir:
        dst = opj(bids_dir, 'QC', 'surfaces')
        os.makedirs(dst, exist_ok=True)
        shutil.copyfile(out, opj(dst, _key(ID, Session) + '_surface_QC.png'))
    return m


# =====================================================================================================================
# per subject: metrics, status, logs
# =====================================================================================================================
def _metric_rows(data_path, modality, folder):
    rows = []
    for f in sorted(glob.glob(opj(data_path, modality, folder, '*_QC_values.json'))):
        v = _load_json(f)
        if isinstance(v, dict):
            rows.append(dict(run=opb(f)[:-len('_QC_values.json')], **_flatten(v)))
    return rows


def _specificity_rows(data_path):
    f = opj(data_path, 'func', 'QC', 'all_runs_classification.csv')
    if not opi(f):
        return []
    import pandas as pd
    try:
        return pd.read_csv(f).to_dict('records')
    except Exception:
        return []


def _pet_motion_rows(data_path):
    rows = []
    for mod in ('pet', 'ct'):
        for f in sorted(glob.glob(opj(data_path, mod, 'orig', 'preprocessing', 'forMC', '*_motion.tsv'))):
            import pandas as pd
            try:
                df = pd.read_csv(f, sep='\t')
            except (OSError, ValueError) as e:
                print(f'WARNING: unreadable motion file {f}: {e}')
                continue
            row = dict(file=opb(f), frames=len(df))
            for c in df.columns:
                v = pd.to_numeric(df[c], errors='coerce')
                if v.notna().any() and ('displ' in c.lower() or c.lower().startswith(('t', 'r'))):
                    row[c + '_max'] = float(np.nanmax(np.abs(v)))
            rows.append(row)
    return rows


def _log_summary(data_path):
    rows = []
    for mod in ('anat', 'func', 'pet', 'ct'):
        for f in sorted(glob.glob(opj(data_path, mod, '* session *.txt'))):
            if 'launcher_parameters' in opb(f):
                continue
            try:
                lines = open(f, errors='replace').read().splitlines()
            except OSError:
                continue
            errs = [l for l in lines if ERROR_RE.search(l)]
            warns = [l for l in lines if not ERROR_RE.search(l) and WARN_RE.search(l)]
            rows.append(dict(modality=mod, log=opb(f), errors=len(errs), warnings=len(warns),
                             last_error=(errs[-1].strip()[:200] if errs else '')))
    return rows


def status(data_path, ID):
    """Presence of the main outputs (cheap file checks; the launchers' resume logic is in Tools/ednix_checkpoint)."""
    a = opj(data_path, 'anat')
    vol, msk = opj(a, 'native', 'volumes'), opj(a, 'native', 'volumes', 'masks')
    fs = opj(a, 'freesurfer', ID, 'surf')
    st = {}
    if os.path.isdir(a):
        st.update({
            'anat: skull-strip (SS)': bool(glob.glob(opj(vol, ID + '_space-acpc_desc-SS_*.nii.gz'))),
            'anat: registration to template': bool(glob.glob(opj(a, 'matrices', 'native_to_*_Final_0GenericAffine.mat'))),
            'anat: atlases in native space': bool(glob.glob(opj(vol, 'labels', ID + '_seg-*_dseg.nii.gz'))),
            'anat: final mask': opi(opj(msk, ID + '_space-acpc_mask.nii.gz')),
            'anat: manual mask(s)': bool(glob.glob(opj(msk, '*final_mask*')) + glob.glob(opj(msk, '*backtonative_mask*'))),
            'surf: white': all(opi(opj(fs, h + '.white')) for h in HEMIS),
            'surf: pial': all(opi(opj(fs, h + '.pial')) for h in HEMIS),
            'surf: sphere.reg (optional)': all(opi(opj(fs, h + '.sphere.reg')) for h in HEMIS),
            'surf: Workbench spec': opi(opj(a, 'native', 'surfaces', 'Native_resol', ID + '_native_LR.wb.spec')),
            'anat: QC (step 16)': bool(glob.glob(opj(a, 'QC_anat', '*_QC_summary.png'))),
        })
    fdir = opj(data_path, 'func')
    if os.path.isdir(fdir):
        raw = [f for f in glob.glob(opj(fdir, '*_bold.nii*')) if opb(f).startswith('sub-')]
        qcs = glob.glob(opj(fdir, 'QC', '*_QC_values.json'))
        st.update({
            'func: raw runs': len(raw),
            'func: runs with QC values (step 12)': len(qcs),
            'func: FC specificity (step 13)': opi(opj(fdir, 'QC', 'all_runs_classification.csv')),
        })
    for mod in ('pet', 'ct'):
        d = opj(data_path, mod)
        if os.path.isdir(d):
            st[f'{mod}: registration record'] = bool(glob.glob(opj(d, 'orig', 'matrices', '*.json')))
            st[f'{mod}: QC images'] = bool(glob.glob(opj(d, 'QC', '*.png')))
    return st


def _subject_images(data_path):
    out = {}
    for label, sub in (('anat/QC_anat', opj('anat', 'QC_anat')), ('func/QC', opj('func', 'QC')),
                       ('pet/QC', opj('pet', 'QC')), ('ct/QC', opj('ct', 'QC'))):
        imgs = sorted(f for f in glob.glob(opj(data_path, sub, '*')) if f.lower().endswith(IMG_EXT))
        if imgs:
            out[label] = imgs
    return out


def _adaptive_registration(data_path):
    """Decisions of the SyNadaptive registration (modalities/anat/registration_adaptive.py), if it was used."""
    rows = []
    for f in sorted(glob.glob(opj(data_path, 'anat', '**', '*adaptive_QC.json'), recursive=True)):
        d = _load_json(f) or {}
        chosen = d.get(d.get('chosen', ''), {})
        rows.append(dict(file=opb(f), chosen=d.get('chosen'), tissue_dice=chosen.get('tissue_dice'),
                         mask_dice=chosen.get('mask_dice'), folding=chosen.get('folding'),
                         decision=d.get('decision')))
    return rows


def collect_subject(bids_dir, ID, Session, data_path, surface=True):
    """Every QC element of one subject / session (dict), recomputing the surface QC if it is missing."""
    surf_json = opj(data_path, 'anat', 'QC_anat', f'{ID}_surface_QC.json')
    surf = _load_json(surf_json) or {}
    if surface and not surf and os.path.isdir(opj(data_path, 'anat', 'freesurfer', ID, 'surf')):
        try:
            surf = surface_qc_subject(data_path, ID, Session, bids_dir)
        except Exception as e:
            surf = {'error': f'surface QC failed: {e}'}
    return dict(ID=ID, Session=Session, data_path=data_path, key=_key(ID, Session),
                status=status(data_path, ID), logs=_log_summary(data_path),
                anat=_metric_rows(data_path, 'anat', 'QC_anat'), func=_metric_rows(data_path, 'func', 'QC'),
                specificity=_specificity_rows(data_path), pet_motion=_pet_motion_rows(data_path),
                surface=surf, images=_subject_images(data_path),
                registration=_adaptive_registration(data_path))


def subject_report(bids_dir, info, check_images):
    """Write QC/subjects/<key>.html for one subject (info from collect_subject). Returns the path."""
    out_dir = opj(bids_dir, 'QC', 'subjects')
    os.makedirs(out_dir, exist_ok=True)
    out = opj(out_dir, info['key'] + '.html')
    b = ['<p><a href="../QC_all_subjects.html">← all subjects</a> · <a href="../QC_dataset.html">dataset</a></p>']
    n_err = sum(r['errors'] for r in info['logs'])
    n_warn = sum(r['warnings'] for r in info['logs'])
    b.append(f'<p><span class="badge {"err" if n_err else "fine"}">{n_err} error line(s) in the logs</span>'
             f'<span class="badge {"warn" if n_warn else "fine"}">{n_warn} warning line(s)</span></p>')
    b.append('<h2>Processing status</h2>')
    b.append(_table([dict(output=k, done=v) for k, v in info['status'].items()], ['output', 'done']))
    if info['logs']:
        b.append('<h3>Logs</h3>' + _table(info['logs'], ['modality', 'log', 'errors', 'warnings', 'last_error']))
    if info['surface']:
        b.append('<h2>Surfaces</h2>')
        b.append(_table([dict(metric=k, value=v) for k, v in info['surface'].items()], ['metric', 'value']))
        png = opj(info['data_path'], 'anat', 'QC_anat', f"{info['ID']}_surface_QC.png")
        if opi(png):
            b.append(_images_html([png], out_dir, big=True))
    for title, rows in (('Anatomical metrics (step 16)', info['anat']), ('Functional metrics per run (step 12)',
                                                                           info['func']),
                        ('FC specificity per run (step 13)', info['specificity']),
                        ('PET / CT motion', info['pet_motion']),
                        ('Adaptive registration (SyNadaptive)', info.get('registration', []))):
        if rows:
            cols = list(dict.fromkeys(c for r in rows for c in r))
            if len(rows) == 1 and len(cols) > 12:
                b.append(f'<h2>{_esc(title)}</h2>' + _table([dict(metric=c, value=rows[0][c]) for c in cols],
                                                           ['metric', 'value']))
            else:
                b.append(f'<h2>{_esc(title)}</h2>' + _table(rows, cols))
    b.append('<h2>Images</h2>')
    for label, imgs in info['images'].items():
        b.append(f'<h3>{_esc(label)}</h3>' + _images_html(imgs, out_dir))
    for check, imgs in sorted(check_images.items()):
        b.append(f'<h3>QC/{_esc(check)}</h3>' + _images_html(imgs, out_dir))
    with open(out, 'w') as fh:
        fh.write(_page(f'EDNiX QC — {info["key"]}', ''.join(b), opb(os.path.normpath(bids_dir))))
    return out


# =====================================================================================================================
# all subjects
# =====================================================================================================================
def _summary_row(info):
    r = dict(subject=info['ID'], session=info['Session'])
    r['log_errors'] = sum(x['errors'] for x in info['logs'])
    for k, v in info['status'].items():
        if k.startswith(('anat: SS', 'anat: skull', 'anat: registration', 'surf: pial', 'func: runs with')):
            r[k] = v
    for m in KEY_METRICS['surface']:
        r['surf.' + m] = info['surface'].get(m) if isinstance(info['surface'], dict) else None
    runs = info['func']
    for m in KEY_METRICS['func']:
        vals = [x.get(m) for x in runs if isinstance(x.get(m), (int, float))]
        r['func.' + m] = float(np.mean(vals)) if vals else None
    reg = info.get('registration') or []
    if reg:
        r['reg.chosen'] = ' / '.join(str(x['chosen']) for x in reg)
        r['reg.tissue_dice'] = reg[-1]['tissue_dice']
    if info['specificity']:
        cats = [str(x.get('category', x.get('classification', ''))) for x in info['specificity']]
        r['func.specific_runs'] = f"{sum('specific' == c.lower() for c in cats)}/{len(cats)}"
    if info['anat']:
        for c, v in info['anat'][0].items():
            if isinstance(v, (int, float)) and any(s in c.lower() for s in ('cnr', 'cjv', 'snr_gm', 'snr_wm',
                                                                          'snr_total', 'fwhm')):
                r['anat.' + c] = v
    return r


def group_report(bids_dir, infos, check_images):
    out = opj(bids_dir, 'QC', 'QC_all_subjects.html')
    rows = [_summary_row(i) for i in infos]
    cols = list(dict.fromkeys(c for r in rows for c in r))
    outliers = set()
    for c in cols:
        if c in ('subject', 'session') or c.endswith('_runs'):
            continue
        z = _robust_z([r.get(c) for r in rows])
        outliers |= {(i, c) for i in range(len(rows)) if np.isfinite(z[i]) and abs(z[i]) > Z_THRESH}
    for i, r in enumerate(rows):
        if r.get('log_errors'):
            outliers.add((i, 'log_errors'))
    links = {i: f"subjects/{infos[i]['key']}.html" for i in range(len(infos))}
    b = ['<p><a href="QC_dataset.html">dataset summary</a> · <a href="EDNiX_BIDS_summary.xlsx">Excel</a></p>',
         f'<p class="muted">{len(rows)} subject/session(s). Red: |robust z| > {Z_THRESH} (median / MAD over the '
         f'dataset) or errors in the logs. Click a subject for its recap.</p>',
         _table(rows, cols, outliers, links)]
    for check, imgs in sorted(check_images.items()):
        b.append(f'<h2>QC/{_esc(check)} ({len(imgs)})</h2>' + _images_html(imgs, opd(out)))
    for label in ('anat/QC_anat', 'func/QC', 'pet/QC', 'ct/QC'):
        imgs = [p for i in infos for p in i['images'].get(label, [])]
        if imgs:
            b.append(f'<h2>{label} ({len(imgs)})</h2>' + _images_html(imgs, opd(out)))
    with open(out, 'w') as fh:
        fh.write(_page('EDNiX QC — all subjects', ''.join(b), opb(os.path.normpath(bids_dir))))
    return out, rows


# =====================================================================================================================
# dataset (BIDS) level: acquisition parameters from the raw files
# =====================================================================================================================
TR_RANGE = (0.1, 10.0)        # s; outside = unit error in the sidecar (e.g. ms written in a field in s)
TR_TOLERANCE = 0.05


def _sidecar(img, bids_dir):
    """BIDS inheritance: <bids>/<suffix>.json, <bids>/task-<x>_<suffix>.json, then the file's own sidecar (wins)."""
    name = re.sub(r'\.nii(\.gz)?$', '', opb(img))
    suffix = name.split('_')[-1]
    meta = {}
    task = re.search(r'task-([A-Za-z0-9]+)', name)
    for f in [opj(bids_dir, suffix + '.json')] + ([opj(bids_dir, f'task-{task.group(1)}_{suffix}.json')] if task
                                                   else []) + [opj(opd(img), name + '.json')]:
        d = _load_json(f) if opi(f) else None
        if isinstance(d, dict):
            meta.update(d)
    return meta


def _header(img):
    import nibabel as nib
    try:
        h = nib.load(img).header
    except Exception:
        return {}
    shape = tuple(int(x) for x in h.get_data_shape())
    zooms = tuple(float(x) for x in h.get_zooms())
    out = dict(matrix='x'.join(map(str, shape[:3])), voxel_mm='x'.join(f'{z:.3g}' for z in zooms[:3]))
    if len(shape) > 3:
        out['volumes'] = shape[3]
        tr = zooms[3] if len(zooms) > 3 else None
        try:
            unit = h.get_xyzt_units()[1]
        except Exception:
            unit = 'sec'
        if tr is not None:
            tr = tr / 1000.0 if unit == 'msec' else (tr / 1e6 if unit == 'usec' else tr)
        out['TR_header_s'] = tr
    return out


def _entity(name, key):
    m = re.search(rf'(?:^|_){key}-([A-Za-z0-9]+)', name)
    return m.group(1) if m else ''


def acquisition_rows(bids_dir):
    """One row per raw image (anat, func, fmap, pet, ct) of the BIDS tree (EDNiX outputs are not read)."""
    rows = []
    for ID, S, dp in sessions(bids_dir):
        for mod in ('anat', 'func', 'fmap', 'pet', 'ct'):
            for img in sorted(glob.glob(opj(dp, mod, 'sub-*.nii*'))):
                name = re.sub(r'\.nii(\.gz)?$', '', opb(img))
                meta = _sidecar(img, bids_dir)
                r = dict(subject=ID, session=S, modality=mod, file=opb(img), suffix=name.split('_')[-1],
                         task=_entity(name, 'task'), run=_entity(name, 'run'), acq=_entity(name, 'acq'))
                r.update(_header(img))
                for k in ('MagneticFieldStrength', 'Manufacturer', 'ManufacturersModelName', 'RepetitionTime',
                          'EchoTime', 'FlipAngle', 'PhaseEncodingDirection', 'TotalReadoutTime', 'EffectiveEchoSpacing',
                          'TracerName', 'InjectedRadioactivity', 'InjectedRadioactivityUnits', 'FrameDuration',
                          'ReconMethodName'):
                    if k in meta:
                        v = meta[k]
                        r[k] = (len(v) if k == 'FrameDuration' and isinstance(v, list) else
                                v if isinstance(v, (int, float, str)) else str(v))
                r['SliceTiming'] = 'SliceTiming' in meta
                if mod == 'func':
                    tr_s, tr_h = r.get('RepetitionTime'), r.get('TR_header_s')
                    flag = []
                    if not isinstance(tr_s, (int, float)):
                        flag.append('no sidecar TR')
                    elif not TR_RANGE[0] <= tr_s <= TR_RANGE[1]:
                        flag.append('sidecar TR out of range')
                    if isinstance(tr_s, (int, float)) and tr_h and abs(tr_s - tr_h) > TR_TOLERANCE * tr_h:
                        flag.append('sidecar vs header differ')
                    r['TR_check'] = '; '.join(flag) or 'ok'
                    tr = tr_s if isinstance(tr_s, (int, float)) and TR_RANGE[0] <= tr_s <= TR_RANGE[1] else tr_h
                    r['TR_used_s'] = tr
                    if tr and r.get('volumes'):
                        r['duration_min'] = round(tr * r['volumes'] / 60.0, 2)
                rows.append(r)
    return rows


def _participants(bids_dir):
    f = opj(bids_dir, 'participants.tsv')
    if not opi(f):
        return None
    import pandas as pd
    try:
        return pd.read_csv(f, sep='\t')
    except Exception:
        return None


def _hist_png(values, title, xlabel, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    v = np.asarray([x for x in values if isinstance(x, (int, float)) and np.isfinite(x)], dtype=float)
    if len(v) < 2:
        return None
    fig, ax = plt.subplots(figsize=(4, 2.6))
    ax.hist(v, bins=min(30, max(5, len(v) // 2)), color='0.45')
    ax.axvline(np.median(v), color='k', lw=1)
    ax.set_title(f'{title} (n={len(v)}, median {np.median(v):.3g})', fontsize=9)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.tick_params(labelsize=7)
    fig.tight_layout()
    os.makedirs(opd(out), exist_ok=True)
    fig.savefig(out, dpi=110)
    plt.close(fig)
    return out


def dataset_report(bids_dir, infos, summary_rows):
    """QC/EDNiX_BIDS_summary.xlsx and QC/QC_dataset.html. Returns (xlsx, html)."""
    import pandas as pd
    qc = opj(bids_dir, 'QC')
    acq = pd.DataFrame(acquisition_rows(bids_dir))
    part = _participants(bids_dir)
    ss = sessions(bids_dir)

    # subjects sheet: one row per subject / session with the number of raw images per modality
    subj = pd.DataFrame([dict(subject=i, session=s) for i, s, _ in ss])
    if len(acq):
        counts = acq.groupby(['subject', 'session', 'modality']).size().unstack(fill_value=0).reset_index()
        subj = subj.merge(counts, on=['subject', 'session'], how='left').fillna(0)
    if part is not None and 'participant_id' in part.columns:
        p = part.copy()
        p['subject'] = p['participant_id'].astype(str).str.replace('^sub-', '', regex=True)
        subj = subj.merge(p.drop(columns=['participant_id']), on='subject', how='left')

    # overview: per modality and suffix, n subjects / sessions / images and parameter ranges
    over = []
    if len(acq):
        for (mod, suf), g in acq.groupby(['modality', 'suffix']):
            r = dict(modality=mod, suffix=suf, n_subjects=g['subject'].nunique(),
                     n_sessions=g[['subject', 'session']].drop_duplicates().shape[0], n_images=len(g))
            for c in ('MagneticFieldStrength', 'Manufacturer', 'voxel_mm', 'matrix'):
                if c in g:
                    r[c] = ' / '.join(map(str, g[c].dropna().unique()[:4]))
            for c in ('TR_used_s', 'EchoTime', 'volumes', 'duration_min'):
                if c in g and g[c].notna().any():
                    v = pd.to_numeric(g[c], errors='coerce').dropna()
                    if len(v):
                        r[c] = f'{v.median():.4g} [{v.min():.4g}-{v.max():.4g}]' if v.min() != v.max() else f'{v.iloc[0]:.4g}'
            if 'TR_check' in g:
                r['TR_flagged'] = int((g['TR_check'] != 'ok').sum())
            over.append(r)
    over = pd.DataFrame(over)

    status_df = pd.DataFrame([dict(subject=i['ID'], session=i['Session'], **i['status']) for i in infos])
    func_df = pd.DataFrame([dict(subject=i['ID'], session=i['Session'], **r) for i in infos for r in i['func']])
    anat_df = pd.DataFrame([dict(subject=i['ID'], session=i['Session'], **r) for i in infos for r in i['anat']])
    surf_df = pd.DataFrame([dict(subject=i['ID'], session=i['Session'], **i['surface']) for i in infos
                            if isinstance(i['surface'], dict) and i['surface']])
    spec_df = pd.DataFrame([dict(subject=i['ID'], session=i['Session'], **r) for i in infos for r in i['specificity']])
    logs_df = pd.DataFrame([dict(subject=i['ID'], session=i['Session'], **r) for i in infos for r in i['logs']])

    xlsx = opj(qc, 'EDNiX_BIDS_summary.xlsx')
    os.makedirs(qc, exist_ok=True)
    with pd.ExcelWriter(xlsx) as w:
        for name, df in (('overview', over), ('subjects', subj), ('acquisition', acq), ('status', status_df),
                         ('QC_summary', pd.DataFrame(summary_rows)), ('anat_QC', anat_df), ('surface_QC', surf_df),
                         ('func_QC', func_df), ('FC_specificity', spec_df), ('logs', logs_df)):
            (df if len(df) else pd.DataFrame({'info': ['no data']})).to_excel(w, sheet_name=name, index=False)

    # html
    figs = opj(qc, 'dataset')
    plots = []
    if len(acq) and 'TR_used_s' in acq:
        plots.append(_hist_png(acq['TR_used_s'].tolist(), 'TR', 's', opj(figs, 'TR.png')))
    if len(func_df):
        for c, xl in (('mean_fd', 'mm'), ('censor_fraction', 'fraction'), ('TSNR', '')):
            if c in func_df:
                plots.append(_hist_png(pd.to_numeric(func_df[c], errors='coerce').tolist(), c, xl,
                                       opj(figs, c + '.png')))
    if len(surf_df) and 'thick_median_mm' in surf_df:
        plots.append(_hist_png(pd.to_numeric(surf_df['thick_median_mm'], errors='coerce').tolist(),
                               'median cortical thickness', 'real mm', opj(figs, 'thickness.png')))
    plots = [p for p in plots if p]
    html_out = opj(qc, 'QC_dataset.html')
    b = ['<p><a href="QC_all_subjects.html">all subjects</a> · <a href="EDNiX_BIDS_summary.xlsx">Excel</a></p>',
         f'<p>{subj["subject"].nunique() if len(subj) else 0} subject(s), {len(ss)} subject/session(s), '
         f'{len(acq)} raw image(s).</p>']
    if len(over):
        b.append('<h2>Acquisitions</h2>' + _table(over.to_dict('records'), list(over.columns)))
    if len(acq) and 'TR_check' in acq:
        bad = acq[(acq['modality'] == 'func') & (acq['TR_check'] != 'ok')]
        if len(bad):
            b.append(f'<h3>TR flagged ({len(bad)})</h3>' +
                     _table(bad[['subject', 'session', 'file', 'RepetitionTime', 'TR_header_s', 'TR_check']]
                            .to_dict('records'), ['subject', 'session', 'file', 'RepetitionTime', 'TR_header_s',
                                                  'TR_check']))
    if len(status_df):
        done = status_df.drop(columns=['subject', 'session'])
        rows = []
        for c in done.columns:
            v = done[c]
            n = int((v.astype(float) > 0).sum()) if v.dtype != object else int(v.astype(bool).sum())
            rows.append(dict(output=c, sessions_done=f'{n}/{len(v)}'))
        b.append('<h2>Processing status</h2>' + _table(rows, ['output', 'sessions_done']))
    if plots:
        b.append('<h2>Distributions</h2>' + _images_html(plots, qc))
    with open(html_out, 'w') as fh:
        fh.write(_page('EDNiX QC — dataset', ''.join(b), opb(os.path.normpath(bids_dir))))
    return xlsx, html_out


# =====================================================================================================================
# entry points
# =====================================================================================================================
def build(bids_dir, subjects=None, surface=True, verbose=True):
    """
    Build the three QC levels of one BIDS dataset. subjects: optional list of (ID, Session) to restrict the
    recaps (the dataset pages always cover every subject that has outputs). Returns a dict of the written paths.
    """
    bids_dir = os.path.normpath(bids_dir)
    ss = sessions(bids_dir)
    check_images = _qc_check_images(bids_dir)
    flat = [p for v in check_images.values() for p in v]
    per_subject = assign_images(flat, ss)
    infos, written = [], []
    for ID, S, dp in ss:
        try:
            info = collect_subject(bids_dir, ID, S, dp, surface=surface)
        except Exception:
            if verbose:
                print(f'QC_report: {_key(ID, S)} failed:\n{traceback.format_exc()}')
            continue
        infos.append(info)
        if subjects is None or (ID, S) in subjects or (str(ID), str(S)) in subjects:
            checks = {}
            for p in per_subject.get((ID, S), []):
                checks.setdefault(os.path.relpath(opd(p), opj(bids_dir, 'QC')), []).append(p)
            written.append(subject_report(bids_dir, info, checks))
    group, rows = group_report(bids_dir, infos, check_images)
    xlsx, ds = dataset_report(bids_dir, infos, rows)
    if verbose:
        print(f'QC_report: {len(written)} subject recap(s) in {opj(bids_dir, "QC", "subjects")}\n'
              f'           all subjects: {group}\n           dataset: {ds}\n           Excel: {xlsx}')
    return dict(subjects=written, group=group, dataset=ds, xlsx=xlsx)


def study_summary(bids_dirs, out_xlsx):
    """One Excel file across datasets: the overview and acquisition sheets of every BIDS folder, stacked."""
    import pandas as pd
    acq_all, over_all = [], []
    for b in bids_dirs:
        name = opb(os.path.normpath(b))
        acq = pd.DataFrame(acquisition_rows(b))
        if not len(acq):
            continue
        acq.insert(0, 'dataset', name)
        acq_all.append(acq)
        for (mod, suf), g in acq.groupby(['modality', 'suffix']):
            r = dict(dataset=name, modality=mod, suffix=suf, n_subjects=g['subject'].nunique(),
                     n_sessions=g[['subject', 'session']].drop_duplicates().shape[0], n_images=len(g))
            for c in ('TR_used_s', 'volumes', 'duration_min'):
                if c in g:
                    v = pd.to_numeric(g[c], errors='coerce').dropna()
                    if len(v):
                        r[c + '_median'] = float(v.median())
                        r[c + '_min'] = float(v.min())
                        r[c + '_max'] = float(v.max())
            for c in ('voxel_mm', 'MagneticFieldStrength'):
                if c in g:
                    r[c] = ' / '.join(map(str, g[c].dropna().unique()[:4]))
            if 'TR_check' in g:
                r['TR_flagged'] = int((g['TR_check'] != 'ok').sum())
            over_all.append(r)
    os.makedirs(opd(out_xlsx) or '.', exist_ok=True)
    with pd.ExcelWriter(out_xlsx) as w:
        pd.DataFrame(over_all or [{'info': 'no data'}]).to_excel(w, sheet_name='overview', index=False)
        (pd.concat(acq_all, ignore_index=True) if acq_all else pd.DataFrame({'info': ['no data']})).to_excel(
            w, sheet_name='acquisition', index=False)
    return out_xlsx


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description='EDNiX QC reports (subject recaps, all subjects, dataset + Excel)')
    ap.add_argument('bids_dir', nargs='+')
    ap.add_argument('--no-surface', action='store_true', help='do not compute missing surface QC figures')
    a = ap.parse_args()
    for bd in a.bids_dir:
        build(bd, surface=not a.no_surface)
