#!/usr/bin/env python3
"""
Test of the EDNiX surface reconstruction (steps 10-12) on the Atlases_library templates, with several ways of
preparing the volumes given to FreeSurfer. Stand-alone: it only READS the library and EDNiX code, and writes
everything under --out (default /scratch2/EDNiX/surface_tests/library). The pipeline is not modified.

Variants (same aseg, same FreeSurfer commands as EDNiX: FS_surf.Wcreate then FS_surf.Pcreate):
  V0  current EDNiX: preFS.prepa_img + smallbrain header rescale + down-sampling of factor/10 (30 for the
      library templates, as in atlases/template.py, i.e. x3 per axis; 20 for subjects in _loop3.py, x2).
  V1  V0 without the extra down-sampling: crop to the brain (+8 voxels) and header rescale only; the grid is
      only coarsened if a dimension still exceeds 256 voxels (Mouse 0.04 mm, Mouselemur 0.09 mm templates).
  V2  V1 + other synthetic brain.mgz / wm.mgz / filled.mgz:
        - WM (and fill) fixed at 110, the value FreeSurfer expects after mri_normalize (as in Precon_all,
          Benn et al., bioRxiv 2025.04.16.649072);
        - GM and CSF keep their within-class intensity variance (1st-99th percentile) but in bands that no
          longer overlap the WM: CSF 10-35, GM 55-85 (V0: CSF 1-35, GM 40-90, WM 80-110 -> GM/WM overlap);
        - hippocampus and amygdala (FS 17/18/53/54) are filled with the subcortex (wm.mgz 250, filled =
          hemisphere value, brain = 110) instead of being cortical GM, so the white surface passes over them
          and the pial surface does not wrap them;
        - brain.mgz smoothed with sigma = 0.5 voxel (sub-voxel boundary instead of a voxel staircase).
  V0m V0 + pial search / thickness limit raised (mris_make_surfaces -max): the default 5 mm is applied in
      rescaled FreeSurfer mm, i.e. 0.87 / 1.41 / 1.26 real mm for the Mouse / Rat / Mouselemur templates,
      below their real isocortical thickness. -max = 3 real mm (>= 5).
  V2m V2 intensities on the V0 grid (EDNiX down-sampling kept) + the same -max.
  V0n V0 intensities on the NATIVE template grid, with the header rescaled so that 1 native voxel = 1 FS mm
      (instead of the smallbrain volume rule): no down-sampling, and FreeSurfer stays in its standard 1 mm
      mode (V1 failed because voxels < 1 FS mm switch mris_make_surfaces to its high-resolution mode). The
      grid is only coarsened if the cropped brain still exceeds 256 voxels (Mouse, Mouselemur). Same idea as
      the HCP-NHP ScaleFactor (Autio et al. 2020). -max = 3 real mm (5 for species with voxels >= 0.3 mm).
  V2n V2 intensities on the V0n grid, same -max.
  V3m V2m with CONSTANT class intensities (CSF 30, GM 75, WM/subcortex 110, smoothed sigma = 1 voxel):
      no within-class texture left in brain.mgz, so the only edges FreeSurfer can follow are the aseg
      boundaries. Tests whether the false sulci/gyri of V2/V2n come from intracortical intensity texture.
  V3n V3m intensities on the V0n grid (native resolution); run for every species.
  summary.csv also gives the fidelity of each surface to the aseg it was built from: median distance (real
  mm) of white vertices to the filled.mgz boundary and of pial vertices to the outer GM boundary of aseg.mgz,
  and the fraction of vertices within 1 voxel of it (medial-wall vertices lower it equally for all variants).
  Sphere registration / aparc (FS_surf.Wconfig) are skipped: only smoothwm/inflated/curvature, which the
  pial step uses, are computed.

Outputs (per species / variant):  <out>/<species>/<variant>/freesurfer/<species>/{mri,surf,label}
  <out>/summary.csv                     status, runtime, vertices, cortical thickness in REAL mm
  <out>/QC_<species>.png                white (blue) / pial (red) on orig.mgz, rows = variants, same relative
                                        slices, + thickness histograms in real mm

Run:  python3 test_surface_variants_library.py                                  (all species found)
      python3 test_surface_variants_library.py --species Mouse Marmoset --variants V0 V2 --jobs 4
      python3 test_surface_variants_library.py --qc-only                      (figures + table from outputs)
"""
import argparse
import csv
import json
import math
import os
import shutil
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

MAIN_PATH_DEFAULT = '/home/cgarin/PycharmProjects/EDNiX/'   # trailing '/' as in the launchers:
# EDNiX builds dirname(MAIN_PATH)/Atlases_library and /Singularity_library, i.e. folders INSIDE EDNiX/
OUT_DEFAULT = '/scratch2/EDNiX/surface_tests/library'
SPECIES_DEFAULT = ['Mouse', 'Rat', 'Mouselemur', 'Marmoset', 'Macaque', 'Cat', 'Dog', 'Pig', 'Chimpanzee', 'Human']
SMALLBRAIN_LIBRARY = ['Mouse', 'Rat', 'Mouselemur', 'Bat']      # as in recontruct_surface_based_on_aseg.py
M_SPECIES = SMALLBRAIN_LIBRARY + ['Marmoset']   # V*m: smallbrain rule forced (Marmoset template < 10 cm3)
FACTOR_LIBRARY = 30                                             # as in atlases/template.py
MAX_DIM = 256
VARIANTS = ['V0', 'V1', 'V2', 'V0m', 'V2m', 'V3m', 'V0n', 'V2n', 'V3n']
NATIVE_GRID = ('V0n', 'V2n', 'V3n')             # 1 native voxel = 1 FS mm
LEGACY_GRID = ('V0', 'V0m', 'V2m', 'V3m')        # variants that keep the EDNiX smallbrain down-sampling
MAX_REAL_THICKNESS = 3.0                 # mm, real units: above any small-species isocortex
HEMIS = ['lh', 'rh']
CSF_CODES = [4, 5, 43, 44, 14, 15, 24, 72, 30, 62, 31, 63]     # preFS.CSF_edit

opj, opb, opd, ope, opi = os.path.join, os.path.basename, os.path.dirname, os.path.exists, os.path.isfile


def ednix_imports(main_path):
    if main_path.rstrip('/') not in sys.path:
        sys.path.insert(0, main_path.rstrip('/'))
    from Tools import run_cmd, Load_EDNiX_requirement, get_orientation
    from atlases import templatefeat
    from modalities.anat.freesurfer import preFS, smallbrain, FS_surf
    return dict(run_cmd=run_cmd, req=Load_EDNiX_requirement, get_orientation=get_orientation,
                templatefeat=templatefeat, preFS=preFS, smallbrain=smallbrain, FS_surf=FS_surf)


# ------------------------------------------------------------------------------------------------------------
def library_paths(E, species, lib):
    fs_tools = opj(lib, 'freesurfer_refs')
    out = E['templatefeat'].get(species, lib, fs_tools, opj(lib, 'BALSA'), 'EDNiX', '', '', 'T1w', 'anat',
                                [['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]], doBALSAspace=True)
    BASE_SS, Aseg_ref = out[9], out[18]
    return fs_tools, BASE_SS, Aseg_ref


def prepa_v2(E, Subname, Ref_file, aseg_file, prep_dir, const=False):
    """V2 (const=False) or V3 (const=True) versions of desc-norm_T1w (brain), seg-wm (wm) and seg-filled
    (filled), written to prep_dir.
    V2: CSF 10-35 and GM 55-85 keep the within-class intensity variation of the template, WM 110, sigma 0.5 vox.
    V3: every class at one value (CSF 30, GM 75, WM + subcortex + hippocampus/amygdala 110), sigma 1 vox:
        the image edges are the aseg boundaries only."""
    import ants
    P = E['preFS']
    hippo_amyg = [17, 18, 53, 54]
    BRAIN = ants.image_read(Ref_file)
    aseg = ants.image_read(aseg_file)
    b = BRAIN.numpy().astype(np.float32)
    a = np.round(aseg.numpy()).astype(np.int32)

    def band(codes, lo, hi):
        m = np.isin(a, codes)
        o = np.zeros(a.shape, np.float32)
        if const or m.sum() < 10:
            o[m] = (lo + hi) / 2.
            return o
        v = b[m].astype(np.float64)
        p1, p99 = np.percentile(v, [1, 99])
        o[m] = lo + np.clip((v - p1) / max(p99 - p1, 1e-6), 0, 1) * (hi - lo)
        return o

    gm_codes = [c for c in P.Cereb_GM_edit if c not in hippo_amyg]
    if const:
        norm = band(P.CSF_edit, 30., 30.) + band(gm_codes, 75., 75.)
        norm[np.isin(a, P.WM_edit + P.aseg_edit + P.wm_aseg_edit + hippo_amyg)] = 110.
        sigma = 1.0
    else:
        norm = band(P.CSF_edit, 10., 35.) + band(gm_codes, 55., 85.)
        norm[np.isin(a, P.WM_edit + hippo_amyg)] = 110.
        sigma = 0.5
    norm_img = ants.smooth_image(BRAIN.new_image_like(norm), sigma, sigma_in_physical_coordinates=False)
    ants.image_write(norm_img, opj(prep_dir, Subname + '_desc-norm_T1w.nii.gz'), ri=False)

    wm = np.zeros(a.shape, np.float32)
    wm[np.isin(a, P.wm_aseg_edit)] = 110.
    wm[np.isin(a, P.aseg_edit + hippo_amyg)] = 250.
    ants.image_write(BRAIN.new_image_like(wm), opj(prep_dir, Subname + '_seg-wm_dseg.nii.gz'), ri=False)

    # filled: same rule as preFS.prepa_img (thalamus centroid decides which side is 127/255), + hippo/amygdala
    hd = ants.image_header_info(Ref_file)
    midsection = math.ceil(int(hd['dimensions'][0]) / 2)
    cent = ants.get_centroids(ants.mask_image(aseg, aseg, 10))
    pos = cent[0][0] * hd['direction'][0][0]
    orig = hd['origin'][0] * hd['direction'][0][0]
    dist = np.abs(pos - orig) / hd['spacing'][0]
    filled = np.zeros(a.shape, np.float32)
    filled[np.isin(a, P.filled_edit[0] + [17, 18])] = 127
    filled[np.isin(a, P.filled_edit[1] + [53, 54])] = 255
    bil = np.isin(a, P.filled_edit[2])
    left, right = bil.copy(), bil.copy()
    left[midsection:, :, :] = False
    right[:midsection, :, :] = False
    filled[left] = 127 if dist < midsection else 255
    filled[right] = 255 if dist < midsection else 127
    ants.image_write(BRAIN.new_image_like(filled), opj(prep_dir, Subname + '_seg-filled_dseg.nii.gz'), ri=False)


def to_freesurfer(E, species, files, FS_dir, change_hd, scaling, legacy_target, variant, diary, sing_fs):
    """Writes orig/brain/wm/aseg/filled.mgz in <FS_dir>/<species>/mri.
    V0: EDNiX behaviour (header rescale, then resampling to legacy_target with smoothing of brain).
    V1/V2: crop to the aseg bounding box (+8 voxels), header rescale, resampling only if a dimension
    would still exceed MAX_DIM voxels."""
    import ants
    names = ['orig.mgz', 'brain.mgz', 'wm.mgz', 'aseg.mgz', 'filled.mgz']
    interps = [3, 3, 1, 1, 1]
    mri = opj(FS_dir, species, 'mri')
    prep = opj(opd(FS_dir), 'resampled')
    os.makedirs(prep, exist_ok=True)
    imgs = [ants.image_read(f) for f in files]
    real_mm_per_fs = float(imgs[0].spacing[0]) / float(scaling[0]) if change_hd else 1.0
    if variant not in LEGACY_GRID:
        nz = np.argwhere(imgs[3].numpy() > 0)
        lo = np.maximum(nz.min(0) - 8, 0)
        hi = np.minimum(nz.max(0) + 9, np.asarray(imgs[3].shape))
        imgs = [ants.crop_indices(im, lo.tolist(), hi.tolist()) for im in imgs]
    ref = imgs[0]
    X, S = ref.origin, ref.spacing
    target = legacy_target
    if variant in NATIVE_GRID:
        # real voxel size kept (coarsened only to fit MAX_DIM), then written with a 1 mm header
        real = max(float(S[0]), math.ceil(max(ref.shape) * float(S[0]) / MAX_DIM * 1000) / 1000)
        real_mm_per_fs = real
        change_hd, scaling = 1, [1., 1., 1.]
        target = None
        if real > float(S[0]):
            imgs = [ants.resample_image(im, (real,) * 3, False, it) for im, it in zip(imgs, interps)]
            for im in imgs:
                im[im < 0] = 0
        ref = imgs[0]
        X, S = ref.origin, ref.spacing
    elif variant not in LEGACY_GRID:
        sp_now = float(scaling[0]) if change_hd else float(S[0])
        target = None if max(ref.shape) <= MAX_DIM else math.ceil(max(ref.shape) * sp_now / MAX_DIM * 1000) / 1000
    new_files = []
    for f, img, interp in zip(files, imgs, interps):
        if change_hd:
            ants.set_spacing(img, scaling)
            ants.set_origin(img, [X[i] * scaling[i] / S[i] for i in range(3)])
        if target is not None:
            img = ants.resample_image(img, (target, target, target), False, interp)
            img[img < 0] = 0
            if variant in LEGACY_GRID and 'desc-norm' in f:
                img = ants.smooth_image(img, target)
        nf = opj(prep, opb(f))
        ants.image_write(img, nf)
        new_files.append(nf)
    _, _, fwd, _ = E['get_orientation'].use_ants(new_files[0])
    for nf, n in zip(new_files, names):
        E['run_cmd'].run(sing_fs + 'mri_convert ' + fwd + nf + ' ' + opj(mri, n), diary)
    return ants.image_read(new_files[0]).shape, real_mm_per_fs


def light_config(E, FS_dir, sp, diary, sing_fs, export_fs):
    """First half of FS_surf.Wconfig (smoothwm, inflated, curvatures); no sphere registration / aparc."""
    rc = E['run_cmd']
    s = opj(FS_dir, sp, 'surf')
    for h in HEMIS:
        rc.run(sing_fs + 'mris_smooth -n 3 -nw ' + opj(s, h + '.white') + ' ' + opj(s, h + '.smoothwm'), diary)
        rc.run(sing_fs + 'mris_inflate ' + opj(s, h + '.smoothwm') + ' ' + opj(s, h + '.inflated'), diary)
        rc.run(sing_fs + 'mris_curvature -w ' + opj(s, h + '.white'), diary)
        rc.run(sing_fs + 'mris_curvature -thresh .999 -n -a 5 -w -distances 10 5 ' + opj(s, h + '.inflated'), diary)
        rc.do(export_fs + 'mris_curvature_stats -m --writeCurvatureFiles -G -o ' +
              opj(FS_dir, sp, 'stats', h + '.curv.stats') + ' -F smoothwm ' + sp + ' ' + h + ' curv sulc', diary)


def run_one(job):
    """One species x variant. Every stage is logged; a failure stops this job only."""
    sp, variant, a = job['species'], job['variant'], job['args']
    t0 = time.time()
    root = opj(a['out'], sp, variant)
    FS_dir = opj(root, 'freesurfer')
    prep_dir = opj(root, 'prep')
    res = dict(species=sp, variant=variant, status='ok', stage='', seconds=0, real_mm_per_fs_mm=1.0,
               grid='', note='')
    if opi(opj(root, 'done.json')) and not a['force']:
        return json.loads(open(opj(root, 'done.json')).read())
    if ope(root):
        shutil.rmtree(root)
    for d in ('mri', 'surf', 'stats', 'label', 'scripts'):
        os.makedirs(opj(FS_dir, sp, d), exist_ok=True)
    os.makedirs(prep_dir, exist_ok=True)
    diary = opj(root, 'diary.txt')
    open(diary, 'w').close()
    try:
        E = ednix_imports(a['main_path'])
        E['run_cmd'].STRICT = True
        import ants
        res['stage'] = 'paths'
        fs_tools, BASE_SS, Aseg_ref = library_paths(E, sp, a['lib'])
        _, _, sing_fs, _, _, _, _, _ = E['req'].load_requirement(a['main_path'], a['lib'], a['out'], 'yes')
        _, _, _, export_fs = E['req'].FS(fs_tools, FS_dir, diary, sing_fs)

        res['stage'] = 'prepa'
        if variant.startswith('V2') or variant.startswith('V3'):
            prepa_v2(E, sp, BASE_SS, Aseg_ref, prep_dir, const=variant.startswith('V3'))
        else:
            E['preFS'].prepa_img(sp, BASE_SS, prep_dir, Aseg_ref, prep_dir, '', diary)
        files = [BASE_SS, opj(prep_dir, sp + '_desc-norm_T1w.nii.gz'), opj(prep_dir, sp + '_seg-wm_dseg.nii.gz'),
                 Aseg_ref, opj(prep_dir, sp + '_seg-filled_dseg.nii.gz')]

        res['stage'] = 'smallbrain'
        change_hd, scaling, legacy_target = 0, [1., 1., 1.], None
        hd = ants.image_header_info(BASE_SS)
        if sp in SMALLBRAIN_LIBRARY or (variant in ('V0m', 'V2m', 'V3m') and sp in M_SPECIES):
            change_hd, scaling, resamp, new_size, _ = E['smallbrain'].get(sp, prep_dir, BASE_SS, '', FACTOR_LIBRARY,
                                                                         diary)
        if change_hd:
            res['real_mm_per_fs_mm'] = float(hd['spacing'][0]) / float(scaling[0])
            legacy_target = float(new_size) if resamp == 1 else None
        res['stage'] = 'mri_convert'
        shape, res['real_mm_per_fs_mm'] = to_freesurfer(E, sp, files, FS_dir, change_hd, scaling, legacy_target,
                                                        variant, diary, sing_fs)
        res['grid'] = 'x'.join(str(int(s)) for s in shape)
        res['stage'] = 'white (FS_surf.Wcreate)'
        E['FS_surf'].Wcreate(FS_dir, sp, 'lr', diary, sing_fs, export_fs)
        res['stage'] = 'curvature'
        light_config(E, FS_dir, sp, diary, sing_fs, export_fs)
        res['stage'] = 'pial (FS_surf.Pcreate)'
        if variant[-1] in 'mn':
            # FS_surf.Pcreate + '-max': mris_make_surfaces searches the pial surface and caps the thickness at
            # 5 mm by default, i.e. 5 x real_mm_per_fs_mm real mm after the smallbrain header rescale.
            lim = MAX_REAL_THICKNESS if res['real_mm_per_fs_mm'] < 0.3 else 5.0
            mx = max(5.0, math.ceil(lim / res['real_mm_per_fs_mm']))
            res['note'] = f'pial -max {mx:g} FS mm'
            for h in HEMIS:
                E['run_cmd'].do(export_fs + 'mris_make_surfaces -white NOWRITE -aseg aseg -orig white -noaparc '
                                f'-mgz -max {mx:g} -T1 brain ' + sp + ' ' + h, diary)
        else:
            E['FS_surf'].Pcreate(sp, 'lr', diary, export_fs)
        for h in HEMIS:
            if not opi(opj(FS_dir, sp, 'surf', h + '.pial')) and not opi(opj(FS_dir, sp, 'surf', h + '.pial.T1')):
                raise RuntimeError(h + '.pial not written')
        res['stage'] = 'done'
    except Exception as e:
        res['status'] = 'failed'
        res['note'] = f'{type(e).__name__}: {str(e)[:300]}'
        with open(diary, 'a') as fh:
            fh.write('\n' + traceback.format_exc())
    res['seconds'] = round(time.time() - t0)
    with open(opj(root, 'done.json'), 'w') as fh:
        json.dump(res, fh, indent=1)
    return res


# ------------------------------------------------------------------------------------------------------------
def surf_path(FS_dir, sp, h, kind):
    s = opj(FS_dir, sp, 'surf')
    if kind == 'pial':
        for n in (h + '.pial', h + '.pial.T1'):
            if opi(opj(s, n)):
                return opj(s, n)
        return None
    p = opj(s, h + '.' + kind)
    return p if opi(p) else None


def thickness_real(FS_dir, sp, k):
    """Cortical thickness (FreeSurfer ?h.thickness, inside ?h.cortex.label if present) converted to real mm."""
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


def plane_segments(v, f, axis, c):
    """Segments where the triangle mesh (v in voxel coords) crosses the plane v[:, axis] == c."""
    d = v[:, axis] - c
    s = d[f] > 0
    cross = (s.sum(1) == 1) | (s.sum(1) == 2)
    f = f[cross]
    pts = []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        a_, b_ = f[:, i], f[:, j]
        m = (d[a_] > 0) != (d[b_] > 0)
        t = np.where(m, d[a_] / np.where(m, d[a_] - d[b_], 1), np.nan)
        pts.append(v[a_] + t[:, None] * (v[b_] - v[a_]))
    P = np.stack(pts, 1)                                    # faces x 3 edges x 3
    ok = ~np.isnan(P[:, :, 0])
    segs = np.stack([P[r][ok[r]][:2] for r in range(len(P)) if ok[r].sum() == 2]) if len(P) else np.zeros((0, 2, 3))
    keep = [k for k in range(3) if k != axis]
    return segs[:, :, keep] if len(segs) else np.zeros((0, 2, 2))


def qc_species(sp, out, rows):
    import nibabel as nib
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    rows = [r for r in rows if r['species'] == sp]
    if not rows:
        return
    planes = [('coronal 30 %', 2, 0.30), ('coronal 50 %', 2, 0.50), ('coronal 70 %', 2, 0.70), ('axial 55 %', 1, 0.55)]
    fig, axes = plt.subplots(len(rows), len(planes) + 1, figsize=(3.2 * (len(planes) + 1), 3.2 * len(rows)),
                             squeeze=False)
    for i, r in enumerate(rows):
        FS_dir = opj(out, sp, r['variant'], 'freesurfer')
        orig_f = opj(FS_dir, sp, 'mri', 'orig.mgz')
        label = f"{r['variant']}  {r['status']}\ngrid {r['grid']}"
        if not opi(orig_f):
            for ax in axes[i]:
                ax.axis('off')
            axes[i, 0].text(0.5, 0.5, label + '\n' + r['note'][:80], ha='center', va='center', fontsize=7,
                            transform=axes[i, 0].transAxes)
            continue
        img = nib.load(orig_f)
        vol = np.asarray(img.dataobj).astype(np.float32)
        inv = np.linalg.inv(img.header.get_vox2ras_tkr())
        nz = np.argwhere(vol > 0)
        lo, hi = nz.min(0), nz.max(0)
        meshes = {}
        for kind in ('white', 'pial'):
            for h in HEMIS:
                p = surf_path(FS_dir, sp, h, kind)
                if p:
                    vv, ff = nib.freesurfer.read_geometry(p)
                    vv = (inv @ np.c_[vv, np.ones(len(vv))].T).T[:, :3]
                    meshes[(kind, h)] = (vv, ff)
        for j, (title, axis, frac) in enumerate(planes):
            ax = axes[i, j]
            c = int(round(lo[axis] + frac * (hi[axis] - lo[axis])))
            sl = np.take(vol, c, axis=axis)
            ax.imshow(sl.T, cmap='gray', origin='upper', interpolation='nearest')
            for (kind, h), (vv, ff) in meshes.items():
                seg = plane_segments(vv, ff, axis, c)
                if len(seg):
                    ax.add_collection(LineCollection(seg, colors='tab:blue' if kind == 'white' else 'tab:red',
                                                     linewidths=0.6))
            keep = [k for k in range(3) if k != axis]
            ax.set_xlim(lo[keep[0]] - 3, hi[keep[0]] + 3)
            ax.set_ylim(hi[keep[1]] + 3, lo[keep[1]] - 3)
            ax.set_xticks([]); ax.set_yticks([])
            if i == 0:
                ax.set_title(title, fontsize=9)
        axes[i, 0].set_ylabel(label, fontsize=8)
        ax = axes[i, -1]
        t = thickness_real(FS_dir, sp, r['real_mm_per_fs_mm'])
        if len(t):
            ax.hist(t, bins=60, color='tab:gray')
            ax.axvline(np.median(t), color='k', lw=1)
            ax.set_title(f'median {np.median(t):.2f} mm (real)', fontsize=8)
        ax.set_xlabel('thickness (real mm)', fontsize=8)
        ax.tick_params(labelsize=7)
    fig.suptitle(f'{sp} template: white (blue) / pial (red) on orig.mgz; rows = variants', fontsize=11)
    fig.tight_layout()
    fig.savefig(opj(out, f'QC_{sp}.png'), dpi=130)
    plt.close(fig)


def aseg_fidelity(FS_dir, sp, k):
    """Distance of the white / pial vertices to the boundaries they should follow, in real mm:
    white -> boundary of filled.mgz > 0; pial -> outer boundary of the aseg brain (aseg > 0 without CSF codes,
    holes filled). Returns {white|pial: (median |d| mm, fraction of vertices with |d| <= 1 voxel)}."""
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
        sd = ndi.distance_transform_edt(~m) - ndi.distance_transform_edt(m)    # voxels, 0 on the boundary
        d = []
        for h in HEMIS:
            p = surf_path(FS_dir, sp, h, kind)
            if p:
                v = nib.freesurfer.read_geometry(p)[0]
                ijk = (inv @ np.c_[v, np.ones(len(v))].T)[:3]
                d.append(np.abs(ndi.map_coordinates(sd.astype(np.float32), ijk, order=1, mode='nearest')))
        if d:
            d = np.concatenate(d)
            out[kind] = (float(np.median(d)) * vox * k, float(np.mean(d <= 1.0)))
    return out


def write_summary(out, rows):
    keys = ['species', 'variant', 'status', 'stage', 'seconds', 'grid', 'real_mm_per_fs_mm', 'n_vertices_white',
            'thick_median_mm', 'thick_p05_mm', 'thick_p95_mm', 'white_dist_med_mm', 'white_within_1vox',
            'pial_dist_med_mm', 'pial_within_1vox', 'note']
    import nibabel as nib
    for r in rows:
        FS_dir = opj(out, r['species'], r['variant'], 'freesurfer')
        n = 0
        for h in HEMIS:
            p = surf_path(FS_dir, r['species'], h, 'white')
            if p:
                n += len(nib.freesurfer.read_geometry(p)[0])
        r['n_vertices_white'] = n
        t = thickness_real(FS_dir, r['species'], r['real_mm_per_fs_mm'])
        for k, q in (('thick_median_mm', 50), ('thick_p05_mm', 5), ('thick_p95_mm', 95)):
            r[k] = round(float(np.percentile(t, q)), 3) if len(t) else ''
        try:
            fid = aseg_fidelity(FS_dir, r['species'], r['real_mm_per_fs_mm'])
        except Exception as e:
            print(f"aseg fidelity {r['species']} {r['variant']} failed: {e}")
            fid = {}
        for kind in ('white', 'pial'):
            md, fr = fid.get(kind, ('', ''))
            r[kind + '_dist_med_mm'] = round(md, 3) if md != '' else ''
            r[kind + '_within_1vox'] = round(fr, 3) if fr != '' else ''
    with open(opj(out, 'summary.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, keys, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    print(f"\n{'species':<12}{'var':<5}{'status':<8}{'grid':<14}{'min':>6}{'thick med':>10}{'p05-p95':>13}"
          f"{'pial<=1vox':>11}  note")
    for r in rows:
        print(f"{r['species']:<12}{r['variant']:<5}{r['status']:<8}{r['grid']:<14}{r['seconds'] // 60:>6}"
              f"{str(r['thick_median_mm']):>10}{(str(r['thick_p05_mm']) + '-' + str(r['thick_p95_mm'])):>13}"
              f"{str(r['pial_within_1vox']):>11}  {r['note'][:60]}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--main-path', default=MAIN_PATH_DEFAULT, help='EDNiX code folder (Atlases_library is next to it)')
    ap.add_argument('--out', default=OUT_DEFAULT)
    ap.add_argument('--species', nargs='+', default=SPECIES_DEFAULT)
    ap.add_argument('--variants', nargs='+', default=VARIANTS, choices=VARIANTS)
    ap.add_argument('--jobs', type=int, default=6, help='species x variant runs in parallel')
    ap.add_argument('--force', action='store_true', help='redo runs that already have done.json')
    ap.add_argument('--qc-only', action='store_true')
    args = ap.parse_args()

    args.main_path = args.main_path.rstrip('/') + '/'
    lib = opj(opd(args.main_path), 'Atlases_library')       # same rule as atlases/template.py
    os.makedirs(args.out, exist_ok=True)
    E = ednix_imports(args.main_path)
    species = []
    for sp in args.species:
        try:
            _, ss, aseg = library_paths(E, sp, lib)
            if opi(ss) and opi(aseg):
                species.append(sp)
            else:
                print(f'skip {sp}: missing {ss if not opi(ss) else aseg}')
        except Exception as e:
            print(f'skip {sp}: {e}')
    jobs = []
    for sp in species:
        for v in args.variants:
            if v == 'V1' and sp not in SMALLBRAIN_LIBRARY:
                continue
            if v in ('V0m', 'V2m', 'V3m') and sp not in M_SPECIES:
                continue        # no smallbrain rescale / down-sampling for these species
            if v in ('V0n', 'V2n') and sp not in SMALLBRAIN_LIBRARY and sp not in ('Marmoset', 'Macaque'):
                continue        # V3n runs for every species (cross-species comparison of the same rule)        # no rescale / down-sampling in V0 for these species: V1 would be V0 + crop
            jobs.append(dict(species=sp, variant=v, args=dict(main_path=args.main_path, lib=lib, out=args.out,
                                                                force=args.force)))
    print(f'{len(jobs)} runs: ' + ', '.join(f"{j['species']}/{j['variant']}" for j in jobs))

    rows = []
    if args.qc_only:
        for j in jobs:
            f = opj(args.out, j['species'], j['variant'], 'done.json')
            if opi(f):
                rows.append(json.loads(open(f).read()))
    else:
        with ProcessPoolExecutor(max(1, args.jobs)) as ex:
            futs = {ex.submit(run_one, j): j for j in jobs}
            for fu in as_completed(futs):
                r = fu.result()
                rows.append(r)
                print(f"[{len(rows)}/{len(jobs)}] {r['species']} {r['variant']}: {r['status']} "
                      f"({r['seconds'] // 60} min) {r['note'][:120]}", flush=True)
    order = {v: i for i, v in enumerate(VARIANTS)}
    rows.sort(key=lambda r: (species.index(r['species']) if r['species'] in species else 99, order[r['variant']]))
    write_summary(args.out, rows)
    for sp in species:
        try:
            qc_species(sp, args.out, rows)
        except Exception as e:
            print(f'QC figure {sp} failed: {e}')
    print(f'\noutputs: {args.out}/summary.csv, {args.out}/QC_<species>.png')


if __name__ == '__main__':
    main()
