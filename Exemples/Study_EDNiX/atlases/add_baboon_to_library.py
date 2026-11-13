"""
Add the Baboon (Papio anubis) to the EDNiX atlas library, with the EDNiX atlases of the Macaque registered onto the
baboon template (interim library entry, to be replaced when a baboon-specific EDNiX atlas is published).

Inputs
  baboon template (WIP folder, already in the old library format, 0.6 mm, RAS):
      <SRC>/modified/native/volumes/Baboon_space-acpc_desc-SS_T1w.nii.gz          brain (registration target)
      <SRC>/modified/native/volumes/Baboon_space_acpc_desc-template_T1w.nii.gz    head  (renamed space-acpc)
      <SRC>/modified/native/volumes/masks/Baboon_mask.nii.gz
      <SRC>/modified/native/volumes/labels/Baboon_seg-Agaronyan2022_dseg.nii.gz   baboon atlas, copied as is
  macaque library entry: <LIBRARY>/atlas/mammals/primates/Macaque/EDNiX (SS T1w, mask, seg-4FS, EDNIxCSC, EDNIxCSCLR)

Output (same layout and names as the Macaque entry, one resolution, no fs_LR 32k / BALSA space)
  <LIBRARY>/atlas/mammals/primates/Baboon/
      label_code/                       LUTs of the propagated atlases (copied from Macaque) + Agaronyan2022
      EDNiX/volumes/                    Baboon_space-acpc_desc-{SS,template}_T1w.nii.gz
      EDNiX/volumes/labels/             Baboon_seg-{4FS,EDNIxCSC,EDNIxCSCLR,Agaronyan2022}_dseg.nii.gz
      EDNiX/volumes/masks/              Baboon_mask + White/Vent/Cerebral_Gray/Whole_Brain_Gray (+erod) + dilat
                                        (preFS.msk_RS, from the propagated seg-4FS); infra/supra layers at step surfaces
      EDNiX/freesurfer/, surfaces/Native_resol/, preprocessing/   atlases/template.modif (surface method V2h)
      EDNiX/registration_from_Macaque/  transforms, warped macaque, registration metrics, provenance json
      EDNiX/QC/                         label overlay + surface QC png

Registration (ANTsPy): Macaque SS T1w -> Baboon SS T1w, both masked and intensity-normalised; initial isotropic scaling
(cube root of the brain-mask volume ratio) about the centres of mass, then rigid + affine (Mattes MI) + SyN (CC, radius
4) = ANTs 'SyNCC' preset, restricted to the dilated baboon mask. Labels are warped with genericLabel interpolation and
cut to the baboon brain mask.

Steps (--step, default all, each one resumable):  prepare  register  labels  masks  surfaces  qc
  (+ 'wb': Workbench conversion only, on existing FreeSurfer surfaces)
  python add_baboon_to_library.py                       # everything
  python add_baboon_to_library.py --plan                # paths and checks only
  python add_baboon_to_library.py --step qc             # one step
  python add_baboon_to_library.py --force               # recompute the registration even if it exists
"""
import os
import sys
import json
import glob
import shutil
import argparse
import datetime
import importlib.util

import numpy as np
import nibabel as nib
import ants

opj = os.path.join
opd = os.path.dirname
opb = os.path.basename
opi = os.path.isfile

# ---------------------------------------------------------------------------------------------------------------------
# settings
# ---------------------------------------------------------------------------------------------------------------------
MAIN_PATH = '/home/cgarin/PycharmProjects/EDNiX/'
LIBRARY = opj(opd(MAIN_PATH), 'Atlases_library')     # same rule as atlases/template.modif (V2 since 2026-11-11 swap)
SRC = '/srv/benhalab2/CASCAD/EDNiX_project/Atlas/1_Anatomical_atlases/WIP/mammals/primates/Baboon'
GROUP = opj('atlas', 'mammals', 'primates')
SPECIES = 'Baboon'
REF_SPECIES = 'Macaque'
REFERENCE = 'EDNiX'

# atlases propagated from the macaque; EXTRA_MACAQUE_ATLASES can add macaque-specific ones (e.g. 'CHARM', 'SARM',
# 'D99v2'): they are then copied with their LUTs, but they are NOT baboon atlases (keep them out of any paper claim).
EDNIX_ATLASES = ['4FS', 'EDNIxCSC', 'EDNIxCSCLR']
EXTRA_MACAQUE_ATLASES = []
LUT_ALWAYS = ['FreeSurfer_ITK.txt', 'FreeSurfer_label.txt', 'AsegStatsLUT.txt']

N_THREADS = 8
SEED = 1


def paths(main_path=MAIN_PATH, library=LIBRARY, src=SRC):
    mac = opj(library, GROUP, REF_SPECIES)
    bab = opj(library, GROUP, SPECIES)
    src_vol = opj(src, 'modified', 'native', 'volumes')
    out_vol = opj(bab, REFERENCE, 'volumes')
    return dict(
        main_path=main_path, library=library, src=src,
        mac=mac, mac_vol=opj(mac, REFERENCE, 'volumes'),
        mac_SS=opj(mac, REFERENCE, 'volumes', f'{REF_SPECIES}_space-acpc_desc-SS_T1w.nii.gz'),
        mac_mask=opj(mac, REFERENCE, 'volumes', 'masks', f'{REF_SPECIES}_mask.nii.gz'),
        mac_labels=opj(mac, REFERENCE, 'volumes', 'labels'),
        mac_code=opj(mac, 'label_code'),
        src_SS=opj(src_vol, f'{SPECIES}_space-acpc_desc-SS_T1w.nii.gz'),
        src_template=opj(src_vol, f'{SPECIES}_space_acpc_desc-template_T1w.nii.gz'),
        src_mask=opj(src_vol, 'masks', f'{SPECIES}_mask.nii.gz'),
        src_agaronyan=opj(src_vol, 'labels', f'{SPECIES}_seg-Agaronyan2022_dseg.nii.gz'),
        src_agaronyan_code=opj(src, 'atlas', 'agaronyan'),
        bab=bab, ref=opj(bab, REFERENCE), vol=out_vol,
        labels=opj(out_vol, 'labels'), masks=opj(out_vol, 'masks'), code=opj(bab, 'label_code'),
        SS=opj(out_vol, f'{SPECIES}_space-acpc_desc-SS_T1w.nii.gz'),
        template=opj(out_vol, f'{SPECIES}_space-acpc_desc-template_T1w.nii.gz'),
        mask=opj(out_vol, 'masks', f'{SPECIES}_mask.nii.gz'),
        reg=opj(bab, REFERENCE, 'registration_from_' + REF_SPECIES),
        qc=opj(bab, REFERENCE, 'QC'),
        diary=opj(bab, REFERENCE, f'add_{SPECIES}_diary.txt'))


def log(P, msg):
    line = f'[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] {msg}'
    print(line, flush=True)
    os.makedirs(opd(P['diary']), exist_ok=True)
    with open(P['diary'], 'a') as f:
        f.write(line + '\n')


def write_json(path, d):
    with open(path, 'w') as f:
        json.dump(d, f, indent=2)


# ---------------------------------------------------------------------------------------------------------------------
# checks
# ---------------------------------------------------------------------------------------------------------------------
def check_inputs(P):
    need = ['mac_SS', 'mac_mask', 'src_SS', 'src_template', 'src_mask']
    missing = [f'{k}: {P[k]}' for k in need if not opi(P[k])]
    for a in EDNIX_ATLASES + EXTRA_MACAQUE_ATLASES:
        f = opj(P['mac_labels'], f'{REF_SPECIES}_seg-{a}_dseg.nii.gz')
        if not opi(f):
            missing.append(f'macaque atlas {a}: {f}')
    if missing:
        raise FileNotFoundError('missing input(s):\n  ' + '\n  '.join(missing))
    # templatefeat.get takes the shallowest folder named <species> under <library>/atlas: refuse a second one
    others = [d for d, dirs, _ in os.walk(opj(P['library'], 'atlas'))
              if opb(d) == SPECIES and 'freesurfer' not in d.split(os.sep) and os.path.normpath(d) != os.path.normpath(P['bab'])]
    if others:
        raise RuntimeError(f'other folder(s) named {SPECIES} in the library (templatefeat would be ambiguous): {others}')
    a = nib.load(P['src_SS'])
    for k in ('src_template', 'src_mask'):
        b = nib.load(P[k])
        if b.shape[:3] != a.shape[:3] or not np.allclose(b.affine, a.affine, atol=1e-3):
            raise ValueError(f'{P[k]} is not on the grid of {P["src_SS"]}')
    if opi(P['src_agaronyan']):
        b = nib.load(P['src_agaronyan'])
        if b.shape[:3] != a.shape[:3] or not np.allclose(b.affine, a.affine, atol=1e-3):
            raise ValueError(f'{P["src_agaronyan"]} is not on the grid of {P["src_SS"]}')


def plan(P):
    for k, v in P.items():
        print(f'{k:20s} {v}' + ('' if not k.startswith(('src_', 'mac_')) else ('   OK' if os.path.exists(v) else '   MISSING')))
    check_inputs(P)
    print(f'\natlases propagated: {EDNIX_ATLASES + EXTRA_MACAQUE_ATLASES}\ninputs OK.')


# ---------------------------------------------------------------------------------------------------------------------
# steps
# ---------------------------------------------------------------------------------------------------------------------
def step_prepare(P):
    """library folders, baboon volumes renamed into EDNiX nomenclature, binary brain mask"""
    check_inputs(P)
    for k in ('vol', 'labels', 'masks', 'code', 'reg', 'qc'):
        os.makedirs(P[k], exist_ok=True)
    os.makedirs(opj(P['ref'], 'surfaces', 'Native_resol'), exist_ok=True)
    ref = nib.load(P['src_SS'])
    for src, dst in ((P['src_SS'], P['SS']), (P['src_template'], P['template'])):
        im = nib.load(src)
        nib.save(nib.Nifti1Image(np.asarray(im.dataobj, dtype=np.float32), ref.affine, ref.header), dst)
    m = (np.asarray(nib.load(P['src_mask']).dataobj) > 0.5).astype(np.float32)
    nib.save(nib.Nifti1Image(m, ref.affine, ref.header), P['mask'])
    write_json(P['SS'].replace('.nii.gz', '.json'),
               {'Sources': [P['src_SS']], 'Description': 'Baboon template (WIP library folder), skull-stripped.'})
    write_json(P['template'].replace('.nii.gz', '.json'),
               {'Sources': [P['src_template']], 'Description': 'Baboon template (WIP library folder), head.'})
    write_json(P['mask'].replace('.nii.gz', '.json'),
               {'Sources': [P['src_mask']], 'Description': 'Baboon brain mask (binarised at 0.5).'})
    log(P, f'prepare: volumes and mask written in {P["vol"]} (brain {m.sum() * np.prod(ref.header.get_zooms()[:3]) / 1000:.1f} cm3)')


def _norm(img, mask):
    """masked, clipped to the 0.5-99.5th percentiles of the brain, scaled to [0, 1]"""
    a = img.numpy()
    m = mask.numpy() > 0
    lo, hi = np.percentile(a[m], [0.5, 99.5])
    return img.new_image_like(np.where(m, np.clip((a - lo) / max(hi - lo, 1e-6), 0, 1), 0).astype(np.float32))


def _initial_scaling(fixed_mask, moving_mask, out):
    """ANTs affine (fixed -> moving points) = isotropic scaling by the cube root of the brain-volume ratio about the
    fixed centre of mass, plus the centre-of-mass translation."""
    def com_and_vol(m):
        a = m.numpy() > 0
        idx = np.argwhere(a).mean(0)
        phys = np.array(m.origin) + np.array(m.direction) @ (np.array(m.spacing) * idx)   # ITK physical (LPS)
        return phys, a.sum() * np.prod(m.spacing)
    cf, vf = com_and_vol(fixed_mask)
    cm, vm = com_and_vol(moving_mask)
    s = (vm / vf) ** (1. / 3.)           # < 1: a baboon point maps onto a smaller macaque brain
    tx = ants.create_ants_transform(transform_type='AffineTransform', dimension=3,
                                    matrix=np.eye(3) * s, center=list(cf), translation=list(cm - cf))
    ants.write_transform(tx, out)
    return float(vf / vm), float(s)


def step_register(P, force=False):
    os.environ['ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS'] = str(N_THREADS)
    prefix = opj(P['reg'], f'{REF_SPECIES}_to_{SPECIES}_')
    affine, warp, iwarp = prefix + '0GenericAffine.mat', prefix + '1Warp.nii.gz', prefix + '1InverseWarp.nii.gz'
    if not force and all(opi(f) for f in (affine, warp, iwarp)):
        log(P, 'register: transforms exist, skipped (--force to recompute)')
        return
    fixed_mask = ants.image_read(P['mask'])
    moving_mask = ants.threshold_image(ants.image_read(P['mac_mask']), 0.5, 1e9)
    fixed = _norm(ants.image_read(P['SS']), fixed_mask)
    moving = _norm(ants.image_read(P['mac_SS']), moving_mask)
    init = prefix + 'init_scaling.mat'
    vol_ratio, s = _initial_scaling(fixed_mask, moving_mask, init)
    log(P, f'register: brain volume ratio {SPECIES}/{REF_SPECIES} = {vol_ratio:.2f} (initial scaling {1 / s:.3f})')
    reg_mask = ants.morphology(fixed_mask, operation='dilate', radius=3, mtype='binary', shape='ball')
    reg = ants.registration(fixed=fixed, moving=moving, type_of_transform='SyNCC', initial_transform=[init],
                            mask=reg_mask, aff_metric='mattes', outprefix=prefix, random_seed=SEED, verbose=False)
    fwd = reg['fwdtransforms']
    inv = reg['invtransforms']
    # ANTs collapses the initial transform into the output affine (-z 1): keep exactly two forward transforms
    if len(fwd) != 2:
        raise RuntimeError(f'unexpected forward transforms: {fwd}')
    for f, dst in zip(fwd, [warp, affine]):
        if os.path.abspath(f) != os.path.abspath(dst):
            shutil.copyfile(f, dst)
    iw = [f for f in inv if 'InverseWarp' in f][0]
    if os.path.abspath(iw) != os.path.abspath(iwarp):
        shutil.copyfile(iw, iwarp)
    warped = ants.apply_transforms(fixed, ants.image_read(P['mac_SS']), [warp, affine], interpolator='lanczosWindowedSinc')
    ants.image_write(warped, opj(P['reg'], f'{REF_SPECIES}_space-{SPECIES}_desc-SS_T1w.nii.gz'))
    wmask = ants.apply_transforms(fixed, moving_mask, [warp, affine], interpolator='genericLabel')
    ants.image_write(wmask, opj(P['reg'], f'{REF_SPECIES}_space-{SPECIES}_mask.nii.gz'))
    m = metrics(P)
    m.update(volume_ratio=vol_ratio, initial_scaling=1 / s, transform='SyNCC', interpolation_labels='genericLabel',
             fixed=P['SS'], moving=P['mac_SS'], seed=SEED, antspy=ants.__version__,
             date=datetime.datetime.now().isoformat(timespec='seconds'))
    write_json(opj(P['reg'], 'registration_metrics.json'), m)
    log(P, f'register: done, mask Dice {m["mask_dice"]:.3f}, NCC in brain {m["ncc_in_brain"]:.3f}')


def metrics(P):
    fm = nib.load(P['mask']).get_fdata() > 0.5
    wm = nib.load(opj(P['reg'], f'{REF_SPECIES}_space-{SPECIES}_mask.nii.gz')).get_fdata() > 0.5
    a = nib.load(P['SS']).get_fdata()[fm]
    b = nib.load(opj(P['reg'], f'{REF_SPECIES}_space-{SPECIES}_desc-SS_T1w.nii.gz')).get_fdata()[fm]
    ncc = float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else float('nan')
    return dict(mask_dice=float(2 * (fm & wm).sum() / (fm.sum() + wm.sum())), ncc_in_brain=ncc)


def _warp_labels(P, src, dst, fixed, transforms, brain):
    img = ants.image_read(src)
    src_dtype = nib.load(src).get_data_dtype()
    vols = ants.ndimage_to_list(img) if img.dimension == 4 else [img]
    out = []
    for v in vols:
        w = ants.apply_transforms(fixed, v, transforms, interpolator='genericLabel').numpy()
        out.append(np.where(brain, np.rint(w), 0))
    arr = np.stack(out, -1) if len(out) > 1 else out[0]
    dtype = src_dtype if (np.issubdtype(src_dtype, np.integer) and src_dtype.itemsize <= 4) else np.dtype(np.int32)
    ref = nib.load(P['SS'])
    hdr = ref.header.copy()
    hdr.set_data_dtype(dtype)
    im = nib.Nifti1Image(arr.astype(dtype), ref.affine, hdr)
    im.header.set_xyzt_units('mm')
    nib.save(im, dst)
    return arr


def step_labels(P):
    fixed = ants.image_read(P['SS'])
    prefix = opj(P['reg'], f'{REF_SPECIES}_to_{SPECIES}_')
    transforms = [prefix + '1Warp.nii.gz', prefix + '0GenericAffine.mat']
    brain = nib.load(P['mask']).get_fdata() > 0.5
    summary = {}
    for a in EDNIX_ATLASES + EXTRA_MACAQUE_ATLASES:
        src = opj(P['mac_labels'], f'{REF_SPECIES}_seg-{a}_dseg.nii.gz')
        dst = opj(P['labels'], f'{SPECIES}_seg-{a}_dseg.nii.gz')
        arr = _warp_labels(P, src, dst, fixed, transforms, brain)
        src_lab = np.asarray(nib.load(src).dataobj)
        first = (lambda x: x[..., 0] if x.ndim == 4 else x)
        lost = sorted(set(np.unique(first(src_lab)).astype(int)) - set(np.unique(first(arr)).astype(int)))
        summary[a] = dict(n_labels_macaque=int(len(np.unique(first(src_lab))) - 1),
                          n_labels_baboon=int(len(np.unique(first(arr))) - 1), labels_lost=lost)
        write_json(dst.replace('.nii.gz', '.json'),
                   {'Sources': [src, P['SS']] + transforms,
                    'Description': f'{REF_SPECIES} EDNiX atlas {a} registered onto the {SPECIES} template (ANTs SyNCC, '
                                   f'genericLabel interpolation, cut to the {SPECIES} brain mask). Interim atlas: '
                                   f'macaque anatomy, not a {SPECIES}-specific delineation.'})
        log(P, f'labels: {a} -> {opb(dst)} ({summary[a]["n_labels_baboon"]}/{summary[a]["n_labels_macaque"]} labels'
               + (f', lost: {lost}' if lost else '') + ')')
    # 4FS coverage of the brain mask
    aseg = np.asarray(nib.load(opj(P['labels'], f'{SPECIES}_seg-4FS_dseg.nii.gz')).dataobj)
    summary['4FS_unlabelled_fraction_of_brain_mask'] = float(((aseg == 0) & brain).sum() / brain.sum())

    # baboon-specific atlas already in the WIP folder: copied as is
    if opi(P['src_agaronyan']):
        dst = opj(P['labels'], f'{SPECIES}_seg-Agaronyan2022_dseg.nii.gz')
        shutil.copyfile(P['src_agaronyan'], dst)
        js = P['src_agaronyan'].replace('.nii.gz', '.json')
        if opi(js):
            shutil.copyfile(js, dst.replace('.nii.gz', '.json'))
        for f in glob.glob(opj(P['src_agaronyan_code'], 'Agaronyan2022_*.txt')):
            shutil.copyfile(f, opj(P['code'], opb(f)))
        log(P, 'labels: Agaronyan2022 copied (baboon atlas of the WIP folder)')

    # label codes of the propagated atlases
    n = 0
    for f in sorted(os.listdir(P['mac_code'])):
        if f.startswith('.~lock'):
            continue
        if f in LUT_ALWAYS or any(f.startswith(a + '_') for a in EDNIX_ATLASES + EXTRA_MACAQUE_ATLASES):
            shutil.copyfile(opj(P['mac_code'], f), opj(P['code'], f))
            n += 1
    log(P, f'labels: {n} label_code files copied from {REF_SPECIES}')
    write_json(opj(P['reg'], 'labels_summary.json'), summary)


def step_masks(P):
    sys.path.insert(1, P['main_path'])
    from modalities.anat.freesurfer import preFS
    aseg = opj(P['labels'], f'{SPECIES}_seg-4FS_dseg.nii.gz')
    preFS.msk_RS(SPECIES, P['mask'], aseg, P['diary'], 'aseg')
    log(P, 'masks: White / Vent / Cerebral_Gray / Whole_Brain_Gray (+ erod) and dilat written')


def _load_build_library(P):
    f = opj(P['main_path'], 'Exemples', 'Study_EDNiX', 'atlases', 'build_library_V2.py')
    spec = importlib.util.spec_from_file_location('build_library_V2', f)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def step_surfaces(P, wb_only=False):
    """FreeSurfer + Workbench surfaces of the template, as for every library species (build_library_V2 settings,
    surface method V2h). No BALSA / fs_LR space exists for the baboon in templatefeat: native resolution only.
    wb_only=True keeps the FreeSurfer surfaces and redoes only the FreeSurfer -> Workbench conversion (FS_2_WB_short),
    as fix_library_wb_surfaces.py does."""
    bl = _load_build_library(P)
    if SPECIES in bl.SMALLBRAIN:
        raise RuntimeError('Baboon must not be in SMALLBRAIN')
    # FS_2_WB_short writes into surfaces/Native_resol without creating it (it exists for the V1 species)
    os.makedirs(opj(P['ref'], 'surfaces', 'Native_resol'), exist_ok=True)
    if wb_only:
        if not opi(opj(P['ref'], 'freesurfer', SPECIES, 'surf', 'rh.mid')):
            raise FileNotFoundError('FreeSurfer surfaces incomplete (no rh.mid): run --step surfaces')
        os.environ.setdefault('MPLBACKEND', 'Agg')
        sys.path.insert(1, P['main_path'])
        from atlases import template
        from modalities.anat.freesurfer import preFS
        preFS.SURFACE_METHOD = bl.SURFACE['method']
        log(P, 'wb: FreeSurfer -> Workbench conversion only (template.modif wb_only=True)')
        template.modif(bl.SURFACE['doBALSAspace'], bl.SMALLBRAIN, SPECIES, REFERENCE, bl.ATLAS_FOLLOWERS, P['main_path'],
                       P['bab'], P['diary'], proj=bl.SURFACE['proj'], doflat=bl.SURFACE['doflat'],
                       library_path=P['library'], smallbrain_factor=bl.SURFACE['smallbrain_factor'],
                       sphere=bl.SURFACE['sphere'], wb_only=True)
    else:
        log(P, f'surfaces: atlases/template.modif, SURFACE={bl.SURFACE}')
        minutes = bl.build_species(SPECIES, REFERENCE, P['bab'], P['library'], bl.SURFACE) / 60
        log(P, f'surfaces: done in {minutes:.0f} min')
    for h in ('l', 'r'):
        f = opj(P['ref'], 'surfaces', 'Native_resol', f'{SPECIES}.{h}.midthickness.surf.gii')
        if not opi(f):
            raise RuntimeError(f'Workbench surfaces missing ({f}): see {P["diary"]} and the build diary')
    try:
        q = bl.species_qc(SPECIES, REFERENCE, P['bab'], P['qc'])
        write_json(opj(P['qc'], 'surface_QC.json'), q if isinstance(q, dict) else {'result': str(q)})
    except Exception as e:
        log(P, f'surfaces: surface QC failed ({type(e).__name__}: {e})')


def step_qc(P):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    T1 = nib.as_closest_canonical(nib.load(P['SS'])).get_fdata()
    wmac = nib.as_closest_canonical(nib.load(opj(P['reg'], f'{REF_SPECIES}_space-{SPECIES}_desc-SS_T1w.nii.gz'))).get_fdata()
    aseg = np.asarray(nib.as_closest_canonical(nib.load(opj(P['labels'], f'{SPECIES}_seg-4FS_dseg.nii.gz'))).dataobj)
    csc = nib.as_closest_canonical(nib.load(opj(P['labels'], f'{SPECIES}_seg-EDNIxCSCLR_dseg.nii.gz')))
    csc = np.asarray(csc.dataobj)
    csc = csc[..., -1] if csc.ndim == 4 else csc
    brain = nib.as_closest_canonical(nib.load(P['mask'])).get_fdata() > 0.5
    ijk = np.argwhere(brain)
    lo, hi = ijk.min(0), ijk.max(0)
    rng = np.random.default_rng(0)
    cmap = ListedColormap(rng.uniform(0.15, 1, (512, 3)))

    def sl(a, ax, f):
        i = int(lo[ax] + f * (hi[ax] - lo[ax]))
        s = np.take(a, i, axis=ax)
        s = s[lo[1 if ax == 0 else 0]:hi[1 if ax == 0 else 0] + 1, lo[1 if ax == 2 else 2]:hi[1 if ax == 2 else 2] + 1]
        return np.rot90(s)

    vmax = np.percentile(T1[brain], 99.5)
    vmax2 = np.percentile(wmac[brain], 99.5)
    tissue = {'cortex': (3, 42), 'WM': (2, 41)}
    planes = [(0, 0.5, 'sagittal 50%'), (1, 0.35, 'coronal 35%'), (1, 0.6, 'coronal 60%'), (2, 0.55, 'axial 55%')]
    fig, axs = plt.subplots(3, len(planes), figsize=(4.2 * len(planes), 12.5), facecolor='k')
    for c, (ax, f, name) in enumerate(planes):
        t = sl(T1, ax, f)
        a = sl(aseg, ax, f)
        axs[0, c].imshow(t, cmap='gray', vmin=0, vmax=vmax)
        for (lab, col) in zip(tissue.values(), ('#ffd400', '#00e5ff')):
            axs[0, c].contour(np.isin(a, lab).astype(float), [0.5], colors=col, linewidths=0.7)
        axs[1, c].imshow(sl(wmac, ax, f), cmap='gray', vmin=0, vmax=vmax2)
        axs[1, c].contour(np.isin(a, tissue['WM']).astype(float), [0.5], colors='#00e5ff', linewidths=0.7)
        axs[2, c].imshow(t, cmap='gray', vmin=0, vmax=vmax)
        lab = sl(csc, ax, f).astype(float)
        axs[2, c].imshow(np.ma.masked_where(lab == 0, lab % 512), cmap=cmap, vmin=0, vmax=511, alpha=0.55,
                         interpolation='nearest')
        axs[0, c].set_title(name, color='w')
    for r, txt in enumerate([f'{SPECIES} T1w + seg-4FS (yellow cortex, cyan WM)',
                             f'{REF_SPECIES} T1w warped + {SPECIES} WM contour',
                             f'{SPECIES} T1w + EDNIxCSCLR (last level)']):
        axs[r, 0].set_ylabel(txt, color='w', fontsize=9)
    for a in axs.ravel():
        a.set_xticks([]); a.set_yticks([])
        for s in a.spines.values():
            s.set_visible(False)
    m = {}
    if opi(opj(P['reg'], 'registration_metrics.json')):
        with open(opj(P['reg'], 'registration_metrics.json')) as fh:
            m = json.load(fh)
    fig.suptitle(f'{REF_SPECIES} EDNiX atlases on the {SPECIES} template - mask Dice {m.get("mask_dice", float("nan")):.3f}, '
                 f'NCC {m.get("ncc_in_brain", float("nan")):.3f}', color='w')
    fig.tight_layout()
    out = opj(P['qc'], f'{SPECIES}_from_{REF_SPECIES}_labels_QC.png')
    fig.savefig(out, dpi=110, facecolor='k')
    plt.close(fig)
    log(P, f'qc: {out}')


def provenance(P):
    write_json(opj(P['ref'], f'{SPECIES}_EDNiX_provenance.json'), {
        'species': SPECIES, 'reference': REFERENCE,
        'status': f'interim: {REF_SPECIES} EDNiX atlases registered onto the {SPECIES} template; replace by a '
                  f'{SPECIES}-specific atlas when available',
        'template_source': P['src'], 'macaque_source': P['mac'],
        'atlases_propagated': EDNIX_ATLASES + EXTRA_MACAQUE_ATLASES,
        'registration': 'ANTsPy SyNCC (rigid+affine Mattes, SyN CC r=4) after volume-ratio scaling, dilated baboon mask',
        'script': os.path.abspath(__file__), 'date': datetime.datetime.now().isoformat(timespec='seconds')})


STEPS = ['prepare', 'register', 'labels', 'masks', 'surfaces', 'qc']

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--step', nargs='+', choices=STEPS + ['wb', 'all'], default=['all'],
                    help="'wb' = Workbench conversion only, on existing FreeSurfer surfaces")
    ap.add_argument('--plan', action='store_true', help='print paths, check inputs, run nothing')
    ap.add_argument('--force', action='store_true', help='recompute the registration')
    ap.add_argument('--main-path', default=MAIN_PATH)
    ap.add_argument('--library', default=None, help='default: <main-path>/../Atlases_library')
    ap.add_argument('--src', default=SRC)
    a = ap.parse_args()
    P = paths(a.main_path, a.library or opj(opd(a.main_path), 'Atlases_library'), a.src)
    if a.plan:
        plan(P)
        sys.exit(0)
    steps = STEPS if 'all' in a.step else [s for s in STEPS + ['wb'] if s in a.step]
    for s in steps:
        log(P, f'===== step {s}')
        if s == 'prepare':
            step_prepare(P)
        elif s == 'register':
            step_register(P, a.force)
        elif s == 'labels':
            step_labels(P)
        elif s == 'masks':
            step_masks(P)
        elif s == 'surfaces':
            step_surfaces(P)
        elif s == 'wb':
            step_surfaces(P, wb_only=True)
        elif s == 'qc':
            step_qc(P)
    provenance(P)
    log(P, 'finished: ' + ' '.join(steps))
