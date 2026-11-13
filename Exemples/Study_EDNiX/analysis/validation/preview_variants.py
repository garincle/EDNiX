"""
Preliminary check of calibrate_registration_variants.py WHILE it runs (read-only on its files).

For every subject folder: progress (which registrations exist), label Dice vs recon-all aseg of every variant already
computable (A, B from the cache; C, D when their warp exists; E always, DE when D exists), and a QC figure
<OUT_DIR>/QC_preview/<subject>_preview.png (one row per available variant). Nothing of the main run is modified.

    cd ~/Downloads && python preview_variants.py
(calibrate_registration_variants.py must be in the same folder: its functions are reused.)
"""
import os
import sys
import glob
import importlib.util

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('variants', os.path.join(HERE, 'calibrate_registration_variants.py'))
V = importlib.util.module_from_spec(spec)
spec.loader.exec_module(V)
opj, opi = os.path.join, os.path.isfile


def template():
    """Template files without rewriting the derived priors the running job is using."""
    hits = sorted(glob.glob(opj(V.LIBRARY, 'atlas', '**', V.SPECIES, V.REFERENCE, 'volumes',
                                f'{V.SPECIES}_space-acpc_desc-SS_T1w.nii.gz'), recursive=True),
                  key=lambda p: (p.count(os.sep), p))
    vol = os.path.dirname([h for h in hits if f'{os.sep}freesurfer{os.sep}' not in h][0])
    f = dict(t1=opj(vol, f'{V.SPECIES}_space-acpc_desc-SS_T1w.nii.gz'),
             seg=opj(vol, 'labels', f'{V.SPECIES}_seg-4FS_dseg.nii.gz'))
    for k in ('CSF', 'GM', 'WM'):
        lib = opj(vol, 'priors', f'{V.SPECIES}_label-{k}_probseg.nii.gz')
        der = opj(V.OUT_DIR, 'template_priors', f'{V.SPECIES}_label-{k}_probseg_from-seg4FS.nii.gz')
        f[k] = lib if opi(lib) else der
        if not opi(f[k]):
            sys.exit(f'No prior for {k} yet ({der}): wait until the main run has started.')
    return f


def preview(sub_dir, tpl):
    import ants
    fsid = os.path.basename(sub_dir.rstrip('/'))
    c = opj(V.CACHE_DIR, fsid)
    aseg = opj(c, 'fs_aseg_in_ednix_native.nii.gz')
    t1s = sorted(t for t in glob.glob(opj(V.BIDS_DIR, f'sub-{fsid[4:].split("_ses-")[0]}',
                                          f'ses-{fsid.split("_ses-")[1]}', 'anat', 'native', 'volumes',
                                          '*_space-acpc_desc-SS_*T1w.nii.gz')) if 'desc-SS-step' not in t)
    if not (opi(aseg) and t1s and opi(opj(c, 'A_1Warp.nii.gz'))):
        return fsid, {}, 'A or recon-all aseg not ready'
    moving = ants.image_read(t1s[0])
    mmask = moving.numpy() > 0
    seg_t = ants.image_read(tpl['seg'])
    pri = [ants.image_read(tpl[k]) for k in ('CSF', 'GM', 'WM')]
    ref = ants.image_read(aseg).numpy().astype(int)

    aA, wA, iA = (opj(c, 'A_' + s) for s in ('0GenericAffine.mat', '1Warp.nii.gz', '1InverseWarp.nii.gz'))
    chains = {'A': ([aA, iA], [True, False])}
    if opi(opj(c, 'B_1InverseWarp.nii.gz')):
        chains['B'] = ([opj(c, 'B_0GenericAffine.mat'), opj(c, 'B_1InverseWarp.nii.gz')], [True, False])
    if opi(opj(sub_dir, 'C_1InverseWarp.nii.gz')):
        chains['C'] = ([aA, iA, opj(sub_dir, 'C_1InverseWarp.nii.gz')], [True, False, False])
    if opi(opj(sub_dir, 'D_1InverseWarp.nii.gz')):
        chains['D'] = ([opj(sub_dir, 'D_0GenericAffine.mat'), opj(sub_dir, 'D_1InverseWarp.nii.gz')], [True, False])

    def to_native(img, ch, interp):
        return ants.apply_transforms(moving, img, ch[0], whichtoinvert=ch[1], interpolator=interp).numpy()

    segs = {v: to_native(seg_t, ch, 'nearestNeighbor').astype(int) for v, ch in chains.items()}
    for base, name in (('A', 'E'), ('D', 'DE')):
        if base in chains:
            post = V.gmm_posteriors(moving.numpy(), [to_native(p, chains[base], 'linear') for p in pri], mmask)
            segs[name] = V.relabel_ribbon(segs[base], post, mmask)
    order = [v for v in V.VARIANTS if v in segs]
    res = {v: V.label_dice(ref, segs[v]) for v in order}
    out = V.qc_figure(fsid, t1s[0], ref, {v: segs[v] for v in order}, res,
                      opj(V.OUT_DIR, 'QC_preview', f'{fsid}_preview.png'), note=' (preview)')
    return fsid, res, out


def main():
    tpl = template()
    dirs = sorted(glob.glob(opj(V.OUT_DIR, 'sub-*/')))
    print('progress (cache A/B, run C/D/result):')
    for d in dirs:
        s = os.path.basename(d.rstrip('/'))
        c = opj(V.CACHE_DIR, s)
        flags = [opi(opj(c, 'A_1Warp.nii.gz')), opi(opj(c, 'B_1Warp.nii.gz')), opi(opj(d, 'C_1Warp.nii.gz')),
                 opi(opj(d, 'D_1Warp.nii.gz')), opi(opj(d, 'result.json'))]
        print(f'  {s}: ' + ' '.join(f'{n}:{"ok" if f else "-"}' for n, f in zip(('A', 'B', 'C', 'D', 'done'), flags)))
    allres = {}
    print('\nlabel Dice vs recon-all aseg (mean of 12 structures | cortex):')
    for d in dirs:
        try:
            fsid, res, out = preview(d, tpl)
        except Exception as e:  # warp still being written by the main run
            fsid, res, out = os.path.basename(d.rstrip('/')), {}, f'skipped ({type(e).__name__}: {e})'
        if not res:
            print(f'  {fsid}: {out}')
            continue
        allres[fsid] = res
        print(f'  {fsid}: ' + ' | '.join(f"{v} {np.nanmean(list(r.values())):.3f}/{r['cortex']:.3f}"
                                         for v, r in res.items()) + f'\n     -> {out}')
    if allres:
        print('\nmean over subjects (delta vs A):')
        for v in V.VARIANTS:
            vals = [(np.nanmean(list(r[v].values())), np.nanmean(list(r['A'].values()))) for r in allres.values()
                    if v in r]
            if vals:
                x = np.array(vals)
                print(f'  {v:3s} n={len(x)}  mean {x[:, 0].mean():.3f}  delta {np.mean(x[:, 0] - x[:, 1]):+.3f}  '
                      f'better than A in {(x[:, 0] > x[:, 1]).sum()}/{len(x)}')


if __name__ == '__main__':
    main()
