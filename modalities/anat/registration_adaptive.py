"""
'SyNadaptive': anat -> template registration that adapts to the anatomical variability of the subject.

Rationale
---------
Inter-individual variability of the human brain (sulcal pattern, ventricle size, atrophy) is larger than in the
non-human species EDNiX was tuned on. A single fixed SyN setting is either too stiff for the atypical brains or too
compliant (folding, noise-driven warps) for the typical ones. SyNadaptive therefore:

  A. runs the standard EDNiX 'SyNCC' (rigid -> affine -> SyN[0.15,3,0], CC radius 4; Avants et al., 2008), initialised
     with the translation computed by norm2template.norm. When A passes the checks, outputs are IDENTICAL to SyNCC.
  B. only if A fails the quality checks, refines A with a more compliant SyN (SyN[0.2,2,0], CC radius 4) restricted
     to the brain masks at every stage (cost-function masking, Brett et al., 2001), starting from A's affine.
     B replaces A only if the tissue overlap improves by more than `min_gain` without folding of the warp.

Quality checks (computed in the fixed / template space, inside the template brain mask):
  - mask_dice   : Dice between the warped subject brain mask and the template mask,
  - tissue_dice : mean Dice of a 3-class k-means segmentation (CSF/GM/WM-like) of the warped subject vs the template,
                  a label-free proxy of anatomical overlap (Klein et al., 2009 evaluate registration by label overlap),
  - folding     : fraction of voxels with Jacobian determinant <= 0 (non-diffeomorphic voxels).

Every decision and metric is written to <outprefix>adaptive_QC.json, so the choice is traceable per subject.

The files written are the standard ANTs names expected downstream
(<outprefix>0GenericAffine.mat, <outprefix>1Warp.nii.gz, <outprefix>1InverseWarp.nii.gz): the rest of the pipeline does
not need to know which candidate was kept. ANTs collapses the initial translation into 0GenericAffine.mat
(collapse-output-transforms), for A as for B.

References: Avants et al., Med Image Anal 2008; Klein et al., NeuroImage 2009; Brett et al., NeuroImage 2001.
"""
import os
import json
import shutil
import numpy as np
import ants

opj = os.path.join
opi = os.path.isfile

NAME = 'SyNadaptive'

# Default thresholds (tunable through register(..., thresholds=...)).
THRESHOLDS = dict(tissue_dice=0.85, mask_dice=0.95, max_folding=1e-4, min_gain=0.01)

# Candidate B: more compliant SyN, masked at every stage.
REFINE = dict(syn_metric='CC', syn_sampling=4, grad_step=0.2, flow_sigma=2, total_sigma=0,
              reg_iterations=(100, 70, 50, 20))


def ants_type(type_of_transform):
    """ANTsPy has no 'SyNadaptive': the places that call ants.registration directly use its first candidate."""
    return 'SyNCC' if type_of_transform == NAME else type_of_transform


def _mask(img, mask):
    if mask is not None:
        return ants.threshold_image(mask, 0.5, 1e9, 1, 0)
    return ants.get_mask(img)


def _dice(a, b, label=1):
    a = a.numpy() == label
    b = b.numpy() == label
    s = a.sum() + b.sum()
    return float(2 * (a & b).sum() / s) if s else np.nan


def _tissue_dice(fixed, warped, mask, k=3):
    """Mean Dice of k-means classes (sorted by intensity, so class i means the same tissue in both images)."""
    seg_f = ants.kmeans_segmentation(fixed, k, kmask=mask)['segmentation']
    seg_w = ants.kmeans_segmentation(warped, k, kmask=mask)['segmentation']
    return float(np.nanmean([_dice(seg_f, seg_w, c) for c in range(1, k + 1)]))


def _folding(fixed, warp_file, mask):
    jac = ants.create_jacobian_determinant_image(fixed, warp_file, do_log=False)
    m = mask.numpy() > 0
    return float((jac.numpy()[m] <= 0).mean()) if m.any() else np.nan


def evaluate(fixed, moving, fwdtransforms, fixed_mask, moving_mask):
    """Quality metrics of a registration moving -> fixed, inside the fixed brain mask."""
    fmask = _mask(fixed, fixed_mask)
    mmask = _mask(moving, moving_mask)
    warped = ants.apply_transforms(fixed=fixed, moving=moving, transformlist=fwdtransforms, interpolator='linear')
    wmask = ants.apply_transforms(fixed=fixed, moving=mmask, transformlist=fwdtransforms,
                                  interpolator='nearestNeighbor')
    warp = [t for t in fwdtransforms if t.endswith('.nii.gz')]
    return dict(mask_dice=_dice(fmask, ants.threshold_image(wmask, 0.5, 1e9, 1, 0)),
                tissue_dice=_tissue_dice(fixed, warped, fmask),
                folding=_folding(fixed, warp[0], fmask) if warp else 0.0)


def passes(m, th):
    return m['tissue_dice'] >= th['tissue_dice'] and m['mask_dice'] >= th['mask_dice'] \
        and m['folding'] <= th['max_folding']


def _result(prefix):
    aff, warp, inv = prefix + '0GenericAffine.mat', prefix + '1Warp.nii.gz', prefix + '1InverseWarp.nii.gz'
    return {'fwdtransforms': [warp, aff], 'invtransforms': [aff, inv]}


def refine(fixed, moving, affine, outprefix, fixed_mask=None, moving_mask=None, params=None):
    """
    Candidate B: compliant SyN from `affine`, restricted to the brain masks at every stage. When no mask is given
    (EDNiX step 6 registers skull-stripped images without masks), the non-zero voxels of the images are used.
    Returns the ANTsPy-like dict of its files (<outprefix>0GenericAffine.mat, 1Warp, 1InverseWarp).
    """
    fm = _mask(fixed, fixed_mask)
    mm = _mask(moving, moving_mask)
    ants.registration(fixed=fixed, moving=moving, type_of_transform='SyNOnly', initial_transform=affine,
                      mask=fm, moving_mask=mm, mask_all_stages=True, outprefix=outprefix, **(params or REFINE))
    return _result(outprefix)


def register(fixed, moving, initial_transform, outprefix, fixed_mask=None, moving_mask=None,
             thresholds=None, log=print):
    """
    fixed, moving     : ANTsImage (fixed = template / reference, moving = subject), as in norm2template.norm
    initial_transform : list of transforms (the translation of norm2template.norm)
    outprefix         : prefix of the standard output files
    Returns the ANTsPy-like dict {'fwdtransforms': [warp, affine], 'invtransforms': [affine, inverse warp]}.
    """
    th = dict(THRESHOLDS, **(thresholds or {}))
    report = dict(method=NAME, thresholds=th, refine_parameters={k: list(v) if isinstance(v, tuple) else v
                                                                  for k, v in REFINE.items()})

    # A: the standard EDNiX SyNCC, written directly under the standard names.
    ants.registration(fixed=fixed, moving=moving, type_of_transform='SyNCC', initial_transform=initial_transform,
                      mask=fixed_mask, moving_mask=moving_mask, outprefix=outprefix)
    res = _result(outprefix)
    mA = evaluate(fixed, moving, res['fwdtransforms'], fixed_mask, moving_mask)
    report['A_SyNCC'] = mA
    report['chosen'] = 'A_SyNCC'
    log(f'INFO: {NAME} candidate A (SyNCC): ' + ', '.join(f'{k}={v:.4f}' for k, v in mA.items()))

    if passes(mA, th):
        report['decision'] = 'A passes the checks: kept (identical to SyNCC).'
    else:
        resB = refine(fixed, moving, res['fwdtransforms'][1], outprefix + 'refine_', fixed_mask, moving_mask)
        mB = evaluate(fixed, moving, resB['fwdtransforms'], fixed_mask, moving_mask)
        report['B_refine'] = mB
        log(f'INFO: {NAME} candidate B (compliant masked SyN): ' + ', '.join(f'{k}={v:.4f}' for k, v in mB.items()))
        if mB['tissue_dice'] - mA['tissue_dice'] > th['min_gain'] and mB['folding'] <= th['max_folding']:
            for a, b in zip(res['fwdtransforms'] + [res['invtransforms'][1]],
                            resB['fwdtransforms'] + [resB['invtransforms'][1]]):
                shutil.copyfile(b, a)
            report['chosen'] = 'B_refine'
            report['decision'] = 'A failed the checks; B improves the tissue overlap: B kept.'
        else:
            report['decision'] = 'A failed the checks; B does not improve enough: A kept (check this subject).'
        for f in resB['fwdtransforms'] + [resB['invtransforms'][1]]:
            if opi(f):
                os.remove(f)

    log(f'INFO: {NAME}: {report["decision"]}')
    with open(outprefix + 'adaptive_QC.json', 'w') as f:
        json.dump(report, f, indent=2)
    return res
