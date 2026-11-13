"""prepa_v2 copied verbatim from Exemples/Study_EDNiX/analysis/validation/test_surface_variants_library.py
(branch surface-variants-test, PR #84), the implementation validated by the surface tests. Reference for
test_surface_v2h.py only."""
import math
import os

import numpy as np

opj = os.path.join


def prepa_v2(E, Subname, Ref_file, aseg_file, prep_dir, const=False, hippo_wm=True, cerebrum_only=False):
    """V2 (const=False) or V3 (const=True) versions of desc-norm_T1w (brain), seg-wm (wm) and seg-filled
    (filled), written to prep_dir.
    V2: CSF 10-35 and GM 55-85 keep the within-class intensity variation of the template, WM 110, sigma 0.5 vox.
    V3: every class at one value (CSF 30, GM 75, WM + subcortex + hippocampus/amygdala 110), sigma 1 vox:
        the image edges are the aseg boundaries only.
    hippo_wm=False (V2h): hippocampus/amygdala left in the GM class and out of wm/filled, as in EDNiX (V0).
    cerebrum_only=True (V2c): cerebellum (7/8/46/47) and brainstem (16) set to 0 in brain.mgz. Otherwise the
        cerebellar WM sits at 110 (preFS.WM_edit) right under the ventral occipito-temporal cortex and the
        white/pial surfaces can snag on it."""
    import ants
    P = E['preFS']
    hippo_amyg = [17, 18, 53, 54] if hippo_wm else []
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
    if cerebrum_only:
        norm[np.isin(a, [7, 8, 46, 47, 16])] = 0.
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
    filled[np.isin(a, P.filled_edit[0] + ([17, 18] if hippo_wm else []))] = 127
    filled[np.isin(a, P.filled_edit[1] + ([53, 54] if hippo_wm else []))] = 255
    bil = np.isin(a, P.filled_edit[2])
    left, right = bil.copy(), bil.copy()
    left[midsection:, :, :] = False
    right[:midsection, :, :] = False
    filled[left] = 127 if dist < midsection else 255
    filled[right] = 255 if dist < midsection else 127
    ants.image_write(BRAIN.new_image_like(filled), opj(prep_dir, Subname + '_seg-filled_dseg.nii.gz'), ri=False)
