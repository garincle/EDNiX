"""
energetic_cost.py -- two cross-species tests of the "energetic cost of signalling" hypothesis
(Castrillon et al., Sci. Adv. 9:eadi7632, 2023; Tomasi et al., PNAS 2013).

Castrillon et al. regress the cerebral glucose metabolism of the cortex on its degree of functional connectivity
(dFC: sum of the positive correlations, p < 0.001, of the BOLD series of a voxel with every other cortical voxel;
r ~ 0.4-0.5 voxel-wise within subject, P_smash < 0.024). The residual of that regression is the "energetic cost" of
signalling; it is higher in the cortex that expanded most between chimpanzee and human (r = 0.30, P_smash = 0.047).

  1. coupling_analysis : is the FDG ~ degree coupling present in each species?  (voxel-wise, spatial-autocorrelation-aware null)
  2. expansion_analysis: is the human-macaque difference in FDG larger where the human region is relatively larger?
                         (region-wise; the relative expansion is measured on our own morphometry, no external data)

Limits to state in the paper: bratio is a RATIO to the whole-brain mean (not CMRglc); the macaque PET is n = 2 under
propofol; the surrogate maps below preserve the amplitude spectrum of the map (phase randomisation), which is an
approximation of BrainSMASH (Burt et al., 2020) with edge effects at the mask border.
"""
import numpy as np
import pandas as pd
import nibabel as nib
from scipy import ndimage, stats


# ===================================================================================================== degree
def r_threshold(n_t, p=0.001):
    """Smallest |r| with two-sided P < p for n_t time points (t test, n_t - 2 d.f.)."""
    t = stats.t.isf(p / 2.0, n_t - 2)
    return float(t / np.sqrt(n_t - 2 + t * t))


def degree_centrality(ts, p=0.001, chunk=2000):
    """ts: (n_voxels, n_t). Positive degree of every voxel: sum of the correlations r_ij > r_threshold(n_t, p) with the
    other voxels (Castrillon 2023, Tomasi 2013: 3dDegreeCentrality). Computed by blocks (no full voxel x voxel matrix)."""
    x = np.asarray(ts, np.float32)
    x = x - x.mean(1, keepdims=True)
    sd = np.sqrt((x ** 2).sum(1, keepdims=True))
    ok = sd[:, 0] > 0
    x[ok] /= sd[ok]
    x[~ok] = 0
    thr = np.float32(r_threshold(x.shape[1], p))
    deg = np.zeros(x.shape[0], np.float64)
    for a in range(0, x.shape[0], chunk):
        r = x[a:a + chunk] @ x.T
        r[r <= thr] = 0.0
        idx = np.arange(a, min(a + chunk, x.shape[0]))
        r[np.arange(len(idx)), idx] = 0.0            # no self-correlation
        deg[a:a + chunk] = r.sum(1, dtype=np.float64)
    deg[~ok] = np.nan
    return deg


def degree_map(bold_path, mask=None, p=0.001):
    """Degree map of one run: 3D array of the z-scored (within the mask) positive degree, NaN outside the mask. mask:
    boolean 3D array (cortex) in the grid of the run; None = every voxel with a non-constant series."""
    img = nib.load(bold_path)
    data = np.asarray(img.dataobj, np.float32)
    if data.ndim != 4:
        raise ValueError(f'{bold_path}: expected a 4D image, got {data.shape}')
    m = np.std(data, axis=3) > 0
    if mask is not None:
        if mask.shape != m.shape:
            raise ValueError(f'mask {mask.shape} and {bold_path} {m.shape} are not on the same grid')
        m &= mask
    deg = degree_centrality(data[m], p=p)
    out = np.full(m.shape, np.nan, np.float32)
    z = (deg - np.nanmean(deg)) / np.nanstd(deg)
    out[m] = z
    return out, img.affine


def group_degree_map(maps, min_fraction=0.8):
    """Mean of the per-run z maps; voxels valid in fewer than min_fraction of the runs are NaN."""
    st = np.stack(maps).astype(np.float32)
    n = np.isfinite(st).sum(0)
    with np.errstate(all='ignore'):
        mean = np.nanmean(st, axis=0)
    mean[n < min_fraction * len(maps)] = np.nan
    return mean


# =============================================================================================== resampling
def resample_to_grid(data, src_affine, ref_shape, ref_affine, order=1):
    """Resample the 3D array `data` (grid src_affine) onto the grid (ref_shape, ref_affine); NaN outside the source.
    Both affines map voxel indices to the same world (template) space."""
    ijk = np.indices(ref_shape).reshape(3, -1).astype(np.float64)
    world = ref_affine[:3, :3] @ ijk + ref_affine[:3, 3:4]
    src = np.linalg.inv(src_affine)
    sij = src[:3, :3] @ world + src[:3, 3:4]
    filled = np.where(np.isfinite(data), data, 0.0)
    valid = np.isfinite(data).astype(np.float64)
    v = ndimage.map_coordinates(filled, sij, order=order, mode='constant', cval=0.0)
    w = ndimage.map_coordinates(valid, sij, order=order, mode='constant', cval=0.0)
    with np.errstate(all='ignore'):
        out = np.where(w > 0.5, v / w, np.nan)
    return out.reshape(ref_shape).astype(np.float32)


def smooth_masked(data, fwhm_mm, affine):
    """Gaussian smoothing restricted to the finite voxels (normalised convolution): no leakage of the NaN background."""
    vox = np.sqrt((affine[:3, :3] ** 2).sum(0))
    sigma = fwhm_mm / 2.3548 / vox
    ok = np.isfinite(data)
    num = ndimage.gaussian_filter(np.where(ok, data, 0.0).astype(np.float64), sigma)
    den = ndimage.gaussian_filter(ok.astype(np.float64), sigma)
    with np.errstate(all='ignore'):
        out = np.where(ok & (den > 1e-6), num / den, np.nan)
    return out.astype(np.float32)


# ============================================================================================ spatial null
def phase_surrogates(data, mask, n, seed=0):
    """n surrogate maps with the amplitude spectrum of `data` (inside `mask`, demeaned) and random phases, masked
    again: they keep the spatial autocorrelation of the map but not its alignment with another map."""
    rng = np.random.default_rng(seed)
    x = np.where(mask, data - np.nanmean(data[mask]), 0.0)
    amp = np.abs(np.fft.rfftn(x))
    for _ in range(n):
        ph = np.angle(np.fft.rfftn(rng.standard_normal(x.shape)))
        s = np.fft.irfftn(amp * np.exp(1j * ph), s=x.shape)
        yield s[mask]


def coupling_analysis(deg, pet, mask, n_surr=1000, seed=0, method='spearman'):
    """Correlation over the voxels of `mask` between the degree map and the PET map (same grid), with
    P (two-sided) from phase-randomised surrogates of the degree map and the residual map (PET - linear fit on degree,
    "energetic cost"). Returns (dict of statistics, residual map)."""
    corr = (lambda a, b: stats.spearmanr(a, b).statistic) if method == 'spearman' else \
           (lambda a, b: float(np.corrcoef(a, b)[0, 1]))
    m = mask & np.isfinite(deg) & np.isfinite(pet)
    a, b = deg[m], pet[m]
    r = corr(a, b)
    null = np.array([corr(s, b) for s in phase_surrogates(np.where(m, deg, np.nan), m, n_surr, seed)])
    p = (1 + np.sum(np.abs(null) >= abs(r))) / (1 + len(null))
    slope, icpt = np.polyfit(a, b, 1)
    res = np.full(deg.shape, np.nan, np.float32)
    res[m] = b - (slope * a + icpt)
    return dict(n_voxels=int(m.sum()), r=float(r), p_surrogate=float(p), slope=float(slope),
                null_sd=float(null.std()), n_surrogates=int(n_surr), method=method), res


def regional_means(data, labels, names):
    """Mean of `data` per label of the integer 3D array `labels` (same grid); names {label id: region}; regions with
    several labels (L + R) are pooled."""
    out = {}
    for i, name in names.items():
        v = data[(labels == i) & np.isfinite(data)]
        if v.size:
            out.setdefault(name, []).append(v)
    return pd.Series({k: np.concatenate(v).mean() for k, v in out.items()}, name='value')


# ============================================================================================== expansion
def expansion_table(volume, fdg, species=('Human', 'Macaque'), regions=None):
    """volume, fdg: long tables (species, subject, region, value) after name harmonisation (volume = % of brain volume,
    fdg = bratio). One row per region present in both modalities and both species:
      log2_expansion = log2(volume_species[0] / volume_species[1])      (relative expansion of the region)
      log2_fdg_ratio = log2(bratio_species[0] / bratio_species[1])      (species difference of the PET ratio)."""
    a, b = species

    def mean(t):
        t = t[t['species'].isin(species)]
        s = t.groupby(['species', 'subject', 'region'])['value'].mean().reset_index()
        return s.groupby(['species', 'region'])['value'].mean().unstack('species')

    v, f = mean(volume), mean(fdg)
    df = pd.DataFrame({'volume_a': v[a], 'volume_b': v[b], 'fdg_a': f[a], 'fdg_b': f[b]}).dropna()
    if regions is not None:
        df = df.loc[[r for r in regions if r in df.index]]
    df['log2_expansion'] = np.log2(df['volume_a'] / df['volume_b'])
    df['log2_fdg_ratio'] = np.log2(df['fdg_a'] / df['fdg_b'])
    return df.rename(columns={'volume_a': f'volume_{a}', 'volume_b': f'volume_{b}',
                              'fdg_a': f'bratio_{a}', 'fdg_b': f'bratio_{b}'})


def expansion_analysis(tab, cortical=None, n=10000, seed=0):
    """Spearman r (permutation P, bootstrap CI over regions) between the relative expansion and the species difference
    in FDG, over all regions of `tab` and, if `cortical` (list of region names) is given, over the cortical regions
    only (the propofol effect on the macaque is mostly cortical, so the two are reported). Holm over the rows."""
    from Exemples.Study_EDNiX.analysis.library import multimodal_profiles as MP
    rng = np.random.default_rng(seed)
    rows = []
    sets = [('all regions', list(tab.index))]
    if cortical is not None:
        sets.append(('cortical regions', [r for r in tab.index if r in set(cortical)]))
    for name, regs in sets:
        t = tab.loc[regs]
        r, p, lo, hi = MP._perm_boot(t['log2_expansion'].values, t['log2_fdg_ratio'].values, rng, n)
        rows.append(dict(regions=name, n_regions=len(t), r=r, p=p, ci_low=lo, ci_high=hi))
    out = pd.DataFrame(rows)
    out['p_holm'] = MP._holm(np.nan_to_num(out['p'].values, nan=1.0))
    return out
