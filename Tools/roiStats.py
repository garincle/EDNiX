"""
Regional statistics of an image over a label image (atlas), in numpy.

Suggested by the PET engineer as a wrapper of AFNI 3dROIstats; rewritten in numpy: 3dROIstats sorts
the labels and drops the absent ones, which the wrapper had to undo column by column (with a bug:
`img.unique()` called on a file name), whereas here every value is attached to its label explicitly.
"""
import numpy as np
import pandas as pd
import nibabel as nib


def _labels_volume(atlas_file, level=0, ref=None):
    img = nib.load(atlas_file)
    lab = np.asarray(img.dataobj)
    if lab.ndim > 3:
        if level >= lab.shape[3]:
            raise ValueError(f'{atlas_file} has {lab.shape[3]} levels; level {level} requested')
        lab = lab[..., level]
    lab = np.rint(lab).astype(np.int64)
    if ref is not None:
        r = nib.load(ref)
        if r.shape[:3] != lab.shape or not np.allclose(r.affine, img.affine, atol=1e-3):
            raise ValueError(f'{atlas_file} is not on the grid of {ref}')
    return lab


def regional_table(labels_vol, label_df, voxel_mm3, labels=None):
    """label, region, n_voxels, volume_mm3 of every label of label_df (columns label_id, region_name),
    or of every label present when label_df is None; absent labels get n_voxels = 0."""
    present, counts = np.unique(labels_vol[labels_vol > 0], return_counts=True)
    count = dict(zip(present.tolist(), counts.tolist()))
    if labels is None:
        labels = label_df['label_id'].tolist() if label_df is not None else present.tolist()
    names = dict(zip(label_df['label_id'], label_df['region_name'])) if label_df is not None else {}
    return pd.DataFrame({'label': labels,
                         'region': [names.get(l, f'label_{l}') for l in labels],
                         'n_voxels': [count.get(l, 0) for l in labels],
                         'volume_mm3': [count.get(l, 0) * voxel_mm3 for l in labels]})


def stats(data3d, labels_vol, label_df, voxel_mm3, extra_mask=None):
    """Mean, SD, median of data3d in each region (NaN for empty regions)."""
    tab = regional_table(labels_vol, label_df, voxel_mm3)
    valid = np.isfinite(data3d)
    if extra_mask is not None:
        valid &= extra_mask
    means, sds, meds, ns = [], [], [], []
    for l in tab['label']:
        v = data3d[(labels_vol == l) & valid]
        ns.append(v.size)
        means.append(float(v.mean()) if v.size else np.nan)
        sds.append(float(v.std(ddof=1)) if v.size > 1 else np.nan)
        meds.append(float(np.median(v)) if v.size else np.nan)
    tab['n_valid'] = ns
    tab['mean'] = means
    tab['sd'] = sds
    tab['median'] = meds
    return tab


def tacs(data4d, labels_vol, label_df, voxel_mm3):
    """Time-activity curves: DataFrame frames x labels (present labels only), with attrs names / volumes."""
    tab = regional_table(labels_vol, label_df, voxel_mm3)
    tab = tab[tab['n_voxels'] > 0].reset_index(drop=True)
    flat = labels_vol.reshape(-1)
    d = data4d.reshape(-1, data4d.shape[-1])
    cols = {}
    with np.errstate(invalid='ignore'):
        for l in tab['label']:
            # NaN voxels (outside the field of view after resampling) ignored, as in stats()
            cols[int(l)] = np.nanmean(d[flat == l], axis=0)
    df = pd.DataFrame(cols)
    df.attrs['names'] = tab['region'].tolist()
    df.attrs['volumes_mm3'] = tab['volume_mm3'].tolist()
    return df


def get(img, atlas_file, label_df, level=0, ref_check=True):
    """Regional statistics of a 3D image file over level `level` of an atlas file on the same grid."""
    im = nib.load(img)
    data = np.asarray(im.dataobj, dtype=np.float32)
    if data.ndim == 4 and data.shape[3] == 1:
        data = data[..., 0]
    if data.ndim > 3:
        raise ValueError(f'{img} is 4D: use tacs()')
    lab = _labels_volume(atlas_file, level, img if ref_check else None)
    return stats(data, lab, label_df, float(np.prod(im.header.get_zooms()[:3])))
