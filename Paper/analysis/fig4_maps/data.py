"""Inputs of the Fig. 4 map analyses: run table, Fisher-z run maps, scores, bins,
regional means and correlation matrices of one finished selection (its 07_state,
read through paper_stats.Inputs(fmri_root=...).fmri()).

Group map of a set of runs = voxelwise MEAN of their Fisher-z run maps, i.e. the
profiled sub-brick of the pipeline's 3dttest++ maps (STAT_SUFFIX
'ttest-stat_fisher_cmap' = sub-brick 0, unthresholded). The run maps are the
inputs of the all-runs 3dttest++ (input_runs/ next to that map: same grid, same
group mask group_mask_rs.nii.gz), so every group is profiled exactly as the figure
maps are: same function (profiles.compute_region_profile), atlas volume and atlas
level. The mean map is handed over in memory; no NIfTI, no AFNI, no ClustSim.
"""
import contextlib
import io
import os
import re

import nibabel as nib
import numpy as np
import pandas as pd

# ── run keys ────────────────────────────────────────────────────────────────────────────
_NULL = {"", "none", "nan", "na", "all", "null"}


def _tok(x, prefix):
    """Normalised BIDS token: prefix removed, lower case, leading zeros of a number
    removed ('01' = '1'), absent values ('', None, 'NA', 'ALL') -> None."""
    x = str(x).strip()
    x = x[len(prefix):] if x.lower().startswith(prefix) else x
    if x.lower() in _NULL:
        return None
    x = x.lower()
    return (x.lstrip("0") or "0") if x.isdigit() else x


def run_key(sub, ses, run):
    return _tok(sub, "sub-"), _tok(ses, "ses-"), _tok(run, "run-")


def pair_runs(vec_keys, map_keys):
    """Match every run vector to one run map. vec_keys {row: (sub, ses, run)}, map_keys
    iterable of (sub, ses, run), both normalised by run_key. Exact match first;
    otherwise the unique map of the same subject whose session and run agree wherever
    both are given (a token absent on one side, e.g. a dataset without ses- folders,
    does not block the match). A map is never used twice.
    Returns ({row: map key}, [rows unmatched], [map keys unmatched], [rows ambiguous])."""
    maps = list(dict.fromkeys(map_keys))
    free = set(maps)
    out, amb = {}, []
    for row, k in vec_keys.items():
        if k in free:
            out[row] = k
            free.discard(k)
    for row, (s, se, r) in vec_keys.items():
        if row in out:
            continue
        cand = [m for m in maps if m in free and m[0] == s
                and (se is None or m[1] is None or m[1] == se)
                and (r is None or m[2] is None or m[2] == r)]
        if len(cand) == 1:
            out[row] = cand[0]
            free.discard(cand[0])
        elif len(cand) > 1:
            amb.append(row)
    miss = [row for row in vec_keys if row not in out]
    return out, miss, [m for m in maps if m in free], amb


def input_runs(map_path):
    """{run_key: path} of the Fisher-z run maps that went into the 3dttest++ of map_path
    (files '<NN>_sub-<s>_ses-<ses>_run-<r>.nii.gz' of input_runs/)."""
    d = os.path.join(os.path.dirname(map_path), "input_runs")
    if not os.path.isdir(d):
        raise FileNotFoundError(f"no input_runs/ next to {map_path}")
    out = {}
    for f in sorted(os.listdir(d)):
        m = re.match(r"^\d+_sub-(.+?)_ses-(.+?)_run-(.+?)\.nii\.gz$", f)
        if m:
            out[run_key(*m.groups())] = os.path.join(d, f)
    return out


# ── run table ───────────────────────────────────────────────────────────────────────────
class RunTable:
    """Every analysed run of the selection: its P80/P20 vector over c['common']
    (normalised as in the selection, ZSCORE), dataset index, selected flag (final
    round) and normalised key. Requires UNIT = 'run'."""

    def __init__(self, c):
        from Paper.analysis.fmri_fingerprints import profile_stats as fa
        cfg = c["cfg"]
        value_col, norm = cfg["VALUE_COL"], cfg["ZSCORE"]
        sol, com = c["species_of_label"], c["common"]
        self.labels = sorted(c["profiles_by_label"])
        vecs, ds, selected, keys = [], [], [], []
        for i, lab in enumerate(self.labels):
            kept = set(c["selected_units"].get(lab, {}))
            mem = {tuple(str(x) for x in k): v for k, v in c["members_by_label"].get(lab, {}).items()}
            for key, prof in c["profiles_by_label"][lab].items():
                v = fa._ordered_common_values(prof, com, value_col)
                if not np.isfinite(v).any():
                    continue
                runs = mem.get(tuple(str(x) for x in key), [key])
                if len(runs) != 1:
                    raise RuntimeError(f"{lab}: unit {key} has {len(runs)} runs -- the Fig. 4 "
                                       f"map analyses need UNIT='run'")
                r = runs[0]
                vecs.append(fa._normalize_values(v, norm))
                ds.append(i)
                selected.append(key in kept)
                keys.append(run_key(r[0], r[1], r[2]))
        if not vecs:
            raise RuntimeError("no run profile in the selection state")
        self.vecs, self.ds = np.vstack(vecs), np.array(ds)
        self.selected, self.keys = np.array(selected), keys
        self.sp_of_ds = np.array([sol.get(l, l) for l in self.labels])
        self.sets = np.unique(self.ds)

    def rows(self, d):
        return np.flatnonzero(self.ds == d)


def scores(c, rt):
    """r of every run vector with its dataset's round-final selection criterion
    (c['final_criterion'], the target the runs were selected against) on the regions
    of the selection mask c['region_mask']."""
    from Paper.analysis.fmri_fingerprints import profile_stats as fa
    com, vc = c["common"], c["cfg"]["VALUE_COL"]
    mask = np.asarray(c.get("region_mask", np.ones(len(com), bool)), bool)
    if mask.shape != (len(com),):
        mask = np.ones(len(com), bool)
    out = np.full(len(rt.ds), np.nan)
    for d in rt.sets:
        crit = c["final_criterion"].get(rt.labels[d])
        if crit is None:
            raise RuntimeError(f"{rt.labels[d]}: no final selection criterion in the state")
        t = fa._ordered_common_values(crit, com, vc)
        for i in rt.rows(d):
            ok = mask & np.isfinite(rt.vecs[i]) & np.isfinite(t)
            if ok.sum() >= 3:
                out[i] = np.corrcoef(rt.vecs[i][ok], t[ok])[0, 1]
    return out, int(mask.sum())


def regions_of_common(common):
    """Homologous regions of c['common'] ('region|q80' -> 'region'), in order."""
    return list(dict.fromkeys(str(x).rpartition("|")[0].lower().strip() or str(x).lower().strip()
                              for x in common))


# ── Fisher-z run maps of one dataset ────────────────────────────────────────────────────
class DatasetMaps:
    """Fisher-z run maps of one dataset (input_runs/ of its all-runs 3dttest++), kept in
    memory on the voxels of the group mask that carry an atlas label."""

    def __init__(self, c, label, map_all):
        from Statistics.Group_fMRI.run_selection import profiles as sri
        from Statistics.Group_fMRI.run_selection import units as fu
        from Paper.analysis.fmri_fingerprints import profile_stats as fa
        self.sri, self.fa = sri, fa
        self.label, self.map_all = label, map_all
        self.level = c["cfg"]["ATLAS_LEVEL"]
        self.value_col, self.common = c["cfg"]["VALUE_COL"], c["common"]
        cb = {b["label"]: b for b in c["bids_configs"]}[label]
        run_keys = [rk for rks in c["members_by_label"].get(label, {}).values() for rk in rks]
        self.atlas = (fu.find_matching_atlas_volume(sri, cb, run_keys, map_all, verbose=False,
                                                    label=label)
                      or sri.pick_representative_atlas_volume(cb["bids_dir"]))
        self.atlas_img = nib.load(self.atlas)
        self.atlas_img.get_fdata()                       # cached by nibabel for later calls
        self.label_df = sri.parse_itk_snap_label_file(cb["label_path"])
        ref = nib.load(map_all)
        self.shape, self.affine = ref.shape[:3], ref.affine
        mpath = os.path.join(os.path.dirname(map_all), "group_mask_rs.nii.gz")
        if os.path.isfile(mpath):
            mask = np.squeeze(nib.load(mpath).get_fdata()) > 0
            self.mask_src = "group_mask_rs.nii.gz"
        else:                                            # 3dttest++ writes 0 outside -mask
            mask = np.squeeze(ref.get_fdata()) != 0
            self.mask_src = "non-zero voxels of the all-runs map"
        atl = np.squeeze(self.atlas_img.get_fdata())
        atl_any = (atl > 0).any(axis=3) if atl.ndim == 4 else atl > 0
        if mask.shape != self.shape or atl_any.shape != self.shape:
            raise ValueError(f"{label}: grid mismatch map {self.shape} / mask {mask.shape} / "
                             f"atlas {atl_any.shape}")
        self.vox = np.flatnonzero((mask & atl_any).ravel())
        self.paths = input_runs(map_all)
        self.rows = {}

    def load(self, keys):
        for k in keys:
            if k in self.rows:
                continue
            img = nib.load(self.paths[k])
            if img.shape[:3] != self.shape:
                raise ValueError(f"{self.label}: {self.paths[k]} grid {img.shape} != {self.shape}")
            self.rows[k] = np.squeeze(img.get_fdata(dtype=np.float32)).ravel()[self.vox]

    def mean_image(self, keys):
        keys = list(keys)
        self.load(keys)
        m = np.mean(np.vstack([self.rows[k] for k in keys]), axis=0, dtype=np.float64)
        vol = np.zeros(int(np.prod(self.shape)), np.float32)
        vol[self.vox] = m
        return nib.Nifti1Image(vol.reshape(self.shape), self.affine)

    def profile(self, keys):
        """Unnormalised P80/P20 vector over c['common'] of the mean map of `keys`."""
        with contextlib.redirect_stdout(io.StringIO()):     # one [INFO] line per call
            prof = self.sri.compute_region_profile(self.mean_image(keys), self.atlas_img,
                                                   self.label_df, atlas_level=self.level)
        return self.fa._ordered_common_values(prof, self.common, self.value_col)

    def region_means(self, regions):
        """{run key: mean Fisher z of each region (left and right pooled, group-mask
        voxels)}. Linear, so the regional means of a mean map are the mean of these."""
        atl = np.squeeze(self.atlas_img.get_fdata())
        if atl.ndim == 4:
            atl = atl[..., (self.level or atl.shape[-1]) - 1]
        lab = atl.ravel()[self.vox].astype(int)
        base = dict(zip(self.label_df["label_id"].astype(int),
                        self.label_df["base_name"].astype(str).str.lower().str.strip()))
        reg = np.array([regions.index(base[x]) if base.get(x) in regions else -1 for x in lab])
        missing = [r for i, r in enumerate(regions) if not np.any(reg == i)]
        if missing:
            raise RuntimeError(f"{self.label}: no voxel for region(s) {missing} at atlas level "
                               f"{self.level}")
        W = np.zeros((len(self.vox), len(regions)), np.float32)
        W[reg >= 0, reg[reg >= 0]] = 1
        W /= W.sum(axis=0, keepdims=True)
        keys = list(self.paths)
        out = {}
        for j in range(0, len(keys), 20):
            ks = keys[j:j + 20]
            self.load(ks)
            X = np.vstack([self.rows[k] for k in ks])
            for k, v in zip(ks, np.nan_to_num(X) @ W):
                out[k] = v.astype(float)
        return out


# ── correlation matrices ────────────────────────────────────────────────────────────────
def matrix_runs(c, cb, regions, log=print):
    """{run_key: {edge: mean Fisher z}} of one dataset, from the per-run region x region
    correlation matrices of EDNiX step 10
    (<bids>/sub-*/[ses-*/]func/acpc-func/Stats/Correl_matrix/EDNIxCSCLR/correlation/
    EDNIxCSCLR_<ATLAS_LEVEL - 1>_run_<i>_correlation_matrix.csv), read with the reader
    of the matrix selection (sensitivity_analysis/fmri_selection/matrices.py). Nodes =
    `regions`; edges L/R pooled ('a <> b |intra', '|inter', 'a |homotopic'), value =
    mean arctanh(r) of the hemispheric versions."""
    from Statistics.Group_fMRI.run_selection import profiles as sri
    from Paper.analysis.sensitivity_analysis.fmri_selection import matrices as mx
    cfg = c["cfg"]
    mp = mx.MatrixProfiles(sri, c["bids_configs"], node_regions=regions, value_col="median_r")
    label_df = sri.parse_itk_snap_label_file(cb["label_path"])
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        profs, _ = mp.compute_all_run_profiles(cb["species"], cb["bids_dir"], label_df,
                                               cfg.get("CLASSIFICATION_PATH"),
                                               category=cfg.get("CATEGORY", "All"),
                                               atlas_level=cfg["ATLAS_LEVEL"])
    for ln in buf.getvalue().splitlines():
        if "WARN" in ln:
            log("    " + ln.strip())
    out = {}
    for key, p in profs.items():
        v = pd.to_numeric(p["median_r"], errors="coerce").clip(-0.999999, 0.999999)
        zz = pd.Series(np.arctanh(v.to_numpy()), index=p["base_name"].astype(str))
        out[run_key(key[0], key[1], key[2])] = zz.groupby(level=0).mean().to_dict()
    return out
