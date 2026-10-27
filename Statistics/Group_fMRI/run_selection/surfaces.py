"""
Surface renders of the group maps (all runs / kept / not kept) on one colour scale.

Moved VERBATIM from Exemples/Study_EDNiX/analysis/Fig5/fig5_paper.py (2026-10-27);
fig5_paper imports these names back.
"""
import os
import math
import re
import glob
import shutil
import itertools

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats as st



opj = os.path.join


def _oversampled(map_path, stat="ttest-stat_fisher_cmap_oversampled"):
    hits = sorted(glob.glob(opj(os.path.dirname(map_path), f"*{stat}.nii.gz")))
    return hits[0] if hits else map_path


def surface_maps_for_5a(all_units_index, map_index_final):
    """{label: {group: map path}} for the three groups of 5a / S5a: every run
    ('all'), kept by the final round ('kept_final'), the same number not kept
    ('not_kept_final')."""
    maps = map_paths_by_dataset(all_units_index, map_index_final)
    groups3 = ("all", "kept_final", "not_kept_final")
    maps = {lab: {g: p for g, p in gs.items() if g in groups3} for lab, gs in maps.items()}
    return {lab: gs for lab, gs in maps.items() if gs}


def map_paths_by_dataset(all_units_index, map_index_final):
    out = {}
    for idx, groups in ((all_units_index, ("all",)),
                        (map_index_final, ("kept_final", "not_kept_final"))):
        if idx is None or idx.empty:
            continue
        for r in idx[idx["group"].isin(groups) & idx["found"]].itertuples():
            out.setdefault(r.label, {})[r.group] = r.map_path
    return out


def render_common_scale(sri, out_dir, maps, species_of_label, pct=90.0):
    """Re-render every group map of a dataset (all runs, kept, not kept) with
    ONE fixed symmetric colour range: +/- the pct-th percentile of |value|
    over the non-zero voxels of all of that dataset's maps. Same wb_command
    call as SBA_Region_Profile_Index.rerun_balanced_by_run, except
    colorbar_mode='fixed' instead of a per-map 'abs_pct'.
    Returns {(label, group): png path} and {label: range}."""
    import nibabel as nib
    from Tools import niftotoWBsurface
    fig5 = sri.fig5
    pngs, ranges = {}, {}
    for lab, groups in maps.items():
        try:
            src = {g: _oversampled(p) for g, p in groups.items() if p}
            vals = []
            for p in src.values():
                d = np.asarray(nib.load(p).get_fdata()).ravel()
                vals.append(np.abs(d[np.isfinite(d) & (d != 0)]))
            vals = np.concatenate(vals) if vals else np.array([])
            if not vals.size:
                continue
            rng_ = float(np.percentile(vals, pct))
            ranges[lab] = rng_
            species = species_of_label[lab]
            scfg = fig5.build_species_cfg(species)
            atlas_ednix = opj(fig5.ATLAS_LIB, scfg["species_fragment"], "EDNiX")
            surface_cfg = dict(
                species_name=species, atlas_lib_path=fig5.ATLAS_LIB,
                species_fragment=scfg["species_fragment"], sing_wb=fig5.SING_WB,
                surface_type="inflated", palette_name="RY-BC-BL", pos_max=rng_,
                neg_max=-rng_, colorbar_mode="fixed",
                scene_template=opj(atlas_ednix, "surfaces", "Native_resol",
                                   f"{species}_thickness.scene"),
                scene_label="New Scene 1", display_zero=False)
            for k in ("projection", "delete_surf_file"):
                if scfg.get(k) is not None:
                    surface_cfg[k] = scfg[k]
            for g, p in src.items():
                d = opj(out_dir, f"{lab}__{g}")
                os.makedirs(d, exist_ok=True)
                _, out_pngs = niftotoWBsurface.nifti_to_surface_wb(
                    nifti_path=p, output_dir=d, map_prefix=f"{lab}__{g}", **surface_cfg)
                if out_pngs:
                    pngs[(lab, g)] = out_pngs[0]
            print(f"  [5a] {lab}: {sorted(src)} rendered with a common range +/-{rng_:.3g}")
        except Exception as e:
            print(f"  [5a] [WARN] {lab}: common-scale render failed ({e}) -- per-map renders "
                  f"kept for this dataset")
    return pngs, ranges
