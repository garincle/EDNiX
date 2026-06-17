"""
ednix_pgls_fc_scatter.py  ?  v4
================================

Two outputs
-----------
A) Morphometry PGLS  (surface / thickness / volume)
   ? All pairwise region PGLS at atlas level 3  ?  run_pgls_all_region_pairs
   ? Per-region allometric scaling vs body weight, brain weight, total neurons
     ?  plot_scaling_with_neurons  (4 panels: body / brain / neurons / internal)

B) FC normalized fingerprint  (after your 2020 radar/similarity script)
   ? Per species: mean connectivity matrix across all subjects
   ? Fingerprint = for each ROI k, take its row, remove self-connection,
     normalize min-max  (exactly as in your old script)
   ? All-pairs cosine similarity + permutation test (100k shuffles, Weibull fit)
   ? Output: heatmap of similarity (species pairs × ROIs) + p-value heatmap

Usage
-----
    python ednix_pgls_fc_scatter.py

Prerequisites
-------------
  - Morphometry data via collect_multi_species
  - Per-species mean connectivity matrices via _collect_corr_per_bids
    (called internally; reads BIDS directories directly)
"""

import os
import sys
import math
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as tck
import seaborn as sns
from scipy import stats as _stats
from PIL import Image

# ?? project imports ???????????????????????????????????????????????????????????
sys.path.insert(0, "/home/cgarin/PycharmProjects/EDNiX")

from Exemples.Study_EDNiX.analysis.ednix_bids_tools import (
    PAPER_RC, collect_multi_species,
    get_atlas_label_path,
    find_species_path,
    extract_corr_matrix_paths,
    load_corr_matrix,
)
from Statistics.Evo.PGLS import (
    load_newick_vcv,
    run_pgls_all_region_pairs,
    plot_pgls_region_pair,
    pgls_on_ax,
    _species_colors,
    _phylo_sort,
)

opj = os.path.join

# ?????????????????????????????????????????????????????????????????????????????
# CONFIGURATION
# ?????????????????????????????????????????????????????????????????????????????

NEWICK_FILE    = "/scratch2/EDNiX/Liste des espèces.nwk"
ATLAS_LIB      = "/home/cgarin/PycharmProjects/EDNiX/Atlases_library"
ATLAS_NAME     = "EDNIxCSC"
ATLAS_LEVEL    = 3
HERCULANO_XLSX = "/scratch2/EDNiX/Herculano-Houzel_mammalian_brain_dataset_complete.xlsx"

OUT_DIR          = "/scratch2/EDNiX/results/multispecies_analysis"
FIG_DIR          = opj(OUT_DIR, "figures")
SCALING_DIR      = opj(FIG_DIR, "cross_species", "brain_scaling_pgls")
PAIRS_DIR        = opj(FIG_DIR, "pgls", "pairs_lvl3")
FINGERPRINT_DIR  = opj(FIG_DIR, "cross_species", "fc_fingerprint")

N_PERMUTATIONS = 100   # permutation test (as in your original script)

SPECIES_BIDS = {
    "Rat":        "/scratch2/EDNiX/Rat/BIDS_Grandjean",
    "Mouse":      "/scratch2/EDNiX/Mouse/BIDS_Grandjean",
    "Dog":        "/scratch2/EDNiX/Dog/BIDS_k9",
    "Marmoset":   "/scratch2/EDNiX/Marmoset/BIDS_NIH_MBM",
    "Mouselemur": "/scratch2/EDNiX/Mouselemur/BIDS_Garin",
}
SPECIES_MULTI_BIDS = {
    "Macaque": ["/scratch2/EDNiX/Macaque/BIDS_BenHamed",
                "/scratch2/EDNiX/Macaque/BIDS_Cdt_Garin"],
    "Human":   ["/scratch2/EDNiX/Human/ds004513-download"],
}

REGIONS_OF_INTEREST = ['Auditory cortex (Superior temporal )', 'Insula and others in lateral sulcus', 'Middle Temporal, Inferior temporal , Temporal pole  (MIPT)', 'Motor and premotor', 'Olfactory cortex', 'Orbital PFC (oPFC)', 'Orbital frontal cortex (oFC)', 'Periarchicortex', 'Posterior medial cortex (PMC)', 'Posterior parietal cortex', 'Prefrontal cortex (PFC)', 'Somatosensory cortex', 'Ventral areas of the temporal lobe (vent Temp)', 'Visual pre and extra striate cortex', 'Visual striate cortex']


# ?????????????????????????????????????????????????????????????????????????????
# HERCULANO-HOUZEL LOADER
# ?????????????????????????????????????????????????????????????????????????????

_HH_SEARCH = {
    "Mouse":      "musculus",
    "Rat":        "norvegicus",
    "Mouselemur": "Microcebus",
    "Marmoset":   "jacchus",
    "Macaque":    "mulatta",
    "Human":      "sapiens",
}

def _parse_num(x):
    if pd.isna(x): return np.nan
    s = str(x).split("±")[0].split("+-")[0].replace(",", "").strip()
    try:    return float(s)
    except: return np.nan

def load_herculano(xlsx_path):
    out = {}
    try:
        df5 = pd.read_excel(xlsx_path, sheet_name="T5 \u2013 Whole Brain",
                            header=1, skiprows=[0])
        df5.columns = ["Species","Order","n","Body_mass_g","Brain_mass_g",
                       "Neurons","Other_cells","pct_Neurons","Source"]
        for sp, sub in _HH_SEARCH.items():
            row = df5[df5["Species"].str.contains(sub, case=False, na=False)]
            if not row.empty:
                r = row.iloc[0]
                out[sp] = dict(
                    body_g  = _parse_num(r["Body_mass_g"]),
                    brain_g = _parse_num(r["Brain_mass_g"]),
                    neurons = _parse_num(str(r["Neurons"]).replace(",", "")),
                )
    except Exception as e:
        warnings.warn(f"load_herculano T5: {e}")

    try:
        df_c = pd.read_excel(xlsx_path, sheet_name="Carnivora (2017)", header=None)
        for _, row in df_c.iterrows():
            if "Canis" in str(row[0]):
                bkg = _parse_num(row[2])
                n_s = str(row[7]).replace("\u00d7" + "10\u2076","e6").replace("x10^6","e6").replace(",","")
                try:    neurons = float(n_s)
                except: neurons = np.nan
                out["Dog"] = dict(body_g=bkg*1000 if np.isfinite(bkg) else np.nan,
                                  brain_g=_parse_num(row[3]), neurons=neurons)
                break
    except Exception as e:
        warnings.warn(f"load_herculano Carnivora: {e}")

    if "Bat" not in out:
        out["Bat"] = dict(body_g=20.0, brain_g=0.80, neurons=6.0e7)

    for sp, v in sorted(out.items()):
        print(f"  [HH] {sp:12s}  body={v['body_g']}g  "
              f"brain={v['brain_g']}g  neurons={v['neurons']:.2e}")
    return out


# ?????????????????????????????????????????????????????????????????????????????
# A) MORPHOMETRY  ? allometric scaling with neurons panel
# ?????????????????????????????????????????????????????????????????????????????

def plot_scaling_with_neurons(df_morph, metric_col, metric_label,
                               output_path, vcv_df, y_region, hh,
                               atlas_level=3, hemisphere="bilateral"):
    """
    4-panel PGLS allometric scaling for one region:
      A) body weight  B) brain weight  C) total neurons  D) Somatosensory (internal)
    """
    df = df_morph.copy()
    if "atlas_level" in df.columns:
        df = df[df["atlas_level"] == atlas_level]
    if "hemisphere" in df.columns and hemisphere != "all":
        h_df = df[df["hemisphere"] == hemisphere]
        df   = h_df if not h_df.empty else df[df["hemisphere"] == "bilateral"]

    df_y = df[df["region"] == y_region].groupby("species")[metric_col].mean()
    if df_y.empty:
        warnings.warn(f"plot_scaling_with_neurons: no data for {y_region}")
        return None

    sp_list   = list(df_y.index)
    y_vals    = df_y.values
    sp_colors = _species_colors(sp_list)

    body_arr  = np.array([hh[s]["body_g"]  if s in hh else np.nan for s in sp_list])
    brain_arr = np.array([hh[s]["brain_g"] if s in hh else np.nan for s in sp_list])
    neur_arr  = np.array([hh[s]["neurons"] if s in hh else np.nan for s in sp_list])

    x_region = "Somatosensory cortex"
    df_x = df[df["region"] == x_region].groupby("species")[metric_col].mean()
    xreg = np.array([float(df_x.loc[s]) if s in df_x.index else np.nan
                     for s in sp_list])

    panels = [(body_arr,  "Body weight (g)",             True),
              (brain_arr, "Brain weight (g)",             True),
              (neur_arr,  "Total neurons",                True)]
    if np.isfinite(xreg).sum() >= 3:
        panels.append((xreg, f"Somatosensory ({metric_label})", True))

    panels = [(xv, xl, lx) for xv, xl, lx in panels
              if np.isfinite(xv).sum() >= 3]
    if not panels:
        return None

    n_panels = len(panels)
    with plt.rc_context(PAPER_RC):
        fig, axes = plt.subplots(1, n_panels, figsize=(n_panels * 4.5, 4.5))
        if n_panels == 1: axes = [axes]

        for pi, (ax, (xv, xlabel, log_x)) in enumerate(zip(axes, panels)):
            for sp, x, y in zip(sp_list, xv, y_vals):
                if not (np.isfinite(x) and np.isfinite(y)): continue
                ax.scatter(x, y, color=sp_colors[sp], s=55, zorder=5,
                           edgecolors="k", linewidths=0.5)
                ax.annotate(sp, (x, y), textcoords="offset points",
                            xytext=(4, 3), fontsize=6, color=sp_colors[sp])
            pgls_on_ax(ax, xv, y_vals, vcv_df, sp_list,
                       log_x=log_x, log_y=True)
            ax.set_xlabel(xlabel, fontsize=9)
            ax.set_ylabel(f"{y_region}\n({metric_label})" if pi == 0 else "", fontsize=8)
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(tck.LogFormatterSciNotation())
            ax.text(-0.14, 1.04, chr(65+pi), transform=ax.transAxes,
                    fontsize=13, fontweight="bold", va="bottom", ha="left")

        handles = [mpatches.Patch(facecolor=sp_colors[sp], label=sp)
                   for sp in _phylo_sort(sp_list)]
        plt.tight_layout()
        fig.legend(handles=handles, loc="lower center",
                   bbox_to_anchor=(0.5, -0.06),
                   ncol=min(len(handles), 7),
                   fontsize=8, frameon=False, title="Species")
        fig.subplots_adjust(bottom=0.18)
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight", dpi=200)
        plt.close(fig)
        print(f"  [plot] {output_path}")
    return output_path


def run_morphometry_scaling(df_morph, metric_col, metric_label,
                             mod_name, vcv_df, hh, out_dir, atlas_level=3):
    os.makedirs(out_dir, exist_ok=True)
    available = df_morph["region"].unique() if df_morph is not None else []
    regions   = [r for r in REGIONS_OF_INTEREST if r in available]
    if not regions:
        print(f"  [{mod_name}] no regions ? skipped"); return

    individual = []
    for region in regions:
        safe = (region.replace(" ","_").replace("(","").replace(")","")
                      .replace("/","-").replace(",","")[:40])
        out = opj(out_dir, f"scaling_{mod_name}_{safe}.png")
        p = plot_scaling_with_neurons(df_morph, metric_col, metric_label,
                                      out, vcv_df, region, hh,
                                      atlas_level=atlas_level)
        if p: individual.append(p)

    assemble_vertical(individual, opj(out_dir, f"scaling_{mod_name}_composite.png"))


# ?????????????????????????????????????????????????????????????????????????????
# B) FC FINGERPRINT  ? exactly as in your 2020 script
#    collect per-species mean matrices ? normalize ? cosine similarity
#    + permutation test (Weibull fit) ? heatmap
# ?????????????????????????????????????????????????????????????????????????????

def _collect_species_mean_matrices(species_config):
    """
    For each species in species_config, load all subject matrices at
    ATLAS_LEVEL (LR version) and compute the per-species mean matrix.

    Returns
    -------
    dict: {species: {"mean": np.ndarray, "rois": list[str], "n": int}}
    """
    result = {}
    for sp, cfg in species_config.items():
        bids_dirs = cfg.get("bids_dirs", [])
        lk = cfg.get("list_to_keep",   [])
        lr = cfg.get("list_to_remove", [])

        all_mats  = []
        rois_ref  = None

        for bd in bids_dirs:
            recs = extract_corr_matrix_paths(
                bd, ATLAS_NAME, ATLAS_LEVEL, True, lk, lr)

            # one mean per subject (average over runs)
            sub_runs = {}
            for rec in recs:
                try:
                    rois, mat = load_corr_matrix(rec["path"])
                    if rois_ref is None: rois_ref = rois
                    key = (rec["subject"], rec["session"])
                    sub_runs.setdefault(key, []).append((rois, mat))
                except Exception as e:
                    print(f"  [WARN] {rec['path']}: {e}")

            for (sub, ses), run_list in sub_runs.items():
                # align to common ROIs across runs
                common = set(run_list[0][0])
                for r, _ in run_list[1:]: common &= set(r)
                ref = [x for x in run_list[0][0] if x in common]
                mats = []
                for r_rois, r_mat in run_list:
                    idx = [list(r_rois).index(x) for x in ref]
                    mats.append(r_mat[np.ix_(idx, idx)])
                all_mats.append((ref, np.nanmean(np.stack(mats, 0), axis=0)))

        if not all_mats:
            print(f"  [{sp}] no matrices found")
            continue

        # align to common ROIs across all subjects of this species
        common_sp = set(all_mats[0][0])
        for r, _ in all_mats[1:]: common_sp &= set(r)
        ref_sp = [x for x in all_mats[0][0] if x in common_sp]

        aligned = []
        for r_rois, r_mat in all_mats:
            idx = [list(r_rois).index(x) for x in ref_sp]
            aligned.append(r_mat[np.ix_(idx, idx)])

        stack = np.stack(aligned, axis=0)
        result[sp] = {"mean": np.nanmean(stack, axis=0),
                      "all":  stack,
                      "rois": ref_sp,
                      "n":    len(aligned)}
        print(f"  [{sp}] n={len(aligned)} subjects  "
              f"rois={len(ref_sp)}")

    # restrict to ROIs common across ALL species
    if not result: return result
    roi_sets = [set(v["rois"]) for v in result.values()]
    common_all = roi_sets[0].intersection(*roi_sets[1:])
    ref_order  = result[next(iter(result))]["rois"]
    common_rois = [r for r in ref_order if r in common_all]
    print(f"\n  Common ROIs across all species: {len(common_rois)}")

    for sp, d in result.items():
        idx = [d["rois"].index(r) for r in common_rois]
        d["mean"] = d["mean"][np.ix_(idx, idx)]
        d["all"]  = d["all"][:, idx, :][:, :, idx]
        d["rois"] = common_rois

    return result


def _fingerprint_normalize(mat, k):
    """
    For ROI k: take its row, remove self-connection (index k),
    normalize min-max.  Exactly as in the 2020 script.
    """
    row = mat[k].copy()
    row = np.delete(row, k)
    mn, mx = row.min(), row.max()
    if mx - mn < 1e-10:
        return np.zeros_like(row)
    return (row - mn) / (mx - mn)


def _cosine_similarity(x1, x2):
    """Cosine similarity between two vectors (same as in 2020 script)."""
    denom = math.sqrt(sum(a**2 for a in x1) * sum(b**2 for b in x2))
    if denom < 1e-12: return 0.0
    return sum(a * b for a, b in zip(x1, x2)) / denom


def _permutation_test_weibull(matA, matB, k, n_perm=100_000, rng=None):
    """
    Permutation test for cosine similarity of two groups (per ROI k).
    Exactly mirrors the 2020 script logic:
      - concatenate individual fingerprints, shuffle, split, compute similarity
      - fit Weibull, estimate p-value
    Returns (observed_similarity, p_value).
    """
    rng = rng or np.random.default_rng(42)

    # individual fingerprints per subject per ROI k
    fpA = np.array([_fingerprint_normalize(mat, k) for mat in matA])
    fpB = np.array([_fingerprint_normalize(mat, k) for mat in matB])

    # observed similarity of group means
    obs = _cosine_similarity(fpA.mean(0).tolist(), fpB.mean(0).tolist())

    # permutation distribution
    concat = np.concatenate([fpA, fpB], axis=0)
    nA     = len(fpA)
    sr     = []
    for _ in range(n_perm):
        perm = rng.permutation(len(concat))
        c    = concat[perm]
        mA   = c[:nA].mean(0).tolist()
        mB   = c[nA:].mean(0).tolist()
        sr.append(_cosine_similarity(mA, mB))

    # Weibull fit ? p-value (as in original script)
    shape, loc, scale = _stats.weibull_min.fit(sr)
    ej = ((obs - loc) / scale) ** shape * math.log10(math.e)
    p  = 10 ** (-ej)

    return obs, float(p)


def compute_fingerprint_similarity(species_data, n_perm=N_PERMUTATIONS):
    """
    For every pair of species and every ROI, compute:
      - cosine similarity of mean normalized fingerprints
      - permutation-test p-value (Weibull fit)

    Returns
    -------
    df_sim : DataFrame  [species_pair × ROI]  ? cosine similarity
    df_p   : DataFrame  [species_pair × ROI]  ? p-value
    """
    species_order = _phylo_sort(list(species_data.keys()))
    rois   = species_data[species_order[0]]["rois"]
    n_roi  = len(rois)
    pairs  = [(s1, s2) for i, s1 in enumerate(species_order)
                        for s2 in species_order[i+1:]]
    pair_labels = [f"{s1} vs {s2}" for s1, s2 in pairs]

    sim_mat = np.full((len(pairs), n_roi), np.nan)
    p_mat   = np.full((len(pairs), n_roi), np.nan)

    rng = np.random.default_rng(42)

    for pi, (s1, s2) in enumerate(pairs):
        matA = species_data[s1]["all"]   # (n_subj, n_roi, n_roi)
        matB = species_data[s2]["all"]
        print(f"  [{s1} vs {s2}]  nA={len(matA)}  nB={len(matB)}")
        for ki in range(n_roi):
            sim, p = _permutation_test_weibull(matA, matB, ki,
                                               n_perm=n_perm, rng=rng)
            sim_mat[pi, ki] = sim
            p_mat[pi, ki]   = p

    df_sim = pd.DataFrame(sim_mat, index=pair_labels, columns=rois)
    df_p   = pd.DataFrame(p_mat,   index=pair_labels, columns=rois)
    return df_sim, df_p


def plot_fingerprint_heatmap(df_sim, df_p, out_dir):
    """
    Two heatmaps side by side (or separate):
      Left  : cosine similarity  (plasma_r, vmin=0.8 vmax=1.0  as in original)
      Right : 1-p  (plasma_r, vmin=0 vmax=0.6  as in original)
    Also saves the DataFrames as CSV.
    """
    os.makedirs(out_dir, exist_ok=True)

    # shorten ROI labels for display
    short_rois = [r.replace("L_","").replace("R_","")
                   .replace("_(","(")[:22]
                  for r in df_sim.columns]

    for df, fname, title, vmin, vmax, cbar_lbl in [
        (df_sim, "fingerprint_similarity.png",
         "Cosine similarity of FC fingerprints", 0.80, 1.00, "Similarity"),
        (df_p.applymap(lambda x: 1 - x),
         "fingerprint_1minusp.png",
         "Fingerprint similarity  (1 ? p,  permutation test)", 0.00, 0.60,
         "1 ? p"),
    ]:
        fig, ax = plt.subplots(figsize=(max(14, len(df.columns)*0.22 + 3),
                                        max(4,  len(df.index)*0.55 + 1.5)))
        df_plot = df.copy()
        df_plot.columns = short_rois
        sns.heatmap(df_plot, ax=ax, annot=True, fmt=".2f",
                    annot_kws={"size": 5.5},
                    cmap="plasma_r", vmin=vmin, vmax=vmax,
                    cbar_kws={"label": cbar_lbl, "shrink": 0.6},
                    linewidths=0.3, linecolor="#eee")
        ax.set_title(title, fontweight="bold", fontsize=11, pad=10)
        ax.set_xlabel("ROI", fontsize=9)
        ax.set_ylabel("Species pair", fontsize=9)
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=6)
        plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)
        plt.tight_layout()
        out = opj(out_dir, fname)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  [plot] {out}")

    # also save a combined side-by-side figure
    fig, axes = plt.subplots(1, 2,
                             figsize=(max(28, len(df_sim.columns)*0.44 + 6),
                                      max(4,  len(df_sim.index)*0.55 + 2)))
    for ax, df, title, vmin, vmax, cbar_lbl in zip(
        axes,
        [df_sim, df_p.applymap(lambda x: 1-x)],
        ["Cosine similarity", "1 ? p  (permutation test)"],
        [0.80, 0.00], [1.00, 0.60],
        ["Similarity", "1 ? p"],
    ):
        df_plot = df.copy(); df_plot.columns = short_rois
        sns.heatmap(df_plot, ax=ax, annot=True, fmt=".2f",
                    annot_kws={"size": 5},
                    cmap="plasma_r", vmin=vmin, vmax=vmax,
                    cbar_kws={"label": cbar_lbl, "shrink": 0.6},
                    linewidths=0.3, linecolor="#eee")
        ax.set_title(title, fontweight="bold", fontsize=11, pad=8)
        ax.set_xlabel("ROI", fontsize=9)
        ax.set_ylabel("")
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=6)
        plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)

    plt.tight_layout()
    out = opj(out_dir, "fingerprint_combined.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [plot] {out}")

    # CSV exports
    df_sim.to_csv(opj(out_dir, "fingerprint_similarity.csv"))
    df_p.to_csv(opj(out_dir,   "fingerprint_pvalue.csv"))
    print(f"  [csv]  fingerprint_similarity.csv  /  fingerprint_pvalue.csv")


# ?????????????????????????????????????????????????????????????????????????????
# HELPERS
# ?????????????????????????????????????????????????????????????????????????????

def assemble_vertical(png_paths, output_path):
    imgs = [Image.open(p).convert("RGB")
            for p in png_paths if os.path.isfile(p)]
    if not imgs:
        warnings.warn(f"assemble_vertical: nothing to assemble ? {output_path}")
        return
    max_w   = max(im.width  for im in imgs)
    total_h = sum(im.height for im in imgs)
    canvas  = Image.new("RGB", (max_w, total_h), (255, 255, 255))
    y = 0
    for im in imgs:
        canvas.paste(im, (0, y)); y += im.height
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    canvas.save(output_path, dpi=(220, 220))
    print(f"  [composite] {output_path}  ({len(imgs)} panels)")


def _build_species_config():
    species_config     = {}
    atlas_label_paths  = {}
    atlas_fragments    = {}
    for sp, bd in SPECIES_BIDS.items():
        species_config[sp] = {"bids_dirs": [bd],
                               "list_to_keep": [], "list_to_remove": []}
    for sp, bds in SPECIES_MULTI_BIDS.items():
        species_config[sp] = {"bids_dirs": bds,
                               "list_to_keep": [], "list_to_remove": []}
    for sp in species_config:
        try:
            frag = find_species_path(ATLAS_LIB, sp)
            atlas_fragments[sp]   = frag
            atlas_label_paths[sp] = get_atlas_label_path(
                ATLAS_LIB, frag, ATLAS_NAME, prefer_statslut=False)
        except Exception as e:
            print(f"  [WARN] {sp}: {e}")
    return species_config, atlas_label_paths, atlas_fragments


# ?????????????????????????????????????????????????????????????????????????????
# MAIN
# ?????????????????????????????????????????????????????????????????????????????

def main():
    print("=" * 65)
    print("  EDNiX PGLS + FC fingerprint")
    print("=" * 65)

    # 0. VCV + Herculano-Houzel
    print("\nLoading VCV ...")
    vcv_df, _ = load_newick_vcv(NEWICK_FILE)
    print(f"  VCV: {vcv_df.shape}")

    print("\nLoading Herculano-Houzel data ...")
    hh = load_herculano(HERCULANO_XLSX)

    # build species config once
    species_config, atlas_label_paths, atlas_fragments = _build_species_config()

    # ?? A: Morphometry ????????????????????????????????????????????????????
    print("\n" + "=" * 65)
    print("  A ? Morphometry PGLS")
    print("=" * 65)

    print("\nCollecting morphometry data ...")
    data = collect_multi_species(
        species_config,
        regions_of_interest     = REGIONS_OF_INTEREST,
        extract                 = ("surface", "volume", "thickness"),
        atlas_name              = ATLAS_NAME,
        atlas_label_paths       = atlas_label_paths,
        atlas_library_root      = ATLAS_LIB,
        species_atlas_fragments = atlas_fragments,
    )
    df_surface   = data.get("surface")
    df_volume    = data.get("volume")
    df_thickness = data.get("thickness")

    for name, df in [("surface", df_surface),
                     ("volume",  df_volume),
                     ("thickness", df_thickness)]:
        n = len(df) if df is not None and not df.empty else 0
        print(f"  {name}: {n} rows")

    # A1 ? all-pairs PGLS
    print("\n  [A1] All-pairs region PGLS at level 3 ...")
    os.makedirs(PAIRS_DIR, exist_ok=True)

    for df, metric_col, metric_label, mod in [
        (df_surface,   "surface_area_mm2", "Surface area (mm²)", "surface"),
        (df_volume,    "volume_mm3",        "Volume (mm³)",       "volume"),
        (df_thickness, "thickness_mm",      "Thickness (mm)",     "thickness"),
    ]:
        if df is None or df.empty:
            print(f"  [{mod}] no data ? skipped"); continue
        df_res = run_pgls_all_region_pairs(
            df_morph    = df,
            metric_col  = metric_col,
            vcv_df      = vcv_df,
            regions     = REGIONS_OF_INTEREST,
            atlas_level = ATLAS_LEVEL,
            hemisphere  = "bilateral",
            min_species = 4,
            log_x=True, log_y=True,
        )
        print(f"  [{mod}] {len(df_res)} pairs  "
              f"({df_res['significant_fdr05'].sum()} FDR<0.05)")
        if df_res.empty: continue

        mod_pairs_dir = opj(PAIRS_DIR, mod)
        os.makedirs(mod_pairs_dir, exist_ok=True)
        df_res.to_excel(opj(mod_pairs_dir, f"pgls_pairs_{mod}.xlsx"),
                        index=False)
        for _, row in df_res.iterrows():
            rx = row["region_x"].replace(" ","_").replace("/","-")[:25]
            ry = row["region_y"].replace(" ","_").replace("/","-")[:25]
            plot_pgls_region_pair(
                result_row   = row,
                df_morph     = df,
                metric_col   = metric_col,
                metric_label = metric_label,
                vcv_df       = vcv_df,
                output_path  = opj(mod_pairs_dir, f"{rx}__x__{ry}.png"),
                atlas_level  = ATLAS_LEVEL,
                hemisphere   = "bilateral",
            )

    # A2 ? allometric scaling with neurons
    print("\n  [A2] Allometric scaling (body / brain / neurons) ...")
    os.makedirs(SCALING_DIR, exist_ok=True)
    for df, metric_col, metric_label, mod in [
        (df_surface,   "surface_area_mm2", "Surface area (mm²)", "surface"),
        (df_volume,    "volume_mm3",        "Volume (mm³)",       "volume"),
        (df_thickness, "thickness_mm",      "Thickness (mm)",     "thickness"),
    ]:
        if df is None or df.empty: continue
        run_morphometry_scaling(df, metric_col, metric_label, mod,
                                vcv_df, hh,
                                out_dir=opj(SCALING_DIR, mod),
                                atlas_level=ATLAS_LEVEL)

    # ?? B: FC fingerprint ?????????????????????????????????????????????????
    print("\n" + "=" * 65)
    print("  B ? FC normalized fingerprint")
    print("=" * 65)

    print("\nCollecting per-species mean connectivity matrices ...")
    species_data = _collect_species_mean_matrices(species_config)

    if len(species_data) < 2:
        print("  Not enough species with matrices ? fingerprint skipped")
    else:
        print(f"\nRunning pairwise similarity + permutation test "
              f"({N_PERMUTATIONS:,} shuffles) ...")
        df_sim, df_p = compute_fingerprint_similarity(
            species_data, n_perm=N_PERMUTATIONS)

        os.makedirs(FINGERPRINT_DIR, exist_ok=True)
        df_sim.to_csv(opj(FINGERPRINT_DIR, "fingerprint_similarity.csv"))
        df_p.to_csv(opj(FINGERPRINT_DIR,   "fingerprint_pvalue.csv"))
        plot_fingerprint_heatmap(df_sim, df_p, FINGERPRINT_DIR)

    print(f"\nDone.\n"
          f"  Pairs      : {PAIRS_DIR}/\n"
          f"  Scaling    : {SCALING_DIR}/\n"
          f"  Fingerprint: {FINGERPRINT_DIR}/")


if __name__ == "__main__":
    main()