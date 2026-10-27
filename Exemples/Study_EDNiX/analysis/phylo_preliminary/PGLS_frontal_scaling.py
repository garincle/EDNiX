"""
PGLS -- frontal lobe scaling (preliminary phylogenetic analysis, not in the paper).

Former second half of Figure_2_anat.py. Seven species are too few for a stable
estimate of Pagel's lambda (Freckleton, Harvey & Pagel 2002, Am Nat 160:712), so
these results are kept as exploratory.

  A) frontal lobe (atlas level 2) vs body weight, brain weight, neuron count
     (Herculano-Houzel compilation) -- one figure per modality, 3 panels
  B) frontal vs occipital lobe -- 1 figure, 3 panels (surface, volume, thickness)

PGLS with lambda estimated by ML on the tree NEWICK_FILE; species means of the
bilateral, one-value-per-subject tables of paper_data (same input as Figure 2).
Output: <PHYLO_DIR>/PGLS_frontal_scaling/{figures/cross_species/brain_scaling, stats}
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as tck

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_data as pdat
from Exemples.Study_EDNiX.analysis import paper_style as ps
from Statistics.Evo.PGLS import (
    load_newick_vcv,
    pgls_regression,
    pgls_on_ax,
    PAPER_RC,
    _phylo_sort,
    _species_colors,
    _pgls_on_ax_with_prediction
)

opj = os.path.join
OUT_DIR = opj(pdat.PHYLO_DIR, 'PGLS_frontal_scaling')
NEWICK_FILE = pdat.NEWICK_FILE
HERCULANO_XLSX = pdat.HERCULANO_XLSX
SPECIES_COLORS = ps.SPECIES_COLORS

# bilateral, one value per subject, every atlas level (the helpers select level 2)
_morpho = pdat.load_morphometry(level=None)
df_surface_plot = _morpho['surface']
df_volume_plot = _morpho['volume']
df_thickness_plot = _morpho['thickness']

# ── Output directories ────────────────────────────────────────────────────────
pgls_dir = opj(OUT_DIR, 'figures', 'cross_species', 'brain_scaling')
pgls_stats_dir = opj(OUT_DIR, 'stats')
os.makedirs(pgls_dir, exist_ok=True)
os.makedirs(pgls_stats_dir, exist_ok=True)

print("\n" + "=" * 64)
print("  PGLS — Evolutionary comparative analyses")
print("=" * 64)

# ── Load phylogeny ──────────────────────────────────────────────────────────
vcv_df, _ = load_newick_vcv(NEWICK_FILE)


# ── Helpers ─────────────────────────────────────────────────────────────────
def _species_mean(df, metric_col, region, atlas_level, hemisphere="bilateral"):
    """Return {species: mean} for a region."""
    if df is None or df.empty:
        return {}
    d = df.copy()
    if "atlas_level" in d.columns:
        d = d[d["atlas_level"] == atlas_level]
    if "hemisphere" in d.columns and hemisphere != "all":
        h = d[d["hemisphere"] == hemisphere]
        d = h if not h.empty else d[d["hemisphere"] == "bilateral"]
    d = d[d["region"] == region]
    return {} if d.empty else d.groupby("species")[metric_col].mean().to_dict()


def _load_herculano_from_excel(xlsx_path):
    """
    Load body weight, brain weight, neuron count from Herculano-Houzel Excel file.
    Uses T5 – Whole Brain sheet primarily.
    """
    # EDNiX species to Excel species mapping
    SPECIES_MAP = {
        "Mouse": ["Mus musculus"],
        "Rat": ["Rattus norvegicus"],
        "Mouselemur": ["Microcebus murinus"],
        "Marmoset": ["Callithrix jacchus"],
        "Macaque": ["Macaca mulatta", "Macaca fascicularis", "Macaca radiata"],
        "Human": ["Homo sapiens"],
        "Dog": ["Canis familiaris"],
        "Cat": ["Felis catus"],
        "Bat": ["Pteropus vampyrus"],
        "Pig": ["Sus scrofa domesticus"],
    }

    body_g, brain_g, neuron_M = {}, {}, {}

    if not os.path.exists(xlsx_path):
        print(f"  [Herculano] File not found: {xlsx_path}")
        return body_g, brain_g, neuron_M

    # ── Read T5 – Whole Brain ──────────────────────────────────────────────────
    try:
        df_raw = pd.read_excel(xlsx_path, sheet_name="T5 – Whole Brain", header=None)

        header_row = None
        for idx, row in df_raw.iterrows():
            row_str = ' '.join([str(x).strip() for x in row.values if pd.notna(x)])
            if 'Species' in row_str and 'Body mass' in row_str and 'Brain mass' in row_str:
                header_row = idx
                break

        if header_row is not None:
            df_t5 = pd.read_excel(xlsx_path, sheet_name="T5 – Whole Brain", header=header_row)
            df_t5.columns = [str(c).replace('\n', ' ').strip() for c in df_t5.columns]

            # Print columns for debugging
            print(f"  [Herculano] T5 columns: {df_t5.columns.tolist()}")

            sp_col = None
            body_col = None
            brain_col = None
            neuron_col = None

            for col in df_t5.columns:
                col_clean = col.lower().strip()
                if 'species' in col_clean:
                    sp_col = col
                elif 'body mass' in col_clean:
                    body_col = col
                elif 'brain mass' in col_clean:
                    brain_col = col
                elif 'neurons' in col_clean and 'mg' not in col_clean and '%' not in col_clean:
                    neuron_col = col

            print(f"  [Herculano] neuron_col = '{neuron_col}'")

            if sp_col and body_col and brain_col and neuron_col:
                for idx, row in df_t5.iterrows():
                    sp_name = str(row[sp_col]).strip().lower() if pd.notna(row[sp_col]) else ""
                    if not sp_name:
                        continue

                    for ednix_name, excel_names in SPECIES_MAP.items():
                        if any(excel_name.lower() in sp_name for excel_name in excel_names):
                            # Body mass
                            if body_col and pd.notna(row[body_col]) and ednix_name not in body_g:
                                try:
                                    val_str = str(row[body_col]).replace(',', '').strip()
                                    if '±' in val_str:
                                        val_str = val_str.split('±')[0]
                                    val = float(val_str)
                                    body_g[ednix_name] = val
                                except:
                                    pass

                            # Brain mass
                            if brain_col and pd.notna(row[brain_col]) and ednix_name not in brain_g:
                                try:
                                    val_str = str(row[brain_col]).replace(',', '').strip()
                                    if '±' in val_str:
                                        val_str = val_str.split('±')[0]
                                    val = float(val_str)
                                    brain_g[ednix_name] = val
                                except:
                                    pass

                            # Raw neuron count - THIS IS THE KEY
                            if neuron_col and pd.notna(row[neuron_col]) and ednix_name not in neuron_M:
                                try:
                                    val_str = str(row[neuron_col]).replace(',', '').strip()
                                    print(f"  [DEBUG] {ednix_name}: neuron_str = '{val_str}'")
                                    if '±' in val_str:
                                        val_str = val_str.split('±')[0]
                                    # Remove ×10⁶ if present
                                    if '×10⁶' in val_str:
                                        val_str = val_str.replace('×10⁶', '')
                                    val = float(val_str)
                                    # If value is > 1e6, it's raw count, convert to millions
                                    if val > 1e6:
                                        val = val / 1e6
                                    neuron_M[ednix_name] = val
                                except Exception as e:
                                    print(f"  [DEBUG] Failed to parse {ednix_name}: {e}")
                            break

    except Exception as e:
        print(f"  [Herculano] Error reading T5: {e}")
        import traceback
        traceback.print_exc()

    # ── Read Carnivora sheet for Dog/Cat ──────────────────────────────────────
    try:
        df_carn = pd.read_excel(xlsx_path, sheet_name="Carnivora (2017)", header=1)
        df_carn.columns = [str(c).replace('\n', ' ').strip() for c in df_carn.columns]

        sp_col = None
        body_col = None
        brain_col = None
        neuron_col = None

        for col in df_carn.columns:
            col_clean = col.lower().strip()
            if 'species' in col_clean:
                sp_col = col
            elif 'body mass' in col_clean:
                body_col = col
            elif 'brain mass' in col_clean:
                brain_col = col
            elif 'n whole brain' in col_clean:
                neuron_col = col

        if sp_col and body_col and brain_col and neuron_col:
            for idx, row in df_carn.iterrows():
                sp_name = str(row[sp_col]).strip().lower() if pd.notna(row[sp_col]) else ""
                if 'dog' in sp_name and 'Dog' not in body_g:
                    body_g['Dog'] = float(str(row[body_col]).replace(',', '').strip()) * 1000
                    brain_g['Dog'] = float(str(row[brain_col]).replace(',', '').strip())
                    val_str = str(row[neuron_col]).replace('×10⁶', '').replace(',', '').strip()
                    neuron_M['Dog'] = float(val_str)
                elif 'cat' in sp_name and 'Cat' not in body_g:
                    body_g['Cat'] = float(str(row[body_col]).replace(',', '').strip()) * 1000
                    brain_g['Cat'] = float(str(row[brain_col]).replace(',', '').strip())
                    val_str = str(row[neuron_col]).replace('×10⁶', '').replace(',', '').strip()
                    neuron_M['Cat'] = float(val_str)
    except Exception as e:
        print(f"  [Herculano] Carnivora sheet error: {e}")

    # ── Read Chiroptera sheet for Bat ─────────────────────────────────────────
    try:
        df_bat = pd.read_excel(xlsx_path, sheet_name="Chiroptera - Bats (2020)", header=1)
        df_bat.columns = [str(c).replace('\n', ' ').strip() for c in df_bat.columns]

        sp_col = None
        body_col = None
        brain_col = None
        neuron_col = None

        for col in df_bat.columns:
            col_clean = col.lower().strip()
            if 'species' in col_clean:
                sp_col = col
            elif 'body mass' in col_clean:
                body_col = col
            elif 'brain mass' in col_clean:
                brain_col = col
            elif 'n whole brain' in col_clean:
                neuron_col = col

        if sp_col and body_col and brain_col and neuron_col:
            for idx, row in df_bat.iterrows():
                sp_name = str(row[sp_col]).strip().lower() if pd.notna(row[sp_col]) else ""
                if 'pteropus vampyrus' in sp_name and 'Bat' not in body_g:
                    body_g['Bat'] = float(str(row[body_col]).replace(',', '').strip())
                    brain_g['Bat'] = float(str(row[brain_col]).replace(',', '').strip())
                    val_str = str(row[neuron_col]).replace('×10⁶', '').replace(',', '').strip()
                    neuron_M['Bat'] = float(val_str)
    except Exception as e:
        print(f"  [Herculano] Chiroptera sheet error: {e}")

    # ── Report ─────────────────────────────────────────────────────────────────
    all_species = ["Mouse", "Rat", "Mouselemur", "Marmoset", "Macaque", "Human", "Dog", "Cat", "Bat", "Pig"]

    print(f"\n  [Herculano] Loaded from Excel (raw neuron counts in millions):")
    print(f"    body_g:    {len(body_g)} species")
    print(f"    brain_g:   {len(brain_g)} species")
    print(f"    neuron_M:  {len(neuron_M)} species")

    missing_body = [s for s in all_species if s not in body_g]
    missing_brain = [s for s in all_species if s not in brain_g]
    missing_neuron = [s for s in all_species if s not in neuron_M]

    if missing_body:
        print(f"\n  ⚠️  MISSING body_g: {missing_body}")
    if missing_brain:
        print(f"  ⚠️  MISSING brain_g: {missing_brain}")
    if missing_neuron:
        print(f"  ⚠️  MISSING neuron_M: {missing_neuron}")

    print("\n  [Herculano] Values (whole brain neurons in millions):")
    print(f"    {'Species':12s} {'Body (g)':>10s} {'Brain (g)':>10s} {'Neurons (M)':>12s}")
    print(f"    {'-' * 12} {'-' * 10} {'-' * 10} {'-' * 12}")
    for s in all_species:
        b = body_g.get(s, 'MISSING')
        br = brain_g.get(s, 'MISSING')
        n = neuron_M.get(s, 'MISSING')
        b_str = f"{b:.0f}" if isinstance(b, (int, float)) else str(b)
        br_str = f"{br:.2f}" if isinstance(br, (int, float)) else str(br)
        n_str = f"{n:.0f}" if isinstance(n, (int, float)) else str(n)
        print(f"    {s:12s} {b_str:>10s} {br_str:>10s} {n_str:>12s}")

    return body_g, brain_g, neuron_M

# ── Load Herculano-Houzel data from Excel ──────────────────────────────────
body_g, brain_g, neuron_M = _load_herculano_from_excel(HERCULANO_XLSX)

# ── Modalities ──────────────────────────────────────────────────────────────
MODALITIES = [
    ("surface", df_surface_plot, "surface_area_mm2", "Surface area (mm²)"),
    ("volume", df_volume_plot, "volume_mm3", "Volume (mm³)"),
    ("thickness", df_thickness_plot, "thickness_mm", "Thickness (mm)"),
]

# ── Results storage ─────────────────────────────────────────────────────────
all_results = []

# ═══════════════════════════════════════════════════════════════════════════
# A) Frontal lobe scaling — 3 separate figures (surface, volume, thickness)
#    Each figure: 1x3 panels (body, brain, neuron)
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "─" * 64)
print("  A) Frontal lobe scaling — body / brain / neuron (Fig. 2B)")
print("─" * 64)

predictors = [
    ("body_weight", "Body weight (g)", body_g),
    ("brain_weight", "Brain weight (g)", brain_g),
    ("neuron_count", "Neuron count (×10⁶)", neuron_M),
]

for mod, df_mod, metric_col, metric_label in MODALITIES:
    if df_mod is None or df_mod.empty:
        print(f"  [{mod}] no data — skipped")
        continue

    # Get Frontal lobe means
    means = _species_mean(df_mod, metric_col, "Frontal_lobe",
                          atlas_level=2, hemisphere="bilateral")
    if len(means) < 3:
        print(f"  [{mod}] <3 species with Frontal data — skipped")
        continue

    sp_list = sorted(means.keys())
    y_vals = np.array([means[s] for s in sp_list])
    print(f"\n  [{mod}] n={len(sp_list)} species: {sp_list}")

    # Create figure: 1 row, 3 columns
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # Use the SPECIES_COLORS defined in the main script
    sp_col = SPECIES_COLORS

    for col, (pred_name, x_label, xdata) in enumerate(predictors):
        ax = axes[col]
        x_vals = np.array([xdata.get(s, np.nan) for s in sp_list])

        # Filter finite values
        valid = [(s, xv, yv) for s, xv, yv in zip(sp_list, x_vals, y_vals)
                 if np.isfinite(xv) and np.isfinite(yv)]
        if len(valid) < 3:
            ax.text(0.5, 0.5, "Insufficient data", ha="center", va="center",
                    transform=ax.transAxes, fontsize=10)
            ax.set_xscale("log")
            ax.set_yscale("log")
            continue

        sp_sub, x_sub, y_sub = zip(*valid)
        x_arr = np.array(x_sub)
        y_arr = np.array(y_sub)

        # Scatter
        for s, xv, yv in zip(sp_sub, x_arr, y_arr):
            ax.scatter(xv, yv, color=sp_col.get(s, "#888"), s=55, zorder=5,
                       edgecolors="k", linewidths=0.5)

        # PGLS line with PREDICTION intervals (wider than confidence intervals)
        res = pgls_regression(x_arr, y_arr, vcv_df, list(sp_sub),
                              log_x=True, log_y=True)
        if res:
            # Use the modified pgls_on_ax with prediction intervals
            _pgls_on_ax_with_prediction(ax, res, log_x=True, log_y=True)

            all_results.append(dict(
                modality=mod, predictor=pred_name, **res,
                region="Frontal_lobe", atlas_level=2
            ))

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(x_label, fontsize=10)
        ax.xaxis.set_major_formatter(tck.LogFormatterSciNotation())
        ax.yaxis.set_major_formatter(tck.LogFormatterSciNotation())

        # Only first column gets y-label
        if col == 0:
            ax.set_ylabel(f"Frontal lobe {metric_label}", fontsize=10)

    # Legend
    handles = [mpatches.Patch(facecolor=sp_col.get(s, "#888"), label=s)
               for s in _phylo_sort(sp_list)]
    fig.legend(handles=handles, fontsize=7,
               bbox_to_anchor=(0.5, -0.02), loc="upper center", ncol=5)

    plt.tight_layout()
    out_png = opj(pgls_dir, f"Fig2B_frontal_scaling_{mod}.png")
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [plot] {out_png}")

# ═══════════════════════════════════════════════════════════════════════════
# B) Frontal vs Occipital — 1 figure with 3 panels (surface, volume, thickness)
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "─" * 64)
print("  B) Frontal vs Occipital (Fig. 2C)")
print("─" * 64)

# Collect data for all modalities
occipital_data = {}
for mod, df_mod, metric_col, metric_label in MODALITIES:
    if df_mod is None or df_mod.empty:
        continue

    means_frontal = _species_mean(df_mod, metric_col, "Frontal_lobe",
                                  atlas_level=2, hemisphere="bilateral")
    means_occipital = _species_mean(df_mod, metric_col, "Occipital_lobe",
                                    atlas_level=2, hemisphere="bilateral")

    common_sp = sorted(set(means_frontal) & set(means_occipital) & set(vcv_df.index))
    if len(common_sp) < 3:
        print(f"  [{mod}] <3 species with both Frontal & Occipital — skipped")
        continue

    occipital_data[mod] = {
        'x': means_occipital,
        'y': means_frontal,
        'species': common_sp,
        'label': metric_label
    }

if occipital_data:
    # Create 1x3 grid
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    sp_col = SPECIES_COLORS

    for idx, (mod, data) in enumerate(occipital_data.items()):
        ax = axes[idx]
        sp_list = data['species']
        x_occ = np.array([data['x'][s] for s in sp_list])
        y_frontal = np.array([data['y'][s] for s in sp_list])

        # Scatter
        for s, xv, yv in zip(sp_list, x_occ, y_frontal):
            ax.scatter(xv, yv, color=sp_col.get(s, "#888"), s=55, zorder=5,
                       edgecolors="k", linewidths=0.5)

        # PGLS with prediction intervals
        res = pgls_regression(x_occ, y_frontal, vcv_df, sp_list,
                              log_x=True, log_y=True)
        if res:
            _pgls_on_ax_with_prediction(ax, res, log_x=True, log_y=True)
            all_results.append(dict(
                modality=mod, predictor="occipital_lobe", **res,
                region_x="Occipital_lobe", region_y="Frontal_lobe", atlas_level=2
            ))

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(f"Occipital lobe {data['label']}", fontsize=10)
        ax.xaxis.set_major_formatter(tck.LogFormatterSciNotation())
        ax.yaxis.set_major_formatter(tck.LogFormatterSciNotation())

        # Only first column gets y-label
        if idx == 0:
            ax.set_ylabel(f"Frontal lobe {data['label']}", fontsize=10)

    # Legend
    first_sp = list(occipital_data.values())[0]['species']
    handles = [mpatches.Patch(facecolor=sp_col.get(s, "#888"), label=s)
               for s in _phylo_sort(first_sp)]
    fig.legend(handles=handles, fontsize=7,
               bbox_to_anchor=(0.5, -0.02), loc="upper center", ncol=5)

    plt.tight_layout()
    out_png = opj(pgls_dir, "Fig2C_frontal_vs_occipital_all_modalities.png")
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [plot] {out_png}")

# ═══════════════════════════════════════════════════════════════════════════
# Export stats
# ═══════════════════════════════════════════════════════════════════════════

if all_results:
    df_results = pd.DataFrame(all_results)
    out_csv = opj(pgls_stats_dir, "pgls_frontal_scaling.csv")
    df_results.to_csv(out_csv, index=False)

    # Print summary table
    print("\n" + "─" * 64)
    print("  SUMMARY — PGLS results")
    print("─" * 64)
    summary_cols = ["modality", "predictor", "slope", "r_pgls", "p", "lambda_pagel", "n_species"]
    summary = df_results[summary_cols]
    print(summary.to_string(index=False))
    print(f"\n  [stats] {out_csv}")

print("\n  PGLS done.")
print(f"  Figures → {pgls_dir}")
print(f"  Stats   → {pgls_stats_dir}")
