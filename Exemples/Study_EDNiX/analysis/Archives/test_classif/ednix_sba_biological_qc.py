"""
EDNiX — Biological SBA QC Metrics
=====================================

Three anatomically-grounded QC metrics computed from the ALREADY-PRODUCED
retrosplenial seed-based-analysis (SBA) correlation maps, as a complement
to the matrix-based dual-threshold classification (ednix_scoring_core.py).

No raw 4D BOLD needed — everything is derived from:
  - the seed's r-map / Fisher-z map already in Stats/SBA/{seed}/
  - the seed's own anatomical mask (eroded version preferred)
  - the pipeline's own `higher_threshold.nii.gz` significance map
  - the atlas 4D volume (index 2 = "niveau 3") + its LUT, to build the
    somatosensory ROI mask

Metrics
-------
1. Intra-somato coherence: dispersion of r-values inside the somatosensory
   mask (median, std, fraction negative). A genuine anticorrelated network
   should give consistently negative, low-dispersion values; noise gives
   scattered values with no consistent sign.

2. Spatial contiguity: using higher_threshold.nii.gz as the pipeline's own
   significance map, find the connected component touching the seed's
   anatomical mask, and compute its Dice overlap with that mask. A focal,
   biologically plausible map gives high Dice; a diffuse/scattered map
   gives low Dice (large "significant" volume barely overlapping the seed).

3. Anticorrelation distance: median(r) inside the seed's own mask (near
   self-correlation, ~1) minus median(r) inside the somatosensory mask.
   A clean, well-differentiated network gives a large positive distance;
   a noisy acquisition gives a small/near-zero distance.

Assumptions to verify on your side
-----------------------------------
- The level-index-2 atlas volume and the SBA correlation maps share the
  same affine/voxel grid (both are in "templates/EDNiX" space) — a
  resample fallback is included but untested against real data.
- `higher_threshold.nii.gz` — CONFIRMED: produced via
  `nilearn.image.threshold_img(fisher_z_map, percentile_threshold)`, which
  keeps the original Fisher-z value where |z| exceeds the percentile cutoff
  and zeros it elsewhere. So `data != 0` correctly identifies the
  significant voxels — no change needed.
- Retrosplenial seed folder is found by globbing "*retrosplenial*" and
  preferring the R_ hemisphere version (matching your example path).
  Override with `seed_name_pattern` if a study uses the left hemisphere.

Usage
-----
    import ednix_sba_biological_qc as bio

    lut = bio.parse_lut("/path/to/EDNIxCSCLR_label.txt")
    feats = bio.compute_biological_qc_features(
        sba_folder="/scratch2/EDNiX/Dog/BIDS_Boch_K9/sub-18/ses-1/func/templates/EDNiX/Stats/SBA/R_BA_26_29_30__retrosplenial",
        atlas_4d_path="/scratch2/EDNiX/Dog/BIDS_Boch_K9/sub-18/ses-1/func/templates/EDNiX/labels/18_seg-EDNIxCSCLR_dseg.nii.gz",
        lut=lut)

Batch over a whole species/BIDS dir:
    df = bio.batch_compute_sba_qc("/scratch2/EDNiX", species="Dog", bids_dir="BIDS_Boch_K9",
                                    lut_path="/path/to/EDNIxCSCLR_label.txt")
    bio.export_to_xlsx(df, "sba_biological_qc_Dog_BochK9.xlsx")
"""

import os
import re
import glob
import json
import warnings

import numpy as np
import pandas as pd
import nibabel as nib
from scipy import ndimage

opj = os.path.join
ope = os.path.exists

LEVEL_INDEX = 2  # 0-based index into the 4D atlas ("niveau 3" per Clement)

# =============================================================================
# LUT PARSING
# =============================================================================
def load_3d_image(path):
    """Charge une image NIfTI et la retourne en 3D (squeeze la dimension temporelle si présente)."""
    img = nib.load(path)
    data = img.get_fdata()
    if data.ndim == 4 and data.shape[-1] == 1:
        data = np.squeeze(data, axis=-1)
        img = nib.Nifti1Image(data, img.affine, img.header)
    return img

def parse_lut(lut_path):
    """
    Parse the EDNIxCSC(LR)_label.txt format:
        RegionName
        label R G B A
        RegionName2
        label2 R G B A
        ...
    Returns {region_name: label_value}.
    """
    with open(lut_path) as f:
        lines = [l.strip() for l in f if l.strip()]
    lut = {}
    i = 0
    while i < len(lines) - 1:
        name = lines[i]
        parts = lines[i + 1].split()
        try:
            label_val = int(parts[0])
            lut[name] = label_val
            i += 2
        except (ValueError, IndexError):
            i += 1
    return lut


def find_labels(lut, name_contains, hemisphere=None):
    """hemisphere: None, 'L', or 'R'. Returns list of (name, label_value)."""
    matches = []
    for name, val in lut.items():
        if hemisphere == "L" and not name.startswith("L_"):
            continue
        if hemisphere == "R" and not name.startswith("R_"):
            continue
        if name_contains.lower() in name.lower():
            matches.append((name, val))
    return matches

# =============================================================================
# MASK / MAP LOADING
# =============================================================================

def get_atlas_level_volume(atlas_4d_path, level_index=LEVEL_INDEX):
    img = nib.load(atlas_4d_path)
    data = img.get_fdata()
    if data.ndim != 4:
        raise ValueError(f"Expected a 4D atlas at {atlas_4d_path}, got shape {data.shape}")
    return data[..., level_index].astype(int), img.affine


def build_roi_mask(atlas_level_vol, lut, name_contains, hemisphere=None):
    """Bilateral (or hemisphere-restricted) boolean mask for a named ROI."""
    matches = find_labels(lut, name_contains, hemisphere=hemisphere)
    if not matches:
        raise ValueError(f"No LUT labels matched '{name_contains}' (hemisphere={hemisphere})")
    values = [v for _, v in matches]
    return np.isin(atlas_level_vol, values), matches


def _resample_if_needed(source_img, reference_shape, reference_affine, label="volume"):
    """If source_img's grid doesn't match the reference, resample (nearest
    neighbour) onto it. Only triggers if shapes/affines actually differ."""
    if source_img.shape == reference_shape and np.allclose(source_img.affine, reference_affine, atol=1e-3):
        return source_img.get_fdata()
    warnings.warn(f"[resample] {label} grid does not match reference — resampling (nearest neighbour). "
                  f"Verify this is expected for your data.")
    from nilearn.image import resample_img
    resampled = resample_img(source_img, target_affine=reference_affine,
                              target_shape=reference_shape, interpolation="nearest")
    return resampled.get_fdata()


def find_sba_folder(subject_session_root, seed_name_pattern="*retrosplenial*"):
    """subject_session_root: .../func/templates/EDNiX/Stats/SBA/
    Returns the matched seed folder path (prefers R_ hemisphere if both exist)."""
    candidates = sorted(glob.glob(opj(subject_session_root, seed_name_pattern)))
    candidates = [c for c in candidates if os.path.isdir(c)]
    if not candidates:
        raise FileNotFoundError(f"No SBA seed folder matching {seed_name_pattern} under {subject_session_root}")
    r_candidates = [c for c in candidates if os.path.basename(c).startswith("R_")]
    return (r_candidates or candidates)[0]


def load_correlation_maps_per_run(sba_folder):
    """
    Retourne une liste de dict pour chaque fichier *_correlations.nii.gz trouvé.
    Chaque dict contient : r_data, z_data, run_id, run_idx, affine, shape.
    """
    corr_files = glob.glob(opj(sba_folder, "*_correlations.nii.gz"))
    # Exclure les fichiers fish (ils ont le même nom de base)
    corr_files = [f for f in corr_files if "_fish" not in f]
    if not corr_files:
        raise FileNotFoundError(f"Aucun fichier *_correlations.nii.gz dans {sba_folder}")

    # Trier pour avoir un ordre déterministe
    corr_files = sorted(corr_files)

    runs = []
    for corr_path in corr_files:
        fname = os.path.basename(corr_path)
        run_match = re.search(r'run-(\d+)', fname)
        if run_match:
            run_id_str = run_match.group(1)
            run_idx = int(run_id_str) - 1  # 0‑based pour correspondre à la matrice
        else:
            # Pas de run dans le nom : on suppose que c'est le premier (et unique) run
            run_idx = 0
            run_id_str = "0"  # ou None, mais on garde une chaîne pour l'affichage

        # Chargement 3D pour la carte r
        r_img = load_3d_image(corr_path)
        r_data = r_img.get_fdata()

        # Chargement 3D pour la carte Fisher-z si elle existe
        z_path = corr_path.replace("_correlations.nii.gz", "_correlations_fish.nii.gz")
        if os.path.exists(z_path):
            z_img = load_3d_image(z_path)
            z_data = z_img.get_fdata()
        else:
            z_data = np.arctanh(np.clip(r_data, -0.999, 0.999))

        runs.append({
            'run_id': run_id_str,          # chaîne (ex: "1", "2", ou "0")
            'run_idx': run_idx,            # entier 0‑based
            'r_data': r_data,
            'z_data': z_data,
            'affine': r_img.affine,
            'shape': r_img.shape
        })
    return runs

def load_seed_mask(sba_folder, prefer_eroded=True):
    """Loads the seed's own anatomical mask (not the correlation map)."""
    exclude = ("correlations", "higher_threshold", "rsp")
    all_masks = [p for p in glob.glob(opj(sba_folder, "*.nii.gz"))
                 if not any(x in os.path.basename(p) for x in exclude)]
    eroded = [p for p in all_masks if "eroded" in p]
    plain = [p for p in all_masks if "eroded" not in p]
    chosen = (eroded if prefer_eroded and eroded else plain or eroded)
    if not chosen:
        raise FileNotFoundError(f"No seed mask file found in {sba_folder} (candidates: {all_masks})")
    # Chargement 3D
    img = load_3d_image(chosen[0])   # <-- CORRECTION
    return (img.get_fdata() > 0), img, chosen[0]


def load_higher_threshold_mask(sba_folder, reference_affine, reference_shape):
    """
    Loads higher_threshold.nii.gz — produced by
    nilearn.image.threshold_img(fisher_z_map, percentile_threshold), so
    nonzero voxels are exactly the significant ones (original Fisher-z
    value kept above threshold, zeroed below).
    """
    path = opj(sba_folder, "higher_threshold.nii.gz")
    if not ope(path):
        raise FileNotFoundError(f"higher_threshold.nii.gz not found in {sba_folder}")
    img = nib.load(path)
    data = _resample_if_needed(img, reference_shape, reference_affine, label="higher_threshold")
    binary = data != 0
    return binary

# =============================================================================
# METRIC 1 — INTRA-SOMATO COHERENCE
# =============================================================================

def metric_intra_somato_coherence(r_map, somato_mask):
    vals = r_map[somato_mask]
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return dict(median_r_somato=np.nan, std_r_somato=np.nan,
                    iqr_r_somato=np.nan, frac_negative_somato=np.nan, n_voxels_somato=0)
    return dict(
        median_r_somato=float(np.median(vals)),
        std_r_somato=float(np.std(vals)),
        iqr_r_somato=float(np.percentile(vals, 75) - np.percentile(vals, 25)),
        frac_negative_somato=float(np.mean(vals < 0)),
        n_voxels_somato=int(vals.size),
    )

# =============================================================================
# METRIC 2 — SPATIAL CONTIGUITY (Dice of significant cluster vs seed mask)
# =============================================================================

def metric_spatial_contiguity(sig_mask, seed_anat_mask):
    labeled, n_clusters = ndimage.label(sig_mask)
    if n_clusters == 0:
        return dict(contiguity_dice=0.0, sig_cluster_voxels=0, seed_voxels=int(seed_anat_mask.sum()))

    overlap_labels = np.unique(labeled[seed_anat_mask & sig_mask])
    overlap_labels = overlap_labels[overlap_labels != 0]
    if overlap_labels.size == 0:
        # no significant cluster touches the seed at all
        return dict(contiguity_dice=0.0, sig_cluster_voxels=0, seed_voxels=int(seed_anat_mask.sum()))

    cluster_mask = np.isin(labeled, overlap_labels)
    intersection = np.sum(cluster_mask & seed_anat_mask)
    dice = 2.0 * intersection / (cluster_mask.sum() + seed_anat_mask.sum())
    return dict(
        contiguity_dice=float(dice),
        sig_cluster_voxels=int(cluster_mask.sum()),
        seed_voxels=int(seed_anat_mask.sum()),
    )

# =============================================================================
# METRIC 3 — ANTICORRELATION DISTANCE
# =============================================================================

def metric_anticorrelation_distance(r_map, seed_anat_mask, somato_mask):
    seed_vals = r_map[seed_anat_mask]
    seed_vals = seed_vals[np.isfinite(seed_vals)]
    somato_vals = r_map[somato_mask]
    somato_vals = somato_vals[np.isfinite(somato_vals)]
    median_seed = float(np.median(seed_vals)) if seed_vals.size else np.nan
    median_somato = float(np.median(somato_vals)) if somato_vals.size else np.nan
    distance = (median_seed - median_somato) if np.isfinite(median_seed) and np.isfinite(median_somato) else np.nan
    return dict(median_r_seed=median_seed, median_r_somato_dup=median_somato, anticorr_distance=distance)

# =============================================================================
# ORCHESTRATOR — ONE RUN
# =============================================================================
def compute_biological_qc_features(sba_folder, atlas_4d_path, lut,
                                   level_index=LEVEL_INDEX, somato_hemisphere=None):
    """Calcule les métriques pour chaque run du dossier SBA. Retourne une liste de dict."""
    # Charger le volume atlas et le masque somatosensoriel (une seule fois)
    level_vol, atlas_affine = get_atlas_level_volume(atlas_4d_path, level_index=level_index)
    somato_mask, somato_labels_used = build_roi_mask(level_vol, lut, "Somatosensory_cortex",
                                                     hemisphere=somato_hemisphere)

    # Charger les cartes de corrélation pour chaque run
    run_data_list = load_correlation_maps_per_run(sba_folder)

    # Charger le masque anatomique de la graine (une seule fois, commun à tous les runs)
    seed_mask, seed_img, seed_mask_path = load_seed_mask(sba_folder)
    sba_shape = seed_img.shape
    sba_affine = seed_img.affine

    # Resampler le masque somato sur la grille SBA si nécessaire
    if somato_mask.shape != sba_shape or not np.allclose(atlas_affine, sba_affine, atol=1e-3):
        somato_img = nib.Nifti1Image(somato_mask.astype(np.uint8), atlas_affine)
        somato_mask = _resample_if_needed(somato_img, sba_shape, sba_affine, label="somato_mask") > 0

    results = []
    for run_info in run_data_list:
        run_id = run_info['run_id']
        r_data = run_info['r_data']
        z_data = run_info['z_data']

        # --- Recalcul du seuil à partir de la carte Fisher-z (méthode percentile 99) ---
        z_vals = z_data[np.isfinite(z_data)]
        if z_vals.size == 0:
            warnings.warn(f"run {run_id}: carte Fisher-z vide ou NaN")
            continue
        threshold = np.percentile(np.abs(z_vals), 99)  # seuil absolu
        sig_mask = np.abs(z_data) >= threshold  # binaire, significatif

        # --- Métriques ---
        feats = {}
        feats.update(metric_intra_somato_coherence(r_data, somato_mask))
        feats.update(metric_spatial_contiguity(sig_mask, seed_mask))
        feats.update(metric_anticorrelation_distance(r_data, seed_mask, somato_mask))
        feats["run"] = run_id
        feats["run_idx"] = run_info['run_idx']
        feats["sba_folder"] = sba_folder
        feats["seed_mask_path"] = seed_mask_path
        feats["somato_labels_used"] = ", ".join(n for n, _ in somato_labels_used)
        results.append(feats)

    return results

# =============================================================================
# BATCH OVER SPECIES / BIDS DIR
# =============================================================================

def _find_atlas_path(subject_session_root_templates_ednix, subject_id, atlas_name="EDNIxCSCLR"):
    """subject_session_root_templates_ednix: .../func/templates/EDNiX/
    Looks for labels/{subject_id}_seg-{atlas_name}_dseg.nii.gz"""
    candidates = glob.glob(opj(subject_session_root_templates_ednix, "labels",
                               f"*_seg-{atlas_name}_dseg.nii.gz"))
    if len(candidates) != 1:
        raise FileNotFoundError(
            f"Expected exactly one {atlas_name} dseg.nii.gz under "
            f"{opj(subject_session_root_templates_ednix, 'labels')}, found {candidates}")
    return candidates[0]

def batch_compute_sba_qc(bids_root, lut_path, species=None, bids_dir=None,
                         atlas_name="EDNIxCSCLR", seed_name_pattern="*retrosplenial*",
                         verbose=True):
    """
    Walk {bids_root}/{species}/{bids_dir}/sub-*/... and compute QC metrics
    for every subject/session/run that has an SBA retrosplenial folder + atlas.
    """
    lut = parse_lut(lut_path)

    species_glob = species or "*"
    bids_glob = bids_dir or "*"
    pattern = opj(bids_root, species_glob, bids_glob, "sub-*", "**", "func", "templates", "EDNiX")
    template_roots = sorted(set(glob.glob(pattern, recursive=True)))   # guaranteed to be a list

    if not template_roots:
        if verbose:
            print(f"[batch_compute_sba_qc] No template roots found for {species}/{bids_dir}")
        return pd.DataFrame()

    rows = []
    for troot in template_roots:
        rel = os.path.relpath(troot, bids_root)
        parts = rel.split(os.sep)
        if len(parts) < 3:
            continue   # skip malformed paths
        sp, bd = parts[0], parts[1]
        # subject and session
        sub_part = parts[2] if len(parts) > 2 else ""
        subject = sub_part.replace("sub-", "", 1) if sub_part.startswith("sub-") else ""
        ses_match = re.search(r"ses-([^/\\]+)", rel)
        session = ses_match.group(1) if ses_match else ""

        try:
            atlas_path = _find_atlas_path(troot, subject, atlas_name=atlas_name)
            sba_root = opj(troot, "Stats", "SBA")
            sba_folder = find_sba_folder(sba_root, seed_name_pattern)
            feats_list = compute_biological_qc_features(sba_folder, atlas_path, lut)
        except Exception as e:
            warnings.warn(f"[batch_compute_sba_qc] skipped {sp}/{bd}/sub-{subject}/ses-{session}: {e}")
            continue

        for feats in feats_list:
            row = {"species": sp, "bids_dir": bd, "subject": subject, "session": session}
            row.update(feats)
            rows.append(row)
            if verbose:
                print(f"  [ok] {sp}/{bd}/sub-{subject}/ses-{session} run={feats.get('run')}  "
                      f"dice={feats.get('contiguity_dice', 0.0):.3f}  "
                      f"anticorr_dist={feats.get('anticorr_distance', 0.0):.3f}")

    df = pd.DataFrame(rows)
    if verbose:
        print(f"\n[batch_compute_sba_qc] {len(df)} runs processed successfully")
    return df


def export_to_xlsx(df, output_path, verbose=True):
    os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)
    cols = ["species", "bids_dir", "subject", "session",
            "median_r_somato", "std_r_somato", "iqr_r_somato", "frac_negative_somato",
            "contiguity_dice", "sig_cluster_voxels", "seed_voxels",
            "median_r_seed", "anticorr_distance",
            "n_voxels_somato", "somato_labels_used", "sba_folder"]
    cols = [c for c in cols if c in df.columns]
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        df[cols].to_excel(writer, sheet_name="SBA_Biological_QC", index=False)
        if "species" in df.columns:
            summary = df.groupby("species")[
                ["median_r_somato", "contiguity_dice", "anticorr_distance"]
            ].agg(["mean", "std", "count"])
            summary.to_excel(writer, sheet_name="Species_Summary")
    if verbose:
        print(f"[export] wrote {output_path}")
