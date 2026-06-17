"""
EDNiX BIDS Analysis Tools  —  v7
=================================

ARCHITECTURE
------------
  §1  IMPORTS & CONSTANTS
  §2  UTILITIES            — path helpers, name normalisation, shell wrappers
  §3  ATLAS / LABEL        — parse StatsLUT / label files, legend xlsx, hierarchy
  §4  PATH FINDERS         — scan BIDS tree for xlsx / NIfTI / JSON files
  §5  STAGE 0              — wb_command surface area + thickness extraction
                             Key design: bulk extraction (8 calls/level/subject,
                             not 4*N_regions calls). Uses cifti-math-free approach:
                             dense scalars separated once per hemisphere, then
                             numpy masking per label ID — works in wb 1.5.0.
  §6  STAGE 1              — read xlsx / NIfTI / JSON → tidy DataFrames
  §7  STAGE 2              — export DataFrames → Excel / CSV
  §8  DIAGNOSTICS
  §9  PLOTTING HELPERS

KEY CHANGES v7
--------------
• Bulk extraction: separate dense maps ONCE per hemisphere per level, then
  use numpy to extract all region values — ~100× fewer wb_command calls.
• Label IDs always sourced from the dlabel itself (cifti-label-export-table),
  never from StatsLUT, so IDs are always correct regardless of atlas version.
• StatsLUT / label file used only for: (a) volume extraction label lookup,
  (b) hierarchy level assignment from legend xlsx, (c) hemisphere inference.
• Legend xlsx drives atlas_level assignment: each region is assigned to the
  level that matches its hierarchy depth (level1..level4).
• No cifti-label-to-roi -name calls (broken in wb 1.5.0 for special chars).
• No legend name remapping in Stage 0 (caused l_Allocortex mismatches).
• Bilateral regions (no L_/R_ prefix) handled: sum areas L+R, mean thickness
  pooling L+R vertices.
• Works with any atlas whose StatsLUT follows:  ID  NAME  R  G  B  A
  and whose dlabel uses either L_/R_ or l_/r_ hemisphere prefixes.
"""

# ═══════════════════════════════════════════════════════════════════════════════
# §1  IMPORTS & CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════════

import os
import re
import glob
import json
import math
import hashlib
import traceback
import warnings
import subprocess
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats

try:
    from joblib import Parallel, delayed
    _JOBLIB = True
except ImportError:
    _JOBLIB = False
    warnings.warn("joblib not found — parallel processing disabled.")

opj = os.path.join

N_ATLAS_LEVELS = 4
_TIMEOUT_SHORT = 120
_TIMEOUT_HEAVY = 600

PAPER_RC = {
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.linewidth":     1.2,
    "xtick.major.width":  1.2,
    "ytick.major.width":  1.2,
    "font.family":        "sans-serif",
    "font.size":          11,
    "axes.titlesize":     12,
    "axes.labelsize":     11,
    "legend.frameon":     False,
    "figure.dpi":         150,
    "figure.facecolor":   "white",
    "axes.facecolor":     "white",
}

PALETTE = [
    "#0072B2", "#E69F00", "#009E73", "#CC79A7",
    "#56B4E9", "#D55E00", "#F0E442", "#000000",
]

AWAKE_MARKER  = "o"
ANESTH_MARKER = "^"
_AWAKE_KEYWORDS = {"awake", "eveil", "éveil", "no anesthesia", "unanesthetized"}


# ═══════════════════════════════════════════════════════════════════════════════
# §2  UTILITIES
# ═══════════════════════════════════════════════════════════════════════════════

def _linux_path(p):
    return str(p).replace("\\", "/")


def _parse_bids_entity(path, entity):
    for part in Path(path).parts:
        if part.startswith(entity + "-"):
            return part.split("-", 1)[1]
    return None


def _subject_session_filter(all_subs, all_ses, list_to_keep, list_to_remove):
    pairs = list(zip(all_subs, all_ses))
    if list_to_keep:
        pairs = [p for p in pairs if p in list_to_keep]
    if list_to_remove:
        pairs = [p for p in pairs if p not in list_to_remove]
    return pairs


def _spco_timed(cmd, timeout=_TIMEOUT_SHORT):
    """Run shell command; return stdout or None on failure/timeout."""
    try:
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                           timeout=timeout)
        # Filter out the LD_PRELOAD vgl noise before checking failure
        real_stderr = "\n".join(
            l for l in r.stderr.splitlines()
            if "libdlfaker" not in l and "libvglfaker" not in l
        ).strip()
        if r.returncode != 0:
            warnings.warn(f"Command failed: {cmd}\n{real_stderr}")
            return None
        return r.stdout
    except subprocess.TimeoutExpired:
        warnings.warn(f"Command timed out after {timeout}s: {cmd}")
        return None


def _spco(cmd):
    return _spco_timed(cmd, _TIMEOUT_SHORT)


def _normalize_name(s):
    """
    Canonical form for fuzzy region matching.
    Strips L_/R_/l_/r_ prefix, lower-cases, collapses whitespace/underscores,
    removes punctuation.
    """
    s = re.sub(r"^[lLrR]_", "", str(s))
    s = s.lower()
    s = re.sub(r"[\s_]+", "_", s)
    s = re.sub(r"[(),.\-/\\]", "", s)
    return re.sub(r"_+", "_", s).strip("_")


def _names_match(a, b):
    """
    Region name match: True if normalised names are equal.

    Exact equality after normalisation (strip hemisphere prefix, lowercase,
    collapse punctuation). Substring matching intentionally NOT used —
    it causes false positives: e.g. 'isocortex' would match
    'orbital_proisocortex_and_preiallocortex' as a substring, adding
    spurious ~91mm3 dots to the Isocortex panel.

    'L_Isocortex' ~ 'Isocortex'  -> True  (prefix stripped)
    'Periallocortex ' ~ 'Periallocortex' -> True  (whitespace)
    'Orbital_proisocortex...' ~ 'Isocortex' -> False  (no longer a false positive)
    """
    return _normalize_name(a) == _normalize_name(b)


def _filter_regions(all_regions, patterns):
    if not patterns:
        return list(all_regions)
    return [r for r in all_regions
            if any(_names_match(r, p) for p in patterns)]


def _df_region_mask(df_col, patterns):
    if not patterns:
        return pd.Series([True] * len(df_col), index=df_col.index)
    return df_col.apply(lambda r: any(_names_match(r, p) for p in patterns))


def _scan_bids(bids_dir, sub_glob="sub-*", ses_glob="ses-*"):
    for sub_dir in sorted(glob.glob(opj(bids_dir, sub_glob))):
        sub = _parse_bids_entity(sub_dir, "sub")
        if sub is None:
            continue
        ses_dirs = sorted(glob.glob(opj(sub_dir, ses_glob)))
        if not ses_dirs:
            yield sub, None, sub_dir, sub_dir
        else:
            for sd in ses_dirs:
                ses = _parse_bids_entity(sd, "ses")
                yield sub, ses, sub_dir, sd


def _iter_subjects(bids_dir, list_to_keep, list_to_remove):
    for sub, ses, _, ses_dir in _scan_bids(bids_dir):
        ses_label = ses or "1"
        if list_to_keep   and (sub, ses_label) not in list_to_keep:   continue
        if list_to_remove and (sub, ses_label) in list_to_remove:     continue
        native_dir = opj(ses_dir, "anat", "native", "surfaces", "Native_resol")
        if not os.path.exists(native_dir):
            continue
        yield sub, ses_label, ses_dir, native_dir


def _parse_run(filename):
    m = re.search(r"run-(\w+)", os.path.basename(filename))
    return m.group(1) if m else None


def _read_gifti_data(filepath):
    import nibabel as nib
    img = nib.load(filepath)
    if hasattr(img, "darrays") and img.darrays:
        return img.darrays[0].data.ravel()
    return img.get_fdata().ravel()


def _bids_label(bids_dir):
    return os.path.basename(str(bids_dir).rstrip("/"))


def _bids_offsets(bids_dirs, spread=0.22):
    n = len(bids_dirs)
    if n <= 1:
        return {b: 0.0 for b in bids_dirs}
    offs = np.linspace(-spread / 2, spread / 2, n)
    return {b: float(o) for b, o in zip(bids_dirs, offs)}


def _anesthesia_marker(anesth_str):
    if not isinstance(anesth_str, str):
        return AWAKE_MARKER
    s = anesth_str.strip().lower()
    if not s or s in ("nan", "none") or any(k in s for k in _AWAKE_KEYWORDS):
        return AWAKE_MARKER
    return ANESTH_MARKER


def _phylo_sort(species_list, phylo_order=None):
    if phylo_order is None:
        phylo_order = ["Bat", "Rat", "Mouse", "Mouselemur",
                       "Marmoset", "Macaque", "Human", "Dog", "Pig", "Cat"]
    known   = [s for s in phylo_order if s in species_list]
    unknown = sorted(s for s in species_list if s not in phylo_order)
    return known + unknown


def _save_failed(failed, out_dir, modality):
    if not failed:
        return
    print(f"\n{'!' * 60}")
    print(f"  [{modality}] {len(failed)} subject(s) FAILED:")
    for e in failed:
        print(f"    • {e['subject']} ses={e['session']}  {e['reason']}")
    print(f"{'!' * 60}\n")
    os.makedirs(out_dir, exist_ok=True)
    log = opj(out_dir, f"failed_subjects_{modality}.txt")
    with open(log, "w") as fh:
        fh.write(f"# Failed — {modality}\n")
        for e in failed:
            fh.write(f"{e['subject']}\t{e['session']}\t{e['reason']}\n")
    print(f"  [{modality}] failure log → {log}")


# ═══════════════════════════════════════════════════════════════════════════════
# §3  ATLAS / LABEL HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

def parse_label_file(label_path):
    """
    Parse an EDNiX atlas label / StatsLUT file.

    Supported formats:
      StatsLUT:  ID  NAME  R  G  B  A      (one entry per line)
      label:     NAME\\nID R G B A\\n...   (alternating name / colour lines)

    Returns DataFrame with columns:
        label_id, region_name, base_name, hemisphere, R, G, B, A
    where hemisphere is 'left' | 'right' | 'bilateral' and
    base_name has the L_/R_ prefix stripped.
    """
    with open(label_path, encoding="utf-8-sig") as f:
        lines = [l.strip() for l in f if l.strip()]

    def _hemi(name):
        if   re.match(r"^[Ll]_", name): return "left",     re.sub(r"^[Ll]_", "", name)
        elif re.match(r"^[Rr]_", name): return "right",    re.sub(r"^[Rr]_", "", name)
        else:                            return "bilateral", name

    # Detect format: StatsLUT if first non-empty token of first line is a digit
    rows = []
    first_parts = lines[0].split()
    is_statslut = len(first_parts) >= 2 and first_parts[0].isdigit()

    if is_statslut:
        for line in lines:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                lid = int(parts[0])
            except ValueError:
                continue
            name = parts[1]
            if name.lower() in ("unknown", "???"):
                continue
            try:
                r, g, b, a = int(parts[2]), int(parts[3]), int(parts[4]), int(parts[5])
            except (IndexError, ValueError):
                r = g = b = a = 0
            hemi, base = _hemi(name)
            rows.append(dict(label_id=lid, region_name=name, base_name=base,
                             hemisphere=hemi, R=r, G=g, B=b, A=a))
    else:
        # Alternating name / "ID R G B A" lines
        i = 0
        while i < len(lines) - 1:
            name_line  = lines[i]
            color_line = lines[i + 1]
            parts = color_line.split()
            if len(parts) >= 5 and parts[0].isdigit():
                lid = int(parts[0])
                if lid != 0 and name_line not in ("???", "Unknown"):
                    r, g, b, a = int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])
                    hemi, base = _hemi(name_line)
                    rows.append(dict(label_id=lid, region_name=name_line,
                                     base_name=base, hemisphere=hemi,
                                     R=r, G=g, B=b, A=a))
                i += 2
            else:
                i += 1

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df[df["label_id"] != 0].reset_index(drop=True)
    return df


def parse_legend_xlsx(legend_path, sheet="Legend_2023"):
    """
    Parse the EDNiX legend Excel file.

    Returns dict:
        {
          'level1': [region_name, ...],
          'level2': [...],
          'level3': [...],
          'level4': [...],
          'hierarchy': {level1_name: {level2_name: {level3_name: [level4_names]}}},
          'region_to_level': {base_name: int}   # lowest (finest) level for each name
        }

    region_to_level maps each region's base_name (no hemisphere prefix) to the
    hierarchy level at which it first appears (1=coarsest). This is used by
    Stage 0 to decide which dlabel map_number to query for each region.
    """
    try:
        df = pd.read_excel(legend_path, sheet_name=sheet)
    except Exception as e:
        warnings.warn(f"parse_legend_xlsx: cannot read {legend_path}: {e}")
        return {}

    def _valid(v):
        return not (pd.isna(v) or str(v).strip() in ("", "0", "NA"))

    levels = {}
    region_to_level = {}

    for lvl in range(1, 5):
        col = f"NEWLVL{lvl}"
        if col not in df.columns:
            continue
        names = sorted(str(v).strip() for v in df[col].unique() if _valid(v))
        levels[lvl] = names
        for n in names:
            base = re.sub(r"^[lLrR]_", "", n).strip()
            if base not in region_to_level:
                region_to_level[base] = lvl

    # Build hierarchy dict
    hierarchy = {}
    for _, row in df.iterrows():
        l1 = row.get("NEWLVL1")
        if not _valid(l1):
            continue
        l1 = str(l1).strip()
        hierarchy.setdefault(l1, {})
        l2 = row.get("NEWLVL2")
        if not _valid(l2): continue
        l2 = str(l2).strip()
        hierarchy[l1].setdefault(l2, {})
        l3 = row.get("NEWLVL3")
        if not _valid(l3): continue
        l3 = str(l3).strip()
        hierarchy[l1][l2].setdefault(l3, [])
        l4 = row.get("NEWLVL4")
        if _valid(l4):
            l4 = str(l4).strip()
            if l4 not in hierarchy[l1][l2][l3]:
                hierarchy[l1][l2][l3].append(l4)

    return {
        "level1":          levels.get(1, []),
        "level2":          levels.get(2, []),
        "level3":          levels.get(3, []),
        "level4":          levels.get(4, []),
        "hierarchy":       hierarchy,
        "region_to_level": region_to_level,
    }


# Keep old name for backward compatibility
def extract_regions_from_legend(file_path):
    return parse_legend_xlsx(file_path)


def get_atlas_label_path(atlas_library_root, species_path_fragment, atlas_name,
                          prefer_statslut=False):
    label_dir = opj(atlas_library_root, species_path_fragment, "label_code")
    priority = (
        [f"{atlas_name}_StatsLUT.txt", f"{atlas_name}_label.txt"]
        if prefer_statslut
        else [f"{atlas_name}_label.txt", f"{atlas_name}_StatsLUT.txt"]
    )
    for fn in priority:
        c = opj(label_dir, fn)
        if os.path.exists(c):
            return _linux_path(c)
    raise FileNotFoundError(f"No label file for atlas '{atlas_name}' in {label_dir}")


def find_species_path(atlas_library_root, species_name):
    atlas_root = opj(atlas_library_root, "atlas")
    matches = []
    for dirpath, dirnames, _ in os.walk(atlas_root):
        if os.path.basename(dirpath).lower() == species_name.lower():
            matches.append(_linux_path(os.path.relpath(dirpath, atlas_library_root)))
            dirnames.clear()
    if not matches:
        raise FileNotFoundError(f"No species folder '{species_name}' under {atlas_root}")
    if len(matches) > 1:
        raise ValueError(f"Ambiguous species '{species_name}':\n" +
                         "\n".join(f"  {m}" for m in matches))
    return matches[0]


def build_atlas_config(label_path, legend_path=None):
    """
    Build a complete atlas configuration dict from a StatsLUT/label file
    and optionally a legend xlsx.

    Returns:
        {
          'label_df':        DataFrame from parse_label_file(),
          'id_to_info':      {label_id: {region_name, base_name, hemisphere}},
          'name_to_ids':     {base_name: {'L': id_or_None, 'R': id_or_None,
                                          'bilateral': id_or_None}},
          'legend':          dict from parse_legend_xlsx() or {},
          'region_to_level': {base_name: atlas_level_int}  # from legend or default 1
        }

    This is the single source of truth passed to Stage 0 and Stage 1.
    Works for any atlas — EDNiX or custom — as long as the StatsLUT format
    is:  ID  NAME  R  G  B  A
    and names optionally carry L_/R_ or l_/r_ hemisphere prefixes.
    """
    label_df = parse_label_file(label_path)

    id_to_info = {}
    name_to_ids = {}   # base_name -> {'L': id, 'R': id, 'bilateral': id}

    for _, row in label_df.iterrows():
        lid  = int(row["label_id"])
        hemi = row["hemisphere"]       # 'left' | 'right' | 'bilateral'
        base = row["base_name"]
        id_to_info[lid] = {
            "region_name": row["region_name"],
            "base_name":   base,
            "hemisphere":  hemi,
        }
        entry = name_to_ids.setdefault(base, {"L": None, "R": None, "bilateral": None})
        if   hemi == "left":      entry["L"]         = lid
        elif hemi == "right":     entry["R"]         = lid
        else:                     entry["bilateral"] = lid

    legend = parse_legend_xlsx(legend_path) if legend_path else {}
    region_to_level = legend.get("region_to_level", {})
    # Default: any region not in legend gets level 1
    for base in name_to_ids:
        if base not in region_to_level:
            region_to_level[base] = 1

    return {
        "label_df":        label_df,
        "id_to_info":      id_to_info,
        "name_to_ids":     name_to_ids,
        "legend":          legend,
        "region_to_level": region_to_level,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# §4  PATH FINDERS
# ═══════════════════════════════════════════════════════════════════════════════

def extract_surface_paths(bids_dir, list_to_keep=None, list_to_remove=None,
                           surface_pattern="surface.xlsx"):
    list_to_keep   = list_to_keep   or []
    list_to_remove = list_to_remove or []
    pattern = opj(bids_dir, "sub-*", "ses-*", "anat", "native", "surfaces",
                  "Native_resol", surface_pattern)
    subs, sess, paths = [], [], []
    for f in sorted(glob.glob(pattern)):
        f   = _linux_path(f)
        sub = _parse_bids_entity(f, "sub")
        ses = _parse_bids_entity(f, "ses") or "1"
        subs.append(sub); sess.append(ses); paths.append(f)
    keep = _subject_session_filter(subs, sess, list_to_keep, list_to_remove)
    out  = {"subject": [], "session": [], "surface_path": []}
    for sub, ses, p in zip(subs, sess, paths):
        if (sub, ses) in keep:
            out["subject"].append(sub); out["session"].append(ses)
            out["surface_path"].append(p)
    print(f"  [surface]   {len(out['surface_path'])} files  → {bids_dir}")
    return out


def extract_thickness_paths(bids_dir, list_to_keep=None, list_to_remove=None,
                             thickness_pattern="thickness.xlsx"):
    list_to_keep   = list_to_keep   or []
    list_to_remove = list_to_remove or []
    pattern = opj(bids_dir, "sub-*", "ses-*", "anat", "native", "surfaces",
                  "**", thickness_pattern)
    subs, sess, paths = [], [], []
    for f in sorted(glob.glob(pattern, recursive=True)):
        f   = _linux_path(f)
        sub = _parse_bids_entity(f, "sub")
        ses = _parse_bids_entity(f, "ses") or "1"
        subs.append(sub); sess.append(ses); paths.append(f)
    keep = _subject_session_filter(subs, sess, list_to_keep, list_to_remove)
    out  = {"subject": [], "session": [], "thickness_path": []}
    for sub, ses, p in zip(subs, sess, paths):
        if (sub, ses) in keep:
            out["subject"].append(sub); out["session"].append(ses)
            out["thickness_path"].append(p)
    print(f"  [thickness] {len(out['thickness_path'])} files  → {bids_dir}")
    return out


def extract_volume_paths(bids_dir, atlas_name, list_to_keep=None,
                          list_to_remove=None):
    list_to_keep   = list_to_keep   or []
    list_to_remove = list_to_remove or []
    pattern = opj(bids_dir, "sub-*", "ses-*", "anat", "native", "volumes",
                  "labels", f"*_seg-{atlas_name}_dseg.nii.gz")
    subs, sess, paths = [], [], []
    for f in sorted(glob.glob(pattern)):
        f   = _linux_path(f)
        sub = _parse_bids_entity(f, "sub")
        ses = _parse_bids_entity(f, "ses") or "1"
        subs.append(sub); sess.append(ses); paths.append(f)
    keep = _subject_session_filter(subs, sess, list_to_keep, list_to_remove)
    out  = {"subject": [], "session": [], "volume_seg_path": []}
    for sub, ses, p in zip(subs, sess, paths):
        if (sub, ses) in keep:
            out["subject"].append(sub); out["session"].append(ses)
            out["volume_seg_path"].append(p)
    print(f"  [volume]    {len(out['volume_seg_path'])} seg files  → {bids_dir}")
    return out


def extract_qc_paths(bids_dir, list_to_keep=None, list_to_remove=None, fit_kind='correlation',
                     func_qc_suffix="*_QC_values.json",
                     anat_qc_suffix="*_QC_values.json",
                     full_results_suffix="*_full_results.json",
                     atlas_name="EDNIxCSC", atlas_level=2, use_lr=True):
    """
    Extract QC file paths for all subjects/sessions in a BIDS directory.

    Now also finds the best correlation matrix CSV per session
    (EDNIxCSCLR_{atlas_level}_run_*_matrix.csv) for homotopic specificity
    computation in process_qc.  The atlas_name/level/use_lr parameters
    control which matrix file is picked (defaults match the standard analysis).
    """
    list_to_keep   = list_to_keep   or []
    list_to_remove = list_to_remove or []
    sub_ses = set()
    for sub_dir in glob.glob(opj(bids_dir, "sub-*")):
        sub      = _parse_bids_entity(sub_dir, "sub")
        ses_dirs = glob.glob(opj(sub_dir, "ses-*"))
        if ses_dirs:
            for sd in ses_dirs:
                ses = _parse_bids_entity(sd, "ses")
                if sub and ses: sub_ses.add((sub, ses))
        elif sub:
            sub_ses.add((sub, "1"))
    all_subs = [p[0] for p in sub_ses]
    all_ses  = [p[1] for p in sub_ses]
    keep     = _subject_session_filter(all_subs, all_ses, list_to_keep, list_to_remove)

    # Corr matrix pattern
    suffix_lr   = "LR" if use_lr else ""
    corr_pattern = f"{atlas_name}{suffix_lr}_{atlas_level}_run_*_matrix.csv"

    out = {"subject": [], "session": [],
           "func_qc_path": [], "anat_qc_path": [],
           "full_results_path": [], "corr_matrix_path": []}
    for sub, ses in sorted(keep):
        ses_dir_c = opj(bids_dir, f"sub-{sub}", f"ses-{ses}")
        ses_dir   = ses_dir_c if os.path.exists(ses_dir_c) else opj(bids_dir, f"sub-{sub}")
        qc_dir    = opj(ses_dir, "func", "QC")
        fk_re = re.escape(fit_kind)
        corr_dir  = opj(ses_dir, "func", "acpc-func", "Stats", "Correl_matrix", f"{atlas_name}{suffix_lr}", fk_re)
        func_m    = glob.glob(opj(qc_dir, func_qc_suffix))
        full_m    = glob.glob(opj(qc_dir, full_results_suffix))
        anat_m    = glob.glob(opj(ses_dir, "anat",
                                  "QC_anat", anat_qc_suffix))
        # Best corr matrix: prefer run_0, exclude check_fit
        corr_m    = sorted([f for f in glob.glob(opj(corr_dir, corr_pattern))
                            if not any(x in os.path.basename(f)
                                       for x in ("check_fit","flattened","pval","tstat"))])
        out["subject"].append(sub); out["session"].append(ses)
        out["func_qc_path"].append(_linux_path(func_m[0]) if func_m else None)
        out["anat_qc_path"].append(_linux_path(anat_m[0]) if anat_m else None)
        out["full_results_path"].append(_linux_path(full_m[0]) if full_m else None)
        out["corr_matrix_path"].append(_linux_path(corr_m[0]) if corr_m else None)
    n_f = sum(p is not None for p in out["func_qc_path"])
    n_a = sum(p is not None for p in out["anat_qc_path"])
    n_c = sum(p is not None for p in out["corr_matrix_path"])
    print(f"  [QC]        {n_f} func / {n_a} anat / {n_c} corr_matrix  → {bids_dir}")
    return out



def extract_bold_paths(
    bids_dir,
    list_to_keep=None,
    list_to_remove=None,
    bold_pattern="*_space-template_desc-fMRI_residual.nii.gz",
    task="rest",
):
    list_to_keep   = list_to_keep   or []
    list_to_remove = list_to_remove or []

    pattern = opj(bids_dir, "sub-*", "ses-*", "func", "**", bold_pattern)
    files   = glob.glob(pattern, recursive=True)
    if not files:
        pattern = opj(bids_dir, "sub-*", "func", "**", bold_pattern)
        files   = glob.glob(pattern, recursive=True)

    subjects, sessions, runs, paths = [], [], [], []
    for f in sorted(files):
        f   = _linux_path(f)
        sub = _parse_bids_entity(f, "sub")
        ses = _parse_bids_entity(f, "ses") or "1"
        if sub is None:
            continue
        bn = os.path.basename(f)
        if task and f"task-{task}" not in bn and "task-" in bn:
            continue
        run = _parse_run(f)
        subjects.append(sub); sessions.append(ses)
        runs.append(run); paths.append(f)

    keep = _subject_session_filter(subjects, sessions, list_to_keep, list_to_remove)
    out  = {"subject": [], "session": [], "run": [], "bold_path": []}
    for sub, ses, run, p in zip(subjects, sessions, runs, paths):
        if (sub, ses) in keep:
            out["subject"].append(sub); out["session"].append(ses)
            out["run"].append(run); out["bold_path"].append(p)

    print(f"  [bold]      {len(out['bold_path'])} files  ? {bids_dir}")
    return out


def extract_corr_matrix_paths(bids_dir, atlas_name="EDNIxCSC", atlas_level=3, fit_kind='correlation',
                               use_lr=False, list_to_keep=None, list_to_remove=None):
    fk_re = re.escape(fit_kind)
    list_to_keep   = list_to_keep   or []
    list_to_remove = list_to_remove or []
    suffix       = "LR" if use_lr else ""
    pattern_name = f"{atlas_name}{suffix}_{atlas_level}_run_*_matrix.csv"
    records = []
    for f in sorted(glob.glob(
        opj(bids_dir, "sub-*", "ses-*", "func", "acpc-func",
            "Stats", "Correl_matrix", f"{atlas_name}{suffix}", fk_re, pattern_name)
    )):
        bn = os.path.basename(f)
        if any(x in bn for x in ("check_fit", "flattened", "pval", "tstat")):
            continue
        f   = _linux_path(f)
        sub = _parse_bids_entity(f, "sub")
        ses = _parse_bids_entity(f, "ses") or "1"
        if not sub: continue
        if list_to_keep   and (sub, ses) not in list_to_keep:   continue
        if list_to_remove and (sub, ses) in list_to_remove:     continue
        m   = re.search(r"run_(\d+)_matrix", bn)
        run = int(m.group(1)) if m else 0
        records.append(dict(subject=sub, session=ses, run=run, path=f))
    print(f"  [corr_matrix] {len(records)} files  → {bids_dir}")
    return records


def load_corr_matrix(path):
    df  = pd.read_csv(path, index_col=0)
    idx = list(df.index)
    cols = list(df.columns)
    def _looks_like_rois(names):
        return any(not str(n).lstrip("-").replace(".", "").isdigit() for n in names)
    if _looks_like_rois(idx):
        return idx, df.values.astype(float)
    elif _looks_like_rois(cols):
        return cols, df.T.values.astype(float)
    return idx, df.values.astype(float)

# ═══════════════════════════════════════════════════════════════════════════════
# HOMOTOPIC BILATERAL SPECIFICITY
# ═══════════════════════════════════════════════════════════════════════════════

def compute_homotopic_specificity(mat, roi_names,
                                   thresh_homotopic=0.1,
                                   thresh_delta=0.05):
    """
    Atlas-agnostic homotopic bilateral specificity.

    Finds all L_X / R_X pairs present in roi_names, then computes:

      homotopic_mean   = mean r(L_X, R_X)   across all matched pairs X
      cross_mean       = mean r(L_X, R_Y) for all Y ≠ X  (and R_X, L_Y)
      specificity_index = homotopic_mean − cross_mean

    This replaces the primate-centric sensory/mPFC circuit used in
    analyze_specificity().  Works for any species and any atlas that
    uses L_ / R_ name prefixes (EDNIxCSC, Grandjean, NIH-MBM …).

    Classification (same four categories as before):
      Specific   : specificity_index > thresh_delta AND homotopic_mean > thresh_homotopic
      Unspecific : homotopic_mean > thresh_homotopic but specificity_index ≤ thresh_delta
      No         : homotopic_mean ≤ thresh_homotopic
      Spurious   : specificity_index < 0  (cross-hemisphere > homotopic → noise/artefact)

    Parameters
    ----------
    mat          : (N,N) correlation matrix (numpy array)
    roi_names    : list of N ROI name strings
    thresh_homotopic : minimum homotopic_mean to leave the "No" category (default 0.1)
    thresh_delta     : minimum specificity_index for "Specific" (default 0.05)

    Returns
    -------
    dict with keys:
      homotopic_mean, cross_mean, specificity_index, sp_category,
      sp_specific_correlation (= homotopic_mean),
      sp_nonspecific_correlation (= cross_mean),
      n_homotopic_pairs
    """
    n = len(roi_names)
    if mat.shape != (n, n):
        return {"sp_category": "No", "n_homotopic_pairs": 0,
                "homotopic_mean": np.nan, "cross_mean": np.nan,
                "specificity_index": np.nan,
                "sp_specific_correlation": np.nan,
                "sp_nonspecific_correlation": np.nan}

    def _base(name):
        s = re.sub(r"^[LlRr]_", "", str(name))
        return s.lower().strip()

    name_to_idx = {name: i for i, name in enumerate(roi_names)}
    l_names = [n for n in roi_names if n.startswith(("L_", "l_"))]
    r_names = [n for n in roi_names if n.startswith(("R_", "r_"))]

    # Homotopic pairs: L_X ↔ R_X (exact base-name match)
    homotopic_pairs = []
    for lname in l_names:
        base  = _base(lname)
        rname = next((r for r in r_names if _base(r) == base), None)
        if rname is not None:
            homotopic_pairs.append((name_to_idx[lname], name_to_idx[rname]))

    if not homotopic_pairs:
        return {"sp_category": "No", "n_homotopic_pairs": 0,
                "homotopic_mean": np.nan, "cross_mean": np.nan,
                "specificity_index": np.nan,
                "sp_specific_correlation": np.nan,
                "sp_nonspecific_correlation": np.nan}

    # Homotopic values
    homo_vals = [float(mat[i, j]) for i, j in homotopic_pairs
                 if not np.isnan(mat[i, j])]
    homotopic_mean = float(np.nanmean(homo_vals)) if homo_vals else np.nan

    # Cross-network: L_X ↔ R_Y where X ≠ Y (and mirror)
    homo_set = set(homotopic_pairs) | {(j, i) for i, j in homotopic_pairs}
    l_idx = [i for i, _ in homotopic_pairs]
    r_idx = [j for _, j in homotopic_pairs]
    cross_vals = []
    for li in l_idx:
        for rj in r_idx:
            if (li, rj) not in homo_set and not np.isnan(mat[li, rj]):
                cross_vals.append(float(mat[li, rj]))
    for ri in r_idx:
        for lj in l_idx:
            if (ri, lj) not in homo_set and not np.isnan(mat[ri, lj]):
                cross_vals.append(float(mat[ri, lj]))

    cross_mean = float(np.nanmean(cross_vals)) if cross_vals else np.nan
    if not np.isnan(homotopic_mean) and not np.isnan(cross_mean):
        specificity_index = homotopic_mean - cross_mean
    else:
        specificity_index = np.nan

    # Classification
    if np.isnan(homotopic_mean) or np.isnan(specificity_index):
        cat = "No"
    elif specificity_index < 0:
        cat = "Spurious"
    elif homotopic_mean <= thresh_homotopic:
        cat = "No"
    elif specificity_index <= thresh_delta:
        cat = "Unspecific"
    else:
        cat = "Specific"

    return {
        "homotopic_mean":             homotopic_mean,
        "cross_mean":                 cross_mean,
        "specificity_index":          float(specificity_index)
                                      if not np.isnan(specificity_index) else np.nan,
        "sp_category":                cat,
        "sp_specific_correlation":    homotopic_mean,
        "sp_nonspecific_correlation": cross_mean,
        "n_homotopic_pairs":          len(homotopic_pairs),
    }

# ═══════════════════════════════════════════════════════════════════════════════
# §5  STAGE 0 — BULK wb_command EXTRACTION
# ═══════════════════════════════════════════════════════════════════════════════
#
# DESIGN
# ──────
# Old approach: one cifti-label-to-roi call per region per hemisphere per level
#   → 4 × N_regions × N_levels wb_command calls per subject.
#
# New approach: per atlas level per subject:
#   1. Read label table from dlabel (cifti-label-export-table)  → 1 call, cached
#   2. cifti-separate the full dlabel into LEFT + RIGHT dense scalars    → 2 calls
#   3. Load those scalars with nibabel (fast, in Python)
#   4. For each region: numpy mask (data == label_id) → area/thickness
#      → 0 additional wb_command calls
#   5. surface-vertex-areas on pial surfaces (one per hemisphere)        → 2 calls
#
# Total: 5 wb_command calls per level per subject (all regions simultaneously)
# vs 4 × N_regions calls in the old design.
#
# Label IDs are always taken from the dlabel itself (cifti-label-export-table),
# never from the StatsLUT, so they are always correct regardless of atlas version.
# The StatsLUT / atlas_config is used only for:
#   • hemisphere inference per region
#   • level assignment (from legend xlsx)
#   • volume extraction
#

# ── Per-process cache: (seg_file, map_number) → [(dlabel_name, label_id), ...] ──
_DLABEL_LABEL_CACHE = {}


def _dlabel_get_labels(sing_wb, seg_file, map_number=1):
    """
    Export label table from dlabel for one map.
    Returns list of (dlabel_name_str, label_id_int) — background excluded.
    Cached per (seg_file, map_number).
    """
    key = (seg_file, map_number)
    if key in _DLABEL_LABEL_CACHE:
        return _DLABEL_LABEL_CACHE[key]

    seg_hash = hashlib.md5(seg_file.encode()).hexdigest()[:8]
    tmp = f"/tmp/ednix_ltable_{seg_hash}_m{map_number}.txt"
    _spco(f'{sing_wb} wb_command -cifti-label-export-table '
          f'"{seg_file}" {map_number} "{tmp}"')

    results = []
    if os.path.exists(tmp):
        with open(tmp) as f:
            lines = [l.strip() for l in f if l.strip()]
        os.remove(tmp)
        i = 0
        while i < len(lines):
            name_line = lines[i]
            if name_line in ("???", "Unknown") or name_line.startswith("0 "):
                i += 1; continue
            if i + 1 < len(lines):
                parts = lines[i + 1].split()
                if parts and parts[0].isdigit():
                    lid = int(parts[0])
                    if lid != 0:
                        results.append((name_line, lid))
                    i += 2; continue
            i += 1
        print(f"  [dlabel_labels] map={map_number} → {len(results)} labels "
              f"({os.path.basename(seg_file)})")
    else:
        warnings.warn(f"  [dlabel_labels] export-table failed: {seg_file} map {map_number}")

    _DLABEL_LABEL_CACHE[key] = results
    return results


def _find_surf_gii(native_dir, sub, hemi_char, surf_type="pial"):
    h = hemi_char.lower()
    candidates = [
        opj(native_dir, f"{sub}.{h}.{surf_type}.native.surf.gii"),
        opj(native_dir, f"{sub}.{h}.{surf_type}.surf.gii"),
        opj(native_dir, f"sub-{sub}.{h}.{surf_type}.native.surf.gii"),
        opj(native_dir, f"sub-{sub}.{h}.{surf_type}.surf.gii"),
    ]
    for c in candidates:
        if os.path.exists(c): return c
    matches = glob.glob(opj(native_dir, f"*.{h}.{surf_type}*.surf.gii"))
    return sorted(matches)[0] if matches else None


def _get_seg_file(native_dir, sub, atlas_name="EDNIxCSC"):
    for prefix in [sub, f"sub-{sub}"]:
        seg_files = glob.glob(opj(native_dir, f"{prefix}.{atlas_name}*.dlabel.nii"))
        if seg_files:
            # Prefer non-LR version if both exist
            non_lr = [s for s in seg_files if "LR" not in os.path.basename(s)]
            return non_lr[0] if non_lr else seg_files[0]
    warnings.warn(f"No dlabel segmentation for {sub} in {native_dir}")
    return None


def _separate_dense_to_metric(sing_wb, cifti_f, hemi_side, out_f, timeout=_TIMEOUT_HEAVY):
    """
    wb_command -cifti-separate cifti_f COLUMN -metric CORTEX_{LEFT|RIGHT} out_f
    Returns True on success.
    """
    if os.path.exists(out_f):
        os.remove(out_f)
    r = _spco_timed(
        f'{sing_wb} wb_command -cifti-separate "{cifti_f}" COLUMN'
        f' -metric CORTEX_{hemi_side} "{out_f}"',
        timeout=timeout)
    return r is not None and os.path.exists(out_f)


def _surface_vertex_areas(sing_wb, surf_f, out_f):
    """
    wb_command -surface-vertex-areas surf_f out_f
    Returns True on success.
    """
    if os.path.exists(out_f):
        os.remove(out_f)
    r = _spco_timed(f'{sing_wb} wb_command -surface-vertex-areas "{surf_f}" "{out_f}"')
    return r is not None and os.path.exists(out_f)


def _bulk_extract_surface(sing_wb, seg_file, map_number, native_dir, sub,
                           dlabel_labels, regions_of_interest, atlas_config,
                           roi_dir):
    """
    Extract surface area (mm²) for all requested regions at one atlas level.

    Strategy
    --------
    1. cifti-separate the dlabel (map_number) into L and R dense metric files.
    2. Load the dense metric data with nibabel → integer label arrays.
    3. Compute vertex areas for L and R pial surfaces.
    4. For each requested region, mask by label_id and sum areas.

    Returns dict {base_name: area_mm2}.
    """
    os.makedirs(roi_dir, exist_ok=True)
    safe_map = f"m{map_number}"

    # ── Separate dlabel map to L and R metric ─────────────────────────────────
    label_l = opj(roi_dir, f"dense_labels_{safe_map}_L.func.gii")
    label_r = opj(roi_dir, f"dense_labels_{safe_map}_R.func.gii")
    ok_l = _separate_dense_to_metric(sing_wb, seg_file, "LEFT",  label_l)
    ok_r = _separate_dense_to_metric(sing_wb, seg_file, "RIGHT", label_r)

    if not ok_l and not ok_r:
        warnings.warn(f"  [surface_bulk] {sub}: cifti-separate failed both hemis map {map_number}")
        return {}

    # ── Vertex areas ──────────────────────────────────────────────────────────
    area_l = area_r = None
    data_l = data_r = None

    if ok_l:
        surf_l  = _find_surf_gii(native_dir, sub, "l", "pial")
        if surf_l:
            area_f_l = opj(roi_dir, f"pial_area_L.func.gii")
            if _surface_vertex_areas(sing_wb, surf_l, area_f_l):
                area_l = _read_gifti_data(area_f_l)
                data_l = _read_gifti_data(label_l).astype(np.int32)
    if ok_r:
        surf_r  = _find_surf_gii(native_dir, sub, "r", "pial")
        if surf_r:
            area_f_r = opj(roi_dir, f"pial_area_R.func.gii")
            if _surface_vertex_areas(sing_wb, surf_r, area_f_r):
                area_r = _read_gifti_data(area_f_r)
                data_r = _read_gifti_data(label_r).astype(np.int32)

    if data_l is None and data_r is None:
        warnings.warn(f"  [surface_bulk] {sub}: no vertex data map {map_number}")
        return {}

    # ── Build label_id lookup from dlabel labels ───────────────────────────────
    # dlabel_labels: [(dlabel_name, label_id), ...]
    # We need to find for each requested base_name the L and R label_ids
    # by matching against atlas_config (which knows hemisphere per base_name)
    # If atlas_config not available, infer from name prefix.
    name_to_ids_dlabel = _build_dlabel_name_map(dlabel_labels, atlas_config)

    # ── Filter to requested regions ───────────────────────────────────────────
    all_base_names = list(name_to_ids_dlabel.keys())
    if regions_of_interest:
        all_base_names = _filter_regions(all_base_names, regions_of_interest)

    results = {}
    for base in all_base_names:
        ids = name_to_ids_dlabel[base]
        total = 0.0
        found = False
        if data_l is not None and ids.get("L") is not None:
            mask = (data_l == ids["L"])
            if mask.any():
                total += float(np.sum(area_l[mask]))
                found = True
        if data_r is not None and ids.get("R") is not None:
            mask = (data_r == ids["R"])
            if mask.any():
                total += float(np.sum(area_r[mask]))
                found = True
        if found:
            results[base] = total

    return results


def _bulk_extract_thickness(sing_wb, seg_file, dscalar_f, map_number,
                              native_dir, sub, dlabel_labels,
                              regions_of_interest, atlas_config, roi_dir):
    """
    Extract mean cortical thickness (mm) for all requested regions at one level.

    Strategy: same as _bulk_extract_surface but reads thickness dscalar instead
    of vertex areas. Mean taken over all vertices in the ROI mask.
    """
    os.makedirs(roi_dir, exist_ok=True)
    safe_map = f"m{map_number}"

    label_l = opj(roi_dir, f"dense_labels_{safe_map}_L.func.gii")
    label_r = opj(roi_dir, f"dense_labels_{safe_map}_R.func.gii")
    ok_l = _separate_dense_to_metric(sing_wb, seg_file, "LEFT",  label_l)
    ok_r = _separate_dense_to_metric(sing_wb, seg_file, "RIGHT", label_r)

    if not ok_l and not ok_r:
        warnings.warn(f"  [thick_bulk] {sub}: cifti-separate failed map {map_number}")
        return {}

    thick_l = opj(roi_dir, f"thickness_{safe_map}_L.func.gii")
    thick_r = opj(roi_dir, f"thickness_{safe_map}_R.func.gii")
    ok_tl = ok_l and _separate_dense_to_metric(sing_wb, dscalar_f, "LEFT",  thick_l, _TIMEOUT_HEAVY)
    ok_tr = ok_r and _separate_dense_to_metric(sing_wb, dscalar_f, "RIGHT", thick_r, _TIMEOUT_HEAVY)

    data_l = data_r = thick_data_l = thick_data_r = None
    if ok_l and os.path.exists(label_l):
        data_l = _read_gifti_data(label_l).astype(np.int32)
    if ok_r and os.path.exists(label_r):
        data_r = _read_gifti_data(label_r).astype(np.int32)
    if ok_tl and os.path.exists(thick_l):
        thick_data_l = _read_gifti_data(thick_l).astype(np.float32)
    if ok_tr and os.path.exists(thick_r):
        thick_data_r = _read_gifti_data(thick_r).astype(np.float32)

    if data_l is None and data_r is None:
        return {}

    name_to_ids_dlabel = _build_dlabel_name_map(dlabel_labels, atlas_config)
    all_base_names = list(name_to_ids_dlabel.keys())
    if regions_of_interest:
        all_base_names = _filter_regions(all_base_names, regions_of_interest)

    results = {}
    for base in all_base_names:
        ids = name_to_ids_dlabel[base]
        all_verts = []
        if data_l is not None and thick_data_l is not None and ids.get("L") is not None:
            mask = (data_l == ids["L"])
            if mask.any():
                all_verts.append(thick_data_l[mask])
        if data_r is not None and thick_data_r is not None and ids.get("R") is not None:
            mask = (data_r == ids["R"])
            if mask.any():
                all_verts.append(thick_data_r[mask])
        if all_verts:
            combined = np.concatenate(all_verts)
            combined = combined[np.isfinite(combined) & (combined > 0)]
            if len(combined):
                results[base] = float(np.mean(combined))

    return results


def _build_dlabel_name_map(dlabel_labels, atlas_config):
    """
    Map base_name → {'L': dlabel_id_L, 'R': dlabel_id_R} using:
      - dlabel_labels: [(dlabel_name, dlabel_id), ...] from cifti-label-export-table
      - atlas_config: provides hemisphere inference if available

    Works even without atlas_config (falls back to name prefix).
    The dlabel IDs are ALWAYS from the dlabel file itself, never from StatsLUT.
    """
    name_to_ids = {}   # base_name -> {'L': id, 'R': id}

    ac_name_to_hemi = {}
    if atlas_config:
        for base, entry in atlas_config.get("name_to_ids", {}).items():
            if entry["L"] is not None:   ac_name_to_hemi[base] = {"L", "bilateral"}
            if entry["R"] is not None:   ac_name_to_hemi[base] = {"R", "bilateral"}
            if entry["bilateral"] is not None: ac_name_to_hemi[base] = {"bilateral"}

    for dlabel_name, dlabel_id in dlabel_labels:
        # Infer hemisphere from name prefix
        rn = dlabel_name
        if   re.match(r"^[Ll]_", rn): hemi_key = "L"; base = re.sub(r"^[Ll]_", "", rn)
        elif re.match(r"^[Rr]_", rn): hemi_key = "R"; base = re.sub(r"^[Rr]_", "", rn)
        else:                           hemi_key = "bilateral"; base = rn

        entry = name_to_ids.setdefault(base, {"L": None, "R": None})
        if   hemi_key == "L":          entry["L"] = dlabel_id
        elif hemi_key == "R":          entry["R"] = dlabel_id
        else:
            # Bilateral: assign same id to both sides so mask works on either
            if entry["L"] is None:  entry["L"] = dlabel_id
            if entry["R"] is None:  entry["R"] = dlabel_id

    return name_to_ids


def _wb_compute_thickness_cifti(sing_wb, sub, ses_dir, native_dir, overwrite=True):
    """Compute per-vertex cortical thickness dscalar. Returns path or None."""
    dscalar_out = opj(native_dir, f"{sub}.thicknesswb.dscalar.nii")
    if os.path.exists(dscalar_out) and not overwrite:
        return dscalar_out

    palette = (" MODE_AUTO_SCALE_PERCENTAGE -pos-percent 4 96"
               " -interpolate true -palette-name videen_style"
               " -disp-pos true -disp-neg false -disp-zero false")
    shape_files = {}

    for h in ("l", "r"):
        pial_f  = _find_surf_gii(native_dir, sub, h, "pial")
        white_f = _find_surf_gii(native_dir, sub, h, "white")
        if pial_f is None or white_f is None:
            warnings.warn(f"  [thickness_wb] {sub}: surface missing hemi={h}")
            return None

        thick_shape = opj(native_dir, f"{sub}.{h}.thicknesswb.shape.gii")
        if os.path.exists(thick_shape): os.remove(thick_shape)

        if (_spco_timed(
                f'{sing_wb} wb_command -surface-to-surface-3d-distance'
                f' "{pial_f}" "{white_f}" "{thick_shape}"',
                _TIMEOUT_HEAVY) is None or not os.path.exists(thick_shape)):
            warnings.warn(f"  [thickness_wb] {sub}: 3d-distance failed hemi={h}")
            return None

        roi_f = opj(native_dir, f"{sub}.{h}.roi.shape.gii")
        if os.path.exists(roi_f):
            tmp = thick_shape + ".tmp.shape.gii"
            if os.path.exists(tmp): os.remove(tmp)
            if (_spco_timed(
                    f'{sing_wb} wb_command -metric-math "roi * thickness"'
                    f' "{tmp}" -var roi "{roi_f}" -var thickness "{thick_shape}"')
                    is not None and os.path.exists(tmp)):
                os.replace(tmp, thick_shape)
        _spco_timed(f'{sing_wb} wb_command -metric-palette "{thick_shape}"{palette}')
        shape_files[h] = thick_shape

    if len(shape_files) < 2:
        return None

    if os.path.exists(dscalar_out): os.remove(dscalar_out)
    roi_l = opj(native_dir, f"{sub}.l.roi.shape.gii")
    roi_r = opj(native_dir, f"{sub}.r.roi.shape.gii")
    cmd4  = (f'{sing_wb} wb_command -cifti-create-dense-scalar "{dscalar_out}"'
             f' -left-metric "{shape_files["l"]}"')
    if os.path.exists(roi_l): cmd4 += f' -roi-left "{roi_l}"'
    cmd4 += f' -right-metric "{shape_files["r"]}"'
    if os.path.exists(roi_r): cmd4 += f' -roi-right "{roi_r}"'

    if _spco_timed(cmd4, _TIMEOUT_HEAVY) is None or not os.path.exists(dscalar_out):
        return None
    _spco_timed(f'{sing_wb} wb_command -set-map-names "{dscalar_out}"'
                f' -map 1 {sub}_thickness_from_wb')
    return dscalar_out


def _save_morphometry_xlsx(results_by_level, ses_dir, modality):
    """
    Write surface.xlsx or thickness.xlsx.
    results_by_level: {atlas_level: {base_name: value}}
    Columns: one per atlas level (Surface_Area_lvl1 … or Thickness_mm_lvl1 …)
    """
    if not results_by_level:
        return
    out_dir = opj(ses_dir, "anat", "native", "surfaces", "Native_resol")
    os.makedirs(out_dir, exist_ok=True)
    fname  = "surface.xlsx" if modality == "surface" else "thickness.xlsx"
    col_pfx = "Surface_Area" if modality == "surface" else "Thickness_mm"
    out    = opj(out_dir, fname)
    if os.path.exists(out): os.remove(out)

    frames = {}
    for lvl in sorted(results_by_level):
        col = f"{col_pfx}_lvl{lvl}"
        frames[col] = pd.Series(results_by_level[lvl], name=col)

    if not frames:
        return
    df = pd.concat(frames.values(), axis=1)
    df.index.name = "Region"
    df.to_excel(out)
    non_nan = {c: int(df[c].notna().sum()) for c in df.columns}
    print(f"  [{modality} xlsx] {out}  ({len(df)} regions, non-NaN: {non_nan})")


def _build_records(results_by_level, sub, ses_label, value_key):
    """Convert {atlas_level: {base_name: value}} → list of record dicts."""
    records = []
    for atlas_level, results in results_by_level.items():
        for base_name, val in results.items():
            records.append(dict(
                subject=sub, session=ses_label,
                region=base_name,
                region_full=base_name,
                hemisphere="bilateral",    # bulk extraction returns base names
                atlas_level=atlas_level,
                **{value_key: val},
            ))
    return records


def _process_one_subject(sing_wb, sub, ses_label, ses_dir, native_dir,
                          atlas_name, n_atlas_levels, regions_of_interest,
                          overwrite_surface, overwrite_thickness,
                          atlas_config=None):
    """
    Stage 0 worker: bulk extract surface area + thickness for one subject.

    For each atlas level:
      • 2 cifti-separate calls for the dlabel → label metric files
      • 2 surface-vertex-areas calls (once, reused for thickness too)
      • All regions extracted via numpy masking (0 extra wb_command calls)
      • 2 more cifti-separate calls for the thickness dscalar

    Total: ~5 wb_command calls per level (all regions at once).
    """
    try:
        seg_file = _get_seg_file(native_dir, sub, atlas_name)
        if seg_file is None:
            return {}, {}, {"subject": sub, "session": ses_label,
                            "reason": "no dlabel seg file"}

        surf_by_level  = {}
        thick_by_level = {}
        roi_dir = opj(native_dir, "ROIs", "bulk")

        # Compute thickness dscalar once (shared across all levels)
        dscalar_f = None
        if overwrite_thickness or not os.path.exists(
                opj(native_dir, f"{sub}.thicknesswb.dscalar.nii")):
            dscalar_f = _wb_compute_thickness_cifti(
                sing_wb, sub, ses_dir, native_dir, overwrite=overwrite_thickness)
        else:
            dscalar_f = opj(native_dir, f"{sub}.thicknesswb.dscalar.nii")

        for atlas_level in range(1, n_atlas_levels + 1):
            # Get label table for this map from the dlabel
            dlabel_labels = _dlabel_get_labels(sing_wb, seg_file, atlas_level)
            if not dlabel_labels:
                print(f"  [stage0] {sub} ses-{ses_label} map {atlas_level}: "
                      f"no labels, skipping")
                continue

            lvl_roi_dir = opj(roi_dir, f"lvl{atlas_level}")

            # Surface extraction
            if overwrite_surface or not os.path.exists(
                    opj(native_dir, "Native_resol", "surface.xlsx")):
                surf = _bulk_extract_surface(
                    sing_wb, seg_file, atlas_level, native_dir, sub,
                    dlabel_labels, regions_of_interest, atlas_config, lvl_roi_dir)
                if surf:
                    surf_by_level[atlas_level] = surf

            # Thickness extraction
            if dscalar_f and (overwrite_thickness or not os.path.exists(
                    opj(native_dir, "Native_resol", "thickness.xlsx"))):
                thick = _bulk_extract_thickness(
                    sing_wb, seg_file, dscalar_f, atlas_level, native_dir, sub,
                    dlabel_labels, regions_of_interest, atlas_config, lvl_roi_dir)
                if thick:
                    thick_by_level[atlas_level] = thick

        print(f"  [stage0] {sub} ses-{ses_label}: "
              f"surface levels={sorted(surf_by_level)}, "
              f"thickness levels={sorted(thick_by_level)}")
        return surf_by_level, thick_by_level, None

    except Exception as exc:
        reason = f"{type(exc).__name__}: {exc}"
        warnings.warn(f"  [stage0] {sub} ses-{ses_label} FAILED → {reason}\n"
                      + traceback.format_exc())
        return {}, {}, {"subject": sub, "session": ses_label, "reason": reason}


def run_wb_extraction(sing_wb, all_bids, regions_to_process=None,
                       overwrite_thickness=True, overwrite_surface=True,
                       n_atlas_levels=N_ATLAS_LEVELS, atlas_name='EDNiX',
                       atlas_label_paths=None, legend_paths=None):
    """
    STAGE 0 — top-level entry point.

    Parameters
    ----------
    sing_wb             : str  — Singularity wb_command prefix.
    all_bids            : list of (species, bids_dir).
    regions_to_process  : list of region name patterns (fuzzy), or None for all.
    overwrite_thickness : bool — recompute even if dscalar exists.
    overwrite_surface   : bool — recompute even if xlsx exists.
    n_atlas_levels      : int  — number of dlabel maps to process (default 4).
    atlas_label_paths   : dict {species: label_file_path} — StatsLUT or label.txt.
    legend_paths        : dict {species: legend_xlsx_path} — for hierarchy levels.

    Side effects: writes surface.xlsx / thickness.xlsx per subject/session.
    """
    atlas_label_paths = atlas_label_paths or {}
    legend_paths      = legend_paths      or {}

    for species, bids_dir in all_bids:
        bids_lbl = os.path.basename(bids_dir.rstrip("/"))
        print(f"\n  [STAGE 0] {species} → {bids_lbl}")

        # Build atlas_config once per species
        label_path  = atlas_label_paths.get(species)
        legend_path = legend_paths.get(species)
        atlas_config = None
        if label_path and os.path.exists(label_path):
            try:
                atlas_config = build_atlas_config(label_path, legend_path)
                print(f"  [atlas_config] {species}: "
                      f"{len(atlas_config['name_to_ids'])} regions loaded")
            except Exception as e:
                warnings.warn(f"  [atlas_config] {species}: {e}")

        subjects = list(_iter_subjects(bids_dir, [], []))
        def _worker(sub, ses_label, ses_dir, native_dir):
            surf_by_level, thick_by_level, failure = _process_one_subject(
                sing_wb, sub, ses_label, ses_dir, native_dir,
                atlas_name, n_atlas_levels, regions_to_process,
                overwrite_surface, overwrite_thickness, atlas_config)
            if surf_by_level:
                _save_morphometry_xlsx(surf_by_level, ses_dir, "surface")
            if thick_by_level:
                _save_morphometry_xlsx(thick_by_level, ses_dir, "thickness")
            return failure

        if _JOBLIB and len(subjects) > 1:
            failures = Parallel(n_jobs=-1, backend="loky", verbose=5)(
                delayed(_worker)(sub, ses_label, ses_dir, native_dir)
                for sub, ses_label, ses_dir, native_dir in subjects)
        else:
            failures = [_worker(sub, ses_label, ses_dir, native_dir)
                        for sub, ses_label, ses_dir, native_dir in subjects]

        _save_failed([f for f in failures if f], out_dir=bids_dir,
                     modality="stage0")


def _infer_atlas_name(bids_dir):
    """Guess atlas name from dlabel files found in the BIDS directory."""
    pattern = opj(bids_dir, "sub-*", "ses-*", "anat", "native", "surfaces",
                  "Native_resol", "*.dlabel.nii")
    files = glob.glob(pattern)
    if not files:
        return "EDNIxCSC"
    bn = os.path.basename(files[0])
    # e.g. 300700.EDNIxCSC.dlabel.nii → EDNIxCSC
    m = re.search(r"\.([^.]+)\.dlabel\.nii$", bn)
    return m.group(1) if m else "EDNIxCSC"


# Backward-compatible wrappers
def extract_and_process_surfaces_from_rois(*args, **kwargs):
    warnings.warn("extract_and_process_surfaces_from_rois is deprecated in v7. "
                  "Use run_wb_extraction with overwrite_surface=False.", DeprecationWarning)


def extract_and_process_thickness_from_rois(*args, **kwargs):
    warnings.warn("extract_and_process_thickness_from_rois is deprecated in v7. "
                  "Use run_wb_extraction with overwrite_thickness=False.", DeprecationWarning)


# ═══════════════════════════════════════════════════════════════════════════════
# §6  STAGE 1 — READ xlsx / NIfTI / JSON → TIDY DATAFRAMES
# ═══════════════════════════════════════════════════════════════════════════════

def _read_one_surface(sub, ses, path, regions_of_interest):
    records = []
    if not path or not os.path.exists(path):
        return records
    try:
        df = pd.read_excel(path, index_col=0)
        if df.empty:
            warnings.warn(f"surface.xlsx is empty: {path}")
            return records
        surf_cols = sorted(
            [c for c in df.columns if re.search(r"lvl\d+", c, re.I)
             and ("surface" in c.lower() or "area" in c.lower())],
            key=lambda c: int(re.search(r"lvl(\d+)", c, re.I).group(1)))
        if not surf_cols:
            surf_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
        seen = {}
        for col in surf_cols:
            lvl_m = re.search(r"lvl(\d+)", col, re.I)
            lvl   = int(lvl_m.group(1)) if lvl_m else 1
            for region in df.index:
                rs  = str(region)
                val = df.loc[region, col]
                if pd.isna(val) or rs in seen:
                    continue
                if regions_of_interest and not any(
                        _names_match(rs, p) for p in regions_of_interest):
                    continue
                if   rs.lower().startswith("l_"): hemi, base = "left",     rs[2:]
                elif rs.lower().startswith("r_"): hemi, base = "right",    rs[2:]
                else:                              hemi, base = "bilateral", rs
                seen[rs] = dict(subject=sub, session=ses,
                                region_full=rs, region=base,
                                hemisphere=hemi, atlas_level=lvl,
                                surface_area_mm2=float(val))
        records = list(seen.values())
    except Exception as e:
        warnings.warn(f"Error reading surface {path}: {e}")
    return records


def _read_one_thickness(sub, ses, path, regions_of_interest):
    records = []
    if not path or not os.path.exists(path):
        return records
    try:
        df = pd.read_excel(path, index_col=0)
        if df.empty:
            warnings.warn(f"thickness.xlsx is empty: {path}")
            return records
        thick_cols = sorted(
            [c for c in df.columns if re.search(r"lvl\d+", c, re.I)
             and ("thick" in c.lower() or "mean" in c.lower())],
            key=lambda c: int(re.search(r"lvl(\d+)", c, re.I).group(1)))
        if not thick_cols:
            thick_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
        seen = {}
        for col in thick_cols:
            lvl_m = re.search(r"lvl(\d+)", col, re.I)
            lvl   = int(lvl_m.group(1)) if lvl_m else 1
            for region in df.index:
                rs  = str(region)
                val = df.loc[region, col]
                if pd.isna(val) or rs in seen:
                    continue
                if regions_of_interest and not any(
                        _names_match(rs, p) for p in regions_of_interest):
                    continue
                if   rs.lower().startswith("l_"): hemi, base = "left",     rs[2:]
                elif rs.lower().startswith("r_"): hemi, base = "right",    rs[2:]
                else:                              hemi, base = "bilateral", rs
                seen[rs] = dict(subject=sub, session=ses,
                                region_full=rs, region=base,
                                hemisphere=hemi, atlas_level=lvl,
                                thickness_mm=float(val))
        records = list(seen.values())
    except Exception as e:
        warnings.warn(f"Error reading thickness {path}: {e}")
    return records


def process_surfaces(surface_paths_dict, regions_of_interest=None, **kwargs):
    triples = list(zip(surface_paths_dict["subject"],
                       surface_paths_dict["session"],
                       surface_paths_dict["surface_path"]))
    if _JOBLIB:
        results = Parallel(n_jobs=-1)(
            delayed(_read_one_surface)(s, e, p, regions_of_interest)
            for s, e, p in triples)
    else:
        results = [_read_one_surface(s, e, p, regions_of_interest)
                   for s, e, p in triples]
    records = [r for sub in results for r in sub]
    df = pd.DataFrame(records)
    if df.empty:
        warnings.warn("process_surfaces: empty DataFrame")
    return df


def process_thickness(thickness_paths_dict, regions_of_interest=None):
    triples = list(zip(thickness_paths_dict["subject"],
                       thickness_paths_dict["session"],
                       thickness_paths_dict["thickness_path"]))
    if _JOBLIB:
        results = Parallel(n_jobs=-1)(
            delayed(_read_one_thickness)(s, e, p, regions_of_interest)
            for s, e, p in triples)
    else:
        results = [_read_one_thickness(s, e, p, regions_of_interest)
                   for s, e, p in triples]
    records = [r for sub in results for r in sub]
    df = pd.DataFrame(records)
    if df.empty:
        warnings.warn("process_thickness: empty DataFrame")
    return df


def process_volumes(volume_paths_dict, label_path, regions_of_interest=None,
                    voxel_volume_mm3=None):
    try:
        import nibabel as nib
    except ImportError:
        raise ImportError("nibabel is required")
    label_df     = parse_label_file(label_path)
    label_lookup = {int(r["label_id"]): r for _, r in label_df.iterrows()}
    records = []
    for sub, ses, seg_path in zip(volume_paths_dict["subject"],
                                   volume_paths_dict["session"],
                                   volume_paths_dict["volume_seg_path"]):
        if not seg_path or not os.path.exists(seg_path):
            warnings.warn(f"Seg not found: {seg_path}"); continue
        try:
            img    = nib.load(seg_path)
            data   = np.asarray(img.dataobj, dtype=np.int32)
            vox_mm = (float(voxel_volume_mm3) if voxel_volume_mm3
                      else float(np.prod(img.header.get_zooms()[:3])))
            ids, counts = np.unique(data[data > 0], return_counts=True)
            for lid, cnt in zip(ids, counts):
                if lid not in label_lookup: continue
                row  = label_lookup[lid]
                base = row["base_name"]
                if regions_of_interest and not any(
                        _names_match(base, p) for p in regions_of_interest):
                    continue
                records.append(dict(subject=sub, session=ses,
                                    label_id=int(lid),
                                    region_name=row["region_name"],
                                    region=base,
                                    hemisphere=row["hemisphere"],
                                    voxel_count=int(cnt),
                                    volume_mm3=float(cnt * vox_mm)))
        except Exception as e:
            warnings.warn(f"Error processing {seg_path}: {e}")
    return pd.DataFrame(records)


def process_volumes_and_save(bids_dir, atlas_name, label_path,
                              list_to_keep=None, list_to_remove=None,
                              regions_of_interest=None, overwrite=False):
    vol_paths = extract_volume_paths(bids_dir, atlas_name, list_to_keep, list_to_remove)
    df_all    = process_volumes(vol_paths, label_path, regions_of_interest)
    for (sub, ses), grp in df_all.groupby(["subject", "session"]):
        ses_dir = opj(bids_dir, f"sub-{sub}", f"ses-{ses}")
        if not os.path.exists(ses_dir):
            ses_dir = opj(bids_dir, f"sub-{sub}")
        out_dir  = opj(ses_dir, "anat", "native", "volumes", "labels")
        os.makedirs(out_dir, exist_ok=True)
        out_xlsx = opj(out_dir, f"volumes_{atlas_name}.xlsx")
        if not os.path.exists(out_xlsx) or overwrite:
            grp.set_index("region_name")[
                ["region", "hemisphere", "voxel_count", "volume_mm3"]
            ].to_excel(out_xlsx)
    return df_all


def _load_json_with_nan(path):
    with open(path) as f:
        raw = f.read()
    cleaned = re.sub(r"\bNaN\b", "null", raw)
    data    = json.loads(cleaned)
    def _restore(v):
        if v is None:           return np.nan
        if isinstance(v, dict): return {k: _restore(vv) for k, vv in v.items()}
        if isinstance(v, list): return [_restore(i) for i in v]
        return v
    return _restore(data)


def _flatten(d, prefix="", sep="_"):
    out = {}
    for k, v in d.items():
        key = f"{prefix}{sep}{k}" if prefix else k
        if isinstance(v, dict):
            out.update(_flatten(v, key, sep))
        elif isinstance(v, list):
            if len(v) == 1:
                out[key] = v[0]
            else:
                for i, item in enumerate(v):
                    if isinstance(item, dict): out.update(_flatten(item, f"{key}_{i}", sep))
                    else: out[f"{key}_{i}"] = item
        else:
            out[key] = v
    return out


def _parse_full_results(path):
    """
    Parse full_results.json for network metrics and hemisphere stats.
    Specificity columns (sp_category, sp_specific_correlation, etc.) are
    intentionally NOT populated here — they are computed from the correlation
    matrix directly by process_qc() using compute_homotopic_specificity().
    The JSON-based specificity (sensory/mPFC circuit) is kept only as
    sp_category_legacy for reference.
    """
    try:
        data = _load_json_with_nan(path)
    except Exception as e:
        warnings.warn(f"full_results parse error {path}: {e}"); return {}
    out = {}
    nm  = data.get("network_metrics", {})
    for jk, ck in [("Mean_Correlation","net_mean_correlation"),
                   ("Std_Correlation","net_std_correlation"),
                   ("Top_Eigenvalue","net_top_eigenvalue"),
                   ("Eigenvalue_Ratio","net_eigenvalue_ratio"),
                   ("Silhouette_Score","net_silhouette"),
                   ("Davies_Bouldin","net_davies_bouldin")]:
        if jk in nm: out[ck] = nm[jk]
    hr = data.get("hemisphere_results", {})
    for k in ("intra_left_mean","intra_right_mean","inter_mean",
              "p_intra_vs_inter","p_left_vs_right"):
        if k in hr: out[f"net_{k}"] = hr[k]
    if "intra_left_mean" in hr and "intra_right_mean" in hr:
        out["net_intra_mean"] = (hr["intra_left_mean"] + hr["intra_right_mean"]) / 2
    # Legacy circuit-based specificity (primate-centric; kept for reference only)
    sr = data.get("specificity_results", {}).get("target_specificity", {})
    if "Category" in sr:
        out["sp_category_legacy"] = str(sr["Category"]).strip()
    return out


def process_qc(qc_paths_dict,
               thresh_homotopic=0.1, thresh_delta=0.05):
    """
    Build QC DataFrame from per-session QC JSON files and correlation matrices.

    Homotopic bilateral specificity is computed directly from each session's
    correlation matrix CSV using compute_homotopic_specificity().  This
    replaces the circuit-based (primate-centric) specificity from full_results.json.

    New columns added to df_qc:
      sp_category             : Specific / Unspecific / No / Spurious
      sp_specific_correlation : mean homotopic (L_X ↔ R_X) correlation
      sp_nonspecific_correlation: mean cross-network correlation
      sp_specificity_index    : homotopic_mean − cross_mean
      net_homotopic_mean      : same as sp_specific_correlation (alias)
      net_cross_mean          : same as sp_nonspecific_correlation (alias)
      n_homotopic_pairs       : number of L/R matched pairs used
      sp_category_legacy      : original JSON-based category (for reference)

    Parameters
    ----------
    thresh_homotopic : r threshold below which homotopic_mean → "No" category
    thresh_delta     : minimum specificity_index for "Specific" category
    """
    n_sub = len(qc_paths_dict["subject"])
    records = []
    for i in range(n_sub):
        sub    = qc_paths_dict["subject"][i]
        ses    = qc_paths_dict["session"][i]
        func_p = qc_paths_dict["func_qc_path"][i]
        anat_p = qc_paths_dict["anat_qc_path"][i]
        full_p = qc_paths_dict.get("full_results_path", [None]*n_sub)[i]
        corr_p = qc_paths_dict.get("corr_matrix_path",  [None]*n_sub)[i]

        row = {"subject": sub, "session": ses}

        # Functional QC JSON
        if func_p and os.path.exists(func_p):
            try: row.update({f"func_{k}": v
                             for k, v in _flatten(_load_json_with_nan(func_p)).items()})
            except Exception as e: warnings.warn(f"func QC error {func_p}: {e}")

        # Anatomical QC JSON
        if anat_p and os.path.exists(anat_p):
            try: row.update({f"anat_{k}": v
                             for k, v in _flatten(_load_json_with_nan(anat_p)).items()})
            except Exception as e: warnings.warn(f"anat QC error {anat_p}: {e}")

        # Network metrics from full_results.json (no sp_category here)
        if full_p and os.path.exists(full_p):
            row.update(_parse_full_results(full_p))

        # Homotopic specificity from correlation matrix CSV
        if corr_p and os.path.exists(corr_p):
            try:
                rois, mat = load_corr_matrix(corr_p)
                hspec = compute_homotopic_specificity(
                    mat, rois,
                    thresh_homotopic=thresh_homotopic,
                    thresh_delta=thresh_delta)
                row["sp_category"]              = hspec["sp_category"]
                row["sp_specific_correlation"]  = hspec["homotopic_mean"]
                row["sp_nonspecific_correlation"]= hspec["cross_mean"]
                row["sp_specificity_index"]     = hspec["specificity_index"]
                row["net_homotopic_mean"]        = hspec["homotopic_mean"]
                row["net_cross_mean"]            = hspec["cross_mean"]
                row["n_homotopic_pairs"]         = hspec["n_homotopic_pairs"]
            except Exception as e:
                warnings.warn(f"homotopic specificity error {corr_p}: {e}")
        elif "sp_category_legacy" in row:
            # Fallback: use legacy JSON category if no corr matrix available
            row["sp_category"] = row["sp_category_legacy"]

        records.append(row)

    df      = pd.DataFrame(records)
    id_cols = [c for c in ("subject","session","species","bids_dir") if c in df.columns]
    data_cols = [c for c in df.columns if c not in id_cols]
    if data_cols: df = df.dropna(subset=data_cols, how="all")

    # Summary
    if "sp_category" in df.columns:
        cat_counts = df["sp_category"].value_counts()
        print(f"  [process_qc] sp_category: {dict(cat_counts)}")

    return df


def collect_multi_species(species_config, regions_of_interest=None,
                           extract=("surface","volume","thickness","qc"),
                           atlas_name="EDNIxCSC", atlas_label_paths=None, fit_kind='correlation',
                           atlas_library_root=None, species_atlas_fragments=None):
    atlas_label_paths       = atlas_label_paths       or {}
    species_atlas_fragments = species_atlas_fragments or {}
    accumulators   = {m: [] for m in ("surface","volume","thickness","qc")}
    bold_paths_all = {}

    for species, cfg in species_config.items():
        bids_dirs = cfg.get("bids_dirs", [])
        print(bids_dirs)
        if isinstance(bids_dirs, str): bids_dirs = [bids_dirs]
        lk = cfg.get("list_to_keep",   [])
        lr = cfg.get("list_to_remove", [])
        print(f"\n{'='*64}\n  Species: {species}  ({len(bids_dirs)} BIDS dir(s))\n{'='*64}")

        label_path = atlas_label_paths.get(species)
        if label_path is None and atlas_library_root and species in species_atlas_fragments:
            try:
                label_path = get_atlas_label_path(
                    atlas_library_root, species_atlas_fragments[species], atlas_name)
            except FileNotFoundError as e:
                warnings.warn(str(e))

        bold_paths_all[species] = []

        for bids_dir in bids_dirs:
            bids_dir = _linux_path(bids_dir)
            print(bids_dir)
            bids_lbl = os.path.basename(bids_dir)
            print(bids_lbl)
            print(f"\n  ↳ BIDS: {bids_dir}")

            if "surface" in extract:
                sp = extract_surface_paths(bids_dir, lk, lr)
                df = process_surfaces(sp, regions_of_interest)
                if not df.empty:
                    df.insert(0, "bids_dir", bids_lbl)
                    df.insert(0, "species",  species)
                    accumulators["surface"].append(df)

            if "thickness" in extract:
                tp = extract_thickness_paths(bids_dir, lk, lr)
                df = process_thickness(tp, regions_of_interest)
                if not df.empty:
                    df.insert(0, "bids_dir", bids_lbl)
                    df.insert(0, "species",  species)
                    accumulators["thickness"].append(df)

            if "volume" in extract and label_path:
                df = process_volumes_and_save(
                    bids_dir, atlas_name, label_path, lk, lr,
                    regions_of_interest, overwrite=False)
                if not df.empty:
                    df.insert(0, "bids_dir", bids_lbl)
                    df.insert(0, "species",  species)
                    accumulators["volume"].append(df)
            elif "volume" in extract:
                warnings.warn(f"No label file for {species} → volume skipped.")

            if "qc" in extract:
                qp = extract_qc_paths(bids_dir, lk, lr, fit_kind)
                df = process_qc(qp)
                if not df.empty:
                    df.insert(0, "bids_dir", bids_lbl)
                    print("insert " + str(bids_lbl))
                    df.insert(0, "species",  species)
                    accumulators["qc"].append(df)

    result = {
        k: pd.concat(v, ignore_index=True) if v else pd.DataFrame()
        for k, v in accumulators.items()
    }
    result["bold_paths"] = bold_paths_all
    return result


def collect_corr_matrices(species_config, atlas_name="EDNIxCSC", atlas_level=3, fit_kind='correlation',
                           use_lr=False, save_csv_dir=None):
    from scipy import stats as _stats
    species_stacks, species_rois = {}, {}
    for species, cfg in species_config.items():
        bids_dirs = cfg.get("bids_dirs", [])
        if isinstance(bids_dirs, str): bids_dirs = [bids_dirs]
        lk, lr = cfg.get("list_to_keep", []), cfg.get("list_to_remove", [])
        sub_ses_runs = {}; rois_ref = None; n_loaded = n_skipped = 0
        for bids_dir in bids_dirs:
            for rec in extract_corr_matrix_paths(bids_dir, atlas_name, atlas_level, fit_kind, use_lr, lk, lr):
                try:
                    rois, mat = load_corr_matrix(rec["path"])
                    valid = mat[~np.isnan(mat)]
                    if valid.size == 0 or np.all(valid == 0): n_skipped += 1; continue
                    if rois_ref is None: rois_ref = rois
                    sub_ses_runs.setdefault((rec["subject"], rec["session"]), []).append((rois, mat))
                    n_loaded += 1
                except Exception as e:
                    print(f"    [ERROR] {rec['path']}: {e}"); n_skipped += 1
        print(f"  [{species}] loaded={n_loaded}  skipped={n_skipped}")
        if not sub_ses_runs: continue
        subject_means = []
        for (sub, ses), run_list in sorted(sub_ses_runs.items()):
            common = set(run_list[0][0]).intersection(*[set(r) for r, _ in run_list[1:]])
            ref    = [r for r in run_list[0][0] if r in common]
            if not ref: continue
            mats = []
            for rr, mm in run_list:
                try:
                    idx = [list(rr).index(r) for r in ref]
                    mats.append(mm[np.ix_(idx, idx)])
                except ValueError: pass
            if mats: subject_means.append((ref, np.nanmean(np.stack(mats, 0), 0)))
        if not subject_means: continue
        all_roi_sets = [set(r) for r, _ in subject_means]
        common_sp    = all_roi_sets[0].intersection(*all_roi_sets[1:])
        rois_ref     = [r for r in subject_means[0][0] if r in common_sp]
        aligned      = []
        for rr, mm in subject_means:
            idx = [list(rr).index(r) for r in rois_ref if r in rr]
            if len(idx) == len(rois_ref): aligned.append(mm[np.ix_(idx, idx)])
        if not aligned: continue
        species_stacks[species] = aligned; species_rois[species] = rois_ref
    if not species_stacks: return {}
    # Each species keeps its OWN ROI set — no cross-species intersection.
    # Different species use different atlases with different ROI names;
    # intersecting across species would give an empty set (e.g. Rat ∩ Human = {}).
    result = {}
    for species, mats in species_stacks.items():
        rois  = species_rois[species]
        stack = np.stack(mats, 0)
        n     = stack.shape[0]
        mean  = np.nanmean(stack, 0); var = np.nanvar(stack, 0)
        _, p  = (_stats.ttest_1samp(stack, 0, axis=0, nan_policy="omit")
                 if n >= 2 else (None, np.full(mean.shape, np.nan)))
        result[species] = dict(mean=mean, var=var, pval=p, rois=rois, n=n)
        if save_csv_dir:
            os.makedirs(save_csv_dir, exist_ok=True)
            for mat, nm in [(mean,"mean"),(var,"var"),(p,"pval")]:
                pd.DataFrame(mat, index=rois, columns=rois).to_csv(
                    opj(save_csv_dir, f"{species}_corr_{nm}.csv"))
    return result


def collect_corr_matrices_per_bids(species_config, atlas_name="EDNIxCSC", fit_kind='correlation',
                                    atlas_level=3, use_lr=False, save_csv_dir=None):
    from scipy import stats as _stats
    bids_stacks = {}
    for species, cfg in species_config.items():
        bids_dirs = cfg.get("bids_dirs", [])
        if isinstance(bids_dirs, str): bids_dirs = [bids_dirs]
        lk, lr = cfg.get("list_to_keep", []), cfg.get("list_to_remove", [])
        for bids_dir in bids_dirs:
            bids_lbl = _bids_label(bids_dir)
            recs = extract_corr_matrix_paths(bids_dir, atlas_name, atlas_level, fit_kind, use_lr, lk, lr)
            sub_ses_runs = {}; rois_ref = None; n_loaded = n_skipped = 0
            for rec in recs:
                try:
                    rois, mat = load_corr_matrix(rec["path"])
                    valid = mat[~np.isnan(mat)]
                    if valid.size == 0 or np.all(valid == 0): n_skipped += 1; continue
                    if rois_ref is None: rois_ref = rois
                    sub_ses_runs.setdefault((rec["subject"], rec["session"]), []).append((rois, mat))
                    n_loaded += 1
                except Exception as e:
                    print(f"    [ERROR] {rec['path']}: {e}"); n_skipped += 1
            print(f"  [{species}/{bids_lbl}] loaded={n_loaded} skipped={n_skipped}")
            if not sub_ses_runs: continue
            subject_means = []
            for (sub, ses), run_list in sorted(sub_ses_runs.items()):
                common = set(run_list[0][0]).intersection(*[set(r) for r, _ in run_list[1:]])
                ref    = [r for r in run_list[0][0] if r in common]
                if not ref: continue
                mats = []
                for rr, mm in run_list:
                    try:
                        idx = [list(rr).index(r) for r in ref]
                        mats.append(mm[np.ix_(idx, idx)])
                    except ValueError: pass
                if mats: subject_means.append((ref, np.nanmean(np.stack(mats, 0), 0)))
            if not subject_means: continue
            all_roi_sets = [set(r) for r, _ in subject_means]
            common_sp    = all_roi_sets[0].intersection(*all_roi_sets[1:])
            rois_ref     = [r for r in subject_means[0][0] if r in common_sp]
            aligned      = []
            for rr, mm in subject_means:
                idx = [list(rr).index(r) for r in rois_ref if r in rr]
                if len(idx) == len(rois_ref): aligned.append(mm[np.ix_(idx, idx)])
            if not aligned: continue
            bids_stacks[bids_lbl] = {"species": species, "mats": aligned, "rois": rois_ref}
    result = {}
    for bids_lbl, d in bids_stacks.items():
        stack = np.stack(d["mats"], 0); n = stack.shape[0]
        mean  = np.nanmean(stack, 0); var = np.nanvar(stack, 0)
        _, p  = (_stats.ttest_1samp(stack, 0, axis=0, nan_policy="omit")
                 if n >= 2 else (None, np.full(mean.shape, np.nan)))
        result[bids_lbl] = dict(mean=mean, var=var, pval=p,
                                rois=d["rois"], n=n, species=d["species"])
        if save_csv_dir:
            os.makedirs(save_csv_dir, exist_ok=True)
            for mat, nm in [(mean,"mean"),(var,"var"),(p,"pval")]:
                pd.DataFrame(mat, index=d["rois"], columns=d["rois"]).to_csv(
                    opj(save_csv_dir, f"{bids_lbl}_corr_{nm}.csv"))
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# §7  STAGE 2 — EXPORT DATAFRAMES → EXCEL / CSV
# ═══════════════════════════════════════════════════════════════════════════════

def export_morphometry_xlsx(data_dict, output_path):
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for mod in ("surface","thickness","volume","qc"):
            df = data_dict.get(mod)
            if df is None or (isinstance(df, pd.DataFrame) and df.empty): continue
            df.to_excel(writer, sheet_name=mod[:31], index=False)
            print(f"  [export] sheet '{mod}' → {len(df):,} rows")
    print(f"\n✓ Morphometry data → {output_path}")
    return output_path


def export_to_excel(data_dict, output_path, modalities=None):
    return export_morphometry_xlsx(
        {k: v for k, v in data_dict.items()
         if isinstance(v, pd.DataFrame) and (not modalities or k in modalities)},
        output_path)


def export_to_csv(data_dict, output_dir, modalities=None):
    modalities = modalities or [k for k, v in data_dict.items()
                                if isinstance(v, pd.DataFrame) and not v.empty]
    os.makedirs(output_dir, exist_ok=True)
    out = {}
    for mod in modalities:
        df = data_dict.get(mod)
        if df is None or (isinstance(df, pd.DataFrame) and df.empty): continue
        p = opj(output_dir, f"{mod}.csv")
        df.to_csv(p, index=False)
        out[mod] = p
        print(f"  [export_csv] {p}")
    return out


def export_summary_stats(data_dict, output_path, groupby=("species",), modalities=None):
    metric_cols = {"surface":"surface_area_mm2","volume":"volume_mm3","thickness":"thickness_mm"}
    modalities  = modalities or list(metric_cols.keys())
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for mod in modalities:
            df  = data_dict.get(mod)
            if df is None or (isinstance(df, pd.DataFrame) and df.empty): continue
            col = metric_cols.get(mod)
            if col not in df.columns: continue
            grp_cols = [g for g in groupby if g in df.columns]
            for extra in ("atlas_level","hemisphere"):
                if extra in df.columns: grp_cols.append(extra)
            grp_cols.append("region")
            summary = (df.groupby(grp_cols)[col]
                       .agg(n="count",mean="mean",std="std",
                            median="median",min="min",max="max")
                       .reset_index())
            summary.to_excel(writer, sheet_name=f"{mod}_summary"[:31], index=False)
    print(f"✓ Summary stats → {output_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# §8  DIAGNOSTICS
# ═══════════════════════════════════════════════════════════════════════════════

def diagnose_surface_xlsx(bids_dir, n_subjects=3):
    print(f"\n{'='*70}\n  DIAGNOSE surface.xlsx  —  {bids_dir}\n{'='*70}")
    for count, (sub, ses, _, ses_dir) in enumerate(_scan_bids(bids_dir)):
        if count >= n_subjects: break
        native_dir = opj(ses_dir, "anat", "native", "surfaces", "Native_resol")
        xlsx       = opj(native_dir, "surface.xlsx")
        print(f"\n  sub-{sub}  ses-{ses or '1'}")
        print(f"    native_dir : {os.path.exists(native_dir)}")
        print(f"    xlsx exists: {os.path.exists(xlsx)}")
        if os.path.exists(xlsx):
            try:
                df = pd.read_excel(xlsx, index_col=0)
                print(f"    shape      : {df.shape}")
                for col in df.columns:
                    print(f"      {col}: {df[col].notna().sum()}/{len(df)} non-NaN")
            except Exception as e:
                print(f"    ERROR: {e}")
        seg = _get_seg_file(native_dir, sub)
        print(f"    seg file   : {seg}")
        for h in ("l","r"):
            print(f"    {h}.pial : {_find_surf_gii(native_dir, sub, h, 'pial')}")
    print(f"\n{'='*70}\n")


def diagnose_thickness_xlsx(bids_dir, n_subjects=3):
    print(f"\n{'='*70}\n  DIAGNOSE thickness.xlsx  —  {bids_dir}\n{'='*70}")
    for count, (sub, ses, _, ses_dir) in enumerate(_scan_bids(bids_dir)):
        if count >= n_subjects: break
        native_dir = opj(ses_dir, "anat", "native", "surfaces", "Native_resol")
        xlsx       = opj(native_dir, "thickness.xlsx")
        print(f"\n  sub-{sub}  ses-{ses or '1'}")
        print(f"    xlsx exists: {os.path.exists(xlsx)}")
        if os.path.exists(xlsx):
            try:
                df = pd.read_excel(xlsx, index_col=0)
                print(f"    shape      : {df.shape}")
                for col in df.columns:
                    print(f"      {col}: {df[col].notna().sum()}/{len(df)} non-NaN")
            except Exception as e:
                print(f"    ERROR: {e}")
    print(f"\n{'='*70}\n")


def diagnose_atlas_config(label_path, legend_path=None, n_show=15):
    print(f"\n{'='*70}\n  DIAGNOSE atlas_config\n{'='*70}")
    cfg = build_atlas_config(label_path, legend_path)
    df  = cfg["label_df"]
    print(f"  Label file rows : {len(df)}")
    print(f"  Unique base names: {len(cfg['name_to_ids'])}")
    print(f"  Hemisphere distribution:")
    print(df["hemisphere"].value_counts().to_string(header=False))
    if legend_path:
        rl = cfg["region_to_level"]
        from collections import Counter
        print(f"  Region-to-level distribution: {dict(Counter(rl.values()))}")
    print(f"\n  First {n_show} entries:")
    for base, ids in list(cfg["name_to_ids"].items())[:n_show]:
        lvl = cfg["region_to_level"].get(base, "?")
        print(f"    {base:50s}  L={ids['L']}  R={ids['R']}  "
              f"bil={ids['bilateral']}  level={lvl}")
    print(f"\n{'='*70}\n")


# ═══════════════════════════════════════════════════════════════════════════════
# §9  PLOTTING HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

def _fmt_y1(ax):
    import matplotlib.ticker as ticker
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))


def _violin_strip_quartiles(ax, data, x_pos, color, width=0.4):
    data = np.asarray(data, dtype=float)
    data = data[~np.isnan(data)]
    if len(data) == 0: return
    if np.std(data) == 0 or len(data) < 2:
        ax.hlines(data[0], x_pos-width/2, x_pos+width/2, colors=color, linewidths=2, alpha=0.6)
        return
    from scipy.stats import gaussian_kde
    try:
        kde = gaussian_kde(data, bw_method="scott")
    except Exception: return
    y_g  = np.linspace(data.min(), data.max(), 200)
    dens = kde(y_g); dens = dens / dens.max() * (width / 2)
    ax.fill_betweenx(y_g, x_pos-dens, x_pos+dens, color=color, alpha=0.25, linewidth=0)
    ax.plot(x_pos-dens, y_g, color=color, lw=0.8, alpha=0.5)
    ax.plot(x_pos+dens, y_g, color=color, lw=0.8, alpha=0.5)
    for q, lw in zip(np.percentile(data, [25,50,75]), [1,2,1]):
        dq = float(kde(q)) / kde(y_g).max() * (width/2)
        ax.hlines(q, x_pos-dq, x_pos+dq, colors=color, linewidths=lw, zorder=3)


def _scatter_with_anesth(ax, x_vals, y_vals, subjects, sessions, anesth_map, color, rng):
    from collections import defaultdict
    if anesth_map and subjects is not None:
        grp = defaultdict(lambda: ([], []))
        for xv, yv, sub, ses in zip(x_vals, y_vals, subjects, sessions):
            m = _anesthesia_marker(
                anesth_map.get((str(sub), str(ses)),
                               anesth_map.get((str(sub), "1"), "")))
            grp[m][0].append(xv); grp[m][1].append(yv)
        for m, (xs, ys) in grp.items():
            ax.scatter(xs, ys, color=color, marker=m, s=22, zorder=5,
                       alpha=0.82, linewidths=0.5, edgecolors="k")
    else:
        ax.scatter(x_vals, y_vals, color=color, marker=AWAKE_MARKER, s=22,
                   zorder=5, alpha=0.82, linewidths=0.5, edgecolors="k")


def _filter_atlas_level(df, atlas_level):
    if df is None or df.empty: return df
    if "atlas_level" in df.columns: return df[df["atlas_level"] == atlas_level]
    return df


def _hemis_for(df):
    if df is None or df.empty or "hemisphere" not in df.columns: return ["bilateral"]
    return sorted(df["hemisphere"].unique())


def _anesthesia_legend_handles():
    import matplotlib.lines as mlines
    return [
        mlines.Line2D([], [], marker=AWAKE_MARKER,  color="grey",
                      linestyle="None", markersize=7, label="Awake"),
        mlines.Line2D([], [], marker=ANESTH_MARKER, color="grey",
                      linestyle="None", markersize=7, label="Anaesthetised"),
    ]


def plot_morphometry_intra_bids(df, metric_col, metric_label, regions, output_path,
                                 species=None, hemisphere=None, bids_label=None,
                                 atlas_level=1, bids_col=None, figsize=None,
                                 anesth_map=None, rng=None):
    _rng = rng or np.random.default_rng(42)
    with plt.rc_context(PAPER_RC):
        plot_df = df.copy()
        if hemisphere:
            filt = plot_df[plot_df["hemisphere"] == hemisphere] if "hemisphere" in plot_df.columns else plot_df
            plot_df = filt if not filt.empty else plot_df
        if "atlas_level" in plot_df.columns:
            plot_df = plot_df[plot_df["atlas_level"] == atlas_level]
        plot_df = plot_df[_df_region_mask(plot_df["region"], regions)]
        if plot_df.empty:
            warnings.warn("plot_morphometry_intra_bids: no data"); return None
        n    = len(regions)
        w, h = figsize or (max(5, n * 2.4), 5)
        fig, ax = plt.subplots(figsize=(w, h))
        for i, region in enumerate(regions):
            rdf  = plot_df[_df_region_mask(plot_df["region"], [region])]
            vals = rdf[metric_col].dropna().values
            _violin_strip_quartiles(ax, vals, i, PALETTE[0])
            subs = rdf.loc[rdf[metric_col].notna(), "subject"].values if "subject" in rdf.columns else None
            sess = rdf.loc[rdf[metric_col].notna(), "session"].values if "session" in rdf.columns else None
            jitter = _rng.uniform(-0.12, 0.12, len(vals))
            _scatter_with_anesth(ax, i+jitter, vals, subs, sess, anesth_map, PALETTE[0], _rng)
        ax.set_xticks(range(n))
        ax.set_xticklabels(regions, rotation=35, ha="right")
        ax.set_ylabel(metric_label); _fmt_y1(ax)
        title = metric_label
        if species:    title += f" — {species}"
        if bids_label: title += f" ({bids_label})"
        ax.set_title(title, fontweight="bold", pad=8)
        ax.set_xlim(-0.6, n - 0.4)
        if anesth_map:
            ax.legend(handles=_anesthesia_legend_handles(), loc="upper right",
                      fontsize=8, frameon=False)
        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight"); plt.close(fig)
        print(f"  [plot] {output_path}")
    return output_path


def plot_qc_dashboard(qc_df, output_path, func_metrics=None, anat_metrics=None,
                      figsize=(14, 8), anesth_map=None):
    with plt.rc_context(PAPER_RC):
        auto_func = ["func_avg_snr_gray","func_TSNR_0","func_mean_fd","func_gcor"]
        auto_anat = ["anat_template_correlation","anat_cortical_contrast"]
        func_metrics = func_metrics or [m for m in auto_func if m in qc_df.columns]
        anat_metrics = anat_metrics or [m for m in auto_anat if m in qc_df.columns]
        all_m = [m for m in func_metrics + anat_metrics if m in qc_df.columns]
        if not all_m: warnings.warn("No QC metrics found"); return None
        has_sp = "species" in qc_df.columns
        colors = ({s: PALETTE[i % len(PALETTE)] for i, s in
                   enumerate(sorted(qc_df["species"].unique()))} if has_sp else None)
        ncols = min(4, len(all_m)); nrows = math.ceil(len(all_m) / ncols)
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        axes = np.array(axes).flatten()
        for i, metric in enumerate(all_m):
            ax   = axes[i]
            data = qc_df[["species", metric]].dropna() if has_sp else qc_df[[metric]].dropna()
            if has_sp:
                sns.boxplot(data=data, x="species", y=metric, palette=colors,
                            width=0.45, linewidth=1.1, ax=ax,
                            flierprops=dict(marker="o", markersize=3, alpha=0.5))
                sns.stripplot(data=data, x="species", y=metric, palette=colors,
                              size=4, alpha=0.65, jitter=True, ax=ax)
                ax.tick_params(axis="x", rotation=30)
            else:
                ax.boxplot(data[metric].values, widths=0.4)
            ax.set_title(metric.replace("func_","").replace("anat_",""), fontweight="bold")
            ax.set_xlabel(""); ax.set_ylabel(""); _fmt_y1(ax)
        for j in range(i+1, len(axes)): axes[j].set_visible(False)
        fig.suptitle("QC Dashboard", fontsize=14, fontweight="bold")
        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight"); plt.close(fig)
        print(f"  [plot] {output_path}")
    return output_path