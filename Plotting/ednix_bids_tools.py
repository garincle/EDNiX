"""
EDNiX BIDS Analysis Tools  ?  v3  (patched)
============================================
Patch summary (atlas-level management fixes)
--------------------------------------------

FIX 1 ? _read_one_surface / _read_one_thickness
    All Surface_Area_lvlN / Thickness_mm_lvlN columns are extracted.
    When the SAME base-region name appears at multiple atlas levels,
    the LOWEST level wins (lvl1 > lvl2 > lvl3 > lvl4).
    The resulting DataFrame always has exactly one row per
    (subject, session, region_full, hemisphere) tuple.

FIX 2 ? _save_surface_xlsx / _save_thickness_xlsx
    The old .xlsx is DELETED before writing a fresh one.
    No more column-merge logic; each call produces a clean file.
    The file now contains one column per atlas level that was
    successfully extracted (Surface_Area_lvl1 ? lvl4 / Thickness_mm_lvl1 ? lvl4).

FIX 3 ? process_volumes / export_summary_stats
    Volumes have no atlas_level column.
    export_summary_stats no longer tries to group by atlas_level for volumes.
    A safe helper _group_cols() picks grouping columns that actually exist
    in the DataFrame being summarised.

No other logic was changed.
"""

import os
import re
import glob
import json
import math
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
    warnings.warn("joblib not found ? surface/thickness reading will be sequential. "
                  "Install with: pip install joblib")

opj = os.path.join
opn = os.path.normpath

# Number of atlas hierarchy levels stored in the dlabel cifti
N_ATLAS_LEVELS = 4


# ??????????????????????????????????????????????????????????????????????????????
# PART 0  ?  HELPERS & LABEL PARSING
# ??????????????????????????????????????????????????????????????????????????????

def _linux_path(p: str) -> str:
    return str(p).replace("\\", "/")


def _parse_bids_entity(path: str, entity: str):
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


def _spco(cmd: str):
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            warnings.warn(f"Command failed: {cmd}\n{result.stderr.strip()}")
            return None
        return result.stdout
    except subprocess.TimeoutExpired:
        warnings.warn("Command timed out...")
        return None


# ?????????????????????????????????????????????????????????????????????????????
# Label-file parsing
# ?????????????????????????????????????????????????????????????????????????????

def parse_label_file(label_path: str) -> pd.DataFrame:
    with open(label_path, encoding='utf-8-sig') as f:
        lines = [l.strip() for l in f if l.strip()]

    def _parse_hemi(name):
        if name.startswith("L_"):
            return "left",  name[2:]
        elif name.startswith("R_"):
            return "right", name[2:]
        else:
            return "bilateral", name

    first_parts = lines[0].split()
    is_statslut = len(first_parts) >= 6 and first_parts[0].isdigit()
    rows = []

    if is_statslut:
        for line in lines:
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                label_id = int(parts[0])
            except ValueError:
                continue
            name = parts[1]
            try:
                r, g, b, a = int(parts[2]), int(parts[3]), int(parts[4]), int(parts[5])
            except (ValueError, IndexError):
                r = g = b = a = 0
            hemi, base = _parse_hemi(name)
            rows.append(dict(region_name=name, base_region=base,
                             label_id=label_id, hemisphere=hemi,
                             R=r, G=g, B=b, A=a))
    else:
        i = 0
        while i < len(lines) - 1:
            name_line  = lines[i]
            color_line = lines[i + 1]
            parts = color_line.split()
            if len(parts) == 5 and parts[0].isdigit():
                label_id = int(parts[0])
                r, g, b, a = int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])
                hemi, base = _parse_hemi(name_line)
                rows.append(dict(region_name=name_line, base_region=base,
                                 label_id=label_id, hemisphere=hemi,
                                 R=r, G=g, B=b, A=a))
                i += 2
            else:
                i += 1

    return pd.DataFrame(rows)


def get_atlas_label_path(
    atlas_library_root: str,
    species_path_fragment: str,
    atlas_name: str,
    prefer_statslut: bool = False,
) -> str:
    label_dir = opj(atlas_library_root, species_path_fragment, "label_code")

    if prefer_statslut:
        priority = [f"{atlas_name}_StatsLUT.txt", f"{atlas_name}_label.txt"]
    else:
        priority = [f"{atlas_name}_label.txt", f"{atlas_name}_StatsLUT.txt"]

    for filename in priority:
        candidate = opj(label_dir, filename)
        if os.path.exists(candidate):
            return _linux_path(candidate)

    candidates = glob.glob(opj(label_dir, f"{atlas_name}_label.txt"))
    if candidates:
        return _linux_path(sorted(candidates)[0])

    raise FileNotFoundError(
        f"No label file for atlas '{atlas_name}' in {label_dir}\n"
        f"Looked for: {priority}"
    )


def find_species_path(atlas_library_root: str, species_name: str) -> str:
    atlas_root = opj(atlas_library_root, "atlas")
    matches = []

    for dirpath, dirnames, filenames in os.walk(atlas_root):
        folder_name = os.path.basename(dirpath)
        if folder_name.lower() == species_name.lower():
            rel = os.path.relpath(dirpath, atlas_library_root)
            matches.append(_linux_path(rel))
            dirnames.clear()

    if not matches:
        raise FileNotFoundError(
            f"No species folder named '{species_name}' found under {atlas_root}"
        )
    if len(matches) > 1:
        raise ValueError(
            f"Ambiguous species name '{species_name}' ? found multiple matches:\n"
            + "\n".join(f"  {m}" for m in matches)
        )
    return matches[0]


# ??????????????????????????????????????????????????????????????????????????????
# PART 1  ?  PATH EXTRACTORS
# ??????????????????????????????????????????????????????????????????????????????

def _scan_bids(bids_dir, sub_glob, ses_glob="ses-*"):
    for sub_dir in sorted(glob.glob(opj(bids_dir, sub_glob))):
        sub = _parse_bids_entity(sub_dir, "sub")
        if sub is None:
            continue
        ses_dirs = sorted(glob.glob(opj(sub_dir, ses_glob)))
        if not ses_dirs:
            yield sub, None, sub_dir, sub_dir
        else:
            for ses_dir in ses_dirs:
                ses = _parse_bids_entity(ses_dir, "ses")
                yield sub, ses, sub_dir, ses_dir


def _parse_run(filename: str):
    m = re.search(r"run-(\w+)", os.path.basename(filename))
    return m.group(1) if m else None


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


def build_func_index(bold_paths_dict, sba_seed_name=None, mean_img_pattern="*mean*.nii.gz"):
    records = []
    for sub, ses, run, bold_p in zip(
        bold_paths_dict["subject"], bold_paths_dict["session"],
        bold_paths_dict["run"],    bold_paths_dict["bold_path"],
    ):
        bn   = os.path.basename(bold_p)
        root = bn
        for suffix in (
            "_space-template_desc-fMRI_residual.nii.gz",
            "_space-template_desc-fMRI_residual.nii",
            ".nii.gz", ".nii",
        ):
            if root.endswith(suffix):
                root = root[: -len(suffix)]; break

        bold_dir    = os.path.dirname(bold_p)
        func_dir    = _func_dir_from_bold(bold_p)
        postpro_dir = bold_dir

        mean_img = None
        if mean_img_pattern:
            candidates = glob.glob(opj(postpro_dir, mean_img_pattern))
            shared = [c for c in candidates if sub in c]
            mean_img = _linux_path(shared[0]) if shared else (
                _linux_path(candidates[0]) if candidates else None)

        sba_fish = None
        if sba_seed_name:
            sba_dir       = opj(postpro_dir, "10_Results", "SBA", sba_seed_name)
            sba_fish_cand = opj(sba_dir, f"{root}_correlations_fish.nii.gz")
            sba_fish      = _linux_path(sba_fish_cand) if os.path.exists(sba_fish_cand) else None

        run_idx       = run if run is not None else "01"
        records.append(dict(
            subject=sub, session=ses, run=run, bold_path=bold_p,
            bold_root=root, func_dir=func_dir, postpro_dir=postpro_dir,
            mean_img_path=mean_img, sba_fish_path=sba_fish,
            sba_fish_exists=os.path.exists(sba_fish) if sba_fish else False,
            run_label=f"run_{run_idx}", session_label=f"Sess_{ses}",
        ))
    return pd.DataFrame(records)


def _func_dir_from_bold(bold_path: str) -> str:
    p = Path(bold_path)
    for parent in p.parents:
        if parent.name == "func":
            return _linux_path(str(parent))
    return _linux_path(str(p.parents[3]))


def extract_surface_paths(bids_dir, list_to_keep=None, list_to_remove=None,
                           surface_pattern="surface.xlsx"):
    list_to_keep   = list_to_keep   or []
    list_to_remove = list_to_remove or []
    pattern = opj(bids_dir, "sub-*", "ses-*", "anat", "native", "surfaces",
                  "Native_resol", surface_pattern)
    subjects, sessions, paths = [], [], []
    for f in sorted(glob.glob(pattern)):
        f   = _linux_path(f)
        sub = _parse_bids_entity(f, "sub")
        ses = _parse_bids_entity(f, "ses") or "1"
        subjects.append(sub); sessions.append(ses); paths.append(f)

    keep = _subject_session_filter(subjects, sessions, list_to_keep, list_to_remove)
    out  = {"subject": [], "session": [], "surface_path": []}
    for sub, ses, p in zip(subjects, sessions, paths):
        if (sub, ses) in keep:
            out["subject"].append(sub); out["session"].append(ses)
            out["surface_path"].append(p)

    print(f"  [surface]   {len(out['surface_path'])} files  ? {bids_dir}")
    return out


def extract_thickness_paths(bids_dir, list_to_keep=None, list_to_remove=None,
                             thickness_pattern="thickness.xlsx"):
    list_to_keep   = list_to_keep   or []
    list_to_remove = list_to_remove or []
    pattern = opj(bids_dir, "sub-*", "ses-*", "anat", "native", "surfaces",
                  "**", thickness_pattern)
    subjects, sessions, paths = [], [], []
    for f in sorted(glob.glob(pattern, recursive=True)):
        f   = _linux_path(f)
        sub = _parse_bids_entity(f, "sub")
        ses = _parse_bids_entity(f, "ses") or "1"
        subjects.append(sub); sessions.append(ses); paths.append(f)

    keep = _subject_session_filter(subjects, sessions, list_to_keep, list_to_remove)
    out  = {"subject": [], "session": [], "thickness_path": []}
    for sub, ses, p in zip(subjects, sessions, paths):
        if (sub, ses) in keep:
            out["subject"].append(sub); out["session"].append(ses)
            out["thickness_path"].append(p)

    print(f"  [thickness] {len(out['thickness_path'])} files  ? {bids_dir}")
    return out


def extract_volume_paths(bids_dir, atlas_name, list_to_keep=None, list_to_remove=None):
    list_to_keep   = list_to_keep   or []
    list_to_remove = list_to_remove or []
    pattern = opj(bids_dir, "sub-*", "ses-*", "anat", "native", "volumes",
                  "labels", f"*_seg-{atlas_name}_dseg.nii.gz")
    subjects, sessions, paths = [], [], []
    for f in sorted(glob.glob(pattern)):
        f   = _linux_path(f)
        sub = _parse_bids_entity(f, "sub")
        ses = _parse_bids_entity(f, "ses") or "1"
        subjects.append(sub); sessions.append(ses); paths.append(f)

    keep = _subject_session_filter(subjects, sessions, list_to_keep, list_to_remove)
    out  = {"subject": [], "session": [], "volume_seg_path": []}
    for sub, ses, p in zip(subjects, sessions, paths):
        if (sub, ses) in keep:
            out["subject"].append(sub); out["session"].append(ses)
            out["volume_seg_path"].append(p)

    print(f"  [volume]    {len(out['volume_seg_path'])} seg files  ? {bids_dir}  (atlas={atlas_name})")
    return out


def extract_qc_paths(bids_dir, list_to_keep=None, list_to_remove=None,
                     func_qc_suffix="*_QC_values.json",
                     anat_qc_suffix="*_QC_values.json",
                     full_results_suffix="*_full_results.json"):
    list_to_keep   = list_to_keep   or []
    list_to_remove = list_to_remove or []

    sub_ses = set()
    for sub_dir in glob.glob(opj(bids_dir, "sub-*")):
        sub      = _parse_bids_entity(sub_dir, "sub")
        ses_dirs = glob.glob(opj(sub_dir, "ses-*"))
        if ses_dirs:
            for sd in ses_dirs:
                ses = _parse_bids_entity(sd, "ses")
                if sub and ses:
                    sub_ses.add((sub, ses))
        elif sub:
            sub_ses.add((sub, "1"))

    all_subs = [p[0] for p in sub_ses]
    all_ses  = [p[1] for p in sub_ses]
    keep     = _subject_session_filter(all_subs, all_ses, list_to_keep, list_to_remove)

    out = {"subject": [], "session": [],
           "func_qc_path": [], "anat_qc_path": [], "full_results_path": []}
    for sub, ses in sorted(keep):
        ses_dir_candidate = opj(bids_dir, f"sub-{sub}", f"ses-{ses}")
        ses_dir = ses_dir_candidate if os.path.exists(ses_dir_candidate) \
                  else opj(bids_dir, f"sub-{sub}")
        qc_dir  = opj(ses_dir, "func", "QC")

        func_m = glob.glob(opj(qc_dir, func_qc_suffix))
        func_p = _linux_path(func_m[0]) if func_m else None

        full_m = glob.glob(opj(qc_dir, full_results_suffix))
        full_p = _linux_path(full_m[0]) if full_m else None

        anat_p_cand = glob.glob(opj(ses_dir, "anat", "native", "volumes",
                                    "QC_anat", anat_qc_suffix))
        anat_p = _linux_path(anat_p_cand[0]) if anat_p_cand else None

        out["subject"].append(sub); out["session"].append(ses)
        out["func_qc_path"].append(func_p)
        out["anat_qc_path"].append(anat_p)
        out["full_results_path"].append(full_p)

    n_f = sum(p is not None for p in out["func_qc_path"])
    n_a = sum(p is not None for p in out["anat_qc_path"])
    n_r = sum(p is not None for p in out["full_results_path"])
    print(f"  [QC]        {n_f} func / {n_a} anat / {n_r} full_results  ? {bids_dir}")
    return out


# ??????????????????????????????????????????????????????????????????????????????
# PART 2  ?  DATA PROCESSORS
# ??????????????????????????????????????????????????????????????????????????????

# ?????????????????????????????????????????????????????????????????????????????
# Surface area
# ?????????????????????????????????????????????????????????????????????????????

def _read_one_surface(sub, ses, path, regions_of_interest):
    """
    Read one surface.xlsx and return tidy record dicts.

    All Surface_Area_lvlN columns are read.  When the same region appears
    at more than one atlas level, the LOWEST level wins (lvl1 > lvl2 > ?).
    One output record per unique (region_full, hemisphere) pair.
    """
    records = []
    if not path or not os.path.exists(path):
        return records
    try:
        df = pd.read_excel(path, index_col=0)

        # Collect all level columns, sorted ascending so lvl1 is processed first
        surf_cols = sorted(
            [c for c in df.columns if re.search(r'lvl\d+', c, re.IGNORECASE)
             and ("surface" in c.lower() or "area" in c.lower())],
            key=lambda c: int(re.search(r'lvl(\d+)', c, re.IGNORECASE).group(1))
        )
        # Fallback: any column if none match the pattern
        if not surf_cols:
            surf_cols = list(df.columns)

        # Build one record per region; lowest level wins
        seen: dict[str, dict] = {}   # region_full -> record

        for col in surf_cols:
            lvl_match = re.search(r'lvl(\d+)', col, re.IGNORECASE)
            lvl = int(lvl_match.group(1)) if lvl_match else 1

            for region in df.index:
                rs = str(region)
                if regions_of_interest and not any(r in rs for r in regions_of_interest):
                    continue
                val = df.loc[region, col]
                if pd.isna(val):
                    continue

                # Only store if this region has not been seen yet (lvl1 wins)
                if rs in seen:
                    continue

                if rs.lower().startswith("l_"):
                    hemi, base = "left",  rs[2:]
                elif rs.lower().startswith("r_"):
                    hemi, base = "right", rs[2:]
                else:
                    hemi, base = "bilateral", rs

                seen[rs] = dict(
                    subject=sub, session=ses,
                    region_full=rs, region=base,
                    hemisphere=hemi, atlas_level=lvl,
                    surface_area_mm2=float(val),
                )

        records = list(seen.values())
    except Exception as e:
        warnings.warn(f"Error reading surface {path}: {e}")
    return records


def process_surfaces(surface_paths_dict, regions_of_interest=None, atlas_name="EDNIxCSC"):
    """
    Build tidy DataFrame from pre-computed surface.xlsx files.
    One row per (subject, session, region_full): lowest atlas level wins.
    Uses parallel reading when joblib is available.
    """
    triples = list(zip(
        surface_paths_dict["subject"],
        surface_paths_dict["session"],
        surface_paths_dict["surface_path"],
    ))

    if _JOBLIB:
        results = Parallel(n_jobs=-1)(
            delayed(_read_one_surface)(sub, ses, path, regions_of_interest)
            for sub, ses, path in triples
        )
    else:
        results = [_read_one_surface(sub, ses, path, regions_of_interest)
                   for sub, ses, path in triples]

    records = [r for sublist in results for r in sublist]
    return pd.DataFrame(records)


def _read_one_thickness(sub, ses, path, regions_of_interest):
    """
    Read one thickness.xlsx and return tidy record dicts.

    All Thickness_mm_lvlN columns are read.  When the same region appears
    at more than one atlas level, the LOWEST level wins (lvl1 > lvl2 > ?).
    One output record per unique (region_full, hemisphere) pair.
    """
    records = []
    if not path or not os.path.exists(path):
        return records
    try:
        df = pd.read_excel(path, index_col=0)

        # Collect all level columns, sorted ascending so lvl1 is processed first
        thick_cols = sorted(
            [c for c in df.columns if re.search(r'lvl\d+', c, re.IGNORECASE)
             and ("thick" in c.lower() or "mean" in c.lower())],
            key=lambda c: int(re.search(r'lvl(\d+)', c, re.IGNORECASE).group(1))
        )
        if not thick_cols:
            thick_cols = list(df.columns)

        seen: dict[str, dict] = {}   # region_full -> record

        for col in thick_cols:
            lvl_match = re.search(r'lvl(\d+)', col, re.IGNORECASE)
            lvl = int(lvl_match.group(1)) if lvl_match else 1

            for region in df.index:
                rs = str(region)
                if regions_of_interest and not any(r in rs for r in regions_of_interest):
                    continue
                val = df.loc[region, col]
                if pd.isna(val):
                    continue

                if rs in seen:
                    continue   # lowest level already recorded

                if rs.lower().startswith("l_"):
                    hemi, base = "left",  rs[2:]
                elif rs.lower().startswith("r_"):
                    hemi, base = "right", rs[2:]
                else:
                    hemi, base = "bilateral", rs

                seen[rs] = dict(
                    subject=sub, session=ses,
                    region_full=rs, region=base,
                    hemisphere=hemi, atlas_level=lvl,
                    thickness_mm=float(val),
                )

        records = list(seen.values())
    except Exception as e:
        warnings.warn(f"Error reading thickness {path}: {e}")
    return records


def process_thickness(thickness_paths_dict, regions_of_interest=None):
    """
    Build tidy DataFrame from pre-computed thickness.xlsx files.
    One row per (subject, session, region_full): lowest atlas level wins.
    Uses parallel reading when joblib is available.
    """
    triples = list(zip(
        thickness_paths_dict["subject"],
        thickness_paths_dict["session"],
        thickness_paths_dict["thickness_path"],
    ))

    if _JOBLIB:
        results = Parallel(n_jobs=-1)(
            delayed(_read_one_thickness)(sub, ses, path, regions_of_interest)
            for sub, ses, path in triples
        )
    else:
        results = [_read_one_thickness(sub, ses, path, regions_of_interest)
                   for sub, ses, path in triples]

    records = [r for sublist in results for r in sublist]
    return pd.DataFrame(records)


# ?????????????????????????????????????????????????????????????????????????????
# FIX 2 ? xlsx persistence: delete old file, write fresh
# ?????????????????????????????????????????????????????????????????????????????

def _save_surface_xlsx(results_by_level: dict[int, dict], ses_dir: str):
    """
    Persist surface results to Native_resol/surface.xlsx.

    Parameters
    ----------
    results_by_level : {atlas_level: {region_name: value_mm2}}
        All levels to write, collected by the caller before this function
        is invoked.  The file is written exactly once per subject/session.

    The old file is deleted first so no stale columns survive.
    Columns: Surface_Area_lvl1, Surface_Area_lvl2, ? (one per level present).
    """
    if not results_by_level:
        return
    out_dir = opj(ses_dir, "anat", "native", "surfaces", "Native_resol")
    os.makedirs(out_dir, exist_ok=True)
    out = opj(out_dir, "surface.xlsx")

    # Delete old file unconditionally
    if os.path.exists(out):
        os.remove(out)

    frames = {}
    for lvl in sorted(results_by_level):
        col = f"Surface_Area_lvl{lvl}"
        frames[col] = pd.Series(results_by_level[lvl], name=col)

    df = pd.concat(frames.values(), axis=1)
    df.index.name = "Region"
    df.to_excel(out)


def _save_thickness_xlsx(results_by_level: dict[int, dict], ses_dir: str):
    """
    Persist thickness results to Native_resol/thickness.xlsx.

    Parameters
    ----------
    results_by_level : {atlas_level: {region_name: value_mm}}

    The old file is deleted first so no stale columns survive.
    Columns: Thickness_mm_lvl1, Thickness_mm_lvl2, ? (one per level present).
    """
    if not results_by_level:
        return
    out_dir = opj(ses_dir, "anat", "native", "surfaces", "Native_resol")
    os.makedirs(out_dir, exist_ok=True)
    out = opj(out_dir, "thickness.xlsx")

    if os.path.exists(out):
        os.remove(out)

    frames = {}
    for lvl in sorted(results_by_level):
        col = f"Thickness_mm_lvl{lvl}"
        frames[col] = pd.Series(results_by_level[lvl], name=col)

    df = pd.concat(frames.values(), axis=1)
    df.index.name = "Region"
    df.to_excel(out)


# ?????????????????????????????????????????????????????????????????????????????
# wb_command helpers (unchanged except _save_* call-sites below)
# ?????????????????????????????????????????????????????????????????????????????

def _wb_get_region_names(sing_wb, seg_file, map_number=1):
    import hashlib
    seg_hash = hashlib.md5(seg_file.encode()).hexdigest()[:8]
    tmp = f"/tmp/ednix_label_table_{seg_hash}_m{map_number}.txt"
    cmd = f'{sing_wb} wb_command -cifti-label-export-table "{seg_file}" {map_number} "{tmp}"'
    _spco(cmd)
    names = []
    if os.path.exists(tmp):
        with open(tmp) as f:
            for line in f:
                line = line.strip()
                if line and (line.startswith("l_") or line.startswith("r_")
                             or line.startswith("L_") or line.startswith("R_")):
                    names.append(line)
        os.remove(tmp)
    return names


def _filter_regions(all_regions, patterns):
    out = []
    for pat in patterns:
        if pat.startswith(("l_", "r_", "L_", "R_")):
            if pat in all_regions:
                out.append(pat)
        else:
            out.extend(r for r in all_regions if pat in r)
    return out


_LEGACY_SUFFIXES = (
    "_rois.dscalar.nii",
    "_rois.shape.gii",
    "_midthickness_shape.gii",)


def _cleanup_legacy_roi_files(parent_roi_dir: str, region_name: str):
    dirs_to_check = [
        parent_roi_dir,
        opj(parent_roi_dir, "surfaces"),
        opj(parent_roi_dir, "thickness"),
    ]
    for d in dirs_to_check:
        if not os.path.isdir(d):
            continue
        for suffix in _LEGACY_SUFFIXES:
            legacy = opj(d, region_name + suffix)
            if os.path.exists(legacy):
                os.remove(legacy)


def _wb_extract_one_surface(sing_wb, seg_file, region_name, subject_id, ses_dir, map_number=1):
    native_dir     = opj(ses_dir, "anat", "native", "surfaces", "Native_resol")
    parent_roi_dir = opj(native_dir, "ROIs")
    roi_dir        = opj(parent_roi_dir, "surfaces")
    os.makedirs(roi_dir, exist_ok=True)

    _cleanup_legacy_roi_files(parent_roi_dir, region_name)

    prefix = region_name[0].lower()
    H_SIDE = "LEFT" if prefix == "l" else "RIGHT"
    surf   = opj(native_dir, f"{subject_id}.{prefix}.pial.surf.gii")
    if not os.path.exists(surf):
        return None

    roi_f   = opj(roi_dir, f"{region_name}_m{map_number}_roi.dscalar.nii")
    shape_f = opj(roi_dir, f"{region_name}_m{map_number}_roi.shape.gii")
    area_f  = opj(roi_dir, f"{region_name}_m{map_number}_pial_area.gii")

    steps = [
        (roi_f,
         f'{sing_wb} wb_command -cifti-label-to-roi "{seg_file}" "{roi_f}"'
         f' -map {map_number} -name "{region_name}"'),
        (shape_f,
         f'{sing_wb} wb_command -cifti-separate "{roi_f}" COLUMN'
         f' -metric CORTEX_{H_SIDE} "{shape_f}"'),
        (area_f,
         f'{sing_wb} wb_command -surface-vertex-areas "{surf}" "{area_f}"'),
    ]
    for out_path, cmd in steps:
        if os.path.exists(out_path):
            os.remove(out_path)
        if _spco_timed(cmd, _TIMEOUT_SHORT) is None or not os.path.exists(out_path):
            return None

    data = _spco_timed(
        f'{sing_wb} wb_command -metric-stats "{area_f}" -reduce SUM -roi "{shape_f}"',
        _TIMEOUT_SHORT,
    )
    try:
        return float(data.strip()) if data else None
    except (ValueError, AttributeError):
        return None


def _average_hemispheres(results: dict) -> dict:
    averaged, done = {}, set()
    for name, val in results.items():
        if name in done:
            continue
        base = name[2:]
        opp  = ("r_" if name.startswith("l_") else "l_") + base
        if opp in results:
            averaged[base] = (val + results[opp]) / 2
            done.update([name, opp])
        else:
            averaged[name] = val
            warnings.warn(f"Only one hemisphere for {name}")
            done.add(name)
    return averaged


# ?????????????????????????????????????????????????????????????????????????????
# Thickness
# ?????????????????????????????????????????????????????????????????????????????

import traceback

try:
    from joblib import Parallel, delayed
    _JOBLIB = True
except ImportError:
    _JOBLIB = False

_TIMEOUT_SHORT = 120
_TIMEOUT_HEAVY = 600


def _spco_timed(cmd: str, timeout: int = _TIMEOUT_SHORT):
    import subprocess
    try:
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, timeout=timeout)
        if result.returncode != 0:
            warnings.warn(f"Command failed: {cmd}\n{result.stderr.strip()}")
            return None
        return result.stdout
    except subprocess.TimeoutExpired:
        warnings.warn(f"Command timed out after {timeout}s: {cmd}")
        return None


def _wb_compute_thickness_cifti(sing_wb, sub, ses_dir, native_dir, overwrite=True):
    dscalar_out = opj(native_dir, f"{sub}.thicknesswb.dscalar.nii")

    if os.path.exists(dscalar_out) and not overwrite:
        shape_files_exist = all(
            os.path.exists(opj(native_dir, f"{sub}.{h}.thicknesswb.shape.gii"))
            for h in ("l", "r")
        )
        if shape_files_exist:
            print(f"  [thickness_wb] {sub}: reusing cached dscalar")
            return dscalar_out
        else:
            warnings.warn(
                f"  [thickness_wb] {sub}: cached dscalar found but shape files missing "
                f"? recomputing")

    palette = (
        ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 4 96'
        ' -interpolate true -palette-name videen_style'
        ' -disp-pos true -disp-neg false -disp-zero false'
    )

    shape_files = {}

    for h_short in ("l", "r"):
        pial_f  = opj(native_dir, f"{sub}.{h_short}.pial.native.surf.gii")
        white_f = opj(native_dir, f"{sub}.{h_short}.white.native.surf.gii")
        roi_f   = opj(native_dir, f"{sub}.{h_short}.roi.shape.gii")

        if not os.path.exists(pial_f):
            pial_f  = opj(native_dir, f"{sub}.{h_short}.pial.surf.gii")
        if not os.path.exists(white_f):
            white_f = opj(native_dir, f"{sub}.{h_short}.white.surf.gii")

        if not os.path.exists(pial_f) or not os.path.exists(white_f):
            warnings.warn(
                f"  [thickness_wb] {sub}: surface not found hemi={h_short} "
                f"(pial={os.path.exists(pial_f)} white={os.path.exists(white_f)})")
            return None

        thick_shape = opj(native_dir, f"{sub}.{h_short}.thicknesswb.shape.gii")

        if os.path.exists(thick_shape):
            os.remove(thick_shape)

        cmd1 = (f'{sing_wb} wb_command -surface-to-surface-3d-distance'
                f' "{pial_f}" "{white_f}" "{thick_shape}"')
        if _spco_timed(cmd1, _TIMEOUT_HEAVY) is None or not os.path.exists(thick_shape):
            warnings.warn(f"  [thickness_wb] {sub}: 3d-distance failed hemi={h_short}")
            return None

        if os.path.exists(roi_f):
            thick_shape_tmp = thick_shape + ".tmp.shape.gii"
            if os.path.exists(thick_shape_tmp):
                os.remove(thick_shape_tmp)

            cmd2 = (f'{sing_wb} wb_command -metric-math "roi * thickness"'
                    f' "{thick_shape_tmp}"'
                    f' -var roi "{roi_f}"'
                    f' -var thickness "{thick_shape}"')

            if _spco_timed(cmd2, _TIMEOUT_SHORT) is not None and os.path.exists(thick_shape_tmp):
                os.replace(thick_shape_tmp, thick_shape)
                print(f"  [thickness_wb] {sub}: ROI mask applied hemi={h_short}")
            else:
                warnings.warn(
                    f"  [thickness_wb] {sub}: ROI mask failed hemi={h_short} "
                    f"? keeping unmasked thickness")
                if os.path.exists(thick_shape_tmp):
                    os.remove(thick_shape_tmp)
        else:
            warnings.warn(
                f"  [thickness_wb] {sub}: no ROI file found hemi={h_short} "
                f"({roi_f}) ? thickness will NOT be masked")

        cmd3 = (f'{sing_wb} wb_command -metric-palette'
                f' "{thick_shape}"{palette}')
        _spco_timed(cmd3, _TIMEOUT_SHORT)

        shape_files[h_short] = thick_shape

    if len(shape_files) < 2:
        warnings.warn(f"  [thickness_wb] {sub}: missing hemisphere shape files")
        return None

    if os.path.exists(dscalar_out):
        os.remove(dscalar_out)

    roi_l = opj(native_dir, f"{sub}.l.roi.shape.gii")
    roi_r = opj(native_dir, f"{sub}.r.roi.shape.gii")

    cmd4 = (f'{sing_wb} wb_command -cifti-create-dense-scalar'
            f' "{dscalar_out}"'
            f' -left-metric  "{shape_files["l"]}"')
    if os.path.exists(roi_l):
        cmd4 += f' -roi-left "{roi_l}"'
    cmd4 += f' -right-metric "{shape_files["r"]}"'
    if os.path.exists(roi_r):
        cmd4 += f' -roi-right "{roi_r}"'

    if _spco_timed(cmd4, _TIMEOUT_HEAVY) is None or not os.path.exists(dscalar_out):
        warnings.warn(f"  [thickness_wb] {sub}: cifti-create-dense-scalar failed")
        return None

    cmd5 = (f'{sing_wb} wb_command -set-map-names "{dscalar_out}"'
            f' -map 1 {sub}_thickness_from_wb')
    _spco_timed(cmd5, _TIMEOUT_SHORT)

    cmd6 = (f'{sing_wb} wb_command -cifti-palette "{dscalar_out}"'
            f'{palette} "{dscalar_out}"')
    _spco_timed(cmd6, _TIMEOUT_HEAVY)

    print(f"  [thickness_wb] {sub}: dscalar ? {os.path.basename(dscalar_out)}")
    return dscalar_out


def _wb_extract_one_thickness_from_cifti(
        sing_wb, dscalar_f, seg_file, region_name, ses_dir, map_number=1):
    native_dir     = opj(ses_dir, "anat", "native", "surfaces", "Native_resol")
    parent_roi_dir = opj(native_dir, "ROIs")
    roi_dir        = opj(parent_roi_dir, "thickness")
    os.makedirs(roi_dir, exist_ok=True)

    _cleanup_legacy_roi_files(parent_roi_dir, region_name)

    prefix  = region_name[0].lower()
    H_SIDE  = "LEFT" if prefix == "l" else "RIGHT"

    roi_cifti   = opj(roi_dir, f"{region_name}_m{map_number}_roi.dscalar.nii")
    roi_shape   = opj(roi_dir, f"{region_name}_m{map_number}_roi.shape.gii")
    thick_shape = opj(roi_dir, f"{region_name}_m{map_number}_thick.shape.gii")

    steps = [
        (roi_cifti,
         f'{sing_wb} wb_command -cifti-label-to-roi "{seg_file}" "{roi_cifti}"'
         f' -map {map_number} -name "{region_name}"'),
        (roi_shape,
         f'{sing_wb} wb_command -cifti-separate "{roi_cifti}" COLUMN'
         f' -metric CORTEX_{H_SIDE} "{roi_shape}"'),
        (thick_shape,
         f'{sing_wb} wb_command -cifti-separate "{dscalar_f}" COLUMN'
         f' -metric CORTEX_{H_SIDE} "{thick_shape}"'),
    ]
    for out_path, cmd in steps:
        if os.path.exists(out_path):
            os.remove(out_path)
        if _spco_timed(cmd, _TIMEOUT_SHORT) is None or not os.path.exists(out_path):
            return None

    data = _spco_timed(
        f'{sing_wb} wb_command -metric-stats "{thick_shape}"'
        f' -reduce MEAN -roi "{roi_shape}"',
        _TIMEOUT_SHORT,
    )
    try:
        return float(data.strip()) if data else None
    except (ValueError, AttributeError):
        return None


# ?????????????????????????????????????????????????????????????????????????????
# Shared helpers for parallel subject processing
# ?????????????????????????????????????????????????????????????????????????????

def _save_failed(failed: list, out_dir: str, modality: str):
    if not failed:
        return
    print(f"\n{'!'*60}")
    print(f"  [{modality}] {len(failed)} subject(s) FAILED:")
    for entry in failed:
        print(f"    ? {entry['subject']}  (ses={entry['session']})  ? {entry['reason']}")
    print(f"{'!'*60}\n")

    os.makedirs(out_dir, exist_ok=True)
    log_path = opj(out_dir, f"failed_subjects_{modality}.txt")
    with open(log_path, "w") as fh:
        fh.write(f"# Failed subjects ? {modality}\n")
        for entry in failed:
            fh.write(f"{entry['subject']}\t{entry['session']}\t{entry['reason']}\n")
    print(f"  [{modality}] failure log ? {log_path}")


def _build_records(results: dict, sub: str, ses_label: str,
                   atlas_level: int, value_key: str) -> list:
    records = []
    for region_name, val in results.items():
        if region_name.lower().startswith("l_"):
            hemi, base = "left",  region_name[2:]
        elif region_name.lower().startswith("r_"):
            hemi, base = "right", region_name[2:]
        else:
            hemi, base = "bilateral", region_name
        records.append(dict(
            subject=sub, session=ses_label,
            region_full=region_name, region=base,
            hemisphere=hemi, atlas_level=atlas_level,
            **{value_key: val},
        ))
    return records


def _iter_subjects(bids_dir, list_to_keep, list_to_remove):
    for sub, ses, _, ses_dir in _scan_bids(bids_dir, "sub-*"):
        ses_label = ses or "1"
        if list_to_keep   and (sub, ses_label) not in list_to_keep:   continue
        if list_to_remove and (sub, ses_label) in list_to_remove:     continue
        native_dir = opj(ses_dir, "anat", "native", "surfaces", "Native_resol")
        if not os.path.exists(native_dir):
            continue
        yield sub, ses_label, ses_dir, native_dir


def _get_seg_file(native_dir: str, sub: str):
    seg_files = glob.glob(opj(native_dir, f"{sub}.EDNIxCSC*.dlabel.nii"))
    if not seg_files:
        warnings.warn(f"No dlabel segmentation for {sub}")
        return None
    return next((s for s in seg_files if "LR" not in s), seg_files[0])


# ?????????????????????????????????????????????????????????????????????????????
# Per-subject workers  ?  accumulate ALL levels, then write xlsx ONCE
# ?????????????????????????????????????????????????????????????????????????????

def _process_one_subject_surface(
    sing_wb, sub, ses_label, ses_dir, native_dir,
    regions_to_process, should_average_hemispheres,
    n_atlas_levels, overwrite,
):
    """
    Extract surface area for every atlas level of one subject.

    All levels are accumulated into results_by_level before the xlsx is
    written so the file is created exactly once (not once per level).
    """
    try:
        seg_file = _get_seg_file(native_dir, sub)
        if seg_file is None:
            return [], {"subject": sub, "session": ses_label, "reason": "no dlabel seg file"}

        all_records:    list       = []
        results_by_level: dict     = {}   # {atlas_level: {region: value}}

        for atlas_level in range(1, n_atlas_levels + 1):
            all_regions = _wb_get_region_names(sing_wb, seg_file, atlas_level)
            if not all_regions:
                continue
            regions = (_filter_regions(all_regions, regions_to_process)
                       if regions_to_process else all_regions)

            if _JOBLIB and len(regions) >= 8:
                areas = Parallel(n_jobs=-1, prefer="threads")(
                    delayed(_wb_extract_one_surface)(
                        sing_wb, seg_file, r, sub, ses_dir, atlas_level)
                    for r in regions
                )
                results = {r: a for r, a in zip(regions, areas) if a is not None}
            else:
                results = {
                    r: a for r in regions
                    if (a := _wb_extract_one_surface(
                        sing_wb, seg_file, r, sub, ses_dir, atlas_level)) is not None
                }

            if should_average_hemispheres:
                results = _average_hemispheres(results)

            results_by_level[atlas_level] = results
            all_records.extend(
                _build_records(results, sub, ses_label, atlas_level, "surface_area_mm2"))

        # Write xlsx once with all levels (old file deleted inside)
        _save_surface_xlsx(results_by_level, ses_dir)
        return all_records, None

    except Exception as exc:
        reason = f"{type(exc).__name__}: {exc}"
        warnings.warn(f"  [surface_wb] {sub} ses-{ses_label} FAILED ? {reason}\n"
                      + traceback.format_exc())
        return [], {"subject": sub, "session": ses_label, "reason": reason}


def _process_one_subject_thickness(
    sing_wb, sub, ses_label, ses_dir, native_dir,
    regions_to_process, n_atlas_levels, overwrite,
):
    """
    Compute & extract thickness for every atlas level of one subject.

    All levels are accumulated into results_by_level before the xlsx is
    written so the file is created exactly once (not once per level).
    """
    try:
        seg_file = _get_seg_file(native_dir, sub)
        if seg_file is None:
            return [], {"subject": sub, "session": ses_label, "reason": "no dlabel seg file"}

        dscalar_f = _wb_compute_thickness_cifti(
            sing_wb, sub, ses_dir, native_dir, overwrite=overwrite)
        if dscalar_f is None:
            reason = "could not compute thickness dscalar (check pial/white surfaces)"
            warnings.warn(f"  [thickness_wb] {sub}: {reason}")
            return [], {"subject": sub, "session": ses_label, "reason": reason}

        all_records:    list  = []
        results_by_level: dict = {}   # {atlas_level: {region: value}}

        for atlas_level in range(1, n_atlas_levels + 1):
            all_regions = _wb_get_region_names(sing_wb, seg_file, atlas_level)
            if not all_regions:
                continue
            regions = (_filter_regions(all_regions, regions_to_process)
                       if regions_to_process else all_regions)

            if _JOBLIB and len(regions) >= 8:
                thicks = Parallel(n_jobs=-1, prefer="threads")(
                    delayed(_wb_extract_one_thickness_from_cifti)(
                        sing_wb, dscalar_f, seg_file, r, ses_dir, atlas_level)
                    for r in regions
                )
                results = {r: t for r, t in zip(regions, thicks) if t is not None}
            else:
                results = {
                    r: t for r in regions
                    if (t := _wb_extract_one_thickness_from_cifti(
                        sing_wb, dscalar_f, seg_file, r, ses_dir, atlas_level)) is not None
                }

            results_by_level[atlas_level] = results
            all_records.extend(
                _build_records(results, sub, ses_label, atlas_level, "thickness_mm"))

        # Write xlsx once with all levels (old file deleted inside)
        _save_thickness_xlsx(results_by_level, ses_dir)
        return all_records, None

    except Exception as exc:
        reason = f"{type(exc).__name__}: {exc}"
        warnings.warn(f"  [thickness_wb] {sub} ses-{ses_label} FAILED ? {reason}\n"
                      + traceback.format_exc())
        return [], {"subject": sub, "session": ses_label, "reason": reason}


# ?????????????????????????????????????????????????????????????????????????????
# Public extraction API
# ?????????????????????????????????????????????????????????????????????????????

def extract_and_process_surfaces_wb(
    sing_wb,
    bids_dir,
    list_to_keep=None,
    list_to_remove=None,
    regions_to_process=None,
    should_average_hemispheres=False,
    n_atlas_levels=N_ATLAS_LEVELS,
    overwrite=True,
):
    subjects = list(_iter_subjects(bids_dir, list_to_keep or [], list_to_remove or []))

    if _JOBLIB and len(subjects) > 1:
        results = Parallel(n_jobs=-1, backend="loky", verbose=5)(
            delayed(_process_one_subject_surface)(
                sing_wb, sub, ses_label, ses_dir, native_dir,
                regions_to_process, should_average_hemispheres,
                n_atlas_levels, overwrite,
            )
            for sub, ses_label, ses_dir, native_dir in subjects
        )
    else:
        results = [
            _process_one_subject_surface(
                sing_wb, sub, ses_label, ses_dir, native_dir,
                regions_to_process, should_average_hemispheres,
                n_atlas_levels, overwrite,
            )
            for sub, ses_label, ses_dir, native_dir in subjects
        ]

    all_records, failed = [], []
    for records, failure in results:
        all_records.extend(records)
        if failure:
            failed.append(failure)

    _save_failed(failed, out_dir=bids_dir, modality="surface")
    return pd.DataFrame(all_records)


def extract_and_process_thickness_wb(
    sing_wb,
    bids_dir,
    list_to_keep=None,
    list_to_remove=None,
    regions_to_process=None,
    n_atlas_levels=N_ATLAS_LEVELS,
    overwrite=True,
):
    subjects = list(_iter_subjects(bids_dir, list_to_keep or [], list_to_remove or []))

    if _JOBLIB and len(subjects) > 1:
        results = Parallel(n_jobs=-1, backend="loky", verbose=5)(
            delayed(_process_one_subject_thickness)(
                sing_wb, sub, ses_label, ses_dir, native_dir,
                regions_to_process, n_atlas_levels, overwrite,
            )
            for sub, ses_label, ses_dir, native_dir in subjects
        )
    else:
        results = [
            _process_one_subject_thickness(
                sing_wb, sub, ses_label, ses_dir, native_dir,
                regions_to_process, n_atlas_levels, overwrite,
            )
            for sub, ses_label, ses_dir, native_dir in subjects
        ]

    all_records, failed = [], []
    for records, failure in results:
        all_records.extend(records)
        if failure:
            failed.append(failure)

    _save_failed(failed, out_dir=bids_dir, modality="thickness")
    return pd.DataFrame(all_records)


# ??????????????????????????????????????????????????????????????????????????????
# Volumes  ?  FIX 3: no atlas_level column
# ??????????????????????????????????????????????????????????????????????????????

def process_volumes(volume_paths_dict, label_path, regions_of_interest=None,
                    voxel_volume_mm3=None):
    """
    Volumes have no atlas-level concept.
    The returned DataFrame does NOT contain an atlas_level column.
    """
    try:
        import nibabel as nib
    except ImportError:
        raise ImportError("nibabel is required: pip install nibabel")

    label_df     = parse_label_file(label_path)
    label_lookup = {row["label_id"]: row for _, row in label_df.iterrows()}

    records = []
    for sub, ses, seg_path in zip(
        volume_paths_dict["subject"],
        volume_paths_dict["session"],
        volume_paths_dict["volume_seg_path"],
    ):
        if not seg_path or not os.path.exists(seg_path):
            warnings.warn(f"Segmentation not found: {seg_path}"); continue
        try:
            img  = nib.load(seg_path)
            data = np.asarray(img.dataobj, dtype=np.int32)
            if voxel_volume_mm3 is None:
                zooms  = img.header.get_zooms()[:3]
                vox_mm = float(np.prod(zooms))
            else:
                vox_mm = float(voxel_volume_mm3)

            unique_ids, counts = np.unique(data[data > 0], return_counts=True)
            for lid, cnt in zip(unique_ids, counts):
                if lid not in label_lookup:
                    continue
                row  = label_lookup[lid]
                base = row["base_region"]
                if regions_of_interest and not any(r in base for r in regions_of_interest):
                    continue
                records.append(dict(
                    subject=sub, session=ses,
                    label_id=int(lid),
                    region_name=row["region_name"],
                    region=base,
                    hemisphere=row["hemisphere"],
                    voxel_count=int(cnt),
                    volume_mm3=float(cnt * vox_mm),
                    # NOTE: no atlas_level ? volumes don't have hierarchy levels
                ))
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
            pivot = grp.set_index("region_name")[["region", "hemisphere",
                                                   "voxel_count", "volume_mm3"]]
            pivot.to_excel(out_xlsx)

    return df_all


# ?????????????????????????????????????????????????????????????????????????????
# QC  (unchanged)
# ?????????????????????????????????????????????????????????????????????????????

def _load_json_with_nan(path):
    with open(path) as f:
        raw = f.read()
    cleaned = re.sub(r"\bNaN\b", "null", raw)
    data    = json.loads(cleaned)

    def _restore(v):
        if v is None:            return np.nan
        if isinstance(v, dict):  return {k: _restore(vv) for k, vv in v.items()}
        if isinstance(v, list):  return [_restore(i) for i in v]
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
                    if isinstance(item, dict):
                        out.update(_flatten(item, f"{key}_{i}", sep))
                    else:
                        out[f"{key}_{i}"] = item
        else:
            out[key] = v
    return out


def _parse_full_results(path):
    try:
        data = _load_json_with_nan(path)
    except Exception as e:
        warnings.warn(f"full_results parse error {path}: {e}")
        return {}

    out = {}

    nm = data.get("network_metrics", {})
    key_map = {
        "Mean_Correlation":    "net_mean_correlation",
        "Std_Correlation":     "net_std_correlation",
        "Top_Eigenvalue":      "net_top_eigenvalue",
        "Eigenvalue_Ratio":    "net_eigenvalue_ratio",
        "Silhouette_Score":    "net_silhouette",
        "Davies_Bouldin":      "net_davies_bouldin",
        "Calinski_Harabasz":   "net_calinski_harabasz",
        "Median_Correlation":  "net_median_correlation",
    }
    for json_k, col_k in key_map.items():
        if json_k in nm:
            out[col_k] = nm[json_k]

    hr = data.get("hemisphere_results", {})
    for k in ("intra_left_mean", "intra_right_mean", "inter_mean",
              "p_intra_vs_inter", "p_left_vs_right",
              "t_intra_vs_inter", "t_left_vs_right"):
        if k in hr:
            out[f"net_{k}"] = hr[k]

    if "intra_left_mean" in hr and "intra_right_mean" in hr:
        out["net_intra_mean"] = (hr["intra_left_mean"] + hr["intra_right_mean"]) / 2

    for side in ("intra_left", "intra_right", "inter"):
        vals = hr.get(side, [])
        if isinstance(vals, list) and len(vals) > 0:
            arr = [v for v in vals
                   if v is not None and not (isinstance(v, float) and np.isnan(v))]
            if arr:
                out[f"net_{side}_pct_pos"] = 100.0 * sum(v > 0 for v in arr) / len(arr)
                out[f"net_{side}_pct_neg"] = 100.0 * sum(v < 0 for v in arr) / len(arr)

    sr = data.get("specificity_results", {})
    ts = sr.get("target_specificity", {})
    for k in ("Specific_Correlation", "NonSpecific_Correlation"):
        if k in ts:
            out[f"sp_{k.lower()}"] = ts[k]
    if "Category" in ts:
        cat = str(ts["Category"]).strip()
        out["sp_category"] = cat
        enc = {"Specific": 1.0, "Unspecific": 0.5, "Spurious": 0.25, "No": 0.0}
        out["sp_specific"] = enc.get(cat, np.nan)

    return out


def process_qc(qc_paths_dict):
    records  = []
    iters = zip(
        qc_paths_dict["subject"],
        qc_paths_dict["session"],
        qc_paths_dict["func_qc_path"],
        qc_paths_dict["anat_qc_path"],
        qc_paths_dict.get("full_results_path", [None] * len(qc_paths_dict["subject"])),
    )

    for sub, ses, func_p, anat_p, full_p in iters:
        row = {"subject": sub, "session": ses}

        if func_p and os.path.exists(func_p):
            try:
                row.update({f"func_{k}": v
                             for k, v in _flatten(_load_json_with_nan(func_p)).items()})
            except Exception as e:
                warnings.warn(f"func QC error {func_p}: {e}")

        if anat_p and os.path.exists(anat_p):
            try:
                row.update({f"anat_{k}": v
                             for k, v in _flatten(_load_json_with_nan(anat_p)).items()})
            except Exception as e:
                warnings.warn(f"anat QC error {anat_p}: {e}")

        if full_p and os.path.exists(full_p):
            row.update(_parse_full_results(full_p))

        records.append(row)

    df = pd.DataFrame(records)
    id_cols   = [c for c in ('subject', 'session', 'species', 'bids_dir') if c in df.columns]
    data_cols = [c for c in df.columns if c not in id_cols]
    if data_cols:
        df = df.dropna(subset=data_cols, how='all')
    return df


# ?????????????????????????????????????????????????????????????????????????????
# Correlation matrix helpers (unchanged)
# ?????????????????????????????????????????????????????????????????????????????

def extract_corr_matrix_paths(
    bids_dir,
    atlas_name="EDNIxCSC",
    atlas_level=3,
    use_lr=False,
    list_to_keep=None,
    list_to_remove=None,
):
    list_to_keep   = list_to_keep   or []
    list_to_remove = list_to_remove or []

    suffix       = "LR" if use_lr else ""
    pattern_name = f"{atlas_name}{suffix}_{atlas_level}_run_*_matrix.csv"

    records = []
    for f in sorted(glob.glob(
        opj(bids_dir, "sub-*", "ses-*", "func", "acpc-func",
            "Stats", "Correl_matrix", pattern_name)
    )):
        bn = os.path.basename(f)
        if any(x in bn for x in ('check_fit', 'flattened', 'pval', 'tstat')):
            continue
        f   = _linux_path(f)
        sub = _parse_bids_entity(f, "sub")
        ses = _parse_bids_entity(f, "ses") or "1"
        if not sub:
            continue
        if list_to_keep   and (sub, ses) not in list_to_keep:   continue
        if list_to_remove and (sub, ses) in list_to_remove:     continue
        m   = re.search(r"run_(\d+)_matrix", os.path.basename(f))
        run = int(m.group(1)) if m else 0
        records.append(dict(subject=sub, session=ses, run=run, path=f))

    print(f"  [corr_matrix] {len(records)} files  ? {bids_dir}  "
          f"(atlas={atlas_name}{suffix} lvl={atlas_level})")
    return records


def load_corr_matrix(path):
    df = pd.read_csv(path, index_col=0)
    col_names = list(df.columns)
    idx_names = list(df.index)

    def _looks_like_rois(names):
        return any(not str(n).lstrip('-').replace('.', '').isdigit() for n in names)

    if _looks_like_rois(idx_names):
        return idx_names, df.values.astype(float)
    elif _looks_like_rois(col_names):
        return col_names, df.T.values.astype(float)
    else:
        return idx_names, df.values.astype(float)


def collect_corr_matrices(
    species_config,
    atlas_name="EDNIxCSC",
    atlas_level=3,
    use_lr=False,
    save_csv_dir=None,
):
    from scipy import stats as _stats

    species_stacks = {}
    species_rois   = {}

    for species, cfg in species_config.items():
        bids_dirs = cfg.get("bids_dirs", [])
        if isinstance(bids_dirs, str):
            bids_dirs = [bids_dirs]
        lk        = cfg.get("list_to_keep",   [])
        lr_remove = cfg.get("list_to_remove", [])

        sub_ses_runs = {}
        rois_ref     = None
        n_loaded     = 0
        n_skipped    = 0

        for bids_dir in bids_dirs:
            records = extract_corr_matrix_paths(
                bids_dir, atlas_name, atlas_level, use_lr, lk, lr_remove)
            for rec in records:
                try:
                    rois, mat = load_corr_matrix(rec["path"])
                    valid = mat[~np.isnan(mat)]
                    if valid.size == 0 or np.all(valid == 0):
                        print(f"    [SKIP] all-NaN/zero matrix: {rec['path']}")
                        n_skipped += 1
                        continue
                    if rois_ref is None:
                        rois_ref = rois
                    key = (rec["subject"], rec["session"])
                    sub_ses_runs.setdefault(key, []).append((rois, mat))
                    n_loaded += 1
                except Exception as e:
                    print(f"    [ERROR] {rec['path']}: {e}")
                    n_skipped += 1

        print(f"  [{species}] loaded={n_loaded}  skipped={n_skipped}  "
              f"unique subjects={len(sub_ses_runs)}")

        if not sub_ses_runs:
            continue

        subject_means = []
        for (sub, ses), run_list in sorted(sub_ses_runs.items()):
            if not run_list:
                continue
            rois_sets = [set(r) for r, _ in run_list]
            common    = rois_sets[0].intersection(*rois_sets[1:])
            ref_rois  = [r for r in run_list[0][0] if r in common]
            if not ref_rois:
                continue
            mats = []
            for rois_r, mat_r in run_list:
                try:
                    idx = [list(rois_r).index(r) for r in ref_rois]
                    mats.append(mat_r[np.ix_(idx, idx)])
                except ValueError:
                    pass
            if mats:
                subject_means.append((ref_rois, np.nanmean(np.stack(mats, axis=0), axis=0)))

        if not subject_means:
            continue

        all_roi_sets = [set(r) for r, _ in subject_means]
        common_sp    = all_roi_sets[0].intersection(*all_roi_sets[1:])
        rois_ref     = [r for r in subject_means[0][0] if r in common_sp]

        aligned = []
        for rois_s, mat_s in subject_means:
            idx = [list(rois_s).index(r) for r in rois_ref if r in rois_s]
            if len(idx) == len(rois_ref):
                aligned.append(mat_s[np.ix_(idx, idx)])

        if not aligned:
            continue

        species_stacks[species] = aligned
        species_rois[species]   = rois_ref

    if not species_stacks:
        return {}

    common_rois = None
    for sp, rois in species_rois.items():
        rois_set = set(rois)
        common_rois = rois_set if common_rois is None else common_rois & rois_set
    first_rois  = species_rois[next(iter(species_rois))]
    common_rois = [r for r in first_rois if r in common_rois]

    result = {}
    for species, mats in species_stacks.items():
        sp_rois    = species_rois[species]
        local_rois = [r for r in common_rois if r in sp_rois]
        local_idx  = [sp_rois.index(r) for r in local_rois]

        stack = np.stack([m[np.ix_(local_idx, local_idx)] for m in mats], axis=0)
        n     = stack.shape[0]
        mean  = np.nanmean(stack, axis=0)
        var   = np.nanvar(stack, axis=0)

        if n >= 2:
            _, p = _stats.ttest_1samp(stack, 0, axis=0, nan_policy='omit')
        else:
            p = np.full(mean.shape, np.nan)

        result[species] = dict(mean=mean, var=var, pval=p, rois=local_rois, n=n)

        if save_csv_dir:
            os.makedirs(save_csv_dir, exist_ok=True)
            for mat, name in [(mean, 'mean'), (var, 'var'), (p, 'pval')]:
                df_out = pd.DataFrame(mat, index=local_rois, columns=local_rois)
                df_out.to_csv(opj(save_csv_dir, f"{species}_corr_{name}.csv"))

    return result


# ??????????????????????????????????????????????????????????????????????????????
# PART 3  ?  MULTI-BIDS / MULTI-SPECIES COLLECTION
# ??????????????????????????????????????????????????????????????????????????????

def collect_multi_species(
    species_config: dict,
    regions_of_interest=None,
    extract=("surface", "volume", "thickness", "qc"),
    atlas_name="EDNIxCSC",
    atlas_label_paths: dict = None,
    atlas_library_root: str = None,
    species_atlas_fragments: dict = None,
):
    atlas_label_paths       = atlas_label_paths       or {}
    species_atlas_fragments = species_atlas_fragments or {}

    accumulators   = {m: [] for m in ("surface", "volume", "thickness", "qc")}
    bold_paths_all = {}

    for species, cfg in species_config.items():
        bids_dirs = cfg.get("bids_dirs", [])
        if isinstance(bids_dirs, str):
            bids_dirs = [bids_dirs]
        lk = cfg.get("list_to_keep",   [])
        lr = cfg.get("list_to_remove", [])

        print(f"\n{'='*64}")
        print(f"  Species : {species}  ({len(bids_dirs)} BIDS dir(s))")
        print(f"{'='*64}")

        label_path = atlas_label_paths.get(species)
        if label_path is None and atlas_library_root and species in species_atlas_fragments:
            try:
                label_path = get_atlas_label_path(
                    atlas_library_root, species_atlas_fragments[species], atlas_name)
                print(f"  Label file: {label_path}")
            except FileNotFoundError as e:
                warnings.warn(str(e))

        bold_paths_all[species] = []

        for bids_dir in bids_dirs:
            bids_dir = _linux_path(bids_dir)
            print(f"\n  ? BIDS: {bids_dir}")

            if "bold" in extract:
                bp = extract_bold_paths(bids_dir, lk, lr)
                bold_paths_all[species].extend(bp["bold_path"])

            if "surface" in extract:
                sp = extract_surface_paths(bids_dir, lk, lr)
                df = process_surfaces(sp, regions_of_interest, atlas_name=atlas_name)
                if not df.empty:
                    df.insert(0, "bids_dir", os.path.basename(bids_dir))
                    df.insert(0, "species",  species)
                    accumulators["surface"].append(df)

            if "thickness" in extract:
                tp = extract_thickness_paths(bids_dir, lk, lr)
                df = process_thickness(tp, regions_of_interest)
                if not df.empty:
                    df.insert(0, "bids_dir", os.path.basename(bids_dir))
                    df.insert(0, "species",  species)
                    accumulators["thickness"].append(df)

            if "volume" in extract and label_path:
                df = process_volumes_and_save(
                    bids_dir, atlas_name, label_path, lk, lr,
                    regions_of_interest, overwrite=False)
                if not df.empty:
                    df.insert(0, "bids_dir", os.path.basename(bids_dir))
                    df.insert(0, "species",  species)
                    accumulators["volume"].append(df)
            elif "volume" in extract and label_path is None:
                warnings.warn(f"No label file for {species} ? volume extraction skipped.")

            if "qc" in extract:
                qp = extract_qc_paths(bids_dir, lk, lr)
                df = process_qc(qp)
                if not df.empty:
                    df.insert(0, "bids_dir", os.path.basename(bids_dir))
                    df.insert(0, "species",  species)
                    accumulators["qc"].append(df)

    result = {
        k: pd.concat(v, ignore_index=True) if v else pd.DataFrame()
        for k, v in accumulators.items()
    }
    result["bold_paths"] = bold_paths_all
    return result


# ??????????????????????????????????????????????????????????????????????????????
# PART 4  ?  EXPORT
# ??????????????????????????????????????????????????????????????????????????????

def export_to_excel(data_dict, output_path, modalities=None):
    modalities = modalities or [
        k for k, v in data_dict.items()
        if isinstance(v, pd.DataFrame) and not v.empty
    ]
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for mod in modalities:
            df = data_dict.get(mod)
            if df is None or (isinstance(df, pd.DataFrame) and df.empty):
                print(f"  [export] skip empty: {mod}"); continue
            df.to_excel(writer, sheet_name=mod[:31], index=False)
            print(f"  [export] sheet '{mod}' ? {len(df)} rows")

    print(f"\n? Excel: {output_path}")
    return output_path


def export_to_csv(data_dict, output_dir, modalities=None):
    modalities = modalities or [
        k for k, v in data_dict.items()
        if isinstance(v, pd.DataFrame) and not v.empty
    ]
    os.makedirs(output_dir, exist_ok=True)
    out = {}
    for mod in modalities:
        df = data_dict.get(mod)
        if df is None or (isinstance(df, pd.DataFrame) and df.empty):
            continue
        p = opj(output_dir, f"{mod}.csv")
        df.to_csv(p, index=False)
        out[mod] = p
        print(f"  [export] {p}")
    return out


def _group_cols_for_summary(df: pd.DataFrame, base_groups: tuple) -> list:
    """
    Build the groupby column list for export_summary_stats.
    Only includes atlas_level and hemisphere when they actually exist
    in the DataFrame being summarised (volumes don't have atlas_level).
    """
    cols = [g for g in base_groups if g in df.columns]
    if "atlas_level" in df.columns:
        cols.append("atlas_level")
    if "hemisphere" in df.columns:
        cols.append("hemisphere")
    cols.append("region")
    return cols


def export_summary_stats(data_dict, output_path, groupby=("species",), modalities=None):
    """
    Compute descriptive statistics per group and export to Excel.

    atlas_level is included in grouping only when the modality DataFrame
    actually contains that column (surfaces & thickness yes, volumes no).
    """
    metric_cols = {
        "surface":   "surface_area_mm2",
        "volume":    "volume_mm3",
        "thickness": "thickness_mm",
    }
    modalities = modalities or list(metric_cols.keys())
    summaries  = {}
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for mod in modalities:
            df = data_dict.get(mod)
            if df is None or (isinstance(df, pd.DataFrame) and df.empty):
                print(f"  [summary] skip empty: {mod}"); continue
            col = metric_cols.get(mod)
            if col not in df.columns:
                print(f"  [summary] column '{col}' missing in {mod}"); continue

            grp_cols = _group_cols_for_summary(df, groupby)

            summary = (
                df.groupby(grp_cols)[col]
                .agg(n="count", mean="mean", std="std",
                     median="median", min="min", max="max")
                .reset_index()
            )
            sheet = f"{mod}_summary"[:31]
            summary.to_excel(writer, sheet_name=sheet, index=False)
            summaries[mod] = summary
            print(f"  [summary] sheet '{sheet}' ? {len(summary)} rows")

    print(f"? Summary stats: {output_path}")
    return summaries


# ??????????????????????????????????????????????????????????????????????????????
# PART 5  ?  PAPER-QUALITY PLOTS  (unchanged)
# ??????????????????????????????????????????????????????????????????????????????

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

PALETTE = ["#0072B2","#E69F00","#009E73","#CC79A7",
           "#56B4E9","#D55E00","#F0E442","#000000"]


def _species_colors(species_list):
    u = sorted(set(species_list))
    return {s: PALETTE[i % len(PALETTE)] for i, s in enumerate(u)}


def _fmt_y1(ax):
    import matplotlib.ticker as ticker
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))


def _sig_bracket(ax, x1, x2, y, p, h_frac=0.04):
    if   p < 0.001: label = "***"
    elif p < 0.01:  label = "**"
    elif p < 0.05:  label = "*"
    else:           return
    y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
    h = y_range * h_frac
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.1, c="k")
    ax.text((x1+x2)/2, y+h, label, ha="center", va="bottom", fontsize=10)


def _violin_strip_quartiles(ax, data, x_pos, color, width=0.4):
    data = np.asarray(data, dtype=float)
    data = data[~np.isnan(data)]
    if len(data) == 0:
        return

    jitter = np.random.default_rng(42).uniform(-width*0.18, width*0.18, len(data))
    ax.scatter(x_pos + jitter, data, color="white", s=18, zorder=3, linewidths=0)
    ax.scatter(x_pos + jitter, data, color=color,   s=12, zorder=4,
               alpha=0.75, linewidths=0.5, edgecolors="k")
    for q, lw in zip(np.percentile(data, [25, 50, 75]), [1, 2, 1]):
        ax.hlines(q, x_pos - width/2, x_pos + width/2, colors=color, linewidths=lw)

    if np.std(data) == 0 or len(data) < 2:
        ax.hlines(data[0], x_pos - width/2, x_pos + width/2,
                  colors=color, linewidths=2, alpha=0.6)
        return

    from scipy.stats import gaussian_kde
    try:
        kde = gaussian_kde(data, bw_method="scott")
    except Exception:
        return

    y_g  = np.linspace(data.min(), data.max(), 200)
    dens = kde(y_g)
    dens = dens / dens.max() * (width / 2)
    ax.fill_betweenx(y_g, x_pos - dens, x_pos + dens,
                     color=color, alpha=0.35, linewidth=0)
    ax.plot(x_pos - dens, y_g, color=color, lw=0.8, alpha=0.6)
    ax.plot(x_pos + dens, y_g, color=color, lw=0.8, alpha=0.6)
    for q, lw in zip(np.percentile(data, [25, 50, 75]), [1, 2, 1]):
        dens_q = float(kde(q)) / kde(y_g).max() * (width / 2)
        ax.hlines(q, x_pos - dens_q, x_pos + dens_q, colors=color, linewidths=lw)


def plot_morphometry_intra_bids(df, metric_col, metric_label, regions, output_path,
                                 species=None, hemisphere=None, bids_label=None,
                                 atlas_level=1, bids_col=None, figsize=None):
    with plt.rc_context(PAPER_RC):
        plot_df = df.copy()
        if hemisphere and hemisphere != 'bilateral':
            filtered = plot_df[plot_df["hemisphere"] == hemisphere]
            plot_df = filtered if not filtered.empty else plot_df[plot_df["hemisphere"] == 'bilateral']
        elif hemisphere == 'bilateral':
            plot_df = plot_df[plot_df["hemisphere"] == 'bilateral']
        if "atlas_level" in plot_df.columns:
            plot_df = plot_df[plot_df["atlas_level"] == atlas_level]
        plot_df = plot_df[plot_df["region"].isin(regions)]
        if plot_df.empty:
            warnings.warn("plot_morphometry_intra_bids: no data"); return None

        n    = len(regions)
        w, h = figsize or (max(5, n * 2.4), 5)
        fig, ax = plt.subplots(figsize=(w, h))

        color = PALETTE[0]
        np.random.seed(42)
        for i, region in enumerate(regions):
            vals = plot_df[plot_df["region"] == region][metric_col].dropna().values
            _violin_strip_quartiles(ax, vals, i, color)

        ax.set_xticks(range(n))
        ax.set_xticklabels(regions, rotation=35, ha="right")
        ax.set_ylabel(metric_label)
        _fmt_y1(ax)
        title = metric_label
        if species:    title += f" ? {species}"
        if bids_label: title += f" ({bids_label})"
        if hemisphere: title += f" [{hemisphere}]"
        ax.set_title(title, fontweight="bold", pad=8)
        ax.set_xlim(-0.6, n - 0.4)

        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight")
        plt.close(fig)
        print(f"  [plot] {output_path}")
    return output_path


def plot_morphometry_inter_species(df, metric_col, metric_label, regions, output_path,
                                    hemisphere=None, atlas_level=1,
                                    figsize=None, show_stats=True, show_n=True):
    with plt.rc_context(PAPER_RC):
        plot_df = df.copy()
        if hemisphere and hemisphere != 'bilateral':
            filtered = plot_df[plot_df["hemisphere"] == hemisphere]
            plot_df = filtered if not filtered.empty else plot_df[plot_df["hemisphere"] == 'bilateral']
        elif hemisphere == 'bilateral':
            plot_df = plot_df[plot_df["hemisphere"] == 'bilateral']
        if "atlas_level" in plot_df.columns:
            plot_df = plot_df[plot_df["atlas_level"] == atlas_level]
        plot_df = plot_df[plot_df["region"].isin(regions)]
        if plot_df.empty:
            warnings.warn("plot_morphometry_inter_species: no data"); return None

        _PHYLO = ['Bat','Rat','Mouse','Mouse_lemur','Marmoset','Macaque','Human','Dog','Pig','Cat']
        _known = [s for s in _PHYLO if s in plot_df["species"].unique()]
        _unk   = sorted([s for s in plot_df["species"].unique() if s not in _PHYLO])
        species_list = _known + _unk
        colors       = _species_colors(species_list)
        n_reg        = len(regions)
        w, h         = figsize or (max(6, n_reg * 3.6), 5.5)
        fig, axes    = plt.subplots(1, n_reg, figsize=(w, h), sharey=False)
        if n_reg == 1: axes = [axes]

        np.random.seed(42)
        for ax, region in zip(axes, regions):
            rdf = plot_df[plot_df["region"] == region]
            if rdf.empty: ax.set_visible(False); continue

            for xi, sp in enumerate(species_list):
                vals = rdf[rdf["species"] == sp][metric_col].dropna().values
                _violin_strip_quartiles(ax, vals, xi, colors[sp])
                if show_n and len(vals) > 0:
                    ax.text(xi, ax.get_ylim()[0] if ax.get_ylim()[0] != 0 else 0,
                            f"n={len(vals)}", ha="center", va="top",
                            fontsize=8, color="gray")

            if show_stats and len(species_list) >= 2:
                groups = [rdf[rdf["species"] == sp][metric_col].dropna().values
                          for sp in species_list]
                valid = [g for g in groups if len(g) >= 2]
                if len(valid) >= 2:
                    _, kw_p  = stats.kruskal(*valid)
                    y_max    = rdf[metric_col].dropna().max()
                    y_range  = y_max - rdf[metric_col].dropna().min()
                    step     = y_range * 0.14

                    for lvl, (i, j) in enumerate(combinations(range(len(species_list)), 2)):
                        g1 = rdf[rdf["species"] == species_list[i]][metric_col].dropna().values
                        g2 = rdf[rdf["species"] == species_list[j]][metric_col].dropna().values
                        if len(g1) < 2 or len(g2) < 2: continue
                        _, p_mw = stats.mannwhitneyu(g1, g2, alternative="two-sided")
                        _sig_bracket(ax, i, j, y_max + step * (lvl + 1), p_mw)

                    lbl = f"KW p={kw_p:.3f}" if kw_p >= 0.001 else "KW p<0.001"
                    ax.text(0.97, 0.01, lbl, transform=ax.transAxes,
                            ha="right", va="bottom", fontsize=8, color="gray")

            ax.set_xticks(range(len(species_list)))
            ax.set_xticklabels(species_list, rotation=30, ha="right")
            ax.set_title(region, fontweight="bold")
            ax.set_ylabel(metric_label if ax is axes[0] else "")
            ax.set_xlim(-0.6, len(species_list) - 0.4)
            _fmt_y1(ax)

        from matplotlib.patches import Patch
        handles  = [Patch(facecolor=colors[sp], label=sp) for sp in species_list]
        hemi_str = f" [{hemisphere}]" if hemisphere else ""
        fig.suptitle(f"{metric_label} ? inter-species{hemi_str}",
                     fontsize=13, fontweight="bold")
        plt.tight_layout()
        fig.legend(handles=handles, loc="upper right",
                   bbox_to_anchor=(1.0, 1.0), title="Species",
                   fontsize=9, borderaxespad=0.)
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight")
        plt.close(fig)
        print(f"  [plot] {output_path}")
    return output_path


def plot_multi_bids_comparison(df, metric_col, metric_label, regions, output_path,
                                species=None, hemisphere=None, atlas_level=1,
                                bids_col=None, figsize=None):
    with plt.rc_context(PAPER_RC):
        plot_df = df.copy()
        if species:   plot_df = plot_df[plot_df["species"] == species]
        if hemisphere and hemisphere != 'bilateral':
            filtered = plot_df[plot_df["hemisphere"] == hemisphere]
            plot_df = filtered if not filtered.empty else plot_df[plot_df["hemisphere"] == 'bilateral']
        elif hemisphere == 'bilateral':
            plot_df = plot_df[plot_df["hemisphere"] == 'bilateral']
        if "atlas_level" in plot_df.columns:
            plot_df = plot_df[plot_df["atlas_level"] == atlas_level]
        plot_df = plot_df[plot_df["region"].isin(regions)]
        if plot_df.empty:
            warnings.warn("plot_multi_bids_comparison: no data"); return None

        bids_list = sorted(plot_df["bids_dir"].unique())
        colors    = bids_col or {b: PALETTE[i % len(PALETTE)] for i, b in enumerate(bids_list)}
        n_reg     = len(regions)
        w, h      = figsize or (max(6, n_reg * 3.6), 5.5)
        fig, axes = plt.subplots(1, n_reg, figsize=(w, h), sharey=False)
        if n_reg == 1: axes = [axes]

        np.random.seed(42)
        for ax, region in zip(axes, regions):
            rdf = plot_df[plot_df["region"] == region]
            if rdf.empty: ax.set_visible(False); continue
            for xi, bids in enumerate(bids_list):
                vals = rdf[rdf["bids_dir"] == bids][metric_col].dropna().values
                _violin_strip_quartiles(ax, vals, xi, colors[bids])

            ax.set_xticks(range(len(bids_list)))
            ax.set_xticklabels(bids_list, rotation=30, ha="right")
            ax.set_title(region, fontweight="bold")
            ax.set_ylabel(metric_label if ax is axes[0] else "")
            _fmt_y1(ax)

        from matplotlib.patches import Patch
        handles = [Patch(facecolor=colors[b], label=b) for b in bids_list]
        title   = f"{metric_label} ? multi-BIDS"
        if species: title += f" ({species})"
        fig.suptitle(title, fontsize=13, fontweight="bold")
        plt.tight_layout()
        fig.legend(handles=handles, loc="upper right",
                   bbox_to_anchor=(1.0, 1.0), title="BIDS directory",
                   fontsize=9, borderaxespad=0.)
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight")
        plt.close(fig)
        print(f"  [plot] {output_path}")
    return output_path


def plot_qc_dashboard(qc_df, output_path, func_metrics=None, anat_metrics=None,
                      figsize=(14, 8)):
    with plt.rc_context(PAPER_RC):
        auto_func = ["func_avg_snr_gray","func_avg_snr_white","func_cnr",
                     "func_TSNR_0","func_mean_fd","func_mean_dvars",
                     "func_fwhm","func_gcor","func_censor_fraction"]
        auto_anat = ["anat_global_snr","anat_avg_snr_gray","anat_avg_snr_white",
                     "anat_template_correlation","anat_cortical_contrast",
                     "anat_fwhm","anat_noise_estimate"]
        func_metrics = func_metrics or [m for m in auto_func if m in qc_df.columns]
        anat_metrics = anat_metrics or [m for m in auto_anat if m in qc_df.columns]
        all_m = [m for m in func_metrics + anat_metrics if m in qc_df.columns]
        if not all_m:
            warnings.warn("No QC metrics found"); return None

        has_species = "species" in qc_df.columns
        colors  = _species_colors(qc_df["species"].unique()) if has_species else None
        ncols   = min(4, len(all_m))
        nrows   = math.ceil(len(all_m) / ncols)
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        axes    = np.array(axes).flatten()

        for i, metric in enumerate(all_m):
            ax   = axes[i]
            data = qc_df[["species", metric]].dropna() if has_species else qc_df[[metric]].dropna()
            if has_species:
                sns.boxplot(data=data, x="species", y=metric, palette=colors,
                            width=0.45, linewidth=1.1, ax=ax,
                            flierprops=dict(marker="o", markersize=3, alpha=0.5))
                sns.stripplot(data=data, x="species", y=metric, palette=colors,
                              size=4, alpha=0.65, jitter=True, ax=ax)
                ax.tick_params(axis="x", rotation=30)
            else:
                ax.boxplot(data[metric].values, widths=0.4)
            ax.set_title(metric.replace("func_","").replace("anat_",""), fontweight="bold")
            ax.set_xlabel(""); ax.set_ylabel("")
            _fmt_y1(ax)

        for j in range(i+1, len(axes)):
            axes[j].set_visible(False)

        fig.suptitle("QC Dashboard", fontsize=14, fontweight="bold")
        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight")
        plt.close(fig)
        print(f"  [plot] {output_path}")
    return output_path


# ??????????????????????????????????????????????????????????????????????????????
# CONVENIENCE WRAPPER
# ??????????????????????????????????????????????????????????????????????????????

def run_full_pipeline(
    species_config,
    output_dir,
    regions_of_interest=None,
    plot_regions=None,
    hemispheres=("left", "right"),
    extract=("surface", "volume", "thickness", "qc"),
    atlas_name="EDNIxCSC",
    atlas_label_paths=None,
    atlas_library_root=None,
    species_atlas_fragments=None,
):
    os.makedirs(output_dir, exist_ok=True)
    print("\n" + "="*70 + "\nSTEP 1 ? DATA COLLECTION")

    data = collect_multi_species(
        species_config, regions_of_interest=regions_of_interest,
        extract=extract, atlas_name=atlas_name,
        atlas_label_paths=atlas_label_paths,
        atlas_library_root=atlas_library_root,
        species_atlas_fragments=species_atlas_fragments,
    )

    print("\n" + "="*70 + "\nSTEP 2 ? EXPORT")
    excel_path = opj(output_dir, "ednix_results.xlsx")
    export_to_excel(data, excel_path)
    export_to_csv(data, opj(output_dir, "csv"))
    export_summary_stats(data, opj(output_dir, "ednix_summary_stats.xlsx"))

    print("\n" + "="*70 + "\nSTEP 3 ? PLOTS")
    plots_dir = opj(output_dir, "figures")
    os.makedirs(plots_dir, exist_ok=True)
    plots     = []

    metric_map = {
        "surface":   ("surface_area_mm2", "Surface area (mm²)"),
        "volume":    ("volume_mm3",        "Volume (mm³)"),
        "thickness": ("thickness_mm",      "Cortical thickness (mm)"),
    }
    n_species = len(species_config)

    for mod, (mcol, mlabel) in metric_map.items():
        df = data.get(mod)
        if df is None or df.empty:
            continue
        regions = plot_regions or df["region"].value_counts().index[:6].tolist()

        hemis_present = df["hemisphere"].unique() if "hemisphere" in df.columns else ["bilateral"]
        hemi_list = [h for h in hemispheres if h in hemis_present]
        if not hemi_list:
            hemi_list = list(hemis_present)

        for hemi in hemi_list:
            if n_species > 1:
                p = opj(plots_dir, f"{mod}_{hemi}_inter_species.png")
                plot_morphometry_inter_species(df, mcol, mlabel, regions, p, hemisphere=hemi)
                plots.append(p)

            for species in df["species"].unique():
                sdf    = df[df["species"] == species]
                n_bids = sdf["bids_dir"].nunique() if "bids_dir" in sdf.columns else 1

                if n_bids > 1:
                    p = opj(plots_dir, f"{mod}_{hemi}_{species}_multi_bids.png")
                    plot_multi_bids_comparison(sdf, mcol, mlabel, regions, p,
                                               species=species, hemisphere=hemi)
                    plots.append(p)

                p = opj(plots_dir, f"{mod}_{hemi}_{species}_intra.png")
                plot_morphometry_intra_bids(sdf, mcol, mlabel, regions, p,
                                             species=species, hemisphere=hemi)
                plots.append(p)

    if "qc" in extract and not data["qc"].empty:
        p = opj(plots_dir, "QC_dashboard.png")
        plot_qc_dashboard(data["qc"], p)
        plots.append(p)

    print(f"\n? Pipeline complete. Outputs: {output_dir}")
    return {"data": data, "excel_path": excel_path, "plots": plots}