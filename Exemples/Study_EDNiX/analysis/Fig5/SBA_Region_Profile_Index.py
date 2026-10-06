"""
SBA_Region_Profile_Index.py
─────────────────────────────
A new quality index for the retrosplenial-seed SBA (DMN), generalizing the
median_r_somato / std_r_somato approach to EVERY region of the EDNiX atlas
at once, instead of one hand-picked network.

Pipeline
--------
1. Per run: load the retrosplenial-seed raw correlation map (r, not
   Fisher-z) and, for EVERY region in the EDNiX atlas, compute median(r)
   and std(r) over that region's voxels. This is the run's "profile" —
   one point per region.
2. Reference/population curve: the SAME per-region median(r)/std(r)
   profile, computed once from the existing group-level SBA map
   (Overall_Effect_coef_norm / correlation_norm) — "the standard for the
   DMN" per-species.
3. Each run's deviation from that reference curve (spatial-correlation or
   Euclidean distance across the per-region profile, reusing the same
   distance/MAD-threshold machinery as the rest of the pipeline) becomes
   the new quality index — Valid / Not valid per run, per species.
4. Cross-species version: restrict to regions common to every species'
   atlas, pool across species into ONE inter-species reference curve, and
   reclassify every run (any species) by its distance to that shared
   curve instead of its own species' curve.

⚠ Atlas volume path — CONFIRMED, not guessed (2024): the parcellation is
stored PER SUBJECT/SESSION, right next to the SBA maps:
    sub-<subj>/ses-<ses>/func/templates/EDNiX/labels/<subj>_seg-<atlas>_dseg.nii.gz
See find_atlas_volume_path() below. This means, unlike an earlier version
of this module, there is no single "one path per species" to fill in —
each run uses its OWN subject's atlas volume (same-subject registration,
so no cross-subject grid-mismatch risk for that pairing specifically).
For the group-level reference map, pick any one representative subject's
atlas volume via pick_representative_atlas_volume() below — if its grid
doesn't match the group map's, compute_region_profile raises a clear
error rather than silently misaligning voxels (grid mismatches have been
common elsewhere in this pipeline, e.g. Rat 48x48x30 vs 48x82x48).
════════════════════════════════════════════════════════════════════════
"""

import os
import re
import sys
import glob
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib
from scipy.stats import norm as _scipy_norm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Local modules ────────────────────────────────────────────────────────────
try:
    from Exemples.Study_EDNiX.analysis import Figure_5_SBA as fig5
    from Exemples.Study_EDNiX.analysis.Fig5 import _Group_anal_Outlier_QC as qc
except ImportError:
    try:
        _this_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        _this_dir = os.getcwd()
    if _this_dir not in sys.path:
        sys.path.insert(0, _this_dir)
    import Figure_5_SBA as fig5
    import _Group_anal_Outlier_QC as qc

try:
    from Exemples.Study_EDNiX.analysis.ednix_bids_tools import parse_label_file
except ImportError:
    parse_label_file = None  # only needed if you don't pass label_df directly

opj = os.path.join
opb = os.path.basename
ope = os.path.exists

BOLD_SUFFIX = "_space-template_desc-fMRI_residual.nii.gz"
RETROSPLENIAL_SEEDS = ["L_BA_26_29_30__retrosplenial", "R_BA_26_29_30__retrosplenial"]


# =============================================================================
# 0) ITK-SNAP label file parser
# =============================================================================
# ednix_bids_tools.parse_label_file() expects a StatsLUT-style single-line
# "ID NAME R G B A" format, or an alternating name/color-line format —
# neither matches genuine ITK-SNAP label description files (single line:
# "ID R G B A VIS MSH \"NAME\""), and silently produces garbage when fed
# one: confirmed on EDNIxCSCLR_ITK.txt — (9522, 8) rows, 9228 duplicates,
# region_name containing whole raw lines repeated dozens of times instead
# of actual region names. Use parse_itk_snap_label_file() below instead
# for any *_ITK.txt label file.

_ITK_LABEL_LINE_RE = re.compile(
    r'^\s*(?P<id>\d+)\s+(?P<r>\d+)\s+(?P<g>\d+)\s+(?P<b>\d+)\s+'
    r'(?P<a>[\d.]+)\s+(?P<vis>\d+)\s+(?P<msh>\d+)\s+"(?P<name>[^"]*)"\s*$'
)


def parse_itk_snap_label_file(label_path):
    """
    Parse a genuine ITK-SNAP label description file:
        # comment lines start with '#' — skipped, however many there are
        IDX  -R- -G- -B- -A-- VIS MSH  LABEL
        0    0   0   0   0    0   0    "Clear Label"
        1    255 0   0   1    1   1    "L_CSF"
        ...
    (one label per line, unlike the StatsLUT/alternating formats
    ednix_bids_tools.parse_label_file() handles).

    Returns a DataFrame with the SAME columns as parse_label_file():
    label_id, region_name, base_name, hemisphere, R, G, B, A — drop-in
    compatible with the rest of this module. label_id 0 ("Clear
    Label"/background) is dropped, matching parse_label_file's own
    convention. Lines that don't match the expected format are skipped
    rather than raising, so stray header/malformed lines don't crash the
    whole parse — but duplicate label_ids (if any survive) are reported.
    """
    def _hemi(name):
        if re.match(r"^[Ll]_", name):
            return "left", re.sub(r"^[Ll]_", "", name)
        elif re.match(r"^[Rr]_", name):
            return "right", re.sub(r"^[Rr]_", "", name)
        return "bilateral", name

    rows = []
    n_skipped = 0
    with open(label_path, encoding="utf-8-sig") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = _ITK_LABEL_LINE_RE.match(line)
            if not m:
                n_skipped += 1
                continue
            lid = int(m.group("id"))
            if lid == 0:
                continue
            name = m.group("name")
            hemi, base = _hemi(name)
            rows.append(dict(
                label_id=lid, region_name=name, base_name=base, hemisphere=hemi,
                R=int(m.group("r")), G=int(m.group("g")), B=int(m.group("b")),
                A=float(m.group("a")),
            ))

    df = pd.DataFrame(rows, columns=["label_id", "region_name", "base_name",
                                     "hemisphere", "R", "G", "B", "A"])
    if not df.empty:
        n_dup = df["label_id"].duplicated().sum()
        if n_dup:
            print(f"  [WARN] {n_dup} duplicate label_id in {opb(label_path)} after "
                  f"parsing — keeping first occurrence")
            df = df.drop_duplicates(subset="label_id", keep="first").reset_index(drop=True)

    print(f"  [INFO] parsed {len(df)} labels from {opb(label_path)} (ITK-SNAP format)"
          + (f", skipped {n_skipped} unmatched line(s)" if n_skipped else ""))
    return df


# =============================================================================
# 1) Locate the raw (non-Fisher-z) per-run SBA correlation map
# =============================================================================

def _extract_run_id(filename):
    r"""
    Extract the run identifier from a filename — matches
    collect_bold_and_means's OWN confirmed real parsing convention, not
    just a bare "run-<digits>" — some species/BIDS use an
    acquisition-direction-tagged run naming instead: "run-RL-4",
    "run-LR-2", etc. A plain r"run-(\d+)" regex silently fails to match
    these (the letters "RL-"/"LR-" break it right after "run-"),
    falling back to a fixed "NA" for EVERY such run — confirmed real
    symptom: 4 distinct real runs (run-RL-1..4) for the same
    subject/session all collapsed onto the same dict key "NA" in
    profiles, silently overwriting each other (only the last-processed
    run survived) — undercounting real data and producing a "1 subject"
    3dttest++ group after averaging instead of the true multi-run count.

    Returns the run id as a string ("RL-4", "3", ...) — "NA" only if no
    run token is found at all in the filename.
    """
    m = re.search(r"run-(RL|LR)-(\d+)", filename)
    if m:
        return f"{m.group(1)}-{m.group(2)}"
    m = re.search(r"run-(\d+)", filename)
    if m:
        return m.group(1)
    return "NA"


def _bold_root(bold_path):
    """Mirrors _Group_anal_3dTtest.py's _bold_root exactly."""
    bn = opb(bold_path)
    for suffix in (BOLD_SUFFIX, ".nii.gz", ".nii"):
        if bn.endswith(suffix):
            return bn[: -len(suffix)]
    return bn


def sba_raw_corr_path(bold_path, seed_name):
    """
    Path to the RAW (non-Fisher-z) correlation map for one run/seed —
    confirmed real layout:
    .../Stats/SBA/<seed>/<bold_root>_correlations.nii.gz
    (the _fish version used elsewhere in the pipeline is
    <bold_root>_correlations_fish.nii.gz — this is the same folder,
    just without that suffix).
    """
    return opj(
        os.path.dirname(os.path.dirname(bold_path)),
        "Stats", "SBA", seed_name,
        f"{_bold_root(bold_path)}_correlations.nii.gz",
    )


def find_raw_corr_maps_for_run(subject_dir_root, subject, session, seed_name):
    """
    Glob every raw-correlation run map for one subject/session/seed,
    without needing the original bold_path list — matches the confirmed
    filename convention sub-<subj>_ses-<ses>_..._run-<run>_bold_correlations.nii.gz
    directly under Stats/SBA/<seed>/.

    Confirmed real layout ALSO includes a session-token-free convention
    for single-session datasets (e.g. Rat): filenames like
    "sub-<subj>_task-rest_bold_correlations.nii.gz" with NO "_ses-<ses>_"
    segment at all — the same "omit the token when there's only one"
    pattern already confirmed for runs (Dog's "no run-<N> token"
    single-run case). Try the strict ...ses-<ses>... pattern first; if
    that finds nothing, fall back to a session-token-free glob instead
    of assuming every dataset tags session in the filename.
    """
    ses_dir = opj(subject_dir_root, f"sub-{subject}", f"ses-{session}")
    if not ope(ses_dir):
        ses_dir = opj(subject_dir_root, f"sub-{subject}")
    seed_dir = opj(ses_dir, "func", "templates", "EDNiX", "Stats", "SBA", seed_name)
    if not ope(seed_dir):
        return []
    found = sorted(glob.glob(opj(seed_dir, f"sub-{subject}_ses-{session}_*_correlations.nii.gz")))
    if found:
        return found
    return sorted(glob.glob(opj(seed_dir, f"sub-{subject}_*_correlations.nii.gz")))


def find_atlas_volume_path(bids_dir, subject, session, atlas_name="EDNIxCSCLR"):
    """
    CONFIRMED layout: the parcellation volume in template space lives per
    subject/session, next to that subject's SBA maps:
        sub-<subj>/ses-<ses>/func/templates/EDNiX/labels/<subj>_seg-<atlas>_dseg.nii.gz
    (filename uses the bare subject id, not "sub-<subj>", matching e.g.
    "01_seg-EDNIxCSCLR_dseg.nii.gz").
    Returns None if not found rather than raising, so callers can skip/
    warn on a per-run basis instead of crashing the whole batch.
    """
    ses_dir = opj(bids_dir, f"sub-{subject}", f"ses-{session}")
    if not ope(ses_dir):
        ses_dir = opj(bids_dir, f"sub-{subject}")
    path = opj(ses_dir, "func", "templates", "EDNiX", "labels",
              f"{subject}_seg-{atlas_name}_dseg.nii.gz")
    return path if ope(path) else None


def pick_representative_atlas_volume(bids_dir, atlas_name="EDNIxCSCLR"):
    """
    For parcellating the GROUP-level reference map (which isn't tied to
    any one subject): grab the first available subject's atlas volume
    under this bids_dir. Group-level AFNI results and every subject's
    template-space data are supposed to share one common grid, so any
    subject's copy should do — but compute_region_profile() will still
    raise a clear error if this particular one doesn't actually match the
    group map's grid, rather than silently misaligning voxels.
    """
    pattern = opj(bids_dir, "sub-*", "ses-*", "func", "templates", "EDNiX",
                  "labels", f"*_seg-{atlas_name}_dseg.nii.gz")
    matches = sorted(glob.glob(pattern))
    if not matches:
        pattern_no_ses = opj(bids_dir, "sub-*", "func", "templates", "EDNiX",
                             "labels", f"*_seg-{atlas_name}_dseg.nii.gz")
        matches = sorted(glob.glob(pattern_no_ses))
    return matches[0] if matches else None


def find_group_map_path(results_base_dir, species, bids_name, seed_query, category="All",
                         output_suffix="Full", stat_suffix="correlation_norm"):
    """
    Resolve the real group map path for one seed, WITHOUT assuming
    whether the group results folder uses a lateralized name
    (L_BA_.../R_BA_...) or a merged bilateral name (BA_..., no prefix) —
    confirmed real layout varies by BIDS: e.g. Dog's
    Grp_SBA_3dTTEST_AllFull/ contains "BA_26_29_30__retrosplenial"
    (unprefixed — exclude_lateralized=True was used for that run), while
    other BIDS use separate L_.../R_... folders.

    seed_query can be given in EITHER form
    ("L_BA_26_29_30__retrosplenial" or "BA_26_29_30__retrosplenial") —
    matched by substring against whatever subfolder actually exists,
    rather than requiring an exact name match.

    Returns the full path to <seed_folder>/<seed_folder>_<stat_suffix>.nii.gz,
    or None if no matching seed folder was found — prints the list of
    what folders WERE found under the group dir, so a naming mismatch is
    immediately visible instead of a bare FileNotFoundError downstream.
    """
    # Two folder-naming conventions exist for the group dir, and which one
    # is on disk depends on WHICH function wrote it:
    #   "Grp_SBA_3dTTEST_<category><output_suffix>"  (older layout, e.g.
    #       "Grp_SBA_3dTTEST_AllFull")
    #   "Grp_SBA_3dTTEST_<output_suffix>"            (what
    #       rerun_balanced_by_run actually writes -- see its own
    #       results_subdir, which has no category prefix)
    # Confirmed real symptom of only trying the first: every HumanAnchored
    # run wrote ".../All/Grp_SBA_3dTTEST_<suffix>" but the lookup asked for
    # ".../All/Grp_SBA_3dTTEST_All<suffix>", so EVERY species returned None
    # and no distances could be computed at all. Try both, in that order.
    candidate_dirs = [
        opj(results_base_dir, species, bids_name, output_suffix, category,
            f"Grp_SBA_3dTTEST_{category}{output_suffix}"),
        opj(results_base_dir, species, bids_name, output_suffix, category,
            f"Grp_SBA_3dTTEST_{output_suffix}"),
    ]
    group_dir = next((d for d in candidate_dirs if ope(d)), None)
    if group_dir is None:
        print(f"  [WARN] group results dir not found, tried: {candidate_dirs}")
        return None

    available = sorted(d for d in os.listdir(group_dir) if os.path.isdir(opj(group_dir, d)))

    query_m = re.match(r"^([LlRr])_(.+)$", seed_query)
    query_hemi = query_m.group(1).upper() if query_m else None
    query_base = query_m.group(2) if query_m else seed_query

    matches = []
    for d in available:
        d_m = re.match(r"^([LlRr])_(.+)$", d)
        d_hemi = d_m.group(1).upper() if d_m else None
        d_base = d_m.group(2) if d_m else d

        if d_base != query_base:
            continue
        if d_hemi is None or query_hemi is None or d_hemi == query_hemi:
            # merged/bilateral folder (no hemisphere) matches any query;
            # an ambiguous no-hemisphere query matches any folder; a
            # lateralized folder only matches a query for the SAME
            # hemisphere — a lateralized query must NOT silently match
            # the opposite hemisphere's folder.
            matches.append(d)

    if not matches:
        print(f"  [WARN] no seed folder matching {seed_query!r} found under {group_dir} "
              f"— available folders: {available}")
        return None
    exact = [d for d in matches if d == seed_query]
    seed_folder = exact[0] if exact else matches[0]
    if len(matches) > 1 and not exact:
        print(f"  [WARN] multiple seed folders matched {seed_query!r}: {matches} — using {seed_folder}")
    # The stat file inside the seed folder is named
    # "<seed_folder><tag>_<stat_suffix>.nii.gz" where <tag> is whatever
    # output_suffix the writing function embedded -- e.g.
    # "BA_26_29_30__retrosplenial_BA_26_29_30__retrosplenial_HumanAnchored_kcore_ttest-stat_fisher_cmap.nii.gz".
    # Confirmed real symptom of building the name exactly as
    # f"{seed_folder}_{stat_suffix}.nii.gz": that tag was omitted, so the
    # lookup missed every HumanAnchored map even though the folder itself
    # resolved correctly. Glob instead, so any tag matches.
    exact_path = opj(group_dir, seed_folder, f"{seed_folder}_{stat_suffix}.nii.gz")
    if ope(exact_path):
        return exact_path

    pattern = opj(group_dir, seed_folder, f"{seed_folder}*_{stat_suffix}.nii.gz")
    hits = sorted(glob.glob(pattern))
    # Exclude derived variants that are NOT the base stat map itself
    # (e.g. the "_oversampled" siblings written for surface rendering) --
    # matching one of those would silently return a visualization-only
    # volume in place of the statistical result.
    hits = [h for h in hits if "_oversampled" not in opb(h)]
    if not hits:
        print(f"  [WARN] no stat map matching {pattern!r}")
        return None
    if len(hits) > 1:
        print(f"  [WARN] multiple stat maps matched {pattern!r}: "
             f"{[opb(h) for h in hits]} -- using {opb(hits[0])}")
    return hits[0]


# =============================================================================
# 2) Per-region profile: configurable index (median, std, %positive, ...)
# =============================================================================

# Registry of per-region summary statistics — add new ones here and they
# become usable everywhere else in the module (classification, plotting,
# inter-species pooling) just by passing their name as `value_col`,
# without touching any of that downstream code.
_REGION_STATS = {
    "median_r":      lambda vals: float(np.median(vals)),
    "mean_r":        lambda vals: float(np.mean(vals)),
    "std_r":         lambda vals: float(np.std(vals)),
    "iqr_r":         lambda vals: float(np.percentile(vals, 75) - np.percentile(vals, 25)),
    "frac_positive": lambda vals: float(np.mean(vals > 0)),
    "frac_negative": lambda vals: float(np.mean(vals < 0)),
}

DEFAULT_REGION_STATS = ("median_r", "std_r", "frac_positive", "frac_negative")


def _build_label_to_level_map(atlas_data_4d, label_ids):
    """
    For a 4D atlas volume (x, y, z, n_levels), determine — for EACH
    label_id individually — which level slice actually contains it as
    real painted voxels (rather than assuming one global level fits every
    region, which breaks for atlases where ID ranges span multiple
    hierarchy depths in one flat numbering, e.g. IDs 1-13 = coarse,
    21-35 = lobes, 51+ = fine parcellation, all in the same label file).

    Returns (id_to_level: {label_id: level_index_0based}, multi_level_ids:
    list of label_ids that were found painted in more than one level —
    for those, the level with the most voxels for that ID is kept, and
    the ambiguity is reported rather than silently picked).
    """
    label_id_set = {int(x) for x in label_ids}
    n_levels = atlas_data_4d.shape[-1]
    level_value_sets = [set(np.unique(atlas_data_4d[..., lvl]).astype(int)) - {0}
                        for lvl in range(n_levels)]

    id_to_level, multi_level_ids = {}, []
    for lid in label_id_set:
        found = [lvl for lvl in range(n_levels) if lid in level_value_sets[lvl]]
        if not found:
            continue
        if len(found) > 1:
            counts = [(lvl, int((atlas_data_4d[..., lvl] == lid).sum())) for lvl in found]
            counts.sort(key=lambda x: -x[1])
            id_to_level[lid] = counts[0][0]
            multi_level_ids.append(lid)
        else:
            id_to_level[lid] = found[0]
    return id_to_level, multi_level_ids


def compute_region_profile(map_path, atlas_vol_path, label_df, stats=DEFAULT_REGION_STATS,
                            atlas_level=None, expected_range=None, restrict_to_nonzero_map=True):
    """
    For one correlation map and the atlas parcellation volume (SAME
    spatial grid), compute every statistic in `stats` (see _REGION_STATS)
    per atlas region — median(r) and std(r) by default, but also
    %positive/%negative voxels, IQR, etc. All requested stats are
    computed in the same pass over the data, so switching WHICH one
    drives the classification/plot downstream (via `value_col` on the
    other functions in this module) never requires re-extracting voxels.

    label_df: output of parse_itk_snap_label_file() (or
    ednix_bids_tools.parse_label_file()) — columns label_id, region_name,
    base_name, hemisphere.

    atlas_level: if the atlas volume is 4D (multi-level dseg, one 3D
    label volume per hierarchy level — confirmed real layout, e.g. shape
    (37, 63, 27, 4)):
      - None (default): PER-REGION auto-detection — each region is
        looked up individually in whichever level slice actually
        contains its label_id as painted voxels (see
        _build_label_to_level_map). This is the correct behavior when
        different regions genuinely belong to different hierarchy
        depths within one flat label file, which is the confirmed case
        here — a single global level covers only a fraction of regions,
        and will silently DROP any region that doesn't live on that one
        forced level (e.g. forcing level 3 can make the retrosplenial
        seed region itself vanish from the profile if it lives on a
        different level).
      - an int: force EVERY region to that one level instead (1-indexed,
        matching the map numbers used elsewhere in this codebase) —
        only useful if you know all of label_df's regions live on one
        specific level.

    expected_range: if given (e.g. (-1, 1) for a genuine correlation
    map), runs sanity_check_profile() on the result and prints any
    region whose value is outside that range — do NOT set this for
    value_cols computed from an unbounded source map (Z-stats, raw
    beta/coef maps), since those are not bounded and would trip false
    warnings. None (default) skips the check.

    restrict_to_nonzero_map: True (default) — AFNI writes literal 0.0
    (not NaN) for every voxel OUTSIDE the analysis mask used by
    3dttest++/3dLMEr. The atlas region mask alone doesn't know about
    that mask, so a region whose true analysis coverage is small
    relative to its atlas extent gets its median/mean dragged toward 0
    by masked-out zero-fill voxels — confirmed on real data: dozens of
    regions reporting median_r EXACTLY 0.000000, which is not plausible
    for a continuous correlation and un-derestimates real signal.
    Excluding exact-zero voxels before computing stats fixes this,
    matching the "nonzero = valid" convention already used elsewhere in
    this pipeline for these AFNI-produced maps (e.g. the auto-mask
    fallback in Figure_5_SBA_OutlierQC.py). Set False only if you are
    profiling a map where a real value of exactly 0.0 is both possible
    and meaningful to keep.

    Returns DataFrame: region_name, base_name, hemisphere, label_id,
    n_voxels (raw atlas region size), n_voxels_valid (after excluding
    exact-zero/masked-out voxels, if restrict_to_nonzero_map), plus one
    column per requested stat.
    """
    unknown = [s for s in stats if s not in _REGION_STATS]
    if unknown:
        raise ValueError(f"Unknown stat(s) {unknown} — available: {list(_REGION_STATS)}")

    map_data = qc._squeeze_trailing_singleton(nib.load(map_path).get_fdata())

    atlas_img = nib.load(atlas_vol_path)
    atlas_data_raw = qc._squeeze_trailing_singleton(atlas_img.get_fdata())

    is_multilevel = atlas_data_raw.ndim == 4
    spatial_shape = atlas_data_raw.shape[:3] if is_multilevel else atlas_data_raw.shape

    if map_data.shape != spatial_shape:
        raise ValueError(
            f"Grid mismatch: map {map_data.shape} vs atlas volume spatial shape "
            f"{spatial_shape} ({atlas_vol_path}). The atlas volume must be in the "
            f"SAME template space/grid as the SBA maps."
        )

    label_ids = label_df["label_id"].astype(int).tolist()

    if is_multilevel:
        if atlas_level is not None:
            idx = atlas_level - 1
            if not (0 <= idx < atlas_data_raw.shape[-1]):
                raise ValueError(f"atlas_level={atlas_level} out of range — this atlas volume "
                                  f"has {atlas_data_raw.shape[-1]} levels")
            id_to_level = {lid: idx for lid in label_ids}
            print(f"  [INFO] {opb(atlas_vol_path)}: using atlas_level={atlas_level} "
                  f"(forced) for all {len(set(label_ids))} regions")
        else:
            id_to_level, multi = _build_label_to_level_map(atlas_data_raw, label_ids)
            n_labels = len(set(label_ids))
            print(f"  [INFO] {opb(atlas_vol_path)}: per-region level assignment — "
                  f"{len(id_to_level)}/{n_labels} label_ids found painted in the atlas "
                  f"volume (across its {atlas_data_raw.shape[-1]} levels)"
                  + (f"; {len(multi)} appeared in >1 level (kept the one with more voxels)"
                     if multi else ""))
    else:
        id_to_level = None  # atlas_data_raw is already the one 3D label volume

    rows = []
    for _, row in label_df.iterrows():
        lid = int(row["label_id"])
        if is_multilevel:
            lvl = id_to_level.get(lid)
            mask = (atlas_data_raw[..., lvl].astype(int) == lid) if lvl is not None \
                   else np.zeros(spatial_shape, dtype=bool)
        else:
            mask = atlas_data_raw.astype(int) == lid
        n = int(mask.sum())
        out = dict(region_name=row["region_name"], base_name=row["base_name"],
                  hemisphere=row["hemisphere"], label_id=lid, n_voxels=n)
        if n == 0:
            for stat_name in stats:
                out[stat_name] = np.nan
            out["n_voxels_valid"] = 0
        else:
            vals = map_data[mask]
            vals = vals[np.isfinite(vals)]
            if restrict_to_nonzero_map:
                vals = vals[vals != 0]
            out["n_voxels_valid"] = int(vals.size)
            for stat_name in stats:
                out[stat_name] = _REGION_STATS[stat_name](vals) if vals.size else np.nan
        rows.append(out)
    profile = pd.DataFrame(rows)
    if expected_range is not None:
        for stat_name in stats:
            sanity_check_profile(profile, value_col=stat_name, expected_range=expected_range,
                                 label=opb(map_path))
    return profile


def diagnose_unmatched_regions(atlas_vol_path, label_df):
    """
    Return the subset of label_df whose label_id was NOT found painted in
    ANY level of atlas_vol_path — run this to check whether the unmatched
    regions form a sensible pattern (e.g. all fine cortical "BA_*" areas
    that might only exist as surface dlabels, not in this volumetric
    segmentation) rather than a parsing/alignment problem. Prints and
    returns a DataFrame with the same columns as label_df.
    """
    atlas_data = qc._squeeze_trailing_singleton(nib.load(atlas_vol_path).get_fdata())
    label_ids = label_df["label_id"].astype(int).tolist()

    if atlas_data.ndim == 4:
        id_to_level, _ = _build_label_to_level_map(atlas_data, label_ids)
        matched_ids = set(id_to_level.keys())
    else:
        present = set(np.unique(atlas_data).astype(int)) - {0}
        matched_ids = present & set(label_ids)

    missing = label_df[~label_df["label_id"].astype(int).isin(matched_ids)].copy()
    print(f"  {len(missing)}/{len(label_df)} regions from label_df were NOT found "
          f"painted anywhere in {opb(atlas_vol_path)}")
    return missing


def report_level_coverage(atlas_vol_path, label_df):
    """
    For a 4D multi-level atlas volume, report — for EACH of its levels —
    how many of label_df's regions are actually painted there. Use this
    to pick a single atlas_level deliberately (e.g. "level 3 covers
    180/276 of my regions, including the retrosplenial seed and every
    DMN-related region I care about") instead of guessing.
    Returns a DataFrame: level (1-indexed), n_matched, n_total,
    matched_region_names (list).
    """
    atlas_data = qc._squeeze_trailing_singleton(nib.load(atlas_vol_path).get_fdata())
    label_ids = label_df["label_id"].astype(int).tolist()
    id_to_name = dict(zip(label_df["label_id"].astype(int), label_df["region_name"]))

    if atlas_data.ndim != 4:
        present = set(np.unique(atlas_data).astype(int)) - {0}
        matched = [lid for lid in label_ids if lid in present]
        print(f"  Not a multi-level volume (3D) — {len(matched)}/{len(set(label_ids))} "
              f"regions matched")
        return pd.DataFrame([{"level": 1, "n_matched": len(matched), "n_total": len(set(label_ids)),
                              "matched_region_names": [id_to_name[i] for i in matched]}])

    rows = []
    print(f"  Per-level coverage for {opb(atlas_vol_path)}:")
    for lvl in range(atlas_data.shape[-1]):
        present = set(np.unique(atlas_data[..., lvl]).astype(int)) - {0}
        matched = [lid for lid in label_ids if lid in present]
        print(f"    level {lvl + 1}: {len(matched)}/{len(set(label_ids))} regions matched")
        rows.append({"level": lvl + 1, "n_matched": len(matched), "n_total": len(set(label_ids)),
                     "matched_region_names": [id_to_name[i] for i in matched]})
    return pd.DataFrame(rows)


def sanity_check_profile(profile, value_col="median_r", expected_range=(-1, 1), label="profile"):
    """
    Print a warning + the offending region names/label_ids for any region
    whose value_col falls outside expected_range — e.g. a genuine Pearson
    r must be in [-1, 1], so anything beyond that is a red flag (grid
    mismatch, wrong atlas level for that region, or the source map isn't
    actually a correlation map at all — e.g. a Z-stat, which is NOT
    bounded and will trip this check spuriously).

    expected_range=None disables the check entirely — use this when
    value_col holds something unbounded by design (Z-stats, raw AFNI
    beta/coef maps, etc.), since there's no universal "impossible" range
    for those.

    Returns the offending rows (empty DataFrame if none / disabled).
    """
    if expected_range is None:
        return profile.iloc[0:0]
    lo, hi = expected_range
    bad = profile[(profile[value_col] < lo) | (profile[value_col] > hi)]
    if not bad.empty:
        print(f"  [WARN] {label}: {len(bad)} region(s) have {value_col} outside the "
              f"expected range [{lo}, {hi}] — impossible for a bounded correlation. "
              f"Check (a) whether the source map is really a correlation map, not a "
              f"Z-stat/beta/other unbounded quantity, and (b) whether these regions "
              f"are among any 'appeared in >1 level' ambiguous label_ids reported "
              f"during extraction:")
        for _, row in bad.iterrows():
            print(f"      label_id={row['label_id']:<6} {row['region_name']:<40} "
                  f"{value_col}={row[value_col]:.3f}")
    return bad


def copy_run_qc_images(df_class, bids_dir, output_root, clean_before_copy=True):
    """
    Copy each run's SBA preview jpg into <output_root>/kept/ or
    <output_root>/outliers/, based on df_class's is_outlier column — same
    confirmed layout as SBA_Visual_Labeler.py:
        sub-<subj>/ses-<ses>/func/templates/EDNiX/Stats/SBA/<seed>/
        sub-<subj>_ses-<ses>_..._run-<run>_bold_*.jpg

    df_class needs columns: subject, session, run, seed, is_outlier —
    exactly what classify_runs_by_profile_distance returns.

    Destination filenames are prefixed with subject/session/run/seed to
    avoid collisions (two different seeds' jpgs for the same run share a
    basename otherwise — same issue fixed earlier in
    Figure_5_SBA_OutlierQC.copy_subject_qc_images).

    clean_before_copy: True (default) — remove any existing kept/
    outliers subfolders under output_root before copying. Without this,
    os.makedirs(exist_ok=True) never removes prior contents, so a rerun
    (e.g. after changing min_correlation, or after a bugfix) silently
    mixes OLD files left over from a previous run together with the NEW
    ones — confirmed real source of confusion earlier in this project.
    Set False only if you deliberately want to ADD to what's already
    there without touching existing files.
    """
    import shutil as _shutil

    kept_dir = opj(output_root, "kept")
    out_dir = opj(output_root, "outliers")
    if clean_before_copy:
        for d in (kept_dir, out_dir):
            if ope(d):
                _shutil.rmtree(d)
    os.makedirs(kept_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    n_copied, n_missing = 0, 0
    for _, row in df_class.iterrows():
        subj, ses, run, seed = str(row["subject"]), str(row["session"]), str(row["run"]), str(row["seed"])
        dest_dir = out_dir if bool(row["is_outlier"]) else kept_dir

        seed_dir = opj(bids_dir, f"sub-{subj}", f"ses-{ses}", "func", "templates",
                       "EDNiX", "Stats", "SBA", seed)
        # Strict run-number match: the run digits must be followed by a
        # non-digit character, and allow an optional leading zero
        # (run-1 vs run-01) — a plain glob wildcard *run-{run}* wrongly
        # matches run-10, run-11, ... run-19 for run="1" (confirmed on
        # real data: 6 files copied under one run-1 destination prefix).
        # run="NA" means no run-<N> token was found in the source
        # filenames at all — confirmed real case: some datasets have a
        # single run per session with no run token whatsoever
        # (sub-05_ses-1_task-rest_bold_.jpg). There, fall back to
        # matching on subject/session alone — unambiguous since there's
        # only one run to find.
        #
        # Session token can ALSO be omitted entirely in the filename —
        # confirmed real Rat case: "sub-300101_task-rest_bold_.jpg", no
        # "_ses-1_" segment at all (same "omit the token when there's
        # only one" convention already fixed in find_raw_corr_maps_for_run
        # for the correlation maps themselves — this glob needed the
        # identical fallback and didn't have it, confirmed real symptom:
        # "[WARN] no jpg found" for every single Rat subject/run).
        candidates = glob.glob(opj(seed_dir, f"sub-{subj}_ses-{ses}_*.jpg"))
        if not candidates:
            candidates = glob.glob(opj(seed_dir, f"sub-{subj}_*.jpg"))
        if run == "NA":
            found = sorted(candidates)
        else:
            run_re = re.compile(rf"run-0*{re.escape(run)}(?!\d)")
            found = sorted(f for f in candidates if run_re.search(opb(f)))
        if not found:
            print(f"  [WARN] no jpg found for sub-{subj} ses-{ses} run-{run} seed={seed} "
                  f"under {seed_dir}")
            n_missing += 1
            continue
        for f in found:
            dest = opj(dest_dir, f"sub-{subj}_ses-{ses}_run-{run}_{seed}_{opb(f)}")
            _shutil.copy2(f, dest)
            n_copied += 1

    print(f"  {n_copied} image(s) copied, {n_missing} run(s) had no matching jpg "
          f"-> {kept_dir} / {out_dir}")
    return kept_dir, out_dir


def rerun_balanced_kept_vs_outlier(df_class, species, bids_dir, classification_path,
                                    category="All", results_base_dir=None,
                                    output_suffix_kept="Full_KeptBalanced",
                                    output_suffix_outlier="Full_OutlierBalanced",
                                    random_seed=42, skip_analyses=None):
    """
    Using df_class's per-run is_outlier classification, re-run the group
    3dttest++ (via Figure_5_SBA.process_species_category) TWICE — once
    using only "kept" subject-sessions, once using only "outlier"
    subject-sessions — with the SAME number of subject-sessions in both
    groups (the smaller of the two; the larger group is randomly
    subsampled down to match, reproducibly via random_seed) so any
    difference between the two resulting group maps isn't just a
    sample-size artifact.

    Aggregation: df_class is per-run, but process_species_category needs
    subject-SESSION granularity. A (subject, session) counts as
    "outlier" if ANY of its runs is flagged — conservative, matching
    Figure_5_SBA_OutlierQC.rerun_SBA_excluding_outliers's exclusion
    granularity elsewhere in this pipeline.

    Writes two separate group results (output_suffix_kept /
    output_suffix_outlier) — nothing already on disk is touched.

    Returns (success_kept, success_outlier, kept_pairs_used, outlier_pairs_used).
    """
    results_base_dir = results_base_dir or fig5.OUTPUT_BASE

    df_class = df_class.copy()
    df_class["subject"] = df_class["subject"].astype(str)
    df_class["session"] = df_class["session"].astype(str)
    session_outlier = (
        df_class.groupby(["subject", "session"])["is_outlier"]
                .any()
                .reset_index()
    )
    kept_pairs = sorted(map(tuple, session_outlier.loc[
        ~session_outlier["is_outlier"], ["subject", "session"]].values))
    outlier_pairs = sorted(map(tuple, session_outlier.loc[
        session_outlier["is_outlier"], ["subject", "session"]].values))

    n = min(len(kept_pairs), len(outlier_pairs))
    if n == 0:
        raise ValueError(f"One group is empty (kept={len(kept_pairs)}, "
                          f"outlier={len(outlier_pairs)}) — can't balance.")

    rng = np.random.default_rng(random_seed)

    def _subsample(pairs, n):
        if len(pairs) <= n:
            return set(pairs)
        idx = rng.choice(len(pairs), size=n, replace=False)
        return {pairs[i] for i in idx}

    kept_balanced = _subsample(kept_pairs, n)
    outlier_balanced = _subsample(outlier_pairs, n)

    print(f"  Balancing kept/outlier groups to n={n} subject-session(s) each "
          f"(kept had {len(kept_pairs)}, outlier had {len(outlier_pairs)})")
    print(f"  kept (balanced):    {sorted(kept_balanced)}")
    print(f"  outlier (balanced): {sorted(outlier_balanced)}")

    # bids_dir must be the SHORT name (e.g. "BIDS_Zhu_Garin"), not the full
    # path — confirmed real source: load_classification_subjects filters
    # via df_subject['bids_dir'] == bids_dir against the Excel's bids_dir
    # column, which stores only Path(d).name, not the full BIDS path.
    # And on no-match/empty, the real function returns the TUPLE
    # (None, None), not a dict — .get() on that raises AttributeError,
    # so check the type before calling it rather than just truthiness.
    category_subjects = fig5.load_classification_subjects(
        classification_path, species=species, bids_dir=Path(bids_dir).name,
        categories=[category],
    )
    if not isinstance(category_subjects, dict) or not category_subjects.get(category):
        raise RuntimeError(f"Could not reload classification for {species}/{category} "
                           f"(bids_dir={Path(bids_dir).name})")

    def _build_list_to_keep(target_pairs):
        return [
            (sp, bd, subj, ses) for sp, bd, subj, ses in category_subjects[category]
            if (str(subj), str(ses)) in target_pairs
        ]

    list_to_keep_kept = _build_list_to_keep(kept_balanced)
    list_to_keep_outlier = _build_list_to_keep(outlier_balanced)

    if not list_to_keep_kept or not list_to_keep_outlier:
        raise RuntimeError(
            "Balanced group(s) didn't match anything in the reloaded classification "
            f"(kept matched {len(list_to_keep_kept)}, outlier matched "
            f"{len(list_to_keep_outlier)}) — the (subject, session) values from "
            "df_class may not line up with what load_classification_subjects returns."
        )

    print(f"  Re-running SBA -> output_suffix='{output_suffix_kept}' "
          f"({len(list_to_keep_kept)} subject-session(s))")
    success_kept = fig5.process_species_category(
        species, [bids_dir], category, list_to_keep_kept,
        bids_label=Path(bids_dir).name, output_suffix=output_suffix_kept,
        skip_analyses=skip_analyses, results_base_dir=results_base_dir, min_runs=1,
    )

    print(f"  Re-running SBA -> output_suffix='{output_suffix_outlier}' "
          f"({len(list_to_keep_outlier)} subject-session(s))")
    success_outlier = fig5.process_species_category(
        species, [bids_dir], category, list_to_keep_outlier,
        bids_label=Path(bids_dir).name, output_suffix=output_suffix_outlier,
        skip_analyses=skip_analyses, results_base_dir=results_base_dir, min_runs=1,
    )

    return success_kept, success_outlier, kept_balanced, outlier_balanced


def _load_label_df_for_seed_matching(label_path, regions_of_interest):
    """
    _3dttest_EDNiX needs a label_df whose region_name matches
    regions_of_interest — but NOT as a raw substring check.
    _3dttest_EDNiX normalizes region_name through its own _format_seed()
    (replaces " ()/,:;.-" with "_", strips to alnum+underscore) BEFORE
    comparing against regions_of_interest — confirmed real source:
    "L_BA_26/29/30_(retrosplenial)" only becomes
    "L_BA_26_29_30__retrosplenial" (matching REGIONS_OF_INTEREST) AFTER
    that normalization; a raw substring check against the unformatted
    name always fails, regardless of which parser produced it. This
    was the actual bug — not the parser choice.

    WHICH PARSER (fig5.parse_label_file vs parse_itk_snap_label_file)
    still varies by file, so both are tried, but the match test now
    mirrors _format_seed exactly.
    """
    try:
        from Statistics.Group_fMRI._Group_anal_3dTtest import _format_seed
    except ImportError as e:
        raise ImportError(
            f"Could not import _format_seed from Statistics.Group_fMRI._Group_anal_3dTtest "
            f"— needed to normalize region names the same way _3dttest_EDNiX does. "
            f"Original error: {e}"
        )

    def _matches(df):
        if df is None or not isinstance(df, pd.DataFrame) or df.empty \
           or "region_name" not in df.columns:
            return False
        return bool(df["region_name"].apply(
            lambda b: any(r in _format_seed(b) for r in regions_of_interest)).any())

    def _diagnose(df, parser_name):
        if df is None or not isinstance(df, pd.DataFrame):
            print(f"  [DIAG] {parser_name}: returned {type(df).__name__}, not a DataFrame")
            return
        n = len(df)
        sample = df["region_name"].head(8).tolist() if "region_name" in df.columns else "<no region_name column>"
        print(f"  [DIAG] {parser_name}: {n} row(s), sample region_name values: {sample}")

    try:
        df1 = fig5.parse_label_file(label_path)
    except Exception as e:
        print(f"  [WARN] fig5.parse_label_file failed on {opb(label_path)}: {e}")
        df1 = None
    _diagnose(df1, "fig5.parse_label_file")
    if _matches(df1):
        print(f"  [INFO] {opb(label_path)}: using fig5.parse_label_file (StatsLUT format)")
        return df1

    df2 = parse_itk_snap_label_file(label_path)
    _diagnose(df2, "parse_itk_snap_label_file")
    if _matches(df2):
        print(f"  [INFO] {opb(label_path)}: using parse_itk_snap_label_file (ITK-SNAP format)")
        return df2

    raise RuntimeError(
        f"Neither fig5.parse_label_file nor parse_itk_snap_label_file produced a "
        f"label_df matching REGIONS_OF_INTEREST {regions_of_interest} from "
        f"{label_path} (after _format_seed normalization) — see the [DIAG] lines "
        f"above for what each parser actually returned."
    )


def _3dttest_EDNiX_run_level(
    bids_dir, templatehigh, templatelow, oversample_map, mask_func, cut_coords,
    label_df, lower_cutoff, upper_cutoff, MAIN_PATH, alpha, bold_paths, mean_imgs,
    regions_of_interest=None, cluster_size_default=10, opening=1, csim_pthr=0.05,
    vmax_percentile=99, method_mask_func='mask_func_over_Gray', output_suffix=None,
    normalize_target=95, exclude_lateralized=False,
    oversample_for_surface=True, oversample_target=None, unit="run",
):
    """
    Run-level variant of Statistics.Group_fMRI._Group_anal_3dTtest._3dttest_EDNiX.
    Identical in every other respect (same mask building, seed filtering,
    thresholding, plotting, output filenames/paths) — the ONLY difference
    is step 1.

    Confirmed from the real _3dttest_EDNiX source: it ALWAYS averages
    every run belonging to the SAME subject into ONE Fisher-z map (via
    3dMean) BEFORE 3dttest++ ever sees the data — one
    "<subject>_avg_fisher_map.nii.gz" per unique SUBJECT, regardless of
    how many distinct runs bold_paths actually contains — so 3dttest++'s
    degrees of freedom are always the number of SUBJECTS, never the
    number of runs. There is no parameter in the original function to
    disable this.

    This variant skips that per-subject averaging entirely: EACH run's
    own Fisher-z map (resampled to a common master grid) is passed
    directly to 3dttest++ -setA — one input per RUN. Everything
    downstream (thresholding, normalization, plotting, output filenames)
    is unchanged, copied verbatim from the original function.

    ⚠ Statistical caveat this deliberately does not hide: multiple runs
    from the SAME subject are not independent samples (within-subject
    correlation) — 3dttest++'s standard errors here will be
    optimistic/anti-conservative relative to a proper mixed-effects
    model (e.g. 3dLMEr with subject as a random effect) whenever a
    balanced group contains more than one run from the same subject.
    Use this specifically when you want run-level degrees of freedom
    and accept that tradeoff — it is not a silent drop-in replacement
    for the subject-averaged original everywhere else in this pipeline.

    oversample_for_surface: True (default) — after computing the raw
    zmap/correlation stat outputs, also write OVERSAMPLED COPIES of
    just those two (resampled onto oversample_target's grid via
    3dresample) as "..._oversampled.nii.gz" siblings, for visualization
    purposes only. The original zmap/correlation files are never
    touched — thresholding, normalization, and every statistical step
    below this one keep using the ORIGINAL (non-oversampled) data.
    Surface mesh vertices are typically much denser than this
    functional-resolution volume, so volume-to-surface projection of
    the raw volume can look blocky/patchy; oversampling onto a finer
    grid first gives a visually smoother surface result. Set False to
    skip this and only produce the original-resolution files.

    oversample_target: the grid to resample onto. None (default) uses
    `templatehigh` — the same per-species structural template
    (typically the EDNiX library's SS_T1w) this function already loads
    for the "oversample_map" background-image choice — since it's
    already available as an argument and is normally finer than the
    functional grid. Pass a different path only if you want a specific
    grid other than templatehigh.
    """
    from Tools import Load_EDNiX_requirement
    from Statistics.Group_fMRI._mask_utils import build_group_mask
    from Statistics.Group_fMRI._Group_anal_3dTtest import (
        _format_seed, _sba_fish_path, _safe_plot, _normalize_map,
    )

    sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _, sing_synstrip, Unetpath = \
        Load_EDNiX_requirement.load_requirement(MAIN_PATH, templatehigh, bids_dir, 'yes')

    output_results1 = opj(bids_dir, "Results")
    os.makedirs(output_results1, exist_ok=True)
    studytemplatebrain = templatehigh if oversample_map else templatelow

    mask_overlap = build_group_mask(
        mean_imgs=mean_imgs, mask_func=mask_func, output_results1=output_results1,
        sing_afni=sing_afni, studytemplatebrain=studytemplatebrain,
        lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff, opening=opening,
        method_mask_func=method_mask_func,
    )

    all_bold_list = bold_paths["bold_path"]

    folder_name = "Grp_SBA_3dTTEST" + (f"_{output_suffix}" if output_suffix else "")
    output_results = opj(output_results1, folder_name)
    os.makedirs(output_results, exist_ok=True)

    seeds = label_df.copy()
    if regions_of_interest:
        seeds = seeds[seeds["region_name"].apply(
            lambda b: any(r in _format_seed(b) for r in regions_of_interest))]
    if exclude_lateralized:
        seeds = seeds[~seeds["region_name"].apply(
            lambda b: _format_seed(b).startswith(("L_", "R_")))]

    tag = f"_{output_suffix}" if output_suffix else ""
    spgo = subprocess.getoutput
    import shutil as _shutil

    for _, seed_row in seeds.iterrows():
        seed_name = _format_seed(seed_row["region_name"])
        print(f"  Processing seed: {seed_name}")
        output_folder = opj(output_results, seed_name)
        # Wipe the whole working folder before writing — NOT just
        # "TTnew*" (the original narrower cleanup this was copied from).
        # Confirmed real symptom: leftover "<subject>_avg_fisher_map.nii.gz"
        # files from a PREVIOUS run using the subject-averaged
        # _3dttest_EDNiX (before this run-level variant existed) stayed in
        # this folder indefinitely and got copied into the final
        # destination alongside the NEW, correct run-level output,
        # making it look like per-subject averaging was still happening
        # even though 3dttest++ itself was already correctly seeing one
        # input per run (confirmed via "have N volumes" in the AFNI log
        # matching the run count, not the subject count).
        if ope(output_folder):
            _shutil.rmtree(output_folder)
        os.makedirs(output_folder, exist_ok=True)

        # ---- 1) Resolve each RUN's own Fisher-z map -- NOT averaged per
        # subject (this is the whole point of this variant) ----
        all_subj_list = bold_paths.get("subject", [None] * len(all_bold_list))
        all_ses_list = bold_paths.get("session", [None] * len(all_bold_list))
        all_run_list = bold_paths.get("run", [None] * len(all_bold_list))
        run_fish_maps = []  # list of (fish_path, subject, session, run) tuples
        for bold_p, subj_id, ses_id, run_id in zip(all_bold_list, all_subj_list,
                                                   all_ses_list, all_run_list):
            fp = _sba_fish_path(bold_p, seed_name)
            if ope(fp):
                run_fish_maps.append((fp, subj_id, ses_id, run_id))
            else:
                print(f"  [WARN] No Fisher-z map for run: {opb(bold_p)}  seed={seed_name}")

        # ---- 1b) Collapse runs within each UNIT, with AFNI ----
        # unit="run" leaves this untouched (one 3dttest++ input per run).
        # unit="subject"/"session" averages that unit's runs into ONE map
        # with 3dMean first, so 3dttest++ sees one input per ANIMAL. Two
        # reasons this is not cosmetic: (a) runs of the same subject are
        # not independent, so run-level degrees of freedom are
        # anti-conservative; (b) a subject scanned 20 times otherwise
        # contributes 20 times the weight of a subject scanned once, and
        # a dataset with more runs dominates any cross-species comparison
        # for a reason that has nothing to do with biology.
        if unit in ("subject", "session"):
            unit_dir = opj(output_folder, "unit_means")
            os.makedirs(unit_dir, exist_ok=True)
            groups = {}
            for fp, subj_id, ses_id, run_id in run_fish_maps:
                gkey = (str(subj_id),) if unit == "subject" else (str(subj_id), str(ses_id))
                groups.setdefault(gkey, []).append(fp)

            averaged = []
            for gkey, fps in sorted(groups.items()):
                gname = "sub-" + gkey[0] + (f"_ses-{gkey[1]}" if len(gkey) > 1 else "")
                out_path = opj(unit_dir, f"{gname}_avg_fisher_map.nii.gz")
                if len(fps) == 1:
                    _shutil.copy2(fps[0], out_path)
                else:
                    # Runs of one subject can sit on different grids (a
                    # second session is not guaranteed to share the
                    # first's), and 3dMean requires an exact match, so
                    # resample onto this unit's first map before averaging.
                    to_mean = [fps[0]]
                    for i, f in enumerate(fps[1:], start=1):
                        f_rs = opj(unit_dir, f"{gname}_{i:02d}_rs.nii.gz")
                        print(spgo(
                            f"{sing_afni} 3dresample -master {fps[0]} -input {f} "
                            f"-prefix {f_rs} -overwrite -bound_type SLAB"
                        ))
                        to_mean.append(f_rs if ope(f_rs) else f)
                    print(spgo(
                        f"{sing_afni} 3dMean -overwrite -prefix {out_path} "
                        f"{subprocess.list2cmdline(to_mean)}"
                    ))
                if ope(out_path):
                    averaged.append((out_path, gkey[0],
                                    gkey[1] if len(gkey) > 1 else "ALL",
                                    f"MEAN{len(fps)}"))
                else:
                    print(f"  [WARN] 3dMean produced nothing for {gname} -- unit dropped")
            print(f"  {len(run_fish_maps)} run(s) -> {len(averaged)} {unit}-level "
                  f"input(s) for 3dttest++ (runs per {unit}: "
                  f"{sorted(len(v) for v in groups.values())})")
            run_fish_maps = averaged

        if len(run_fish_maps) < 2:
            print(f"  [SKIP] {seed_name}: only {len(run_fish_maps)} {unit}(s) "
                  f"-- 3dttest++ needs at least 2 inputs")
            continue

        # ---- 2) Resample every run's map to a common master grid ----
        master_run, master_subj, master_ses, master_run_id = run_fish_maps[0]
        resampled_runs = [master_run]
        for fm, subj_id, ses_id, run_id in run_fish_maps[1:]:
            fm_rs = fm.replace(".nii.gz", "_rs.nii.gz")
            print(spgo(
                f"{sing_afni} 3dresample -master {master_run} -input {fm} "
                f"-prefix {fm_rs} -overwrite -bound_type SLAB"
            ))
            resampled_runs.append(fm_rs if ope(fm_rs) else fm)

        # Copy every input that actually went into -setA into
        # output_folder/input_runs/, with a clear sub-<subject>_run-<run>
        # name — confirmed real gap: these resampled per-run maps were
        # previously only written as siblings next to their SOURCE files
        # (scattered across each subject's own Stats/SBA/<seed>/ folder),
        # never centralized anywhere the final copied results end up, so
        # there was no easy way to directly confirm (by just looking at
        # the destination folder) which — and how many — individual runs
        # actually fed into 3dttest++'s -setA, as opposed to trusting the
        # AFNI log's "have N volumes" line alone.
        input_runs_dir = opj(output_folder, "input_runs")
        os.makedirs(input_runs_dir, exist_ok=True)
        for i, (rs_path, (_, subj_id, ses_id, run_id)) in enumerate(
                zip(resampled_runs, run_fish_maps), start=1):
            dest_name = f"{i:02d}_sub-{subj_id}_ses-{ses_id}_run-{run_id}.nii.gz"
            _shutil.copy2(rs_path, opj(input_runs_dir, dest_name))
        print(f"  [saved] {len(resampled_runs)} input map(s) at {unit} level "
              f"-> {input_runs_dir}")

        # ---- 3) 3dttest++ -- one input PER RUN, not per subject ----
        # (no separate "TTnew*" cleanup needed here anymore -- the whole
        # output_folder was already wiped fresh above)
        os.chdir(output_folder)

        # Resample the group mask onto the SAME master grid the run maps
        # were just resampled to. Confirmed real failure without this:
        # Rat died with "-mask doesn't match datasets number of voxels"
        # (mask had 15371 voxels) -- build_group_mask produces the mask on
        # its own grid, which is NOT guaranteed to match master_run's, and
        # 3dttest++ requires an exact voxel-count match. -rmode NN keeps
        # the mask binary (linear interpolation would create fractional
        # edge values).
        mask_for_ttest = opj(output_folder, "group_mask_rs.nii.gz")
        print(spgo(
            f"{sing_afni} 3dresample -master {master_run} -input {mask_overlap} "
            f"-prefix {mask_for_ttest} -overwrite -rmode NN"
        ))
        if not ope(mask_for_ttest):
            print(f"  [WARN] mask resample produced nothing -- falling back to the "
                 f"original mask ({mask_overlap}); 3dttest++ may still reject it "
                 f"if the grids differ")
            mask_for_ttest = mask_overlap

        use_clustsim = len(resampled_runs) > 14
        print(spgo(
            f"{sing_afni} "
            f"3dttest++ -setA {subprocess.list2cmdline(resampled_runs)} "
            f"-toz {'-Clustsim ' if use_clustsim else ''}"
            f"-mask {mask_for_ttest}"
        ))

        # ---- 4) Save AFNI outputs as NIfTI (identical to the original) ----
        zmap = opj(output_folder, f"{seed_name}{tag}_ttest-stat_fisher_zmap.nii.gz")
        correlation = opj(output_folder, f"{seed_name}{tag}_ttest-stat_fisher_cmap.nii.gz")

        print(spgo(
            f"{sing_afni} 3dcalc -overwrite "
            f"-a {opj(output_folder, 'TTnew+orig.HEAD[1]')} "
            f"-expr a -prefix {zmap}"
        ))
        print(spgo(
            f"{sing_afni} 3dcalc -overwrite "
            f"-a {opj(output_folder, 'TTnew+orig.HEAD[0]')} "
            f"-expr a -prefix {correlation}"
        ))

        # ---- 4b) Oversampled copies of zmap/cmap, for surface rendering
        # ONLY — never used for thresholding, plotting, or any statistical
        # step below. Surface mesh vertices are typically much denser than
        # this functional-resolution volume, so projecting the volume as-is
        # onto the surface gives a blocky/patchy look; resampling onto the
        # finer structural template grid (templatehigh — the same one this
        # function already loads as `base_ss`/`templatehigh`, e.g. the
        # EDNiX library's per-species SS_T1w) first gives a visually
        # smoother surface result. These are SEPARATE files
        # ("..._oversampled.nii.gz") — the original zmap/correlation files
        # above are never modified or overwritten by this step.
        if oversample_for_surface:
            oversample_grid = oversample_target or templatehigh
            zmap_oversampled = zmap.replace(".nii.gz", "_oversampled.nii.gz")
            corr_oversampled = correlation.replace(".nii.gz", "_oversampled.nii.gz")
            print(spgo(
                f"{sing_afni} 3dresample -master {oversample_grid} -input {zmap} "
                f"-prefix {zmap_oversampled} -overwrite -bound_type SLAB"
            ))
            print(spgo(
                f"{sing_afni} 3dresample -master {oversample_grid} -input {correlation} "
                f"-prefix {corr_oversampled} -overwrite -bound_type SLAB"
            ))

        # ---- 5) Cluster threshold (identical to the original) ----
        if use_clustsim:
            csim_out = spgo(
                f"{sing_afni} "
                f"1d_tool.py -infile {output_folder}/TTnew.CSimA.NN1_2sided.1D "
                f"-csim_show_clustsize -verb 0 -csim_pthr {csim_pthr} -csim_alpha {alpha}"
            )
            try:
                cluster_size = next(
                    int(l.strip()) for l in csim_out.splitlines() if l.strip().isdigit()
                )
            except StopIteration:
                print(f"  [WARN] ClustSim parse failed, using default={cluster_size_default}")
                cluster_size = cluster_size_default
        else:
            cluster_size = cluster_size_default

        # ---- 6) Threshold and prepare images (identical to the original) ----
        zmap_img = nib.load(zmap)
        zmap_data = zmap_img.get_fdata()
        vmax_z = np.percentile(np.abs(zmap_data), vmax_percentile) if np.any(zmap_data) else 5.0
        z_score = _scipy_norm.ppf(1 - alpha / 2)
        print(f"  vmax={vmax_z:.3f}  alpha={alpha}  z_score={z_score:.3f}")

        stat_mask = np.abs(zmap_data) >= z_score

        thr_data = zmap_data.copy()
        thr_data[~stat_mask] = 0
        thr_img = nib.Nifti1Image(thr_data, zmap_img.affine, zmap_img.header)
        thr_path = opj(output_folder, f"{seed_name}{tag}_thresholded_ttest-stat_fisher.nii.gz")
        thr_img.to_filename(thr_path)

        corr_img = nib.load(correlation)
        corr_data = corr_img.get_fdata()
        thr_corr_data = corr_data.copy()
        thr_corr_data[~stat_mask] = 0
        thr_cmap = nib.Nifti1Image(thr_corr_data, corr_img.affine, corr_img.header)
        thr_cmap_path = opj(output_folder, f"{seed_name}{tag}_thresholded_ttest-stat_fisher_cmap.nii.gz")
        thr_cmap.to_filename(thr_cmap_path)

        vmax_corr = np.max(np.abs(corr_data)) if np.any(corr_data) else 1.0

        # ---- 7) Plots (identical to the original) ----
        _safe_plot(
            correlation, output_folder, studytemplatebrain,
            f"{seed_name} mean Fisher-z (unthresholded)",
            "cold_hot" if np.any(corr_data < 0) and np.any(corr_data > 0) else "hot",
            cut_coords, vmax=vmax_corr
        )
        _safe_plot(
            thr_cmap_path, output_folder, studytemplatebrain,
            f"{seed_name} Fisher-z (p<{alpha}, cluster>={cluster_size})",
            "cold_hot" if np.any(thr_corr_data < 0) and np.any(thr_corr_data > 0) else "hot",
            cut_coords, vmax=vmax_corr, threshold=1e-6
        )
        _safe_plot(
            thr_path, output_folder, studytemplatebrain,
            f"{seed_name} t-stat z-map (p<{alpha}, cluster>={cluster_size})",
            "cold_hot" if np.any(thr_data < 0) and np.any(thr_data > 0) else "hot",
            cut_coords, vmax=vmax_z, threshold=z_score
        )

        if normalize_target is not None:
            norm_data = _normalize_map(corr_data, normalize_target)
            if norm_data is not None:
                norm_path = opj(output_folder, f"{seed_name}{tag}_zstat_norm.nii.gz")
                nib.save(nib.Nifti1Image(norm_data, corr_img.affine), norm_path)
                nz = norm_data[norm_data != 0]
                cmap_norm = "cold_hot" if (np.any(nz < 0) and np.any(nz > 0)) else "hot"
                vmax_norm = min(1.0, np.percentile(np.abs(nz), vmax_percentile)) if nz.size else 1.0
                _safe_plot(norm_path, output_folder, studytemplatebrain,
                          f"{seed_name} (Normalized, med={normalize_target:.2f})",
                          cmap_norm, cut_coords, vmax=vmax_norm)

        if normalize_target is not None:
            norm_data = _normalize_map(corr_data, normalize_target)
            if norm_data is not None:
                norm_path = opj(output_folder, f"{seed_name}{tag}_correlation_norm.nii.gz")
                nib.save(nib.Nifti1Image(norm_data, corr_img.affine), norm_path)
                nz = norm_data[norm_data != 0]
                cmap_norm = "cold_hot" if (np.any(nz < 0) and np.any(nz > 0)) else "hot"
                vmax_norm = min(1.0, np.percentile(np.abs(nz), vmax_percentile)) if nz.size else 1.0
                _safe_plot(
                    norm_path, output_folder, studytemplatebrain,
                    f"{seed_name} Normalized Fisher-z corr (p{normalize_target:.0f})",
                    cmap_norm, cut_coords, vmax=vmax_norm
                )

        plt.close("all")


def rerun_balanced_by_run(df_class, species, bids_dir, seed=None, results_base_dir=None,
                          output_suffix_kept="Full_KeptRunBalanced",
                          output_suffix_outlier="Full_OutlierRunBalanced",
                          category="All", random_seed=42, exclude_lateralized=True,
                          render_surface=True, balance=True, unit="subject"):
    """
    Like rerun_balanced_kept_vs_outlier, but balances and reruns at the
    RUN level instead of subject-session level.

    seed: which seed/region to actually run 3dttest++ for. REQUIRED to
    be the exact seed you care about — previously this used
    fig5.REGIONS_OF_INTEREST (the module's FULL list) instead, which
    made _3dttest_EDNiX also try processing every OTHER region in that
    list (e.g. Posterior_medial_cortex__PMC alongside the retrosplenial
    seed) even though that other region was never the actual SBA seed
    for these subjects — confirmed real symptom: "[WARN] No Fisher-z
    maps: subject=X seed=L_Posterior_medial_cortex__PMC" for every
    subject. None (default) auto-derives it from df_class's own 'seed'
    column (every row should share one) — printed either way so it's
    always visible which seed is actually being reprocessed.

    exclude_lateralized: True (default) — matches the confirmed real
    convention for how the original group results were produced
    (process_species_category's CLI default). _3dttest_EDNiX's own seed
    filter is a SUBSTRING match, not exact: querying the unprefixed
    "BA_26_29_30__retrosplenial" also matches "L_BA_..." and "R_BA_..."
    label_df rows (the unprefixed name is literally a substring of
    both) — so without exclude_lateralized=True, an unprefixed seed
    query silently ALSO reprocesses both lateralized regions alongside
    the merged one. Confirmed real symptom: asking for
    "BA_26_29_30__retrosplenial" produced separate
    L_BA_26_29_30__retrosplenial/ and R_BA_26_29_30__retrosplenial/
    output folders too. Set this False only if your seed query itself
    IS a lateralized name ("L_BA_..." / "R_BA_...") — exclude_lateralized
    drops EVERY lateralized row, including the one you asked for.

    ⚠ Architectural limitation of the underlying pipeline (confirmed from
    the real source of process_species_category / collect_bold_and_means):
    list_to_keep there is a list of (subject, session) pairs ONLY — there
    is no run-level filtering hook in that path at all. To actually
    balance/exclude by individual run, this function bypasses
    process_species_category entirely: it calls collect_bold_and_means()
    itself for the union of subject-sessions involved, manually filters
    the returned bold_paths dict down to the exact runs in each balanced
    group (bold_paths already carries subject/session/RUN per file — that
    granularity exists, it's just never exposed as a filter upstream),
    then calls _3dttest_EDNiX directly — replicating the same setup
    process_species_category itself does for analysis 1 (3dttest++) ONLY
    (the two 3dLMEr analyses aren't run here; they don't have a
    per-run-filterable entry point either, and weren't what run-level
    balancing was asked for).

    ⚠ UNVERIFIED beyond load_classification_subjects (confirmed): this
    assumes Figure_5_SBA.py exposes the same helpers/constants shown in
    run_SBA_by_category.py's process_species_category — build_species_cfg,
    ATLAS_LIB, MAIN_PATH, get_atlas_label_path, parse_label_file,
    find_template_low, REGIONS_OF_INTEREST, NORMALIZE_PERCENTILE,
    collect_bold_and_means, SING_WB, and that Tools.niftotoWBsurface is
    importable the same way. Very likely true given everything else
    confirmed about this file so far, but if you hit an AttributeError on
    one of these, that's exactly what it means — tell me which one and
    I'll fix the exact call.

    A subject-session's mean image (mean_imgs) is NOT split by run — that
    quantity is inherently per-session in this pipeline, not per-run, so
    both balanced groups reuse whatever mean images collect_bold_and_means
    finds for the sessions involved.

    Returns (success_kept, success_outlier, kept_runs_used, outlier_runs_used).

    render_surface: True (default, unchanged behavior) — after the
    3dttest++ group stats, also do the wb_command volume-to-surface
    rendering step (niftotoWBsurface.render_stat_maps). Set False to run
    ONLY the statistical test (3dttest++, still full-cost AFNI/ClustSim)
    and skip JUST the surface rendering that follows it — e.g. for
    intermediate steps of an iterative classification where you want the
    real group stats at every step but only want surface reconstruction
    for the FINAL, converged result.
    """
    try:
        from Statistics.Group_fMRI._Group_anal_3dTtest import _format_seed
    except ImportError as e:
        raise ImportError(
            "Could not import _3dttest_EDNiX/_format_seed from "
            f"Statistics.Group_fMRI._Group_anal_3dTtest — required for run-level "
            f"reruns. Original error: {e}"
        )
    import shutil

    if seed is None:
        seeds_in_df = (df_class["seed"].astype(str).unique().tolist()
                      if "seed" in df_class.columns else [])
        if len(seeds_in_df) != 1:
            raise ValueError(
                f"Could not auto-derive a single seed from df_class (found "
                f"{seeds_in_df}) — pass seed=... explicitly."
            )
        seed = seeds_in_df[0]
    print(f"  Restricting this rerun to seed: {seed!r}")

    results_base_dir = results_base_dir or fig5.OUTPUT_BASE

    df_class = df_class.copy()
    df_class["subject"] = df_class["subject"].astype(str)
    df_class["session"] = df_class["session"].astype(str)
    df_class["run"] = df_class["run"].astype(str)

    kept_runs = [tuple(r) for r in
                df_class.loc[~df_class["is_outlier"], ["subject", "session", "run"]].values]
    outlier_runs = [tuple(r) for r in
                    df_class.loc[df_class["is_outlier"], ["subject", "session", "run"]].values]

    def _unit_of(t):
        """(subject, session, run) -> the tuple that identifies its unit."""
        if unit == "subject":
            return (t[0],)
        if unit == "session":
            return (t[0], t[1])
        if unit == "run":
            return t
        raise ValueError(f"Unknown unit {unit!r} -- use 'run', 'session' or 'subject'")

    kept_units = sorted({_unit_of(t) for t in kept_runs})
    outlier_units = sorted({_unit_of(t) for t in outlier_runs})

    if balance:
        # Balancing at UNIT level, not run level: two groups with the
        # same number of runs but 2 animals vs 12 animals are not
        # comparable, and matching run counts hides that.
        n = min(len(kept_units), len(outlier_units))
        if n == 0:
            raise ValueError(f"One group is empty (kept={len(kept_units)}, "
                              f"outlier={len(outlier_units)} {unit}(s)) — can't balance. "
                              f"Pass balance=False if you deliberately want a map over "
                              f"one group only.")
        rng = np.random.default_rng(random_seed)

        def _subsample(items, n):
            if len(items) <= n:
                return set(items)
            idx = rng.choice(len(items), size=n, replace=False)
            return {items[i] for i in idx}

        kept_sel = _subsample(kept_units, n)
        outlier_sel = _subsample(outlier_units, n)
        print(f"  Balancing kept/outlier groups to n={n} {unit}(s) each "
              f"(kept had {len(kept_units)}, outlier had {len(outlier_units)})")
    else:
        # No subsampling at all, and an empty group is allowed -- its
        # 3dttest is simply skipped at the end. This is what a map over
        # EVERY unit needs. The previous workaround (mark one run as a
        # fake outlier so the function accepts the call) silently set
        # n = min(N-1, 1) = 1, so the "all runs" map was computed on a
        # single run and AFNI refused it, leaving an empty folder.
        kept_sel, outlier_sel = set(kept_units), set(outlier_units)
        if not kept_sel and not outlier_sel:
            raise ValueError("Both groups are empty — nothing to run.")
        print(f"  balance=False: no subsampling "
              f"(kept={len(kept_sel)} {unit}(s), outlier={len(outlier_sel)} {unit}(s))")

    kept_balanced = {t for t in kept_runs if _unit_of(t) in kept_sel}
    outlier_balanced = {t for t in outlier_runs if _unit_of(t) in outlier_sel}
    print(f"  kept ({len(kept_sel)} {unit}(s), {len(kept_balanced)} run(s)):    "
          f"{sorted(kept_sel)}")
    print(f"  outlier ({len(outlier_sel)} {unit}(s), {len(outlier_balanced)} run(s)): "
          f"{sorted(outlier_sel)}")

    all_subj_ses = sorted({(s, se) for s, se, _ in kept_balanced | outlier_balanced})
    bold_paths_all, mean_imgs_all = fig5.collect_bold_and_means([bids_dir], all_subj_ses)

    def _filter_bold_paths(target_runs):
        # target_runs holds (subj, ses, run) tuples where run may be "NA"
        # — compute_all_run_profiles's fallback when no run-<N> token
        # exists in the source correlation-map filename. But
        # collect_bold_and_means has its OWN, DIFFERENT fallback for the
        # exact same "no run token in this filename" situation: it
        # hardcodes run='01' (confirmed real source). Same real-world
        # case, two different sentinel strings — so an exact 3-tuple
        # match silently finds nothing (confirmed real symptom: "Matched
        # 0/7 kept run(s) and 0/7 outlier run(s)"). Treat "NA" in the
        # target as "match on (subject, session) alone" instead.
        target_exact = {(s, se, r) for s, se, r in target_runs}
        target_na_subj_ses = {(s, se) for s, se, r in target_runs if r == "NA"}

        out = {"bold_path": [], "subject": [], "session": [], "run": []}
        for bp, subj, ses, run in zip(bold_paths_all["bold_path"], bold_paths_all["subject"],
                                       bold_paths_all["session"], bold_paths_all["run"]):
            key = (str(subj), str(ses), str(run))
            if key in target_exact or (str(subj), str(ses)) in target_na_subj_ses:
                out["bold_path"].append(bp)
                out["subject"].append(subj)
                out["session"].append(ses)
                out["run"].append(run)
        return out

    bold_paths_kept = _filter_bold_paths(kept_balanced)
    bold_paths_outlier = _filter_bold_paths(outlier_balanced)

    n_k, n_o = len(kept_balanced), len(outlier_balanced)
    print(f"  Matched {len(bold_paths_kept['bold_path'])}/{n_k} kept run(s) and "
          f"{len(bold_paths_outlier['bold_path'])}/{n_o} outlier run(s) to real BOLD files")
    if len(bold_paths_kept["bold_path"]) < n_k or len(bold_paths_outlier["bold_path"]) < n_o:
        print(f"  [WARN] not every balanced run was found by collect_bold_and_means — "
              f"the run-number format used here may not exactly match theirs")

    cfg = fig5.build_species_cfg(species)
    atlas_ednix = opj(fig5.ATLAS_LIB, cfg["species_fragment"], "EDNiX")
    base_ss = opj(atlas_ednix, "volumes", f"{species}_space-acpc_desc-SS_T1w.nii.gz")
    mask_func = opj(atlas_ednix, "volumes", "masks", f"{species}_desc-Cerebral_Gray_mask.nii.gz")
    template_low = fig5.find_template_low([bids_dir], species)
    # REVERTED to "EDNIxCSC" (NOT "EDNIxCSCLR") — my earlier change to
    # "EDNIxCSCLR" here was wrong. Confirmed real source:
    # process_species_category (run_SBA_by_category.py) ALWAYS resolves
    # this SPECIFIC label file (the one used for _3dttest_EDNiX's seed
    # matching) via atlas_name="EDNIxCSC", unconditionally, for every
    # species — never "EDNIxCSCLR". Those are two DIFFERENT atlas
    # artifacts for two different purposes: "EDNIxCSCLR" (with "LR") is
    # the volumetric 4D dseg + its ITK label file used elsewhere in this
    # module for per-region MASKING (find_atlas_volume_path,
    # pick_representative_atlas_volume — correct to keep "EDNIxCSCLR"
    # there); "EDNIxCSC" (no "LR") is what _3dttest_EDNiX itself uses to
    # decide WHICH SEED to correlate against. Using "EDNIxCSCLR" here
    # confirmed real symptom: for species whose EDNIxCSCLR_label.txt has
    # ONLY lateralized L_.../R_... rows (no merged entry at all for a
    # region), an unprefixed seed query can only ever match those
    # lateralized rows — producing separate L_.../R_... 3dttest++ output
    # folders even when the ORIGINAL group results (and REGIONS_OF_INTEREST)
    # use the merged, unprefixed convention — and then render_stat_maps
    # can't find a literal "<seed>" folder that was never created,
    # hence "[WARN] Result dir not found, skipping: .../BA_..._retrosplenial".
    label_path = fig5.get_atlas_label_path(fig5.ATLAS_LIB, cfg["species_fragment"], "EDNIxCSC")
    # Which parser actually works varies BY FILE, not just by species —
    # confirmed: EDNIxCSC_label.txt (no "LR", used here, matching
    # process_species_category's own default) parses fine with
    # fig5.parse_label_file (genuine StatsLUT format) but produced 0
    # matches with parse_itk_snap_label_file; the SEPARATE
    # EDNIxCSCLR_ITK.txt file used elsewhere in this module for the 4D
    # volumetric region masking is genuine ITK-SNAP format and needs the
    # opposite parser. Auto-detect instead of hardcoding either one.
    label_df = _load_label_df_for_seed_matching(label_path, fig5.REGIONS_OF_INTEREST)

    # Pre-flight check: replicate _3dttest_EDNiX's own seed filter here —
    # including its _format_seed() normalization — so a bad label_df is
    # caught with a clear message BEFORE launching AFNI, instead of it
    # silently producing nothing. Checked against THIS call's actual
    # restriction ([seed]), not fig5.REGIONS_OF_INTEREST's full list —
    # that broader list would misleadingly show matches (e.g. PMC) that
    # _3dttest_EDNiX will never actually be asked to process here.
    matched_seeds = label_df[label_df["region_name"].apply(
        lambda b: any(r in _format_seed(b) for r in [seed]))]
    if matched_seeds.empty:
        raise RuntimeError(
            f"0/{len(label_df)} regions in label_df match seed={seed!r} — "
            f"_3dttest_EDNiX would silently produce nothing. Check label_path "
            f"({label_path})."
        )

    lateralized = matched_seeds[matched_seeds["region_name"].apply(
        lambda b: _format_seed(b).startswith(("L_", "R_")))]
    merged = matched_seeds[~matched_seeds.index.isin(lateralized.index)]
    print(f"  {len(matched_seeds)} region(s) in label_df matched seed={seed!r}: "
          f"{matched_seeds['region_name'].tolist()} "
          f"({len(merged)} merged/bilateral, {len(lateralized)} lateralized)")

    if exclude_lateralized and merged.empty and not lateralized.empty:
        print(f"  [WARN] exclude_lateralized=True but label_df has ONLY lateralized "
              f"matches for seed={seed!r} — after exclude_lateralized filtering, "
              f"_3dttest_EDNiX will have 0 seeds left and will silently produce "
              f"nothing. Try exclude_lateralized=False for this species/label_path.")
    elif not exclude_lateralized and not lateralized.empty and not merged.empty:
        print(f"  [INFO] both merged and lateralized rows matched seed={seed!r} — "
              f"with exclude_lateralized=False, _3dttest_EDNiX will process BOTH "
              f"the merged region AND each hemisphere separately (multiple output "
              f"folders). Set exclude_lateralized=True if you only want the "
              f"merged one.")

    # Same surface_cfg process_species_category itself builds — needed so
    # the group result gets rendered/visualized, not just left as raw
    # NIfTI files (which is exactly what was missing: files existed but
    # nothing was viewable in them without this render step).
    # ⚠ UNVERIFIED beyond what's shown in run_SBA_by_category.py: assumes
    # Figure_5_SBA.py also has a SING_WB constant and Tools.niftotoWBsurface
    # is importable the same way — tell me if either name is wrong.
    try:
        from Tools import niftotoWBsurface
    except ImportError as e:
        raise ImportError(
            f"Could not import niftotoWBsurface from Tools — needed to render "
            f"the group result into viewable images. Original error: {e}"
        )
    sing_wb = getattr(fig5, "SING_WB", None)
    if sing_wb is None:
        raise AttributeError(
            "Figure_5_SBA module has no SING_WB constant — needed to render "
            "surface visualizations. Tell me what it's actually called there "
            "(or its value) and I'll fix this."
        )
    scene_template = opj(atlas_ednix, "surfaces", "Native_resol", f"{species}_thickness.scene")
    surface_cfg = dict(
        species_name=species, atlas_lib_path=fig5.ATLAS_LIB,
        species_fragment=cfg["species_fragment"], sing_wb=sing_wb,
        surface_type="inflated", palette_name="RY-BC-BL", pos_max=90, neg_max=10,
        colorbar_mode="abs_pct", scene_template=scene_template,
        scene_label="New Scene 1", display_zero=False,
    )
    if cfg.get("projection") is not None:
        surface_cfg["projection"] = cfg["projection"]
    if cfg.get("delete_surf_file") is not None:
        surface_cfg["delete_surf_file"] = cfg["delete_surf_file"]

    shared_base = dict(
        templatehigh=base_ss, templatelow=template_low, oversample_map=False,
        mask_func=mask_func, label_df=label_df, MAIN_PATH=fig5.MAIN_PATH,
        # Restricted to THIS call's seed only — NOT fig5.REGIONS_OF_INTEREST
        # (the module's full list, e.g. also including PMC). Passing the
        # full list makes _3dttest_EDNiX try to process every region in
        # it, including seeds that were never actually used to generate
        # this run's Fisher-z correlation maps — confirmed real symptom:
        # "[WARN] No Fisher-z maps: subject=X seed=L_Posterior_medial_cortex__PMC"
        # for every subject, because PMC was never the SBA seed here.
        regions_of_interest=[seed], alpha=cfg["alpha"],
        lower_cutoff=cfg["lower_cutoff"], upper_cutoff=cfg["upper_cutoff"],
        cut_coords=cfg["cut_coords"], normalize_target=fig5.NORMALIZE_PERCENTILE,
        bids_dir=bids_dir, exclude_lateralized=exclude_lateralized,
    )
    if cfg.get("method_mask_func") is not None:
        shared_base["method_mask_func"] = cfg["method_mask_func"]

    bids_name = Path(bids_dir).name
    category_tag = f"{category}"

    def _run_3dttest(bold_paths, output_suffix, label):
        output_dir = opj(results_base_dir, species, bids_name, output_suffix, category)
        os.makedirs(output_dir, exist_ok=True)
        print(f"  Running 3dttest++ for {label} ({len(bold_paths['bold_path'])} runs) "
              f"-> {output_dir}")
        # Pass output_suffix INTO _3dttest_EDNiX_run_level itself, so each
        # group gets its own physically separate internal working folder
        # (bids_dir/Results/Grp_SBA_3dTTEST_<output_suffix>/<seed>/) instead
        # of both groups sequentially reusing the exact same
        # Grp_SBA_3dTTEST/<seed>/ folder.
        #
        # Uses _3dttest_EDNiX_run_level (this module's own variant), NOT
        # the imported _3dttest_EDNiX — the real _3dttest_EDNiX ALWAYS
        # pre-averages every run of a subject into one Fisher-z map
        # before 3dttest++ (confirmed real source: the "3dMean" step,
        # producing one "<subject>_avg_fisher_map.nii.gz" per unique
        # subject), which defeats the entire point of balancing groups
        # at the RUN level — 3dttest++'s degrees of freedom would still
        # reflect subject count, not run count. _3dttest_EDNiX_run_level
        # skips that pre-averaging so each run is its own -setA input.
        _3dttest_EDNiX_run_level(**shared_base, bold_paths=bold_paths, mean_imgs=mean_imgs_all,
                                 output_suffix=output_suffix, unit=unit)
        # Matches _3dttest_EDNiX_run_level's own internal naming exactly:
        # folder_name = "Grp_SBA_3dTTEST" + (f"_{output_suffix}" if output_suffix else "")
        source_dir = opj(bids_dir, "Results", f"Grp_SBA_3dTTEST_{output_suffix}")
        # Drop the redundant category_tag prefix here — category ("All")
        # is already a parent directory level in output_dir, so baking
        # it into the folder name too just makes an already-long name
        # longer for no added information. Now matches source_dir's own
        # naming exactly (both "Grp_SBA_3dTTEST_<output_suffix>").
        results_subdir = f"Grp_SBA_3dTTEST_{output_suffix}"
        target_dir = opj(output_dir, results_subdir)
        if not ope(source_dir):
            print(f"  [WARN] expected source dir not found after _3dttest_EDNiX: {source_dir}")
            return False

        if ope(target_dir):
            shutil.rmtree(target_dir)
        shutil.copytree(source_dir, target_dir)
        print(f"  [saved] {target_dir}")

        # zmap/cmap use their OVERSAMPLED variants for the surface render
        # (see _3dttest_EDNiX_run_level's oversample_for_surface — smoother
        # volume-to-surface projection, since surface mesh vertices are
        # much denser than this functional-resolution volume). correlation_norm
        # is left as-is, unchanged from before.
        if render_surface:
            stat_suffixes = ["ttest-stat_fisher_zmap_oversampled", "ttest-stat_fisher_cmap_oversampled"]
            if fig5.NORMALIZE_PERCENTILE is not None:
                stat_suffixes.append("correlation_norm")
            niftotoWBsurface.render_stat_maps(
                results_dir=output_dir, results_subdir=results_subdir,
                stat_suffixes=stat_suffixes, regions_of_interest=[seed],
                surface_cfg=surface_cfg,
            )
            print(f"  [rendered] {target_dir}")
        else:
            print(f"  [skipped surface] {target_dir} (render_surface=False -- stats only)")
        return True

    if bold_paths_kept["bold_path"]:
        success_kept = _run_3dttest(bold_paths_kept, output_suffix_kept, "kept")
    else:
        print("  [skip] kept group is empty — no 3dttest launched for it")
        success_kept = False
    if bold_paths_outlier["bold_path"]:
        success_outlier = _run_3dttest(bold_paths_outlier, output_suffix_outlier, "outlier")
    else:
        print("  [skip] outlier group is empty — no 3dttest launched for it")
        success_outlier = False

    return success_kept, success_outlier, kept_balanced, outlier_balanced


# =============================================================================
# 3) Reference ("population") curve, from the existing group-level SBA map
# =============================================================================

def build_reference_curve(group_map_path, atlas_vol_path, label_df, stats=DEFAULT_REGION_STATS,
                           atlas_level=None, restrict_to_nonzero_map=True):
    """The group map's own per-region profile — "the standard for the DMN"."""
    return compute_region_profile(group_map_path, atlas_vol_path, label_df, stats=stats,
                                  atlas_level=atlas_level,
                                  restrict_to_nonzero_map=restrict_to_nonzero_map)


# =============================================================================
# 4) Deviation of a run's profile from the reference curve
# =============================================================================

def profile_distance(run_profile, ref_profile, value_col="median_r", metric="correlation"):
    """
    Distance between one run's per-region profile and the reference
    profile, on the SAME region ordering (merged on region_name so a
    missing region in either doesn't misalign the vectors).
    Reuses qc._distance — same metric conventions as the rest of the
    pipeline (correlation / cosine / euclidean).

    Defensive dedup + validate="one_to_one": a region_name that appears
    more than once in either profile (e.g. from a label file with
    duplicate/malformed rows) would otherwise make merge() silently
    produce a cartesian product — inflating n_regions far past the real
    region count instead of erroring. drop_duplicates keeps the first
    occurrence; if you're seeing this warning a lot, the label_df itself
    likely needs fixing at the source.
    """
    run_p = run_profile[["region_name", value_col]]
    ref_p = ref_profile[["region_name", value_col]]
    if run_p["region_name"].duplicated().any() or ref_p["region_name"].duplicated().any():
        print(f"  [WARN] duplicate region_name in profile_distance inputs — "
              f"keeping first occurrence of each (run: {run_p['region_name'].duplicated().sum()} dup, "
              f"ref: {ref_p['region_name'].duplicated().sum()} dup)")
        run_p = run_p.drop_duplicates(subset="region_name")
        ref_p = ref_p.drop_duplicates(subset="region_name")

    merged = run_p.merge(ref_p, on="region_name", suffixes=("_run", "_ref"), validate="one_to_one")
    run_vec = merged[f"{value_col}_run"].values
    ref_vec = merged[f"{value_col}_ref"].values
    finite = np.isfinite(run_vec) & np.isfinite(ref_vec)
    dist, r = qc._distance(run_vec[finite], ref_vec[finite], metric=metric)
    return dist, r, int(finite.sum())


# =============================================================================
# 5) Orchestration — one species/BIDS at a time
# =============================================================================

def compute_all_run_profiles(species, bids_dir, label_df, classification_path,
                              category="All", seeds=RETROSPLENIAL_SEEDS,
                              stats=DEFAULT_REGION_STATS, atlas_name="EDNIxCSCLR",
                              atlas_level=None, restrict_to_nonzero_map=True):
    """
    For one species/BIDS: reload the subject-session list (same
    classification-based approach as elsewhere), find every run's raw SBA
    map for both retrosplenial seeds, and compute its region profile —
    using EACH SUBJECT'S OWN atlas volume (see find_atlas_volume_path),
    since the parcellation is stored per subject/session, not once per
    species.

    Returns a dict {(subject, session, run, seed): profile_df} plus a long
    DataFrame with one row per (run, seed, region).

    ⚠ Confirmed real fix (2024): bids_dir must be the SHORT name (e.g.
    "BIDS_Zhu_Garin"), not the full path — load_classification_subjects
    filters via df_subject['bids_dir'] == bids_dir against the Excel's
    bids_dir column, which stores only the short name. And on a no-match
    (empty dataframe after filtering) the real function returns the
    TUPLE (None, None), not a dict — calling .get() on that raises
    AttributeError, so the type is checked before calling it.
    """
    category_subjects = fig5.load_classification_subjects(
        classification_path, species=species, bids_dir=Path(bids_dir).name,
        categories=[category],
    )
    if not isinstance(category_subjects, dict) or not category_subjects.get(category):
        print(f"  [WARN] 0 subjects from load_classification_subjects for "
              f"{species}/{category} (bids_dir={Path(bids_dir).name})")
        return {}, pd.DataFrame()

    # Use classification ONLY to get the set of approved SUBJECTS — confirmed
    # on real data that load_classification_subjects collapses every subject
    # down to a single (subj, ses) entry regardless of how many sessions
    # actually exist (a Macaque BIDS with subjects having up to 11 sessions
    # each still produced ~1 profile per subject). Real session enumeration
    # is done directly against the BIDS tree below instead of trusting that
    # session field.
    approved_subjects = sorted({str(subj) for sp, bd, subj, ses in category_subjects[category]})

    subs_sessions = []
    for subj in approved_subjects:
        ses_dirs = sorted(glob.glob(opj(bids_dir, f"sub-{subj}", "ses-*")))
        if ses_dirs:
            for d in ses_dirs:
                subs_sessions.append((subj, opb(d).replace("ses-", "")))
        else:
            subs_sessions.append((subj, "1"))  # no ses-* subfolders: single-session layout

    all_rows = []
    profiles = {}
    for subj, ses in subs_sessions:
        atlas_vol_path = find_atlas_volume_path(bids_dir, subj, ses, atlas_name=atlas_name)
        if atlas_vol_path is None:
            print(f"  [WARN] no atlas volume found for sub-{subj} ses-{ses} — skipping")
            continue
        for seed_name in seeds:
            for map_path in find_raw_corr_maps_for_run(bids_dir, subj, ses, seed_name):
                run = _extract_run_id(opb(map_path))
                try:
                    prof = compute_region_profile(map_path, atlas_vol_path, label_df, stats=stats,
                                                  atlas_level=atlas_level,
                                                  restrict_to_nonzero_map=restrict_to_nonzero_map)
                except ValueError as e:
                    print(f"  [WARN] {map_path}: {e}")
                    continue
                key = (subj, ses, run, seed_name)
                profiles[key] = prof
                prof = prof.copy()
                prof["subject"], prof["session"] = subj, ses
                prof["run"], prof["seed"] = run, seed_name
                all_rows.append(prof)

    long_df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    return profiles, long_df


def classify_runs_by_profile_distance(profiles, ref_profile, value_col="median_r",
                                       metric="correlation", outlier_threshold=3.5):
    """
    Compute each run's distance to the reference curve and flag outliers
    (MAD, one-sided "upper" — only a LARGE deviation is bad, same
    reasoning as the rest of the pipeline).
    """
    rows = []
    for (subj, ses, run, seed), prof in profiles.items():
        dist, r, n_regions = profile_distance(prof, ref_profile, value_col=value_col, metric=metric)
        rows.append(dict(subject=subj, session=ses, run=run, seed=seed,
                         distance=dist, pearson_r=r, n_regions=n_regions))
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return qc.flag_outliers(df, column="distance", threshold=outlier_threshold, direction="upper")


def classify_runs_by_absolute_correlation(profiles, ref_profile, value_col="median_r",
                                           min_correlation=0.3):
    """
    Direct, ABSOLUTE acceptability threshold on the spatial (Pearson)
    correlation between each run's profile and the reference — unlike
    classify_runs_by_percentile (always flags a fixed proportion,
    regardless of whether the whole batch is actually good or bad) and
    classify_runs_by_profile_distance's MAD threshold (relative to THIS
    batch's own spread), this means the same thing whether you're
    looking at 5 runs or 500: a run is flagged if its correlation to the
    reference falls below min_correlation.

    Conventional rough guidance for Pearson r magnitude (not a hard
    rule, just an anchor for picking min_correlation): >0.7 strong,
    >0.5 moderate, >0.3 weak-to-moderate, below that generally weak.

    Always uses metric="correlation" internally (that's the whole point
    of thresholding on pearson_r specifically) — value_col still selects
    which per-region stat (median_r, frac_positive, ...) the profiles are
    built from.

    Returns the same columns as classify_runs_by_profile_distance.
    """
    rows = []
    for (subj, ses, run, seed), prof in profiles.items():
        dist, r, n_regions = profile_distance(prof, ref_profile, value_col=value_col,
                                              metric="correlation")
        rows.append(dict(subject=subj, session=ses, run=run, seed=seed,
                         distance=dist, pearson_r=r, n_regions=n_regions))
    df = pd.DataFrame(rows)
    if df.empty:
        return df

    df["is_outlier"] = df["pearson_r"] < min_correlation
    n_out = int(df["is_outlier"].sum())
    print(f"  min_correlation={min_correlation}: {n_out}/{len(df)} runs below threshold "
          f"(pearson_r < {min_correlation}), {len(df) - n_out} kept")
    return df


def classify_runs_by_percentile(profiles, ref_profile, value_col="median_r",
                                 metric="correlation", keep_fraction=0.5):
    """
    Alternative to classify_runs_by_profile_distance's MAD-based cutoff:
    directly keep only the best `keep_fraction` of runs (lowest distance
    to the reference) and flag the rest as outliers.

    Use this instead of pushing outlier_threshold toward 0 to "be more
    selective" — threshold≈0 in the MAD method already selects ~50% of
    the distribution (everything above the median) as outliers, so
    shrinking a small POSITIVE threshold further barely changes anything;
    going more negative would flag more, but that's an unintuitive way to
    control it. This function gives direct, exact control instead:
    keep_fraction=0.5 keeps the best half, 0.1 keeps only the best 10%,
    etc. — no threshold-tuning-by-trial-and-error required.

    Returns the same columns as classify_runs_by_profile_distance
    (subject, session, run, seed, distance, pearson_r, n_regions,
    is_outlier), plus the cutoff distance actually used. modified_z is
    not meaningful for this method and is left out.
    """
    rows = []
    for (subj, ses, run, seed), prof in profiles.items():
        dist, r, n_regions = profile_distance(prof, ref_profile, value_col=value_col, metric=metric)
        rows.append(dict(subject=subj, session=ses, run=run, seed=seed,
                         distance=dist, pearson_r=r, n_regions=n_regions))
    df = pd.DataFrame(rows)
    if df.empty:
        return df

    cutoff = df["distance"].quantile(keep_fraction)
    df["is_outlier"] = df["distance"] > cutoff
    n_out = int(df["is_outlier"].sum())
    print(f"  keep_fraction={keep_fraction}: cutoff distance={cutoff:.4f} -> "
          f"{n_out}/{len(df)} flagged as outliers, {len(df) - n_out} kept")
    return df


# =============================================================================
# 6) Cross-species common-region reference curve
# =============================================================================

def find_common_base_regions(label_dfs_by_species):
    """
    Base region names present in EVERY species' atlas label_df.

    ⚠ Simple name-intersection fallback. Your codebase already has a
    canonical cross-species region mapping (Atlas_nomenclature /
    get_region_to_level, imported elsewhere in ednix_bids_tools.py) — if
    that module exposes a species->species region alias table, use THAT
    for this intersection instead of bare name matching, since atlas
    naming may differ across species even for homologous regions. I don't
    know that module's exact API for this, so this is the safe local
    fallback until you tell me otherwise.
    """
    base_sets = [set(df["base_name"].str.lower().str.strip()) for df in label_dfs_by_species.values()]
    common = set.intersection(*base_sets) if len(base_sets) > 1 else base_sets[0]
    return sorted(common)


def build_interspecies_reference_curve(ref_profiles_by_species, common_base_names, value_col="median_r"):
    """
    Pool the per-species reference curves (see build_reference_curve),
    restricted to common_base_names, into ONE inter-species reference
    curve — the mean of each region's value across species.
    """
    frames = []
    for species, prof in ref_profiles_by_species.items():
        p = prof[prof["base_name"].str.lower().str.strip().isin(common_base_names)].copy()
        p["species"] = species
        frames.append(p)
    pooled = pd.concat(frames, ignore_index=True)
    interspecies = (
        pooled.groupby("base_name")[value_col]
              .mean()
              .reset_index()
              .rename(columns={value_col: value_col})
    )
    # profile_distance merges on "region_name" — reuse that column name
    # here too (base_name doubles as region_name at this pooled level,
    # since hemisphere-specific region_name no longer applies once pooled
    # across species with potentially different L_/R_ naming).
    interspecies = interspecies.rename(columns={"base_name": "region_name"})
    return interspecies


def compute_group_average_profile(profiles, keys, stats=DEFAULT_REGION_STATS):
    """
    Average a SUBSET of run profiles (selected by `keys`, a set of
    (subject, session, run, seed) tuples matching profiles' dict keys)
    into ONE profile DataFrame — mean of each stat, per region, across
    whichever runs are in `keys`. Used to build a "kept average" or
    "outlier average" connectivity profile per BIDS, e.g. for
    kept-vs-outlier comparison or as input to inter-species pooling.

    Returns None if no profile in `profiles` matches `keys` (empty group)
    rather than raising — callers decide how to handle that.
    """
    group_profiles = [prof for key, prof in profiles.items() if key in keys]
    if not group_profiles:
        return None
    combined = pd.concat(group_profiles, ignore_index=True)
    stat_cols = [s for s in stats if s in combined.columns]
    avg = combined.groupby("region_name")[stat_cols].mean().reset_index()
    meta = group_profiles[0][["region_name", "base_name", "hemisphere", "label_id"]].drop_duplicates()
    return meta.merge(avg, on="region_name", how="right")


def compute_group_median_profile(profiles, keys, stats=DEFAULT_REGION_STATS):
    """
    Same as compute_group_average_profile, but MEDIAN instead of mean per
    region across the runs in `keys` — more robust to a handful of
    still-noisy runs remaining in the group (e.g. mid-way through
    classify_runs_iterative_robust's convergence, before the kept set
    has fully stabilized, or if a group's true "clean core" isn't
    perfectly homogeneous). Same return shape (region_name, base_name,
    hemisphere, label_id, plus each stat column) as
    compute_group_average_profile, so both are interchangeable anywhere
    a reference profile is expected.

    Returns None if no profile in `profiles` matches `keys` (empty
    group), same as compute_group_average_profile.
    """
    group_profiles = [prof for key, prof in profiles.items() if key in keys]
    if not group_profiles:
        return None
    combined = pd.concat(group_profiles, ignore_index=True)
    stat_cols = [s for s in stats if s in combined.columns]
    med = combined.groupby("region_name")[stat_cols].median().reset_index()
    meta = group_profiles[0][["region_name", "base_name", "hemisphere", "label_id"]].drop_duplicates()
    return meta.merge(med, on="region_name", how="right")


# =============================================================================
# Reference construction strategies (medoid / weighted / k-core)
# =============================================================================

def compute_pairwise_similarity(profiles, keys=None, value_col="median_r",
                                metric="correlation"):
    """
    Pairwise similarity (pearson_r via profile_distance) between every
    pair of runs in `profiles` (restricted to `keys` if given).

    Shared building block for the medoid/weighted/k-core reference
    strategies below, all of which rest on the same idea: rather than
    averaging every run uniformly (which BLURS the reference when runs
    are poorly aligned or artifact-dominated — averaging only improves
    SNR when the signal is spatially aligned across runs, otherwise it
    degrades it) or hand-picking one "good-looking" subject (an
    arbitrary, hard-to-justify choice), let each run's CENTRALITY among
    all the others decide how much it defines the reference.

    Returns (sim_matrix, key_list): sim_matrix[i, j] is the correlation
    between key_list[i]'s and key_list[j]'s profiles (diagonal = 1.0,
    NaN where a pair couldn't be compared).
    """
    key_list = sorted(keys if keys is not None else profiles.keys(),
                      key=lambda k: tuple(str(x) for x in k))
    n = len(key_list)
    sim = np.full((n, n), np.nan)
    for i in range(n):
        sim[i, i] = 1.0
        for j in range(i + 1, n):
            try:
                _, r, _ = profile_distance(profiles[key_list[i]], profiles[key_list[j]],
                                           value_col=value_col, metric=metric)
            except Exception:
                r = np.nan
            sim[i, j] = sim[j, i] = r
    return sim, key_list


def compute_medoid_reference(profiles, keys=None, value_col="median_r",
                             metric="correlation", verbose=True):
    """
    Use the single MEDOID run as the reference — the run whose MEDIAN
    similarity to every other run is highest, i.e. the most "central"
    real run rather than a synthetic average of all of them.

    Motivation (confirmed by this project's own visual observation that
    a single hand-picked template subject outperformed the group average
    in the noisiest species): averaging misaligned/artifact-dominated
    runs produces a BLURRED reference whose network boundaries are
    washed out, sometimes resembling the true network LESS than the
    single best run does. The medoid captures exactly what makes that
    hand-picked run good — being typical — but selects it
    deterministically and reproducibly, with no manual choice to
    justify.

    Tradeoff vs. an average: keeps the sharpness of a real single run,
    but forgoes all noise-averaging, so it inherits whatever
    idiosyncratic noise that one run happens to carry. See
    compute_kcore_reference for the middle ground.

    Returns (ref_profile, medoid_key, sim_matrix, key_list) — the
    profile itself IS that run's own profile, unmodified.
    """
    sim, key_list = compute_pairwise_similarity(profiles, keys, value_col, metric)
    if not key_list:
        return None, None, sim, key_list
    with np.errstate(invalid="ignore"):
        median_sim = np.nanmedian(sim, axis=1)
    best_i = int(np.nanargmax(median_sim))
    medoid_key = key_list[best_i]
    if verbose:
        print(f"  Medoid reference: {medoid_key} "
             f"(median similarity to others = {median_sim[best_i]:.4f})")
    return profiles[medoid_key].copy(), medoid_key, sim, key_list


def compute_weighted_reference(profiles, keys=None, value_col="median_r",
                               metric="correlation", weight_power=1.0,
                               stats=DEFAULT_REGION_STATS, verbose=True):
    """
    Weighted-average reference: every run contributes, but PROPORTIONALLY
    to its centrality (median similarity to all other runs), instead of
    every run counting equally.

    This is the middle ground between the two failure modes: a uniform
    average lets badly-aligned/artifact-heavy runs blur the reference as
    much as good ones, while a single medoid throws away all
    noise-averaging. Here good runs dominate the reference while poor
    ones contribute little, so you keep the SNR benefit of averaging
    without the dilution.

    weight_power: exponent applied to the (clipped-at-zero) centrality
    weights — 1.0 (default) weights linearly by centrality; higher
    values (2, 3) concentrate the reference more sharply on the most
    central runs, approaching medoid-like behavior in the limit; 0
    reduces to a plain uniform average.

    Returns (ref_profile, weights_series, sim_matrix, key_list) where
    weights_series is indexed by the same keys, for inspection.
    """
    sim, key_list = compute_pairwise_similarity(profiles, keys, value_col, metric)
    if not key_list:
        return None, None, sim, key_list

    with np.errstate(invalid="ignore"):
        centrality = np.nanmedian(sim, axis=1)
    # Negative centrality means "anti-correlated with the group on
    # average" -- such a run should contribute nothing, not negatively
    # (a negative weight would actively push the reference AWAY from the
    # group pattern, which is not what "weight by quality" should mean).
    weights = np.clip(np.nan_to_num(centrality, nan=0.0), 0.0, None) ** weight_power
    if weights.sum() <= 0:
        print("  [WARN] all centrality weights are zero -- falling back to a uniform average")
        weights = np.ones(len(key_list))
    weights = weights / weights.sum()

    frames = []
    for w, key in zip(weights, key_list):
        prof = profiles[key][["region_name"] + [s for s in stats if s in profiles[key].columns]].copy()
        for s in stats:
            if s in prof.columns:
                prof[s] = prof[s] * w
        frames.append(prof)
    combined = pd.concat(frames, ignore_index=True)
    stat_cols = [s for s in stats if s in combined.columns]
    weighted = combined.groupby("region_name")[stat_cols].sum().reset_index()

    meta = profiles[key_list[0]][["region_name", "base_name", "hemisphere", "label_id"]].drop_duplicates()
    ref_profile = meta.merge(weighted, on="region_name", how="right")

    weights_series = pd.Series(weights, index=[tuple(str(x) for x in k) for k in key_list])
    if verbose:
        top = weights_series.sort_values(ascending=False).head(3)
        print(f"  Weighted reference (power={weight_power}): top contributors "
             f"{[(k, round(v, 4)) for k, v in top.items()]}")
    return ref_profile, weights_series, sim, key_list


def compute_kcore_reference(profiles, k=5, keys=None, value_col="median_r",
                            metric="correlation", use_median=True,
                            stats=DEFAULT_REGION_STATS, verbose=True):
    """
    K-core reference: average (or median) of ONLY the medoid run and its
    k-1 most similar neighbours — a small, mutually-consistent core.

    The explicit middle ground between compute_medoid_reference (sharp
    but no noise-averaging, inherits one run's idiosyncrasies) and a
    full group average (maximum noise-averaging but blurred by
    misaligned runs). With k small (3-5), you average enough runs to
    cancel independent noise while keeping every contributor close
    enough to the medoid that alignment-driven blurring stays limited.

    ⚠ Same self-reinforcing caveat as any similarity-based selection:
    picking the runs most similar to the medoid makes the resulting
    reference more similar to those runs by construction. With noise
    that is STRUCTURED rather than independent (very common in animal
    fMRI: anesthesia depth, respiratory/cardiac artifact, coil or
    sequence signatures), a shared artifact can be more mutually
    consistent than the true neural signal, so this can converge on an
    artifact core. Validate the resulting core against something
    EXTERNAL to it (motion/tSNR metrics, anatomical priors, or
    cross-species convergence) rather than trusting internal
    consistency alone.

    Returns (ref_profile, core_keys, sim_matrix, key_list).
    """
    sim, key_list = compute_pairwise_similarity(profiles, keys, value_col, metric)
    if not key_list:
        return None, [], sim, key_list

    with np.errstate(invalid="ignore"):
        median_sim = np.nanmedian(sim, axis=1)
    medoid_i = int(np.nanargmax(median_sim))

    k_eff = max(1, min(k, len(key_list)))
    sims_to_medoid = sim[medoid_i].copy()
    sims_to_medoid[medoid_i] = np.inf  # keep the medoid itself first
    order = np.argsort(-np.nan_to_num(sims_to_medoid, nan=-np.inf))
    core_idx = list(order[:k_eff])
    core_keys = [key_list[i] for i in core_idx]

    if verbose:
        print(f"  K-core reference (k={k_eff}): medoid={key_list[medoid_i]}, "
             f"core={[tuple(str(x) for x in kk) for kk in core_keys]}")
        if k_eff < k:
            print(f"  [INFO] only {len(key_list)} run(s) available -- k reduced from {k} to {k_eff}")

    build = compute_group_median_profile if use_median else compute_group_average_profile
    ref_profile = build(profiles, set(core_keys), stats=stats)
    return ref_profile, core_keys, sim, key_list


def classify_runs_iterative_robust(profiles, value_col="median_r", metric="correlation",
                                   min_correlation=None, keep_fraction=None,
                                   use_median=True, max_iter=10, verbose=True,
                                   on_iteration=None):
    """
    Iterative robust reference classification — an alternative to BOTH a
    group-average reference computed from ALL runs (contaminated by the
    very noisy runs it's meant to help exclude) AND a hand-picked
    individual reference subject (an arbitrary choice with no
    consistent external-QC advantage found — see Fig 11 in this
    project's own analysis).

    Procedure: start from the average/median profile across EVERY run,
    classify Kept/Outlier against it, then RECOMPUTE the reference using
    ONLY the runs kept in that round, reclassify everyone against the
    new (purified) reference, and repeat until the kept SET stops
    changing between two consecutive iterations (or max_iter is hit).
    Each round's reference is a little less contaminated by noisy runs
    than the last, converging toward a "clean core" reference — with no
    subject choice, deterministic, and reproducible from the data alone.

    Stays CHEAP by design: every iteration only re-averages/re-medians
    ALREADY-COMPUTED per-run profiles (compute_group_average_profile /
    compute_group_median_profile — plain pandas groupby) and
    reclassifies — it never reruns any AFNI group analysis or surface
    rendering. If you want a final group map / surface rendering built
    from the CONVERGED kept set specifically, that stays a SEPARATE,
    explicit step (rerun_balanced_by_run / the existing
    RERUN_GROUP_ANALYSIS pipeline) — run it ONCE after this converges,
    never inside the loop.

    profiles: {(subject, session, run, seed): profile_df} from
    compute_all_run_profiles.
    value_col / metric: same meaning as classify_runs_by_absolute_correlation
    / classify_runs_by_percentile.
    min_correlation: if given, every iteration classifies via
    classify_runs_by_absolute_correlation with this fixed threshold.
    keep_fraction: if given instead, every iteration classifies via
    classify_runs_by_percentile with this fraction. Pass EXACTLY ONE of
    min_correlation/keep_fraction.
    use_median: True (default) — build each iteration's reference with
    compute_group_median_profile (robust to whatever noisy runs remain
    in the kept set at that point). False uses
    compute_group_average_profile (plain mean) instead.
    max_iter: safety cap if convergence is never reached (e.g. an
    oscillating assignment) — returns the LAST iteration's result with a
    printed warning if hit, rather than looping forever.
    on_iteration: optional callback called as on_iteration(iteration_index,
    df_class, ref_profile) once per iteration, right after that
    iteration's classification is computed — lets a caller save/plot
    each iteration's intermediate result (e.g. to its own "iter_N"
    folder) WITHOUT this function itself knowing anything about file
    paths or plotting. None (default) skips this entirely — purely
    additive, no behavior change if unused.

    Returns dict(df_class=<final DataFrame>, ref_profile=<final,
    converged reference profile>, n_iterations=<int>, converged=<bool>,
    history=<list of every iteration's df_class, in order, for
    inspecting how the assignment evolved>).
    """
    if (min_correlation is None) == (keep_fraction is None):
        raise ValueError("Pass exactly one of min_correlation or keep_fraction.")

    build_ref = compute_group_median_profile if use_median else compute_group_average_profile

    all_keys = set(profiles.keys())
    ref_profile = build_ref(profiles, all_keys)
    if ref_profile is None:
        raise ValueError("profiles is empty -- nothing to classify.")

    history = []
    prev_kept_set = None

    for it in range(max_iter):
        if min_correlation is not None:
            df_class = classify_runs_by_absolute_correlation(
                profiles, ref_profile, value_col=value_col, min_correlation=min_correlation)
        else:
            df_class = classify_runs_by_percentile(
                profiles, ref_profile, value_col=value_col, metric=metric,
                keep_fraction=keep_fraction)

        kept_set = frozenset(
            (str(row.subject), str(row.session), str(row.run), str(row.seed))
            for row in df_class.itertuples() if not bool(row.is_outlier)
        )
        n_kept, n_total = len(kept_set), len(df_class)
        if verbose:
            print(f"  [iter {it}] {n_kept}/{n_total} kept")
        history.append(df_class)

        if on_iteration is not None:
            on_iteration(it, df_class, ref_profile)

        if prev_kept_set is not None and kept_set == prev_kept_set:
            if verbose:
                print(f"  Converged after {it + 1} iteration(s) -- kept set unchanged.")
            return dict(df_class=df_class, ref_profile=ref_profile,
                       n_iterations=it + 1, converged=True, history=history)

        if n_kept == 0:
            print(f"  [WARN] iteration {it}: 0 runs kept -- cannot recompute a reference "
                 f"from an empty set. Returning the PREVIOUS iteration's result instead.")
            prev_df_class = history[-2] if len(history) > 1 else df_class
            return dict(df_class=prev_df_class, ref_profile=ref_profile,
                       n_iterations=it, converged=False, history=history)

        prev_kept_set = kept_set
        # profiles dict keys are the RAW (subject, session, run, seed)
        # tuples (not necessarily str-cast) -- match kept_set (str-cast,
        # built from df_class) back to profiles' real keys before
        # rebuilding the reference.
        matching_keys = {k for k in profiles if tuple(str(x) for x in k) in kept_set}
        ref_profile = build_ref(profiles, matching_keys)

    print(f"  [WARN] did not converge within {max_iter} iteration(s) -- returning the last one.")
    return dict(df_class=history[-1], ref_profile=ref_profile,
               n_iterations=max_iter, converged=False, history=history)


def classify_runs_progressive_shrink(profiles, keep_fraction_schedule, min_n_kept=8,
                                     value_col="median_r", metric="correlation",
                                     use_median=True, verbose=True, on_step=None):
    """
    Graduated version of classify_runs_iterative_robust: instead of
    purifying the reference at ONE FIXED keep_fraction until it
    converges, steps through a SCHEDULE of increasingly strict
    keep_fraction values (e.g. [0.5, 0.35, 0.25, 0.10]) — at EACH step,
    the reference is rebuilt from the CURRENT kept set (same
    self-purifying mechanism as classify_runs_iterative_robust), then
    everyone is reclassified against it at the NEXT (stricter) fraction.

    Confirmed real risk this guards against: shrinking too far (e.g. all
    the way to 1-2 individuals) recreates the exact "arbitrary individual
    reference" problem this whole approach was built to avoid (no
    consistent external-QC advantage over the group average — see Fig 11
    in this project's own analysis) — AND a shrinking loop can
    self-reinforce by construction (each step keeps whoever is closest to
    the current reference, which then makes the NEXT reference even
    closer to just them), converging toward an increasingly narrow,
    self-consistent subgroup rather than necessarily "the cleanest
    signal." min_n_kept is an ABSOLUTE COUNT floor (not a fraction) so it
    protects small BIDS (e.g. 25 runs) from collapsing to a
    near-individual reference, while still letting large BIDS (e.g. 170+
    runs) shrink much further in real terms.

    profiles: {(subject, session, run, seed): profile_df} from
    compute_all_run_profiles.
    keep_fraction_schedule: list of fractions in the order to try, e.g.
    [0.5, 0.35, 0.25, 0.10] -- typically decreasing (looser to stricter),
    though nothing enforces that ordering if you have another reason.
    min_n_kept: stop BEFORE applying a schedule step that would bring
    the kept count below this absolute number -- returns the last step
    that still met the floor, with a printed reason.
    value_col / metric / use_median: same meaning as
    classify_runs_iterative_robust.
    on_step: optional callback, on_step(step_index, keep_fraction,
    df_class, ref_profile) -- called once per schedule step actually
    applied (same purpose as classify_runs_iterative_robust's
    on_iteration, for saving/plotting each step's intermediate result).

    Returns dict(df_class=<final DataFrame>, ref_profile=<the reference
    that CLASSIFIED df_class -- i.e. one step "behind" the final kept
    set>, final_kept_reference=<ONE MORE average/median pass over the
    final kept set itself -- the "last average" over whoever survived
    all the way down>, keep_fraction=<the LAST fraction actually
    applied>, n_steps=<int>, stopped_reason=<"schedule_exhausted" |
    "min_n_kept_floor">, history=<list of every applied step's
    df_class, in schedule order>).
    """
    if not keep_fraction_schedule:
        raise ValueError("keep_fraction_schedule must have at least one value.")

    build_ref = compute_group_median_profile if use_median else compute_group_average_profile

    all_keys = set(profiles.keys())
    ref_profile = build_ref(profiles, all_keys)
    if ref_profile is None:
        raise ValueError("profiles is empty -- nothing to classify.")

    history = []
    last_df_class, last_ref_profile, last_fraction = None, ref_profile, None
    stopped_reason = "schedule_exhausted"

    for step_i, frac in enumerate(keep_fraction_schedule):
        df_class = classify_runs_by_percentile(
            profiles, ref_profile, value_col=value_col, metric=metric, keep_fraction=frac)

        kept_set = frozenset(
            (str(row.subject), str(row.session), str(row.run), str(row.seed))
            for row in df_class.itertuples() if not bool(row.is_outlier)
        )
        n_kept, n_total = len(kept_set), len(df_class)

        if n_kept < min_n_kept:
            if verbose:
                print(f"  [step {step_i}] keep_fraction={frac}: would keep only {n_kept} "
                     f"(< min_n_kept={min_n_kept}) -- stopping BEFORE this step, keeping "
                     f"the previous one.")
            stopped_reason = "min_n_kept_floor"
            break

        if verbose:
            print(f"  [step {step_i}] keep_fraction={frac}: {n_kept}/{n_total} kept")
        history.append(df_class)
        if on_step is not None:
            on_step(step_i, frac, df_class, ref_profile)

        last_df_class, last_ref_profile, last_fraction = df_class, ref_profile, frac

        # rebuild the reference from THIS step's kept set, for the next
        # (stricter) step -- profiles dict keys may not be str-cast the
        # same way df_class's columns are; match consistently. This is
        # ALSO exactly "final_kept_reference" below if the loop stops
        # right after this (schedule exhausted or next step hits the
        # floor) -- one more average/median pass over whoever is left.
        matching_keys = {k for k in profiles if tuple(str(x) for x in k) in kept_set}
        ref_profile = build_ref(profiles, matching_keys)

    if last_df_class is None:
        raise ValueError(
            f"Even the LOOSEST fraction in keep_fraction_schedule "
            f"({keep_fraction_schedule[0]}) keeps fewer than min_n_kept={min_n_kept} "
            f"runs -- lower min_n_kept, loosen the schedule, or check this BIDS has "
            f"enough runs at all."
        )

    if verbose:
        print(f"  Progressive shrink stopped ({stopped_reason}) after {len(history)} "
             f"step(s) -- final keep_fraction={last_fraction}, "
             f"{int((~last_df_class['is_outlier']).sum())}/{len(last_df_class)} kept.")

    return dict(df_class=last_df_class, ref_profile=last_ref_profile,
               final_kept_reference=ref_profile,
               keep_fraction=last_fraction, n_steps=len(history),
               stopped_reason=stopped_reason, history=history)


def plot_region_boxplots_by_category(profiles, ref_profile, df_class, region_patterns,
                                     value_col="median_r", title=None, output_path=None,
                                     dpi=300, jitter=0.08, random_seed=0):
    """
    Boxplot(s) of value_col for specific named regions (e.g.
    "Somatosensory", "Prefrontal"), split Kept vs Outlier, with
    individual run points overlaid (jittered) and the REFERENCE
    profile's own value for that region marked as a distinctly colored
    diamond (not part of either box — the reference isn't itself a
    classified run).

    region_patterns: list of (label, substring) tuples, e.g.
    [("Somatosensory", "Somatosensory"), ("Prefrontal", "Prefrontal")]
    — substring matched CASE-INSENSITIVELY against ref_profile's
    region_name. Region-name-based matching is used instead of an atlas
    level number, since compute_region_profile's own per-region level
    auto-detection already makes "which level a region lives on" an
    implementation detail, not something callers should need to know.
    If a pattern matches multiple regions (e.g. both hemispheres, or
    several sub-areas of "Prefrontal"), each run's value is the MEAN
    across every matching region for that run — same rule applied to
    the reference.

    profiles: {(subject, session, run, seed): profile_df} from
    compute_all_run_profiles — keys are str-cast and matched against
    df_class's own subject/session/run/seed columns, same convention as
    plot_kept_vs_outlier_average.

    Returns (fig, per_pattern_summary) where per_pattern_summary is
    {pattern_label: dict(kept=[...], outlier=[...], reference=value or
    None, matched_regions=[...])} — the raw numbers behind each panel,
    for a stats test or a supplementary table.
    """
    profiles_str_keys = {tuple(str(x) for x in k): v for k, v in profiles.items()}
    rng = np.random.default_rng(random_seed)

    n_patterns = len(region_patterns)
    fig, axes = plt.subplots(1, n_patterns, figsize=(4.5 * n_patterns, 5), squeeze=False)
    axes = axes[0]

    summary = {}

    for ax, (pat_label, pat_substr) in zip(axes, region_patterns):
        matched_regions = ref_profile[
            ref_profile["region_name"].str.contains(pat_substr, case=False, na=False)
        ]["region_name"].tolist()

        if not matched_regions:
            ax.axis("off")
            ax.set_title(f"{pat_label}\n(no region matched '{pat_substr}')", fontsize=9, color="#999999")
            summary[pat_label] = dict(kept=[], outlier=[], reference=None, matched_regions=[])
            continue

        ref_vals = ref_profile[ref_profile["region_name"].isin(matched_regions)][value_col]
        ref_val = float(ref_vals.mean()) if not ref_vals.empty and ref_vals.notna().any() else None

        kept_vals, outlier_vals = [], []
        for row in df_class.itertuples():
            key = (str(row.subject), str(row.session), str(row.run), str(row.seed))
            prof = profiles_str_keys.get(key)
            if prof is None:
                continue
            vals = prof[prof["region_name"].isin(matched_regions)][value_col].dropna()
            if vals.empty:
                continue
            v = float(vals.mean())
            (outlier_vals if bool(row.is_outlier) else kept_vals).append(v)

        summary[pat_label] = dict(kept=kept_vals, outlier=outlier_vals, reference=ref_val,
                                  matched_regions=matched_regions)

        data = [kept_vals, outlier_vals]
        positions = [1, 2]
        bp = ax.boxplot(data, positions=positions, widths=0.5, showfliers=False,
                        patch_artist=True)
        for patch, color in zip(bp["boxes"], ["#009E73", "#D55E00"]):
            patch.set_facecolor(color)
            patch.set_alpha(0.25)
            patch.set_edgecolor(color)
        for median in bp["medians"]:
            median.set_color("#333333")

        for pos, vals, color in zip(positions, data, ["#009E73", "#D55E00"]):
            if not vals:
                continue
            x_jitter = pos + rng.uniform(-jitter, jitter, size=len(vals))
            ax.scatter(x_jitter, vals, color=color, alpha=0.6, s=22,
                      edgecolors="white", lw=0.4, zorder=3)

        if ref_val is not None:
            ax.scatter([1.5], [ref_val], color="#4472C4", marker="D", s=90,
                      edgecolors="black", lw=0.8, zorder=5, label="reference")
            ax.axhline(ref_val, color="#4472C4", lw=1.0, ls=":", alpha=0.6, zorder=1)

        ax.set_xticks(positions)
        ax.set_xticklabels([f"Kept\n(n={len(kept_vals)})", f"Outlier\n(n={len(outlier_vals)})"],
                           fontsize=8)
        ax.set_ylabel(value_col, fontsize=9)
        n_matched = len(matched_regions)
        ax.set_title(f"{pat_label}" + (f" ({n_matched} regions avg.)" if n_matched > 1 else ""),
                    fontsize=9, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if ref_val is not None:
            ax.legend(fontsize=7, loc="best")

    fig.suptitle(title or f"{value_col} by region, Kept vs Outlier", fontsize=11)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig, summary


def plot_deviation_to_reference(profiles, ref_profile, df_class=None,
                                value_col="median_r", metric="correlation",
                                title=None, output_path=None, dpi=300):
    """
    For EVERY run in `profiles` (not just whichever group a
    classification kept), plot the distribution of its correlation to a
    GIVEN reference profile -- typically
    classify_runs_progressive_shrink's own "final_kept_reference" (the
    one extra average/median pass over the smallest surviving group at
    the end of a shrinking schedule): "how far is every run, including
    the ones excluded along the way, from that final clean-core
    average?" Complements the per-step classification plots by showing
    the FULL spread relative to the converged endpoint, not just a
    kept/outlier split at one particular threshold.

    profiles: {(subject, session, run, seed): profile_df}.
    ref_profile: the reference to measure distance against.
    df_class: optional -- if given, points are colored by THAT
    classification's kept/outlier status (which may reflect an
    earlier/looser threshold than ref_profile itself, e.g. the
    FIRST step's classification, useful to see how the group that
    looked "kept" early on is scattered around the eventual endpoint).
    None (default) plots every run in one color.

    Returns a DataFrame (subject, session, run, seed, pearson_r,
    is_outlier) with every run's correlation to ref_profile, alongside
    creating/saving the figure.
    """
    outlier_lookup = {}
    if df_class is not None:
        for row in df_class.itertuples():
            key = (str(row.subject), str(row.session), str(row.run), str(row.seed))
            outlier_lookup[key] = bool(row.is_outlier)

    rows = []
    for key, prof in profiles.items():
        str_key = tuple(str(x) for x in key)
        try:
            dist, r, n = profile_distance(prof, ref_profile, value_col=value_col, metric=metric)
        except Exception:
            continue
        rows.append(dict(subject=key[0], session=key[1], run=key[2], seed=key[3],
                        pearson_r=r, is_outlier=outlier_lookup.get(str_key)))

    result_df = pd.DataFrame(rows)
    if result_df.empty:
        print("  [plot_deviation_to_reference] nothing to plot -- no runs matched")
        return result_df

    fig, ax = plt.subplots(figsize=(7, 5))
    vals_all = result_df["pearson_r"].dropna()
    bins = np.linspace(vals_all.min(), vals_all.max(), 30) if len(vals_all) > 1 else 10

    if df_class is not None and result_df["is_outlier"].notna().any():
        kept_vals = result_df.loc[result_df["is_outlier"] == False, "pearson_r"].dropna()
        outlier_vals = result_df.loc[result_df["is_outlier"] == True, "pearson_r"].dropna()
        unknown_vals = result_df.loc[result_df["is_outlier"].isna(), "pearson_r"].dropna()
        ax.hist(kept_vals, bins=bins, color="#009E73", alpha=0.7, label=f"Kept (n={len(kept_vals)})")
        ax.hist(outlier_vals, bins=bins, color="#D55E00", alpha=0.7, label=f"Outlier (n={len(outlier_vals)})")
        if len(unknown_vals):
            ax.hist(unknown_vals, bins=bins, color="#888888", alpha=0.5,
                   label=f"Unclassified (n={len(unknown_vals)})")
        ax.legend(fontsize=8)
    else:
        ax.hist(vals_all, bins=bins, color="#4472C4", alpha=0.75)

    ax.set_xlabel("correlation to final reference (pearson_r)")
    ax.set_ylabel("runs")
    ax.set_title(title or "Deviation of every run to the final reference")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return result_df


def plot_kept_vs_outlier_average(profiles, df_class, value_col="median_r",
                                  title="Kept vs outlier — average profile",
                                  output_path=None, dpi=300, drop_empty_regions=True):
    """
    Per-BIDS comparison: the KEPT runs' average profile vs the OUTLIER
    runs' average profile, region by region — shows the systematic shift
    the outlier-exclusion makes, analogous in spirit to the FC-explorer's
    kept-vs-non-kept matrix comparison but for our 1D per-region profiles
    rather than full ROI x ROI matrices (which this pipeline doesn't
    compute).

    df_class needs is_outlier + the same (subject,session,run,seed)
    identity columns as classify_runs_by_profile_distance's output.
    Raises ValueError if either group is empty.
    """
    kept_keys = {(str(r.subject), str(r.session), str(r.run), str(r.seed))
                for r in df_class.itertuples() if not bool(r.is_outlier)}
    out_keys = {(str(r.subject), str(r.session), str(r.run), str(r.seed))
               for r in df_class.itertuples() if bool(r.is_outlier)}
    profiles_str_keys = {tuple(str(x) for x in k): v for k, v in profiles.items()}

    kept_avg = compute_group_average_profile(profiles_str_keys, kept_keys)
    out_avg = compute_group_average_profile(profiles_str_keys, out_keys)
    if kept_avg is None:
        raise ValueError("kept group is empty — nothing to plot")
    if out_avg is None:
        raise ValueError("outlier group is empty — nothing to plot")

    if drop_empty_regions:
        kept_avg = kept_avg[kept_avg[value_col].notna()]
    order = kept_avg.sort_values(value_col)["region_name"].tolist()
    x = np.arange(len(order))
    kept_y = kept_avg.set_index("region_name").reindex(order)[value_col].values
    out_y = out_avg.set_index("region_name").reindex(order)[value_col].values

    fig, ax = plt.subplots(figsize=(max(10, 0.15 * len(order)), 6))
    ax.fill_between(x, kept_y, out_y, color="gray", alpha=0.15, zorder=1)
    ax.plot(x, kept_y, color="#2ca02c", linewidth=2, label="kept (average)", zorder=3)
    ax.plot(x, out_y, color="#d62728", linewidth=2, label="outlier (average)", zorder=3)
    ax.set_xticks(x[::max(1, len(x) // 40)])
    ax.set_xticklabels([order[i] for i in x[::max(1, len(x) // 40)]], rotation=90, fontsize=6)
    ax.set_ylabel(value_col)
    ax.set_xlabel("region (ordered by kept average)")
    ax.set_title(title)
    ax.legend(fontsize=8)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi)
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_interspecies_profile_comparison(profiles_by_bids, common_base_names, value_col="median_r",
                                          title="Inter-species retrosplenial profile comparison",
                                          output_path=None, dpi=300):
    """
    Overlay each BIDS's (already-averaged, e.g. kept-only) profile on the
    SAME x-axis of common regions (base_name, L_/R_ pooled together via
    mean), one line per BIDS — to visually compare how close the
    connectivity profiles are across species/datasets.

    profiles_by_bids: {bids_label: profile_df} — typically each
    compute_group_average_profile()'s output, one per BIDS.
    common_base_names: from find_common_base_regions() — lowercase,
    stripped base names present in every BIDS involved.

    Region order is by the MEAN value across all BIDS (readable sweep,
    same reasoning as plot_region_profile_curve).
    """
    pooled = {}
    for label, prof in profiles_by_bids.items():
        p = prof.copy()
        p["_base"] = p["base_name"].str.lower().str.strip()
        p = p[p["_base"].isin(common_base_names)]
        pooled[label] = p.groupby("_base")[value_col].mean()

    all_df = pd.DataFrame(pooled)
    order = all_df.mean(axis=1).sort_values().index.tolist()
    x = np.arange(len(order))

    fig, ax = plt.subplots(figsize=(max(10, 0.2 * len(order)), 6))
    colors = plt.cm.tab10.colors
    for i, (label, series) in enumerate(pooled.items()):
        y = series.reindex(order).values
        ax.plot(x, y, marker="o", markersize=3, linewidth=1.5,
               label=label, color=colors[i % len(colors)])
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=90, fontsize=6)
    ax.set_ylabel(value_col)
    ax.set_xlabel("common region (base name, ordered by cross-BIDS mean)")
    ax.set_title(title)
    ax.legend(fontsize=8, loc="best")
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi)
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_distance_to_reference_species(profiles_by_bids, common_base_names, reference_label,
                                        value_col="median_r", metric="correlation",
                                        title=None, output_path=None, dpi=300):
    """
    Distance of every OTHER BIDS's pooled profile to reference_label's
    pooled profile (e.g. Human), restricted to common_base_names.
    Reuses qc._distance (same metric conventions as the rest of the
    pipeline). Returns (fig, df) where df has one row per non-reference
    BIDS: bids_label, distance, pearson_r, n_regions — sorted by distance
    ascending (closest to the reference first).
    """
    pooled = {}
    for label, prof in profiles_by_bids.items():
        p = prof.copy()
        p["_base"] = p["base_name"].str.lower().str.strip()
        p = p[p["_base"].isin(common_base_names)]
        pooled[label] = p.groupby("_base")[value_col].mean()

    if reference_label not in pooled:
        raise ValueError(f"reference_label {reference_label!r} not found among "
                         f"{list(pooled)}")
    ref = pooled[reference_label]

    rows = []
    for label, series in pooled.items():
        if label == reference_label:
            continue
        merged = pd.DataFrame({"ref": ref, "this": series}).dropna()
        if merged.empty:
            print(f"  [WARN] {label}: no overlapping regions with {reference_label}, skipped")
            continue
        dist, r = qc._distance(merged["this"].values, merged["ref"].values, metric=metric)
        rows.append({"bids_label": label, "distance": dist, "pearson_r": r,
                     "n_regions": len(merged)})
    df = pd.DataFrame(rows).sort_values("distance").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(max(6, 0.8 * len(df)), 5))
    ax.bar(df["bids_label"], df["distance"], color="#4472c4")
    ax.set_ylabel(f"distance to {reference_label} ({metric})")
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["bids_label"], rotation=45, ha="right")
    ax.set_title(title or f"Distance to {reference_label} reference ({value_col})")
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi)
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig, df


# =============================================================================
# 7) Plot: population curve + individual run overlay
# =============================================================================

def plot_region_profile_curve(ref_profile, run_profiles=None, value_col="median_r",
                               title="DMN region profile", output_path=None, dpi=300,
                               is_outlier_map=None, smooth_window=5, ci_pct=95,
                               drop_empty_regions=True):
    """
    x = regions, ordered by the REFERENCE curve's value (so the shape is
    readable — a monotonic-ish sweep rather than a random-looking atlas
    order). y = value_col.

      - Reference curve: bold black, smoothed with a centered rolling
        mean over `smooth_window` neighboring (in sorted order) regions.
        This is a moving average across the ranked sequence, not a
        continuous-variable smoother — the x-axis is categorical
        (discrete regions), so "smooth" here means "less jagged from
        region to region", not interpolation.
      - Confidence band: `ci_pct`% empirical CI (percentile-based) of
        value_col ACROSS every run in run_profiles, per region, shaded
        gray behind the reference curve — the spread of the population at
        each region, not the reference's own uncertainty.
      - Individual runs: thin lines, blue if kept / red if flagged
        (is_outlier_map: {(subject,session,run,seed): bool}, e.g. from
        classify_runs_by_profile_distance's is_outlier column).
      - drop_empty_regions=True (default): regions where the REFERENCE
        has no voxels (NaN value_col — atlas_level didn't cover that
        region for this map) are excluded from the x-axis entirely,
        instead of showing up as an empty/flat gap.
    """
    if drop_empty_regions:
        n_before = len(ref_profile)
        ref_profile = ref_profile[ref_profile[value_col].notna()]
        n_after = len(ref_profile)
        if n_before != n_after:
            print(f"  [INFO] plot_region_profile_curve: dropped {n_before - n_after} "
                  f"region(s) with no voxels (NaN {value_col}) from the plot")

    ref_sorted = ref_profile.sort_values(value_col).reset_index(drop=True)
    order = ref_sorted["region_name"].tolist()
    x = np.arange(len(order))
    ref_vals = ref_sorted[value_col].values

    def _smooth(y, window):
        if window is None or window <= 1:
            return y
        return pd.Series(y).rolling(window, center=True, min_periods=1).mean().values

    fig, ax = plt.subplots(figsize=(max(10, 0.15 * len(order)), 6))

    if run_profiles:
        stacked = []
        for key, prof in run_profiles.items():
            vals = prof.set_index("region_name").reindex(order)[value_col].values
            stacked.append(vals)
        stacked = np.array(stacked, dtype=float)

        if ci_pct is not None and stacked.shape[0] > 1:
            lo_pct, hi_pct = (100 - ci_pct) / 2, 100 - (100 - ci_pct) / 2
            ci_lo = np.nanpercentile(stacked, lo_pct, axis=0)
            ci_hi = np.nanpercentile(stacked, hi_pct, axis=0)
            ax.fill_between(x, _smooth(ci_lo, smooth_window), _smooth(ci_hi, smooth_window),
                            color="gray", alpha=0.25, zorder=1,
                            label=f"{ci_pct:.0f}% CI across runs")

        for key, vals in zip(run_profiles.keys(), stacked):
            is_out = bool(is_outlier_map.get(key, False)) if is_outlier_map else False
            ax.plot(x, vals, color="#d62728" if is_out else "#4472c4",
                    alpha=0.6 if is_out else 0.12, linewidth=1.3 if is_out else 0.6, zorder=2)

    ref_smoothed = _smooth(ref_vals, smooth_window)
    if smooth_window is not None and smooth_window > 1:
        print(f"  [INFO] {value_col} dynamic range — raw: [{np.nanmin(ref_vals):.4f}, "
              f"{np.nanmax(ref_vals):.4f}] (span={np.nanmax(ref_vals) - np.nanmin(ref_vals):.4f}), "
              f"smoothed (window={smooth_window}): [{np.nanmin(ref_smoothed):.4f}, "
              f"{np.nanmax(ref_smoothed):.4f}] (span={np.nanmax(ref_smoothed) - np.nanmin(ref_smoothed):.4f})"
              f" — if the smoothed span is much smaller, the window is flattening real structure; "
              f"try smooth_window=1 to see the raw curve instead")
        ax.plot(x, ref_vals, color="black", linewidth=1.0, linestyle="--", alpha=0.5,
                label="reference (raw, unsmoothed)", zorder=3)

    ax.plot(x, ref_smoothed, color="black", linewidth=2.5,
            label="reference (group map, smoothed)" if (smooth_window and smooth_window > 1)
                  else "reference (group map)", zorder=3)

    ax.set_xticks(x[::max(1, len(x) // 40)])
    ax.set_xticklabels([order[i] for i in x[::max(1, len(x) // 40)]],
                       rotation=90, fontsize=6)
    ax.set_ylabel(value_col)
    ax.set_xlabel("region (ordered by reference value)")
    ax.set_title(title)
    ax.legend(fontsize=8, loc="upper left")
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=dpi)
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


# =============================================================================
# 8) Optional alternative index: Dice overlap between run and reference
# =============================================================================

def compute_dice_overlap(run_map_path, ref_map_path, mask_path=None, threshold=0.0):
    """
    Fresh, independently-defined Dice metric — NOT assumed to be the same
    as any pre-existing 'contiguity_dice' column you may already have;
    compare the two before relying on this if that one still exists.

    Definition: binarize both maps at `threshold` (default: r > 0, i.e.
    "positively connected voxel"), restricted to `mask_path` if given
    (else wherever the reference map is non-zero), then
        Dice = 2*|A ∩ B| / (|A| + |B|)
    over those binary masks. 1 = identical positive-connectivity
    footprint, 0 = no overlap at all.
    """
    run_data = qc._squeeze_trailing_singleton(nib.load(run_map_path).get_fdata())
    ref_img = nib.load(ref_map_path)
    ref_data = qc._squeeze_trailing_singleton(ref_img.get_fdata())

    if mask_path:
        mask = qc._squeeze_trailing_singleton(nib.load(mask_path).get_fdata()) != 0
    else:
        mask = ref_data != 0

    if run_data.shape != ref_data.shape:
        raise ValueError(f"Grid mismatch: {run_data.shape} vs {ref_data.shape}")

    a = (run_data > threshold) & mask
    b = (ref_data > threshold) & mask
    denom = a.sum() + b.sum()
    if denom == 0:
        return np.nan
    return float(2 * (a & b).sum() / denom)
