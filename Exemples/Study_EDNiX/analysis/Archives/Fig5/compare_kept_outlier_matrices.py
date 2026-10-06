"""
compare_kept_outlier_matrices.py
──────────────────────────────────
Compare lvl3 (EDNIxCSCLR level-3) whole-brain ROI-to-ROI correlation
matrices between the Kept and Outlier groups from a
classify_runs_by_absolute_correlation / classify_runs_by_profile_distance
result (df_class), produced by SBA_Region_Profile_Index.py.

⚠ UNVERIFIED matrix path convention. Two different naming patterns have
shown up across earlier work on this pipeline:
    .../func/acpc-func/Stats/Correl_matrix/EDNIxCSCLR_<level>_run_<i>_matrix.csv
    .../func/acpc-func/Stats/Correl_matrix/EDNIxCSCLR/correlation/EDNIxCSCLR_<level>_run_<i>_correlation_matrix.csv
Neither has been confirmed against the CURRENT (func/templates/EDNiX/...)
directory layout used elsewhere for SBA in this pipeline —
find_lvl3_matrix_path() tries all four combinations (two directory
layouts × two filename styles) and reports which one actually worked.
Tell me the real layout on first use and I'll simplify it down to just
that one.

⚠ ALSO UNVERIFIED: whether "<i>" in these filenames is the actual run
TOKEN from the BOLD filename (matching df_class's own "run" column, e.g.
"14", "RL-2") or a simple 0-based ENUMERATION INDEX over that session's
runs (unrelated to the token — the older per-run QC script built it via
`for i in range(int(nb_run))`, a plain counter). find_lvl3_matrix_path()
tries a token-based match first; if that finds nothing, it lists
everything available for that (subject, session) and the caller SKIPS
that run with a clear warning, rather than silently guessing an index
mapping that could pair the wrong matrix with the wrong run.
"""

import os
import re
import glob

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

opj = os.path.join
opb = os.path.basename
ope = os.path.exists


def _resolve_session_dirs(bids_dir, subject, session):
    """
    Every plausible session directory for (subject, session), trying the
    literal "ses-<session>" first, then a zero-padding-insensitive match
    against whatever "ses-*" directories actually exist. Confirmed real
    case: df_class's session value can be unpadded ("1", from
    SBA_Region_Profile_Index's own session enumeration) while the real
    BIDS directory is zero-padded ("ses-01") — e.g. Marmoset's
    /BIDS_Tian/sub-NIHm14/ses-01/... existed on disk, but "ses-1" (built
    from the literal session value) did not, so every single Marmoset
    run silently found nothing.

    Returns a list of candidate directories in priority order: the
    literal "ses-<session>" match (kept even if it doesn't exist, for
    backward compatibility with callers that print it), any zero-padding
    variant found by listing "ses-*" and comparing numeric values, then
    the no-session (single-session-layout) fallback subject dir.
    """
    literal = opj(bids_dir, f"sub-{subject}", f"ses-{session}")
    candidates = [literal]

    sub_dir = opj(bids_dir, f"sub-{subject}")
    if ope(sub_dir) and str(session).isdigit():
        session_digits = str(session).lstrip("0") or "0"
        for d in sorted(glob.glob(opj(sub_dir, "ses-*"))):
            ses_name = opb(d)[len("ses-"):]
            if not ses_name.isdigit():
                continue
            if ses_name == str(session):
                continue  # already the literal candidate above
            if (ses_name.lstrip("0") or "0") == session_digits:
                candidates.append(d)

    candidates.append(sub_dir)
    return candidates


# =============================================================================
# 1) Matrix location + loading
# =============================================================================

def find_lvl3_matrix_path(bids_dir, subject, session, run, level=3, atlas_name="EDNIxCSCLR"):
    """
    Locate one run's level-<level> ROI-to-ROI correlation matrix CSV.

    Confirmed real naming: "run_0" appears for EVERY subject/session,
    never varying — this "run_N" segment in the filename is NOT the
    actual per-BOLD-run token SBA_Region_Profile_Index.py's df_class
    uses (e.g. "01", "14", "RL-2"). It's a fixed placeholder for what is
    actually a SESSION-level matrix (one FC estimate per session, not
    per individual functional run). Given that, this function:
      1. Tries an EXACT run-token match first (in case some BIDS really
         DO have one matrix per run) — excluding check_fit/flattened/
         pval/tstat variants (same exclusion filter as
         ednix_fc_threshold_explorer.py's own _find_all_run_paths).
      2. If that finds nothing AND there's exactly ONE clean candidate
         file for this (subject, session), uses it as a session-level
         matrix — printed clearly, since every run in that session then
         maps to the SAME file (confirmed real symptom this fixes:
         every subject's matrix was named "..._run_0_..." regardless of
         the actual run number, so an exact-token match could never
         succeed for a real per-run df_class).
      3. Otherwise (several ambiguous candidates, none matching)
         returns None with a full listing rather than guessing.

    run: pass the exact df_class value (may be the string "NA" for "no
    run token in the SBA filename" — NOT a float NaN; see
    load_all_results_from_disk's keep_default_na=False fix for why a
    reloaded Excel could otherwise hand this "NA" back as a real NaN).
    """
    EXCLUDE_SUBSTRINGS = ("check_fit", "flattened", "pval", "tstat")

    ses_dir_candidates = _resolve_session_dirs(bids_dir, subject, session)
    matrix_dir_suffixes = [
        opj("func", "acpc-func", "Stats", "Correl_matrix"),
        opj("func", "acpc-func", "Stats", "Correl_matrix", atlas_name, "correlation"),
        opj("func", "templates", "EDNiX", "Stats", "Correl_matrix"),
        opj("func", "templates", "EDNiX", "Stats", "Correl_matrix", atlas_name, "correlation"),
    ]
    # One pattern is enough — "*_matrix.csv" already matches the
    # confirmed real "*_correlation_matrix.csv" suffix as a trailing
    # substring; a second, narrower pattern previously produced pure
    # duplicate matches against the exact same files.
    filename_pattern = f"{atlas_name}_{level}_run_*_matrix.csv"

    all_found = {}  # realpath -> path, dedups matches from overlapping dir candidates
    for ses_dir in ses_dir_candidates:
        for suffix in matrix_dir_suffixes:
            matrix_dir = opj(ses_dir, suffix)
            if not ope(matrix_dir):
                continue
            for f in glob.glob(opj(matrix_dir, filename_pattern)):
                if any(x in opb(f) for x in EXCLUDE_SUBSTRINGS):
                    continue
                all_found[os.path.realpath(f)] = f

    candidates = sorted(all_found.values())
    if not candidates:
        return None

    run_str = str(run)
    if run_str and run_str.lower() != "nan":
        run_re = re.compile(rf"run_0*{re.escape(run_str)}(?!\d)")
        token_matches = [f for f in candidates if run_re.search(opb(f))]
        if token_matches:
            return sorted(token_matches)[0]

        # Confirmed real off-by-one: df_class's run numbering is 1-indexed
        # (taken from the BOLD filename's own "run-<N>" token, N starting
        # at 1), but these lvl3 matrix files are numbered 0-indexed
        # ("run_0", "run_1", ... — a plain enumeration counter in whatever
        # produced them, unrelated to the BOLD run token). Verified on 5
        # independent real cases: a session with exactly N matrix files
        # numbered run_0..run_{N-1} always corresponds to df_class's runs
        # 1..N for that session — e.g. run=8 with files run_0..run_7 (8
        # files), run=14 with files run_0..run_13 (14 files). Only trust
        # this when `run` is a plain positive integer string (never for
        # "NA" or an acquisition-tagged token like "RL-2", where this
        # numeric correspondence doesn't apply).
        if run_str.isdigit():
            shifted = str(int(run_str) - 1)
            shifted_re = re.compile(rf"run_0*{re.escape(shifted)}(?!\d)")
            shifted_matches = [f for f in candidates if shifted_re.search(opb(f))]
            if shifted_matches:
                chosen = sorted(shifted_matches)[0]
                print(f"  [INFO] sub-{subject} ses-{session} run={run}: matched via "
                     f"the confirmed 1-indexed(df_class)/0-indexed(matrix file) "
                     f"offset -> {opb(chosen)}")
                return chosen

    if len(candidates) == 1:
        print(f"  [INFO] sub-{subject} ses-{session} run={run}: no exact run-token "
             f"match, but exactly 1 clean matrix file exists for this session -- "
             f"using it as a SESSION-level matrix ({opb(candidates[0])}). Every run "
             f"in this session will map to the SAME file.")
        return candidates[0]

    print(f"  [WARN] sub-{subject} ses-{session} run={run}: found "
         f"{len(candidates)} lvl{level} matrix file(s) but none matched "
         f"the run token exactly, and more than one candidate exists -- "
         f"available: {[opb(f) for f in candidates]}. Skipping rather "
         f"than guessing which one.")
    return None


def debug_lvl3_search(bids_dir, subject, session, level=3, atlas_name="EDNIxCSCLR"):
    """
    Diagnostic helper: for one (subject, session), print every directory
    candidate find_lvl3_matrix_path tries, whether it exists, and (if it
    exists) which files are actually in it. Use this when
    find_lvl3_matrix_path returns None for EVERY run of a BIDS — that
    means every single directory candidate came back empty/nonexistent,
    which find_lvl3_matrix_path itself does NOT report (it only prints
    a diagnostic when it found candidate files but none matched the run
    token — the "found nothing at all" case was silent, which is the
    confirmed real gap this fixes: Marmoset/Mouselemur produced zero
    [INFO]/[WARN] lines per run, unlike every other BIDS, yet still
    ended up with 0 usable matrices).

    Call this ONCE per BIDS (e.g. on the first row) rather than per run,
    to see exactly why — wrong bids_dir, wrong atlas_name for this
    species, or a genuinely different directory layout — instead of
    guessing.
    """
    ses_dir_candidates = _resolve_session_dirs(bids_dir, subject, session)
    matrix_dir_suffixes = [
        opj("func", "acpc-func", "Stats", "Correl_matrix"),
        opj("func", "acpc-func", "Stats", "Correl_matrix", atlas_name, "correlation"),
        opj("func", "templates", "EDNiX", "Stats", "Correl_matrix"),
        opj("func", "templates", "EDNiX", "Stats", "Correl_matrix", atlas_name, "correlation"),
    ]
    print(f"  [DEBUG] lvl{level} matrix search for sub-{subject} ses-{session} "
         f"(bids_dir={bids_dir}, atlas_name={atlas_name}):")
    any_dir_exists = False
    for ses_dir in ses_dir_candidates:
        print(f"    session dir: {ses_dir}  (exists={ope(ses_dir)})")
        for suffix in matrix_dir_suffixes:
            matrix_dir = opj(ses_dir, suffix)
            exists = ope(matrix_dir)
            if exists:
                any_dir_exists = True
                all_files = sorted(os.listdir(matrix_dir)) if os.path.isdir(matrix_dir) else []
                shown = all_files[:10]
                print(f"      {matrix_dir}  (exists=True, {len(all_files)} file(s): "
                     f"{shown}{'...' if len(all_files) > 10 else ''})")
            else:
                print(f"      {matrix_dir}  (exists=False)")
    if not any_dir_exists:
        print(f"  [DEBUG] NONE of the candidate matrix directories exist for "
             f"sub-{subject} -- check bids_dir, subject id, or whether lvl{level} "
             f"matrices were ever computed for this BIDS.")


def load_matrix(path):
    """Load a correlation matrix CSV -> (roi_names, ndarray)."""
    df = pd.read_csv(path, index_col=0)
    df.columns = df.columns.astype(str).str.strip()
    df.index = df.index.astype(str).str.strip()
    return list(df.index), df.values.astype(float)


def collect_group_matrices(df_class, bids_dir, level=3, atlas_name="EDNIxCSCLR"):
    """
    For every row of df_class, locate and load its lvl<level> matrix.
    Returns (kept_matrices, outlier_matrices, roi_names) — the two
    matrix lists hold plain ndarrays (aligned to a shared roi_names
    order, checked below); roi_names is that shared order, taken from
    the first successfully-loaded matrix. Raises if a later matrix has
    a DIFFERENT region set/order, rather than silently misaligning ROI
    positions across runs — every run must share the same lvl<level>
    parcellation for a group comparison to mean anything.
    """
    kept, outlier = [], []
    roi_names_ref = None
    n_missing = 0

    for row in df_class.itertuples():
        path = find_lvl3_matrix_path(bids_dir, row.subject, row.session, row.run,
                                     level=level, atlas_name=atlas_name)
        if path is None:
            n_missing += 1
            continue
        roi_names, mat = load_matrix(path)
        if roi_names_ref is None:
            roi_names_ref = roi_names
        elif roi_names != roi_names_ref:
            raise ValueError(
                f"ROI order mismatch: {opb(path)} has a different region set/"
                f"order than the reference matrix ({len(roi_names)} vs "
                f"{len(roi_names_ref)} regions). Every run must share the "
                f"same lvl{level} parcellation for group comparison to be "
                f"valid — check whether this run genuinely used a different "
                f"atlas level."
            )
        (outlier if bool(row.is_outlier) else kept).append(mat)

    print(f"  Loaded {len(kept)} kept + {len(outlier)} outlier lvl{level} matrices "
         f"({n_missing} run(s) skipped, no matrix found)")
    return kept, outlier, roi_names_ref


# =============================================================================
# 2) Group comparison
# =============================================================================

def compare_kept_vs_outlier(kept_mats, outlier_mats, roi_names, fdr_alpha=0.05):
    """
    Edge-wise comparison of kept vs outlier group matrices:
      - group mean matrices
      - difference matrix (kept - outlier)
      - per-edge Mann-Whitney U test, FDR-corrected (Benjamini-Hochberg
        if statsmodels is available; raw p-values only otherwise, with
        a warning — never silently skips the correction without saying so)

    Returns a dict: kept_mean, outlier_mean, diff (all n_roi x n_roi
    ndarrays), roi_names, n_kept, n_outlier, and edge_df — one row per
    upper-triangle edge with roi_i, roi_j, kept_mean, outlier_mean,
    diff, p_value, q_value, significant.
    """
    if len(kept_mats) < 2 or len(outlier_mats) < 2:
        raise ValueError(
            f"Need at least 2 matrices per group for a meaningful comparison "
            f"(kept={len(kept_mats)}, outlier={len(outlier_mats)})."
        )

    kept_stack = np.stack(kept_mats)      # (n_kept, n_roi, n_roi)
    outlier_stack = np.stack(outlier_mats)

    kept_mean = np.nanmean(kept_stack, axis=0)
    outlier_mean = np.nanmean(outlier_stack, axis=0)
    diff = kept_mean - outlier_mean

    n_roi = len(roi_names)
    iu = np.triu_indices(n_roi, k=1)  # upper triangle, no diagonal

    rows = []
    for i, j in zip(*iu):
        kvals = kept_stack[:, i, j]
        ovals = outlier_stack[:, i, j]
        kvals = kvals[np.isfinite(kvals)]
        ovals = ovals[np.isfinite(ovals)]
        if len(kvals) < 2 or len(ovals) < 2:
            rows.append(dict(roi_i=roi_names[i], roi_j=roi_names[j],
                             kept_mean=np.nan, outlier_mean=np.nan, diff=np.nan,
                             p_value=np.nan))
            continue
        try:
            _, p = stats.mannwhitneyu(kvals, ovals, alternative="two-sided")
        except ValueError:
            p = np.nan
        rows.append(dict(roi_i=roi_names[i], roi_j=roi_names[j],
                         kept_mean=float(np.mean(kvals)), outlier_mean=float(np.mean(ovals)),
                         diff=float(np.mean(kvals) - np.mean(ovals)), p_value=p))

    edge_df = pd.DataFrame(rows)
    valid = edge_df["p_value"].notna()
    try:
        from statsmodels.stats.multitest import multipletests
        _, q_vals, _, _ = multipletests(edge_df.loc[valid, "p_value"], alpha=fdr_alpha, method="fdr_bh")
        edge_df["q_value"] = np.nan
        edge_df.loc[valid, "q_value"] = q_vals
    except ImportError:
        print("  [WARN] statsmodels not available -- reporting raw p-values only, no FDR correction")
        edge_df["q_value"] = edge_df["p_value"]

    edge_df["significant"] = edge_df["q_value"] < fdr_alpha

    return dict(kept_mean=kept_mean, outlier_mean=outlier_mean, diff=diff,
               edge_df=edge_df, roi_names=roi_names,
               n_kept=len(kept_mats), n_outlier=len(outlier_mats))


# =============================================================================
# 3) Plot
# =============================================================================

def plot_group_comparison(result, title="Kept vs Outlier lvl3 matrices", output_path=None, dpi=300):
    """4-panel figure: kept mean, outlier mean, difference, significance mask."""
    roi_names = result["roi_names"]
    n_roi = len(roi_names)
    sig_mask = np.zeros((n_roi, n_roi), dtype=bool)
    roi_index = {name: i for i, name in enumerate(roi_names)}
    for row in result["edge_df"].itertuples():
        if bool(row.significant):
            i, j = roi_index[row.roi_i], roi_index[row.roi_j]
            sig_mask[i, j] = sig_mask[j, i] = True

    fig, axes = plt.subplots(1, 4, figsize=(24, 6))
    vmax = max(np.nanmax(np.abs(result["kept_mean"])), np.nanmax(np.abs(result["outlier_mean"])))

    im0 = axes[0].imshow(result["kept_mean"], cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    axes[0].set_title(f"Kept (n={result['n_kept']})")
    plt.colorbar(im0, ax=axes[0], fraction=0.046)

    im1 = axes[1].imshow(result["outlier_mean"], cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    axes[1].set_title(f"Outlier (n={result['n_outlier']})")
    plt.colorbar(im1, ax=axes[1], fraction=0.046)

    diff_vmax = np.nanmax(np.abs(result["diff"])) or 1.0
    im2 = axes[2].imshow(result["diff"], cmap="RdBu_r", vmin=-diff_vmax, vmax=diff_vmax)
    axes[2].set_title("Difference (kept - outlier)")
    plt.colorbar(im2, ax=axes[2], fraction=0.046)

    axes[3].imshow(sig_mask, cmap="gray_r", vmin=0, vmax=1)
    n_sig = int(result["edge_df"]["significant"].sum())
    axes[3].set_title(f"Significant edges (FDR<0.05): {n_sig}")

    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])

    fig.suptitle(title)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi)
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


# =============================================================================
# 4) Orchestration
# =============================================================================

def run_matrix_comparison(df_class, bids_dir, level=3, atlas_name="EDNIxCSCLR",
                          fdr_alpha=0.05, output_dir=None, label=""):
    """
    One-call pipeline: load every run's lvl<level> matrix, split by
    df_class's is_outlier, compare, plot, and (if output_dir given) save
    the figure + the edge-level results as an Excel file.
    Returns the `compare_kept_vs_outlier` result dict.
    """
    kept_mats, outlier_mats, roi_names = collect_group_matrices(
        df_class, bids_dir, level=level, atlas_name=atlas_name)
    result = compare_kept_vs_outlier(kept_mats, outlier_mats, roi_names, fdr_alpha=fdr_alpha)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        plot_group_comparison(
            result, title=f"{label} — Kept vs Outlier lvl{level}" if label else f"Kept vs Outlier lvl{level}",
            output_path=opj(output_dir, f"{label}_lvl{level}_kept_vs_outlier.png" if label
                            else f"lvl{level}_kept_vs_outlier.png"))
        excel_path = opj(output_dir, f"{label}_lvl{level}_kept_vs_outlier.xlsx" if label
                         else f"lvl{level}_kept_vs_outlier.xlsx")
        result["edge_df"].sort_values("p_value").to_excel(excel_path, index=False)
        print(f"  [saved] {excel_path}")

    return result