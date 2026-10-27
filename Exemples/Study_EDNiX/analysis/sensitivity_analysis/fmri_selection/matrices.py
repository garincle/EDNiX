"""
matrices.py
================
The Figure 5 analysis run on whole-brain correlation matrices (EDNiX step 10,
modalities/fMRI/_10_Correl_matrix.py) instead of retrosplenial seed maps --
same selection chain, same consensus, same figures -- and the comparison of
the runs each analysis keeps.

What a "profile" is here
------------------------
Seed analysis: one value per target region (median r of the seed map in it).
Matrix analysis: one value per EDGE (Pearson r between two regional mean time
series). Hemispheres are handled the way the seed analysis handles them (it
averages L and R of a region through `base_name`):
    "a <> b |intra"     L_a-L_b and R_a-R_b          (averaged by profile_stats)
    "a <> b |inter"     L_a-R_b and R_a-L_b          (averaged)
    "a |homotopic"      L_a-R_a
region_name keeps the hemispheric edge (e.g. "L_a <> L_b"), so the k-core
references, which merge on region_name, see every edge separately -- again as
for the seed profiles. Nodes default to the seed run's common regions
(02_selection/common_regions.xlsx), so both analyses cover the same anatomy.

Where the matrices are
----------------------
<bids>/sub-<s>/[ses-<x>/]func/acpc-func/Stats/Correl_matrix/<atlas>/<type>/
    <atlas>_<level index>_run_<i>_<type>_matrix.csv
<level index> is the 0-based volume of the 4D atlas (atlas_data[..., idx] in
_10_Correl_matrix), i.e. ATLAS_LEVEL - 1 in the 1-based convention of
SBA_Region_Profile_Index.compute_region_profile.
<i> is the position of the run in EDNiX's run list (sorted glob of endfmri,
4D files with >= ntimepoint_treshold volumes), NOT the BIDS run number. It is
resolved from <atlas>_<level>_run_<i>_filtered.json, written by the same step
next to the csv, whose "Sources"[1] is the residual time series actually used;
if that json is missing, the run list is rebuilt from the raw BIDS folder.
Every resolution is written to matrix_run_mapping.xlsx -- check it once.
"""

import os
import re
import glob
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as st
from pathlib import Path

import Exemples.Study_EDNiX.analysis.fmri_fingerprints.fingerprints as ff

opj = os.path.join
opb = os.path.basename
EDGE_LINK = ff.EDGE_LINK
RESIDUAL_RE = re.compile(r"_space-acpc-func_desc-fMRI_residual\.nii(\.gz)?$")


# ═══════════════════════════════════════════════════════════════════════════
# Files
# ═══════════════════════════════════════════════════════════════════════════

def session_dir(bids_dir, subj, ses):
    d = opj(bids_dir, f"sub-{subj}", f"ses-{ses}")
    return d if os.path.isdir(d) else opj(bids_dir, f"sub-{subj}")


def matrix_folder(bids_dir, subj, ses, atlas_name, correl_type):
    return opj(session_dir(bids_dir, subj, ses), "func", "acpc-func", "Stats",
               "Correl_matrix", atlas_name, correl_type)


def find_matrix_files(bids_dir, subj, ses, atlas_name, level_idx, correl_type):
    """[(run index i, csv path)] sorted by i."""
    d = matrix_folder(bids_dir, subj, ses, atlas_name, correl_type)
    pat = re.compile(rf"^{re.escape(atlas_name)}_{int(level_idx)}_run_(\d+)_"
                     rf"{re.escape(correl_type)}_matrix\.csv$")
    hits = []
    for f in glob.glob(opj(d, f"{atlas_name}_{int(level_idx)}_run_*_matrix.csv")):
        m = pat.match(opb(f))
        if m:
            hits.append((int(m.group(1)), f))
    return sorted(hits)


def residual_from_json(csv_path, atlas_name, level_idx, i):
    """Residual time series used for matrix i, from the *_filtered.json that
    _10_Correl_matrix writes in <atlas>/ (one level above <type>/)."""
    j = opj(os.path.dirname(os.path.dirname(csv_path)),
            f"{atlas_name}_{int(level_idx)}_run_{int(i)}_filtered.json")
    if not os.path.isfile(j):
        return None
    try:
        with open(j) as f:
            data = json.load(f)
    except (OSError, ValueError):
        return None
    if isinstance(data, list):        # the step dumps a 1-tuple -> JSON list
        data = data[0] if data else {}
    src = data.get("Sources") if isinstance(data, dict) else None
    if isinstance(src, list) and len(src) >= 2 and isinstance(src[1], str):
        return src[1]
    return None


def sorted_bold_list(bids_dir, subj, ses, endfmri, nt_min):
    """EDNiX's run list: sorted(glob(func/endfmri)), 4D with >= nt_min volumes."""
    import nibabel as nib
    func = opj(session_dir(bids_dir, subj, ses), "func")
    keep = []
    for f in sorted(glob.glob(opj(func, endfmri))):
        try:
            shape = nib.load(f).shape
        except Exception:
            continue
        if len(shape) == 4 and shape[3] >= nt_min:
            keep.append(f)
    return keep


def bold_root(path):
    name = opb(path)
    name = RESIDUAL_RE.sub("", name)
    for suf in (".nii.gz", ".nii"):
        if name.endswith(suf):
            name = name[: -len(suf)]
    return name


# ═══════════════════════════════════════════════════════════════════════════
# Matrix -> edge profile
# ═══════════════════════════════════════════════════════════════════════════

def read_matrix(csv_path):
    df = pd.read_csv(csv_path, index_col=0)
    names = [str(n) for n in df.index]
    if list(df.columns.astype(str)) != names:
        raise ValueError(f"{csv_path}: row and column labels differ")
    return names, df.values.astype(float)


def node_info(name, lut):
    """(base_name lower-case, hemisphere) of a matrix row label."""
    if name in lut:
        base, hemi = lut[name]
    else:
        m = re.match(r"^([LlRr])_(.*)$", name)
        base, hemi = (m.group(2), "left" if m.group(1) in "Ll" else "right") if m \
            else (name, "bilateral")
    return str(base).lower().strip(), hemi


def edge_names(a, b, lut):
    """(region_name, base_name, hemisphere code) of the edge a-b."""
    (ba, ha), (bb, hb) = node_info(a, lut), node_info(b, lut)
    lat = ha in ("left", "right") and hb in ("left", "right")
    if lat and ha != hb:
        cls = "homotopic" if ba == bb else "inter"
    else:
        cls = "intra"
    lo, hi = sorted([ba, bb])
    base = f"{lo} |homotopic" if cls == "homotopic" else f"{lo}{EDGE_LINK}{hi} |{cls}"
    n1, n2 = sorted([a, b])
    code = "".join(sorted(h[0].upper() for h in (ha, hb)))
    return f"{n1}{EDGE_LINK}{n2}", base, code


def edge_table(names, lut, keep_nodes=None):
    """Upper-triangle edges among `names` whose two nodes are in keep_nodes
    (base names, lower-case; None = all)."""
    rows = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            if keep_nodes is not None and (node_info(names[i], lut)[0] not in keep_nodes
                                           or node_info(names[j], lut)[0] not in keep_nodes):
                continue
            rn, bn, code = edge_names(names[i], names[j], lut)
            rows.append((rn, bn, code))
    return rows


# ═══════════════════════════════════════════════════════════════════════════
# The stand-in for SBA_Region_Profile_Index
# ═══════════════════════════════════════════════════════════════════════════

class MatrixProfiles:
    """Pass this object instead of SBA_Region_Profile_Index to
    exploration_pipeline.run(). Only compute_all_run_profiles is replaced; every other
    attribute (kcore reference, group median, label parser, classification
    reader) is the real sri one, so the selection is literally the same code.

    Run the pipeline with RUN_ALL_UNITS_MAPS=False, RUN_GROUP_MAPS=False,
    COPY_JPG=False, RENDER_SURFACE=False, RUN_DISTANCE_TO_REF=False: there is
    no AFNI map of a matrix."""

    def __init__(self, sri, bids_configs, atlas_name="EDNIxCSCLR", correl_type="correlation",
                 level_index=None, node_regions=None, value_col="median_r", nt_min=100,
                 default_endfmri="*_task-rest*_bold.nii.gz"):
        self._sri = sri
        self._cfg_by_bids = {c["bids_dir"]: c for c in bids_configs}
        self.atlas_name = atlas_name
        self.correl_type = correl_type
        self.level_index = level_index
        self.node_regions = (None if node_regions is None
                             else {str(n).lower().strip() for n in node_regions})
        self.value_col = value_col
        self.nt_min = nt_min
        self.default_endfmri = default_endfmri
        self.mapping_rows = []

    def __getattr__(self, name):
        return getattr(self._sri, name)

    # -- run index -> BIDS run token ---------------------------------------
    def _resolve_runs(self, cb, subj, ses, hits, level_idx):
        out = {}
        bold_list = None
        for i, csv in hits:
            res = residual_from_json(csv, self.atlas_name, level_idx, i)
            method = "json"
            if res is None:
                if bold_list is None:
                    bold_list = sorted_bold_list(cb["bids_dir"], subj, ses,
                                                 cb.get("endfmri", self.default_endfmri),
                                                 self.nt_min)
                res = bold_list[i] if i < len(bold_list) else None
                method = "sorted_bold"
            root = bold_root(res) if res else None
            token = self._sri._extract_run_id(root) if root else None
            ok_sub = root is not None and f"sub-{subj}" in root
            out[i] = (token, root, method)
            self.mapping_rows.append(dict(
                label=cb["label"], subject=subj, session=ses, matrix_index=i,
                run=token, bold_root=root, task=ff.task_of(root) if root else None,
                method=method, subject_in_name=ok_sub, csv=csv))
        return out

    # -- the replacement ---------------------------------------------------
    def compute_all_run_profiles(self, species, bids_dir, label_df, classification_path,
                                 category="All", seeds=None, atlas_level=None, **_):
        cb = self._cfg_by_bids.get(bids_dir, dict(label=Path(bids_dir).name,
                                                  bids_dir=bids_dir))
        if self.level_index is not None:
            level_idx = int(self.level_index)
        elif atlas_level is not None:
            level_idx = int(atlas_level) - 1
        else:
            raise ValueError("ATLAS_LEVEL is None: set MATRIX_LEVEL_INDEX explicitly")
        seed_name = (seeds or ["matrix"])[0]

        cat = self._sri.fig5.load_classification_subjects(
            classification_path, species=species, bids_dir=Path(bids_dir).name,
            categories=[category])
        if not isinstance(cat, dict) or not cat.get(category):
            print(f"  [WARN] 0 subjects from load_classification_subjects for "
                  f"{species}/{category} (bids_dir={Path(bids_dir).name})")
            return {}, pd.DataFrame()
        subjects = sorted({str(s) for _, _, s, _ in cat[category]})
        subs_sessions = []
        for subj in subjects:
            ses_dirs = sorted(glob.glob(opj(bids_dir, f"sub-{subj}", "ses-*")))
            subs_sessions += ([(subj, opb(d).replace("ses-", "")) for d in ses_dirs]
                              or [(subj, "1")])

        lut = {}
        if label_df is not None and not label_df.empty:
            lut = {str(r.region_name): (r.base_name, r.hemisphere)
                   for r in label_df.itertuples()}

        # pass 1: read every matrix, resolve its run
        raw = []
        for subj, ses in subs_sessions:
            hits = find_matrix_files(bids_dir, subj, ses, self.atlas_name, level_idx,
                                     self.correl_type)
            if not hits:
                continue
            runs = self._resolve_runs(cb, subj, ses, hits, level_idx)
            for i, csv in hits:
                token = runs[i][0]
                if token is None:
                    print(f"  [WARN] {csv}: run index {i} could not be resolved -- skipped")
                    continue
                try:
                    names, M = read_matrix(csv)
                except Exception as e:
                    print(f"  [WARN] {csv}: {e}")
                    continue
                raw.append(((subj, ses, token, seed_name), names, M, csv))
        if not raw:
            print(f"  [WARN] {cb['label']}: no {self.atlas_name}_{level_idx}_run_*_"
                  f"{self.correl_type}_matrix.csv found")
            return {}, pd.DataFrame()

        all_nodes = sorted({n for _, names, _, _ in raw for n in names})
        bases = {node_info(n, lut)[0] for n in all_nodes}
        if self.node_regions is not None:
            hit = bases & self.node_regions
            missing = sorted(self.node_regions - bases)
            print(f"  [{cb['label']}] matrix level index {level_idx}: {len(all_nodes)} nodes, "
                  f"{len(hit)}/{len(self.node_regions)} requested region(s) present"
                  + (f"; absent: {missing}" if missing else ""))
            if len(hit) < max(3, 0.5 * len(self.node_regions)):
                raise ValueError(
                    f"{cb['label']}: only {len(hit)} of the requested regions exist in the "
                    f"level-{level_idx} matrices (nodes there: {sorted(bases)[:12]}...). "
                    f"Wrong level? The seed ATLAS_LEVEL is 1-based, the matrix index "
                    f"0-based -- set MATRIX_LEVEL_INDEX explicitly if your files differ.")
            keep = hit
        else:
            keep = None
            print(f"  [{cb['label']}] matrix level index {level_idx}: {len(all_nodes)} "
                  f"nodes, all used")

        edges = edge_table(all_nodes, lut, keep)
        rn_list = [e[0] for e in edges]
        meta = pd.DataFrame(dict(region_name=rn_list, base_name=[e[1] for e in edges],
                                 hemisphere=[e[2] for e in edges],
                                 label_id=np.arange(1, len(edges) + 1)))
        # pass 2: every run gets the dataset-wide edge list (NaN where absent)
        profiles = {}
        for key, names, M, csv in raw:
            vals = {}
            for i in range(len(names)):
                for j in range(i + 1, len(names)):
                    vals[edge_names(names[i], names[j], lut)[0]] = M[i, j]
            prof = meta.copy()
            prof[self.value_col] = [vals.get(r, np.nan) for r in rn_list]
            if key in profiles:
                # Same (subject, session, run) for two acquisitions, e.g. acq-AP and
                # acq-PA with the same run number. SBA_Region_Profile_Index keys the
                # seed maps the same way and its dict keeps the LAST file in sorted
                # order; keeping the last here too pairs the same acquisition in
                # both analyses (seed vs matrix comparison).
                print(f"  [WARN] {cb['label']}: two matrices map to {key[:3]} (same session "
                      f"and run number: acq/task not in the run key) -- keeping the last "
                      f"({csv}), as the seed analysis does")
            profiles[key] = prof
        n_edges = meta["base_name"].nunique()
        print(f"  [{cb['label']}] {len(profiles)} run matrix(es) -> {len(rn_list)} hemispheric "
              f"edges, {n_edges} edge profile entries after L/R pooling")
        return profiles, pd.DataFrame()

    def mapping_table(self):
        return pd.DataFrame(self.mapping_rows)


def read_node_regions(seed_output_root):
    """The seed run's common regions -- nodes of the matrix analysis."""
    p = opj(seed_output_root, "02_selection", "common_regions.xlsx")
    if not os.path.isfile(p):
        return None
    return pd.read_excel(p)["common_region"].astype(str).str.lower().str.strip().tolist()


# ═══════════════════════════════════════════════════════════════════════════
# Are the selected runs the same?
# ═══════════════════════════════════════════════════════════════════════════

KEY = ["subject", "session", "run"]


def _read_selection(root):
    s = pd.read_excel(opj(root, "02_selection", "selection.xlsx"),
                      dtype={"subject": str, "session": str, "run": str})
    p = opj(root, "02_selection", "selection_by_round.xlsx")
    final = None
    if os.path.isfile(p):
        r = pd.read_excel(p, dtype={"subject": str, "session": str, "run": str})
        last = int(r["round"].max())
        final = r[r["round"] == last][["label"] + KEY].assign(kept_final=True)
        final.attrs["round"] = last
    return s, final


def _overlap(pool, a, b):
    """Hypergeometric test of |a & b| within `pool` (sets of keys)."""
    a, b = a & pool, b & pool
    N, n1, n2, k = len(pool), len(a), len(b), len(a & b)
    exp = n1 * n2 / N if N else np.nan
    p = float(st.hypergeom.sf(k - 1, N, n1, n2)) if N and n1 and n2 else np.nan
    jac = k / len(a | b) if (a | b) else np.nan
    return dict(pool=N, n_seed=n1, n_matrix=n2, overlap=k, expected=exp,
                p_hypergeom=p, jaccard=jac)


def compare_selections(seed_root, matrix_root, out_dir, n_mc=20000, seed=0):
    """Per dataset: Spearman of the round-1 scores of the runs both analyses
    scored, and the overlap of the kept sets (round 1 and final round) against
    the hypergeometric expectation of two independent draws from the shared
    pool. A global p sums the overlaps over datasets (Monte Carlo of
    independent hypergeometrics)."""
    os.makedirs(out_dir, exist_ok=True)
    s1, sf = _read_selection(seed_root)
    m1, mf = _read_selection(matrix_root)
    rows, scatter = [], {}
    rng = np.random.default_rng(seed)
    for lab in sorted(set(s1["label"]) & set(m1["label"])):
        a = s1[s1["label"] == lab]
        b = m1[m1["label"] == lab]
        mg = a.merge(b, on=KEY, suffixes=("_seed", "_matrix"))
        pool = set(map(tuple, mg[KEY].values))
        row = dict(label=lab, species=a["species"].iloc[0],
                   n_seed_units=len(a), n_matrix_units=len(b), n_shared_units=len(mg))
        if len(mg) >= 3:
            rho, p = st.spearmanr(mg["score_seed"], mg["score_matrix"])
            row.update(spearman_rho_round1=float(rho), spearman_p_round1=float(p))
        k1 = _overlap(pool, set(map(tuple, a.loc[a["is_kept"], KEY].values)),
                      set(map(tuple, b.loc[b["is_kept"], KEY].values)))
        row.update({f"{k}_round1": v for k, v in k1.items()})
        if sf is not None and mf is not None:
            kf = _overlap(pool, set(map(tuple, sf.loc[sf["label"] == lab, KEY].values)),
                          set(map(tuple, mf.loc[mf["label"] == lab, KEY].values)))
            row.update({f"{k}_final": v for k, v in kf.items()})
        rows.append(row)
        scatter[lab] = mg
    df = pd.DataFrame(rows)
    glob_rows = []
    for tag in ("round1", "final"):
        if f"overlap_{tag}" not in df:
            continue
        d = df.dropna(subset=[f"pool_{tag}"])
        d = d[(d[f"n_seed_{tag}"] > 0) & (d[f"n_matrix_{tag}"] > 0)]
        obs = int(d[f"overlap_{tag}"].sum())
        sims = np.zeros(n_mc)
        for r in d.itertuples():
            N, n1, n2 = int(getattr(r, f"pool_{tag}")), int(getattr(r, f"n_seed_{tag}")), \
                int(getattr(r, f"n_matrix_{tag}"))
            sims += rng.hypergeometric(n1, N - n1, n2, size=n_mc)
        glob_rows.append(dict(selection=tag, observed_overlap=obs,
                              expected_overlap=float(d[f"expected_{tag}"].sum()),
                              p_mc=(1 + int(np.sum(sims >= obs))) / (n_mc + 1),
                              n_datasets=len(d)))
    g = pd.DataFrame(glob_rows)
    with pd.ExcelWriter(opj(out_dir, "seed_vs_matrix_selection.xlsx")) as xw:
        df.to_excel(xw, sheet_name="per_dataset", index=False)
        g.to_excel(xw, sheet_name="all_datasets", index=False)
    final_round = sf.attrs.get("round") if sf is not None else None
    plot_selection_comparison(df, g, scatter, opj(out_dir, "F19_seed_vs_matrix_selection.png"),
                              final_round)
    ff._legend(
        opj(out_dir, "F19_seed_vs_matrix_selection.txt"),
        "F19 -- Do the seed and the matrix analyses keep the same runs?",
        "(a) per dataset, number of runs kept by BOTH analyses (bar) vs the number "
        "expected if the two selections were independent draws from the runs both "
        "analyses scored (black tick); round 1 and final round. (b) Spearman correlation "
        "of the round-1 scores. (c) round-1 scores, seed vs matrix, kept runs marked.",
        "Pool = runs present in both analyses. Expected overlap = n_seed x n_matrix / "
        "pool; p = hypergeometric upper tail. Global p: Monte Carlo sum of independent "
        "hypergeometric overlaps across datasets.",
        "Both analyses share the run's data quality (motion, global signal, tSNR), and "
        "the matrix contains the retrosplenial edges (a small fraction of all edges), so "
        "some agreement is expected even without a shared biological pattern. Overlap "
        "near the expectation means the selection is specific to the profile it scores -- "
        "the seed selection then cannot be described as 'keeping the good runs' in "
        "general. High overlap means both pick the same runs: then check F10/external QC "
        "to see whether that shared factor is quality rather than biology.")
    return df, g


def plot_selection_comparison(df, g, scatter, out_png, final_round):
    labs = list(df["label"])
    fig = plt.figure(figsize=(17, 9.5))
    gs = fig.add_gridspec(2, 5, height_ratios=[1, 1.05], hspace=0.55, wspace=0.35)
    ax = fig.add_subplot(gs[0, :3])
    x = np.arange(len(labs))
    for off, tag, name in ((-0.2, "round1", "round 1"),
                           (0.2, "final", f"final round ({final_round})")):
        if f"overlap_{tag}" not in df:
            continue
        cols = [ff._color(s) for s in df["species"]]
        ax.bar(x + off, df[f"overlap_{tag}"], width=0.38, color=cols,
               alpha=0.55 if tag == "round1" else 0.95, edgecolor="black", lw=0.5,
               label=name)
        ax.scatter(x + off, df[f"expected_{tag}"], marker="_", s=260, color="black", zorder=4)
        for xi, r in zip(x + off, df.itertuples()):
            p = getattr(r, f"p_hypergeom_{tag}")
            n1 = getattr(r, f"n_seed_{tag}")
            ax.text(xi, getattr(r, f"overlap_{tag}") + 0.15,
                    f"{getattr(r, f'overlap_{tag}'):.0f}/{n1:.0f}  p={ff._fmt_p(p)}",
                    ha="center", va="bottom", rotation=90, fontsize=6.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labs, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("runs kept by both")
    ax.set_ylim(0, ax.get_ylim()[1] * 1.45)
    gl = "; ".join(f"{r.selection}: {r.observed_overlap} vs {r.expected_overlap:.1f} expected, "
                   f"p = {ff._fmt_p(r.p_mc)}" for r in g.itertuples())
    ax.set_title(f"a  Overlap of the kept runs (light = round 1, dark = final; tick = "
                 f"expected if independent)\nall datasets: {gl}", fontsize=9,
                 fontweight="bold", loc="left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax = fig.add_subplot(gs[0, 3:])
    if "spearman_rho_round1" in df:
        ax.barh(x[::-1], df["spearman_rho_round1"], color=[ff._color(s) for s in df["species"]])
        for xi, r in zip(x[::-1], df.itertuples()):
            ax.text(r.spearman_rho_round1 + (0.02 if r.spearman_rho_round1 >= 0 else -0.02),
                    xi, f"n={r.n_shared_units}", va="center", fontsize=7,
                    ha="left" if r.spearman_rho_round1 >= 0 else "right")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(x[::-1])
    ax.set_yticklabels(labs, fontsize=8)
    ax.set_xlim(-1, 1)
    ax.set_xlabel("Spearman rho, round-1 score seed vs matrix")
    ax.set_title("b  Do the two analyses rank runs alike?", fontsize=9, fontweight="bold",
                 loc="left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    sub = gs[1, :].subgridspec(2, 5, hspace=0.6, wspace=0.35)
    for k, lab in enumerate(labs[:10]):
        ax = fig.add_subplot(sub[k // 5, k % 5])
        mg = scatter[lab]
        both = mg["is_kept_seed"] & mg["is_kept_matrix"]
        ax.scatter(mg["score_seed"], mg["score_matrix"], s=8, color="#BBBBBB", lw=0)
        ax.scatter(mg.loc[mg["is_kept_seed"], "score_seed"],
                   mg.loc[mg["is_kept_seed"], "score_matrix"], s=16, facecolor="none",
                   edgecolor="#0072B2", lw=0.8)
        ax.scatter(mg.loc[mg["is_kept_matrix"], "score_seed"],
                   mg.loc[mg["is_kept_matrix"], "score_matrix"], s=16, marker="x",
                   color="#D55E00", lw=0.8)
        ax.set_title(f"{lab} ({int(both.sum())} both)", fontsize=7.5,
                     color=ff._color(mg["species_seed"].iloc[0]) if "species_seed" in mg
                     else "black")
        ax.tick_params(labelsize=6)
        if k % 5 == 0:
            ax.set_ylabel("matrix score", fontsize=7)
        if k // 5 == 1:
            ax.set_xlabel("seed score", fontsize=7)
    fig.text(0.5, 0.47, "c  Round-1 scores (blue ring = kept by seed, orange x = kept by "
             "matrix)", ha="center", fontsize=9, fontweight="bold")
    fig.suptitle("Seed (retrosplenial) vs whole-brain matrix selection", fontsize=11,
                 fontweight="bold")
    ff._save(fig, out_png)


# ═══════════════════════════════════════════════════════════════════════════
# F21 -- mean EDNIxCSCLR matrices of the runs kept by the SEED selection
# ═══════════════════════════════════════════════════════════════════════════

def _node_sort_key(name):
    m = re.match(r"^([LlRr])_(.*)$", str(name))
    hemi, base = (m.group(1).upper(), m.group(2).lower()) if m else ("", str(name).lower())
    order = ff.NODE_ORDER.index(base) if base in ff.NODE_ORDER else len(ff.NODE_ORDER)
    return ({"L": 0, "R": 1}.get(hemi, 2), order, base)


def _node_label(name):
    m = re.match(r"^([LlRr])_(.*)$", str(name))
    if not m:
        return ff.short_region(name)
    return f"{m.group(1).upper()} {ff.short_region(m.group(2).lower())}"


def run_matrix(prof, value_col, nodes):
    """Node x node matrix of one run from its edge profile (region_name = 'A <> B')."""
    idx = {n: i for i, n in enumerate(nodes)}
    M = np.full((len(nodes), len(nodes)), np.nan)
    for name, v in zip(prof["region_name"].astype(str), prof[value_col].astype(float)):
        if EDGE_LINK not in name:
            continue
        a, b = name.split(EDGE_LINK, 1)
        if a in idx and b in idx:
            M[idx[a], idx[b]] = M[idx[b], idx[a]] = v
    return M


def mean_matrix(mats):
    """Fisher-z average of correlation matrices, back to r."""
    with np.errstate(invalid="ignore", divide="ignore"):
        Z = np.arctanh(np.clip(np.stack(mats), -0.999999, 0.999999))
        return np.tanh(np.nanmean(Z, axis=0))


def plot_mean_matrices(seed_root, matrix_runs_by_label, species_of_label, out_dir,
                       value_col="median_r"):
    """Per dataset, mean EDNIxCSCLR matrix (nodes = the seed's common regions,
    left then right hemisphere) over every run, over the runs KEPT by the
    seed's final selection round, and over the runs not kept. One colour scale
    per dataset. Writes F21_mean_matrices_seed_selection.png/.pdf + .txt + .xlsx."""
    os.makedirs(out_dir, exist_ok=True)
    _, sf = _read_selection(seed_root)
    if sf is None:
        print("  [F21] no final-round seed selection found -- skipped")
        return None
    rows, panels = [], []
    for lab in sorted(matrix_runs_by_label):
        profs = matrix_runs_by_label[lab]
        if not profs:
            continue
        kept = set(map(tuple, sf.loc[sf["label"] == lab, KEY].astype(str).values))
        keys = list(profs)
        nodes = set()
        for p in profs.values():
            for nm in p["region_name"].astype(str):
                if EDGE_LINK in nm:
                    nodes |= set(nm.split(EDGE_LINK, 1))
        nodes = sorted(nodes, key=_node_sort_key)
        mats = {k: run_matrix(profs[k], value_col, nodes) for k in keys}
        is_kept = {k: tuple(str(x) for x in k[:3]) in kept for k in keys}
        groups = dict(all=[mats[k] for k in keys],
                      kept=[mats[k] for k in keys if is_kept[k]],
                      not_kept=[mats[k] for k in keys if not is_kept[k]])
        if not groups["kept"]:
            print(f"  [F21] {lab}: no matrix matches the seed's kept runs -- skipped")
            continue
        means = {g: mean_matrix(m) for g, m in groups.items() if m}
        iu = np.triu_indices(len(nodes), 1)
        r_kn = (np.corrcoef(means["kept"][iu], means["not_kept"][iu])[0, 1]
                if "not_kept" in means else np.nan)
        rows.append(dict(label=lab, species=species_of_label.get(lab), n_nodes=len(nodes),
                         n_all=len(groups["all"]), n_kept=len(groups["kept"]),
                         n_not_kept=len(groups["not_kept"]),
                         mean_r_all=float(np.nanmean(means["all"][iu])),
                         mean_r_kept=float(np.nanmean(means["kept"][iu])),
                         mean_r_not_kept=float(np.nanmean(means["not_kept"][iu]))
                         if "not_kept" in means else np.nan,
                         r_kept_vs_not_kept_pattern=float(r_kn)))
        panels.append((lab, nodes, means, {g: len(m) for g, m in groups.items()}))
    if not panels:
        return None
    cols = [("all", "all runs"), ("kept", "selected"), ("not_kept", "not selected")]
    fig, axes = plt.subplots(len(panels), 3, figsize=(13.5, 4.3 * len(panels)), squeeze=False)
    for i, (lab, nodes, means, ns) in enumerate(panels):
        vals = np.concatenate([m[np.triu_indices(len(nodes), 1)] for m in means.values()])
        vmax = float(np.nanpercentile(np.abs(vals), 99)) or 1.0
        im = None
        for j, (g, name) in enumerate(cols):
            ax = axes[i, j]
            if g not in means:
                ax.axis("off")
                continue
            M = means[g].copy()
            np.fill_diagonal(M, np.nan)
            im = ax.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax, interpolation="nearest")
            ax.set_xticks(range(len(nodes)))
            ax.set_yticks(range(len(nodes)))
            show = len(nodes) <= 40
            ax.set_xticklabels([_node_label(n) for n in nodes] if show and i == len(panels) - 1
                               else [], rotation=90, fontsize=5)
            ax.set_yticklabels([_node_label(n) for n in nodes] if show and j == 0 else [],
                               fontsize=5)
            ax.set_title(f"n = {ns[g]} runs", loc="right", fontsize=9)
            if j == 0:
                ax.set_title(lab, loc="left", fontsize=10, fontweight="bold")
            if i == 0:
                ax.text(0.5, 1.12, name, transform=ax.transAxes, ha="center", fontsize=11)
        cb = fig.colorbar(im, ax=axes[i, :].tolist(), shrink=0.8, pad=0.01)
        cb.set_label("mean r (Fisher z average)")
    out = opj(out_dir, "F21_mean_matrices_seed_selection")
    with plt.rc_context({"pdf.fonttype": 42}):
        for ext in ("png", "pdf"):
            fig.savefig(f"{out}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [saved] {out}.png/.pdf")
    tab = pd.DataFrame(rows)
    tab.to_excel(out + ".xlsx", index=False)
    ff._legend(out + ".txt",
               "F21 -- Mean EDNIxCSCLR connectivity matrices of the runs kept by the SEED "
               "selection",
               "Per dataset (rows), the mean region x region correlation matrix (EDNIxCSCLR, "
               "nodes = the seed analysis' common regions, left hemisphere then right) over "
               "every run, over the runs kept by the seed's final selection round, and over the "
               "runs it did not keep. Average of Fisher z, shown back in r. One colour scale "
               "per dataset (+/- 99th percentile of |r|).",
               "Descriptive. In the table (F21_..xlsx): mean off-diagonal r of each matrix and "
               "r_kept_vs_not_kept_pattern, the correlation between the kept and the not-kept "
               "matrices across edges (how much the selection changes the PATTERN, as opposed "
               "to the level).",
               "The selection used the retrosplenial seed profile only. If the kept matrix "
               "shows a sharper whole-brain structure (stronger homotopic and within-network "
               "edges) than the not-kept one, the seed criterion also picks runs whose whole "
               "connectome is cleaner. A higher mean r alone is not that: the global level of "
               "correlation (no global signal regression) raises every edge. Kept and not-kept "
               "groups differ in size (12 vs the rest), so the not-kept mean is less noisy.")
    return tab
