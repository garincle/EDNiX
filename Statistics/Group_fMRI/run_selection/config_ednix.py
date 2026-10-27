"""
Datasets and selection settings of the EDNiX paper (retrosplenial seed, atlas level 3).

Single source of truth: launcher.py (part 1, this package) and the paper analysis
(Exemples/Study_EDNiX/analysis) both import these values, so the selection that the
paper analyses is the one launched here. Values copied unchanged from
Exemples/Study_EDNiX/analysis/Fig5/Launcher_Fig5.py (V15, 2026-10-27).

To add a species: add one dict to BIDS_CONFIGS (label, species, bids_dir, label_path,
seeds, state) and its tree tip to selection.CONFIG["SPECIES_TO_TREE_TIP"].
"""
import pandas as pd

# ── profile of each seed map (changes the selection) ─────────────────────────────
# "percentiles"     : 2 values per region (upper / lower TAIL_Q percentile of the voxels)
# "median_r"        : 1 value per region (median of the voxels; original analysis)
# "fraction_posneg" : 2 values per region (% of voxels in the upper / lower tails)
PROFILE_METHOD = "percentiles"
TAIL_Q = 0.2             # 0.2 / "quintile" (P80 / P20), 0.1 / "decile", 0.25 / "quartile"
_Q_PCT = int(round(100 * float({"decile": 0.1, "deciles": 0.1, "quintile": 0.2, "quintiles": 0.2,
                                "quartile": 0.25, "quartiles": 0.25}.get(TAIL_Q, TAIL_Q))))

LIB = "/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals"
SEED = ["BA_26_29_30__retrosplenial"]

# ── datasets ─────────────────────────────────────────────────────────────────────
# label      : dataset name everywhere (folders, figures, AFNI balance json) -- changing
#              it recomputes its kept / not-kept AFNI maps
# species    : must exist in selection.CONFIG["SPECIES_TO_TREE_TIP"]
# state      : arousal state. String = whole dataset; {"from": "task"} = from task-<x> of
#              the BOLD name; {"from": "subject"} = from the subject name (first pattern
#              found); None = unknown. Anything but "awake" = anaesthetised. "animal"
#              (regex) is removed from the subject name to recover the animal (Mouse).
STATE_HORI = {"from": "task", "map": {"restawake": "awake", "restanest": "anesthetized"}}
STATE_MOUSE = {"from": "subject", "map": [["MEDISO", "medetomidine+isoflurane"],
                                          ["ISO", "isoflurane"],
                                          ["AW", "awake"]],
               # the animal is the suffix: jgrAesAWc11L, jgrAesISOc11L, jgrAesMEDISOc11L
               "animal": r"^jgrAes(?:MEDISO|ISO|AW)"}

BIDS_CONFIGS = [
    dict(label="Dog", species="Dog", bids_dir="/scratch2/EDNiX/Dog/BIDS_Boch_K9",
         label_path=f"{LIB}/carnivorans/Dog/label_code/EDNIxCSCLR_ITK.txt", seeds=SEED,
         state="awake", endfmri="*_task-rest_*.nii.gz"),
    dict(label="Macaque_BenHamed", species="Macaque", bids_dir="/scratch2/EDNiX/Macaque/BIDS_BenHamed",
         label_path=f"{LIB}/primates/Macaque/label_code/EDNIxCSCLR_ITK.txt", seeds=SEED,
         state="awake", endfmri="*_task-rest_*.nii.gz"),
    dict(label="Macaque_ZhuGarin", species="Macaque", bids_dir="/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin",
         label_path=f"{LIB}/primates/Macaque/label_code/EDNIxCSCLR_ITK.txt", seeds=SEED,
         state="isoflurane", endfmri="*_task-rest_*.nii.gz"),
    dict(label="Rat", species="Rat", bids_dir="/scratch2/EDNiX/Rat/BIDS_Grandjean",
         label_path=f"{LIB}/rodents/Rat/label_code/EDNIxCSCLR_ITK.txt", seeds=SEED,
         state="anesthetized", endfmri="*_task-rest_bold.nii.gz"),
    dict(label="Mouse", species="Mouse", bids_dir="/scratch2/EDNiX/Mouse/BIDS_Grandjean2",
         label_path=f"{LIB}/rodents/Mouse/label_code/EDNIxCSCLR_ITK.txt", seeds=SEED,
         state=STATE_MOUSE, endfmri="*_task-rest_*.nii.gz"),
    dict(label="Marmoset_Tian", species="Marmoset", bids_dir="/scratch2/EDNiX/Marmoset/BIDS_Tian",
         label_path=f"{LIB}/primates/Marmoset/label_code/EDNIxCSCLR_ITK.txt", seeds=SEED,
         state="awake", endfmri="*_task-rest_run-RL*.nii.gz"),
    dict(label="Marmoset_Hori", species="Marmoset", bids_dir="/scratch2/EDNiX/Marmoset/BIDS_Hori",
         label_path=f"{LIB}/primates/Marmoset/label_code/EDNIxCSCLR_ITK.txt", seeds=SEED,
         state=STATE_HORI, endfmri="*_task-rest*.nii.gz"),
    dict(label="Mouselemur", species="Mouselemur", bids_dir="/scratch2/EDNiX/Mouselemur/BIDS_Garin",
         label_path=f"{LIB}/primates/Mouselemur/label_code/EDNIxCSCLR_ITK.txt", seeds=SEED,
         state="anesthetized", endfmri="*_task-rest_*.nii.gz"),
    dict(label="Human_Castrillon", species="Human", bids_dir="/scratch2/EDNiX/Human/BIDS_Castrillon",
         label_path=f"{LIB}/primates/Human/label_code/EDNIxCSCLR_ITK.txt", seeds=SEED,
         state="awake", endfmri="*_task-rest_bold.nii.gz"),
    dict(label="Human_Babayan", species="Human", bids_dir="/scratch2/EDNiX/Human/BIDS_Babayan",
         label_path=f"{LIB}/primates/Human/label_code/EDNIxCSCLR_ITK.txt", seeds=SEED,
         state="awake", endfmri="*_task-rest_*.nii.gz"),
]
for c in BIDS_CONFIGS:
    c.update(reference_subject=None, reference_session=None, reference_run=None)

# ── selection settings (keys of selection.CONFIG; absent = default) ──────────────
SELECTION_CONFIG = dict(
    # one output folder per profile method; the AFNI maps already computed with the
    # same inputs (checked in input_runs/) are reused
    OUTPUT_ROOT="/scratch2/EDNiX/results/EDNiX_Figures/Fig_5_V15_"
                + {"median_r": "", "fraction_posneg": f"_posneg{_Q_PCT}", "percentiles": f"_quant{_Q_PCT}",
                   "quantile_posneg": f"_quant{_Q_PCT}"}[PROFILE_METHOD],
    PROFILE_METHOD=PROFILE_METHOD,
    TAIL_Q=TAIL_Q,
    # data and unit of analysis
    UNIT="run",              # "run" | "session" | "subject"
    ATLAS_LEVEL=3,           # EDNiX atlas level of the regions (2 or 3)
    RUN_QC_ENABLE=False,     # run QC first (MAD distance to the dataset medoid)
    # selection
    N_SELECT=10,             # units kept per dataset at every round
    INVARIANCE_MODE="top_fraction",  # regions used by the selection score:
    INVARIANCE_VALUE=0.8,    # the 80 % most invariant across species (1.0 = all)
    SIMILARITY="correlation",  # score = Pearson r between the profile and the target
    ZSCORE="scale",          # profile normalisation: "scale" (RMS, sign kept), True, False
    LEAVE_ONE_OUT=True,      # a dataset's target excludes its whole species
    RUN_ROUND2_SELECTION=True,     # round 2: each dataset re-ranked against its own kept
    RUN_ITERATIVE_SELECTION=True,  # rounds 3..N: cross-species target rebuilt from the kept
    N_ITERATIVE_ROUNDS=10,         # number of the LAST round (final selection)
    SEED_REGION_EXCLUDE="auto",    # atlas region containing the seed: "auto" | names | None
    SEED_REGION_EXCLUDE_SELECTION=False,  # kept in the selection (= V12 selection)
    SEED_REGION_EXCLUDE_ANALYSES=True,    # removed from every analysis
    # AFNI maps
    OVERWRITE_MAPS=True,     # False: a map computed with exactly the same runs is reused
    MAP_TAG=None,            # None = unique suffix derived from OUTPUT_ROOT
    RENDER_SURFACE=True,     # wb_command renders of the group maps
    # time tree: only its file name is checked here (used by the statistics)
    TREE_PATH="/scratch2/EDNiX/Liste des espèces.nwk",
)

# Control of the same-species effect: every dataset is selected towards the SAME
# target (consensus of all species). Python only; common regions = those of the
# main run (COMMON_REGIONS_OVERRIDE, set by the launcher).
CONTROL_SELECTION_CONFIG = dict(
    SELECTION_CONFIG,
    OUTPUT_ROOT=SELECTION_CONFIG["OUTPUT_ROOT"] + "_control_common_target",
    LEAVE_ONE_OUT=False,
    RUN_ALL_UNITS_MAPS=False, RUN_GROUP_MAPS=False, RENDER_SURFACE=False,
    COPY_JPG=False,
)

# ── external QC (tSNR, FD... per run): selection validation and PGLS covariates ──
QC_EXCEL_PATH = "/scratch2/EDNiX/results/EDNiX_Figures/Suppl_Figure_4-5-6_QC/qc_raw_per_individual.xlsx"


def load_qc_by_label(bids_configs=BIDS_CONFIGS, path=QC_EXCEL_PATH):
    """{label: QC rows of that dataset}. EXACT match on the last element of bids_dir
    (str.contains matched "BIDS_Grandjean" (Rat) inside "BIDS_Grandjean2" (Mouse))."""
    df_qc_all = pd.read_excel(path, dtype={"subject": str, "session": str, "run": str})
    qc_bids_name = (df_qc_all["bids_dir"].astype(str).str.strip().str.rstrip("/")
                    .str.split("/").str[-1])
    out = {}
    for cb in bids_configs:
        bids_name = cb["bids_dir"].rstrip("/").split("/")[-1]
        mask = qc_bids_name == bids_name
        if mask.any():
            out[cb["label"]] = df_qc_all[mask].copy()
        else:
            print(f"  [WARN] {cb['label']}: no QC row for '{bids_name}'")
    print(f"  df_qc_by_label: {len(out)}/{len(bids_configs)} dataset(s) with QC")
    return out


def read_common_regions(output_root):
    """Common regions of a finished selection (02_selection/common_regions.xlsx), or None."""
    import os
    p = os.path.join(output_root, "02_selection", "common_regions.xlsx")
    if not os.path.isfile(p):
        return None
    return pd.read_excel(p)["common_region"].astype(str).str.lower().str.strip().tolist()
