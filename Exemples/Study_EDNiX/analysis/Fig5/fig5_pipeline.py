"""
fig5_pipeline.py
================
ONE entry point for the whole Figure 5 analysis. Replaces the scattered
run_*.py drivers: everything is controlled from the CONFIG block below,
and every output lands in one predictable tree under OUTPUT_ROOT.

WHY THIS EXISTS
---------------
Four concrete problems in the previous organisation, all fixed here:

1. THE UNIT OF ANALYSIS IS NO LONGER THE RUN (`UNIT` in CONFIG, default
   "subject"). Everything -- the per-dataset profile, the cross-species
   consensus, the selection, and the 3dttest++ group maps -- now treats
   one ANIMAL as one sample. Before, a subject scanned 20 times counted
   20 times while a subject scanned once counted once, so a dataset's
   profile described its most-scanned animals, N_SELECT "best runs"
   could all come from a single animal, and 3dttest++ spent its degrees
   of freedom on non-independent scans. Runs are still the raw material:
   they are averaged within each subject, in Python for the profiles
   (median) and by AFNI 3dMean for the group maps.

2. STEP 1 ACTUALLY PRODUCES A MAP. It used to call
   `rerun_balanced_by_run` -- a kept-vs-outlier function -- with one
   arbitrary run flagged as a fake outlier just to satisfy it. That
   function then BALANCED the two groups, so min(N-1, 1) = 1: the
   "all runs" map was computed on a single run, AFNI refused it
   ("[SKIP] ...: only 1 run(s)"), and every species got an empty output
   folder while summary.xlsx still claimed 24/39/80/170 runs. Step 1 now
   passes balance=False and keeps every unit, and the `_unused` folder is
   gone.

3. Two different folder-naming conventions existed at once --
   Figure_5_SBA.py wrote `Grp_SBA_3dTTEST_All<suffix>` (category glued
   on) while rerun_balanced_by_run wrote `Grp_SBA_3dTTEST_<suffix>`.
   Lookups guessed wrong and silently reported "group map not found at
   None". Here, ONE convention is used everywhere and the resolved path
   is printed and recorded for every map.

4. Outputs were spread across several roots with opaque names. Here the
   tree is fixed (see OUTPUT LAYOUT), every figure gets a .txt legend,
   and a methods.txt is written with the exact parameters used.

DEPENDS ON
----------
- fig5_analysis.py   (unchanged: it works on generic profile dicts)
- fig5_units.py      (new: the run -> subject aggregation, and the
                      step-1 figure read off the AFNI group maps)
- SBA_Region_Profile_Index.py, patched by patch_SBA_unit_level.py
  (adds `unit=` and `balance=` to rerun_balanced_by_run). Running this
  pipeline against an unpatched module raises a TypeError immediately,
  by design -- silently falling back to run-level would undo point 1.

OUTPUT LAYOUT (under OUTPUT_ROOT/)
----------------------------------
    methods.txt                     parameters + method text for the paper
    00_log.txt                      full console log
    01_all_units/                   group maps from ALL units, per BIDS
        summary.xlsx                what actually entered each map
        map_index.xlsx              every map's resolved path
        F0_species_profiles_afni.png   species profiles read off those maps,
        F0_species_profiles_afni.txt   with the cross-species core shaded
    02_selection/                   the consensus selection itself
        selection.xlsx              per-unit scores, kept flag
        common_regions.xlsx
        region_invariance.xlsx
        unit_counts.xlsx            units available per dataset
        balanced_stats.xlsx         kept vs not-kept, size-matched (descriptive, no p)
        external_qc.xlsx            non-circular QC validation (optional)
    03_group_maps/                  unit-level 3dttest for kept AND not-kept
        map_index.xlsx              every map's resolved path
    04_jpg/                         SBA previews
    05_surfaces/                    every surface .png, renamed unambiguously
        <label>__<group>__<stat>.png
    06_figures/                     all figures + one .txt legend each
"""

import os
import sys
import glob
import shutil
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

import Exemples.Study_EDNiX.analysis.Fig5.fig5_analysis as fa   # pure numpy/pandas/matplotlib
import Exemples.Study_EDNiX.analysis.Fig5.fig5_units as fu      # unit aggregation + step-1 figure
import Exemples.Study_EDNiX.analysis.Fig5.fig5_phylo as fp      # Mantel / Blomberg K / PGLS

# ═══════════════════════════════════════════════════════════════════════════
# CONFIG -- everything controllable lives here
# ═══════════════════════════════════════════════════════════════════════════

CONFIG = dict(
    # --- where everything goes -------------------------------------------
    OUTPUT_ROOT="/scratch2/EDNiX/results/EDNiX_Figures/Fig_5",

    # --- paths formerly read from Figure5_SBA_KeeporNot -------------------
    RESULTS_BASE_DIR=None,  # None (recommandé) = <OUTPUT_ROOT>/03_group_maps,
        # dérivé automatiquement plus bas -- chaque run garde ses propres
        # statmaps AFNI isolées. Mettre un chemin explicite ici les fait
        # toutes atterrir dans un dossier PARTAGÉ entre runs différents ;
        # avec OVERWRITE_MAPS=True, deux lancements avec un SEED/ATLAS_LEVEL
        # identiques mais un SIMILARITY/ZSCORE/N_SELECT différent
        # s'écraseraient silencieusement l'un l'autre au même endroit.
    CLASSIFICATION_PATH=("/scratch2/EDNiX/results/EDNiX_Figures/Figure_4_rsfMRI_QC/"
                         "lvl2/all/classification_results.xlsx"),
    CATEGORY="All",
    STAT_SUFFIX="ttest-stat_fisher_cmap",
    STAT_MAP_USE_THRESHOLDED=False,  # False (correct) : moyenne Fisher-z NON
        # seuillée de 3dttest++. True : version "_thresholded_" (|z|>=1.96,
        # voxels non significatifs mis à 0 puis EXCLUS des profils) -- profil
        # calculé sur les seuls voxels significatifs : biais de sélection
        # (Kriegeskorte et al. 2009), couverture des régions dépendante de N.
    MAP_TAG=None,           # None = tag unique dérivé de OUTPUT_ROOT, ajouté à
        # chaque suffixe AFNI : deux runs en parallèle ne partagent plus le
        # dossier intermédiaire <bids_dir>/Results/Grp_SBA_3dTTEST_<suffix>.
        # "" = pas de tag (ancien nommage).
    REQUIRE_ALL_MAPS=True,  # True : arrêt si une carte all-units manque (la
        # phylogénie sans toutes les espèces n'a pas de sens).

    # --- what is analysed -------------------------------------------------
    SEED="BA_26_29_30__retrosplenial",
    ATLAS_LEVEL=2,          # THE atlas level knob: 2 or 3. Used everywhere.
    VALUE_COL="median_r",
    METRIC="correlation",

    # --- THE UNIT OF ANALYSIS --------------------------------------------
    # "subject" (default): one animal = one sample, everywhere. Runs are
    #     averaged within subject -- median for the Python profiles,
    #     AFNI 3dMean for the group maps.
    # "session": one scanning session = one sample. Use when sessions are
    #     far enough apart to be treated as separate observations and you
    #     accept that a subject with 11 sessions then counts 11 times.
    # "run": the old behaviour. Kept only so the previous results can be
    #     reproduced for comparison; it weights datasets by scan count.
    UNIT="subject",

    # --- selection --------------------------------------------------------
    N_SELECT=5,             # UNITS kept per BIDS (subjects, with UNIT="subject")
    LEAVE_ONE_OUT=True,     # species excluded from the consensus that judges it
    CONSENSUS_SOURCE="all_runs",   # "all_runs" | "reference"
    SIMILARITY="correlation",      # "correlation" | "cosine" | "neg_euclidean"
    ZSCORE=True,            # True/"zscore" | "scale" | False -- normalisation par
        # profil avant consensus. True recentre chaque profil (crée des
        # "négatifs" artificiels) ; "scale" divise par la RMS SANS recentrer
        # (signe réel conservé, amplitudes égalisées -- recommandé si le
        # signe compte) ; False = brut (signe conservé mais consensus en
        # patchwork dominé par les datasets de forte amplitude).
    THRESHOLD_ABS=None,     # None = keep every region (recommended)
    REFERENCE_METHOD="kcore",
    K_REF=5,

    # --- restrict selection to cross-species INVARIANT regions ------------
    INVARIANCE_MODE="top_fraction",   # "top_fraction" | "threshold" | "all"
    INVARIANCE_VALUE=1.0,             # 1.0 = keep every common region (default),
        # 0.5 = top 50% most invariant, 0 = none. Only controls which regions
        # the RUN SELECTION is restricted to; every common region is plotted.

    # --- heavy steps (each independently switchable) ----------------------
    RUN_ALL_UNITS_MAPS=True,    # step 1: group map from ALL units, per BIDS
    RUN_GROUP_MAPS=True,        # step 3: group maps for kept AND not-kept
    RENDER_SURFACE=True,        # surface rendering on top of the group maps
    COPY_JPG=True,
    MAKE_FIGURES=True,

    # --- stats ------------------------------------------------------------
    N_DRAWS_BALANCED=1000,      # subsampling draws for size-matched tests
    N_PERM_PHYLOGENY=10000,

    RUN_QC_ENABLE=True,
    RUN_QC_MAD_THRESHOLD=3.5,  # même défaut que classify_runs_by_profile_distance
    COMMON_REGIONS_REQUIRE_NONNAN=False,  # False: garde une région si présente
        # structurellement dans chaque BIDS, ignore juste les NaN ponctuels
        # (nanmedian/isfinite en aval). True: retire la région PARTOUT dès
        # qu'un seul run, n'importe où, est NaN pour elle (ancien comportement).
    OVERWRITE_MAPS=False,  # False: si le .nii de sortie d'un 3dttest++ existe
        # déjà sur disque (même suffixe), on réutilise son chemin tel quel et
        # on NE relance PAS l'appel AFNI (ClustSim est long et fragile --
        # cf. les crashs NFS/wb_command déjà vus). True: supprime les anciens
        # fichiers de ce suffixe puis recalcule.
    SURFACE_PNG_KEEP=("ttest-stat_fisher_cmap_oversampled",
                      "ttest-stat_fisher_zmap_oversampled"),
        # 05_surfaces/ ne garde que les pngs dont le nom contient une de ces
        # sous-chaînes. () ou None pour tout garder (ancien comportement).

    # --- F6 multi-niveaux (impact de la segmentation) ----------------------
    RUN_MULTILEVEL_F6=False,  # sensibilité seulement : aux niveaux 1, 2 et 4 trop
        # de régions manquent dans au moins une espèce (6, 11 et 10 régions
        # communes, contre 17 au niveau 3) ; le résultat de référence est
        # ATLAS_LEVEL=3.
    MULTILEVEL_F6_LEVELS=(1, 2, 3, 4),  # relit les cartes all-units déjà
        # construites à ces niveaux d'atlas -- pas de recalcul AFNI.

    # --- distance de chaque run à la référence "kept" -----------------------
    RUN_DISTANCE_TO_REF=False,
    DISTANCE_REF_GROUP="kept",  # quelle carte de 03_group_maps sert de réf

    # --- corrélation QC vs distance à la réf (besoin de df_qc_by_label) -----
    QC_METRICS_FOR_CORRELATION=None,  # None = auto-détecte func_/anat_/net_/
        # homotopic_/cross_/sp_/contam_ ; sinon liste explicite de colonnes.
    QC_META_MIN_BIDS=3,   # nb minimum de BIDS avec un r calculé pour tenter
        # le test méta inter-espèces sur une métrique donnée.
    QC_META_ALPHA=0.05,   # seuil de significativité du test méta pour F10.

    # --- chaîne de sélection ----------------------------------------------
    # round 1 : consensus inter-espèces (tous les runs)          -- toujours
    # round 2 : ancré espèce (score vs sa propre moyenne kept1)  -- RUN_ROUND2_SELECTION
    # rounds 3..N : consensus inter-espèces reconstruit à partir
    #   des moyennes kept du tour précédent                      -- RUN_ITERATIVE_SELECTION
    # Le DERNIER tour exécuté ("final round") a ses propres cartes AFNI
    # (03b_group_maps_final/, suffixe Kept<N>/NotKept<N>), ses JPG (04b) et
    # surfaces (05b), sa ligne sur F4/F1b, et alimente F6/F8/F9/F11-F13 si
    # PHYLOGENY_USE_FINAL_ROUND=True.
    RUN_ROUND2_SELECTION=False,
    RUN_ITERATIVE_SELECTION=False,
    N_ITERATIVE_ROUNDS=3,   # numéro du DERNIER tour (tours 1 et 2 inclus) ;
        # ignoré si RUN_ITERATIVE_SELECTION=False.
    PHYLOGENY_USE_FINAL_ROUND=False,  # False (défaut, recommandé) : la phylo
        # utilise les cartes all-units, sans sélection. True : cartes kept du
        # dernier tour -- biais de taille de pool (winner's curse) et, avec
        # LEAVE_ONE_OUT=False, circularité ; à présenter comme sensibilité.
    RUN_REFERENCE_METHOD_COMPARISON=False,

    # --- tests phylogénétiques comparatifs (arbre daté complet) -----------
    TREE_PATH="/scratch2/EDNiX/Liste_des_espèces.nwk",  # arbre daté (Newick) :
        # âges des nœuds pour F6/F7, Mantel, K, PGLS. Le pipeline s'arrête si le
        # fichier est absent (None = anciennes valeurs approximatives codées en
        # dur, à ne pas publier).
    SPECIES_TO_TREE_TIP={"Human": "Homo_sapiens", "Macaque": "Macaca_mulatta",
                         "Marmoset": "Callithrix_jacchus",
                         "Mouselemur": "Microcebus_murinus",
                         "Rat": "Rattus_norvegicus", "Mouse": "Mus_musculus",
                         "Dog": "Canis_lupus"},
    RUN_PHYLO_COMPARATIVE=True,   # Mantel exact + K de Blomberg + PGLS,
        # Python seul (quelques secondes). Principal test phylo, remplace
        # la régression 1-D "distance vs divergence depuis l'humain".
    PGLS_QC_COVARIATES=("func_TSNR", "func_mean_fd"),  # médiane par espèce
        # tirée de df_qc_by_label ; ignoré si pas de table QC.
    SPECIES_COVARIATES={},        # covariables espèce additionnelles, ex.
        # {"brain_volume_mm3": {"Human": ..., "Mouse": ...}} -- vide = rien.  # figure kcore vs uniform_median
        # sur la sélection round 2 -- Python seul, pas d'AFNI, peu coûteux.
)


# ═══════════════════════════════════════════════════════════════════════════
# Small helpers
# ═══════════════════════════════════════════════════════════════════════════

class _Tee:
    """Mirror stdout into 00_log.txt. The debug prints are the main way to
    verify a step did what it claims, so they belong in the saved output,
    not only on screen."""

    def __init__(self, path):
        self.file = open(path, "w", encoding="utf-8")
        self.stdout = sys.stdout

    def write(self, data):
        self.stdout.write(data)
        self.file.write(data)

    def flush(self):
        self.stdout.flush()
        self.file.flush()

    def close(self):
        self.file.close()


def _write_legend(path, title, what, how, reading):
    """Write the .txt legend that accompanies each figure."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{title}\n{'=' * len(title)}\n\n")
        f.write(f"WHAT IS PLOTTED\n---------------\n{what}\n\n")
        f.write(f"HOW IT WAS COMPUTED\n-------------------\n{how}\n\n")
        f.write(f"HOW TO READ IT\n--------------\n{reading}\n")


def _tree(root):
    """Create the fixed output tree and return the paths as a dict."""
    p = dict(
        root=root,
        all_units=os.path.join(root, "01_all_units"),
        selection=os.path.join(root, "02_selection"),
        group_maps=os.path.join(root, "03_group_maps"),
        group_maps_final=os.path.join(root, "03b_group_maps_final"),
        jpg=os.path.join(root, "04_jpg"),
        jpg_final=os.path.join(root, "04b_jpg_final"),
        surfaces=os.path.join(root, "05_surfaces"),
        surfaces_final=os.path.join(root, "05b_surfaces_final"),
        figures=os.path.join(root, "06_figures"),
    )
    for v in p.values():
        os.makedirs(v, exist_ok=True)
    return p


def _norm_sentence(mode):
    if mode is True or mode == "zscore":
        return ("z-scored across regions (profile mean subtracted, divided by its SD) "
                "before combining. This equalises amplitude across datasets but re-centres "
                "every profile on zero, so the sign of a consensus value is relative to "
                "the profile mean, not the sign of the underlying correlation")
    if mode == "scale":
        return ("divided by their root-mean-square across regions before combining, "
                "without centring. This equalises amplitude across datasets, so that a "
                "high-amplitude dataset cannot dictate the consensus region by region, "
                "while preserving the sign of the underlying correlations")
    return ("combined on their raw scale, without normalisation. The sign of the "
            "correlations is preserved, but datasets with larger amplitude weigh more "
            "on the region-wise median")


def _write_methods(cfg, paths, bids_labels, n_common, unit_counts, final_round=None,
                   phylo_source=None, n_species=None):
    """Write methods.txt: what was actually run, in prose that can be
    lifted into a manuscript, plus the caveats that apply to it. Called
    once early (partial) and once at the end of the run (complete)."""
    unit = cfg["UNIT"]
    if unit_counts is not None and not unit_counts.empty:
        n_units_line = (f"{int(unit_counts['n_units'].sum())} {unit}s in total "
                        f"({int(unit_counts['n_units'].min())} to "
                        f"{int(unit_counts['n_units'].max())} per dataset)")
    else:
        n_units_line = "n/a"

    if unit == "run":
        unit_par = (f"Each run was treated as one sample ({n_units_line}). Runs from the "
                    f"same animal are not independent, so group maps (3dttest++, one "
                    f"input per run) and run-level statistics count runs, not animals, "
                    f"as degrees of freedom; the number of runs per animal and per "
                    f"dataset varies by more than an order of magnitude, which weights "
                    f"heavily scanned animals and datasets more.")
    else:
        unit_par = (f"Every step treats one {unit} as one sample ({n_units_line}). Runs "
                    f"were averaged within {unit} before any analysis: by median across "
                    f"regional profiles for the profile-based steps, and by AFNI 3dMean "
                    f"across Fisher-z volumes for the group maps, so that repeated scans "
                    f"of one animal are not treated as independent observations.")

    thr = cfg["STAT_MAP_USE_THRESHOLDED"]
    maps_par = (f"Group maps were computed with AFNI 3dttest++ (one Fisher-z map per "
                f"{unit}). Regional group values were taken from the "
                f"{'THRESHOLDED (|z| >= 1.96, other voxels set to 0 and excluded)' if thr else 'unthresholded'} "
                f"group mean Fisher-z map; within each region the median over in-mask "
                f"voxels was used.")

    chain = [f"Round 1: {unit}s were ranked by {cfg['SIMILARITY']} similarity to the "
             f"cross-species consensus, restricted to the regions most invariant across "
             f"species (median absolute deviation across species; "
             f"{cfg['INVARIANCE_MODE']}={cfg['INVARIANCE_VALUE']}; the same region set is "
             f"used in every round), and the top {cfg['N_SELECT']} per dataset kept."]
    if cfg["RUN_ROUND2_SELECTION"]:
        chain.append("Round 2 (species-anchored): each dataset re-ranked its full pool "
                     "against the reference of its own round-1 kept runs.")
    if final_round is not None and final_round >= 3:
        chain.append(f"Rounds 3-{final_round}: the cross-species consensus was rebuilt "
                     f"from the previous round's kept references and the full pool "
                     f"re-ranked; round-to-round stability is reported (F1e).")
    fr = final_round if final_round is not None else "n/a"
    chain.append(f"The final selection is round {fr}. A fixed number (not a fraction) "
                 f"of {unit}s was kept so that selected groups are comparable across "
                 f"datasets regardless of size.")

    tree_txt = (f"Divergence times were taken from the time-calibrated tree "
                f"{os.path.basename(cfg['TREE_PATH'])} (MRCA age of each species pair; "
                f"tips: {', '.join(f'{k}={v}' for k, v in cfg['SPECIES_TO_TREE_TIP'].items())}); "
                f"cite its source (e.g. TimeTree version and access date). "
                if cfg["TREE_PATH"] else
                "Divergence times were APPROXIMATE hand-entered values (TREE_PATH=None); "
                "replace them with a time-calibrated tree before publication. ")
    phylo_par = (tree_txt + f"Phylogenetic analyses used {phylo_source or 'n/a'}"
                 f"{f' ({n_species} species)' if n_species else ''}. Primary test: exact "
                 f"Mantel test between pairwise functional distance (mean 1 - r over "
                 f"cross-species dataset pairs) and pairwise MRCA age from a "
                 f"time-calibrated tree (all permutations of species labels). "
                 f"Complementary: Blomberg's K per region (exact permutation p, "
                 f"Benjamini-Hochberg across regions) and PGLS of each dataset's "
                 f"leave-one-species-out distance to consensus, with Pagel's lambda "
                 f"estimated by maximum likelihood and also fixed at 1 (Brownian "
                 f"motion), alone and against covariates (log number of {unit}s"
                 f"{', and, when QC tables are given, ' + ', '.join(cfg['PGLS_QC_COVARIATES']) if cfg['PGLS_QC_COVARIATES'] else ''}"
                 f"{', ' + ', '.join(cfg['SPECIES_COVARIATES']) if cfg['SPECIES_COVARIATES'] else ''}"
                 f"). The 1-D regression of distance to consensus on "
                 f"divergence from human (F6) was tested by permuting divergence times "
                 f"({cfg['N_PERM_PHYLOGENY']} shuffles).")
    if cfg["PHYLOGENY_USE_FINAL_ROUND"] and phylo_source and "KEPT" in phylo_source:
        phylo_par += (" CAVEAT: these are kept (selected) maps; selection pressure grows "
                      "with pool size, which is confounded with clade, so this is a "
                      "sensitivity analysis, not the primary test.")

    partial = ("" if final_round is not None else
               "\n*** INCOMPLETE -- written before selection; a complete run overwrites\n"
               "*** this file at the end. If you are reading this, the run did not\n"
               "*** finish: do not use it for a manuscript. ***\n")
    txt = f"""METHODS -- Figure 5
==================={partial}

Data
----
{len(bids_labels)} datasets: {', '.join(bids_labels)}
Seed: {cfg['SEED']}
Atlas level: {cfg['ATLAS_LEVEL']}
Common regions across all datasets: {n_common}
Unit of analysis: {unit} -- {n_units_line}

Unit of analysis
----------------
{unit_par}

Regional profiles and group maps
--------------------------------
For each run, the seed-based correlation map was summarised into a regional
profile: the median correlation ({cfg['VALUE_COL']}) within each atlas region
at level {cfg['ATLAS_LEVEL']}. {maps_par} The set of regions common to every
dataset (regions sampled by every dataset's all-{unit}s group map) was fixed
ONCE ({n_common} regions) and used unchanged for every subsequent step.

Cross-species consensus
-----------------------
Each dataset contributed one profile (source: {cfg['CONSENSUS_SOURCE']}).
Datasets of the same species were first collapsed by median, so a species
represented by several datasets counts once. Profiles were
{_norm_sentence(cfg['ZSCORE'])}. The consensus is the median across species.
{'Leave-one-species-out: for each dataset, the consensus it was scored against excluded its entire species, so no dataset contributed to the criterion that selects its own data.' if cfg['LEAVE_ONE_OUT'] else 'No leave-one-out: every dataset contributed to the consensus it was scored against (circular; use only as a sensitivity check).'}

Selection
---------
{chr(10).join(chain)}

Statistical comparisons
-----------------------
The separation between kept and not-kept {unit}s on the selection score is
reported descriptively (difference of medians, AUC; not-kept group subsampled
to the kept group's size, {cfg['N_DRAWS_BALANCED']} draws), without p-values:
kept {unit}s were chosen on that score, so a test would be circular.
QC metrics were related to distance to the dataset reference within dataset
only: at session level for metrics computed once per session (distance
averaged over that session's runs, one observation per session), at run
level for metrics that vary between runs. Correlations were combined across
datasets by fixed-effect meta-analysis of Fisher-z values (weights n - 3),
with Benjamini-Hochberg correction across metrics.

Phylogeny
---------
{phylo_par}

Parameters (verbatim)
---------------------
"""
    for k, v in cfg.items():
        txt += f"  {k} = {v!r}\n"

    with open(os.path.join(paths["root"], "methods.txt"), "w", encoding="utf-8") as f:
        f.write(txt)
    print(f"  [saved] methods.txt")


def _vec_to_profile(common, vec, value_col):
    """Wrap a consensus vector (common-region order) as a profile
    DataFrame, so it can go anywhere a profile is expected."""
    return pd.DataFrame({"base_name": list(common),
                         value_col: np.asarray(vec, dtype=float)})


LEGACY_CONFIG_KEYS = {"PHYLOGENY_USE_ROUND2": "PHYLOGENY_USE_FINAL_ROUND"}


def _validate_config(cfg, user_config):
    """
    Fails loudly on what used to fail silently. Confirmed real cases:
    N_ITERATIVE_ROUNDS=4 given without RUN_ITERATIVE_SELECTION=True did
    nothing at all; a misspelled key was simply ignored. Returns the list
    of warnings (also printed at the top of 00_log.txt).
    """
    warns = []
    user_config = user_config or {}
    for old, new in LEGACY_CONFIG_KEYS.items():
        if old in user_config:
            if new not in user_config:
                cfg[new] = user_config[old]
            cfg.pop(old, None)
            warns.append(f"{old} is renamed {new} (now means the LAST selection "
                         f"round, not round 2) -- value {user_config[old]!r} applied to {new}")
    unknown = sorted(set(user_config) - set(CONFIG) - set(LEGACY_CONFIG_KEYS))
    if unknown:
        raise KeyError(f"Unknown CONFIG key(s): {unknown} -- typo? Valid keys: "
                       f"{sorted(CONFIG)}")
    if cfg["ZSCORE"] not in fa.VALID_NORMALIZATIONS:
        raise ValueError(f"ZSCORE must be one of {fa.VALID_NORMALIZATIONS}, "
                         f"got {cfg['ZSCORE']!r}")
    if "N_ITERATIVE_ROUNDS" in user_config and not cfg["RUN_ITERATIVE_SELECTION"]:
        warns.append(f"N_ITERATIVE_ROUNDS={cfg['N_ITERATIVE_ROUNDS']} has NO effect: "
                     f"RUN_ITERATIVE_SELECTION is False")
    start = 2 if cfg["RUN_ROUND2_SELECTION"] else 1
    if cfg["RUN_ITERATIVE_SELECTION"] and cfg["N_ITERATIVE_ROUNDS"] <= start:
        warns.append(f"RUN_ITERATIVE_SELECTION=True but N_ITERATIVE_ROUNDS="
                     f"{cfg['N_ITERATIVE_ROUNDS']} <= {start} (last round already reached "
                     f"without the loop) -- the loop adds nothing")
    has_refinement = cfg["RUN_ROUND2_SELECTION"] or (
        cfg["RUN_ITERATIVE_SELECTION"] and cfg["N_ITERATIVE_ROUNDS"] > 1)
    if cfg["PHYLOGENY_USE_FINAL_ROUND"] and not has_refinement:
        warns.append("PHYLOGENY_USE_FINAL_ROUND=True but no round after round 1 is "
                     "enabled -- phylogeny stays on the all-units maps")
    if cfg["PHYLOGENY_USE_FINAL_ROUND"] and not cfg["LEAVE_ONE_OUT"]:
        warns.append("PHYLOGENY_USE_FINAL_ROUND=True with LEAVE_ONE_OUT=False: the maps "
                     "feeding F6/F11/F12/PGLS were selected partly for resembling a "
                     "consensus that CONTAINED their own species -- circular; report "
                     "as a sensitivity analysis only, never as the main phylogenetic test")
    if cfg["PHYLOGENY_USE_FINAL_ROUND"]:
        warns.append("PHYLOGENY_USE_FINAL_ROUND=True: kept maps carry a pool-size "
                     "(winner's-curse) bias confounded with clade -- keep the all-units "
                     "result as the primary phylogenetic test")
    return warns


def _map_suffix(cfg, core):
    """AFNI output suffix: <seed>_<core>_lvl<L>[_<MAP_TAG>]."""
    base = f"{cfg['SEED']}_{core}_lvl{cfg['ATLAS_LEVEL']}"
    return f"{base}_{cfg['MAP_TAG']}" if cfg.get("MAP_TAG") else base


def _norm_id(x):
    """'01' and '1' are the same session/run."""
    return str(x).lstrip("0") or "0"


def _map_input_runs(map_path):
    """(subject, session, run) of every input that went into -setA, read
    from the input_runs/ folder _3dttest_EDNiX writes next to the map
    (files '<NN>_sub-<s>_ses-<ses>_run-<r>.nii.gz'). None if absent."""
    d = os.path.join(os.path.dirname(map_path), "input_runs")
    if not os.path.isdir(d):
        return None
    out = set()
    for f in os.listdir(d):
        m = re.match(r"^\d+_sub-(.+?)_ses-(.+?)_run-(.+?)\.nii\.gz$", f)
        if m:
            out.add((m.group(1), _norm_id(m.group(2)), _norm_id(m.group(3))))
    return out


def _check_map_inputs(map_path, expected_runs):
    """
    Provenance check before REUSING a map (OVERWRITE_MAPS=False): the set
    of runs actually fed to 3dttest++ must equal the runs this pipeline
    run expects. Catches a partial copy, a map written by another run
    launched concurrently on the same dataset, and a map built from an
    earlier, different selection -- all of which used to be reused
    silently. Run "NA" (no run token) matches on (subject, session).
    Returns (ok, reason).
    """
    actual = _map_input_runs(map_path)
    if actual is None:
        return False, "no input_runs/ folder next to the map -- provenance unverifiable"
    exact, by_ses = set(), set()
    for rk in expected_runs:
        sub, ses, run = str(rk[0]), _norm_id(rk[1]), str(rk[2])
        (by_ses.add((sub, ses)) if run == "NA" else exact.add((sub, ses, _norm_id(run))))
    missing = [k for k in exact if k not in actual] + \
              [k for k in by_ses if not any(a[:2] == k for a in actual)]
    extra = [a for a in actual if a not in exact and a[:2] not in by_ses]
    if extra:
        return False, (f"{len(extra)} input run(s) this run did not ask for "
                       f"(e.g. {sorted(extra)[:2]}) -- map from another run/selection")
    if missing:
        return False, (f"{len(missing)}/{len(exact) + len(by_ses)} expected run(s) absent "
                       f"(e.g. {sorted(missing)[:2]}) -- partial or different map")
    return True, f"{len(actual)} input run(s) verified against this run's expectation"


def _expected_balanced_subsample(df_class, unit, random_seed=42):
    """
    The (kept_runs_used, outlier_runs_used) that sri.rerun_balanced_by_run
    WOULD draw for this df_class -- a line-by-line mirror of its balancing
    (string ids, units sorted, n = min(#kept units, #outlier units),
    np.random.default_rng(random_seed), kept subsampled first, then
    outliers). Deterministic, so an existing map is reusable only if its
    saved subsample equals this exactly. If the SBA module's balancing ever
    changes, the mismatch makes maps look stale (recomputed) -- never the
    reverse, so the failure mode is extra compute, not a wrong map.
    """
    df = df_class.copy()
    for c in ("subject", "session", "run"):
        df[c] = df[c].astype(str)
    kept_runs = [tuple(r) for r in df.loc[~df["is_outlier"], ["subject", "session", "run"]].values]
    outlier_runs = [tuple(r) for r in df.loc[df["is_outlier"], ["subject", "session", "run"]].values]
    unit_of = {"subject": lambda t: (t[0],), "session": lambda t: (t[0], t[1]),
               "run": lambda t: t}[unit]
    kept_units = sorted({unit_of(t) for t in kept_runs})
    outlier_units = sorted({unit_of(t) for t in outlier_runs})
    n = min(len(kept_units), len(outlier_units))
    rng = np.random.default_rng(random_seed)

    def _subsample(items, n):
        if len(items) <= n:
            return set(items)
        idx = rng.choice(len(items), size=n, replace=False)
        return {items[i] for i in idx}

    kept_sel = _subsample(kept_units, n)
    outlier_sel = _subsample(outlier_units, n)
    return ({t for t in kept_runs if unit_of(t) in kept_sel},
            {t for t in outlier_runs if unit_of(t) in outlier_sel})


def _kept_pair_reusable(path_kept, path_not, subsample_path, df_class, unit):
    """Kept/not-kept maps are reusable only if (1) the saved balanced
    subsample is EXACTLY what the current selection would draw (a subset
    test is not enough: going from N_SELECT=5 to 6 keeps the old 5 kept
    runs inside the new kept set, so a subset test could accept 5-vs-5
    maps for a 6-vs-6 design) and (2) each map's input_runs/ equals that
    subsample. End-to-end check (fake AFNI): N_SELECT 5 -> 6 marks all 18
    kept/not-kept maps stale and reuses the 9 all-units maps."""
    if not os.path.exists(subsample_path):
        return False, "no balance record (.json) -- cannot verify which runs were used"
    with open(subsample_path) as f:
        saved = json.load(f)
    norm = lambda t: (str(t[0]), _norm_id(t[1]), _norm_id(t[2]))
    exp_kept, exp_out = _expected_balanced_subsample(df_class, unit)
    if ({norm(t) for t in saved["kept_runs_used"]} != {norm(t) for t in exp_kept} or
            {norm(t) for t in saved["outlier_runs_used"]} != {norm(t) for t in exp_out}):
        return False, "saved balanced subsample differs from what the current selection draws"
    ok_k, why_k = _check_map_inputs(path_kept, [tuple(t) for t in saved["kept_runs_used"]])
    ok_n, why_n = _check_map_inputs(path_not, [tuple(t) for t in saved["outlier_runs_used"]])
    if not ok_k:
        return False, f"kept map: {why_k}"
    if not ok_n:
        return False, f"not-kept map: {why_n}"
    return True, "both maps verified: inputs == the balanced subsample the current selection draws"


def _resolve_map(sri, cfg, cfg_b, suffix):
    """
    Resolve one group map's path.

    Confirmed bug this fixes: sri.find_group_map_path globs
    "<seed>*_ttest-stat_fisher_cmap.nii.gz", which ALSO matches
    "<seed>..._thresholded_ttest-stat_fisher_cmap.nii.gz"; after sorting,
    "thresholded" comes before "ttest", so the THRESHOLDED map was used
    for every group-map profile (F0, common regions, F6/F7, F11-F13, F8
    reference, multi-level). Non-significant voxels are 0 there and
    compute_region_profile drops exact zeros, so each region was
    summarized over its significant voxels only -- inflated magnitudes,
    sign bias toward the dominant (positive) clusters, and region coverage
    depending on statistical power, i.e. on the number of runs.
    Here the file is chosen explicitly according to
    cfg["STAT_MAP_USE_THRESHOLDED"] (default False = unthresholded mean).
    """
    bids_name = Path(cfg_b["bids_dir"]).name
    # find_group_map_path prints "multiple stat maps matched ... using
    # <..._thresholded_...>", which is no longer what is used (the choice
    # is made below) and would mislead anyone reading the log -- drop it.
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        path = sri.find_group_map_path(
            cfg["RESULTS_BASE_DIR"], cfg_b["species"], bids_name, cfg["SEED"],
            stat_suffix=cfg["STAT_SUFFIX"], category=cfg["CATEGORY"],
            output_suffix=suffix)
    # Also dropped: "group results dir not found", printed by every pre-check
    # of a map not computed yet (expected); the pipeline prints its own
    # OK/MISSING status after each computation.
    for line in buf.getvalue().splitlines():
        if "multiple stat maps matched" not in line and "group results dir not found" not in line:
            print(line)
    if not path:
        return path, False
    folder = os.path.dirname(path)
    cands = [c for c in glob.glob(os.path.join(folder, f"*{suffix}*_{cfg['STAT_SUFFIX']}.nii.gz"))
             if "_oversampled" not in os.path.basename(c)]
    want_thr = bool(cfg["STAT_MAP_USE_THRESHOLDED"])
    pick = sorted(c for c in cands if ("_thresholded_" in os.path.basename(c)) == want_thr)
    if not pick:
        print(f"    [WARN] {cfg_b['label']}: no {'thresholded' if want_thr else 'UNthresholded'} "
              f"'{cfg['STAT_SUFFIX']}' map for suffix {suffix!r} in {folder} "
              f"(found: {[os.path.basename(c) for c in cands]})")
        return None, False
    if len(pick) > 1:
        print(f"    [WARN] {cfg_b['label']}: {len(pick)} candidate maps, using "
              f"{os.path.basename(pick[0])}")
    return pick[0], os.path.exists(pick[0])


def _delete_old_map_files(existing_path, suffix, label=""):
    """OVERWRITE_MAPS=True support: remove every file next to
    existing_path whose name contains `suffix`, before recomputing.

    3dttest++/wb_command write many files per suffix (HEAD/BRIK, cmap,
    zmap, ClustSim tables, surface pngs) -- calling rerun_balanced_by_run
    again without clearing these first risks a stale file from the OLD
    run (e.g. a surface png) surviving next to new stat maps if the new
    run only rewrites a subset. Only files matching this run's own
    suffix are touched -- nothing else in a shared results directory.
    """
    map_dir = os.path.dirname(existing_path)
    removed = 0
    for f in glob.glob(os.path.join(map_dir, f"*{suffix}*")):
        try:
            if os.path.isdir(f):
                shutil.rmtree(f)
            else:
                os.remove(f)
            removed += 1
        except OSError as e:
            print(f"    [WARN] {label}: could not remove {f}: {e}")
    print(f"    {label}: OVERWRITE_MAPS=True -- removed {removed} old file(s) "
          f"matching '*{suffix}*' in {map_dir}, recomputing")
    return removed


# ═══════════════════════════════════════════════════════════════════════════
# STEP 1 -- group map from ALL units (no selection)
# ═══════════════════════════════════════════════════════════════════════════

def step1_all_units_maps(sri, cfg, paths, bids_configs, members_by_label):
    """Group map per dataset using EVERY unit, via unit-level 3dttest.

    This is the unselected baseline: whatever the selection later does,
    this shows what the dataset looks like with nothing removed.

    Two things are different from the previous version, and both were
    bugs rather than choices:

    - balance=False. rerun_balanced_by_run subsamples both groups to
      min(kept, outlier). The old code marked ONE run as a fake outlier
      to get past the "one group is empty" guard, which set that minimum
      to 1 -- so the map was built from a single run, AFNI skipped it,
      and the folder came out empty for every species while
      summary.xlsx reported the full run count. With balance=False no
      subsampling happens and the empty outlier side is simply not run.

    - every run of every subject goes in, and AFNI averages them within
      subject (unit=cfg["UNIT"]). Nothing is sacrificed to keep a
      kept/outlier API happy.
    """
    failures = {}
    rows, index_rows = [], []
    for cfg_b in bids_configs:
        label = cfg_b["label"]
        members = members_by_label.get(label)
        if not members:
            print(f"\n  [SKIP] {label}: no profiles")
            continue

        run_keys = [rk for rks in members.values() for rk in rks]
        n_units = len(members)
        print(f"\n  {label}: group map from ALL {n_units} {cfg['UNIT']}(s) "
              f"({len(run_keys)} run(s))")

        df_class = pd.DataFrame([
            dict(subject=k[0], session=k[1], run=k[2], seed=k[3], is_outlier=False)
            for k in run_keys])
        suffix = _map_suffix(cfg, f"All{cfg['UNIT'].capitalize()}s")

        existing_path, existing_found = _resolve_map(sri, cfg, cfg_b, suffix)
        reusable = False
        if existing_found and not cfg["OVERWRITE_MAPS"]:
            reusable, why = _check_map_inputs(existing_path, run_keys)
            print(f"    all_{cfg['UNIT']}s on disk: "
                  f"{'REUSED -- ' if reusable else '[STALE] recomputing -- '}{why}")
        if reusable:
            path, found = existing_path, existing_found
        else:
            if existing_found:
                _delete_old_map_files(existing_path, suffix, label=f"{label}/all_{cfg['UNIT']}s")
            try:
                sri.rerun_balanced_by_run(
                    df_class, cfg_b["species"], cfg_b["bids_dir"], seed=cfg["SEED"],
                    results_base_dir=cfg["RESULTS_BASE_DIR"],
                    output_suffix_kept=suffix,
                    output_suffix_outlier=f"{suffix}_EMPTY",   # never launched
                    balance=False, unit=cfg["UNIT"],
                    exclude_lateralized=cfg_b.get("exclude_lateralized", True),
                    render_surface=cfg["RENDER_SURFACE"])
            except TypeError as e:
                raise TypeError(
                    "rerun_balanced_by_run() rejected balance=/unit= -- "
                    "SBA_Region_Profile_Index.py has not been patched. Run "
                    "patch_SBA_unit_level.py against it first."
                ) from e
            except Exception as e:
                msg = str(e)
                print(f"    [ERROR] {label}: {msg[:500]}{' ...' if len(msg) > 500 else ''}")
                failures[label] = msg[:300]
                continue
            path, found = _resolve_map(sri, cfg, cfg_b, suffix)
        if not found:
            failures.setdefault(label, "3dttest++ ran but no map was found on disk")
        print(f"    all_{cfg['UNIT']}s: {'OK  ' if found else 'MISSING '} {path}")
        rows.append(dict(label=label, species=cfg_b["species"],
                         unit=cfg["UNIT"], n_units_in_map=n_units,
                         n_runs_behind_them=len(run_keys),
                         output_suffix=suffix))
        index_rows.append(dict(label=label, species=cfg_b["species"],
                               group="all", output_suffix=suffix,
                               map_path=path or "", found=found))

    if failures:
        report = "\n".join(f"      - {lab}: {why}" for lab, why in sorted(failures.items()))
        msg = (f"STEP 1: no all-{cfg['UNIT']}s group map for {len(failures)} dataset(s):\n"
               f"{report}\n"
               f"    Errors such as a missing 'TTnew.*' file, or files vanishing "
               f"during the copy, mean the intermediate folder "
               f"<bids_dir>/Results/Grp_SBA_3dTTEST_<suffix> was modified by "
               f"ANOTHER process (two runs launched at the same time on the same "
               f"datasets) or by an NFS hiccup. MAP_TAG now makes that folder "
               f"unique per OUTPUT_ROOT; rerun with OVERWRITE_MAPS=True.")
        if cfg["REQUIRE_ALL_MAPS"]:
            raise RuntimeError(msg + "\n    (REQUIRE_ALL_MAPS=True -- stopping before "
                               "any figure is drawn from an incomplete species set)")
        print(f"  [WARN] {msg}")

    df = pd.DataFrame(rows)
    map_index = pd.DataFrame(index_rows)
    if not df.empty:
        df.to_excel(os.path.join(paths["all_units"], "summary.xlsx"), index=False)
        map_index.to_excel(os.path.join(paths["all_units"], "map_index.xlsx"),
                           index=False)
        n_ok = int(map_index["found"].sum())
        print(f"\n  [saved] 01_all_units/summary.xlsx ({len(df)} dataset(s)), "
              f"map_index.xlsx ({n_ok}/{len(map_index)} maps found)")
        if n_ok < len(map_index):
            print(f"  [WARN] {len(map_index) - n_ok} step-1 map(s) missing -- "
                  f"see the 'found' column. A missing map here does NOT "
                  f"invalidate the selection or the consensus (both come from "
                  f"the profiles, not from these maps), but it does mean the "
                  f"F0 figure will be drawn from fewer species.")
    return df, map_index


def step1b_species_profiles_figure(sri, cfg, paths, afni_profiles,
                                   common, region_mask, invariance):
    """The two figures drawn off the all-individuals group maps: F0a with
    one line per BIDS, F0 with one line per species, both with the
    regions that define the cross-species profile shaded.

    Reading the profile back OFF the volumetric group map (rather than
    re-deriving it from the per-unit profiles, which figure F1 already
    does) is the point: it is the only check that the AFNI group
    statistics and the Python profile arithmetic describe the same
    pattern. If F0 and F1 disagree in shape, one of the two is wrong and
    it is worth knowing before either goes in a paper.

    afni_profiles is passed in already computed (the common regions are
    now derived from these same profiles, so they must exist before this
    is called) -- this function no longer reads the maps itself.
    """
    if not afni_profiles:
        print("  [SKIP] F0: no group map could be profiled")
        return None

    # F0a -- one line per BIDS: the reproducibility check. Two datasets of
    # one species come from different labs, scanners and protocols, so
    # their agreement (or lack of it) bounds how much of any cross-species
    # difference can be attributed to acquisition rather than biology.
    png_bids = os.path.join(paths["all_units"], "F0a_bids_profiles_afni.png")
    fu.plot_species_profiles_from_afni(
        afni_profiles, common, region_mask=region_mask, invariance=invariance,
        value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"], group_by="label",
        title=(f"Per-BIDS profiles from the all-{cfg['UNIT']}s AFNI group maps "
               f"(lvl {cfg['ATLAS_LEVEL']}, seed {cfg['SEED']})"),
        output_path=png_bids, analysis_fa=fa)

    png = os.path.join(paths["all_units"], "F0_species_profiles_afni.png")
    fu.plot_species_profiles_from_afni(
        afni_profiles, common, region_mask=region_mask, invariance=invariance,
        value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"], group_by="species",
        title=(f"Species profiles from the all-{cfg['UNIT']}s AFNI group maps "
               f"(lvl {cfg['ATLAS_LEVEL']}, seed {cfg['SEED']})"),
        output_path=png, analysis_fa=fa)

    n_sel = int(np.asarray(region_mask, dtype=bool).sum()) if region_mask is not None else len(common)
    _write_legend(
        os.path.join(paths["all_units"], "F0_species_profiles_afni.txt"),
        "F0/F0a -- Species and per-BIDS profiles from the all-units group maps",
        f"F0a draws one line per BIDS (species colour, line style per dataset) so "
        f"the datasets of one species can be read against each other; F0 collapses "
        f"them into one line per species. "
        f"Top panel: one line per species, over the {len(common)} regions common "
        f"to every dataset, read directly off the step-1 AFNI group maps (every "
        f"{cfg['UNIT']}, no selection). Black = the median across species. The "
        f"green shading and the bold region names mark the "
        f"{n_sel}/{len(common)} regions that define the cross-species profile, "
        f"i.e. the ones the selection was restricted to. Bottom panel: each "
        f"region's invariance across species, green where it was retained.",
        f"Common regions are taken from these group maps themselves: present at "
        f"atlas level {cfg['ATLAS_LEVEL']} in every dataset's atlas AND sampled by "
        f"every group map. "
        f"Each dataset's group map is a 3dttest++ over its {cfg['UNIT']}-level "
        f"Fisher-z maps (runs averaged within {cfg['UNIT']} by AFNI 3dMean). "
        f"The map was then summarised into a regional profile with the same "
        f"{cfg['VALUE_COL']}-per-region routine used for the individual "
        f"profiles, datasets of one species collapsed by median, and each "
        f"profile {fa._norm_label(cfg['ZSCORE']).strip() or 'left unscaled'} before "
        f"combining. Regions are ordered by the consensus value. Invariance and "
        f"the retained set come from the profile-level analysis "
        f"(mode={cfg['INVARIANCE_MODE']}, value={cfg['INVARIANCE_VALUE']}).",
        "Read the green band as the answer to 'which regions is the "
        "cross-species claim actually made on'. Inside it, species lines should "
        "track the black line closely -- that is what 'invariant' means, and it "
        "is the region set the selection scored runs against. Outside it, "
        "divergence between species is expected and is precisely why those "
        "regions were excluded. Compare the shape of this figure with F1: F1 is "
        "computed from the Python profiles, F0 from the AFNI volumes, and they "
        "should agree. If they do not, the discrepancy is in the pipeline, not "
        "in the biology.")
    return png


# ═══════════════════════════════════════════════════════════════════════════
# STEP 3 -- group maps for kept AND not-kept, with resolved-path index
# ═══════════════════════════════════════════════════════════════════════════

def step3_group_maps(sri, cfg, paths, bids_configs, members_by_label, selections):
    """Unit-level 3dttest for the kept group and the not-kept group.

    The selection is made on units, but AFNI needs real BOLD files, so
    each selected unit is expanded back to all of its runs here; the
    averaging over those runs happens inside AFNI (3dMean), and the
    balancing happens over units, not runs.

    Also persists, per label, exactly which runs rerun_balanced_by_run
    actually put in each group (kept_runs_used / outlier_runs_used --
    for "outlier" this is the SIZE-MATCHED SUBSAMPLE that went into the
    not-kept 3dttest++, not every non-kept run) to
    03_group_maps/<label>_balance_subsample.json, so step 4 can tell
    "the not-kept runs that were actually compared against kept" apart
    from "everything else that wasn't selected either way".

    Returns (map_index_df, balance_subsample_by_label).
    """
    cfg_by_label = {c["label"]: c for c in bids_configs}
    rows = []
    balance_subsample_by_label = {}
    for label, keys in selections.items():
        cfg_b = cfg_by_label[label]
        members = members_by_label[label]
        selected_units = {tuple(str(x) for x in k) for k in keys}

        df_class = pd.DataFrame([
            dict(subject=rk[0], session=rk[1], run=rk[2], seed=rk[3],
                 is_outlier=ukey not in selected_units)
            for ukey, run_keys in members.items() for rk in run_keys])

        n_kept_u = len(selected_units)
        n_not_u = len(members) - n_kept_u
        suffix_kept = _map_suffix(cfg, "Kept")
        suffix_not = _map_suffix(cfg, "NotKept")
        print(f"\n  {label}: kept={n_kept_u} {cfg['UNIT']}(s), "
              f"not-kept={n_not_u} {cfg['UNIT']}(s) "
              f"({len(df_class)} run(s) behind them)")

        subsample_path = os.path.join(paths["group_maps"], f"{label}_balance_subsample.json")
        path_kept, found_kept = _resolve_map(sri, cfg, cfg_b, suffix_kept)
        path_not, found_not = _resolve_map(sri, cfg, cfg_b, suffix_not)

        reusable = False
        if found_kept and found_not and not cfg["OVERWRITE_MAPS"]:
            reusable, why = _kept_pair_reusable(path_kept, path_not, subsample_path, df_class, cfg["UNIT"])
            print(f"    kept/not-kept maps on disk: "
                  f"{'REUSED -- ' if reusable else '[STALE] recomputing -- '}{why}")
        if reusable:
            with open(subsample_path) as f:
                saved = json.load(f)
            kept_runs_used = {tuple(t) for t in saved["kept_runs_used"]}
            outlier_runs_used = {tuple(t) for t in saved["outlier_runs_used"]}
        else:
            if found_kept:
                _delete_old_map_files(path_kept, suffix_kept, label=f"{label}/kept")
            if found_not:
                _delete_old_map_files(path_not, suffix_not, label=f"{label}/not_kept")
            try:
                (success_kept, success_outlier,
                 kept_runs_used, outlier_runs_used) = sri.rerun_balanced_by_run(
                    df_class, cfg_b["species"], cfg_b["bids_dir"], seed=cfg["SEED"],
                    results_base_dir=cfg["RESULTS_BASE_DIR"],
                    output_suffix_kept=suffix_kept, output_suffix_outlier=suffix_not,
                    balance=True, unit=cfg["UNIT"],
                    exclude_lateralized=cfg_b.get("exclude_lateralized", True),
                    render_surface=cfg["RENDER_SURFACE"])
            except Exception as e:
                print(f"    [ERROR] 3dttest failed: {e}")
                continue
            with open(subsample_path, "w") as f:
                json.dump(dict(kept_runs_used=sorted(kept_runs_used),
                               outlier_runs_used=sorted(outlier_runs_used)), f, indent=2)
            print(f"    [saved] {os.path.basename(subsample_path)} "
                  f"({len(kept_runs_used)} kept run(s), {len(outlier_runs_used)} "
                  f"not-kept run(s) actually in the 3dttest++)")

        balance_subsample_by_label[label] = dict(
            kept_runs_used=kept_runs_used, outlier_runs_used=outlier_runs_used)

        for group, suffix in (("kept", suffix_kept), ("not_kept", suffix_not)):
            path, found = _resolve_map(sri, cfg, cfg_b, suffix)
            rows.append(dict(label=label, species=cfg_b["species"], group=group,
                             output_suffix=suffix, map_path=path or "",
                             found=found))
            print(f"    {group}: {'OK  ' if found else 'MISSING '} {path}")

    df = pd.DataFrame(rows)
    if not df.empty:
        df.to_excel(os.path.join(paths["group_maps"], "map_index.xlsx"), index=False)
        n_ok = int(df["found"].sum())
        print(f"\n  [saved] 03_group_maps/map_index.xlsx ({n_ok}/{len(df)} maps found)")
        if n_ok < len(df):
            print(f"  [WARN] {len(df) - n_ok} map(s) missing -- see the 'found' column")
    return df, balance_subsample_by_label


# ═══════════════════════════════════════════════════════════════════════════
# STEP 3b -- kept / not-kept maps for the FINAL selection round (2..N).
# Own suffix (Kept<N>/NotKept<N>) and own folder, so round 1's maps and the
# all-units maps are never touched.
# ═══════════════════════════════════════════════════════════════════════════

def step3b_final_group_maps(sri, cfg, paths, bids_configs, members_by_label,
                            selections_final, final_round):
    """
    Same balanced 3dttest++ kept/not-kept logic as step3_group_maps, run
    on the FINAL round's selection (round 2 if only RUN_ROUND2_SELECTION,
    round N_ITERATIVE_ROUNDS if the iterative loop ran). Suffix encodes
    the round number, so maps from different final rounds never collide
    and a map on disk always says which round produced it.
    Groups in the returned index: "kept_final" / "not_kept_final".
    """
    cfg_by_label = {c["label"]: c for c in bids_configs}
    rows = []
    balance_subsample_by_label = {}
    tag = f"R{final_round}"
    for label, keys in selections_final.items():
        cfg_b = cfg_by_label[label]
        members = members_by_label[label]
        selected_units = {tuple(str(x) for x in k) for k in keys}

        df_class = pd.DataFrame([
            dict(subject=rk[0], session=rk[1], run=rk[2], seed=rk[3],
                 is_outlier=ukey not in selected_units)
            for ukey, run_keys in members.items() for rk in run_keys])

        n_kept_u = len(selected_units)
        n_not_u = len(members) - n_kept_u
        suffix_kept = _map_suffix(cfg, f"Kept{final_round}")
        suffix_not = _map_suffix(cfg, f"NotKept{final_round}")
        print(f"\n  {label}: round-{final_round} kept={n_kept_u} {cfg['UNIT']}(s), "
              f"not-kept={n_not_u} {cfg['UNIT']}(s) "
              f"({len(df_class)} run(s) behind them)")

        subsample_path = os.path.join(paths["group_maps_final"],
                                      f"{label}_balance_subsample_{tag}.json")
        path_kept, found_kept = _resolve_map(sri, cfg, cfg_b, suffix_kept)
        path_not, found_not = _resolve_map(sri, cfg, cfg_b, suffix_not)

        reusable = False
        if found_kept and found_not and not cfg["OVERWRITE_MAPS"]:
            reusable, why = _kept_pair_reusable(path_kept, path_not, subsample_path, df_class, cfg["UNIT"])
            print(f"    kept/not-kept maps on disk: "
                  f"{'REUSED -- ' if reusable else '[STALE] recomputing -- '}{why}")
        if reusable:
            with open(subsample_path) as f:
                saved = json.load(f)
            kept_runs_used = {tuple(t) for t in saved["kept_runs_used"]}
            outlier_runs_used = {tuple(t) for t in saved["outlier_runs_used"]}
        else:
            if found_kept:
                _delete_old_map_files(path_kept, suffix_kept, label=f"{label}/kept{tag}")
            if found_not:
                _delete_old_map_files(path_not, suffix_not, label=f"{label}/not_kept{tag}")
            try:
                (success_kept, success_outlier,
                 kept_runs_used, outlier_runs_used) = sri.rerun_balanced_by_run(
                    df_class, cfg_b["species"], cfg_b["bids_dir"], seed=cfg["SEED"],
                    results_base_dir=cfg["RESULTS_BASE_DIR"],
                    output_suffix_kept=suffix_kept, output_suffix_outlier=suffix_not,
                    balance=True, unit=cfg["UNIT"],
                    exclude_lateralized=cfg_b.get("exclude_lateralized", True),
                    render_surface=cfg["RENDER_SURFACE"])
            except Exception as e:
                print(f"    [ERROR] round-{final_round} 3dttest failed: {e}")
                continue
            with open(subsample_path, "w") as f:
                json.dump(dict(kept_runs_used=sorted(kept_runs_used),
                               outlier_runs_used=sorted(outlier_runs_used)), f, indent=2)
            print(f"    [saved] {os.path.basename(subsample_path)}")

        balance_subsample_by_label[label] = dict(
            kept_runs_used=kept_runs_used, outlier_runs_used=outlier_runs_used)

        for group, suffix in (("kept_final", suffix_kept), ("not_kept_final", suffix_not)):
            path, found = _resolve_map(sri, cfg, cfg_b, suffix)
            rows.append(dict(label=label, species=cfg_b["species"], group=group,
                             round=final_round, output_suffix=suffix,
                             map_path=path or "", found=found))
            print(f"    {group} (round {final_round}): "
                  f"{'OK  ' if found else 'MISSING '} {path}")

    df = pd.DataFrame(rows)
    if not df.empty:
        df.to_excel(os.path.join(paths["group_maps_final"], f"map_index_final_{tag}.xlsx"),
                    index=False)
        n_ok = int(df["found"].sum())
        print(f"\n  [saved] 03b_group_maps_final/map_index_final_{tag}.xlsx "
              f"({n_ok}/{len(df)} maps found)")
        if n_ok < len(df):
            print(f"  [WARN] {len(df) - n_ok} final-round map(s) missing -- see 'found'")
    return df, balance_subsample_by_label


def _copy_selection_jpgs(sri, bids_configs, members_by_label, selections,
                         balance_subsample_by_label, jpg_root, step_label="STEP 4"):
    """
    Shared logic behind STEP 4 (round 1) and STEP 4b (final round): copies
    preview jpgs into <jpg_root>/<label>/{kept, outliers, nonselectedoutliers}.
    kept/outliers = the runs actually in the balanced 3dttest++ (from
    balance_subsample_by_label); nonselectedoutliers = everything else
    that was neither kept nor drawn into the balanced not-kept sample.
    """
    cfg_by_label = {c["label"]: c for c in bids_configs}
    for label, keys in selections.items():
        members = members_by_label[label]
        bids_dir = cfg_by_label[label]["bids_dir"]
        all_runs = {rk for run_keys in members.values() for rk in run_keys}
        seed_by_run3 = {tuple(str(x) for x in rk[:3]): rk[3] for rk in all_runs}

        subsample = balance_subsample_by_label.get(label)
        if subsample is None:
            selected_units = {tuple(str(x) for x in k) for k in keys}
            print(f"  [{label}] no {step_label} balance record -- falling back to "
                  f"a plain kept/outliers split")
            df_class = pd.DataFrame([
                dict(subject=rk[0], session=rk[1], run=rk[2], seed=rk[3],
                     is_outlier=ukey not in selected_units)
                for ukey, run_keys in members.items() for rk in run_keys])
            try:
                sri.copy_run_qc_images(df_class, bids_dir,
                                       os.path.join(jpg_root, label))
            except Exception as e:
                print(f"  [WARN] {label}: jpg copy failed: {e}")
            continue

        outlier_runs_used = {tuple(str(x) for x in rk) for rk in subsample["outlier_runs_used"]}
        kept_runs_used = {tuple(str(x) for x in rk) for rk in subsample["kept_runs_used"]}
        all_runs3 = set(seed_by_run3)
        nonselected = all_runs3 - outlier_runs_used - kept_runs_used

        df_used = pd.DataFrame(
            [dict(subject=rk[0], session=rk[1], run=rk[2],
                 seed=seed_by_run3[rk], is_outlier=False) for rk in kept_runs_used] +
            [dict(subject=rk[0], session=rk[1], run=rk[2],
                 seed=seed_by_run3[rk], is_outlier=True) for rk in outlier_runs_used])
        try:
            sri.copy_run_qc_images(df_used, bids_dir, os.path.join(jpg_root, label))
        except Exception as e:
            print(f"  [WARN] {label}: jpg copy (kept/outliers) failed: {e}")

        if nonselected:
            df_rest = pd.DataFrame(
                [dict(subject=rk[0], session=rk[1], run=rk[2],
                     seed=seed_by_run3[rk], is_outlier=True) for rk in nonselected])
            tmp_root = os.path.join(jpg_root, f"__tmp_{label}_nonselected")
            try:
                sri.copy_run_qc_images(df_rest, bids_dir, tmp_root)
                dest = os.path.join(jpg_root, label, "nonselectedoutliers")
                if os.path.exists(dest):
                    shutil.rmtree(dest)
                shutil.move(os.path.join(tmp_root, "outliers"), dest)
            except Exception as e:
                print(f"  [WARN] {label}: jpg copy (nonselectedoutliers) failed: {e}")
            finally:
                if os.path.exists(tmp_root):
                    shutil.rmtree(tmp_root)

        print(f"  [{label}] jpgs: kept={len(kept_runs_used)}, "
             f"outliers={len(outlier_runs_used)}, "
             f"nonselectedoutliers={len(nonselected)}")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 5 -- collect surface PNGs under unambiguous names
# ═══════════════════════════════════════════════════════════════════════════

def step5_collect_surfaces(map_index, paths, cfg):
    """Copy surface .png(s) next to the group maps into 05_surfaces/,
    renamed `<label>__<group>__<stat>.png`.

    Only files whose name contains one of cfg["SURFACE_PNG_KEEP"] are
    copied -- wb_command typically renders several views/thresholds per
    stat, most of which are intermediate; by default this keeps just the
    two oversampled fisher cmap/zmap renders.
    """
    n, n_skipped = 0, 0
    keep = cfg["SURFACE_PNG_KEEP"]
    for row in map_index.itertuples():
        if not row.found or not row.map_path:
            continue
        map_dir = os.path.dirname(row.map_path)
        for png in sorted(glob.glob(os.path.join(map_dir, "*.png"))):
            base = os.path.basename(png)
            if keep and not any(k in base for k in keep):
                n_skipped += 1
                continue
            stat = os.path.splitext(base)[0]
            stat = stat.replace(row.output_suffix, "").strip("_") or "map"
            dest = os.path.join(
                paths["surfaces"], f"{row.label}__{row.group}__{stat}.png")
            shutil.copy2(png, dest)
            n += 1
    print(f"  {n} surface png(s) -> 05_surfaces/ "
          f"(filtered out {n_skipped} not matching {keep})")
    if n == 0:
        print(f"  [INFO] none found -- RENDER_SURFACE was "
              f"{cfg['RENDER_SURFACE']}; surfaces only exist if it was True "
              f"and wb_command succeeded, and their names matched "
              f"SURFACE_PNG_KEEP={keep}")
    return n


# ═══════════════════════════════════════════════════════════════════════════
# STEP 7 -- F6 at several atlas levels, from the SAME already-built
# all-units maps (no AFNI recompute)
# ═══════════════════════════════════════════════════════════════════════════

def step6_multilevel_f6(sri, fa, fu, cfg, paths, bids_configs, all_units_index,
                        members_by_label, levels=(1, 2, 3, 4)):
    """
    Re-summarizes the SAME step-1 all-units cmaps at several atlas
    levels (compute_region_profile just re-groups the same volume by a
    different label hierarchy -- no 3dttest++/ClustSim involved), then
    redoes the leave-one-species-out F6 distance/fit at each level, to
    see whether the phylogenetic signal survives coarser or finer
    segmentation. Returns {level: {"ddf_bids", "ddf_species",
    "test_result", "n_common_regions"}}.
    """
    cfg_by_label = {c["label"]: c for c in bids_configs}
    results = {}
    for level in levels:
        print(f"\n  --- F6 @ atlas_level={level} ---")
        profiles_lvl = {}
        # one map per BIDS: the index may carry several groups (e.g. kept /
        # not-kept) -- a later row would silently overwrite an earlier one
        idx = all_units_index
        if "group" in idx.columns and idx["group"].nunique() > 1:
            keep_group = "kept_final" if "kept_final" in set(idx["group"]) else idx["group"].iloc[0]
            idx = idx[idx["group"] == keep_group]
        for r in idx.itertuples():
            if not r.found or not r.map_path:
                continue
            cb = cfg_by_label.get(r.label)
            if cb is None:
                continue
            try:
                run_keys = []
                if r.label in members_by_label:
                    for rks in members_by_label[r.label].values():
                        run_keys.extend(rks)
                atlas = fu.find_matching_atlas_volume(
                    sri, cb, run_keys, r.map_path, verbose=False,
                    label=f"F6lvl{level}/{r.label}")
                atlas = atlas or sri.pick_representative_atlas_volume(cb["bids_dir"])
                ldf = sri.parse_itk_snap_label_file(cb["label_path"])
                prof = sri.compute_region_profile(r.map_path, atlas, ldf, atlas_level=level)
                profiles_lvl[r.label] = (r.species, prof)
            except Exception as e:
                print(f"    [WARN] {r.label} @ lvl{level}: {e}")

        if not profiles_lvl:
            print(f"    [SKIP] lvl{level}: no profile could be built")
            continue

        common_lvl = fa.get_common_regions_from_profiles(
            {lab: {"map": p} for lab, (sp, p) in profiles_lvl.items()},
            value_col=cfg["VALUE_COL"],
            # ALWAYS True here, unlike the single-level pipeline's
            # COMMON_REGIONS_REQUIRE_NONNAN=False default: at a forced
            # atlas_level, label_df still lists all 276 ITK-SNAP rows
            # regardless of whether a region has any voxels at that
            # specific level (compute_region_profile tries every region
            # at the forced level and returns NaN, not a dropped row, for
            # ones with none there) -- so require_nonnan=False here would
            # count every region as "common" at every level structurally,
            # even ones with no real data at that level, exactly what
            # produced the identical 138-region count across levels 1-4.
            require_nonnan=True, verbose=False)
        print(f"    {len(common_lvl)} common region(s) at level {level}")
        if len(common_lvl) < 3:
            print(f"    [SKIP] lvl{level}: fewer than 3 common regions")
            continue

        ddf_bids_lvl = fa.compute_loo_species_distances(
            profiles_lvl, common_lvl, value_col=cfg["VALUE_COL"],
            zscore=cfg["ZSCORE"], verbose=False)
        if len(ddf_bids_lvl) < 3:
            print(f"    [SKIP] lvl{level}: only {len(ddf_bids_lvl)} usable BIDS")
            continue

        ddf_species_lvl = (ddf_bids_lvl.groupby("species", as_index=False)
                           .agg(distance_to_human=("distance_to_human", "mean"),
                                pearson_r=("pearson_r", "mean")))
        tr_lvl = fa.test_phylogenetic_signal(ddf_species_lvl, exclude_human=False)
        mres_lvl = (fp.run_mantel(profiles_lvl, common_lvl, value_col=cfg["VALUE_COL"],
                                  zscore=cfg["ZSCORE"], label=f" lvl{level}")
                    if cfg["RUN_PHYLO_COMPARATIVE"] else None)
        results[level] = dict(ddf_bids=ddf_bids_lvl, ddf_species=ddf_species_lvl,
                              test_result=tr_lvl, n_common_regions=len(common_lvl),
                              mantel=mres_lvl)
        print(f"    lvl{level}: r={tr_lvl['pearson_r']:.3f}, p={tr_lvl['pearson_p']:.4f}, "
              f"n={tr_lvl['n']} species")

    if not results:
        print("  [SKIP] multilevel F6: no level produced usable results")
        return results

    out_xlsx = os.path.join(paths["group_maps"], "distance_to_consensus_multilevel.xlsx")
    with pd.ExcelWriter(out_xlsx) as xw:
        for level, res in results.items():
            res["ddf_bids"].to_excel(xw, sheet_name=f"lvl{level}_per_bids", index=False)
            res["ddf_species"].to_excel(xw, sheet_name=f"lvl{level}_per_species", index=False)
        mantel_rows = [dict(level=lv, n_regions=r["n_common_regions"],
                            **{f"{meth}_{k}": r["mantel"][f"mantel_{meth}"][k]
                               for meth in ("pearson", "spearman")
                               for k in ("r", "p_greater", "p_two_sided")})
                       for lv, r in results.items() if r.get("mantel")]
        if mantel_rows:
            pd.DataFrame(mantel_rows).to_excel(xw, sheet_name="mantel_by_level", index=False)
    print(f"  [saved] {os.path.basename(out_xlsx)}")

    fa.plot_distance_vs_divergence_multilevel(
        results, output_path=os.path.join(paths["figures"],
                                          "F6_multilevel_segmentation_impact.png"))
    return results


# ═══════════════════════════════════════════════════════════════════════════
# STEP 8 -- reference profile from each BIDS' "kept" (or other) group map,
# for scoring every run's distance to it
# ═══════════════════════════════════════════════════════════════════════════

def step_build_reference_profiles(sri, fu, cfg, bids_configs, map_index, members_by_label,
                                  ref_group="kept"):
    """
    Builds {label: (species, profile_df)} from map_index's `ref_group`
    rows (default "kept"), reusing the same grid-matching atlas picker
    as step 1b/F6 (fu.find_matching_atlas_volume) rather than blindly
    taking the first candidate -- same fix as the Rat grid-mismatch bug.
    """
    cfg_by_label = {c["label"]: c for c in bids_configs}
    out = {}
    ref_rows = map_index[(map_index.group == ref_group) & map_index.found]
    for r in ref_rows.itertuples():
        cb = cfg_by_label.get(r.label)
        if cb is None:
            continue
        try:
            run_keys = []
            if r.label in members_by_label:
                for rks in members_by_label[r.label].values():
                    run_keys.extend(rks)
            atlas = fu.find_matching_atlas_volume(
                sri, cb, run_keys, r.map_path, verbose=True, label=f"ref/{r.label}")
            atlas = atlas or sri.pick_representative_atlas_volume(cb["bids_dir"])
            ldf = sri.parse_itk_snap_label_file(cb["label_path"])
            prof = sri.compute_region_profile(r.map_path, atlas, ldf,
                                              atlas_level=cfg["ATLAS_LEVEL"])
            out[r.label] = (r.species, prof)
            print(f"  [{r.label}] {ref_group} reference profile OK")
        except Exception as e:
            print(f"  [WARN] {r.label}: could not build {ref_group} reference profile: {e}")
    return out


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

def _archive_previous_outputs(paths):
    """
    Move (never delete) the result files of a previous run in the same
    OUTPUT_ROOT into <root>/_previous_runs/<timestamp>/, BEFORE this run
    writes anything. Otherwise a run that stops early (e.g. REQUIRE_ALL_MAPS)
    leaves the previous run's figures, tables, methods.txt and log in place,
    looking current. AFNI maps and their balance records (*.json) stay where
    they are so OVERWRITE_MAPS=False can still reuse them after the
    input_runs/ provenance check. Returns the archive folder or None.
    """
    import datetime
    root = paths["root"]
    cand = []
    for f in ("00_log.txt", "methods.txt"):
        if os.path.isfile(os.path.join(root, f)):
            cand.append(os.path.join(root, f))
    for key in ("figures", "selection"):
        d = paths[key]
        cand += [os.path.join(d, f) for f in os.listdir(d)] if os.path.isdir(d) else []
    for key in ("all_units", "group_maps", "group_maps_final"):
        d = paths[key]
        if os.path.isdir(d):
            cand += [os.path.join(d, f) for f in os.listdir(d)
                     if os.path.isfile(os.path.join(d, f)) and not f.endswith(".json")]
    if not cand:
        return None
    dest_root = os.path.join(root, "_previous_runs",
                             datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
    for src in cand:
        rel = os.path.relpath(src, root)
        dst = os.path.join(dest_root, rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.move(src, dst)
    return dest_root


_DIVERGENCE_DEFAULT = dict(fa.DIVERGENCE_TIME_MYA)   # hand-entered values, at import


def _load_time_tree(cfg, bids_configs):
    """
    Load the time-calibrated Newick tree (cfg["TREE_PATH"]) and make it the
    single source of divergence times: fig5_phylo uses it for every pairwise
    age (Mantel, K, PGLS), and fa.DIVERGENCE_TIME_MYA (F6/F7, multi-level)
    is overwritten with the MRCA age of each species with Human.
    """
    # reset module state first: a tree loaded by an earlier run() in the same
    # Python session must not leak into this one
    fp._TREE_PAIR_AGES, fp._TREE_SOURCE = None, None
    fa.DIVERGENCE_TIME_MYA.clear()
    fa.DIVERGENCE_TIME_MYA.update(_DIVERGENCE_DEFAULT)
    if not cfg["TREE_PATH"]:
        print("  [tree] TREE_PATH=None -- using the hand-entered approximate ages "
              "(fa.DIVERGENCE_TIME_MYA / fp.SPECIES_CLADES); not for publication")
        return
    if not os.path.isfile(cfg["TREE_PATH"]):
        # "è" can be stored as one code point (NFC) or two (NFD) depending on
        # how the file was copied: look for the same name under the other form.
        import unicodedata
        d, name = os.path.split(cfg["TREE_PATH"])
        if os.path.isdir(d or "."):
            want = unicodedata.normalize("NFC", name)
            hit = [f for f in os.listdir(d or ".") if unicodedata.normalize("NFC", f) == want]
            if hit:
                print(f"  [tree] {cfg['TREE_PATH']!r} matched on disk as {hit[0]!r} "
                      f"(different Unicode normalisation of the same name)")
                cfg["TREE_PATH"] = os.path.join(d, hit[0])
    if not os.path.isfile(cfg["TREE_PATH"]):
        raise FileNotFoundError(
            f"TREE_PATH {cfg['TREE_PATH']!r} not found. Point TREE_PATH to the "
            f"time-calibrated Newick file, or set TREE_PATH=None to fall back to "
            f"approximate hand-entered ages.")
    species = sorted({c["species"] for c in bids_configs})
    unmapped = [sp for sp in species if sp not in cfg["SPECIES_TO_TREE_TIP"]]
    if unmapped:
        raise KeyError(f"SPECIES_TO_TREE_TIP has no tree tip for {unmapped}")
    mapping = {sp: cfg["SPECIES_TO_TREE_TIP"][sp] for sp in species}
    fp.use_newick_tree(cfg["TREE_PATH"], mapping)
    old = dict(fa.DIVERGENCE_TIME_MYA)
    for sp in species:
        fa.DIVERGENCE_TIME_MYA[sp] = 0.0 if sp == "Human" else fp.mrca_age("Human", sp)
    for sp in species:
        print(f"  [tree] {sp:11s} ({mapping[sp]}): divergence from Human = "
              f"{fa.DIVERGENCE_TIME_MYA[sp]:7.3f} Mya (hand-entered value was "
              f"{old.get(sp, float('nan'))})")


def run(sri, bids_configs, config=None, df_qc_by_label=None):
    """Run the whole Figure 5 pipeline.

    sri          : SBA_Region_Profile_Index module -- the only external
                   dependency (AFNI / nibabel / atlas work). Must be the
                   version patched by patch_SBA_unit_level.py.
    bids_configs : list of per-dataset dicts (label, species, bids_dir,
                   label_path, seeds, optional exclude_lateralized).
    config       : dict overriding CONFIG entries
    df_qc_by_label : {label: df_qc slice} to enable the external QC check
    """
    cfg = dict(CONFIG)
    if config:
        cfg.update(config)
    config_warnings = _validate_config(cfg, config)
    if cfg["UNIT"] not in fu.VALID_UNITS:
        raise ValueError(f"UNIT must be one of {fu.VALID_UNITS}, got {cfg['UNIT']!r}")

    paths = _tree(cfg["OUTPUT_ROOT"])
    if cfg["MAP_TAG"] is None:
        import hashlib
        cfg["MAP_TAG"] = "r" + hashlib.md5(
            os.path.abspath(cfg["OUTPUT_ROOT"]).encode()).hexdigest()[:6]
    if cfg["RESULTS_BASE_DIR"] is None:
        cfg["RESULTS_BASE_DIR"] = paths["group_maps"]
    # AFNI/singularity fail with "stale file handle" if the process cwd is
    # deleted under them (seen in the real log): start from a directory this
    # run owns. All pipeline paths are absolute, so this changes nothing else.
    archived = _archive_previous_outputs(paths)
    old_cwd = os.getcwd()
    os.chdir(paths["root"])
    tee = _Tee(os.path.join(paths["root"], "00_log.txt"))
    sys.stdout = tee
    if archived:
        print(f"  Previous run's figures/tables/log moved to {archived} "
              f"(AFNI maps and *.json balance records left in place for reuse)")
    print(f"  AFNI statmaps stored under RESULTS_BASE_DIR={cfg['RESULTS_BASE_DIR']!r}, "
          f"suffix tag MAP_TAG={cfg['MAP_TAG']!r}, "
          f"{'THRESHOLDED' if cfg['STAT_MAP_USE_THRESHOLDED'] else 'unthresholded'} "
          f"group maps")
    for w in config_warnings:
        print(f"  [CONFIG WARNING] {w}")
    _load_time_tree(cfg, bids_configs)
    try:
        return _run(sri, cfg, paths, bids_configs, df_qc_by_label)
    finally:
        sys.stdout = tee.stdout
        tee.close()
        try:
            os.chdir(old_cwd)
        except OSError:
            pass
        print(f"\n[log] {os.path.join(paths['root'], '00_log.txt')}")


def _run(sri, cfg, paths, bids_configs, df_qc_by_label):

    unit = cfg["UNIT"]
    print(f"\n{'#' * 70}\n  FIGURE 5 PIPELINE  (atlas level {cfg['ATLAS_LEVEL']}, "
          f"seed {cfg['SEED']}, unit of analysis: {unit})\n{'#' * 70}")

    # ---- profiles: runs, then collapsed to units ------------------------
    print(f"\n{'=' * 70}\n  Per-run profiles, aggregated to {unit} level\n{'=' * 70}")
    profiles_by_label, members_by_label, runs_by_label = {}, {}, {}
    run_qc_by_label = {}
    for c in bids_configs:
        print(f"\n  --- {c['label']} ---")
        try:
            label_df = sri.parse_itk_snap_label_file(c["label_path"])
            profs, _ = sri.compute_all_run_profiles(
                c["species"], c["bids_dir"], label_df, cfg["CLASSIFICATION_PATH"],
                category=cfg["CATEGORY"], seeds=[cfg["SEED"]],
                atlas_level=cfg["ATLAS_LEVEL"])
            print(f"  [{c['label']}] compute_all_run_profiles: {len(profs)} raw run(s)")
            if not profs:
                print(f"  [SKIP] {c['label']}: no profiles")
                continue

            if cfg["RUN_QC_ENABLE"]:
                profs, df_qc_runs = fu.reject_bad_runs(
                    sri, fa, profs, method=cfg["REFERENCE_METHOD"], k=cfg["K_REF"],
                    value_col=cfg["VALUE_COL"], metric=cfg["METRIC"],
                    mad_threshold=cfg["RUN_QC_MAD_THRESHOLD"], label=c["label"])
                run_qc_by_label[c["label"]] = df_qc_runs
                if not profs:
                    print(f"  [SKIP] {c['label']}: 0 run survived run-level QC")
                    continue
            else:
                print(f"  [{c['label']}] RUN_QC_ENABLE=False -- run-level QC skipped")

            unit_profs, members = fu.aggregate_profiles_to_unit(
                sri, profs, unit=unit, label=c["label"])
            print(f"  [{c['label']}] aggregate_profiles_to_unit: "
                  f"{len(profs)} run(s) -> {len(unit_profs)} {unit}(s)")
            if not unit_profs:
                print(f"  [SKIP] {c['label']}: aggregation produced nothing")
                continue
            profiles_by_label[c["label"]] = unit_profs
            members_by_label[c["label"]] = members
            runs_by_label[c["label"]] = profs
        except Exception as e:
            print(f"  [ERROR] {c['label']}: {e}")

    if run_qc_by_label:
        qc_frames = [df.assign(label=lab) for lab, df in run_qc_by_label.items() if not df.empty]
        if qc_frames:
            qc_all = pd.concat(qc_frames, ignore_index=True)
            qc_all.to_excel(os.path.join(paths["selection"], "run_level_qc.xlsx"), index=False)
            print(f"\n  [saved] 02_selection/run_level_qc.xlsx "
                  f"({int(qc_all['is_outlier'].sum())}/{len(qc_all)} run(s) flagged across "
                  f"{qc_all['label'].nunique()} dataset(s))")

    print(f"\n  Summary: {len(profiles_by_label)}/{len(bids_configs)} dataset(s) "
          f"produced usable {unit}-level profiles: {sorted(profiles_by_label)}")
    if not profiles_by_label:
        raise ValueError("No dataset produced profiles.")

    unit_counts, _ = fu.check_unit_counts(profiles_by_label, cfg["N_SELECT"], unit=unit)
    unit_counts.to_excel(os.path.join(paths["selection"], "unit_counts.xlsx"),
                         index=False)

    species_of_label = {c["label"]: c["species"] for c in bids_configs}

    # ---- QC tables: keep only rows that belong to each dataset -----------
    # Confirmed real case: the launcher selected QC rows with
    # str.contains(bids_name), and "BIDS_Grandjean" (Rat) is a substring of
    # "BIDS_Grandjean2" (Mouse), so Rat's QC table also held every Mouse
    # row. The subject join was unaffected (ids differ) but the per-species
    # PGLS covariates (median tSNR/FD) mixed Rat and Mouse. Rows whose
    # subject is not one of this dataset's subjects are dropped here.
    if df_qc_by_label:
        clean = {}
        for label, dq in df_qc_by_label.items():
            if (label not in runs_by_label or dq is None or dq.empty
                    or "subject" not in dq.columns):
                clean[label] = dq
                continue
            subs = {str(k[0]).strip() for k in runs_by_label[label]}
            keep = dq["subject"].astype(str).str.strip().isin(subs)
            if (~keep).any():
                ex = sorted(set(dq.loc[~keep, "subject"].astype(str)))[:3]
                print(f"  [QC] {label}: dropped {int((~keep).sum())}/{len(dq)} QC row(s) "
                      f"whose subject has no analysed run in this dataset (e.g. {ex}): "
                      f"either rows of another dataset matched by the launcher's filter, "
                      f"or subjects excluded upstream (e.g. no atlas volume)")
            clean[label] = dq[keep].copy()
        df_qc_by_label = clean
        missing_qc = sorted(set(runs_by_label) - {l for l, d in clean.items()
                                                  if d is not None and not d.empty})
        if missing_qc:
            print(f"  [QC] [WARN] no QC rows for {missing_qc}")

    # ---- STEP 1: group maps from ALL individuals, FIRST ------------------
    # Moved ahead of the common-region and invariance computations, which
    # now derive from these maps (see get_common_regions_from_group_maps):
    # a region counts as common when it actually exists at ATLAS_LEVEL in
    # every dataset's atlas AND is sampled by every dataset's group map.
    # No selection of any kind has happened at this point -- these maps
    # are built from every run that survived the run-level checks.
    all_units_index = pd.DataFrame()
    afni_profiles = {}
    if cfg["RUN_ALL_UNITS_MAPS"]:
        print(f"\n{'=' * 70}\n  STEP 1: group maps from ALL {unit}s\n{'=' * 70}")
        _, all_units_index = step1_all_units_maps(sri, cfg, paths, bids_configs,
                                                  members_by_label)
        print(f"\n{'=' * 70}\n  STEP 1b: profiles of those group maps\n{'=' * 70}")
        afni_profiles = fu.profiles_from_group_maps(
            sri, all_units_index, bids_configs, None,
            atlas_level=cfg["ATLAS_LEVEL"], members_by_label=members_by_label)
        unprofiled = sorted(set(profiles_by_label) - set(afni_profiles))
        if unprofiled and cfg["REQUIRE_ALL_MAPS"]:
            raise RuntimeError(
                f"STEP 1b: the all-{unit}s group map of {unprofiled} exists but could not "
                f"be profiled (see the [WARN] lines above: unreadable/truncated file or "
                f"atlas grid mismatch). REQUIRE_ALL_MAPS=True -- rerun those datasets "
                f"with OVERWRITE_MAPS=True.")
    else:
        print(f"\n  [SKIP] STEP 1 (RUN_ALL_UNITS_MAPS=False)")

    # ---- common regions -------------------------------------------------
    print(f"\n{'=' * 70}\n  Common regions\n{'=' * 70}")
    if afni_profiles:
        common = fa.get_common_regions_from_group_maps(
            afni_profiles, value_col=cfg["VALUE_COL"])
        common_source = "group maps (atlas level + presence in the z-map)"
    else:
        # No group map available (STEP 1 skipped or every map failed):
        # fall back to the per-run profiles so the pipeline still runs.
        common = fa.get_common_regions_from_profiles(
            profiles_by_label, value_col=cfg["VALUE_COL"],
            require_nonnan=cfg["COMMON_REGIONS_REQUIRE_NONNAN"])
        common_source = (f"per-run profiles, fallback "
                         f"(require_nonnan={cfg['COMMON_REGIONS_REQUIRE_NONNAN']})")
    print(f"  {len(common)} common region(s) -- source: {common_source}")
    if not common:
        raise ValueError("No common region survived -- nothing can be compared. "
                         "Check ATLAS_LEVEL and the group maps' coverage.")
    pd.DataFrame({"common_region": common, "source": common_source}).to_excel(
        os.path.join(paths["selection"], "common_regions.xlsx"), index=False)
    _write_methods(cfg, paths, sorted(profiles_by_label), len(common), unit_counts)

    # ---- invariance -----------------------------------------------------
    # With INVARIANCE_VALUE=1.0 (the default) every common region is kept
    # and the mask is all-True: the invariance score is then purely
    # descriptive -- the bottom panel of F0 still ranks the regions, but
    # nothing is excluded from the run selection. Lower the value to
    # restrict the selection to the most invariant fraction.
    print(f"\n{'=' * 70}\n  Cross-species invariant regions\n{'=' * 70}")
    species_profs = {lab: fa.species_profile(
                        sri, profs, source=cfg["CONSENSUS_SOURCE"],
                        method=cfg["REFERENCE_METHOD"], k=cfg["K_REF"],
                        value_col=cfg["VALUE_COL"], metric=cfg["METRIC"])
                     for lab, profs in profiles_by_label.items()}
    invariance, inv_mat, inv_species = fa.compute_region_invariance(
        species_profs, common, value_col=cfg["VALUE_COL"],
        species_of_label=species_of_label, zscore=cfg["ZSCORE"])
    region_mask = fa.select_invariant_regions(
        invariance, common, mode=cfg["INVARIANCE_MODE"], value=cfg["INVARIANCE_VALUE"])
    print(f"  {int(np.asarray(region_mask, dtype=bool).sum())}/{len(common)} "
          f"region(s) used for the run selection "
          f"(mode={cfg['INVARIANCE_MODE']}, value={cfg['INVARIANCE_VALUE']})")
    pd.DataFrame({"region": common, "invariance": invariance,
                  "used_for_selection": region_mask}).to_excel(
        os.path.join(paths["selection"], "region_invariance.xlsx"), index=False)

    # ---- F0 figures: per BIDS, then per species -------------------------
    if afni_profiles and cfg["MAKE_FIGURES"]:
        print(f"\n{'=' * 70}\n  F0 figures (per BIDS, then per species)\n{'=' * 70}")
        step1b_species_profiles_figure(sri, cfg, paths, afni_profiles,
                                       common, region_mask, invariance)

    # ---- STEP 2: selection ---------------------------------------------
    print(f"\n{'=' * 70}\n  STEP 2: consensus selection ({unit} level)\n{'=' * 70}")
    sel = fa.run_consensus_selection(
        sri, cfg, profiles_by_label, common,
        n_select=cfg["N_SELECT"], leave_one_out=cfg["LEAVE_ONE_OUT"],
        consensus_source=cfg["CONSENSUS_SOURCE"], similarity=cfg["SIMILARITY"],
        method=cfg["REFERENCE_METHOD"], k=cfg["K_REF"], zscore=cfg["ZSCORE"],
        threshold_abs=cfg["THRESHOLD_ABS"], species_of_label=species_of_label,
        region_mask=region_mask,
        value_col=cfg["VALUE_COL"], metric=cfg["METRIC"], output_dir=None)
    selections = sel["selections"]

    rows = []
    for label, ranking in sel["rankings"].items():
        chosen = {tuple(str(x) for x in k) for k in selections.get(label, [])}
        members = members_by_label[label]
        for r in ranking.itertuples():
            key = (str(r.subject), str(r.session), str(r.run), str(r.seed))
            rows.append(dict(label=label, species=species_of_label.get(label),
                             unit=unit, subject=r.subject, session=r.session,
                             run=r.run, n_runs_behind=len(members.get(key, [])),
                             score=r.pearson_r, is_kept=key in chosen))
    pd.DataFrame(rows).to_excel(
        os.path.join(paths["selection"], "selection.xlsx"), index=False)
    print(f"  [saved] 02_selection/selection.xlsx")

    # ---- selection chain: round 2 (species-anchored), rounds 3..N ------
    # all_round_selections[r] = {label: [unit keys]} for EVERY executed
    # round; final_round is the last one. The final round alone gets AFNI
    # maps (STEP 3b), JPGs (4b), surfaces (5b), its F4/F1b line, and feeds
    # the phylogeny when PHYLOGENY_USE_FINAL_ROUND=True.
    all_round_selections = {1: selections}
    # final_criterion[label] = profile_df the FINAL round actually scored
    # that BIDS's runs against (F8b plots correlation to exactly this).
    final_criterion = {lab: _vec_to_profile(common, vec, cfg["VALUE_COL"])
                       for lab, vec in sel["consensus_by_label"].items()}
    selections_round2 = {}
    if cfg["RUN_ROUND2_SELECTION"]:
        print(f"\n{'=' * 70}\n  STEP 2b: round-2 selection (species-anchored)\n{'=' * 70}")
        r2 = fa.run_round2_selection(
            sri, cfg, profiles_by_label, selections, common,
            n_select=cfg["N_SELECT"], region_mask=region_mask)
        selections_round2 = r2["selections"]
        for label in set(selections) - set(selections_round2):
            print(f"  [WARN] {label}: no round-2 selection -- keeps its round-1 selection")
            selections_round2[label] = selections[label]
        all_round_selections[2] = selections_round2
        final_criterion = {lab: (_vec_to_profile(common, r2["criterion_by_label"][lab],
                                                 cfg["VALUE_COL"])
                                 if lab in r2["criterion_by_label"] else final_criterion.get(lab))
                           for lab in selections_round2}
        rows_r2 = []
        for label, ranking in r2["rankings"].items():
            chosen = {tuple(str(x) for x in k) for k in selections_round2.get(label, [])}
            members = members_by_label[label]
            for r in ranking.itertuples():
                key = (str(r.subject), str(r.session), str(r.run), str(r.seed))
                rows_r2.append(dict(label=label, species=species_of_label.get(label),
                                    unit=unit, subject=r.subject, session=r.session,
                                    run=r.run, n_runs_behind=len(members.get(key, [])),
                                    score=r.pearson_r, is_kept=key in chosen))
        if rows_r2:
            pd.DataFrame(rows_r2).to_excel(
                os.path.join(paths["selection"], "selection_round2.xlsx"), index=False)
            print(f"  [saved] 02_selection/selection_round2.xlsx")

        if cfg["RUN_REFERENCE_METHOD_COMPARISON"]:
            print(f"\n{'=' * 70}\n  Effect of REFERENCE_METHOD on round-2 selection\n{'=' * 70}")
            overlap_df = fa.compare_reference_methods_on_round2(
                sri, cfg, profiles_by_label, selections, common,
                n_select=cfg["N_SELECT"], region_mask=region_mask)
            if not overlap_df.empty:
                overlap_df.drop(columns=["kcore_only", "uniform_only"]).to_excel(
                    os.path.join(paths["selection"], "reference_method_overlap.xlsx"),
                    index=False)
                fa.plot_reference_method_overlap(
                    overlap_df, output_path=os.path.join(
                        paths["figures"], "F1d_reference_method_effect.png"))

    start_round = max(all_round_selections)
    if cfg["RUN_ITERATIVE_SELECTION"] and cfg["N_ITERATIVE_ROUNDS"] > start_round:
        print(f"\n{'=' * 70}\n  Iterative selection: rounds {start_round + 1}.."
              f"{cfg['N_ITERATIVE_ROUNDS']} (continuing from round {start_round})\n{'=' * 70}")
        iter_res = fa.run_iterative_selection(
            sri, cfg, profiles_by_label, common, n_select=cfg["N_SELECT"],
            start_selections=all_round_selections[start_round], start_round=start_round,
            n_rounds=cfg["N_ITERATIVE_ROUNDS"], species_of_label=species_of_label,
            region_mask=region_mask)
        for r, rd in iter_res["rounds"].items():
            all_round_selections[r] = rd["selections"]
        last = max(iter_res["rounds"]) if iter_res["rounds"] else None
        if last is not None:
            cons_last = iter_res["rounds"][last]["consensus_by_label"]
            final_criterion = {lab: (_vec_to_profile(common, cons_last[lab], cfg["VALUE_COL"])
                                     if lab in cons_last else final_criterion.get(lab))
                               for lab in all_round_selections[last]}

    final_round = max(all_round_selections)
    selections_final = all_round_selections[final_round]
    print(f"\n  Selection chain: rounds {sorted(all_round_selections)} -- FINAL round = "
          f"{final_round}" + ("" if final_round > 1 else " (no refinement round ran)"))

    # Convergence across the WHOLE chain (round r vs r-1, r >= 2) + per-round table
    if final_round >= 2:
        conv_rows, rows_iter = [], []
        for r in range(2, final_round + 1):
            for label, keys in all_round_selections[r].items():
                cur = {tuple(str(x) for x in k) for k in keys}
                prev = {tuple(str(x) for x in k)
                        for k in all_round_selections[r - 1].get(label, [])}
                conv_rows.append(dict(round=r, label=label, n_select=len(cur),
                                      n_overlap_with_previous=len(cur & prev)))
        for r, sel_r in all_round_selections.items():
            for label, keys in sel_r.items():
                for k in keys:
                    rows_iter.append(dict(round=r, is_final=r == final_round, label=label,
                                          subject=k[0], session=k[1], run=k[2], seed=k[3]))
        conv_df = pd.DataFrame(conv_rows)
        conv_df.to_excel(os.path.join(paths["selection"], "selection_convergence.xlsx"),
                         index=False)
        pd.DataFrame(rows_iter).to_excel(
            os.path.join(paths["selection"], "selection_by_round.xlsx"), index=False)
        fa.plot_iterative_selection_convergence(
            conv_df, output_path=os.path.join(
                paths["figures"], "F1e_selection_convergence.png"))
        print(f"  [saved] 02_selection/selection_by_round.xlsx, selection_convergence.xlsx")

    balanced = fa.balanced_kept_vs_notkept(
        profiles_by_label, selections, common,
        consensus_by_label=sel["consensus_by_label"], value_col=cfg["VALUE_COL"],
        zscore=cfg["ZSCORE"], n_draws=cfg["N_DRAWS_BALANCED"],
        species_of_label=species_of_label)
    if not balanced.empty:
        balanced.to_excel(os.path.join(paths["selection"], "balanced_stats.xlsx"),
                          index=False)

    if df_qc_by_label:
        try:
            from external_qc_validation import merge_qc_with_classification
            # With unit="subject" the classification's session/run are the
            # placeholder "ALL", so the QC merge will fall back to its
            # subject-only join level -- which is the correct granularity
            # here, not a degraded one. Check the printed join_level.
            qc = fa.validate_selection_against_qc(
                profiles_by_label, selections, df_qc_by_label,
                merge_qc_with_classification, n_draws=cfg["N_DRAWS_BALANCED"])
            if not qc.empty:
                qc.to_excel(os.path.join(paths["selection"], "external_qc.xlsx"),
                            index=False)
        except Exception as e:
            print(f"  [WARN] external QC check failed: {e}")

    # ---- STEP 3/4/5 -----------------------------------------------------
    map_index = pd.DataFrame()
    balance_subsample_by_label = {}
    if cfg["RUN_GROUP_MAPS"]:
        print(f"\n{'=' * 70}\n  STEP 3: group maps (kept / not-kept)\n{'=' * 70}")
        map_index, balance_subsample_by_label = step3_group_maps(
            sri, cfg, paths, bids_configs, members_by_label, selections)
    else:
        print(f"\n  [SKIP] STEP 3 (RUN_GROUP_MAPS=False)")

    map_index_final = pd.DataFrame()
    balance_subsample_final = {}
    if cfg["RUN_GROUP_MAPS"] and final_round >= 2:
        print(f"\n{'=' * 70}\n  STEP 3b: final-round group maps "
              f"(round {final_round}: kept / not-kept)\n{'=' * 70}")
        map_index_final, balance_subsample_final = step3b_final_group_maps(
            sri, cfg, paths, bids_configs, members_by_label, selections_final, final_round)

    if cfg["COPY_JPG"]:
        print(f"\n{'=' * 70}\n  STEP 4: SBA preview jpgs "
              f"(kept / outliers=in the 3dttest++ / nonselectedoutliers=rest)\n"
              f"{'=' * 70}")
        _copy_selection_jpgs(sri, bids_configs, members_by_label, selections,
                             balance_subsample_by_label, paths["jpg"], step_label="STEP 4")

    if cfg["COPY_JPG"] and final_round >= 2 and not map_index_final.empty:
        print(f"\n{'=' * 70}\n  STEP 4b: final-round (round {final_round}) "
              f"SBA preview jpgs\n{'=' * 70}")
        _copy_selection_jpgs(sri, bids_configs, members_by_label, selections_final,
                             balance_subsample_final, paths["jpg_final"],
                             step_label="STEP 4b")

    if not map_index.empty:
        print(f"\n{'=' * 70}\n  STEP 5: collect surface renders\n{'=' * 70}")
        step5_collect_surfaces(map_index, paths, cfg)

    if not map_index_final.empty:
        print(f"\n{'=' * 70}\n  STEP 5b: collect final-round (round {final_round}) "
              f"surface renders\n{'=' * 70}")
        step5_collect_surfaces(map_index_final,
                               {**paths, "surfaces": paths["surfaces_final"]}, cfg)

    # ---- phylogeny source (F6, F7, F11-F13, STEP 7) -----------------------
    # Resolved ONCE here so every phylogenetic output uses the same maps.
    # Confirmed bugs this replaces: the switch used to live inside the
    # figures block (overwriting afni_profiles in place), STEP 7 ignored
    # it entirely, and the log/excel still said "all-units, no selection".
    phylo_profiles = afni_profiles
    phylo_index = all_units_index
    phylo_source = f"all-{unit}s group maps (no selection)"
    # Every dataset that has runs must be represented; compare against that,
    # not against afni_profiles (which can itself be incomplete -- confirmed
    # real case: 6/9 all-units maps, while all 9 kept maps existed, and the
    # old comparison wrongly fell back with an empty "missing" list).
    expected = set(profiles_by_label)
    missing_all_units = sorted(expected - set(afni_profiles))
    if missing_all_units:
        print(f"  [PHYLO] [WARN] no all-{unit}s profile for {missing_all_units}")
    if cfg["PHYLOGENY_USE_FINAL_ROUND"] and final_round >= 2:
        kf = (map_index_final[(map_index_final.group == "kept_final") & map_index_final.found]
              if not map_index_final.empty else pd.DataFrame(columns=["label"]))
        missing = sorted(expected - set(kf["label"]))
        if missing:
            print(f"  [PHYLO] [WARN] PHYLOGENY_USE_FINAL_ROUND=True but no round-"
                  f"{final_round} kept map for {missing} -- keeping all-{unit}s maps "
                  f"(mixing sources across species would bias F6)")
        else:
            cand = fu.profiles_from_group_maps(
                sri, kf, bids_configs, None, atlas_level=cfg["ATLAS_LEVEL"],
                members_by_label=members_by_label)
            not_profiled = sorted(expected - set(cand))
            if not_profiled:
                print(f"  [PHYLO] [WARN] could not profile round-{final_round} kept maps "
                      f"for {not_profiled} -- keeping all-{unit}s maps")
            else:
                phylo_profiles, phylo_index = cand, kf
                phylo_source = f"round-{final_round} KEPT group maps (selection applied)"
    if set(phylo_profiles) != expected:
        print(f"  [PHYLO] [WARN] phylogenetic analyses cover {len(phylo_profiles)}/"
              f"{len(expected)} datasets -- missing {sorted(expected - set(phylo_profiles))}")
    print(f"\n  [PHYLO] phylogenetic analyses use: {phylo_source}")

    # ---- STEP 6: figures + legends --------------------------------------
    if cfg["MAKE_FIGURES"]:
        print(f"\n{'=' * 70}\n  STEP 6: figures\n{'=' * 70}")
        fd = paths["figures"]

        fa.plot_consensus_by_species(
            sel["species_profiles"], common, species_of_label,
            value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"],
            title=f"Consensus by species -- all runs (lvl {cfg['ATLAS_LEVEL']}, {unit} level)",
            output_path=os.path.join(fd, "F1_consensus_by_species.png"))
        _write_legend(
            os.path.join(fd, "F1_consensus_by_species.txt"),
            "F1 -- Consensus pattern by species (all runs)",
            f"One line per species, over the {len(common)} regions common to every "
            f"dataset. Black = the consensus (median across species). Built from "
            f"EVERY run (kept and not-kept alike) -- source=all_runs, no selection.",
            f"Built from the {unit}-level profiles (runs averaged within {unit} by "
            f"median). Datasets of the same species were collapsed by median "
            f"first, so a species with several datasets counts once. Profiles were "
            f"{fa._norm_label(cfg['ZSCORE']).strip() or 'left unscaled'} before combining.",
            "Regions are ordered by the consensus value. Species lines tracking the "
            "black line indicates a shared profile shape; a species departing from "
            "it is atypical for that region. F0 shows the same thing computed from "
            "the AFNI group maps instead -- the two should agree. See F1a/F1b for "
            "the same plot built from round-1 / final-round KEPT runs only instead of "
            "every run.")

        # Kept-round-1 and kept-round-2 per-BIDS references, computed once
        # here and reused below by F1a/F1b AND by F4 (was recomputed inline
        # there before; same build_reference_for_profiles calls, now shared).
        kept_ref_by_label, kept_ref_final_by_label = {}, {}
        for label, profs in profiles_by_label.items():
            kept_keys_r1 = {tuple(str(x) for x in k) for k in selections.get(label, [])}
            kept_profs_r1 = {k: v for k, v in profs.items()
                             if tuple(str(x) for x in k) in kept_keys_r1}
            try:
                if len(kept_profs_r1) >= 2:
                    kept_ref_by_label[label], _ = fa.build_reference_for_profiles(
                        sri, kept_profs_r1, method=cfg["REFERENCE_METHOD"],
                        k=min(cfg["K_REF"], len(kept_profs_r1)),
                        value_col=cfg["VALUE_COL"], metric=cfg["METRIC"])
            except Exception as e:
                print(f"  [WARN] {label}: round-1 kept reference failed: {e}")
            if final_round >= 2:
                kept_keys_r2 = {tuple(str(x) for x in k)
                                for k in selections_final.get(label, [])}
                kept_profs_r2 = {k: v for k, v in profs.items()
                                 if tuple(str(x) for x in k) in kept_keys_r2}
                try:
                    if len(kept_profs_r2) >= 2:
                        kept_ref_final_by_label[label], _ = fa.build_reference_for_profiles(
                            sri, kept_profs_r2, method=cfg["REFERENCE_METHOD"],
                            k=min(cfg["K_REF"], len(kept_profs_r2)),
                            value_col=cfg["VALUE_COL"], metric=cfg["METRIC"])
                except Exception as e:
                    print(f"  [WARN] {label}: round-{final_round} kept reference failed: {e}")

        if kept_ref_by_label:
            fa.plot_consensus_by_species(
                kept_ref_by_label, common, species_of_label,
                value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"],
                title=f"Consensus by species -- round-1 kept only "
                      f"(lvl {cfg['ATLAS_LEVEL']}, {unit} level)",
                output_path=os.path.join(fd, "F1a_consensus_by_species_kept1.png"))
            _write_legend(
                os.path.join(fd, "F1a_consensus_by_species_kept1.txt"),
                "F1a -- Consensus pattern by species (round-1 kept only)",
                "Same as F1, but each BIDS's line is now its round-1 KEPT "
                "reference (cross-species-anchored selection) instead of the "
                "all-runs median.",
                f"{cfg['REFERENCE_METHOD']} reference over each BIDS's round-1 "
                f"kept runs only (k={cfg['K_REF']}, metric={cfg['METRIC']}).",
                "Compare against F1: if a species' line moves a lot from F1 to "
                "F1a, round-1 selection changed what that species looks like "
                "here -- either cleaning it up (closer to consensus) or biasing "
                "it toward the consensus by construction, which F6's use of the "
                "all-units map (not this kept map) by default is meant to avoid.")

        if kept_ref_final_by_label:
            fa.plot_consensus_by_species(
                kept_ref_final_by_label, common, species_of_label,
                value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"],
                title=f"Consensus by species -- final round (round {final_round}) kept only "
                      f"(lvl {cfg['ATLAS_LEVEL']}, {unit} level)",
                output_path=os.path.join(fd, "F1b_consensus_by_species_kept_final.png"))
            _write_legend(
                os.path.join(fd, "F1b_consensus_by_species_kept_final.txt"),
                f"F1b -- Consensus pattern by species (final round = round {final_round})",
                f"Same as F1a, but each BIDS's line is its KEPT reference from the LAST "
                f"selection round (round {final_round}; chain: round 1 cross-species"
                f"{', round 2 species-anchored' if 2 in all_round_selections else ''}"
                f"{f', rounds 3-{final_round} cross-species from previous kept' if final_round >= 3 else ''}).",
                f"{cfg['REFERENCE_METHOD']} reference over each BIDS's round-{final_round} "
                f"kept runs only.",
                "Compare against F1a: if F1b and F1a are nearly identical, the refinement "
                "rounds barely changed that species. F1e shows round-to-round stability.")

        fa.plot_consensus_by_bids(
            sel["species_profiles"], common, species_of_label,
            value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"],
            title=f"Consensus by BIDS -- all runs (lvl {cfg['ATLAS_LEVEL']}, {unit} level)",
            output_path=os.path.join(fd, "F1c_consensus_by_bids.png"))
        _write_legend(
            os.path.join(fd, "F1c_consensus_by_bids.txt"),
            "F1c -- Consensus pattern by BIDS (all runs)",
            f"Companion to F1: one line per BIDS (10 here) instead of one per "
            f"species (7 in F1) -- same black consensus line, unchanged.",
            f"Same {unit}-level profiles and z-scoring as F1, just not collapsed "
            f"by species before plotting.",
            "Lets you see whether the 2-BIDS species (Human, Macaque, Marmoset) "
            "actually agree with each other before F1 collapses them into one "
            "line -- two BIDS of the same species tracking each other closely "
            "but diverging from the consensus points at a real species effect; "
            "two BIDS of the same species diverging from EACH OTHER points at a "
            "dataset/protocol effect instead.")

        fa.plot_intra_vs_inter_species_distance(
            profiles_by_label, common, species_of_label,
            value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"],
            title=f"Within- vs across-species similarity (lvl {cfg['ATLAS_LEVEL']})",
            output_path=os.path.join(fd, "F2_within_vs_across_species.png"))
        _write_legend(
            os.path.join(fd, "F2_within_vs_across_species.txt"),
            "F2 -- Within- vs across-species similarity",
            f"Distribution of {unit}-to-{unit} profile correlation for three kinds of "
            f"pairs: same dataset; same species but another dataset; different species.",
            f"All {unit} pairs (subsampled above 20000 pairs). AUC (same vs different "
            f"species) is descriptive. p comes from permuting species labels ACROSS "
            f"DATASETS (each dataset keeps its {unit}s; pair correlations fixed): pairs "
            f"share {unit}s and are not independent, so a Mann-Whitney p on pairs would "
            f"be anti-conservative by orders of magnitude."
            + (f" With UNIT='run', same-dataset pairs include repeated scans of the same "
               f"animal and the same protocol, so they overestimate species similarity."
               if unit == "run" else ""),
            "The informative comparison is 'same species, other dataset' vs 'different "
            "species': it separates species from dataset/protocol and animal identity. "
            "If it is not clearly above the different-species class, species is not the "
            "dominant structure and phylogenetic interpretation is unsupported. Only "
            "species with several datasets contribute to that middle class, so the "
            "permutation test has few distinct configurations and limited power.")

        selected_by_label = {}
        for label, profs in profiles_by_label.items():
            kept_keys = {tuple(str(x) for x in k) for k in selections_final.get(label, [])}
            selected_by_label[label] = {k: v for k, v in profs.items()
                                        if tuple(str(x) for x in k) in kept_keys}
        fa.plot_within_vs_across_all_vs_selected(
            profiles_by_label, selected_by_label, common, species_of_label,
            value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"],
            title=f"Within- vs across-species similarity: all vs final-round selection "
                  f"(lvl {cfg['ATLAS_LEVEL']})",
            selected_name=f"Selected, round {final_round}",
            output_path=os.path.join(fd, "F2b_within_vs_across_all_vs_selected.png"))
        _write_legend(
            os.path.join(fd, "F2b_within_vs_across_all_vs_selected.txt"),
            "F2b -- Within- vs across-species similarity, all runs vs selected",
            f"Same analysis as F2: left, every {unit}; right, the {unit}s kept by the "
            f"final selection round (round {final_round}).",
            "Same pair classes, AUC and dataset-level permutation test as F2, computed "
            "independently on each pool.",
            "The left panel is the unbiased estimate. Selection moves the right panel "
            "in both directions: rounds scored against a cross-species consensus pull "
            "runs toward a shared pattern and should REDUCE species separation, while a "
            "species-anchored round (round 2) pulls each dataset toward its own mean "
            "and can INFLATE it. A separation seen only on the right is therefore not "
            "evidence of species structure.")

        fa.plot_score_distributions(
            sel["rankings"], selections,
            title=f"Consensus score per dataset (lvl {cfg['ATLAS_LEVEL']}, {unit} level)",
            output_path=os.path.join(fd, "F3_score_distributions.png"))
        _write_legend(
            os.path.join(fd, "F3_score_distributions.txt"),
            "F3 -- Selection score distributions",
            f"Per dataset, the distribution of every {unit}'s similarity to the "
            f"consensus, with the selection cutoff marked in orange.",
            f"Similarity metric: {cfg['SIMILARITY']}. Cutoff is the score of the "
            f"{cfg['N_SELECT']}th best {unit}.",
            f"If the kept {unit}s form a separate upper mode, the selection captures "
            f"a real subgroup. If the distribution is smooth across the cutoff, the "
            f"selection is a cut through a continuum and should be described as such. "
            f"Note the x-axis is now much shorter than at run level -- each dataset "
            f"has far fewer {unit}s than runs, so a single {unit} moves the cutoff "
            f"more.")

        own_refs = {}
        for label, profs in profiles_by_label.items():
            try:
                own_refs[label], _ = fa.build_reference_for_profiles(
                    sri, profs, method=cfg["REFERENCE_METHOD"], k=cfg["K_REF"],
                    value_col=cfg["VALUE_COL"], metric=cfg["METRIC"])
            except Exception as e:
                print(f"  [WARN] {label}: reference failed: {e}")

        if own_refs:
            # cons_prof must stay on the RAW value_col scale to be
            # comparable to own_refs and the individual units. The
            # consensus in sel["consensus_by_label"] is z-scored and
            # leave-one-out (correct for the selection SCORING) --
            # rebuild a separate display-only consensus here: all
            # species, no LOO, no z-score.
            display_consensus, _ = fa.build_consensus_pattern(
                species_profs, common, value_col=cfg["VALUE_COL"],
                exclude_label=None, zscore=False, threshold_abs=cfg["THRESHOLD_ABS"],
                species_of_label=species_of_label, verbose=False)
            finite = np.isfinite(display_consensus)
            print(f"  [F4] display consensus (raw scale): "
                  f"{int(finite.sum())}/{len(display_consensus)} region(s) finite, "
                  f"range=[{np.nanmin(display_consensus):.3f}, "
                  f"{np.nanmax(display_consensus):.3f}]")
            sign_df = fa.diagnose_consensus_sign(
                species_profs, common, value_col=cfg["VALUE_COL"],
                species_of_label=species_of_label)
            sign_df.to_excel(os.path.join(
                paths["group_maps"], "consensus_sign_diagnostic.xlsx"), index=False)

            cons_prof = own_refs[next(iter(own_refs))].copy()
            cons_prof[cfg["VALUE_COL"]] = (
                cons_prof["base_name"].str.lower().str.strip()
                .map(dict(zip(common, display_consensus))))

            # Unweighted grand mean: every one of the 10 BIDS's own
            # species_ref counted once, unlike display_consensus which
            # collapses same-species BIDS together first (Human counts
            # once even though 2 BIDS feed it). Answers "what does the
            # average across every dataset we have look like", flat.
            grand_mean_vals = np.nanmedian(
                np.vstack([fa._ordered_common_values(p, common, cfg["VALUE_COL"])
                          for p in own_refs.values()]), axis=0)
            grand_mean_prof = own_refs[next(iter(own_refs))].copy()
            grand_mean_prof[cfg["VALUE_COL"]] = (
                grand_mean_prof["base_name"].str.lower().str.strip()
                .map(dict(zip(common, grand_mean_vals))))

            per = os.path.join(fd, "F4_per_dataset")
            os.makedirs(per, exist_ok=True)
            for label, profs in profiles_by_label.items():
                if label in own_refs:
                    own_vals = own_refs[label][cfg["VALUE_COL"]].values
                    print(f"  [F4] {label}: own_ref range=[{np.nanmin(own_vals):.3f}, "
                          f"{np.nanmax(own_vals):.3f}]  (should be same order of "
                          f"magnitude as the consensus range above)")

                    kept_keys = {tuple(str(x) for x in k) for k in selections.get(label, [])}
                    notkept_profs = {k: v for k, v in profs.items()
                                     if tuple(str(x) for x in k) not in kept_keys}
                    kept_ref = kept_ref_by_label.get(label)
                    notkept_ref = None
                    try:
                        if len(notkept_profs) >= 2:
                            notkept_ref, _ = fa.build_reference_for_profiles(
                                sri, notkept_profs, method=cfg["REFERENCE_METHOD"],
                                k=min(cfg["K_REF"], len(notkept_profs)),
                                value_col=cfg["VALUE_COL"], metric=cfg["METRIC"])
                    except Exception as e:
                        print(f"  [WARN] {label}: not-kept-only reference failed: {e}")

                    kept_ref_r2 = kept_ref_final_by_label.get(label)

                    fa.plot_bids_ref_with_individuals(
                        profs, own_refs[label], cons_prof, common,
                        selected_keys=selections.get(label),
                        kept_ref=kept_ref, notkept_ref=notkept_ref,
                        grand_mean_ref=grand_mean_prof, kept_ref_round2=kept_ref_r2,
                        actual_selection_target=sel["consensus_by_label"].get(label),
                        final_round_label=(f"Kept final round (round {final_round})"
                                           if final_round >= 2 else "Kept final round"),
                        target_label=(f"Round-1 selection target ("
                                      f"{fa._norm_label(cfg['ZSCORE']).strip() or 'raw'}, "
                                      f"{'LOO' if cfg['LEAVE_ONE_OUT'] else 'no LOO'}, "
                                      f"right axis)"),
                        value_col=cfg["VALUE_COL"],
                        title=f"{label} (lvl {cfg['ATLAS_LEVEL']}, {unit} level)",
                        output_path=os.path.join(per, f"{label}.png"))
            _write_legend(
                os.path.join(fd, "F4_per_dataset.txt"),
                "F4 -- Per-dataset reference, consensus and individual units",
                f"One figure per dataset: its own reference (blue), the consensus it "
                f"was scored against (orange dashed), and every individual {unit} "
                f"(green = kept, grey = not kept).",
                f"Dataset reference built with {cfg['REFERENCE_METHOD']} over ALL its "
                f"{unit}s, independent of the selection. Regions ordered by consensus "
                f"value.",
                f"Kept {unit}s should track the orange line more closely than grey "
                f"ones. If they instead look flatter, the metric is rewarding low "
                f"amplitude rather than matching shape -- try "
                f"SIMILARITY='neg_euclidean' with ZSCORE=False.")

        fa.plot_region_invariance(
            invariance, common, mask=region_mask,
            title=f"Region invariance across species (lvl {cfg['ATLAS_LEVEL']})",
            output_path=os.path.join(fd, "F5_region_invariance.png"))
        _write_legend(
            os.path.join(fd, "F5_region_invariance.txt"),
            "F5 -- Cross-species region invariance",
            "Per region, how much species AGREE on its value. Green bars are the "
            "regions the selection was restricted to.",
            f"Invariance = 1 - rescaled median absolute deviation of the per-species "
            f"values, computed on the {unit}-level profiles. Selection used "
            f"mode={cfg['INVARIANCE_MODE']}, value={cfg['INVARIANCE_VALUE']}.",
            "Grey regions vary across species, so matching them would mean matching "
            "cross-species noise. If very few regions are green, the conserved core "
            "is small and any cross-species claim rests on little.")

        # --- phylogenetic regression, from the ALL-INDIVIDUALS maps -------
        # Deliberately NOT the kept maps. The kept units are chosen to
        # maximise correlation with the cross-species consensus, and the
        # pool they are chosen from differs by an order of magnitude
        # between datasets (25 runs for Dog/Human_Castrillon, 171 for
        # Rat). Taking the top N of a larger pool yields a higher
        # expected correlation even with identical biology and identical
        # noise -- an order-statistics effect, the same structure as the
        # winner's curse. Since pool size is confounded with clade
        # (large rodent datasets, small primate ones), a slope fitted on
        # kept maps could reflect how much sorting each species allowed
        # rather than phylogeny. The all-individuals maps carry no
        # selection pressure at all, so every species is measured under
        # the same regime.
        if phylo_profiles:
            print(f"  [F6] using {phylo_source} ({len(phylo_profiles)} dataset(s)), "
                  f"leave-one-SPECIES-out reference per point")
            ddf_bids = fa.compute_loo_species_distances(
                phylo_profiles, common, value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"])
            if not ddf_bids.empty:
                ddf_bids["source"] = phylo_source
            if len(ddf_bids) >= 3:
                # One row per SPECIES for the fit: two BIDS of one species
                # are not independent observations, and fitting on per-BIDS
                # rows would both pseudo-replicate and give double weight to
                # the species that happen to have two datasets.
                ddf = (ddf_bids.groupby("species", as_index=False)
                       .agg(distance_to_human=("distance_to_human", "mean"),
                            pearson_r=("pearson_r", "mean"),
                            n_bids=("label", "count")))
                ddf["source"] = "all_units"

                # Intra- vs inter-species contrast: the direct answer to
                # "are the cross-species differences protocol artefacts?".
                intra = (ddf_bids.groupby("species")["distance_to_human"]
                         .agg(lambda s: float(np.abs(s.max() - s.min())) if len(s) > 1 else np.nan)
                         .dropna())
                d = ddf["distance_to_human"].values
                inter = float(np.nanmean(np.abs(d[:, None] - d[None, :])[np.triu_indices(len(d), 1)])) \
                    if len(d) > 1 else np.nan
                if len(intra):
                    print(f"  [F6] intra-species spread (|BIDS1-BIDS2|): "
                          f"mean={intra.mean():.4f} over {len(intra)} species "
                          f"{dict(intra.round(4))}")
                    print(f"  [F6] inter-species mean pairwise distance: {inter:.4f}")
                    print(f"  [F6] ratio inter/intra = {inter / intra.mean():.2f} "
                          f"(>1 means datasets of a species agree more with each "
                          f"other than species do with each other)")
                pd.DataFrame(dict(metric=["intra_species_spread_mean",
                                          "inter_species_pairwise_mean",
                                          "ratio_inter_over_intra"],
                                  value=[float(intra.mean()) if len(intra) else np.nan,
                                         inter,
                                         float(inter / intra.mean()) if len(intra) else np.nan])
                             ).to_excel(os.path.join(paths["group_maps"],
                                                     "intra_vs_inter_species.xlsx"), index=False)

                with pd.ExcelWriter(os.path.join(paths["group_maps"],
                                                 "distance_to_consensus.xlsx")) as xw:
                    # "distance_to_human" is the historical column name: it is the
                    # distance to the leave-one-species-out consensus. The
                    # divergence time actually used is exported next to it.
                    for df_ in (ddf_bids, ddf):
                        df_["divergence_mya_used"] = df_["species"].map(fa.DIVERGENCE_TIME_MYA)
                        df_["divergence_source"] = cfg["TREE_PATH"] or "hand-entered (approximate)"
                    ddf_bids.to_excel(xw, sheet_name="per_bids", index=False)
                    ddf.to_excel(xw, sheet_name="per_species", index=False)
                print(f"  [F6] distance_to_consensus.xlsx: {len(ddf_bids)} BIDS row(s), "
                      f"{len(ddf)} species row(s) -- fit on species rows only")
                tr = fa.test_phylogenetic_signal(ddf, exclude_human=False)
                fa.plot_distance_vs_divergence(
                    ddf_bids, test_result=tr, exclude_human=False,
                    title=f"Distance to consensus vs divergence (lvl {cfg['ATLAS_LEVEL']})",
                    output_path=os.path.join(fd, "F6_distance_vs_divergence.png"))
                perm = fa.permutation_test_phylogeny(
                    tr["data"], fa.DIVERGENCE_TIME_MYA, n_perm=cfg["N_PERM_PHYLOGENY"])
                if perm is not None:
                    fa.plot_phylogeny_permutation(
                        perm, title=f"Phylogeny permutation (lvl {cfg['ATLAS_LEVEL']})",
                        output_path=os.path.join(fd, "F7_phylogeny_permutation.png"))
                _write_legend(
                    os.path.join(fd, "F6_distance_vs_divergence.txt"),
                    "F6/F7 -- Distance to consensus vs phylogenetic divergence",
                    f"Distance (1 - r) of each dataset's group map ({phylo_source}) to "
                    f"the leave-one-species-out median across species, against "
                    f"divergence time from human. "
                    "One point per BIDS, so the two datasets of a species appear "
                    "separately; the regression is fitted on the per-species means "
                    "only.",
                    (f"Group maps from {unit}-level 3dttest on EVERY {unit}, with no "
                     f"selection -- deliberately not the kept maps, whose selection "
                     f"pressure varies with pool size across datasets and would bias "
                     f"the slope. " if phylo_profiles is afni_profiles else
                     f"SENSITIVITY ANALYSIS: round-{final_round} KEPT maps "
                     f"(PHYLOGENY_USE_FINAL_ROUND=True). Selection pressure varies with "
                     f"pool size, which is confounded with clade (winner's curse); "
                     f"{'LEAVE_ONE_OUT=False also made the selection circular. ' if not cfg['LEAVE_ONE_OUT'] else ''}"
                     f"Compare with the all-units run before interpreting. ")
                    + f"Distance computed on the {len(common)} common "
                    f"regions. Fit on species means (per-BIDS rows are not "
                    f"independent). F7 shuffles divergence times for an exact null.",
                    "A positive slope means more distantly related species differ more "
                    "from the consensus. Check intra_vs_inter_species.xlsx first: if "
                    "the two BIDS of a species sit closer to each other than species "
                    "do to one another (ratio > 1), the cross-species differences are "
                    "not protocol artefacts. With ~7 species one point moves the fit "
                    "a lot -- trust F7's permutation p over the parametric one.")
            else:
                print(f"  [SKIP] phylogeny figures: only {len(ddf_bids)} usable "
                      f"group-map profile(s), need >= 3")

        # ---- F11-F13: comparative phylogenetic tests on the full tree ----
        if cfg["RUN_PHYLO_COMPARATIVE"] and phylo_profiles:
            print(f"\n{'=' * 70}\n  Phylogenetic comparative tests (Mantel / K / PGLS)\n{'=' * 70}")
            fp.check_tree_against_divergence_times(fa.DIVERGENCE_TIME_MYA)
            mres = fp.run_mantel(phylo_profiles, common, value_col=cfg["VALUE_COL"],
                                 zscore=cfg["ZSCORE"])
            if mres is not None:
                with pd.ExcelWriter(os.path.join(paths["group_maps"],
                                                 "phylo_mantel.xlsx")) as xw:
                    pd.DataFrame(mres["D_func"], index=mres["species"],
                                 columns=mres["species"]).to_excel(xw, sheet_name="D_functional")
                    pd.DataFrame(mres["D_phylo"], index=mres["species"],
                                 columns=mres["species"]).to_excel(xw, sheet_name="D_mrca_Mya")
                    mres["bids_pairs"].to_excel(xw, sheet_name="bids_pairs", index=False)
                    pd.DataFrame([{k: v for k, v in m.items() if k != "null"}
                                  for m in (mres["mantel_pearson"], mres["mantel_spearman"])]
                                 ).to_excel(xw, sheet_name="mantel", index=False)
                fp.plot_mantel(mres, title=f"Pairwise functional vs phylogenetic distance "
                                           f"(lvl {cfg['ATLAS_LEVEL']})",
                               output_path=os.path.join(fd, "F11_mantel.png"))
                _write_legend(
                    os.path.join(fd, "F11_mantel.txt"),
                    "F11 -- Mantel test: functional vs phylogenetic distance",
                    "Top: species x species functional distance (mean 1 - r over "
                    f"every cross-species BIDS pair of {phylo_source} profiles) "
                    "and MRCA age from the time-calibrated tree. Bottom-left: every "
                    "species pair (blue, tested) and every BIDS pair (grey); "
                    "within-species BIDS pairs in orange at 0 Mya (reported, not "
                    "tested). Bottom-right: exact permutation null.",
                    f"Mantel correlation between upper triangles; null = all "
                    f"{mres['mantel_pearson']['n_perm']} permutations of species "
                    f"labels (exact). Pearson and Spearman. {len(common)} common "
                    f"regions, profiles z-scored within dataset.",
                    "Primary phylogenetic test: uses the whole tree, not only the "
                    "distance to human, so Rat/Mouse are no longer pseudo-replicates. "
                    "p(>) is the pre-specified directional test (more distant -> "
                    "more different). Within-species BIDS distances as large as "
                    "cross-species ones mean dataset effects rival species effects. "
                    "Mantel tests can be liberal on phylogenetic data (Harmon & "
                    "Glor 2010) -- read with F12/F13.")

            kdf = fp.regionwise_phylo_signal(phylo_profiles, common,
                                             value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"])
            if not kdf.empty:
                kdf.to_excel(os.path.join(paths["group_maps"], "phylo_blomberg_K.xlsx"),
                             index=False)
                fp.plot_regionwise_k(kdf, output_path=os.path.join(fd, "F12_blomberg_K.png"))
                _write_legend(
                    os.path.join(fd, "F12_blomberg_K.txt"),
                    "F12 -- Phylogenetic signal per region (Blomberg's K)",
                    "Per common region, Blomberg's K of the species-level value "
                    "(median over that species' BIDS of the z-scored group-map "
                    "profile). Green = q_FDR < 0.05.",
                    "K from the Brownian covariance of the time-calibrated tree; "
                    "exact permutation p (all tip shuffles); Benjamini-Hochberg "
                    "across regions.",
                    "K = 1: as much signal as Brownian motion; K << 1: less. With 7 "
                    "species, power is low -- a non-significant region is not "
                    "evidence of absence. Useful to see WHICH regions drive F11.")

            ddf_pg = fa.compute_loo_species_distances(
                phylo_profiles, common, value_col=cfg["VALUE_COL"],
                zscore=cfg["ZSCORE"], verbose=False)
            if len(ddf_pg) >= 4:
                sp_df = (ddf_pg.groupby("species", as_index=False)["distance_to_human"]
                         .mean().rename(columns={"distance_to_human": "distance_to_consensus"}))
                n_runs = {sp: sum(len(runs_by_label[l]) for l in runs_by_label
                                  if species_of_label.get(l) == sp)
                          for sp in sp_df["species"]}
                sp_df["log_n_runs"] = sp_df["species"].map(
                    lambda s: np.log(n_runs[s]) if n_runs.get(s) else np.nan)
                covs = ["log_n_runs"]
                qc_cov = fp.species_covariates_from_qc(
                    df_qc_by_label, species_of_label, list(cfg["PGLS_QC_COVARIATES"]))
                if not qc_cov.empty:
                    sp_df = sp_df.merge(qc_cov, on="species", how="left")
                    covs += [c for c in cfg["PGLS_QC_COVARIATES"] if c in sp_df.columns]
                for cname, values in (cfg["SPECIES_COVARIATES"] or {}).items():
                    sp_df[cname] = sp_df["species"].map(values)
                    covs.append(cname)
                pg = fp.pgls_battery(sp_df, "distance_to_consensus", covs)
                if not pg.empty:
                    with pd.ExcelWriter(os.path.join(paths["group_maps"],
                                                     "phylo_pgls.xlsx")) as xw:
                        pg.to_excel(xw, sheet_name="pgls", index=False)
                        sp_df.to_excel(xw, sheet_name="species_data", index=False)
                    print(f"  [saved] phylo_pgls.xlsx")

        print(f"  [saved] figures + legends -> 06_figures/")

    # ---- STEP 7: F6 at several atlas levels ------------------------------
    multilevel_results = {}
    if cfg["RUN_MULTILEVEL_F6"]:
        print(f"\n{'=' * 70}\n  STEP 7: F6 at levels {cfg['MULTILEVEL_F6_LEVELS']}\n{'=' * 70}")
        print(f"  [STEP 7] source: {phylo_source}")
        multilevel_results = step6_multilevel_f6(
            sri, fa, fu, cfg, paths, bids_configs, phylo_index, members_by_label,
            levels=cfg["MULTILEVEL_F6_LEVELS"])

    # ---- STEP 8: per-run distance to the kept reference -------------------
    dist_df = pd.DataFrame()
    corr_df = pd.DataFrame()
    if cfg["RUN_DISTANCE_TO_REF"]:
        # Dispersion is measured around the dataset's reference from the
        # FINAL selection round (its kept map), whatever
        # PHYLOGENY_USE_FINAL_ROUND says -- that flag only concerns the
        # winner's-curse issue of the cross-species phylogenetic tests.
        use_final = (final_round >= 2 and not map_index_final.empty
                     and cfg["DISTANCE_REF_GROUP"] == "kept")
        ref_map_index = map_index_final if use_final else map_index
        ref_group = "kept_final" if use_final else cfg["DISTANCE_REF_GROUP"]
        ref_selections = selections_final if use_final else selections
        ref_name = (f"round-{final_round} kept" if use_final
                    else f"round-1 {cfg['DISTANCE_REF_GROUP']}")
        print(f"\n{'=' * 70}\n  STEP 8: per-run distance to each BIDS' "
              f"{ref_name} group-map reference\n{'=' * 70}")
        ref_profiles = step_build_reference_profiles(
            sri, fu, cfg, bids_configs, ref_map_index, members_by_label,
            ref_group=ref_group)
        dist_df = fa.compute_run_distance_to_reference(
            runs_by_label, ref_profiles, common, value_col=cfg["VALUE_COL"],
            kept_keys_by_label=ref_selections)
        if not dist_df.empty:
            dist_df.to_excel(os.path.join(paths["group_maps"], "run_distance_to_reference.xlsx"),
                            index=False)
            fa.plot_run_distance_distribution(
                dist_df, output_path=os.path.join(
                    paths["figures"], "F8_run_distance_to_reference.png"),
                title=f"Per-run correlation to each BIDS' {ref_name} group-map reference "
                      f"(x = round-{final_round if use_final else 1} kept)")
            fa.plot_run_distance_std_by_bids(
                dist_df, output_path=os.path.join(
                    paths["figures"], "F9_run_distance_std_by_bids.png"))

            # ---- F8b: distance to the FINAL round's selection criterion ----
            # The exact profile the final round scored each BIDS's runs
            # against (round 1: LOO/non-LOO consensus; round 2: the BIDS's
            # own round-1 kept mean; round >= 3: the consensus rebuilt from
            # the previous round's kept means). Crosses MUST sit at the top
            # by construction -- this checks the implementation, it is not
            # a finding. F8 (AFNI kept map) is the independent check.
            ref_profiles_criterion = {
                label: (species_of_label.get(label), ref)
                for label, ref in final_criterion.items() if ref is not None}
            dist_df_b = fa.compute_run_distance_to_reference(
                runs_by_label, ref_profiles_criterion, common,
                value_col=cfg["VALUE_COL"], kept_keys_by_label=selections_final)
            if not dist_df_b.empty:
                dist_df_b.to_excel(os.path.join(
                    paths["group_maps"], "run_distance_to_selection_criterion.xlsx"),
                    index=False)
                fa.plot_run_distance_distribution(
                    dist_df_b, output_path=os.path.join(
                        paths["figures"], "F8b_run_distance_to_selection_criterion.png"),
                    title=f"Per-run correlation to the round-{final_round} SELECTION "
                          f"CRITERION (what the ranking used, not the AFNI map)")

            # ---- STEP 9: correlate against external QC metrics ------------
            if df_qc_by_label:
                print(f"\n{'=' * 70}\n  STEP 9: QC metrics vs distance to reference\n{'=' * 70}")
                corr_df = fa.correlate_qc_with_distance(
                    dist_df, df_qc_by_label, qc_metrics=cfg["QC_METRICS_FOR_CORRELATION"])
                if not corr_df.empty:
                    corr_df.to_excel(os.path.join(
                        paths["group_maps"], "qc_vs_distance_correlation_per_bids.xlsx"),
                        index=False)
                    meta_df = fa.meta_analyze_qc_correlations(
                        corr_df, min_bids=cfg["QC_META_MIN_BIDS"], alpha=cfg["QC_META_ALPHA"])
                    if not meta_df.empty:
                        meta_df.to_excel(os.path.join(
                            paths["group_maps"], "qc_vs_distance_meta_analysis.xlsx"),
                            index=False)
                        fa.plot_qc_distance_correlation(
                            meta_df, alpha=cfg["QC_META_ALPHA"], output_path=os.path.join(
                                paths["figures"], "F10_qc_vs_distance_correlation.png"))
            else:
                print("  [SKIP] STEP 9: no df_qc_by_label given to run() -- "
                      "pass QC data to correlate it against the per-run distances")

    _write_methods(cfg, paths, sorted(profiles_by_label), len(common), unit_counts,
                   final_round=final_round, phylo_source=phylo_source,
                   n_species=len({species_of_label.get(l) for l in phylo_profiles}))
    print(f"\n{'#' * 70}\n  DONE -- {paths['root']}\n{'#' * 70}")
    return dict(config=cfg, paths=paths, common_regions=common,
                multilevel_f6=multilevel_results,
                run_distance_to_ref=dist_df, qc_vs_distance=corr_df,
                profiles_by_label=profiles_by_label,
                members_by_label=members_by_label,
                runs_by_label=runs_by_label,
                unit_counts=unit_counts,
                selection=sel, invariance=invariance, region_mask=region_mask,
                all_units_index=all_units_index, map_index=map_index)
