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
   SBA_group_maps.py (former Figure_5_SBA.py) wrote `Grp_SBA_3dTTEST_All<suffix>` (category glued
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
import Exemples.Study_EDNiX.analysis.Fig5.fig5_fingerprints as ff  # F6ab, fingerprints, F15-F18
import Exemples.Study_EDNiX.analysis.Fig5.fig5_paper as fpp         # panels of Figures 5 and 6
import Exemples.Study_EDNiX.analysis.Fig5.fig5_supp_paper as fsp    # print-size panels, Figs 4-5 / Supp 7-11

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
    PROFILE_METHOD="percentiles",  # comment le profil d'une carte est construit :
        # "percentiles"     : (défaut) 2 valeurs par région, en r : percentile haut et
        #   percentile bas de ses voxels, fixés par TAIL_Q (0.2 -> P80 / P20 ;
        #   0.1 -> P90 / P10 ; 0.25 -> P75 / P25). VALUE_COL devient q20q80_r
        #   (automatique). Ancien nom accepté : "quantile_posneg" ;
        # "median_r"        : 1 valeur par région = médiane des voxels ;
        # "fraction_posneg" : 2 valeurs par région, en % des voxels de la région :
        #     <région>|pos = % dans les TAIL_Q supérieurs de la carte,
        #     <région>|neg = % dans les TAIL_Q inférieurs ;
        #   percentiles de la carte pris sur les voxels des régions communes.
        #   VALUE_COL devient pct_tail20 (automatique).
        #   G et D poolés par nombre de voxels. Tout le reste de l'analyse est
        #   inchangé (normalisation ZSCORE comprise).
    TAIL_Q=0.2,                 # taille des queues : 0.2 ou "quintile" (P80/P20, 20 %),
        # 0.1 ou "decile" (P90/P10, 10 %), 0.25 ou "quartile" (P75/P25, 25 %)
    RUN_PROFILE_RELIABILITY=True,  # F22 : fiabilité des 3 constructions de profil sur
        # tous les runs, avant sélection (discriminabilité, Bridgeford et al. 2021),
        # dans le même run quel que soit PROFILE_METHOD -- critère de choix du profil

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
    RUN_PHYLO_COMPARATIVE=False,  # Mantel exact + K de Blomberg + PGLS (F11-F13).
        # Désactivé par défaut : la phylogénie de la Figure 5 se limite à F6
        # (F6 d'origine + F6ab). Code conservé, remettre True pour le relancer.

    # --- figures d'interprétation (fig5_fingerprints.py) --------------------
    RUN_FINGERPRINTS=True,   # F6ab (deux références), F1/F1c en fingerprints,
        # F15 (fingerprint vs référence), F16 (proximité entre datasets),
        # F17 (avant/après sélection), F18 (état d'éveil). Python seul.
    REFERENCE_SPECIES="Human",   # référence de F6b/F15b/F18
    SEED_REGION_EXCLUDE="auto",  # région du niveau ATLAS_LEVEL qui contient le seed (sa
        # valeur reflète surtout l'autocorrélation du seed, pas le couplage) :
        # "auto" = trouvée dans l'atlas de chaque dataset (voxels du label du seed,
        # niveau ATLAS_LEVEL) ; liste de noms = imposée ; None = aucune (les deux
        # clés suivantes n'ont alors pas d'effet).
    SEED_REGION_EXCLUDE_SELECTION=False,  # True : retirée des profils que la sélection
        # note (standardisation des canaux, invariance, consensus, score) ; False
        # (défaut) : gardée -- la sélection est alors celle d'avant son retrait (V12).
    SEED_REGION_EXCLUDE_ANALYSES=True,    # True (défaut) : retirée de tous les profils
        # utilisés après la sélection (partie 2 : statistiques et figures).
        # Les deux sont des clés de la partie 1 (les profils des deux espaces sont
        # construits avec la sélection) : en changer une demande de relancer la partie 1.
    PROFILE_KIND="seed",     # "seed" (cartes SBA) | "matrix" (fig5_matrices) :
        # change seulement le rendu des fingerprints (polaire vs matrice).

    # --- panneaux du papier (fig5_paper.py) -------------------------------------
    RUN_PAPER_FIGURES=True,  # 06_figures/Fig5, 06_figures/Fig6 (un fichier par panneau,
        # sans titre, n et p seulement) ; toutes les autres figures sont
        # déplacées dans 06_figures/supplementary/ en fin de run.
    BRAIN_MASS_XLSX="/scratch2/EDNiX/Herculano-Houzel_mammalian_brain_dataset_complete.xlsx",
    BRAIN_MASS_G=None,       # {espèce: g} -- remplace le xlsx si donné
    SURFACE_COMMON_SCALE=True,   # 5a : re-rendu wb_command de toutes les cartes d'un
        # dataset (tous les runs, gardés, non gardés) avec UNE échelle de couleur fixe
    SURFACE_COMMON_PCT=90.0,     # échelle = +/- ce percentile de |z| sur ces cartes
    FIG6D=dict(  # 6d (et S6d) : valeur de chaque région, run par run sur les runs
        # sélectionnés, Human + Macaque contre les 5 autres espèces (espèce = unité du
        # test, un test exact par région, non corrigé). Unité = celle de PROFILE_METHOD
        # (percentiles : P80 en 6d, P20 en 6d' ; fraction_posneg : % dans la queue
        # haute / basse ; median_r : r médian) ; les autres indices en supplémentaire.
        regions=("posterior_parietal_cortex", "prefrontal_cortex_(pfc)",
                 "somatosensory_cortex", "visual_striate_cortex", "striatum",
                 "motor_and_premotor"),
        group=("Human", "Macaque")),
    WITHIN_ACROSS_TEST="wilcoxon",  # 5b : test affiché sur la figure, intra- vs
        # inter-espèce. "wilcoxon" = Wilcoxon rank-sum (Mann-Whitney) sur les paires de
        # runs, unilatéral, non corrigé : les paires partagent des runs, le p est donc
        # optimiste. "permutation" = espèces permutées entre datasets (respecte la
        # dépendance). Les deux sont toujours dans le .txt.
    WITHIN_ACROSS_PAIRS="other_dataset",  # 5b / F2 / F2b : paires comparées.
        # "other_dataset" (défaut) : même espèce dans un AUTRE dataset vs espèces
        # différentes -- les paires d'un même dataset (même protocole, même scanner,
        # souvent même animal) ne séparent pas l'espèce du dataset et sont exclues du
        # test (montrées sur F2/F2b). "all" : toutes les paires de même espèce, même
        # dataset compris (5b de la V14, dominé par les paires du même dataset).
    SURFACE_ZOOM={"Human": 1.0, "Macaque": 1.0, "Dog": 1.0, "Marmoset": 0.9,
                  "Mouselemur": 0.9, "Rat": 0.9, "Mouse": 0.9},  # 5a/S5a : taille relative
        # des cerveaux par espèce (1 = remplit la case), après réduction automatique
        # des blancs entre les vues
    SURFACE_DATASETS_PER_ROW=2,  # 5a (gardés | non gardés, en miroir) et S5a (tous |
                                 # gardés | non gardés) : datasets par ligne ; 1 = un par ligne
    COMMON_REGIONS_OVERRIDE=None,  # liste de régions communes imposée (run de
        # contrôle sans cartes AFNI : reprendre celles du run principal)
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
        state=os.path.join(root, "07_state"),
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
    if cfg.get("PROFILE_METHOD") == "percentiles":
        q = 100 * float(cfg["TAIL_Q"])
        profile_par = (
            f"profile with two values per atlas region at level {cfg['ATLAS_LEVEL']}, in the "
            f"units of the map: the {100 - q:g}th and the {q:g}th percentiles of the region's "
            f"valid voxel values ({cfg['VALUE_COL']}), left and right voxels pooled. With the "
            f"median they describe the location, spread and asymmetry of the region's "
            f"distribution; a global offset or gain of the map moves every value alike, so "
            f"the Pearson similarity between profiles ignores it. Each channel was "
            f"standardised over the common regions, so that the Pearson r between two "
            f"profiles is the mean of the two channel correlations.")
    elif cfg.get("PROFILE_METHOD") == "fraction_posneg":
        q = 100 * float(cfg["TAIL_Q"])
        profile_par = (
            f"profile with two values per atlas region at level {cfg['ATLAS_LEVEL']}, both in "
            f"% of the region's valid voxels: the share of its voxels in the top {q:g} % of "
            f"the map ('pos') and in the bottom {q:g} % ('neg') ({cfg['VALUE_COL']}). The two "
            f"percentiles of each map were taken over the voxels of the common regions only "
            f"(both hemispheres), and left and right were pooled by voxel count. The profile "
            f"is rank-based, hence invariant to any monotonic change of the map (global "
            f"offset, gain); without spatial structure every region would score {q:g} % in "
            f"both tails. Each channel was standardised over the common regions, so that the "
            f"Pearson r between two profiles is the mean of the two channel correlations.")
    else:
        profile_par = (f"profile: the median correlation ({cfg['VALUE_COL']}) within each "
                       f"atlas region at level {cfg['ATLAS_LEVEL']}.")
    found = cfg.get("_seed_regions_found") or cfg.get("_seed_regions_excluded") or []
    out_sel = cfg.get("_seed_regions_excluded_selection", cfg.get("_seed_regions_excluded"))
    out_ana = cfg.get("_seed_regions_excluded")
    if found and out_sel and out_ana:
        profile_par += (f" The region that contains the seed ({', '.join(found)}, found from "
                        f"the seed's voxels in each atlas) was removed from every profile: its "
                        f"value reflects the seed's own autocorrelation, not coupling.")
    elif found and out_ana:
        profile_par += (f" The region that contains the seed ({', '.join(found)}, found from "
                        f"the seed's voxels in each atlas) was kept in the profiles scored by "
                        f"the run selection (channel standardisation, invariance ranking, "
                        f"consensus) and removed from every profile used afterwards, because "
                        f"its value reflects the seed's own autocorrelation, not coupling.")
    elif found and out_sel:
        profile_par += (f" The region that contains the seed ({', '.join(found)}) was removed "
                        f"from the profiles scored by the run selection only.")
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
{profile_par} {maps_par} The set of regions common to every
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


def _run_region_indices(sri, runs_by_label, common, cfg, species_of_label):
    """Long table: one row per run x common region, every regional index
    (sri.region_indices) at q = TAIL_Q. Empty when the profiles carry no voxel
    values (matrix profiles)."""
    if cfg.get("PROFILE_KIND", "seed") != "seed" or not hasattr(sri, "region_indices"):
        return pd.DataFrame()
    print(f"\n{'=' * 70}\n  Regional indices of every run (6d, S6d)\n{'=' * 70}")
    rows, n_skip = [], 0
    for lab, profs in runs_by_label.items():
        for key, prof in profs.items():
            if prof is None or "_sketch" not in prof.columns:
                n_skip += 1
                continue
            try:
                t = sri.region_indices(prof, common, q=cfg["TAIL_Q"])
            except Exception as e:
                print(f"  [WARN] {lab} {key}: regional indices failed ({e})")
                continue
            if t.empty:
                continue
            t = t[t["region"].isin(common)]
            k = tuple(str(x) for x in key)
            rows.append(t.assign(label=lab, species=species_of_label.get(lab), subject=k[0],
                                 session=k[1] if len(k) > 1 else "",
                                 run=k[2] if len(k) > 2 else "",
                                 seed=k[3] if len(k) > 3 else ""))
    if n_skip:
        print(f"  [WARN] {n_skip} run profile(s) without voxel values -- no regional indices")
    if not rows:
        return pd.DataFrame()
    tab = pd.concat(rows, ignore_index=True)
    print(f"  {tab[['label', 'subject', 'session', 'run']].drop_duplicates().shape[0]} run(s) "
          f"x {tab['region'].nunique()} region(s), q = {cfg['TAIL_Q']:g}")
    return tab


def _seed_regions_to_exclude(sri, cfg, bids_configs, profiles_by_label, common, paths):
    """Common regions that contain the seed (SEED_REGION_EXCLUDE). 'auto': for
    every dataset, the ATLAS_LEVEL region under the seed's voxels in a
    representative atlas volume; the union is returned. Where they are then
    removed is decided by the caller: SEED_REGION_EXCLUDE_SELECTION (profiles
    scored by the selection) and SEED_REGION_EXCLUDE_ANALYSES (every profile
    used after it). Writes 02_selection/seed_region_exclusion.xlsx."""
    rule = cfg.get("SEED_REGION_EXCLUDE", "auto")
    if not rule or cfg.get("PROFILE_KIND", "seed") != "seed":
        return []
    excl_sel = bool(cfg.get("SEED_REGION_EXCLUDE_SELECTION", True))
    excl_ana = bool(cfg.get("SEED_REGION_EXCLUDE_ANALYSES", True))
    if not isinstance(rule, str):
        wanted = {str(r).lower().strip() for r in rule}
        absent = sorted(wanted - set(common))
        if absent:
            print(f"  [seed region] [WARN] SEED_REGION_EXCLUDE {absent} not among the common "
                  f"regions -- ignored")
        return sorted(wanted & set(common))
    if rule != "auto":
        raise ValueError(f"SEED_REGION_EXCLUDE must be 'auto', a list of regions or None, "
                         f"got {rule!r}")
    rows = []
    for c in bids_configs:
        if c["label"] not in profiles_by_label:
            continue
        try:
            atlas = sri.pick_representative_atlas_volume(c["bids_dir"])
            if not atlas:
                rows.append(dict(label=c["label"], parent=None,
                                 reason=f"no atlas volume under {c['bids_dir']}"))
                continue
            ldf = sri.parse_itk_snap_label_file(c["label_path"])
            res = sri.find_seed_parent_region(atlas, ldf, cfg["SEED"], cfg["ATLAS_LEVEL"])
            rows.append(dict(label=c["label"], atlas=atlas, **{k: (str(v) if k == "per_parent"
                                                                  else v)
                                                              for k, v in res.items()}))
        except Exception as e:
            rows.append(dict(label=c["label"], parent=None, reason=f"failed: {e}"))
    tab = pd.DataFrame(rows)
    found = tab[tab["parent"].notna()] if not tab.empty else tab
    for r in tab.itertuples():
        if getattr(r, "parent", None):
            print(f"  [seed region] {r.label}: seed labels {r.seed_labels} "
                  f"({r.n_seed_voxels} voxels) lie in '{r.parent}' at level "
                  f"{cfg['ATLAS_LEVEL']} ({100 * r.share:.0f} % of them)")
        else:
            print(f"  [seed region] [WARN] {r.label}: {getattr(r, 'reason', '')}")
    parents = sorted(set(found["parent"]) & set(common)) if not found.empty else []
    if not found.empty:
        outside = sorted(set(found["parent"]) - set(common))
        if outside:
            print(f"  [seed region] {outside} contain the seed but are not common regions")
    if tab.empty or found.empty:
        print("  [seed region] [WARN] seed region not found in any atlas -- nothing removed; "
              "set SEED_REGION_EXCLUDE to the region name to remove it")
    is_seed = tab.get("parent", pd.Series(dtype=object)).isin(parents)
    tab.assign(seed_region=is_seed, excluded_from_selection=is_seed & excl_sel,
               excluded_from_analyses=is_seed & excl_ana).to_excel(
        os.path.join(paths["selection"], "seed_region_exclusion.xlsx"), index=False)
    return parents


def _write_profile_reliability_legend(png, rel, tests, cfg):
    lines = []
    for r in rel.sort_values(["label", "method"]).itertuples():
        lines.append(f"{r.label} / {r.method}: discriminability {r.discriminability:.3f} "
                     f"(p vs chance {r.p_discriminability:.3f}); r inter-run "
                     f"{r.r_inter_run:.3f}, same animal {r.r_same_animal:.3f}, different "
                     f"animals {r.r_inter_animal:.3f} (gap {r.gap_animal:+.3f}); "
                     f"{r.n_runs} runs, {r.n_animals} animals, {r.n_values} values")
    for k, d in tests.items():
        for m, t in d.items():
            lines.append(f"PAIRED {k}: {m} - median_r = {t['mean_diff']:+.4f}, higher in "
                         f"{t['n_better']}/{t['n']} datasets, p = {t['p']:.3f} (exact sign-flip, "
                         f"two-sided)")
    _write_legend(
        os.path.splitext(png)[0] + ".txt",
        "F22 -- Which construction of the profile is the most reliable?",
        "Every run of every dataset with several runs per animal, BEFORE any selection. For "
        "each profile method: left, the discriminability (Bridgeford et al. 2021 PLoS Comput "
        "Biol 17:e1009279) = probability that a run is closer (1 - r between profiles) to "
        "another run of the same animal than to a run of another animal, 0.5 = chance; right, "
        "median r between runs of the same animal minus median r between runs of different "
        "animals. One line per dataset, black = mean. Methods: median_r (median of the "
        f"region's voxels); fraction_posneg (% of the region's voxels in the top / bottom "
        f"{100 * cfg['TAIL_Q']:g} % of the map, percentiles over the common regions); "
        f"percentiles ({100 - 100 * cfg['TAIL_Q']:g}th / {100 * cfg['TAIL_Q']:g}th "
        f"percentile of the region's voxels). All three are computed from the same voxels "
        f"of the same runs and normalised as in the analysis (ZSCORE={cfg['ZSCORE']!r}).",
        "p of a dataset: animal labels shuffled across its runs (500 permutations, "
        "one-sided). Method comparison: difference with median_r per dataset, exact sign-flip "
        "over datasets (2^n), two-sided. The criterion was fixed before looking at the "
        "cross-species results and does not involve them.\n" + "\n".join(lines),
        "Discriminability is rank-based, hence comparable between profiles of different "
        "length (16 values for median_r, 32 for the two-channel methods); the r values "
        "themselves are not. A method that is more discriminable measures the animal more "
        "reliably from one run. For Mouse and Marmoset_Hori, runs of one animal span arousal "
        "states, which lowers every method alike.")


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
    methods = ("median_r", "fraction_posneg", "percentiles")
    if cfg["PROFILE_METHOD"] == "quantile_posneg":
        cfg["PROFILE_METHOD"] = "percentiles"
        warns.append("PROFILE_METHOD 'quantile_posneg' is now called 'percentiles' -- same "
                     "profile")
    if cfg["PROFILE_METHOD"] not in methods:
        raise ValueError(f"PROFILE_METHOD must be one of {methods}, "
                         f"got {cfg['PROFILE_METHOD']!r}")
    names = {"decile": 0.1, "deciles": 0.1, "quintile": 0.2, "quintiles": 0.2,
             "quartile": 0.25, "quartiles": 0.25}
    tq = cfg["TAIL_Q"]
    tq = names.get(tq.strip().lower(), tq) if isinstance(tq, str) else tq
    try:
        tq = float(tq)
    except (TypeError, ValueError):
        raise ValueError(f"TAIL_Q must be a number in (0, 0.5) or one of {sorted(names)}, "
                         f"got {cfg['TAIL_Q']!r}") from None
    if not 0.0 < tq < 0.5:
        raise ValueError(f"TAIL_Q must be in (0, 0.5), got {cfg['TAIL_Q']!r}")
    cfg["TAIL_Q"] = tq
    if cfg.get("WITHIN_ACROSS_TEST", "wilcoxon") not in ("wilcoxon", "permutation"):
        raise ValueError(f"WITHIN_ACROSS_TEST must be 'wilcoxon' or 'permutation', got "
                         f"{cfg['WITHIN_ACROSS_TEST']!r}")
    if cfg.get("PROFILE_KIND", "seed") != "seed":
        # matrix profiles are edge values, not voxel maps
        if cfg["PROFILE_METHOD"] != "median_r":
            warns.append(f"PROFILE_METHOD={cfg['PROFILE_METHOD']!r} ignored for PROFILE_KIND="
                         f"{cfg.get('PROFILE_KIND')!r} (edge profiles): median_r used")
            cfg["PROFILE_METHOD"] = "median_r"
        cfg["RUN_PROFILE_RELIABILITY"] = False
    elif cfg["PROFILE_METHOD"] != "median_r":
        q = float(cfg["TAIL_Q"])
        lo, hi = int(round(100 * q)), int(round(100 * (1 - q)))
        cfg["VALUE_COL"] = (f"pct_tail{lo}" if cfg["PROFILE_METHOD"] == "fraction_posneg"
                            else f"q{lo}q{hi}_r")
        if cfg["RUN_MULTILEVEL_F6"]:
            warns.append(f"RUN_MULTILEVEL_F6 switched off: with PROFILE_METHOD="
                         f"{cfg['PROFILE_METHOD']!r} the profiles are converted on the common "
                         f"regions of ATLAS_LEVEL only")
            cfg["RUN_MULTILEVEL_F6"] = False
        if cfg["THRESHOLD_ABS"] is not None:
            warns.append(f"THRESHOLD_ABS={cfg['THRESHOLD_ABS']} is in median_r units -- "
                         f"set to None with PROFILE_METHOD={cfg['PROFILE_METHOD']!r}")
            cfg["THRESHOLD_ABS"] = None
    for k in ("SEED_REGION_EXCLUDE_SELECTION", "SEED_REGION_EXCLUDE_ANALYSES"):
        if not isinstance(cfg[k], (bool, np.bool_)):
            raise ValueError(f"{k} must be True or False, got {cfg[k]!r}")
    if not cfg["SEED_REGION_EXCLUDE"] and (cfg["SEED_REGION_EXCLUDE_SELECTION"]
                                           or cfg["SEED_REGION_EXCLUDE_ANALYSES"]):
        warns.append("SEED_REGION_EXCLUDE=None: no seed region is looked for, "
                     "SEED_REGION_EXCLUDE_SELECTION / _ANALYSES have no effect")
    if cfg["SEED_REGION_EXCLUDE_SELECTION"] and not cfg["SEED_REGION_EXCLUDE_ANALYSES"]:
        warns.append("seed region removed from the selection but KEPT in the analyses: its "
                     "autocorrelation peak then enters every similarity of part 2")
    if cfg.get("WITHIN_ACROSS_PAIRS", "other_dataset") not in ("other_dataset", "all"):
        raise ValueError(f"WITHIN_ACROSS_PAIRS must be 'other_dataset' or 'all', got "
                         f"{cfg['WITHIN_ACROSS_PAIRS']!r}")
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
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

def _archive_previous_outputs(paths, figures_only=False):
    """
    Move (never delete) the result files of a previous run in the same
    OUTPUT_ROOT into <root>/_previous_runs/<timestamp>/, BEFORE this run
    writes anything. Otherwise a run that stops early (e.g. REQUIRE_ALL_MAPS)
    leaves the previous run's figures, tables, methods.txt and log in place,
    looking current. AFNI maps and their balance records (*.json) stay where
    they are so OVERWRITE_MAPS=False can still reuse them after the
    input_runs/ provenance check. The saved part-1 state (07_state) is moved
    too, so a part 1 that stops early cannot leave an older selection behind
    for part 2. figures_only=True (part 2 alone): only 06_figures and the
    part-2 log. Returns the archive folder or None.
    """
    import datetime
    root = paths["root"]
    cand = []
    logs = ("00_log_figures.txt",) if figures_only else ("00_log.txt", "00_log_figures.txt",
                                                          "methods.txt", "Methods_article.txt")
    for f in logs:
        if os.path.isfile(os.path.join(root, f)):
            cand.append(os.path.join(root, f))
    for key in (("figures",) if figures_only else ("figures", "selection", "state")):
        d = paths[key]
        cand += [os.path.join(d, f) for f in os.listdir(d)] if os.path.isdir(d) else []
    if not figures_only:
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


# ═══════════════════════════════════════════════════════════════════════════
# The two parts: SELECTION (part 1) and FIGURES (part 2)
# ═══════════════════════════════════════════════════════════════════════════
# Part 1 = everything that depends on the data and on the selection settings:
# run profiles, run QC, AFNI all-units maps, common regions, regional indices of
# every run, F22 table, selection rounds, AFNI kept / not-kept maps, jpgs,
# surface renders, profiles of every group map. It ends by writing
# <OUTPUT_ROOT>/07_state/selection_state.pkl.
# Part 2 = every statistic and figure, from that file only (no AFNI, no NIfTI
# read except STEP 7, which is off by default): it can be rerun as often as
# needed, with other figure settings, in minutes.
# run(..., part="both") does both in one go (the default).

PARTS = ("both", "selection", "figures")
STATE_FILE = "selection_state.pkl"

# Keys that only change part 2: they can differ between the selection run and a
# later run(part="figures"). Every other key changes the selection (or the
# maps): giving it a different value with part="figures" is an error.
FIGURE_KEYS = frozenset({
    "MAKE_FIGURES", "RUN_FINGERPRINTS", "REFERENCE_SPECIES", "RUN_PAPER_FIGURES",
    "BRAIN_MASS_XLSX", "BRAIN_MASS_G", "SURFACE_ZOOM", "SURFACE_DATASETS_PER_ROW",
    "FIG6D", "WITHIN_ACROSS_TEST", "WITHIN_ACROSS_PAIRS", "N_PERM_PHYLOGENY", "RUN_PHYLO_COMPARATIVE",
    "TREE_PATH", "SPECIES_TO_TREE_TIP", "PGLS_QC_COVARIATES", "SPECIES_COVARIATES",
    "RUN_MULTILEVEL_F6", "MULTILEVEL_F6_LEVELS", "RUN_DISTANCE_TO_REF",
    "DISTANCE_REF_GROUP", "QC_METRICS_FOR_CORRELATION", "QC_META_MIN_BIDS",
    "QC_META_ALPHA", "PHYLOGENY_USE_FINAL_ROUND"})
# derived automatically in part 1: None / absent in the launcher means "as saved"
_DERIVED_KEYS = frozenset({"MAP_TAG", "RESULTS_BASE_DIR", "VALUE_COL"})


def _state_path(root):
    return os.path.join(root, "07_state", STATE_FILE)


def save_state(state):
    """Write the part-1 state (pickle) and a readable summary next to it."""
    import pickle
    import datetime
    path = _state_path(state["cfg"]["OUTPUT_ROOT"])
    os.makedirs(os.path.dirname(path), exist_ok=True)
    state = dict(state, saved=datetime.datetime.now().isoformat(timespec="seconds"),
                 versions=dict(python=sys.version.split()[0], numpy=np.__version__,
                               pandas=pd.__version__))
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        pickle.dump(state, f, protocol=4)
    os.replace(tmp, path)
    with open(os.path.splitext(path)[0] + ".txt", "w", encoding="utf-8") as f:
        f.write(f"Part-1 state of {state['cfg']['OUTPUT_ROOT']}\nsaved {state['saved']} "
                f"(python {state['versions']['python']}, numpy {state['versions']['numpy']}, "
                f"pandas {state['versions']['pandas']})\n\n")
        f.write(f"datasets: {sorted(state['profiles_by_label'])}\n")
        f.write(f"final round: {state['final_round']}\n")
        f.write(f"common regions, analyses ({len(state['anat_common'])}): "
                f"{state['anat_common']}\n")
        ss = state.get("selection_space")
        if ss is not None:
            f.write(f"common regions, selection ({len(ss['anat_common'])}): "
                    f"{ss['anat_common']}\n")
        if state.get("seed_regions"):
            f.write(f"seed region(s): {state['seed_regions']} -- removed from the selection: "
                    f"{bool(state['cfg'].get('SEED_REGION_EXCLUDE_SELECTION'))}, from the "
                    f"analyses: {bool(state['cfg'].get('SEED_REGION_EXCLUDE_ANALYSES'))}\n")
        f.write(f"kept per dataset: " + ", ".join(
            f"{k} {len(v)}" for k, v in sorted(state["selections_final"].items())) + "\n\n")
        f.write("Selection settings (part 2 refuses another value for these):\n")
        for k, v in sorted(state["cfg"].items()):
            if k not in FIGURE_KEYS:
                f.write(f"  {k} = {v!r}\n")
    print(f"  [saved] 07_state/{STATE_FILE} (+ .txt summary) -- part 2 can now be rerun "
          f"alone with run(..., part='figures')")
    return path


def load_state(root):
    import pickle
    path = _state_path(root)
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"{path} not found: run part 1 first (run(..., part='selection') or "
            f"part='both') with OUTPUT_ROOT={root!r}")
    with open(path, "rb") as f:
        return pickle.load(f)


def _merge_figure_config(saved_cfg, user_config):
    """Saved part-1 config + the part-2 keys of `user_config`. A selection key
    with another value than the one the selection was made with raises."""
    saved_cfg = dict(saved_cfg)
    if "SEED_REGION_EXCLUDE_SELECTION" not in saved_cfg:
        # state written before the two switches existed (V14): the seed region
        # was removed everywhere when SEED_REGION_EXCLUDE was set
        excl = bool(saved_cfg.get("SEED_REGION_EXCLUDE"))
        saved_cfg["SEED_REGION_EXCLUDE_SELECTION"] = excl
        saved_cfg["SEED_REGION_EXCLUDE_ANALYSES"] = excl
    cfg = dict(saved_cfg)
    user_config = dict(user_config or {})
    unknown = sorted(set(user_config) - set(CONFIG) - set(LEGACY_CONFIG_KEYS))
    if unknown:
        raise KeyError(f"Unknown CONFIG key(s): {unknown} -- typo? Valid keys: {sorted(CONFIG)}")
    if user_config.get("PROFILE_METHOD") == "quantile_posneg":
        user_config["PROFILE_METHOD"] = "percentiles"
    if isinstance(user_config.get("TAIL_Q"), str):
        names = {"decile": 0.1, "deciles": 0.1, "quintile": 0.2, "quintiles": 0.2,
                 "quartile": 0.25, "quartiles": 0.25}
        user_config["TAIL_Q"] = names.get(user_config["TAIL_Q"].strip().lower(),
                                          user_config["TAIL_Q"])
    changed = []
    for k, v in user_config.items():
        if k in FIGURE_KEYS:
            cfg[k] = v
        elif k in _DERIVED_KEYS and v is None:
            continue
        elif k in LEGACY_CONFIG_KEYS:
            continue
        elif saved_cfg.get(k, CONFIG.get(k)) != v:
            changed.append(f"{k}: selection made with {saved_cfg.get(k, CONFIG.get(k))!r}, "
                           f"now {v!r}")
    if changed:
        raise ValueError(
            "part='figures' reuses the saved selection, but these settings change the "
            "selection or the maps -- rerun part 1 (part='selection' or 'both'), or put them "
            "back:\n  " + "\n  ".join(changed))
    return cfg


def _open_log(paths, name):
    old_cwd = os.getcwd()
    os.chdir(paths["root"])
    tee = _Tee(os.path.join(paths["root"], name))
    sys.stdout = tee
    return tee, old_cwd


def _close_log(tee, old_cwd, paths, name):
    sys.stdout = tee.stdout
    tee.close()
    try:
        os.chdir(old_cwd)
    except OSError:
        pass
    print(f"\n[log] {os.path.join(paths['root'], name)}")


def run(sri, bids_configs, config=None, df_qc_by_label=None, part="both"):
    """Run the Figure 5 pipeline.

    sri          : SBA_Region_Profile_Index module -- the only external
                   dependency (AFNI / nibabel / atlas work). Must be the
                   version patched by patch_SBA_unit_level.py.
    bids_configs : list of per-dataset dicts (label, species, bids_dir,
                   label_path, seeds, optional exclude_lateralized, state).
    config       : dict overriding CONFIG entries
    df_qc_by_label : {label: df_qc slice} to enable the external QC check
    part         : "both" (default) = selection then figures;
                   "selection" = part 1 only (ends with 07_state/selection_state.pkl);
                   "figures" = part 2 only, from that file (only FIGURE_KEYS may
                   differ from the selection run; bids_configs and df_qc_by_label
                   are then taken from the file).
    """
    if part not in PARTS:
        raise ValueError(f"part must be one of {PARTS}, got {part!r}")
    if part == "figures":
        root = (config or {}).get("OUTPUT_ROOT", CONFIG["OUTPUT_ROOT"])
        state = load_state(root)
        cfg = _merge_figure_config(state["cfg"], config)
        paths = _tree(root)
        archived = _archive_previous_outputs(paths, figures_only=True)
        tee, old_cwd = _open_log(paths, "00_log_figures.txt")
        if archived:
            print(f"  Previous figures moved to {archived}")
        print(f"  part 2 (figures) from {_state_path(root)} (saved {state.get('saved')})")
        _load_time_tree(cfg, state["bids_configs"])
        if hasattr(sri, "set_profile_method"):
            sri.set_profile_method(cfg["PROFILE_METHOD"], cfg["TAIL_Q"], keep_sketch=False)
            if cfg["PROFILE_METHOD"] != "median_r":
                sri.set_tail_reference(state["anat_common"])
        try:
            return _run_figures(sri, dict(state, cfg=cfg, paths=paths))
        finally:
            if hasattr(sri, "set_profile_method"):
                sri.set_profile_method("median_r")
            _close_log(tee, old_cwd, paths, "00_log_figures.txt")

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
    tee, old_cwd = _open_log(paths, "00_log.txt")
    if archived:
        print(f"  Previous run's figures/tables/log moved to {archived} "
              f"(AFNI maps and *.json balance records left in place for reuse)")
    print(f"  part = {part!r}")
    print(f"  AFNI statmaps stored under RESULTS_BASE_DIR={cfg['RESULTS_BASE_DIR']!r}, "
          f"suffix tag MAP_TAG={cfg['MAP_TAG']!r}, "
          f"{'THRESHOLDED' if cfg['STAT_MAP_USE_THRESHOLDED'] else 'unthresholded'} "
          f"group maps")
    for w in config_warnings:
        print(f"  [CONFIG WARNING] {w}")
    _load_time_tree(cfg, bids_configs)
    # The profile method is a module-level setting of sri: set it for this run
    # and put it back afterwards, so the next run of the same process (control,
    # matrix) starts clean. The voxel values of every region are kept until the
    # regional indices of every run (6d) and F22 are computed.
    if hasattr(sri, "set_profile_method"):
        sri.set_profile_method(cfg["PROFILE_METHOD"], cfg["TAIL_Q"],
                               keep_sketch=cfg.get("PROFILE_KIND", "seed") == "seed")
    print(f"  Profile method: {cfg['PROFILE_METHOD']}"
          + (f" (q = {cfg['TAIL_Q']:g}, value column {cfg['VALUE_COL']})"
             if cfg["PROFILE_METHOD"] != "median_r" else ""))
    try:
        state = _run_selection(sri, cfg, paths, bids_configs, df_qc_by_label)
        save_state(state)
        if part == "selection":
            _write_methods(cfg, paths, sorted(state["profiles_by_label"]),
                           state["n_regions"], state["unit_counts"],
                           final_round=state["final_round"])
            print(f"\n{'#' * 70}\n  PART 1 DONE -- {paths['root']}\n  figures: "
                  f"run(..., part='figures')\n{'#' * 70}")
            return state
        return _run_figures(sri, state)
    finally:
        if hasattr(sri, "set_profile_method"):
            sri.set_profile_method("median_r")
        _close_log(tee, old_cwd, paths, "00_log.txt")


def _run_selection(sri, cfg, paths, bids_configs, df_qc_by_label):
    """Part 1. Returns the state dict that part 2 (_run_figures) needs.

    Two region spaces when the seed region is treated differently by the
    selection and by the analyses (SEED_REGION_EXCLUDE_SELECTION !=
    SEED_REGION_EXCLUDE_ANALYSES):
      selection space -- the profiles the selection scores (regions
          `sel_regions`, two-channel profiles standardised over them):
          invariance, consensus, scores of every round, F3/F4/F5, balanced
          stats. Saved in state['selection_space'].
      analysis space  -- every profile used after the selection (regions
          `common`, standardised over them): the state's usual keys
          (profiles_by_label, runs_by_label, afni_profiles, group-map
          profiles, common, species_profs, final_criterion...).
    Both come from the same voxels of the same runs; only the region set (and
    with two channels, the per-channel standardisation over it) differs. With
    the same setting for both, the two spaces are one and
    state['selection_space'] is None."""

    unit = cfg["UNIT"]
    two_ch = cfg["PROFILE_METHOD"] != "median_r"
    # Until the common regions are known the profiles are the usual median_r
    # ones (carrying, in addition, each region's voxel values when a
    # two-channel method or F22 needs them): run-level QC and the common-region
    # search use median_r, then every profile is converted (see "two-channel
    # profiles" below).
    pre_vc = "median_r" if two_ch else cfg["VALUE_COL"]
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
                    value_col=pre_vc, metric=cfg["METRIC"],
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
            afni_profiles, value_col=pre_vc)
        common_source = "group maps (atlas level + presence in the z-map)"
    else:
        # No group map available (STEP 1 skipped or every map failed):
        # fall back to the per-run profiles so the pipeline still runs.
        common = fa.get_common_regions_from_profiles(
            profiles_by_label, value_col=pre_vc,
            require_nonnan=cfg["COMMON_REGIONS_REQUIRE_NONNAN"])
        common_source = (f"per-run profiles, fallback "
                         f"(require_nonnan={cfg['COMMON_REGIONS_REQUIRE_NONNAN']})")
    if cfg["COMMON_REGIONS_OVERRIDE"]:
        wanted = [str(c).lower().strip() for c in cfg["COMMON_REGIONS_OVERRIDE"]]
        present = set(fa.get_common_regions_from_profiles(
            profiles_by_label, value_col=pre_vc, require_nonnan=False, verbose=False))
        absent = [c for c in wanted if c not in present]
        if absent:
            print(f"  [WARN] COMMON_REGIONS_OVERRIDE: {absent} absent from the profiles -- dropped")
        common = [c for c in wanted if c in present]
        common_source = "COMMON_REGIONS_OVERRIDE (regions of the main run)"
    print(f"  {len(common)} common region(s) -- source: {common_source}")
    if not common:
        raise ValueError("No common region survived -- nothing can be compared. "
                         "Check ATLAS_LEVEL and the group maps' coverage.")
    # ---- the region that contains the seed ---------------------------------
    # Its value is dominated by the seed's own autocorrelation (a peak whose
    # size depends on resolution, smoothing and on how much of the region the
    # seed fills), not by coupling. SEED_REGION_EXCLUDE finds it; two switches
    # decide where it is removed:
    #   SEED_REGION_EXCLUDE_SELECTION -- from the profiles the selection scores
    #     (per-channel standardisation, invariance ranking, consensus, score);
    #   SEED_REGION_EXCLUDE_ANALYSES  -- from every profile used after the
    #     selection (statistics and figures of part 2).
    # common_regions.xlsx keeps every region with its flags, so the matrix
    # analysis (no seed) keeps it as a node.
    all_common = list(common)
    seed_regions = _seed_regions_to_exclude(sri, cfg, bids_configs, profiles_by_label,
                                            all_common, paths)
    excl_sel = bool(cfg["SEED_REGION_EXCLUDE_SELECTION"]) and bool(seed_regions)
    excl_ana = bool(cfg["SEED_REGION_EXCLUDE_ANALYSES"]) and bool(seed_regions)
    sel_regions = [c for c in all_common if not (excl_sel and c in seed_regions)]
    common = [c for c in all_common if not (excl_ana and c in seed_regions)]
    split = sel_regions != common
    cfg["_seed_regions_found"] = list(seed_regions)
    cfg["_seed_regions_excluded"] = list(seed_regions) if excl_ana else []
    cfg["_seed_regions_excluded_selection"] = list(seed_regions) if excl_sel else []
    if seed_regions:
        print(f"  seed region(s) {seed_regions}: "
              f"{'removed from' if excl_sel else 'KEPT in'} the selection "
              f"({len(sel_regions)} regions), {'removed from' if excl_ana else 'kept in'} "
              f"the analyses ({len(common)} regions)")
    pd.DataFrame({"common_region": all_common, "source": common_source,
                  "seed_region": [c in seed_regions for c in all_common],
                  "used_by_selection": [c in sel_regions for c in all_common],
                  "used_by_analyses": [c in common for c in all_common]}).to_excel(
        os.path.join(paths["selection"], "common_regions.xlsx"), index=False)

    # ---- regional indices of every run (6d) --------------------------------
    # While every profile still carries its voxel values: median, upper / lower
    # percentile, % in the top / bottom tail of the map, share > 0 / < 0, per
    # region, both hemispheres pooled (analysis regions). Part 2 draws 6d for
    # every index from this table, without the voxel values.
    run_index = _run_region_indices(sri, runs_by_label, common, cfg, species_of_label)

    # ---- F22: reliability of the three profile constructions ---------------
    # On EVERY run, before any selection, whatever PROFILE_METHOD: the
    # criterion that decides which profile to use, independent of the
    # cross-species questions (analysis regions).
    reliability = None
    if cfg["RUN_PROFILE_RELIABILITY"]:
        print(f"\n{'=' * 70}\n  F22: profile reliability by method (every run)\n{'=' * 70}")
        try:
            resolvers = {c["label"]: ff.StateResolver(c) for c in bids_configs}
            rel = ff.profile_reliability(sri, runs_by_label, common, resolvers,
                                         cfg["ZSCORE"], q=cfg["TAIL_Q"])
            if not rel.empty:
                tests = {m: ff.compare_methods(rel, m) for m in
                         ("discriminability", "gap_animal", "r_inter_run")}
                with pd.ExcelWriter(os.path.join(paths["selection"],
                                                 "profile_reliability.xlsx")) as xw:
                    rel.to_excel(xw, sheet_name="per_dataset", index=False)
                    pd.DataFrame([dict(measure=k, method=m, **t) for k, d in tests.items()
                                  for m, t in d.items()]).to_excel(
                        xw, sheet_name="method_vs_median_r", index=False)
                reliability = dict(rel=rel, tests=tests)
                for k, d in tests.items():
                    for m, t in d.items():
                        print(f"  {k}: {m} - median_r = {t['mean_diff']:+.3f} "
                              f"({t['n_better']}/{t['n']} datasets higher), p = {t['p']:.3f}")
            else:
                print("  [F22] no dataset with several runs per animal -- skipped")
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"  [F22] [WARN] failed: {e}")
    # the voxel values are no longer needed by median_r
    if not two_ch:
        sri.set_profile_method("median_r", cfg["TAIL_Q"], keep_sketch=False)
        for profs in runs_by_label.values():
            for k in list(profs):
                if "_sketch" in profs[k].columns:
                    profs[k] = profs[k].drop(columns="_sketch")
        for lab in list(profiles_by_label):
            profiles_by_label[lab], members_by_label[lab] = fu.aggregate_profiles_to_unit(
                sri, runs_by_label[lab], unit=unit, label=lab, verbose=False)
        afni_profiles = {lab: (sp_, p.drop(columns="_sketch") if "_sketch" in p else p)
                         for lab, (sp_, p) in afni_profiles.items()}

    # ---- two-channel profiles ---------------------------------------------
    # fraction_posneg takes the map percentiles over the voxels of the common
    # regions, and both two-channel methods standardise each channel over them:
    # only now can the profiles be built. Every profile computed so far is
    # converted; every profile computed from here on (group maps of STEP 3, 3b)
    # comes out converted from sri.compute_region_profile with the LAST tail
    # reference set -- the analysis regions. With two region spaces the
    # selection-space profiles are converted first (reference = sel_regions),
    # then the analysis-space ones (reference = common).
    profiles_sel = profiles_by_label        # one space unless split (below)
    if two_ch:
        print(f"\n{'=' * 70}\n  {cfg['PROFILE_METHOD']} profiles (q = {cfg['TAIL_Q']:g}, "
              f"{len(common)} common regions"
              + (f"; selection: {len(sel_regions)} regions" if split else "")
              + f")\n{'=' * 70}")
        sri.set_profile_method(cfg["PROFILE_METHOD"], cfg["TAIL_Q"], keep_sketch=False)
        if split:
            sri.set_tail_reference(sel_regions)
            runs_sel = {lab: {k: sri.to_channel_profile(p) for k, p in profs.items()}
                        for lab, profs in runs_by_label.items()}
            profiles_sel = {}
            for lab in list(profiles_by_label):
                profiles_sel[lab], _ = fu.aggregate_profiles_to_unit(
                    sri, runs_sel[lab], unit=unit, label=lab, verbose=False)
            del runs_sel
        sri.set_tail_reference(common)
        runs_by_label = {lab: {k: sri.to_channel_profile(p) for k, p in profs.items()}
                         for lab, profs in runs_by_label.items()}
        for lab in list(profiles_by_label):
            profiles_by_label[lab], members_by_label[lab] = fu.aggregate_profiles_to_unit(
                sri, runs_by_label[lab], unit=unit, label=lab, verbose=False)
        afni_profiles = {lab: (sp_, sri.to_channel_profile(p))
                         for lab, (sp_, p) in afni_profiles.items()}
        anat = list(common)
        common = sri.channel_names(anat)
        common_sel = sri.channel_names(sel_regions) if split else common
        n_runs = sum(len(v) for v in runs_by_label.values())
        print(f"  {n_runs} run profile(s) and {len(afni_profiles)} group-map profile(s) "
              f"converted: {len(common)} values per profile ({len(anat)} regions x "
              f"{'/'.join(sri.channel_suffixes())}), column {cfg['VALUE_COL']}"
              + (f"; selection profiles: {len(common_sel)} values ({len(sel_regions)} regions, "
                 f"channels standardised over them)" if split else ""))
    else:
        common_sel = list(sel_regions) if split else common
    anat_common = list(anat) if two_ch else list(common)
    n_regions = len(anat_common)
    _write_methods(cfg, paths, sorted(profiles_by_label), n_regions, unit_counts)

    # ---- invariance (selection space) -----------------------------------
    # With INVARIANCE_VALUE=1.0 (the default) every common region is kept
    # and the mask is all-True: the invariance score is then purely
    # descriptive -- the bottom panel of F0 still ranks the regions, but
    # nothing is excluded from the run selection. Lower the value to
    # restrict the selection to the most invariant fraction.
    print(f"\n{'=' * 70}\n  Cross-species invariant regions\n{'=' * 70}")
    species_profs_sel = {lab: fa.species_profile(
                            sri, profs, source=cfg["CONSENSUS_SOURCE"],
                            method=cfg["REFERENCE_METHOD"], k=cfg["K_REF"],
                            value_col=cfg["VALUE_COL"], metric=cfg["METRIC"])
                         for lab, profs in profiles_sel.items()}
    invariance_sel, inv_mat, inv_species = fa.compute_region_invariance(
        species_profs_sel, common_sel, value_col=cfg["VALUE_COL"],
        species_of_label=species_of_label, zscore=cfg["ZSCORE"])
    region_mask_sel = fa.select_invariant_regions(
        invariance_sel, common_sel, mode=cfg["INVARIANCE_MODE"], value=cfg["INVARIANCE_VALUE"])
    print(f"  {int(np.asarray(region_mask_sel, dtype=bool).sum())}/{len(common_sel)} "
          f"region(s) used for the run selection "
          f"(mode={cfg['INVARIANCE_MODE']}, value={cfg['INVARIANCE_VALUE']})")
    pd.DataFrame({"region": common_sel, "invariance": invariance_sel,
                  "used_for_selection": region_mask_sel}).to_excel(
        os.path.join(paths["selection"], "region_invariance.xlsx"), index=False)
    if split:
        # analysis space: the species profiles recomputed on its profiles; the
        # invariance and the selection mask are those of the selection, carried
        # over by name (they say which values the selection used)
        species_profs = {lab: fa.species_profile(
                            sri, profs, source=cfg["CONSENSUS_SOURCE"],
                            method=cfg["REFERENCE_METHOD"], k=cfg["K_REF"],
                            value_col=cfg["VALUE_COL"], metric=cfg["METRIC"])
                         for lab, profs in profiles_by_label.items()}
        inv_by = dict(zip(common_sel, np.asarray(invariance_sel, float)))
        mask_by = dict(zip(common_sel, np.asarray(region_mask_sel, bool)))
        invariance = np.array([inv_by.get(c, np.nan) for c in common], float)
        region_mask = np.array([bool(mask_by.get(c, False)) for c in common])
    else:
        species_profs, invariance, region_mask = (species_profs_sel, invariance_sel,
                                                  region_mask_sel)

    # ---- STEP 2: selection (selection space) ----------------------------
    print(f"\n{'=' * 70}\n  STEP 2: consensus selection ({unit} level)\n{'=' * 70}")
    sel = fa.run_consensus_selection(
        sri, cfg, profiles_sel, common_sel,
        n_select=cfg["N_SELECT"], leave_one_out=cfg["LEAVE_ONE_OUT"],
        consensus_source=cfg["CONSENSUS_SOURCE"], similarity=cfg["SIMILARITY"],
        method=cfg["REFERENCE_METHOD"], k=cfg["K_REF"], zscore=cfg["ZSCORE"],
        threshold_abs=cfg["THRESHOLD_ABS"], species_of_label=species_of_label,
        region_mask=region_mask_sel,
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
    # final_criterion_sel[label] = profile_df the FINAL round actually scored
    # that BIDS's runs against (selection space).
    final_criterion_sel = {lab: _vec_to_profile(common_sel, vec, cfg["VALUE_COL"])
                           for lab, vec in sel["consensus_by_label"].items()}
    selections_round2 = {}
    overlap_df = pd.DataFrame()
    if cfg["RUN_ROUND2_SELECTION"]:
        print(f"\n{'=' * 70}\n  STEP 2b: round-2 selection (species-anchored)\n{'=' * 70}")
        r2 = fa.run_round2_selection(
            sri, cfg, profiles_sel, selections, common_sel,
            n_select=cfg["N_SELECT"], region_mask=region_mask_sel)
        selections_round2 = r2["selections"]
        for label in set(selections) - set(selections_round2):
            print(f"  [WARN] {label}: no round-2 selection -- keeps its round-1 selection")
            selections_round2[label] = selections[label]
        all_round_selections[2] = selections_round2
        final_criterion_sel = {
            lab: (_vec_to_profile(common_sel, r2["criterion_by_label"][lab], cfg["VALUE_COL"])
                  if lab in r2["criterion_by_label"] else final_criterion_sel.get(lab))
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
                sri, cfg, profiles_sel, selections, common_sel,
                n_select=cfg["N_SELECT"], region_mask=region_mask_sel)
            if not overlap_df.empty:
                overlap_df.drop(columns=["kcore_only", "uniform_only"]).to_excel(
                    os.path.join(paths["selection"], "reference_method_overlap.xlsx"),
                    index=False)

    start_round = max(all_round_selections)
    last_iter = None
    if cfg["RUN_ITERATIVE_SELECTION"] and cfg["N_ITERATIVE_ROUNDS"] > start_round:
        print(f"\n{'=' * 70}\n  Iterative selection: rounds {start_round + 1}.."
              f"{cfg['N_ITERATIVE_ROUNDS']} (continuing from round {start_round})\n{'=' * 70}")
        iter_res = fa.run_iterative_selection(
            sri, cfg, profiles_sel, common_sel, n_select=cfg["N_SELECT"],
            start_selections=all_round_selections[start_round], start_round=start_round,
            n_rounds=cfg["N_ITERATIVE_ROUNDS"], species_of_label=species_of_label,
            region_mask=region_mask_sel)
        for r, rd in iter_res["rounds"].items():
            all_round_selections[r] = rd["selections"]
        last_iter = max(iter_res["rounds"]) if iter_res["rounds"] else None
        if last_iter is not None:
            cons_last = iter_res["rounds"][last_iter]["consensus_by_label"]
            final_criterion_sel = {
                lab: (_vec_to_profile(common_sel, cons_last[lab], cfg["VALUE_COL"])
                      if lab in cons_last else final_criterion_sel.get(lab))
                for lab in all_round_selections[last_iter]}

    final_round = max(all_round_selections)
    selections_final = all_round_selections[final_round]
    print(f"\n  Selection chain: rounds {sorted(all_round_selections)} -- FINAL round = "
          f"{final_round}" + ("" if final_round > 1 else " (no refinement round ran)"))

    # F8b / 6c (analyses): the final round's criterion in the analysis space --
    # same construction, same kept runs at every round, analysis profiles and
    # regions (identical to the selection's own criterion when there is one space)
    if split:
        print(f"\n  final-round criterion rebuilt in the analysis regions ({len(common)} "
              f"values) for F8b / 6c")
        final_criterion = _criterion_in_space(
            sri, cfg, profiles_by_label, common, region_mask, species_profs,
            species_of_label, all_round_selections, round2_ran=bool(selections_round2),
            last_iter_round=last_iter)
    else:
        final_criterion = final_criterion_sel

    # Convergence across the WHOLE chain (round r vs r-1, r >= 2) + per-round table
    conv_df = pd.DataFrame()
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
        print(f"  [saved] 02_selection/selection_by_round.xlsx, selection_convergence.xlsx")

    balanced = fa.balanced_kept_vs_notkept(
        profiles_sel, selections, common_sel,
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
                profiles_sel, selections, df_qc_by_label,
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

    # ---- profiles of every group map (part 2 reads no NIfTI) ---------------
    # Every group of 03_group_maps and 03b_group_maps_final ('kept', 'not_kept',
    # 'kept_final', 'not_kept_final'), plus the all-units maps ('all'): they feed
    # the phylogenetic sensitivity (PHYLOGENY_USE_FINAL_ROUND), F17 and 5d
    # (before / after), STEP 8 (distance to the kept map) -- whichever part 2 uses.
    # Analysis space: the tail reference left set by the conversion is `common`.
    group_map_profiles = {"all": dict(afni_profiles)}
    for idx in (map_index, map_index_final):
        if idx is None or idx.empty:
            continue
        for g in sorted(set(idx["group"])):
            sub = idx[(idx["group"] == g) & idx["found"]]
            if sub.empty:
                continue
            print(f"\n  profiles of the '{g}' group maps ({len(sub)} dataset(s))")
            group_map_profiles[g] = fu.profiles_from_group_maps(
                sri, sub, bids_configs, None, atlas_level=cfg["ATLAS_LEVEL"],
                members_by_label=members_by_label, verbose=False)

    # ---- 5a / S5a: surface renders with one colour range per dataset --------
    surface_pics, surface_ranges = {}, {}
    if (cfg["RENDER_SURFACE"] and cfg.get("SURFACE_COMMON_SCALE", True)
            and cfg.get("PROFILE_KIND", "seed") == "seed"):
        maps_5a = fpp.surface_maps_for_5a(all_units_index, map_index_final)
        if maps_5a:
            print(f"\n{'=' * 70}\n  5a / S5a: surface renders, one colour range per "
                  f"dataset\n{'=' * 70}")
            try:
                surface_pics, surface_ranges = fpp.render_common_scale(
                    sri, os.path.join(paths["surfaces_final"], "common_scale"), maps_5a,
                    species_of_label, pct=cfg.get("SURFACE_COMMON_PCT", 90.0))
            except Exception as e:      # 5a then falls back to the per-map renders
                import traceback
                traceback.print_exc()
                print(f"  [5a] [WARN] common-scale rendering failed ({e}) -- per-map renders "
                      f"will be used")

    selection_space = None
    if split:
        selection_space = dict(
            common=common_sel, anat_common=list(sel_regions), profiles_by_label=profiles_sel,
            species_profs=species_profs_sel, invariance=invariance_sel,
            region_mask=region_mask_sel, final_criterion=final_criterion_sel)
    return dict(
        cfg=cfg, paths=paths, bids_configs=bids_configs, species_of_label=species_of_label,
        profiles_by_label=profiles_by_label, members_by_label=members_by_label,
        runs_by_label=runs_by_label, unit_counts=unit_counts, df_qc_by_label=df_qc_by_label,
        all_units_index=all_units_index, afni_profiles=afni_profiles, common=common,
        anat_common=anat_common, all_common=all_common, common_source=common_source,
        n_regions=n_regions, species_profs=species_profs, invariance=invariance,
        region_mask=region_mask, sel=sel, selections=selections,
        all_round_selections=all_round_selections, final_round=final_round,
        selections_final=selections_final, final_criterion=final_criterion,
        conv_df=conv_df, overlap_df=overlap_df, map_index=map_index,
        map_index_final=map_index_final, group_map_profiles=group_map_profiles,
        surface_pics=surface_pics, surface_ranges=surface_ranges, run_index=run_index,
        reliability=reliability, seed_regions=list(seed_regions),
        selection_space=selection_space)


def _criterion_in_space(sri, cfg, profiles_by_label, common, region_mask, species_profs,
                        species_of_label, all_round_selections, round2_ran, last_iter_round):
    """The final round's selection criterion of every dataset, rebuilt in
    another region space (the analysis space when the seed region is kept for
    the selection and removed from the analyses): same construction and the
    same kept units at every round as the selection actually used, with the
    profiles, regions and selection mask of that space. Mirrors how
    _run_selection updates the criterion round after round: round 1 (LOSO
    consensus of the species profiles), round 2 if it ran (each dataset's own
    round-1 kept reference), then the last iterative round if any (consensus
    rebuilt from the previous round's kept references)."""
    vc = cfg["VALUE_COL"]
    mask = np.asarray(region_mask, dtype=bool)
    _, _, cons1, _ = fa._select_against_species_profiles(
        sri, species_profs, profiles_by_label, common, cfg["N_SELECT"],
        cfg["LEAVE_ONE_OUT"], cfg["ZSCORE"], cfg["THRESHOLD_ABS"], cfg["SIMILARITY"],
        species_of_label, mask, vc, cfg["METRIC"], False)
    crit = {lab: _vec_to_profile(common, v, vc) for lab, v in cons1.items()}
    if round2_ran:
        r2 = fa.run_round2_selection(sri, cfg, profiles_by_label, all_round_selections[1],
                                     common, n_select=cfg["N_SELECT"], region_mask=mask,
                                     verbose=False)
        crit = {lab: (_vec_to_profile(common, r2["criterion_by_label"][lab], vc)
                      if lab in r2["criterion_by_label"] else crit.get(lab))
                for lab in all_round_selections[2]}
    if last_iter_round is not None:
        it = fa.run_iterative_selection(
            sri, cfg, profiles_by_label, common, n_select=cfg["N_SELECT"],
            start_selections=all_round_selections[last_iter_round - 1],
            start_round=last_iter_round - 1, n_rounds=last_iter_round,
            species_of_label=species_of_label, region_mask=mask, verbose=False)
        cons = it["rounds"][last_iter_round]["consensus_by_label"]
        crit = {lab: (_vec_to_profile(common, cons[lab], vc) if lab in cons else crit.get(lab))
                for lab in all_round_selections[last_iter_round]}
    return crit


def _write_methods_article(st, cfg, paths, phylo_source=None):
    """Methods_article.txt: the Methods of Figures 5-6 written as for a
    manuscript, with the numbers of this run (datasets, runs, animals,
    regions, selection, tests), followed by the references cited. Rewritten
    at the end of every part 2."""
    import scipy
    import matplotlib
    unit = cfg["UNIT"]
    labels = sorted(st["profiles_by_label"])
    species_of = st["species_of_label"]
    cfg_by = {c["label"]: c for c in st["bids_configs"]}
    rows = []
    for lab in labels:
        res = ff.StateResolver(cfg_by.get(lab, {}))
        keys = list(st["runs_by_label"].get(lab, {}))
        animals = {res.animal(tuple(str(x) for x in k)) for k in keys}
        try:
            state = ff.dataset_state(res, keys)
        except Exception:
            state = "unknown"
        rows.append((lab, species_of.get(lab), len(keys), len(animals), state))
    n_species = len({r[1] for r in rows})
    n_runs = sum(r[2] for r in rows)
    n_anim = sum(r[3] for r in rows)
    ds_txt = "; ".join(f"{ff.pretty_label(lab)} ({sp}, {nr} runs, {na} animals, {state})"
                       for lab, sp, nr, na, state in rows)
    method = cfg["PROFILE_METHOD"]
    q = float(cfg["TAIL_Q"])
    lo, hi = int(round(100 * q)), int(round(100 * (1 - q)))
    anat = st.get("anat_common") or st["common"]
    seed_found = (st.get("seed_regions") or cfg.get("_seed_regions_found")
                  or cfg.get("_seed_regions_excluded") or [])
    out_ana = bool(cfg.get("_seed_regions_excluded"))
    out_sel = bool(cfg.get("_seed_regions_excluded_selection",
                           cfg.get("_seed_regions_excluded")))
    if seed_found and out_ana and out_sel:
        seed_txt = ("; the region containing the seed (" + ", ".join(seed_found) + "), found "
                    "from the seed voxels in each atlas, was removed from every profile because "
                    "its value reflects the autocorrelation of the seed rather than coupling")
    elif seed_found and out_ana:
        seed_txt = ("; the region containing the seed (" + ", ".join(seed_found) + "), found "
                    "from the seed voxels in each atlas, was kept in the profiles scored by the "
                    "selection and removed from every subsequent analysis, because its value "
                    "reflects the autocorrelation of the seed rather than coupling")
    elif seed_found and out_sel:
        seed_txt = ("; the region containing the seed (" + ", ".join(seed_found) + ") was "
                    "removed from the profiles scored by the selection only")
    else:
        seed_txt = ""
    # the selection score: regions and mask of the selection space
    ss = st.get("selection_space") or dict(common=st["common"], anat_common=anat,
                                           region_mask=st.get("region_mask"))
    n_val = len(ss["common"])
    n_inv = (int(np.asarray(ss["region_mask"], bool).sum())
             if ss.get("region_mask") is not None else n_val)
    sel_anat = ss.get("anat_common") or anat
    inv_txt = (f"{n_inv} of the {n_val} profile values ({len(sel_anat)} regions x 2 channels)"
               if n_val != len(sel_anat) else f"{n_inv} of the {n_val} regions")
    n_kept = {lab: len(v) for lab, v in st["selections_final"].items()}
    if method == "percentiles":
        prof = (f"two values per region: the {hi}th and the {lo}th percentiles of the "
                f"region's voxel values (its most and least coupled {lo} %), left and right "
                f"hemispheres pooled by voxel count. Each of the two channels was standardised "
                f"across the common regions (mean 0, s.d. 1), so that the Pearson correlation "
                f"between two profiles is the mean of the two channel correlations and a global "
                f"offset between channels cannot create similarity")
    elif method == "fraction_posneg":
        prof = (f"two values per region: the percentage of the region's voxels in the top "
                f"{lo} % and in the bottom {lo} % of the map, the map percentiles being taken "
                f"over the voxels of the common regions (both hemispheres); rank-based, hence "
                f"invariant to any monotonic transform of the map. Each channel was standardised "
                f"across the common regions")
    else:
        prof = "the median of the region's voxel values"
    norm = {True: "z-scored across regions", "zscore": "z-scored across regions",
            "scale": "divided by their root mean square across regions (sign preserved)",
            False: "left on their raw scale", None: "left on their raw scale",
            "none": "left on their raw scale"}.get(cfg["ZSCORE"], str(cfg["ZSCORE"]))
    rounds = st["final_round"]
    chain = (f"Round 1 ranked the {unit}s of each dataset by the Pearson correlation of their "
             f"profile with the leave-one-species-out consensus (median across the other "
             f"species of the species profiles; datasets of one species first collapsed by "
             f"median, so that each species counts once)"
             + ("; round 2 re-ranked each dataset against the median of its own round-1 "
                "selection" if cfg["RUN_ROUND2_SELECTION"] else "")
             + ((f"; round 3 rebuilt" if rounds == 3 else f"; rounds 3 to {rounds} each rebuilt")
                + f" the cross-species consensus from the previous round's selected {unit}s and "
                  f"re-ranked the whole pool" if rounds >= 3 else "")
             + f". The {cfg['N_SELECT']} best-ranked {unit}s of each dataset were kept at every "
             f"round (final selection: round {rounds}).")
    versions = (f"Python {sys.version.split()[0]}, NumPy {np.__version__}, pandas "
                f"{pd.__version__}, SciPy {scipy.__version__}, Matplotlib "
                f"{matplotlib.__version__}")
    g6 = " and ".join((cfg.get("FIG6D") or {}).get("group", ("Human", "Macaque")))
    regs6 = ", ".join(ff.short_region(r) for r in (cfg.get("FIG6D") or {}).get(
        "regions", ()))
    pairs5b = ("pairs of runs from the same species but two different datasets (available for "
               "the species with two datasets) and pairs from different species; pairs from the "
               "same dataset, which also share protocol, scanner and usually the animal, were "
               "not compared"
               if cfg.get("WITHIN_ACROSS_PAIRS", "other_dataset") == "other_dataset" else
               "intra-species pairs (same species, same or another dataset) and inter-species "
               "pairs")
    test5b = ("one-sided Wilcoxon rank-sum test on the pairs, not corrected for multiple "
              "comparisons (single primary hypothesis); because every run enters many pairs, "
              "the pairs are not independent and this test is liberal, so a permutation test "
              "of the species labels across datasets (each dataset keeping its runs; 2,000 "
              "permutations; statistic, difference between the intra- and inter-species "
              "medians), which respects that dependence, is reported with it"
              if cfg.get("WITHIN_ACROSS_TEST", "wilcoxon") == "wilcoxon" else
              "permutation test of the species labels across datasets (each dataset keeping its "
              "runs; 2,000 permutations; statistic, difference between the intra- and "
              "inter-species medians), one-sided; a Wilcoxon rank-sum test on the pairs is "
              "reported for comparison")
    txt = f"""METHODS -- Cross-species retrosplenial connectivity profiles (Figures 5 and 6)
(written by fig5_pipeline from the run in {paths['root']}; numbers are those of this run)

Datasets
{len(labels)} resting-state fMRI datasets from {n_species} species ({n_runs} runs, {n_anim} animals) were analysed: {ds_txt}. Every run was preprocessed and registered to its species template with EDNiX, and a seed-based correlation map (Pearson r) of the retrosplenial seed ({cfg['SEED']}) was computed in template space.

Regional profiles
Each run's map was summarised over the regions of the EDNiX atlas at level {cfg['ATLAS_LEVEL']} into a regional profile: {prof}. Only regions present in the atlas of every dataset and sampled by every dataset's group map were kept ({len(anat)} common regions in the analyses){seed_txt}. {'Runs were screened beforehand by their distance to the most central run of their dataset (median absolute deviation threshold ' + str(cfg['RUN_QC_MAD_THRESHOLD']) + ').' if cfg['RUN_QC_ENABLE'] else 'No run was excluded before the selection.'} Profiles were {norm} before any comparison. The unit of the selection and of the group maps was the {unit}.

Consensus-based selection
{chain} The score was computed on the {inv_txt} that were the most invariant across species (top fraction {cfg['INVARIANCE_VALUE']:g}, ranked by the median absolute deviation of the species values); every common region was used everywhere else. Leave-one-species-out ensured that no dataset contributed to the criterion that selected its own runs. Selected per dataset: {', '.join(f'{ff.pretty_label(k)} {v}' for k, v in sorted(n_kept.items()))}.

Group maps and surfaces
For each dataset, group maps were computed with AFNI 3dttest++ (Cox 1996) on the Fisher-z transformed run maps of (i) every {unit}, (ii) the {unit}s of the final selection and (iii) the same number of {unit}s drawn at random among those not selected. They were rendered on each species' inflated cortical surface with Connectome Workbench (wb_command; Marcus et al. 2011), with one symmetric colour range per dataset shared by its three maps (± the {cfg.get('SURFACE_COMMON_PCT', 90):g}th percentile of |z| over their non-zero voxels).

Reliability of the profile constructions
The three constructions of the profile (median; upper and lower percentiles; tail fractions) were compared on every run before the selection, within the datasets holding several runs per animal, by their discriminability (Bridgeford et al. 2021): the probability that a run is closer (1 − r) to another run of the same animal than to a run of a different animal (0.5, chance). Its significance against chance was assessed by permuting the animal labels across runs (500 permutations), and the methods were compared with an exact paired sign-flip test across datasets. The criterion was fixed before, and independently of, the cross-species analyses.

Statistics
Figure 5b. Pearson correlations between the profiles of every pair of {unit}s were split into {pairs5b}, among every {unit} and among the selected {unit}s; the two distributions were compared with a {test5b}.
Figure 5c. The agreement between species profiles was the mean Pearson r over all pairs of species, tested against {10000:,} permutations of the region labels within each species (one-sided); this permutation ignores the spatial autocorrelation between neighbouring regions and is therefore liberal.
Figure 5d. Whether the selection reorders the datasets was tested by the Spearman correlation between the dataset-to-consensus correlations before and after selection ({10000:,} permutations, two-sided).
Figure 6a. Pearson correlations between the species-mean similarity to the leave-one-species-out consensus (or to Human) and divergence time (time-calibrated tree{': ' + os.path.basename(cfg['TREE_PATH']) if cfg.get('TREE_PATH') else ''}) or log10 brain mass (Herculano-Houzel et al. 2015; Jardim-Messeder et al. 2017) were tested exactly over all assignments of the x values to the species (two-sided), with and without Human. Profiles: {phylo_source or 'group maps'}.
Figure 6c. Within the datasets holding both arousal states in the same animals, the per-run correlation with the selection criterion was compared between anaesthetised and awake runs by the mean over animals of the per-animal difference, tested with an exact sign-flip of the animals (2^n assignments).
Figure 6d. For each of six regions ({regs6}), the regional value of every selected run (one panel per index: upper and lower percentile, median, % of voxels in the top and bottom tail of the map, % of voxels with r > 0 and r < 0) was averaged over the runs of each animal, then over the animals of each dataset, then over the datasets of each species, so that neither an animal with more runs nor a dataset with more animals weighs more. The species is the unit of inference -- runs of one species share its brain, protocol and scanner, and testing runs would be pseudoreplication (Hurlbert 1984; Lazic 2010; Aarts et al. 2014). The difference between the mean of {g6} and the mean of the other species was tested exactly over every assignment of 2 of the 7 species to the group (21 assignments; two-sided; smallest possible P = 1/21 = 0.048), one test per region, not corrected for the number of regions. Robustness: the same test with one dataset per species, for every combination of the datasets of the species that have two. Region size: the group difference of each region was compared with the linear fit, across the other common regions, of the group difference against the group difference in log10 relative region size.
Significance marks on the figures: *** P < 0.001, ** P < 0.01, * P < 0.05, n.s. not significant; exact P values, n and statistics are given in the .txt file of each panel. Permutation P values are (1 + number of null statistics at least as extreme) / (1 + number of permutations) except for exact enumerations.

Software
{versions}; AFNI (3dttest++, 3dMean); Connectome Workbench (wb_command); nibabel.

References
Aarts E, Verhage M, Veenvliet JV, Dolan CV, van der Sluis S (2014) A solution to dependency: using multilevel analysis to accommodate nested data. Nat Neurosci 17:491-496.
Bridgeford EW, et al. (2021) Eliminating accidental deviations to minimize generalization error and maximize replicability: applications in connectomics and genomics. PLoS Comput Biol 17:e1009279.
Cox RW (1996) AFNI: software for analysis and visualization of functional magnetic resonance neuroimages. Comput Biomed Res 29:162-173.
Herculano-Houzel S, Catania K, Manger PR, Kaas JH (2015) Mammalian brains are made of these: a dataset of the numbers and densities of neuronal and nonneuronal cells in the brain of glires, primates, scandentia, eulipotyphlans, afrotherians and artiodactyls, and their relationship with body mass. Brain Behav Evol 86:145-163.
Hurlbert SH (1984) Pseudoreplication and the design of ecological field experiments. Ecol Monogr 54:187-211.
Jardim-Messeder D, et al. (2017) Dogs have the most neurons, though not the largest brain: trade-off between body mass and number of neurons in the cerebral cortex of large carnivoran species. Front Neuroanat 11:118.
Lazic SE (2010) The problem of pseudoreplication in neuroscientific studies: is it affecting your analysis? BMC Neurosci 11:5.
Mantel N (1967) The detection of disease clustering and a generalized regression approach. Cancer Res 27:209-220.
Marcus DS, et al. (2011) Informatics and data mining tools and strategies for the Human Connectome Project. Front Neuroinform 5:4.
Nili H, et al. (2014) A toolbox for representational similarity analysis. PLoS Comput Biol 10:e1003553.
"""
    path = os.path.join(paths["root"], "Methods_article.txt")
    with open(path, "w", encoding="utf-8") as f:
        f.write(txt)
    print(f"  [saved] Methods_article.txt")
    return path


def _kept_references(sri, cfg, profiles_by_label, selections, selections_final, final_round):
    """Per dataset, the reference (REFERENCE_METHOD / K_REF) of its round-1
    kept units and of its final-round kept units -- F1a/F1b, F4, fingerprints."""
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
    return kept_ref_by_label, kept_ref_final_by_label


def _f4_per_dataset(sri, cfg, fd, sel, selections, final_round, profiles_by_label, common,
                    species_profs, species_of_label, kept_ref_by_label,
                    kept_ref_final_by_label, space_note=""):
    """F4: one figure per dataset -- its own reference, the consensus it was
    scored against, every unit (kept / not kept), the kept references. Drawn in
    the region space of the SELECTION (it shows the selection's target).
    Returns own_refs ({label: profile}); empty when no reference could be built."""
    own_refs = {}
    for label, profs in profiles_by_label.items():
        try:
            own_refs[label], _ = fa.build_reference_for_profiles(
                sri, profs, method=cfg["REFERENCE_METHOD"], k=cfg["K_REF"],
                value_col=cfg["VALUE_COL"], metric=cfg["METRIC"])
        except Exception as e:
            print(f"  [WARN] {label}: reference failed: {e}")
    if not own_refs:
        return own_refs
    unit = cfg["UNIT"]
    # cons_prof must stay on the RAW value_col scale to be comparable to
    # own_refs and the individual units. The consensus in
    # sel["consensus_by_label"] is z-scored and leave-one-out (correct for the
    # selection SCORING) -- rebuild a separate display-only consensus here: all
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
    cons_prof = own_refs[next(iter(own_refs))].copy()
    cons_prof[cfg["VALUE_COL"]] = (
        cons_prof["base_name"].str.lower().str.strip()
        .map(dict(zip(common, display_consensus))))
    # Unweighted grand mean: every BIDS's own reference counted once, unlike
    # display_consensus which collapses same-species BIDS together first.
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
        if label not in own_refs:
            continue
        own_vals = own_refs[label][cfg["VALUE_COL"]].values
        print(f"  [F4] {label}: own_ref range=[{np.nanmin(own_vals):.3f}, "
              f"{np.nanmax(own_vals):.3f}]  (should be same order of "
              f"magnitude as the consensus range above)")
        kept_keys = {tuple(str(x) for x in k) for k in selections.get(label, [])}
        notkept_profs = {k: v for k, v in profs.items()
                         if tuple(str(x) for x in k) not in kept_keys}
        notkept_ref = None
        try:
            if len(notkept_profs) >= 2:
                notkept_ref, _ = fa.build_reference_for_profiles(
                    sri, notkept_profs, method=cfg["REFERENCE_METHOD"],
                    k=min(cfg["K_REF"], len(notkept_profs)),
                    value_col=cfg["VALUE_COL"], metric=cfg["METRIC"])
        except Exception as e:
            print(f"  [WARN] {label}: not-kept-only reference failed: {e}")
        fa.plot_bids_ref_with_individuals(
            profs, own_refs[label], cons_prof, common,
            selected_keys=selections.get(label),
            kept_ref=kept_ref_by_label.get(label), notkept_ref=notkept_ref,
            grand_mean_ref=grand_mean_prof,
            kept_ref_round2=kept_ref_final_by_label.get(label),
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
        f"(green = kept, grey = not kept).{space_note}",
        f"Dataset reference built with {cfg['REFERENCE_METHOD']} over ALL its "
        f"{unit}s, independent of the selection. Regions ordered by consensus "
        f"value.",
        f"Kept {unit}s should track the orange line more closely than grey "
        f"ones. If they instead look flatter, the metric is rewarding low "
        f"amplitude rather than matching shape -- try "
        f"SIMILARITY='neg_euclidean' with ZSCORE=False.")
    return own_refs


def _paper_out(paths, name):
    """06_figures/Paper/<name>: print-size panels laid out by analysis/Figure_4_*.py and
    analysis/Supplementary/Supp_Fig_07-11*.py."""
    d = os.path.join(paths["figures"], "Paper")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, name)


def _paper(fn, *a, **k):
    """A print-size panel failing must not stop the pipeline."""
    try:
        return fn(*a, **k)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"  [paper panel] [WARN] {getattr(fn, '__name__', fn)} failed: {e}")
        return None


def _run_figures(sri, st):
    """Part 2: every statistic and figure, from the part-1 state `st` only.

    Two region spaces when the seed region is kept for the selection and
    removed from the analyses (state['selection_space'] set): F3, F4 and F5
    describe the selection and are drawn in its regions; every other output
    uses the analysis regions (the state's usual keys)."""
    cfg, paths, bids_configs = st["cfg"], st["paths"], st["bids_configs"]
    unit = cfg["UNIT"]
    species_of_label = st["species_of_label"]
    profiles_by_label, members_by_label = st["profiles_by_label"], st["members_by_label"]
    runs_by_label, unit_counts = st["runs_by_label"], st["unit_counts"]
    df_qc_by_label = st["df_qc_by_label"]
    all_units_index, afni_profiles = st["all_units_index"], st["afni_profiles"]
    common, n_regions = st["common"], st["n_regions"]
    species_profs, invariance, region_mask = (st["species_profs"], st["invariance"],
                                              st["region_mask"])
    sel, selections = st["sel"], st["selections"]
    all_round_selections, final_round = st["all_round_selections"], st["final_round"]
    selections_final, final_criterion = st["selections_final"], st["final_criterion"]
    map_index, map_index_final = st["map_index"], st["map_index_final"]
    group_map_profiles = st["group_map_profiles"]
    # ---- region space of the selection vs of the analyses -------------------
    split = st.get("selection_space") is not None
    ss = st["selection_space"] if split else dict(
        common=common, profiles_by_label=profiles_by_label, species_profs=species_profs,
        invariance=invariance, region_mask=region_mask, final_criterion=final_criterion)
    # species profiles of every analysis (F1, F1c, fingerprints, 5c): those of the
    # analysis space (with one space, the selection's own, as before)
    ana_species_profiles = species_profs if split else sel["species_profiles"]
    seed_txt = ", ".join(st.get("seed_regions") or cfg.get("_seed_regions_found") or [])
    sel_note = (f" Drawn in the regions of the SELECTION ({len(ss['common'])} values): the "
                f"region containing the seed ({seed_txt}) is "
                f"{'included' if not cfg.get('SEED_REGION_EXCLUDE_SELECTION') else 'excluded'}"
                f", as it was for the selection; the analyses use {len(common)} values."
                if split else "")
    # selections are unit keys; every run-level step needs the RUN keys behind
    # them (with UNIT="subject" a unit key (subject, 'ALL', 'ALL', seed) never
    # equals a run key, so kept flags of runs were all False before)
    kept_runs_final = fpp.selected_run_keys(selections_final, members_by_label)
    kept_runs_r1 = fpp.selected_run_keys(selections, members_by_label)
    print(f"\n{'#' * 70}\n  PART 2: statistics and figures\n{'#' * 70}")
    if split:
        print(f"  region spaces: selection {len(ss['common'])} values (seed region "
              f"{seed_txt} {'excluded' if cfg.get('SEED_REGION_EXCLUDE_SELECTION') else 'kept'})"
              f", analyses {len(common)} values (seed region "
              f"{'excluded' if cfg.get('SEED_REGION_EXCLUDE_ANALYSES') else 'kept'}); F3, F4, "
              f"F5 in the selection's regions, everything else in the analysis regions")

    # ---- figures computed in part 1, drawn here ------------------------------
    if cfg["MAKE_FIGURES"]:
        if afni_profiles:
            print(f"\n{'=' * 70}\n  F0 figures (per BIDS, then per species)\n{'=' * 70}")
            step1b_species_profiles_figure(sri, cfg, paths, afni_profiles,
                                           common, region_mask, invariance)
        if st["overlap_df"] is not None and not st["overlap_df"].empty:
            fa.plot_reference_method_overlap(
                st["overlap_df"], output_path=os.path.join(
                    paths["figures"], "F1d_reference_method_effect.png"))
        if st["conv_df"] is not None and not st["conv_df"].empty:
            fa.plot_iterative_selection_convergence(
                st["conv_df"], output_path=os.path.join(
                    paths["figures"], "F1e_selection_convergence.png"))
            if cfg.get("RUN_PAPER_FIGURES", True):
                _paper(fsp.f1e_convergence, st["conv_df"], _paper_out(paths, "F1e_convergence"))
        rel = st.get("reliability")
        if rel and not rel["rel"].empty:
            png = os.path.join(paths["figures"], "F22_profile_reliability_by_method.png")
            ff.plot_profile_reliability(rel["rel"], species_of_label, png,
                                        current_method=cfg["PROFILE_METHOD"])
            _write_profile_reliability_legend(png, rel["rel"], rel["tests"], cfg)

    # ---- phylogeny source (F6, F7, F11-F13, STEP 7) -----------------------
    # Resolved ONCE here so every phylogenetic output uses the same maps.
    phylo_profiles = afni_profiles
    phylo_index = all_units_index
    phylo_source = f"all-{unit}s group maps (no selection)"
    # Every dataset that has runs must be represented; compare against that,
    # not against afni_profiles (which can itself be incomplete).
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
            cand = {lab: v for lab, v in group_map_profiles.get("kept_final", {}).items()
                    if lab in expected}
            not_profiled = sorted(expected - set(cand))
            if not_profiled:
                print(f"  [PHYLO] [WARN] could not profile round-{final_round} kept maps "
                      f"for {not_profiled} -- keeping all-{unit}s maps")
            else:
                phylo_profiles, phylo_index = cand, kf
                phylo_source = f"round-{final_round} KEPT group maps (selection applied)"
    phylo_kind = "afni_kept_final" if phylo_profiles is not afni_profiles else "afni_all"
    if not phylo_profiles:
        # No AFNI group map at all (RUN_ALL_UNITS_MAPS=False, e.g. the matrix
        # analysis of fig5_matrices): fall back to the median over every unit
        # of each dataset -- still no selection, only a different estimator.
        phylo_profiles = {
            lab: (species_of_label.get(lab), sri.compute_group_median_profile(profs, set(profs)))
            for lab, profs in profiles_by_label.items()}
        phylo_source = f"median of all {unit} profiles (no AFNI group map, no selection)"
        phylo_kind = "python_all"
    if set(phylo_profiles) != expected:
        print(f"  [PHYLO] [WARN] phylogenetic analyses cover {len(phylo_profiles)}/"
              f"{len(expected)} datasets -- missing {sorted(expected - set(phylo_profiles))}")
    print(f"\n  [PHYLO] phylogenetic analyses use: {phylo_source}")

    # ---- STEP 6: figures + legends --------------------------------------
    fingerprints = {}
    kept_ref_by_label, kept_ref_final_by_label = {}, {}
    if cfg["MAKE_FIGURES"]:
        print(f"\n{'=' * 70}\n  STEP 6: figures\n{'=' * 70}")
        fd = paths["figures"]

        fa.plot_consensus_by_species(
            ana_species_profiles, common, species_of_label,
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

        # Kept-round-1 and kept-final per-BIDS references (analysis space), reused
        # by F1a/F1b and the fingerprints; F4 has its own, in the selection space.
        kept_ref_by_label, kept_ref_final_by_label = _kept_references(
            sri, cfg, profiles_by_label, selections, selections_final, final_round)

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
            ana_species_profiles, common, species_of_label,
            value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"],
            title=f"Consensus by BIDS -- all runs (lvl {cfg['ATLAS_LEVEL']}, {unit} level)",
            output_path=os.path.join(fd, "F1c_consensus_by_bids.png"))
        if cfg.get("RUN_PAPER_FIGURES", True):
            _paper(fsp.f1c_consensus_by_dataset, ana_species_profiles, common,
                   species_of_label, cfg["VALUE_COL"], cfg["ZSCORE"],
                   _paper_out(paths, "F1c_consensus_by_dataset"))
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

        pairs_def = cfg.get("WITHIN_ACROSS_PAIRS", "other_dataset")
        fa.plot_intra_vs_inter_species_distance(
            profiles_by_label, common, species_of_label,
            value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"],
            title=f"Within- vs across-species similarity (lvl {cfg['ATLAS_LEVEL']})",
            output_path=os.path.join(fd, "F2_within_vs_across_species.png"),
            intra=pairs_def)
        compared = ("'same species, other dataset' vs 'different species' (same-dataset "
                    "pairs shown, not compared)" if pairs_def == "other_dataset" else
                    "every same-species pair (same dataset included) vs 'different species'")
        _write_legend(
            os.path.join(fd, "F2_within_vs_across_species.txt"),
            "F2 -- Within- vs across-species similarity",
            f"Distribution of {unit}-to-{unit} profile correlation for three kinds of "
            f"pairs: same dataset; same species but another dataset; different species.",
            f"All {unit} pairs (subsampled above 20000 pairs). Classes compared "
            f"(WITHIN_ACROSS_PAIRS={pairs_def!r}): {compared}. AUC is descriptive. p comes "
            f"from permuting species labels ACROSS DATASETS (each dataset keeps its {unit}s; "
            f"pair correlations fixed): pairs share {unit}s and are not independent, so a "
            f"Mann-Whitney p on pairs would be anti-conservative by orders of magnitude."
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
            output_path=os.path.join(fd, "F2b_within_vs_across_all_vs_selected.png"),
            intra=pairs_def)
        _write_legend(
            os.path.join(fd, "F2b_within_vs_across_all_vs_selected.txt"),
            "F2b -- Within- vs across-species similarity, all runs vs selected",
            f"Same analysis as F2: left, every {unit}; right, the {unit}s kept by the "
            f"final selection round (round {final_round}). Paper panel 5b draws the two "
            f"compared classes of this figure with the same pairs, bins and statistics.",
            f"Same pair classes, AUC and dataset-level permutation test as F2 (compared: "
            f"{compared}), computed independently on each pool.",
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
        if cfg.get("RUN_PAPER_FIGURES", True):
            _paper(fsp.f3_score_distributions, sel["rankings"], selections,
                   _paper_out(paths, "F3_score_distributions"))
        _write_legend(
            os.path.join(fd, "F3_score_distributions.txt"),
            "F3 -- Selection score distributions",
            f"Per dataset, the distribution of every {unit}'s similarity to the "
            f"consensus, with the selection cutoff marked in orange.{sel_note}",
            f"Similarity metric: {cfg['SIMILARITY']}. Cutoff is the score of the "
            f"{cfg['N_SELECT']}th best {unit}.",
            f"If the kept {unit}s form a separate upper mode, the selection captures "
            f"a real subgroup. If the distribution is smooth across the cutoff, the "
            f"selection is a cut through a continuum and should be described as such. "
            f"Note the x-axis is now much shorter than at run level -- each dataset "
            f"has far fewer {unit}s than runs, so a single {unit} moves the cutoff "
            f"more.")

        # F4 shows the selection's own target: drawn in the selection's regions
        if split:
            f4_kept, f4_kept_final = _kept_references(
                sri, cfg, ss["profiles_by_label"], selections, selections_final, final_round)
        else:
            f4_kept, f4_kept_final = kept_ref_by_label, kept_ref_final_by_label
        own_refs = _f4_per_dataset(
            sri, cfg, fd, sel, selections, final_round, ss["profiles_by_label"], ss["common"],
            ss["species_profs"], species_of_label, f4_kept, f4_kept_final,
            space_note=sel_note)
        if own_refs:
            sign_df = fa.diagnose_consensus_sign(
                species_profs, common, value_col=cfg["VALUE_COL"],
                species_of_label=species_of_label)
            sign_df.to_excel(os.path.join(
                paths["group_maps"], "consensus_sign_diagnostic.xlsx"), index=False)

        fa.plot_region_invariance(
            ss["invariance"], ss["common"], mask=ss["region_mask"],
            title=f"Region invariance across species (lvl {cfg['ATLAS_LEVEL']})",
            output_path=os.path.join(fd, "F5_region_invariance.png"))
        if cfg.get("RUN_PAPER_FIGURES", True):
            _paper(fsp.f5_region_invariance, ss["invariance"], ss["common"], ss["region_mask"],
                   _paper_out(paths, "F5_region_invariance"))
        _write_legend(
            os.path.join(fd, "F5_region_invariance.txt"),
            "F5 -- Cross-species region invariance",
            "Per region, how much species AGREE on its value. Green bars are the "
            f"regions the selection was restricted to.{sel_note}",
            f"Invariance = 1 - rescaled median absolute deviation of the per-species "
            f"values, computed on the {unit}-level profiles. Selection used "
            f"mode={cfg['INVARIANCE_MODE']}, value={cfg['INVARIANCE_VALUE']}.",
            "Grey regions vary across species, so matching them would mean matching "
            "cross-species noise. If very few regions are green, the conserved core "
            "is small and any cross-species claim rests on little.")

        # --- phylogenetic regression (F6/F7) --------------------------------
        if phylo_profiles:
            print(f"  [F6] using {phylo_source} ({len(phylo_profiles)} dataset(s)), "
                  f"leave-one-SPECIES-out reference per point")
            ddf_bids = fa.compute_loo_species_distances(
                phylo_profiles, common, value_col=cfg["VALUE_COL"], zscore=cfg["ZSCORE"])
            if not ddf_bids.empty:
                ddf_bids["source"] = phylo_source
            if len(ddf_bids) >= 3:
                # One row per SPECIES for the fit: two BIDS of one species
                # are not independent observations.
                ddf = (ddf_bids.groupby("species", as_index=False)
                       .agg(distance_to_human=("distance_to_human", "mean"),
                            pearson_r=("pearson_r", "mean"),
                            n_bids=("label", "count")))
                ddf["source"] = "all_units"

                # Intra- vs inter-species contrast
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
                    # distance to the leave-one-species-out consensus.
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
                     f"the slope. " if phylo_kind == "afni_all" else
                     f"Median over EVERY {unit} profile of each dataset, no selection "
                     f"(no AFNI group map in this run). " if phylo_kind == "python_all" else
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

        # ---- F6ab, fingerprints, F15-F18 (fig5_fingerprints) ------------
        if cfg["RUN_FINGERPRINTS"] and phylo_profiles:
            print(f"\n{'=' * 70}\n  Fingerprints, two-reference F6, proximity, "
                  f"before/after, arousal state\n{'=' * 70}")
            # Before/after needs the SAME estimator on both sides: AFNI
            # all-units vs AFNI kept-final maps when every dataset has both,
            # otherwise the median of all vs kept unit profiles.
            kept_final_keys = {lab: {tuple(str(x) for x in k)
                                     for k in selections_final.get(lab, [])}
                               for lab in profiles_by_label}
            py_all = {lab: (species_of_label.get(lab),
                            sri.compute_group_median_profile(profs, set(profs)))
                      for lab, profs in profiles_by_label.items()}
            py_sel = {}
            for lab, profs in profiles_by_label.items():
                keys = {k for k in profs if tuple(str(x) for x in k) in kept_final_keys[lab]}
                if keys:
                    py_sel[lab] = (species_of_label.get(lab),
                                   sri.compute_group_median_profile(profs, keys))
            ba_all, ba_sel = py_all, py_sel
            ba_source = (f"median of {unit} profiles: all vs round-{final_round} kept")
            if (afni_profiles and set(afni_profiles) == expected and final_round >= 2
                    and not map_index_final.empty):
                kf = map_index_final[(map_index_final.group == "kept_final")
                                     & map_index_final.found]
                if set(kf["label"]) >= expected:
                    afni_sel = {lab: v for lab, v in
                                group_map_profiles.get("kept_final", {}).items()
                                if lab in expected}
                    if set(afni_sel) == expected:
                        ba_all, ba_sel = afni_profiles, afni_sel
                        ba_source = (f"AFNI group maps: all {unit}s vs round-{final_round} "
                                     f"kept")
            fingerprints = ff.make_fingerprint_figures(
                group_all=phylo_profiles, group_all_source=phylo_source,
                ba_all=ba_all, ba_sel=ba_sel, ba_source=ba_source,
                f1_all=ana_species_profiles,
                f1_sel=kept_ref_final_by_label or kept_ref_by_label,
                common=common, region_mask=region_mask, species_of_label=species_of_label,
                bids_configs=bids_configs, runs_by_label=runs_by_label,
                selections_final=kept_runs_final, final_round=final_round, cfg=cfg,
                fig_dir=fd, table_dir=paths["group_maps"],
                divergence=dict(fa.DIVERGENCE_TIME_MYA),
                ref_species=cfg["REFERENCE_SPECIES"], n_perm=cfg["N_PERM_PHYLOGENY"])

            if cfg["RUN_PAPER_FIGURES"]:
                print(f"\n{'=' * 70}\n  Paper panels: Figures 5 and 6\n{'=' * 70}")
                try:
                    fingerprints["paper"] = fpp.make_paper_figures(
                        sri=sri, fu=fu, cfg=cfg, paths=paths, bids_configs=bids_configs,
                        species_of_label=species_of_label, common=common,
                        region_mask=region_mask, final_round=final_round,
                        profiles_by_label=profiles_by_label, runs_by_label=runs_by_label,
                        members_by_label=members_by_label, selections_final=selections_final,
                        rankings=sel["rankings"], selections_r1=selections,
                        f1_all=ana_species_profiles,
                        f1_sel=kept_ref_final_by_label or kept_ref_by_label,
                        group_all=phylo_profiles, group_all_source=phylo_source,
                        fingerprints=fingerprints, all_units_index=all_units_index,
                        map_index_final=map_index_final, phylo_index=phylo_index,
                        df_qc_by_label=df_qc_by_label,
                        divergence=dict(fa.DIVERGENCE_TIME_MYA),
                        surface_pics=st["surface_pics"], surface_ranges=st["surface_ranges"],
                        run_index=st["run_index"], anat_common=st["anat_common"])
                except Exception as e:
                    import traceback
                    traceback.print_exc()
                    print(f"  [WARN] paper panels failed: {e}")

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
        # PHYLOGENY_USE_FINAL_ROUND says.
        use_final = (final_round >= 2 and not map_index_final.empty
                     and cfg["DISTANCE_REF_GROUP"] == "kept")
        ref_group = "kept_final" if use_final else cfg["DISTANCE_REF_GROUP"]
        ref_selections = kept_runs_final if use_final else kept_runs_r1
        ref_name = (f"round-{final_round} kept" if use_final
                    else f"round-1 {cfg['DISTANCE_REF_GROUP']}")
        print(f"\n{'=' * 70}\n  STEP 8: per-run distance to each BIDS' "
              f"{ref_name} group-map reference\n{'=' * 70}")
        # profiles of these maps were read in part 1 (07_state)
        ref_profiles = dict(group_map_profiles.get(ref_group, {}))
        if not ref_profiles:
            print(f"  [STEP 8] [WARN] no '{ref_group}' group-map profile in the saved state "
                  f"(groups: {sorted(group_map_profiles)})")
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
            if cfg.get("RUN_PAPER_FIGURES", True):
                _paper(fsp.f9_run_distance_sd, dist_df, _paper_out(paths, "F9_run_distance_sd"))

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
                        if cfg.get("RUN_PAPER_FIGURES", True):
                            _paper(fsp.f10_qc_vs_distance, meta_df,
                                   _paper_out(paths, "F10_qc_vs_distance"),
                                   alpha=cfg["QC_META_ALPHA"])
            else:
                print("  [SKIP] STEP 9: no df_qc_by_label given to run() -- "
                      "pass QC data to correlate it against the per-run distances")

    if cfg["RUN_DISTANCE_TO_REF"]:
        # ---- F8b: distance to the FINAL round's selection criterion ----
        # The profile the final round scored each BIDS's runs against (round 1:
        # LOO/non-LOO consensus; round 2: the BIDS's own round-1 kept mean;
        # round >= 3: the consensus rebuilt from the previous round's kept
        # means). With two region spaces it is the same criterion rebuilt in the
        # analysis regions (_criterion_in_space, part 1).
        crit_note = None
        if split:
            crit_note = (f"The criterion is rebuilt in the analysed regions ({len(common)} "
                         f"values; seed region {seed_txt} removed), with the same construction "
                         f"and the same kept runs at every round as the selection, which "
                         f"scored {len(ss['common'])} values (seed region included); the r "
                         f"is therefore close to, not identical with, the selection score.")
        ref_profiles_criterion = {
            label: (species_of_label.get(label), ref)
            for label, ref in final_criterion.items() if ref is not None}
        dist_df_b = fa.compute_run_distance_to_reference(
            runs_by_label, ref_profiles_criterion, common,
            value_col=cfg["VALUE_COL"], kept_keys_by_label=kept_runs_final)
        if not dist_df_b.empty:
            dist_df_b.to_excel(os.path.join(
                paths["group_maps"], "run_distance_to_selection_criterion.xlsx"),
                index=False)
            # F8b with the arousal state of every run (circles awake,
            # triangles anaesthetised; Mouse and Marmoset_Hori split and
            # tested within animal); the same figure is paper panel 6c.
            try:
                fpp.fig6c_f8b(paths, dist_df_b, bids_configs, final_round,
                              n_perm=cfg["N_PERM_PHYLOGENY"],
                              paper=cfg.get("RUN_PAPER_FIGURES", True),
                              criterion_note=crit_note)
            except Exception as e:
                import traceback
                traceback.print_exc()
                print(f"  [F8b] [WARN] state version failed ({e}) -- plain F8b drawn")
                fa.plot_run_distance_distribution(
                    dist_df_b, output_path=os.path.join(
                        paths["figures"], "F8b_run_distance_to_selection_criterion.png"),
                    title=f"Per-run correlation to the round-{final_round} SELECTION "
                          f"CRITERION (what the ranking used, not the AFNI map)")

    if cfg["MAKE_FIGURES"] and cfg["RUN_PAPER_FIGURES"]:
        fpp.organise_figures(paths["figures"])
    _write_methods(cfg, paths, sorted(profiles_by_label), n_regions, unit_counts,
                   final_round=final_round, phylo_source=phylo_source,
                   n_species=len({species_of_label.get(l) for l in phylo_profiles}))
    try:
        _write_methods_article(st, cfg, paths, phylo_source=phylo_source)
    except Exception as e:
        print(f"  [WARN] Methods_article.txt failed: {e}")
    print(f"\n{'#' * 70}\n  DONE -- {paths['root']}\n{'#' * 70}")
    return dict(config=cfg, paths=paths, common_regions=common,
                multilevel_f6=multilevel_results,
                run_distance_to_ref=dist_df, qc_vs_distance=corr_df,
                profiles_by_label=profiles_by_label,
                members_by_label=members_by_label,
                runs_by_label=runs_by_label,
                unit_counts=unit_counts,
                selection=sel, invariance=invariance, region_mask=region_mask,
                all_units_index=all_units_index, map_index=map_index,
                all_round_selections=all_round_selections, final_round=final_round,
                selections_final=selections_final, phylo_profiles=phylo_profiles,
                phylo_source=phylo_source, fingerprints=fingerprints,
                run_index=st["run_index"])
