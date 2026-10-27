"""
Run selection, part 1 of the former Fig5 pipeline: run profiles -> unit (subject)
profiles -> AFNI group maps of all units -> common regions -> cross-species
consensus selection (iterative, leave-one-species-out) -> group maps of kept and
not-kept units -> 07_state/selection_state.pkl.

Moved VERBATIM from Exemples/Study_EDNiX/analysis/Fig5/fig5_pipeline.py (2026-10-27):
CONFIG, the steps and _run_selection are the same code with the same defaults; only
the module imports point to this package. Entry point: run_selection() below, or
launcher.py. The paper statistics (part 2) read the state file this writes.
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

from Statistics.Group_fMRI.run_selection import consensus as fa     # selection, invariance, balance
from Statistics.Group_fMRI.run_selection import units as fu         # run -> unit aggregation, run QC
from Statistics.Group_fMRI.run_selection import reliability as ff   # F22 profile reliability, arousal state
from Statistics.Group_fMRI.run_selection import surfaces as fpp     # surface renders on one scale


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
    # subject list: classification of the former network-QC figure (archived); None or a
    # missing file with CATEGORY="All" = every subject with a func/ folder in the BIDS
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
            from Statistics.Group_fMRI.run_selection.external_qc import merge_qc_with_classification
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


# ═══════════════════════════════════════════════════════════════════════════
# Entry point (part 1 only)
# ═══════════════════════════════════════════════════════════════════════════
def _normalise_tree_path(cfg):
    """TREE_PATH as the former pipeline stored it in the state: "è" can be stored as
    one code point (NFC) or two (NFD); the file name found on disk is kept. The tree
    itself is not read here (divergence times are only used by the statistics)."""
    if not cfg.get("TREE_PATH") or os.path.isfile(cfg["TREE_PATH"]):
        return
    import unicodedata
    d, name = os.path.split(cfg["TREE_PATH"])
    if os.path.isdir(d or "."):
        want = unicodedata.normalize("NFC", name)
        hit = [f for f in os.listdir(d or ".") if unicodedata.normalize("NFC", f) == want]
        if hit:
            cfg["TREE_PATH"] = os.path.join(d, hit[0])


def run_selection(sri, bids_configs, config=None, df_qc_by_label=None):
    """Part 1 alone: same steps, same order and same outputs as the former
    fig5_pipeline.run(..., part="selection"). Ends with
    <OUTPUT_ROOT>/07_state/selection_state.pkl (+ .txt summary) and methods.txt.

    sri          : Statistics.Group_fMRI.run_selection.profiles
    bids_configs : list of per-dataset dicts (label, species, bids_dir, label_path,
                   seeds, state; see launcher.py)
    config       : dict overriding CONFIG entries
    df_qc_by_label : {label: QC table} for the external QC validation (optional)
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
    archived = _archive_previous_outputs(paths)
    tee, old_cwd = _open_log(paths, "00_log.txt")
    if archived:
        print(f"  Previous run's figures/tables/log moved to {archived} "
              f"(AFNI maps and *.json balance records left in place for reuse)")
    print("  part = 'selection'")
    print(f"  AFNI statmaps stored under RESULTS_BASE_DIR={cfg['RESULTS_BASE_DIR']!r}, "
          f"suffix tag MAP_TAG={cfg['MAP_TAG']!r}, "
          f"{'THRESHOLDED' if cfg['STAT_MAP_USE_THRESHOLDED'] else 'unthresholded'} "
          f"group maps")
    for w in config_warnings:
        print(f"  [CONFIG WARNING] {w}")
    _normalise_tree_path(cfg)
    if hasattr(sri, "set_profile_method"):
        sri.set_profile_method(cfg["PROFILE_METHOD"], cfg["TAIL_Q"],
                               keep_sketch=cfg.get("PROFILE_KIND", "seed") == "seed")
    print(f"  Profile method: {cfg['PROFILE_METHOD']}"
          + (f" (q = {cfg['TAIL_Q']:g}, value column {cfg['VALUE_COL']})"
             if cfg["PROFILE_METHOD"] != "median_r" else ""))
    try:
        state = _run_selection(sri, cfg, paths, bids_configs, df_qc_by_label)
        save_state(state)
        _write_methods(cfg, paths, sorted(state["profiles_by_label"]),
                       state["n_regions"], state["unit_counts"],
                       final_round=state["final_round"])
        print(f"\n{'#' * 70}\n  PART 1 DONE -- {paths['root']}\n{'#' * 70}")
        return state
    finally:
        if hasattr(sri, "set_profile_method"):
            sri.set_profile_method("median_r")
        _close_log(tee, old_cwd, paths, "00_log.txt")
