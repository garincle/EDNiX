"""
Launcher_Fig5.py -- lance l'analyse des Figures 5 et 6 (EDNiX, seed rétrosplénial).

Tout ce qui se règle est dans ce fichier ; les valeurs par défaut de chaque clé et
leur explication détaillée sont dans fig5_pipeline.CONFIG. Le pipeline a deux parties :

  PARTIE 1 -- sélection (longue : profils de tous les runs, cartes AFNI 3dttest++,
              rendus de surface). Finit par écrire <OUTPUT_ROOT>/07_state/selection_state.pkl
              (+ selection_state.txt, lisible : datasets, régions, runs gardés, réglages).
  PARTIE 2 -- statistiques et figures, relues depuis ce fichier seul : pas d'AFNI, pas de
              NIfTI. Se relance autant de fois qu'on veut (quelques minutes), avec d'autres
              réglages de figures (clés de fig5_pipeline.FIGURE_KEYS : FIG6D, WITHIN_ACROSS_TEST,
              SURFACE_ZOOM, REFERENCE_SPECIES, ...). Si une clé de SÉLECTION change entre les
              deux parties (N_SELECT, PROFILE_METHOD, SEED, ...), la partie 2 refuse de tourner
              et liste les clés en cause : il faut relancer la partie 1.

Sorties (sous OUTPUT_ROOT) : 01_all_units, 02_selection, 03_group_maps, 03b_group_maps_final,
04_jpg, 05_surfaces, 05b_surfaces_final, 06_figures (Fig5, Fig6, supplementary,
Figure_legends.txt), 07_state, methods.txt, Methods_article.txt, 00_log.txt
(00_log_figures.txt pour une partie 2 seule). À chaque lancement, les figures et tables du
lancement précédent sont déplacées dans _previous_runs/<date>/ (jamais effacées).
"""
import os

import pandas as pd

from Exemples.Study_EDNiX.analysis.Fig5 import SBA_Region_Profile_Index as sri
from Exemples.Study_EDNiX.analysis.Fig5.fig5_pipeline import run
import Exemples.Study_EDNiX.analysis.Fig5.fig5_matrices as fm
import Exemples.Study_EDNiX.analysis.Fig5.fig5_paper as fpp

# ═══════════════════════════════════════════════════════════════════════════
# 1. CE QUI EST LANCÉ
# ═══════════════════════════════════════════════════════════════════════════
RUN_PART = "both"        # "both"      : partie 1 (sélection) puis partie 2 (figures)
                         # "selection" : partie 1 seule
                         # "figures"   : partie 2 seule, depuis 07_state/ d'une partie 1 déjà
                         #               faite avec les mêmes réglages de sélection
RUN_SEED = True          # analyse principale : cartes SBA rétrospléniales -> Fig 5, Fig 6
RUN_CONTROL = True       # contrôle de 5b : même sélection, mais une cible commune à tous les
                         # datasets (LEAVE_ONE_OUT=False) ; Python seul, pas d'AFNI
RUN_MATRIX = False        # même analyse sur les matrices EDNIxCSCLR (supplémentaire)
RUN_COMPARISON = False    # les runs sélectionnés par la seed et par les matrices sont-ils les
                         # mêmes ? (F19, supplémentaire)
# Figure 5e (dernière étape, automatique si RUN_PART != "selection") : la sélection fait-elle
# ressortir l'espèce, au-delà de sa construction ? Lit 07_state du run principal et du run de
# contrôle (RUN_CONTROL) ; sans run de contrôle, 5e est dessinée pour le run principal seul.

# ═══════════════════════════════════════════════════════════════════════════
# 2. CONSTRUCTION DU PROFIL DE CHAQUE CARTE SEED (change la sélection)
# ═══════════════════════════════════════════════════════════════════════════
# "percentiles"     : (défaut) 2 valeurs par région, en r : percentile haut et percentile bas
#                     des voxels de la région (TAIL_Q = 0.2 -> P80 et P20). Ancien nom :
#                     "quantile_posneg" (accepté).
# "median_r"        : 1 valeur par région = médiane des voxels (analyse d'origine)
# "fraction_posneg" : 2 valeurs par région, en % des voxels de la région, dans les TAIL_Q
#                     supérieurs (|pos) et inférieurs (|neg) de la carte
# Le reste de l'analyse est identique. Les matrices EDNIxCSCLR gardent median_r.
# Chaque run seed calcule aussi F22 (fiabilité des 3 profils sur tous les runs, avant
# sélection) : c'est le critère pour choisir PROFILE_METHOD.
# 6d est dessinée pour CHAQUE indice régional (06_figures/Fig6/6d_<indice>) : P80, P20,
# r médian, % dans la queue haute / basse de la carte, % de voxels r > 0 / r < 0.
PROFILE_METHOD = "percentiles"
TAIL_Q = 0.2             # taille des queues : 0.2 ou "quintile" (P80 / P20), 0.1 ou "decile"
                         # (P90 / P10), 0.25 ou "quartile" (P75 / P25)

# suffixe du dossier de sortie : même valeur pour 0.2 et "quintile"
_Q_PCT = int(round(100 * float({"decile": 0.1, "deciles": 0.1, "quintile": 0.2, "quintiles": 0.2,
                                "quartile": 0.25, "quartiles": 0.25}.get(TAIL_Q, TAIL_Q))))

LIB = "/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals"
SEED = ["BA_26_29_30__retrosplenial"]

# ═══════════════════════════════════════════════════════════════════════════
# 3. DATASETS
# ═══════════════════════════════════════════════════════════════════════════
# label     : nom du dataset partout (dossiers, figures, json d'équilibrage AFNI) -- le
#             changer fait recalculer ses cartes AFNI kept / not kept
# species   : espèce (doit exister dans SPECIES_TO_TREE_TIP de fig5_pipeline.CONFIG)
# bids_dir  : dossier BIDS EDNiX ; label_path : fichier de labels ITK-SNAP de l'atlas
# state     : état d'éveil. Chaîne = tout le dataset ; {"from": "task"} = d'après task-<x>
#             du nom du BOLD ; {"from": "subject"} = d'après le nom du sujet (premier motif
#             trouvé : MEDISO avant ISO) ; None = inconnu. Tout état autre que "awake" =
#             anesthésié. "animal" (regex) : retirée du nom du sujet pour retrouver l'animal
#             quand l'état est codé dans le nom (Mouse).
# endfmri   : motif des BOLD du dataset. Non utilisé ici : Macaque Tremblay (propofol).
STATE_HORI = {"from": "task", "map": {"restawake": "awake", "restanest": "anesthetized"}}
STATE_MOUSE = {"from": "subject", "map": [["MEDISO", "medetomidine+isoflurane"],
                                          ["ISO", "isoflurane"],
                                          ["AW", "awake"]],
               # l'animal est le suffixe : jgrAesAWc11L, jgrAesISOc11L, jgrAesMEDISOc11L
               # = même souris c11L dans 3 états -> contraste intra-animal
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

# ═══════════════════════════════════════════════════════════════════════════
# 4. RÉGLAGES (clés de fig5_pipeline.CONFIG ; une clé absente garde sa valeur par défaut)
# ═══════════════════════════════════════════════════════════════════════════
BASE_CONFIG = dict(
    # --- sorties ---------------------------------------------------------------------
    # Un dossier par méthode de profil : les résultats d'une méthode n'écrasent pas ceux
    # d'une autre. Même dossier possible : les cartes AFNI déjà calculées (même suffixe,
    # mêmes runs en entrée, vérifiés) sont réutilisées.
    # Dossier de V12 : sa sélection est celle que redonne SEED_REGION_EXCLUDE_SELECTION=False,
    # donc ses cartes AFNI (all-units, kept / not kept des tours 1 et final) sont réutilisées
    # (vérifiées dans input_runs/) ; seules celles de Human_Babayan (ex-BIDS_Babayan) sont
    # recalculées. Les figures et tables de V12 vont dans _previous_runs/<date>/.
    OUTPUT_ROOT="/scratch2/EDNiX/results/EDNiX_Figures/Fig_5_V15_"
                + {"median_r": "", "fraction_posneg": f"_posneg{_Q_PCT}", "percentiles": f"_quant{_Q_PCT}",
                   "quantile_posneg": f"_quant{_Q_PCT}"}[PROFILE_METHOD],
    PROFILE_METHOD=PROFILE_METHOD,   # voir section 2
    TAIL_Q=TAIL_Q,                   # voir section 2

    # --- données et unité d'analyse -------------------------------------------------
    UNIT="run",              # unité de la sélection et des cartes AFNI : "run" | "session" |
                             # "subject" (un animal = un échantillon partout)
    ATLAS_LEVEL=3,           # niveau d'atlas EDNiX des régions (2 ou 3)
    RUN_QC_ENABLE=False,     # QC des runs avant tout (distance MAD au médoïde du dataset)

    # --- sélection ---------------------------------------------------------------------
    N_SELECT=10,             # unités gardées par dataset à chaque tour
    INVARIANCE_MODE="top_fraction",  # régions utilisées pour le score de sélection :
    INVARIANCE_VALUE=0.8,    # les 80 % les plus invariantes entre espèces (1.0 = toutes)
    SIMILARITY="correlation",  # score = r de Pearson entre le profil et la cible
    ZSCORE="scale",          # normalisation des profils : "scale" (RMS, signe gardé), True
                             # (z-score) ou False (brut)
    LEAVE_ONE_OUT=True,      # la cible d'un dataset exclut toute son espèce (pas de circularité)
    RUN_ROUND2_SELECTION=True,     # tour 2 : chaque dataset reclassé contre ses propres gardés
    RUN_ITERATIVE_SELECTION=True,  # tours 3..N : cible inter-espèces reconstruite des gardés
    N_ITERATIVE_ROUNDS=10,         # numéro du DERNIER tour (sélection finale)
    SEED_REGION_EXCLUDE="auto",    # région de l'atlas qui contient le seed (sa valeur reflète
                                   # l'autocorrélation du seed) : "auto" = trouvée dans
                                   # l'atlas ; liste de noms ; None = aucune
    SEED_REGION_EXCLUDE_SELECTION=False,  # False : gardée dans la sélection (= sélection V12)
    SEED_REGION_EXCLUDE_ANALYSES=True,    # True : retirée de toutes les analyses et figures

    # --- phylogénie : F6 seulement (F6 d'origine + F6ab) ------------------------------
    PHYLOGENY_USE_FINAL_ROUND=True,   # (figures) cartes sélectionnées ; la version sans
        # sélection est dans 5d / F17 et fingerprint_tables.xlsx
    TREE_PATH="/scratch2/EDNiX/Liste des espèces.nwk",  # (figures) arbre daté (Newick)
    RUN_MULTILEVEL_F6=False,          # (figures) F6 aux autres niveaux d'atlas
    RUN_PHYLO_COMPARATIVE=True,      # (figures) Mantel / K de Blomberg / PGLS (F11-F13)

    # --- figures d'interprétation (F1/F1c fingerprints, F6ab, F15-F18) ---------------
    RUN_FINGERPRINTS=True,   # (figures)
    REFERENCE_SPECIES="Human",  # (figures) référence de F6b / F15b / F18

    # --- panneaux du papier : 06_figures/Fig5, Fig6, supplementary --------------------
    RUN_PAPER_FIGURES=True,  # (figures) taille d'impression (89 / 120 / 183 mm), 5-7 pt,
        # PDF à polices TrueType, PNG 600 dpi ; légendes dans chaque .txt et dans
        # 06_figures/Figure_legends.txt
    BRAIN_MASS_XLSX="/scratch2/EDNiX/Herculano-Houzel_mammalian_brain_dataset_complete.xlsx",
    WITHIN_ACROSS_PAIRS="other_dataset",  # (figures) 5b / F2 / F2b, paires comparées :
        # "other_dataset" = même espèce dans un AUTRE dataset vs espèces différentes (les
        # paires d'un même dataset partagent protocole, scanner et souvent l'animal : exclues
        # du test) ; "all" = toutes les paires de même espèce (5b de V14)
    WITHIN_ACROSS_TEST="wilcoxon",  # (figures) 5b, test affiché : "wilcoxon" (rank-sum sur
        # les paires, unilatéral, non corrigé ; optimiste car les paires partagent des runs)
        # ou "permutation" (espèces permutées entre datasets) ; les deux sont dans le .txt
    FIG6D=dict(              # (figures) 6d, un panneau par indice : une valeur par région et par run
        # sélectionné ; Human + Macaque contre les 5 autres espèces, un test exact
        # bilatéral par région (espèce = unité, 21 affectations, p >= 0,048), non corrigé
        regions=("posterior_parietal_cortex", "prefrontal_cortex_(pfc)",
                 "somatosensory_cortex", "visual_striate_cortex", "striatum",
                 "motor_and_premotor"),
        group=("Human", "Macaque")),
    SURFACE_ZOOM={"Human": 1.0, "Macaque": 1.0, "Dog": 1.0, "Marmoset": 0.9,
                  "Mouselemur": 0.9, "Rat": 0.9, "Mouse": 0.9},  # (figures) 5a / S5a :
        # taille relative des cerveaux par espèce (1 = remplit la case)
    SURFACE_DATASETS_PER_ROW=2,      # (figures) 5a / S5a : datasets par ligne

    # --- analyses complémentaires ----------------------------------------------------
    RUN_DISTANCE_TO_REF=True,  # (figures) F8 / F8b / 6c : r de chaque run avec la carte
        # gardée et avec le critère de sélection, par état d'éveil

    # --- cartes AFNI (partie 1) --------------------------------------------------------
    OVERWRITE_MAPS=True,    # False : une carte déjà calculée avec exactement les mêmes runs
                             # est réutilisée (vérifié dans input_runs/) ; True : recalcul
    MAP_TAG=None,            # None = suffixe unique dérivé de OUTPUT_ROOT
    RENDER_SURFACE=True,     # rendus wb_command des cartes de groupe (5a, S5a)
)

# Même analyse sur les matrices : pas de carte AFNI, rendu matriciel des fingerprints.
MATRIX_ATLAS = "EDNIxCSCLR"
MATRIX_TYPE = "correlation"
MATRIX_CONFIG = dict(
    BASE_CONFIG,
    OUTPUT_ROOT=BASE_CONFIG["OUTPUT_ROOT"] + f"_matrix_{MATRIX_ATLAS}",
    SEED=f"{MATRIX_ATLAS}_{MATRIX_TYPE}_matrix",
    PROFILE_KIND="matrix",
    RUN_ALL_UNITS_MAPS=False, RUN_GROUP_MAPS=False, RENDER_SURFACE=False,
    COPY_JPG=False, RUN_DISTANCE_TO_REF=False,
)
COMPARISON_DIR = BASE_CONFIG["OUTPUT_ROOT"] + "_seed_vs_matrix"

# Contrôle de 5b : chaque dataset est sélectionné vers la MÊME cible (consensus de
# toutes les espèces). Deux datasets d'une même espèce n'ont alors plus de cible
# commune qui leur soit propre : si l'écart « même espèce > espèces différentes »
# persiste, il n'est pas mécanique. Python seul (pas d'AFNI) ; régions communes =
# celles du run principal.
CONTROL_CONFIG = dict(
    BASE_CONFIG,
    OUTPUT_ROOT=BASE_CONFIG["OUTPUT_ROOT"] + "_control_common_target",
    LEAVE_ONE_OUT=False,
    PHYLOGENY_USE_FINAL_ROUND=False,
    RUN_ALL_UNITS_MAPS=False, RUN_GROUP_MAPS=False, RENDER_SURFACE=False,
    COPY_JPG=False, RUN_DISTANCE_TO_REF=False,
)

# ═══════════════════════════════════════════════════════════════════════════
# 5. QC EXTERNE (tSNR, FD, ... par run ; PGLS et validation externe de la sélection)
# ═══════════════════════════════════════════════════════════════════════════
QC_EXCEL_PATH = "/scratch2/EDNiX/results/EDNiX_Figures/Suppl_Figure_4-5-6_QC/qc_raw_per_individual.xlsx"
df_qc_by_label = {}
if RUN_PART != "figures":           # la partie 2 seule relit la table QC de 07_state
    df_qc_all = pd.read_excel(QC_EXCEL_PATH, dtype={"subject": str, "session": str, "run": str})
    # Correspondance EXACTE sur le dernier élément du chemin : str.contains faisait
    # matcher "BIDS_Grandjean" (Rat) dans "BIDS_Grandjean2" (Mouse).
    qc_bids_name = (df_qc_all["bids_dir"].astype(str).str.strip().str.rstrip("/")
                    .str.split("/").str[-1])
    for cb in BIDS_CONFIGS:
        bids_name = cb["bids_dir"].rstrip("/").split("/")[-1]
        mask = qc_bids_name == bids_name
        if mask.any():
            df_qc_by_label[cb["label"]] = df_qc_all[mask].copy()
        else:
            print(f"  [WARN] {cb['label']}: aucune ligne QC trouvée pour '{bids_name}'")
    print(f"  df_qc_by_label: {len(df_qc_by_label)}/{len(BIDS_CONFIGS)} dataset(s) avec QC")

# ═══════════════════════════════════════════════════════════════════════════
# 6. LANCEMENT
# ═══════════════════════════════════════════════════════════════════════════
# ── 1. analyse seed ─────────────────────────────────────────────────────────
if RUN_SEED:
    res = run(sri, BIDS_CONFIGS, config=BASE_CONFIG, df_qc_by_label=df_qc_by_label,
              part=RUN_PART)

# ── 2. contrôle de 5b : cible commune ─────────────────────────────────────────
if RUN_CONTROL:
    regions = fm.read_node_regions(BASE_CONFIG["OUTPUT_ROOT"])
    if regions is None:
        print("  [CONTROL] pas de common_regions.xlsx du run principal -- contrôle sauté "
              "(lancer RUN_SEED d'abord)")
    else:
        res_control = run(sri, BIDS_CONFIGS,
                          config=dict(CONTROL_CONFIG, COMMON_REGIONS_OVERRIDE=regions),
                          df_qc_by_label=df_qc_by_label, part=RUN_PART)

# ── 3. même analyse sur les matrices ──────────────────────────────────────────
if RUN_MATRIX:
    # noeuds = régions communes de l'analyse seed (même anatomie des deux côtés)
    nodes = fm.read_node_regions(BASE_CONFIG["OUTPUT_ROOT"])
    if nodes is None:
        print("  [MATRIX] pas de common_regions.xlsx côté seed -- tous les noeuds du "
              "niveau d'atlas sont utilisés")
    msri = fm.MatrixProfiles(sri, BIDS_CONFIGS, atlas_name=MATRIX_ATLAS,
                             correl_type=MATRIX_TYPE, node_regions=nodes,
                             value_col=BASE_CONFIG.get("VALUE_COL", "median_r"))
    res_matrix = run(msri, BIDS_CONFIGS, config=MATRIX_CONFIG, df_qc_by_label=df_qc_by_label,
                     part=RUN_PART)
    if RUN_PART != "figures":
        mapping = msri.mapping_table()
        mapping.to_excel(os.path.join(MATRIX_CONFIG["OUTPUT_ROOT"], "02_selection",
                                      "matrix_run_mapping.xlsx"), index=False)
        print(f"  [MATRIX] matrix_run_mapping.xlsx : {len(mapping)} matrice(s), méthode "
              f"{mapping['method'].value_counts().to_dict() if len(mapping) else {}}")
    # F21 : matrices EDNIxCSCLR moyennes des runs gardés par la sélection SEED,
    # non gardés et tous (dans le dossier de comparaison seed vs matrices)
    if RUN_PART != "selection":
        fm.plot_mean_matrices(BASE_CONFIG["OUTPUT_ROOT"], res_matrix["runs_by_label"],
                              {c["label"]: c["species"] for c in BIDS_CONFIGS}, COMPARISON_DIR,
                              value_col=BASE_CONFIG.get("VALUE_COL", "median_r"))

# ── 4. les runs sélectionnés sont-ils les mêmes ? ─────────────────────────────
if RUN_COMPARISON and RUN_PART != "selection":
    per_dataset, overall = fm.compare_selections(BASE_CONFIG["OUTPUT_ROOT"],
                                                 MATRIX_CONFIG["OUTPUT_ROOT"], COMPARISON_DIR)
    print(overall.to_string(index=False))

# ── 5. Figure 5e : la sélection fait-elle ressortir l'espèce, au-delà de sa construction ? ──
# Unité = la paire de datasets ; gain = écart(sélectionnés) - écart(tous les runs), testé sur
# toutes les répartitions des espèces entre datasets ; run principal et run de contrôle.
# Écrit 06_figures/Fig5/5e_selection_species_gain (.png/.pdf/.txt/.xlsx) et complète le .txt de 5b.
if RUN_PART != "selection" and (RUN_SEED or RUN_CONTROL):
    fpp.fig5e_selection_species_gain(BASE_CONFIG["OUTPUT_ROOT"], CONTROL_CONFIG["OUTPUT_ROOT"])
