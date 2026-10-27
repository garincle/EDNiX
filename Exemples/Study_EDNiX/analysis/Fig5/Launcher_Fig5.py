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

from Statistics.Group_fMRI.run_selection import profiles as sri
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
# 2-3. PROFIL, DATASETS ET RÉGLAGES DE SÉLECTION
# ═══════════════════════════════════════════════════════════════════════════
# Une seule source : Statistics/Group_fMRI/run_selection/config_ednix.py (méthode de profil,
# TAIL_Q, datasets, réglages de la sélection). La partie 1 se lance aussi seule depuis
# Statistics/Group_fMRI/run_selection/launcher.py.
from Statistics.Group_fMRI.run_selection.config_ednix import (   # noqa: E402
    PROFILE_METHOD, TAIL_Q, _Q_PCT, LIB, SEED, STATE_HORI, STATE_MOUSE, BIDS_CONFIGS,
    SELECTION_CONFIG, QC_EXCEL_PATH)

# ═══════════════════════════════════════════════════════════════════════════
# 4. RÉGLAGES (clés de fig5_pipeline.CONFIG ; une clé absente garde sa valeur par défaut)
# ═══════════════════════════════════════════════════════════════════════════
BASE_CONFIG = dict(
    SELECTION_CONFIG,   # sorties, profil, unité, sélection, cartes AFNI : config_ednix.py
    # Réglages des figures (partie 2) ci-dessous.

    # --- phylogénie : F6 seulement (F6 d'origine + F6ab) ------------------------------
    PHYLOGENY_USE_FINAL_ROUND=True,   # (figures) cartes sélectionnées ; la version sans
        # sélection est dans 5d / F17 et fingerprint_tables.xlsx
    # TREE_PATH (arbre daté, Newick) : config_ednix.SELECTION_CONFIG
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
# QC_EXCEL_PATH : config_ednix.py
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
