# Figures 5 et 6 — dossier `Fig5/`

Point d'entrée : `Launcher_Fig5.py`. Chaque option y est décrite ; les valeurs par défaut et leur explication détaillée sont dans `fig5_pipeline.CONFIG`.

## Deux parties : sélection, puis figures (`RUN_PART`)

| `RUN_PART` | ce qui tourne |
|---|---|
| `"both"` (défaut) | partie 1 puis partie 2 |
| `"selection"` | **partie 1** : profils de tous les runs, QC, cartes AFNI all-units, régions communes et région du seed, indices régionaux de chaque run (6d), table F22, tours de sélection, cartes AFNI kept / not kept, jpg, rendus de surface (5a), profils de toutes les cartes de groupe. Finit par `07_state/selection_state.pkl` (+ `selection_state.txt` lisible : datasets, régions, runs gardés, réglages de sélection) |
| `"figures"` | **partie 2** seule : toutes les statistiques et figures, relues depuis `07_state/` — pas d'AFNI, pas de NIfTI (sauf `RUN_MULTILEVEL_F6`, désactivé). Relançable autant de fois qu'on veut |

- Seules les clés de `fig5_pipeline.FIGURE_KEYS` peuvent changer entre les deux parties (`FIG6D`, `WITHIN_ACROSS_TEST`, `SURFACE_ZOOM`, `SURFACE_DATASETS_PER_ROW`, `REFERENCE_SPECIES`, `PHYLOGENY_USE_FINAL_ROUND`, `TREE_PATH`, `BRAIN_MASS_XLSX`, `RUN_DISTANCE_TO_REF`, …). Si une clé de sélection a une autre valeur (`N_SELECT`, `PROFILE_METHOD`, `TAIL_Q`, `SEED`, …), la partie 2 refuse de tourner et liste les clés en cause.
- Une partie 1 déplace l'ancien `07_state/` dans `_previous_runs/` dès le départ : une partie 1 interrompue ne laisse jamais une ancienne sélection à la partie 2.
- Une partie 2 seule ne déplace que `06_figures/` (et son log `00_log_figures.txt`).
- Les runs de contrôle et de matrices suivent `RUN_PART` (chacun a son `07_state/`).

Étapes du launcher : `RUN_SEED` (Fig 5, Fig 6), `RUN_CONTROL` (contrôle de 5b, cible commune, Python seul), `RUN_MATRIX` (matrices EDNIxCSCLR, supplémentaire, + F21), `RUN_COMPARISON` (F19).

## Construction du profil : `PROFILE_METHOD`, `TAIL_Q`

| valeur | profil d'une carte seed |
|---|---|
| `"percentiles"` (**défaut** ; ancien nom `"quantile_posneg"` accepté) | 2 valeurs par région, **en r** : percentile haut et percentile bas de ses voxels. `TAIL_Q = 0.2` ou `"quintile"` → P80 / P20 ; `0.1` ou `"decile"` → P90 / P10 ; `0.25` ou `"quartile"` → P75 / P25. Colonne `q20q80_r` |
| `"median_r"` | 1 valeur par région : la médiane de ses voxels |
| `"fraction_posneg"` | 2 valeurs par région, **en % des voxels de la région**, dans les `TAIL_Q` supérieurs (`|pos`) et inférieurs (`|neg`) de la carte (percentiles pris sur les voxels des régions communes). Colonne `pct_tail20` |

- **F22** (`supplementary/F22_profile_reliability_by_method`, `02_selection/profile_reliability.xlsx`) : dans chaque run seed, les trois profils sont construits sur tous les runs avant sélection, et comparés par leur discriminabilité (Bridgeford et al. 2021). C'est le critère fixé d'avance pour choisir `PROFILE_METHOD`.
- **Gauche / droite** poolés par nombre de voxels ; **chaque canal standardisé** sur les régions communes (sans ça, le décalage P80 > P20 dominerait le Pearson ; avec, le r entre deux profils est la moyenne des r des deux canaux). La valeur brute est gardée dans `<colonne>_raw`.
- **Exactitude** : jusqu'à 2000 voxels, tous les voxels de la région ; au-delà, ses 2001 percentiles (erreur < 0,05 point).
- **Dossier** : `OUTPUT_ROOT` reçoit le suffixe `_quant20` (percentiles) ou `_posneg20` (fraction_posneg), pour ne pas écraser les résultats d'une autre méthode.
- **Figures** : chaque rayon des fingerprints (5c, F1, F18...) porte son canal, comme en V12 : `PFC P80` / `PFC P20` (ou `PFC +` / `PFC −`), tous les canaux hauts puis tous les canaux bas. En V14, le canal n'était écrit qu'une fois par moitié de radar : le profil avait bien ses deux canaux, seul l'affichage avait changé.
- **Matrices** EDNIxCSCLR : toujours `median_r`.

## Région du seed (`SEED_REGION_EXCLUDE`, `SEED_REGION_EXCLUDE_SELECTION`, `SEED_REGION_EXCLUDE_ANALYSES`)

La région du niveau `ATLAS_LEVEL` qui contient le seed a une valeur dominée par l'autocorrélation du seed, pas par un couplage. `SEED_REGION_EXCLUDE` la trouve : `"auto"` (défaut) par les voxels du label du seed dans l'atlas de chaque dataset (periarchicortex, label 29, qui contient aussi BA 28, 35/36, présubiculum, BA 33, BA 25) ; une liste de noms l'impose ; `None` = aucune. Détail dans `02_selection/seed_region_exclusion.xlsx`. Deux interrupteurs disent où elle est retirée :

| clé | défaut | effet |
|---|---|---|
| `SEED_REGION_EXCLUDE_SELECTION` | `False` | `True` : retirée des profils que la sélection note. `False` : gardée, ce qui redonne exactement la sélection de V12 |
| `SEED_REGION_EXCLUDE_ANALYSES` | `True` | `True` : retirée de tous les profils utilisés après la sélection (partie 2) |

- **Pourquoi la sélection change avec cette région** : en V12 elle était déjà hors du score (le masque d'invariance écartait `periarchicortex|q80` et `|q20`), mais chaque canal est standardisé sur les régions communes : son pic (P80 très haut) entrait dans la moyenne et l'écart-type du canal P80, donc dans les valeurs standardisées de toutes les autres régions, et dans le classement d'invariance.
- **Deux espaces de régions** quand les deux clés diffèrent. Les deux viennent des mêmes voxels des mêmes runs ; seuls l'ensemble de régions et, avec deux canaux, la standardisation sur cet ensemble changent.
  - L'espace de **sélection** sert aux profils notés, à l'invariance, au consensus, aux scores de tous les tours, à `selection*.xlsx`, `region_invariance.xlsx` et `balanced_stats.xlsx`. Il sert aussi à **F3, F4 et F5**, qui décrivent la sélection.
  - L'espace d'**analyse** sert à tout le reste : F0–F2, 5b–6d, F6–F22, Mantel/K/PGLS, profils des cartes de groupe et 6d.
- **6c / F8b** : le critère du tour final est **reconstruit dans les régions analysées**, avec la même construction et les mêmes runs gardés à chaque tour (`_criterion_in_space`). Le .txt le signale. Quand il n'y a qu'un espace, ce critère reconstruit est identique au critère de la sélection (vérifié).
- `common_regions.xlsx` : toutes les régions, avec `seed_region`, `used_by_selection` et `used_by_analyses`. `07_state/selection_state.pkl` garde l'espace de sélection dans `selection_space`.
- Ces deux clés sont des clés de **partie 1** : en changer une impose de relancer la partie 1. Un état écrit par V14, où la région était retirée partout, est lu comme `True` / `True`.

## Panneaux du papier (qualité publication)

- **Taille d'impression** : largeurs d'une page Nature (89 mm, 120 mm, 183 mm), texte 5–7 pt, traits 0,5–1 pt, police Arial (ou Liberation Sans, son clone métrique), PDF à polices TrueType (texte éditable), PNG 600 dpi ; 5a et S5a aussi en `.svg`.
- **Sur les figures** : seulement les tests, en **étoiles** (*** P < 0,001, ** P < 0,01, * P < 0,05, n.s.) avec une barre entre les éléments comparés (5b, 6c, 6d), ou l'effet suivi des étoiles (5c : r moyen entre espèces ; 5d : ρ ; 6a : r). **Aucun n ni p exact sur les figures** : ils sont dans le `.txt`.
- **Chaque `.txt`** commence par `LEGEND` (la légende de la figure, style article), puis question, ce qui est tracé, test, résultat complet (n, effets, p exacts), lecture. Toutes les légendes sont aussi réunies dans **`06_figures/Figure_legends.txt`**.
- **`Methods_article.txt`** (à la racine) : les Méthodes des Figures 5–6 rédigées comme pour l'article, avec les chiffres du run (datasets, runs, animaux, régions, sélection, tests) et les références. `methods.txt` garde les paramètres verbatim.

| panneau | = supplémentaire | question | test (sur la figure : étoiles) |
|---|---|---|---|
| **5a** `5a_surface_maps_selected_vs_excluded` | S5a | (illustration) cartes de surface, selected au centre, excluded sur les côtés, une échelle de couleur par dataset (rendu en partie 1) | — |
| **5b** `5b_intra_vs_inter_species` | F2b | Les profils de deux runs d'une même espèce se ressemblent-ils plus que ceux de deux espèces différentes (hypothèse principale) ? `WITHIN_ACROSS_PAIRS="other_dataset"` (défaut) : paires **même espèce, autre dataset** contre paires **espèces différentes** ; les paires d'un même dataset (même protocole, même scanner, souvent même animal) sont exclues du test et restent visibles sur F2b. `"all"` : toutes les paires de même espèce (5b de V14, dominé par les paires du même dataset). Mêmes paires, classes, bins et statistiques que F2b | `WITHIN_ACROSS_TEST="wilcoxon"` : Wilcoxon rank-sum unilatéral sur les paires, non corrigé. **Les paires partagent des runs, donc ne sont pas indépendantes : ce p est optimiste** (Mantel 1967 ; Nili et al. 2014). Le test par permutation des espèces entre datasets, qui respecte cette dépendance, est toujours dans le `.txt` (`"permutation"` le met sur la figure) |
| **5e** `5e_selection_species_gain` | — | La sélection fait-elle ressortir l'espèce, au-delà de ce qu'elle fait à toutes les paires et de sa construction ? (a) matrice dataset × dataset du r médian entre runs, tous les runs / sélectionnés, paires de même espèce encadrées ; (b, c) une ligne par paire de datasets, tous → sélectionnés, sélection principale et sélection de contrôle (cible commune) ; (d) gain observé et sa distribution nulle | unité = **paire de datasets** (45, dont 3 de même espèce) ; écart = moyenne des paires de même espèce − moyenne des autres ; **gain = écart(sélectionnés) − écart(tous)** = interaction sélection × espèce ; nul exact sur toutes les répartitions des espèces entre datasets qui gardent la structure observée (3150) ; unilatéral. Dessinée en dernière étape du launcher (lit `07_state` des deux runs) ; résultat ajouté au `.txt` de 5b |
| **5c** `5c_fingerprint_species` | F1_fingerprint | Le profil rétrosplénial est-il partagé entre espèces ? | r moyen entre espèces vs 10 000 permutations des régions (libéral : ignore l'autocorrélation spatiale) |
| **5d** `5d_rank_before_after_selection` | F17 a | La sélection réordonne-t-elle les datasets ? | Spearman avant / après, 10 000 permutations |
| **6a** `6a_distance_divergence_brain_mass` | F6ab + masse cérébrale | Similarité au consensus selon la divergence, à Human, selon la masse cérébrale | Pearson sur les moyennes d'espèce, p exact (7! ou 6!), avec et sans Human |
| **6b** `6b_fingerprint_vs_consensus` | F15a | Quelles régions portent l'écart de chaque dataset ? | descriptif. Mise en page de V12 (même taille de radar, 5 par ligne), **sans nom de région** (même ordre que 5c), espaces resserrés |
| **6c** `6c_arousal_state_F8b` | F8b | Les runs anesthésiés sont-ils moins corrélés au critère de sélection que les runs éveillés des mêmes animaux ? | sign-flip exact des différences par animal (barre entre éveillé et anesthésié) ; n animaux dans le `.txt` |
| **6d** : un panneau par indice dans `Fig6/` : `6d_upper_P80`, `6d_lower_P20`, `6d_median_r`, `6d_top20_fraction`, `6d_bottom20_fraction`, `6d_positive_fraction` (% de voxels r > 0), `6d_negative_fraction` (% r < 0) | — | La valeur de six régions (pariétal postérieur, PFC, somatosensoriel, V1, striatum, moteur) diffère-t-elle chez Human + Macaque par rapport aux 5 autres espèces, au-delà de ce que prédit la taille de la région ? | un test exact bilatéral par région, **non corrigé** ; espèce = unité (21 affectations, p ≥ 0,048) ; valeurs run → animal → dataset → espèce sur les runs sélectionnés. Droite : contrôle de taille (ajustement sur les autres régions ± 2 SD) |

**5e : pourquoi la paire de datasets.** Les centaines de paires de runs d'une paire de datasets partagent les deux datasets (protocole, scanner, animaux, sélection) : elles valent une observation. Le Wilcoxon de 5b, sur les paires de runs, en compte des milliers et donne un p beaucoup trop petit. La sélection rapproche toutes les paires (cible inter-espèces commune), donc l'écart après sélection seul ne dit rien ; c'est le gain qui compte. Dans la sélection principale, deux datasets d'une même espèce sont classés contre la même cible (consensus sans leur espèce), ce qui les rapproche par construction ; le run de contrôle supprime cette cible propre : un gain qui persiste dans le contrôle n'est pas mécanique. Avec 3 paires de même espèce, la puissance est faible.

**6d pour chaque indice**, automatiquement, l'indice de `PROFILE_METHOD` en premier ; il n'y a plus de S6d. Tous viennent de la table des indices régionaux de chaque run, calculée en partie 1 sur les voxels (`sri.region_indices`), donc la partie 2 n'a pas besoin des voxels. Réglage : `FIG6D` (`regions`, `group`). Une carte de run est plus bruitée qu'une carte de groupe ; le contraste de groupe confond taille du cerveau, gyrification, état d'éveil et acquisition ; Human et Macaque sont les deux seuls catarrhiniens échantillonnés (Marmoset et Mouselemur sont dans l'autre groupe).

**Source des panneaux de distance (6a, 6b)** : les cartes de groupe de l'analyse phylogénétique (cartes sélectionnées si `PHYLOGENY_USE_FINAL_ROUND=True`). La version sans sélection est en 5d et dans `fingerprint_tables.xlsx`.

**Masse cérébrale** : classeur Herculano-Houzel (`BRAIN_MASS_XLSX`) ; `BRAIN_MASS_G={espèce: g}` le remplace.

## Autres figures (supplémentaires)

**F20** : la sélection réduit-elle la variabilité entre runs, entre sessions ou entre animaux ? **F21** (`<OUTPUT_ROOT>_seed_vs_matrix/`) : matrices EDNIxCSCLR moyennes tous / gardés / non gardés d'après la sélection seed. **Diagnostic des datasets** : `03_group_maps/dataset_diagnostics.xlsx`.

## Fichiers

| fichier | rôle |
|---|---|
| `Launcher_Fig5.py` | réglages et lancement (toutes les options décrites) |
| `fig5_pipeline.py` | `run(..., part=)`, `CONFIG`, `FIGURE_KEYS` ; `_run_selection` (partie 1), `_run_figures` (partie 2), `save_state` / `load_state`, `Methods_article.txt` |
| `SBA_Region_Profile_Index.py` | profils régionaux ; `set_profile_method`, `to_tail_profile`, `to_quantile_profile`, `region_indices` |
| `fig5_paper.py` | panneaux 5a–6d, S5a, S6d, F8b par état, légendes (`Figure_legends.txt`), rangement |
| `fig5_fingerprints.py` | F6ab, fingerprints, F15–F18, F20, F22, F8b par état ; mode `paper=True` à la taille d'impression |
| `fig5_analysis.py` | sélection, F2 / F2b (5b), style d'impression partagé (`paper_style`, `stars`, `bracket`, `repel_labels`, `save_paper`) |
| `fig5_matrices.py` | adaptateur matrices, F19, F21 |
| `tests/` | tests unitaires et de bout en bout (`python -m pytest -q tests` depuis `Fig5/`) |

## États d'éveil (`BIDS_CONFIGS[...]["state"]`)

| dataset | état |
|---|---|
| Human_Castrillon, Human_Babayan, Macaque_BenHamed, Marmoset_Tian, Dog | éveillé |
| Macaque_ZhuGarin | isoflurane |
| Marmoset_Hori | par run : `task-restawake` / `task-restanest` (l'état est aussi la session) |
| Mouse | par sujet : `MEDISO` → médétomidine + isoflurane, `ISO` → isoflurane, `AW` → éveillé ; l'animal est le suffixe (`jgrAesISOc11L` → `c11l`) |
| Mouselemur, Rat | anesthésié |

## Analyse matrices (supplémentaire)

- Fichiers lus : `EDNIxCSCLR_<ATLAS_LEVEL−1>_run_<i>_correlation_matrix.csv`, `<i>` résolu en run BIDS via `*_filtered.json` ; vérifier `02_selection/matrix_run_mapping.xlsx`.
- Profil : arêtes entre les régions communes de la seed (`|intra`, `|inter`, `|homotopic`).
- Doublons acq-AP / acq-PA de même numéro de run (Babayan) : le dernier fichier est gardé, comme dans l'analyse seed ; inclure `acq` dans l'identifiant de run reste à faire.
