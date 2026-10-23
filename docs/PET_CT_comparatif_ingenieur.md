# PET / CT : journal comparatif, version de l'ingénieur vs version corrigée

Comparaison entre le code PET/CT **tel que reçu** (commit `91ebd01`, « last version, not tested ») et la version corrigée de la branche `pet-ct-import` (commit `9352616`). Le détail ligne à ligne des erreurs est dans `docs/PET_CT_review_journal.md` (en anglais) ; ce document donne la vue d'ensemble, les mesures et ce qui reste à valider. Reproduire : `git diff 91ebd01 HEAD -- modalities/pet modalities/ct`.

## 1. Verdict

- **Conception** : celle de l'ingénieur est conservée (traitement dans l'espace PET natif, atlas et masques ramenés dans cet espace, arbre d'appels `pet3b.svg`, scénarios `pet4.svg`). Les modules gardent leurs noms et leurs rôles.
- **Exécution** : sa version ne pouvait pas tourner (voir §2). Elle n'avait jamais été exécutée.
- **Science** : même après correction des erreurs bloquantes, une vingtaine de choix auraient donné des résultats faux sans message d'erreur (TAC de référence, SUV, décroissance, chaînes de transformations, recalage du mouvement, correction de volume partiel…).
- **Limite de ma version** : elle est testée sur fantômes analytiques (outils AFNI et Turku émulés, ANTs réel), **pas sur des données réelles**. La liste des points à vérifier sur les premières données est au §7.

## 2. Mesures objectives sur le code reçu

Mesures faites sur les deux commits, avec le même outillage (pyflakes, import des modules, contrôle AST de la liaison de chaque appel à la signature de la fonction appelée).

| mesure | ingénieur (`91ebd01`) | version corrigée (`9352616`) |
|---|---|---|
| modules PET/CT/Tools/atlases importables | **11 sur 29** (18 échecs) | **27 sur 27** (les 2 autres ont été supprimés volontairement) |
| noms indéfinis (pyflakes) | **18** | 0 |
| erreur de syntaxe | **1** (`ct/_0_Pipeline_launcher.py:42`, `:` manquant) | 0 |
| appels dont les arguments ne se lient pas à la signature | non évaluable (le contrôle s'arrête dès le premier import en échec) | 0 sur 369 appels vérifiés |
| imports vers des modules inexistants | `anat.freesurfer`, `Tools.filename`, `Tools.setTransfo`, `tools`, `roi_signal` (import relatif implicite) | 0 |
| lancement par `python Launcher_...py` | échoue avant le pipeline (import d'un launcher à l'import, `baseline` indéfini) | corrigé (imports après `sys.path.insert`) |
| tests automatiques PET/CT | **0** | **36** (fantôme de tête analytique) |
| lignes de code `modalities/pet` + `modalities/ct` | 3 470 | 4 234 (+22 %) |
| documentation | aucune page PET | `docs/pet_pipeline.md` (184 lignes), références de paramètres générées PET et CT, journal |
| interface graphique | non câblée | onglets PET et CT, étapes, QC, rapports |

Les 18 échecs d'import de l'ingénieur (première erreur de chaque module) : 13 imports de modules inexistants (`anat` 2, `Tools.filename` 8, `roi_signal` 1, `tools` 2), 2 `os.path.exist` (`sifFile`, `getsif`), 1 erreur de syntaxe (launcher CT), 1 `DataFrame.subject` à l'import du launcher PET, et 1 `tkinter` absent de mon environnement (importé dès le chargement ; la version corrigée ne l'importe qu'à l'ouverture d'une fenêtre).

Le nombre d'environ 40 erreurs bloquantes cité dans le journal détaillé est un décompte à la lecture ; les chiffres ci-dessus sont ceux d'outils automatiques et sont des minorants (pyflakes ne voit ni les mauvais noms de fichiers, ni les signatures dont l'appelant est dynamique).

## 3. Architecture et registrations

| thème | ingénieur | corrigé |
|---|---|---|
| lien avec l'anatomie | transformations natif → template écrites en dur, sous deux noms différents (le second jamais produit par EDNiX) | chaînes de `loop3/backtonative.get`, comme le pipeline fonctionnel : template de référence, template d'étude et session de référence (longitudinal) |
| composition PET → template | `[shift, acpc]` + natif → template, même liste et mêmes drapeaux dans les deux sens | le shift est déjà dans la matrice AC-PC, il n'est plus recomposé ; chaînes aller et retour distinctes (`transfo.py`) |
| sans anatomie | `print('Sorry not ready yet....')` | scénarios 2v1, 2v2, 3v1, 3v2 |
| scénarios de `pet4.svg` | 1v2a seul pensé | 1v1, 1v2a, 1v2b, 2v1, 2v2, 3v1, 3v2 (2v2 et 3v2 jusqu'au template d'étude, comme indiqué sur le schéma) |
| fiche de recalage | transformations recalculées à chaque étape depuis des noms de fichiers | un enregistrement JSON par session (`orig/matrices/…_desc-registration_pet.json`) lu par les étapes 7 à 10 |
| granularité | un recalage par run | un par session, sur la moyenne après correction entre runs |
| PET et anatomie dans le même espace (TEP-IRM, TEP/TDM) | non traité | transformation du monde de la réorientation suivie par le PET et le CT (`3drefit -orient` change le signe de l'origine, cf. AFNI `3drefit.c`) + ancrage sur la position de tête de `same_space_window` |
| recalages sur template | `SyN` / `Affine` ANTsPy (un seul étage affine au début : rotation de 12° manquée de 7 mm sur le fantôme) | `SyNRA` par défaut (< 1 mm) |

## 4. Prétraitement et mouvement

| thème | ingénieur | corrigé |
|---|---|---|
| lecture des images brutes | cherchées dans `pet/raw/` avec un autre nom : rien trouvé | chemins BIDS réels (`pet/`, `sub-<ID>_ses-…`) |
| réorientation | test inversé : seules les images à en-tête **faux** étaient réorientées | corrigé ; interpolation `-linear` au lieu de `-wsinc5` (dépassements des concentrations) |
| moyennes temporelles | `nii = nii / n` dans la boucle ; moyennes non pondérées de trames de durées différentes | moyennes pondérées par la durée (moyenne temporelle de la concentration) |
| correction du mouvement, trames ignorées | `if i > nb_void` : une trame de plus laissée non corrigée | `>=` ; les trames initiales reçoivent la transformation de la première trame recalée |
| métrique de recalage des trames | ANTsPy `Rigid` par défaut, sans masque : 0,4 à 2 mm de mouvement parasite sur un fantôme immobile | information mutuelle dans un masque de tête, 3 niveaux jusqu'à la pleine résolution, tous les voxels : 0,05–0,1 mm parasite, 0,2–0,4 mm d'erreur sur le mouvement relatif |
| référence | moyenne des trames non corrigées (floue), une passe | deux passes |
| entre runs | un run dynamique rééchantillonné deux fois (étape 4 puis 5) | transformations composées, un seul rééchantillonnage |
| trame sans signal (correctif du 22/10) | n/a | tolérée : reçoit la matrice de la trame précédente, avertissement ; erreur explicite si la référence ou le masque de recadrage est vide |
| contrôle qualité | aucun | `*_motion.tsv` et `.png` (translations, rotations, déplacement moyen) |

## 5. Quantification (le plus important pour une publication)

| thème | ingénieur | corrigé |
|---|---|---|
| TAC de référence | moyenne des moyennes régionales (une petite région pèse autant qu'une grande) | moyenne sur l'union des voxels |
| TAC de référence à partir d'un masque | **somme** de l'activité | moyenne de la concentration |
| TAC régionales et référence (étape 8) | calculées sur l'image sans crâne, voxels hors masque comptés comme zéros : BPnd' biaisé, `BPnd' = (a/s)(BPnd+1) − 1` | calculées sur l'image entière ; modèles voxel à voxel sur l'image sans crâne |
| temps des trames | comptés depuis `TimeZero` (BIDS autorise `TimeZero` = début du scan) | temps post-injection (`InjectionStart`) : sinon `SUV_window=[40,60]` tombe 70–90 min après l'injection et le DVR simulé passe de 2,5 à 2,53–2,65 |
| SUV | unité supposée Bq/mL, dose non décroissante jusqu'à l'heure de correction, aucune correction de décroissance pour les données non corrigées | unité lue (`Units`), dose corrigée de la décroissance, facteur exact moyenné sur la trame ; pas de SUV si l'isotope est inconnu (erreur de 31 % sinon) |
| glycémie | SUV multiplié par la glycémie en g/L, sans valeur de référence | `SUVglu = SUV × glycémie / 5,0 mmol/L`, gardé à part |
| `bfmsrtm -i=<isotope>` | donné pour toutes les images (option réservée aux données **non** corrigées de la décroissance : biais) | donné seulement si `ImageDecayCorrected` est faux |
| DFT / SIF | volumes en nombre de voxels (pas en mm³), temps tronqués à l'entier, unité toujours `kBq/mL` | mm³, temps exacts, unité de l'image ; NaN ignorés (un voxel NaN annulait une TAC entière) |
| `white+csf` | produit des masques (intersection vide) | union |
| volume partiel (PETPVC) | PSF = taille du voxel ; classes tissulaires binarisées, cervelet et bords érodés dans « reste » | PSF = résolution de l'image (`pvc_fwhm`) ; 4 classes à fractions sommant à 1 (GM, WM hors GM, LCR, reste) ; méthodes sans image de sortie (GTM, STC) refusées |
| isotope | `'^18^Fluorine'.split('^')` donne `-18` | période (0018,1075) puis libellé DICOM, ou `TracerRadionuclide` BIDS |
| Turku | décompression sur place (un échec laissait le pipeline sans `.nii.gz`) ; `imgdv -n=` avec `split('.')[0]` | copies temporaires, géométrie réécrite ; `-k2=` via `DV_k2ref` ; `desc-STRM` renommé `desc-SRTM` (Lammertsma & Hume 1996) |

## 6. Interface, documentation, tests, utilisation

- **Point d'entrée** : `preprocess_pet` (nom attendu par l'interface et le générateur de documentation). L'ancien `preprocess_data` lève une `TypeError` : un appel positionnel à 53 arguments se serait lié aux mauvais paramètres sans erreur.
- **Étapes interactives conservées** (freeview, fenêtre tkinter pour la région de référence) ; tkinter n'est importé qu'à l'ouverture d'une fenêtre.
- **Tests** : 36 tests sur fantôme. Précision mesurée : corrélation inter-trames 0,97 → > 0,995 pour des mouvements de 3–5 mm / 4–6° ; pose PET → anatomie à moins d'un demi-voxel PET ; SUV analytique ; temps post-injection ; ancrage ; fractions tissulaires ; reprises des étapes tardives ; tous les scénarios exécutés de bout en bout.
- **Suite complète** : les mêmes 23 échecs qu'avant le travail PET (21 `test_post_ttt`, `Launcher_func_Dog`, `test_qc_scan_table_and_path_safety`), sans rapport avec le PET.
- **Propositions de l'ingénieur hors PET/CT** : `Tools/getpath.py` (`ct()`) intégré sous forme de dictionnaires ; `FS_freeview.check_coreg` intégré (sans le renommage `FSview` → `check_seg`, qui casse `_loop3.py`) ; `Tools/roiStats.py` et `atlases/getlabels.py` réécrits ; `Tools/extract_filename.py` et `Tools/get_orientation.py` **non pris** (le second réintroduit des miroirs gauche-droite corrigés en septembre). Aucun fichier IRMf modifié.

## 7. Ce que cette comparaison ne démontre pas

1. **Aucun test sur données réelles.** Les fantômes valident la logique et la géométrie, pas le comportement sur tes images (bruit, artefacts, tables, unités réelles des sidecars).
2. **Outils AFNI et Turku émulés** dans les tests : les options ont été vérifiées contre leur documentation, pas exécutées.
3. **Premier retour réel** : le recalage des trames a levé « Total Mass of the image was zero » sur `sub-PD09659` (trame sans signal après recadrage), cas non couvert par les fantômes ; corrigé et testé depuis. D'autres cas de ce type sont probables sur les premiers sujets.
4. **À contrôler en priorité** : QC de l'image PET dans l'espace de l'anatomie pour le premier sujet de chaque scénario ; `orientation='RPS'` du launcher macaque ; unités et `InjectionStart` de CERMEP ; validité des régions de référence par espèce (littérature) ; graphiques de mouvement des premiers runs dynamiques.
5. **Code plus long (+22 %)** : la surface à relire est plus grande ; le gain vient des vérifications, des enregistrements de recalage et de la documentation.

## Références

- Lammertsma AA, Hume SP (1996) Simplified reference tissue model for PET receptor studies. NeuroImage 4:153-158.
- Logan J et al. (1996) Distribution volume ratios without blood sampling from graphical analysis of PET data. J Cereb Blood Flow Metab 16:834-840.
- Boellaard R et al. (2015) FDG PET/CT: EANM procedure guidelines for tumour imaging, v2.0. Eur J Nucl Med Mol Imaging 42:328-354 (SUV, décroissance).
- Thomas BA et al. (2016) PETPVC: a toolbox for performing partial volume correction techniques in PET. Phys Med Biol 61:7975-7993.
- Avants BB et al. (2011) An open source multivariate framework for n-tissue segmentation… NeuroImage 54:2033-2044 (ANTs, SyN).
- Mattes D et al. (2003) PET-CT image registration in the chest using free-form deformations. IEEE Trans Med Imaging 22:120-128 (information mutuelle).
- Norgaard M et al. (2022) PET-BIDS, an extension to the brain imaging data structure for positron emission tomography. Sci Data 9:65.
- Innis RB et al. (2007) Consensus nomenclature for in vivo imaging of reversibly binding radioligands. J Cereb Blood Flow Metab 27:1533-1539 (BPnd).
