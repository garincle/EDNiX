# CT pipeline (for the PET pipeline): parameters of `preprocess_ct`

Generated from `modalities/ct/_0_Pipeline_launcher.py` by `docs/generate_parameter_reference.py`; do not edit by hand.

## Steps (`Skip_step`)

Identifiers accepted by `Skip_step`: the listed steps are skipped.

| step | what it does |
|---|---|
| `1` | Reorientation of the header (as the PET of the same PET/CT acquisition) |
| `2` | CT study template (creat_study_template_ct, no anatomy) |
| `3` | Registration to the anatomy / CT template / study template, record for the PET pipeline, QC |

## Parameters

| parameter | default | description |
|---|---|---|
| `Skip_step` | required | Steps NOT to run. |
| `MAIN_PATH` | required | Folder of the EDNiX installation, with its trailing '/' as in the example launchers: Atlases_library and Singularity_library are looked for in os.path.dirname(MAIN_PATH) (modalities/anat/set_launcher.py). |
| `bids_dir` | required | Path to the BIDS dataset. |
| `species` | required | Species of the dataset, as in the anatomical pipeline. |
| `allinfo_study_c` | required | BIDS table of the study (Tools.load_bids.Load_BIDS_to_pandas with modalities=['ct']). |
| `endct` | required | Glob pattern of the raw CT images in sub-*/ses-*/ct, e.g. '*_ct.nii.gz'. |
| `animalPosition` | required | Position of the animal and position declared to the scanner, as in the anatomical pipeline; use the same values as the PET of the same PET/CT acquisition. |
| `humanPosition` | required | Position of the animal and position declared to the scanner, as in the anatomical pipeline; use the same values as the PET of the same PET/CT acquisition. |
| `orientation` | required | Orientation code forced on the CT header; '' to derive it from the positions. |
| `T_CT` | required | Anatomical image the CT is registered to ('T1w' or 'T2w'). |
| `type_norm` | required | Anatomical contrast of the normalisation of the anatomical pipeline (same value as there). |
| `creat_study_template` | required | Whether the anatomical pipeline used a study template (same value as there). |
| `coregistration_longitudinal` | required | Whether the anatomical pipeline used a longitudinal registration (same value as there). |
| `IhaveanANAT` | `True` | The subjects have an anatomical MRI processed by the anatomical pipeline (scenario 1v2b). |
| `anat_ct_same_space` | `False` | CT and anatomy acquired in the same scanner space: the CT follows the AC-PC matrix of the anatomy. |
| `template_CT` | `''` | CT template in the space of the reference template (scenario 2v1, no anatomy). |
| `creat_study_template_ct` | `False` | Build a CT study template (step 2) and register the CTs to it (scenario 2v2, no anatomy, no CT template). |
| `studyT_transform` | `'SyNRA'` | ANTs transform of the study template registrations: 'SyNRA' (default: rigid + affine + SyN stages). ANTsPy's 'SyN' and 'Affine' start with a single affine stage, which missed a 12-degree rotation by about 7 mm on the test phantom (< 1 mm with 'SyNRA' or the linear 'TRSAA'). |
| `type_of_transform_template` | `'SyNRA'` | ANTs transform of the CT -> CT template registration (2v1): 'SyNRA' (default: rigid + affine + SyN stages). ANTsPy's 'SyN' and 'Affine' start with a single affine stage, which missed a 12-degree rotation by about 7 mm on the test phantom (< 1 mm with 'SyNRA' or the linear 'TRSAA'). |
| `coreg_metric` | `'mattes'` | Metric of the registrations ('mattes': mutual information). |
| `CT_header_wrong` | `False` | The reconstruction wrote a wrong orientation: the orientation code of the raw anatomy of the same scanner space is copied before the reorientation (with anat_ct_same_space). |
| `deoblique` | `'WARP'` | 'WARP' resamples an oblique CT on a cardinal grid (3dWarp, as the anatomy); 'header' drops the obliquity. |
| `do_manual` | `False` | Manual correction of the CT -> anatomy registration in freeview. |
| `coreg_transform` | `'Rigid'` | ANTsPy type_of_transform of the CT -> anatomy registration ('Rigid': the CT and the MRI of one subject differ by a rigid motion). |
| `coreg_init` | `'auto'` | Initialisation of the CT -> anatomy registration: 'auto' (centres of mass of the CT head without the bed and of the anatomy, centres of the fields of view, AC-PC matrix of the anatomy, identity, search over the rotations; the start with the highest mutual information covering the brain is kept, scores in ct/orig/matrices/*candidates.json), 'com', 'fov', 'acpc', 'identity' or 'search'. |
| `ct_brain_refine` | `False` | Second pass of the CT -> anatomy registration (1v2b) in the brain: CT in the brain window (0-80 HU, bone set to 0) within the brain of the anatomy brought into the CT, from the head registration. The grey / white matter and the ventricles drive it instead of the skull. False by default (no gain on the phantoms, where the head registration is already < 0.5 mm); to try when the QC shows the brain shifted inside the skull. |
| `Need_check_orient` | `False` | Open freeview to check the orientation and the registration. |
| `Need_check_coreg` | `False` | Open freeview to check the orientation and the registration. |
| `reference` | `'EDNiX'` | Reference template ('EDNiX'). |
| `atlas_followers` | `[['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]]` | Atlases of the anatomical pipeline, as there. |
| `use_erode_WM_func_masks` | `True` | Same values as in the PET pipeline (masks of the anatomy). |
| `use_erode_V_func_masks` | `True` | Same values as in the PET pipeline (masks of the anatomy). |
| `list_to_keep` | `[]` | (subject, session) pairs to process / exclude. |
| `list_to_remove` | `[]` | (subject, session) pairs to process / exclude. |
| `overwrite_option` | `True` | Not used: kept for the launchers of the other pipelines (the PET / CT steps always rewrite their outputs). |
| `doBALSAspace` | `False` | Same value as in the anatomical pipeline. |

`preprocess_ct` also accepts other keyword arguments, which it ignores: a misspelled parameter name is therefore reported as a warning, not an error. The interface and the unit tests check the example launchers for this.
