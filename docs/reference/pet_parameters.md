# PET pipeline: parameters of `preprocess_pet`

Generated from `modalities/pet/_0_Pipeline_launcher.py` by `docs/generate_parameter_reference.py`; do not edit by hand.

## Steps (`Skip_step`)

Identifiers accepted by `Skip_step`: the listed steps are skipped.

| step | what it does |
|---|---|
| `1` | Frame timing, SIF and mean image of the raw runs; reorientation of the header (as the anatomy) |
| `2` | PET study template (creat_study_template_pet) |
| `3` | Crop drawn on the study template in freeview (box used by the motion correction) |
| `4` | Motion correction within each dynamic run (two passes, rigid) |
| `5` | Motion correction between the runs of the session; mean image of the session |
| `6` | Registration to the anatomy / CT / template (scenarios of pet4.svg), quality control; atlases and masks in the PET space |
| `bridge` | Scenario 2v2: PET study template in the CT study template, reference template registered to it, atlases to every session (ct_template_bridge) |
| `7` | Parametric images (SUV, SUVglu, ratio to the brain mean, time blocks) and regional values / TACs |
| `8` | Reference-region kinetic models: SRTM (bfmsrtm, imgbfbp, imgsrtm), Logan DVR (imgdv) |
| `pvc` | Partial-volume correction (PETPVC, needs pvc_fwhm) |
| `9` | PET during the fMRI runs (simultaneous PET-MR, do_pet_to_func) |
| `10` | Export to the AC-PC and template spaces (+ skull-stripped) |

## Parameters

| parameter | default | description |
|---|---|---|
| `Skip_step` | required | Steps NOT to run (numbers above, or 'pvc'). |
| `MAIN_PATH` | required | Folder of the EDNiX installation, with its trailing '/' as in the example launchers: Atlases_library and Singularity_library are looked for in os.path.dirname(MAIN_PATH) (modalities/anat/set_launcher.py). |
| `bids_dir` | required | Path to the BIDS dataset. |
| `species` | required | Species of the dataset ('Macaque', 'Human'...), as in the anatomical pipeline. |
| `allinfo_study_c` | required | BIDS table of the study (Tools.load_bids.Load_BIDS_to_pandas with modalities=['pet']). |
| `endpet` | required | Glob pattern of the raw PET images in sub-*/ses-*/pet, e.g. '*_trc-FDG_*pet.nii.gz'. |
| `tracer` | required | Tracer name, as in the BIDS trc- entity (FDG, RACLO, PE2I...); names the study template and the session files. |
| `animalPosition` | required | Position of the animal and position declared to the scanner, as in the anatomical pipeline (['']: humanlike / PatientPosition of the sidecar). One value, or one per run. |
| `humanPosition` | required | Position of the animal and position declared to the scanner, as in the anatomical pipeline (['']: humanlike / PatientPosition of the sidecar). One value, or one per run. |
| `orientation` | required | Orientation code forced on the PET header (as `orientation` of the anatomical pipeline); '' to derive it from the positions. |
| `T_PET` | required | Anatomical image the PET is registered to ('T1w' or 'T2w'). |
| `type_norm` | required | Anatomical contrast of the normalisation of the anatomical pipeline (same value as there). |
| `creat_study_template` | required | Declaration, not an action: whether the ANATOMICAL pipeline used an anatomical study template (same value as in the anatomical launcher). Tells where the anatomy -> template transforms are; builds nothing for the PET (see creat_study_template_pet). |
| `coregistration_longitudinal` | required | Declaration, not an action: whether the ANATOMICAL pipeline registered the sessions of a subject to its reference session (same value as in the anatomical launcher). The PET sessions are not registered to each other: each one goes to the anatomy of its session. |
| `IhaveanANAT` | `True` | The subjects have an anatomical MRI processed by the anatomical pipeline (scenario 1 of pet4.svg). |
| `anat_pet_same_space` | `False` | PET and anatomy acquired in the same session and scanner space (hybrid PET-MR, scenario 1v1): the PET is reoriented like the anatomy and follows its AC-PC matrix, without registration. |
| `IhaveaCT` | `False` | A CT processed by the CT pipeline (modalities/ct) is used: PET -> CT -> anatomy (1v2b) or CT template (2v1, 2v2). |
| `ct_pet_same_space` | `True` | PET and CT share the scanner space (PET/CT scanner): no PET -> CT registration. |
| `template_PET` | `''` | PET template in the space of the reference template (scenario 3v1, no anatomy, no CT). |
| `creat_study_template_pet` | `False` | Action: build a PET study template (step 2) from the session means. Used for the crop of crop_for_MoCo='auto' (drawn once on it, step 3, and brought to every session) and, without anatomy, CT or PET template, as the common space (3v2, no atlas). It is not registered to the anatomy: with an anatomy the sessions are registered to their anatomy as without it. False: a crop of a study template left on disk by an earlier run is not used. |
| `studyT_transform` | `'SyNRA'` | ANTs transform of the PET study template registrations: 'SyNRA' (default: rigid + affine + SyN stages). ANTsPy's 'SyN' and 'Affine' start with a single affine stage, which missed a 12-degree rotation by about 7 mm on the test phantom (< 1 mm with 'SyNRA' or the linear 'TRSAA'). |
| `isotope` | `'unknown'` | Radionuclide ('F-18', 'C-11'...); 'unknown': TracerRadionuclide of the sidecar, else from the tracer name. |
| `PET_header_wrong` | `False` | The reconstruction wrote a wrong orientation: the orientation code of the raw anatomy (or CT) of the same scanner space is copied before the reorientation. |
| `deoblique` | `'WARP'` | 'WARP' resamples an oblique PET on a cardinal grid (3dWarp, as the anatomy); 'header' drops the obliquity. |
| `dcmfolder` | `''` | DICOM folder of the PET and file extension, to read the isotope and scan start for the SIF (optional). |
| `endfiledcm` | `'dcm'` | DICOM folder of the PET and file extension, to read the isotope and scan start for the SIF (optional). |
| `nb_void` | `0` | Number of first frames not registered by the motion correction (they get the transform of the first registered frame). |
| `crop_for_MoCo` | `'auto'` | Head box used by the motion correction and the registration: 'auto' (crop of the PET study template brought to every session when creat_study_template_pet=True, else drawn on the session), 'manual' (drawn on every session, freeview), 'none' (whole field of view). A session crop that no longer matches the reoriented image of its session (step 1 run again with another orientation) is drawn again. The box is shown in pet/QC/<session>_desc-cropbox_pet.png (and <BIDS>/QC/pet_crop), with a warning when it misses the activity. |
| `redo_crop` | `False` | Draw the crops again: the crop of the PET study template (step 3) and the crop of every session (<session>_desc-crop_pet.nii.gz in pet/orig/preprocessing/forMC, drawn in freeview at the first of steps 4, 5, 6 that runs). They are kept from one run to the next otherwise; a session crop that no longer matches the reoriented image of its session (step 1 run again with another orientation) is drawn again. |
| `MoCo_transform` | `'Rigid'` | ANTs transform of the motion correction ('Rigid'). |
| `coreg_metric` | `'mattes'` | Metric of the rigid registrations (motion, PET -> anatomy / CT): 'mattes' (mutual information). |
| `type_of_transform_template` | `'SyNRA'` | ANTs transform of the PET -> PET template registration (scenario 3v1): 'SyNRA' (default: rigid + affine + SyN stages). ANTsPy's 'SyN' and 'Affine' start with a single affine stage, which missed a 12-degree rotation by about 7 mm on the test phantom (< 1 mm with 'SyNRA' or the linear 'TRSAA'). |
| `do_manual` | `False` | Manual correction of the PET -> anatomy registration in freeview (scenario 1v2a). |
| `same_space_window` | `None` | PET in the scanner space of the anatomy (1v1) or of the CT (ct_pet_same_space): [start, end] (min after the injection) of the anatomy / CT acquisition. The motion-corrected PET (average head position) is brought back to the mean head position of the frames of that window; None = the first 10 min of the first run from its first registered frame (CT acquired before the PET). |
| `coreg_transform` | `'Rigid'` | ANTsPy type_of_transform of the PET -> anatomy / CT registration: 'Rigid' (default; a PET and an MRI of the same subject differ by a rigid motion). 'Similarity' / 'Affine' add scale parameters that can absorb a misregistration (PET brain shrunk into the MRI brain); non-linear transforms ('SyNRA'...) are not physiological for one subject. A warning is written for both. |
| `coreg_init` | `'auto'` | Initialisation of the PET -> anatomy / CT registration: 'auto' (default: centres of mass of the masked images, centres of the fields of view ('FOV to FOV'), AC-PC matrix of the anatomy (the PET is near the native anatomy), identity (headers), search over the rotations (antsAffineInitializer: large rotations); the start with the highest mutual information that covers the brain is kept, scores in pet/orig/matrices/*candidates.json), or one of 'com', 'fov', 'acpc', 'identity', 'search'. |
| `pet_coreg_mask` | `'none'` | Mask of the PET for the registration (after the crop box of the motion correction): 'none' (default), 'threshold' (voxels above pet_mask_percentile of the positive voxels, largest connected component, holes filled, dilated pet_mask_dilate_mm), 'anat' (second pass: brain of the anatomy brought into the PET by a first registration, dilated pet_mask_dilate_mm, then registration of the PET brain only; removes the extracerebral uptake from the metric), or a mask file in the PET space. |
| `pet_mask_percentile` | `30.0` | Percentile of the positive PET voxels for pet_coreg_mask='threshold' (30). |
| `pet_mask_dilate_mm` | `3.0` | Dilation (mm) of the PET mask ('threshold', 'anat') and of the brain of the anatomy in the second pass (3). |
| `coreg_target` | `'anat'` | Fixed image of the PET -> anatomy registration (1v2a): 'anat' (default, AC-PC head image), 'brain' (skull-stripped anatomy) or 'synthetic' (FDG-like image of the anatomy: grey matter 1, white matter 0.35, CSF 0.05, smoothed to synthetic_fwhm; Ashburner & Friston 1997). FDG only. |
| `synthetic_fwhm` | `None` | FWHM (mm) of the synthetic target; None = twice the largest PET voxel. |
| `ct_template_bridge` | `True` | Scenario 2v2 (CT study template, no anatomy): step 'bridge' links the CT study template to the reference template through a PET study template (modalities/pet/bridge.py): steps 4-6 for every session, mean of the session means in the CT study template, skull-stripped reference template registered to it (multi-start similarity then affine, mutual information), brain mask and atlases brought back to every session, then steps 7-10 with the atlases. True by default; False keeps 2v2 without atlas. FDG-like tracers (brain contrast). |
| `bridge_transform` | `'Affine'` | Reference template -> PET study template: 'Affine' (default) or 'SyNRA' (affine then SyN). |
| `SUV` | `True` | Compute the SUV (body weight and injected radioactivity from the sidecar / sessions.tsv / participants.tsv) and, with the blood glucose, the glucose-corrected SUV. |
| `bratio` | `True` | Compute the ratio of the PET to its mean over the brain mask. |
| `SUV_window` | `None` | [start, end] in minutes after the injection (BIDS frame times minus InjectionStart) of the frames averaged for the SUV and ratio images of a dynamic run (frames whose middle is in the window); None = all frames. |
| `Do_cut` | `False` | Also write the mean of every time block of `duration` minutes (divided by its brain mean). |
| `duration` | `10` | Duration (min) of the time blocks of Do_cut (from the first frame; a frame belongs to the block that contains its middle). |
| `srtm` | `False` | Reference-region kinetic models of the dynamic runs (step 8): SRTM (regional bfmsrtm, voxel-wise imgbfbp, imgsrtm) and Logan DVR (imgdv) with the Turku PET Centre tools. |
| `reference_signal` | `'unknown'` | Reference region (atlases/getlabels.py): {'type': 'atlas', 'atlas': 'EDNIxCSCLR', 'level': 0, 'labels': [7, 1007] or a name pattern or 'default', 'name': 'Cer'}, {'type': 'mask', 'name': 'white'}, or 'unknown' for the interactive choice (saved for the next subjects). |
| `DV_time` | `0` | Start t* (min after the injection) of the linear part of the Logan plot (imgdv); 0 = no DVR image. |
| `DV_k2ref` | `None` | Population average k2' (min-1) of the reference region for the Logan plot (imgdv -k2=); None = the C_ref/k2' term is omitted (simplified reference Logan, valid when C_ref / C_T is constant after DV_time). |
| `kinetic_end` | `None` | End (min after the injection) of the kinetic fits; None = all frames. |
| `pvc_fwhm` | `None` | Spatial resolution of the reconstructed images (PSF FWHM, mm: scanner, reconstruction and post-reconstruction filter, sqrt(FWHM_scanner^2 + FWHM_filter^2) for Gaussian kernels): partial-volume correction (PETPVC) of the static images; None = no correction. |
| `pvc_method` | `'RBV'` | PETPVC method that writes an image ('RBV', 'MG', 'LABBE+RBV', 'IY'...; see PVC_METHODS). |
| `do_pet_to_func` | `False` | Simultaneous PET-MR: mean PET during every fMRI run of the session (step 9). |
| `endfmri` | `'*_bold.nii.gz'` | Glob pattern of the raw fMRI runs in sub-*/ses-*/func (step 9). |
| `T1_eq` | `0` | Number of first fMRI volumes removed by the functional pipeline (step 9). |
| `export_resolution` | `'pet'` | Grid of the exported images: 'pet' (target space at the PET voxel size) or 'target' (anatomy / template grid). |
| `export_dynamic` | `False` | Also export the motion-corrected dynamic images (large files). |
| `selected_atlases` | `'all'` | Atlases of the regional analysis: 'all' (every level of every atlas of atlas_followers) or [[name, level or 'all'], ...]. |
| `atlas_followers` | `[['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]]` | Atlases of the anatomical pipeline ([[names], [types], [levels], [...]]), as there. |
| `reference` | `'EDNiX'` | Reference template ('EDNiX'). |
| `use_erode_WM_func_masks` | `True` | Use the eroded white-matter / ventricle masks (reference region 'white', 'csf'). |
| `use_erode_V_func_masks` | `True` | Use the eroded white-matter / ventricle masks (reference region 'white', 'csf'). |
| `list_to_keep` | `[]` | (subject, session) pairs to process / exclude. |
| `list_to_remove` | `[]` | (subject, session) pairs to process / exclude. |
| `overwrite_option` | `True` | Not used: kept for the launchers of the other pipelines (the PET / CT steps always rewrite their outputs). |
| `Need_check_orient` | `False` | Open freeview to check the orientation, the motion correction and the registrations. |
| `Need_check_MC` | `False` | Open freeview to check the orientation, the motion correction and the registrations. |
| `Need_check_coreg` | `False` | Open freeview to check the orientation, the motion correction and the registrations. |
| `doBALSAspace` | `False` | Same value as in the anatomical pipeline. |

`preprocess_pet` also accepts other keyword arguments, which it ignores: a misspelled parameter name is therefore reported as a warning, not an error. The interface and the unit tests check the example launchers for this.
