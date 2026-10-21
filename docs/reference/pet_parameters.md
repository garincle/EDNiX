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
| `6` | Registration to the anatomy / CT / template (scenarios of pet4.svg); atlases and masks in the PET space |
| `7` | Parametric images (SUV, SUVglu, ratio to the brain mean, time blocks) and regional values / TACs |
| `8` | Reference-region kinetic models: SRTM (bfmsrtm, imgbfbp, imgsrtm), Logan DVR (imgdv) |
| `pvc` | Partial-volume correction (PETPVC, needs pvc_fwhm) |
| `9` | PET during the fMRI runs (simultaneous PET-MR, do_pet_to_func) |
| `10` | Export to the AC-PC and template spaces (+ skull-stripped), quality control |

## Parameters

| parameter | default | description |
|---|---|---|
| `Skip_step` | required | Steps NOT to run (numbers above, or 'pvc'). |
| `MAIN_PATH` | required | Folder of the EDNiX installation (the one that contains Atlases_library). |
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
| `creat_study_template` | required | Whether the anatomical pipeline used a study template (same value as there; anatomy -> template bridge). |
| `coregistration_longitudinal` | required | Whether the anatomical pipeline used a longitudinal registration (same value as there). |
| `IhaveanANAT` | `True` | The subjects have an anatomical MRI processed by the anatomical pipeline (scenario 1 of pet4.svg). |
| `anat_pet_same_space` | `False` | PET and anatomy acquired in the same session and scanner space (hybrid PET-MR, scenario 1v1): the PET is reoriented like the anatomy and follows its AC-PC matrix, without registration. |
| `IhaveaCT` | `False` | A CT processed by the CT pipeline (modalities/ct) is used: PET -> CT -> anatomy (1v2b) or CT template (2v1, 2v2). |
| `ct_pet_same_space` | `True` | PET and CT share the scanner space (PET/CT scanner): no PET -> CT registration. |
| `template_PET` | `''` | PET template in the space of the reference template (scenario 3v1, no anatomy, no CT). |
| `creat_study_template_pet` | `False` | Build a PET study template (step 2): used to crop the images for the motion correction and, without anatomy, CT or PET template, as the common space (3v2, no atlas). |
| `studyT_transform` | `'SyN'` | ANTs transform of the PET study template registrations ('SyN'). |
| `isotope` | `'unknown'` | Radionuclide ('F-18', 'C-11'...); 'unknown': TracerRadionuclide of the sidecar, else from the tracer name. |
| `PET_header_wrong` | `False` | The reconstruction wrote a wrong orientation: the orientation code of the raw anatomy (or CT) of the same scanner space is copied before the reorientation. |
| `deoblique` | `'WARP'` | 'WARP' resamples an oblique PET on a cardinal grid (3dWarp, as the anatomy); 'header' drops the obliquity. |
| `dcmfolder` | `''` | DICOM folder of the PET and file extension, to read the isotope and scan start for the SIF (optional). |
| `endfiledcm` | `'dcm'` | DICOM folder of the PET and file extension, to read the isotope and scan start for the SIF (optional). |
| `nb_void` | `0` | Number of first frames not registered by the motion correction (they get the transform of the first registered frame). |
| `crop_for_MoCo` | `'auto'` | Head box used by the motion correction: 'auto' (crop of the study template if it exists, else drawn on the session), 'manual' (drawn on every session, freeview), 'none' (whole field of view). |
| `redo_crop` | `False` | Draw the crops again (they are kept from one run to the next otherwise). |
| `MoCo_transform` | `'Rigid'` | ANTs transform of the motion correction ('Rigid'). |
| `coreg_metric` | `'mattes'` | Metric of the rigid registrations (motion, PET -> anatomy / CT): 'mattes' (mutual information). |
| `type_of_transform_template` | `'SyN'` | ANTs transform of the PET -> PET template registration (scenario 3v1), e.g. 'SyN'. |
| `do_manual` | `False` | Manual correction of the PET -> anatomy registration in freeview (scenario 1v2a). |
| `SUV` | `True` | Compute the SUV (body weight and injected radioactivity from the sidecar / sessions.tsv / participants.tsv) and, with the blood glucose, the glucose-corrected SUV. |
| `bratio` | `True` | Compute the ratio of the PET to its mean over the brain mask. |
| `SUV_window` | `None` | [start, end] in minutes of the frames averaged for the SUV and ratio images of a dynamic run (frames whose middle is in the window); None = all frames. |
| `Do_cut` | `False` | Also write the mean of every time block of `duration` minutes (divided by its brain mean). |
| `duration` | `10` | Duration (min) of the time blocks of Do_cut. |
| `srtm` | `False` | Reference-region kinetic models of the dynamic runs (step 8): SRTM (regional bfmsrtm, voxel-wise imgbfbp, imgsrtm) and Logan DVR (imgdv) with the Turku PET Centre tools. |
| `reference_signal` | `'unknown'` | Reference region (atlases/getlabels.py): {'type': 'atlas', 'atlas': 'EDNIxCSCLR', 'level': 0, 'labels': [7, 1007] or a name pattern or 'default', 'name': 'Cer'}, {'type': 'mask', 'name': 'white'}, or 'unknown' for the interactive choice (saved for the next subjects). |
| `DV_time` | `0` | Start (min) of the linear part of the Logan plot (imgdv); 0 = no DVR image. |
| `kinetic_end` | `None` | End (min) of the kinetic fits; None = all frames. |
| `pvc_fwhm` | `None` | Spatial resolution of the scanner (PSF FWHM, mm): partial-volume correction (PETPVC) of the static images; None = no correction. |
| `pvc_method` | `'RBV'` | PETPVC method ('RBV', 'GTM', 'MG'...). |
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
| `overwrite_option` | `True` | Overwrite the AFNI outputs. |
| `Need_check_orient` | `False` | Open freeview to check the orientation, the motion correction and the registrations. |
| `Need_check_MC` | `False` | Open freeview to check the orientation, the motion correction and the registrations. |
| `Need_check_coreg` | `False` | Open freeview to check the orientation, the motion correction and the registrations. |
| `doBALSAspace` | `False` | Same value as in the anatomical pipeline. |

`preprocess_pet` also accepts other keyword arguments, which it ignores: a misspelled parameter name is therefore reported as a warning, not an error. The interface and the unit tests check the example launchers for this.
