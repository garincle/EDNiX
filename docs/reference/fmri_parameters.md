# Functional pipeline: parameters of `preprocess_fmri`

Generated from `modalities/fMRI/_0_Pipeline_launcher.py` by `docs/generate_parameter_reference.py`; do not edit by hand.

## Steps (`Skip_step`)

Identifiers accepted by `Skip_step`: the listed steps are skipped.

| step | what it does |
|---|---|
| `1` | Native-space preprocessing: slice timing, despiking, motion correction, FD censoring |
| `2` | Run-to-reference coregistration; distortion field from the field maps (applied at step 5) |
| `3` | Functional brain mask |
| `itk_1` | Manual check of the functional masks in ITK-SNAP |
| `5` | Anatomy to functional registration; one-shot resampling of every volume (motion, distortion, registration) |
| `itk_2` | Manual check of the functional data in ITK-SNAP |
| `6` | ICA (MELODIC) cleaning, interactive component selection |
| `7` | Nuisance regression, temporal filtering, smoothing |
| `8` | Functional data to anatomical space |
| `9` | Functional data to template space |
| `10` | Correlation matrices between atlas regions |
| `11` | Seed-based analysis |
| `12` | Functional quality control |
| `13` | Connectivity quality control (FC specificity) |
| `itk_3` | Manual check of the final data in ITK-SNAP |
| `Clean` | Delete intermediate files |

## Parameters

| parameter | default | description |
|---|---|---|
| `Skip_step` | required | Steps to skip in the pipeline (integers or descriptive labels). |
| `MAIN_PATH` | required | Base path for output and intermediate results. |
| `bids_dir` | required | Path to BIDS-formatted dataset. |
| `species` | required | Species label ('Human', 'Rat', 'Mouse', etc.). |
| `allinfo_study_c` | required | DataFrame containing BIDS subject/session information. |
| `endfmri` | required | File suffixes for fMRI, JSON, and derived maps. |
| `endjson` | required | File suffixes for fMRI, JSON, and derived maps. |
| `endmap` | required | File suffixes for fMRI, JSON, and derived maps. |
| `resting_or_task` | required | Type of fMRI acquisition ('resting' or 'task'). |
| `animalPosition` | required | Orientation information for animals or humans during acquisition. |
| `humanPosition` | required | Orientation information for animals or humans during acquisition. |
| `orientation` | required | Target orientation for preprocessing (e.g., 'LPI'). |
| `Slice_timing_info` | required | Slice acquisition timing information. |
| `TfMRI` | required | Suffix of the anat use to register functional MRI data ('T2w', etc.). |
| `type_norm` | required | Anatomical type_norm contrast used originally for normalization ('T1w' or 'T2w', should be the same as in type_norm anat). |
| `creat_study_template` | required | Whether study-specific anat template has been generated in the anat regitration (should be the same as in creat_study_template anat) |
| `anat_func_same_space` | required | Flag to indicate if anat and functional images share EXACTLY the same space. |
| `coregistration_longitudinal` | required | Whether to apply longitudinal coregistration across sessions. (should be the same as in coregistration_longitudinal anat) |
| `Method_mask_func` | required | Method for generating functional masks. |
| `extra_erode` | `0` | Number of voxels by which the edge of the functional brain mask is eroded, removing dark edge voxels of the mean image; 0 disables it. |
| `do_anat_to_func` | `True` | Perform anat-to-functional registration (default: True). |
| `folderforTemplate_Anat` | `''` | Path to anat templates. |
| `IhaveanANAT` | `True` | Whether anat images are available (default: True). |
| `ntimepoint_treshold` | `100` | Minimum number of timepoints required to include a run. |
| `REF_int` | `0` | Reference volume index for registration. |
| `T1_eq` | `5` | Number of T1-equilibration volumes to remove. |
| `correction_direction` | `'Auto'` | Phase-encoding (unwarping) direction for FSL FUGUE: 'x', 'x-', 'y', 'y-', 'z' or 'z-'; 'Auto' reads PhaseEncodingDirection from the JSON sidecar. Also restricts the non-linear func-to-anat registration (do_anat_to_func) to that axis. |
| `overwrite_option` | `True` | Overwrite existing files (default: True). |
| `DwellT` | `'Auto'` | Effective echo spacing in seconds, used by FUGUE; 'Auto' reads EffectiveEchoSpacing, or TotalReadoutTime / (ReconMatrixPE - 1), from the JSON sidecar. |
| `SED` | `'Auto'` | Slice-encoding direction; 'Auto' reads SliceEncodingDirection from the JSON sidecar (logged; the registration restriction uses the phase-encoding direction). |
| `TR` | `'Auto'` | Repetition time in seconds; 'Auto' reads RepetitionTime from the JSON sidecar (or estimates it from SliceTiming, with a warning). |
| `TRT` | `'Auto'` | Total readout time in seconds, used by TOPUP; 'Auto' reads TotalReadoutTime from the JSON sidecar. |
| `nb_ICA_run` | `20` | Number of ICA components for denoising. |
| `ICA_cleaning` | `'Skip'` | ICA cleaning method ('Skip', 'MELODIC', etc.). |
| `costAllin` | `'lpa'` | Registration cost function for 3dAllineate alignment ('lpa', etc.). ## will be removed |
| `doWARPonfunc` | `False` | How oblique functional images are de-obliqued: 'WARP' resamples them (AFNI 3dWarp -deoblique), 'header' keeps the data and only rewrites the header; other values (default False) leave oblique images as they are. |
| `registration_type` | `'origin'` | Parameter set of the func-to-anat registration when do_anat_to_func is True: 'origin' (default), 'fast' (fewer iterations) or 'V2'. The transform itself is type_of_transform (a non-linear one, e.g. 'SyN', is needed for a distortion-like correction). |
| `type_of_transform` | `'BOLDAffine'` | Type of transformation for alignment ('BOLDAffine', 'Rigid', 'SyN', etc.). |
| `n_for_ANTS` | `'Lanczos'` | Interpolation method for ANTs registration ('Lanczos', 'Linear', etc.). |
| `aff_metric_ants` | `'meansquares'` | Affine registration metrics for ANTs. |
| `aff_metric_ants_Transl` | `'mattes'` | Affine registration metrics for ANTs. |
| `dilate_mask` | `0` | Number of voxels to dilate functional masks. |
| `list_to_keep` | `[]` | Subject/session tuples to include or exclude. |
| `list_to_remove` | `[]` | Subject/session tuples to include or exclude. |
| `atlas_followers` | `[[], [], [], []]` | Multi-level atlas selection for processing. |
| `reference` | `'EDNiX'` | Reference dataset/session ('EDNiX' by default). |
| `post_treatment_method` | `'Grandjean'` | Method for post-processing (default: 'Grandjean'). |
| `band` | `'0.01 0.1'` | Temporal filter band (default: '0.01 0.1'). |
| `blur` | `0` | Gaussian smoothing kernel in mm. |
| `do_not_correct_signal` | `False` | Skip temporal or spatial signal corrections. |
| `extract_WM` | `True` | Nuisance signals of the white matter (Wmask), the ventricles (Vmask) and the brain (global signal). They select the regressors of the 'AFNI' method only. 'Grandjean' always uses the mean signal of the white-matter and ventricle masks (when the masks exist) and 'fMRIPrep_style' always uses aCompCor of both; extract_GS adds the global signal to 'AFNI' and 'fMRIPrep_style'. |
| `extract_Vc` | `False` | Nuisance signals of the white matter (Wmask), the ventricles (Vmask) and the brain (global signal). They select the regressors of the 'AFNI' method only. 'Grandjean' always uses the mean signal of the white-matter and ventricle masks (when the masks exist) and 'fMRIPrep_style' always uses aCompCor of both; extract_GS adds the global signal to 'AFNI' and 'fMRIPrep_style'. |
| `extract_GS` | `False` | Nuisance signals of the white matter (Wmask), the ventricles (Vmask) and the brain (global signal). They select the regressors of the 'AFNI' method only. 'Grandjean' always uses the mean signal of the white-matter and ventricle masks (when the masks exist) and 'fMRIPrep_style' always uses aCompCor of both; extract_GS adds the global signal to 'AFNI' and 'fMRIPrep_style'. |
| `use_erode_WM_func_masks` | `True` | Apply erosion to WM and ventricle masks. |
| `use_erode_V_func_masks` | `True` | Apply erosion to WM and ventricle masks. |
| `normalize` | `'Skip'` | Normalization method or 'Skip'. |
| `acompcor_mode` | `'fixed'` | 'fixed' (default): acompcor_n_fixed components per tissue (Behzadi et al. 2007); 'parallel': the components that stand above chance (parallel analysis against a phase-randomised null), at most acompcor_max; 'adaptive': as many components as needed to explain acompcor_variance of the variance, at most acompcor_max (the cap is usually reached: it then acts as a fixed number). Only for post_treatment_method='fMRIPrep_style'. The log reports, for every tissue, the number found by parallel analysis. |
| `acompcor_n_fixed` | `5` | Number of aCompCor components per tissue in 'fixed' mode (default 5). |
| `acompcor_variance` | `0.50` | Fraction of the variance to explain in 'adaptive' mode (default 0.50). |
| `acompcor_max` | `15` | Cap on the number of aCompCor components per tissue in 'adaptive' and 'parallel' modes (default 15). |
| `motion_mask` | `'None'` | Brain mask for the motion estimation of step 1. 'None' (default): the whole image is used (ANTs rigid, 'BOLDRigid'). 'previous': the functional mask written by step 3 of an earlier execution (run the pipeline once, then again from step 1). Or the path of a mask in the frame of the reoriented functional series. Needed when tissue that does not move with the brain (skull, muscle, skin) dominates the image: the log of step 1 compares the two estimations. |
| `motion_mask_erode` | `3` | Voxels removed from the border of the motion mask (default 3): it must stay inside the brain whatever the head does, so use at least the largest head displacement in voxels. |
| `motion_method` | `'ants'` | Program that estimates the head motion in step 1: 'ants' (default, ANTs rigid), '3dvolreg' (AFNI) or 'mcflirt' (FSL). The transforms of 3dvolreg and MCFLIRT are converted into the ANTs matrices used by step 5 (one interpolation per volume with the registration and distortion correction), after a check that the converted matrices reproduce the registered volumes written by the program (the run stops otherwise). Needs the python package nitransforms. The brain mask (motion_mask) is used with 'ants' only. |
| `selected_atlases_matrix` | `'all'` | Atlas selection for matrix-based analyses ('all', 'default', etc.). |
| `wanted_level_matrix` | `'all'` | Atlas selection for matrix-based analyses ('all', 'default', etc.). |
| `matrix_fit` | `['covariance', 'correlation', 'partial correlation']` | Connectivity estimators (nilearn ConnectivityMeasure), e.g. 'correlation', 'partial correlation'. |
| `selected_atlases_SBA` | `'default'` | Atlas and file settings for seed-based analyses. |
| `panda_files_SBA` | `'default'` | Atlas and file settings for seed-based analyses. |
| `SBAspace` | `['func', 'anat', 'atlas']` | Space(s) for SBA outputs ('func', 'anat', 'atlas'). |
| `erod_seed` | `True` | Erode seed masks for SBA. |
| `smoothSBA` | `False` | Smooth SBA results. |
| `intra_thresh` | `0.240` | Connectivity-specificity thresholds of step 13 (Grandjean et al. 2020): a run is 'Specific' when homotopic connectivity exceeds intra_thresh and exceeds the cross-network connectivity by more than delta_thresh; 'Unspecific' when it exceeds intra_thresh but not by delta_thresh; 'Spurious' when only the cross-network connectivity exceeds intra_thresh; 'No' otherwise. |
| `delta_thresh` | `0.05` | Connectivity-specificity thresholds of step 13 (Grandjean et al. 2020): a run is 'Specific' when homotopic connectivity exceeds intra_thresh and exceeds the cross-network connectivity by more than delta_thresh; 'Unspecific' when it exceeds intra_thresh but not by delta_thresh; 'Spurious' when only the cross-network connectivity exceeds intra_thresh; 'No' otherwise. |
| `oversample_map` | `False` | Oversample maps to higher resolution. |
| `use_cortical_mask_func` | `False` | Restrict analyses to cortical regions. |
| `n_cut` | `10` | Number of voxels to exclude at edges. |
| `threshold_val` | `10` | Intensity threshold for masking. |
| `doBALSAspace` | `False` | Human, chimpanzee and macaque: also bring the subject to the BALSA template space (S1200, CY29, MY19), resample its surfaces to fs_LR, and for the macaque build the BALSA atlas in native space. |

`preprocess_fmri` also accepts other keyword arguments, which it ignores: a misspelled parameter name is therefore reported as a warning, not an error. The interface and the unit tests check the example launchers for this.
