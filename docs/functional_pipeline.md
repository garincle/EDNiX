# Functional pipeline

Entry point: `modalities/fMRI/_0_Pipeline_launcher.preprocess_fmri`. Every parameter is listed in the [functional parameter reference](reference/fmri_parameters.md). Run the anatomical pipeline of the dataset first.

## Steps

| step | module | what it does |
|---|---|---|
| `1` | `_1_fMRI_preTTT_in_fMRIspace` | removal of the first volumes (`T1_eq`), despiking, slice timing, motion correction (ANTs), framewise displacement and censoring, mean image, bias correction |
| `2` | `_2_coregistration_to_norm` | field maps (TOPUP / FUGUE) and distortion field, orientation fix, registration of each run to the reference run |
| `3` | `_3_mask_fMRI` | functional brain mask (see [skullstripping.md](skullstripping.md)) |
| `itk_1` | | manual check of the masks in ITK-SNAP |
| `5` | `_5_anat_to_fMRI` | anatomy-to-function registration; one-shot resampling of every volume |
| `itk_2` | | manual check in ITK-SNAP |
| `6` | `_6_Melodic` | optional ICA cleaning (MELODIC, interactive component selection) |
| `7` | `_7_post_TTT` | nuisance regression, band-pass filtering, smoothing, normalisation |
| `8` | `_8_fMRI_to_anat` | functional data to anatomical space |
| `9` | `_9_coregistration_to_template_space` | functional data, masks and atlases in template space |
| `10` | `_10_Correl_matrix` | correlation matrices between atlas regions |
| `11` | `_11_Seed_base_many_regionsatlas` | seed-based correlation maps |
| `12` | `_12_fMRI_QC` | quality-control metrics (FD, DVARS, tSNR, GCOR...) |
| `13` | `_14_fMRI_QC_matrix` | connectivity specificity |
| `itk_3` | | manual check of the final data |
| `Clean` | | deletion of intermediate files |

## Timing

`TR`, `DwellT`, `SED` and `TRT` default to `'Auto'`: they are read from the JSON sidecar of the first run. `TR` is then always used as a number in seconds; when `RepetitionTime` is missing, it is estimated from `SliceTiming` (with a warning to check). A value given in the launcher (number or text) overrides the sidecar.

## Motion and censoring

Motion is estimated with ANTs (`BOLDRigid`) on each run. When the orientation fix changes the image frame, the motion matrices are re-estimated in the reoriented frame, so that rotations are applied about the right axes.

Framewise displacement (FD) converts rotations into millimetres on a sphere of the same volume as the brain, and the censoring threshold is scaled with the same radius: 0.5 mm for a 1200 cm³ human brain, about 0.2 mm for a macaque and 0.04 mm for a mouse. Censored volumes are excluded from the regression and from the normalisation. No cross-species threshold is established in the literature: report the proportion of censored volumes per species.

## Distortion correction

When field maps exist (`recordings` is `'new'`, `'old'` or `'2_mapdir'`):

1. TOPUP estimates the field and FUGUE corrects the run mean, saving its voxel shift map; FUGUE needs the effective echo spacing (`DwellT`), read from `EffectiveEchoSpacing` in the JSON sidecar, or computed as `TotalReadoutTime / (ReconMatrixPE - 1)`; give it in the launcher (seconds) if neither exists;
2. `modalities/fMRI/sdc_field.py` turns the shift map into a displacement field; its sign is checked by reproducing FUGUE's corrected image, and the field is not used if it does not match or if the distortion stays below 0.1 voxel;
3. the field is applied to **every volume** in the one-shot resampling of step 5 (reference registration, distortion correction and motion of that volume, in a single interpolation).

With one field map per run (`'new'`), each run uses its own field; with a shared field map (`'old'`, `'2_mapdir'`), the reference run's field is used for all runs. The log states for each run whether the field was applied (`SDC ... status ok`). `EDNIX_SDC=0` disables the correction. Without field maps (`'very_old'`), `do_anat_to_func=True` with a non-linear `type_of_transform` (for example `'SyN'`) registers the functional mean to the anatomy with the deformation restricted to the phase-encoding axis (weights 1 along it, 0.1 on the other axes); with the default `'BOLDAffine'` the registration is affine only.

## Denoising (`post_treatment_method`)

`band` is either a high-pass cutoff (`'0.01'`) or a band (`'0.01 0.1'`), in Hz.

| method | what it does |
|---|---|
| `'AFNI'` | design matrix built with `3dDeconvolve` (motion and derivatives, tissue signals, band-pass as sine / cosine regressors, polynomials), projected with `3dTproject` together with censoring and smoothing |
| `'fMRIPrep_style'` | 24-parameter motion model, aCompCor on white matter and CSF (see below), optional global signal, projected with `3dTproject` with the temporal filter; the regressors are saved in `*_confounds_fMRIPrep_style.tsv` |
| `'Grandjean'` | motion and mean white-matter and CSF signals, with the temporal filter; the mean signal is added back unless `normalize` is set |

**Which regressors each method uses.** `extract_WM`, `extract_Vc` and `extract_GS` select the regressors of the `'AFNI'` method. `'Grandjean'` always uses the mean white-matter and ventricle signals (when the masks exist), whatever `extract_WM` and `extract_Vc` say, and `'fMRIPrep_style'` always uses aCompCor of both tissues; `extract_GS` adds the global signal to `'AFNI'` and `'fMRIPrep_style'`. `extract_exterior_CSF` is not available (it must be `False`): no step produces its mask, and the signals are extracted from skull-stripped data, which are zero outside the brain.

**Low-pass.** The `band` is applied in full by the three methods: `'0.01 0.1'` removes the frequencies below 0.01 Hz and above 0.1 Hz (`3dTproject -passband`, or the equivalent regressors for `'AFNI'`); a single value `'0.01'` is a high-pass only. fMRIPrep itself does not low-pass, so a band-pass is a choice of the lab, not of fMRIPrep: for `'fMRIPrep_style'`, give a single value to follow it. `3dTproject` reads the TR in the header of its input while `1dBport` receives it on the command line, so the header TR is checked against the TR of the analysis (`'Grandjean'` and `'fMRIPrep_style'`): a difference stops the run, because the filter would remove the wrong frequencies (`3drefit -TR <seconds> <file>`).

**aCompCor (`acompcor_mode`).** `'fixed'` (default) takes `acompcor_n_fixed` components (5) per tissue (Behzadi et al. 2007). `'adaptive'` takes as many as needed for `acompcor_variance` (50 %) of the variance, at most `acompcor_max` (15): with spatially independent noise the variance is spread over tens of components, so the cap is almost always reached and the mode amounts to a fixed 15 (the log says so). `'parallel'` keeps the components that stand above chance (parallel analysis against phase-randomised surrogates of the tissue time series, 99th percentile), at most `acompcor_max`. In simulations (300 volumes, three shared sources, independent AR(1) noise), the variance criterion needed 30 to 94 components, parallel analysis found exactly the three sources in 68 % of the masks with white noise and 20 % with AR(1) = 0.5 noise (overestimating otherwise, by up to 8), and kept components by mistake in 1 % (white) to 18 % (AR(1) = 0.5) of pure-noise masks: the null is not well calibrated for autocorrelated noise, so it is a diagnostic to read, not a reference. Whatever the mode, the log reports, for each tissue, the number of components found by parallel analysis; compare it with the number used on one subject per species.

Tissue masks for the regressors are eroded by default (`use_erode_WM_func_masks`, `use_erode_V_func_masks`). `normalize` scales the residuals: `'zscore'` over the uncensored volumes, `'psc'` relative to the mean of the data before regression, or `'Skip'`.

**Small tissue masks.** A white-matter or ventricle mask with fewer than 10 voxels at the functional resolution (small brains, eroded masks) is not used: its mean signal would be noise. The log gives the voxel count of each mask, and `*_tissue_regressors.json` (next to the residual) records which tissue regressors each run has. Subjects and species can therefore have different models: report it.

**Several species with one band.** The degrees of freedom left by a band-pass depend on the duration of the acquisition, not on the number of volumes: the 0.01–0.1 Hz band holds about 2 x 0.09 x duration independent samples (54 for 5 min, 108 for 10 min), so a correlation from 5 min of data has a standard error (Fisher z) of about 0.14 against 0.10 for 10 min. With a Grandjean model (6 motion parameters, 2 tissue signals, polynomial), a run of 300 volumes at TR = 1 s keeps 46 degrees of freedom (15 % of the volumes), 300 volumes at TR = 2 s keep 100. Match the duration of the runs across species, or analyse them with the same number of volumes, before comparing connectivity values.

**Degrees of freedom.** For each run, the log gives the uncensored volumes minus the regressors minus the regressors of the temporal filter. A band-pass removes a large part of the spectrum (with 0.01–0.1 Hz, about 80 % of it at TR = 1 s and 65 % at TR = 2 s): with short runs and many regressors, nothing is left. The run then stops with an error; below 20 % of the volumes, a warning is given. Use fewer aCompCor components, a high-pass only, or longer runs.

## Connectivity

- **Matrices** (step 10): regional time series of the selected atlases and levels (`selected_atlases_matrix`, `wanted_level_matrix`), with regions eroded by up to 2 voxels; estimators given by `matrix_fit` (covariance, correlation, partial correlation by default).
- **Seed-based maps** (step 11): seeds given by atlas / level (`selected_atlases_SBA`) and label IDs (`panda_files_SBA`), maps in functional, anatomical and template spaces (`SBAspace`), Fisher z-transformed.
- **Specificity** (step 13): each run is classified as `Specific`, `Unspecific`, `Spurious` or `No` connectivity following Grandjean et al. (2020), from homotopic (specific) and cross-network (non-specific) connectivity.

The atlas tab of the interface picks seeds by clicking regions in 3D and writes them into the launcher.

## Main outputs

- `func/acpc-func/postprocessed_rs/`: denoised time series in AC-PC space at functional resolution.
- `func/templates/<reference>/postprocessed_rs/`: the same in template space, with correlation matrices and seed-based maps.
- `func/raw/matrices/`: motion matrices, FD and censoring files, distortion fields (`*_desc-sdc_1Warp.nii.gz` with a JSON report).
