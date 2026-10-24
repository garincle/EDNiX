# Changelog

## Unit tests green again, October 2026

The `tests` workflow had been failing on master since 6 October (28 tests on Python 3.11).
- `extract_exterior_CSF` had been removed from `signal_regression` / `preprocess_fmri`: the AFNI emulator of the tests still passed it (21 `test_post_ttt` failures), `Launcher_func_Dog.py` still passed it (removed), and its test now checks that the option no longer exists (it used to be refused with an error).
- `preprocess_pet`: the `coreg_transform` entry of the docstring was written as code, so the parameter reference and the interface did not see it.
- Interface QC test updated to the current `qc_scan` paths (relative to the BIDS folder).
- PET tests: `ants.config` exists only from antspyx 0.6; the hand-built options of the end-to-end test lacked `coreg_transform`.
- **[results]** CT -> T1w registration (scenario 1v2b): the quality score of `rigid_to` assumes a bright brain (PET); for a CT a correct mutual-information registration scored ~0.2 and the `meansquares` retry, not valid between a CT and a T1w, was kept (11 mm error on the test phantom; probably the CT / anatomy misalignment seen on CERMEP). The CT call now keeps the mutual-information registration. Stopgap: the registration rework (PR #15) replaces this logic.

## PET group analyses, October 2026

New `Statistics/Group_PET` (documentation: [docs/group_analysis.md](docs/group_analysis.md#pet-group-analyses); example `Exemples/Study_EDNiX/analysis/Launcher_PET_group_CERMEP.py`; tests `tests/unit/test_pet_group.py`):

- `pet_group.run_group_analysis`: voxel-wise one-sample, two-sample, regression and paired designs with covariates on the PET images exported to a template space; normalisation to the brain mean, to a reference region (same syntax as `reference_signal`) or none, with the normalising value tested between groups; analysis mask with coverage and the relative threshold of SPM (excludes CSF and background); smoothing restricted to the mask; cluster FWE by Freedman-Lane permutation (default) or AFNI 3dClustSim on the residual ACF; cluster table with atlas regions.
- `pet_roi_group.run_roi_analysis`: the earlier ROI scripts (drop regions with fewer than 10 PET voxels in any subject, division by a reference region, tables as Excel) on the regional tables of the PET space, with GLM and FDR per region.
- `pet_surface.render_pet_surfaces`: maps projected on the surfaces of the species with the SBA surface function, ribbon-constrained projection by default.
## Cross-species regional profiles, October 2026

- `Plotting/Plot_interspecies_regional_profiles.py`: gathers the ROI tables of the PET pipeline (`ROI_analysis/*_seg-<atlas>-lvl<k>_desc-<desc>_ROIs.tsv`) for several species and compares them on the regions of the EDNiX atlas (shared label table, so regions are matched by name). One subject = one sample (runs averaged), hemispheres merged by valid-voxel weighting or kept separate, optional ratio to reference regions or within-subject z-score. Regions kept when present in every species with >= 10 voxels per label in every subject (rule of the former `3_correct_atlas_for_small_V.py`); optional size threshold in FWHM units against partial volume. Per region, each species vs the reference species: permutation test (or Welch / Mann-Whitney), Hedges g with 95% CI, Benjamini-Hochberg FDR. Profile similarity: Pearson / Spearman r with subject-bootstrap CI, split-half reliabilities and noise ceiling. Figures: polar fingerprint + polar effect sizes, grouped bars with every subject, region-mean scatter; each with a .txt legend, plus an Excel workbook and a methods text.
- Example `Exemples/Study_EDNiX/analysis/PET_interspecies_profiles.py` (CERMEP human vs BIDS_Tremblay macaque, FDG). Tests `tests/unit/test_interspecies_profiles.py`.

## PET and CT pipelines, October 2026

Import of the PET engineer's version (`modalities/pet`, `modalities/ct`), reviewed and corrected: the detail of every change is in [docs/PET_CT_review_journal.md](docs/PET_CT_review_journal.md), the methods in [docs/pet_pipeline.md](docs/pet_pipeline.md).

- `modalities.pet._0_Pipeline_launcher.preprocess_pet` and `modalities.ct._0_Pipeline_launcher.preprocess_ct`: processing in the PET space, registration scenarios 1v1, 1v2a, 1v2b, 2v1, 2v2, 3v1, 3v2 (`Illustration/pet4.svg`), output folders of `Illustration/pet2b.svg`, bridge with the anatomy through `backtonative.get` as the functional pipeline, registration records.
- Motion correction within (two passes) and between runs; SUV (units, decay reference), SUVglu, ratio to the brain mean, time blocks; regional values and TACs; SRTM / Logan with the Turku tools (decay correction checked, reference TAC over the union of the voxels); partial-volume correction (PETPVC, scanner PSF); PET during fMRI runs; export to the AC-PC and template spaces.
- `Tools/getpath.py`: `pet()` and new `ct()` return the folders of `pet2b.svg` as dictionaries (only the PET / CT code used `pet()`).
- `Tools/roiStats.py` (numpy), `atlases/getlabels.py`, `FS_freeview.check_coreg` (from the engineer's suggestions; his `get_orientation.py` not taken: it reintroduced the left-right mirrors corrected in September; his `extract_filename.py` helpers not taken: unused).
- Interface: CT tab, PET / CT steps and options, PET / CT logs and QC. Example launchers `Launcher_pet` (BIDS_Tremblay, CERMEP) and `Launcher_ct` (CERMEP). Tests `tests/unit/test_pet.py`.
- Deep verification (22/10/2026, journal section 6) **[results]**: PET / CT of the same scanner space follow the world transform of the reorientation of the anatomy / CT (`3drefit -orient` re-signs the origins) and its head position (`same_space_window`); motion correction with a head mask for the metric (relative motion within 0.2-0.4 mm instead of 1.5-2 mm on a phantom); `SyNRA` by default for the template registrations; frame times counted from the injection; partial-volume correction with complete tissue classes and fractions; TACs of the kinetic models on the non-skull-stripped image; `DV_k2ref`; anatomy of another session and longitudinal reference session; `preprocess_data` refused.

## Code review, September 2026

48 corrections from a full code review. The PET module is unchanged (under development).
Entries marked **[results]** change pipeline outputs: re-run the affected steps before reusing numbers.

### Execution and robustness
- `Tools/run_cmd.py`: every command's return code is checked and failures are logged; `EDNIX_STRICT=1` makes the first failure stop the pipeline. Arguments are split with `shlex` (quotes, repeated spaces). `do()` reports the return code; `msg()` has a default colour.
- `nilearn.input_data` replaced by `nilearn.maskers` (`_10`, `_11`, `_12`): the functional pipeline imports with current nilearn.
- Broken imports fixed in `Load_BIDS_data_for_analysis`, `WB_matrices`, `FS_modif`, `anat/diary.py`.
- `Tools/Clean_BIDS.py` no longer deletes files when imported (code under `if __name__ == '__main__'`).
- `preprocess_fmri` warns about unknown keyword arguments (they were silently ignored); the confounds-length warning no longer crashes.
- Atlas label helpers moved to `Tools/atlas_labels.py`: the core no longer imports study scripts from `Exemples/`. `atlas4func` raises a clear error for `'all'` with a non-EDNiX reference.
- `Skullstrip_func`: commands go through `run_cmd` (failures visible); a bare `-dilate_input` that made `3dmask_tool` fail is removed.
- `Load_subject_with_BIDS`: excluded sessions (`list_to_remove`) no longer enter the study template or the longitudinal reference; unmatched filter entries are reported. `_3_make_template` reads the session from the path. **[results]**
- `_12_fMRI_QC`: invalid escape sequence fixed.

### Geometry and registration
- `Tools/check_nii.py`: grids are compared on dimensions, spacing, origin and direction (images with the same size but a shifted origin were combined voxel-wise); `comphd` requires origins within one voxel; `keep_header` refuses mismatched grids. **[results]**
- `Tools/get_orientation.py`: orientation mapping computed geometrically; 5 animal/scanner position combinations encoded a left-right mirror. **[results]**
- `transfoparams` / `load_transfo_parameters`: one independent block per step with the names the pipeline looks up (`SS1`, `SS2`, `SS3`, `between` were never found and fell back to a rigid registration). **[results]**
- Motion correction: stale transform files are removed before a run; motion matrices are re-estimated in the frame of the reoriented data when the orientation fix changed it (they were applied about the wrong axes). **[results]**
- `transfo_for_func`: one-shot resampling splits the series in memory (volume order no longer depends on file-name sorting).
- `do_anat_to_func=True` runs the non-linear functional-to-anatomical registration again (the block was commented out); `_8` finds the warp file. **[results]**
- Anatomical images resampled to the functional grid and study-template inputs use cubic interpolation (was nearest neighbour). **[results]**
- Jenkinson FD uses the brain radius of each subject (was 80 mm for every species); the affine matrix file stores `[R | t]`. **[results]**
- `backtonative`: correct transform name in the longitudinal linear branch.
- Field maps are read from their real location (`list_map`). **[results]**
- `templatefeat`: deterministic choice of the species atlas folder.
- `acpcalign` (3dAllineate option): linear interpolation during the optimisation. **[results]**

### Denoising and connectivity
- FD censoring threshold scaled with brain size (0.5 mm for a 1200 cm³ brain) instead of twice the voxel size. **[results]**
- `_7_post_TTT`: percent signal change uses the pre-regression mean; z-scoring uses uncensored volumes only. **[results]**
- `_10_Correl_matrix`: label-wise erosion really erodes up to 2 voxels; the resampling fallback works and checks labels. **[results]**
- `_14_fMRI_QC_matrix`: FC-specificity categories follow Grandjean et al. (2020). **[results]**
- `_6_Melodic`: selected components passed to `fsl_regfilt` 1-based (the previous component was removed); an empty selection copies the data. **[results]**
- `_adaptive_n4_correction`: B-spline spacing scaled with brain size (small brains were barely corrected). **[results]**

### Group statistics
- `_Group_anal_3dLMEr_SBA`, `_Group_anal_3dLMEr_Mirror`, `_Group_anal_3dTtest`: the cluster-size threshold is applied (it was only printed); `3dClustSim` uses the ACF of the model residuals; voxel threshold 0.001 by default; 3dLMEr sub-bricks found by label. **[results]**
- `_mask_utils`: the group mask keeps only voxels covered by every subject. **[results]**
- `EDNiX_longitudinal`: exact ML random-intercept model, likelihood-ratio test of the age effect (the minimum coefficient p-value was used), prediction basis with the fitted knots, full covariance for the band, correct marginal R². **[results]**
- `PGLS`: GLS confidence band and prediction interval, likelihood-ratio tests of lambda against 0 and 1, robust Benjamini-Hochberg. **[results]**
- `Select_cpmnt_number`, `Select_parameters_DicL`: component matching insensitive to sign. **[results]**
- `_Group_anal__func_DicLearn`: first image used as reference grid.

### Anatomy
- `FS_2_WB`, `WB_sphere_reg`: thickness, sulcal depth and curvature rescaled for brains enlarged before FreeSurfer. **[results]**
- `_16_anat_QC_SNR`: tissue SNR, CNR and CJV (the previous noise estimate was not noise). **[results]**
- `WB_vol_atlas`: hemisphere label maps combined without adding label values. **[results]**
- `myelin`: T2w resampled with continuous interpolation for the T1w/T2w ratio. **[results]**
- `_1_correct_orient`: N4/N3 bias correction restricted to a head mask. **[results]**
- `Skullstrip_method`: subject-specific branches removed (they never ran); `correctaseg` dilation radius in voxels.

### Study figures (`Exemples/Study_EDNiX/analysis`)
- Figure 2: cross-modality correlations on log values with FDR. **[results]**
- Figure 3: paired test for association vs primary cortex; summary text computed from the results (was hard-coded). **[results]**
- Supplementary figures 4–6: partial-correlation degrees of freedom; fixed-effect meta-analysis across datasets. **[results]**

### Distortion correction of the time series
- The FUGUE correction was applied to the run mean only. FUGUE now also saves its voxel shift map; `modalities/fMRI/sdc_field.py` turns it into an ANTs displacement field, whose sign is checked by reproducing FUGUE's own unwarped image (the field is not used if it does not match, or if the distortion is below 0.1 voxel). The field is applied to every volume in the one-shot resampling (`run_inRef`, then SDC, then motion) and to the run-to-reference registration. With one field map per run (`new`), each run uses only its own field (a run whose field was rejected is not corrected with another one); with a shared field map (`old`, `2_mapdir`), the reference run's field is used for every run. Fields of a previous execution are deleted before they are rebuilt. `EDNIX_SDC=0` disables it. **[results]**
- With exactly `[affine.mat, field]`, ANTsPy inverts the affine when `whichtoinvert` is not given: every transform list that can contain a field now passes `whichtoinvert` explicitly.
- Without field maps, `do_anat_to_func=True` (non-linear registration restricted to the phase-encoding axis) remains the available correction.

### Timing parameters of the functional pipeline
- The TR is resolved once, as a float in seconds (`resolve_TR`). When it was given in the launcher (as text, as the documentation asked), it stayed text: the fMRIPrep-style DCT regressors computed `2 * n * TR / period`, Python repeated the text and the division raised a `TypeError`, caught by the denoising step, which then produced no residual. When `RepetitionTime` is missing, the estimate from `SliceTiming` was multiplied by 1000 (a value in milliseconds reported as seconds) and counted simultaneous multiband slices as separate steps. **[results]**
- MELODIC received the launcher value of `TR` (`'Auto'` by default) instead of the resolved TR. **[results]**
- The dwell time given to FUGUE (`--dwell`, the effective echo spacing) is read from `EffectiveEchoSpacing`, or computed as `TotalReadoutTime / (ReconMatrixPE - 1)` (`resolve_dwell`). It was read from `DwellTime`, which in BIDS is the receiver dwell time per readout sample (microseconds), making the correction about a hundred times too small; the second fallback divided a non-BIDS `TotalReadOutTimeEPI` by the number of slices. `DwellTime` is still used when its value is a plausible echo spacing (non-standard sidecars), with a warning. **[results]**
- The non-linear functional-to-anatomical registration (`do_anat_to_func`) received the launcher value of `SED` (`'Auto'` by default), so its deformation was never restricted; and `SED` is the slice-encoding direction, while EPI distortions happen along the phase-encoding direction. The deformation is now restricted to the physical axis of the phase-encoding direction (from `PhaseEncodingDirection` or `correction_direction`). **[results]** (only with a non-linear `type_of_transform`)
- Docstrings of `preprocess_fmri` corrected for `correction_direction`, `DwellT`, `SED`, `TR`, `TRT`, `doWARPonfunc`, `registration_type`, `intra_thresh` and `delta_thresh` (the parameter reference is generated from them).

### Denoising (step 7, `_7_post_TTT`)
- New version: high-pass-only or band-pass filtering (`band` = `'hp'` or `'hp lp'`), fMRIPrep-style 24-parameter motion model with fixed or adaptive aCompCor, regressor table saved for each run.
- The band is validated (one or two frequencies, low-pass above high-pass; a low-pass at or above Nyquist is ignored).
- A degrees-of-freedom budget is computed for every run (uncensored volumes minus regressors minus filter regressors): the run stops if none are left, with a warning below 20 %. With a 0.01–0.1 Hz band-pass, the filter alone removes 60–80 % of the degrees of freedom; the fMRIPrep-style design exceeded the number of volumes for 300 volumes at TR = 1 s. **[results]**
- fMRIPrep style: no cosine drift regressors (the temporal filter already removes those frequencies); aCompCor on voxel time series high-passed with the cosine basis, as in fMRIPrep; global signal = mean in the brain mask (it was never added: its mask did not exist). **[results]**
- Global-signal regressor: `maskDilat.nii.gz` does not exist; the brain mask of step 5 (`mask_ref.nii.gz`) is used. With the AFNI method, requesting the global signal made every run fail. **[results]**
- Exterior CSF: its mask (`exterior_ligne.nii.gz`) is not produced by the pipeline; a missing mask is reported instead of being replaced by the ventricle mask (which duplicated the ventricle regressor).
- Grandjean: `normalize` was ignored; z-score / PSC are applied to the residual (without the mean restoration, which would make PSC wrong).
- Tissue masks with fewer than 10 voxels are skipped with a warning in the three methods (`'Grandjean'` and aCompCor had no guard; `'AFNI'` had one for the ventricles only): small brains and eroded masks can leave a handful of voxels, and the mean of a handful of voxels is not a tissue signal. The voxel counts and the decision are written to `*_tissue_regressors.json` for each run, since the model then differs between subjects and species. **[results]** for small brains.
- Low-pass: the three methods apply the whole `band` (high and low edge); the tests now check it on the spectrum of the residual (the emulated `3dTproject` and `1dBport` apply the band, and the tests fail when the low-pass is removed), whereas they only compared the command line before. `3dTproject` takes the TR from the header of its input for `-passband` / `-stopband` while `1dBport` receives it on the command line: a header TR different from the TR of the analysis is now an error for `'Grandjean'` and `'fMRIPrep_style'` (the filter would remove the wrong frequencies). **[results]** if your headers had a wrong TR.
- aCompCor: `acompcor_mode`, `acompcor_n_fixed`, `acompcor_variance`, `acompcor_max` are now parameters of `preprocess_fmri` (they could not be set). The default mode is `'fixed'` (5 components per tissue, Behzadi et al. 2007) instead of `'adaptive'`, whose 50 % variance target is out of reach with spatially independent noise (30 to 94 components in simulations): the cap of 15 always decided, so it was a fixed 15 (the log now says so). New mode `'parallel'` (components above chance, phase-randomised null) and, in every mode, a log line with the number found by parallel analysis. Only `post_treatment_method='fMRIPrep_style'` is concerned (no example launcher uses it). **[results]**
- `extract_exterior_CSF=True` stops the step with an error: no step produces `exterior_ligne.nii.gz` and the signals are extracted from skull-stripped data (zero outside the brain); it used to skip the regressor with a warning. The behaviour of `extract_WM` / `extract_Vc` per method is documented (`'Grandjean'` and `'fMRIPrep_style'` ignore them).
- Custom confounds table (`*_confounds_correct.tsv`): fMRIPrep tables have NaN in the first row of the derivative columns, which `3dTproject` and `3dDeconvolve` refuse; NaN are replaced by 0 with a warning, text columns are ignored, and the number of rows is checked for the AFNI method as well. **[results]**
- Errors: every command goes through `run_cmd` (3dDeconvolve and 1dcat bypassed the log); a failing run is logged with its traceback and listed at the end of the step, and stops the pipeline with `EDNIX_STRICT=1` (the previous fallback wrote an unused `*_3dDeconvolve_failed` copy, so failures only appeared later as a missing residual). Regressor files of a previous execution are not reused.

### Tests
- `tests/unit` (pytest): every module compiles and imports; orientation, grids, registration parameters, normalisation, erosion, one-shot resampling, SDC, longitudinal model (type-I error), cluster extent, component matching, PGLS, example launchers (no unknown parameter), interface. Run with `pytest -n auto`.
- `.github/workflows/tests.yml` runs them on every push.
- `tests/e2e/run_smoke.py`: every example launcher whose data is on the machine, on one subject with shortened functional runs, with `EDNIX_STRICT=1`; report in `smoke_report.md`. `.github/workflows/smoke.yml` runs it on a self-hosted runner of the lab server, on demand.
- The marmoset functional launchers passed `registration_fast`, which `preprocess_fmri` does not know (it was ignored): removed.

### New
- `GUI/`: desktop interface: 3D atlases (volume and surface layers, side by side, selections linked across layers by region name or label ID), launcher editor for the anatomical, functional and PET pipelines with parameter descriptions, background runs with the current step, quality-control viewer (`<BIDS>/QC` images and tables with outliers highlighted), remote use through an SSH tunnel (`GUI/open_remote.sh`). See `GUI/README.md`.
- `requirements.txt`: Python dependencies.
- Documentation: new `README.md` (overview, installation, quick start, interface) and `docs/` (installation, launchers, BIDS, anatomical and functional pipelines, skull-stripping, atlases, group statistics, quality control, tests); parameter reference generated from the code by `docs/generate_parameter_reference.py`. Corrected on the way: the second manual anatomical mask is `<ID>_final_mask_2.nii.gz`, and `MAIN_PATH` must end with `/`.
