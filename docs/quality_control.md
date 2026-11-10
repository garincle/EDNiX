# Quality control

## What the pipeline produces

**Images** in `<BIDS>/QC/<check>/`, one per subject, session or run, for example:

| folder | check |
|---|---|
| `skullstrip_step1` | first anatomical skull-strip |
| `align_rigid_to_template` | AC-PC / rigid alignment to the template |
| `fMRI_runMean_in_REF` | each run registered to the reference run |
| `check_mask_fMRI` | functional brain mask |
| `mask_to_fMRI_orig`, `meanIMG_in_anat` | functional-anatomical registration |
| `fMRI_MEAN_in_REF_FUNC_ACPC` | functional mean in AC-PC space |
| `pet_in_acpc`, `pet_in_template` | PET session mean in AC-PC / template space, brain mask of that space as contour |
| `ct_in_acpc`, `ct_in_template` | CT in AC-PC / template space |

**Metrics**

- Anatomy (step 16, `_16_anat_QC_SNR`): tissue signal-to-noise ratio (mean / standard deviation within a tissue), grey-white contrast-to-noise ratio, coefficient of joint variation, smoothness, per atlas region and hemisphere.
- Functional (step 12, `_12_fMRI_QC`): framewise displacement (Jenkinson, with the brain radius of the subject), fraction of censored volumes, DVARS, outlier counts, grey- and white-matter SNR, contrast, ghost ratio, global correlation, smoothness, and the agreement of the functional and anatomical images after registration (correlation, mutual information, Jaccard index of the masks).
- Connectivity (step 13, `_14_fMRI_QC_matrix`): specificity category of each run (Grandjean et al., 2020).
- PET: motion of every frame (`pet/orig/preprocessing/forMC/*_motion.tsv`, `.png`: translations, rotations, mean displacement of a sphere of the radius of the head), time-activity curves and kinetic fits (`pet/orig/postprocessing/ROI_analysis/pictures`), and the subject values used by the SUV with their source (`pet/Prepro_Settings.json`).

**Visual checks** in ITK-SNAP (`itk_1`, `itk_2`, `itk_3` steps): masks, registrations and atlases overlaid on the images, with the possibility to correct masks by hand ([skullstripping.md](skullstripping.md)).

**Logs**: every command, its output and return code are in the subject's log. With `EDNIX_STRICT=1` the first failing command stops the run; without it, failures are shown in red.

**Surface QC** (step 16, `Tools/QC_report.surface_qc`): white, pial and midthickness contours on three planes over the brain, thickness distribution in real millimetres (small brains: FreeSurfer millimetres divided by the header scale factor) and the distance of the white / pial surfaces to the aseg boundaries. `anat/QC_anat/<ID>_surface_QC.png` and `.json`, copied to `<BIDS>/QC/surfaces/`.

**Adaptive registration** (`SyNadaptive`): candidate kept and its overlap metrics, `<prefix>_adaptive_QC.json`.

## Reports at three levels

`Tools/QC_report.build(bids_dir)` (run by `run_all_ednix.py`, or `python -m Tools.QC_report <BIDS> [<BIDS> ...]`) gathers everything above, without recomputing it, in:

| file | content |
|---|---|
| `QC/subjects/sub-<ID>_ses-<S>.html` | one subject / session: processing status of every step, log errors and warnings, anatomical, functional, specificity, PET motion and surface metrics, all its QC images |
| `QC/QC_all_subjects.html` | one row per subject / session with the key metrics, values more than 3.5 robust SD from the median highlighted, one thumbnail grid per QC check |
| `QC/QC_dataset.html` | the dataset: modalities, subjects, acquisition parameters, TR flags, completion of each step, distributions |
| `QC/EDNiX_BIDS_summary.xlsx` | sheets `overview`, `subjects`, `acquisition` (TR from the sidecar and the header, `TR_check` when they differ by more than 5 % or fall outside 0.1–10 s, volumes, duration, voxel size), `status`, `QC_summary`, `anat_QC`, `surface_QC`, `func_QC`, `FC_specificity`, `logs` |

`QC_report.study_summary([bids_1, bids_2, ...], out.xlsx)` stacks the overview and acquisition sheets of several datasets, with a `dataset` column.

## Reviewing quality control

The **Quality control** tab of the interface lists the QC folders of a dataset, shows the images as a grid (subject filter, full-screen view with the arrow keys) and the QC tables, where values more than 3.5 robust standard deviations (median absolute deviation) from the median are highlighted.

## Recommendations

- Look at every skull-strip and registration image of a new species or acquisition before running the group analyses.
- Report, per species, the proportion of censored volumes and the temporal SNR: they differ strongly between awake and anaesthetised animals and between field strengths.
- For distortion correction, check in the log that the field was applied (`SDC ... status ok`) and compare the corrected mean image with the anatomy.
