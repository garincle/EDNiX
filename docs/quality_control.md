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

**Metrics**

- Anatomy (step 16, `_16_anat_QC_SNR`): tissue signal-to-noise ratio (mean / standard deviation within a tissue), grey-white contrast-to-noise ratio, coefficient of joint variation, smoothness, per atlas region and hemisphere.
- Functional (step 12, `_12_fMRI_QC`): framewise displacement (Jenkinson, with the brain radius of the subject), fraction of censored volumes, DVARS, outlier counts, grey- and white-matter SNR, contrast, ghost ratio, global correlation, smoothness, and the agreement of the functional and anatomical images after registration (correlation, mutual information, Jaccard index of the masks).
- Connectivity (step 13, `_14_fMRI_QC_matrix`): specificity category of each run (Grandjean et al., 2020).

**Visual checks** in ITK-SNAP (`itk_1`, `itk_2`, `itk_3` steps): masks, registrations and atlases overlaid on the images, with the possibility to correct masks by hand ([skullstripping.md](skullstripping.md)).

**Logs**: every command, its output and return code are in the subject's log. With `EDNIX_STRICT=1` the first failing command stops the run; without it, failures are shown in red.

## Reviewing quality control

The **Quality control** tab of the interface lists the QC folders of a dataset, shows the images as a grid (subject filter, full-screen view with the arrow keys) and the QC tables, where values more than 3.5 robust standard deviations (median absolute deviation) from the median are highlighted.

## Recommendations

- Look at every skull-strip and registration image of a new species or acquisition before running the group analyses.
- Report, per species, the proportion of censored volumes and the temporal SNR: they differ strongly between awake and anaesthetised animals and between field strengths.
- For distortion correction, check in the log that the field was applied (`SDC ... status ok`) and compare the corrected mean image with the anatomy.
