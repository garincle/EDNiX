# Group and comparative statistics

The group analyses are in `Statistics/`. They read the outputs of the functional and anatomical pipelines of several subjects.

## Voxel-wise functional statistics

| module | analysis |
|---|---|
| `Group_fMRI/_Group_anal_3dLMEr_SBA.py` | linear mixed-effects models (AFNI `3dLMEr`) on seed-based maps: repeated sessions, groups, covariates, general linear tests |
| `Group_fMRI/_Group_anal_3dTtest.py` | one-sample tests (AFNI `3dttest++`) on seed-based maps |
| `Group_fMRI/_Group_anal_3dLMEr_Mirror.py` | left versus right hemisphere comparisons of mirrored maps |
| `Group_fMRI/_Group_anal__func_DicLearn.py` | dictionary learning (nilearn) of group functional networks; `Select_cpmnt_number.py` and `Select_parameters_DicL.py` choose the number of components and the parameters from split-half stability |

**Group mask.** `_mask_utils.build_group_mask` builds the mask from the mean functional images of the subjects, restricted to grey matter or with white matter excluded (`method_mask_func`), and keeps only voxels covered by every subject (or uses a provided mask).

**Cluster correction.** Maps are thresholded voxel-wise at `csim_pthr` (0.001 by default) and clusters smaller than the size given by AFNI `3dClustSim` for the family-wise `alpha` are removed. `3dClustSim` uses the spatial autocorrelation (ACF) estimated on the residuals of the group model (`3dFWHMx -acf`); for `3dttest++` with more than 14 subjects, its built-in permutation option (`-Clustsim`) is used. The voxel threshold, cluster size and ACF parameters are printed in the log and shown in the figure titles. `3dLMEr` sub-bricks are found by the labels of the general linear tests.

**Mirror analyses** flip one hemisphere onto the other around `midline_x`: use a symmetric template, otherwise homologous regions do not overlap.

## PET group analyses

`Group_PET/` analyses the images of the PET pipeline (`docs/pet_pipeline.md`) of several subjects. Examples: `Exemples/Study_EDNiX/analysis/Launcher_PET_group_CERMEP.py` (human, static FDG) and `Launcher_PET_group_Tremblay.py` (macaque, dynamic FDG), both on the brain ratio (`desc='bratio'`, `normalization='none'`).

**Static and dynamic runs.** Step 7 names its images `desc-<x>` for a static run and `desc-<x>-mean` (mean of the frames of `SUV_window`) for a dynamic run. `desc='SUV'` and `desc='SUV-mean'` are interchangeable in the voxel-wise and regional analyses: each session uses the one it has, and the choice is printed. When the file names carry no `ses-` entity, the session is the name of the `ses-<label>` folder.

**Subject tables.** Covariates and groups are read from `participants.tsv` (ID column `participant_id`, or `subject`, `sub`, `id`...) and `sub-<ID>/sub-<ID>_sessions.tsv`. They are optional: with `covariates=None` and no group column to read (one-sample, paired on sessions, or two-sample with `group_assignment={'A': [IDs], 'B': [IDs]}`) the tables are not read.

| module | analysis |
|---|---|
| `Group_PET/pet_group.py` | voxel-wise GLM of the images exported to a template space (step 10): `run_group_analysis` |
| `Group_PET/pet_roi_group.py` | the same models on the regional values computed in the PET space (step 7): `run_roi_analysis` |
| `Group_PET/pet_surface.py` | projection of the maps on the surfaces of the species, with the function of the SBA figures (`Tools/niftotoWBsurface.py`) |
| `Group_PET/_glm.py` | design matrices, OLS, Freedman-Lane permutations, FWE and FDR |

**Designs.** `'one-sample'`, `'two-sample'` (`group_col`, `groups=[A, B]`, contrast A - B), `'regression'` (effect of `test_covariate`), `'paired'` (difference between two sessions of the same subjects, `paired_levels`), with `covariates` read from `participants.tsv` and `sub-<ID>_sessions.tsv` (numeric: centred; text: categorical). Runs of a session are averaged; sessions of a subject too (`aggregate='subject'`), because the observations of the model must be independent.

**Intensity normalisation** (`normalization`), computed on each image before smoothing:

- `'global'`: proportional scaling to the mean of the brain (Fox et al. 1988; SPM);
- a reference region (SUVR), with the syntax of `reference_signal` of the PET pipeline: `{'type': 'atlas', 'atlas': 'EDNIxCSCLR', 'level': 0, 'labels': 'default', 'name': 'Cer'}` (cerebellum), label numbers, a name pattern, or `{'type': 'mask', 'file': <mask in the template space>}` (e.g. pons, Minoshima et al. 1995);
- `'none'`: SUV.

Global normalisation turns a large or strong focal decrease into apparent increases elsewhere, because the decrease lowers the global mean (Borghammer et al. 2008; Yakushev et al. 2009; in one simulated study of the tests' phantom, a 15 % focal decrease produced a significant spurious increase in the cerebellum with `'global'`, none with the cerebellum reference). The normalising value of every image is written in `inputs.tsv`, and for two groups it is compared between them (`normalising_value_test` in the summary, warning when p < 0.05): choose a reference region not affected by the condition studied, for each species and condition (Dukart et al. 2010; López-González et al. 2020).

**Mask and smoothing.** Analysis mask = brain mask of the template (and `tissue_mask`, e.g. the grey-matter mask) x voxels covered by every image x relative threshold (`rel_threshold=0.8`: voxels below 0.8 x the brain mean in any image are excluded, the SPM default for PET): it removes the ventricles, the CSF and the background. Smoothing (`fwhm`, mm; twice the voxel size by default) is restricted to this mask (normalised convolution, Knutsson & Westin 1993), so that the signal outside the brain does not spread into the cortex. Scale `fwhm` with the brain: 8-12 mm in humans for FDG (Signorini et al. 1999), a few mm in macaques. Images exported at the PET resolution of different scanners are resampled to the grid of the first one.

**Inference.** Voxel-wise two-sided threshold `voxel_p` (0.001), cluster-level family-wise error `alpha`:

- `cluster_method='permutation'` (default): Freedman-Lane permutations (sign flips for one-sample and paired designs), maximum cluster size and maximum |t| over the mask (Nichols & Holmes 2002; Winkler et al. 2014; Eklund et al. 2016 for the validity of permutation against parametric cluster inference). Exhaustive when there are fewer than `n_perm` permutations; the smallest attainable p is 1 / number of permutations (with 5 subjects per group, 252 distinct relabellings). Positive and negative clusters are formed separately (face connectivity). Also written: voxel-level FWE p (`_logp-FWE-voxel`). On null smooth data, the cluster FWE rate was 0.07 and the voxel FWE rate 0.01 for alpha = 0.05 (100 simulations, 6 vs 6, 200 permutations each; `tests/unit/test_pet_group.py`);
- `cluster_method='clustsim'`: AFNI `3dClustSim` with the ACF of the residuals of the model (bi-sided table), as the fMRI analyses.

**Outputs** (`<BIDS>/Results/Grp_PET_<design>_desc-<desc>_norm-<norm>[_<suffix>]/`): `mask_group`, `inputs_4D` (normalised, smoothed images) and `inputs.tsv`, group means and SD, `<contrast>_cope / _tstat / _zstat / _zstat_thresholded / _cope_thresholded`, percent difference (two groups), `_clusters.tsv` (size, peak z and coordinates, region of the peak and main regions of the cluster, cluster FWE p), `_summary.json`, PNG slices. `pet_surface.render_pet_surfaces` projects the z maps, the percent difference and the group mean on the surfaces of the reference template of the atlas library (ribbon-constrained projection by default: only the voxels between the white and pial surfaces enter a vertex); use it with maps in that template space.

**Regional analysis.** `run_roi_analysis(bids_dir, atlas, level, desc, design, ...)` reads the regional tables of the PET space (`ROI_analysis/<base>_seg-<atlas>-lvl<k>_desc-<desc>_ROIs.tsv`), drops the regions with fewer than `min_voxels` PET voxels (10) in at least one session (rule of the earlier EDNiX scripts; partial volume and registration error dominate smaller regions), normalises to a reference region (voxel-weighted mean of its labels) or to the voxel-weighted mean of all regions (`'global'`), and fits the same models per region with Benjamini-Hochberg FDR. Tables `values_raw` / `values_normalised` (TSV and Excel, one row per session), `regions_voxel_counts.tsv`, `<contrast>_regions.tsv` (effect, 95 % CI, t, p, FDR, group means), and a figure of the effects.

## Longitudinal models

`Group_fMRI/EDNiX_longitudinal.py` models a regional measure against age with a random intercept per subject:

- `fit_gamm`: cubic B-spline of age (4 degrees of freedom), maximum-likelihood fit; the age effect is tested with a likelihood-ratio test against the model without age. The fitted curve and its confidence band come from the fitted spline basis and the full covariance of the coefficients.
- `fit_lme`: linear age effect, with the marginal R² of Nakagawa and Schielzeth.

## Phylogenetic comparative analyses

`Evo/PGLS.py` relates regional measures across species while accounting for their phylogeny (a Newick tree):

- phylogenetic generalised least squares with Pagel's λ estimated by maximum likelihood;
- likelihood-ratio tests of λ against 0 (no phylogenetic signal) and 1 (Brownian motion): with about ten species, λ is weakly identified, report these tests with it;
- confidence band of the regression line and prediction interval for a new species;
- all region pairs with Benjamini-Hochberg correction (`run_pgls_all_region_pairs`).

Species means ignore within-species variability; with individual data, a phylogenetic mixed model is more powerful.

## References (PET group analyses)

- Anderson MJ, Robinson J (2001) Permutation tests for linear models. Aust N Z J Stat 43:75-88.
- Borghammer P, Jonsdottir KY, Cumming P, et al. (2008) Normalization in PET group comparison studies: the importance of a valid reference region. NeuroImage 40:529-540.
- Dukart J, Mueller K, Horstmann A, et al. (2010) Differential effects of global and cerebellar normalization on detection and differentiation of dementia in FDG-PET studies. NeuroImage 49:1490-1495.
- Eklund A, Nichols TE, Knutsson H (2016) Cluster failure: why fMRI inferences for spatial extent have inflated false-positive rates. PNAS 113:7900-7905.
- Fox PT, Mintun MA, Reiman EM, Raichle ME (1988) Enhanced detection of focal brain responses using intersubject averaging and change-distribution analysis of subtracted PET images. J Cereb Blood Flow Metab 8:642-653.
- Freedman D, Lane D (1983) A nonstochastic interpretation of reported significance levels. J Bus Econ Stat 1:292-298.
- Knutsson H, Westin CF (1993) Normalized and differential convolution. Proc IEEE CVPR 515-523.
- López-González FJ, Silva-Rodríguez J, Paredes-Pacheco J, et al. (2020) Intensity normalization methods in brain FDG-PET quantification. NeuroImage 222:117229.
- Minoshima S, Frey KA, Foster NL, Kuhl DE (1995) Preserved pontine glucose metabolism in Alzheimer disease: a reference region for functional brain image (PET) analysis. J Comput Assist Tomogr 19:541-547.
- Nichols TE, Holmes AP (2002) Nonparametric permutation tests for functional neuroimaging: a primer with examples. Hum Brain Mapp 15:1-25.
- Phipson B, Smyth GK (2010) Permutation p-values should never be zero. Stat Appl Genet Mol Biol 9:39.
- Signorini M, Paulesu E, Friston K, et al. (1999) Rapid assessment of regional cerebral metabolic abnormalities in single subjects with quantitative and nonquantitative [18F]FDG PET: a clinical validation of statistical parametric mapping. NeuroImage 9:63-80.
- Winkler AM, Ridgway GR, Webster MA, Smith SM, Nichols TE (2014) Permutation inference for the general linear model. NeuroImage 92:381-397.
- Yakushev I, Hammers A, Fellgiebel A, et al. (2009) SPM-based count normalization provides excellent discrimination of mild Alzheimer's disease and amnestic mild cognitive impairment from healthy aging. NeuroImage 44:43-50.
