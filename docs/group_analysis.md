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
