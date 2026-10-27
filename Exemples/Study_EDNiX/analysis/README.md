# EDNiX paper -- analysis scripts

One script per paper figure at this level; everything else is filed by category.
All paper figures are written to `/scratch2/EDNiX/results/EDNiX_Figures/Paper/<Figure>/`
(`paper_data.PAPER_DIR`), as PDF (fonts embedded) + PNG 600 dpi, Nature widths (89 / 120 / 183 mm).

## Paper figures

| Figure | Script | Input |
|---|---|---|
| Fig. 1 | `Figure_1_atlases.py` | atlases |
| Fig. 2 | `Figure_2_morphometry.py` | morphometry, atlas level 2, bilateral, log10 |
| Fig. 3 | `Figure_3_longitudinal.py` | BIDS_Zhu_Garin (macaque), BIDS_Park (human) |
| Fig. 4 | `Figure_4_seed_maps.py` | panel 5a of the Fig5 pipeline (seed maps) |
| Fig. 5 | `Figure_5_fingerprints.py` | panels of the Fig5 pipeline (F15b, 5b, F16a) |
| Supp. 4 | `Supplementary/Supp_Fig_04_anat_QC.py` | anat QC |
| Supp. 5 | `Supplementary/Supp_Fig_05_func_QC.py` | func QC |
| Supp. 6 | `Supplementary/Supp_Fig_06_human_assoc_primary.py` | BIDS_Park |
| Supp. 7 | `Supplementary/Supp_Fig_07_selection.py` | F1e, F1c, F5 |
| Supp. 8 | `Supplementary/Supp_Fig_08_scores_arousal.py` | F3, 6c |
| Supp. 9 | `Supplementary/Supp_Fig_09_dataset_proximity.py` | F16b |
| Supp. 10 | `Supplementary/Supp_Fig_10_variability_levels.py` | F20 |
| Supp. 11 | `Supplementary/Supp_Fig_11_run_distance_qc.py` | F9, F10 |

Figures 4, 5 and Supp. 7-11 only assemble panels (the "Fig5" pipeline and its panel names
5a, 5b... predate the paper numbering): run `Fig5/Launcher_Fig5.py` first
(`RUN_PART='figures'` is enough once the matrices exist); the paper panels are written to
`<FIG5_ROOT>/06_figures/Paper/`.

## Shared modules

- `paper_style.py` -- colours (species, one fixed colour per dataset, lobes), species order, region order
  (atlas / Figure 1 sunburst), Nature sizes, panel letters, `compose_images`.
- `paper_data.py` -- datasets, paths, loading with cache, one value per subject, bilateral values,
  robust outliers.
- `ednix_bids_tools.py`, `EDNiX_figures.py`, `_0_extract_data_for_analysis.py` -- engines (unchanged).
- `Fig5/` -- fMRI selection / fingerprint pipeline (`Launcher_Fig5.py`), group seed maps
  (`SBA_group_maps.py`), paper panels (`fig5_supp_paper.py`).

## Not in the paper

| Folder | Content |
|---|---|
| `sensitivity_analysis/morphometry/` | `morphometry_all_levels.py`: levels 2-4, raw/log/normalised, L/R, cross-modality correlations (former `Figure_2_anat.py`); `outlier_rates.py`: regional outlier rates (former Supp. 6-8, removed from the paper) |
| `sensitivity_analysis/fmri_qc/` | `qc_all_metrics.py` (former `Supp_Figure_4-6_QC.py`), `fmri_acquisition_parameters.py`, `check_fd_inputs.py` (why mean FD is missing) |
| `sensitivity_analysis/longitudinal/` | `longitudinal_all_outputs.py`: every Fig. 3 output besides the paper figure |
| `phylo_preliminary/` | `PGLS_frontal_scaling.py` (former second half of `Figure_2_anat.py`) |
| `multimodal_PET/` | Figure XX (multimodal radar), PET group launchers, PET inter-species profiles |
| `Archives/` | former Figure 4 (network QC) and older exploratory scripts |
