# EDNiX paper -- analysis scripts

This folder holds the paper and nothing else at its top level: statistics, figures and one
runner. Every other analysis is filed by category below.

    python run_paper.py                    # statistics + figures, every figure
    python run_paper.py figure4 supp10     # some figures (stable1 = Supp. Table 1)
    python run_paper.py --stats-only       # numbers only
    python run_paper.py --refresh          # re-read the EDNiX outputs (no cache)

Before Figures 4-5 / Supp. 8-12: the run selection must exist
(`python Statistics/Group_fMRI/run_selection/launcher.py` from the repository root; it writes
`<OUTPUT_ROOT>/07_state`, read by the paper, `OUTPUT_ROOT` from `config_ednix.SELECTION_CONFIG`).

Outputs: `/scratch2/EDNiX/results/EDNiX_Figures/Paper/<Figure>/` (`paper_data.PAPER_DIR`):
`<Figure>.pdf|png|svg` (fonts embedded, PNG 600 dpi, Nature widths 89 / 120 / 183 mm),
`<Figure>_stats.xlsx` + `<Figure>_stats.txt` (every number quoted in the text and legend),
`panels/` (single panels of multi-panel figures).

## Paper files

| File | Role |
|---|---|
| `run_paper.py` | runner: `paper_stats.run` then `paper_figures.run` |
| `paper_stats.py` | every statistic, one function per figure, in figure order (`SECTIONS`); `Inputs` = lazy, cached inputs |
| `paper_figures.py` | every figure, drawn from the results of `paper_stats` (nothing computed there) |
| `paper_style.py` | colours (species; a dataset = its species colour + its own marker), lobes, species / region order, Nature sizes, panel letters, `compose_images` |
| `paper_data.py` | datasets, paths, loading with cache, one value per subject, bilateral values |

## Figure -> functions

| Figure | Content | Statistics / drawing engine |
|---|---|---|
| Fig. 1 | atlas hierarchy (sunburst) | `library/atlas_sunburst.py` |
| Fig. 2 | cortical morphometry, atlas level 2, bilateral, log10 | `library/morphometry_figure.py` |
| Fig. 3 | longitudinal morphometry (BIDS_Zhu_Garin, BIDS_Park) | `library/longitudinal.py` |
| Fig. 4 | the run selection: a common profile + each dataset vs its leave-one-species-out consensus, b run-to-run r (5b), c dataset similarity matrix r (F16a) | `fmri_fingerprints/` |
| Fig. 5 | its result: seed maps, selected vs excluded runs (5a) | `fmri_fingerprints/paper_panels.py` |
| Fig. 6 | multimodal profile, Human / Macaque (former Figure XX) | `library/multimodal_figure.py`, `library/multimodal_profiles.py` |
| Supp. 1 | schematic: spatial registration of the anatomical module (example session of BIDS_Zhu_Garin, else synthetic) | `library/schematic_anat.py` |
| Supp. 2 | schematic: spatial registration of the fMRI module | `library/schematic_func.py` |
| Supp. 3 | schematic: FDG-PET registration, scenarios 1v1 / 1v2a / 1v2b with the CT path (BIDS_Tremblay, BIDS_Merida) | `library/schematic_pet_ct.py` |
| Supp. 4 | acquisition parameters reported for every fMRI dataset (at most 4: TR, run duration, voxel volume, B0 / TE) | `library/acquisition_table.py` |
| Supp. 5 / 6 | anatomical / functional QC | `library/qc_figures.py` |
| Supp. 7 | human association vs primary cortex (BIDS_Park) | `library/human_assoc_primary.py` |
| Supp. 8 | a convergence (F1e), b consensus by dataset (F1c), c region invariance (F5) | `fmri_fingerprints/supp_panels.py` |
| Supp. 9 | a selection scores (F3), b r with the criterion by arousal state (6c) | `fmri_fingerprints/` |
| Supp. 10 | proximity rows (F16b) | `fmri_fingerprints/supp_panels.py` |
| Supp. 11 | selection vs variability levels (F20) | `fmri_fingerprints/` |
| Supp. 12 | a SD of run r with the kept map (F9), b QC vs that r (F10) | `fmri_fingerprints/supp_panels.py` |
| Supp. Table 1 | acquisition parameters per dataset (BIDS sidecars + NIfTI headers; a fMRI, b anatomy) | `library/acquisition_table.py` |

The panel codes (F15b, 5b, ...) are those of the former Fig5 pipeline, kept so that its outputs
and the paper can be matched.

## Folders

| Folder | Content |
|---|---|
| `library/` | statistics and drawing engines of Figures 1-3, 6, Supp. 1-7 and Supp. Table 1, called by the paper files. `pipeline_schematic.py` is the drawing engine of Supp. 1-3 (frames from a processed session, else from a synthetic brain built on the MNI ICBM152 template of nilearn, tagged 'synthetic'); each `schematic_*.py` also runs alone with `--synthetic` / `--out DIR` |
| `fmri_fingerprints/` | fMRI statistics and panels from the saved run selection: `context.py` (loads `07_state`), `settings.py`, `figure4_panels.py`, `fingerprints.py`, `profile_stats.py`, `paper_panels.py`, `supp_panels.py`, `phylo.py` |
| `validation/` | `paper_vs_pipeline.py`: checks on synthetic data that the paper numbers equal those of the former Fig5 pipeline |
| `sensitivity_analysis/morphometry/` | `morphometry_all_levels.py`: levels 2-4, raw/log/normalised, L/R, cross-modality correlations; `outlier_rates.py`: regional outlier rates (former Supp. 6-8) |
| `sensitivity_analysis/fmri_qc/` | `qc_all_metrics.py`, `fmri_acquisition_parameters.py`, `check_fd_inputs.py` |
| `sensitivity_analysis/longitudinal/` | `longitudinal_all_outputs.py`: every Fig. 3 output besides the paper figure |
| `sensitivity_analysis/fmri_selection/` | exploration of the run selection (every F-panel, matrices, seed vs matrix): `Launcher_exploration.py` |
| `phylo_preliminary/` | `PGLS_frontal_scaling.py` |
| `multimodal_PET/` | PET group launchers (CERMEP, Tremblay), PET inter-species profiles |
| `Archives/` | former Figure 4 (network QC) and older exploratory scripts |

Engines kept unchanged at this level: `ednix_bids_tools.py`, `EDNiX_figures.py`,
`_0_extract_data_for_analysis.py`.

The run selection itself (AFNI group maps, consensus selection, `07_state`) is a standalone
module at the repository root: `Statistics/Group_fMRI/run_selection/` (`launcher.py`,
`config_ednix.py`, `compare_states.py`).
