# EDNiX paper -- analysis scripts

This folder holds the paper and nothing else at its top level: statistics, figures and one
runner. Every other analysis is filed by category below.

    python run_paper.py                    # statistics + figures, every figure
    python run_paper.py figure4 supp12     # some figures (stable1 = Supp. Table 1)
    python run_paper.py --stats-only       # numbers only
    python run_paper.py --refresh          # re-read the EDNiX outputs (no cache)

Before Figures 4-5 / Supp. 10-14: the run selection must exist
(`python Paper/run_selection/launcher.py` from the repository root; it writes
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
| Fig. 1 | atlas surfaces (Workbench scenes) + hierarchy (sunburst) | `library/figure1_surfaces.py`, `library/atlas_sunburst.py` |
| Fig. 2 | cortical morphometry, atlas level 2, bilateral, log10 | `library/morphometry_figure.py` |
| Fig. 3 | longitudinal morphometry (BIDS_Zhu_Garin, BIDS_Park) | `library/longitudinal.py` |
| Fig. 4 | the run selection: a common profile + each dataset vs its leave-one-species-out consensus, b selection vs random runs (noise, pair classes), c dataset similarity matrix r (F16a), d phylogenetic gradient, e species identity; run-to-run r (5b) in the tables only | `fmri_fingerprints/figure4_panels.py`, `fmri_fingerprints/selection_validation.py` |
| Fig. 4 (new layout, 2026-11-27) | on the 3D group maps: a phylogenetic signal of the all-runs maps (K_mult, Mantel, exact 7! P), b common + dataset profiles (all runs by default, `--profiles selected`), c similarity raw / residual; harmonisation figure (amplitude, within/between-species ratio + size null, seed maps and correlation matrices); Supp: selected vs random, noise bins, surfaces. `python -m Paper.analysis.fig4_maps.run <ROOT> <OUT>`; night run per invariance value: `Paper/run_selection/run_final.py 0.8` | `fig4_maps/` |
| Fig. 5 | its result: seed maps, selected vs excluded runs (5a) | `fmri_fingerprints/paper_panels.py` |
| Fig. 6 | multimodal profile: a-c FDG bratio (human, macaque), b / d reproducibility within species, e regional bratio vs the pooled macaques, f macaque - human, g multimodal fingerprints | `library/multimodal_figure.py`, `library/multimodal_profiles.py`, `library/pet_datasets.py` |
| Supp. 1 | schematic: spatial registration of the anatomical module (example session of BIDS_Zhu_Garin, else synthetic) | `library/schematic_anat.py` |
| Supp. 2 | schematic: spatial registration of the fMRI module | `library/schematic_func.py` |
| Supp. 3 | schematic: FDG-PET registration, scenarios 1v1 / 1v2a / 1v2b with the CT path (BIDS_Tremblay, BIDS_Merida) | `library/schematic_pet_ct.py` |
| Supp. 4 | acquisition parameters reported for every fMRI dataset (at most 4: TR, run duration, voxel volume, B0 / TE) | `library/acquisition_table.py` |
| Supp. 5 / 6 | anatomical / functional QC | `library/qc_figures.py` |
| Supp. 7 | EDNiX vs FreeSurfer recon-all, human (BIDS_Park); standalone script, not a `run_paper` section (writes `Supp_Fig_07/`) | `validation/supp_freesurfer_comparison.py` |
| Supp. 8 | longitudinal vs cross-sectional processing of the Fig. 3 datasets (BIDS_Park vs BIDS_Park_nolong, BIDS_Zhu_Garin vs its non-longitudinal copy): residual SD, slope SE, slope bias, paired subject bootstrap | `library/longitudinal_ablation.py` |
| Supp. 9 | human association vs primary cortex (BIDS_Park) | `library/human_assoc_primary.py` |
| Supp. 10 | a convergence (F1e), b consensus by dataset (F1c), c region invariance (F5) | `fmri_fingerprints/supp_panels.py` |
| Supp. 11 | a selection vs variability levels (F20), b raw species gap (selected vs random) | `fmri_fingerprints/` |
| Supp. 12 | proximity rows (F16b) | `fmri_fingerprints/supp_panels.py` |
| Supp. 13 | a selection scores (F3), b r with the criterion by arousal state (6c), c kept units per state | `fmri_fingerprints/` |
| Supp. 14 | a SD of run r with the kept map (F9), b QC vs that r (F10, random-effects meta-analysis) | `fmri_fingerprints/supp_panels.py` |
| Supp. 15 | similarity of the regional FDG profiles between pairs of subjects (same subject / other session, within dataset, between human datasets, human vs macaque) | `library/multimodal_figure.py` |
| Supp. Table 1 | acquisition parameters per dataset (BIDS sidecars + NIfTI headers; a fMRI, b anatomy) | `library/acquisition_table.py` |

Numbering of the V4 manuscript (supplementary figures renumbered by first citation, PR #88).
The section names of `run_paper.py` are the figure numbers (`supp8` -> `Supp_Fig_08/`); the
authoritative list is `paper_stats.SECTIONS` / `FOLDER` and the docstring of `paper_figures.py`.

The panel codes (F15b, 5b, ...) are those of the former Fig5 pipeline, kept so that its outputs
and the paper can be matched.

## Folders

| Folder | Content |
|---|---|
| `library/` | statistics and drawing engines of Figures 1-3, 6, Supp. 1-6, 8, 9, 15 and Supp. Table 1, called by the paper files. `pipeline_schematic.py` is the drawing engine of Supp. 1-3 (frames from a processed session, else from a synthetic brain built on the MNI ICBM152 template of nilearn, tagged 'synthetic'); each `schematic_*.py` also runs alone with `--synthetic` / `--out DIR` |
| `fmri_fingerprints/` | fMRI statistics and panels from the saved run selection: `context.py` (loads `07_state`), `settings.py`, `figure4_panels.py`, `fingerprints.py`, `profile_stats.py`, `paper_panels.py`, `supp_panels.py`, `phylo.py` |
| `validation/` | `paper_vs_pipeline.py`: checks on synthetic data that the paper numbers equal those of the former Fig5 pipeline. `paper_claims.py` + `paper_claims_registry.py`: traceability of the manuscript (see below) |
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

## Traceability: manuscript number -> stats cell -> code line -> data

    python run_paper.py                                         # (re)compute every stats file
    python validation_morphometry.py                            # Fig. 2 validation numbers
    python paper_extra_stats.py                                 # Fig. 2 scaling, Fig. 6c
    python validation/paper_claims.py /path/to/EDNiX_PaperV3_Claude_v11.docx

`validation/paper_claims_registry.py` lists each number of the text (the words just before
it, the stats file / sheet / row / column that holds it, the function and line that writes
it, the datasets it must come from) and the Methods statements to check against the code.
`paper_claims.py` reads the .docx (tracked insertions included), finds each number, compares
it with the stats file at the precision written, locates the code line (AST), flags stats
older than the code, and lists every number of the Results and legends that no entry
covers. Report: `<PAPER_DIR>/Traceability/paper_claims_report.xlsx|md`.
Status: OK, MISMATCH, PENDING ([XX] in the text: the computed value is given), NOT_IN_TEXT
(sentence rewritten: update the prefix), NO_SOURCE (stats not computed yet).
Covered now: Fig. 2 validation, Fig. 3, Fig. 4 b-e, Fig. 6, Supp. 7, 8, 9, 11b, 13c. To add: Fig. 4a, 5, Supp. 10, 12, 14, 15.
