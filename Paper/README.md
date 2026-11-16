# Paper: reproducing the EDNiX article

This folder is not needed to use EDNiX. It holds what produced the figures, tables and statistics of the article, kept for reproducibility.

| Folder | Content |
|---|---|
| [`analysis/`](analysis/README.md) | every figure, table and statistic: `python Paper/analysis/run_paper.py` (index figure -> code in `analysis/README.md`), sensitivity analyses, validation against FreeSurfer, checks of the numbers quoted in the text (`analysis/validation/paper_claims.py`) |
| `run_selection/` | configuration of the published datasets for the run selection (`config_ednix.py`), its launcher (`launcher.py`), the SBA group maps of the paper (`sba_group_maps.py`) and the equivalence checks done when the code was moved (`validation/`). The method itself is `Statistics/Group_fMRI/run_selection/` |
| `library_build/` | scripts that built the atlas library distributed with EDNiX from the original species atlases (`build_library_V2.py` and the one-off conversions), and the Workbench scenes of Figure 1 |
| `rebuild_all_ednix.py` | rebuilds the library, resets and reprocesses every example dataset (maintainers) |
| `tests/` | tests of the figure code on synthetic data (run with the rest of the suite) |

## Order

1. Preprocess every dataset with the example launchers (`Exemples/Study_EDNiX/Launcher_*`, same parameters as in the article), for example with `python Exemples/Study_EDNiX/run_all_ednix.py`, `EDNIX_DATA` pointing to the folder of the datasets (see [docs/launchers.md](../docs/launchers.md)).
2. Run selection: `python Paper/run_selection/launcher.py`.
3. Statistics and figures: `python Paper/analysis/run_paper.py`.

## Paths

The scripts of this folder keep the paths of the authors' server (`/scratch2/EDNiX/...` for the data and results, the EDNiX folder for the library); they are gathered in `analysis/paper_data.py` and `run_selection/config_ednix.py`. Change them there to rerun the analyses elsewhere. The datasets are public (sources in the article, Supplementary Table 1); the atlas library and the containers are on Zenodo (see [docs/installation.md](../docs/installation.md)).
