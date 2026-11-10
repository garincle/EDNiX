# Launchers

A launcher is a short Python script that describes one dataset (paths, species, subjects) and calls a pipeline entry point:

| Pipeline | Entry point | Examples |
|---|---|---|
| Anatomy | `modalities.anat._0_Pipeline_launcher.preprocess_anat` | [`Exemples/Study_EDNiX/Launcher_anat`](../Exemples/Study_EDNiX/Launcher_anat) |
| Functional MRI | `modalities.fMRI._0_Pipeline_launcher.preprocess_fmri` | [`Exemples/Study_EDNiX/Launcher_func`](../Exemples/Study_EDNiX/Launcher_func) |
| PET | `modalities.pet._0_Pipeline_launcher.preprocess_pet` | [`Exemples/Study_EDNiX/Launcher_pet`](../Exemples/Study_EDNiX/Launcher_pet) |
| CT (for the PET) | `modalities.ct._0_Pipeline_launcher.preprocess_ct` | [`Exemples/Study_EDNiX/Launcher_ct`](../Exemples/Study_EDNiX/Launcher_ct) |

The example launchers cover bat, dog, human, mouse lemur, macaque, marmoset, mouse and rat datasets. Start from the one closest to your species and acquisition.

## Minimal edits

1. `MAIN_PATH`: your EDNiX folder, **ending with `/`**.
2. `bids_dir`: your BIDS dataset.
3. `species`: a species folder name of the atlas library (`Mouse`, `Rat`, `Bat`, `Mouselemur`, `Marmoset`, `Macaque`, `Dog`, `Human`...).
4. `list_to_keep` / `list_to_remove`: `(subject, session)` pairs, see [bids.md](bids.md).
5. `Skip_step`: steps not to run (see the step tables of the [anatomical](reference/anat_parameters.md), [functional](reference/fmri_parameters.md), [PET](reference/pet_parameters.md) and [CT](reference/ct_parameters.md) references).

Every parameter is described in the parameter references ([anatomy](reference/anat_parameters.md), [functional](reference/fmri_parameters.md), [PET](reference/pet_parameters.md), [CT](reference/ct_parameters.md)), and in the launcher editor of the interface.

## Order of the runs

For a dataset, run the anatomical launcher first: the functional and PET pipelines reuse the anatomical transforms, masks and atlases (unless `IhaveanANAT=False`). A PET registered through a CT needs the CT launcher before the PET launcher ([pet_pipeline.md](pet_pipeline.md)).

## Running

```bash
EDNIX_STRICT=1 python my_func_launcher.py
```

- Every external command goes through `Tools/run_cmd.py`, which logs it with its output and return code in the subject's log (`diary`).
- Without `EDNIX_STRICT`, a failing command is logged in red and the pipeline continues; with `EDNIX_STRICT=1`, it stops with an error. Use strict mode for new datasets and for validation runs.
- `EDNIX_SDC=0` disables the distortion correction of the functional time series (see [functional_pipeline.md](functional_pipeline.md)).

## Resuming and skipping

`Skip_step` lists the steps to skip. The example launchers often skip the steps already done on their data: when you start a new dataset, begin with an empty or minimal `Skip_step` (the interactive ITK-SNAP steps `itk_1`, `itk_2`, `itk_3` can be skipped for unattended runs).

## Cleaning before a relaunch

Every example launcher has a reset block before `Skip_step`:

```python
RESET_BEFORE_RUN = False   # True: run the cleaner of this modality on the selected subjects first
RESET_DRY_RUN = True       # True: list what would be removed; False: remove, after ONE 'YES' for the whole batch
```

| modality | cleaner | removed | kept |
|---|---|---|---|
| anat | `modalities/anat/_000_Data_Reset` | everything EDNiX wrote under `anat/` (`RESET_SCOPE='surfaces'`: only `freesurfer/` and `native/surfaces/`, to rerun steps 10–15) | raw BIDS files, `*_final_mask*` (itk_1/itk_2), `*_backtonative_mask*` (itk_3), `*_desc-manual*`, `segmentation.mgz`, `*.scene` |
| func | `modalities/fMRI/_000_Data_Reset` | EDNiX outputs under `func/` | raw BIDS files, `func/raw/masks` |
| PET / CT | `modalities/pet/_000_Data_Reset`, `modalities/ct/_000_Data_Reset` | EDNiX outputs under `pet/` / `ct/` (`RESET_KEEP_REGISTRATION=True`: keeps `orig/matrices`) | raw BIDS files, manual crops `*_desc-crop_*` |

Files are removed one by one (manual files stay where the pipeline looks for them), then empty folders are pruned. Every action is written to the subject's diary.

## Full rebuild in one command: `rebuild_all_ednix.py`

`Exemples/Study_EDNiX/rebuild_all_ednix.py` chains everything with one confirmation: it previews the library build, the BIDS reset (dry run, manual masks and crops kept) and the jobs, asks one `YES`, then detaches (survives closing the terminal) and logs to `run_reports/rebuild_<date>.log`. It builds `Atlases_library_V2` (`atlases/build_library_V2.py`), stops if a species fails or if the library does not match its manifest (nothing in the BIDS folders is touched before that), swaps the libraries (`Atlases_library` → `Atlases_library_V1`, `Atlases_library_V2` → `Atlases_library`), resets the BIDS folders and runs `run_all_ednix.py`. If `Atlases_library` is already the V2 library, the build is skipped (resume).

```bash
python Exemples/Study_EDNiX/rebuild_all_ednix.py --plan     # preview only
python Exemples/Study_EDNiX/rebuild_all_ednix.py            # preview, YES, then everything in the background
```

## Relaunching everything: `run_all_ednix.py`

`Exemples/Study_EDNiX/run_all_ednix.py` runs the example launchers as a dependency graph: per BIDS folder, anat before CT, CT before PET, anat before func; then the QC reports of every dataset and a study-level Excel summary.

```bash
python Exemples/Study_EDNiX/run_all_ednix.py --plan                 # jobs, dependencies, core budget; runs nothing
python Exemples/Study_EDNiX/run_all_ednix.py --only Human --phases anat,func
python Exemples/Study_EDNiX/run_all_ednix.py --reset                # one preview + one YES for all the resets, then run
```

- Launchers are run as they are (their own `Skip_step` and parameters). A launcher is split by subject into up to `MAX_SPLIT` workers of `THREADS_PER_WORKER` threads within the cores of the machine (minus `CORE_MARGIN`); the longest jobs start first. Study-template runs are not split.
- A job whose prerequisite failed is not started; the others continue. One log per job and a `report.txt` (with the end of each failed log) are in `run_reports/<run>/` of the EDNiX folder.
- Workers run with `EDNIX_SKIP_RESET=1`: the resets of the launchers are skipped there (they are done once, interactively, by `--reset`).

## Things to know

- `preprocess_fmri`, `preprocess_pet` and `preprocess_ct` accept unknown keyword arguments and ignore them (they print a warning): check the spelling of parameters. The interface (`Check`) and the unit tests detect them.
- `Launcher_func_Marmoset_Hori.py` and `Launcher_func_marmoset_MBM.py` have `RESET_BEFORE_RUN = True, RESET_DRY_RUN = False`: they delete the previous results of the dataset each time they run (after one `YES`). The interface warns about it.
- `list_transfo` / `transfo_message` control the registration parameters of each step (`align`, `between`, `coreg`, `SS1`, `SS2`, `SS3`, `stdyT`); `modalities/anat/load_transfo_parameters.build_transfos` builds them with defaults.

## With the interface

The **Launchers** tabs (Anatomy, Functional, PET, CT) open any launcher as a form: each parameter with its description, `Skip_step` as a checklist of steps, a check for unknown parameters, saving (the previous version is kept) and running in the background with live logs. See [GUI/README.md](../GUI/README.md).
