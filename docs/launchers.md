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

## Things to know

- `preprocess_fmri`, `preprocess_pet` and `preprocess_ct` accept unknown keyword arguments and ignore them (they print a warning): check the spelling of parameters. The interface (`Check`) and the unit tests detect them.
- Some example functional launchers call `_000_Data_Reset.reset_all(..., dry_run=False)`, which deletes previous results of the dataset each time they run. The interface warns about it.
- `list_transfo` / `transfo_message` control the registration parameters of each step (`align`, `between`, `coreg`, `SS1`, `SS2`, `SS3`, `stdyT`); `modalities/anat/load_transfo_parameters.build_transfos` builds them with defaults.

## With the interface

The **Launchers** tabs (Anatomy, Functional, PET, CT) open any launcher as a form: each parameter with its description, `Skip_step` as a checklist of steps, a check for unknown parameters, saving (the previous version is kept) and running in the background with live logs. See [GUI/README.md](../GUI/README.md).
