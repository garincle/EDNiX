# Installation

## Requirements

- Linux (the pipeline runs its neuroimaging programs in Singularity / Apptainer containers).
- Python 3.10 or later (tested with 3.11 and 3.12).
- Singularity or Apptainer.
- Disk space for the atlas library and the containers (several GB), plus the outputs of your datasets.

## 1. Code and Python packages

```bash
git clone https://github.com/garincle/EDNiX.git
cd EDNiX
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

A conda environment works as well. `requirements.txt` lists the Python packages without pinned versions; once your environment is validated, record it:

```bash
pip freeze > requirements-lock.txt
```

## 2. Containers

The pipeline calls AFNI, ANTs, FSL, FreeSurfer, Connectome Workbench, ITK-SNAP and SynthStrip through Singularity images (`Tools/Load_EDNiX_requirement.py`). They are expected in:

```
<EDNiX folder>/Singularity_library/Singularity/
    afni_ub24_latest.sif
    fsl_6.0.5.1-cuda9.1.sif
    freesurfer_NHP.sif
    connectome_workbench_1.5.0-freesurfer-update.sif
    itksnap_5.0.9.sif
    synthstrip.1.5.sif
    tpcclib.sif                      (PET)
```

The U-Net skull-stripping models are read from the same folder. FreeSurfer needs its licence file in the container environment.

> `MAIN_PATH` must end with `/` in the launchers: the container folder is derived from it with `os.path.dirname`, which returns the parent folder when the trailing slash is missing.

## 3. Atlas library

Templates, atlases, label tables and reference surfaces live in `<EDNiX folder>/Atlases_library/`. Its layout is described in [atlases.md](atlases.md).

## 4. Check the installation

```bash
pip install pytest pytest-xdist
pytest -n auto
```

All tests should pass (a few minutes). Then run one example launcher on one subject with `EDNIX_STRICT=1`, or use the end-to-end smoke test ([testing.md](testing.md)).

## 5. Interface (optional)

```bash
bash GUI/install_shortcut.sh
```

creates an **EDNiX** icon on the desktop and in the applications menu. See [GUI/README.md](../GUI/README.md).

## Local data folders

The following folders sit inside the EDNiX folder but are not part of the repository (they are listed in `.gitignore`): `Atlases_library/`, `Singularity_library/`, `Tool_library/`. When you update EDNiX from a new copy of the code, keep these folders.
