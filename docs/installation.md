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
bash Tools/setup_venv.sh            # creates EDNiX/venvEDNiX with the exact versions of requirements-lock.txt
source venvEDNiX/bin/activate
```

| file | content |
|---|---|
| `requirements.txt` | every third-party package imported by the code, with compatible ranges |
| `requirements-lock.txt` | exact versions of all packages (direct and transitive) of the validated environment |

Options of `Tools/setup_venv.sh`: `--cpu` (PyTorch without CUDA), `--dev` (pytest), `--venv DIR`,
`--python python3.11`, `--unlocked` (newest versions allowed by `requirements.txt`, to test an upgrade).
Without `python3.X-venv` (no `ensurepip`, common on Debian/Ubuntu without root access), the script falls back to
`virtualenv` if available, else to a venv created `--without-pip` with pip installed from its PyPI wheel.
The script checks the installation (`pip check`, import of every package, `Tools/check_requirements.py`)
and writes the installed versions to `venvEDNiX/ednix-freeze.txt`.

`venvEDNiX` at the root of the EDNiX folder is where the interface (`GUI/launch_EDNiX_GUI.sh`) looks for Python; `run_all_ednix.py` uses the interpreter that runs it.
The U-Net skull-stripping (`muSkullStrip.py`, methods `muSkullStrip_*`) is run with the `python3` found on
the `PATH`: activate the venv before running a launcher, so that it finds `torch`.

After a change of environment that has been validated (tests + smoke test), update the lock from it:

```bash
bash Tools/setup_venv.sh --freeze --venv venvEDNiX
```

`tests/unit/test_requirements.py` fails when the code imports a package missing from either file.

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
    tpcclib.sif                      (PET kinetic models: Turku PET Centre tools)
```

The U-Net skull-stripping models are read from the same folder. FreeSurfer needs its licence file in the container environment.

Programs outside pip and the containers:

| program | used by | required |
|---|---|---|
| Singularity / Apptainer | every neuroimaging step | yes (or the programs installed natively, `singularity='no'`) |
| NHP-BrainExtraction (`muSkullStrip.py` + models) in `Singularity_library/Singularity/` | skull-strip methods `muSkullStrip_*` | when these methods are used |
| `petpvc` (PETPVC) on the `PATH` | PET partial-volume correction, `pvc_engine='petpvc'` (or `'auto'`) | no: `pvc_engine='python'` gives RBV and MG without it |
| FreeSurfer licence file | FreeSurfer container | yes |
| `git` | provenance written by the atlas library builder | no |

`pydicom` (in the venv) reads the isotope and scan start from the DICOM files (`dcmfolder`).

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
