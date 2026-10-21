# PET / CT import: review journal (October 2026)

Review and correction of the PET and CT code delivered by the PET engineer ("last version, not tested"), and of his suggested changes to `Tools/`, `atlases/` and `modalities/anat/`. Git history of the branch `pet-ct-import`:

| commit | content |
|---|---|
| `25ed18e` | EDNiX before the import (old `modalities/pet` removed) |
| `91ebd01` | engineer's files **as received** (reference point: `git diff 91ebd01 -- modalities/pet` shows every change) |
| next commits | corrected version, tests, interface, documentation |

**Summary.** The engineer's version could not run: about 40 blocking errors (imports of modules that do not exist, functions called with the wrong arguments, undefined variables, a syntax error in the CT launcher), and several scientific errors in what would have run (reference TAC, decay handling of the kinetic fits, transform chains, motion correction). Its **design** is sound and is kept: processing in the PET space with atlases and masks brought into it (`Illustration/pet2b.svg`), the call tree (`pet3b.svg`) and the registration scenarios (`pet4.svg`). The modules were therefore rewritten one by one with the same names and roles; three helper modules were added (`petutils.py`, `transfo.py`, `turkutools/_nifti.py`), two were merged into others (`baseline_def.py` -> `atlases/getlabels.py`, `getsif.py` -> `sifFile.py`).

Decisions taken with the user (21/10/2026): every scenario of `pet4.svg` made operational, with the limits marked there (2v2 / 3v2 up to the study template only); hybrid PET-MR (1v1) handled **as the functional pipeline** (`anat_func_same_space`): the PET is reoriented like the anatomy and follows its AC-PC matrix; the **interactive steps are kept** (freeview crops and checks, tkinter window for the reference region); example launchers for the macaques of BIDS_Tremblay and the human CERMEP dataset.

Legend: **[blocking]** the code stopped there; **[results]** the code would have run with a wrong result; **[design]** change of behaviour or structure.

## 1. Bridge with the anatomy (the core change)

| engineer | correction |
|---|---|
| native -> template transforms hard-coded, under two different names (`native_to_<T>_Final_*` in `coreg2anat`, `anat2<T>_Final_*` in `norm2template`, the second never written by EDNiX) **[blocking]** | chains of `modalities/anat/loop3/backtonative.get`, as the functional pipeline: reference template, study template and reference session (longitudinal) are handled |
| PET -> template = `[shift, acpc]` + native -> template, with the same `[affine, inverse warp]` list and flags for both directions **[results]** | the AC-PC matrix of the anatomical pipeline already contains the shift (it is the initial transform of the AC-PC registration, collapsed by ANTs): `shift` is not composed again. Forward (`[warp, affine]`) and inverse (`[affine^-1, inverse warp]`) chains are distinct (`transfo.py`) |
| no anatomy: `print('Sorry not ready yet....')` | scenarios 2v1, 2v2, 3v1, 3v2 implemented (`coreg2anat.register`) |
| transforms recomputed by every step from file names | **registration record** (`orig/matrices/<session>_desc-registration_pet.json`): target spaces, reference images, masks, both chains; read by steps 7-10 and by the PET when it goes through the CT |
| anatomy looked up with `chooseanat.retrieve4pet` (does not exist) **[blocking]** | `chooseanat.retrieve`, with the same fallback on the other sessions of the subject as the functional pipeline |
| one registration per run | one per session, on the mean image of the session after the between-runs correction (all runs share that space) **[design]** |

## 2. File by file

### `modalities/pet/_0_Pipeline_launcher.py`
- imported `Exemples.Study_EDNiX.Launcher_pet.Launcher_fdg_CERMEP` (runs a launcher at import, `baseline` not defined) **[blocking]**; called the modules `loop1(...)` and `loop2(...)` as functions **[blocking]**; undefined variables `anat_raw`, `Need_check_orient`, `ref_CT`, `templateT1w`, `templatemask`, `atlas_file`, `atlas_label`, `studyTfile` (if step 2 skipped), `IhaveaCT` passed to `loop2` **[blocking]**.
- `participants.tsv`: `ID.split('-')[1]` (EDNiX IDs have no `sub-` prefix) **[blocking]**, non-BIDS columns `id` / `session_name`, `.values[0]` on an empty selection. Now `petutils.subject_values`: PET sidecar (`InjectedRadioactivity` + units), `sessions.tsv`, `participants.tsv` (BIDS columns, and the engineer's `weight_kg`, `tracor_qte_MBq`, `blood_concentration`), values and sources written in `Prepro_Settings.json`.
- reference region: `getlabels.gettype(..., '')` then `pd.read_csv('')` **[blocking]**, `pd.read_csv(header=[...])` misused, atlas JSON `atlas_config_V3.json`; now `atlases/getlabels.resolve` (launcher value or interactive choice, saved for the study).
- step messages all said `skip step 1`; ~15 parameters never used (`endjson`, `Method_mask_pet`, `doWARPonfunc`, `n_for_ANTS`...) removed.
- entry point renamed `preprocess_pet` (the interface and the documentation generator look for it; `preprocess_data` kept as an alias), full docstring (every parameter documented in the interface).

### `loop1.py` (step 1)
- raw images searched in `pet/raw/` with `ID + '_trc-'...` (BIDS raw images are in `pet/`, named `sub-<ID>_ses-<ses>_...`): nothing was found **[blocking]**.
- frame-timing JSON written with `opj(dir, img.replace(...))` on an absolute path: written **inside the raw BIDS folder** **[results]**; now `pet/raw/<run>_desc-frametiming_pet.json`.
- `REF_image` undefined when neither the anatomy nor the CT shares the PET space **[blocking]**; `if not pethd:` reoriented only the images with a *wrong* header (inverted test): correct images were never reoriented and step 4 found nothing **[blocking]**.
- `FrameReferenceTime` required (not a BIDS field): the frame middles are used when it is absent.
- SIF never created (`getsif.py` was a stub); now created at step 1 for the dynamic runs when `srtm=True`.

### `reorient.py` (shared with the CT)
- `animalPosition == 'na'` while EDNiX uses `''` / `'humanlike'`: `getreal` received `''` **[blocking]**; `orientation` (forced code, as in the anatomy) not supported **[design]**.
- 3dWarp `-wsinc5` interpolation: over/undershoots of the activity concentrations; now `-linear` **[results]**.
- `PET_header_wrong`: the orientation code of the anatomy was copied, but the PET was then not reoriented like the anatomy, so it did not reach the space of the reoriented anatomy **[results]**; now copy then the same reorientation.
- mean images weighted by the frame durations (time average of the concentration).
- `extract_filename.change` did not exist in EDNiX.

### `MoCo.py` (steps 4, 5)
- `apply_transforms(fixed=mcdir, ...)`: a folder as reference image **[blocking]**.
- `if i > nb_void`: one frame more than asked left out; these frames were not corrected. Now `>=`, and they receive the transform of the first registered frame (nearest in time) **[results]**.
- single pass on the mean of uncorrected frames (blurred reference): two passes **[results]**.
- `static()` called with 4 arguments for 8 **[blocking]**; 3dMean of 4D runs of different lengths **[blocking]**; crop file name built two different ways (file not found) **[blocking]**; 4D images resampled without `imagetype=3` **[blocking]**; outputs written in `orig/` and read from `orig/preprocessing/` **[blocking]**; SIF copied without checking that it exists **[blocking]**.
- between runs, a dynamic run was resampled twice (step 4, then step 5): the per-frame and between-runs transforms are now composed and the reoriented frames resampled once **[results]**.
- added: motion parameters and mean displacement per frame (`forMC/*_motion.tsv`, `.png`).

### `studyT.py` (steps 2, 3)
- `os.mkdir(..., exist_ok=True)` **[blocking]**; `ants.build_template` and `ants.registration` given file names instead of images **[blocking]**.
- images not scaled before averaging (a high dose weighed more): scaled by their 99th percentile **[results]**.
- the manual crop file was overwritten at every run (crop lost): created once, `redo_crop` to draw it again.
- automatic crop: `[affine, inverse warp]` applied without inverting the affine (wrong direction) and template intensities thresholded at 0.1 instead of the drawn box **[results]**; now the box of the crop through the 'from' chain.

### `coreg2anat.py` (step 6)
- `from anat.freesurfer import FS_freeview` **[blocking]**; `transPET = transPET + T_PET` inside the function (UnboundLocalError) **[blocking]**; `os.mkdirs` **[blocking]**; initial transform read from another folder than the one it was written to **[blocking]**; `full_IMG` undefined for a static image **[blocking]**; `Tf = T1; Tf.append(...)` modified `T1` for the later steps **[results]**; the atlas dictionary passed where a file was expected **[blocking]**.
- masks: genericLabel then 3dresample NN (two resamplings); now one resampling (masks: linear + threshold 0.5, atlases: genericLabel, every level of a 4D atlas).
- PET -> anatomy: rigid, mutual information, restricted to the dilated brain mask of the anatomy, registered to the AC-PC anatomy directly (no composition with the raw space).
- manual correction (`do_manual`, freeview `-transform-volume`) kept, re-estimated from the edited image starting from the current registration.
- QC: PET in every target space with the brain mask as contour (`pet/QC`, `<BIDS>/QC/pet_in_acpc`, `pet_in_template`).

### `norm2template.py` (step 10)
- `from Tools import filename, setTransfo` (do not exist) **[blocking]**; transform names never written by EDNiX **[blocking]**.
- export through the record (one resampling per image), `-SS` with the mask of the target space, export grid at the PET resolution by default (`export_resolution`), atlases and masks of the template spaces on the same grid.

### `preproPET.py`, `paramImg.py`, `cut_image.py` (step 7)
- implicit relative imports (`import roi_signal, ...`) and `Tools.filename` **[blocking]**; brain mask searched as `*desc-coreg_mask.nii.gz`, never written **[blocking]**.
- time means: `nii = nii / len(...)` **inside** the loop (wrong average) **[results]**; unweighted averages of frames of different durations **[results]**: now duration-weighted, with `SUV_window` for dynamic acquisitions.
- SUV: image units assumed Bq/mL (now read from `Units`; an image already in SUV is kept), injected dose not decayed to the decay-correction time of the image, no decay correction of non-decay-corrected data (static or dynamic) **[results]**; injected radioactivity read for every run (two injections in one session).
- glucose: the SUV was multiplied by the glycaemia in g/L without reference value **[results]**; now `SUVglu = SUV x glucose / 5.0 mmol/L`, kept apart from the SUV.
- time blocks: frames assigned with indices shifted by one, divided by a "brain mean" that included the zeros outside the brain (`ants.mask_image(...).mean()`) **[results]**.
- `srtm`: `run_cmd(cmd, ...)` (calling the module) **[blocking]**; `'-expr "a*b"'` glued to the prefix **[blocking]**; `white+csf` = product of the masks (intersection: empty) **[results]**, now their union; `selectmsk['white+csf']` was a folder.

### `roi_signal.py`, `dftFile.py`, `Tools/roiStats.py` (steps 7, 8)
- reference TAC = average of the regional means of the chosen labels: a small region weighed as much as a large one **[results]**; now the mean over the union of the voxels.
- reference from a mask: **sum** of the activity instead of the mean concentration **[results]**.
- DFT volumes = voxel counts instead of mm³; frame times truncated to integers; unit always `kBq/mL` **[results]**; now volumes in mm³, exact times, unit of the image.
- `bfmsrtm -i=<isotope>` for every image: this option is for **non-decay-corrected** data (the tool then models the decay) and biases the fit of decay-corrected images **[results]**; now given only when `ImageDecayCorrected` is false, and the voxel-wise tools (which require decay-corrected data) are then skipped.
- 3dROIstats wrapper (`Tools/roiStats.py`: `from tools import`, `img.unique()` on a file name **[blocking]**, reordering of the labels by hand) rewritten in numpy: every value is attached to its label; TAC tables carry the frame times.
- the voxel-wise models run on the skull-stripped dynamic image (fewer voxels, no fit in the background).

### `turkutools/` (`turkusif`, `turkubfbp`, `turkustrm`, `turkudv`)
- `from Tools import filename` **[blocking]**; inputs decompressed in place (`gzip -d` / `gzip -f`): a failure in between left the pipeline without its `.nii.gz`; now temporary copies (`_nifti.TurkuWorkspace`) inside the BIDS folder (bound in the container), and outputs written back as `.nii.gz` with the geometry of the input (voxel order restored when the tool writes a valid header).
- `imgdv -n=` built with `path.split('.')[0]` (paths with dots); figures `range(16)` crashed with fewer than 16 regions.
- `desc-STRM` renamed `desc-SRTM` (simplified reference tissue model, Lammertsma & Hume 1996) **[design]**.

### `sifFile.py`, `dcmFile.py`
- `os.path.exist` (AttributeError at import) and `from motilities...` **[blocking]**.
- isotope from the DICOM code meaning: `'^18^Fluorine'.split('^')` gave `'-18'` **[results]**; now from the radionuclide half-life (0018,1075), then the code meaning; BIDS `TracerRadionuclide` used first.
- SIF frame times no longer truncated; scan start from DICOM or BIDS `TimeZero` + `ScanStart`.

### `pvc.py`
- the point-spread function was the **voxel size** **[results]**: it is the resolution of the scanner and reconstruction, now a launcher parameter (`pvc_fwhm`); segmentation from the grey / white / CSF masks in the PET space; wired as step `pvc`.

### `cut_func.py`, `run_cut.py` (step 9, PET during fMRI)
- `from pet import ...` **[blocking]**; `Settings.json` of the functional pipeline (no longer written) **[blocking]**; clock times compared with `time.mktime` on strings; whole frames taken on both sides. Now overlap-weighted mean of the frames during the run (BOLD `AcquisitionTime`, `RepetitionTime`, `T1_eq`; PET `TimeZero`), in the PET space, exported at step 10.

### `baseline_def.py`, `getsif.py` (removed)
- `baseline_def`: search of the reference labels across the levels of an atlas JSON; it wrote `_test.nii.gz` and returned another name. The reference region is now its own mask (labels of any level), so no modified atlas is needed; the name search is `getlabels.find_labels`.
- `getsif`: stub with undefined names.

### CT (`modalities/ct`)
- `_0_Pipeline_launcher.py`: **syntax error** (no `:` after the signature) **[blocking]**; `chooseanat.retrieve4pet` **[blocking]**; only the case with an anatomy. Now `preprocess_ct`: reorientation, CT study template (2v2), registration to the anatomy (rigid, mutual information on the CT clipped to [-1000, 2000] HU, brain mask dilated by 10 mm so that the skull counts), to a CT template (2v1) or to the study template; registration record read by the PET pipeline.
- `coreg2anat.py`: sidecars written under the wrong names, `FS_freeview.check_coreg` called with 5 and 6 arguments, no record.
- `reorient.py`: copy of the PET module, replaced by a call to it (same orientation for the PET and the CT of one PET/CT).

## 3. Engineer's suggestions outside `modalities/pet` and `modalities/ct`

| file | decision |
|---|---|
| `Tools/getpath.py`: `ct()` | **integrated**, with `pet()`, as dictionaries following the folders of `pet2b.svg` (`raw`, `orig/{preprocessing/forMC, matrices, labels, masks, postprocessing/ROI_analysis/pictures}`, `acpc`, `templates/<name>`, `QC`) instead of 22-element tuples copied from `func` (`acpc-func` / `acpc-anat`). Only the PET / CT code uses them. |
| `Tools/extract_filename.py`: `change`, `move`, `get_subname`, `get_runname`, `get_typename`, `get_endname`, `get_desc` | `get_desc` **integrated** (entity matched on its key, empty string when absent); the others **not integrated**: unused by the corrected code and fragile (`get_runname` assumes three entities; `change` / `move` return `None` when the file does not exist yet). |
| `Tools/get_orientation.py` | **not integrated**: it is the old version with lookup tables, whose 5 left-right mirrors were corrected in September 2026 (geometric mapping, tested), and `OPPOSITES.get(humanPosition, animalPosition, default)` raises a TypeError. |
| `Tools/roiStats.py` | **rewritten** (numpy), see above. |
| `modalities/anat/freeview/FS_freeview.py` | **partly integrated**, in `modalities/anat/freesurfer/FS_freeview.py` (the `freeview` folder does not exist): `check_coreg` added as a top-level function. Not taken: the renaming `FSview` -> `check_seg` (breaks `_loop3.py`) and the nesting of `check_coreg` inside `check_seg` (indentation; `NameError` when `space=''`). |
| `atlases/getlabels.py` | **rewritten**, interactive windows kept: tkinter imported only when a window opens (the module imports on a server), EDNiX label files (`label_code/<atlas>_label.txt`), cerebellum table corrected (missing commas in `['BAL2013','BAL2020''BAL2023''BAL2024']`, `REF` instead of `ref` for INIA19), non-interactive resolution of `reference_signal`, saved choice. |

## 4. Interface, documentation, launchers, tests

- Interface: **CT** tab; steps of `preprocess_pet` and `preprocess_ct` in the `Skip_step` checklist; choices for the new options; the run progress reads the PET / CT step banners (`Working on step N (function: pet_...)`); PET / CT diaries in the log review; `pet/QC` and `ct/QC` in the quality-control tab.
- Documentation: `docs/pet_pipeline.md` (methods, scenarios, inputs, outputs, interactive steps, QC, limitations, references), generated references `docs/reference/pet_parameters.md` and `ct_parameters.md`, updates of README, launchers, BIDS, installation, QC and tests pages. The three diagrams are in `Illustration/`.
- Example launchers: `Launcher_pet/Launcher_fdg_Macaque_imagina.py` (BIDS_Tremblay, 1v2a; `orientation='RPS'` kept from the previous launcher, to check on one subject), `Launcher_pet/Launcher_fdg_CERMEP.py` (Human, 1v2a by default, 1v1 / 1v2b explained), `Launcher_ct/Launcher_ct_CERMEP.py`.
- Tests: `tests/unit/test_pet.py` (22 tests, analytic head phantom, AFNI and Turku tools emulated, ANTs real): helpers, SUV and decay, DFT and reference TAC, motion correction (inter-frame correlation 0.97 -> > 0.995 for motions of 3-5 mm / 4-6 degrees), PET -> anatomy registration (pose within half a PET voxel), chains, export alignment, and runs of `preprocess_pet` / `preprocess_ct` in the scenarios 1v1, 1v2a, 1v2b, 2v1, 3v1, 3v2. The PET and CT modules are now in the import test; the PET launchers in the launcher checks.
- Unchanged, already failing before this work (not related to the PET): 21 tests of `test_post_ttt.py` (TypeError in `_7_post_TTT`), `Launcher_func_Dog.py` (parameter `extract_exterior_CSF` unknown to `preprocess_fmri`), `test_qc_scan_table_and_path_safety`.

## 5. To check on the first real data

1. **1v1** (hybrid PET-MR): the PET follows the reorientation of the anatomy as the functional images do; check the QC image of the first subject (assumption: the AFNI header operations move images of different grids identically).
2. **Turku tools**: one parametric image (BP) against the anatomy (geometry written back from the input); the SIF date is 01/01/1970 without DICOM (display only).
3. **BIDS_Tremblay**: `orientation='RPS'` (from the previous PET launcher; the anatomy of this dataset uses `'LSP'`): run step 1 with `Need_check_orient=True` on one subject.
4. **CERMEP**: folder and name of the CT (`ct/*_ct.nii.gz` assumed), units of the PET (`Units`: Bq/mL or SUV) and presence of `InjectedRadioactivity` / `body_weight`; with the NIfTI version coregistered to the PET by the authors, `anat_pet_same_space=True` is exact.
5. **Reference regions** of the receptor tracers (RACLO, PE2I, BUPRE...): validity per species (literature).
