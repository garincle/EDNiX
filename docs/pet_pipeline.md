# PET and CT pipelines

The PET pipeline (`modalities/pet`, entry point `preprocess_pet`) processes static and dynamic PET images of a BIDS study. Its design follows three diagrams of the PET engineer (folder `Illustration/`):

- `pet4.svg`: how the PET reaches the template space depending on the images available (scenarios below);
- `pet2b.svg`: the output folders;
- `pet3b.svg`: the call tree of the modules.

**Principle.** Every computation is done in the space of the PET (`orig`): motion correction, parametric images, regional time-activity curves and kinetic models. The atlases and masks are brought *into* the PET space (one resampling of a label image), and only the final parametric images are resampled to the AC-PC space of the anatomy and to the template spaces. The dynamic data are therefore never resampled more than once after the motion correction, and the regional curves are computed on the measured voxels.

The CT pipeline (`modalities/ct`, entry point `preprocess_ct`) registers a CT to the anatomy or to a CT template; the PET pipeline then goes through the CT (`IhaveaCT=True`).

## Scenarios

| | Images | PET -> template | Atlas in PET space |
|---|---|---|---|
| **1v1** | anatomical MRI, same session and scanner space (hybrid PET-MR) | PET reoriented like the anatomy, then the AC-PC matrix of the anatomy, then anatomy -> template | yes (atlases of the anatomy) |
| **1v2a** | anatomical MRI, other session | rigid PET -> anatomy (mutual information, dilated brain mask), then anatomy -> template | yes |
| **1v2b** | anatomical MRI and CT (PET/CT scanner) | PET -> CT (same space: nothing) -> anatomy (CT pipeline) -> template | yes |
| **2v1** | CT and a CT template, no MRI | PET -> CT -> CT template (affine + SyN) | yes (template atlases) |
| **2v2** | CT, no CT template | PET -> CT -> CT study template | no |
| **3v1** | PET template in the space of the reference template | PET -> PET template (affine + SyN) | yes (template atlases) |
| **3v2** | nothing else | PET -> PET study template | no |

The scenario is chosen from `IhaveanANAT`, `anat_pet_same_space`, `IhaveaCT` (and the record of the CT pipeline), `template_PET` and `creat_study_template_pet`; it is written in the diary and in the registration record. In 2v2 and 3v2 the study template is not linked to the reference template ("possible in theory, not available yet" in `pet4.svg`): without an atlas or a brain mask, only the SUV images are computed, in the PET and study-template spaces.

With an anatomy, the AC-PC chain is composed with the anatomy -> template chains of `modalities/anat/loop3/backtonative.get`, the bridge used by the functional pipeline: reference template, study template and reference session (longitudinal), with the same `creat_study_template` and `coregistration_longitudinal` as the anatomical launcher. One registration is estimated per session, on the mean image of the session after the between-runs motion correction.

**1v1, same space.** As for functional images acquired in the space of the anatomy (`anat_func_same_space`), the PET header is reoriented with the same operations as the anatomy (3dWarp -deoblique, 3drefit -orient with `orientation` or the positions), and the PET is assumed to be in the space of the reoriented anatomy. This assumes that these header operations move two images of different grids in the same way; check the first subject (`Need_check_coreg=True`, or the QC image `pet/QC/<session>_space-acpc_pet.png`). When in doubt, use 1v2a, which also works for images already coregistered.

## Running

1. Anatomical pipeline of the dataset (except scenarios 2 and 3).
2. CT pipeline (`Exemples/Study_EDNiX/Launcher_ct`), if a CT is used.
3. PET pipeline (`Exemples/Study_EDNiX/Launcher_pet`).

```bash
EDNIX_STRICT=1 python Exemples/Study_EDNiX/Launcher_pet/Launcher_fdg_CERMEP.py
```

| step | PET pipeline |
|---|---|
| 1 | frame timing, SIF and mean image of the raw runs (`pet/raw`); reorientation of the header |
| 2 | PET study template (`creat_study_template_pet`), every session registered to it |
| 3 | crop drawn once on the study template in freeview |
| 4 | motion correction within each dynamic run |
| 5 | motion correction between the runs of the session; mean image of the session |
| 6 | registration (scenario), atlases / masks / anatomy in the PET space |
| 7 | parametric images and regional values |
| 8 | reference-region kinetic models (`srtm=True`) |
| `pvc` | partial-volume correction (`pvc_fwhm`) |
| 9 | PET during the fMRI runs (`do_pet_to_func=True`, simultaneous PET-MR) |
| 10 | export to the AC-PC and template spaces, quality control |

| step | CT pipeline |
|---|---|
| 1 | reorientation of the header (same positions as the PET of the PET/CT) |
| 2 | CT study template (scenario 2v2) |
| 3 | registration, record for the PET pipeline, CT in the AC-PC / template space, anatomy in the CT space, QC |

Every parameter is described in [reference/pet_parameters.md](reference/pet_parameters.md) and [reference/ct_parameters.md](reference/ct_parameters.md).

## Inputs

BIDS-PET layout: `sub-<ID>/ses-<ses>/pet/sub-<ID>_ses-<ses>_trc-<tracer>[_run-<n>]_pet.nii.gz` and its sidecar (`endpet` selects the files). Fields read from the sidecar:

| field | use |
|---|---|
| `FrameTimesStart`, `FrameDuration` (s) | frame timing (required for a dynamic image); `FrameReferenceTime` used for the mid-frame times when present |
| `Units` | conversion to Bq/mL for the SUV (`Bq/mL`, `kBq/mL`, `kBq/cc`, `nCi/mL`...); an image in `SUV` / `g/mL` is kept as it is |
| `InjectedRadioactivity`, `InjectedRadioactivityUnits` | SUV |
| `ImageDecayCorrected`, `ImageDecayCorrectionTime`, `InjectionStart` | decay reference of the SUV; non-decay-corrected data for the kinetic tools |
| `TracerRadionuclide` | isotope (else `isotope`, else the tracer table) |
| `TimeZero`, `ScanStart` | SIF title line; PET-fMRI timing |
| `PatientPosition` | `humanPosition` when left empty |

Body weight and blood glucose: `sub-<ID>/sub-<ID>_sessions.tsv`, then `participants.tsv` (columns `body_weight`, `weight` or `weight_kg` in kg; `blood_glucose` in mmol/L or `blood_concentration` in g/L; the injected radioactivity can also be given as `injected_radioactivity` or `tracor_qte_MBq`, in MBq). The values found and their source are written in `pet/Prepro_Settings.json`.

CT: `sub-<ID>/ses-<ses>/ct/*_ct.nii.gz` (`endct`), in Hounsfield units.

## Methods

**Reorientation (step 1).** As the anatomy: an oblique image is resampled on a cardinal grid with linear interpolation (`deoblique='WARP'`; no overshoot of the activity concentrations), the axes are relabelled with `orientation` or `get_orientation.getreal(humanPosition, animalPosition)`. `PET_header_wrong=True` (wrong orientation written by the reconstruction) copies the orientation code of the raw anatomy or CT of the same scanner space before this reorientation.

**Motion correction within a run (step 4).** Rigid registration (ANTs, mutual information) of every frame to the duration-weighted mean of the frames from `nb_void` on, in two passes: the second registers to the mean of the frames corrected by the first, a mean of misaligned frames being a blurred reference (as SPM Realign or AFNI 3dvolreg with a mean base). The registration uses the frames cropped to the head box (crop drawn once on the study template and brought to every session, or drawn on the session; `crop_for_MoCo`), the transform is applied to the full frame with linear interpolation. The first `nb_void` frames, too low in counts and with another distribution of the tracer, receive the transform of the first registered frame (nearest in time). Motion parameters and the mean displacement of a sphere of the radius of the head are written in `orig/preprocessing/forMC/*_motion.tsv` and `.png`.

**Between runs (step 5).** Every run (its mean) is registered rigidly to the mean of the runs and resampled on the grid of the first run; for a dynamic run the per-frame and between-runs transforms are composed, so that the reoriented frames are resampled only once. The session mean is then the reference of the registration.

**Atlases and masks in the PET space (step 6).** Atlases: generic-label interpolation (every level of a 4D atlas); masks: linear interpolation thresholded at 0.5 (voxels at least half covered). With an anatomy, the atlases are those of the subject (`anat/native/volumes/labels`, including manual corrections); without, those of the template.

**Parametric images (step 7).** For a dynamic run, the duration-weighted mean of the frames whose middle falls in `SUV_window` (minutes) gives the activity image (`desc-activity-mean`; the dynamic SUV is also written).

- SUV (body weight, g/mL) = C (Bq/mL) x body weight (g) / injected radioactivity (Bq), with a density of 1 g/mL (Boellaard 2009; Boellaard et al. 2015; Kinahan & Fletcher 2010). The injected radioactivity is decayed to `ImageDecayCorrectionTime` when it differs from `InjectionStart`; non-decay-corrected dynamic data are decay-corrected to the injection frame by frame (frame-averaged factor). Without this information the image is assumed decay-corrected to the injection (written in the diary).
- SUVglu = SUV x glucose / 5.0 mmol/L, when the blood glucose is known: SUV normalised to a reference glycaemia, whose value differs between protocols (5.0 mmol/L here, 100 mg/dL = 5.55 mmol/L in others; constant `GLUCOSE_REF_MMOL_L` of `paramImg.py`). The glucose correction is optional in the clinical guidelines (Boellaard et al. 2015).
- bratio = image / its mean over the brain mask (global normalisation).
- Time blocks (`Do_cut`): duration-weighted mean of the frames of every block of `duration` minutes (frame assigned to the block containing its middle), divided by its brain mean.
- Skull-stripped versions (`-SS`) and regional tables (mean, SD, median, voxels, volume) for every atlas and level; regional time-activity curves of the dynamic runs, with figures.

The regional values are computed on PET voxels: they include partial-volume effects (spatial resolution of the scanner of the order of the size of small structures, especially in small brains). Use `pvc_fwhm` to correct them, and interpret small regions with caution.

**Kinetic models (step 8, `srtm=True`).** Reference-region models with the Turku PET Centre tools (tpcclib, container `tpcclib.sif`):

- regional BPnd: simplified reference tissue model (SRTM, Lammertsma & Hume 1996) fitted with basis functions (`bfmsrtm`, Gunn et al. 1997), frame weights from the counts estimated by `imgweigh` / `dftweigh` (Mazoyer et al. 1986);
- voxel-wise BPnd: `imgbfbp` (basis functions; `desc-BP`, and `desc-BP-err` where theta3 reached a limit) and `imgsrtm` (`desc-SRTM`);
- voxel-wise DVR: Logan plot with a reference region (`imgdv`, Logan et al. 1996) from `DV_time` minutes.

The voxel-wise tools need decay-corrected data: they are not run otherwise (`bfmsrtm` then receives the isotope, `-i=`). The reference TAC is the mean concentration over all the voxels of the reference region (atlas labels, or a tissue mask: `white`, `csf`, `white+csf`, `gray`, `brain`). The reference region is given in `reference_signal`, or chosen in a window (`'unknown'`) and saved in `<BIDS>/sty_template/pet/reference_region_trc-<tracer>.json` for the other subjects. The validity of a reference region (no specific binding, same non-displaceable uptake as the target regions) depends on the tracer and the species: check it in the literature for each tracer (e.g. cerebellar grey matter for [11C]raclopride in primates).

**Partial-volume correction (`pvc`).** PETPVC (Thomas et al. 2016) with the grey matter, white matter and CSF masks in the PET space; `pvc_fwhm` is the spatial resolution of the scanner and reconstruction (point-spread function, FWHM in mm), not the voxel size. Needs `nipype` and the `petpvc` program.

**PET during fMRI (step 9).** For a simultaneous PET-MR session, the mean PET during every fMRI run (frames weighted by their overlap with the run, from the clock times `TimeZero` of the PET and `AcquisitionTime` of the BOLD, minus the `T1_eq` first volumes), divided by its brain mean.

**Export (step 10).** Every 3D image of `orig/postprocessing` is resampled once (chain of the record) to the AC-PC space and to each template space, with its skull-stripped version; `export_resolution='pet'` (default) resamples the target space at the PET voxel size, `'target'` keeps the grid of the anatomy / template. The atlases and masks of the template spaces are written on the same grid.

## Outputs

```
sub-<ID>/ses-<ses>/pet/
    <ID> session <ses><date>.txt         diary          Prepro_Settings.json
    raw/        <run>_desc-frametiming_pet.json, <run>_desc-mean_pet.nii.gz, <run>_pet.sif
    orig/
        preprocessing/   <run>_desc-reorient[-mean]_pet, _desc-MoCo-within[-mean]_pet, _desc-MoCo-between[-mean]_pet
                         <session>_desc-MoCo-between-mean_pet (session mean), <session>_space-orig_T1w (anatomy in PET space)
            forMC/       crop box, per-frame matrices (_mcNNN_0GenericAffine.mat), motion tables, PET in each target space
        matrices/        <session>_desc-registration_pet.json (record), <session>_from-pet-to-T1w_0GenericAffine.mat ...
        labels/          <session>_space-orig_seg-<atlas>_dseg.nii.gz
        masks/           <session>_space-orig[_desc-gray|white|csf]_mask.nii.gz
        postprocessing/  _desc-activity-mean, -SUV[-mean], -SUVglu, -bratio[-mean], -bratio-period..., -BP, -SRTM, -DVR (+ -SS)
            ROI_analysis/    _seg-<atlas>-lvl<k>_vol_ROIs.tsv, _desc-<desc>_ROIs.tsv, _desc-Kinetic_ROIs.tsv / .dft, reference .dft, BPnd .res
                pictures/    time-activity curves, fits
    acpc/                    parametric images in the AC-PC space (+ -SS)
    templates/<name>/        parametric images in each template space (+ -SS), labels/, masks/
    QC/                      PET over the brain mask of every target space
sub-<ID>/ses-<ses>/ct/       same layout (orig/preprocessing, orig/matrices/<base>_desc-registration_ct.json, acpc, templates, QC)
<BIDS>/sty_template/pet/     study template, its crop, reference region
<BIDS>/QC/pet_in_acpc/, pet_in_template/, ct_in_acpc/ ...   QC images of all subjects
```

The registration record (`orig/matrices/*_desc-registration_*.json`) lists, for every target space, its reference image, mask and the two chains of transforms (to and from the PET or CT space): it is what the later steps, the CT -> PET link and group analyses use.

## Interactive steps

As in the engineer's version: freeview for the orientation (`Need_check_orient`), the motion correction (`Need_check_MC`), the registrations (`Need_check_coreg`), the manual correction of the PET -> anatomy registration (`do_manual`: move the PET with *Tools > Transform Volume*, save it under the same name, the registration is estimated again from it), the crops (`crop_for_MoCo`, saved and reused; `redo_crop=True` to draw them again), and a window for the reference region (`reference_signal='unknown'`, saved). For unattended runs (interface in the background, server), use `crop_for_MoCo='none'`, give `reference_signal`, and leave the checks off.

## Quality control

- `pet/QC/<session>_space-<target>_pet.png`: PET in each target space with the brain mask of that space as contour (also in `<BIDS>/QC/pet_in_acpc` and `pet_in_template`, shown by the interface);
- `orig/preprocessing/forMC/*_motion.png`: motion of every frame;
- `ROI_analysis/pictures`: time-activity curves and kinetic fits;
- diary: scenario, isotope, subject values and sources of the SUV, warnings (assumed decay correction, units).

## Limitations, to check on the first datasets

- 1v1: consistency of the reorientation of two images of different grids (see above).
- Turku tools: the outputs are written back with the geometry of the input image (voxel order restored when the tool writes a valid header); check one parametric image against the anatomy.
- The SIF date is 01/01/1970 without DICOM (BIDS has no date); tpcclib only displays it.
- Scenarios 2v2 / 3v2: no link to the reference template, no atlas.
- MoCo of very early frames: they inherit the transform of the first registered frame (motion before that frame is not corrected).

## References

- Boellaard R (2009) Standards for PET image acquisition and quantitative data analysis. J Nucl Med 50 Suppl 1:11S-20S.
- Boellaard R et al. (2015) FDG PET/CT: EANM procedure guidelines for tumour imaging: version 2.0. Eur J Nucl Med Mol Imaging 42:328-354.
- Gunn RN, Lammertsma AA, Hume SP, Cunningham VJ (1997) Parametric imaging of ligand-receptor binding in PET using a simplified reference region model. NeuroImage 6:279-287.
- Kinahan PE, Fletcher JW (2010) Positron emission tomography-computed tomography standardized uptake values in clinical practice and assessing response to therapy. Semin Ultrasound CT MR 31:496-505.
- Lammertsma AA, Hume SP (1996) Simplified reference tissue model for PET receptor studies. NeuroImage 4:153-158.
- Logan J, Fowler JS, Volkow ND, Wang GJ, Ding YS, Alexoff DL (1996) Distribution volume ratios without blood sampling from graphical analysis of PET data. J Cereb Blood Flow Metab 16:834-840.
- Mazoyer BM, Huesman RH, Budinger TF, Knittel BL (1986) Dynamic PET data analysis. J Comput Assist Tomogr 10:645-653.
- Mérida I et al. (2021) CERMEP-IDB-MRXFDG: a database of 37 normal adult human brain [18F]FDG PET, T1 and FLAIR MRI, and CT images available for research. EJNMMI Res 11:91, doi:10.1186/s13550-021-00830-6.
- Norgaard M et al. (2022) PET-BIDS, an extension to the brain imaging data structure for positron emission tomography. Sci Data 9:65.
- Thomas BA et al. (2016) PETPVC: a toolbox for performing partial volume correction techniques in positron emission tomography. Phys Med Biol 61:7975-7993.
- Turku PET Centre, tpcclib documentation: https://www.turkupetcentre.net/petanalysis/tpcclib/doc/
