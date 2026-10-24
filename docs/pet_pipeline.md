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
| **1v1** | anatomical MRI, same session and scanner space (hybrid PET-MR) | world transform of the reorientation of the anatomy applied to the PET, head position of the anatomy acquisition, AC-PC matrix of the anatomy, then anatomy -> template | yes (atlases of the anatomy) |
| **1v2a** | anatomical MRI, other session | rigid PET -> anatomy (mutual information, dilated brain mask), then anatomy -> template | yes |
| **1v2b** | anatomical MRI and CT (PET/CT scanner) | PET -> CT (same space: world transform of the CT and head position of the CT acquisition; otherwise rigid, head mask) -> anatomy (CT pipeline) -> template | yes |
| **2v1** | CT and a CT template, no MRI | PET -> CT -> CT template (`SyNRA`: rigid + affine + SyN) | yes (template atlases) |
| **2v2** | CT, no CT template | PET -> CT -> CT study template; group step `bridge`: CT study template -> PET study template -> reference template (`ct_template_bridge=True`) | yes (template atlases, through the bridge) |
| **3v1** | PET template in the space of the reference template | PET -> PET template (`SyNRA`) | yes (template atlases) |
| **3v2** | nothing else | PET -> PET study template | no |

The scenario is chosen from `IhaveanANAT`, `anat_pet_same_space`, `IhaveaCT` (and the record of the CT pipeline), `template_PET` and `creat_study_template_pet`; it is written in the diary and in the registration record. In 3v2 (and in 2v2 with `ct_template_bridge=False`) the study template is not linked to the reference template: without an atlas or a brain mask, only the SUV images are computed, in the PET and study-template spaces. In 2v2 the group step `bridge` (below) provides that link.

With an anatomy, the AC-PC chain is composed with the anatomy -> template chains of `modalities/anat/loop3/backtonative.get`, the bridge used by the functional pipeline: reference template, study template and reference session (longitudinal), with the same `creat_study_template` and `coregistration_longitudinal` as the anatomical launcher. One registration is estimated per session, on the mean image of the session after the between-runs motion correction.

**Same scanner space (1v1; PET/CT with `ct_pet_same_space=True`).** Two conditions make the PET follow the anatomy (or the CT) without registration:

1. *Header.* `3drefit -orient` keeps the magnitude of every origin and sets its sign from the new orientation code (AFNI `3drefit.c`: `xxorg = (ORIENT_sign == '+') ? -|xxorg| : |xxorg|`). Two images of one scanner space with different fields of view are therefore relabelled consistently only when their first voxels lie on the same side of the coordinate origin; otherwise each is shifted by twice its own first-voxel coordinate (7 mm between a 72-mm and an 80-mm field of view in the tests, 20-40 mm with a table offset). The PET is not relabelled on its own: after the obliquity correction it receives the world transform W = A_after . A_before^-1 that the reorientation applied to the reference (raw and reoriented T1w of `type_norm` of the anatomical pipeline; the CT's own transform, recorded in its sidecar). The positions and `orientation` of the PET launcher are then not used. Keep `deoblique='WARP'` (as the anatomical pipeline).
2. *Head position.* The motion correction puts the session at the mean head position of the frames (and of the runs). The record brings it back to the position of the frames acquired with the anatomy / CT: mean (`MoCo.mean_rigid`) of the motion transforms of the frames of `same_space_window` ([start, end], min after the injection; default: the first 10 min of the first run, CT acquired before the PET), composed with the between-runs transform of that run.

Check the first subject anyway (`Need_check_coreg=True`, or the QC image `pet/QC/<session>_space-acpc_pet.png`). When in doubt, use 1v2a, which also works for images already coregistered. The functional pipeline (`anat_func_same_space`) and the T1w / T2w pair of the anatomical pipeline are relabelled independently and have the same exposure to off-centre fields of view.

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
| 6 | registration (scenario), quality control images, atlases / masks / tissue fractions / anatomy in the PET space |
| `bridge` | 2v2 only, once for the group after step 6 of every session: PET study template in the CT study template space, reference template -> PET study template, reference target and atlases added to every session |
| 7 | parametric images and regional values |
| 8 | reference-region kinetic models (`srtm=True`) |
| `pvc` | partial-volume correction (`pvc_fwhm`) |
| 9 | PET during the fMRI runs (`do_pet_to_func=True`, simultaneous PET-MR) |
| 10 | export to the AC-PC and template spaces |

| step | CT pipeline |
|---|---|
| 1 | reorientation of the header (same positions as the PET of the PET/CT) |
| 2 | CT study template (scenario 2v2), built on the CTs restricted to the head (`desc-head_ct`: > -500 HU, largest component, outside set to -1000 HU, cropped) |
| 3 | registration, record for the PET pipeline, CT in the AC-PC / template space, anatomy in the CT space, QC |

Every parameter is described in [reference/pet_parameters.md](reference/pet_parameters.md) and [reference/ct_parameters.md](reference/ct_parameters.md).

## Inputs

BIDS-PET layout: `sub-<ID>/ses-<ses>/pet/sub-<ID>_ses-<ses>_trc-<tracer>[_run-<n>]_pet.nii.gz` and its sidecar (`endpet` selects the files). Fields read from the sidecar:

| field | use |
|---|---|
| `FrameTimesStart`, `FrameDuration` (s) | frame timing (required for a dynamic image; overlapping frames refused); `FrameReferenceTime` used for the mid-frame times when it lies within its frame |
| `InjectionStart` (s, relative to `TimeZero`) | origin of all the times of the pipeline: frame times are converted to times after the injection (BIDS allows `TimeZero` = scan start); `SUV_window`, `duration`, `DV_time`, `kinetic_end`, `same_space_window` and the DFT / SIF times are post-injection |
| `Units` | conversion to Bq/mL for the SUV (`Bq/mL`, `kBq/mL`, `kBq/cc`, `kBq/cm3`, `nCi/mL`, DICOM `BQML`...); an image in `SUV` / `g/mL` is kept as it is |
| `InjectedRadioactivity`, `InjectedRadioactivityUnits` | SUV |
| `ImageDecayCorrected`, `ImageDecayCorrectionTime` | decay reference of the SUV; non-decay-corrected data are decay-corrected at step 7 (also used by the kinetic tools) |
| `TracerRadionuclide` | isotope (else `isotope`, else the tracer table); `F18`, `18F`, `[18F]`, `Fluorine-18` are recognised, an unrecognised `isotope` stops the run |
| `TimeZero` | SIF title line (injection clock time); PET-fMRI timing |
| `PatientPosition`, `PatientWeight` | `humanPosition` when left empty; body weight when not in the tsv files |
| `ReconFilterType`, `ReconFilterSize` | reported with the partial-volume correction (`pvc_fwhm` must include the filter) |

Body weight and blood glucose: `sub-<ID>/sub-<ID>_sessions.tsv`, then `participants.tsv` (columns `body_weight`, `weight` or `weight_kg` in kg; `blood_glucose` in mmol/L or `blood_concentration` in g/L; the injected radioactivity can also be given as `injected_radioactivity` or `tracor_qte_MBq`, in MBq). The values found and their source are written in `pet/Prepro_Settings.json`.

CT: `sub-<ID>/ses-<ses>/ct/*_ct.nii.gz` (`endct`), in Hounsfield units.

## Methods

**Reorientation (step 1).** As the anatomy: an oblique image is resampled on a cardinal grid with linear interpolation (`deoblique='WARP'`; no overshoot of the activity concentrations), the axes are relabelled with `orientation` or `get_orientation.getreal(humanPosition, animalPosition)`; in the same scanner space as the anatomy / CT, the PET receives their world transform instead (see above). `PET_header_wrong=True` (wrong orientation written by the reconstruction) copies the orientation code of the raw anatomy or CT of the same scanner space before this reorientation. A static image stored as 4D with one frame is written as 3D.

**Motion correction within a run (step 4).** Rigid registration (ANTs, mutual information) of every frame to the duration-weighted mean of the frames from `nb_void` on, in two passes: the second registers to the mean of the frames corrected by the first, a mean of misaligned frames being a blurred reference (as SPM Realign or AFNI 3dvolreg with a mean base). The metric is computed in a head mask of the reference (voxels above 10 % of its 99th percentile, dilated by 2 voxels), at three levels ending at full resolution with every voxel sampled (`MoCo.MOCO_REG`). Without the mask the background dominates the mutual-information histogram: on a still phantom the frames received 0.4-2 mm of spurious motion, and the relative motion between frames was recovered within 1.5-2 mm (mean; max 3-4 mm) instead of 0.2-0.4 mm (max < 1 mm) with it (`tests/unit/test_pet.py::test_within_run_motion_estimates_are_accurate`). The registration uses the frames cropped to the head box (crop drawn once on the study template and brought to every session, or drawn on the session; `crop_for_MoCo`), the transform is applied to the full frame with linear interpolation. The first `nb_void` frames, too low in counts and with another distribution of the tracer, receive the transform of the first registered frame (nearest in time). Motion parameters and the mean displacement of a sphere of the radius of the head are written in `orig/preprocessing/forMC/*_motion.tsv` and `.png` (copied to `pet/QC`). `MoCo_transform` must be rigid (`Rigid`, `QuickRigid`, `DenseRigid`, `BOLDRigid`, `Translation`).

**Between runs (step 5).** Every run (its mean) is registered rigidly to the mean of the runs (two passes, same settings) and resampled on the grid of the first run; for a dynamic run the per-frame and between-runs transforms are composed, so that the reoriented frames are resampled only once. The session mean is then the reference of the registration. When step 4 is skipped and was never run, the frames are used without within-run correction (warning).

**Registrations (step 6).** PET -> anatomy (1v2a), PET -> CT (other space) and CT -> anatomy go through `coreg2anat.rigid_to`:

- *Multi-start initialisation* (`coreg_init='auto'`). ANTs starts a registration from the centres of mass of the two whole images; with the torso or the heart in the PET field of view, or with a CT in HU (the air, -1000 HU, is a negative mass larger than the head), that start is centimetres away and mutual information converges to a wrong optimum (85 mm error on the phantom with the torso, `test_multistart_registration_of_a_pet_with_the_torso_in_the_field_of_view`). Five starts are tried: centres of mass of the masked images (moments; Alpert et al. 1990, AFNI `-cmass`), centres of the two fields of view ('FOV to FOV'), the AC-PC matrix of the anatomy (the PET or CT is assumed near the native anatomy), the headers (identity) and a search over the rotations (`search`: ANTs `antsAffineInitializer`, rotations by 15 degrees within +/- 90 degrees on images resampled to 3 voxels of the anatomy, reduced to its rotation; for a head flexed on the CT head holder or an animal in another position: an 80-degree pitch is missed by the other starts and recovered within 0.5 mm on the phantom, `test_ct_large_rotation_search_and_brain_window_pass`). Each is followed by a translation then the rigid registration; the start kept has the highest Mattes mutual information (Wells et al. 1996; Maes et al. 1997) among those whose moving image covers 90 % of the brain (MI can be increased by sliding the object out of the mask). Multi-start searches against the local optima of MI are standard (Jenkinson & Smith 2001, FLIRT; AFNI `3dAllineate -twobest`). Scores: `orig/matrices/*candidates.json` and the diary. The CT is shifted to >= 0 and restricted to its head (> -500 HU, largest component: the bed is removed); mutual information is the reference metric for CT-MRI (West et al. 1997).
- *CT brain window* (`ct_brain_refine`, CT pipeline, off by default). On [-1000, 2000] HU the parenchyma (CSF ~5, white matter ~30, grey matter ~40 HU) falls in one or two of the 32 histogram bins of the Mattes metric: the skull and the scalp drive the CT -> T1w registration. The optional second pass registers the CT in the brain window (smoothed 1 mm against the low-dose noise, 0-80 HU mapped to 0-1, bone above 200 HU set to 0 as the cortical bone of a T1w) within the brain of the anatomy brought into the CT and dilated by 3 mm, from the head registration (brain window of head CT: Muschelli et al. 2015). Neutral on the phantoms (0.1 mm before and after); to try when the QC shows the brain shifted inside a well-aligned skull.
- *Before registering a CT* (geometry that no rigid transform corrects): a CT acquired with a gantry tilt must be resampled at the DICOM conversion (dcm2niix does it, Li et al. 2016; `GantryDetectorTilt` in the sidecar), otherwise the image is sheared; overlay the reoriented CT and the T1w (`Need_check_orient=True`): a flipped axis or a rotation beyond the search range comes from the headers, not from the metric. Rigid CT-MRI registration by mutual information is accurate to about 1-2 mm in the human head (West et al. 1997; Studholme et al. 1999).
- *Mask of the PET* (`pet_coreg_mask`, as `3dAllineate -source_mask` in the previous PET scripts). `'threshold'`: voxels above `pet_mask_percentile` of the positive voxels, largest connected component, holes filled, dilated `pet_mask_dilate_mm`. `'anat'`: after a first registration, the brain mask of the anatomy is brought into the PET, dilated `pet_mask_dilate_mm`, and the registration is estimated again on the PET brain only, within the same dilated brain of the anatomy: the extracerebral uptake (temporal muscles, Harderian / lacrimal glands, eye muscles, nasal mucosa: strong in FDG, large relative to the brain in macaques and rodents) no longer enters the metric. A mask file in the PET space can also be given.
- *Target* (`coreg_target`). `'anat'` (AC-PC head image, default), `'brain'` (skull-stripped), or `'synthetic'`: an FDG-like image of the anatomy (grey matter 1, white matter 0.35, CSF 0.05, outside 0; grey/white ratio about 3, Phelps et al. 1979) smoothed to the PET resolution, the segmentation-based PET-MRI registration of Ashburner & Friston (1997). FDG only.
- `coreg_transform='Rigid'`: a PET and an MRI of one subject differ by a rigid motion; `Similarity` / `Affine` can absorb a misregistration with a scale (PET brain shrunk into the MRI brain) and non-linear transforms are not physiological for one subject (warnings).

On the phantoms the four starts, the masks and the synthetic target all reach < 0.4 mm when the start is good; the phantoms do not reproduce the T1 contrast of the scalp nor the partial volume of a 2-mm macaque cortex, so the choice between `'anat'` masks and the `'synthetic'` target must be made on the QC of real data. A cortex whose FDG signal spreads symmetrically into the CSF and the white matter is partial volume (PET FWHM larger than the cortical thickness: see the `pvc` step), not misregistration; a misregistration shifts the signal into the CSF on one side and into the white matter on the other.

**Group step `bridge` (scenario 2v2, `ct_template_bridge=True`).** A CT has almost no contrast inside the brain, so the CT study template cannot be registered to a T1 template; the FDG grey / white matter contrast can (the chain of the dog notebook of the previous PET scripts, `modalities/pet/bridge.py`). The pipeline runs steps 1-6 for every session, then: (1) PET study template = mean of the session means brought to the CT study template through their own PET -> CT -> CT study template chains, each scaled by its 99th percentile (dose). A mean and not `ants.build_template`, whose additional warps would not be part of the session chains (the notebook built the PET template again with its own warps and did not propagate them to the subjects). (2) Skull-stripped reference template (moving, its dilated brain mask as moving mask) -> PET study template: multi-start similarity (`rigid_to`, same starts as above, no coverage condition since the T1 brain covers only part of the head-wide PET template), then affine initialised by it; `bridge_transform='SyNRA'` adds a SyN stage (Avants et al. 2008) in the dilated FDG mask of the PET template. An affine started from a translation stopped in a sheared local optimum on the phantom (6 mm error, against 1.2 mm with the cascade, the Translation / Rigid / Similarity / Affine order of ANTs `TRSAA`). (3) Brain mask of the reference template in the PET study template, ratio image, PET study template in the reference space (QC). (4) Every 2v2 session receives the reference target (session -> CT study template -> reference) and the atlases / masks of the reference template in its PET space; steps 7-10 then run as in 2v1. Record: `<BIDS>/sty_template/pet/matrices/studyTemplate_to_<reference>_desc-bridge_pet.json`; `Skip_step=['bridge']` reuses it. On the phantom: session -> reference chain within 2 mm of the truth and Dice > 0.7 of a template atlas region in the PET space (`test_ct_study_template_without_anatomy_2v2`). The bridge is a between-subject registration of templates of different modalities: its accuracy is that of the T1 / FDG-mean contrast, check the QC image first.

PET -> PET template and CT -> CT template: whole head, `SyNRA` by default (rigid + affine + SyN). ANTsPy's `SyN` and `Affine` start with a single affine stage, which missed a 12-degree rotation by about 7 mm on the test phantom (< 1 mm with `SyNRA` or the linear `TRSAA`), and a brain mask on these whole-head registrations increased the error (1.4-4.6 mm). The files of a previous registration with the same prefix are removed first (a warp left by an earlier non-linear run would otherwise be reused).

**Atlases and masks in the PET space (step 6).** Atlases: generic-label interpolation (every level of a 4D atlas); masks: linear interpolation thresholded at 0.5 (voxels at least half covered); tissue fractions (`_label-GM|WM|brain_probseg`, for the partial-volume correction): whole-brain grey matter, non-eroded white matter and brain mask averaged over the extent of a PET voxel (Gaussian of the variance of a box of the PET voxel size) before the linear resampling. With an anatomy, the atlases are those of the subject (`anat/native/volumes/labels`, including manual corrections); without, those of the template.

**Parametric images (step 7).** For a dynamic run, the duration-weighted mean of the frames whose middle falls in `SUV_window` (minutes after the injection) gives the activity image (`desc-activity-mean`; the dynamic SUV is also written); a static image is written as `desc-activity` (exported, and the only result without SUV data and without mask).

- SUV (body weight, g/mL) = C (Bq/mL) x body weight (g) / injected radioactivity (Bq), with a density of 1 g/mL (Boellaard 2009; Boellaard et al. 2015; Kinahan & Fletcher 2010). The injected radioactivity is decayed to `ImageDecayCorrectionTime` when it differs from `InjectionStart` (the SUV is not computed when the isotope is then unknown); non-decay-corrected data, static or dynamic, are decay-corrected to the injection frame by frame (exact frame-averaged factor λΔt / (1 - e^-λΔt) · e^λ(t_start - t_inj)). Without this information the image is assumed decay-corrected to the injection (written in the diary). Tested analytically: a static, non-decay-corrected frame 50-60 min after the injection gives the SUV of the white matter within 2 %.
- SUVglu = SUV x glucose / 5.0 mmol/L, when the blood glucose is known: SUV normalised to a reference glycaemia, whose value differs between protocols (5.0 mmol/L here, 100 mg/dL = 5.55 mmol/L in others; constant `GLUCOSE_REF_MMOL_L` of `paramImg.py`). The glucose correction is optional in the clinical guidelines (Boellaard et al. 2015).
- bratio = image / its mean over the brain mask (global normalisation).
- Time blocks (`Do_cut`): duration-weighted mean of the frames of every block of `duration` minutes from the first frame (frame assigned to the block containing its middle: the label `period10to20min` is the block, the frames may start before it), divided by its brain mean.
- Skull-stripped versions (`-SS`) and regional tables (mean, SD, median, voxels, volume) for every atlas and level; regional time-activity curves of the dynamic runs, with figures.

The regional values are computed on PET voxels: they include partial-volume effects (spatial resolution of the scanner of the order of the size of small structures, especially in small brains). Use `pvc_fwhm` to correct them, and interpret small regions with caution.

**Kinetic models (step 8, `srtm=True`).** Reference-region models with the Turku PET Centre tools (tpcclib, container `tpcclib.sif`):

- regional BPnd: simplified reference tissue model (SRTM, Lammertsma & Hume 1996) fitted with basis functions (`bfmsrtm`, Gunn et al. 1997), frame weights from the counts estimated by `imgweigh` / `dftweigh` (Mazoyer et al. 1986);
- voxel-wise BPnd: `imgbfbp` (basis functions; `desc-BP`, and `desc-BP-err` where theta3 reached a limit) and `imgsrtm` (`desc-SRTM`);
- voxel-wise DVR: Logan plot with a reference region (`imgdv`, Logan et al. 1996) from `DV_time` minutes after the injection, with the population k2' of the reference region when `DV_k2ref` is given (`-k2=`), otherwise the simplified form without the C_ref/k2' term (valid when C_ref / C_T is constant after t*).

The frame times given to the tools are post-injection (`InjectionStart`). The voxel-wise tools need decay-corrected data: data that are not are decay-corrected at step 7 and that image is used; without it the voxel-wise tools are not run (`bfmsrtm` then receives the isotope, `-i=`). The regional and reference TACs are computed on the motion-corrected image (not skull-stripped: label voxels outside the PET-space brain mask would count as zeros and dilute a TAC, BPnd' = (a/s)(BPnd + 1) - 1), the voxel-wise models on its skull-stripped version. The reference TAC is the mean concentration over all the voxels of the reference region (atlas labels, or a tissue mask: `white`, `csf`, `white+csf`, `gray`, `brain`); NaN voxels are ignored. Prefer label numbers to a name pattern for the reference region (`'Cerebell'` also matches the cerebellar white matter of some atlases). The reference region is given in `reference_signal`, or chosen in a window (`'unknown'`) and saved in `<BIDS>/sty_template/pet/reference_region_trc-<tracer>.json` for the other subjects. The validity of a reference region (no specific binding, same non-displaceable uptake as the target regions) depends on the tracer and the species: check it in the literature for each tracer (e.g. cerebellar grey matter for [11C]raclopride in primates).

**Partial-volume correction (`pvc`).** PETPVC (Thomas et al. 2016) of the first static result (SUV-mean, SUV, activity-mean or activity) with four tissue classes summing to one in every voxel: whole-brain grey matter, white matter outside it, CSF = rest of the brain mask, and the rest, from the tissue fractions of step 6 (the eroded white-matter and ventricle masks used for the reference regions would leave tissue in the 'rest' class). `pvc_fwhm` is the spatial resolution of the reconstructed image (point-spread function of the scanner and reconstruction combined with the post-reconstruction filter, FWHM in mm; sqrt(a^2 + b^2) for Gaussian kernels), not the voxel size. `pvc_method`: a PETPVC method that writes an image (`RBV`, `MG`, `LABBE+RBV`, `IY`...; GTM writes a table and is refused). Needs `nipype` and the `petpvc` program.

**PET during fMRI (step 9).** For a simultaneous PET-MR session, the mean PET during every fMRI run (frames weighted by their overlap with the run, from the clock times `TimeZero` of the PET and `AcquisitionTime` of the BOLD, minus the `T1_eq` first volumes; a run after midnight is handled), divided by its brain mean.

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
        masks/           <session>_space-orig[_desc-gray|white|csf]_mask.nii.gz, <session>_space-orig_label-GM|WM|brain_probseg.nii.gz
        postprocessing/  _desc-activity[-mean], -SUV[-mean], -SUVglu, -bratio[-mean], -bratio-period..., -BP, -SRTM, -DVR (+ -SS)
            ROI_analysis/    _seg-<atlas>-lvl<k>_vol_ROIs.tsv, _desc-<desc>_ROIs.tsv, _desc-Kinetic_ROIs.tsv / .dft, reference .dft, BPnd .res
                pictures/    time-activity curves, fits
    acpc/                    parametric images in the AC-PC space (+ -SS)
    templates/<name>/        parametric images in each template space (+ -SS), labels/, masks/
    QC/                      PET over the brain mask of every target space (head of a study template), motion plots
sub-<ID>/ses-<ses>/ct/       same layout (orig/preprocessing, orig/matrices/<base>_desc-registration_ct.json, acpc, templates, QC)
<BIDS>/sty_template/pet/     study template, its crop, reference region
<BIDS>/QC/pet_in_acpc/, pet_in_template/, ct_in_acpc/ ...   QC images of all subjects
```

Group analyses of these outputs (voxel-wise in a template space, regional, surfaces): [group_analysis.md](group_analysis.md#pet-group-analyses).

The registration record (`orig/matrices/*_desc-registration_*.json`) lists, for every target space, its reference image, mask and the two chains of transforms (to and from the PET or CT space): it is what the later steps, the CT -> PET link and group analyses use.

## Interactive steps

As in the engineer's version: freeview for the orientation (`Need_check_orient`), the motion correction (`Need_check_MC`), the registrations (`Need_check_coreg`), the manual correction of the PET -> anatomy registration (`do_manual`: move the PET with *Tools > Transform Volume*, save it under the same name, the registration is estimated again from it), the crops (`crop_for_MoCo`, saved and reused; `redo_crop=True` to draw them again), and a window for the reference region (`reference_signal='unknown'`, saved). For unattended runs (interface in the background, server), use `crop_for_MoCo='none'`, give `reference_signal`, and leave the checks off.

## Quality control

- `pet/QC/<session>_space-<target>_pet.png` (step 6): PET in each target space with the brain mask of that space as contour, or the head of a study template (also in `<BIDS>/QC/pet_in_acpc` and `pet_in_template`, shown by the interface); the same for the CT (`ct/QC`, `ct_in_*`);
- `pet/QC/*_motion.png` (step 4, also in `orig/preprocessing/forMC`): motion of every frame;
- diary: a warning when the session mean is newer than the registration record (steps 4-5 run again without step 6);
- `ROI_analysis/pictures`: time-activity curves and kinetic fits;
- diary: scenario, isotope, subject values and sources of the SUV, warnings (assumed decay correction, units).

## Limitations, to check on the first datasets

- Same scanner space: the anatomy / CT must not have been moved by hand after its reorientation, and `same_space_window` must match its acquisition.
- Turku tools: the outputs are written back with the geometry of the input image (voxel order restored when the tool writes a valid header); check one parametric image against the anatomy.
- The SIF date is 01/01/1970 without DICOM (BIDS has no date); tpcclib only displays it.
- Scenario 3v2: no link to the reference template, no atlas. Scenario 2v2: the bridge is a linear (or `SyNRA`) registration of a T1 template on an FDG group mean; check `<study template>_space-<reference>_pet.png` before using the atlases.
- MoCo of very early frames: they inherit the transform of the first registered frame (motion before that frame is not corrected).
- `imgdv` / `bfmsrtm` unit handling: the DFT carries the image units (`Units`, Bq/mL when absent) and the NIfTI image none; tpcclib then assumes both are in the same units (they are: the TACs come from the image).
- Registration accuracies above are on noise-free or Gaussian-noise phantoms with 2-mm PET voxels; real data (lower counts, tracer kinetics changing the contrast between frames) will do worse: check the motion plots and QC images.

## References

- Alpert NM, Bradshaw JF, Kennedy D, Correia JA (1990) The principal axes transformation: a method for image registration. J Nucl Med 31:1717-1722.
- Ashburner J, Friston KJ (1997) Multimodal image coregistration and partitioning: a unified framework. NeuroImage 6:209-217.
- Avants BB, Epstein CL, Grossman M, Gee JC (2008) Symmetric diffeomorphic image registration with cross-correlation: evaluating automated labeling of elderly and neurodegenerative brain. Med Image Anal 12:26-41.
- Boellaard R (2009) Standards for PET image acquisition and quantitative data analysis. J Nucl Med 50 Suppl 1:11S-20S.
- Boellaard R et al. (2015) FDG PET/CT: EANM procedure guidelines for tumour imaging: version 2.0. Eur J Nucl Med Mol Imaging 42:328-354.
- Gunn RN, Lammertsma AA, Hume SP, Cunningham VJ (1997) Parametric imaging of ligand-receptor binding in PET using a simplified reference region model. NeuroImage 6:279-287.
- Jenkinson M, Smith S (2001) A global optimisation method for robust affine registration of brain images. Med Image Anal 5:143-156.
- Kinahan PE, Fletcher JW (2010) Positron emission tomography-computed tomography standardized uptake values in clinical practice and assessing response to therapy. Semin Ultrasound CT MR 31:496-505.
- Li X, Morgan PS, Ashburner J, Smith J, Rorden C (2016) The first step for neuroimaging data analysis: DICOM to NIfTI conversion. J Neurosci Methods 264:47-56.
- Lammertsma AA, Hume SP (1996) Simplified reference tissue model for PET receptor studies. NeuroImage 4:153-158.
- Logan J, Fowler JS, Volkow ND, Wang GJ, Ding YS, Alexoff DL (1996) Distribution volume ratios without blood sampling from graphical analysis of PET data. J Cereb Blood Flow Metab 16:834-840.
- Maes F, Collignon A, Vandermeulen D, Marchal G, Suetens P (1997) Multimodality image registration by maximization of mutual information. IEEE Trans Med Imaging 16:187-198.
- Mazoyer BM, Huesman RH, Budinger TF, Knittel BL (1986) Dynamic PET data analysis. J Comput Assist Tomogr 10:645-653.
- Mérida I et al. (2021) CERMEP-IDB-MRXFDG: a database of 37 normal adult human brain [18F]FDG PET, T1 and FLAIR MRI, and CT images available for research. EJNMMI Res 11:91, doi:10.1186/s13550-021-00830-6.
- Muschelli J, Ullman NL, Mould WA, Vespa P, Hanley DF, Crainiceanu CM (2015) Validated automatic brain extraction of head CT images. NeuroImage 114:379-385.
- Norgaard M et al. (2022) PET-BIDS, an extension to the brain imaging data structure for positron emission tomography. Sci Data 9:65.
- Phelps ME, Huang SC, Hoffman EJ, Selin C, Sokoloff L, Kuhl DE (1979) Tomographic measurement of local cerebral glucose metabolic rate in humans with (F-18)2-fluoro-2-deoxy-D-glucose: validation of method. Ann Neurol 6:371-388.
- Studholme C, Hill DLG, Hawkes DJ (1999) An overlap invariant entropy measure of 3D medical image alignment. Pattern Recognit 32:71-86.
- Thomas BA et al. (2016) PETPVC: a toolbox for performing partial volume correction techniques in positron emission tomography. Phys Med Biol 61:7975-7993.
- Turku PET Centre, tpcclib documentation: https://www.turkupetcentre.net/petanalysis/tpcclib/doc/
- Wells WM, Viola P, Atsumi H, Nakajima S, Kikinis R (1996) Multi-modal volume registration by maximization of mutual information. Med Image Anal 1:35-51.
- West J et al. (1997) Comparison and evaluation of retrospective intermodality brain image registration techniques. J Comput Assist Tomogr 21:554-566.
