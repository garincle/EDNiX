# Anatomical pipeline: parameters of `preprocess_anat`

Generated from `modalities/anat/_0_Pipeline_launcher.py` by `docs/generate_parameter_reference.py`; do not edit by hand.

## Steps (`Skip_step`)

Identifiers accepted by `Skip_step`: the listed steps are skipped.

| step | what it does |
|---|---|
| `1` | Orientation fix, averaging of repeated scans, bias-field correction |
| `2` | Cropping, denoising, AC-PC alignment and first skull-strip |
| `itk_1` | Manual check of the first skull-strip in ITK-SNAP |
| `5` | Individual brain template (second skull-strip) |
| `itk_2` | Manual check of the second skull-strip in ITK-SNAP |
| `6` | Registration to the study / reference template |
| `8` | Atlases and masks brought back to native space |
| `itk_3` | Manual check of the native-space atlases in ITK-SNAP |
| `9` | FreeSurfer preparation (masks, T2w to T1w) |
| `10` | FreeSurfer input images (rescaling of small brains) |
| `11` | FreeSurfer white surface |
| `12` | FreeSurfer pial surface |
| `flat_map` | Flat maps |
| `14` | Atlases projected on the surfaces |
| `15` | FreeSurfer to Workbench conversion |
| `16` | Anatomical quality control (SNR, CNR...) |
| `Clean` | Delete intermediate files |

## Parameters

| parameter | default | description |
|---|---|---|
| `Skip_step` | required | Steps of the pipeline to skip (integers or descriptive strings) |
| `MAIN_PATH` | required | Base path for EDNiX installation and resources |
| `bids_dir` | required | Path to BIDS-formatted dataset |
| `BIDStype` | required | Format identifier for BIDS dataset structure |
| `species` | required | Subject species ('Human', 'Rat', 'Mouse', etc. as in Atlases_library/atlas) |
| `allinfo_study_c` | required | DataFrame with subject/session metadata |
| `list_to_keep` | required | List of (subject, session) tuples to include |
| `list_to_remove` | required | List of (subject, session) tuples to exclude |
| `type_norm` | required | Anatomical contrast for normalization ('T1w' or 'T2w') |
| `otheranat` | required | Optional secondary anatomical contrast |
| `masking_img` | required | Image contrast for masking/skull-stripping |
| `orientation` | required | Target image orientation ('LPI', etc.) |
| `animalPosition` | required | Animal positioning during scanning |
| `humanPosition` | required | Human positioning |
| `coregistration_longitudinal` | required | Perform longitudinal coregistration |
| `creat_study_template` | required | Flag to create study-specific template |
| `which_on` | required | Template session selection: 'all' or 'max' |
| `brain_skullstrip_1` | required | Coarse skull-stripping method |
| `brain_skullstrip_2` | required | Fine skull-stripping method |
| `template_skullstrip` | required | Skull-stripping method for study/session template |
| `list_transfo` | required | Parameters for rigid and nonlinear transformations |
| `Align_img_to_template` | required | Coregistration tool/method |
| `MNIBcorrect_indiv` | required | Bias field correction method ('N4', 'N3', or '') |
| `fMRImasks` | required | Type of fMRI masks to generate ('aseg' or 'custom') |
| `reference` | `'EDNiX'` | Reference dataset or session |
| `do_fMRImasks` | `True` | Generate fMRI masks |
| `atlas_followers` | `[[], [], [], []]` | Multi-level atlas-following regions |
| `addatlas` | `''` | Additional atlas to include |
| `transfo_message` | `'do_as_I_said'` | Transformation control ('do_it_for_me' or 'do_as_I_said') |
| `force_myelin_same_space` | `False` | Maintain myelin maps in same space |
| `check_visualy_final_mask` | `False` | Open QC viewer for final mask |
| `check_visualy_each_img` | `False` | Open QC viewer for each anatomical image |
| `overwrite_option` | `True` | Overwrite existing files |
| `preftool` | `'ITK'` | Preferred QC visualization tool ('ITK' or 'freeview') |
| `doBALSAspace` | `False` | Macaque only: also build the volumetric atlas of the BALSA (Yerkes) surface parcellation in native space. |
