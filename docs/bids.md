# BIDS data and subject selection

EDNiX reads datasets organised in [BIDS](https://bids.neuroimaging.io/) (`sub-<ID>/ses-<session>/anat|func|fmap/...`) and writes its outputs next to the raw data of each subject and session.

## Reading the dataset

The launchers build a table of the files to process with `Tools/load_bids.Load_BIDS_to_pandas`, a light BIDS parser (no database):

```python
from Tools import load_bids, Load_subject_with_BIDS

allinfo_study_c = load_bids.Load_BIDS_to_pandas(
    bids_dir, modalities=['anat'], suffixes=['T1w'], extensions=['.nii.gz'])
Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)   # prints the (subject, session) pairs
```

Datasets without session folders are read with session `'1'`.

## Selecting subjects and sessions

`list_to_keep` and `list_to_remove` take `(subject, session)` pairs, without the `sub-` / `ses-` prefixes:

```python
list_to_keep = []                    # empty: every subject
list_to_remove = [('05', '2')]       # exclude subject 05, session 2
```

If `list_to_keep` is not empty, only its pairs are processed (and `list_to_remove` is not used). Identifiers are compared as text, so `('05', 2)` also matches; pairs that match nothing in the dataset are reported in the log.

`Tools/Load_subject_with_BIDS.load_data_bids` returns two groups of lists:

- `all_ID`, `all_Session`, `all_data_path`: every selected (subject, session);
- `all_ID_max`, `all_Session_max`, `all_data_path_max`: the **latest non-excluded session** of each subject, used to build the study template (`which_on='max'`) and as reference for longitudinal registration.

Sessions excluded with `list_to_remove` never enter these "max" lists. `list_to_keep` does not restrict them, because it is also used to process subjects in batches (`Tools/ednix_checkpoint.py`): the study template keeps using every non-excluded subject.

## Paths

`Tools/Load_subject_with_BIDS.linux_path` converts Windows paths (`C:\data\...`) into the Linux form used inside the containers (`/mnt/c/data/...`).

## Where outputs go

For each subject and session (`<BIDS>/sub-<ID>/ses-<session>/`):

```
anat/
    native/volumes/            images in the subject's AC-PC space (labels/, masks/)
    native/surfaces/           Workbench surfaces (Native_resol, fsaverage_LR_32k)
    templates/<reference>/     images and atlases in template space
    templates/studyTemplate/   (study template option)  templates/refSession/ (longitudinal option)
    freesurfer/                FreeSurfer subject folder
    matrices/                  transforms
    preprocessing/             intermediate files
func/
    raw/                       native functional space: preprocessing/, masks/, matrices/
    acpc-func/                 functional resolution, AC-PC space: postprocessed_rs/ (resting state)
    acpc-anat/                 anatomical resolution, AC-PC space
    templates/<reference>/     template space: postprocessed_rs/, labels/, masks/
fmap/preprocessing/            field-map processing (TOPUP / FUGUE)
```

Dataset-level quality-control images and tables are written in `<BIDS>/QC/` (see [quality_control.md](quality_control.md)).
