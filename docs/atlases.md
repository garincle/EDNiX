# Atlases and the atlas library

## Atlas library layout

The templates, atlases and reference surfaces of every species are in `Atlases_library/` inside the EDNiX folder:

```
Atlases_library/atlas/<group>/<subgroup>/<Species>/       e.g. atlas/mammals/primates/Marmoset
    label_code/
        <atlas>_label.txt          region names and colours (Workbench format)
        <atlas>_StatsLUT.txt       same, FreeSurfer colour-table format
        <atlas>_ITK.txt            same, ITK-SNAP format
        <atlas>_l.ctab, _r.ctab    FreeSurfer colour tables per hemisphere
    <reference>/                   e.g. EDNiX
        volumes/
            <Species>_space-acpc_desc-template_T1w.nii.gz   template (T1w, T2w)
            <Species>_space-acpc_desc-SS_T1w.nii.gz         skull-stripped template (T1w, T2w)
            labels/<Species>_seg-<atlas>_dseg.nii.gz        atlases (EDNIxCSC, EDNIxCSCLR, species atlases...)
            masks/<Species>_desc-<tissue>_mask.nii.gz       tissue masks (Cerebral_Gray, White, Vent,
                                                            Whole_Brain_Gray, eroded versions, cortical layers)
        surfaces/Native_resol/
            <Species>.<l|r>.<surface>.surf.gii              white, pial, midthickness, inflated, very_inflated, flat
            <Species>.<l|r>.<atlas>.native.label.gii        atlases on the surface
            <Species>.<atlas>.dlabel.nii                    CIFTI versions
        freesurfer/<Species>/                               FreeSurfer subject of the template
```

`species` in the launchers is the `<Species>` folder name, `reference` the `<reference>` folder name.

## EDNiX atlases

- **EDNIxCSC**: a cross-species parcellation with **4 levels** of detail, stored as the 4 volumes of one 4D image (level 0: coarsest, level 3: finest). Region IDs are unique across levels, so one label table describes every level.
- **EDNIxCSCLR**: the same regions split by hemisphere; the right-hemisphere ID is the left-hemisphere ID plus 1000.

Species atlases (for example a published macaque or marmoset parcellation) are listed in `atlas_followers` and processed in the same way.

## Label tables

Region names come from `label_code/<atlas>_label.txt`, a Workbench label file (one line with the name, one line with `ID R G B A`). `_StatsLUT.txt` (FreeSurfer: `ID name R G B A`) and `_ITK.txt` (ITK-SNAP: `ID R G B A visibility mesh "name"`) contain the same information. Names prefixed with `L_` / `R_` (or `lh`, `rh`) are read as left / right regions.

## Where atlases are used

- **Anatomy**: atlases are registered to every subject (`anat/native/volumes/labels/`), projected on its surfaces, and used for regional volumes, thickness and surface area.
- **Functional**: correlation matrices (`selected_atlases_matrix`, `wanted_level_matrix`) and seeds (`selected_atlases_SBA`, `panda_files_SBA`); `atlases/atlas4func.py` prepares them.
- **Interface**: the **Atlases** tab reads the same library: volume atlases (surfaces computed at the first display and cached in `~/.ednix_gui/cache`) and cortical surfaces with their labels, several species side by side, and seeds sent to the functional launcher.

## Adding an atlas

1. Put the atlas volume in `volumes/labels/` as `<Species>_seg-<atlas>_dseg.nii.gz`, in the space of the template.
2. Write its label table in `label_code/<atlas>_label.txt` (and the other formats if you use FreeSurfer or ITK-SNAP with it).
3. Add it to `atlas_followers` in your anatomical launcher; `atlases/templatefeat.py` finds it from its name.

## How the EDNiX atlases were built

The scripts that build the EDNiX atlases from the original species atlases are in [`Exemples/Study_EDNiX/atlases`](../Exemples/Study_EDNiX/atlases):

- `Atlas_nomenclature.py`: the region hierarchy of the EDNiX atlases (regions of each level and their parent regions), the reference for region names;
- `Transform_orig_atlas_to_EDNiX.py`: translation of the labels of an original atlas into EDNiX labels (a table of old and new IDs applied with AFNI `3dcalc`) and resampling to the template;
- `Update_EDNiX_labels.py`: label tables of the atlases (volume and, optionally, surface) for a list of species;
- `Change_EDNiX_spaceV2_to_acpc.py`: transfer of atlases between template spaces;
- `recontruct_surface_based_on_aseg.py`: reconstruction of template surfaces from a segmentation.

They document how the library was made; you do not need them to run the pipeline.

