# Skull-stripping

Two orchestrators choose the method from its name: `modalities/anat/skullstrip/Skullstrip_method.py` for anatomical images (`brain_skullstrip_1`, `brain_skullstrip_2`, `template_skullstrip`) and `modalities/fMRI/Skullstrip_func.py` for the functional mean image (`Method_mask_func`). Every mask has a JSON sidecar.

## Anatomical methods

| family | names |
|---|---|
| AFNI `3dSkullStrip` | `3dSkullStrip`, `3dSkullStrip_noDilate`, `3dSkullStrip_macaque`, `3dSkullStrip_marmoset`, `3dSkullStrip_dog_macaque`, `3dSkullStrip_Rat`, `3dSkullStrip_rat`, `3dSkullStrip_rat_no_dil` |
| Intensity histogram (nilearn, AFNI clean-up) | `Custum_Macaque`, `Custum_Baboon`, `Custum_mouse`, `Custum_bat`, `Custum_dog`, `custum_rat`, `Custum_Chimp`, `CustumNilearn_<low>_<high>`, `CustumNilearnExcludeZeros_<low>_<high>`, `CustumThreshold_<percentile>` |
| Sammba histogram masks | `sammba_mouse`, `sammba_dog`, `sammba_rat`, `Vol_sammba_<volume>` |
| FSL BET | `_bet<F>` (for example `_bet070`), `bet2_ANTS` |
| Template projection (ANTs, AFNI 3dQwarp) | `Custum_ANTS_NL`, `Custum_ANTS_Garin`, `CustumANTSGarin_<dil>`, `QWARP`, `Custum_QWARP`, `Custum_QWARPT2`, `Custum_QWARPT2_dil` |
| Deep learning | `synthstrip`, `muSkullStrip_Human`, `muSkullStrip_macaque`, `muSkullStrip_cross_species`, `muSkullStrip_cross_species_dil`, `Custum_Macaque2` (macaque U-Net with a 3-voxel dilation) |
| Other | `NoSkullStrip` (positive voxels), `Manual` |

The template-projection methods register the subject to the template with the ANTs parameters of the `SS1`, `SS2` or `SS3` block of `list_transfo`.

## Functional methods

| family | names |
|---|---|
| AFNI `3dSkullStrip` | `3dSkullStrip`, `3dSkullStrip_monkey`, `3dSkullStrip_monkeynodil`, `3dSkullStrip_dog`, `3dSkullStrip_marmoset` |
| FSL BET | `bet2`, `_bet<F>` |
| Intensity | `CustomNilearn_<low>_<high>`, `CustomNilearnExcludeZeros_<low>_<high>`, `CustomThreshold_<percentile>`, `Vol_sammba_<volume>` |
| Registration to the anatomy | `Custom_ANTS_NL`, `Custom_ANTS_Garin` |
| Deep learning | `muSkullStrip_cross_species` |
| Other | `NoSkullStrip`, `Manual` |

`dilate_mask` and `extra_erode` adjust the functional mask afterwards.

## Manual masks

- `Manual` opens ITK-SNAP (or Freeview, `preftool`) on the image and an initial mask; save your edit as the file indicated in the log (`*_desc-manual.nii.gz`), it then becomes the final mask.
- To impose your own anatomical mask, save it in `anat/native/volumes/masks/` of the subject as `<ID>_final_mask.nii.gz` (first skull-strip) or `<ID>_final_mask_2.nii.gz` (second skull-strip); it is used instead of any method.
- For the functional pipeline, `<ID>_final_mask.nii.gz` in `func/raw/masks/` is resampled to the mean image and used as the functional mask.
- `manualaseg` corrects FreeSurfer's aseg in ITK-SNAP or Freeview; the corrected segmentation is smoothed (`correctaseg`) and used from then on.

The `itk_*` steps open ITK-SNAP on the masks and registrations for a visual check.
