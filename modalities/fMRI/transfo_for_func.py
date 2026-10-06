import os
import glob
import json
import subprocess
import shutil
from os.path import join as opj
import ants
import nibabel as nib
import numpy as np

# SPECIFIC FUNCTION FOR MOTION CORRECTION + TRANSFORMATIONS
def apply_motion_correction_and_transforms(input_4d, reference_image, output_4d,
                                           mat_files_pattern, additional_transformations,
                                           TR, sing_afni, temp_dir=None, interpolator='linear',
                                           metadata_sources=None, description=None):
    """
    One-shot resampling of a 4D series: for volume j, the transform list is
    additional_transformations + [motion_j] (ANTs point-mapping order), applied
    in a single interpolation onto reference_image.

    Volumes are split in memory (ants.ndimage_to_list), so volume order cannot
    depend on file-name sorting; motion matrices are sorted on their
    zero-padded index (img_00000_...). `sing_afni` and `temp_dir` are kept for
    signature compatibility and are not used.
    """
    mat_files = sorted(glob.glob(mat_files_pattern))
    if not mat_files:
        raise ValueError(f"No .mat file found with pattern: {mat_files_pattern}")

    img4d = ants.image_read(input_4d)
    if img4d.dimension != 4:
        raise ValueError(f"{input_4d} is not a 4D image")
    volumes = ants.ndimage_to_list(img4d)
    if len(volumes) != len(mat_files):
        raise ValueError(f"{len(volumes)} volumes but {len(mat_files)} motion matrices for {input_4d}")

    ref_img = ants.image_read(reference_image)
    out_vols = []
    for vol, mat_file in zip(volumes, mat_files):
        transform_list = list(additional_transformations) + [mat_file]
        # explicit whichtoinvert: ANTsPy inverts the first transform of a 2-element
        # [affine.mat, field] list when whichtoinvert is None
        out_vols.append(ants.apply_transforms(fixed=ref_img, moving=vol,
                                              transformlist=transform_list,
                                              whichtoinvert=[False] * len(transform_list),
                                              interpolator=interpolator).numpy().astype(np.float32))

    # geometry of the output grid = reference_image (written through ANTs to get the NIfTI affine)
    tmp_ref = output_4d + '_refgrid_tmp.nii.gz'
    ants.image_write(ref_img, tmp_ref)
    ref_nib = nib.load(tmp_ref)
    affine, header = ref_nib.affine, ref_nib.header.copy()
    os.remove(tmp_ref)

    data_4d = np.stack(out_vols, axis=-1)
    header.set_data_dtype(np.float32)
    out_img = nib.Nifti1Image(data_4d, affine=affine, header=header)
    out_img.header.set_xyzt_units('mm', 'sec')
    zooms = list(out_img.header.get_zooms()[:3]) + [float(TR)]
    out_img.header.set_zooms(zooms)
    out_path = output_4d + '.nii.gz'
    nib.save(out_img, out_path)

    if metadata_sources is None:
        metadata_sources = [input_4d, reference_image] + mat_files + list(additional_transformations)
    metadata = {
        "Sources": metadata_sources,
        "Description": description or f"Motion correction + {len(additional_transformations)} transformations",
        "Command": [f"ants.apply_transforms per volume, {len(additional_transformations) + 1} transforms, "
                    f"single interpolation ({interpolator})"],
        "Parameters": {
            "interpolator": interpolator,
            "TR": TR,
            "motion_matrices": len(mat_files),
            "additional_transformations": len(additional_transformations)
        }
    }
    with open(output_4d + '.json', "w") as outfile:
        json.dump(metadata, outfile, indent=3)

    print(f"Motion correction + transformations completed: {out_path}")
    return output_4d
