import os
import nilearn
import nibabel as nib
import subprocess
from nilearn.masking import compute_epi_mask

opj = os.path.join
ope = os.path.exists
spgo = subprocess.getoutput


def build_group_mask(
    mean_imgs,
    mask_func,
    output_results1,
    sing_afni,
    studytemplatebrain,
    lower_cutoff,
    upper_cutoff,
    opening=1,
    method_mask_func='mask_func_over_Gray',
):
    """
    Build the group brain mask and save it to output_results1.

    Parameters
    ----------
    mean_imgs        : list of mean functional image paths
    mask_func        : path to the provided gray/white matter mask
    output_results1  : output directory (Results/)
    sing_afni        : singularity AFNI call prefix
    studytemplatebrain : template brain path (used for 'onlyprovidedmask')
    lower_cutoff     : lower cutoff for compute_epi_mask
    upper_cutoff     : upper cutoff for compute_epi_mask
    opening          : morphological opening radius for compute_epi_mask (default 1)
    method_mask_func : 'mask_func_over_Gray'   ? EPI mask AND provided mask
                       'mask_func_minus_White' ? EPI mask minus white matter mask
                       'onlyprovidedmask'      ? use provided mask directly

    Returns
    -------
    mask_overlap : path to the final overlap mask nii.gz
    """
    mask_path    = opj(output_results1, "mask_mean_func.nii.gz")
    mask_orig    = opj(output_results1, "mask_mean_func_orig.nii.gz")
    mask_overlap = opj(output_results1, "mask_mean_func_overlapp.nii.gz")

    if method_mask_func in ('mask_func_over_Gray', 'mask_func_minus_White'):
        mean_imgs_rs = nilearn.image.concat_imgs(
            mean_imgs, ensure_ndim=None, memory=None, memory_level=0,
            auto_resample=True, verbose=0,
        )
        mask_img = compute_epi_mask(
            mean_imgs_rs,
            lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff,
            connected=True, opening=opening,
            exclude_zeros=True, ensure_finite=True,
        )
        mask_img.to_filename(mask_path)

        print(spgo(
            f"{sing_afni} 3dmask_tool -overwrite "
            f"-prefix {mask_path} -input {mask_path} -fill_holes"
        ))
        print(spgo(
            f"{sing_afni} 3dresample -master {mask_path} "
            f"-prefix {mask_orig} -input {mask_func} -overwrite -bound_type SLAB"
        ))

        if method_mask_func == 'mask_func_over_Gray':
            print(spgo(
                f"{sing_afni} 3dcalc -a {mask_path} -b {mask_orig} "
                f"-expr 'a*b' -prefix {mask_overlap} -overwrite"
            ))
        else:  # mask_func_minus_White
            print(spgo(
                f"{sing_afni} 3dcalc -a {mask_path} -b {mask_orig} "
                f"-expr 'ispositive(a-b)' -prefix {mask_overlap} -overwrite"
            ))

        print(spgo(
            f"{sing_afni} 3dmask_tool -overwrite "
            f"-prefix {mask_overlap} -input {mask_overlap} -fill_holes"
        ))

    elif method_mask_func == 'onlyprovidedmask':
        print(spgo(
            f"{sing_afni} 3dresample -master {studytemplatebrain} "
            f"-prefix {mask_overlap} -input {mask_func} -overwrite -bound_type SLAB"
        ))

    else:
        raise ValueError(
            f"Unknown method_mask_func='{method_mask_func}'. "
            f"Choose from: 'mask_func_over_Gray', 'mask_func_minus_White', 'onlyprovidedmask'."
        )

    return mask_overlap