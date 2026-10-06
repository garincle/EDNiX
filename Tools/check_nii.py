#import
import os
import numpy as np
import ants
import nibabel as nib
from nilearn.image import resample_to_img
from nilearn.image import new_img_like
opj = os.path.join

from Tools import run_cmd

def comphd(img1,img2):
    """
    True if both images share dimensions, spacing (1e-3), direction (1e-3) and
    an origin closer than one voxel diagonal. The previous version first compared
    np.ceil(origin), which accepted offsets up to ~1 mm, i.e. ~10 voxels at 0.1 mm.
    """
    hd1 = ants.image_header_info(img1)
    hd2 = ants.image_header_info(img2)
    test = []
    test.append(tuple(hd1['dimensions'][0:3]) == tuple(hd2['dimensions'][0:3]))
    test.append(bool(np.allclose(hd1['spacing'][0:3], hd2['spacing'][0:3], atol=1e-3)))
    dist_origin = float(np.linalg.norm(np.asarray(hd1['origin'][0:3]) - np.asarray(hd2['origin'][0:3])))
    voxel_diag = float(np.linalg.norm(np.asarray(hd1['spacing'][0:3])))
    test.append(dist_origin < voxel_diag)
    d1 = np.asarray(hd1['direction'])[0:3, 0:3]
    d2 = np.asarray(hd2['direction'])[0:3, 0:3]
    same_dir = bool(np.allclose(d1, d2, atol=1e-3))
    if not same_dir:
        print('beware that the geometry seems to differ between the two files')
        print(d1)
        print(d2)
    test.append(same_dir)
    return all(test)


def compare_headers(source, target, fields_to_compare=None, tolerance=1e-6, check_obliquity=True):
    """
    Compare ANTs image headers. Integer fields (dimensions, numpixels, pixeltype)
    must match exactly; real-valued fields (spacing, origin, direction) are
    compared with an absolute tolerance. Returns {field: (source_value, target_value)}.
    """
    hd1 = ants.image_header_info(source)
    hd2 = ants.image_header_info(target)
    if fields_to_compare is None:
        fields_to_compare = ['dimensions', 'spacing', 'direction', 'pixeltype']
    differences = {}
    for field in fields_to_compare:
        val1, val2 = hd1.get(field), hd2.get(field)
        if field in ('spacing', 'origin', 'direction'):
            v1 = np.asarray(val1, dtype=float)
            v2 = np.asarray(val2, dtype=float)
            if field == 'direction':
                v1, v2 = v1[0:3, 0:3], v2[0:3, 0:3]
            else:
                v1, v2 = v1[0:3], v2[0:3]
            if v1.shape != v2.shape or not np.allclose(v1, v2, atol=tolerance):
                differences[field] = (v1.tolist(), v2.tolist())
        elif field == 'dimensions':
            if tuple(val1[0:3]) != tuple(val2[0:3]):
                differences[field] = (str(val1), str(val2))
        elif val1 != val2:
            differences[field] = (val1, val2)

    if check_obliquity:
        def max_tilt_degrees(img):
            direction = np.asarray(ants.image_header_info(img)['direction'])[0:3, 0:3]
            angles = [np.degrees(np.arccos(np.clip(abs(np.dot(direction[:, i], np.eye(3)[:, i])), -1, 1)))
                      for i in range(3)]
            return np.max(np.abs(angles))
        tilt1, tilt2 = max_tilt_degrees(source), max_tilt_degrees(target)
        if not np.isclose(tilt1, tilt2, atol=0.1):
            differences['obliquity_deg_from_plumb'] = (round(tilt1, 3), round(tilt2, 3))
    return differences


def resamp(source,target,imgtype,path_code,labelname,diary_file,sing_wb):
    """
    Resample `source` onto the grid of `target` (in place) when their spatial
    grids differ. Origin and direction are now part of the comparison: grids
    with identical dimensions/spacing but a different origin were previously
    treated as identical and combined voxel-wise without resampling.
    """
    fields_to_compare = ['dimensions', 'spacing', 'origin', 'direction']

    if imgtype == 'label' or imgtype == 'msk':
        interp = 'nearest'
    elif imgtype == 'pet' or imgtype == 'bold':
        interp = 'linear'
    else:
        interp = 'continuous'

    differences = compare_headers(source, target, fields_to_compare=fields_to_compare,
                                  tolerance=1e-4, check_obliquity=False)

    if differences:
        run_cmd.msg("Differences found in the following fields:", diary_file, 'OKGREEN')
        for field, values in differences.items():
            run_cmd.msg(str(field), diary_file, 'OKGREEN')
            run_cmd.msg('Image 1:' + str(values[0]), diary_file, 'OKGREEN')
            run_cmd.msg('Image 2:' + str(values[1]), diary_file, 'OKGREEN')

        dummy = resample_to_img(source, target, interpolation=interp)
        dummy.to_filename(source)

        if imgtype == 'label':
            cmd = (sing_wb + 'wb_command -volume-label-import' + ' ' + source +
                   ' ' + opj(path_code, labelname + '_label.txt') + ' ' + source + ' -drop-unused-labels')
            run_cmd.run(cmd, diary_file)


def resamp_no_check(source,target,imgtype,path_code,labelname,diary_file,sing_wb):

    if imgtype == 'label' or imgtype == 'msk':
        interp = 'nearest'
    elif imgtype == 'pet' or imgtype == 'bold':
        interp = 'linear'
    else:
        interp = 'continuous'

    dummy = resample_to_img(source, target, interpolation=interp)
    dummy.to_filename(source)
    extracted_data = nib.load(source).get_fdata()
    labeled_img2 = new_img_like(target, extracted_data, copy_header=True)
    labeled_img2.to_filename(source)

    if imgtype == 'label':
        cmd = (sing_wb + 'wb_command -volume-label-import' + ' ' + source +
               ' ' + opj(path_code, labelname + '_label.txt') + ' ' + source + ' -drop-unused-labels')
        run_cmd.run(cmd, diary_file)

def normalize_anat_type(name: str) -> str:
    """Normalize anat type to 'T1w' or 'T2w'."""
    if not name:  # empty string or None
        return ''
    if 'T1w' in name:
        return 'T1w'
    if 'T2w' in name:
        return 'T2w'
    if 'T1' in name:
        return 'T1w'
    if 'T2' in name:
        return 'T2w'
    raise ValueError(
        f"Cannot infer T1/T2 suffix from filename '{name}'. "
        "Expected a chunk containing 'T1' or 'T2'.")

def keep_header(input_file, reference_image):
    """
    Overwrite the header of input_file with the header of reference_image,
    WITHOUT resampling. Only valid when both images share the same voxel grid;
    raises otherwise (the previous version silently relabelled mismatched arrays).
    """
    src = nib.load(input_file)
    ref = nib.load(reference_image) if isinstance(reference_image, str) else reference_image
    if src.shape[:3] != ref.shape[:3]:
        raise ValueError(f'keep_header: spatial shapes differ {src.shape[:3]} vs {ref.shape[:3]} '
                         f'({input_file} / {reference_image}); resample instead.')
    zooms_src = np.asarray(src.header.get_zooms()[:3])
    zooms_ref = np.asarray(ref.header.get_zooms()[:3])
    if not np.allclose(zooms_src, zooms_ref, rtol=1e-3):
        raise ValueError(f'keep_header: voxel sizes differ {zooms_src} vs {zooms_ref}; resample instead.')
    extracted_data = src.get_fdata()
    labeled_img2 = new_img_like(ref, extracted_data, copy_header=True)
    labeled_img2.to_filename(input_file)
    return input_file
