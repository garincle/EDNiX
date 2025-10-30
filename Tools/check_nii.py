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
    # check quality
    hd1 = ants.image_headevr_info(img1)
    hd2 = ants.image_header_info(img2)
    test = []
    test.append(hd1['dimensions'][0:3] == hd2['dimensions'][0:3])
    test.append(np.around(hd1['spacing'],3) == np.around(hd2['spacing'],3))

    ori = all(np.ceil(hd1['origin']) == np.ceil(hd2['origin']))
    if ori == False:
        dist1 = np.sqrt(
            np.sqrt(np.power((hd1['origin'][0] - hd2['origin'][0]),2) +
            np.power((hd1['origin'][1] - hd2['origin'][1]),2) +
            np.power((hd1['origin'][2] - hd2['origin'][2]),2)))
        dist2 = np.sqrt(
            np.sqrt(np.power(hd1['spacing'][0],2) +
            np.power(hd1['spacing'][1],2) +
            np.power(hd1['spacing'][2],2)))
        if dist2 - dist1 > 0:
            ori = True
    test.append(ori)
    for i in range(3):
        res = all(np.around(hd1['direction'][i], 3) == np.around(hd2['direction'][i], 3))
        if res == False:
            print('beware that the geometry seems to differ between the two files')
            print(hd1['direction'][i])
            print(hd2['direction'][i])
        test.append(res)
    var = all(test)
    return var


def compare_headers(source, target, fields_to_compare=None, tolerance=1e-6, check_obliquity=True):
    """
    Compare ANTs image headers and optionally check obliquity.

    Args:
        source (ANTsImage): Source image.
        target (ANTsImage): Target image.
        fields_to_compare (list): Header fields to compare. Defaults to key geometric fields.
        tolerance (float): Tolerance for numeric comparisons.
        check_obliquity (bool): If True, compares maximum tilt from RAS axes.

    Returns:
        dict: Differences found (field: (source_value, target_value)).
    """
    hd1 = ants.image_header_info(source)
    hd2 = ants.image_header_info(target)

    # Default fields to compare (expandable)
    if fields_to_compare is None:
        fields_to_compare = [
            'dimensions', 'spacing', 'direction',
            'pixeltype', 'numpixels'  # Added basic metadata
        ]

    differences = {}
    # --- Core header comparison ---
    for field in fields_to_compare:
        val1, val2 = hd1.get(field), hd2.get(field)

        if isinstance(val1, (tuple)):
            if val1 != val2:
                differences[field] = (str(val1), str(val2))

        # Handle arrays (e.g., direction, spacing)
        elif isinstance(val1, (np.ndarray)):
            val1, val2 = np.array(val1), np.array(val2)
            if not np.allclose(val1, val2, atol=tolerance):
                differences[field] = (val1.tolist(), val2.tolist())  # Convert for readability
        # Handle scalars (e.g., numpixels)
        elif val1 != val2:
            differences[field] = (val1, val2)

    # --- Optional: Obliquity check ---
    if check_obliquity:
        def max_tilt_degrees(img):
            img1 = ants.image_read(img)
            direction = np.array(img1.direction)
            angles = [
                np.degrees(np.arccos(np.clip(np.dot(direction[:, i], np.eye(3)[:, i]), -1, 1)))
                for i in range(3)
            ]
            return np.max(np.abs(angles))

        tilt1, tilt2 = max_tilt_degrees(source), max_tilt_degrees(target)
        if not np.isclose(tilt1, tilt2, atol=0.1):
            differences['obliquity_deg_from_plumb'] = (round(tilt1, 3), round(tilt2, 3))

    return differences

def resamp(source,target,imgtype,path_code,labelname,diary_file,sing_wb):

    # Compare specific fields
    fields_to_compare = ['dimensions', 'spacing']
    tolerance = 0.0000006

    if imgtype == 'label' or imgtype == 'msk':
        interp = 'nearest'
    elif imgtype == 'pet' or imgtype == 'bold':
        interp = 'linear'
    else:
        interp = 'continuous'

    differences = compare_headers(source, target, fields_to_compare=None, tolerance=1e-6, check_obliquity=True)

    if differences:
        # Print differences
        run_cmd.msg("Differences found in the following fields:", diary_file, 'OKGREEN')
        for field, values in differences.items():
            run_cmd.msg(str(field), diary_file, 'OKGREEN')
            run_cmd.msg('Image 1:' + str(values[0]), diary_file, 'OKGREEN')
            run_cmd.msg('Image 2:' + str(values[1]), diary_file, 'OKGREEN')

        dummy = resample_to_img(source, target, interpolation=interp)
        dummy.to_filename(source)
        extracted_data = nib.load(source).get_fdata()
        labeled_img2 = new_img_like(target, extracted_data, copy_header=True)
        labeled_img2.to_filename(source)

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
        cmd = (sing_wb + 'wb_command -volume-label-iimgtypemport' + ' ' + source +
               ' ' + opj(path_code, labelname + '_label.txt') + ' ' + source + ' -drop-unused-labels')
        run_cmd.run(cmd, diary_file)

def normalize_anat_type(name: str) -> str:
    """Normalize anatomical type to 'T1w' or 'T2w'."""
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
    Ensures header conservation when transforming extracted data to NIfTI format.

    Args:
        input_file (str): Input file path
        reference_image: Reference image to copy the header from

    Returns:
        str: Path to the generated file
    """
    # Load extracted data
    extracted_data = nib.load(input_file).get_fdata()

    # Create new image with reference header
    labeled_img2 = new_img_like(reference_image, extracted_data, copy_header=True)

    # Save result (overwrites original file)
    labeled_img2.to_filename(input_file)

    return input_file