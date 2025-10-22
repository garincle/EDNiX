import nibabel as nib
from Tools import run_cmd, get_orientation
import numpy as np
import json

def make_iso(img, img_out, diary_file, sing_afni, type, overwrite):
    # Load the image directly
    img_nib = nib.load(img)
    # Get voxel sizes
    delta_x, delta_y, delta_z = [str(round(abs(x), 10)) for x in img_nib.header.get_zooms()[:3]]
    size = np.min([float(delta_x), float(delta_y), float(delta_z)])
    # Usage
    orient = get_orientation.get_orientation_nibabel(img)

    if type == 'seg':
        rmode = 'NN'
    elif type == 'anat':
        rmode = 'Cu'

    command = (sing_afni + '3dresample' + overwrite +
               ' -orient ' + orient +
               ' -prefix ' + img_out +
               ' -dxyz ' + str(size) + ' ' + str(size) + ' ' + str(size) +
               ' -rmode ' + rmode + ' -input ' + img)

    run_cmd.run(command, diary_file)
    dictionary = {"Sources": img,
                  "Description": 'make '+ str(size) + 'mm isotropric image with mininum detected voxel size:'}
    json_object = json.dumps(dictionary, indent=2)
    with open(img_out.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)