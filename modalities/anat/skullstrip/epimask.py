import os
import json

from nilearn.masking import compute_epi_mask

opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists

from Tools import run_cmd

def run(input_for_msk,output_for_mask,cutoff,exzeros,option2,sing_afni,sing_fs,diary_file):

    # convert to float
    command = (sing_fs + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz')
    run_cmd.run(command, diary_file)

    mask_img = compute_epi_mask(input_for_msk[:-7] + '_float.nii.gz',
                                lower_cutoff=cutoff[0], upper_cutoff=cutoff[1], connected=True, opening=3,
                                exclude_zeros=exzeros, ensure_finite=True)
    mask_img.to_filename(output_for_mask)

    command = (sing_afni + '3dmask_tool -overwrite -prefix ' + output_for_mask +
               ' -input ' + output_for_mask + option2)
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": input_for_msk,
                  "Description": 'Brain mask (compute_epi_mask from Nilearn).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(output_for_mask[:-7] + '.json', "w") as outfile:
        outfile.write(json_object)