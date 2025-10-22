
import json

from Tools import run_cmd
from anatomical.skullstrip import Histrogram_mask_EMB

def run(input_for_msk,output_for_mask,setmasker,sing_fs,diary_file):

    # convert to float
    command = (sing_fs + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz')
    run_cmd.run(command, diary_file)

    nichols_masker = Histrogram_mask_EMB.HistogramMask()
    nichols_masker.inputs.in_file = input_for_msk[:-7] + '_float.nii.gz'

    nichols_masker.inputs.volume_threshold    = setmasker[0]
    nichols_masker.inputs.upper_cutoff        = setmasker[1]
    nichols_masker.inputs.lower_cutoff        = setmasker[2]
    nichols_masker.inputs.intensity_threshold = setmasker[3]
    nichols_masker.inputs.opening             = setmasker[4]
    nichols_masker.inputs.closing             = setmasker[5]

    nichols_masker.inputs.dilation_size = (1, 2, 3)
    nichols_masker.inputs.connected = True
    nichols_masker.inputs.out_file = output_for_mask

    res = nichols_masker.run()  # doctest: +SKIP

    dictionary = {"Sources": input_for_msk,
                  "Description": 'Brain mask (Histrogram_mask_EMB from sammba).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(output_for_mask[:-7] + '.json', "w") as outfile:
        outfile.write(json_object)
