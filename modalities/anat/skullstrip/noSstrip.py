
import json
from Tools import run_cmd

def run(input_for_msk,output_for_mask,sing_afni,diary_file):

    command = (sing_afni + '3dcalc -a ' + input_for_msk +
               ' -expr "step(a)" -prefix ' + output_for_mask + ' -overwrite')
    run_cmd.do(command, diary_file)

    dictionary = {"Sources": input_for_msk,
                  "Description": 'Brain mask (positive voxels).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(output_for_mask[:-7] + '.json', "w") as outfile:
        outfile.write(json_object)