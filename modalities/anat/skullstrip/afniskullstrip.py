
import json

from Tools import run_cmd

def run(input_for_msk,output_for_mask,option1,option2,sing_afni,diary_file):

    command = (sing_afni + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite -input ' + input_for_msk +
               option1)
    run_cmd.run(command, diary_file)

    command = (sing_afni + '3dmask_tool -overwrite -prefix ' + output_for_mask +
               ' -input ' + output_for_mask + option2)
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": input_for_msk,
                  "Description": 'Brain mask (3dSkullStrip from AFNI).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(output_for_mask[:-7] + '.json', "w") as outfile:
        outfile.write(json_object)

