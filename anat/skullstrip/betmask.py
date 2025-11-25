
import json
from Tools import run_cmd


def do(input_for_msk,output_for_mask,brain_skullstrip,sing_fsl,diary_file):

    # Extract the last two digits to use as the -f value
    # and Create the approximate brain mask using bet2
    command = (sing_fsl + 'bet2 ' + input_for_msk + ' ' + output_for_mask[:-12] +
               ' -f ' + brain_skullstrip[-3:] + ' -m -n')
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": input_for_msk,
                  "Description": 'Brain mask (bet from FSL).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(output_for_mask[:-12] + '.json', "w") as outfile:
        outfile.write(json_object)