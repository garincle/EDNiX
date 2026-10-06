import os
import json

opj = os.path.join
opd = os.path.dirname

from Tools import run_cmd
from Tools.extract_filename import extract_filename

def do(input_for_msk,output_for_mask,sing_synstrip,diary_file):
    process_dir = opd(output_for_mask)

    command = (sing_synstrip + '-o ' + opj(process_dir, extract_filename(input_for_msk) + 'skullstriped.nii.gz') +
               ' -m ' + output_for_mask + ' -i ' + input_for_msk + ' --no-csf')
    print(command)
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": input_for_msk,
                  "Description": 'Brain mask (synthstrip).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(output_for_mask[:-7] + '.json', "w") as outfile:
        outfile.write(json_object)