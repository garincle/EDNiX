import os
import shutil
import json

opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists

from Tools import run_cmd
from Tools.extract_filename import extract_filename

def do(input_for_msk,output_for_mask,model,Unetpath,option2,sing_afni,diary_file):

    Unet_dir = opj(Unetpath, 'NHP-BrainExtraction', 'UNet_Model')

    process_dir = opd(output_for_mask)
    print(input_for_msk)
    command = ('python3 ' + opj(Unet_dir, 'muSkullStrip.py') +
               ' -in ' + input_for_msk +
               ' -model ' + opj(Unet_dir, 'models', model) +
               ' -out ' + process_dir)
    print(command)
    run_cmd.run(command, diary_file)

    shutil.copyfile(opj(process_dir, extract_filename(input_for_msk) + '_pre_mask.nii.gz'), output_for_mask)

    command = (sing_afni + '3dmask_tool -overwrite -prefix ' + output_for_mask +
               ' -input ' + output_for_mask + option2)
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": input_for_msk,
                  "Description": 'Brain mask (U-Net).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(output_for_mask[:-7] + '.json', "w") as outfile:
        outfile.write(json_object)