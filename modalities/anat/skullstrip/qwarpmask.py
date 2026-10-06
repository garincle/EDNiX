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

def run(input_for_msk,output_for_mask,process_dir,masks_dir,masking_img,BASE_SS_coregistr,BASE_SS_mask,correl,maxlev,option2,sing_afni,diary_file):

    command = (sing_afni + '3dQwarp -overwrite' + correl + ' -iwarp' +
               ' -base ' + BASE_SS_coregistr +
               ' -prefix ' + opj(process_dir,
                                 extract_filename(input_for_msk) + 'template_to_' + masking_img + '_AFNIQ.nii.gz') +
               ' -source ' + input_for_msk + ' -maxlev '+ str(maxlev) + ' -resample')
    run_cmd.run(command, diary_file)

    command = (sing_afni + '3dNwarpApply -nwarp ' + opj(process_dir, extract_filename(
        input_for_msk) + 'template_to_' + masking_img + '_AFNIQ_WARPINV.nii.gz') +
               ' -source ' + BASE_SS_mask + ' -master ' + input_for_msk +
               ' -interp NN -prefix ' + opj(masks_dir, extract_filename(input_for_msk) + 'template_brainmask.nii.gz') +
               ' -overwrite')
    run_cmd.run(command, diary_file)

    #Ex_Mask = opj(masks_dir, extract_filename(input_for_msk) + 'mask_tmp' + masking_img + '.nii.gz')
    if option2== '':
        command = (sing_afni + '3dmask_tool -overwrite -fill_holes -prefix ' + output_for_mask +
                   ' -input ' + opj(masks_dir, extract_filename(input_for_msk) + 'template_brainmask.nii.gz'))
    else:
        command = (sing_afni + '3dmask_tool -overwrite -fill_holes -prefix ' + output_for_mask +
                   ' -input ' + opj(masks_dir, extract_filename(input_for_msk) + 'template_brainmask.nii.gz') +
                   option2)
    run_cmd.run(command, diary_file)

    #shutil.copyfile(Ex_Mask, output_for_mask)
    dictionary = {"Sources": [input_for_msk,
                              BASE_SS_coregistr,
                              BASE_SS_mask],
                  "Description": 'Brain mask (3dQwarp from AFNI).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(output_for_mask[:-7] + '.json', "w") as outfile:
        outfile.write(json_object)
