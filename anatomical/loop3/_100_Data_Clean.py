# import
import os

opj = os.path.join
opi = os.path.isfile

from Tools import run_cmd
from Tools import getpath



def clean(dir_prepro, volumes_dir, masks_dir, ID, diary_file):
    nl = 'Run anatomical._100_Data_Clean.clean'
    run_cmd.msg(nl, diary_file, 'HEADER')


    # PRESERVE final_mask - DO NOT REMOVE
    final_mask = opj(masks_dir, ID + '_final_mask.nii.gz')

    list_to_remove = [opj(dir_prepro, ID + '_space-acpc_desc-64_T1w'),
                      opj(dir_prepro, ID + '_space-acpc_desc-cropped_T1w'),
                      opj(dir_prepro, ID + '_brain_for_Align_Center_T1w'),
                      opj(volumes_dir, ID + '_space-acpc_desc-SS-step1_T1w'),
                      opj(volumes_dir, ID + '_space-acpc_desc-SS-step2_T1w'),
                      opj(dir_prepro, ID + '_space-acpc_desc-64_T2w'),
                      opj(dir_prepro, ID + '_space-acpc_desc-cropped_T2w'),
                      opj(dir_prepro, ID + '_brain_for_Align_Center_T2w'),
                      opj(volumes_dir, ID + '_space-acpc_desc-SS-step1_T2w'),
                      opj(volumes_dir, ID + '_space-acpc_desc-SS-step2_T2w'),
                      opj(masks_dir, ID + '_desc-step1_mask'),
                      opj(masks_dir, ID + '_desc-step2_mask')
                      ]

    # Log what we're preserving
    if opi(final_mask) or opi(final_mask.replace('.nii.gz', '.json')):
        nl = f'PRESERVING: {final_mask} and associated files'
        run_cmd.msg(nl, diary_file, 'OKGREEN')

    for i in range(len(list_to_remove)):
        for z in ['.nii.gz', '.json']:
            file_path = list_to_remove[i] + z
            # Skip if this is the final_mask we want to preserve
            if file_path == final_mask or file_path == final_mask.replace('.nii.gz', '.json'):
                nl = 'INFO: ' + file_path + ' preserved (final_mask).'
                run_cmd.msg(nl, diary_file, 'OKBLUE')
                continue

            if opi(file_path):
                os.remove(file_path)
                nl = 'INFO: ' + file_path + ' removed.'
                run_cmd.msg(nl, diary_file, 'OKGREEN')
            else:
                nl = 'INFO: ' + file_path + ' not found.'