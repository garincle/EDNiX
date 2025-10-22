# import
import os

opj = os.path.join
opi = os.path.isfile

from Tools import run_cmd
from Tools import getpath

#################################################################################################
#### get rid of the unnecessary step files
#################################################################################################

def clean(all_ID, all_Session, all_data_path,diary_file):

    nl = 'Run anatomical._100_Data_Clean.clean'
    run_cmd.msg(nl, diary_file, 'HEADER')

    for ID, Session, data_path in zip(all_ID, all_Session, all_data_path):

        _, _, _, dir_prepro, _, volumes_dir, labels_dir, masks_dir = getpath.anat(data_path,'', '', False, False, 'native')

        list_to_remove = [opj(dir_prepro, ID + '_space-acpc_desc-64_T1w'),
                          opj(dir_prepro, ID + '_space-acpc_desc-cropped_T1w'),
                          opj(dir_prepro, ID + '_brain_for_Align_Center_T1w'),
                          opj(volumes_dir,ID + '_space-acpc_desc-SS-step1_T1w'),
                          opj(dir_prepro, ID + '_space-acpc_desc-64_T2w'),
                          opj(dir_prepro, ID + '_space-acpc_desc-cropped_T2w'),
                          opj(dir_prepro, ID + '_brain_for_Align_Center_T2w'),
                          opj(volumes_dir,ID + '_space-acpc_desc-SS-step1_T2w'),
                          opj(masks_dir,  ID + '_desc-step1_mask'),
                          opj(masks_dir,  ID + '_desc-step2_mask')
                          ]


        for i in range(len(list_to_remove)):

            for z in ['.nii.gz','.json']:
                if opi(list_to_remove[i] + z):
                    os.remove(list_to_remove[i] + z)
                    nl = 'INFO: ' + list_to_remove[i] + z + ' removed.'
                    run_cmd.msg(nl, diary_file, 'OKGREEN')
                else:
                    nl = 'INFO: ' + list_to_remove[i] + z + '  not found.'
                    run_cmd.msg(nl, diary_file, 'WARNING')
