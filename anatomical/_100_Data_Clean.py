import shutil
import subprocess
import os
import datetime

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
# Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output

#################################################################################################
####Seed base analysis
#################################################################################################
def clean(all_ID, all_Session, all_data_path,diary_file):
    ct = datetime.datetime.now()
    nl = 'Run anatomical._100_Data_Clean.clean'
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{nl}')

    for ID, Session, data_path in zip(all_ID, all_Session, all_data_path):

            # The anatomy
            path_anat     = opj(data_path, 'anat')
            dir_native    = opj(path_anat, 'native')
            dir_prepro    = opj(dir_native, '01_preprocess')
            wb_native_dir = opj(dir_native, '02_Wb')
            volumes_dir   = opj(wb_native_dir, 'volumes')
            labels_dir    = opj(volumes_dir, 'labels')
            masks_dir     = opj(volumes_dir, 'masks')

            list_to_remove = [opj(dir_prepro, 'mask_tmpT2w.nii.gz'),
                              opj(dir_prepro, ID + '_acpc_64T2w.nii.gz'),
                              opj(dir_prepro, ID + '_brain_for_Align_CenterT2w.nii.gz'),
                              opj(dir_prepro, ID + '_acpc_64_orig_3dAllineateT2w.nii.gz'),
                              opj(masks_dir, ID + 'T2w_mask_2.nii.gz'),
                              opj(masks_dir, ID + 'T2w_mask_1_rsp.nii.gz'),
                              opj(masks_dir, ID + 'T2w_mask_1.nii.gz'),
                              opj(masks_dir, ID + '_mask_rsp_T2w.nii.gz'),
                              opj(labels_dir, ID + 'T2w_mask_2.nii.gz'),
                              opj(labels_dir, ID + 'T2w_mask_1_rsp.nii.gz'),
                              opj(labels_dir, ID + 'T2w_mask_1.nii.gz'),
                              opj(labels_dir, ID + '_mask_rsp_T2w.nii.gz'),
                              opj(dir_prepro, 'mask_tmpT2w.json'),
                              opj(dir_prepro, ID + '_acpc_64T2w.json'),
                              opj(dir_prepro, ID + '_brain_for_Align_CenterT2w.json'),
                              opj(dir_prepro, ID + '_acpc_64_orig_3dAllineateT2w.json'),
                              opj(masks_dir, ID + 'T2w_mask_2.json'),
                              opj(masks_dir, ID + 'T2w_mask_1_rsp.json'),
                              opj(masks_dir, ID + 'T2w_mask_1.json'),
                              opj(masks_dir, ID + '_mask_rsp_T2w.json'),
                              opj(labels_dir, ID + 'T2w_mask_2.json'),
                              opj(labels_dir, ID + 'T2w_mask_1_rsp.json'),
                              opj(labels_dir, ID + 'T2w_mask_1.json'),
                              opj(labels_dir, ID + '_mask_rsp_T2w.json')
                              ]

            for remove_data in list_to_remove:
                if ope(remove_data):
                    os.remove(remove_data)
                    print(bcolors.OKGREEN + 'INFO: ' + remove_data + ' removed !' + bcolors.ENDC)
                else:
                    print(bcolors.OKGREEN + 'INFO: ' + remove_data + ' not found' + bcolors.ENDC)

    diary.write(f'\n')
    diary.close()
