#import
import os

opj = os.path.join
ope = os.path.exists

from Tools import run_cmd
from Tools import diaryfile
from Tools import getpath

from anat.loop1 import _1_correct_orient
from anat.loop1 import _2_clean_anat

def run(ID, Session, data_path, path_rawanat,BIDStype, listTimage,otheranat, listTimage_orig, type_norm, orientation,IgotbothT1T2, force_myelin_same_space, Align_img_to_template,list_transfo,masking_img,
        brain_skullstrip_1, anat_ref_path, BASE_SS_coregistr, BASE_SS_mask, BASE_SS,check_visualy_each_img,check_visualy_final_mask, overwrite,bids_dir,Skip_step,
        MNIBcorrect_indiv,animalPosition,humanPosition,sing_afni, sing_fsl, sing_fs, sing_itk,sing_wb, sing_synstrip, Unetpath,preftool):

    NL1 = '############################ work on subject: ' + str(ID) + ' Session ' + str(Session) + ' BLOCK 1 ##################################'
    run_cmd.printcolor(NL1, 'HEADER')

    DIR = os.getcwd()
    run_cmd.printcolor('INFO: Working path : ' + DIR, 'OKGREEN')

    # create the folder's architecture
    path_anat, dir_transfo, _, dir_prepro, dir_native, volumes_dir, labels_dir, masks_dir = getpath.anat(data_path, '', '', False, False, 'native')

    if ope(dir_prepro) == False: os.makedirs(dir_prepro)
    if ope(dir_transfo) == False: os.makedirs(dir_transfo)
    if ope(labels_dir) == False: os.makedirs(labels_dir)
    if ope(masks_dir) == False: os.makedirs(masks_dir)
    if ope(opj(bids_dir, 'QC')) == False: os.makedirs(opj(bids_dir, 'QC'))

    diary_file = diaryfile.create(opj(path_anat, str(ID) + ' session ' + str(Session)), NL1)


    if 1 in Skip_step:
        run_cmd.msg('INFO: skip step ' + str(1), diary_file, 'OKGREEN')
    else:
        _1_correct_orient.correct_orient(BIDStype,listTimage,path_anat,path_anat,ID, Session, otheranat, listTimage_orig, type_norm, MNIBcorrect_indiv,orientation, dir_prepro,dir_transfo,IgotbothT1T2,force_myelin_same_space,list_transfo,overwrite,
                   animalPosition,humanPosition,sing_afni,sing_fs, sing_wb,diary_file)


    if 2 in Skip_step:
        run_cmd.msg('INFO: skip step ' + str(2), diary_file, 'OKGREEN')
    else:
        _2_clean_anat.clean_anat(Align_img_to_template, bids_dir, listTimage,
                                                  list_transfo, ID, Session,
                                                  type_norm, data_path, masking_img, brain_skullstrip_1,
                                                  BASE_SS_coregistr,BASE_SS_mask, BASE_SS, IgotbothT1T2, check_visualy_each_img,
                                                  check_visualy_final_mask,overwrite, anat_ref_path,
                                                  sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb,sing_synstrip,
                                   Unetpath, diary_file,preftool)



