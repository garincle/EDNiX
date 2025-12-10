#import
import os
import json
import shutil
import glob

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd
from Tools import diaryfile
from Tools import getpath

from anat.studytemplate import _3_make_template
from anat.studytemplate import _4_skullstrip_template
from anat.studytemplate import _7_stdyT_to_AtlasT

def create(study_template_atlas_folder,Skip_step,which_on, all_ID_max, all_Session_max, all_data_path_max, all_ID,
           all_Session, all_data_path,type_norm,BASE_SS_coregistr,
           BASE_SS_mask, list_transfo,
           template_skullstrip,preftool,check_visualy_final_mask,
           sing_afni, sing_fsl, sing_fs, sing_itk, sing_synstrip,reference,BALSAname, Unetpath):


    # second loop # first loop : (contains steps  3 and 4: reorientation and skull stripping)

    studyacpc_dir,studyprepro_dir,studytransfo_dir, studyvolume_dir, studylabels_dir, studymasks_dir,_,_,_,_,_,_,_,_ = getpath.stytemplate(
        study_template_atlas_folder,reference,BALSAname)

    anat_img           = opj(studyprepro_dir, 'warped_3_adjusted_mean.nii.gz')
    stdy_template      = opj(studyprepro_dir, 'studyTemplate_' + type_norm + '.nii.gz')
    stdy_template_SS = opj(studyprepro_dir, 'studyTemplate_SS_' + type_norm + '.nii.gz')
    stdy_template_mask = opj(studyprepro_dir, 'studyTemplate_mask.nii.gz')
    endname            = '_'.join(['studyTemplate', 'final', 'mask.nii.gz'])

    if not ope(study_template_atlas_folder): os.mkdir(study_template_atlas_folder)

    diary_file = diaryfile.create(opj(study_template_atlas_folder, 'Study_template_BLOCK1'),
                                  'Create the study template')

    if 3 in Skip_step:
        run_cmd.msg('INFO: skip step ' + str(3), diary_file, 'OKGREEN')
    else:
        _3_make_template.make_template(which_on, all_ID_max, all_Session_max, all_data_path_max, all_ID,
                                                        all_Session, all_data_path, type_norm, study_template_atlas_folder,
                                                        sing_afni, diary_file)

    if '4' in Skip_step:
        run_cmd.msg('INFO: skip step ' + str('4'), diary_file, 'OKGREEN')
    else:
        refnb = 0
        for i, j in enumerate(list_transfo):
            if list_transfo[i]["name"] == 'SS3':
                refnb = i
        _4_skullstrip_template.skullstrip_T(stdy_template_SS, anat_img, stdy_template_mask,endname,type_norm, BASE_SS_coregistr, BASE_SS_mask, list_transfo[refnb]["type_of_transform"],  list_transfo[refnb]["affmetric"],
                                                             study_template_atlas_folder, template_skullstrip, preftool,
                                                             check_visualy_final_mask, sing_afni, sing_fsl, sing_fs, sing_itk, sing_synstrip, Unetpath, diary_file)


def use(study_template_atlas_folder,Skip_step,list_transfo, list_atlases,
        BASE_SS,BASE_mask, BASE_atlas_folder,species,stdy_template,fMRImasks,reference,BALSAname,path_label_code,
        sing_afni,sing_wb,which_on,all_data_path_max,
        all_data_path,listTimage,creat_study_template,coregistration_longitudinal,refimagename):

    (studyacpc_dir,dir_prepro, studytransfo_dir, studyvolume_dir, studylabels_dir, studymasks_dir,
     wb_template_dir, wb_template_vol, wb_template_labels, wb_template_masks,
     wb_balsa_dir, wb_balsa_vol, wb_balsa_labels, wb_balsa_masks)= getpath.stytemplate(study_template_atlas_folder,reference,BALSAname)

    diary_file = diaryfile.create(opj(study_template_atlas_folder, 'Study_template_BLOCK2'), 'ADD the atlases')

    #######################################################################
    ###############     finalise the study template     ###################
    #######################################################################


    if which_on == 'max':  # all or max
        all_data_path_temp = all_data_path_max
    elif which_on == 'all':
        all_data_path_temp = all_data_path

    # average the available anat files

    second_template = []
    for Timage in listTimage:

        for data_path in all_data_path_temp:
            (_, _, _, _, _, _, _, _, _, wb_studytemplate_vol, _, _) = getpath.anat(data_path, reference, '',
                                                                                   creat_study_template,
                                                                                   coregistration_longitudinal,
                                                                                   'studyTemplate')


            name = opj(wb_studytemplate_vol, refimagename + 'desc-SS_' + Timage + '.nii.gz')
            print(name)
            if ope(name):
                second_template.append(name)
            else :
                run_cmd.msg('the coregistration '+ Timage + ' image within ' + name + ' does not exist !',diary_file, 'WARNING')

        # Check if all paths exist
        all_paths_exist = all(opi(path) for path in second_template)

        if all_paths_exist:
            TEMP = ' '.join(second_template)
            if not ope(studyvolume_dir):
                os.mkdir(studyvolume_dir)

            command = (sing_afni + '3dMean -overwrite -prefix ' + opj(studyvolume_dir,
                                                           'studyTemplate_' + Timage + '.nii.gz') +
                       ' ' + TEMP)
            run_cmd.run(command, diary_file)

            dictionary = {"Sources": TEMP,
                          "Description": 'Mean.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(studyvolume_dir, 'studyTemplate_' + Timage + '.json'), "w") as outfile:
                outfile.write(json_object)

        if not ope(opj(dir_prepro, 'studyTemplate_' + Timage + '.nii.gz')):
            shutil.copy(opj(studyvolume_dir, 'studyTemplate_' + Timage + '.nii.gz'),
                        opj(dir_prepro, 'studyTemplate_' + Timage + '.nii.gz'))
            shutil.copy(opj(studyvolume_dir, 'studyTemplate_' + Timage + '.json'),
                        opj(dir_prepro, 'studyTemplate_' + Timage + '.json'))

            ################### redefine new atlases variable!!!
    if 7 in Skip_step:
        run_cmd.msg('INFO: skip step ' + str(7), diary_file, 'OKGREEN')
    else:
        nl = 'Run anat._7_stdyT_to_AtlasT.stdyT_to_AtlasT'
        run_cmd.msg(nl, diary_file, 'HEADER')

        refnb = 0
        for i, j in enumerate(list_transfo):
            if list_transfo[i]["name"] == 'stdyT':
                refnb = i

        _7_stdyT_to_AtlasT.stdyT_to_AtlasT(list_transfo[refnb]["affmetricT"], list_atlases, BASE_SS, BASE_atlas_folder,BASE_mask, species, fMRImasks,
                        list_transfo[refnb]["interpol"], list_transfo[refnb]["affmetric"], study_template_atlas_folder, stdy_template,
                        list_transfo[refnb]["type_of_transform"], path_label_code,
                        diary_file, reference, sing_wb, sing_afni)

        # Handle the atlas files
        listatlas = sorted(glob.glob(opj(studylabels_dir, 'studyTemplate_space-acpc_*_dseg.**')))
        for filepath in listatlas:
            name = opb(filepath).split('_')
            newname = '_'.join([name[0]] + name[2:])  # Removes 'space-acpc'
            newpath = opj(opd(filepath), newname)

            if os.path.exists(newpath):
                os.remove(newpath)  # Optional: remove existing file if overwriting is safe
            shutil.move(filepath, newpath)

        # Handle the mask files
        listmask = sorted(glob.glob(opj(studymasks_dir, 'studyTemplate_space-acpc_mask.**')))
        for filepath in listmask:
            name = opb(filepath).split('_')
            newname = '_'.join([name[0]] + name[2:])  # Removes 'space-acpc'
            newpath = opj(opd(filepath), newname)

            if os.path.exists(newpath):
                os.remove(newpath)
            shutil.move(filepath, newpath)

    ###### re-define the variable: study template atlas, etc should now be the new template !!
    BASE_SS           = stdy_template
    BASE_mask         = opj(studymasks_dir, 'studyTemplate_mask.nii.gz')
    Aseg_ref          = opj(studylabels_dir, 'studyTemplate_seg-4FS_dseg.nii.gz')
    return BASE_SS, BASE_mask, Aseg_ref