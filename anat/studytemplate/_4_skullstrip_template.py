#################################################
########    creat brain image of animal  ########
#################################################
#co-register linear to the new sty template...
##############################################################################
####### CREATE THE STUDY TEMPLATE (IF YOU WANT ON) ###########################
##############################################################################

import os
import json

opj = os.path.join
ope = os.path.exists

from Tools import run_cmd
from Tools import getpath

from anat.skullstrip import Skullstrip_method

def skullstrip_T(stdy_template_SS, stdy_template,input4mask,output4msk,endname,type_norm,BASE_SS_coregistr, BASE_SS_mask, type_of_transform, aff_metric_ants,
                 study_template_atlas_folder, template_skullstrip,preftool,
                 check_visualy_final_mask,sing_afni, sing_fsl, sing_fs, sing_itk, sing_synstrip, Unetpath, diary_file):

    nl = 'Run anat._4_skullstrip_template.skullstrip_T'
    run_cmd.msg(nl, diary_file,'HEADER')

    (studyacpc_dir, studyprepro_dir, studytransfo_dir, studyvolume_dir, studylabels_dir, studymasks_dir,
     _, _, _, _, _, _, _, _) = getpath.stytemplate(study_template_atlas_folder, '', '')

    if not ope(studyvolume_dir):
        os.makedirs(studyvolume_dir)
    if not ope(studymasks_dir):
        os.makedirs(studymasks_dir)

    stdy_template_mask = Skullstrip_method.Skullstrip_method(template_skullstrip,endname,input4mask,output4msk,
                                                                                   type_norm,studyprepro_dir, studymasks_dir, BASE_SS_coregistr, BASE_SS_mask,
                                                                                   type_of_transform, '', aff_metric_ants,
                                                                                   check_visualy_final_mask, sing_afni, sing_fsl, sing_fs, sing_itk, sing_synstrip, Unetpath, diary_file,preftool)

    nl = 'Run anat._4_skullstrip_template.skullstrip_T, after Skullstrip_method'
    run_cmd.msg(nl, diary_file,'HEADER')

    ## extract brain
    command = (sing_afni + '3dcalc -overwrite -a ' + stdy_template_mask + ' -b ' + input4mask +
               ' -expr "(a*b)" -prefix ' + stdy_template_SS)
    run_cmd.do(command, diary_file)

    dictionary = {"Sources": [stdy_template_mask,
                              input4mask],
                  "Description": 'Skull stripping.', }
    json_object = json.dumps(dictionary, indent=2)
    with open(stdy_template_SS.replace('.nii.gz','.json'), "w") as outfile:
        outfile.write(json_object)

    ## extract brain
    command = (sing_afni + '3dcalc -overwrite -a ' + input4mask +
               ' -expr "a" -prefix ' + stdy_template)
    run_cmd.do(command, diary_file)

    dictionary = {"Sources": [input4mask],
                  "Description": 'Template orig', }
    json_object = json.dumps(dictionary, indent=2)
    with open(stdy_template.replace('.nii.gz','.json'), "w") as outfile:
        outfile.write(json_object)



