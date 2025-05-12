import os
import subprocess
import glob
import shutil
from fonctions import extract_filename
import json
from nilearn.image import resample_to_img
from nilearn.image import math_img

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput


##########################################
###########define orientation#############
##########################################

def correct_orient(BIDStype, 
                   listTimage, 
                   path_anat, 
                   ID, 
                   Session, 
                   otheranat, 
                   type_norm, 
                   deoblique_exeption1, 
                   deoblique_exeption2, 
                   deoblique, 
                   orientation, 
                   dir_prepro, 
                   IgotbothT1T2, 
                   overwrite):

    for Timage in listTimage:
        # get variables from json
        if BIDStype == 1:
            list_anat = sorted(glob.glob(opj(path_anat, 'sub-' + ID + '_ses-' + str(Session) + '_run-*' + Timage + '.nii.gz')))
        if BIDStype == 2:
            list_anat = sorted(glob.glob(opj(path_anat, 'sub-' + ID + '_' + Timage + '.nii.gz')))

        print(opj(path_anat, 'sub-' + ID + '_ses-' + str(Session) + '_run-*' + Timage + '.nii.gz'))
        print(list_anat)

        if len(list_anat)==1:
            print('test for one run')
            #list_anat = [opj(path_anat, 'sub-' + ID + '_ses-' + str(Session) + '_' + Timage + '.nii.gz')]
            print(ope(list_anat[0]))
            #shutil.copyfile(list_anat[0], opj(path_anat, ID + 'template_indiv' + Timage + '.nii.gz'))
            command = SINGULARITY + 'run' + AFNI_sif  + '3dcalc -a ' + list_anat[0] + \
            ' -prefix ' + opj(path_anat, ID + 'template_indiv' + Timage + '.nii.gz') + ' -expr "a"' + overwrite
            spco([command], shell=True)

            ###################################
            ###########remove norm anat #########
            ###################################
        elif len(list_anat)>1:
            for i, anat in enumerate(list_anat):
                read_json = opj(opd(list_anat[i]), extract_filename(list_anat[i]) + '.json')
                ff = open(read_json)
                anat_T1 = json.load(ff)
                ImageType = anat_T1["ImageType"]

                # Compare the contents of the two files
                if 'NORM' in ImageType:
                    print('Removing ' + extract_filename(list_anat[i]) + ' as it is a NORM')
                    list_anat.pop(i)
                    # After removing the file, no need to increment i since the next file takes the place of the removed one

            if len(list_anat) > 1:
                ANAT= ' '.join(list_anat)
                command = SINGULARITY + 'run' + FS_sif  + 'mri_robust_template --mov ' + ANAT + \
                ' --template ' + opj(path_anat, ID + 'template_indiv' + Timage + '.nii.gz') + \
                ' --satit --inittp 1 --fixtp --noit --iscale --average 0'
                spco([command], shell=True)
            else:
                #shutil.copyfile(list_anat[0], opj(path_anat, ID + 'template_indiv' + Timage + '.nii.gz'))
                command = SINGULARITY + 'run' + AFNI_sif  + '3dcalc -a ' + list_anat[0] + \
                ' -prefix ' + opj(path_anat, ID + 'template_indiv' + Timage + '.nii.gz') + ' -expr "a"' + overwrite
                spco([command], shell=True)

        #command = '3drefit -space ORIG ' + opj(path_anat, ID + 'template_indiv' + Timage + '.nii.gz')
        #spco([command], shell=True)
        ####################################################
        ####### creat indiv T1/T2 template if necessary ####
        ####################################################

    if IgotbothT1T2 == True:
        #spco(['3dresample', '-master', opj(dir_prepro, ID + '_acpc_cropped' + type_norm + '.nii.gz'), '-prefix', opj(volumes_dir, ID + '_' + otheranat + '_template_RSPL.nii.gz'), '-input', opj(dir_prepro, ID + '_' + otheranat + '.nii.gz'), '-overwrite'])
        caca2 = resample_to_img(opj(path_anat, ID + 'template_indiv' + otheranat + '.nii.gz'),
                                opj(path_anat, ID + 'template_indiv' + type_norm + '.nii.gz'),
                                interpolation='nearest')
        caca2.to_filename(opj(path_anat, ID + 'template_indiv' + otheranat + 'rsp_TO_Tnorm.nii.gz'))

        command = SINGULARITY + 'run' + AFNI_sif + '3dcalc -a ' +  opj(path_anat, ID + 'template_indiv' + type_norm + '.nii.gz') +  \
        ' -b ' + opj(path_anat, ID + 'template_indiv' + otheranat + 'rsp_TO_Tnorm.nii.gz') + \
        ' -expr b*(step(a)) -prefix ' + opj(path_anat, ID + 'template_indiv' + otheranat + 'rsp_TO_Tnorm.nii.gz') + \
        '-overwrite'
        spco([command], shell=True)
        
        
        log_img = math_img("np.where((img1 > 1) & (img2 > 1), img2, 0)", 
                           img1=opj(path_anat, ID + 'template_indiv' + type_norm + '.nii.gz'),
                           img2=opj(path_anat, ID + 'template_indiv' + otheranat + 'rsp_TO_Tnorm.nii.gz'))
        log_img.to_filename(opj(path_anat, ID + 'template_indiv' + otheranat + 'rsp_T1T2.nii.gz'))

        log_img = math_img("np.where((img1 > 1) & (img2 > 1), img1, 0)", 
                           img1=opj(path_anat, ID + 'template_indiv' + type_norm + '.nii.gz'), 
                           img2=opj(path_anat, ID + 'template_indiv' + otheranat + 'rsp_TO_Tnorm.nii.gz'))
        log_img.to_filename(opj(path_anat, ID + 'template_indiv' + type_norm + 'rsp_T1T2.nii.gz'))

        command = SINGULARITY + 'run' + AFNI_sif + '3dcalc -overwrite -a ' + opj(path_anat, ID + 'template_indiv' + type_norm + 'rsp_T1T2.nii.gz') + \
                  ' -b ' + opj(path_anat, ID + 'template_indiv' + otheranat + 'rsp_T1T2.nii.gz') + \
                  ' -expr "a/b" -datum float  -prefix ' + opj(path_anat, ID + 'template_indiv' + type_norm + '_' + otheranat + '.nii.gz')
        spco([command], shell=True)

    if IgotbothT1T2 == True:
        listTimage2 = list(listTimage)
        listTimage2.append(str(type_norm) + '_' + str(otheranat))
        print(listTimage2)
    else:
        listTimage2 = list(listTimage)
    print(listTimage2)
    print(listTimage)



    ###define if the two image have the same oblique
    if IgotbothT1T2 == True:
        command = SINGULARITY + 'run' + AFNI_sif + '3dinfo -same_obl ' + \
        opj(path_anat, ID + 'template_indiv' + otheranat + '.nii.gz') + ' ' + \
        opj(path_anat, ID + 'template_indiv' + type_norm + '.nii.gz')
        nx = spco([command], shell=True)
    else:
        nx="b'1\\n1\\n'" #code to say "same oblique"...


    debolique_spe = ID + 'ses-' + str(Session)
    if debolique_spe in deoblique_exeption1:
        deoblique_1='exeption1'
        print(deoblique_1)
    elif debolique_spe in deoblique_exeption2:
        deoblique_1='exeption2'
    else:
        deoblique_1=deoblique
        print(deoblique_1)

    for Timage in listTimage2:
        list_anat = sorted(glob.glob(opj(path_anat, '*_' + Timage + '.nii.gz')))
        ## add option for if oblique = 0 !! XX
        if deoblique_1=='header':
            command = SINGULARITY + 'run' + AFNI_sif + '3drefit -deoblique -NN' + overwrite + ' -orient ' + orientation + \
            ' ' + opj(path_anat, ID + 'template_indiv' + Timage + '.nii.gz')
            spco([command], shell=True)

            shutil.copyfile(opj(path_anat, ID + 'template_indiv' + Timage + '.nii.gz'), 
                            opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz'))
        
        #work
        elif deoblique_1=='WARP':
            if str(nx)=="b'1\\n1\\n'":
                print('same oblique')
                command = SINGULARITY + 'run' + AFNI_sif + '3drefit' + overwrite + ' -orient ' + orientation + \
                ' ' + opj(path_anat, ID + 'template_indiv' + Timage + '.nii.gz')
                spco([command], shell=True)
                # reorient the fiedls according to the json file
                command = SINGULARITY + 'run' + AFNI_sif + '3dWarp' + overwrite + \
                ' -deoblique -NN -prefix ' +  opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz') + \
                ' ' + opj(path_anat, ID + 'template_indiv' + Timage + '.nii.gz')
                spco([command], shell=True)

            else:
                print('oblique different')
                # reorient the fiedls according to the json file
                command = SINGULARITY + 'run' + AFNI_sif + '3dWarp' + overwrite + \
                ' -deoblique -NN -prefix ' +  opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz') + \
                ' ' + opj(path_anat, ID + 'template_indiv' + Timage + '.nii.gz')
                spco([command], shell=True)
                
                command = SINGULARITY + 'run' + AFNI_sif + '3drefit' + overwrite + ' -keepcen -orient ' + orientation + \
                ' -duporigin' + ' ' + opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz') + \
                ' ' + opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz')
                spco([command], shell=True)

        elif deoblique_1=='exeption1': #do nothing
            print('exeption1')
            shutil.copyfile(opj(path_anat, ID + 'template_indiv' + Timage + '.nii.gz'), 
                            opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz'))

        elif deoblique_1=='exeption2': #re-alineate
            print('exeption2') ### add something else not usefull!! XX
            command = '3dWarp' + overwrite + ' -deoblique -NN -prefix ' +  opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz') + \
            ' ' + opj(path_anat, ID + 'template_indiv' + Timage + '.nii.gz')
            spco([command], shell=True)
            print('need to realinate or use just one')

    for Timage in listTimage:
        # B1 Non uniform correction (MINC tools via freesurfer)
        command = SINGULARITY + 'run' + FS_sif + 'mri_nu_correct.mni --i ' + opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz') + \
        ' --o ' + opj(dir_prepro, ID + '_anat_reorient_NU' + Timage + '.nii.gz') + \
        ' --distance 60 --proto-iters 150 --stop 0.000'
        spco([command], shell=True)