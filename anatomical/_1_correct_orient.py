import os
import subprocess
import shutil
import glob
from fonctions.extract_filename import extract_filename
import json
from nilearn.image import resample_to_img
from nilearn.image import math_img
import ants
import numpy as np
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
                   deoblique, 
                   orientation, 
                   dir_prepro, 
                   IgotbothT1T2, 
                   overwrite,
                   s_bind,afni_sif,fs_sif,
                   diary_file):

    ct = datetime.datetime.now()
    nl = 'Run  anatomical._1_correct_orient.correct_orient'
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{nl}')

    for Timage in listTimage:
        # get variables from json
        if BIDStype == 1:
            list_anat = sorted(glob.glob(opj(path_anat, 'sub-' + ID + '_ses-' + str(Session) + '_*' + Timage + '.nii*')))
        elif BIDStype == 2:
            list_anat = sorted(glob.glob(opj(path_anat, 'sub-' + ID + '_' + Timage + '.nii*')))
        else:
            nl = 'WARNING: BIDStype was not 1 or 2 so it must be define by a string that you provided'
            print(bcolors.WARNING + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')
            list_anat = sorted(glob.glob(opj(path_anat, BIDStype.format(ID=ID, Session=Session, Timage=Timage))))

        if len(list_anat)==1:
            nl = 'INFO: We found only one anat images for this session'
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')
            nl = 'list of the anat found is:' + str(list_anat) + bcolors.ENDC
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

            command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + list_anat[0] + \
            ' -prefix ' + opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz') + ' -expr "a"' + overwrite
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": list_anat[0],
                          "Description": 'Nothing done.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_prepro, ID + 'template_indiv' + Timage + '.json'), "w") as outfile:
                outfile.write(json_object)

            ###################################
            ###########remove norm anat #########
            ###################################
        elif len(list_anat)>1:
            nl = 'INFO: We found ' + str(len(list_anat)) + ' anat images for this session'
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')
            nl = 'list of the anat found is:' + str(list_anat)
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

            for i, anat in enumerate(list_anat):
                if ope(opj(opd(list_anat[i]), extract_filename(list_anat[i]) + '.json')):
                    read_json = opj(opd(list_anat[i]), extract_filename(list_anat[i]) + '.json')
                else:
                    nl = "WARNING: No .json associated to the the anat image !! you might want to check that!"
                    print(bcolors.WARNING + nl  + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    ff = open(read_json)
                    anat_T1 = json.load(ff)
                    try:
                        ImageType = anat_T1["ImageType"]
                        # Compare the contents of the two files
                        if 'NORM' in ImageType:
                            nl = 'Removing ' + extract_filename(list_anat[i]) + ' as it is a NORM'
                            print(bcolors.WARNING + nl + bcolors.ENDC)
                            diary.write(f'\n{nl}')

                            list_anat.pop(i)
                            # After removing the file, no need to increment i since the next file takes the place of the removed one
                    except:
                        nl = "WARNING No ImageType, could not check if NORM image was included in your BIDS, you might want to check that!"
                        print(bcolors.WARNING + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')

            if len(list_anat) > 1:
                ANAT= ' '.join(list_anat)
                command = 'singularity run' + s_bind + fs_sif + 'mri_robust_template --mov ' + ANAT + \
                ' --template ' + opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz') + \
                ' --satit --inittp 1 --fixtp --noit --iscale --average 0'
                nl= spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
                dictionary = {"Sources": ANAT,
                              "Description": 'Co-registation and Average (mri_robust_template from Freesurfer)) .', }
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(dir_prepro, ID + 'template_indiv' + Timage + '.json'), "w") as outfile:
                    outfile.write(json_object)
            else:
                #shutil.copyfile(list_anat[0], opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz'))
                command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + list_anat[0] + \
                ' -prefix ' + opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz') + ' -expr "a"' + overwrite
                nl = spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
                dictionary = {"Sources": list_anat[0],
                              "Description": 'Nothing done.', }
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(dir_prepro, ID + 'template_indiv' + Timage + '.json'), "w") as outfile:
                    outfile.write(json_object)

        else:
            nl = 'HELP: anat is detected with BIDStype=' + str(BIDStype) + ' the command: '
            print(bcolors.WARNING + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')
            if BIDStype == 1:
                list_anat = "glob.glob(" + path_anat + '/sub-' + ID + '_ses-' + str(Session) + '_*' + Timage + '.nii*"'
            elif BIDStype == 2:
                list_anat = '"glob.glob(' + path_anat + '/sub-' + ID + '_' + Timage + '.nii*"'
            else:
                list_anat = "glob.glob(" + str(opj(path_anat, BIDStype.format(ID=ID, Session=Session, Timage=Timage)))
            nl = str(list_anat) + ' but failed'
            diary.write(nl)
            raise Exception(bcolors.FAIL + nl + bcolors.ENDC)


        ####################################################
        ####### create indiv T1/T2 template if necessary ####
        ####################################################

    if IgotbothT1T2 == True:
        caca2 = resample_to_img(opj(dir_prepro, ID + 'template_indiv' + otheranat + '.nii.gz'),
                                opj(dir_prepro, ID + 'template_indiv' + type_norm + '.nii.gz'),
                                interpolation='nearest')
        caca2.to_filename(opj(dir_prepro, ID + 'template_indiv' + otheranat + 'rsp_TO_Tnorm.nii.gz'))
        dictionary = {"Sources": [opj(dir_prepro, ID + 'template_indiv' + otheranat + '.nii.gz'),
                                  opj(dir_prepro, ID + 'template_indiv' + type_norm + '.nii.gz')],
                      "Description": 'Resample .', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_prepro, ID + 'template_indiv' + otheranat + 'rsp_TO_Tnorm.json'), "w") as outfile:
            outfile.write(json_object)

        command = 'singularity run' + s_bind + afni_sif + \
                  '3dcopy ' + opj(dir_prepro, ID + 'template_indiv' + type_norm + '.nii.gz') + ' ' + \
                  opj(dir_prepro, ID + 'template_indiv' + type_norm + 'rsp_TO_Tnorm.nii.gz') + overwrite
        nl= spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": opj(dir_prepro, ID + 'template_indiv' + type_norm + '.nii.gz'),
                      "Description": 'Copy .', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_prepro, ID + 'template_indiv' + type_norm  + 'rsp_TO_Tnorm.json'), "w") as outfile:
            outfile.write(json_object)

        for Timage in listTimage:
            BRAIN = ants.image_read(opj(dir_prepro, ID + 'template_indiv' + Timage + 'rsp_TO_Tnorm.nii.gz'))
            image_shape = BRAIN.shape
            if len(image_shape) == 3:
                # The image is 3D, no need for averaging
                print("The input image is 3D. No averaging needed.")
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
            elif len(image_shape) == 4:
                # The image is 4D (multi-echo or time series), average across the 4th dimension
                nl = f"The input image is 4D with {image_shape[3]} volumes. Averaging across the 4th dimension."
                print(nl)
                diary.write(f'\n{nl}')

                BRAINnumpy = BRAIN.numpy()
                averaged_array = np.mean(BRAINnumpy, axis=3)
                averaged_image = ants.from_numpy(
                    averaged_array,
                    origin=BRAIN.origin[:3],
                    spacing=BRAIN.spacing[:3],
                    direction=np.eye(3)  # Set direction matrix to identity to avoid dimensionality issues
                )
                nl = 'Averaged 4D image to 3D'
                print(nl)
                diary.write(f'\n{nl}')

                ants.image_write(averaged_image, opj(dir_prepro, ID + 'template_indiv' + Timage + 'rsp_TO_Tnorm.nii.gz'), ri=False)

                dictionary = {"Sources": opj(dir_prepro, ID + 'template_indiv' + Timage + 'rsp_TO_Tnorm.nii.gz'),
                              "Description": '4D averaging .', }
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(dir_prepro, ID + 'template_indiv' + Timage + 'rsp_TO_Tnorm.json'), "w") as outfile:
                    outfile.write(json_object)

        #### if T2 or T1 not mentionned in the name then it's a mess... XXX
        for Timage in listTimage:
            if 'T1' in Timage:
                T1 = Timage
            elif 'T2' in Timage:
                T2 = Timage

        command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' +  opj(dir_prepro, ID + 'template_indiv' + T1 + 'rsp_TO_Tnorm.nii.gz') +  \
        ' -b ' + opj(dir_prepro, ID + 'template_indiv' + T2 + 'rsp_TO_Tnorm.nii.gz') + \
        ' -expr "b*(step(a))" -prefix ' + opj(dir_prepro, ID + 'template_indiv' + T2 + 'rsp_TO_Tnorm.nii.gz') + \
        ' -overwrite'
        nl= spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": [opj(dir_prepro, ID + 'template_indiv' + T1 + 'rsp_TO_Tnorm.nii.gz'),
                                  opj(dir_prepro, ID + 'template_indiv' + T2 + 'rsp_TO_Tnorm.nii.gz')],
                      "Description": 'T2 signal "correction" .', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_prepro, ID + 'template_indiv' + T2+ 'rsp_TO_Tnorm.json'), "w") as outfile:
            outfile.write(json_object)
        
        
        log_img = math_img("np.where((img1 > 1) & (img2 > 1), img2, 0)", 
                           img1=opj(dir_prepro, ID + 'template_indiv' + T1 + 'rsp_TO_Tnorm.nii.gz'),
                           img2=opj(dir_prepro, ID + 'template_indiv' + T2 + 'rsp_TO_Tnorm.nii.gz'))
        log_img.to_filename(opj(dir_prepro, ID + 'template_indiv' + otheranat + 'rsp_T1T2.nii.gz'))

        dictionary = {"Sources": [opj(dir_prepro, ID + 'template_indiv' + T1 + 'rsp_TO_Tnorm.nii.gz'),
                                  opj(dir_prepro, ID + 'template_indiv' + T2 + 'rsp_TO_Tnorm.nii.gz')],
                      "Description": 'T2 signal "selection" .', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_prepro, ID + 'template_indiv' + otheranat + 'rsp_T1T2.json'), "w") as outfile:
            outfile.write(json_object)

        log_img = math_img("np.where((img1 > 1) & (img2 > 1), img1, 0)", 
                           img1=opj(dir_prepro, ID + 'template_indiv' + T1 + 'rsp_TO_Tnorm.nii.gz'),
                           img2=opj(dir_prepro, ID + 'template_indiv' + T2 + 'rsp_TO_Tnorm.nii.gz'))
        log_img.to_filename(opj(dir_prepro, ID + 'template_indiv' + T1 + 'rsp_T1T2.nii.gz'))

        dictionary = {"Sources": [opj(dir_prepro, ID + 'template_indiv' + T1 + 'rsp_TO_Tnorm.nii.gz'),
                                  opj(dir_prepro, ID + 'template_indiv' + T2 + 'rsp_TO_Tnorm.nii.gz')],
                      "Description": 'T1 signal "selection" .', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_prepro, ID + 'template_indiv' + T1 + 'rsp_T1T2.json'), "w") as outfile:
            outfile.write(json_object)

        command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + opj(dir_prepro, ID + 'template_indiv' + T1 + 'rsp_T1T2.nii.gz') + \
                  ' -b ' + opj(dir_prepro, ID + 'template_indiv' + T2 + 'rsp_T1T2.nii.gz') + \
                  ' -expr "a/b" -datum float  -prefix ' + opj(dir_prepro, ID + 'template_indiv' + T1 + '_' + T2 + '.nii.gz')
        nl= spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": [opj(dir_prepro, ID + 'template_indiv' + T1 + 'rsp_T1T2.nii.gz'),
                                  opj(dir_prepro, ID + 'template_indiv' + T2 + 'rsp_T1T2.nii.gz')],
                      "Description": 'T1 divided by T2 .', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_prepro, ID + 'template_indiv' + T1 + '_' + T2 + '.json'), "w") as outfile:
            outfile.write(json_object)


    if IgotbothT1T2 == True:
        listTimage2 = list(listTimage)
        listTimage2.append(str(type_norm) + '_' + str(otheranat))
    else:
        listTimage2 = list(listTimage)

    ###define if the two image have the same oblique
    if IgotbothT1T2 == True:
        command = 'singularity run' + s_bind + afni_sif + '3dinfo -same_obl ' + \
        opj(dir_prepro, ID + 'template_indiv' + otheranat + '.nii.gz') + ' ' + \
        opj(dir_prepro, ID + 'template_indiv' + type_norm + '.nii.gz')
        arg = spgo(command).split('\n')
        nx = arg[-1]
    else:
        nx="b'1\\n1\\n'" #code to say "same oblique"...

    for Timage in listTimage2:
        ## add option for if oblique = 0 !! XX
        if deoblique=='header':
            command = 'singularity run' + s_bind + afni_sif + '3drefit -deoblique' + overwrite + ' -orient ' + orientation + \
            ' ' + opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz')
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

            shutil.copyfile(opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz'), 
                            opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz'))

            dictionary = {"Sources": opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz'),
                          "Description": 'Correction of the fields orientation .', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.json'), "w") as outfile:
                outfile.write(json_object)

        if deoblique=='header_WO_deob':
            command = 'singularity run' + s_bind + afni_sif + '3drefit ' + overwrite + ' -orient ' + orientation + \
            ' ' + opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz')
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

            shutil.copyfile(opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz'),
                            opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz'))

            dictionary = {"Sources": opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz'),
                          "Description": 'Correction of the fields orientation .', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.json'), "w") as outfile:
                outfile.write(json_object)

        elif deoblique == 'WARPbaboon':
            command = 'singularity run' + s_bind + afni_sif + '3drefit' + overwrite + ' -duporigin ' + list_anat[0] + ' -orient ' + orientation + \
                      ' ' + opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz')
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

            # reorient the fields according to the json file
            command = 'singularity run' + s_bind + afni_sif + '3dWarp ' + overwrite + \
                      ' -deoblique -wsinc5 -prefix ' + opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz') + \
                      ' ' + opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz')
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

            dictionary = {"Sources": opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz'),
                          "Description": 'Correction of the fields orientation .', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.json'), "w") as outfile:
                outfile.write(json_object)

        #work
        elif deoblique=='WARP':
            if str(nx)=="b'1\\n1\\n'":
                command = 'singularity run' + s_bind + afni_sif + '3drefit' + overwrite + ' -orient ' + orientation + \
                ' ' + opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz')
                nl= spgo(command)
                diary.write(f'\n{nl}')
                print(nl)

                # reorient the fields according to the json file
                command = 'singularity run' + s_bind + afni_sif + '3dWarp' + overwrite + \
                ' -deoblique -NN -prefix ' +  opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz') + \
                ' ' + opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz')
                nl = spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
                dictionary = {"Sources": opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz'),
                              "Description": 'Correction of the fields orientation .', }
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.json'), "w") as outfile:
                    outfile.write(json_object)

                #command = 'singularity run' + s_bind + afni_sif + '3drefit' + overwrite + ' -deoblique -orient ' + orientation + \
                #' ' + opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz')
                #nl= spgo(command)
                #diary.write(f'\n{nl}')
                #print(nl)

            else:
                # reorient the fields according to the json file
                command = 'singularity run' + s_bind + afni_sif + '3dWarp' + overwrite + \
                ' -deoblique -NN -prefix ' +  opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz') + \
                ' ' + opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz')
                nl= spgo(command)
                diary.write(f'\n{nl}')
                print(nl)

                command = 'singularity run' + s_bind + afni_sif + '3drefit' + overwrite + ' -keepcen -orient ' + orientation + \
                ' -duporigin' + ' ' + opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz') + \
                ' ' + opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz')
                nl= spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
                dictionary = {"Sources": opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz'),
                              "Description": 'Correction of the fields orientation .', }
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.json'), "w") as outfile:
                    outfile.write(json_object)

        elif deoblique=='no_deoblique': #do nothing
            nl='exeption1'
            print(nl)
            diary.write(f'\n{nl}')
            shutil.copyfile(opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz'), 
                            opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz'))
            dictionary = {"Sources": opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz'),
                          "Description": 'Copy .', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.json'), "w") as outfile:
                outfile.write(json_object)

        elif deoblique == 'deob_WO_orient':  # do nothing
            command = 'singularity run' + s_bind + afni_sif + '3drefit -deoblique' + overwrite + \
            ' ' + opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz')
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

            shutil.copyfile(opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz'),
                            opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz'))

            dictionary = {"Sources": opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz'),
                          "Description": 'Correction of the fields orientation .', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.json'), "w") as outfile:
                outfile.write(json_object)

        elif deoblique=='WARP_without_3drefit': #re-alineate
            nl = 'exeption2'
            print(nl) ### add something else not usefull!! XX
            diary.write(f'\n{nl}')
            command = 'singularity run' + s_bind + afni_sif + '3dWarp' + overwrite + ' -deoblique -NN -prefix ' +  opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz') + \
            ' ' + opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz')
            nl= spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            nl = 'need to realign or use just one'
            print(nl)
            diary.write(f'\n{nl}')
            dictionary = {"Sources": opj(dir_prepro, ID + 'template_indiv' + Timage + '.nii.gz'),
                         "Description": 'deoblique .', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.json'), "w") as outfile:
                outfile.write(json_object)

    for Timage in listTimage:
        ###N4BiasFieldCorrection
        BRAIN = ants.image_read(opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz'))
        image_shape = BRAIN.shape

        # Apply N4 bias field correction
        N4 = ants.n4_bias_field_correction(BRAIN, shrink_factor=2, convergence={'iters': [50, 50, 30, 20], 'tol': 1e-7})
        nl = 'N4 Bias Field Correction done'
        print(nl)
        diary.write(f'\n{nl}')

        if len(image_shape) == 3:
            # The image is 3D, no need for averaging
            nl = "The input image is 3D. No averaging needed."
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

            averaged_image = N4
        elif len(image_shape) == 4:
            # The image is 4D (multi-echo or time series), average across the 4th dimension

            nl = f"The input image is 4D with {image_shape[3]} volumes. Averaging across the 4th dimension."
            print(nl)
            diary.write(f'\n{nl}')

            # Convert the 4D image to a NumPy array
            n4_array = N4.numpy()

            # Average across the 4th dimension
            averaged_array = np.mean(n4_array, axis=3)

            # Convert the averaged NumPy array back to an ANTsImage
            # Use metadata from the original image but adjust for 3D
            averaged_image = ants.from_numpy(
                averaged_array,
                origin=N4.origin[:3],
                spacing=N4.spacing[:3],
                direction=np.eye(3)  # Set direction matrix to identity to avoid dimensionality issues
            )
            nl = 'Averaged 4D image to 3D'
            print(nl)
            diary.write(f'\n{nl}')

        else:
            nl = f"Unsupported image dimensions: {image_shape}"
            diary.write(f'\n{nl}')
            raise ValueError(bcolors.FAIL + nl + bcolors.ENDC)

        # Write the output image
        ants.image_write(averaged_image, opj(dir_prepro, ID + '_anat_reorient_NU' + Timage + '.nii.gz'), ri=False)
        dictionary = {"Sources": opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz'),
                      "Description": 'B1 bias field non uniformity correction (N4 from ANTspy).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_prepro, ID + '_anat_reorient_NU' + Timage + '.json'), "w") as outfile:
            outfile.write(json_object)

    diary.write(f'\n')
    diary.close()


