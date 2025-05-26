import os
import subprocess
from nilearn.image import resample_to_img
from fonctions.extract_filename import extract_filename
from nilearn.image.image import mean_img
import fonctions.Skullstrip_func
import datetime
import json
import re
import nibabel as nib

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

opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

spco = subprocess.check_output
spgo = subprocess.getoutput

def Refimg_to_meanfMRI(anat_func_same_space, BASE_SS_coregistr,TfMRI , dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, use_master_for_Allineate,
                       dir_fMRI_Refth_RS_prepro3, RS, nb_run, REF_int, ID, dir_prepro, brainmask, V_mask, W_mask, G_mask, dilate_mask,
                       costAllin, anat_subject, Method_mask_func, overwrite, type_of_transform, aff_metric_ants,
                       s_bind,afni_sif,fs_sif, fsl_sif, itk_sif,diary_file):

    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  Working on step ' + str(3) + '(function: _3_mask_fMRI).  ##'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    ##### create new variable for template space (we will need to store and downsample template img to func resolution)
    if ope(dir_fMRI_Refth_RS_prepro3) == False:

        os.makedirs(dir_fMRI_Refth_RS_prepro3)
        os.makedirs(opj(dir_fMRI_Refth_RS_prepro3,'tmp'))

    ### create a list of the image to be corrected

    root_RS_ref = extract_filename(RS[REF_int])

    MEAN_im_list = ' ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS_ref + '_xdtr_mean_deob.nii.gz') #image "reference" (to be created)
    MEAN_im_list_1    = [opj(dir_fMRI_Refth_RS_prepro1,root_RS_ref + '_xdtr_mean_deob.nii.gz')] #image "reference" (to be created)

    for r in range(0, int(nb_run)):
        root_RS = extract_filename(RS[r])
        if not RS[r] == RS[REF_int]: # do not process ref...
            # add the co-registered mean image to the list
            MEAN_im_list = MEAN_im_list + ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref.nii.gz') #add images in the same space
            MEAN_im_list_1.append(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref.nii.gz'))

    ################################################################# create a mean image to use for the anat to func and recenter

    ###### average all func data and clean the image #####

    command = 'singularity run' + s_bind + afni_sif + '3dTcat' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'residual_motion.nii.gz') + MEAN_im_list
    nl = spgo(command)
    diary.write(f'\n{nl}')
    print(nl)
    dictionary = {"Sources": MEAN_im_list,
                  "Description": '4D concatenation (3dTcat,AFNI).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro1,'residual_motion.json'), "w") as outfile:
        outfile.write(json_object)

    #################################### production of Mean image ####################################

    mean_haxby = mean_img(MEAN_im_list_1)
    mean_haxby.to_filename(opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz'))
    dictionary = {"Sources": MEAN_im_list_1,
                  "Description": 'mean image (mean_image, nilearn).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.json'), "w") as outfile:
        outfile.write(json_object)

    mean_haxby.to_filename(opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_test.nii.gz'))
    dictionary = {"Sources": MEAN_im_list_1,
                  "Description": 'mean image (mean_image, nilearn).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_test.json'), "w") as outfile:
        outfile.write(json_object)

    # Load the image directly
    mean_img_path = os.path.join(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz')
    img = nib.load(mean_img_path)
    # Get voxel sizes
    delta_x, delta_y, delta_z = [str(round(abs(x), 10)) for x in img.header.get_zooms()[:3]]

    def get_orientation_nibabel(nifti_path):
        """Get 3-letter orientation code using NiBabel."""
        img = nib.load(nifti_path)
        aff = img.affine

        # Extract orientation from affine matrix
        ornt = nib.orientations.io_orientation(aff)
        codes = nib.orientations.ornt2axcodes(ornt)
        orient_code = ''.join(codes)

        # Validate (should already be valid from NiBabel)
        if not re.fullmatch(r'^[RLAPSI]{3}$', orient_code):
            raise ValueError(f"Invalid orientation: {orient_code}")
        return orient_code

    # Usage
    orient_meanimg = get_orientation_nibabel(mean_img_path)
    print(f"Orientation: {orient_meanimg}")

    command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + BASE_SS_coregistr +  \
    ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + ' -expr "a"'
    nl = spgo(command)
    diary.write(f'\n{nl}')
    print(nl)

    command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
    ' -orient ' + orient_meanimg + \
    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + \
    ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + \
    ' -rmode Cu -input  ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz')
    nl = spgo(command)
    diary.write(f'\n{nl}')
    print(nl)
    dictionary = {"Sources": [BASE_SS_coregistr,
                              opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')],
                  "Description": 'Resampling (3dresample, AFNI).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro3, 'BASE_SS_fMRI.json'), "w") as outfile:
        outfile.write(json_object)

    # create anat space dir
    if ope(dir_fMRI_Refth_RS_prepro2) == False:
        os.makedirs(dir_fMRI_Refth_RS_prepro2)
        os.makedirs(opj(dir_fMRI_Refth_RS_prepro2,'tmp'))

    ############################### ############################### ############################### 
    ##                         resample the masks for signal extraction                         ###
    ############################### ############################### ###############################

    ##### first you need to re-create a dilate ref anat image (works better for co-registration anat to fMRI)
    # dilate a little bit MORE the "maskDilat"

    if dilate_mask != 0:
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + \
        ' -input ' + brainmask + ' -fill_holes -dilate_input ' + str(dilate_mask)
    else:
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + \
        ' -input ' + brainmask + ' -fill_holes'

    nl = spgo(command)
    diary.write(f'\n{nl}')
    print(nl)
    dictionary = {"Sources": brainmask,
                  "Description": 'dilation of a factor ' + str(dilate_mask) + ' (3dmask_tool, nilearn).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro2, 'maskDilatanat.json'), "w") as outfile:
        outfile.write(json_object)

    if anat_func_same_space == True:
        # Original matrix path
        orig_mat = opj(dir_prepro, f"{ID}_brain_for_Align_Center.1D")
        # Inverse matrix path
        mvt_shft = opj(dir_prepro, f"{ID}_brain_for_Align_Center_inv.1D")

        command = (f'singularity run {s_bind} {afni_sif} '
            f'cat_matvec -ONELINE {orig_mat} '
            f'> {opj(dir_fMRI_Refth_RS_prepro2, "_brain_for_Align_Center.1D")}')
        nl = spgo(command)
        print(command)
        diary.write(f'\n{nl}')
        print(nl)

        # Generate inverse matrix with proper formatting
        command = f"""singularity run {s_bind} {afni_sif} \
        cat_matvec -ONELINE {orig_mat} -I > {mvt_shft}"""
        nl = spgo(command)
        print(command)
        diary.write(f'\n{nl}')
        print(nl)

        # move the atlases to the space before the AFNI shift
        for input1, output2 in zip([anat_subject, brainmask, opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz'), V_mask, W_mask, G_mask],
            [opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz'),
             opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz'),
            opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz')]):

            if ope(input1):
                ##### apply the recenter fmri
                command = 'singularity run' + s_bind + afni_sif + '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + output2 + \
                          ' ' + input1 + ' -overwrite'
                nl = spgo(command)
                diary.write(f'\n{nl}')
                print(nl)

                if use_master_for_Allineate==True:
                    ####  may be added if failed!! don't know how to handle that or why it fail!
                    command = 'singularity run' + s_bind + afni_sif + '3dAllineate -final NN' + overwrite + ' -overwrite -1Dmatrix_apply ' + mvt_shft + \
                    ' -prefix ' + output2 + \
                    ' -master ' + opj(dir_prepro, ID + '_mprage_reorient' + TfMRI + '.nii.gz') + \
                    ' -input  ' + output2
                    nl = spgo(command)
                    diary.write(f'\n{nl}')
                    print(nl)
                else:
                    command = 'singularity run' + s_bind + afni_sif + '3dAllineate -final NN' + overwrite + ' -overwrite -1Dmatrix_apply ' + mvt_shft + \
                              ' -prefix ' + output2 + \
                              ' -input  ' + output2
                    nl = spgo(command)
                    diary.write(f'\n{nl}')
                print(nl)

                ### to reapply the original obliquity
                caca2 = resample_to_img(output2, opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz'), interpolation='nearest')
                caca2.to_filename(output2)

                dictionary = {"Sources": [input1,
                                          opj(dir_prepro, ID + '_mprage_reorient' + TfMRI + '.nii.gz'),
                                          opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz')],
                              "Description": 'Co-registration and resampling (3dAllineate, AFNI; resample_to_img, nilearn).', }
                json_object = json.dumps(dictionary, indent=2)
                with open(output2[:-7] + '.json', "w") as outfile:
                    outfile.write(json_object)

            else:
                nl = 'WARNING:' + str(input1) + ' not found!!! this may be because you have not provided an aseg file,' + \
                     ' then no extraction of WM or Ventricles or GM will be possible... pls check that!'
                print(bcolors.WARNING + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

        # skullstrip the anat
        command = ('singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') +
                   ' -b ' + opj(dir_fMRI_Refth_RS_prepro2, 'orig_anat_for_fMRI.nii.gz') +
                  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.nii.gz') + ' -expr "a*b"')
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro2, 'orig_anat_for_fMRI.nii.gz'),
                                  opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz')],
                      "Description": 'Skull stripping (3dcalc, AFNI).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.json'), "w") as outfile:
            outfile.write(json_object)

    else:
        for input1, output2 in zip([anat_subject, brainmask, opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz'), V_mask, W_mask, G_mask],
            [opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz'), 
            opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz')]):

            if ope(input1):
                if input1 == anat_subject:
                    command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
                    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + \
                    ' -master ' + anat_subject + \
                    ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz')
                    nl = spgo(command)
                    diary.write(f'\n{nl}')
                    print(nl)

                    # skullstrip the anat
                    command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + \
                              ' -b ' + anat_subject + \
                              ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz') + ' -expr "a*b"'
                    nl = spgo(command)
                    diary.write(f'\n{nl}')
                    print(nl)

                    command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
                    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.nii.gz') + \
                    ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                    ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz')
                    nl = spgo(command)
                    diary.write(f'\n{nl}')
                    print(nl)
                    dictionary = {"Sources": [anat_subject,
                                              opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz'),
                                              opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')],
                                  "Description": 'Skull stripping and resampling (3dcalc and 3dresample, AFNI).', }
                    json_object = json.dumps(dictionary, indent=2)
                    with open(opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.json'), "w") as outfile:
                        outfile.write(json_object)

                else:
                    command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
                    ' -prefix ' + output2 + \
                    ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                    ' -input  ' + input1
                    nl = spgo(command)
                    diary.write(f'\n{nl}')
                    print(nl)
                    dictionary = {"Sources": [input1,
                                              opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')],
                                  "Description": 'resampling (3dresample, AFNI).', }
                    json_object = json.dumps(dictionary, indent=2)
                    with open(output2[:-7] + '.json', "w") as outfile:
                        outfile.write(json_object)
            else:
                print(bcolors.WARNING + 'WARNING:'
                      + str(input1) + ' not found!!! this may be because you have not provided an aseg file, then no '
                                      'extraction of WM or Ventricles or GM will be possible... please check that!' + bcolors.ENDC)

    nl1 = "INFO: brain_skullstrip method is " + Method_mask_func
    nl2 = 'INFO: looking for manual segmentation named:' + opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz') + '...'
    print(bcolors.OKGREEN + nl1 + bcolors.ENDC)
    diary.write(f'\n{nl1}')
    print(bcolors.OKGREEN + nl2 + bcolors.ENDC)
    diary.write(f'\n{nl2}')

    #### explore if manual_mask.nii.gz exists?
    if ope(opj(dir_fMRI_Refth_RS_prepro1,'manual_mask.nii.gz')):

        nl1 = 'WARNING: We found a final mask to skullstrip the functional image !!! no Skullstrip will be calculated!'
        nl2 = 'INFO: please delete' +opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz') + ' if you want retry to create a skulstripp images'
        print(bcolors.WARNING + nl1 + bcolors.ENDC)
        diary.write(f'\n{nl1}')
        print(bcolors.OKGREEN + nl2 + bcolors.ENDC)
        diary.write(f'\n{nl2}')

        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1,'manual_mask.nii.gz') + \
        ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1,'maskDilat_Allineate_in_func.nii.gz') + ' -expr "a"'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1,'manual_mask.nii.gz'),
                      "Description": 'Copy.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1,'maskDilat_Allineate_in_func.json'), "w") as outfile:
            outfile.write(json_object)

    #### if not, create an fMRI mask
    else:
        nl1 = "INFO: no manual mask found "
        print(bcolors.OKGREEN + nl1 + bcolors.ENDC)
        if anat_func_same_space:
            command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro2, 'maskDilat.nii.gz') + \
                      ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.nii.gz') + ' -expr "a"'
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro2, 'maskDilat.nii.gz'),
                          "Description": 'Copy.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.json'), "w") as outfile:
                outfile.write(json_object)

            nl = 'you are using the mask from the anat img'
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')
            diary.write(f'\n')
            diary.close()
        else:
            diary.write(f'\n')
            diary.close()
            fonctions.Skullstrip_func.Skullstrip_func(Method_mask_func, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                                                      overwrite, costAllin, type_of_transform,
                                                      aff_metric_ants, s_bind, afni_sif, fsl_sif, fs_sif, itk_sif, diary_file)
