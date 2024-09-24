import os
import numpy as np
import subprocess
import ants
from nilearn import plotting
from nilearn.image import resample_to_img
import anatomical.Skullstrip_method
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

#################################################
########    creat brain image of animal  ########
#################################################


def clean_anat(Align_img_to_template, cost3dAllineate, bids_dir, listTimage, type_of_transform, ID, aff_metric_ants, Session, otheranat, type_norm, dir_prepro, masking_img, do_manual_crop,
    brain_skullstrip_1, brain_skullstrip_2, masks_dir, volumes_dir, dir_transfo, BASE_SS_coregistr, BASE_SS_mask, BASE_SS, IgotbothT1T2, check_visualy_each_img, check_visualy_final_mask, template_skullstrip, study_template_atlas_forlder, overwrite,
               s_bind,afni_sif,fsl_sif,fs_sif, itk_sif, strip_sif):

    step_skullstrip = 1
    output_for_mask =  anatomical.Skullstrip_method.Skullstrip_method(step_skullstrip, template_skullstrip, study_template_atlas_forlder, masking_img, brain_skullstrip_1, brain_skullstrip_2, masks_dir, volumes_dir, dir_prepro, type_norm, BASE_SS_coregistr, BASE_SS_mask,
    type_of_transform, ID, aff_metric_ants, check_visualy_final_mask, s_bind, afni_sif, fsl_sif, fs_sif, itk_sif, strip_sif)

    #### Apply masking to other anat images

    for Timage in listTimage:
        maskRS = opj(masks_dir, ID + Timage + '_mask_1_rsp.nii.gz')
        caca2 = resample_to_img(output_for_mask, opj(dir_prepro, ID + '_anat_reorient_NU' + Timage + '.nii.gz'), interpolation='nearest')
        caca2.to_filename(maskRS)
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + opj(dir_prepro, ID + '_anat_reorient_NU' + Timage + '.nii.gz') + \
        ' -b ' + maskRS + \
        ' -expr "a*b" -prefix ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + Timage + '.nii.gz')
        spco([command], shell=True)

        #QC
        if not os.path.exists(bids_dir + '/QC/skullstrip_step1'):
            os.mkdir(bids_dir + '/QC/skullstrip_step1')
        try:
            display = plotting.plot_anat(opj(dir_prepro, ID + '_anat_reorient_NU' + Timage + '.nii.gz'), display_mode='mosaic')
            display.add_contours(opj(dir_prepro, ID + '_brain_for_Align_Center' + Timage + '.nii.gz'),
            linewidths=.3, colors=['red'])
            display.savefig(bids_dir + '/QC/skullstrip_step1/' + ID + '_' + str(Session) + '_' + Timage + '_skullstriped.png')
            # Don't forget to close the display
            display.close()
        except:
            display = plotting.plot_anat(opj(dir_prepro, ID + '_brain_for_Align_Center' + Timage + '.nii.gz'), display_mode='x', cut_coords=10)
            display.savefig(bids_dir + '/QC/skullstrip_step1/' + ID + '_' + str(Session) + '_' + Timage + '_skullstriped.png')
            # Don't forget to close the display
            display.close()

        ####################################################################################
        ########################## transfo rigid to atlas template (BASE_SS)   #############
        ####################################################################################
        command = 'singularity run' + s_bind + afni_sif + '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + Timage + '.nii.gz') + ' ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + Timage + '.nii.gz') + ' -overwrite'
        spco([command], shell=True)

        command = 'singularity run' + s_bind + afni_sif + '3dresample -master ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + Timage + '.nii.gz') + ' -prefix ' + opj(masks_dir, ID + Timage + '_mask_1_rsp.nii.gz') + \
                  ' -input ' + opj(masks_dir, ID + Timage + '_mask_1_rsp.nii.gz') + ' -overwrite'
        spco([command], shell=True)

    #### align anat img to atlas template
    if Align_img_to_template == '3dAllineate':
        command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -warp shift_rotate -cmass -base ' + BASE_SS + \
        ' -cost ' + cost3dAllineate + ' -interp NN' + \
        ' -source_mask ' + opj(masks_dir, ID + type_norm + '_mask_1_rsp.nii.gz') + \
        ' -prefix ' + opj(dir_prepro, ID + '_acpc_64_orig_3dAllineate' + type_norm + '.nii.gz') + \
        ' -source ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + type_norm + '.nii.gz') + ' -1Dmatrix_save ' + \
        opj(dir_prepro,ID + '_brain_for_Align_Center.1D')
        spco(command, shell=True)

    elif Align_img_to_template == '@Align_Centers':
        current_working_directory = os.getcwd()
        os.chdir(dir_prepro)
        command = 'singularity run' + s_bind + afni_sif + '@Align_Centers -base ' + BASE_SS + \
        ' -dset ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + type_norm + '.nii.gz') + ' -cm -prefix ' + ID + '_brain_for_Align_Center.nii.gz' + overwrite
        spco([command], shell=True)
        os.chdir(str(current_working_directory))

        command = 'singularity run' + s_bind + afni_sif + '3dcopy ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + type_norm + '.nii.gz') + ' ' + opj(dir_prepro, ID + '_acpc_64_orig_3dAllineate' + type_norm + '.nii.gz') + overwrite
        spco([command], shell=True)

    elif Align_img_to_template == 'No':
        command = 'singularity run' + s_bind + afni_sif + '3dcopy ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + type_norm + '.nii.gz') + ' ' + opj(dir_prepro, ID + '_acpc_64_orig_3dAllineate' + type_norm + '.nii.gz') + overwrite
        spco([command], shell=True)
        def create_1D_matrix(filename, matrix):
            if matrix.shape != (3, 4):
                raise Exception(bcolors.FAIL + 'ERROR: Matrix must be 3x4 for a valid AFNI transformation.' + bcolors.ENDC)
            # Flatten the matrix to 1D and save it
            flat_matrix = matrix.flatten()
            np.savetxt(filename, [flat_matrix], fmt='%.6f')

        # Example usage
        filename = '/srv/projects/easymribrain/data/MRI/Bat/BIDS_bat/sub-1/ses-1/anat/native/01_preprocess/1_brain_for_Align_Center.1D'
        matrix = np.array([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0]])
        create_1D_matrix(filename, matrix)
        print(bcolors.OKGREEN + f"INFO: .1D matrix file saved as {filename}" + bcolors.ENDC)
    else:
        raise Exception(bcolors.FAIL + 'ERROR: Align_img_to_template need to be define as string (3dAllineate, @Align_Centers, No)' + bcolors.ENDC)


    command = 'export AFNI_NIFTI_TYPE_WARN=NO'
    spco([command], shell=True)

    command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -di ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + type_norm + '.nii.gz')
    dummy = spgo(command).split('\n')
    delta_x = str(round(abs(float(dummy[-1])), 2))
    command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -dj ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + type_norm + '.nii.gz')
    dummy = spgo(command).split('\n')
    delta_y = str(round(abs(float(dummy[-1])), 2))
    command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -dk ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + type_norm + '.nii.gz')
    dummy = spgo(command).split('\n')
    delta_z = str(round(abs(float(dummy[-1])), 2))

    command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
              ' -prefix ' + opj(dir_prepro, ID + '_acpc_64_orig_3dAllineate' + type_norm + '.nii.gz') + \
              ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
              ' -input  ' + opj(dir_prepro, ID + '_acpc_64_orig_3dAllineate' + type_norm + '.nii.gz')
    spco([command], shell=True)

    '''
    spco(['singularity run' + s_bind + afni_sif + '3dresample', '-master', opj(dir_prepro, ID + '_brain_for_Align_Center' + type_norm + '.nii.gz'), '-prefix', opj(dir_prepro, ID + '_acpc_64_orig_3dAllineate' + type_norm + '.nii.gz'),
          '-input', opj(dir_prepro, ID + '_acpc_64_orig_3dAllineate' + type_norm + '.nii.gz'), '-overwrite'])

    caca2 = resample_to_img(opj(dir_prepro, ID + '_acpc_64_orig_3dAllineate' + type_norm + '.nii.gz'),
                            opj(dir_prepro, ID + '_brain_for_Align_Center' + type_norm + '.nii.gz'),
                            interpolation='nearest')
    
    caca2.to_filename(opj(dir_prepro, ID + '_acpc_64_orig_3dAllineate' + type_norm + '.nii.gz'))
    '''

        ##########################################
        ########## correct anat images ###########
        ##########################################

    for Timage in listTimage:
        ####### optional manual acpc center #######
        if check_visualy_each_img == True:
            #Check the result and do manually the finest correction save the file as "ID_acpc_tmp.nii.gz"
            command = 'singularity run' + s_bind + fs_sif + 'freeview -v ' + BASE_SS + ' ' + opj(dir_prepro, ID + '_acpc_64' + Timage + '.nii.gz') + ':opacity=0.6:visible=1'
            spco([command], shell=True)

    if IgotbothT1T2 == True:
        listTimage2 = list(listTimage)
        listTimage2.append(str(type_norm) + '_' + str(otheranat))
    else:
        listTimage2 = list(listTimage)

    ##### Apply transfo to other anat images
    for enu, Timage in enumerate(listTimage2):
        if enu<1:
            command = 'singularity run' + s_bind + afni_sif + '3dAllineate -overwrite -interp NN -1Dmatrix_apply ' + opj(dir_prepro,ID + '_brain_for_Align_Center.1D') + \
            ' -prefix ' + opj(dir_prepro, ID + '_acpc_64' + Timage + '.nii.gz') + \
            ' -master ' + opj(dir_prepro, ID + '_acpc_64_orig_3dAllineate' + type_norm + '.nii.gz') + \
            ' -input  ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + Timage + '.nii.gz')
            spco([command], shell=True)
        else:
            command = 'singularity run' + s_bind + afni_sif + '3dAllineate -overwrite -interp NN -1Dmatrix_apply ' + opj(dir_prepro,ID + '_brain_for_Align_Center.1D') + \
            ' -prefix ' + opj(dir_prepro, ID + '_acpc_64' + str(type_norm) + '_' + str(otheranat) + '.nii.gz') + \
            ' -master ' + opj(dir_prepro, ID + '_acpc_64_orig_3dAllineate' + type_norm + '.nii.gz') + \
            ' -input  ' + opj(dir_prepro, ID + '_mprage_reorient' + str(type_norm) + '_' + str(otheranat) + '.nii.gz')
            spco([command], shell=True)

        ####### optional crop !! save as "ID_acpc_cropped.nii.gz" #######
        if do_manual_crop == True:
            command = 'singularity run' + s_bind + fs_sif + 'freeview -v ' + opj(dir_prepro, ID + '_acpc' + Timage + '.nii.gz')
        else:
            command = 'singularity run' + s_bind + afni_sif + '3dAutobox' + overwrite + ' -input ' + opj(dir_prepro, ID + '_acpc_64' + Timage + '.nii.gz') + ' -prefix ' + opj(dir_prepro, ID + '_acpc_cropped' + Timage + '.nii.gz') + ' -noclust  -overwrite'
            spco(command, shell=True)

        ######################################################################
        ####### B0 correction ####### ==> creat INDIV template !!!!!!!!!!!####
        ######################################################################

    for Timage in listTimage:
        IMG = ants.image_read(opj(dir_prepro,ID + '_acpc_cropped' + Timage + '.nii.gz'))
        IMG = ants.denoise_image(IMG,r=3,noise_model='Gaussian')
        ants.image_write(IMG, opj(volumes_dir, ID + '_' + Timage + '_template.nii.gz'), ri=False)
        if not os.path.exists(bids_dir + '/QC/align_rigid_to_template'):
            os.mkdir(bids_dir + '/QC/align_rigid_to_template')
        try:
            display = plotting.plot_anat(opj(volumes_dir, ID + '_' + Timage + '_template.nii.gz'), display_mode='mosaic', dim=4)
            display.add_contours(BASE_SS_coregistr,
            linewidths=.2, colors=['red'])
            display.savefig(bids_dir + '/QC/align_rigid_to_template/' + ID + '_' + str(Session) + '_' + Timage + '_align_rigid_to_template.png')
            # Don't forget to close the display
            display.close()
        except:
            display = plotting.plot_anat(opj(volumes_dir, ID + '_' + Timage + '_template.nii.gz'), display_mode='mosaic', dim=4)
            display.savefig(bids_dir + '/QC/align_rigid_to_template/' + ID + '_' + str(Session) + '_' + Timage + '_align_rigid_to_template.png')
            # Don't forget to close the display
            display.close()

        ####### check visualy indiv template #######
        if check_visualy_each_img == True:
            # now look at the coordinate of the brain and ajust bet2 according to them
            command = 'singularity run' + s_bind + fs_sif + 'freeview -v ' + opj(volumes_dir, ID + '_' + Timage + '_template.nii.gz')
            spco([command], shell=True)





