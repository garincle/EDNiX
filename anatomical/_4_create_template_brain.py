#################################################
########    creat brain image of animal  ########
#################################################
#co-register linear to the new sty template...
##############################################################################
####### CREATE THE STUDY TEMPLATE (IF YOU WANT ON) ###########################
##############################################################################
import os
import nibabel as nib
import subprocess
from nilearn import plotting
import anatomical.Skullstrip_method

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

def create_indiv_template_brain(dir_prepro, ID, Session, listTimage, volumes_dir, masking_img, brain_skullstrip_1, brain_skullstrip_2, masks_dir, type_norm, n_for_ANTS, dir_transfo, BASE_SS_coregistr, BASE_SS_mask, otheranat,
    check_visualy_final_mask, useT1T2_for_coregis, bids_dir, overwrite, BASE_bet, s_bind,afni_sif,fsl_sif,fs_sif, itk_sif):

    if ope(opj(masks_dir, ID + '_finalmask.nii.gz')):
        output_for_mask = opj(masks_dir, ID + '_finalmask.nii.gz')

    else:
        step_skullstrip = 2
        brain_skullstrip = 'brain_skullstrip_2'
        output_for_mask =  anatomical.Skullstrip_method.Skullstrip_method(step_skullstrip, brain_skullstrip, masking_img, brain_skullstrip_1, brain_skullstrip_2, masks_dir, volumes_dir, dir_prepro, type_norm, n_for_ANTS, dir_transfo, BASE_SS_coregistr, BASE_SS_mask,
        otheranat, ID, Session, check_visualy_final_mask, overwrite, BASE_bet, s_bind,afni_sif,fsl_sif,fs_sif, itk_sif)

        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes'
        spco(command, shell=True)

        output_for_mask = opj(masks_dir, ID + masking_img + '_mask_2.nii.gz')

    ####################################
    ########### choose Ref_file ########
    ####################################
    if useT1T2_for_coregis == True:
        Ref_file = opj(volumes_dir,ID + type_norm + '_' + otheranat + '_brain.nii.gz')
        ###creat brain of the subject template
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + opj(volumes_dir, ID + '_' + type_norm + '_' + otheranat + '_template.nii.gz') + \
        ' -b ' + output_for_mask + \
        ' -prefix ' + Ref_file + ' -expr "a*b"'
        spco([command], shell=True)
    else:
    # Organization of the folders
        Ref_file = opj(volumes_dir,ID + type_norm + '_brain.nii.gz')


    #################################### signal correction 
    for Timage in listTimage:
        command = 'singularity run' + s_bind + afni_sif + '3dresample -master ' + opj(volumes_dir, ID + '_' + Timage + '_template.nii.gz') + \
            ' -prefix ' + opj(masks_dir, ID + '_mask_rsp_' + Timage + '.nii.gz') + \
            ' -input ' + output_for_mask + ' -overwrite'
        spco([command], shell=True)

        ### masking !!!!
        ###creat brain of the subject template
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + opj(volumes_dir, ID + '_' + Timage + '_template.nii.gz') + \
                ' -b ' + opj(masks_dir, ID + '_mask_rsp_' + Timage + '.nii.gz') + \
                ' -expr "a*b" -prefix ' + opj(volumes_dir,ID + Timage + '_brain.nii.gz')
        spco([command], shell=True)

    for Timage in listTimage:
        ####plot the Reffile
        try:
            display = plotting.plot_anat(opj(dir_prepro, ID + '_acpc_cropped' + Timage + '.nii.gz'),
                                         display_mode='mosaic', dim=4)
            display.add_contours(Ref_file,
                                 linewidths=.2, colors=['red'])
            display.savefig(bids_dir + '/QC/' + ID + '_' + str(Session) + '_' + Timage + 'Ref_file.png')
            # Don't forget to close the display
            display.close()
        except:
            display = plotting.plot_anat(Ref_file,
                                         display_mode='mosaic', dim=4)
            display.savefig(bids_dir + '/QC/' + ID + '_' + str(Session) + '_' + Timage + 'Ref_file.png')
            # Don't forget to close the display
            display.close()

    ###N4BiasFieldCorrection
    try:
        BRAIN = ants.image_read(Ref_file)
        MSK   = ants.image_read(output_for_mask)
        N4    = ants.n4_bias_field_correction(BRAIN, mask=MSK,
                                           shrink_factor=4,
                                           convergence={'iters': [50, 50, 50, 50], 'tol': 1e-07},
                                           spline_param=200)
        ants.image_write(N4, opj(dir_prepro,ID + '_' + type_norm + '_brain_N4.nii.gz'), ri=False)

        try:
            # Load the NIfTI image
            img = nib.load(file_path)
            # Get the image data
            data = img.get_fdata()
            # Check if all the data is zero
            return np.all(data == 0)
        except:
            print("WARNING: N4BiasFieldCorrection failed, we can continue it is not the end of the world =)")
            command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + Ref_file + \
                      ' -expr "a" -prefix ' + opj(dir_prepro, ID + '_' + type_norm + '_brain_N4.nii.gz')
            spco([command], shell=True)

        if np.all(data == 0)== True:
            print("WARNING: N4BiasFieldCorrection failed, we can continue it is not the end of the world =)")
            command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + Ref_file + \
                      ' -expr "a" -prefix ' + opj(dir_prepro, ID + '_' + type_norm + '_brain_N4.nii.gz')
            spco([command], shell=True)

    except:
        print("WARNING: N4BiasFieldCorrection failed, we can continue it is not the end of the world =)")
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + Ref_file + \
                ' -expr "a" -prefix ' + opj(dir_prepro,ID + '_' + type_norm + '_brain_N4.nii.gz')
        spco([command], shell=True)