#import
import os
import subprocess
import sys


#########################################
########### Subject loader################
##########################################

MAIN_PATH   = opj('/','srv','projects','easymribrain')
sys.path.append(opj(MAIN_PATH + 'Code','EasyMRI_brain-master'))

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
ops = os.path.splitext
spco = subprocess.check_output
spgo = subprocess.getoutput

def QuantifT1T2(all_ID, all_Session, all_data_path, IgotbothT1T2, list_atlases, otheranat, type_norm, max_sessionlist, overwrite_option,s_bind,afni_sif):

    ###########################################################################################################################################################
    ############################################################## start the proces ###########################################################################
    ###########################################################################################################################################################
    ######### define other usefull paramater automatically (do no touch)#########

    if overwrite_option == True:
        overwrite = ' -overwrite'
    else:
        overwrite = ''

    listTimage = []
    if IgotbothT1T2 == True:
        listTimage = [otheranat, type_norm]
    else:
        listTimage = [type_norm]
    ####################################################################################
    ########################## Start the pipeline !!!!!!!!!!!!!!!!!!!!!!   #############
    ####################################################################################


    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, max_sessionlist):
        print('############################################################################################################################################################################' +
        '###################################################### work on ' + str(ID) + ' Session ' + str(Session) + '###############################################################' +
        '############################################################################################################################################################################')

        # The anatomy
        path_anat    = opj(data_path,'anat/')
        dir_native    = opj(path_anat,'native')
        dir_prepro    = opj(dir_native,'01_preprocess')
        wb_native_dir = opj(dir_native, '02_Wb')
        wb_T1T2_dir = opj(dir_native, '03_T1T2')


        volumes_dir = opj(wb_native_dir, 'volumes')
        masks_dir = opj(volumes_dir, 'masks')
        labels_dir = opj(volumes_dir, 'labels')

        #creat path
        if ope(wb_T1T2_dir) == False:
            os.makedirs(wb_T1T2_dir)

        mvt_shft = opj(dir_prepro,ID + '_brain_for_Align_Center_inv.1D')
        if not ope(mvt_shft):
            command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + 'cat_matvec ' + opj(dir_prepro,ID + '_brain_for_Align_Center.1D') + \
            ' -I | tail -n +3 > ' + mvt_shft
            spco([command], shell=True)

        brainmask = opj(masks_dir, 'brain_mask_in_anat_DC.nii.gz')

        for atlas in list_atlases:
            ##### apply the recenter fmri
            command = 'singularity run' + s_bind + afni_sif + '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(wb_T1T2_dir, opb(atlas)) + ' ' + opj(labels_dir, type_norm + opb(atlas)) + ' -overwrite'
            spco([command], shell=True)

            command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -overwrite -final NN -1Dmatrix_apply ' + mvt_shft + \
            ' -prefix ' + opj(wb_T1T2_dir, opb(atlas)) + \
            ' -input  ' + opj(wb_T1T2_dir, opb(atlas))
            spco([command], shell=True)

            spco(['singularity run' + s_bind + afni_sif + '3dresample', '-master', opj(dir_prepro, ID + '_mprage_reorient' + type_norm + '.nii.gz'), '-prefix',
                  opj(wb_T1T2_dir, opb(atlas)), '-input',
                  opj(wb_T1T2_dir, opb(atlas)), '-overwrite'])
        '''
        ##### apply the recenter fmri
        command = 'singularity run' + s_bind + afni_sif + '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(wb_T1T2_dir, opb(brainmask)) + ' ' + brainmask + ' -overwrite'
        spco([command], shell=True)

        command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -overwrite -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(wb_T1T2_dir, opb(brainmask)) + \
        ' -input  ' + opj(wb_T1T2_dir, opb(brainmask))
        spco([command], shell=True)


        ### creat a mask of the T1
        spco(['singularity run' + s_bind + afni_sif + '3dresample', '-master', opj(dir_prepro, ID + '_mprage_reorient' + type_norm + '.nii.gz'), '-prefix',
              opj(wb_T1T2_dir, opb(brainmask)), '-input',
              opj(wb_T1T2_dir, opb(brainmask)), '-overwrite'])
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(dir_prepro, ID + '_mprage_reorient' + type_norm + '.nii.gz') + ' -b ' +  opj(wb_T1T2_dir, 'brain_mask_in_anat_DC.nii.gz') + \
                  ' -prefix ' + opj(wb_T1T2_dir, ID + '_' + type_norm + '_brain_RSPL.nii.gz') + ' -expr "a*b"'
        spco([command], shell=True)

        ### creat a mask of the T2
        spco(['singularity run' + s_bind + afni_sif + '3dresample', '-master', opj(dir_prepro, ID + '_mprage_reorient' + type_norm + '.nii.gz'), '-prefix',
              opj(wb_T1T2_dir, ID + '_' + otheranat + '_brain_RSPL.nii.gz'), '-input',
              opj(dir_prepro, ID + '_mprage_reorient' + otheranat + '.nii.gz'), '-overwrite'])
        spco(['singularity run' + s_bind + afni_sif + '3dcalc', '-a', opj(wb_T1T2_dir, ID + '_' + type_norm + '_brain_RSPL.nii.gz'), '-b',
              opj(wb_T1T2_dir, ID + '_' + otheranat + '_brain_RSPL.nii.gz'), '-expr', 'b*(step(a))', '-prefix',
              opj(wb_T1T2_dir, ID + '_' + otheranat + '_brain_RSPL.nii.gz'), '-overwrite'])

        from nilearn.image import math_img
        log_img = math_img("np.where((img1 > 1) & (img2 > 1), img2, 0)",
                           img1=opj(wb_T1T2_dir, ID + '_' + type_norm + '_brain_RSPL.nii.gz'),
                           img2=opj(wb_T1T2_dir, ID + '_' + otheranat + '_brain_RSPL.nii.gz'))
        log_img.to_filename(opj(wb_T1T2_dir, ID + '_' + otheranat + '_brain_RSPL.nii.gz'))

        log_img = math_img("np.where((img1 > 1) & (img2 > 1), img1, 0)",
                           img1=opj(wb_T1T2_dir, ID + '_' + type_norm + '_brain_RSPL.nii.gz'),
                           img2=opj(wb_T1T2_dir, ID + '_' + otheranat + '_brain_RSPL.nii.gz'))
        log_img.to_filename(opj(wb_T1T2_dir, ID + '_' + type_norm + '_brain_RSPL.nii.gz'))

        command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + opj(wb_T1T2_dir, ID + '_' + type_norm + '_brain_RSPL.nii.gz') + \
                  ' -b ' + opj(wb_T1T2_dir, ID + '_' + otheranat + '_brain_RSPL.nii.gz') + \
                  ' -expr "a/b" -datum float  -prefix ' + opj(wb_T1T2_dir,
                                                              ID + '_' + type_norm + '_' + otheranat + '_brain_RSPL.nii.gz')
        spco([command], shell=True)
        '''

