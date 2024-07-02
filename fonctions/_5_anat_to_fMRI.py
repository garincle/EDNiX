import os
import subprocess
import glob
import shutil
import sys
import nilearn
from nilearn.image import resample_to_img

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

from fonctions.extract_filename import extract_filename

def Refimg_to_meanfMRI(SED, anat_func_same_space, BASE_SS_coregistr, TfMRI, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
    dir_fMRI_Refth_RS_prepro3, RS, nb_run, REF_int, ID, dir_prepro, n_for_ANTS, brainmask, V_mask, W_mask, G_mask, dilate_mask,
                       list_atlases, labels_dir, costAllin, anat_subject, IhaveanANAT, doMaskingfMRI, do_anat_to_func, Method_mask_func, lower_cutoff, upper_cutoff, overwrite):

    command = '3dinfo -di ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    delta_x = str(abs(round(float(spgo([command])[-8:]), 10)))
    command = '3dinfo -dj ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    delta_y= str(abs(round(float(spgo([command])[-8:]), 10)))
    command = '3dinfo -dk ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    delta_z = str(abs(round(float(spgo([command])[-8:]), 10)))

    if ope(opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz')):
        command = '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz') + \
                  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.nii.gz') + ' -expr "a"'
        spco([command], shell=True)

    #### take opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    command = '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat_Allineate_in_func.nii.gz') + ' -b ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
    ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz') + ' -expr "a*b"'
    spco([command], shell=True)

    ####################################################################################
    ########################## use template and transfo to anat (average indiv anat)  ##
    ####################################################################################
    #average indiv anat =>coreg_ref
    # mean of all func SS => opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS.nii.gz')

    mvt_shft_INV_ANTs = ' -t ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_1InverseWarp.nii.gz') + \
                        ' -t ' + str([opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_0GenericAffine.mat'), 1])

    mvt_shft_ANTs = ' -t ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_1Warp.nii.gz') + \
                    ' -t ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_0GenericAffine.mat')

    if anat_func_same_space == True and do_anat_to_func==False:

        command = '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz') + \
                  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, 'FakeWarp.nii.gz') + ' -expr "a*0"'
        spco([command], shell=True)

        command = '3dTcat ' + overwrite +  ' ' + opj(dir_fMRI_Refth_RS_prepro1, 'FakeWarp.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1, 'FakeWarp.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1, 'FakeWarp.nii.gz') + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_1InverseWarp.nii.gz')
        spco([command], shell=True)
        command = '3dTcat ' + overwrite +  ' ' + opj(dir_fMRI_Refth_RS_prepro1, 'FakeWarp.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1, 'FakeWarp.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1, 'FakeWarp.nii.gz') + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_1Warp.nii.gz')
        spco([command], shell=True)

        import ants
        tx = ants.new_ants_transform(dimension=3)
        tx.set_parameters([1., 0., 0., 0., 1., 0., 0., 0., 1., 0., 0., 0])
        ants.write_transform(tx, opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_0GenericAffine.mat'))

    else:

        if 'i' in SED:
            restrict = '1x0.1x0.1'
        elif 'j' in SED:
            restrict = '0.1x1x0.1'
        elif 'k' in SED:
            restrict = '0.1x0.1x1'

        command = 'antsRegistration -d 3 --float 0 --verbose 1 -w [0.01,0.99] -n HammingWindowedSinc' + \
                  ' -o [' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_') + ',' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped.nii.gz') + ']' + \
                  ' --use-histogram-matching 1' + \
                  ' -t Translation[0.1]' + \
                  ' -r [' + opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.nii.gz') + ',' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz') + ',1]' + \
                  ' -m MI[' + opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.nii.gz') + ',' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz') + ',1,32,Regular,0.25]' + \
                  ' -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
                  ' -t Rigid[0.1]' + \
                  ' -m MI[' + opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.nii.gz') + ',' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz') + ',1,32,Regular,0.25]' + \
                  ' -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
                  ' -t Affine[0.1]' + \
                  ' -m MI[' + opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.nii.gz') + ',' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz') + ',1,32,Regular,0.25]' + \
                  ' -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
                  ' -t Syn[0.1,3,0]' + \
                  ' -m MI[' + opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.nii.gz') + ',' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz') + ',1,32,Regular,0.20]' + \
                  ' -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
                  ' --restrict-deformation ' + restrict
        spco(command, shell=True)

    ########## ########## help for antsRegistration ########## ########## ########## ##########
    # -f -s -c X1 x X2 x X3 x X4 from large mvt to small si 4 prend toujours 4
    ## erreur accepté a la fin de la coregistration
    # -w [0.01,0.99] threshold
    # function pour -u pour normaliser les images entre elles
    #  -n Linear ==> transfo comme NN
    ########## ########## ########## ########## ########## ########## ########## ##########

    #don't work for two different 1d matrices... so lets do it speparatly....
    for input1, output2 in zip([opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz'), 
        opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz')], 
        [opj(dir_fMRI_Refth_RS_prepro1,'mask_ref.nii.gz'), opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz'), opj(dir_fMRI_Refth_RS_prepro1,'Vmask.nii.gz'), 
        opj(dir_fMRI_Refth_RS_prepro1,'Wmask.nii.gz'), opj(dir_fMRI_Refth_RS_prepro1,'Gmask.nii.gz')]):

        # mask
        command = 'antsApplyTransforms -d 3 -e 3 -i ' + input1 + ' -n ' + n_for_ANTS + \
        ' -r ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
        ' -o ' + output2 + \
        ' --interpolation nearestNeighbor'  + \
        mvt_shft_INV_ANTs
        spco([command], shell=True)

        command = '3dmask_tool' + overwrite + ' -prefix ' + output2 + \
        ' -input ' + output2 + ' -fill_holes'
        spco(command, shell=True)

        command = '3dclust -NN1 10 -prefix ' + output2 + output2
        spco(command, shell=True)

    ## in func space resample to func 
    command = 'antsApplyTransforms -d 3 -e 3 -i ' + opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz') + ' -n ' + n_for_ANTS + \
    ' -r ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
    ' -o ' + opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI.nii.gz') + \
    mvt_shft_INV_ANTs
    spco([command], shell=True)

    #### apply to all atlases
    for atlas in list_atlases:
        ## in anat space resample to func 

        if anat_func_same_space == True:

            mvt_shft = opj(dir_prepro, ID + '_brain_for_Align_Center_inv.1D')
            command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + ' ' + opj(labels_dir, TfMRI + opb(atlas)) + ' -overwrite'
            spco([command], shell=True)

            command = '3dAllineate' + overwrite + ' -final NN -overwrite -1Dmatrix_apply ' + mvt_shft + \
            ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + \
            ' -master ' + opj(dir_prepro, ID + '_mprage_reorient' + TfMRI + '.nii.gz') + \
            ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas))
            spco([command], shell=True)

            #command = '3drefit -atrcopy ' + opj(labels_dir, TfMRI + opb(atlas)) + ' IJK_TO_DICOM_REAL ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas))
            #spco([command], shell=True)
            caca = resample_to_img(opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)), opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz'),  interpolation='nearest')
            caca.to_filename(opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)))

        elif IhaveanANAT == False:
            command = '3dresample' + overwrite + \
            ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + \
            ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
            ' -input  ' + atlas
            spco([command], shell=True)
        else:
            command = '3dresample' + overwrite + \
            ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + \
            ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
            ' -input  ' + opj(labels_dir, TfMRI + opb(atlas))
            spco([command], shell=True)

        ## in func space resample to func 
        command = 'antsApplyTransforms -d 3 -e 3 -i ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + ' -n ' + n_for_ANTS + \
        ' -r ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
        ' -o ' + opj(dir_fMRI_Refth_RS_prepro1, opb(atlas)) + \
        ' --interpolation nearestNeighbor' + \
        mvt_shft_INV_ANTs
        spco([command], shell=True)


    command = '3dinfo -di ' + anat_subject
    delta_x1 = str(abs(round(float(spgo([command])[-8:]), 10)))
    command = '3dinfo -dj ' + anat_subject
    delta_y1= str(abs(round(float(spgo([command])[-8:]), 10)))
    command = '3dinfo -dk ' + anat_subject
    delta_z1 = str(abs(round(float(spgo([command])[-8:]), 10)))

    command = '3dresample' + overwrite + \
    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS_anat_resolution.nii.gz') + \
    ' -dxyz ' + delta_x1 + ' ' + delta_y1 + ' ' + delta_z1 + ' ' + \
    ' -rmode Cu -input  ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    spco([command], shell=True)


    #### creat a nice anat in func space
    if anat_func_same_space == True:
        anatstd = opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_plot.nii.gz')
        ##### apply the recenter fmri
        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + anatstd + ' ' + anat_subject + ' -overwrite'
        spco([command], shell=True)

        command = '3dAllineate' + overwrite + ' -overwrite -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + anatstd + \
        ' -input  ' + anatstd + \
        ' -master ' + opj(dir_prepro, ID + '_mprage_reorient' + TfMRI + '.nii.gz')

        spco([command], shell=True)

    else:
        anatstd = anat_subject
        command = '3dcalc' + overwrite + ' -a ' + anatstd + \
        ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_plot.nii.gz') + ' -expr "a"'
        spco([command], shell=True)


    ## in func space resample to func 
    command = 'antsApplyTransforms -d 3 -e 3 -i ' + anatstd + ' -n ' + n_for_ANTS + \
    ' -r ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS_anat_resolution.nii.gz') + \
    ' -o ' + opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI_anat_resolution.nii.gz') + \
    mvt_shft_INV_ANTs
    spco([command], shell=True)


    ###finaly mask the func with mask
    for i in range(0, int(nb_run)):
        root_RS = extract_filename(RS[i])
        ### take opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
        
        command = '3dcalc' + overwrite + ' -a ' +   opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz') + ' -b ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref.nii.gz') + \
        ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref_RcT_masked.nii.gz') + ' -expr "a*b"'
        spco([command], shell=True)
