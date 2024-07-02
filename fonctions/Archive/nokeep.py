    mvt_shft = opj(dir_prepro,ID + '_brain_for_Align_Center_inv.1D')
    if os.path.isfile(mvt_shft) == False and os.path.isfile(opj(dir_prepro,ID + '_brain_for_Align_Center.1D')) == True:
        print('dd')
        command = 'cat_matvec ' + opj(dir_prepro,ID + '_brain_for_Align_Center.1D') + \
        ' -I > ' + mvt_shft
        spco([command], shell=True)



        #####XXX might not work... to test

        # mask
        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz') + ' ' + brainmask + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz')
        spco([command], shell=True)  

        # dilate mask
        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat_orig.nii.gz') + ' ' + brainmask_D + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat_orig.nii.gz')
        spco([command], shell=True)  

        # white matter mask
        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz') + ' ' + W_mask + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz')
        spco([command], shell=True)  

        # ventricule mask
        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz') + ' ' + V_mask + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz')
        spco([command], shell=True) 

        ####skullstrip the indiv pre-template (and clean it)
        # B0 correction
        command = 'DenoiseImage' + ' -d 3 -i ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + TfMRI + '.nii.gz') + ' -n Gaussian -r 3 ' + \
            ' -o ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + TfMRI + '_NU_denoise.nii.gz')
        spco([command], shell=True)