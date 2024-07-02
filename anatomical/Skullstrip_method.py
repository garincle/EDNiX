###################################
###      Skullstrip method      ###
###################################
import os
import subprocess
import glob
import shutil
import sys


#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput


def Skullstrip_method(step_skullstrip, brain_skullstrip, masking_img, brain_skullstrip_1, brain_skullstrip_2, masks_dir, volumes_dir, dir_prepro, type_norm, n_for_ANTS, dir_transfo, BASE_SS_coregistr, BASE_SS_mask, 
    otheranat, ID, Session, check_visualy_final_mask, overwrite):

    if step_skullstrip == 1:
        masking_img = masking_img
        input_for_msk = opj(dir_prepro, ID + '_anat_reorient_NU' + masking_img + '.nii.gz')
        output_for_mask = opj(masks_dir, ID + masking_img + '_mask_1.nii.gz')
        brain_skullstrip = brain_skullstrip_1

    elif step_skullstrip == 2:
        masking_img = type_norm
        input_for_msk = opj(volumes_dir, ID + '_' + masking_img + '_template.nii.gz')
        output_for_mask = opj(masks_dir, ID + masking_img + '_mask_2.nii.gz')
        brain_skullstrip = brain_skullstrip_2


    else:
        print('no step_skullstrip ??')

    if brain_skullstrip =='bet2_ANTS':
        #####creat an approximate brain mask
        command = 'bet2 ' + input_for_msk + ' ' + opj(masks_dir, ID + '_bet' + masking_img) + \
        ' -f 0.01 -m'
        spco([command], shell=True)

        #####creat a precise brain mask
        command = 'antsRegistration -d 3 --float 0 --verbose 1 -u 1 -w [0.05,0.95]  -n ' + n_for_ANTS + \
        ' -o [' + opj(dir_transfo,'template_to_' + masking_img + '_SyN_') + ',' + opj(dir_prepro,'template_to_' + masking_img + '_SyN.nii.gz') + ']'+ \
        ' -t Affine[0.1] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
        ' -m MI[' + opj(masks_dir, ID + '_bet' + masking_img + '.nii.gz') + ',' + BASE_SS_coregistr + ',1,32,Regular,0.2]' + \
        ' -t Syn[0.1,3,0] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
        ' -m CC[' + opj(masks_dir, ID + '_bet' + masking_img + '.nii.gz') + ',' + BASE_SS_coregistr + ',1,4,Regular,0.2]'
        spco([command], shell=True)

        command = 'antsApplyTransforms -d 3 -i ' + BASE_SS_mask + \
            ' -r ' + opj(masks_dir, ID + '_bet' + masking_img + '.nii.gz') + \
            ' -o ' + opj(dir_prepro, masking_img + 'template_brainmask.nii.gz') + \
            ' -t ' + opj(dir_transfo,'template_to_' + masking_img + '_SyN_1Warp.nii.gz') + \
            ' -t ' + opj(dir_transfo,'template_to_' + masking_img + '_SyN_0GenericAffine.mat') + \
            ' -n gaussian'
        spco([command], shell=True)

        Ex_Mask    = opj(dir_prepro,'mask_tmp' + masking_img + '.nii.gz')

        command = 'ThresholdImage 3 ' + opj(dir_prepro, masking_img + 'template_brainmask.nii.gz') + \
        ' ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + ' 0.5 1 1 0'
        spco([command], shell=True)

        command = 'ImageMath 3 ' + Ex_Mask + ' MD ' + opj(dir_prepro, masking_img + 'template_brainmask.nii.gz') + ' 2'
        spco([command], shell=True)

        command = 'ImageMath 3 ' + Ex_Mask + ' GetLargestComponent ' + Ex_Mask
        spco([command], shell=True)

        shutil.copyfile(Ex_Mask,output_for_mask)

        if check_visualy_each_img == 'YES':
            command = 'itksnap -g ' + Ref_file + \
            ' -s ' + opj(masks_dir, ID + masking_img + '_mask.nii.gz')
            spco([command], shell=True)

    elif brain_skullstrip =='Custum_ANTS_NL':

        #####creat a precise brain mask

        command = 'antsRegistration -d 3 --float 0 --verbose 1 -u 1 -w [0.05,0.95]  -n ' + n_for_ANTS + \
        ' -o [' + opj(dir_transfo,'template_to_' + masking_img + '_SyN_') + ',' + opj(dir_prepro,'template_to_' + masking_img + '_SyN.nii.gz') + ']'+ \
        ' -t Affine[0.1] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
        ' -m MI[' + input_for_msk + ',' + BASE_SS_coregistr + ',1,32,Regular,0.2]' + \
        ' -t Syn[0.1,3,0] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
        ' -m CC[' + input_for_msk + ',' + BASE_SS_coregistr + ',1,4,Regular,0.2]'
        spco([command], shell=True)

        command = 'antsApplyTransforms -d 3 -i ' + BASE_SS_mask + \
            ' -r ' + input_for_msk + \
            ' -o ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + \
            ' -t ' + opj(dir_transfo,'template_to_' + masking_img + '_SyN_1Warp.nii.gz') + \
            ' -t ' + opj(dir_transfo,'template_to_' + masking_img + '_SyN_0GenericAffine.mat') + \
            ' -n gaussian'
        spco([command], shell=True)

        Ex_Mask    = opj(dir_prepro,'mask_tmp' + masking_img + '.nii.gz')

        command = 'ThresholdImage 3 ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + \
        ' ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + ' 0.5 1 1 0'
        spco([command], shell=True)
        command = 'ImageMath 3 ' + Ex_Mask + ' MD ' + opj(dir_prepro, masking_img + 'template_brainmask.nii.gz') + ' 2'
        spco([command], shell=True)
        command = 'ImageMath 3 ' + Ex_Mask + ' GetLargestComponent ' + Ex_Mask
        spco([command], shell=True)

        shutil.copyfile(Ex_Mask, output_for_mask)

        if check_visualy_final_mask == True:
            command = 'itksnap -g ' + input_for_msk + \
            ' -s ' + output_for_mask
            spco([command], shell=True)

    elif brain_skullstrip =='Custum_ANTS':

        ####remplace by applytransfo
        command = 'antsRegistration -d 3 --float 0 --verbose 1 -u 1 -w [0.05,0.95] -n Linear' + \
        ' -o [' + opj(dir_transfo,'template_to_' + masking_img + '_SyN_') + ',' + opj(dir_prepro,'template_to_' + masking_img + '_SyN.nii.gz') + ']'+ \
        ' -t Affine[0.1] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
        ' -m MI[' + input_for_msk + ',' + BASE_SS_coregistr + ',1,32,Regular,0.2]'
        spco([command], shell=True)

        command = 'antsApplyTransforms -d 3 -i ' + BASE_SS_mask + \
            ' -r ' + input_for_msk + \
            ' -o ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + \
            ' -t ' + opj(dir_transfo,'template_to_' + masking_img + '_SyN_0GenericAffine.mat') + \
            ' -n gaussian'
        spco([command], shell=True)

        Ex_Mask    = opj(dir_prepro,'mask_tmp' + masking_img + '.nii.gz')

        command = 'ThresholdImage 3 ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + \
        ' ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + ' 0.5 1 1 0'
        spco([command], shell=True)
        command = 'ImageMath 3 ' + Ex_Mask + ' MD ' + opj(dir_prepro, masking_img + 'template_brainmask.nii.gz') + ' 2'
        spco([command], shell=True)
        command = 'ImageMath 3 ' + Ex_Mask + ' GetLargestComponent ' + Ex_Mask
        spco([command], shell=True)

        shutil.copyfile(Ex_Mask, output_for_mask)

        if check_visualy_final_mask == True:
            command = 'itksnap -g ' + input_for_msk + \
            ' -s ' + output_for_mask
            spco([command], shell=True)

    elif brain_skullstrip == 'MachinL':
        spco(['python3', '/home/cgarin/Documents/0_Clement/CODE/CODE/NHP-BrainExtraction-master/UNet_Model/muSkullStrip.py',
        '-in', input_for_msk, '-model', '/home/cgarin/Documents/1_Macaque_MRI/4_SSwarper_muSkull/result/model-20-epoch', '-out', 
        masks_dir])
        shutil.copyfile(opj(masks_dir, opb(ID + '_' + masking_img + '_template_pre_mask.nii.gz')), output_for_mask)

    elif brain_skullstrip == 'antsBrainExtraction':
        spco(['antsBrainExtraction.sh', '-d', '3', 
            '-a', input_for_msk, 
            '-e', BASE_SS_coregistr, '-m', BASE_SS_mask, '-o', volumes_dir + '/'])
        shutil.copyfile(opj(volumes_dir, 'BrainExtractionMask.nii.gz'), output_for_mask) 

    elif brain_skullstrip == '3dSkullStrip':

        spco(['3dSkullStrip', '-prefix', output_for_mask, '-overwrite',
            '-input', input_for_msk, '-blur_fwhm', '2', '-orig_vol', '-mask_vol', '-use_skull', '-monkey'])

        command = '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 2'
        spco(command, shell=True)

    elif brain_skullstrip == 'Custum_1':
        import nibabel as nib
        import nilearn
        loadimg = nib.load(opj(volumes_dir, ID + '_' + otheranat + '_template.nii.gz')).get_data() 
        loadimgsort85 =  np.percentile(np.abs(loadimg)[np.abs(loadimg)>0], 40)
        mask_imag = nilearn.image.threshold_img(opj(volumes_dir, ID + '_' + otheranat + '_template.nii.gz'), loadimgsort85, cluster_threshold=10)
        mask_imag.to_filename(output_for_mask)
        command = '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 1'
        spco(command, shell=True)
        spco(['3dresample', '-master', input_for_msk, '-prefix', output_for_mask, '-input', output_for_mask, '-overwrite', '-bound_type', 'SLAB'])

    elif brain_skullstrip == 'Custum_Baboon':
        from nilearn.masking import compute_epi_mask
        #convert to float
        command = 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
        spco([command], shell=True)
        mask_img = compute_epi_mask(input_for_msk[:-7] + '_float.nii.gz', lower_cutoff=0.85, upper_cutoff=0.95, connected=True, opening=3,
            exclude_zeros=False, ensure_finite=True)
        mask_img.to_filename(output_for_mask)
        command = '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 15'
        spco(command, shell=True)

    elif brain_skullstrip == 'Custum_Macaque':
        from nilearn.masking import compute_epi_mask
        #convert to float
        command = 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
        spco([command], shell=True)
        mask_img = compute_epi_mask(input_for_msk[:-7] + '_float.nii.gz', lower_cutoff=0.75, upper_cutoff=0.90, connected=True, opening=3,
            exclude_zeros=False, ensure_finite=True)
        mask_img.to_filename(output_for_mask)
        command = '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 10'
        spco(command, shell=True)


    elif brain_skullstrip == 'Custum_mouse':
        from nilearn.masking import compute_epi_mask
        #convert to float
        command = 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
        spco([command], shell=True)
        mask_img = compute_epi_mask(input_for_msk[:-7] + '_float.nii.gz', lower_cutoff=0.2, upper_cutoff=0.80, connected=True, opening=3,
            exclude_zeros=False, ensure_finite=True)
        mask_img.to_filename(output_for_mask)
        command = '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 3 -1'
        spco(command, shell=True)

    elif brain_skullstrip == 'Custum_dog':
        from nilearn.masking import compute_epi_mask
        #convert to float
        command = 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
        spco([command], shell=True)
        mask_img = compute_epi_mask(input_for_msk[:-7] + '_float.nii.gz', lower_cutoff=0.85, upper_cutoff=0.95, connected=True, opening=3,
            exclude_zeros=False, ensure_finite=True)
        mask_img.to_filename(output_for_mask)
        command = '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 23'
        spco(command, shell=True)

    elif brain_skullstrip =='bet2':
        #####creat an approximate brain mask
        command = 'bet2 ' + input_for_msk + ' ' + opj(masks_dir, ID + '_bet' + masking_img + '.nii.gz') + \
        ' -f 0.70'
        spco([command], shell=True)
        spco(['3dcalc', '-a', opj(masks_dir, ID + '_bet' + masking_img + '.nii.gz'), '-expr', 'step(a)', '-prefix', output_for_mask, '-overwrite'])

    elif brain_skullstrip == 'Custum T1/T2':
        spco(['3dresample', '-master', input_for_msk, '-prefix', opj(volumes_dir, ID + '_' + otheranat + '_template_RSPL.nii.gz'), '-input', opj(volumes_dir, ID + '_' + otheranat + '_template.nii.gz'), '-overwrite'])
        spco(['3dcalc', '-a', input_for_msk, '-b', opj(volumes_dir, ID + '_' + otheranat + '_template_RSPL.nii.gz'), '-expr', 'b*(step(a))', '-prefix', opj(volumes_dir, ID + '_' + otheranat + '_template_RSPL.nii.gz'), '-overwrite'])
        from nilearn.image import math_img
        log_img = math_img("np.where((img1 > 1) & (img2 > 1), img2, 0)", img1=input_for_msk, img2=opj(volumes_dir, ID + '_' + otheranat + '_template_RSPL.nii.gz'))
        log_img.to_filename(opj(volumes_dir, ID + '_' + otheranat + '_template_RSPL.nii.gz'))
        log_img = math_img("np.where((img1 > 1) & (img2 > 1), img1, 0)", img1=input_for_msk, img2=opj(volumes_dir, ID + '_' + otheranat + '_template_RSPL.nii.gz'))
        log_img.to_filename(opj(volumes_dir, ID + '_' + masking_img + '_template_RSPL.nii.gz'))
        command = '3dcalc -overwrite -a ' + opj(volumes_dir, ID + '_' + masking_img + '_template_RSPL.nii.gz') + \
        ' -b ' + opj(volumes_dir, ID + '_' + otheranat + '_template_RSPL.nii.gz') + \
        ' -expr "a/b" -datum float  -prefix ' + opj(volumes_dir, ID + '_' + masking_img + '_' + otheranat + '_template.nii.gz')
        spco([command], shell=True)
        spco(['antsBrainExtraction.sh', '-d', '3', 
        '-a',opj(volumes_dir, ID + '_' + masking_img + '_' + otheranat + '_template.nii.gz'), 
        '-e', BASE_SS_coregistr, '-m', BASE_SS_mask, 
        '-o', volumes_dir + '/'])
        mask_img = compute_epi_mask(opj(volumes_dir, ID + '_' + masking_img + '_' + otheranat + '_template.nii.gz'), lower_cutoff=0.10, upper_cutoff=0.60, connected=True, opening=3)
        mask_img.to_filename(output_for_mask)

    elif brain_skullstrip =='QWARP':

        command = '3dQwarp -overwrite -iwarp' + \
        ' -base ' + BASE_SS_coregistr + \
        ' -prefix ' + opj(dir_prepro,'template_to_' + masking_img + '_AFNIQ.nii.gz') + \
        ' -source ' + input_for_msk + ' -maxlev 5 -resample'
        spco(command, shell=True)

        spco(['3dNwarpApply', '-nwarp', opj(dir_prepro,'template_to_' + masking_img + '_AFNIQ_WARPINV.nii.gz'), 
        '-source', BASE_SS_mask, '-master', input_for_msk, '-interp', 'NN',
        '-prefix', opj(dir_prepro,masking_img + 'template_brainmask.nii.gz'), '-overwrite'])

        Ex_Mask    = opj(dir_prepro,'mask_tmp' + masking_img + '.nii.gz')

        command = '3dmask_tool -overwrite -prefix ' + Ex_Mask + \
        ' -input ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + ' -fill_holes' # -dilate_input 2'
        spco(command, shell=True)

        shutil.copyfile(Ex_Mask,output_for_mask)

        if check_visualy_final_mask == True:
            command = 'itksnap -g ' + input_for_msk + \
            ' -s ' + output_for_mask
            spco([command], shell=True)

    elif brain_skullstrip == '3dSkullStrip_Rat':
        if ID in ['301502','302101','302105','302106','301603', '300908','301500','301501','301503','301504','301505','301508','301509','302107','302108','300600']:
            spco(['3dSkullStrip', '-prefix', output_for_mask, '-overwrite',
                '-input', input_for_msk, '-orig_vol', '-mask_vol', '-rat'])

            command = '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes -dilate_input 15'
            spco(command, shell=True)
        else:
            spco(['3dSkullStrip', '-prefix', output_for_mask, '-overwrite',
                '-input', input_for_msk, '-orig_vol', '-mask_vol', '-surface_coil', '-rat'])

            command = '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes -dilate_input 4'
            spco(command, shell=True)

    elif brain_skullstrip == 'NoSkullStrip':
        spco(['3dcalc', '-a', input_for_msk, '-expr', 'step(a)', '-prefix', output_for_mask, '-overwrite'])

    elif brain_skullstrip =='sammba_rat':
        from sammba import interfaces
        from sammba.interfaces import HistogramMask

        command = 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
        spco([command], shell=True)

        nichols_masker = HistogramMask()
        nichols_masker.inputs.in_file = input_for_msk[:-7] + '_float.nii.gz'
        nichols_masker.inputs.volume_threshold = 2500
        #nichols_masker.inputs.upper_cutoff = 0.2
        #nichols_masker.inputs.lower_cutoff = 0.8
        #nichols_masker.inputs.intensity_threshold = 500
        #nichols_masker.inputs.opening = 2
        #nichols_masker.inputs.closing = 10
        nichols_masker.inputs.dilation_size = (1, 2, 3)
        nichols_masker.inputs.connected = True
        nichols_masker.inputs.out_file = output_for_mask
        res = nichols_masker.run()  # doctest: +SKIP

        if check_visualy_final_mask == True:
            command = 'itksnap -g ' + input_for_msk + \
                      ' -s ' + output_for_mask
            spco([command], shell=True)

    elif brain_skullstrip =='custum_rat':
            from nilearn.masking import compute_epi_mask
            #convert to float
            command = 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
            spco([command], shell=True)
            mask_img = compute_epi_mask(input_for_msk[:-7] + '_float.nii.gz', lower_cutoff=0.8, upper_cutoff=0.85, connected=True, opening=3,
                                        exclude_zeros=True, ensure_finite=True)
            mask_img.to_filename(output_for_mask)
            command = '3dmask_tool -overwrite -prefix ' + output_for_mask + \
                      ' -input ' + output_for_mask + ' -fill_holes -dilate_input 3'
            spco(command, shell=True)

        ################################################
        ###### in use for study don't touch ############
        ################################################

    elif brain_skullstrip == 'Custum_Macaque2':
        if ID == 'Pickle' and Session == 5:
            spco(['python3', '/home/cgarin/Documents/0_Clement/CODE/CODE/NHP-BrainExtraction-master/UNet_Model/muSkullStrip.py',
            '-in', input_for_msk, '-model', '/home/cgarin/Documents/1_Macaque_MRI/4_SSwarper_muSkull/result/model-20-epoch', '-out', 
            masks_dir])
            shutil.copyfile(opj(masks_dir, opb(ID + '_anat_reorient_NU' + masking_img + '_pre_mask.nii.gz')), output_for_mask)
            command = '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes -dilate_input 10' #10 for Pickel session 5 and Trinity Session 6
            spco(command, shell=True)

        elif ID == 'Trinity' and Session == 6:
            spco(['python3', '/home/cgarin/Documents/0_Clement/CODE/CODE/NHP-BrainExtraction-master/UNet_Model/muSkullStrip.py',
            '-in', input_for_msk, '-model', '/home/cgarin/Documents/1_Macaque_MRI/4_SSwarper_muSkull/result/model-20-epoch', '-out', 
            masks_dir])
            shutil.copyfile(opj(masks_dir, opb(ID + '_anat_reorient_NU' + masking_img + '_pre_mask.nii.gz')), output_for_mask)
            command = '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes -dilate_input 10' #10 for Pickel session 5 and Trinity Session 6
            spco(command, shell=True)

        else:
            spco(['python3', '/home/cgarin/Documents/0_Clement/CODE/CODE/NHP-BrainExtraction-master/UNet_Model/muSkullStrip.py',
            '-in', input_for_msk, '-model', '/home/cgarin/Documents/1_Macaque_MRI/4_SSwarper_muSkull/result/model-20-epoch', '-out', 
            masks_dir])
            shutil.copyfile(opj(masks_dir, opb(ID + '_anat_reorient_NU' + masking_img + '_pre_mask.nii.gz')), output_for_mask)
            command = '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes -dilate_input 3' #10 for Pickel session 5 and Trinity Session 6
            spco(command, shell=True)


    elif brain_skullstrip =='Custum_QWARP':

        if ID in ['Oliver'] and Session in [3] or ID in ['Roshan'] and Session in [3] or ID in ['Quantum'] and Session in [3] or ID in ['Roshan'] and Session in [6]:

            command = '3dQwarp -overwrite -iwarp' + \
            ' -base ' + BASE_SS_coregistr + \
            ' -prefix ' + opj(dir_prepro,'template_to_' + masking_img + '_AFNIQ.nii.gz') + \
            ' -source ' + input_for_msk + ' -maxlev 5 -lpa -resample'
            spco(command, shell=True)

            spco(['3dNwarpApply', '-nwarp', opj(dir_prepro,'template_to_' + masking_img + '_AFNIQ_WARPINV.nii.gz'), 
            '-source', BASE_SS_mask, '-master', input_for_msk, '-interp', 'NN',
            '-prefix', opj(dir_prepro,masking_img + 'template_brainmask.nii.gz'), '-overwrite'])

            Ex_Mask    = opj(dir_prepro,'mask_tmp' + masking_img + '.nii.gz')

            command = '3dmask_tool -overwrite -prefix ' + Ex_Mask + \
            ' -input ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + ' -fill_holes' # -dilate_input 2'
            spco(command, shell=True)

            shutil.copyfile(Ex_Mask,output_for_mask)

            if check_visualy_final_mask == True:
                command = 'itksnap -g ' + input_for_msk + \
                ' -s ' + output_for_mask
                spco([command], shell=True)

        else:

            command = '3dQwarp -overwrite -iwarp' + \
            ' -base ' + BASE_SS_coregistr + \
            ' -prefix ' + opj(dir_prepro,'template_to_' + masking_img + '_AFNIQ.nii.gz') + \
            ' -source ' + input_for_msk + ' -maxlev 5 -resample'
            spco(command, shell=True)

            spco(['3dNwarpApply', '-nwarp', opj(dir_prepro,'template_to_' + masking_img + '_AFNIQ_WARPINV.nii.gz'), 
            '-source', BASE_SS_mask, '-master', input_for_msk, '-interp', 'NN',
            '-prefix', opj(dir_prepro,masking_img + 'template_brainmask.nii.gz'), '-overwrite'])

            Ex_Mask    = opj(dir_prepro,'mask_tmp' + masking_img + '.nii.gz')

            command = '3dmask_tool -overwrite -prefix ' + Ex_Mask + \
            ' -input ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + ' -fill_holes' # -dilate_input 2'
            spco(command, shell=True)

            shutil.copyfile(Ex_Mask,output_for_mask)

            if check_visualy_final_mask == True:
                command = 'itksnap -g ' + input_for_msk + \
                ' -s ' + output_for_mask
                spco([command], shell=True)


    elif brain_skullstrip =='Custum_QWARPT2':

        command = '3dQwarp -overwrite -lpa -iwarp' + \
        ' -base ' + BASE_SS_coregistr + \
        ' -prefix ' + opj(dir_prepro,'template_to_' + masking_img + '_AFNIQ.nii.gz') + \
        ' -source ' + input_for_msk + ' -maxlev 3 -resample'
        spco(command, shell=True)

        spco(['3dNwarpApply', '-nwarp', opj(dir_prepro,'template_to_' + masking_img + '_AFNIQ_WARPINV.nii.gz'), 
        '-source', BASE_SS_mask, '-master', input_for_msk, '-interp', 'NN',
        '-prefix', opj(dir_prepro,masking_img + 'template_brainmask.nii.gz'), '-overwrite'])

        Ex_Mask    = opj(dir_prepro,'mask_tmp' + masking_img + '.nii.gz')

        command = '3dmask_tool -overwrite -prefix ' + Ex_Mask + \
        ' -input ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + ' -fill_holes'
        spco(command, shell=True)

        shutil.copyfile(Ex_Mask,output_for_mask)

        if check_visualy_final_mask == True:
            command = 'itksnap -g ' + input_for_msk + \
            ' -s ' + output_for_mask
            spco([command], shell=True)


    else:
        print("ERROR in brain_skullstrip name???")

    return(output_for_mask)