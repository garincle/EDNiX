###################################
###      Skullstrip method      ###
###################################
import os
import subprocess
import shutil
import nibabel as nib
import ants
import math
import nilearn
from nilearn.masking import compute_epi_mask
from nilearn.image import math_img
import anatomical.Histrogram_mask_EMB

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput


def Skullstrip_method(step_skullstrip, brain_skullstrip, masking_img, brain_skullstrip_1, brain_skullstrip_2, masks_dir, volumes_dir, dir_prepro, type_norm, n_for_ANTS, dir_transfo, BASE_SS_coregistr, BASE_SS_mask,
    otheranat, ID, Session, check_visualy_final_mask, overwrite, BASE_bet, s_bind,afni_sif,fsl_sif,fs_sif):

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
        #
        hd_IMG = ants.image_header_info(input_for_msk)
        cmd = 'singularity run' + s_bind + fsl_sif + 'bet2 ' + input_for_msk + ' ' + opj(masks_dir, ID + '_bet' + masking_img) + \
              ' -f 0.40 -c ' + str(math.ceil(int(hd_IMG['dimensions'][0]) / 2)) + ' ' + str(math.ceil(int(hd_IMG['dimensions'][1]) / 2)) + \
              ' ' + str(math.ceil(int(hd_IMG['dimensions'][2]) / 2)) + ' -m'
        spco([cmd], shell=True)

        #####creat a precise brain mask

        IMG      = ants.image_read(input_for_msk)
        REF_BET  = ants.image_read(BASE_bet)
        tmp_bet  = ants.image_read(opj(masks_dir, ID + '_bet' + masking_img + '.nii.gz'))
        REF_MASK = ants.image_read(BASE_SS_mask)

        mTx = ants.registration(fixed=tmp_bet, moving=REF_BET,
                                 type_of_transform='SyNRA',
                                outprefix=opj(dir_transfo,'template_to_' + masking_img + '_SyN_'),
                                 grad_step=0.1,
                                 flow_sigma=3,
                                 total_sigma=0,
                                 aff_sampling=32,
                                 aff_random_sampling_rate=0.2,
                                 syn_sampling=32,
                                 aff_iterations=(1000, 500, 250, 100),
                                 aff_shrink_factors=(8, 4, 2, 1),
                                 aff_smoothing_sigmas=(3, 2, 1, 0),
                                 reg_iterations=(1000, 500, 250, 100),
                                 reg_smoothing_sigmas=(3, 2, 1, 0),
                                 reg_shrink_factors=(8, 4, 2, 1),
                                 verbose=True)


        tmp_mask1 = ants.apply_transforms(fixed=tmp_bet, moving=REF_MASK,
                                          transformlist=mTx['fwdtransforms'], interpolator='nearestNeighbor')

        tmp_mask1 = ants.threshold_image(tmp_mask1, 0.5, 1, 1, 0, True)
        tmp_mask1 = ants.morphology(tmp_mask1, operation='dilate', radius=2, mtype='binary', shape='ball')
        tmp_mask1 = ants.iMath(tmp_mask1, operation='GetLargestComponent')

        seg_tmp = ants.atropos(a=IMG, m='[0.1,1x1x1]', c='[3,0]', i='kmeans[3]', x=tmp_mask1)

        #  Clean up:
        seg_tmp2 = ants.iMath(seg_tmp['segmentation'], 'Pad', 10)
        W = ants.threshold_image(seg_tmp2, 3, 3, 1, 0, True)
        W = ants.iMath(W, operation='GetLargestComponent')
        W = W * 3
        G = ants.threshold_image(seg_tmp2, 2, 2, 1, 0, True)
        G = ants.iMath(G, operation='GetLargestComponent')
        TMP1 = ants.iMath(G, 'FillHoles', 2)
        G = G * TMP1
        C = ants.threshold_image(seg_tmp2, 1, 1, 1, 0, True)
        TMP2 = ants.morphology(C, operation='erode', radius=10, mtype='binary', shape='ball')
        G[G == 0] = TMP2[G == 0]
        G = G * 2
        seg_tmp2 = W
        seg_tmp2[W == 0] = G[W == 0]

        #  clean the Brainmask
        tmp_mask3 = ants.threshold_image(seg_tmp2, 3, 3, 1, 0, True)
        TMP3 = ants.threshold_image(seg_tmp2, 2, 2, 1, 0, True)
        tmp_mask3[tmp_mask3 == 0] = TMP3[tmp_mask3 == 0]
        tmp_mask3 = ants.morphology(tmp_mask3, operation='erode', radius=2, mtype='binary', shape='ball')
        tmp_mask3 = ants.iMath(tmp_mask3, operation='GetLargestComponent')
        tmp_mask3 = ants.morphology(tmp_mask3, operation='dilate', radius=4, mtype='binary', shape='ball')
        tmp_mask3 = ants.iMath(tmp_mask3, 'FillHoles', 2)

        # tmp_mask1 = ants.iMath(tmp_mask1,'Pad',10)
        # tmp_mask3[tmp_mask3==0]=tmp_mask1[tmp_mask3==0]
        tmp_mask3 = ants.morphology(tmp_mask3, operation='dilate', radius=5, mtype='binary', shape='ball')
        tmp_mask3 = ants.morphology(tmp_mask3, operation='erode', radius=5, mtype='binary', shape='ball')
        tmp_mask3 = ants.iMath(tmp_mask3, 'Pad', -10)

        ants.image_write(tmp_mask3, output_for_mask, ri=False)

        if check_visualy_final_mask == 'YES':
            command = 'singularity run' + s_bind + fs_sif + 'freeview - v' + input_for_msk + \
                      ' ' + output_for_mask + ': colormap = heat:opacity = 0.5:visible = 1'
            spco([command], shell=True)
            NEW_mask = ants.image_read(output_for_mask)
            NEW_mask = ants.iMath(NEW_mask, operation='GetLargestComponent')
            NEW_mask = ants.smooth_image(NEW_mask, 0.7)
            NEW_mask = ants.threshold_image(NEW_mask, 0.5, 1, 1, 0, True)
            ants.image_write(NEW_mask, output_for_mask, ri=False)




    elif brain_skullstrip =='Custum_ANTS_NL':

        IMG = ants.image_read(input_for_msk)
        mTx = ants.registration(fixed=IMG, moving=BASE_SS_coregistr,
                                type_of_transform='SyNCC',
                                outprefix=opj(dir_transfo, 'template_to_' + masking_img + '_SyN_'),verbose=True)

        REF_MASK = ants.image_read(BASE_SS_mask)
        tmp_mask1 = ants.apply_transforms(fixed=IMG, moving=REF_MASK,
                                          transformlist=mTx['fwdtransforms'], interpolator='nearestNeighbor')

        tmp_mask1 = ants.threshold_image(tmp_mask1, 0.5, 1, 1, 0, True)
        tmp_mask1 = ants.morphology(tmp_mask1, operation='dilate', radius=2, mtype='binary', shape='ball')
        tmp_mask1 = ants.iMath(tmp_mask1, operation='GetLargestComponent')

        seg_tmp = ants.atropos(a=IMG, m='[0.1,1x1x1]', c='[3,0]', i='kmeans[3]', x=tmp_mask1)

        #  Clean up:
        seg_tmp2 = ants.iMath(seg_tmp['segmentation'], 'Pad', 10)
        W = ants.threshold_image(seg_tmp2, 3, 3, 1, 0, True)
        W = ants.iMath(W, operation='GetLargestComponent')
        W = W * 3
        G = ants.threshold_image(seg_tmp2, 2, 2, 1, 0, True)
        G = ants.iMath(G, operation='GetLargestComponent')
        TMP1 = ants.iMath(G, 'FillHoles', 2)
        G = G * TMP1
        C = ants.threshold_image(seg_tmp2, 1, 1, 1, 0, True)
        TMP2 = ants.morphology(C, operation='erode', radius=10, mtype='binary', shape='ball')
        G[G == 0] = TMP2[G == 0]
        G = G * 2
        seg_tmp2 = W
        seg_tmp2[W == 0] = G[W == 0]

        #  clean the Brainmask
        tmp_mask3 = ants.threshold_image(seg_tmp2, 3, 3, 1, 0, True)
        TMP3 = ants.threshold_image(seg_tmp2, 2, 2, 1, 0, True)
        tmp_mask3[tmp_mask3 == 0] = TMP3[tmp_mask3 == 0]
        tmp_mask3 = ants.morphology(tmp_mask3, operation='erode', radius=2, mtype='binary', shape='ball')
        tmp_mask3 = ants.iMath(tmp_mask3, operation='GetLargestComponent')
        tmp_mask3 = ants.morphology(tmp_mask3, operation='dilate', radius=4, mtype='binary', shape='ball')
        tmp_mask3 = ants.iMath(tmp_mask3, 'FillHoles', 2)

        # tmp_mask1 = ants.iMath(tmp_mask1,'Pad',10)
        # tmp_mask3[tmp_mask3==0]=tmp_mask1[tmp_mask3==0]
        tmp_mask3 = ants.morphology(tmp_mask3, operation='dilate', radius=5, mtype='binary', shape='ball')
        tmp_mask3 = ants.morphology(tmp_mask3, operation='erode', radius=5, mtype='binary', shape='ball')
        tmp_mask3 = ants.iMath(tmp_mask3, 'Pad', -10)

        ants.image_write(tmp_mask3, output_for_mask, ri=False)

        if check_visualy_final_mask == 'YES':
            command = 'singularity run' + s_bind + fs_sif + 'freeview - v' + input_for_msk + \
                      ' ' + output_for_mask + ': colormap = heat:opacity = 0.5:visible = 1'
            spco([command], shell=True)
            NEW_mask = ants.image_read(output_for_mask)
            NEW_mask = ants.iMath(NEW_mask, operation='GetLargestComponent')
            NEW_mask = ants.smooth_image(NEW_mask, 0.7)
            NEW_mask = ants.threshold_image(NEW_mask, 0.5, 1, 1, 0, True)
            ants.image_write(NEW_mask, output_for_mask, ri=False)



    elif brain_skullstrip =='Custum_ANTS':

        IMG = ants.image_read(input_for_msk)
        mTx = ants.registration(fixed=IMG, moving=BASE_SS_coregistr,
                                type_of_transform='Affine',
                                outprefix=opj(dir_transfo, 'template_to_' + masking_img + '_SyN_'), verbose=True)

        REF_MASK = ants.image_read(BASE_SS_mask)
        tmp_mask1 = ants.apply_transforms(fixed=IMG, moving=REF_MASK,
                                          transformlist=mTx['fwdtransforms'], interpolator='nearestNeighbor')

        tmp_mask1 = ants.threshold_image(tmp_mask1, 0.5, 1, 1, 0, True)
        tmp_mask1 = ants.morphology(tmp_mask1, operation='dilate', radius=2, mtype='binary', shape='ball')
        tmp_mask1 = ants.iMath(tmp_mask1, operation='GetLargestComponent')

        seg_tmp = ants.atropos(a=IMG, m='[0.1,1x1x1]', c='[3,0]', i='kmeans[3]', x=tmp_mask1)

        #  Clean up:
        seg_tmp2 = ants.iMath(seg_tmp['segmentation'], 'Pad', 10)
        W = ants.threshold_image(seg_tmp2, 3, 3, 1, 0, True)
        W = ants.iMath(W, operation='GetLargestComponent')
        W = W * 3
        G = ants.threshold_image(seg_tmp2, 2, 2, 1, 0, True)
        G = ants.iMath(G, operation='GetLargestComponent')
        TMP1 = ants.iMath(G, 'FillHoles', 2)
        G = G * TMP1
        C = ants.threshold_image(seg_tmp2, 1, 1, 1, 0, True)
        TMP2 = ants.morphology(C, operation='erode', radius=10, mtype='binary', shape='ball')
        G[G == 0] = TMP2[G == 0]
        G = G * 2
        seg_tmp2 = W
        seg_tmp2[W == 0] = G[W == 0]

        #  clean the Brainmask
        tmp_mask3 = ants.threshold_image(seg_tmp2, 3, 3, 1, 0, True)
        TMP3 = ants.threshold_image(seg_tmp2, 2, 2, 1, 0, True)
        tmp_mask3[tmp_mask3 == 0] = TMP3[tmp_mask3 == 0]
        tmp_mask3 = ants.morphology(tmp_mask3, operation='erode', radius=2, mtype='binary', shape='ball')
        tmp_mask3 = ants.iMath(tmp_mask3, operation='GetLargestComponent')
        tmp_mask3 = ants.morphology(tmp_mask3, operation='dilate', radius=4, mtype='binary', shape='ball')
        tmp_mask3 = ants.iMath(tmp_mask3, 'FillHoles', 2)

        # tmp_mask1 = ants.iMath(tmp_mask1,'Pad',10)
        # tmp_mask3[tmp_mask3==0]=tmp_mask1[tmp_mask3==0]
        tmp_mask3 = ants.morphology(tmp_mask3, operation='dilate', radius=5, mtype='binary', shape='ball')
        tmp_mask3 = ants.morphology(tmp_mask3, operation='erode', radius=5, mtype='binary', shape='ball')
        tmp_mask3 = ants.iMath(tmp_mask3, 'Pad', -10)

        ants.image_write(tmp_mask3, output_for_mask, ri=False)

        if check_visualy_final_mask == 'YES':
            command = 'singularity run' + s_bind + fs_sif + 'freeview - v' + input_for_msk + \
                      ' ' + output_for_mask + ': colormap = heat:opacity = 0.5:visible = 1'
            spco([command], shell=True)
            NEW_mask = ants.image_read(output_for_mask)
            NEW_mask = ants.iMath(NEW_mask, operation='GetLargestComponent')
            NEW_mask = ants.smooth_image(NEW_mask, 0.7)
            NEW_mask = ants.threshold_image(NEW_mask, 0.5, 1, 1, 0, True)
            ants.image_write(NEW_mask, output_for_mask, ri=False)


    elif brain_skullstrip == '3dSkullStrip':

        command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite ' + \
            '-input ' + input_for_msk + ' -blur_fwhm 2 -orig_vol -mask_vol -use_skull -monkey'
        spco(command, shell=True)

        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 2'
        spco(command, shell=True)

    elif brain_skullstrip == 'Custum_1':
        loadimg = nib.load(opj(volumes_dir, ID + '_' + otheranat + '_template.nii.gz')).get_data() 
        loadimgsort85 =  np.percentile(np.abs(loadimg)[np.abs(loadimg)>0], 40)
        mask_imag = nilearn.image.threshold_img(opj(volumes_dir, ID + '_' + otheranat + '_template.nii.gz'), loadimgsort85, cluster_threshold=10)
        mask_imag.to_filename(output_for_mask)
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 1'
        spco(command, shell=True)
        command = 'singularity run' + s_bind + afni_sif + '3dresample -master ' + input_for_msk + ' -prefix' + output_for_mask + ' -input ' + output_for_mask + ' -overwrite -bound_type SLAB'
        spco(command, shell=True)

    elif brain_skullstrip == 'Custum_Baboon':
        #convert to float
        command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
        spco([command], shell=True)
        mask_img = compute_epi_mask(input_for_msk[:-7] + '_float.nii.gz', lower_cutoff=0.85, upper_cutoff=0.95, connected=True, opening=3,
            exclude_zeros=False, ensure_finite=True)
        mask_img.to_filename(output_for_mask)
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 15'
        spco(command, shell=True)

    elif brain_skullstrip == 'Custum_Macaque':
        #convert to float
        command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
        spco([command], shell=True)
        mask_img = compute_epi_mask(input_for_msk[:-7] + '_float.nii.gz', lower_cutoff=0.2, upper_cutoff=0.90, connected=True, opening=3,
            exclude_zeros=False, ensure_finite=True)
        mask_img.to_filename(output_for_mask)
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 8 -2'
        spco(command, shell=True)


    elif brain_skullstrip == 'Custum_mouse':
        #convert to float
        command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
        spco([command], shell=True)
        mask_img = compute_epi_mask(input_for_msk[:-7] + '_float.nii.gz', lower_cutoff=0.2, upper_cutoff=0.80, connected=True, opening=3,
            exclude_zeros=False, ensure_finite=True)
        mask_img.to_filename(output_for_mask)
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 3 -1'
        spco(command, shell=True)

    elif brain_skullstrip == 'Custum_dog':
        #convert to float
        command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
        spco([command], shell=True)
        mask_img = compute_epi_mask(input_for_msk[:-7] + '_float.nii.gz', lower_cutoff=0.85, upper_cutoff=0.95, connected=True, opening=3,
            exclude_zeros=False, ensure_finite=True)
        mask_img.to_filename(output_for_mask)
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
        ' -input ' + output_for_mask + ' -fill_holes -dilate_input 23'
        spco(command, shell=True)

    elif brain_skullstrip =='bet2':
        #####creat an approximate brain mask
        command = 'singularity run' + s_bind + fsl_sif + 'bet2 ' + input_for_msk + ' ' + opj(masks_dir, ID + '_bet' + masking_img + '.nii.gz') + \
        ' -f 0.70'
        spco([command], shell=True)
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + opj(masks_dir, ID + '_bet' + masking_img + '.nii.gz') + ' -expr "step(a)" -prefix ' + output_for_mask + ' -overwrite'
        spco(command, shell=True)

    elif brain_skullstrip == 'Custum T1/T2':
        command = 'singularity run' + s_bind + afni_sif + '3dresample -master ' + input_for_msk + ' -prefix ' + opj(volumes_dir, ID + '_' + otheranat + '_template_RSPL.nii.gz') + ' -input ' + opj(volumes_dir, ID + '_' + otheranat + '_template.nii.gz') + ' -overwrite'
        spco([command], shell=True)
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + input_for_msk + ' -b ' + opj(volumes_dir, ID + '_' + otheranat + '_template_RSPL.nii.gz') + ' -expr b*(step(a)) -prefix ' + opj(volumes_dir, ID + '_' + otheranat + '_template_RSPL.nii.gz') + ' -overwrite'
        spco([command], shell=True)
        log_img = math_img("np.where((img1 > 1) & (img2 > 1), img2, 0)", img1=input_for_msk, img2=opj(volumes_dir, ID + '_' + otheranat + '_template_RSPL.nii.gz'))
        log_img.to_filename(opj(volumes_dir, ID + '_' + otheranat + '_template_RSPL.nii.gz'))
        log_img = math_img("np.where((img1 > 1) & (img2 > 1), img1, 0)", img1=input_for_msk, img2=opj(volumes_dir, ID + '_' + otheranat + '_template_RSPL.nii.gz'))
        log_img.to_filename(opj(volumes_dir, ID + '_' + masking_img + '_template_RSPL.nii.gz'))
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + opj(volumes_dir, ID + '_' + masking_img + '_template_RSPL.nii.gz') + \
        ' -b ' + opj(volumes_dir, ID + '_' + otheranat + '_template_RSPL.nii.gz') + \
        ' -expr "a/b" -datum float  -prefix ' + opj(volumes_dir, ID + '_' + masking_img + '_' + otheranat + '_template.nii.gz')
        spco([command], shell=True)

        IMG = ants.image_read(opj(volumes_dir, ID + '_' + masking_img + '_' + otheranat + '_template.nii.gz'))
        mTx = ants.registration(fixed=IMG, moving=BASE_SS_coregistr,
                                type_of_transform='SyNCC',
                                outprefix=opj(dir_transfo, 'template_to_' + masking_img + '_SyN_'), verbose=True)

        REF_MASK = ants.image_read(BASE_SS_mask)
        tmp_mask1 = ants.apply_transforms(fixed=IMG, moving=REF_MASK,
                                          transformlist=mTx['fwdtransforms'], interpolator='nearestNeighbor')

        tmp_mask1 = ants.threshold_image(tmp_mask1, 0.5, 1, 1, 0, True)
        tmp_mask1 = ants.morphology(tmp_mask1, operation='dilate', radius=2, mtype='binary', shape='ball')
        tmp_mask1 = ants.iMath(tmp_mask1, operation='GetLargestComponent')

        seg_tmp = ants.atropos(a=IMG, m='[0.1,1x1x1]', c='[3,0]', i='kmeans[3]', x=tmp_mask1)

        #  Clean up:
        seg_tmp2 = ants.iMath(seg_tmp['segmentation'], 'Pad', 10)
        W = ants.threshold_image(seg_tmp2, 3, 3, 1, 0, True)
        W = ants.iMath(W, operation='GetLargestComponent')
        W = W * 3
        G = ants.threshold_image(seg_tmp2, 2, 2, 1, 0, True)
        G = ants.iMath(G, operation='GetLargestComponent')
        TMP1 = ants.iMath(G, 'FillHoles', 2)
        G = G * TMP1
        C = ants.threshold_image(seg_tmp2, 1, 1, 1, 0, True)
        TMP2 = ants.morphology(C, operation='erode', radius=10, mtype='binary', shape='ball')
        G[G == 0] = TMP2[G == 0]
        G = G * 2
        seg_tmp2 = W
        seg_tmp2[W == 0] = G[W == 0]

        #  clean the Brainmask
        tmp_mask3 = ants.threshold_image(seg_tmp2, 3, 3, 1, 0, True)
        TMP3 = ants.threshold_image(seg_tmp2, 2, 2, 1, 0, True)
        tmp_mask3[tmp_mask3 == 0] = TMP3[tmp_mask3 == 0]
        tmp_mask3 = ants.morphology(tmp_mask3, operation='erode', radius=2, mtype='binary', shape='ball')
        tmp_mask3 = ants.iMath(tmp_mask3, operation='GetLargestComponent')
        tmp_mask3 = ants.morphology(tmp_mask3, operation='dilate', radius=4, mtype='binary', shape='ball')
        tmp_mask3 = ants.iMath(tmp_mask3, 'FillHoles', 2)

        # tmp_mask1 = ants.iMath(tmp_mask1,'Pad',10)
        # tmp_mask3[tmp_mask3==0]=tmp_mask1[tmp_mask3==0]
        tmp_mask3 = ants.morphology(tmp_mask3, operation='dilate', radius=5, mtype='binary', shape='ball')
        tmp_mask3 = ants.morphology(tmp_mask3, operation='erode', radius=5, mtype='binary', shape='ball')
        tmp_mask3 = ants.iMath(tmp_mask3, 'Pad', -10)

        ants.image_write(tmp_mask3, opj(volumes_dir, ID + '_' + masking_img + '_' + otheranat + '_template_msk.nii.gz'), ri=False)

        mask_img = compute_epi_mask(opj(volumes_dir, ID + '_' + masking_img + '_' + otheranat + '_template.nii.gz'), lower_cutoff=0.10, upper_cutoff=0.60, connected=True, opening=3)
        mask_img.to_filename(output_for_mask)

    elif brain_skullstrip =='QWARP':

        command = 'singularity run' + s_bind + afni_sif + '3dQwarp -overwrite -iwarp' + \
        ' -base ' + BASE_SS_coregistr + \
        ' -prefix ' + opj(dir_prepro,'template_to_' + masking_img + '_AFNIQ.nii.gz') + \
        ' -source ' + input_for_msk + ' -maxlev 5 -resample'
        spco(command, shell=True)

        command = 'singularity run' + s_bind + afni_sif + '3dNwarpApply -nwarp ' + opj(dir_prepro,'template_to_' + masking_img + '_AFNIQ_WARPINV.nii.gz') + \
        ' -source ' + BASE_SS_mask + ' -master ' + input_for_msk + ' -interp NN' + \
        ' -prefix ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + ' -overwrite'
        spco(command, shell=True)

        Ex_Mask    = opj(dir_prepro,'mask_tmp' + masking_img + '.nii.gz')

        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + Ex_Mask + \
        ' -input ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + ' -fill_holes' # -dilate_input 2'
        spco(command, shell=True)

        shutil.copyfile(Ex_Mask, output_for_mask)

        if check_visualy_final_mask == True:
            command = 'singularity run' + s_bind + fs_sif + 'freeview - v' + input_for_msk + \
                      ' ' + output_for_mask + ': colormap = heat:opacity = 0.5:visible = 1'
            spco([command], shell=True)

    elif brain_skullstrip == '3dSkullStrip_Rat':
        if ID in ['301502','302101','302105','302106','301603', '300908','301500','301501','301503','301504','301505','301508','301509','302107','302108','300600']:
            command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite -input ' + input_for_msk + ' -orig_vol -mask_vol -rat'
            spco([command], shell=True)

            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes -dilate_input 15'
            spco(command, shell=True)
        else:
            command = 'singularity run' + s_bind + afni_sif + '3dSkullStrip -prefix ' + output_for_mask + ' -overwrite -input ' + input_for_msk + ' -orig_vol -mask_vol -surface_coil -rat'
            spco([command], shell=True)

            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
            ' -input ' + output_for_mask + ' -fill_holes -dilate_input 4'
            spco(command, shell=True)

    elif brain_skullstrip == 'NoSkullStrip':
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + input_for_msk + \
                  ' -expr "step(a)" -prefix ' + output_for_mask + ' -overwrite'
        spco(command, shell=True)

    elif brain_skullstrip =='sammba_rat':
        command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
        spco([command], shell=True)

        nichols_masker = Histrogram_mask_EMB()
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
            command = 'singularity run' + s_bind + fs_sif + 'freeview - v' + input_for_msk + \
                      ' ' + output_for_mask + ': colormap = heat:opacity = 0.5:visible = 1'
            spco([command], shell=True)

    elif brain_skullstrip =='custum_rat':
            #convert to float
            command = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + input_for_msk + ' ' + input_for_msk[:-7] + '_float.nii.gz'
            spco([command], shell=True)
            mask_img = compute_epi_mask(input_for_msk[:-7] + '_float.nii.gz', lower_cutoff=0.8, upper_cutoff=0.85, connected=True, opening=3,
                                        exclude_zeros=True, ensure_finite=True)
            mask_img.to_filename(output_for_mask)
            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + output_for_mask + \
                      ' -input ' + output_for_mask + ' -fill_holes -dilate_input 3'
            spco(command, shell=True)

        ################################################
        ###### in use for study don't touch ############
        ################################################

    elif brain_skullstrip =='Custum_QWARP':

        if ID in ['Oliver'] and Session in [3] or ID in ['Roshan'] and Session in [3] or ID in ['Quantum'] and Session in [3] or ID in ['Roshan'] and Session in [6]:

            command = 'singularity run' + s_bind + afni_sif + '3dQwarp -overwrite -iwarp' + \
            ' -base ' + BASE_SS_coregistr + \
            ' -prefix ' + opj(dir_prepro,'template_to_' + masking_img + '_AFNIQ.nii.gz') + \
            ' -source ' + input_for_msk + ' -maxlev 5 -lpa -resample'
            spco(command, shell=True)

            command = 'singularity run' + s_bind + afni_sif + '3dNwarpApply -nwarp ' + opj(dir_prepro,'template_to_' + masking_img + '_AFNIQ_WARPINV.nii.gz') + \
            ' -source ' + BASE_SS_mask + ' -master ' + input_for_msk + ' -interp NN' + \
            ' -prefix ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + ' -overwrite'
            spco(command, shell=True)

            Ex_Mask    = opj(dir_prepro,'mask_tmp' + masking_img + '.nii.gz')

            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + Ex_Mask + \
            ' -input ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + ' -fill_holes' # -dilate_input 2'
            spco(command, shell=True)

            shutil.copyfile(Ex_Mask,output_for_mask)

            if check_visualy_final_mask == True:
                command = 'singularity run' + s_bind + fs_sif + 'freeview - v' + input_for_msk + \
                ' ' + output_for_mask + ': colormap = heat:opacity = 0.5:visible = 1'
                spco([command], shell=True)

        else:

            command = 'singularity run' + s_bind + afni_sif + '3dQwarp -overwrite -iwarp' + \
            ' -base ' + BASE_SS_coregistr + \
            ' -prefix ' + opj(dir_prepro,'template_to_' + masking_img + '_AFNIQ.nii.gz') + \
            ' -source ' + input_for_msk + ' -maxlev 5 -resample'
            spco(command, shell=True)

            command = 'singularity run' + s_bind + afni_sif + '3dNwarpApply -nwarp ' + opj(dir_prepro,'template_to_' + masking_img + '_AFNIQ_WARPINV.nii.gz') + \
            ' -source ' + BASE_SS_mask + ' -master ' + input_for_msk + ' -interp NN' + \
            ' -prefix ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + ' -overwrite'
            spco(command, shell=True)

            Ex_Mask    = opj(dir_prepro,'mask_tmp' + masking_img + '.nii.gz')

            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + Ex_Mask + \
            ' -input ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + ' -fill_holes' # -dilate_input 2'
            spco(command, shell=True)

            shutil.copyfile(Ex_Mask,output_for_mask)

            if check_visualy_final_mask == True:
                command = 'singularity run' + s_bind + fs_sif + 'freeview - v' + input_for_msk + \
                ' ' + output_for_mask + ': colormap = heat:opacity = 0.5:visible = 1'
                spco([command], shell=True)


    elif brain_skullstrip =='Custum_QWARPT2':

        command = 'singularity run' + s_bind + afni_sif + '3dQwarp -overwrite -lpa -iwarp' + \
        ' -base ' + BASE_SS_coregistr + \
        ' -prefix ' + opj(dir_prepro,'template_to_' + masking_img + '_AFNIQ.nii.gz') + \
        ' -source ' + input_for_msk + ' -maxlev 3 -resample'
        spco(command, shell=True)

        command = 'singularity run' + s_bind + afni_sif + '3dNwarpApply -nwarp ' + opj(dir_prepro,'template_to_' + masking_img + '_AFNIQ_WARPINV.nii.gz') + \
        ' -source ' + BASE_SS_mask + ' -master ' + input_for_msk + ' -interp NN' + \
        ' -prefix ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + ' -overwrite'
        spco(command, shell=True)

        Ex_Mask    = opj(dir_prepro,'mask_tmp' + masking_img + '.nii.gz')

        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + Ex_Mask + \
        ' -input ' + opj(dir_prepro,masking_img + 'template_brainmask.nii.gz') + ' -fill_holes'
        spco(command, shell=True)

        shutil.copyfile(Ex_Mask,output_for_mask)

        if check_visualy_final_mask == True:
            command = 'singularity run' + s_bind + fs_sif + 'freeview - v' + input_for_msk + \
            ' ' + output_for_mask + ': colormap = heat:opacity = 0.5:visible = 1'
            spco([command], shell=True)


    else:
        print("ERROR in brain_skullstrip name???")

    return(output_for_mask)