import os
import subprocess
from nilearn.image import resample_to_img
from fonctions.extract_filename import extract_filename
import ants

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput



def Refimg_to_meanfMRI(SED, anat_func_same_space, TfMRI, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, RS, nb_run, ID, dir_prepro, n_for_ANTS, list_atlases, labels_dir, anat_subject, IhaveanANAT, do_anat_to_func, type_of_transform, overwrite, s_bind, afni_sif):

    command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -di ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    ADI = spgo(command).split('\n')
    delta_x = str(abs(round(float(ADI[-1]), 10)))
    command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -dj ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    ADJ = spgo(command).split('\n')
    delta_y= str(abs(round(float(ADJ[-1]), 10)))
    command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -dk ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    ADK = spgo(command).split('\n')
    delta_z = str(abs(round(float(ADK[-1]), 10)))

    if ope(opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz')):
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz') + \
                  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.nii.gz') + ' -expr "a"'
        spco([command], shell=True)

    #### take opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat_Allineate_in_func.nii.gz') + ' -b ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
    ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz') + ' -expr "a*b"'
    spco([command], shell=True)

    ####################################################################################
    ########################## use template and transfo to anat (average indiv anat)  ##
    ####################################################################################
    #average indiv anat =>coreg_ref
    # mean of all func SS => opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS.nii.gz')


    if anat_func_same_space == True and do_anat_to_func == False:
        print('No anat to func step required')
    else:

        if 'i' in SED:
            restrict = (1,0.1,0.1)
        elif 'j' in SED:
            restrict = (0.1,1,0.1)
        elif 'k' in SED:
            restrict = (0.1,0.1,1)
        elif 'None' in SED:
            restrict = (1, 1, 1)

        MEAN = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz'))
        ANAT = ants.image_read(opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.nii.gz'))
        mtx1 = ants.registration(fixed=ANAT, moving=MEAN,type_of_transform='Translation', outprefix=opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_shift_'))
        MEAN_tr = ants.apply_transforms(fixed=ANAT, moving=MEAN,transformlist=mtx1['fwdtransforms'],interpolator=n_for_ANTS)
        ants.image_write(MEAN_tr, opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_shift.nii.gz'), ri=False)
        mTx2 = ants.registration(fixed=ANAT, moving=MEAN,
                                 type_of_transform=type_of_transform,
                                 initial_transform=mtx1['fwdtransforms'],
                                 outprefix=opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_'),
                                 grad_step=0.1,flow_sigma=3,total_sigma=0,aff_sampling=32,
                                 aff_random_sampling_rate=0.2,
                                 syn_sampling=32,
                                 aff_iterations=(1000, 500, 250, 100),
                                 aff_shrink_factors=(8, 4, 2, 1),
                                 aff_smoothing_sigmas=(3, 2, 1, 0),
                                 reg_iterations=(1000, 500, 250, 100),
                                 reg_smoothing_sigmas=(3, 2, 1, 0),
                                 reg_shrink_factors=(8, 4, 2, 1),
                                 verbose=True,
                                 restrict_transformation=restrict)
        MEAN_tr = ants.apply_transforms(fixed=ANAT, moving=MEAN_tr,transformlist=mTx2['fwdtransforms'],interpolator=n_for_ANTS)
        ants.image_write(MEAN_tr, opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped.nii.gz'), ri=False)



    mvt_shft_INV_ANTs = []
    w2inv_inv = []
    for elem1, elem2 in zip([#opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_shift_0GenericAffine.mat'),
                             opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_0GenericAffine.mat'),
                             opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_1Warp.nii.gz')],
                            [True,True]):
        if ope(elem1):
            mvt_shft_INV_ANTs.append(elem1)
            w2inv_inv.append(elem2)

    print(mvt_shft_INV_ANTs)
    print(w2inv_inv)

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
        if ope(input1):
            # mask
            MEAN = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz'))
            IMG = ants.image_read(input1)
            TRANS = ants.apply_transforms(fixed=MEAN, moving=IMG,
                                          transformlist=mvt_shft_INV_ANTs,
                                          interpolator='nearestNeighbor',
                                          which2invert=w2inv_inv)
            ants.image_write(TRANS,output2,ri=False)


            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool' + overwrite + ' -prefix ' + output2 + \
            ' -input ' + output2 + ' -fill_holes'
            spco(command, shell=True)

            command = 'singularity run' + s_bind + afni_sif + '3dclust -NN1 10 -prefix ' + output2 + output2
            spco(command, shell=True)
        else:
            print('WARNING:' + str(input1) + ' not found!!! this may be because you have not provided an aseg file, then no extraction of WM or Ventricules or GM will be possible... pls check that!!!!')

    ## in func space resample to func
    BRAIN = ants.image_read(opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz'))
    TRANS = ants.apply_transforms(fixed=MEAN, moving=BRAIN,
                                  transformlist=mvt_shft_INV_ANTs,
                                  interpolator=n_for_ANTS,
                                  which2invert=w2inv_inv)
    ants.image_write(TRANS, opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI.nii.gz'), ri=False)

    #### apply to all atlases

    if len(list_atlases) > 0:
        for atlas in list_atlases:
            ## in anat space resample to func

            if anat_func_same_space == True:

                mvt_shft = opj(dir_prepro, ID + '_brain_for_Align_Center_inv.1D')
                command = 'singularity run' + s_bind + afni_sif + '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + ' ' + opj(labels_dir, TfMRI + opb(atlas)) + ' -overwrite'
                spco([command], shell=True)

                command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -final NN -overwrite -1Dmatrix_apply ' + mvt_shft + \
                ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + \
                ' -master ' + opj(dir_prepro, ID + '_mprage_reorient' + TfMRI + '.nii.gz') + \
                ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas))
                spco([command], shell=True)

                #command = '3drefit -atrcopy ' + opj(labels_dir, TfMRI + opb(atlas)) + ' IJK_TO_DICOM_REAL ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas))
                #spco([command], shell=True)
                caca = resample_to_img(opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)), opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz'),  interpolation='nearest')
                caca.to_filename(opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)))

            elif IhaveanANAT == False:
                command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
                ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + \
                ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                ' -input  ' + atlas
                spco([command], shell=True)
            else:
                command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
                ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + \
                ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                ' -input  ' + opj(labels_dir, TfMRI + opb(atlas))
                spco([command], shell=True)

            ## in func space resample to func

            ATLAS = ants.image_read(opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)))
            TRANS = ants.apply_transforms(fixed=MEAN, moving=ATLAS,
                                          transformlist=mvt_shft_INV_ANTs,
                                          interpolator='nearestNeighbor',
                                          which2invert=w2inv_inv)
            ants.image_write(TRANS, opj(dir_fMRI_Refth_RS_prepro1, opb(atlas)), ri=False)
    else:
        print('WARNING: list_atlases is empty!')

    command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -di ' + anat_subject
    ADI = spgo(command).split('\n')[-1]
    delta_x1 = str(abs(round(float(ADI), 10)))
    command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -dj ' + anat_subject
    ADJ = spgo(command).split('\n')
    delta_y1= str(abs(round(float(ADJ[-1]), 10)))
    command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -dk ' + anat_subject
    ADK = spgo(command).split('\n')
    delta_z1 = str(abs(round(float(ADK[-1]), 10)))

    command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS_anat_resolution.nii.gz') + \
    ' -dxyz ' + delta_x1 + ' ' + delta_y1 + ' ' + delta_z1 + ' ' + \
    ' -rmode Cu -input  ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    spco([command], shell=True)


    #### creat a nice anat in func space
    if anat_func_same_space == True:
        mvt_shft = opj(dir_prepro, ID + '_brain_for_Align_Center_inv.1D')
        anatstd = opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_plot.nii.gz')
        ##### apply the recenter fmri
        command = 'singularity run' + s_bind + afni_sif + '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + anatstd + ' ' + anat_subject + ' -overwrite'
        spco([command], shell=True)

        command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -overwrite -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + anatstd + \
        ' -input  ' + anatstd + \
        ' -master ' + opj(dir_prepro, ID + '_mprage_reorient' + TfMRI + '.nii.gz')

        spco([command], shell=True)

    else:
        anatstd = anat_subject
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + anatstd + \
        ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_plot.nii.gz') + ' -expr "a"'
        spco([command], shell=True)


    ## in func space resample to func 

    ANAT = ants.image_read(anatstd)
    MEAN_RES = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS_anat_resolution.nii.gz'))
    TRANS = ants.apply_transforms(fixed=MEAN_RES, moving=ANAT,
                                  transformlist=mvt_shft_INV_ANTs,
                                  interpolator=n_for_ANTS,
                                  which2invert=w2inv_inv)
    ants.image_write(TRANS, opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI_anat_resolution.nii.gz'), ri=False)


    ###finaly mask the func with mask
    for i in range(0, int(nb_run)):
        root_RS = extract_filename(RS[i])
        ### take opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
        
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' +   opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz') + \
                  ' -b ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref.nii.gz') + \
                  ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref_RcT_masked.nii.gz') + ' -expr "a*b"'
        spco([command], shell=True)
