import os
import subprocess
from nilearn.image import resample_to_img
from fonctions.extract_filename import extract_filename
import ants
import datetime
import json
from fonctions.plot_QC_func import plot_qc
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
spco = subprocess.check_output
spgo = subprocess.getoutput

def Refimg_to_meanfMRI(REF_int, SED, anat_func_same_space, TfMRI, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, RS, nb_run, ID, dir_prepro,
                       n_for_ANTS, aff_metric_ants, list_atlases, labels_dir, anat_subject, IhaveanANAT, do_anat_to_func, type_of_transform,
                       overwrite, s_bind, afni_sif,diary_file):

    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  Working on step ' + str(5) + '(function: _5_anat_to_fMRI).  ##'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

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
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"ADD_from_Sources": opj(dir_fMRI_Refth_RS_prepro1, 'manual_mask.nii.gz'),
                      "ADD_Description": 'Copy.', },
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.json'), "a") as outfile:
            outfile.write(json_object)

    #### apply skull stripping
    command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat_Allineate_in_func.nii.gz') + \
              ' -b ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
    ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz') + ' -expr "a*b"'
    nl = spgo(command)
    diary.write(f'\n{nl}')
    print(nl)
    dictionary = {"ADD_Sources": [opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz'),
                              opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.nii.gz')],
                  "ADD_Description": 'Skull strippng (3dcalc, AFNI).', },
    # Path to the JSON file
    json_file_path = opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.json')
    # Load existing JSON data
    with open(json_file_path, "r") as infile:
        existing_data = json.load(infile)
    # Update the existing data with the new dictionary
    existing_data.update(dictionary)
    # Save the updated content back to the file
    with open(json_file_path, "w") as outfile:
        json.dump(existing_data, outfile, indent=2)

    ####################################################################################
    ########################## use template and transfo to anat (average indiv anat)  ##
    ####################################################################################

    if anat_func_same_space == True and do_anat_to_func == False:
        nl = 'No anat to func step required'
        print(bcolors.OKGREEN + nl + bcolors.ENDC)
        diary.write(f'\n{nl}')

        command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz' ) \
        + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped.nii.gz') + ' -expr "a"' + overwrite
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

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
        dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz'),
                                  opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.nii.gz')],
                      "Description": 'Co-registration (translation, ANTspy).', },
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_shift.json'), "w") as outfile:
            outfile.write(json_object)

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
                                 aff_metric=aff_metric_ants,
                                 restrict_transformation=restrict)

        transfo_concat = \
            [opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_1Warp.nii.gz'),
             opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_0GenericAffine.mat')]

        MEAN_trAff = ants.apply_transforms(fixed=ANAT, moving=MEAN,transformlist=transfo_concat, interpolator='nearestNeighbor')
        ants.image_write(MEAN_trAff, opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped.nii.gz'), ri=False)

        dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz'),
                                  opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.nii.gz'),
                                  mtx1['fwdtransforms']],
                      "Description": 'Co-registration (Non linear, ANTspy).', },
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped.json'), "w") as outfile:
            outfile.write(json_object)


    if do_anat_to_func == True:
        mvt_shft_INV_ANTs = []
        w2inv_inv = []
        for elem1, elem2 in zip([  # opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_shift_0GenericAffine.mat'),
            opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_0GenericAffine.mat'),
            opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_1InverseWarp.nii.gz')],
                [True, False]):
            if ope(elem1):
                mvt_shft_INV_ANTs.append(elem1)
                w2inv_inv.append(elem2)
    elif do_anat_to_func == False and anat_func_same_space == True:
        mvt_shft_INV_ANTs = []
        w2inv_inv = []
    else:
        nl = 'ERROR: If Anat and Func are not in the same space you need to perform that trasnformation (do_anat_to_func = True)'
        diary.write(f'\n{nl}')
        raise Exception(bcolors.FAIL + nl + bcolors.ENDC)

    # doesn't work for two different 1d matrices... so let's do it separately....

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
                                          whichtoinvert=w2inv_inv)
            ants.image_write(TRANS,output2,ri=False)

            command = 'singularity run' + s_bind + afni_sif + '3dmask_tool' + overwrite + ' -prefix ' + output2 + \
            ' -input ' + output2 + ' -fill_holes'
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

            command = 'singularity run' + s_bind + afni_sif + '3dclust -NN1 10 -prefix ' + output2 + output2
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

            dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz'),
                                      input1],
                          "Description": ['1. Normalization (nearestNeighbo,ANTspy).',
                                          '2. fill holes (3dmask_tool, AFNI',
                                          '3. make sure there is enough voxels (3dclust, AFNI)'], }
            json_object = json.dumps(dictionary, indent=2)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        else:
            nl = 'WARNING:' + str(input1) + ' not found!!! this may be because you have not provided an aseg file, ' + \
                 ' then no extraction of WM or Ventricles or GM will be possible... pls check that!'
            print(bcolors.WARNING + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

    ## in func space resample to func
    BRAIN = ants.image_read(opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz'))
    TRANS = ants.apply_transforms(fixed=MEAN, moving=BRAIN,
                                  transformlist=mvt_shft_INV_ANTs,
                                  interpolator=n_for_ANTS,
                                  whichtoinvert=w2inv_inv)
    ants.image_write(TRANS, opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI.nii.gz'), ri=False)
    dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz'),
                              opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz')],
                  "Description": 'Normalization (ANTspy).', },
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI.json'), "w") as outfile:
        outfile.write(json_object)


    #### apply to each available atlas
    if len(list_atlases) > 0:
        for atlas in list_atlases:
            ## in anat space resample to funcdo_anat_to_func
            if anat_func_same_space == True:
                mvt_shft = opj(dir_prepro, ID + '_brain_for_Align_Center_inv.1D')
                command = 'singularity run' + s_bind + afni_sif + '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + \
                          opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + ' ' + opj(labels_dir, TfMRI + opb(atlas)) + ' -overwrite'
                nl = spgo(command)
                diary.write(f'\n{nl}')
                print(nl)

                command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -final NN -overwrite -1Dmatrix_apply ' + mvt_shft + \
                ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + \
                ' -master ' + opj(dir_prepro, ID + '_mprage_reorient' + TfMRI + '.nii.gz') + \
                ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas))
                nl = spgo(command)
                diary.write(f'\n{nl}')
                print(nl)

                caca = resample_to_img(opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)),
                                       opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz'),  interpolation='nearest')
                caca.to_filename(opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)))

                dictionary = {"Sources": [opj(labels_dir, TfMRI + opb(atlas)),
                                          opj(dir_prepro, ID + '_mprage_reorient' + TfMRI + '.nii.gz'),
                                          opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')],
                              "Description": ['Normalization (3dAllineate, AFNI).',
                                              'Resampling (resample_to_img, nilearn)']},
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)[:-7] + '.json'), "w") as outfile:
                    outfile.write(json_object)



            elif IhaveanANAT == False:
                command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
                ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + \
                ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                ' -input  ' + atlas
                nl = spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
                dictionary = {"Sources": [atlas,
                                          opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz')],
                              "Description": 'Resampling (3dresample, AFNI)'},
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)[:-7] + '.json'), "w") as outfile:
                    outfile.write(json_object)

            else:
                command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
                ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + \
                ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                ' -input  ' + opj(labels_dir, TfMRI + opb(atlas))
                nl = spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
                dictionary = {"Sources": [opj(labels_dir, TfMRI + opb(atlas)),
                                          opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz')],
                              "Description": 'Resampling (3dresample, AFNI)'},
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)[:-7] + '.json'), "w") as outfile:
                    outfile.write(json_object)

            ## in func space resample to func
            ATLAS = ants.image_read(opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)))
            TRANS = ants.apply_transforms(fixed=MEAN, moving=ATLAS,
                                          transformlist=mvt_shft_INV_ANTs,
                                          interpolator='nearestNeighbor',
                                          whichtoinvert=w2inv_inv)
            ants.image_write(TRANS, opj(dir_fMRI_Refth_RS_prepro1, opb(atlas)), ri=False)

            dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)),
                                      opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz')],
                          "Description": 'Normalization (nearestNeighbor, AFNI)'},
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, opb(atlas)[:-7] + '.json'), "w") as outfile:
                outfile.write(json_object)

    else:
        nl = 'WARNING: list_atlases is empty!'
        print(bcolors.WARNING + 'WARNING: list_atlases is empty!' + bcolors.ENDC)
        diary.write(f'\n{nl}')


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
    nl = spgo(command)
    diary.write(f'\n{nl}')
    print(nl)
    dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz'),
                              anat_subject],
                  "Description": 'resampling (3dresample, AFNI)'},
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS_anat_resolution.json'), "w") as outfile:
        outfile.write(json_object)

    #### create a nice anat in func space
    if anat_func_same_space == True:

        mvt_shft = opj(dir_prepro, ID + '_brain_for_Align_Center_inv.1D')
        anatstd = opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_plot.nii.gz')

        ##### apply the recenter fmri
        command = 'singularity run' + s_bind + afni_sif + '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + anatstd + ' ' + anat_subject + ' -overwrite'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -overwrite -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + anatstd + \
        ' -input  ' + anatstd + \
        ' -master ' + opj(dir_prepro, ID + '_mprage_reorient' + TfMRI + '.nii.gz')
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        dictionary = {"Sources": [anat_subject,
                                  anatstd,
                                  opj(dir_prepro, ID + '_mprage_reorient' + TfMRI + '.nii.gz')],
                      "Description": 'Normalization (3dAllineate, AFNI).'},
        json_object = json.dumps(dictionary, indent=2)
        with open(anatstd[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

    else:
        anatstd = anat_subject
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + anatstd + \
        ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_plot.nii.gz') + ' -expr "a"'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": anatstd,
                      "Description": 'Copy.'},
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_plot.json'), "w") as outfile:
            outfile.write(json_object)


    ## in func space resample to func
    ANAT = ants.image_read(anatstd)
    MEAN_RES = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS_anat_resolution.nii.gz'))
    TRANS = ants.apply_transforms(fixed=MEAN_RES, moving=ANAT,
                                  transformlist=mvt_shft_INV_ANTs,
                                  interpolator='nearestNeighbor',
                                  whichtoinvert=w2inv_inv)
    ants.image_write(TRANS, opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI_anat_resolution.nii.gz'), ri=False)
    dictionary = {"Sources": [anatstd,
                              opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS_anat_resolution.nii.gz')],
                  "Description": 'Normalization (nearestNeighbore, ANTSpy).'},
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI_anat_resolution.json'), "w") as outfile:
        outfile.write(json_object)


    ### finally mask the func with mask
    for i in range(0, int(nb_run)):
        root_RS = extract_filename(RS[i])
        root_RS_ref = extract_filename(RS[REF_int])

        ### take opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
        
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' +   opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz') + \
                  ' -b ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref.nii.gz') + \
                  ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref_RcT_masked.nii.gz') + ' -expr "a*b"'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref.nii.gz'),
                                  opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz')],
                      "Description": 'Skull stripping (3dcalc, AFNI).'},
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref_RcT_masked.json'), "w") as outfile:
            outfile.write(json_object)


        if not root_RS == root_RS_ref:  # do not process ref...
            #### send mask the bold "not in norm", essentially for QC
            mvt_shft_ANTs_func_to_norm = [opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp_0GenericAffine.mat'),
                                          opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp_1InverseWarp.nii.gz')]
            w2inv_fwd_func_to_norm = [True, False]

            mask = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz'))
            MEAN = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_deob.nii.gz'))
            TRANS = ants.apply_transforms(fixed=MEAN, moving=mask,
                                          transformlist=mvt_shft_ANTs_func_to_norm,
                                          interpolator='nearestNeighbor',
                                          whichtoinvert=w2inv_fwd_func_to_norm)
            ants.image_write(TRANS, opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_mask_final_in_fMRI_orig.nii.gz'), ri=False)

            dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz'),
                                      opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_deob.nii.gz')],
                          "Description": 'Normalization (ANTspy).'},
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_mask_final_in_fMRI_orig.json'), "w") as outfile:
                outfile.write(json_object)

        else:
            command = 'singularity run' + s_bind + afni_sif + '3dcopy ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz') + \
                      ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_mask_final_in_fMRI_orig.nii.gz') + overwrite
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat.nii.gz'),
                          "Description": 'Copy.'},
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_mask_final_in_fMRI_orig.json'), "w") as outfile:
                outfile.write(json_object)

        diary.write(f'\n')

        bids_dir = opd(opd(opd(opd(opd(dir_fMRI_Refth_RS_prepro1)))))
        if not ope(opj(bids_dir + '/QC/')):
            os.mkdir(opj(bids_dir + '/QC/'))

        bids_dir = opd(opd(opd(opd(opd(dir_fMRI_Refth_RS_prepro1)))))
        if not ope(opj(bids_dir + '/QC/','mask_to_fMRI_orig')):
            os.mkdir(opj(bids_dir + '/QC/','mask_to_fMRI_orig'))

        ####plot the QC
        plot_qc(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_deob.nii.gz'), opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_mask_final_in_fMRI_orig.nii.gz'), opj(bids_dir, 'QC', 'mask_to_fMRI_orig', root_RS + '_mask_final_in_fMRI_orig.png'))

    diary.close()