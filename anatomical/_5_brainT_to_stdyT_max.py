import os
import subprocess
from nilearn import plotting
import ants

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

def brainT_to_T_max(dir_prepro, ID, Session, listTimage, n_for_ANTS, dir_transfo, type_norm, BASE_SS_coregistr, Ref_file, volumes_dir, transfo_concat_inv,w2inv_inv, creat_study_template, which_on, all_ID_max, max_session, all_data_path_max, all_ID, all_Session, all_data_path, study_template_atlas_forlder, otheranat, bids_dir, overwrite):
    #################################################################### coregistration with ANTs ####################################################################

    print('Coregistration ready baby')
    print('3')
    print('2')
    print('1')
    print('prouut')

    ######Coregistration!!!! (template to anat)

    IMG = ants.image_read(Ref_file)
    REF = ants.image_read(BASE_SS_coregistr)

    mTx = ants.registration(fixed=IMG, moving=REF,
                            type_of_transform='SyNRA',
                            outprefix=opj(dir_transfo, 'template_to_' + type_norm + '_SyN_final_max_'),
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
    TRANS = ants.apply_transforms(fixed=IMG, moving=REF,
                                  transformlist=mTx['fwdtransforms'], interpolator=n_for_ANTS)
    ants.image_write(TRANS, opj(dir_prepro, 'template_to_' + type_norm + '_SyN_final_max.nii.gz'), ri=False)


    for Timage in listTimage:

        img_ref = opj(volumes_dir,ID + Timage + '_brain.nii.gz')
        print(img_ref)
        print(Timage)

        IMG = ants.image_read(img_ref)
        TRANS = ants.apply_transforms(fixed=REF, moving=IMG,
                                      transformlist=transfo_concat_inv,
                                      interpolator=n_for_ANTS,whichtoinvert=w2inv_inv)
        ants.image_write(TRANS, opj(dir_prepro, 'template_to_' + Timage + '_to_template_test_matrix_max.nii.gz'), ri=False)


        ####plot the Reffile
        try:
            display = plotting.plot_anat(BASE_SS_coregistr, threshold='auto',
                                         display_mode='mosaic', dim=4)
            display.add_contours(opj(dir_prepro,'template_to_' + Timage + '_SyN_final_test_matrix_max.nii.gz'),
                                 linewidths=.2, colors=['red'])
            display.savefig(bids_dir + '/QC/' + ID + '_' + str(Session) + '_' + Timage + '_SyN_final_test_matrix_max.png')
            # Don't forget to close the display
            display.close()
        except:
            display = plotting.plot_anat(opj(dir_prepro,'template_to_' + Timage + '_SyN_final_test_matrix_max.nii.gz'), threshold='auto',
                                         display_mode='mosaic', dim=4)
            display.savefig(bids_dir + '/QC/' + ID + '_' + str(Session) + '_' + Timage + '_SyN_final_test_matrix_max.png')
            # Don't forget to close the display
            display.close()


    if creat_study_template == True:
    ###creat brain of the subject template
        # average the available T1 files
        if which_on == 'max': # all or max
            all_ID_temp = all_ID_max
            all_Session_temp = max_session
            all_data_path_temp =all_data_path_max

        elif which_on == 'all':
            all_ID_temp = all_ID
            all_Session_temp = all_Session
            all_data_path_temp = all_data_path

        second_template = []
        for data_path in all_data_path_temp:
            path_anat    = opj(data_path,'anat/')
            dir_native    = opj(path_anat,'native')
            dir_prepro    = opj(dir_native,'01_preprocess')
            second_template.append(opj(dir_prepro, Timage + '_to_template_test_matrix_max.nii.gz'))

        TEMP= ' '.join(second_template)
        if not ope(opj(study_template_atlas_forlder, 'studytemplate2_' + otheranat)): os.mkdir(opj(study_template_atlas_forlder, 'studytemplate2_' + otheranat))
        #### to change for 3dmean!!!!!!!
        #command = 'mri_robust_template  ' + TEMP + ' --template ' + opj(study_template_atlas_forlder, 'studytemplate2_' + otheranat, 'study_template.nii.gz') + ' --satit --inittp 1 --fixtp --noit --iscale --average 0'
        #spco([command], shell=True)
