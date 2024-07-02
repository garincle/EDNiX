import os
import subprocess
import glob
import shutil
import sys
import nilearn
from nilearn import plotting


#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

def brainT_to_T(dir_prepro, ID, Session, listTimage, n_for_ANTS, dir_transfo, type_norm, BASE_SS_coregistr, Ref_file, volumes_dir, transfo_concat_inv, creat_sutdy_template, which_on, all_ID_max, max_session, all_data_path_max, all_ID, all_Session, all_data_path, study_template_atlas_forlder, otheranat, bids_dir, overwrite):

    #################################################################### coregistration with ANTs ####################################################################

    print('Coregistration ready baby')
    print('3')
    print('2')
    print('1')
    print('prouut')

    ######Coregistration!!!! (template to anat)
    command = 'antsRegistration -d 3 --float 0 --verbose 1 -n ' + n_for_ANTS +\
        ' -o [' + opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_') + ',' + opj(dir_prepro,'template_to_' + type_norm + '_SyN_final.nii.gz') + ']' + \
        ' -t Affine[0.1] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
        ' -m MI[' + Ref_file + ',' + BASE_SS_coregistr + ',1,32,Regular,0.2]' + \
        ' -t Syn[0.1,3,0] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
        ' -m CC[' + Ref_file + ',' + BASE_SS_coregistr + ',1,4,Regular,0.2]'
    spco([command], shell=True)


    for Timage in listTimage:

        img_ref = opj(volumes_dir,ID + Timage + '_brain.nii.gz')
        print(img_ref)
        print(Timage)

        command = 'antsApplyTransforms -d 3 -i ' + img_ref + \
        ' -r ' + BASE_SS_coregistr + \
        ' -o ' + opj(dir_prepro,'template_to_' + Timage + '_SyN_final_test_matrix.nii.gz') + \
        transfo_concat_inv + \
        ' -n ' + n_for_ANTS
        spco([command], shell=True)

        ####plot the Reffile
        try:
            display = plotting.plot_anat(BASE_SS_coregistr, threshold='auto',
                                         display_mode='mosaic', dim=4)
            display.add_contours(opj(dir_prepro,'template_to_' + Timage + '_SyN_final_test_matrix.nii.gz'),
                                 linewidths=.2, colors=['red'])
            display.savefig(bids_dir + '/QC/' + ID + '_' + str(Session) + '_' + Timage + '_SyN_final_test_matrix.png')
            # Don't forget to close the display
            display.close()
        except:
            display = plotting.plot_anat(opj(dir_prepro,'template_to_' + Timage + '_SyN_final_test_matrix.nii.gz'), threshold='auto',
                                         display_mode='mosaic', dim=4)
            display.savefig(bids_dir + '/QC/' + ID + '_' + str(Session) + '_' + Timage + '_SyN_final_test_matrix.png')
            # Don't forget to close the display
            display.close()