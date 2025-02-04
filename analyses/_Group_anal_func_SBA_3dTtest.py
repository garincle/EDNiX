import nilearn
from nilearn import plotting
import glob
import subprocess
import os
import numpy as np
import nibabel as nib
from nilearn.masking import compute_epi_mask
import matplotlib.pyplot as plt

#################################################################################################
####LOADER YUNG LEMUR
#################################################################################################
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

#################################################################################################
####Seed base analysis
#################################################################################################
def LME_EDNiX(BASE_SS, oversample_map, mask_func, folder_atlases, cut_coordsX, cut_coordsY, panda_files, selected_atlases,
              cut_coordsZ,  cortical_mask_func, bids_dir, mean_imgs, lower_cutoff, upper_cutoff, s_bind, afni_sif, treshold_or_stat,
              seed_to_voxel_correlations_all_fish):

















    output_results1 = opj(bids_dir, 'Results')
    if not os.path.exists(output_results1): os.mkdir(output_results1)

    if oversample_map == True:
        studytemplatebrain = BASE_SS
    else:
        studytemplatebrain = opj(folder_atlases, 'BASE_SS_fMRI.nii.gz')

    mean_imgs_rs = nilearn.image.concat_imgs(mean_imgs, ensure_ndim=None, memory=None, memory_level=0,
                                             auto_resample=True, verbose=0)
    mask_img = compute_epi_mask(mean_imgs_rs,
                                lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff,
                                connected=True, opening=1,
                                exclude_zeros=True, ensure_finite=True)

    mask_img.to_filename(opj(output_results1, 'mask_mean_func.nii.gz'))

    command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + opj(output_results1, 'mask_mean_func.nii.gz') + \
              ' -input ' + opj(output_results1, 'mask_mean_func.nii.gz') + ' -fill_holes'
    nl = spgo(command)
    print(nl)

    command = 'singularity run' + s_bind + afni_sif + '3dresample -master ' + mask_func + \
              + ' -prefix ' + opj(output_results1, 'mask_mean_func.nii.gz') + ' -input ' + \
              opj(output_results1, 'mask_mean_func.nii.gz') + ' -overwrite -bound_type SLAB'
    nl = spgo(command)
    print(nl)

    command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + opj(output_results1, 'mask_mean_func.nii.gz') + \
              ' -b ' + mask_func + \
              ' -expr "a*b" -prefix ' + opj(output_results1, 'mask_mean_func_overlapp.nii.gz') + ' -overwrite'
    nl = spgo(command)
    print(nl)

    for panda_file, atlas in zip(panda_files, selected_atlases):

        output_results = opj(output_results1, 'Grp_SBA')
        if not os.path.exists(opj(output_results1, 'Grp_SBA')):
            os.mkdir(opj(output_results1, 'Grp_SBA'))
        if not os.path.exists(output_results):
            os.mkdir(output_results)

        ##########################################################################
        def format_seed_name(seed_name):
            # Define characters to be replaced with underscores
            replace_chars = [' ', '(', ')', ',', '/', ':', ';', '.', '-']

            # Replace each character in replace_chars with underscores
            for char in replace_chars:
                seed_name = seed_name.replace(char, '_')

            # Remove any other non-alphanumeric characters except for underscores
            formatted_name = ''.join(char for char in seed_name if char.isalnum() or char == '_')

            # Remove any leading or trailing underscores (optional)
            formatted_name = formatted_name.strip('_')
            return formatted_name

        for column, row in panda_file.T.items():
            Seed_name = row['region']
            Seed_name = format_seed_name(Seed_name)

            output_folder = opj(output_results, Seed_name)
            if not os.path.exists(output_folder):
                os.mkdir(output_folder)

            output_folder = opj(output_results1, 'Grp_SBA' + str(Seed_name))
            if not os.path.exists(output_folder): os.mkdir(output_folder)

            output_z1_list = []
            #if dolvl1==True:

            for all_dinv_statmap in seed_to_voxel_correlations_all_fish:
                #################stat test
                os.chdir(output_folder)
                if not os.path.exists(output_folder + 'ttest-stat_fisher.nii.gz'):

                    command = 'singularity run' + s_bind + afni_sif + '3dttest++ -setA ' + subprocess.list2cmdline(all_dinv_statmap) + '-toz -mask ' + \
                          opj(output_results1, 'mask_mean_func_overlapp.nii.gz') + ' -prefix ' + output_folder + 'ttest-stat_fisher.nii.gz'
                    nl = spgo(command)
                    print(nl)

                    command = 'singularity run' + s_bind + afni_sif + '1d_tool.py -infile ' + \
                                        output_folder + 'ttest-stat_fisher.CSimA.NN1_2sided.1D -csim_show_clustsize -verb 0 -csim_pthr ' + \
                                        tresh
                    nl = spgo(command)
                    print(nl)
                    csize = float(str(nl)[2:][:-3])
                    print(csize)

                z_map = output_folder + 'ttest-stat_fisher.nii.gz'
                z_map_data = nib.load(z_map).get_data()

                extracted_datacorect = []
                for n in list(range(0, 2)):
                    z_map_data_img = (z_map_data[:, :, :, :, n])
                    maxim = z_map_data_img.max()
                    labeled_img = nilearn.image.new_img_like(z_map, z_map_data_img, copy_header=True)
                    labeled_img.to_filename(output_folder + 'Seed' + str(n) + str(Seed_name) + '_all_chimp.nii.gz')

                output_z1 = output_folder + 'Seed' + str(n) + str(Seed_name) + '.nii.gz'
                output_z1_list.append(output_z1)
            else:
                output_z1_list = seed_to_voxel_correlations_all_fish

            #################stat test lvl 2
            if not os.path.exists(output_folder + 'ttest-stat_fisher.nii.gz'):
                os.chdir
                spco(['3dttest++', '-setA', subprocess.list2cmdline(output_z1_list), '-toz', '-mask', cortical_mask_func,
                '-prefix', output_folder + 'ttest-stat_fisher.nii.gz'])
                csize = float(str(spco(['1d_tool.py', '-infile', output_folder + 'ttest-stat_fisher.CSimA.NN1_2sided.1D', '-csim_show_clustsize', '-verb', '0', '-csim_pthr', tresh]))[2:][:-3])
                print(csize)

            z_map = output_folder + 'ttest-stat_fisher.nii.gz'
            z_map_data = nib.load(z_map).get_data()

            extracted_datacorect = []
            for n in list(range(0, 2)):
                z_map_data_img = (z_map_data[:, :, :, :, n])
                maxim = z_map_data_img.max()
                labeled_img = nilearn.image.new_img_like(z_map, z_map_data_img, copy_header=True)
                labeled_img.to_filename(output_folder + 'Seed' + str(n) + '_all_macaque.nii.gz')

            output_z = output_folder + 'Seed' + str(n) + '_all_macaque.nii.gz'

            output_z1_list.append(output_z)

            ####remove a percentage of the zmap
            threshold_val = 99
            loadimg = nib.load(output_z).get_data()
            loadimgsort99 =  np.percentile(np.abs(loadimg)[np.abs(loadimg)>0], threshold_val)

            if treshold_or_stat == 'stat':
                loadimgsort = '1.65'
            elif treshold_or_stat == int:
                loadimg = nib.load(output_z).get_data()
                loadimgsort =  np.percentile(np.abs(loadimg)[np.abs(loadimg)>0], treshold_or_stat)
            else:
                print('ERROR treshold_or_stat must be an int or the string "stat"')

            if size_or_stat == 'stat':
                csize_final = csize
            elif treshold_or_stat == int:
                csize_final = size_or_stat
            else:
                print('ERROR treshold_or_stat must be an int or the string "stat"')


            mask_imag = nilearn.image.threshold_img(output_z, loadimgsort, cluster_threshold=csize_final)
            mask_imag.to_filename(output_folder + 'treshold_or_stat.nii.gz')
            thresholded_map1 = output_folder + 'treshold_or_stat.nii.gz'

            display = plotting.plot_stat_map(thresholded_map1, dim=0, threshold=loadimgsort, vmax=loadimgsort99,
                                             colorbar=True, bg_img=studytemplatebrain,
                                             display_mode='mosaic',
                                             cut_coords=(len(cut_coordsY), len(cut_coordsY), len(cut_coordsY)))
            display.savefig(output_folder + '/treshold_or_stat_y_mosaic.jpg')
            display.close()
            plt.close('all')






