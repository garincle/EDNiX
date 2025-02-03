import nilearn
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Set the non-interactive backend
from nilearn import image
from nilearn import plotting
import subprocess
import os
from nilearn.input_data import NiftiMasker
import nibabel as nib
from nilearn.input_data import NiftiLabelsMasker
from nilearn.image import resample_to_img
from fonctions.extract_filename import extract_filename
import scipy.ndimage as ndimage
import ants
import datetime
import json
import matplotlib.pyplot as plt

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

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

spco = subprocess.check_output
spgo = subprocess.getoutput

#################################################################################################
####Seed base analysis
################################################################################################# 
def SBA(SBAspace, BASE_SS_coregistr, erod_seed, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
    dir_fMRI_Refth_RS_prepro3, RS, nb_run, selected_atlases, panda_files, oversample_map, use_cortical_mask_func,
        cut_coordsX, cut_coordsY, cut_coordsZ, threshold_val, s_bind, afni_sif, diary_file):

    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  Working on step ' + str(11) + '(function: _11_Seed_base_many_regionsatlas).  ##'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    for space in SBAspace:
        for i in range(0, int(nb_run)):
            root_RS = extract_filename(RS[i])
            if space=='func':
                direction_results = dir_fMRI_Refth_RS_prepro1
                func_filename     = opj(direction_results, root_RS + '_residual.nii.gz')

                if oversample_map == True:
                    ### not possible yet
                    studytemplatebrain = opj(direction_results, 'Ref_anat_in_fMRI_anat_resolution.nii.gz')
                else:
                    studytemplatebrain = opj(direction_results, 'Ref_anat_in_fMRI.nii.gz')

                if use_cortical_mask_func == True:
                    cortical_mask_func = opj(direction_results,'Gmask.nii.gz')
                else:
                    cortical_mask_func = opj(direction_results,'mask_ref.nii.gz')

            elif space=='anat':
                direction_results = dir_fMRI_Refth_RS_prepro2
                func_filename     = opj(direction_results, root_RS + '_residual_in_anat.nii.gz')

                if oversample_map == True:
                    ## to test
                    studytemplatebrain = opj(direction_results,'orig_anat_for_plot.nii.gz')
                else:
                    studytemplatebrain = opj(direction_results,'anat_rsp_in_func.nii.gz')

                if use_cortical_mask_func == True:
                    cortical_mask_func = opj(direction_results,'Gmask.nii.gz')
                else:
                    cortical_mask_func = opj(direction_results,'mask_ref.nii.gz')

            elif space == 'atlas':
                direction_results = dir_fMRI_Refth_RS_prepro3
                func_filename     = opj(direction_results, root_RS + '_residual_in_template.nii.gz')

                if oversample_map == True:
                    studytemplatebrain = BASE_SS_coregistr
                else:
                    studytemplatebrain = opj(direction_results,'BASE_SS_fMRI.nii.gz')

                if use_cortical_mask_func == True:
                    cortical_mask_func = opj(direction_results,'Gmask.nii.gz')
                else:
                    cortical_mask_func = opj(direction_results,'mask_brain.nii.gz') # should change the name as mask_ref.nii.gz
            else:
                nl = 'WARNING: will not perform ' + str(direction_results) + ' space because SBAspace is ' + str(SBAspace)
                print(bcolors.WARNING + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')


            if ope(func_filename):
                # Get the global signal within the brain (or gray matter) mask
                resampled_cortical_mask_func = resample_to_img(cortical_mask_func, func_filename,
                                                               interpolation='nearest')
                resampled_cortical_mask_func.to_filename(opj(direction_results, 'cortical_mask_funcrsp.nii.gz'))

                extracted_data2 = nib.load(opj(direction_results, 'cortical_mask_funcrsp.nii.gz')).get_fdata()
                labeled_img2 = image.new_img_like(func_filename, extracted_data2, copy_header=True)
                labeled_img2.to_filename(opj(direction_results, 'cortical_mask_funcrsp.nii.gz'))
                dictionary = {"Sources": [cortical_mask_func,
                                          func_filename],
                              "Description": ' resampling (nilearn)', },
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(direction_results, 'cortical_mask_funcrsp.json'), "w") as outfile:
                    outfile.write(json_object)

                brain_masker = NiftiMasker(standardize="zscore", smoothing_fwhm=None,
                                           memory_level=0, verbose=1,
                                           mask_img=opj(direction_results, 'cortical_mask_funcrsp.nii.gz'))
                brain_time_series = brain_masker.fit_transform(func_filename)

                for panda_file, atlas in zip(panda_files, selected_atlases):

                    atlas_filename = opj(direction_results, atlas)
                    output_results = opj(direction_results, '10_Results', 'SBA')

                    if not os.path.exists(opj(direction_results, '10_Results')):
                        os.mkdir(opj(direction_results, '10_Results'))
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
                        Seed_label = row['label']

                        output_folder = opj(output_results, Seed_name)
                        if not os.path.exists(output_folder):
                            os.mkdir(output_folder)

                        atlas_path = opj(output_folder, Seed_name + '.nii.gz')

                        command = ('singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + atlas_filename  +
                                   ' -expr "ispositive(a)*(iszero(ispositive((a-' + str(Seed_label) + ')^2)))" -prefix '
                                   + atlas_path)
                        nl = spgo(command)
                        diary.write(f'\n{nl}')
                        print(nl)
                        dictionary = {"Sources": atlas_filename,
                                      "Description": 'ispositive(a)*(iszero(ispositive((a-' + str(Seed_label) + ')^2))) (3dcalc, AFNI).'},
                        json_object = json.dumps(dictionary, indent=2)
                        with open(atlas_path[:-7] +  '.json',"w") as outfile:
                            outfile.write(json_object)

                        # Load the atlas
                        atlas_img = ants.image_read(atlas_path)

                        # Get the unique region labels in the atlas
                        labels = np.unique(atlas_img.numpy())  # Assumes label 0 is background

                        # Check if all the values in the image are zero
                        if np.all(labels == 0):
                            nl = "WARNING: The NIfTI image is empty (all voxel values are zero)."
                            print(nl)
                            diary.write(f'\n{nl}')

                            diary_file_WARNING = opj(opd(opd(dir_fMRI_Refth_RS_prepro1)), 'SEED_' + str(Seed_name) + '_EMPTY_WARNING.txt')
                            if not opi(diary_file_WARNING):
                                diary_file_WARNING_file = open(diary_file_WARNING, "w")
                                diary_file_WARNING_file.write(f'\n{nl}')
                            else:
                                diary_file_WARNING_file = open(diary_file_WARNING, "a")
                                diary_file_WARNING_file.write(f'\n{nl}')
                            diary_file_WARNING_file.close()

                        else:
                            ## to prevent as much as possible the occurrence of spatial overlaps between the seeds
                            # that would create statistical dependencies between the voxels
                            # (because of interpolations during every normalization steps and the spatial smoothing)
                            # we added this additional "security" :

                            #####        Erosion of the seed, if possible        #####
                            # Erosion function
                            def erode_region(region, iterations=1):
                                binary_mask = region > 0  # Create a binary mask of the region
                                eroded_mask = ndimage.binary_erosion(binary_mask, iterations=iterations)
                                return eroded_mask.astype(np.uint8)

                            # Create an empty array to store the final results
                            final_result = np.zeros_like(atlas_img.numpy())

                            # Iterate over each label in the atlas, skipping background (0)

                            for label in labels:
                                if label == 0:
                                    continue  # Skip background
                                nl = f"Processing label: {label}"
                                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                                diary.write(f'\n{nl}')

                                # Isolate the current region (region is 1 where label matches, else 0)
                                region_mask = (atlas_img.numpy() == label).astype(np.uint8)
                                if erod_seed == True:

                                    # Try erosion by 1 voxel
                                    eroded_region_1 = erode_region(region_mask, iterations=1)
                                    if np.sum(eroded_region_1) > 10:  # Check if the region still has enough voxels
                                        nl = f"Label {label}: Eroded by 1 voxel, still has voxels."
                                        print(bcolors.OKGREEN + nl + bcolors.ENDC)
                                        diary.write(f'\n{nl}')
                                        final_result[eroded_region_1 > 0] = label
                                    else:
                                        nl = f"Label {label}: Eroding by 1 voxel removes too many voxels, keeping original region."
                                        print(bcolors.OKGREEN + nl + bcolors.ENDC)
                                        diary.write(f'\n{nl}')
                                        final_result[region_mask > 0] = label  # Keep the original mask
                                else:
                                    final_result[region_mask > 0] = label  # Keep the original mask

                            # Save the final eroded atlas
                            output_path = opj(output_folder, Seed_name + 'eroded.nii.gz')

                            final_atlas = ants.from_numpy(final_result, spacing=atlas_img.spacing,
                                                          origin=atlas_img.origin, direction=atlas_img.direction)
                            ants.image_write(final_atlas, output_path)

                            dictionary = {"Sources": atlas_filename,
                                          "Description": 'ispositive(a)*(iszero(ispositive((a-' + str(
                                              Seed_label) + ')^2))) (3dcalc, AFNI).'},
                            json_object = json.dumps(dictionary, indent=2)
                            with open(output_path[:-7] + '.json', "w") as outfile:
                                outfile.write(json_object)

                            nl = f"Eroded atlas saved to: {output_path}"
                            print(bcolors.OKGREEN + nl + bcolors.ENDC)
                            diary.write(f'\n{nl}')

                            # Get the global signal within each seed
                            labels_img = resample_to_img(output_path, func_filename, interpolation='nearest')
                            labels_img.to_filename(opj(output_folder, Seed_name + 'rsp.nii.gz'))

                            extracted_data2 = nib.load(opj(output_folder,Seed_name + 'rsp.nii.gz')).get_fdata()
                            labeled_img2    = image.new_img_like(func_filename,
                                extracted_data2, copy_header=True)
                            labeled_img2.to_filename(opj(output_folder,Seed_name + 'rsp.nii.gz'))
                            dictionary = {"Sources": [output_path,
                                                      func_filename],
                                          "Description": ' resampling (nilearn)',},
                            json_object = json.dumps(dictionary, indent=2)
                            with open(opj(output_folder,Seed_name + 'rsp.json'), "w") as outfile:
                                outfile.write(json_object)

                            seed_masker = NiftiLabelsMasker(labels_img=opj(output_folder,Seed_name + 'rsp.nii.gz'),
                                                            standardize='zscore', resampling_target= 'data', smoothing_fwhm=None,
                                                            memory_level=0, verbose=1)
                            seed_time_serie = seed_masker.fit_transform(func_filename)

                            ###########################################################################################
                            ##               Perform the seed-to-voxel correlation analysis                         ##
                            ###########################################################################################
                            # calculate the correlations
                            seed_to_voxel_correlations = (np.dot(brain_time_series.T, seed_time_serie) /
                                                          seed_time_serie.shape[0])
                            seed_to_voxel_correlations_img = brain_masker.inverse_transform(
                                seed_to_voxel_correlations.T)
                            seed_to_voxel_correlations_img.to_filename(opj(output_folder, root_RS + '_correlations.nii.gz'))
                            seed_to_voxel_correlations_img = nib.load(opj(output_folder, root_RS + '_correlations.nii.gz')).get_fdata()
                            seed_to_voxel_correlations_img = image.new_img_like(func_filename,seed_to_voxel_correlations_img, copy_header=True)
                            seed_to_voxel_correlations_img.to_filename(opj(output_folder, root_RS + '_correlations.nii.gz'))

                            dictionary = {"Sources": [opj(output_folder,Seed_name + 'rsp.nii.gz'),
                                                      opj(output_folder, Seed_name + 'cortical_mask_funcrsp.nii.gz')],
                                          "Description": ' Signal Correlation (nilearn)', },
                            json_object = json.dumps(dictionary, indent=2)
                            with open(opj(output_folder, root_RS + '_correlations.json'), "w") as outfile:
                                outfile.write(json_object)


                            # Fisher-z transformation

                            seed_to_voxel_correlations_fisher_z = np.arctanh(seed_to_voxel_correlations)
                            nl = "Seed-to-voxel correlation Fisher-z transformed: min = %.3f; max = %.3f" % (seed_to_voxel_correlations_fisher_z.min(),
                                                                                                             seed_to_voxel_correlations_fisher_z.max())
                            print(bcolors.OKGREEN + nl + bcolors.ENDC)
                            diary.write(f'\n{nl}')

                            seed_to_voxel_correlations_img_fish = brain_masker.inverse_transform(
                                seed_to_voxel_correlations_fisher_z.T)
                            seed_to_voxel_correlations_img_fish.to_filename(opj(output_folder, root_RS + '_correlations_fish.nii.gz'))
                            seed_to_voxel_correlations_img_fish = nib.load(opj(output_folder, root_RS + '_correlations_fish.nii.gz')).get_fdata()
                            seed_to_voxel_correlations_img_fish = image.new_img_like(opj(output_folder, root_RS + '_correlations.nii.gz'),seed_to_voxel_correlations_img_fish, copy_header=True)
                            seed_to_voxel_correlations_img_fish.to_filename(opj(output_folder, root_RS + '_correlations_fish.nii.gz'))

                            dictionary = {"Sources": opj(output_folder, root_RS + '_correlations.nii.gz'),
                                          "Description": ' Fisher transformation (numpy)', },
                            json_object = json.dumps(dictionary, indent=2)
                            with open(opj(output_folder, root_RS + '_correlations_fish.json'), "w") as outfile:
                                outfile.write(json_object)

                            correlation_img = ants.image_read(opj(output_folder, root_RS + '_correlations_fish.nii.gz'))
                            correlation_val = np.unique(correlation_img.numpy())  # Assumes label 0 is background

                            # Check if all the values in the image are zero
                            if np.all(correlation_val == 0):
                                nl = "WARNING: The NIfTI image is empty (all voxel values are zero)."
                                print(nl)
                                diary.write(f'\n{nl}')

                            else:

                                #### Remove an arbitrary percentage of the Zscore map
                                thresholded_map = opj(output_folder, 'higher_threshold.nii.gz')
                                threshold_val99 = 99
                                loadimg         = nib.load(opj(output_folder, root_RS + '_correlations_fish.nii.gz')).get_fdata()

                                loadimgsort99 =  np.percentile(np.abs(loadimg)[np.abs(loadimg)>0], threshold_val99)
                                custom_thresh =  np.percentile(np.abs(loadimg)[np.abs(loadimg)>0], threshold_val)

                                mask_imag = nilearn.image.threshold_img(opj(output_folder, root_RS + '_correlations_fish.nii.gz'), custom_thresh)
                                mask_imag.to_filename(thresholded_map)


                                labels_img = resample_to_img(thresholded_map,studytemplatebrain, interpolation='nearest')
                                labels_img.to_filename(thresholded_map)
                                extracted_data2 = nib.load(thresholded_map).get_fdata()
                                labeled_img2    = image.new_img_like(studytemplatebrain,extracted_data2, copy_header=True)
                                labeled_img2.to_filename(thresholded_map)

                                dictionary = {"Sources": opj(output_folder, root_RS + '_correlations_fish_custom_thresh.nii.gz'),
                                              "Description": ' remove the ' + str(custom_thresh) + ' percent', },
                                json_object = json.dumps(dictionary, indent=2)
                                with open(thresholded_map[:-7] + '.json', "w") as outfile:
                                    outfile.write(json_object)

                                # PLot the results
                                if direction_results == dir_fMRI_Refth_RS_prepro1:
                                    display = plotting.plot_stat_map(thresholded_map, threshold=custom_thresh, vmax=loadimgsort99,
                                        colorbar=True, bg_img=studytemplatebrain, display_mode='mosaic', cut_coords=(len(cut_coordsY), len(cut_coordsY), len(cut_coordsY)))
                                    display.savefig(opj(output_folder, root_RS + '_.jpg'))
                                    display.close()
                                    plt.close('all')

                                elif direction_results == dir_fMRI_Refth_RS_prepro2:
                                    display = plotting.plot_stat_map(thresholded_map, threshold=custom_thresh, vmax=loadimgsort99,
                                        colorbar=True, bg_img=studytemplatebrain, display_mode='mosaic', cut_coords=(len(cut_coordsY), len(cut_coordsY), len(cut_coordsY)))
                                    display.savefig(opj(output_folder, root_RS + '_.jpg'))
                                    display.close()
                                    plt.close('all')

                                elif direction_results == dir_fMRI_Refth_RS_prepro3:
                                    display = plotting.plot_stat_map(thresholded_map, threshold=custom_thresh, vmax=loadimgsort99,
                                        colorbar=True, bg_img=studytemplatebrain, display_mode='mosaic', cut_coords=(len(cut_coordsY), len(cut_coordsY), len(cut_coordsY)))
                                    display.savefig(opj(output_folder, root_RS + '_.jpg'))
                                    display.close()
                                    plt.close('all')
            else:
                nl = 'WARNING: ' + str(func_filename) + ' not found!!'
                print(bcolors.WARNING + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

    diary.write(f'\n')
    diary.close()

