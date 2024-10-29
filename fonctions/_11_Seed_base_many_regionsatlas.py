import nilearn
import numpy as np
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
spco = subprocess.check_output
opd = os.path.dirname
ope = os.path.exists
#################################################################################################
####Seed base analysis
################################################################################################# 
def SBA(volumes_dir, BASE_SS_coregistr, TfMRI, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, 
    dir_fMRI_Refth_RS_prepro3, RS, nb_run, ID, selected_atlases, panda_files, oversample_map, use_cortical_mask_func, cut_coordsX, cut_coordsY, cut_coordsZ, threshold_val, overwrite,s_bind,afni_sif):

    for panda_file, atlas in zip(panda_files, selected_atlases):
        for i in range(0, int(nb_run)):
            root_RS = extract_filename(RS[i])
            for direction_results in [dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3]:
                if direction_results == dir_fMRI_Refth_RS_prepro1:
                    if oversample_map == True:
                        ###not possible yet
                        studytemplatebrain = opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI_anat_resolution.nii.gz')
                    else:
                        studytemplatebrain = opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI.nii.gz')
                    func_filename = opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz')
                    if use_cortical_mask_func == True:
                        cortical_mask_func = opj(dir_fMRI_Refth_RS_prepro1,'Gmask.nii.gz')
                    else:
                        cortical_mask_func = opj(dir_fMRI_Refth_RS_prepro1,'mask_ref.nii.gz')

                    atlas_filename = opj(dir_fMRI_Refth_RS_prepro1, atlas)
                    output_results = opj(dir_fMRI_Refth_RS_prepro1, '10_Results')
                    if not os.path.exists(output_results): os.mkdir(output_results)
                    output_results = opj(dir_fMRI_Refth_RS_prepro1, '10_Results/SBA')
                    if not os.path.exists(output_results): os.mkdir(output_results)

                if direction_results == dir_fMRI_Refth_RS_prepro2:
                    if oversample_map == True:
                        ###to test
                        studytemplatebrain = opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_plot.nii.gz')
                    else:
                        studytemplatebrain = opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz')

                    func_filename = opj(dir_fMRI_Refth_RS_prepro2, root_RS + '_residual_in_anat.nii.gz')
                    if use_cortical_mask_func == True:
                        cortical_mask_func = opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz')
                    else:
                        cortical_mask_func = opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz')

                    atlas_filename = opj(dir_fMRI_Refth_RS_prepro2, atlas)
                    output_results = opj(dir_fMRI_Refth_RS_prepro2, '10_Results')
                    if not os.path.exists(output_results): os.mkdir(output_results)
                    output_results = opj(dir_fMRI_Refth_RS_prepro2, '10_Results/SBA')
                    if not os.path.exists(output_results): os.mkdir(output_results)


                if direction_results == dir_fMRI_Refth_RS_prepro3:
                    if oversample_map == True:
                        studytemplatebrain = BASE_SS_coregistr
                    else:
                        studytemplatebrain = opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz')

                    func_filename = opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz')
                    if use_cortical_mask_func == True:
                        cortical_mask_func = opj(dir_fMRI_Refth_RS_prepro3,'Gmask.nii.gz')
                    else:
                        cortical_mask_func = opj(dir_fMRI_Refth_RS_prepro3,'mask_brain.nii.gz')
                    atlas_filename = opj(dir_fMRI_Refth_RS_prepro3, atlas)
                    output_results = opj(dir_fMRI_Refth_RS_prepro3, '10_Results')
                    if not os.path.exists(output_results): os.mkdir(output_results)
                    output_results = opj(dir_fMRI_Refth_RS_prepro3, '10_Results/SBA')
                    if not os.path.exists(output_results): os.mkdir(output_results)

                if ope(func_filename):
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

                    for colomn, row in panda_file.T.items():
                        Seed_name = row['region']
                        Seed_name = format_seed_name(Seed_name)
                        Seed_label = row['label']

                        output_folder = opj(output_results, Seed_name + '/')
                        if not os.path.exists(output_folder): os.mkdir(output_folder)

                        command = ('singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + atlas_filename  +
                                   ' -expr "ispositive(a)*(iszero(ispositive((a-' + str(Seed_label) + ')^2)))" -prefix '
                                   + output_folder + '/' + Seed_name + '.nii.gz')
                        spco([command], shell=True)

                        # Load the NIfTI image
                        nifti_image = nib.load(output_folder + '/' + Seed_name + '.nii.gz')
                        # Get the image data as a NumPy array
                        image_data = nifti_image.get_fdata()
                        # Check if all the values in the image are zero
                        if np.all(image_data == 0):
                            print(bcolors.WARNING + "WARNING: The NIfTI image is empty (all voxel values are zero)." + bcolors.ENDC)

                        else:
                            ################## EROD the seed if possible ##################
                            # Load your atlas image
                            atlas_path = output_folder + '/' + Seed_name + '.nii.gz'
                            atlas_img = ants.image_read(atlas_path)
                            # Get the unique region labels in the atlas
                            labels = np.unique(atlas_img.numpy())  # Assumes label 0 is background

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
                                print(f"Processing label: {label}")
                                # Isolate the current region (region is 1 where label matches, else 0)
                                region_mask = (atlas_img.numpy() == label).astype(np.uint8)
                                # Try erosion by 1 voxel
                                eroded_region_1 = erode_region(region_mask, iterations=1)
                                if np.sum(eroded_region_1) > 10:  # Check if the region still has voxels
                                    print(f"Label {label}: Eroded by 1 voxel, still has voxels.")
                                    final_result[eroded_region_1 > 0] = label
                                else:
                                    print(f"Label {label}: Eroding by 1 voxel removes all voxels, keeping original region.")
                                    final_result[region_mask > 0] = label  # Keep the original mask
                            # Save the final eroded atlas
                            print(final_result)
                            final_atlas = ants.from_numpy(final_result, spacing=atlas_img.spacing, origin=atlas_img.origin, direction=atlas_img.direction)
                            output_path = output_folder + '/' + Seed_name + 'eroded.nii.gz'
                            ants.image_write(final_atlas, output_path)
                            print(f"Eroded atlas saved to: {output_path}")

                            labels_img = resample_to_img(output_path, func_filename, interpolation='nearest')
                            labels_img.to_filename(output_folder + '/' + Seed_name + 'rsp.nii.gz')
                            # Plot the generated mask using the mask_img_ attribute
                            extracted_data2 = nib.load(output_folder + '/' + Seed_name + 'rsp.nii.gz').get_fdata()
                            labeled_img2 = image.new_img_like(func_filename,
                                extracted_data2, copy_header=True)
                            labeled_img2.to_filename(output_folder + '/' + Seed_name + 'rsp.nii.gz')
                            seed_masker = NiftiLabelsMasker(labels_img=output_folder + '/' + Seed_name + 'rsp.nii.gz', standardize='zscore', resampling_target= 'data', smoothing_fwhm=None,
                                memory_level=1, verbose=1)

                            ##########################################################################
                            seed_time_serie = seed_masker.fit_transform(func_filename)
                            resampled_cortical_mask_func = resample_to_img(cortical_mask_func, func_filename, interpolation='nearest')
                            resampled_cortical_mask_func.to_filename(output_folder + '/' + Seed_name + 'cortical_mask_funcrsp.nii.gz')
                            # Plot the generated mask using the mask_img_ attribute
                            extracted_data2 = nib.load(output_folder + '/' + Seed_name + 'cortical_mask_funcrsp.nii.gz').get_fdata()
                            labeled_img2 = image.new_img_like(func_filename,
                                extracted_data2, copy_header=True)
                            labeled_img2.to_filename(output_folder + '/' + Seed_name + 'cortical_mask_funcrsp.nii.gz')

                            ##########################################################################
                            brain_masker = NiftiMasker(standardize='zscore', smoothing_fwhm=None,
                                memory_level=1, verbose=1, mask_img=output_folder + '/' + Seed_name + 'cortical_mask_funcrsp.nii.gz')

                            ##########################################################################

                            brain_time_series = brain_masker.fit_transform(func_filename)

                            ##########################################################################
                            # Performing the seed-to-voxel correlation analysis
                            seed_to_voxel_correlations = (np.dot(brain_time_series.T, seed_time_serie) /
                                                          seed_time_serie.shape[0])

                            ################################################

                            seed_to_voxel_correlations_img = brain_masker.inverse_transform(
                                seed_to_voxel_correlations.T)
                            seed_to_voxel_correlations_img.to_filename(output_folder + '/' + root_RS + '_correlations.nii.gz')

                            ##########################################################################
                            # Fisher-z transformation and save nifti
                            seed_to_voxel_correlations_fisher_z = np.arctanh(seed_to_voxel_correlations)
                            print(bcolors.OKGREEN + "Seed-to-voxel correlation Fisher-z transformed: min = %.3f; max = %.3f"
                                  % (seed_to_voxel_correlations_fisher_z.min(),
                                     seed_to_voxel_correlations_fisher_z.max()) + bcolors.ENDC)

                            seed_to_voxel_correlations_img_fish = brain_masker.inverse_transform(
                                seed_to_voxel_correlations_fisher_z.T)
                            seed_to_voxel_correlations_img_fish.to_filename(output_folder + '/' + root_RS + '_correlations_fish.nii.gz')

                            ##########################################################################
                            # Fisher-z transformation and save nifti
                            ####remove a percentage of the zmap
                            threshold_val99 = 99
                            loadimg = nib.load(output_folder + '/' + root_RS + '_correlations_fish.nii.gz').get_fdata()
                            loadimgsort99 =  np.percentile(np.abs(loadimg)[np.abs(loadimg)>0], threshold_val99)

                            loadimg = nib.load(output_folder + '/' + root_RS + '_correlations_fish.nii.gz').get_fdata()
                            custom_thresh =  np.percentile(np.abs(loadimg)[np.abs(loadimg)>0], threshold_val)

                            mask_imag = nilearn.image.threshold_img(output_folder + '/' + root_RS + '_correlations.nii.gz', custom_thresh)
                            mask_imag.to_filename(output_folder + 'higher_thresold.nii.gz')

                            labels_img = resample_to_img(output_folder + 'higher_thresold.nii.gz', studytemplatebrain, interpolation='nearest')
                            labels_img.to_filename(output_folder + 'higher_thresold_res.nii.gz')
                            extracted_data2 = nib.load(output_folder + 'higher_thresold_res.nii.gz').get_fdata()
                            labeled_img2 = image.new_img_like(studytemplatebrain,
                                extracted_data2, copy_header=True)
                            labeled_img2.to_filename( output_folder + 'higher_thresold_res.nii.gz')

                            thresholded_map1 = output_folder + 'higher_thresold_res.nii.gz'

                            if direction_results == dir_fMRI_Refth_RS_prepro1:
                                display = plotting.plot_stat_map(thresholded_map1, threshold=custom_thresh, vmax=loadimgsort99,
                                    colorbar=True, bg_img=studytemplatebrain, display_mode='mosaic', cut_coords=(len(cut_coordsY), len(cut_coordsY), len(cut_coordsY)))
                                display.savefig(output_folder + '/' + root_RS + '_.jpg')
                                display.close()

                            elif direction_results == dir_fMRI_Refth_RS_prepro2:
                                display = plotting.plot_stat_map(thresholded_map1, threshold=custom_thresh, vmax=loadimgsort99,
                                    colorbar=True, bg_img=studytemplatebrain, display_mode='mosaic', cut_coords=(len(cut_coordsY), len(cut_coordsY), len(cut_coordsY)))
                                display.savefig(output_folder + '/' + root_RS + '_.jpg')
                                display.close()

                            elif direction_results == dir_fMRI_Refth_RS_prepro3:
                                display = plotting.plot_stat_map(thresholded_map1, threshold=custom_thresh, vmax=loadimgsort99,
                                    colorbar=True, bg_img=studytemplatebrain, display_mode='x', cut_coords=cut_coordsX)
                                display.savefig(output_folder + '/' + root_RS + '_x_.jpg')
                                display.close()

                                display = plotting.plot_stat_map(thresholded_map1, threshold=custom_thresh, vmax=loadimgsort99,
                                    colorbar=True, bg_img=studytemplatebrain, display_mode='y', cut_coords=cut_coordsY)
                                display.savefig(output_folder + '/' + root_RS + '_y_.jpg')
                                display.close()

                                display = plotting.plot_stat_map(thresholded_map1, threshold=custom_thresh, vmax=loadimgsort99,
                                    colorbar=True, bg_img=studytemplatebrain, display_mode='z', cut_coords=cut_coordsZ)
                                display.savefig(output_folder + '/' + root_RS + '_z_.jpg')
                                display.close()

                                display = plotting.plot_stat_map(thresholded_map1, threshold=custom_thresh, vmax=loadimgsort99,
                                    colorbar=True, bg_img=studytemplatebrain, display_mode='mosaic', cut_coords=(len(cut_coordsY), len(cut_coordsY), len(cut_coordsY)))
                                display.savefig(output_folder + '/' + root_RS + '_.jpg')
                                display.close()
                else:
                    print(bcolors.WARNING + 'WARNING: ' + str(func_filename) + ' not found!!' + bcolors.ENDC)

