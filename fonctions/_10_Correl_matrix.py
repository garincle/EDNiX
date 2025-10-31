import os
import pandas as pd
import nilearn
from nilearn import plotting
from nilearn.input_data import NiftiLabelsMasker
from nilearn.connectome import ConnectivityMeasure
import matplotlib.pyplot as plt
import numpy as np
import nibabel as nib
import scipy.ndimage as ndimage
import ants
import json

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists

from Tools import run_cmd
from fonctions.extract_filename import extract_filename


#################################################################################################
####Seed base analysis
#################################################################################################
def correl_matrix(dir_prepro_orig_postprocessed, RS, nb_run, selected_atlases_matrix, segmentation_name_list,
                  ID, Session, TR_val, dir_prepro_orig_labels, dir_prepro_orig,
                  sing_afni,diary_file):

    nl = '##  Working on step ' + str(10) + '(function: _10_Correl_matrix).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    for panda_file, atlas in zip(segmentation_name_list, selected_atlases_matrix):
        for i in range(int(nb_run)):

            runname = '_'.join([atlas[0],atlas[1], 'run', str(i)])

            root_RS = extract_filename(RS[i])
            func_filename = opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_residual.nii.gz')

            if ope(func_filename):
                atlas_filename = opj(dir_prepro_orig_labels, ID + '_seg-' + atlas[0] + '_dseg.nii.gz')

                output_results = opj(dir_prepro_orig, 'Stats')
                if not ope(output_results):
                    os.mkdir(output_results)
                output_results = opj(dir_prepro_orig, 'Stats', 'Correl_matrix')
                if not ope(output_results):
                    os.mkdir(output_results)

                nl = 'working on matrix for ' + str(panda_file['region'])
                run_cmd.msg(nl, diary_file, 'OKGREEN')

                # Load the NIfTI images
                img1 = nib.load(atlas_filename)
                img2 = nib.load(func_filename)

                # Extract headers
                header1 = img1.header
                header2 = img2.header

                # Compare specific fields
                fields_to_compare = ['dim', 'pixdim']
                tolerance = 0.00006
                differences = {}
                for field in fields_to_compare:
                    value1 = header1[field][:4]
                    value2 = header2[field][:4]
                    if not np.allclose(value1, value2, atol=tolerance):
                        differences[field] = (value1, value2)

                # Print differences
                if differences:
                    nl = "INFO: Differences between images:"
                    run_cmd.msg(nl, diary_file, 'OKGREEN')

                    for field, values in differences.items():

                        nl = f"{field}:"
                        run_cmd.msg(nl, diary_file, 'OKGREEN')
                        nl = f"  Image 1: {values[0]}"
                        run_cmd.msg(nl, diary_file, 'OKGREEN')
                        nl = f"  Image 2: {values[1]}"
                        run_cmd.msg(nl, diary_file, 'OKGREEN')

                        dummy = nilearn.image.resample_to_img(atlas_filename, func_filename, interpolation='nearest')
                        dummy.to_filename(atlas_filename)

                        extracted_data = nib.load(atlas_filename).get_fdata()
                        labeled_img2 = nilearn.image.new_img_like(func_filename, extracted_data, copy_header=True)
                        labeled_img2.to_filename(atlas_filename)

                        dictionary = {"Sources": [atlas_filename,
                                                  func_filename],
                                      "Description": ' Resampling (resample_to_img, nilearn).'},
                        json_object = json.dumps(dictionary, indent=2)
                        with open(atlas_filename[:-7] + '.json', "w") as outfile:
                            outfile.write(json_object)

                else:
                    nl = "INFO: No differences found in the specified fields."
                    run_cmd.msg(nl, diary_file, 'OKGREEN')

                # Load your brain atlas data
                atlas_img = nib.load(atlas_filename)
                atlas_data = atlas_img.get_fdata()[:,:,:,atlas[1]]

                # Flatten the atlas data for easy comparison
                atlas_flat = np.unique(atlas_data.flatten().astype(int))
                atlas_flat[atlas_flat != 0]

                missing_labels = set(panda_file['label']) - set(atlas_flat)
                if missing_labels:
                    run_cmd.msg(nl, diary_file, 'WARNING')

                # Step 1: Filter the labels
                filtered_labels = panda_file[panda_file['label'].isin(atlas_flat)]
                # Sort panda_file to match atlas_flat order
                filtered_labels = filtered_labels.sort_values('label')
                # Step 2: Create a dictionary for quick lookup of regions by label
                label_to_region = filtered_labels.set_index('label')['region'].to_dict()
                # Step 3: Order the regions according to atlas_flat
                atlas_filtered_list = [label_to_region[label] for label in atlas_flat if label in label_to_region]

                nl =  "INFO: label that you provided"
                run_cmd.msg(nl, diary_file, 'OKGREEN')
                nl =  str(panda_file)
                run_cmd.msg(nl, diary_file, 'OKGREEN')
                nl =  "INFO: label found in the atlas after extraction of the signal :" + str(atlas)
                run_cmd.msg(nl, diary_file, 'OKGREEN')
                nl =   str(atlas_flat)
                run_cmd.msg(nl, diary_file, 'OKGREEN')
                nl =   "INFO: if we combine them (analyse similarity we ware able to keep:"
                run_cmd.msg(nl, diary_file, 'OKGREEN')
                nl =   str(filtered_labels)
                run_cmd.msg(nl, diary_file, 'OKGREEN')
                nl =   "INFO: corresponding to the labels"
                run_cmd.msg(nl, diary_file, 'OKGREEN')
                nl =   str(atlas_filtered_list)
                run_cmd.msg(nl, diary_file, 'OKGREEN')
                nl =   "INFO: if they are different from the one you provided, you should investigate why!: region too small " + \
                       "for fMRI? error in the labeling? error in the coregistration?, this is not always " + \
                       "a big dill but you should know what happend!!"
                run_cmd.msg(nl, diary_file, 'OKGREEN')

                ### create a new list of label with zero if not in the new label list (for 3dcalc)
                filtered_labels_list = list(filtered_labels['label'])
                list_old_label = []
                for lab in atlas_flat:
                    if lab in np.array(filtered_labels_list):
                        list_old_label.append(lab)
                    else:
                        list_old_label.append(0)

                string_build_atlas = str('')

                for new_label, old_label in zip(list_old_label, list(atlas_flat)):
                    string_build_atlas = string_build_atlas + '(' + str(int(new_label)) + \
                                         '*(equals(a,' + str(int(old_label)) + ')))+'
                string_build_atlas2 = "'" + string_build_atlas[:-1] + "'"

                command = (sing_afni + '3dcalc' + ' -a ' + atlas_filename + '[' + str(atlas[1]) + '] -expr ' + string_build_atlas2 +
                           ' -prefix ' + opj(output_results, '_'.join([runname,'filtered.nii.gz'])) +
                           ' -overwrite')
                run_cmd.do(command, diary_file)

                dictionary = {"Sources": [atlas_filename,
                                          func_filename],
                              "Description": string_build_atlas2 + ' (3dcalc, AFNI).'},
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(output_results, '_'.join([runname,'filtered.json'])), "w") as outfile:
                    outfile.write(json_object)

                ################## EROD the seed if possible ##################
                # Load your atlas image
                atlas_path = opj(output_results, '_'.join([runname,'filtered.nii.gz']))
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
                labels_final = []
                for label in labels:
                    if label == 0:
                        continue  # Skip background

                    nl = f"Processing label: {label}"
                    run_cmd.msg(nl, diary_file, 'ENDC')
                    labels_final.append(label)

                    # Isolate the current region (region is 1 where label matches, else 0)
                    region_mask = (atlas_img.numpy() == label).astype(np.uint8)
                    # Try erosion by 1 voxel
                    eroded_region_1 = erode_region(region_mask, iterations=1)
                    if np.sum(eroded_region_1) > 0:  # Check if the region still has voxels
                        nl = f"Label {label}: Eroded by 1 voxel, still has voxels."
                        run_cmd.msg(nl, diary_file, 'ENDC')
                        final_result[eroded_region_1 > 0] = label
                    else:
                        nl = f"Label {label}: Eroding by 1 voxel removes all voxels, keeping original region."
                        run_cmd.msg(nl, diary_file, 'ENDC')
                        final_result[region_mask > 0] = label  # Keep the original mask
                    # Try erosion by 2 voxels
                    eroded_region_2 = erode_region(region_mask, iterations=2)

                    if np.sum(eroded_region_2) > 0:  # Check if the region still has voxels
                        nl = f"Label {label}: Eroded by 2 voxels, still has voxels."
                        run_cmd.msg(nl, diary_file, 'ENDC')
                        final_result[eroded_region_2 > 0] = label
                    else:
                        nl = f"Label {label}: Eroding by 2 voxels removes all voxels, keeping 1 voxel erosion or original region."
                        run_cmd.msg(nl, diary_file, 'ENDC')

                        # No need to change, we keep the previously applied erosion or original region.
                # Save the final eroded atlas
                nl = final_result
                run_cmd.msg(nl, diary_file, 'ENDC')

                final_atlas = ants.from_numpy(final_result, spacing=atlas_img.spacing, origin=atlas_img.origin,
                                              direction=atlas_img.direction)
                output_path = opj(output_results, '_'.join([runname,'filtered','eroded.nii.gz']))
                ants.image_write(final_atlas, output_path)
                dictionary = {"Sources": atlas_filename,
                              "Description": 'Erosion (numpy and ANTspy).'},
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(output_results, '_'.join([runname,'filtered','eroded.json'])), "w") as outfile:
                    outfile.write(json_object)

                nl = f"Eroded atlas saved to: {output_path}"
                run_cmd.msg(nl, diary_file, 'ENDC')

                ##########################################################################
                NAD_masker = NiftiLabelsMasker(labels_img=opj(output_results, '_'.join([runname,'filtered','eroded.nii.gz'])),
                                                detrend=False,
                                                smoothing_fwhm=None,
                                                low_pass=None,
                                                high_pass=None,
                                                t_r=TR_val,
                                                standardize='zscore_sample',
                                                memory=None, verbose=5)
                try:
                    time_series = NAD_masker.fit_transform(func_filename)
                except:
                    dummy = nilearn.image.resample_to_img(opj(output_results, atlas + '_filtered.nii.gz'), func_filename, interpolation='nearest')
                    dummy.to_filename(opj(output_results, atlas + '_filtered.nii.gz'))
                    extracted_data = nib.load(opj(output_results, atlas + '_filtered.nii.gz')).get_fdata()
                    labeled_img2 = nilearn.image.new_img_like(func_filename, extracted_data, copy_header=True)
                    labeled_img2.to_filename(opj(output_results, atlas + '_filtered.nii.gz'))

                    NAD_masker = NiftiLabelsMasker(
                        labels_img=opj(output_results, '_'.join([runname,'filtered.nii.gz'])), # must be a mistake ......
                        detrend=False,
                        smoothing_fwhm=None,
                        low_pass=None,
                        high_pass=None,
                        t_r=TR_val,
                        standardize= 'zscore_sample',
                        memory=None, verbose=5)

                    time_series = NAD_masker.fit_transform(func_filename)

                correlation_measure = ConnectivityMeasure(kind="correlation", standardize='zscore_sample',)
                correlation_matrix = correlation_measure.fit_transform([time_series])[0]

                # Plot the correlation matrix
                # Mask the main diagonal for visualization:
                np.fill_diagonal(correlation_matrix, 0)
                # The labels we have start with the background (0), hence we skip the
                # first label
                # matrices are ordered for block-like representation
                nl = str(correlation_matrix.shape)
                run_cmd.msg(nl, diary_file, 'OKGREEN')
                nl = str(len(panda_file['region']))
                run_cmd.msg(nl, diary_file, 'OKGREEN')
                nl = str(correlation_matrix)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

                plotting.plot_matrix(
                    correlation_matrix,
                    labels=atlas_filtered_list,
                    vmax=0.8,
                    vmin=-0.8)
                plt.savefig(opj(output_results,'_'.join([runname,'matrix.png'])))

                # Convert correlation matrix to a DataFrame
                corr_df = pd.DataFrame(correlation_matrix, index=atlas_filtered_list, columns=atlas_filtered_list)
                corr_df.to_csv(opj(output_results,'_'.join([runname,'matrix.csv'])))

                # Correctly construct new_df from lists
                new_df = pd.DataFrame({
                    'labels_final': labels_final,  # List of labels (e.g., [1, 2, 3])
                    'atlas_filtered_list': atlas_filtered_list  # List of regions (e.g., ['A', 'B', 'C'])
                })

                # Merge with panda_file to verify regions
                merged_df = pd.merge(
                    new_df,
                    panda_file[['label', 'region']],
                    left_on='labels_final',
                    right_on='label',
                    how='left'  # Keep all rows from new_df, add matches from panda_file
                )

                # Check for mismatches (optional)
                merged_df['is_consistent'] = (
                        merged_df['atlas_filtered_list'] == merged_df['region']
                )
                print("Mismatches:\n", merged_df[~merged_df['is_consistent']])

                # Save the merged DataFrame
                merged_df.to_csv(opj(output_results,'_'.join([runname,'check','fit','matrix.csv'])),index=False)

                connection = []
                value  = []
                # Iterate through the matrix to populate the list
                for num in range(len(atlas_filtered_list)):
                    for num2 in range(num + 1, len(atlas_filtered_list)):
                        connection.append(f"{atlas_filtered_list[num]} to {atlas_filtered_list[num2]}")
                        value.append(correlation_matrix[num, num2])

                nl = str(connection)
                run_cmd.msg(nl, diary_file, 'OKGREEN')
                nl = str(value)
                run_cmd.msg(nl, diary_file, 'OKGREEN')
                nl = str(len(connection))
                run_cmd.msg(nl, diary_file, 'OKGREEN')
                nl = str(len(value))
                run_cmd.msg(nl, diary_file, 'OKGREEN')

                # Create DataFrame from the list
                flattened_df = pd.DataFrame([value], columns=connection)
                flattened_df['ID'] = ID
                flattened_df['Session'] = Session

                # Display the flattened DataFrame
                nl = str(flattened_df)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

                # Optionally, save the flattened DataFrame to a CSV file
                flattened_df.to_csv(opj(output_results,'_'.join([runname,'flattened','correlation','matrix.csv'])), index=False)

            else:
                nl = 'WARNING: ' + str(func_filename) + ' not found!!'
                run_cmd.msg(nl, diary_file, 'WARNING')

