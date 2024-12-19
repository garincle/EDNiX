import nilearn
import subprocess
import pandas as pd
import os
from nilearn.input_data import NiftiLabelsMasker
import matplotlib.pyplot as plt
import numpy as np
from nilearn import plotting
import nibabel as nib
from nilearn.connectome import ConnectivityMeasure
from fonctions.extract_filename import extract_filename
import scipy.ndimage as ndimage
import ants
import datetime
import json

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


#################################################################################################
####Seed base analysis
#################################################################################################
def correl_matrix(dir_fMRI_Refth_RS_prepro1, RS, nb_run, selected_atlases_matrix, segmentation_name_list,
                  ID, Session,
                  bids_dir,s_bind,afni_sif,diary_file):

    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  Working on step ' + str(10) + '(function: _10_Correl_matrix).  ##'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')


    out_results = opj(bids_dir, 'Results')
    out_results_V = opj(out_results, 'fMRI_matrix')

    if not os.path.exists(out_results):
        os.mkdir(out_results)

    if not os.path.exists(out_results_V):
        os.mkdir(out_results_V)

    for panda_file, atlas in zip(segmentation_name_list, selected_atlases_matrix):
        for i in range(0, int(nb_run)):

            root_RS = extract_filename(RS[i])
            func_filename = opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz')

            if ope(func_filename):
                atlas = extract_filename(atlas)
                atlas_filename = opj(dir_fMRI_Refth_RS_prepro1, atlas + '.nii.gz')

                output_results = opj(dir_fMRI_Refth_RS_prepro1, '10_Results')
                if not os.path.exists(output_results):
                    os.mkdir(output_results)
                output_results = opj(dir_fMRI_Refth_RS_prepro1, '10_Results', 'correl_matrix')
                if not os.path.exists(output_results):
                    os.mkdir(output_results)

                nl = 'working on matrix for ' + str(panda_file['region'])
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

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
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    for field, values in differences.items():

                        nl = f"{field}:"
                        print(bcolors.OKGREEN + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')
                        nl = f"  Image 1: {values[0]}"
                        print(bcolors.OKGREEN + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')
                        nl = f"  Image 2: {values[1]}"
                        print(bcolors.OKGREEN + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')

                        caca = nilearn.image.resample_to_img(atlas_filename, func_filename, interpolation='nearest')
                        caca.to_filename(atlas_filename)

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
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                # Load your brain atlas data
                atlas_img = nib.load(atlas_filename)
                atlas_data = atlas_img.get_fdata()

                # Flatten the atlas data for easy comparison
                atlas_flat = np.unique(atlas_data.flatten().astype(int))
                atlas_flat[atlas_flat != 0]

                # Check if labels in the segmentation DataFrame are in the atlas data
                filtered_labels = panda_file[panda_file['label'].isin(atlas_flat)]
                atlas_filtered_list = list(filtered_labels['region'])

                nl =  "INFO: label that you provided"
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                nl =  str(panda_file)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                nl =  "INFO: label found in the atlas after extraction of the signal :" + str(atlas)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                nl =   str(atlas_flat)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                nl =   "INFO: if we combine them (analyse similarity we ware able to keep:"
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                nl =   str(filtered_labels)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                nl =   "INFO: corresponding to the labels"
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                nl =   str(atlas_filtered_list)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                nl =   "INFO: if they are different from the one you provided, you should investigate why!: region too small " + \
                       "for fMRI? error in the labeling? error in the coregistration?, this is not always " + \
                       "a big dill but you should know what happend!!"
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

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

                command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + atlas_filename + ' -expr ' + string_build_atlas2 + \
                          ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, atlas + '_run_' + str(i) + '_filtered.nii.gz') + ' -overwrite'
                nl = spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
                dictionary = {"Sources": [atlas_filename,
                                          func_filename],
                              "Description": string_build_atlas2 + ' (3dcalc, AFNI).'},
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(dir_fMRI_Refth_RS_prepro1, atlas + '_run_' + str(i) + '_filtered.json'), "w") as outfile:
                    outfile.write(json_object)

                ################## EROD the seed if possible ##################
                # Load your atlas image
                atlas_path = opj(dir_fMRI_Refth_RS_prepro1, atlas + '_run_' + str(i) + '_filtered.nii.gz')
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
                    nl = f"Processing label: {label}"
                    print(nl)
                    diary.write(f'\n{nl}')
                    # Isolate the current region (region is 1 where label matches, else 0)
                    region_mask = (atlas_img.numpy() == label).astype(np.uint8)
                    # Try erosion by 1 voxel
                    eroded_region_1 = erode_region(region_mask, iterations=1)
                    if np.sum(eroded_region_1) > 0:  # Check if the region still has voxels
                        nl = f"Label {label}: Eroded by 1 voxel, still has voxels."
                        print(nl)
                        diary.write(f'\n{nl}')
                        final_result[eroded_region_1 > 0] = label
                    else:
                        nl = f"Label {label}: Eroding by 1 voxel removes all voxels, keeping original region."
                        print(nl)
                        diary.write(f'\n{nl}')
                        final_result[region_mask > 0] = label  # Keep the original mask
                    # Try erosion by 2 voxels
                    eroded_region_2 = erode_region(region_mask, iterations=2)

                    if np.sum(eroded_region_2) > 0:  # Check if the region still has voxels
                        nl = f"Label {label}: Eroded by 2 voxels, still has voxels."
                        print(nl)
                        diary.write(f'\n{nl}')
                        final_result[eroded_region_2 > 0] = label
                    else:
                        nl = f"Label {label}: Eroding by 2 voxels removes all voxels, keeping 1 voxel erosion or original region."
                        print(nl)
                        diary.write(f'\n{nl}')

                        # No need to change, we keep the previously applied erosion or original region.
                # Save the final eroded atlas
                nl = final_result
                print(nl)
                diary.write(f'\n{nl}')
                final_atlas = ants.from_numpy(final_result, spacing=atlas_img.spacing, origin=atlas_img.origin,
                                              direction=atlas_img.direction)
                output_path = opj(dir_fMRI_Refth_RS_prepro1, atlas + '_run_' + str(i) + '_filtered_eroded.nii.gz')
                ants.image_write(final_atlas, output_path)
                dictionary = {"Sources": atlas_filename,
                              "Description": 'Erosion (numpy and ANTspy).'},
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(dir_fMRI_Refth_RS_prepro1, atlas + '_run_' + str(i) + '_filtered_eroded.json'), "w") as outfile:
                    outfile.write(json_object)

                nl = f"Eroded atlas saved to: {output_path}"
                print(nl)
                diary.write(f'\n{nl}')

                ##########################################################################
                NAD_masker = NiftiLabelsMasker(labels_img=opj(dir_fMRI_Refth_RS_prepro1, atlas + '_run_' + str(i) + '_filtered_eroded.nii.gz'),
                                                detrend=False,
                                                smoothing_fwhm=None,
                                                standardize=False,
                                                low_pass=None,
                                                high_pass=None,
                                                t_r=None,
                                                memory=None, verbose=5)
                try:
                    time_series = NAD_masker.fit_transform(func_filename)
                except:

                    caca = nilearn.image.resample_to_img(opj(dir_fMRI_Refth_RS_prepro1, atlas + '_filtered.nii.gz'), func_filename, interpolation='nearest')
                    caca.to_filename(opj(dir_fMRI_Refth_RS_prepro1, atlas + '_filtered.nii.gz'))
                    extracted_data = nib.load(opj(dir_fMRI_Refth_RS_prepro1, atlas + '_filtered.nii.gz')).get_fdata()
                    labeled_img2 = nilearn.image.new_img_like(func_filename, extracted_data, copy_header=True)
                    labeled_img2.to_filename(opj(dir_fMRI_Refth_RS_prepro1, atlas + '_filtered.nii.gz'))

                    NAD_masker = NiftiLabelsMasker(
                        labels_img=opj(dir_fMRI_Refth_RS_prepro1, atlas + '_run_' + str(i) + '_filtered.nii.gz'), # must be a mistake ......
                        detrend=False,
                        smoothing_fwhm=None,
                        standardize=False,
                        low_pass=None,
                        high_pass=None,
                        t_r=None,
                        memory=None, verbose=5)

                    time_series = NAD_masker.fit_transform(func_filename)

                correlation_measure = ConnectivityMeasure(kind="partial correlation")
                correlation_matrix = correlation_measure.fit_transform([time_series])[0]

                # Plot the correlation matrix

                # Make a large figure
                # Mask the main diagonal for visualization:
                np.fill_diagonal(correlation_matrix, 0)
                # The labels we have start with the background (0), hence we skip the
                # first label
                # matrices are ordered for block-like representation
                nl = str(correlation_matrix.shape)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                nl = str(len(panda_file['region']))
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                nl = str(correlation_matrix)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

                plotting.plot_matrix(
                    correlation_matrix,
                    figure=(10, 8),
                    labels=atlas_filtered_list,
                    vmax=0.8,
                    vmin=-0.8)
                plt.savefig(opj(output_results,atlas + '_run_' + str(i) + 'matrix.png'))

                # Convert correlation matrix to a DataFrame
                corr_df = pd.DataFrame(correlation_matrix, index=atlas_filtered_list, columns=atlas_filtered_list)

                connection = []
                value  = []
                # Iterate through the matrix to populate the list
                for i in range(len(atlas_filtered_list)):
                    for j in range(i + 1, len(atlas_filtered_list)):
                        connection.append(f"{atlas_filtered_list[i]} to {atlas_filtered_list[j]}")
                        value.append(correlation_matrix[i, j])

                nl = str(connection)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                nl = str(value)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                nl = str(len(connection))
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                nl = str(len(value))
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

                # Create DataFrame from the list
                flattened_df = pd.DataFrame([value], columns=connection)
                flattened_df['ID'] = ID
                flattened_df['Session'] = Session

                # Display the flattened DataFrame
                nl = str(flattened_df)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

                # Optionally, save the flattened DataFrame to a CSV file
                flattened_df.to_csv(output_results + '/' + atlas + '_run_' + str(i) + '_flattened_correlation_matrix.csv', index=False)

            else:
                nl =  'WARNING: ' + str(func_filename) + ' not found!!'
                print(bcolors.WARNING + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

    diary.write(f'\n')
    diary.close()
