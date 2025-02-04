import subprocess
import os
import numpy as np
import nibabel as nib
from sklearn.cluster import KMeans
import shutil
from fonctions.extract_filename import extract_filename
import datetime

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
opi = os.path.isfile
#################################################################################################
####Seed base analysis
#################################################################################################
def fMRI_QC_SBA(SBAspace, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
    dir_fMRI_Refth_RS_prepro3, RS, nb_run, selected_atlases, panda_files,diary_file):

    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  Working on step ' + str(13) + '(function: _13_fMRI_QC_SBA).  ##'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')


    for panda_file, atlas in zip(panda_files, selected_atlases):
        for space in SBAspace:
            if space == 'func':
                direction_results = dir_fMRI_Refth_RS_prepro1
                output_results = opj(dir_fMRI_Refth_RS_prepro1, '10_Results')
                if not ope(output_results): os.mkdir(output_results)
                output_results = opj(dir_fMRI_Refth_RS_prepro1, '10_Results','SBA')
                if not ope(output_results): os.mkdir(output_results)

            if space == 'anat':
                direction_results = dir_fMRI_Refth_RS_prepro2
                output_results = opj(dir_fMRI_Refth_RS_prepro2, '10_Results')
                if not ope(output_results): os.mkdir(output_results)
                output_results = opj(dir_fMRI_Refth_RS_prepro2, '10_Results','SBA')
                if not ope(output_results): os.mkdir(output_results)

            if space == 'atlas':
                direction_results = dir_fMRI_Refth_RS_prepro3
                output_results = opj(dir_fMRI_Refth_RS_prepro3, '10_Results')
                if not ope(output_results): os.mkdir(output_results)
                output_results = opj(dir_fMRI_Refth_RS_prepro3, '10_Results','SBA')
                if not ope(output_results): os.mkdir(output_results)
            else:
                nl = 'WARNING: will not perform ' + str(direction_results) + ' space because SBAspace is ' + str(SBAspace)
                print(bcolors.WARNING + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

            if not ope(output_results): os.mkdir(output_results)
            output_results_result = opj(output_results, 'fMRI_QC_SBA')
            print(output_results_result)
            diary.write(f'\n{output_results_result}')

            if ope(output_results_result):
                shutil.rmtree(output_results_result)
            os.mkdir(output_results_result)

            for i in range(0, int(nb_run)):
                root_RS = extract_filename(RS[i])
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

                    if ope(opj(output_folder,root_RS + '_correlations_fish.nii.gz')):
                        image = nib.load(opj(output_folder,root_RS + '_correlations_fish.nii.gz')).get_fdata()
                        # Reshape the image data to a 2D array where each row is a voxel with its intensity values
                        # Remove the unnecessary dimension if it exists
                        if image.ndim == 4 and image.shape[-1] == 1:
                            image = image[:, :, :, 0]

                        voxels = image.reshape(-1, 1)
                        # Replace inf and -inf with NaN
                        voxels = np.where(np.isinf(voxels), np.nan, voxels)
                        # Compute the mean of the non-NaN values in each column
                        mean_values = np.nanmean(voxels, axis=0)
                        # Replace NaN values with the mean of each column
                        voxels = np.where(np.isnan(voxels), mean_values, voxels)

                        # Apply K-means clustering
                        num_clusters = 5
                        kmeans = KMeans(n_clusters=num_clusters, random_state=42).fit(voxels)
                        labels = kmeans.labels_
                        clustered_image = labels.reshape(image.shape)

                        # Function to calculate local clustering capacity in 3D
                        def local_clustering_capacity_3d(clustered_img, window_size=3):
                            depth, height, width = clustered_img.shape
                            capacity = np.zeros((depth, height, width))
                            offset = window_size // 2

                            for d in range(offset, depth - offset):
                                for i in range(offset, height - offset):
                                    for j in range(offset, width - offset):
                                        local_window = clustered_img[d - offset:d + offset + 1, i - offset:i + offset + 1,
                                                       j - offset:j + offset + 1]
                                        central_label = clustered_img[d, i, j]
                                        capacity[d, i, j] = np.sum(local_window == central_label) / (window_size ** 3)

                            return capacity

                        # Calculate local clustering capacity
                        capacity = local_clustering_capacity_3d(clustered_image, window_size=3)

                        # Calculate global clustering index as the average local clustering capacity
                        global_clustering_index = np.mean(capacity)

                        # Open the file in 'r+' mode (read and write) and add the new line
                        with open(opj(output_results,'fMRI_QC_SBA',root_RS + 'QC_result.txt'), 'w') as file:
                            file.seek(0, 2)  # Move the pointer to the end of the file
                            if file.tell() > 0:  # Ensure the file is not empty
                                file.seek(file.tell() - 1)  # Move to the last character
                                last_char = file.read(1)
                                if last_char != '\n':
                                    file.write('\n')  # Add a newline if the last character is not a newline
                            # Write the new content
                            file.write("Global Clustering Index SBA " + Seed_name + " :" + str(global_clustering_index) + '\n')
                            nl = "Global Clustering Index SBA " + Seed_name + " :" + str(global_clustering_index)
                            print(bcolors.OKGREEN + nl + bcolors.ENDC)
                            diary.write(f'\n{nl}')

                    else:
                        nl = "WARNING: correlation file of the seed does not exists, please check that"
                        print(bcolors.WARNING + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')

    diary.write(f'\n')
    diary.close()






