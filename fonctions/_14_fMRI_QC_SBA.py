import subprocess
import os
import numpy as np
import nibabel as nib
from sklearn.cluster import KMeans
from fonctions.extract_filename import extract_filename

# Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
spco = subprocess.check_output
opd = os.path.dirname
ope = os.path.exists

#################################################################################################
####Seed base analysis
#################################################################################################
def fMRI_QC_SBA(Seed_name, BASE_SS_coregistr, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
    dir_fMRI_Refth_RS_prepro3, RS, nb_run, selected_atlases, panda_files, oversample_map, use_cortical_mask_func):

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

                ##########################################################################

                output_folder = opj(output_results, Seed_name + '/')
                if not os.path.exists(output_folder): os.mkdir(output_folder)

                if os.path.exists(output_folder + '/' + root_RS + '_correlations_fish.nii.gz'):

                    image = nib.load(output_folder + '/' + root_RS + '_correlations_fish.nii.gz').get_data()
                    # Reshape the image data to a 2D array where each row is a voxel with its intensity values
                    # Remove the unnecessary dimension if it exists
                    if image.ndim == 4 and image.shape[-1] == 1:
                        image = image[:, :, :, 0]


                    voxels = image.reshape(-1, 1)

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
                    print("Global Clustering Index:", global_clustering_index)


                    '''
    
                    from scipy.spatial import distance
                    from scipy.stats import zscore
    
                    def calculate_morans_i(image):
                        # Flatten the image and get the coordinates
                        flat_image = image.flatten()
                        z_image = zscore(flat_image)
                        coordinates = np.indices(image.shape).reshape(3, -1).T
    
                        # Calculate spatial weights using inverse distance
                        dist_matrix = distance.cdist(coordinates, coordinates, metric='euclidean')
                        np.fill_diagonal(dist_matrix, np.inf)
                        weight_matrix = 1 / dist_matrix
                        weight_matrix[dist_matrix == np.inf] = 0
    
                        # Calculate Moran's I
                        N = flat_image.size
                        W = weight_matrix.sum()
                        I = (N / W) * (np.sum(weight_matrix * np.outer(z_image, z_image)) / np.sum(z_image ** 2))
    
                        return I
    
                    # Calculate Moran's I for the clustered image
                    morans_i = calculate_morans_i(clustered_image)
                    print("Moran's I:", morans_i)
    '''

                    # Open the file in append mode and add the new line
                    with open(opd(output_results) + '/fMRI_QC/'+ root_RS + 'QC_result.txt', 'a+') as file:
                        # Move the pointer to the end of the file
                        file.seek(0, 2)
                        # Check if the file already ends with a newline
                        if file.tell() > 0:
                            file.seek(-1, 2)
                            last_char = file.read(1)
                            if last_char != '\n':
                                file.write('\n')
                        file.write("Global Clustering Index:", global_clustering_index)

                else:
                    print("correlation file of the seed does not exists, please check that")






