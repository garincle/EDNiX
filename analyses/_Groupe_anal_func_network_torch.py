import nilearn
import numpy as np
from nilearn.image import iter_img
from nilearn import plotting
import glob
import subprocess
import os
import torch
from torch.nn.functional import softmax
import nibabel as nib
import numpy.ma as ma
from nilearn.masking import compute_epi_mask
import matplotlib.pyplot as plt

# Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

#################################################################################################
#### Dynamic and Statistical Network Analysis with PyTorch
#################################################################################################
def dicstat_torch(BASE_SS, oversample_map, mask_func, folder_atlases, cut_coordsX, cut_coordsY, alpha, component_list,
                  cut_coordsZ, bids_dir, images_dir, mean_imgs, min_size, lower_cutoff, upper_cutoff, afni_sif, s_bind):

    #### Preparing Output Directories
    output_results1 = opj(bids_dir, 'Results')
    if not os.path.exists(output_results1):
        os.mkdir(output_results1)

    if oversample_map:
        studytemplatebrain = BASE_SS
    else:
        studytemplatebrain = opj(folder_atlases, 'BASE_SS_fMRI.nii.gz')

    mean_imgs_rs = nilearn.image.concat_imgs(mean_imgs, ensure_ndim=None, memory=None, memory_level=0, auto_resample=True, verbose=0)
    mask_img = compute_epi_mask(mean_imgs_rs, lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff,
                                connected=True, opening=1, exclude_zeros=True, ensure_finite=True)
    mask_img.to_filename(opj(output_results1, 'mask_mean_func.nii.gz'))

    command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + opj(output_results1, 'mask_mean_func.nii.gz') + \
              ' -input ' + opj(output_results1, 'mask_mean_func.nii.gz') + ' -fill_holes'
    spgo(command)

    command = 'singularity run' + s_bind + afni_sif + '3dresample -master ' + mask_func + ' -prefix ' + opj(output_results1, 'mask_mean_func.nii.gz') + \
              ' -input ' + opj(output_results1, 'mask_mean_func.nii.gz') + ' -overwrite -bound_type SLAB'
    spgo(command)

    command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + opj(output_results1, 'mask_mean_func.nii.gz') + \
              ' -b ' + mask_func + ' -expr "a*b" -prefix ' + opj(output_results1, 'mask_mean_func_overlapp.nii.gz') + ' -overwrite'
    spgo(command)

    #### Dynamic and Statistical Analysis using PyTorch
    dynamic_result_dir = opj(output_results1, 'dynamic_analysis')
    if not os.path.exists(dynamic_result_dir):
        os.mkdir(dynamic_result_dir)

    for image in iter_img(mean_imgs):
        image_data = image.get_fdata()
        tensor_data = torch.tensor(image_data, dtype=torch.float32)

        # Compute Dynamic Network (Mean Activation)
        mean_activation = torch.mean(tensor_data, dim=(0, 1, 2))
        dynamic_network = softmax(mean_activation, dim=0).numpy()

        dynamic_network_img = nilearn.image.new_img_like(image, dynamic_network, copy_header=True)
        dynamic_network_img.to_filename(opj(dynamic_result_dir, 'dynamic_network.nii.gz'))

        display = plotting.plot_stat_map(dynamic_network_img, bg_img=studytemplatebrain, display_mode='z', threshold=0.2)
        display.savefig(opj(dynamic_result_dir, 'dynamic_network_plot.jpg'))
        display.close()
        plt.close('all')

    #### Statistical Network Extraction
    for component in component_list:
        result_dir = opj(output_results1, 'stat_network_' + str(component))
        if not os.path.exists(result_dir):
            os.mkdir(result_dir)

        print(component)
        tensor_data = []

        # Load all images into tensors
        for image in iter_img(images_dir):
            image_data = image.get_fdata()
            tensor_data.append(torch.tensor(image_data, dtype=torch.float32))

        tensor_data = torch.stack(tensor_data)

        # Compute Statistical Components using PCA
        mean_tensor = torch.mean(tensor_data, dim=0)
        u, s, v = torch.pca_lowrank(tensor_data, q=component)
        components = torch.matmul(tensor_data, v[:, :component])

        for i in range(component):
            component_data = components[..., i].numpy()

            # Save and visualize component
            component_img = nilearn.image.new_img_like(mean_imgs_rs, component_data, copy_header=True)
            component_img.to_filename(opj(result_dir, f'component_{i + 1}.nii.gz'))

            display = plotting.plot_stat_map(component_img, bg_img=studytemplatebrain, display_mode='mosaic', threshold=0.2)
            display.savefig(opj(result_dir, f'component_{i + 1}_plot.jpg'))
            display.close()
            plt.close('all')

    print("Dynamic and Statistical Analysis Completed.")
