import nilearn
import numpy as np
from nilearn.image import iter_img
from nilearn import plotting
import glob
import subprocess
import os
####DL
from nilearn.decomposition import DictLearning
from nilearn.image import threshold_img
from nilearn import regions
import nibabel as nib
import numpy.ma as ma
from nilearn.masking import compute_epi_mask
import matplotlib.pyplot as plt

#Path to the excels files and data structure
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
def dicstat(BASE_SS, oversample_map, mask_func, folder_atlases, cut_coordsX, cut_coordsY, alpha, component_list,
              cut_coordsZ, bids_dir, images_dir, mean_imgs, min_size, lower_cutoff, upper_cutoff, afni_sif, s_bind):

    #### DL analysis and functional atlas building
    output_results1 = opj(bids_dir, 'Results')
    if not os.path.exists(output_results1): os.mkdir(output_results1)

    if oversample_map == True:
        studytemplatebrain = BASE_SS
    else:
        studytemplatebrain = opj(folder_atlases, 'BASE_SS_fMRI.nii.gz')

    mean_imgs_rs = nilearn.image.concat_imgs(mean_imgs, ensure_ndim=None, memory=None, memory_level=0, auto_resample=True, verbose=0)
    mask_img = compute_epi_mask(mean_imgs_rs,
                                lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff,
                                connected=True, opening=1,
                                exclude_zeros=True, ensure_finite=True)
    mask_img.to_filename(opj(output_results1, 'mask_mean_func.nii.gz'))

    command = 'singularity run' + s_bind + afni_sif + '3dmask_tool -overwrite -prefix ' + opj(output_results1, 'mask_mean_func.nii.gz') + \
    ' -input ' + opj(output_results1, 'mask_mean_func.nii.gz') + ' -fill_holes'
    nl = spgo(command)
    print(nl)

    command = 'singularity run' + s_bind + afni_sif + '3dresample -master ' + mask_func + ' -prefix ' + opj(output_results1, 'mask_mean_func.nii.gz') + ' -input ' + opj(output_results1, 'mask_mean_func.nii.gz') + ' -overwrite -bound_type SLAB'
    nl = spgo(command)
    print(nl)

    command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + opj(output_results1, 'mask_mean_func.nii.gz') + \
                ' -b ' + mask_func + \
              ' -expr "a*b" -prefix ' + opj(output_results1, 'mask_mean_func_overlapp.nii.gz') + ' -overwrite'
    nl = spgo(command)
    print(nl)


    for component in component_list:
        result_dir = opj(output_results1, 'dicL'+ str(component))
        if not os.path.exists(result_dir): os.mkdir(result_dir)

        print(component)
        dict_learning = DictLearning(mask= opj(output_results1, 'mask_mean_func_overlapp.nii.gz'),n_components=component, alpha=alpha,
                                     verbose=10, random_state=0, n_jobs=-3, smoothing_fwhm=False, detrend=False)

        components_imgs = []
        dict_learning.fit(images_dir)
        print('[Example] Saving results')
        # Decomposition dict_learning embeds their own masker
        masker = dict_learning.masker_
        # Drop output maps to a Nifti   file
        components_img = masker.inverse_transform(dict_learning.components_)
        components_img.to_filename(result_dir + '/DL' + str(component) + 'cpts_DicL.nii.gz')

        #load images
        Dl_i = result_dir + '/DL' + str(component) + 'cpts_DicL.nii.gz'

        # Two types of strategies can be used from this threshold function
        # Type 1: strategy used will be based on scoreatpercentile

        for i, cur_img in enumerate(iter_img(Dl_i)):
            tmap_filename = (result_dir + '/network' + str(i) + 'dl.nii.gz')
            cur_img.to_filename(tmap_filename)

            display = plotting.plot_stat_map(tmap_filename,
                                             colorbar=True, bg_img=studytemplatebrain,
                                             display_mode='mosaic',
                                             cut_coords=(len(cut_coordsY), len(cut_coordsY), len(cut_coordsY)))
            display.savefig(result_dir + '/groupttest_mosaic_' + str(i) + '_all.jpg')
            display.close()
            plt.close('all')

            ## transform network into masks
            extracted_data = nib.load(tmap_filename).get_fdata()
            labelnetwork = np.where(extracted_data!=0, 1, 0)
            labeled_img = nilearn.image.new_img_like(tmap_filename, labelnetwork, copy_header=True)
            labeled_img.to_filename(result_dir + '/network_mask' + str(i) + 'dl.nii.gz')

            # Plot the generated mask using the mask_img_ attribute
            thres_connect_img =  regions.connected_label_regions(result_dir + '/network_mask' + str(i) + 'dl.nii.gz', connect_diag=False, min_size=min_size)
            extracted_data2 = thres_connect_img.get_fdata()
            roi_sup0masklab = extracted_data2>0

            #with stat values
            extracted_data = nib.load(tmap_filename).get_fdata()
            labelnetwork = np.where(roi_sup0masklab, extracted_data, 0)
            labeled_img2 = nilearn.image.new_img_like(result_dir + '/network_mask' + str(i) + 'dl.nii.gz', labelnetwork, copy_header=True)
            labeled_img2.to_filename(result_dir + '/network_thresh_mask_size_statmap' + str(i) + 'dl.nii.gz')


        #slelect arg max to concatenate regions per image into an atlas(choose components to remove(if it's needed (not here))
        #img as base to produce a standart  img (can be change)
        networks =  glob.glob(result_dir + '/network_thresh_mask_size_statmap*')
        networks_concate = nilearn.image.concat_imgs(networks)
        networks_concate.to_filename(result_dir + '/concate_network_clean.nii.gz')
        netw6_img = nib.load(result_dir + '/concate_network_clean.nii.gz')
        netw6_extracted_network = netw6_img.get_fdata()
        #choose components to remove

        extracted_datacorect = []
        for n in list(range(0, component)):
            thresh_indexmask = ma.masked_not_equal(netw6_extracted_network[:, :, :, n], 0)
            applymask = np.where(thresh_indexmask, netw6_extracted_network[:, :, :, n], 0)
            extracted_datacorect.append(applymask)

        extracted_datacorect = np.stack(extracted_datacorect, axis=-1)
        #arg maximum when overlapp
        labeled_data = np.argmax(np.abs(extracted_datacorect), axis=-1)+1
        labeled_data[np.max(np.abs(extracted_datacorect), axis=-1)==0]=0


        #save img
        labeled_img = nilearn.image.new_img_like(netw6_img, labeled_data, copy_header=True)
        if not os.path.exists(result_dir + 'atlas/'): os.mkdir(result_dir + 'atlas/')
        labeled_img.to_filename(result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat.nii.gz')

        display = plotting.plot_roi(
            result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat.nii.gz',
            view_type="contours",
            bg_img=studytemplatebrain,
            display_mode='mosaic',
            cut_coords=(len(cut_coordsY), len(cut_coordsY), len(cut_coordsY)))
        display.savefig(result_dir + '/func_atlas.jpg')
        display.close()
        plt.close('all')

        from nilearn.regions import connected_label_regions
        region_labels_min_size = connected_label_regions(result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat.nii.gz',
         min_size=100, connect_diag=False)
        region_labels_min_size.to_filename(result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat_break.nii.gz')

        img2 = nib.load(result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat_break.nii.gz')
        #save img
        extracted_data2 = img2.get_fdata()
        labeled_img2 = nilearn.image.new_img_like(result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat.nii.gz', extracted_data2, copy_header=True)
        labeled_img2.to_filename(result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat_break.nii.gz')
