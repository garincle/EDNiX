import nilearn
import numpy as np
from nilearn.image import iter_img
from nilearn import plotting
import glob
import subprocess
import os
####DL
from nilearn.decomposition import DictLearning
from nilearn import regions
import nibabel as nib
import numpy.ma as ma
from nilearn.masking import compute_epi_mask
import matplotlib.pyplot as plt
import ants
import Tools.Load_EDNiX_requirement
# Plot the scores
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

def dicstat(oversample_map, mask_func, cut_coords, alpha_dic, component_list, oversample_dictionary,
              bids_dir, images_dir, mean_imgs, min_size, lower_cutoff, upper_cutoff, MAIN_PATH, FS_dir, templatelow, templatehigh, TR, smoothing):

    s_path, afni_sif, fsl_sif, fs_sif, itk_sif, wb_sif, strip_sif, s_bind = Tools.Load_EDNiX_requirement.load_requirement(
        MAIN_PATH, bids_dir, FS_dir)

    #### DL analysis and functional atlas building
    output_results1 = opj(bids_dir, 'Results')
    if not os.path.exists(output_results1): os.mkdir(output_results1)

    if oversample_map == True:
        studytemplatebrain = templatehigh
    else:
        studytemplatebrain = templatelow

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

    # Resample to match the mask function
    command = f"singularity run {s_bind} {afni_sif} 3dresample -master {opj(output_results1, 'mask_mean_func.nii.gz')} -prefix {opj(output_results1, 'mask_mean_func_orig.nii.gz')} " \
              f"-input {mask_func} -overwrite -bound_type SLAB"
    nl = spgo(command)
    print(nl)

    command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + opj(output_results1, 'mask_mean_func.nii.gz') + \
                ' -b ' + opj(output_results1, 'mask_mean_func_orig.nii.gz') + \
              ' -expr "a*b" -prefix ' + opj(output_results1, 'mask_mean_func_overlapp.nii.gz') + ' -overwrite'
    nl = spgo(command)
    print(nl)

    if oversample_dictionary == True:
        command = f"singularity run {s_bind} {afni_sif} 3dresample -master {studytemplatebrain} -prefix {opj(output_results1, 'mask_mean_func_overlapp.nii.gz')} " \
                  f"-input {opj(output_results1, 'mask_mean_func_overlapp.nii.gz')} -overwrite -bound_type SLAB"
        nl = spgo(command)
        print(nl)

    lowresanat = ants.image_read(templatelow)  # Low-resolution atlas
    anat = ants.image_read(templatehigh)  # High-resolution anatomical image
    # Perform affine + SyN nonlinear registration
    reg = ants.registration(
        fixed=anat,  # High-resolution anatomical image (target)
        moving=lowresanat,  # Low-resolution atlas (source)
        type_of_transform="SyN")  # Symmetric normalization (nonlinear)

    # Apply transformation with nearest-neighbor interpolation to preserve labels
    registered_atlas = ants.apply_transforms(
        fixed=anat,
        moving=lowresanat,
        transformlist=reg["fwdtransforms"],
        invert_transform_flags=[False],  # Apply inverse transformation
        interpolator="nearestNeighbor")  # Critical for label preservation
    # Save output
    registered_atlas.to_filename( opj(output_results1, 'low_to_highR_template.nii.gz'))

    for component in component_list:
        result_dir = opj(output_results1, 'dicL'+ str(component))
        if not os.path.exists(result_dir): os.mkdir(result_dir)

        if oversample_dictionary == True:
            dict_learning = DictLearning(mask= opj(output_results1, 'mask_mean_func_overlapp.nii.gz'),n_components=component, alpha=alpha_dic, n_epochs=1,
                                         verbose=10, standardize="zscore_sample", random_state=0, n_jobs=1, smoothing_fwhm=smoothing, detrend=False, t_r=TR)
        else:
            dict_learning = DictLearning(mask=opj(output_results1, 'mask_mean_func_overlapp.nii.gz'),
                                         n_components=component, alpha=alpha_dic, batch_size=5, standardize="zscore_sample", n_epochs=1,
                                         verbose=10, random_state=0, n_jobs=1, smoothing_fwhm=smoothing, detrend=False,
                                         t_r=TR)
        dict_learning.fit(images_dir)
        print('[Example] Saving results')
        # Decomposition dict_learning embeds their own masker
        masker = dict_learning.masker_
        # Drop output maps to a Nifti   file
        components_img = masker.inverse_transform(dict_learning.components_)
        components_img.to_filename(result_dir + '/DL' + str(component) + 'cpts_DicL.nii.gz')
        #load images
        Dl_i = result_dir + '/DL' + str(component) + 'cpts_DicL.nii.gz'

        scores = dict_learning.score(images_dir, per_component=True)

        plt.figure(figsize=(4, 4), constrained_layout=True)
        positions = np.arange(len(scores))
        plt.barh(positions, scores)
        plt.ylabel("Component #", size=12)
        plt.xlabel("Explained variance", size=12)
        plt.yticks(np.arange(20))
        plt.gca().xaxis.set_major_formatter(FormatStrFormatter("%.3f"))
        plt.savefig(result_dir + '/explained_var.jpg')
        plt.close('all')

        for i, cur_img in enumerate(iter_img(Dl_i)):
            tmap_filename = (result_dir + '/network' + str(i) + 'dl.nii.gz')
            cur_img.to_filename(tmap_filename)

            display = plotting.plot_stat_map(tmap_filename, symmetric_cbar=False,
                                             colorbar=True, bg_img=studytemplatebrain,
                                             display_mode='mosaic',
                                             cut_coords=cut_coords)
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

        from scipy.ndimage import label, generic_filter

        def mode_filter(values):
            """ Returns the most frequent nonzero value in the neighborhood """
            unique, counts = np.unique(values[values > 0], return_counts=True)
            return unique[np.argmax(counts)] if unique.size > 0 else 0

        def reassign_small_regions(labeled_data, statistical_scores, min_size=50):
            """ Reassigns small regions to surrounding larger labels while preserving zero-value areas. """

            # Create a mask of voxels where statistical scores are 0
            zero_mask = (statistical_scores == 0)

            # Copy labeled data to modify
            new_labeled_data = labeled_data.copy()

            unique_labels = np.unique(labeled_data)
            for lbl in unique_labels:
                if lbl == 0:  # Skip background
                    continue

                # Identify connected components within the current label
                mask = (labeled_data == lbl)
                labeled_mask, num_features = label(mask)

                for i in range(1, num_features + 1):
                    region_size = np.sum(labeled_mask == i)
                    if region_size < min_size:
                        # Small region: reassign based on local majority vote
                        small_region_mask = (labeled_mask == i)
                        new_labeled_data[small_region_mask] = 0  # Temporarily remove small region

            # Apply mode filtering to fill in gaps with surrounding majority label
            new_labeled_data = generic_filter(new_labeled_data, mode_filter, size=5)

            # Restore original zero-value voxels
            new_labeled_data[zero_mask] = 0

            return new_labeled_data

        # Apply the function, using the original statistical scores as a reference
        labeled_data_corrected = reassign_small_regions(labeled_data, np.max(np.abs(extracted_datacorect), axis=-1),
                                                        min_size=50)
        # save img
        labeled_img = nilearn.image.new_img_like(netw6_img, labeled_data_corrected, copy_header=True)
        if not os.path.exists(result_dir + 'atlas/'): os.mkdir(result_dir + 'atlas/')
        labeled_img.to_filename(result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat_smoothed.nii.gz')

        # Resample to match the mask function
        command = f"singularity run {s_bind} {afni_sif} 3dresample -master {mask_func} -prefix {result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat_anatR.nii.gz'} " \
                  f"-input {result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat.nii.gz'} -overwrite -bound_type SLAB"
        nl = spgo(command)
        print(nl)

        display = plotting.plot_roi(
            result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat.nii.gz',
            view_type="contours",
            bg_img=studytemplatebrain,
            display_mode='mosaic',
            cut_coords=cut_coords)
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

        atlas = ants.image_read(result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat.nii.gz')
        # Apply transformation with nearest-neighbor interpolation to preserve labels
        registered_atlas = ants.apply_transforms(
            fixed=anat,
            moving=atlas,
            transformlist=reg["fwdtransforms"],
            invert_transform_flags=[False],  # Apply inverse transformation
            interpolator="nearestNeighbor" ) # Critical for label preservation
        # Save output
        registered_atlas.to_filename(result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat_anatR.nii.gz')

        command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat_anatR.nii.gz' + \
                ' -b ' + templatehigh  + \
                ' -expr "a*step(b)" -prefix ' + result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat_anatR.nii.gz' + ' -overwrite'
        nl = spgo(command)
        print(nl)


        atlas = ants.image_read(result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat_smoothed.nii.gz')
        # Apply transformation with nearest-neighbor interpolation to preserve labels
        registered_atlas = ants.apply_transforms(
            fixed=anat,
            moving=atlas,
            transformlist=reg["fwdtransforms"],
            invert_transform_flags=[False],  # Apply inverse transformation
            interpolator="nearestNeighbor" ) # Critical for label preservation
        # Save output
        registered_atlas.to_filename(result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat_smoothed_anatR.nii.gz')

        command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat_smoothed_anatR.nii.gz' + \
                ' -b ' + templatehigh  + \
                ' -expr "a*step(b)" -prefix ' + result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat_smoothed_anatR.nii.gz' + ' -overwrite'
        nl = spgo(command)
        print(nl)