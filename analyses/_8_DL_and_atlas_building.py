import nilearn
import numpy as np
from nilearn import image
from nilearn.image import iter_img
from nilearn import plotting
import glob
import subprocess
import nipype
import pandas as pd
import os
from datetime import datetime as dt
from datetime import timedelta as td
####DL
from nilearn.decomposition import DictLearning
from nilearn.decomposition import CanICA
from nilearn.input_data import NiftiMasker
from nilearn.image import threshold_img
from nilearn import regions
import nibabel as nib
import numpy.ma as ma
from nilearn.image import math_img

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
spco = subprocess.check_output
opd = os.path.dirname

import nilearn
import numpy as np
from nilearn import image
from nilearn.image import iter_img
from nilearn import plotting
import glob
import subprocess
import nipype
import pandas as pd
import os
from datetime import datetime as dt
from datetime import timedelta as td
####DL
from nilearn.decomposition import DictLearning
from nilearn.decomposition import CanICA
from nilearn.input_data import NiftiMasker
from nilearn.image import threshold_img
from nilearn import regions
import nibabel as nib
import numpy.ma as ma
from nilearn.image import math_img
from nilearn.input_data import NiftiLabelsMasker
from nilearn.image import resample_to_img

# Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
spco = subprocess.check_output
opd = os.path.dirname

from fonctions.extract_filename import extract_filename


#################################################################################################
####Seed base analysis
#################################################################################################
def dicstat(BASE_SS, oversample_map, use_cortical_mask_func, cut_coordsX, cut_coordsY,
              cut_coordsZ, bids_dir, images_dir, min_size):

    direction_results = opd(images_dir[0])
    #### DL analysis and functional atlas building

    if oversample_map == True:
        studytemplatebrain = BASE_SS
    else:
        studytemplatebrain = opj(direction_results, 'BASE_SS_fMRI.nii.gz')

    if use_cortical_mask_func == True:
        cortical_mask_func = opj(direction_results, 'Gmask.nii.gz')
    else:
        cortical_mask_func = opj(direction_results, 'mask_brain.nii.gz')

    for component in [7, 17]:
        output_results1 = opj(bids_dir, 'Results')
        if not os.path.exists(output_results1): os.mkdir(output_results1)
        result_dir = opj(output_results1, 'dicL'+ str(component))
        if not os.path.exists(result_dir): os.mkdir(result_dir)

        print(component)
        dict_learning = DictLearning(mask= cortical_mask_func,n_components=component, n_epochs=10,  alpha=9,
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
                                             colorbar=True, bg_img=studytemplatebrain, display_mode='x',
                                             cut_coords=cut_coordsX)
            display.savefig(result_dir + '/groupttest_x_' + str(i) + 'all.jpg')
            display.close()

            display = plotting.plot_stat_map(tmap_filename,
                                             colorbar=True, bg_img=studytemplatebrain, display_mode='y',
                                             cut_coords=cut_coordsY)
            display.savefig(result_dir + '/groupttest_y_' + str(i) + 'all.jpg')
            display.close()

            display = plotting.plot_stat_map(tmap_filename,
                                             colorbar=True, bg_img=studytemplatebrain, display_mode='z',
                                             cut_coords=cut_coordsZ)
            display.savefig(result_dir + '/groupttest_z_' + str(i) + 'all.jpg')
            display.close()

            display = plotting.plot_stat_map(tmap_filename,
                                             colorbar=True, bg_img=studytemplatebrain,
                                             display_mode='mosaic',
                                             cut_coords=(len(cut_coordsY), len(cut_coordsY), len(cut_coordsY)))
            display.savefig(result_dir + '/groupttest_mosaic_' + str(i) + '_all.jpg')
            display.close()

            extracted_data = nib.load(tmap_filename).get_data()
            labelnetwork = np.where(extracted_data>0, 1, 0)
            labeled_img = nilearn.image.new_img_like(tmap_filename, labelnetwork, copy_header=True)
            labeled_img.to_filename(result_dir + '/network_mask' + str(i) + 'dl.nii.gz')

            thres_img = threshold_img(tmap_filename, threshold='98%', mask_img=result_dir + '/network_mask' + str(i) + 'dl.nii.gz')
            thres_img.to_filename(result_dir + '/network_thresh' + str(i) + 'dl.nii.gz')

            extracted_data = nib.load(result_dir + '/network_thresh' + str(i) + 'dl.nii.gz').get_data()
            labelnetwork = np.where(extracted_data>0, 1, 0)
            labeled_img = nilearn.image.new_img_like(tmap_filename, labelnetwork, copy_header=True)
            labeled_img.to_filename(result_dir + '/network_mask_tresh' + str(i) + 'dl.nii.gz')

            # Plot the generated mask using the mask_img_ attribute
            thres_connect_img =  regions.connected_label_regions(result_dir + '/network_mask' + str(i) + 'dl.nii.gz', connect_diag=False, min_size=min_size)
            extracted_data2 = thres_connect_img.get_data()
            roi_sup0masklab = extracted_data2>0
            extracted_data = nib.load(result_dir + '/network_mask_tresh' + str(i) + 'dl.nii.gz').get_data()
            labelnetwork = np.where(roi_sup0masklab, extracted_data, 0)
            labeled_img2 = nilearn.image.new_img_like(thres_img, labelnetwork, copy_header=True)
            #with stat values
            extracted_data = nib.load(tmap_filename).get_data()
            labelnetwork = np.where(roi_sup0masklab, extracted_data, 0)
            labeled_img2 = nilearn.image.new_img_like(thres_img, labelnetwork, copy_header=True)
            labeled_img2.to_filename(result_dir + '/network_thresh_mask_size_statmap' + str(i) + 'dl.nii.gz')


        #slelect arg max to concatenate regions per image into an atlas(choose components to remove(if it's needed (not here))
        #img as base to produce a standart  img (can be change)
        networks =  glob.glob(result_dir + '/network_thresh_mask_size_statmap*')
        networks_concate = nilearn.image.concat_imgs(networks)
        networks_concate.to_filename(result_dir + '/concate_network_clean.nii.gz')
        netw6_img = nib.load(result_dir + '/concate_network_clean.nii.gz')
        netw6_extracted_network = netw6_img.get_data()
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



        from nilearn.regions import connected_label_regions
        region_labels_min_size = connected_label_regions(result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat.nii.gz',
         min_size=100, connect_diag=False)
        region_labels_min_size.to_filename(result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat_break.nii.gz')

        img2 = nib.load(result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat_break.nii.gz')
        #save img
        extracted_data2 = img2.get_data()
        labeled_img2 = nilearn.image.new_img_like(result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat.nii.gz', extracted_data2, copy_header=True)
        labeled_img2.to_filename(result_dir + 'atlas/dict_learning_' + str(component) + 'compos_concat_break.nii.gz')
