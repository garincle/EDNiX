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
def tteststat(BASE_SS, selected_atlases, panda_files, oversample_map, use_cortical_mask_func, cut_coordsX, cut_coordsY, cut_coordsZ, threshold_val, bids_dir, all_ID):

    for panda_file, atlas in zip(panda_files, selected_atlases):
        for colomn, row in panda_file.T.iteritems():
            Seed_name = row['region']
            Seed_label = row['label']

            for ID in pd.unique(all_ID):
                print(ID)
                # get useful informations
                List_correlations = sorted(glob.glob(bids_dir + '/sub-' + str(ID) + '/**/func/01_prepro/03_atlas_space/10_Results/SBA/' + Seed_name + '/**_task-rest_bold_correlations_fish.nii.gz'))


                if len(List_correlations)>1:
                    direction_results = opd(opd(opd(opd(List_correlations[0]))))
                    if oversample_map == True:
                        studytemplatebrain = BASE_SS
                    else:
                        studytemplatebrain = opj(direction_results,'BASE_SS_fMRI.nii.gz')

                    if use_cortical_mask_func == True:
                        cortical_mask_func = opj(direction_results,'Gmask.nii.gz')
                    else:
                        cortical_mask_func = opj(direction_results,'mask_brain.nii.gz')

                    output_results1 = opj(opd(opd(opd(opd((direction_results))))), 'Results')
                    if not os.path.exists(output_results1): os.mkdir(output_results1)
                    output_results2 = opj(output_results1, 'SBA')
                    if not os.path.exists(output_results2): os.mkdir(output_results2)
                    output_results = opj(output_results2, Seed_name)
                    if not os.path.exists(output_results): os.mkdir(output_results)

                    command = '3dresample -overwrite' + \
                              ' -prefix ' + opj(direction_results, 'mask_stat.nii.gz') + \
                              ' -master ' + List_correlations[0] + \
                              ' -rmode Cu -input  ' + cortical_mask_func
                    spco([command], shell=True)
                    if os.path.exists(output_results + '/ttest-stat_fisher.nii.gz'):
                        os.remove(output_results + '/ttest-stat_fisher.nii.gz')
                    # '-mask', opj(direction_results,'mask_stat.nii.gz')
                    spco(['3dttest++', '-setA', subprocess.list2cmdline(List_correlations), '-toz', '-prefix', output_results + '/ttest-stat_fisher.nii.gz'])

                    z_map = output_results + '/ttest-stat_fisher.nii.gz'
                    z_map_data = nib.load(z_map).get_data()

                    extracted_datacorect = []
                    for n in list(range(0, 2)):
                        z_map_data_img = (z_map_data[:, :, :, :, n])
                        maxim = z_map_data_img.max()
                        labeled_img = nilearn.image.new_img_like(z_map, z_map_data_img, copy_header=True)
                        labeled_img.to_filename(output_results + '/' + str(n) + 'ttest-stat_fisher.nii.gz')


            List_correlations_all = sorted(glob.glob(bids_dir + '/**/Results/SBA/' + Seed_name + '/**1ttest-stat_fisher.nii.gz'))

            output_results1 = opj(bids_dir, 'Results')
            if not os.path.exists(output_results1): os.mkdir(output_results1)
            output_results2 = opj(output_results1, 'SBA')
            if not os.path.exists(output_results2): os.mkdir(output_results2)
            output_results = opj(output_results2, Seed_name)
            if not os.path.exists(output_results): os.mkdir(output_results)
            # '-mask', opj(direction_results, 'mask_stat.nii.gz'),
            if os.path.exists(output_results + '/ttest-stat_fisher.nii.gz'):
                os.remove(output_results + '/ttest-stat_fisher.nii.gz')
            spco(['3dttest++', '-setA', subprocess.list2cmdline(List_correlations_all), '-toz', '-prefix',
                  output_results + '/ttest-stat_fisher.nii.gz'])

            z_map = output_results + '/ttest-stat_fisher.nii.gz'
            z_map_data = nib.load(z_map).get_data()

            extracted_datacorect = []
            for n in list(range(0, 2)):
                z_map_data_img = (z_map_data[:, :, :, :, n])
                maxim = z_map_data_img.max()
                labeled_img = nilearn.image.new_img_like(z_map, z_map_data_img, copy_header=True)
                labeled_img.to_filename(output_results + '/' + str(n) + 'ttest-stat_fisher.nii.gz')

            ##########################################################################

            ####remove a percentage of the zmap
            threshold_val99 = 99
            loadimg = nib.load(output_results + '/1ttest-stat_fisher.nii.gz').get_data()
            loadimgsort99 = np.percentile(np.abs(loadimg)[np.abs(loadimg) > 0], threshold_val99)

            display = plotting.plot_stat_map(output_results + '/1ttest-stat_fisher.nii.gz', threshold=threshold_val, vmax=loadimgsort99,
                                             colorbar=True, bg_img=studytemplatebrain, display_mode='x',
                                             cut_coords=cut_coordsX)
            display.savefig(output_results + '/groupttest_x_.jpg')
            display.close()
    
            display = plotting.plot_stat_map(output_results + '/1ttest-stat_fisher.nii.gz', threshold=threshold_val, vmax=loadimgsort99,
                                             colorbar=True, bg_img=studytemplatebrain, display_mode='y',
                                             cut_coords=cut_coordsY)
            display.savefig(output_results + '/groupttest_y_.jpg')
            display.close()
    
            display = plotting.plot_stat_map(output_results + '/1ttest-stat_fisher.nii.gz', threshold=threshold_val, vmax=loadimgsort99,
                                             colorbar=True, bg_img=studytemplatebrain, display_mode='z',
                                             cut_coords=cut_coordsZ)
            display.savefig(output_results + '/groupttest_z_.jpg')
            display.close()
    
            display = plotting.plot_stat_map(output_results + '/1ttest-stat_fisher.nii.gz', threshold=threshold_val, vmax=loadimgsort99,
                                             colorbar=True, bg_img=studytemplatebrain,
                                             display_mode='mosaic', cut_coords=(len(cut_coordsY), len(cut_coordsY), len(cut_coordsY)))
            display.savefig(output_results + '/groupttest_mosaic__.jpg')
            display.close()




            #### try all at ones!

            List_correlations_all = sorted(glob.glob(bids_dir + '/sub-**/**/func/01_prepro/03_atlas_space/10_Results/SBA/' + Seed_name + '/**_task-rest_bold_correlations_fish.nii.gz'))
            output_results1 = opj(bids_dir, 'Results')
            if not os.path.exists(output_results1): os.mkdir(output_results1)
            output_results2 = opj(output_results1, 'SBA')
            if not os.path.exists(output_results2): os.mkdir(output_results2)
            output_results = opj(output_results2, Seed_name)
            if not os.path.exists(output_results): os.mkdir(output_results)
            # '-mask', opj(direction_results, 'mask_stat.nii.gz'),
            if os.path.exists(output_results + '/ttest-stat_fisherall.nii.gz'):
                os.remove(output_results + '/ttest-stat_fisherall.nii.gz')
            spco(['3dttest++', '-setA', subprocess.list2cmdline(List_correlations_all), '-toz', '-prefix',
                  output_results + '/ttest-stat_fisherall.nii.gz'])

            z_map = output_results + '/ttest-stat_fisherall.nii.gz'
            z_map_data = nib.load(z_map).get_data()

            extracted_datacorect = []
            for n in list(range(0, 2)):
                z_map_data_img = (z_map_data[:, :, :, :, n])
                maxim = z_map_data_img.max()
                labeled_img = nilearn.image.new_img_like(z_map, z_map_data_img, copy_header=True)
                labeled_img.to_filename(output_results + '/' + str(n) + 'ttest-stat_fisherall.nii.gz')

            ##########################################################################

            ####remove a percentage of the zmap
            threshold_val99 = 99
            loadimg = nib.load(output_results + '/1ttest-stat_fisherall.nii.gz').get_data()
            loadimgsort99 = np.percentile(np.abs(loadimg)[np.abs(loadimg) > 0], threshold_val99)

            display = plotting.plot_stat_map(output_results + '/1ttest-stat_fisherall.nii.gz', threshold=threshold_val,
                                             vmax=loadimgsort99,
                                             colorbar=True, bg_img=studytemplatebrain, display_mode='x',
                                             cut_coords=cut_coordsX)
            display.savefig(output_results + '/groupttest_x_all.jpg')
            display.close()

            display = plotting.plot_stat_map(output_results + '/1ttest-stat_fisherall.nii.gz', threshold=threshold_val,
                                             vmax=loadimgsort99,
                                             colorbar=True, bg_img=studytemplatebrain, display_mode='y',
                                             cut_coords=cut_coordsY)
            display.savefig(output_results + '/groupttest_y_all.jpg')
            display.close()

            display = plotting.plot_stat_map(output_results + '/1ttest-stat_fisherall.nii.gz', threshold=threshold_val,
                                             vmax=loadimgsort99,
                                             colorbar=True, bg_img=studytemplatebrain, display_mode='z',
                                             cut_coords=cut_coordsZ)
            display.savefig(output_results + '/groupttest_z_all.jpg')
            display.close()

            display = plotting.plot_stat_map(output_results + '/1ttest-stat_fisherall.nii.gz', threshold=threshold_val,
                                             vmax=loadimgsort99,
                                             colorbar=True, bg_img=studytemplatebrain,
                                             display_mode='mosaic',
                                             cut_coords=(len(cut_coordsY), len(cut_coordsY), len(cut_coordsY)))
            display.savefig(output_results + '/groupttest_mosaic__all.jpg')
            display.close()