import nilearn
import numpy as np
from nilearn.image import iter_img
import glob
import subprocess
import pandas as pd
import os
####DL
from nilearn.decomposition import DictLearning
from nilearn.image import threshold_img
from nilearn import regions
import nibabel as nib
import numpy.ma as ma


#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
spco = subprocess.check_output
opd = os.path.dirname



#############load data
#mask and stdy template brain
# refs
#/02_WB/volumes/

#############load data
ml_excel = opn('/home/cgarin/Documents/8_Multispecies/fMRI/Chimpanzee/DataLoader/Data_loader.xlsx')


animalinfo = pd.read_excel(ml_excel, sheet_name='animalinfo')
MRIsessions = pd.read_excel(ml_excel, sheet_name='MRIsessions')
allinfo = pd.merge(animalinfo, MRIsessions)


path_ATLAS  = opj('/home/cgarin/Documents/4_Interspecies_fMRI/1_Atlases/Chimpanzee/native/')
Major_fMRI_Refth = opj('/home/cgarin/Documents/8_Multispecies/fMRI/')
study_fMRI_Refth = opj(Major_fMRI_Refth,study)

if not os.path.exists(opn('/home/cgarin/Documents/8_Multispecies/fMRI/Chimpanzee/Results')):
    os.mkdir(opn('/home/cgarin/Documents/8_Multispecies/fMRI/Chimpanzee/Results'))


'''
output_dir_process = opn('/home/cgarin/Documents/8_Multispecies/fMRI/Chimpanzee/Results/1_Seed_base/')
template_filename = opj(path_ATLAS,'0_orig/template.nii.gz')

atlas_filename = '/home/cgarin/Documents/8_Multispecies/13_Atlas_project/New_atlas_PET/Chimpanzee/PMC_func_reso.nii.gz'
Aseg_ref    = '/home/cgarin/Documents/4_Interspecies_fMRI/1_Atlases/Chimpanzee/native/02_WB/volumes/labels/atlas_forSEG_final.nii.gz'
BASE_mask   = opj(path_ATLAS,'0_orig/brain_mask.nii.gz')


####creat mask dilate func
Cortex_edit  = ',3,42,10,49,11,12,13,17,18,26,28,50,51,52,53,54,58,60'
command = '3dcalc -a ' + Aseg_ref + ' -expr "amongst(a' + Cortex_edit + ')" -overwrite -prefix ' + opj(path_ATLAS,'02_WB', 'volumes', 'labels', 'cortical_mask_anat.nii.gz')
spco([command], shell=True)

#dilate MORE the "maskDilat"
command = '3dmask_tool -prefix ' + opj(path_ATLAS,'02_WB', 'volumes', 'labels', 'cortical_mask_anat.nii.gz') + \
' -input ' + opj(path_ATLAS,'02_WB', 'volumes', 'labels', 'cortical_mask_anat.nii.gz') + ' -overwrite -dilate_input 4'
spco(command, shell=True)

command = '3dcalc -a ' + opj(path_ATLAS,'02_WB', 'volumes', 'labels', 'cortical_mask_anat.nii.gz') + ' -b ' + BASE_mask + \
' -prefix ' +  opj(path_ATLAS,'02_WB', 'volumes', 'labels', 'cortical_mask_anat.nii.gz') + ' -expr "a*b" -overwrite'
spco([command], shell=True)

command = '3dresample -master ' + '/home/cgarin/Documents/8_Multispecies/fMRI/Chimpanzee/BIDS/sub-Ranette/ses-1/func/02_residual/02_ICA_Norm/sub-Ranette_ses-1_task-rest_run-01_epi_norm_final.nii.gz' + \
' -prefix ' +  opj(path_ATLAS,'02_WB', 'volumes', 'labels', 'cortical_mask_func_ctx.nii.gz') + \
' -input ' + opj(path_ATLAS,'02_WB', 'volumes', 'labels', 'cortical_mask_anat.nii.gz') + ' -overwrite'
spco([command], shell=True)

func_filenames = []
for ID, Session, DICOMdir in zip(allinfo.ID, allinfo.Session, allinfo.DICOMdir):
    # Organization of the folders
    data_fMRI_Refth = opj(Major_fMRI_Refth,study,'BIDS','sub-' + ID,'ses-' + str(Session))
    dir_fMRI_Refth_RS          = opj(data_fMRI_Refth,'func')
    dir_fMRI_Refth_RS_residual = opj(dir_fMRI_Refth_RS,'02_residual')
    dir_RS_norm          = opj(dir_fMRI_Refth_RS_residual,'02_ICA_Norm')
    dir_fMRI_Refth_RS_prepro   = opj(dir_fMRI_Refth_RS,'01_prepro')
    dir_fMRI_Refth_RS_prepro3  = opj(dir_fMRI_Refth_RS_prepro,'03_atlas_space')

    # get useful informations
    list_RS = sorted(glob.glob(opj(dir_fMRI_Refth_RS, '*_task-rest_*epi.nii.gz')))
    nb_run  = len(list_RS)
    RS      = [os.path.basename(i) for i in list_RS]

    for i in range(0, int(nb_run)):
        func_filenames.append(opj(dir_fMRI_Refth_RS_prepro3, RS[i].replace('.nii.gz','_residual_norm_RsemP.nii.gz')))

        command = '3dresample -master ' + '/home/cgarin/Documents/8_Multispecies/fMRI/Chimpanzee/BIDS/sub-Ranette/ses-1/func/02_residual/02_ICA_Norm/sub-Ranette_ses-1_task-rest_run-01_epi_norm_final.nii.gz' + \
        ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro3, RS[i].replace('.nii.gz','_residual_norm_RsemP.nii.gz')) + \
        ' -input ' + opj(dir_fMRI_Refth_RS_prepro3, RS[i].replace('.nii.gz','_residual_norm.nii.gz')) + ' -overwrite'
        spco([command], shell=True)


import nilearn
from nilearn.masking import compute_multi_epi_mask
mask_img = nilearn.masking.compute_multi_epi_mask(func_filenames, lower_cutoff=0.05, upper_cutoff=0.95, \
    connected=True, opening=2, threshold=0.8, target_affine=None, target_shape=None, exclude_zeros=False, n_jobs=1, memory=None, verbose=0)
mask_img.to_filename(opj(path_ATLAS,'02_WB', 'volumes', 'labels', 'mask_img_nilearn.nii.gz'))

command = '3dcalc -a ' + opj(path_ATLAS,'02_WB', 'volumes', 'labels', 'cortical_mask_func_ctx.nii.gz') + \
' -b ' + opj(path_ATLAS,'02_WB', 'volumes', 'labels', 'mask_img_nilearn.nii.gz') + \
' -prefix ' +  opj(path_ATLAS,'02_WB', 'volumes', 'labels', 'cortical_mask_func.nii.gz') + ' -expr "a*b" -overwrite'
spco([command], shell=True)
'''


cortical_mask_func = opj(path_ATLAS,'02_WB', 'volumes', 'labels', 'cortical_mask_func.nii.gz')

#################################################################################################
####Seed base analysis
################################################################################################# 


func_filenames = []
for ID, Session, DICOMdir in zip(allinfo.ID, allinfo.Session, allinfo.DICOMdir):


    # Organization of the folders
    data_fMRI_Refth = opj(Major_fMRI_Refth,study,'BIDS','sub-' + ID,'ses-' + str(Session))
    dir_fMRI_Refth_RS          = opj(data_fMRI_Refth,'func')
    dir_fMRI_Refth_RS_residual = opj(dir_fMRI_Refth_RS,'02_residual')
    dir_RS_norm          = opj(dir_fMRI_Refth_RS_residual,'02_ICA_Norm')

    # get useful informations
    list_RS = sorted(glob.glob(opj(dir_fMRI_Refth_RS, '*_task-rest_*epi.nii.gz')))
    nb_run  = len(list_RS)
    RS      = [os.path.basename(i) for i in list_RS]

    for i in range(0, int(nb_run)):
        '''
        command = '3dresample -master ' + '/home/cgarin/Documents/8_Multispecies/fMRI/Chimpanzee/BIDS/sub-Ranette/ses-1/func/02_residual/02_ICA_Norm/sub-Ranette_ses-1_task-rest_run-01_epi_norm_final.nii.gz' + \
        ' -prefix ' +  opj(dir_RS_norm, RS[i].replace('.nii.gz','_norm_final_RsemP.nii.gz')) + \
        ' -input ' + opj(dir_RS_norm, RS[i].replace('.nii.gz','_norm_final.nii.gz'))
        spco([command], shell=True)
        '''

        func_filenames.append(opj(dir_RS_norm, RS[i].replace('.nii.gz','_norm_final_clean.nii.gz')))
        print(os.path.exists(opj(dir_RS_norm, RS[i].replace('.nii.gz','_norm_final_clean.nii.gz'))))


#### DL analysis and functional atlas building

for component in [7, 17]:

    result_dir = '/home/cgarin/Documents/8_Multispecies/fMRI/Chimpanzee/Results/2_DL/' + str(component) + 'component_preproc/'
    if not os.path.exists(result_dir): os.mkdir(result_dir)

    print(component)
    dict_learning = DictLearning(mask=cortical_mask_func, n_components=component, n_epochs=10,  alpha=10, 
                                 verbose=10, random_state=0, n_jobs=-3, smoothing_fwhm=False, detrend=False)

    components_imgs = []
    dict_learning.fit(func_filenames)
    print('[Example] Saving results')
    # Decomposition dict_learning embeds their own masker
    masker = dict_learning.masker_
    # Drop output maps to a Nifti   file
    components_img = masker.inverse_transform(dict_learning.components_)
    components_img.to_filename(result_dir + '/DL' + str(component) + 'cpts_DicL_yung.nii.gz')


    #load images
    Dl_i = result_dir + '/DL' + str(component) + 'cpts_DicL_yung.nii.gz'


    # Two types of strategies can be used from this threshold function
    # Type 1: strategy used will be based on scoreatpercentile

    for i, cur_img in enumerate(iter_img(Dl_i)):
        tmap_filename = (result_dir + '/network' + str(i) + 'dl.nii.gz')
        cur_img.to_filename(tmap_filename)

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
        thres_connect_img =  regions.connected_label_regions(result_dir + '/network_mask' + str(i) + 'dl.nii.gz', connect_diag=False, min_size=300)
        extracted_data2 = thres_connect_img.get_data()
        roi_sup0masklab = extracted_data2>0
        extracted_data = nib.load(result_dir + '/network_mask_tresh' + str(i) + 'dl.nii.gz').get_data()
        labelnetwork = np.where(roi_sup0masklab, extracted_data, 0)
        labeled_img2 = nilearn.image.new_img_like(thres_img, labelnetwork, copy_header=True)
        #with stat values
        extracted_data = nib.load(tmap_filename).get_data()
        labelnetwork = np.where(roi_sup0masklab, extracted_data, 0)
        labeled_img2 = nilearn.image.new_img_like(thres_img, labelnetwork, copy_header=True )
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
