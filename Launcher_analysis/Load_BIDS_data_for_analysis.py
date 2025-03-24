import os
import subprocess
import numpy as np
import pandas as pd
import glob
import nibabel as nib
from fonctions.extract_filename import extract_filename
#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput


def load_data(bids_dir, df, ntimepoint_treshold, list_to_keep, list_to_remove, endfmri):
    allinfo_study_c = df[(df['suffix'] == 'bold') & (df['extension'] == '.nii.gz')]
    allinfo_study_c.rename(columns={'session': 'Session'}, inplace=True)
    allinfo_study_c.rename(columns={'subject': 'ID'}, inplace=True)
    allinfo_study_c.rename(columns={'path': 'DICOMdir'}, inplace=True)
    allinfo_study_c_formax = allinfo_study_c.copy()

    ############################################################## NOTHING TO DO HERE ##############################################################

    #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
    #### #### now, we need to creat the list of variable to feed to function "_0_Pipeline_launcher" #### ####
    #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

    #########creat lists of indiv and usefull variables
    # just a list of all datapath
    all_data_path = []
    # just a list of all subject ID
    all_ID = []
    # just a list of all Session
    all_Session = []

    ### for longitudinal datasets, you may want to know which session is the last, this is the way to handle it
    all_data_path_max = []
    all_ID_max = []
    max_session = []
    max_sessionlist = []

    ### this is just to list all subject in a ID + 'ses-' + str(Session) way (usefull for the variable "deoblique_exeption")
    animal_ID = []

    #let's add all the the string to those lists
    for ID in pd.unique(allinfo_study_c_formax.ID):
        list_session = allinfo_study_c_formax.loc[allinfo_study_c_formax['ID'] == ID].Session.dropna()
        listereverse = list(list_session)
        listereverse.reverse()
        if len(list(list_session))>1:
            max_session.append(np.array(listereverse).max())
        else:
            max_session.append(np.array(listereverse))

        for Session in listereverse:
            print('session numuber ' + str(Session))

            # Organization of the folders
            data_path = opj(bids_dir,'sub-' + ID,'ses-' + str(Session))
            all_data_path.append(data_path)
            all_Session.append(Session)
            all_ID.append(ID)
            if len(list(list_session)) > 1:
                max_sessionlist.append(np.array(listereverse).max())
            else:
                max_sessionlist.append(np.array(listereverse))
            animal_ID.append(ID + 'ses-' + str(Session))

    for ID, Session in zip(pd.unique(allinfo_study_c_formax.ID), max_session):
        # Organization of the folders
        data_path = opj(bids_dir,'sub-' + ID,'ses-' + str(Session))
        all_data_path_max.append(data_path)
        all_ID_max.append(ID)

    ##############################################################  TO DO !!!! ##############################################################

    ######## sometime you aleady analysed some subjects, and you don't want to do it again,
    # this can not be done in the previous step, as if you are working with the longitudinal dataset it can mess up your analysis
    ######## select animals that have been analyzed already

    if len(list_to_keep) >0 and len(list_to_remove)==0:
        removelist = []
        ######### select the indiv you want to remove !!!
        for num, (ID, Session, data_path, max_ses) in enumerate(zip(all_ID, all_Session, all_data_path, max_sessionlist)):
            if ID in list_to_remove:
                removelist.append(num)

        #### apply this to the already created lists
        all_ID =  [item for i, item in enumerate(all_ID) if i in removelist]
        all_Session =  [item for i, item in enumerate(all_Session) if i in removelist]
        all_data_path =  [item for i, item in enumerate(all_data_path) if i in removelist]
        max_sessionlist =  [item for i, item in enumerate(max_sessionlist) if i in removelist]


    elif len(list_to_keep) == 0 and len(list_to_remove) > 0:
        removelist = []
        ######### select the indiv you want to remove !!!
        for num, (ID, Session, data_path, max_ses) in enumerate(zip(all_ID, all_Session, all_data_path, max_sessionlist)):
            if ID in list_to_remove:
                removelist.append(num)

        #### apply this to the already created lists
        all_ID =  [item for i, item in enumerate(all_ID) if i not in removelist]
        all_Session =  [item for i, item in enumerate(all_Session) if i not in removelist]
        all_data_path =  [item for i, item in enumerate(all_data_path) if i not in removelist]
        max_sessionlist =  [item for i, item in enumerate(max_sessionlist) if i not in removelist]

    elif len(list_to_keep) == 0 and len(list_to_remove) == 0:
        print('all subjects are included for the analysis')

    else:
        raise ValueError('one of list_to_keep or list_to_remove need to be empty')


    images_dir = []
    mean_imgs = []
    seed_to_voxel_correlations_all_fish = []
    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, max_sessionlist):

        dir_fMRI_Refth_RS = opj(data_path, 'func')
        dir_fMRI_Refth_RS_prepro = opj(dir_fMRI_Refth_RS, '01_prepro')
        dir_fMRI_Refth_RS_prepro1 = opj(dir_fMRI_Refth_RS_prepro, '01_funcspace')
        dir_fMRI_Refth_RS_prepro2 = opj(dir_fMRI_Refth_RS_prepro, '02_anatspace')
        dir_fMRI_Refth_RS_prepro3 = opj(dir_fMRI_Refth_RS_prepro, '03_atlas_space')

        # Check the func runs
        list_RS = sorted(glob.glob(opj(dir_fMRI_Refth_RS, endfmri)))
        RS = [os.path.basename(i) for i in list_RS]

        if len(list_RS) == 0:
            nl = 'ERROR : No func image found, we are look for an image define such as opj(dir_fMRI_Refth_RS, endfmri) and here it is ' + str(
                opj(dir_fMRI_Refth_RS, endfmri)) + ' I would check how you define "endfmri"'
            raise ValueError(nl)
        list_RS_list = list_RS.copy()
        list_pop_index = []

        for imageF in list_RS_list:
            # Load the fMRI NIfTI image
            fmri_image = nib.load(imageF)
            # Get the shape of the image (x, y, z, t)
            image_shape = fmri_image.shape
            # Check the number of time points (4th dimension)
            if len(image_shape) == 4:
                ntimepoint = image_shape[3]  # The 4th dimension represents time

                if int(ntimepoint) < ntimepoint_treshold:
                    index_of_imageF = list_RS.index(imageF)
                    list_RS.pop(index_of_imageF)
                    list_pop_index.append(index_of_imageF)

        nb_run = len(list_RS)
        # Setup for distortion correction using Fieldmaps

        if ope(opj(dir_fMRI_Refth_RS_prepro3, 'Mean_Image_RcT_SS_in_template.nii.gz')):
            mean_imgs.append(opj(dir_fMRI_Refth_RS_prepro3, 'Mean_Image_RcT_SS_in_template.nii.gz'))

        for i in range(0, int(nb_run)):
            root_RS = extract_filename(RS[i])
            ### here we want to select all the image preTTT in a common space
            if ope(opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz')):
                images_dir.append(opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz'))

        templatelow = opj(dir_fMRI_Refth_RS_prepro3, 'BASE_SS_fMRI.nii.gz')

    print("Remaining all_data_path:", all_data_path)
    return(images_dir, all_ID, all_Session, all_data_path, max_sessionlist, mean_imgs, templatelow)