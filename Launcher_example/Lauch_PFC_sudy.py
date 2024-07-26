#import
import os
import subprocess
import numpy as np
import pandas as pd
import sys


#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output


MAIN_PATH   = opj('/','srv','projects','easymribrain')

SINGULARITY = 'singularity -- bind ' + opj(MAIN_PATH)
s_path      =  opj(MAIN_PATH + 'code' + 'singularity')
AFNI_sif    =  ' ' + opj(s_path + 'afni_make_build_AFNI_23.1.10.sif ')
FSL_sif     =  ' ' + opj(s_path + 'fsl_6.0.5.1-cuda9.1.sif ')
FS_sif      =  ' ' + opj(s_path + 'freesurfer_7.4.1.sif ')
WB_sif      =  ' ' + opj(s_path + 'connectome_workbench_1.5.0-freesurfer-update.sif ')
ITKs        =  ' ' + opj(s_path + 'itksnap_5.0.9.sif ')


################################### if re-use this script auto: ####################################################
##### Trinity Session 6 (T1 or T2) and Unity Session take the other T1. Otherwise: anat image not in same space ####
########################################v###########################################################################

##########################################
########### Subject loader################
##########################################

#https://bids-standard.github.io/pybids/reports/index.html

sys.path.append(opj(MAIN_PATH + 'Code' + 'EasyMRI_Brain-master'))     # to check
import anatomical._0_Pipeline_launcher

from sammba import io_conversions, registration

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
ops = os.path.splitext
spco = subprocess.check_output

# Freesurfer set up
FS_dir = opj(MAIN_PATH + 'FS_Dir_tmp')

command = 'export SUBJECTS_DIR=' + FS_dir
spco([command], shell=True)


###where to store the BIDS data?
study = 'Macaque'
bids_dir = opj('/srv/projects/easymribrain/Data'+ study + 'BIDS_NIH') # to check


##########################################
########### Subject loader################
##########################################


mldir = opn('/home/cgarin/Documents/1_Macaque_MRI/')                # to check
ml_excel = opn(opj(mldir, '5_Data_launcher/Garin_macaque.xlsx'))    # to check

animalinfo = pd.read_excel(ml_excel, sheet_name='animalinfo')
MRIsessions = pd.read_excel(ml_excel, sheet_name='MRIsessions')

#### create column of PV folder names

PVdirs_base = [opb(opn(PVdir_base)) for PVdir_base in MRIsessions.Pvdir]
PVdirs_base = pd.DataFrame({'PVdir_base': PVdirs_base})


#### create column of PV folder names
# concatenate above two columns to MRIsessions, merge with animalinfo and
# calculate age at scan
MRIsessions = pd.concat([MRIsessions, PVdirs_base], axis=1)
allinfo = pd.merge(animalinfo, MRIsessions, on='ID')
age = allinfo.DOS - allinfo.DOB
age = age.dt.days

allinfo = pd.concat([allinfo,
                     pd.DataFrame({'age': age})],
                     axis=1)

'''
##############Select the "good" monkey for the study
filter1 = allinfo["Study"].isin(["pFC_AS", "pFC_ODR"]) & allinfo["Trained_yung"].isin(["YES"])
filter2 = allinfo["ID"].isin(["Trinity"]) & allinfo["Other_sequence"].isin(["SE-fMRI_LR_1.9_150vol_SE-fMRI_LR_1.9_150vol_10311296.nii"])
#filter6 = allinfo["ID"].isin(["Quantum"]) & allinfo["Session"].isin([2])
#filter7 = allinfo["ID"].isin(["Pickle"]) & allinfo["Session"].isin([2])
filter3 = allinfo["Session"].isin([1])
filter4 = allinfo["Session"].isin([12])
filter5 = allinfo["Session"].isin([13])

allinfo_study_c = allinfo[filter1]
allinfo_study_c = allinfo_study_c[~filter2]
allinfo_study_c = allinfo_study_c[~filter3]
allinfo_study_c = allinfo_study_c[~filter4]
allinfo_study_c = allinfo_study_c[~filter5]
#allinfo_study_c = allinfo_study_c[~filter6]
#allinfo_study_c = allinfo_study_c[~filter7]

'''

##############Select all the monkey of the study 
filter1 = allinfo["Study"].isin(["pFC_AS", "pFC_ODR"]) & allinfo["Trained_yung"].isin(["YES"])
filter2 = allinfo["ID"].isin(["Trinity"]) & allinfo["Other_sequence"].isin(["SE-fMRI_LR_1.9_150vol_SE-fMRI_LR_1.9_150vol_10311296.nii"])
filter4 = allinfo["Session"].isin([12])
filter5 = allinfo["Session"].isin([13])
filter6 = allinfo["ID"].isin(["Quantum", "Oliver"])
filter7 = allinfo["Session"].isin([95])

allinfo_study_c_formax = allinfo[filter1]
allinfo_study_c_formax = allinfo_study_c_formax[~filter2]
allinfo_study_c_formax = allinfo_study_c_formax[~filter4]
allinfo_study_c_formax = allinfo_study_c_formax[~filter5]
allinfo_study_c_formax = allinfo_study_c_formax[filter6]
allinfo_study_c_formax = allinfo_study_c_formax[~filter7]

#########creat lists of indiv and usefull variables
all_data_path = []
all_ID = []
all_Session = []
all_data_path_max = []
all_ID_max = []
max_session = []
max_sessionlist = []

for ID in pd.unique(allinfo_study_c_formax.ID):
    list_session = allinfo_study_c_formax.loc[allinfo_study_c_formax['ID'] == ID].Session.dropna()
    listereverse = list(list_session)
    listereverse.reverse()
    if 95 in listereverse:
        listereverse.remove(95)

    max_session.append(np.array(listereverse).max())

    for Session in listereverse:
        print('session numuber ' + str(Session))

        # Organization of the folders
        data_path = opj(bids_dir,'sub-' + ID,'ses-' + str(Session))
        all_data_path.append(data_path)
        all_Session.append(Session)
        all_ID.append(ID)
        max_sessionlist.append(np.array(listereverse).max())


for ID, Session in zip(pd.unique(allinfo_study_c_formax.ID), max_session):
    # Organization of the folders
    data_path = opj(bids_dir,'sub-' + ID,'ses-' + str(Session))
    all_data_path_max.append(data_path)
    all_ID_max.append(ID)


######## select animals that have not been analyzed yet
removelist = []


######### select the indiv you want to analyse!!!
for num, (ID, Session, data_path, max_ses) in enumerate(zip(all_ID, all_Session, all_data_path, max_sessionlist)):
    if ID in ["Quantum"] and Session in [1,2,3,4,5,6]:
        removelist.append(num)
    elif ID in ["Pickle"] and Session in [1,2,3,4,5,6,7,8,9,10]:
        removelist.append(num)
    elif ID in ["Sonic"] and Session in [1,2,3,4,5,6,7,8,9,10]:
        removelist.append(num)
    elif ID in ["Oliver"] and Session in [1,2,3,4,5,6,7,8,9,10]:
        removelist.append(num)
    elif ID in ["Trinity"] and Session in [1,2,3,4,5,6,7,8,9,10]:
        removelist.append(num)
    elif ID in ["Unity"] and Session in [1,2,3,4,5,6,7,8,9,10]:
        removelist.append(num)
    elif ID in ["Viking"] and Session in [1,2,3,4,5,6,7,8,9,10]:
        removelist.append(num)
    elif ID in ["Roshan"] and Session in [1,2,3,4,5,6,7,8,9,10]:
        removelist.append(num)
'''
######### select the indiv you want to analyse!!!
for num, (ID, Session, data_path, max_ses) in enumerate(zip(all_ID, all_Session, all_data_path, max_sessionlist)):
    if ID in ["Roshan", "Pickle", "Oliver", "Quantum", "Unity", "Viking", "Sonic"]:
        removelist.append(num)
    elif ID in ["Trinity"]:
        removelist.append(num)
'''



all_ID =  [item for i, item in enumerate(all_ID) if i not in removelist]
all_Session =  [item for i, item in enumerate(all_Session) if i not in removelist]
all_data_path =  [item for i, item in enumerate(all_data_path) if i not in removelist]
max_sessionlist =  [item for i, item in enumerate(max_sessionlist) if i not in removelist]


deoblique_exeption1 = ['Roshanses-9']
deoblique_exeption2 = ['Trinityses-6', 'Unityses-5']

#if BIDStype == 1:
#    list_anat = sorted(glob.glob(opj(path_anat, 'sub-' + ID + '_ses-' + str(Session) + '_run-*' + Timage + '.nii.gz')))
#if BIDStype == 2:
#    list_anat = sorted(glob.glob(opj(path_anat, 'sub-' + ID + '_' + Timage + '.nii.gz')))
BIDStype = 1

##########################################
################To Do ####################
##########################################
###########specific steps#################
##########################################
coregistration_longitudinal = True #YES or NO
check_visualy_each_img = False #YES or NO
do_manual_crop = False #YES or NO

#### only use this option if T1 and T2 are in the same space!!!!! 
IgotbothT1T2 = True #YES or NO

check_visualy_final_mask = False #YES or NO
deoblique='WARP' #header or WARP
n_for_ANTS='HammingWindowedSinc'
overwrite_option = True #YES or NO

####Choose to normalize using T1 or T2
type_norm = 'T1' # T1 or T2
otheranat = 'T2FLAIR' #NA if none
###masking
#ruf XXX!!!!!!
###img use for masking in Skullstrip 1 'maybe this need to be change'!!!!!! because Skullstrip 2 is in auto equal to type_norm.... not sure that it will not creat problem in the futur
masking_img = 'T1'

brain_skullstrip_1 ='Custum_Macaque2' # bet2_ANTS or MachinL

#precise
brain_skullstrip_2 ='Custum_QWARP' # bet2_ANTS or MachinL

useT1T2_for_coregis = False

do_fMRImasks = True
#######################################################################
######################### study template??? ###########################
#######################################################################

#creat_sutdy_template with type_norm img
creat_sutdy_template = True

#folder where you want to store the stdy template
study_template_atlas_forlder = opj(MAIN_PATH + study + 'Study_template_test')
#then 
dir_out = study_template_atlas_forlder + '/Atlases_ref_in_stdy_template/'

#do you want to use all the data or only the last one of each subject (for longitudinal co-registration)
which_on = 'max' # all or max

###use type_norm or otheranat for atlas template to study template co-registration
Atemplate_to_Stemplate = 'T1'

template_skullstrip = 'MachinL'

stdy_template_mask = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template_mask.nii.gz')
stdy_template = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template.nii.gz')

do_surfacewith = 'T1' #'T1' 'T1andT2'

    ##########################################
    ###########define orientation#############
    ##########################################

###question
####WITH deoblique='WARP'
##orig LPI
orientation = 'RAI'
####WITH deoblique='header'
#orientation = 'LSP'


####### attention!! change LPS -r based on what you can observe =====  RAI => LSP ;  LIP => LAS ; LSP => LPS (xxxchange RL?)
#RAI
#first 
#second LR
#third

### file for standardize space
FS_buckner40 = opj('/home/cgarin/Documents/8_Multispecies/fMRI/')

    ##########################################
    ##### define atlases and tempates ########
    ##########################################
diratlas_orig = '/home/cgarin/Documents/1_Macaque_MRI/2_Atlas/NMT_v2_modif/NTM2V6/'
list_atlases = [diratlas_orig + 'atlaslvl1.nii.gz',
diratlas_orig + 'atlaslvl2.nii.gz',
diratlas_orig + 'atlaslvl3.nii.gz',
diratlas_orig + 'atlaslvl4.nii.gz',
diratlas_orig + 'CHARM_1_in_NMT_v2.0_sym.nii.gz',
diratlas_orig + 'CHARM_2_in_NMT_v2.0_sym.nii.gz',
diratlas_orig + 'CHARM_3_in_NMT_v2.0_sym.nii.gz',
diratlas_orig + 'CHARM_4_in_NMT_v2.0_sym.nii.gz',
diratlas_orig + 'CHARM_5_in_NMT_v2.0_sym.nii.gz',
diratlas_orig + 'CHARM_6_in_NMT_v2.0_sym.nii.gz',
diratlas_orig + 'atlaslvl1_LR.nii.gz',
diratlas_orig + 'atlaslvl2_LR.nii.gz',
diratlas_orig + 'atlaslvl3_LR.nii.gz',
diratlas_orig + 'atlaslvl4_LR.nii.gz',
diratlas_orig + 'CHARM_1_in_NMT_v2.0_sym_LR.nii.gz',
diratlas_orig + 'CHARM_2_in_NMT_v2.0_sym_LR.nii.gz',
diratlas_orig + 'CHARM_3_in_NMT_v2.0_sym_LR.nii.gz',
diratlas_orig + 'CHARM_4_in_NMT_v2.0_sym_LR.nii.gz',
diratlas_orig + 'CHARM_5_in_NMT_v2.0_sym_LR.nii.gz',
diratlas_orig + 'CHARM_6_in_NMT_v2.0_sym_LR.nii.gz',
diratlas_orig + 'atlas_forSEG_final_LR.nii.gz']

BASE_SS     = diratlas_orig + 'NMT_SS_v6.nii.gz'
BASE_mask   = diratlas_orig + 'NMT_brainmask_v6.nii.gz'

####atlases files
Aseg_ref   = diratlas_orig + 'Yerkes_aseg_final.nii.gz'
Aseg_refLR = diratlas_orig + 'atlas_forSEG_final_LR.nii.gz'

###question of Aseg_refLR
define_center = '0'

### if it doesn't exists let's make it!!! but you need an Aseg_ref!!
if not ope(Aseg_refLR):
    command = '3dcalc -a ' + Aseg_ref + ' -expr "step(ispositive(x-' + define_center + ')*a)" -prefix ' + opd(Aseg_refLR) + '/Aseg_ref_L.nii.gz'
    spco([command], shell=True)

    command = '3dcalc -a ' + opd(Aseg_refLR) + '/Aseg_ref_L.nii.gz -b ' + Aseg_ref + ' -expr "ifelse(a, a*255,step(b)*127)" -prefix ' + Aseg_refLR
    spco([command], shell=True)


cost3dAllineate = 'hel'

#### for 14 ####
list_atlases_2 = [diratlas_orig + 'aseg.nii.gz',
diratlas_orig + 'atlaslvl1.nii.gz',
diratlas_orig + 'atlaslvl2.nii.gz',
diratlas_orig + 'atlaslvl3.nii.gz',
diratlas_orig + 'atlaslvl4.nii.gz',
diratlas_orig + 'CHARM_1_in_NMT_v2.0_sym.nii.gz',
diratlas_orig + 'CHARM_2_in_NMT_v2.0_sym.nii.gz',
diratlas_orig + 'CHARM_3_in_NMT_v2.0_sym.nii.gz',
diratlas_orig + 'CHARM_4_in_NMT_v2.0_sym.nii.gz',
diratlas_orig + 'CHARM_5_in_NMT_v2.0_sym.nii.gz',
diratlas_orig + 'CHARM_6_in_NMT_v2.0_sym.nii.gz']

FreeSlabel_ctab_list = ['/home/cgarin/Documents/8_Multispecies/13_Atlas_project/LUT_files/FreeSurferColorLUT.txt',
'/home/cgarin/Documents/8_Multispecies/13_Atlas_project/LUT_files/Multispecies_LUT.txt',
'/home/cgarin/Documents/8_Multispecies/13_Atlas_project/LUT_files/Multispecies_LUT.txt',
'/home/cgarin/Documents/8_Multispecies/13_Atlas_project/LUT_files/Multispecies_LUT.txt',
'/home/cgarin/Documents/8_Multispecies/13_Atlas_project/LUT_files/Multispecies_LUT.txt',
'/home/cgarin/Documents/8_Multispecies/13_Atlas_project/LUT_files/CHARMSARM/CHARMSARM_1_StatsLUT.txt',
'/home/cgarin/Documents/8_Multispecies/13_Atlas_project/LUT_files/CHARMSARM/CHARMSARM_2_StatsLUT.txt',
'/home/cgarin/Documents/8_Multispecies/13_Atlas_project/LUT_files/CHARMSARM/CHARMSARM_3_StatsLUT.txt',
'/home/cgarin/Documents/8_Multispecies/13_Atlas_project/LUT_files/CHARMSARM/CHARMSARM_4_StatsLUT.txt',
'/home/cgarin/Documents/8_Multispecies/13_Atlas_project/LUT_files/CHARMSARM/CHARMSARM_5_StatsLUT.txt',
'/home/cgarin/Documents/8_Multispecies/13_Atlas_project/LUT_files/CHARMSARM/CHARMSARM_6_StatsLUT.txt']


######### define other usefull paramater automatically (do no touch)#########

Hmin     = ['l','r']

Skip_step = [1,2,4,5,6,7,8,9,10,11,12,13,14,15]

Lut_file = '/home/cgarin/Documents/8_Multispecies/13_Atlas_project/LUT_files/AsegStatsLUT.txt'

species = 'Macaque'

anatomical._0_Pipeline_launcher.preprocess_anat(BIDStype, deoblique_exeption1, deoblique_exeption2, deoblique, BASE_mask, coregistration_longitudinal, creat_study_template,
    orientation, masking_img, brain_skullstrip_1, brain_skullstrip_2, n_for_ANTS, Skip_step, check_visualy_each_img, do_manual_crop, do_fMRImasks,
    BASE_SS,BASE_bet, which_on, all_ID_max, max_session, all_data_path_max, all_ID, all_Session, all_data_path, study_template_atlas_forlder, template_skullstrip,
    IgotbothT1T2, list_atlases, Aseg_ref, Aseg_refLR, dir_out, FS_dir, do_surfacewith, Atemplate_to_Stemplate,
    FS_buckner40_TIF,FS_buckner40_GCS, Hmin, Lut_file, otheranat, type_norm, max_sessionlist, bids_dir, check_visualy_final_mask, useT1T2_for_coregis, FreeSlabel_ctab_list, list_atlases_2, cost3dAllineate, Align_img_to_template,
    species, overwrite_option,MAIN_PATH)