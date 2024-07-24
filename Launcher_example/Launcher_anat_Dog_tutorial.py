#import
import os
import subprocess
import glob
import numpy as np
import pandas as pd
import sys

from bids import BIDSLayout
from bids.reports import BIDSReport

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output


MAIN_PATH = opj('/','srv','projects','easymribrain')
# Freesurfer set up
FS_dir    = opj(MAIN_PATH,'FS_Dir_tmp')
sys.path.append(opj(MAIN_PATH,'code','EasyMRI_brain-master'))
import anatomical._0_Pipeline_launcher


################################### if re-use this script auto: ####################################################
##### Trinity Session 6 (T1 or T2) and Unity Session take the other T1. Otherwise: anat image not in same space ####
########################################v###########################################################################

##########################################
########### Subject loader################
##########################################

#https://bids-standard.github.io/pybids/reports/index.html

###where to store the BIDS data?
bids_dir = opj(MAIN_PATH,'data','MRI','Dog','BIDS_k9')
study = 'Dog'

##########################################
########### Subject loader################
##########################################

layout= BIDSLayout(bids_dir,  validate=True)
print(layout)
subject = layout.get_subjects()
print(subject)
tasks = layout.get_tasks()
print(tasks)

###report
report = BIDSReport(layout)
#counter = report.generate()
#main_report = counter.most_common()[0][0]
#print(main_report)


# Ask get() to return the ids of subjects that have T1w files #return_type='filename
T1 = layout.get(return_type='filename', target='subject', suffix='T1w', extension='nii.gz')
print(T1)
###question


# Ask get() to return the ids of subjects that have T2w files
#T2 = layout.get(return_type='filename', target='subject', suffix='T2w', extension='nii.gz')
#print(T2)

# Ask get() to return the ids of subjects that have T1w files
Bold = layout.get(return_type='filename', target='subject', suffix='epi', extension='nii.gz')

# Ask get() to return the ids of subjects that have T1w files
#topup_dir = layout.get(return_type='filename', target='subject', suffix='epi', extension='nii.gz')

# Convert the layout to a pandas dataframe
df = layout.to_df()
df.head()


allinfo_study_c = df[(df['suffix'] == 'epi') & (df['extension'] == '.nii.gz')]

#### Session number
list_of_ones = [1] * len(allinfo_study_c)


allinfo_study_c['Session'] = list_of_ones
allinfo_study_c.rename(columns={'subject': 'ID'}, inplace=True)
allinfo_study_c.rename(columns={'path': 'DICOMdir'}, inplace=True)

folders = glob.glob(opj(bids_dir,'sub*'))

##############################################################  TO DO !! ##############################################################

#### Create a pandas sheet for the dataset (I like it, it help to know what you are about to process
allinfo_study_c = df[(df['suffix'] == 'T1') & (df['extension'] == '.nii.gz')]
list_of_ones = [1] * len(allinfo_study_c)
allinfo_study_c['Session'] = list_of_ones
allinfo_study_c.rename(columns={'subject': 'ID'}, inplace=True)
allinfo_study_c.rename(columns={'path': 'DICOMdir'}, inplace=True)

filter1 = allinfo_study_c["ID"].isin([])

allinfo_study_c_formax = allinfo_study_c.copy()



##############Select all the monkey of the study
### equal to allinfo_study_c, espcially if not longitudinal  and you have not selected specific subjects

#########creat lists of indiv and usefull variables
all_data_path = []
all_ID = []
all_Session = []
all_data_path_max = []
all_ID_max = []
max_session = []
max_sessionlist = []
animal_ID = []

for ID in pd.unique(allinfo_study_c_formax.ID):
    list_session = allinfo_study_c_formax.loc[allinfo_study_c_formax['ID'] == ID].Session.dropna()
    listereverse = list(list_session)
    listereverse.reverse()
    max_session.append(np.array(listereverse).max())

    for Session in listereverse:
        print('session numuber ' + str(Session))

        # Organization of the folders
        data_path = opj(bids_dir,'sub-' + ID,'ses-' + str(Session))
        all_data_path.append(data_path)
        all_Session.append(Session)
        all_ID.append(ID)
        max_sessionlist.append(np.array(listereverse).max())
        animal_ID.append(ID + 'ses-' + str(Session))


for ID, Session in zip(pd.unique(allinfo_study_c_formax.ID), max_session):
    # Organization of the folders
    data_path = opj(bids_dir,'sub-' + ID,'ses-' + str(Session))
    all_data_path_max.append(data_path)
    all_ID_max.append(ID)


######## select animals that have not been analyzed yet
removelist = []
''' 
######### select the indiv you want to analyse!!!
for num, (ID, Session, data_path, max_ses) in enumerate(zip(all_ID, all_Session, all_data_path, max_sessionlist)):
    if ID in []:
        removelist.append(num)
'''
all_ID =  [item for i, item in enumerate(all_ID) if i not in removelist]
all_Session =  [item for i, item in enumerate(all_Session) if i not in removelist]
all_data_path =  [item for i, item in enumerate(all_data_path) if i not in removelist]
max_sessionlist =  [item for i, item in enumerate(max_sessionlist) if i not in removelist]


###### needs to be ID + 'ses-' + str(Session)
deoblique_exeption1 = []
deoblique_exeption2 = animal_ID

##########################################
################To Do ####################
##########################################
###########specific steps#################
##########################################
coregistration_longitudinal = False #YES or NO
check_visualy_each_img = False #YES or NO
do_manual_crop = False #YES or NO

#### only use this option if T1 and T2 are in the same space!!!!!
IgotbothT1T2 = False #YES or NO

check_visualy_final_mask = False #YES or NO
deoblique='WARP' #header or WARP
n_for_ANTS='hammingWindowedSinc'
overwrite_option = True #YES or NO

####Choose to normalize using T1 or T2
type_norm = 'T2w' # T1 or T2
otheranat = '' #NA if none
###masking
#ruf XXX!!!!!!
###img use for masking in Skullstrip 1 'maybe this need to be change'!!!!!! because Skullstrip 2 is in auto equal to type_norm.... not sure that it will not creat problem in the futur
masking_img = 'T2w'

brain_skullstrip_1 ='3dSkullStrip_Rat' # bet2_ANTS or MachinL

#precise
brain_skullstrip_2 ='Custum_QWARPT2' # bet2_ANTS or MachinL

useT1T2_for_coregis = False

do_fMRImasks = True
Align_img_to_template = '@Align_Centers' #3dAllineate or No or @Align_Centers
cost3dAllineate = 'lpa'
'''
ls   *OR*  leastsq         = Least Squares [Pearson Correlation]
mi   *OR*  mutualinfo      = Mutual Information [H(b)+H(s)-H(b,s)]
crM  *OR*  corratio_mul    = Correlation Ratio (Symmetrized*)
nmi  *OR*  norm_mutualinfo = Normalized MI [H(b,s)/(H(b)+H(s))]
hel  *OR*  hellinger       = Hellinger metric
crA  *OR*  corratio_add    = Correlation Ratio (Symmetrized+)
crU  *OR*  corratio_uns    = Correlation Ratio (Unsym)
lpc  *OR*  localPcorSigned = Local Pearson Correlation Signed
lpa  *OR*  localPcorAbs    = Local Pearson Correlation Abs
lpc+ *OR*  localPcor+Others= Local Pearson Signed + Others
lpa+ *OR*  localPcorAbs+Others= Local Pearson Abs + Others
'''
#######################################################################
######################### study template??? ###########################
#######################################################################

#creat_study_template with type_norm img
creat_study_template = False

#folder where you want to store the stdy template
study_template_atlas_forlder = ''
#then
dir_out = study_template_atlas_forlder + ''

#do you want to use all the data or only the last one of each subject (for longitud inal co-registration)
which_on = '' # all or max

###use type_norm or otheranat for atlas template to study template co-registration
Atemplate_to_Stemplate = ''

template_skullstrip = ''

stdy_template_mask = ''
stdy_template = ''

do_surfacewith = 'T2w' #'T1' 'T1andT2'

    ##########################################
    ###########define orientation#############
    ##########################################

###question
####WITH deoblique='WARP'
##orig LPI
orientation = 'RIP'
####WITH deoblique='header'
#orientation = 'LSP'

#if BIDStype == 1:
#    list_anat = sorted(glob.glob(opj(path_anat, 'sub-' + ID + '_ses-' + str(Session) + '_run-*' + Timage + '.nii.gz')))
#if BIDStype == 2:
#    list_anat = sorted(glob.glob(opj(path_anat, 'sub-' + ID + '_' + Timage + '.nii.gz')))
BIDStype = 2

####### attention!! change LPS -r based on what you can observe =====  RAI => LSP ;  LIP => LAS ; LSP => LPS (xxxchange RL?)
#RAI
#first
#second LR
#third

### file for standardize space
FS_buckner40_TIF = opj(FS_dir,'MacaqueYerkes19')
FS_buckner40_GCS = opj(FS_dir,'MacaqueYerkes19')

    ##########################################
    ##### define atlases and tempates ########
    ##########################################
diratlas_orig = opj(MAIN_PATH,'data','Atlas','13_Atlas_project','New_atlas_Dual','RatWHS') # sting # sting
list_atlases = [opj(diratlas_orig, 'atlaslvl1.nii.gz'),
opj(diratlas_orig, 'atlaslvl2.nii.gz'),
opj(diratlas_orig, 'atlaslvl3.nii.gz'),
opj(diratlas_orig, 'atlaslvl4.nii.gz'),
opj(diratlas_orig, 'atlaslvl1_LR.nii.gz'),
opj(diratlas_orig, 'atlaslvl2_LR.nii.gz'),
opj(diratlas_orig, 'atlaslvl3_LR.nii.gz'),
opj(diratlas_orig, 'atlaslvl4_LR.nii.gz')]

BASE_SS     = opj(MAIN_PATH,'data','Atlas','13_Atlas_project','0_Atlas_modify','Atlas','RatWHS','templateT2.nii.gz') # sting
BASE_mask   = opj(MAIN_PATH,'data','Atlas','13_Atlas_project','New_atlas_Dual','RatWHS','BrainMask.nii.gz') # sting
BASE_bet =''
GM_mask     = opj(MAIN_PATH,'data','Atlas','13_Atlas_project','New_atlas_Dual','RatWHS','CorticalMask.nii.gz') # sting

####atlases files
Aseg_ref    = opj(MAIN_PATH,'data','Atlas','13_Atlas_project','New_atlas_Dual','RatWHS','atlas_forSEG_final.nii.gz')
Aseg_refLR  = opj(MAIN_PATH,'data','Atlas','13_Atlas_project','New_atlas_Dual','RatWHS','atlas_forSEG_final_LR.nii.gz')

###question of Aseg_refLR
define_center = '0'


### if it doesn't exists let's make it!!! but you need an Aseg_ref!!
if not ope(Aseg_refLR):
    command = '3dcalc -a ' + Aseg_ref + ' -expr "step(ispositive(x-' + define_center + ')*a)" -prefix ' + opd(Aseg_refLR) + '/Aseg_ref_L.nii.gz'
    spco([command], shell=True)

    command = '3dcalc -a ' + opd(Aseg_refLR) + '/Aseg_ref_L.nii.gz -b ' + Aseg_ref + ' -expr "ifelse(a, a*255,step(b)*127)" -prefix ' + Aseg_refLR
    spco([command], shell=True)



#### for 14 ####
list_atlases_2 = [opj(diratlas_orig, 'atlaslvl1.nii.gz'),
opj(diratlas_orig, 'atlaslvl2.nii.gz'),
opj(diratlas_orig, 'atlaslvl3.nii.gz'),
opj(diratlas_orig, 'atlaslvl4.nii.gz')]

FreeSlabel_ctab_list = [opj(MAIN_PATH,'data','Atlas','13_Atlas_project','LUT_files','Multispecies_LUT.txt'),
opj(MAIN_PATH,'data','Atlas','13_Atlas_project','LUT_files','Multispecies_LUT.txt'),
opj(MAIN_PATH,'data','Atlas','13_Atlas_project','LUT_files','Multispecies_LUT.txt'),
opj(MAIN_PATH,'data','Atlas','13_Atlas_project','LUT_files','Multispecies_LUT.txt')]


######### define other usefull paramater automatically (do no touch)#########
Hmin     = ['l','r']
### Block1: step 1,2,3
### Block2: step 4,5
### Block3: step 6 (template)
### Block4: step 7,8 (altases, masks, fmri masks)
### Block5: step 9, 10, 11, 12, 13, 14, 15 (surfaces)

Skip_step = [10,11,12,13,14,15]

Lut_file = opj(MAIN_PATH,'data','Atlas','13_Atlas_project','LUT_files','Multispecies_LUT.txt')
species = 'RatWHS'


anatomical._0_Pipeline_launcher.preprocess_anat(BIDStype, deoblique_exeption1, deoblique_exeption2, deoblique, BASE_mask, coregistration_longitudinal, creat_study_template,
    orientation, masking_img, brain_skullstrip_1, brain_skullstrip_2, n_for_ANTS, Skip_step, check_visualy_each_img, do_manual_crop, do_fMRImasks,
    BASE_SS,BASE_bet, which_on, all_ID_max, max_session, all_data_path_max, all_ID, all_Session, all_data_path, study_template_atlas_forlder, template_skullstrip,
    IgotbothT1T2, list_atlases, Aseg_ref, Aseg_refLR, dir_out, FS_dir, diratlas_orig, do_surfacewith, Atemplate_to_Stemplate,
    FS_buckner40_TIF,FS_buckner40_GCS, Hmin, Lut_file, otheranat, type_norm, max_sessionlist, bids_dir, check_visualy_final_mask, useT1T2_for_coregis, FreeSlabel_ctab_list, list_atlases_2, cost3dAllineate, Align_img_to_template,
    species, overwrite_option,MAIN_PATH)
