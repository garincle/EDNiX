#import
import os
import subprocess
import glob
import json
import shutil
import numpy as np
import time
import datetime
import ants
import pydicom as pdcm


opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile
spco = subprocess.check_output
spr  = subprocess.run
spgo = subprocess.getoutput


def orga(study,Subject_Dir,Sess_to_process,situation,
         reference,tracor,
         Bq_inj,
         weight,
         REF_name):
  
  # Where to find whatever you need 
  if situation == 'home':
    Major_path = opj('/','media','simon','LaCie','travail','simon','current_works','UCBL','rsFMRI_Tremblay')    
    path_ATLAS = opj('/','media','simon','LaCie','travail','simon','fMRI','MRI_atlas','atlas','macaque')
  
  elif situation == 'work':
    Major_path  = opj('/','home','clavagni','Documents','UCBL','rsFMRI_Tremblay')
    path_ATLAS  = opj('/','home','clavagni','Bureau','MRI_atlas','atlas','macaque')
  
  elif situation == 'gpu-05-01':
    Major_path  = opj('/','srv','projects')
    path_ATLAS  = opj(Major_path,'tremblaybgtemplate','Atlas','NHP','macaque')
    
    # singularity set up
    BIND = ' --bind ' + opj('/','scratch','in_process/') + ',' + opj(Major_path,study) + ',' + path_ATLAS
    s_path      =  opj(Major_path,study,'BIDS','code','singularity')
    AFNI_sif    =  ' ' + opj(s_path , 'afni_make_build_AFNI_23.1.10.sif') + ' '
    FSL_sif     =  ' ' + opj(s_path , 'fsl_6.0.5.1-cuda9.1.sif') + ' '
    FS_sif      =  ' ' + opj(s_path , 'freesurfer_NHP.sif') + ' '
    ITKs        =  ' ' + opj(s_path , 'itksnap_5.0.9.sif') + ' '
    
  
  ## Organization 
  
  # references 
  # need to be simplified
  if reference == 'ballanger':
    REF_fasci    = opj(path_ATLAS,'fascicularis','Wb','Native')
    ref_template = opj(REF_fasci,'volumes','fascicularis_T1_2022_brain.nii.gz')
    ref_mask     = opj(REF_fasci,'volumes','masks','fascicularis_2022_brainmask.nii.gz')
    gray         = opj(REF_fasci,'volumes','masks','Gray_3.nii.gz')
    ref_Wmask    = opj(REF_fasci,'volumes','masks','wmask_erod.nii.gz')
    ref_Vmask    = opj(REF_fasci,'volumes','masks','Vmask_erod.nii.gz')
    ref_csfmask  = opj(REF_fasci,'volumes','masks','ant_csf_mask.nii.gz')
    atlas        = opj(REF_fasci,'volumes','labels','atlasFascicularis_2024.nii.gz')
    atlas_ctab   = opj(REF_fasci,'volumes','labels','atlasFascicularis_2024_StatsLUT.txt')
    
    template_name   = 'BAL2013'
    template_folder = 'fascicularis'
    
    if REF_name == 'Cer':
      REF       = [29,30]
    
    
  elif reference == 'MY19':
    REF_MY19     = opj(path_ATLAS,'Yerkes','macaque19')
    ref_template = opj(REF_MY19,'atlas','MY19_T1w_0.5mm_brain_modif.nii.gz')
    ref_mask     = opj(REF_MY19,'atlas','NMT_mask2.nii.gz')
    gray         = opj(REF_MY19,'atlas','Gray.nii.gz')
    ref_Wmask    = opj(REF_MY19,'atlas','MY19_Wmask_erod.nii.gz')
    ref_Vmask    = opj(REF_MY19,'atlas','MY19_Vmask_erod.nii.gz')
    ref_csfmask  = opj(REF_MY19,'atlas','volumes','masks','ant_csf_mask.nii.gz') # need to do it
    atlas        = opj(REF_MY19,'atlas','D99_atlasv2_Y19.nii.gz')
    atlas_ctab   = opj(path_ATLAS,'NMT_v2_modif','atlas','D99v2_StatsLUT.txt')
    template_name   = 'MY19'
    template_folder = 'MY19'
    
  elif reference  == 'MEBRAINS':
    template_name = 'MEBRAINS'
    print('not ready yet')
  
  elif reference  == 'NMT':
    template_name = 'NMT'
    print('not ready yet')
  


  # Datas

  N = Subject_Dir.split('-')
  Subname = N[0] + '-' + N[1]

  Sess = Sess_to_process.split('-')
  type_sess = Sess[1].lower()
  loca_sess = Sess[1].split('+') 
  
  pet_dir = opj(Major_path,study,'derivatives', Subject_Dir,Sess_to_process,'pet')
  
  if study == 'eating_circuit':
    
    data_path_raw   = opj(Major_path, study, 'sourcedata', Subject_Dir, Sess_to_process)
    
    if type_ses == 'control':
      cond ='_CTL'
      pet_raw_dcm = opj(Major_path,study,'sourcedata' ,'dicoms',N[1] + '_' + N[2], type_ses ,tracor + '_new_pet', N[1] + '_' + S[3] + '_' + S[1] + cond)
    elif type_ses == 'flx':
      cond ='_FLX'
      pet_raw_dcm = opj(Major_path,study,'sourcedata' ,'dicoms',N[1] + '_' + N[2], S[2] + '_' + S[3]  ,tracor + '_new_pet', N[1] + '_' + S[4] + '_' + S[1] + cond)
    elif type_ses == 'way':
      cond ='_WAY'
      if Nb_Scond<6:
        pet_raw_dcm = opj(Major_path,study,'sourcedata' ,'dicoms',N[1] + '_' + N[2], S[2] + '_' + S[3]  ,tracor + '_new_pet', N[1] + '_' + S[4] + '_' + S[1] + cond)
      else:
        pet_raw_dcm = opj(Major_path,study,'sourcedata' ,'dicoms',N[1] + '_' + N[2], S[2] + '_' + S[3]  ,tracor + '_new_pet', N[1] + '_' + S[5] + '_' + S[1] + cond)
    elif type_ses == 'sb':
      cond ='_SB'
    
    pet_raw_dcm = opj(Major_path,study,'sourcedata' ,'dicoms',N[1] + '_' + N[2], S[2] + '_' + S[3]  ,tracor + '_new_pet', N[1] + '_' + S[4] + '_' + S[1] + cond)
    dir_dcm = sorted(os.listdir(pet_raw_dcm))
    dir_dcm = opj(pet_raw_dcm,dir_dcm[-1])
    
    
  elif study == 'imagina':
    data_path_raw   = opj(Major_path, study, 'rawdata', Subject_Dir, Sess_to_process)
    
    if type_sess[:-1] == 'control':
      cond = Sess[0] + '-baseline_' + type_sess[-1]  
    else: 
      if len(loca_sess)>1:
        cond = Sess[0] + '-' + loca_sess[0] + '_' + loca_sess[1][:-1] + '_' + loca_sess[1][-1] 
    
    pet_raw_dcm = opj(Major_path,study,'sourcedata', Subname, cond ,'DICOM','PET_FDG')
    dir_dcm = opj(pet_raw_dcm,'PET_DATA_dyn')
  
  
  pet_raw_dir = opj(data_path_raw,'pet')

  # anatomy
  path_anat     = opj(opd(pet_dir),'anat')
  dir_transfo   = opj(path_anat,'matrices')
  dir_prepro    = opj(path_anat,'preprocess')
  volumes_dir   = opj(path_anat,'Wb','volumes')
  masks_dir     = opj(volumes_dir,'masks')
  labels_dir    = opj(volumes_dir,'labels')

  Ref_native_brain      = opj(volumes_dir,Subname + '_space-acpc_desc-SS_T1w.nii.gz')
  Ref_native_template   = opj(volumes_dir,Subname + '_space-acpc_desc-template_T1w.nii.gz')
  
  Ref_native_template_mask     = opj(masks_dir,  Subname + '_mask.nii.gz')
  
  Ref_native_raw              = opj(dir_prepro, Subname  + '_desc-preproc-n3Bias_T1w.nii.gz')
  Ref_native_acpc             = opj(dir_prepro, Subname  + '_space-acpc_T1w.nii.gz')
  
  
  Ref_native_template_mask_D  = opj(masks_dir, Subname   + '_desc-dilat_mask.nii.gz')
  Ref_native_template_atlas   = opj(labels_dir, Subname  + '_label-'+ template_name + '.nii.gz')
  
  
  REF_Brain    = opj(path_anat,'templates',template_folder ,'volumes', Subname + '_space-'+ template_name + '_desc-SS_T1w.nii.gz')
    
  if ope(pet_dir) == False:
    os.makedirs(pet_dir)
    
  # Create (or implement) a diary : one per day
  date_file = datetime.date.today()
  ct        = datetime.datetime.now()
  diary_name = Subject_Dir + '_' + str(date_file) + '.txt'
  diary = open(opj(pet_dir,diary_name), "w")
  
  if opi(diary_name) == False:   
    diary.write('Process of the pet data')
  diary.write(f'\n{ct}')
  
  print(dir_dcm) 

  # get useful informations

  Study_n = 'study_01234'
  
  ratio = Bq_inj/weight *1000
  
  list_dcm = sorted(glob.glob(opj(dir_dcm,'*.ima')))
  ds1 = pdcm.dcmread(opj(dir_dcm,list_dcm[0]))
  
  Iso_code = ds1[0x54,0x16][0][0x54,0x300][0][0x8,0x104].value
  ISO      = Iso_code.split('^')
  isotope  = ISO[0] + '-' + ISO[1]
  
  Day_scan =  ds1[0x8,0x20].value
  acq_scan =  ds1[0x8,0x31].value
  
  START = Day_scan[6:8] + '/' + Day_scan[4:6] + '/' + Day_scan[0:4] + \
  ' ' + acq_scan[0:2] + ':' + acq_scan[2:4] + ':' + acq_scan[4:6]
  
  
  # extracting some usefull information ...............................................................................................................

  f = open(opj(pet_raw_dir,Subname + '_trc-'+ tracor + '_pet.json'))
  info_pet = json.load(f)
  if info_pet["PatientPosition"] == 'FFS':
    orientation = 'RSA'  # RIA
  elif info_pet["PatientPosition"] == 'HFS':
    orientation = 'LIP'
  print(orientation)
  FS = info_pet["FrameTimesStart"]
  FT = info_pet["FrameReferenceTime"]
  FD = info_pet["FrameDuration"]
  
  
  # create the sif file :
  # Weighting the Radioactivity concentration according to signal intensity, physical decay, frame duration  and voxel volume

  SIF = open(opj(pet_raw_dir,Subname + '.txt'),"w")
  SIF.write(START + ' ' + str(len(FT)) + ' 4 1 ' + Study_n + ' ' + isotope)
  for j in range(len(FT)):
    SIF.write('\n' + str(FS[j]) + ' ' + str(int(FS[j]) + int(FD[j])) + ' 0 0')
  SIF.close()
  shutil.move(opj(pet_raw_dir,Subname + '.txt'), opj(pet_raw_dir,Subname + '_trc-'+ tracor + '_pet.sif'))
  
  cmd = 'gzip -d ' + opj(pet_raw_dir,Subname+ '_trc-'+ tracor + '_pet.nii.gz')
  nl = spgo(cmd); diary.write(f'\n{nl}'); print(nl)
  
  cmd = 'imgweigh ' + opj(pet_raw_dir,Subname + '_trc-'+ tracor + '_pet.nii') + ' ' + opj(pet_raw_dir,Subname+ '_trc-' + tracor + '_pet.sif')
  nl = spgo(cmd); diary.write(f'\n{nl}'); print(nl)
  
  cmd = 'gzip ' + opj(pet_raw_dir,Subname+ '_trc-'+ tracor + '_pet.nii')
  nl = spgo(cmd); diary.write(f'\n{nl}'); print(nl)
  
  
  '''
  dictionary = {"study": study,
                "Subject": Subject_Dir,
                "Session": Sess_to_process,
                "Template": reference,
                "Tracor": Realign_nb+1,
                "Slice order": tshift,
                "Remove the first n volumes ": T1_eq,
                "Spatial blurring": blur,
                "manual correction of brain mask": will[Need_check],
                "ICA ": ['nb per run: ' + str(nb_ICA_run),
                         'nb when concatenated: ' + str(nb_ICA),
                         'ICA in native space at the end of the preprocessing: ' +  will[ICA_native_end],
                         'manual cleaning in native space: ' +  will[ICA_native_clean],
                         'ICA in template space at the end of the preprocessing: ' +  will[ICA_ref_end],
                         'manual cleaning in native space: ' +  will[ICA_ref_clean]],
                "Field distortion correction method": [Fmap,
                                                       G_fmap],}
  json_object = json.dumps(dictionary, indent=2)
  with open(opj(dir_path_RS,'Settings.json'), "w") as outfile:
    outfile.write(json_object)
  '''  
  
  diary.write(f'\n')
  diary.close()
  
  
  return [isotope,Study_n,ratio,pet_dir,Subname,orientation,diary_name,dir_transfo,REF,
          Ref_native_brain,Ref_native_template,Ref_native_template_mask,Ref_native_raw,Ref_native_acpc,template_folder,template_name,ref_mask,
          Ref_native_template_mask_D,Ref_native_template_atlas,REF_Brain,pet_raw_dir,FT,FD]