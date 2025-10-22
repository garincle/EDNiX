#import
import os
import glob
import json
import datetime
import numpy as np
import fnmatch

opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile


from atlases import templatefeat
import study_struct
import dcmFile
import sifFile

AFNI_sif  = 'afni_ub24_latest.sif'
FSL_sif   = 'fsl_6.0.5.1-cuda9.1.sif'
FS_sif    = 'freesurfer_NHP.sif'
WB_sif    = 'connectome_workbench_1.5.0.sif'
ITK_sif   = 'itksnap_5.0.9.sif'
turku_sif = 'tpcclib.sif'


def orga(study,Major_path,path_ATLAS,singularity,Subject_Dir,Sess_to_process,Sess_anat,anat_ref,specie,
         template,atlas_name,
         tracor,im_type,to_remove,Do_cut,duration,
         Bq_inj,
         weight,
         blood_sugar,
         REF_name,
         nbInj):


  if singularity == 'no':
    sing_afni = ''
    sing_fsl = ''
    sing_fs = ''
    sing_wb = ''
    sing_itk = ''
    sing_turku = ''

  elif singularity == 'yes':
    working_path = opj('/', 'scratch', 'in_process')
    s_path = opj(Major_path, study, 'BIDS', 'code', 'singularity')
    sing_func1 = 'singularity run --bind '  + ','.join([working_path, opj(Major_path, study), path_ATLAS])
    sing_func2 = 'singularity exec --bind ' + ','.join([working_path, opj(Major_path, study), path_ATLAS])
    sing_afni  = sing_func1 + ' ' + opj(s_path, AFNI_sif) + ' '
    sing_fsl   = sing_func1 + ' ' + opj(s_path, FSL_sif) + ' '
    sing_fs    = sing_func1 + ' ' + opj(s_path, FS_sif) + ' '
    sing_itk   = sing_func1 + ' ' + opj(s_path, ITK_sif) + ' '
    sing_wb    = sing_func1 + ' ' + opj(s_path, WB_sif) + ' '
    sing_turku = sing_func2 + ' ' + opj(s_path, turku_sif) + ' '

  
  ## Organization 
  
  # references 

  for dirpath, dirnames, filenames in os.walk(path_ATLAS):
    folder_name = fnmatch.filter(dirnames, specie)
    if folder_name:
      if not opb(dirpath) == 'freesurfer':
        BASE_ref = opj(dirpath, folder_name[0])


  # Datas
  N = Subject_Dir.split('-')
  Subname = '-'.join([N[0], N[1]])
  
  pet_dir    = opj(Major_path,study,'BIDS','derivatives', Subject_Dir,Sess_to_process,'pet')
  pet_prepro = opj(pet_dir,'preprocessing')
  pet_raw_dir,dir_dcm = study_struct.path(Major_path,study,Subject_Dir,Sess_to_process,tracor,nbInj)

  # anatomy
  if Sess_anat == '':
    path_anat     = opj(opd(pet_dir),'anat')
  else:
    path_anat     = opj(Sess_anat,'anat')
  
  dir_transfo   = opj(path_anat,'matrices')
  dir_prepro    = opj(path_anat,'preprocessing')
  volumes_dir   = opj(path_anat,'native','volumes')
  masks_dir     = opj(volumes_dir,'masks')

  T1wacpc        = opj(volumes_dir,Subname + '_space-acpc_desc-template_' + anat_ref + '.nii.gz')
  T1wacpcmask    = opj(masks_dir,  Subname + '_mask.nii.gz')
  T1worig        = opj(dir_prepro, Subname + '_desc-preproc-n3Bias_' + anat_ref + '.nii.gz')
  T1wacpcfull    = opj(dir_prepro, Subname  + '_space-acpc_' + anat_ref + '.nii.gz')
  T1wacpcmaskD   = opj(masks_dir, Subname   + '_desc-dilat_mask.nii.gz')
  templateT1w    = opj(path_anat,'templates',template ,'volumes', Subname + '_space-'+ template + '_desc-SS_' + anat_ref + '.nii.gz')
  

  [_, templatemask, REF,atlas, atlas_label,surf_dir,LRside,roi] = templatefeat.get(specie, path_ATLAS, '',template,
  atlas_name, REF_name,opj(path_anat, 'templates',template),BASE_ref, 'pet')
  

  if ope(pet_dir) == False:
    os.makedirs(pet_dir)
    if ope(pet_prepro) == False:
      os.makedirs(pet_prepro)
    
  # Create (or implement) a diary : one per day
  date_file = datetime.date.today()
  ct        = datetime.datetime.now()
  diary_name = opj(pet_dir,Subject_Dir + '_' + str(date_file) + '.txt')

  if opi(diary_name) == False:
    diary = open(diary_name, "w")
    diary.write('Process of the pet data')
  else :
    diary = open(diary_name, "a")
  diary.write(f'\n{ct}')
  
  print(dir_dcm) 

  # get useful information

  if not Bq_inj == '' and not weight == '': 
    if not blood_sugar == '':
      ratio = (Bq_inj/weight *1000)/blood_sugar
    else:
      blood_sugar = 'na'
      ratio = Bq_inj / weight * 1000
    maxval = '2.5'
  else:
    Bq_inj = 'na'
    weight = 'na'
    ratio  = 1
    blood_sugar = 'na'
    maxval = '10'
  

  # Extract some useful information ...............................................................................................................
  # From the dicom
  list_dcm = sorted(glob.glob(opj(dir_dcm,'*.ima')))
  isotope,START = dcmFile.extract(opj(dir_dcm,list_dcm[0]))


  list_img = sorted(glob.glob(opj(pet_raw_dir,Subname + '_trc-' + tracor + '*_pet.nii.gz')))
  print(list_img)

  frame_stop = []
  steps      = []
  orientation=[]
  
  for i in range(len(list_img)):
    file = opj(pet_raw_dir,list_img[i])

    # From the json
    f = open(list_img[i].replace('.nii.gz','.json'))
    info_pet = json.load(f)
    if info_pet["PatientPosition"] == 'FFS':
      orientation.append('RIA')  # RIA or RSA
    elif info_pet["PatientPosition"] == 'HFS':
      orientation.append('LIP')
    print(orientation[i])
    FS=[]
    FT=[]
    FD=[]

    if im_type[i] =='dynamic':
      FS.append(info_pet["FrameTimesStart"])
      FT.append(info_pet["FrameReferenceTime"])
      FD.append(info_pet["FrameDuration"])
      
      # create the sif file if needed:
      sifFile.create(file,START,FT[i],FS[i],FD[i],isotope,diary_name,sing_turku)
      
      if Do_cut[i] == 0:
        frame_stop.append('na')
        steps.append('na')

      elif Do_cut[i] == 1:
        f = open(file.replace('.nii.gz','.sif'), "r")
        matrix = f.readlines()
        matrix = matrix[1:]
        
        starts = np.zeros(len(matrix))
        stops  = np.zeros(len(matrix))
        Bq     = np.zeros(len(matrix))
        for q in range(len(matrix)):
          starts[q] = int(matrix[q].split(' ')[0])
          stops[q]  = int(matrix[q].split(' ')[1])
          Bq[q]     = int(matrix[q].split(' ')[2])
      
        steps.append(np.arange(starts[0], stops[-1], duration * 60)) # blocks
        frame_stop.append(stops)

  dictionary = {"study": study,
                "Subject":  Subname,
                "Session": Sess_to_process,
                "Template": template,
                "atlas": atlas_name,
                "Orientation": orientation,
                "Tracor": tracor,
                "isotope": isotope,
                "Nb of recordings" : str(len(list_img)),
                "date first recording":  START,
                "type of image": im_type,
                "nb of ignored frames for motion correction" : str(to_remove),
                "Duration of images" : str(duration),
                "Dose injected (in MBq)" : str(Bq_inj),
                "Subject Weight (in Kg)" : str(weight),
                "Blood_sugar (in g/L)" : str(blood_sugar),
                "Ratio " : str(ratio),
                }
  json_object = json.dumps(dictionary, indent=2)
  with open(opj(pet_dir,'Prepro_Settings.json'), "w") as outfile:
    outfile.write(json_object)


  return [Subname,pet_raw_dir,pet_dir,pet_prepro,list_img,
  dir_transfo,T1wacpc,T1wacpcmask,T1wacpcmaskD,T1wacpcfull,T1worig,templateT1w,templatemask,atlas,atlas_label,surf_dir,
  LRside,roi,isotope,orientation,FT,FD,REF,ratio,steps,frame_stop,maxval,
  diary_name,sing_afni,sing_fsl,sing_fs,sing_wb,sing_itk,sing_turku]