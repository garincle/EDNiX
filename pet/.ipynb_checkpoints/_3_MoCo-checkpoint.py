#import
import os
import subprocess
import json
import shutil
import datetime
import ants


opj = os.path.join
ope = os.path.exists
spgo = subprocess.getoutput

def doit(pet_dir,Subname,tracor,nb_void,diary_name):
  
  
  ct = datetime.datetime.now()
  diary = open(opj(pet_dir,diary_name), "a")
  diary.write(f'\n{ct}')
  nl = 'Motion correction'
  diary.write(f'\n{nl}'); print(nl)
  
  if ope(opj(pet_dir,'tmp')) == False:
    os.makedirs(opj(pet_dir,'tmp'))

  Name = Subname + '_trc-'+ tracor
  
  shutil.move(opj(pet_dir,Name + '_desc-reorient-mean_pet.nii.gz'),opj(pet_dir,'tmp', Name+ '_desc-mean-crop_pet.nii.gz'))

  cmd = 'freeview -v ' + opj(pet_dir,'tmp',Name + '_desc-mean-crop_pet.nii.gz') 
  nl = spgo(cmd); diary.write(f'\n{nl}'); print(nl)

  ref_tmp = ants.image_read(opj(pet_dir,'tmp',Name + '_desc-mean-crop_pet.nii.gz'))

  cmd = '3dresample -master ' + opj(pet_dir,'tmp',Name + '_desc-mean-crop_pet.nii.gz') + \
  ' -prefix ' + opj(pet_dir,Name + '_desc-crop_pet.nii.gz') + ' -input ' + opj(pet_dir,Name + '_desc-reorient_pet.nii.gz')
  nl = spgo(cmd); diary.write(f'\n{nl}'); print(nl)
  
    cmd = '3dTstat -mean -prefix ' + opj(pet_dir,Name + '_desc-reorient-mean_pet.nii.gz') + ' ' + opj(pet_dir,Name + '_desc-reorient_pet.nii.gz')
  nl = spgo(cmd); diary.write(f'\n{nl}'); print(nl)

  img_tmp       = ants.image_read(opj(pet_dir,Name + '_desc-crop_pet.nii.gz'))
  tmp_unmerged  = ants.ndimage_to_list(img_tmp)

  To_be_corrected  = ants.image_read(opj(pet_dir,Name + '_desc-reorient_pet.nii.gz'))
  ref_mean         = ants.image_read(opj(pet_dir,Name + '_desc-reorient-mean_pet.nii.gz'))
  img_unmerged = ants.ndimage_to_list(To_be_corrected)

  motion_corrected = list()
  for i in range(len(img_unmerged)): 
      if i > nb_void:
          mc = ants.registration(fixed=ref_tmp, moving=tmp_unmerged[i], type_of_transform = 'Rigid',
                               outprefix=opj(pet_dir,'tmp',Name + '_mc_'+ str(i) + '_'))
          moved = ants.apply_transforms(fixed=ref_mean, moving=img_unmerged[i],transformlist=opj(pet_dir,'tmp',Name + '_mc_'+ str(i) + '_0GenericAffine.mat')
                                      ,interpolator='linear',whichtoinvert=[True])
        
      else:
        moved=img_unmerged[i]
      motion_corrected.append(moved)
    
  motCorr = ants.list_to_ndimage(To_be_corrected, motion_corrected)
  ants.image_write(motCorr, opj(pet_dir,Name + '_desc-MoCo_pet.nii.gz'))

  diary.write(f'\n')
  diary.close()