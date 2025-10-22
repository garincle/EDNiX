#import
import os
import subprocess
import json
import shutil
import datetime


opj = os.path.join
spgo = subprocess.getoutput

def doit(pet_dir,pet_raw_dir,Subname,tracor,orientation,Ref_native_raw,Need_check,diary_name):
  
  
  ct = datetime.datetime.now()
  diary = open(opj(pet_dir,diary_name), "a")
  diary.write(f'\n{ct}')
  nl = 'REORIENT THE FILE'
  diary.write(f'\n{nl}'); print(nl)
  
  Name = Subname + '_trc-'+ tracor
  
  cmd = '3dTstat -mean -prefix ' + opj(pet_dir,Name   + '_desc-mean_pet.nii.gz') + ' ' + opj(pet_raw_dir,Name  + '_pet.nii.gz'),
  nl = spgo(cmd)
  diary.write(f'\n{nl}')
  print(nl)
  
  
  shutil.copyfile(opj(pet_dir,Name  + '_desc-mean_pet.nii.gz'),opj(pet_dir,Name  + '_desc-reorient-mean_pet.nii.gz'))
  
  cmd = '3drefit -orient ' +  orientation + ' ' + opj(pet_dir,Name  + '_desc-reorient-mean_pet.nii.gz')
  nl = spgo(cmd)
  diary.write(f'\n{nl}')
  print(nl)
  
  shutil.copyfile(opj(pet_raw_dir,Name  + '_pet.nii.gz'),opj(pet_dir,Name  + '_desc-reorient_pet.nii.gz'))
  
  cmd = '3drefit -orient ' +  orientation + ' ' +   opj(pet_dir,Name  + '_desc-reorient_pet.nii.gz')
  nl = spgo(cmd)
  diary.write(f'\n{nl}')
  print(nl)
  
  print('better check that the orientation is correct')
  if Need_check ==1:
    cmd = 'freeview -v ' + opj(pet_dir,Name  + '_desc-reorient-mean_pet.nii.gz') + \
    ' ' + Ref_native_raw + ':colormap=heat:opacity=0.5:visible=1'
    nl = spgo(cmd)
    diary.write(f'\n')
    diary.write('check of the image')
  
  diary.write(f'\n')
  diary.close()