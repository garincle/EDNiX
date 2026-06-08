#import
import os
import sys
import json
import shutil

opj = os.path.join
opd = os.path.dirname
opb = os.path.basename
ope = os.path.exists
opi = os.path.isfile

sys.path.insert(1, '../tools')
import run_cmd

def doit(new_dir,img,orientation,Ref_native_raw,im_type,Need_check,diary_name,sing_afni,sing_fs):

  N = opb(img).split('_')
  E = -1
  for i in range(len(N)):
    a = N[i].split('-')
    if a[0] == 'desc':
      E = i
  if E == -1:
    E = N.index('pet.nii.gz')
  Name = '_'.join(N[0:E])

  nl = 'REORIENT THE FILE'
  run_cmd.msg(nl, diary_name)

  if im_type == 'dynamic':

    cmd = (sing_afni + '3dTstat -mean -prefix ' + opj(new_dir, Name + '_desc-mean_pet.nii.gz') + ' ' + img)
    run_cmd.run(cmd,diary_name)

    dictionary = {"Sources":  img,
                  "Description": 'Mean image (3dTstat, AFNI).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(new_dir,Name + '_desc-mean_pet.json'), "w") as outfile:
      outfile.write(json_object)

    shutil.copyfile(opj(new_dir,Name + '_desc-mean_pet.nii.gz'),
                    opj(new_dir, Name + '_desc-reorient-mean_pet.nii.gz'))
  
    cmd = (sing_afni + '3drefit -orient ' +  orientation + ' ' + opj(new_dir, Name + '_desc-reorient-mean_pet.nii.gz'))
    run_cmd.run(cmd,diary_name)

    dictionary = {"Sources": opj(new_dir,Name + '_desc-mean_pet.nii.gz'),
                  "Description": 'header reorientation according to : ' + orientation + ' (3drefit, AFNI).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(new_dir, Name + '_desc-reorient-mean_pet.json'), "w") as outfile:
      outfile.write(json_object)
    to_see = opj(new_dir, Name + '_desc-reorient-mean_pet.nii.gz')
  else:

    to_see = opj(new_dir, Name + '_desc-reorient_pet.nii.gz')

  shutil.copyfile(img, opj(new_dir, Name + '_desc-reorient_pet.nii.gz'))
  if opi(img.replace('.nii.gz','.sif')):
    shutil.copyfile(img.replace('.nii.gz','.sif'), opj(new_dir, Name + '_desc-reorient_pet.sif'))
  
  cmd = (sing_afni + '3drefit -orient ' +  orientation + ' ' +   opj(new_dir, Name + '_desc-reorient_pet.nii.gz'))
  run_cmd.run(cmd,diary_name)


  dictionary = {"Sources": img,
                "Description": 'Header reorientation according to : ' + orientation + ' (3drefit, AFNI).', }
  json_object = json.dumps(dictionary, indent=2)
  with open(opj(new_dir, Name + '_desc-reorient_pet.json'), "w") as outfile:
    outfile.write(json_object)
  
  print('better check that the orientation is correct')
  if Need_check ==1:
    nl = 'check of the orientation'
    run_cmd.msg(nl, diary_name)
    cmd = (sing_fs + 'freeview -v ' + to_see + ' ' + Ref_native_raw + ':colormap=heat:opacity=0.5:visible=1')
    run_cmd.do(cmd, diary_name)

