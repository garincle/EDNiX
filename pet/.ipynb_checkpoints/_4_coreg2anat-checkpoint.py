#import
import os
import subprocess
import glob
import json
import shutil
import ants
import math
import numpy as np
import pandas as pd
import datetime
import matplotlib.pyplot as plt
import pydicom as pdcm


opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile
spco = subprocess.check_output
spgo = subprocess.getoutput


def doit():
  
    ct = datetime.datetime.now()
  diary = open(opj(pet_dir,diary_name), "a")
  diary.write(f'\n{ct}')
  nl = 'Motion correction'
  diary.write(f'\n{nl}'); print(nl)
  
  nl = 'Co-registration with the T1w MRI anatomy'
  diary.write(f'\n{nl}'); print(nl)
  cmd = '3dTstat -mean -prefix ' + opj(pet_dir,Subname + '_desc-MoCo-mean_pet.nii.gz') + ' ' + opj(pet_dir,Subname + '_desc-MoCo_pet.nii.gz')
  nl = spgo(cmd); diary.write(f'\n{nl}'); print(nl)

  cmd='3dresample -master ' + Ref_native_raw + ' -rmode Cu -prefix ' + opj(pet_dir,'tmp','mean_pet2T1w.nii.gz') + \
  ' -input ' +  opj(pet_dir,Subname + '_desc-MoCo-mean_pet.nii.gz')
  nl = spgo(cmd); diary.write(f'\n{nl}'); print(nl)

shutil.copyfile(opj(pet_dir,'tmp','mean_pet2T1w.nii.gz'),opj(pet_dir,'tmp','mean_pet2T1w_tmp.nii.gz'))
cmd = 'freeview -v ' + opj(pet_dir,'tmp','mean_pet2T1w_tmp.nii.gz') + \
' ' + Ref_native_raw + ':colormap=heat:opacity=0.5:visible=1'
nl = spgo(cmd); diary.write(f'\n{nl}'); print(nl)


mPet_old  = ants.image_read(opj(pet_dir,'tmp','mean_pet2T1w.nii.gz'))
mPet_new  = ants.image_read(opj(pet_dir,'tmp','mean_pet2T1w_tmp.nii.gz'))
mTx_Pet2T1w_tmp = ants.registration(fixed=mPet_new, moving=mPet_old, type_of_transform = 'Rigid', outprefix=opj(pet_dir,'tmp','mean_pet2T1w_new'))

brain = ants.image_read(Ref_native_raw)
mTx_Pet2T1w = ants.registration(fixed=brain, moving=mPet_old, type_of_transform = 'Rigid', 
                                initial_transform=opj(pet_dir,'tmp','mean_pet2T1w_new0GenericAffine.mat'),
                                outprefix=opj(pet_dir,Subname + '_from-Pet_to-T1w_'))
# Transform into the T1w MRI anatomy ..........................................................................................................
moved = ants.apply_transforms(fixed=brain, moving=mPet_old ,transformlist=mTx_Pet2T1w['fwdtransforms'],interpolator='hammingWindowedSinc')
ants.image_write(moved, opj(pet_dir,'tmp',Subname + '_mean_Pet2T1w.nii.gz'), ri=False)


# Transform the ROIs defined by the atlas to fit the PET data .................................................................................
Mean       = ants.image_read(opj(pet_dir,Subname + '_desc-MoCo-mean_pet.nii.gz'))
brain_acpc = ants.image_read(Ref_native_template)
atlas      = ants.image_read(Ref_native_template_atlas)

moved = ants.apply_transforms(fixed=Mean, moving=atlas,
                              transformlist=[opj(pet_dir,Subname + '_from-Pet_to-T1w_0GenericAffine.mat'),
                                             opj(dir_transfo,'shift_0GenericAffine.mat'),
                                             opj(dir_transfo,'acpc_0GenericAffine.mat')],
                              interpolator='nearestNeighbor',whichtoinvert=[True,True,True])
ants.image_write(moved, opj(pet_dir,'tmp',Subname + '_atlas_tmp.nii.gz'), ri=False)

cmd = '3dresample -master ' + opj(pet_dir,Subname + '_desc-reorient-mean_pet.nii.gz') + ' -prefix ' + opj(pet_dir,Subname + '_desc-coreg_dseg.nii.gz') + \
' -input ' +  opj(pet_dir,'tmp',Subname + '_atlas_tmp.nii.gz')
nl = spgo(cmd); diary.write(f'\n{nl}'); print(nl)

moved = ants.apply_transforms(fixed=Mean, moving=brain_acpc,
                              transformlist=[opj(pet_dir,Subname + '_from-Pet_to-T1w_0GenericAffine.mat'),
                                             opj(dir_transfo,'shift_0GenericAffine.mat'),
                                             opj(dir_transfo,'acpc_0GenericAffine.mat')],
                              interpolator='nearestNeighbor',whichtoinvert=[True,True,True])
ants.image_write(moved, opj(pet_dir,'tmp','T1w_pet_tmp.nii.gz'), ri=False)

cmd = '3dresample -master ' + opj(pet_dir,Subname + '_desc-reorient-mean_pet.nii.gz') + ' -prefix ' + opj(pet_dir, Subname + '_desc-coreg_T1w.nii.gz') + \
' -input ' +  opj(pet_dir,'tmp','T1w_pet_tmp.nii.gz')
nl = spgo(cmd); diary.write(f'\n{nl}'); print(nl)