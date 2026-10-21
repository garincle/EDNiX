#import
import os
import glob
import json
import shutil
import ants

opj = os.path.join
opd = os.path.dirname
opb = os.path.basename
ope = os.path.exists
opi = os.path.isfile

antsaff     = '0GenericAffine.mat'
antswarp    = '1Warp.nii.gz'
antsinvwarp = '1InverseWarp.nii.gz'

resamp_interp = '-rmode Cu'
petinterp     = 'linear'

from Tools import run_cmd,extract_filename
from modalities.pet import studyT

def dynamic(img,cropped_3D,nb_void,Need_check,diary_name,sing_afni,sing_fs):
    
    nl = 'Motion correction'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    pet_dir = opd(img)
    mcdir   = opd(cropped_3D)

    fname,descname = extract_filename.get_desc(opb(img),'desc')

    cropped_4D   = opj(mcdir, '_'.join([fname, 'desc-crop', 'pet.nii.gz']))
    
    mc4D   = opj(pet_dir, '_'.join([fname, 'desc-MoCo-within','pet.nii.gz']))
    mc3D   = opj(pet_dir, '_'.join([fname, 'desc-MoCo-within-mean', 'pet.nii.gz']))

    
    cmd = (sing_afni + '3dresample ' + resamp_interp + ' -overwrite -master ' + cropped_3D +' -prefix ' + cropped_4D + ' -input ' + img)
    run_cmd.run(cmd,diary_name)

    dictionary = {"Sources": [img,
                              cropped_3D],
                  "Description": 'Cropped image (3dresample, AFNI).', }

    json_object = json.dumps(dictionary, indent=2)
    with open(cropped_4D.replace('.nii.gz','.json'), "w") as outfile:
        outfile.write(json_object)

    ref_img       = ants.image_read(cropped_3D)
    img_tmp       = ants.image_read(cropped_4D)
    tmp_unmerged  = ants.ndimage_to_list(img_tmp)

    To_be_corrected  = ants.image_read(img)
    img_unmerged     = ants.ndimage_to_list(To_be_corrected)

    motion_corrected = list()
    for i in range(len(img_unmerged)):
        if i > nb_void:
            ants.registration(fixed=tmp_unmerged[i], moving=ref_img,
                                   type_of_transform='Rigid',
                                   outprefix=opj(mcdir,'_'.join([fname,'mc'+ str(i),''])))

            moved = ants.apply_transforms(fixed=mcdir, moving=img_unmerged[i],
                                          transformlist=opj(mcdir,'_'.join([fname,'mc'+ str(i),antsaff])),
                                          whichtoinvert=[True])
        else:
            moved=img_unmerged[i]
        motion_corrected.append(moved)
    
    motCorr = ants.list_to_ndimage(To_be_corrected, motion_corrected)
    ants.image_write(motCorr, mc4D, ri=False)
    shutil.copyfile(img.replace('.nii.gz','.sif'), mc4D.replace('.nii.gz','.sif'))

    dictionary = {"Sources": [img,
                              cropped_4D,
                              cropped_3D],
                  "Description": 'Motion correction (Antspy).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(mc4D.replace('.nii.gz','.json'), "w") as outfile:
        outfile.write(json_object)

    cmd = (sing_afni  + '3dTstat -overwrite -mean -prefix ' + mc3D + ' ' + mc4D)
    run_cmd.run(cmd,diary_name)

    dictionary = {"Sources": mc4D,
                  "Description": 'Mean image (3dTstat, AFNI).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(mc3D.replace('.nii.gz','.json'), "w") as outfile:
        outfile.write(json_object)

    if Need_check==1:
        cmd = (sing_fs + 'freeview -v ' + mc4D + ' -subtitle "Check-residual-head-motion ! "')
        run_cmd.do(cmd,diary_name)


def static(list_img,subname,tracor,dir_prepro_orig,dir_prepro_orig_process,diary_name,sing_afni,sing_fs):


    tmp_dir = opj(dir_prepro_orig_process,'forMC')

    mean_img      = opj(dir_prepro_orig_process, '_'.join([subname, 'trc-' + tracor  ,'desc-mean', 'pet.nii.gz']))
    crop3D        = opj(tmp_dir, '_'.join([subname, 'trc-' + tracor , 'desc-mean-crop', 'pet.nii.gz']))
    MoCo_mean_img = opj(dir_prepro_orig_process, '_'.join([subname, 'trc-' + tracor , 'desc-MoCo-between-mean', 'pet.nii.gz']))


    nl = 'Between-images motion correction'
    run_cmd.msg(nl, diary_name,'OKGREEN')
    
    os.makedirs(tmp_dir,exist_ok=True)

    # Do the mean image
    mean_list = ' '.join(list_img)
    print(mean_list)

    cmd = (sing_afni + '3dMean -overwrite -float -prefix ' + mean_img + ' ' + mean_list)
    run_cmd.run(cmd, diary_name)

    ref_mean = ants.image_read(mean_img)

    # Cropped the mean_image
    studyT.manual_crop(tmp_dir,mean_img,'desc-mean-crop',diary_name,sing_fs)

    ref_tmp  = ants.image_read(crop3D)
  
    for s in range(len(list_img)):

        Name_file = extract_filename.get_runname(opb(list_img[s]))

        cropfile = opj(tmp_dir, '_'.join([Name_file,'desc-crop', 'pet.nii.gz']))
        MCfile   = opj(dir_prepro_orig, '_'.join([Name_file,  'desc-MoCo-between', 'pet.nii.gz']))
        allfile  = '_'.join(['*',  'desc-MoCo-between', 'pet.nii.gz'])

        cmd = (sing_afni + '3dresample -overwrite -master ' + crop3D +
               ' -prefix ' + cropfile +
               ' -input '  + list_img[s])
        run_cmd.run(cmd,diary_name)

        dictionary = {"Sources": [list_img[s],
                                  crop3D],
                      "Description": 'Cropped image (3dresample, AFNI).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(cropfile.replace('.nii.gz','.json'), "w") as outfile:
            outfile.write(json_object)

        img_tmp          = ants.image_read(cropfile)
        To_be_corrected  = ants.image_read(list_img[s])
        ants.registration(fixed=ref_tmp, moving=img_tmp, type_of_transform='Translation',
                          outprefix=opj(tmp_dir, Name_file + '_'))
        ants.registration(fixed=ref_tmp, moving=img_tmp, type_of_transform = 'Rigid',
                          initial_transform=opj(tmp_dir, '_'.join([Name_file,antsaff])),
                          outprefix=opj(tmp_dir, Name_file + '_'))

        moved = ants.apply_transforms(fixed=ref_mean, moving=To_be_corrected,
                                      transformlist=opj(tmp_dir, '_'.join([Name_file,antsaff])),
                                      interpolator=petinterp,
                                      whichtoinvert=[False])

        ants.image_write(moved, MCfile, ri=False)
        if opi(list_img[s].replace('.nii.gz','.sif')):
            shutil.copyfile(list_img[s].replace('.nii.gz','.sif'), MCfile.replace('.nii.gz','.sif'))

        dictionary = {"Sources": [list_img[s],
                                  cropfile,
                                  crop3D],
                      "Description": 'Motion correction (Antspy).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(MCfile.replace('.nii.gz','.json'), "w") as outfile:
            outfile.write(json_object)
    
    list_name = sorted(glob.glob(opj(dir_prepro_orig,allfile)))
    todolist = ' '.join(list_name)
  
    cmd = (sing_afni + '3dMean -overwrite -float -prefix ' + MoCo_mean_img + ' ' + todolist)
    run_cmd.run(cmd,diary_name)

    dictionary = {"Sources": list_name,
                  "Description": 'Mean image (3dMean, AFNI).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(MoCo_mean_img.replace('.nii.gz','.json'), "w") as outfile:
        outfile.write(json_object)



