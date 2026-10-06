#import
import os
import sys
import glob
import json
import shutil
import ants

opj = os.path.join
opd = os.path.dirname
opb = os.path.basename
ope = os.path.exists
opi = os.path.isfile

sys.path.insert(1, '../tools')
import run_cmd

def dynamic(img,ref_im,nb_void,Need_check,diary_name,sing_afni,sing_fs):

    pet_prepro = opd(img)
    pet_dir = opd(pet_prepro)
    N = opb(img).split('_')
    E = -1
    for i in range(len(N)):
        a = N[i].split('-')
        if a[0] == 'desc':
            E = i
    if E == -1:
        E = N.index('pet.nii.gz')
    Name = '_'.join(N[0:E])

    nl = 'Motion correction'
    run_cmd.msg(nl, diary_name)
  
    if ope(opj(pet_prepro,'tmp')) == False:
        os.makedirs(opj(pet_prepro,'tmp'))

    if ref_im == '':
        a = N[E].split('-')

        cmd = (sing_afni + '3dTstat -overwrite -mean -prefix ' + opj(pet_prepro, Name + '_desc-' + a[1] + '-mean_pet.nii.gz') +
               ' ' + img + '[' + str(nb_void) + '..$]')
        run_cmd.run(cmd,diary_name)

        dictionary = {"Sources":  img + '[' + str(nb_void) + '..$]',
                      "Description": 'Mean image (3dTstat, AFNI).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(pet_prepro, Name + '_desc-' + a[1] + '-mean_pet.json'), "w") as outfile:
            outfile.write(json_object)

        shutil.copy(opj(pet_prepro, Name + '_desc-' + a[1] + '-mean_pet.nii.gz'),
                    opj(pet_prepro, 'tmp', Name + '_desc-mean-crop_pet.nii.gz'))

        ref_cropp = opj(pet_prepro, 'tmp', Name + '_desc-mean-crop_pet.nii.gz')
        ref_mean = ants.image_read(ref_cropp)
    else:
        ref_mean = ants.image_read(ref_im)
        shutil.copy(ref_im,
                    opj(pet_prepro, 'tmp', Name + '_desc-mean-crop_pet.nii.gz'))
        ref_cropp = ref_im


    cmd = (sing_fs + 'freeview -v ' + opj(pet_prepro,'tmp',Name + '_desc-mean-crop_pet.nii.gz') + ' -subtitle "CROP-THE-IMAGE ! "')
    run_cmd.do(cmd,diary_name)

    dictionary = {"Sources": ref_cropp,
                "Description": 'Manually Cropped image (Freeview, Freesurfer).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(pet_prepro,'tmp',Name + '_desc-mean-crop_pet.json'), "w") as outfile:
        outfile.write(json_object)


    ref_tmp = ants.image_read(opj(pet_prepro,'tmp',Name + '_desc-mean-crop_pet.nii.gz'))

    cmd = (sing_afni + '3dresample -overwrite -master ' + opj(pet_prepro,'tmp',Name + '_desc-mean-crop_pet.nii.gz') +
           ' -prefix ' + opj(pet_prepro,'tmp',Name + '_desc-crop_pet.nii.gz') +
           ' -input ' + img)
    run_cmd.run(cmd,diary_name)

    dictionary = {"Sources": [img,
                              opj(pet_prepro, 'tmp', Name + '_desc-mean-crop_pet.nii.gz')],
                  "Description": 'Cropped image (3dresample, AFNI).', }

    json_object = json.dumps(dictionary, indent=2)
    with open(opj(pet_prepro,'tmp',Name + '_desc-crop_pet.json'), "w") as outfile:
        outfile.write(json_object)

    img_tmp       = ants.image_read(opj(pet_prepro,'tmp',Name + '_desc-crop_pet.nii.gz'))
    tmp_unmerged  = ants.ndimage_to_list(img_tmp)

    To_be_corrected  = ants.image_read(img)
    img_unmerged     = ants.ndimage_to_list(To_be_corrected)

    motion_corrected = list()
    for i in range(len(img_unmerged)):
        if i > nb_void:
            mc = ants.registration(fixed=tmp_unmerged[i], moving=ref_tmp,
                                   type_of_transform='Rigid',
                                   outprefix=opj(pet_prepro,'tmp',Name + '_mc_'+ str(i) + '_'))

            moved = ants.apply_transforms(fixed=ref_mean, moving=img_unmerged[i],
                                          transformlist=opj(pet_prepro,'tmp',Name + '_mc_'+ str(i) + '_0GenericAffine.mat'),
                                          whichtoinvert=[True])
        
        else:
            moved=img_unmerged[i]
        motion_corrected.append(moved)
    
    motCorr = ants.list_to_ndimage(To_be_corrected, motion_corrected)
    ants.image_write(motCorr, opj(pet_dir, Name + '_desc-MoCo_pet.nii.gz'), ri=False)
    shutil.copyfile(img.replace('.nii.gz','.sif'), opj(pet_dir, Name + '_desc-MoCo_pet.sif'))

    dictionary = {"Sources": [img,
                              opj(pet_prepro,'tmp',Name + '_desc-crop_pet.nii.gz'),
                              opj(pet_prepro, 'tmp', Name + '_desc-mean-crop_pet.nii.gz')],
                  "Description": 'Motion correction (Antspy).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(pet_dir, Name + '_desc-MoCo_pet.json'), "w") as outfile:
        outfile.write(json_object)

    cmd = (sing_afni  + '3dTstat -overwrite -mean -prefix ' + opj(pet_dir, Name + '_desc-MoCo-mean_pet.nii.gz') +
           ' ' + opj(pet_dir, Name + '_desc-MoCo_pet.nii.gz'))
    run_cmd.run(cmd,diary_name)

    dictionary = {"Sources": opj(pet_prepro, Name + '_desc-MoCo_pet.nii.gz'),
                  "Description": 'Mean image (3dTstat, AFNI).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(pet_dir, Name + '_desc-MoCo-mean_pet.json'), "w") as outfile:
        outfile.write(json_object)

    if Need_check==1:
        cmd = (sing_fs + 'freeview -v ' + opj(pet_dir, Name + '_desc-MoCo_pet.nii.gz') + ' -subtitle "Check-residual-head-motion ! "')
        run_cmd.do(cmd,diary_name)


def static(list_img,diary_name,sing_afni,sing_fs):

    pet_prepro = opd(list_img[0])
    pet_dir = opd(pet_prepro)
    N = opb(list_img[0]).split('_')
    E = -1
    R = -1
    for i in range(len(N)):
        a = N[i].split('-')
        if a[0] == 'desc':
            E = i
        elif a[0] =='run':
            R = i

    if E == -1:
        E = N.index('pet.nii.gz')
    
    Name = '_'.join(N[0:R])
    
    suffix = N[E].split('-')[-1]

    nl = 'Between-images motion correction'
    run_cmd.msg(nl, diary_name)
    
    if ope(opj(pet_prepro,'tmp')) == False:
        os.makedirs(opj(pet_prepro,'tmp'))

    mean_list = ' '.join(list_img)
    print(mean_list)
    print(Name)
    print(suffix)
    cmd = (sing_afni + '3dMean -overwrite -prefix ' + opj(pet_prepro, Name + '_desc-' + suffix + 'mean_pet.nii.gz') + ' ' + mean_list)
    run_cmd.run(cmd, diary_name)

    ref_mean = ants.image_read(opj(pet_prepro, Name + '_desc-' + suffix + 'mean_pet.nii.gz'))

    shutil.move(opj(pet_prepro, Name + '_desc-' + suffix + 'mean_pet.nii.gz'),
                opj(pet_prepro,'tmp', Name + '_desc-mean-crop_pet.nii.gz'))
  
    # Crop the head. Make sure to include the ears and the nose  but remove the lower jaw
    nl = 'Manual cropping around the head\n'
    run_cmd.msg(nl, diary_name)
    cmd = (sing_fs + 'freeview -v ' + opj(pet_prepro,'tmp', Name + '_desc-mean-crop_pet.nii.gz') + ' -subtitle "Crop The Head ! "')
    run_cmd.do(cmd,diary_name)

    dictionary = {"Sources": opj(pet_prepro, Name + '_desc-' + suffix + 'mean_pet.nii.gz'),
                  "Description": 'Manually Cropped image (Freeview, Freesurfer).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(pet_prepro, Name + '_desc-' + suffix + 'mean_pet.json'), "w") as outfile:
        outfile.write(json_object)

    ref_tmp  = ants.image_read(opj(pet_prepro,'tmp', Name + '_desc-mean-crop_pet.nii.gz'))
  
    for s in range(len(list_img)):
        N = opb(list_img[s]).split('_')
        Name_file = '_'.join(N[0:E])
    
        cmd = (sing_afni + '3dresample -overwrite -master ' + opj(pet_prepro,'tmp', Name + '_desc-mean-crop_pet.nii.gz') +
               ' -prefix ' + opj(pet_prepro,'tmp', Name_file + '_desc-crop_pet.nii.gz') +
               ' -input '  + list_img[s])
        run_cmd.run(cmd,diary_name)

        dictionary = {"Sources": [list_img[s],
                                  opj(pet_prepro,'tmp', Name + '_desc-mean-crop_pet.nii.gz')],
                      "Description": 'Cropped image (3dresample, AFNI).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(pet_prepro,'tmp', Name_file + '_desc-crop_pet.json'), "w") as outfile:
            outfile.write(json_object)

        img_tmp          = ants.image_read(opj(pet_prepro,'tmp', Name_file + '_desc-crop_pet.nii.gz'))
        To_be_corrected  = ants.image_read(list_img[s])

    
        mc = ants.registration(fixed=ref_tmp, moving=img_tmp, type_of_transform = 'Rigid',
                               outprefix=opj(pet_prepro,'tmp', Name_file + '_'))

        moved = ants.apply_transforms(fixed=ref_mean, moving=To_be_corrected,
                                      transformlist=opj(pet_prepro,'tmp', Name_file + '_0GenericAffine.mat'),
                                      interpolator='linear',
                                      whichtoinvert=[False])

        ants.image_write(moved, opj(pet_dir,Name_file + '_desc-MoCo_pet.nii.gz'), ri=False)
        if opi(list_img[s].replace('.nii.gz','.sif')):
            shutil.copyfile(list_img[s].replace('.nii.gz','.sif'), opj(pet_dir, Name_file + '_desc-MoCo_pet.sif'))

        dictionary = {"Sources": [list_img[s],
                                  opj(pet_prepro,'tmp', Name_file + '_desc-crop_pet.nii.gz'),
                                  opj(pet_prepro, Name + '_desc-' + suffix + 'mean_pet.nii.gz')],
                      "Description": 'Motion correction (Antspy).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(pet_dir,Name_file + '_desc-MoCo_pet_pet.json'), "w") as outfile:
            outfile.write(json_object)
    
    list_name = sorted(glob.glob(opj(pet_dir,'*_desc-MoCo_pet.nii.gz')))
    todolist = ' '.join(list_name)
  
    cmd = (sing_afni + '3dMean -overwrite -prefix ' + opj(pet_dir, Name + '_desc-MoCo-mean_pet.nii.gz') + ' ' + todolist)
    run_cmd.run(cmd,diary_name)

    dictionary = {"Sources": list_name,
                  "Description": 'Mean image (3dMean, AFNI).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(pet_dir, Name + '_desc-MoCo-mean_pet.json'), "w") as outfile:
        outfile.write(json_object)


def redo(img1,img2,dyn_img,start,diary_name,sing_afni,sing_fs):

    pet_prepro = opd(img1)
    pet_dir = opd(pet_prepro)
    N = opb(img1).split('_')
    E = -1
    for i in range(len(N)):
        a = N[i].split('-')
        if a[0] == 'desc':
            E = i
    if E == -1:
        E = N.index('pet.nii.gz')
    Name = '_'.join(N[0:E])

    nl = 'Motion correction... again'
    run_cmd.msg(nl, diary_name)
    
    if ope(opj(pet_prepro,'tmp')) == False:
        os.makedirs(opj(pet_prepro,'tmp'))

    shutil.copy(img1, opj(pet_prepro,'tmp', 'crop_ref.nii.gz'))

    # Crop the head. Make sure to include the ears and the nose  but remove the lower jaw
    nl = 'Manual cropping around the head\n'
    run_cmd.msg(nl, diary_name)
    cmd = (sing_fs + 'freeview -v ' + opj(pet_prepro,'tmp', 'crop_ref.nii.gz') + ' -subtitle "Crop The Head ! "')
    run_cmd.do(cmd,diary_name)

    ref_1  = ants.image_read(opj(pet_prepro,'tmp', 'crop_ref.nii.gz'))

    cmd = (sing_afni + '3dresample -overwrite -master ' + opj(pet_prepro,'tmp', 'crop_ref.nii.gz') +
           ' -prefix ' + opj(pet_prepro,'tmp', 'crop_2.nii.gz') +
           ' -input '  + img2)
    run_cmd.run(cmd,diary_name)

    ref_2 = ants.image_read(opj(pet_prepro,'tmp', 'crop_2.nii.gz'))

    mc1 = ants.registration(fixed=ref_1, moving=ref_2, type_of_transform='Rigid',
                            outprefix=opj(pet_prepro, 'tmp', 'section2_to_section1_'))

    IM1 = ants.image_read(img1)
    IM2   = ants.image_read(img2)

    ref_mean = ants.apply_transforms(fixed=IM1, moving=IM2,transformlist=mc1['fwdtransforms'],
                                interpolator='linear')

    cmd = (sing_afni + '3dresample -overwrite -master ' + opj(pet_prepro,'tmp', 'crop_ref.nii.gz') +
           ' -prefix ' + opj(pet_prepro,'tmp', Name + '_desc-crop_pet.nii.gz') +
           ' -input '  + dyn_img)
    run_cmd.run(cmd,diary_name)

    img_tmp = ants.image_read(opj(pet_prepro,'tmp', Name + '_desc-crop_pet.nii.gz'))
    tmp_unmerged = ants.ndimage_to_list(img_tmp)

    To_be_corrected = ants.image_read(dyn_img)
    img_unmerged    = ants.ndimage_to_list(To_be_corrected)

    motion_corrected = list()
    for i in range(len(img_unmerged)):
        if i > start[0] & i<start[1]:
            mc = ants.registration(fixed=ref_1, moving=tmp_unmerged[i],
                                   type_of_transform='Rigid',
                                   outprefix=opj(pet_prepro, 'tmp', Name + '_mc_' + str(i) + '_'))

            moved = ants.apply_transforms(fixed=ref_mean, moving=img_unmerged[i],
                                          transformlist=opj(pet_prepro, 'tmp', Name + '_mc_' + str(i) + '_0GenericAffine.mat'),
                                          interpolator='linear',
                                          whichtoinvert=[False])
        elif i > start[1]:
            mc = ants.registration(fixed=ref_1, moving=tmp_unmerged[i],
                                   type_of_transform='Rigid',initial_transform=mc1['fwdtransforms'],
                                   outprefix=opj(pet_prepro, 'tmp', Name + '_mc_' + str(i) + '_'))

            moved = ants.apply_transforms(fixed=ref_mean, moving=img_unmerged[i],
                                          transformlist=opj(pet_prepro, 'tmp', Name + '_mc_' + str(i) + '_0GenericAffine.mat'),
                                          interpolator='linear',
                                          whichtoinvert=[False])
        else:
            moved = img_unmerged[i]

    motion_corrected.append(moved)
    motCorr = ants.list_to_ndimage(To_be_corrected, motion_corrected)
    ants.image_write(motCorr, opj(pet_dir, Name + '_desc-MoCo_pet.nii.gz'), ri=False)

    dictionary = {"Sources": [dyn_img,
                              opj(pet_prepro,'tmp', Name + '_desc-crop_pet.nii.gz'),
                              img1],
                  "Description": 'Motion correction (Antspy).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(pet_dir, Name + '_desc-MoCo_pet.json'), "w") as outfile:
        outfile.write(json_object)


    cmd = (sing_afni + '3dTstat -overwrite -mean -prefix ' + opj(pet_dir, Name + '_desc-MoCo-mean_pet.nii.gz') +
           ' ' + opj(pet_prepro, Name + '_desc-MoCo_pet.nii.gz'))
    run_cmd.run(cmd,diary_name)

    dictionary = {"Sources": opj(pet_prepro,'tmp', Name + '_desc-crop_pet.nii.gz'),
                  "Description": 'Mean image (3dTstat, AFNI).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(pet_dir, Name + '_desc-MoCo-mean_pet.json'), "w") as outfile:
        outfile.write(json_object)
