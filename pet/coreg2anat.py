#import
import os
import sys
import json
import ants

opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

sys.path.insert(1, '../tools')
import run_cmd

def doit(IMG,
         Ref_native_raw,Ref_native_template,Ref_native_mask,template_name,atlas_file,atlas_name,dir_transfo,
         Need_check,diary_name,sing_afni,sing_fs):

    pet_dir = opd(IMG)
    pet_prepro=opj(pet_dir,'preprocessing')
    N = opb(IMG).split('_')
    for i in range(len(N)):
        a = N[i].split('-')
        if a[0] == 'desc':
            E = i
    Name = '_'.join(N[0:E])

    # ..........................................................................................................

    nl = 'Co-registration with the T1w MRI anatomy'
    run_cmd.msg(nl, diary_name)

    cmd=(sing_afni + '3dresample -master ' + Ref_native_raw + ' -rmode Cu -prefix ' + opj(pet_prepro,'tmp','pet2T1w.nii.gz') +
       ' -input ' +  IMG)
    run_cmd.run(cmd,diary_name)

    dictionary = {"Sources": [IMG,
                            Ref_native_raw],
                "Description": 'Resampling pet image to T1w (3dresample (AFNI)).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(pet_prepro,'tmp','pet2T1w.json'), "w") as outfile:
        outfile.write(json_object)

    mPet  = ants.image_read(opj(pet_prepro,'tmp','pet2T1w.nii.gz'))
    brain = ants.image_read(Ref_native_raw)

    mTx_Pet2T1w_trans = ants.registration(fixed=brain, moving=mPet,
                                      type_of_transform = 'Translation',
                                      outprefix=opj(pet_prepro,'tmp','translation'))

    source = opj(pet_prepro,'tmp','pet2T1w.nii.gz')

    if Need_check == 1:

        moved = ants.apply_transforms(fixed=brain, moving=mPet
                                    ,transformlist=mTx_Pet2T1w_trans['fwdtransforms'],
                                    interpolator='hammingWindowedSinc')
        ants.image_write(moved, opj(pet_prepro,'tmp','pet2T1w_tmp.nii.gz'), ri=False)

        cmd = (sing_fs + 'freeview -v ' + Ref_native_raw + ':colormap=grayscale:opacity=1:visible=1 ' + 
        opj(pet_prepro,'tmp','pet2T1w_tmp.nii.gz') + ':colormap=pet:opacity=0.5:visible=1')
        run_cmd.do(cmd,diary_name)

        mPetcorrect = ants.image_read(opj(pet_prepro, 'tmp', 'pet2T1w_tmp.nii.gz'))
        source = opj(pet_prepro, 'tmp', 'pet2T1w_tmp.nii.gz')
        
        mTx_Pet2T1w_correct = ants.registration(fixed=mPetcorrect, moving=mPet, type_of_transform = 'Rigid',
                                  initial_transform=opj(pet_prepro,'tmp','translation0GenericAffine.mat'),
                                  outprefix=opj(pet_prepro,'tmp','translation'))



  
    mTx_Pet2T1w = ants.registration(fixed=brain, moving=mPet, type_of_transform = 'Rigid',
                                  initial_transform=opj(pet_prepro,'tmp','translation0GenericAffine.mat'),
                                  outprefix=opj(pet_dir,Name + '_from-Pet-to-T1w_'))


    # -----------------------------------------------------------------------------------------------------------------

    # from Pet to T1w
    nl = 'Transform the pet image into the T1w MRI anatomy'
    run_cmd.msg(nl, diary_name)

    moved = ants.apply_transforms(fixed=brain, moving=mPet ,
                                transformlist=mTx_Pet2T1w['fwdtransforms'],
                                interpolator='hammingWindowedSinc')
    ants.image_write(moved, opj(pet_dir, Name + '_desc-coreg_pet.nii.gz'), ri=False)

    dictionary = {"Sources": [source,
                            Ref_native_raw],
                "Description": 'Co-registration of the PET image with the T1w image (Antspy).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(pet_dir, Name + '_desc-coreg_pet.json'), "w") as outfile:
        outfile.write(json_object)


    # -----------------------------------------------------------------------------------------------------------------

    Mean = ants.image_read(IMG)

    if opi(opj(dir_transfo,'shift_0GenericAffine.mat')):
        T1 = [opj(pet_dir,Name + '_from-Pet-to-T1w_0GenericAffine.mat'),
            opj(dir_transfo,'shift_0GenericAffine.mat'),
            opj(dir_transfo,'acpc_0GenericAffine.mat')]
        w1 = [True,True,True]
    else :
        T1 = [opj(pet_dir, Name + '_from-Pet-to-T1w_0GenericAffine.mat'),
            opj(dir_transfo, 'acpc_0GenericAffine.mat')]
        w1 = [True, True]

    T2 = [opj(dir_transfo,'anat2' + template_name + '_SyN_0GenericAffine.mat'),
          opj(dir_transfo,'anat2' + template_name + '_SyN_1InverseWarp.nii.gz')]
    w2 = [True,False]


    # from Atlas template to pet
    nl = 'Transform the ROIs defined by the atlas to fit the PET data'
    run_cmd.msg(nl, diary_name)

    atlas = ants.image_read(atlas_file)
    moved = ants.apply_transforms(fixed=Mean, moving=atlas,
                                transformlist= T1 + T2,
                                interpolator='genericLabel',whichtoinvert= w1 + w2)
    ants.image_write(moved, opj(pet_prepro,'tmp',Name + '_atlas_tmp.nii.gz'), ri=False)

    cmd = (sing_afni + '3dresample -master ' + IMG + ' -prefix ' + opj(pet_dir, Name + '_seg-'+ atlas_name + '_dseg.nii.gz') +
           ' -input ' +  opj(pet_prepro,'tmp',Name + '_atlas_tmp.nii.gz'))
    run_cmd.run(cmd,diary_name)

    dictionary = {"Sources": [source,
                            atlas_file],
                "Description": 'Co-registration of the labelled image with the PET image (Antspy).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(pet_dir, Name + '_seg-' + atlas_name + '_dseg.json'), "w") as outfile:
        outfile.write(json_object)


    # from acpc space to pet

    for old_img,interp, tmp,name,in zip([Ref_native_template,Ref_native_mask],
                                        ['hammingWindowedSinc','nearestNeighbor'],
                                        ['T1w','brain mask'],
                                        ['coreg_T1w','coreg_mask']):

        nl = 'Transform the ' + tmp + ' image to fit the PET data'
        run_cmd.msg(nl, diary_name)

        img = ants.image_read(old_img)

        moved = ants.apply_transforms(fixed=Mean, moving=img,transformlist= T1,interpolator=interp,whichtoinvert=w1)
        ants.image_write(moved, opj(pet_prepro,'tmp','tmp.nii.gz'), ri=False)

        cmd = (sing_afni + '3dresample -master ' + IMG + ' -prefix ' + opj(pet_dir, Name + '_desc-' + name + '.nii.gz') +
            ' -input ' +  opj(pet_prepro,'tmp','tmp.nii.gz'))
        run_cmd.run(cmd,diary_name)

        dictionary = {"Sources": [source,
                                  old_img],
                    "Description": 'Co-registration of the ' + tmp + ' image with the PET image (Antspy).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(pet_dir, Name + '_desc-' + name + '.json'), "w") as outfile:
            outfile.write(json_object)

