#import
import os
import json
import ants
import shutil

opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd,extract_filename
from anat.freesurfer import FS_freeview

T1interp    = 'hammingWindowedSinc'
Boldinterp  = 'lanczosWindowedSinc'
labelinterp = 'genericLabel'

antsaff     = '0GenericAffine.mat'
antswarp    = '1Warp.nii.gz'
antsinvwarp = '1InverseWarp.nii.gz'

resamp_interp = 'Cu'

tmp_trans = 'translation'
transPET  = 'from-Pet-to-'
transCT   = 'from-ct-to-'
acpcname  = 'acpc'
shiftname = 'shift'


def doit(IMG,Ref_native_raw,T_PET,Ref_native_template,Ref_native_mask,Ref_native_whitemask,Ref_nativegray_mask,template_name,
         atlas_file,atlas,dir_transfo,dir_prepro_orig_process, dir_prepro_orig_labels, dir_prepro_orig_masks,
         dir_prepro_orig_matrices,anat_pet_same_space,do_manual,IhaveanANAT,IhaveaCT,dir_prepro_orig_ct_matrices,dir_prepro_acpc,
         dir_prepro_acpc_process,dir_prepro_acpc_label,dir_prepro_acpc_masks,dir_prepro_template,
         dir_prepro_template_labels, dir_prepro_template_masks,BASE_template,
         dir_prepro_template_process,BASE_mask, BASE_Gmask, BASE_Wmask,
         Need_check_coreg,diary_name,sing_afni,sing_fs):

    transPET = transPET + T_PET
    transCT  = transCT + T_PET

    tmp_dir = opj(dir_prepro_orig_process, 'forMC')
    fname, _ = extract_filename.get_desc(opb(IMG), 'desc')

    # pet image in T1w space (only for visualisation)
    step1        = opj(tmp_dir, 'pet2T1w.nii.gz')
    petanat      = opj(dir_prepro_orig_process, '_'.join([fname, 'space-anat', 'pet.nii.gz']))
    petacpc      = opj(dir_prepro_acpc_process, '_'.join([fname, 'space-acpc', 'pet.nii.gz']))
    pet_template = opj(dir_prepro_template_process, '_'.join([fname, 'space-' + template_name, 'pet.nii.gz']))

    imghd = ants.image_header_info(IMG)
    if imghd['nDimensions'] == 4:
        full_IMG = IMG
        IMG = IMG.replace('_pet.nii.gz','-mean_pet.nii.gz')
        full_petanat = petanat
        full_petacpc = petacpc
        full_pet_template = pet_template
        petanat = opj(dir_prepro_orig_process, '_'.join([fname, 'space-anat','desc-mean' ,'pet.nii.gz']))
        petacpc = opj(dir_prepro_acpc_process, '_'.join([fname, 'space-acpc', 'desc-mean','pet.nii.gz']))
        pet_template = opj(dir_prepro_template_process, '_'.join([fname, 'space-' + template_name,'desc-mean', 'pet.nii.gz']))

    crop_pet = opj(tmp_dir, '_'.join([fname, 'desc-mean-crop', 'pet.nii.gz']))


    # template atlas to pet image
    atlas2         = opj(dir_prepro_orig_labels, '_'.join([fname, 'seg-' + atlas["name"], 'dseg.nii.gz']))
    atlas_acpc     = opj(dir_prepro_acpc_label, '_'.join([fname, 'seg-' + atlas["name"], 'dseg.nii.gz']))
    atlas_template = opj(dir_prepro_template_labels, '_'.join([fname, 'seg-' + atlas["name"], 'dseg.nii.gz']))

    # anat image in pet space (only for visualisation)
    anat2petfile  = opj(dir_prepro_orig_process, '_'.join([fname, 'space-orig',T_PET + '.nii.gz']))

    mask2petfile  = opj(dir_prepro_orig_masks,   '_'.join([fname,'space-orig','mask.nii.gz']))
    white2petfile = opj(dir_prepro_orig_masks,   '_'.join([fname, 'space-orig','desc-white', 'mask.nii.gz']))
    gray2petfile  = opj(dir_prepro_orig_masks,   '_'.join([fname, 'space-orig','desc-gray', 'mask.nii.gz']))

    mask2petfile_acpc = opj(dir_prepro_acpc_masks, '_'.join([fname, 'space-acpc', 'mask.nii.gz']))
    white2petfile_acpc = opj(dir_prepro_acpc_masks, '_'.join([fname, 'space-acpc', 'desc-white', 'mask.nii.gz']))
    gray2petfile_acpc = opj(dir_prepro_acpc_masks, '_'.join([fname, 'space-acpc', 'desc-gray', 'mask.nii.gz']))

    mask2petfile_template = opj(dir_prepro_template_masks, '_'.join([fname, 'space-' + template_name, 'mask.nii.gz']))
    white2petfile_template = opj(dir_prepro_template_masks, '_'.join([fname, 'space-' + template_name, 'desc-white', 'mask.nii.gz']))
    gray2petfile_template = opj(dir_prepro_template_masks, '_'.join([fname, 'space-' + template_name, 'desc-gray', 'mask.nii.gz']))



    # get the "native to template" transformation matrices

    orderT1 = [opj(dir_transfo, '_'.join([shiftname, antsaff])),
               opj(dir_transfo, '_'.join([acpcname, antsaff]))]
    orderW1 = [True, True]

    T1 = []
    w1 = []
    for i in range(len(orderT1)):
        if opi(orderT1[i]):
            T1.append(orderT1[i])
            w1.append(orderW1[i])

    orderT2= [opj(dir_transfo, '_'.join(['native','to',template_name,'Final',antsaff])),
              opj(dir_transfo, '_'.join(['native','to',template_name,'Final',antsinvwarp]))]
    orderW2 = [True, False]

    T2 = []
    w2 = []
    for i in range(len(orderT2)):
        if opi(orderT2[i]):
            T2.append(orderT2[i])
            w2.append(orderW2[i])

    if IhaveaCT:
        pet_to_anat = opj(dir_prepro_orig_ct_matrices, '_'.join([fname, transCT, antsaff]))
    else:
        pet_to_anat = opj(dir_prepro_orig_matrices, '_'.join([fname, transPET, antsaff]))

    list_file = [IMG, Ref_native_raw, Ref_native_template, Ref_native_mask, atlas_file,
                 opj(dir_transfo, '_'.join([acpcname, antsaff])), T2[0]]

    for fneeded in list_file:
        if not opi(fneeded):
            raise FileNotFoundError(f"Input file not found: {fneeded}")


    # ..........................................................................................................

    nl = 'Co-registration with the MRI anatomy'
    run_cmd.msg(nl, diary_name,'OkGREEN')

    # from Pet to T1w

    if not opi(crop_pet):
        os.mkdirs(tmp_dir,exist_ok = True)

        shutil.copyfile(IMG, crop_pet)
        cmd = (sing_fs + 'freeview -v ' + IMG + ' -subtitle "CROP-THE-IMAGE ! "')
        run_cmd.do(cmd, diary_name)

        dictionary = {"Sources": IMG,
                      "Description": 'Manually Cropped image (Freeview, Freesurfer).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(crop_pet.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)



    nl = 'Transform the pet image into the MRI anatomy'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    # Crop the PET image around the skull
    cmd = (sing_afni + '3dresample -overwrite -master ' + crop_pet + ' -rmode ' + resamp_interp +
           ' -prefix ' + step1 + ' -input ' + IMG)
    run_cmd.run(cmd, diary_name)

    # Resample to anat resolution
    cmd=(sing_afni + '3dresample -overwrite -master ' + Ref_native_raw + ' -rmode ' + resamp_interp +
         ' -prefix ' + step1 + ' -input ' +  step1)
    run_cmd.run(cmd,diary_name)

    dictionary = {"Sources": [IMG,
                            Ref_native_raw],
                "Description": 'Resampling pet image to T1w (3dresample (AFNI)).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(step1.replace('.nii.gz','.json'), "w") as outfile:
        outfile.write(json_object)


    Pet       = ants.image_read(IMG)
    brain = ants.image_read(Ref_native_raw)




    if IhaveanANAT:
        if anat_pet_same_space:
            shutil.copyfile(step1,petanat)
            shutil.copyfile(step1.replace('.nii.gz','.json'),
                            petanat.replace('.nii.gz','.json'))

            # Anatomical file into pet file
            cmd = (sing_afni + '3dresample -overwrite -master ' + IMG + ' -rmode ' + resamp_interp +
                   ' -prefix ' + anat2petfile + ' -input ' + Ref_native_raw)
            run_cmd.run(cmd, diary_name)

            dictionary = {"Sources": [IMG,
                                      Ref_native_raw],
                          "Description": 'Resampling pet image to T1w (3dresample (AFNI)).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(step1.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)

        else:
            if not IhaveaCT:
                ants.registration(fixed=brain, moving=Pet,
                                  type_of_transform='Translation',
                                  outprefix=opj(dir_prepro_orig_matrices, tmp_trans + '_'))
                ants.registration(fixed=brain, moving=Pet, type_of_transform='Rigid',
                                  initial_transform=opj(tmp_dir, '_'.join([tmp_trans, antsaff])),
                                  outprefix=opj(dir_prepro_orig_matrices, '_'.join([fname, transPET, ''])))

                if do_manual:
                    mPet = ants.apply_transforms(fixed=brain, moving=Pet, transformlist=[pet_to_anat], interpolator=Boldinterp)
                    refacpc = ants.image_read(Ref_native_template)

                    T1b = []
                    for i in range(len(orderT1)):
                        if opi(orderT1[i]):
                            T1b.append(orderT1[i])

                    moved = ants.apply_transforms(fixed=refacpc, moving=mPet,transformlist=T1b,
                                                  interpolator=Boldinterp)
                    ants.image_write(moved, opj(dir_prepro_orig_matrices, 'dummy.nii.gz'))

                    FS_freeview.check_coreg(opj(dir_prepro_orig_matrices, 'dummy.nii.gz'), Ref_native_template, 'acpc','', 'pet',
                                            sing_fs)

                    dummy = ants.image_read(opj(dir_prepro_orig_matrices, 'dummy.nii.gz'))

                    mv = ants.apply_transforms(fixed=mPet, moving=dummy,
                                               transformlist=[opj(dir_transfo, '_'.join([acpcname, antsaff]))],
                                               interpolator=Boldinterp, whichtoinvert=[True])

                    ants.registration(fixed=mv, moving=Pet,
                                      type_of_transform='Translation',
                                      outprefix=opj(dir_prepro_orig_matrices, tmp_trans + '_'))
                    ants.registration(fixed=mv, moving=mPet,
                                      type_of_transform='Rigid',
                                      initial_transform=opj(dir_prepro_orig_matrices, '_'.join([tmp_trans, antsaff])),
                                      outprefix=opj(dir_prepro_orig_matrices, '_'.join([fname, transPET, ''])))

                    os.remove(opj(dir_prepro_orig_matrices, 'dummy.nii.gz'))


            moved = ants.apply_transforms(fixed=brain, moving=Pet, transformlist=[pet_to_anat], interpolator=Boldinterp)
            ants.image_write(moved, petanat)
            dictionary = {"Sources": [IMG,
                                      Ref_native_raw],
                          "Description": 'Co-registration of the PET image with the anatomical image (Antspy).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(petanat.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)

            moved = ants.apply_transforms(fixed=Pet, moving=brain, transformlist=[pet_to_anat], interpolator=Boldinterp,
                                          whichtoinvert=[True])
            ants.image_write(moved, anat2petfile)
            dictionary = {"Sources": [IMG,
                                      Ref_native_raw],
                          "Description": 'Co-registration of the anatomical image with the PET image (Antspy).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(anat2petfile.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)


        # Atlas and masks (from template space and acpc space):
        for old_img, tmp, file in zip([atlas_file, Ref_native_mask, Ref_native_whitemask, Ref_nativegray_mask],
                                      ['atlas_file', 'brain_mask', 'white_mask', 'gray_mask'],
                                      [atlas2, mask2petfile, white2petfile, gray2petfile]):
            if opi(old_img):
                nl = 'Transform the ' + tmp + ' image to fit the PET data'
                run_cmd.msg(nl, diary_name, 'OKGREEN')

                dummy_img = opj(dir_prepro_orig_matrices, 'tmp.nii.gz')
                img = ants.image_read(old_img)

                imd = 0
                if img.dimension > 3:
                    imd = 3

                Tf = []
                wf = []

                if opi(pet_to_anat):
                    Tf.append(pet_to_anat)
                    wf.append(True)

                Tf += T1
                wf += w1
                if old_img == atlas_file:
                    Tf += T2
                    wf += w2

                moved = ants.apply_transforms(fixed=Pet, moving=img, transformlist=Tf, interpolator=labelinterp,
                                              whichtoinvert=wf, imagetype=imd)
                ants.image_write(moved, dummy_img, ri=False)

                cmd = (sing_afni + '3dresample -rmode NN -overwrite -master ' + IMG +
                       ' -prefix ' + file + ' -input ' + dummy_img)
                run_cmd.run(cmd, diary_name)

                dictionary = {"Sources": [IMG,
                                          old_img],
                              "Description": 'Co-registration of the ' + tmp + ' image with the PET image (Antspy).', }
                json_object = json.dumps(dictionary, indent=2)
                with open(file.replace('.nii.gz', '.json'), "w") as outfile:
                    outfile.write(json_object)
                os.remove(dummy_img)


            if Need_check_coreg:
                FS_freeview.check_coreg(petanat, Ref_native_raw, 'acpc', '', 'pet', sing_fs)



            new_dir = [dir_prepro_acpc,dir_prepro_acpc_process,dir_prepro_acpc_label,dir_prepro_acpc_masks,dir_prepro_template, dir_prepro_template_labels, dir_prepro_template_masks,
         dir_prepro_template_process]
            for dir in new_dir:
                os.makedirs(dir, exist_ok=True)

            # acpc space
            acpc_image = ants.image_read(Ref_native_template)
            # pet image
            Tf = T1
            if opi(pet_to_anat):
                Tf.append(pet_to_anat)
            moved_ref = ants.apply_transforms(fixed=acpc_image, moving=Pet, transformlist=Tf, interpolator='Linear')
            ants.image_write(moved_ref, petacpc)
            dictionary = {"Sources": [IMG,
                                      Ref_native_template],
                          "Description": 'Co-registration of the PET image with the anatomical image (Antspy).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(petacpc.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)

            if full_IMG:
                full_Pet = ants.image_read(full_IMG)
                moved = ants.apply_transforms(fixed=acpc_image, moving= full_Pet, transformlist=Tf, interpolator='Linear',imagetype=3)
                ants.image_write(moved, full_petacpc)
                dictionary = {"Sources": [full_IMG,
                                          Ref_native_template],
                              "Description": 'Co-registration of the PET image with the anatomical image (Antspy).', }
                json_object = json.dumps(dictionary, indent=2)
                with open(full_petacpc.replace('.nii.gz', '.json'), "w") as outfile:
                    outfile.write(json_object)

            # labelled images
            atlas_ref = ants.image_read(atlas_file)
            moved_atlas = ants.apply_transforms(fixed=moved_ref, moving=atlas_ref , transformlist=T2, interpolator=labelinterp)
            ants.image_write(moved_atlas, atlas_acpc)
            dictionary = {"Sources": [petacpc,
                                      atlas_file],
                          "Description": 'Co-registration of the atlas image with the pet image (Antspy).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(atlas_acpc.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)


            for img, file in zip([Ref_native_mask,Ref_native_whitemask,Ref_nativegray_mask],
                                 [mask2petfile_acpc, white2petfile_acpc,
                                  gray2petfile_acpc]):
                cmd = (sing_afni + '3dresample -overwrite -master ' + petacpc + ' -rmode NN ' +
                       ' -prefix ' + file + ' -input ' + img)
                run_cmd.run(cmd, diary_name)

                dictionary = {"Sources": [petacpc,
                                          img],
                              "Description": 'Resampled to pet image (3dresample (AFNI)).', }
                json_object = json.dumps(dictionary, indent=2)
                with open(file.replace('.nii.gz', '.json'), "w") as outfile:
                    outfile.write(json_object)



            # template space
            T_image = ants.image_read(BASE_template)
            # pet image
            Tf = T2 + T1
            if opi(pet_to_anat):
                Tf.append(pet_to_anat)
            
            moved = ants.apply_transforms(fixed=T_image, moving=Pet, transformlist=Tf, interpolator='Linear')
            ants.image_write(moved, pet_template)
            dictionary = {"Sources": [IMG,
                                      BASE_template],
                          "Description": 'Co-registration of the PET image with the anatomical image (Antspy).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(pet_template.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)
            if full_IMG:
                full_Pet = ants.image_read(full_IMG)
                moved = ants.apply_transforms(fixed=T_image, moving= full_Pet, transformlist=Tf, interpolator='Linear',imagetype=3)
                ants.image_write(moved, full_pet_template)
                dictionary = {"Sources": [full_IMG,
                                          BASE_template],
                              "Description": 'Co-registration of the PET image with the anatomical image (Antspy).', }
                json_object = json.dumps(dictionary, indent=2)
                with open(pet_template.replace('.nii.gz', '.json'), "w") as outfile:
                    outfile.write(json_object)


            # labelled images
            for img, file in zip([atlas_file,BASE_mask, BASE_Gmask, BASE_Wmask],
                                 [atlas_template,mask2petfile_template,white2petfile_template,gray2petfile_template]):
                cmd = (sing_afni + '3dresample -overwrite -master ' + pet_template + ' -rmode NN ' +
                       ' -prefix ' + file + ' -input ' + img)
                run_cmd.run(cmd, diary_name)

                dictionary = {"Sources": [pet_template,
                                          img],
                              "Description": 'Resampled to pet image (3dresample (AFNI)).', }
                json_object = json.dumps(dictionary, indent=2)
                with open(file.replace('.nii.gz', '.json'), "w") as outfile:
                    outfile.write(json_object)


        else:
            print('Sorry not ready yet....')




