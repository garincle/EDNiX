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

def pet(img_file,template_name,
         Ref_native_template_mask,Ref_native_template,
         REF_Brain,ref_mask,
         maxval,
         dir_transfo,
         diary_name,sing_wb):
    
    E=-1
    Dir_pet = opd(img_file)
    N = opb(img_file).split('_')
    
    E = [i for i, elem in enumerate(N) if 'desc' in elem]
    Name = '_'.join(N[0:E[0]])
    suffix = N[E[0]]
    

    nl = 'Normalize to template'
    run_cmd.msg(nl, diary_name)

    Palette = (' MODE_USER_SCALE -thresholding THRESHOLD_TYPE_OFF THRESHOLD_TEST_SHOW_INSIDE 0 2.5' +
              ' -pos-user 0 ' + maxval + ' -interpolate true -palette-name ROY-BIG -disp-pos true -disp-neg false -disp-zero false')


    img = ants.image_read(img_file)
    im_head = ants.image_header_info(img_file)
    
    if len(im_head['dimensions'])>3:
      dim4 = 3
    else:
      dim4 = 0

    Tpet = opj(Dir_pet, Name + '_from-Pet-to-T1w_0GenericAffine.mat')
    
    if not opi(Tpet):
        R = [i for i, elem in enumerate(N) if 'run' in elem]
        name1= '_'.join(N[0:R[0]])
        Tpet = opj(Dir_pet, name1 + '_from-Pet-to-T1w_0GenericAffine.mat')


    if opi(opj(dir_transfo, 'shift_0GenericAffine.mat')):
        T1 = [opj(dir_transfo, 'acpc_0GenericAffine.mat'),
              opj(dir_transfo, 'shift_0GenericAffine.mat'),
              Tpet]
    else:
        T1 = [opj(dir_transfo, 'acpc_0GenericAffine.mat'),
              Tpet]

    T2 = [opj(dir_transfo,'anat2' + template_name + '_SyN_1Warp.nii.gz'),
          opj(dir_transfo,'anat2' + template_name + '_SyN_0GenericAffine.mat')]


    for brain_file, msk_file,transfolist,space_name,folder in zip ([Ref_native_template,REF_Brain],
                                                                   [Ref_native_template_mask,ref_mask],
                                                                   [T1,T2 + T1],
                                                                   ['acpc', template_name],
                                                                   [opj(Dir_pet), opj(Dir_pet,'templates')]):

        brain = ants.image_read(brain_file)
        mask = ants.image_read(msk_file)
        foldername = opj(folder,space_name)

        if not ope(foldername):
            os.makedirs(foldername)

        moved = ants.apply_transforms(fixed=brain, moving=img,
                                      transformlist=transfolist,interpolator='linear',imagetype=dim4)
        ants.image_write(moved, opj(foldername, Name + '_space-' + space_name + '_' + suffix + '_pet.nii.gz'),ri=False)

        cmd = (sing_wb + 'wb_command -volume-palette' +
           ' ' + opj(foldername, Name + '_space-' + space_name + '_' + suffix + '_pet.nii.gz') + Palette)
        run_cmd.run(cmd,diary_name)

        dictionary = {"Sources": [img_file,
                                  brain_file],
                      "Description": 'Co-registration of the PET image with a template image (Antspy).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(foldername, Name + '_space-' + space_name + '_' + suffix + '_pet.json'), "w") as outfile:
            outfile.write(json_object)

        if len(im_head['dimensions'])==4:
            moved_unmerged = ants.ndimage_to_list(moved)
            moved_SS_list = list()
            # beware : here you need to reset orientation, spacing and origin
            # image = ants.make_image((*mask_templ.shape,int(im_head['dimensions'][3])))
            for i in range(len(moved_unmerged)):
                temp = moved_unmerged[i] * mask
                moved_SS_list.append(temp)

            moved_SS = ants.list_to_ndimage(moved, moved_SS_list)
        else:
            moved_SS = moved * mask

        ants.image_write(moved_SS,opj(foldername, Name + '_space-' + space_name + '_' + suffix + '-SS_pet.nii.gz'),ri=False)
        cmd = (sing_wb + 'wb_command -volume-palette' +
               ' ' + opj(foldername, Name + '_space-' + space_name + '_' + suffix + '-SS_pet.nii.gz') + Palette)
        run_cmd.run(cmd,diary_name)

        dictionary = {"Sources": [opj(foldername, Name + '_space-' + space_name + '_' + suffix + '_pet.nii.gz'),
                                  msk_file],
                      "Description": 'Skull-stripping (Antspy).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(foldername, Name + '_space-' + space_name + '_' + suffix + '-SS_pet.json'), "w") as outfile:
            outfile.write(json_object)

    nl ='Done'
    run_cmd.msg(nl, diary_name)
