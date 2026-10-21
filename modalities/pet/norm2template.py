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

from Tools import run_cmd,filename,setTransfo


T1interp    = 'hammingWindowedSinc'
Boldinterp  = 'lanczosWindowedSinc'
labelinterp = 'genericLabel'
petinterp   = 'linear'

antsaff     = '0GenericAffine.mat'
antswarp    = '1Warp.nii.gz'
antsinvwarp = '1InverseWarp.nii.gz'

tounwarp     = '_'.join(['bold','to','unwarp'])
tounwarptype = 'SynCC'
resamp_interp = 'Cu'

tmp_trans = 'translation'
transPET  = 'from-Pet-to-T1w'

shiftname = 'shift'
acpcname = 'acpc'

def pet(img_file,template_name,
         Ref_native_template_mask,Ref_native_template,
         REF_Brain,ref_mask,
         maxval,
         dir_transfo,
         diary_name,sing_wb,same_space):
    

    Dir_pet = opd(img_file)

    fname,suffix = filename.get_desc(opb(img_file),'desc')

    if opi(opj(dir_transfo, '_'.join(['anat2' + template_name, 'Final', antswarp]))):
        if opi(opj(dir_transfo, '_'.join(['anat2' + template_name, 'Final', antsaff]))):
            T2 = [opj(dir_transfo, '_'.join(['anat2' + template_name, 'Final', antswarp])),
                  opj(dir_transfo, '_'.join(['anat2' + template_name, 'Final', antsaff]))]
        else:
            T2 = [opj(dir_transfo, '_'.join(['anat2' + template_name, 'Final', antswarp]))]

    else:
        print('sorry impossible to find ' + opj(dir_transfo, '_'.join(['anat2' + template_name, 'Final', antsinvwarp])))

    nl = 'Normalize to template'
    run_cmd.msg(nl, diary_name,'OKGREEN')
    if str(maxval)=='':
        val = '2.5'
    else:
        val = str(maxval)    

    Palette = (' MODE_USER_SCALE -thresholding THRESHOLD_TYPE_OFF THRESHOLD_TEST_SHOW_INSIDE 0 2.5' +
              ' -pos-user 0 ' + val + ' -interpolate true -palette-name ROY-BIG -disp-pos true -disp-neg false -disp-zero false')

    for fneeded in [img_file,Ref_native_template,REF_Brain,ref_mask,
                    opj(dir_transfo,'_'.join([acpcname,antsaff])),T2[0]]:
        if not opi(fneeded):
            raise FileNotFoundError(f"Input file not found: {fneeded}")


    img = ants.image_read(img_file)
    im_head = ants.image_header_info(img_file)
    
    if len(im_head['dimensions'])>3:
      dim4 = 3
    else:
      dim4 = 0

    if not same_space:
        Tpet = opj(Dir_pet, '_'.join([fname,transPET,antsaff]))

        if not opi(Tpet):
            name1, _ = filename.get_desc(opb(img_file), 'run')
            Tpet = opj(Dir_pet, '_'.join([name1 ,transPET,antsaff]))


        if opi(opj(dir_transfo,'_'.join([shiftname,antsaff]))):
            T1 = [opj(dir_transfo,'_'.join([acpcname,antsaff])),
                  opj(dir_transfo,'_'.join([shiftname,antsaff])),
                  Tpet]
        else:
            T1 = [opj(dir_transfo,'_'.join([acpcname,antsaff])),
                  Tpet]

    else:
        if opi(opj(dir_transfo,'_'.join([shiftname,antsaff]))):
            T1 = [opj(dir_transfo,'_'.join([acpcname,antsaff])),
                  opj(dir_transfo,'_'.join([shiftname,antsaff]))]
        else:
            T1 = [opj(dir_transfo,'_'.join([acpcname,antsaff]))]


    for brain_file, msk_file,transfolist,space_name,folder in zip ([Ref_native_template,REF_Brain],
                                                                   [Ref_native_template_mask,ref_mask],
                                                                   [T1,T2 + T1],
                                                                   ['acpc', template_name],
                                                                   [opj(Dir_pet), opj(Dir_pet,'templates')]):

        brain      = ants.image_read(brain_file)
        mask       = ants.image_read(msk_file)
        foldername = opj(folder,space_name)

        file1 = opj(foldername, '_'.join([fname,'space-' + space_name,suffix,'pet.nii.gz']))
        file2 = opj(foldername, '_'.join([fname, 'space-' + space_name, suffix + '-SS', 'pet.nii.gz']))

        if not ope(foldername):
            os.makedirs(foldername)

        moved = ants.apply_transforms(fixed=brain, moving=img,transformlist=transfolist,interpolator=petinterp,imagetype=dim4)
        ants.image_write(moved, file1)

        cmd = (sing_wb + 'wb_command -volume-palette ' + file1 + Palette)
        run_cmd.run(cmd,diary_name)

        dictionary = {"Sources": [img_file,
                                  brain_file],
                      "Description": 'Co-registration of the PET image with a template image (Antspy).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(file1.replace('.nii.gz','.json'), "w") as outfile:
            outfile.write(json_object)

        if len(im_head['dimensions'])==4:
            moved_unmerged = ants.ndimage_to_list(moved)
            moved_SS_list = list()

            for i in range(len(moved_unmerged)):
                temp = moved_unmerged[i] * mask
                moved_SS_list.append(temp)

            moved_SS = ants.list_to_ndimage(moved, moved_SS_list)
        else:
            moved_SS = moved * mask

        ants.image_write(moved_SS,file2)
        cmd = (sing_wb + 'wb_command -volume-palette ' + file2 + Palette)
        run_cmd.run(cmd,diary_name)

        dictionary = {"Sources": [file1,
                                  msk_file],
                      "Description": 'Skull-stripping (Antspy).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(file2.replace('.nii.gz','.json'), "w") as outfile:
            outfile.write(json_object)

    nl ='Done'
    run_cmd.msg(nl, diary_name,'OKGREEN')
