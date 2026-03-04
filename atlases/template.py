#import
import os
import json
import numpy as np
import ants

opj = os.path.join
opi = os.path.isfile
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists

from Tools import run_cmd, Load_EDNiX_requirement,get_orientation
from atlases import templatefeat
from anat.freesurfer import smallbrain, preFS, FS_surf,  FS_2_WB_short, FS_flat


Hmin    = ['l','r']
species =  'Pig'
reference = 'EDNiX'
atlas_followers=[['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]]
MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
bids_dir = '/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/ungulates/Pig/'
diary_file = '/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/ungulates/Pig/EDNiX/diary.txt'

def modif(species,reference,atlas_followers,MAIN_PATH,bids_dir,diary_file,proj='ribbon',doflat=0):


    # reference templates
    if reference == '':
        reference = 'EDNiX'

    reftemplate_path = opj(opd(MAIN_PATH), "Atlases_library")
    fs_tools         = opj(reftemplate_path, 'freesurfer_refs')
    path_BALSA       = opj(reftemplate_path, 'BALSA')

    ### singularity set up
    sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _, sing_synstrip, Unetpath = Load_EDNiX_requirement.load_requirement(
        MAIN_PATH, reftemplate_path, bids_dir, 'yes')


    (FS_refs,template_dir,reference,
                balsa_folder,BALSAname,balsa_brainT1,
                BASE_folder,BASE_atlas_folder, BASE_template,BASE_SS, BASE_mask, BASE_Gmask, BASE_Wmask, BASE_Vmask,
                CSF, GM, WM, Aseg_ref,
                list_atlas, path_label_code) = templatefeat.get(species, reftemplate_path, fs_tools, path_BALSA, reference,
                                                     '', '','T1w','anat', atlas_followers)

    FS_dir     = opj(template_dir,reference,'freesurfer')
    if not ope(opj(FS_dir,species)):
        os.makedirs(FS_dir)
        os.makedirs(opj(FS_dir,species,'mri'))
        os.makedirs(opj(FS_dir,species,'surf'))
        os.makedirs(opj(FS_dir,species,'stats'))
        os.makedirs(opj(FS_dir,species, 'label'))
        os.makedirs(opj(FS_dir,species,  'scripts'))
    dir_prepro = opj(template_dir,reference,'preprocessing')
    if not ope(dir_prepro):
        os.makedirs(dir_prepro)

    cmd_tksurfer, cmd_flatten, cmd_mris, export_fs = Load_EDNiX_requirement.FS(fs_tools, FS_dir, diary_file, sing_fs)

    #loop for specie

    # modifications of the atlas to T1 already done

    # a little help for segmentation
    preFS.prepa_img(species, BASE_SS, dir_prepro, Aseg_ref, BASE_atlas_folder, '', diary_file)

    list1 = [BASE_SS,
             '_'.join([species, 'desc-norm', 'T1w']),
             '_'.join([species, 'seg-wm', 'dseg']),
             '_'.join([species, 'seg-4FS', 'dseg']),
             '_'.join([species, 'seg-filled', 'dseg'])]

    list2 = ['orig.mgz',
             'brain.mgz',
             'wm.mgz',
             'aseg.mgz',
             'filled.mgz']

    if species in ['Mouse','Rat','Mouselemur','Bat']:
        change_hd, scaling, resamp, new_size, DATfile = smallbrain.get(species, dir_prepro, species, '', diary_file)

        # set images to resolution 1mm

        brain_img = ants.image_read(list[0])

        X = brain_img.origin
        S = brain_img.spacing

        new_orig = np.zeros(3)
        new_orig[0] = X[0] * scaling[0] / S[0]
        new_orig[1] = X[1] * scaling[1] / S[1]
        new_orig[2] = X[2] * scaling[2] / S[2]

        for new_img, dir_file, interp in zip(list,
                                             ['',dir_prepro,BASE_atlas_folder,BASE_atlas_folder,BASE_atlas_folder],
                                             [3, 3, 1, 1, 1]):
            print(new_img)
            if new_img == list[0]:
                resamp_img = ants.image_read(new_img)
            else:
                resamp_img = ants.image_read(opj(dir_file, new_img + '.nii.gz'))
            ants.set_spacing(resamp_img, scaling)
            ants.set_origin(resamp_img, [new_orig[0], new_orig[1], new_orig[2]])
            if resamp == 1:
                resamp_img = ants.resample_image(resamp_img, (new_size, new_size, new_size), False, interp)
                resamp_img[resamp_img < 0] = 0
                if 'desc-norm' in new_img:
                    resamp_img = ants.smooth_image(resamp_img, new_size)
                descript = 'reformat image with a voxel size of ' + str(new_size) + 'and a resampling of ' + str(
                    scaling)
            else:
                descript = 'reformat image with a resampling of ' + str(scaling)
            if new_img == list[0]:
                newname = opj(dir_prepro, opb(new_img).replace('.nii.gz', '_resamp-' + str(scaling[0]) + '.nii.gz'))
                print(newname)
            else:
                newname = opj(dir_prepro, new_img + '_resamp-' + str(scaling[0]) + '.nii.gz')
                print(newname)
            ants.image_write(resamp_img, newname)
            dictionary = {"Sources": [opj(dir_file, new_img),
                                      list[0]],
                          "Description": descript + ' (Antspy)'}
            json_object = json.dumps(dictionary, indent=2)
            with open(newname.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)

        nl = 'Done'
        run_cmd.msg(nl, diary_file, 'OKGREEN')

        Ref = '_resamp-' + str(scaling)
        list1[0] = opj(dir_prepro, opb(list1[0].replace('.nii.gz', Ref + '.nii.gz')))
        labels_dir = dir_prepro
    else:
        Ref=''
        labels_dir = BASE_atlas_folder




    list1[1] = opj(dir_prepro, list1[1] + Ref + '.nii.gz')
    for i in range(2, 5):
        list1[i] = opj(labels_dir, list1[i] + Ref + '.nii.gz')

    [_, _, fwdFS_cmd, _] = get_orientation.use_ants(list1[0])

    # set the data and folder for freesurfer
    for i in range(len(list1)):
        if opi(list1[i])==True:
            cmd = sing_fs + 'mri_convert ' + fwdFS_cmd + list1[i] + ' ' + opj(FS_dir,species, 'mri', list2[i])
            run_cmd.run(cmd,diary_file)

    nl = 'Done'
    run_cmd.msg(nl, diary_file,'OKGREEN')

    FS_surf.Wcreate(FS_dir, species, 'lr', diary_file, sing_fs, export_fs)

    FS_surf.Wconfig(FS_dir, FS_refs, species, 'lr', diary_file, sing_fs, export_fs)

    if doflat ==0:
        run_cmd.msg('INFO: skip step ' + str('flat_map'), diary_file, 'OKGREEN')
    else:
        # Flat maps
        #### 4.4.1 Left hemisphere
        FS_flat.map(FS_dir, species, 'left', cmd_tksurfer, cmd_flatten, diary_file, sing_fs, export_fs)
        #### 4.4.2 Right hemisphere
        FS_flat.map(FS_dir, species, 'right', cmd_tksurfer, cmd_flatten, diary_file, sing_fs, export_fs)

    FS_surf.Pcreate(species, 'lr', diary_file, export_fs)

    cmd = export_fs + 'mris_volmask --aseg_name aseg --label_left_white 2 --label_left_ribbon 3 --label_right_white 41 --label_right_ribbon 42 --save_ribbon ' + species
    run_cmd.do(cmd, diary_file)

    os.remove(opj(FS_dir, species, 'mri', 'lh.ribbon.mgz'))
    os.remove(opj(FS_dir, species, 'mri', 'rh.ribbon.mgz'))

    for H in range(2):

        # 2) Midthickness surface
        nl = 'Creation of the Midthickness surface'
        run_cmd.msg(nl, diary_file, 'OKGREEN')

        cmd = sing_fs + 'mris_expand -thickness ' + opj(FS_dir,species,'surf', Hmin[H] + 'h.white') + ' 0.5 ' + opj(FS_dir,species,'surf',Hmin[H] + 'h.mid')
        run_cmd.run(cmd,diary_file)

    FS_2_WB_short.WB_prep(cmd_mris, FS_dir, FS_refs, species, species, BASE_SS, DATfile, reference, list_atlas, balsa_folder,
            BALSAname,
            path_label_code, template_dir, proj, diary_file, sing_fs, sing_wb, export_fs)