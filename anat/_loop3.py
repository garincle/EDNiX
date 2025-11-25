#import
import os
import json
import nibabel as nib

opj = os.path.join
opi = os.path.isfile
opb = os.path.basename

from Tools import run_cmd, diaryfile, getpath, make_iso_img, Load_EDNiX_requirement
from atlases import templatefeat
from anat.loop3 import backtonative, _16_anat_QC_SNR, _100_Data_Clean, _200_Data_QC
from anat import transfo_T2toT1w
from anat.freesurfer import smallbrain, preFS, FS_surf, FS_finalise, FS_2_WB, FS_flat, FS_freeview

def run(all_ID, all_Session, all_data_path,ID, Session, data_path, max_ses,creat_study_template,study_template_atlas_folder,BASE_SS,BASE_mask,BASE_atlas_folder,
        coregistration_longitudinal,bids_dir,type_norm,Skip_step,check_visualy_each_img,do_fMRImasks, n_for_ANTS,
        fs_tools,sing_afni, sing_itk,sing_fs,sing_wb,fMRImasks,list_atlas,listTimage,species,
        path_label_code,FS_refs,reference,balsa_folder,BALSAname,balsa_brainT1,do_Nu,addatlas, list_transfo):


    nl = '###################################################### work on subject: ' + str(ID) + ' Session ' + str(
        Session) + ' BLOCK 3 ###############################################################'
    run_cmd.printcolor(nl, 'HEADER')

    #ID = 'sub-' + ID + '_ses-' + str(Session)

    # The anatomy
    (path_anat, dir_transfo, FS_dir, dir_prepro, dir_native,
     volumes_dir, labels_dir, masks_dir) = getpath.anat(data_path,reference,BALSAname,creat_study_template,
                                                        coregistration_longitudinal,'native')

    diary_file = diaryfile.create(opj(path_anat, str(ID) + ' session ' + str(Session)),'BLOCK3')
    cmd_tksurfer, cmd_flatten, cmd_mris, export_fs = Load_EDNiX_requirement.FS(fs_tools, FS_dir, diary_file, sing_fs)
    targetsuffix = 'space-acpc_desc-SS'


    if 8 in Skip_step:
        run_cmd.msg('INFO: skip step ' + str(8), diary_file, 'OKGREEN')
    else:
        
        ####################################################################################
        ########################## Coregistration template to anat #########################
        ####################################################################################
        
        info = backtonative.get(ID,data_path,bids_dir,Session,max_ses,targetsuffix,type_norm,BASE_SS, BASE_mask,BASE_atlas_folder,study_template_atlas_folder,
            creat_study_template,coregistration_longitudinal,reference,'Final',species,diary_file)

        targetsuffix= 'space-acpc_desc-SS-step2'
        print(info)
        backtonative.apply(ID,volumes_dir,masks_dir,labels_dir,bids_dir,info,listTimage,targetsuffix,n_for_ANTS,list_atlas,path_label_code,type_norm,diary_file,sing_afni,sing_wb)

    # ................................................................................................................
            ########################## Building fMRI masks for EPI analysis ##############################

    targetsuffix = '_'.join(['space-acpc','desc-SS'])
    brain_mask   = opj(masks_dir, '_'.join([ID,'space-acpc','mask']) + '.nii.gz')
    brain_mask_iso = opj(masks_dir, '_'.join([ID, 'space-acpc', 'res-iso', 'mask']) + '.nii.gz')
    aseg_img     = opj(labels_dir,'_'.join([ID,'seg-4FS','dseg']) + '.nii.gz')
    aseg_img_iso    = opj(labels_dir, '_'.join([ID, 'seg-4FS','res-iso', 'dseg']) + '.nii.gz')

    Ref_file     = ''
    T2_img       = ''
    T1_suffix_T  = ''
    T2_img_iso   = ''
    Ref_file_iso = ''

    for suffix_T in listTimage:
        if 'T1' in suffix_T:
            T1_suffix_T = suffix_T
            if type_norm == T1_suffix_T:
                Ref_file     = opj(volumes_dir, '_'.join([ID,targetsuffix,T1_suffix_T]) + '.nii.gz')
                Ref_file_iso = opj(volumes_dir, '_'.join([ID,targetsuffix,'res-iso',T1_suffix_T]) + '.nii.gz')
        elif 'T2' in suffix_T:
            T2_suffix_T = suffix_T
            T2_img     = opj(volumes_dir, '_'.join([ID,targetsuffix,T2_suffix_T]) + '.nii.gz')
            T2_img_iso = opj(volumes_dir, '_'.join([ID,targetsuffix,'res-iso',T2_suffix_T]) + '.nii.gz')
            if type_norm == T2_suffix_T:
                Ref_file     = T2_img
                Ref_file_iso = T2_img_iso
        else :
            nl = 'T1 or T2 is not in suffix name of your BIDS anat images, this is a requirement for EDNiX!'
            raise Exception(run_cmd.error(nl, diary_file))
    print('Ref_file is ' + str(Ref_file))

    iso = 1
    for img, img_out, type in zip([brain_mask, aseg_img, Ref_file, T2_img], [brain_mask_iso, aseg_img_iso, Ref_file_iso, T2_img_iso], ['seg', 'seg', 'anat', 'anat']):
        if opi(img):
            img_nib = nib.load(img)
            # Get voxel sizes
            delta_x, delta_y, delta_z = [str(round(abs(x), 10)) for x in img_nib.header.get_zooms()[:3]]
            if float(delta_x) == float(delta_y) == float(delta_z):
                print("Image is isotropic:", delta_x)
            else:
                print("Image is not isotropic:", delta_x, delta_y, delta_z)
                make_iso_img.make_iso(img, img_out, diary_file, sing_afni, type, ' -overwrite')
                iso = 0

    if iso == 0:
        brain_mask = brain_mask_iso
        aseg_img   = aseg_img_iso
        Ref_file   = Ref_file_iso
        T2_img     = T2_img_iso

    if 9 in Skip_step:
        run_cmd.msg('INFO: skip step ' + str(9), diary_file, 'OKGREEN')
    else:
        if do_fMRImasks == True:
            preFS.msk_RS(ID, brain_mask, aseg_img, diary_file,fMRImasks)

            #_9_do_fMRImasks.do_fMRImasks(masks_dir, labels_dir, type_norm, fMRImasks, overwrite,sing_afni, diary_file)

    # ...................................................................................................................
                    ########################## White Surface construction ##############################

    if 10 in Skip_step:
        run_cmd.msg('INFO: skip step ' + str(10), diary_file, 'OKGREEN')
    else:
        # modifications of the atlas to T1
        if 'T2' in type_norm:
            name = '_'.join(opb(Ref_file).split('_')[0:-1])
            if T1_suffix_T == '':
                print('creat fake T1 to build surfaces')
                transfo_T2toT1w.prepa(Ref_file, brain_mask, diary_file)
                new_Ref_file = opj(volumes_dir, name + '_T1w.nii.gz')
            else:
                new_Ref_file = opj(volumes_dir, name + '_' + T1_suffix_T + '.nii.gz')

        else:
            print('nothing to do, so far so good')
            new_Ref_file = Ref_file

        # a little help for segmentation
        preFS.prepa_img(ID, new_Ref_file, dir_prepro, aseg_img,labels_dir, '', diary_file)

        # need to lie on the brain volume for small brain animals
        change_hd,scaling,resamp,new_size,DATfile = smallbrain.get(ID, dir_prepro, new_Ref_file, '', diary_file,sing_fs)

        if iso == 0:
            list1 = [new_Ref_file,
                    '_'.join([ID, 'desc-norm', 'T1w']),
                    '_'.join([ID, 'seg-wm',    'dseg']),
                    '_'.join([ID, 'seg-4FS','res-iso','dseg']),
                    '_'.join([ID, 'seg-filled','dseg'])]
        else:
            list1 = [new_Ref_file,
                     '_'.join([ID, 'desc-norm', 'T1w']),
                     '_'.join([ID, 'seg-wm',    'dseg']),
                     '_'.join([ID, 'seg-4FS',   'dseg']),
                     '_'.join([ID, 'seg-filled','dseg'])]
        if change_hd == 1:
            preFS.resamp(list1, data_path,reference,BALSAname, scaling, new_size, resamp, diary_file)

        # set the data and folder for freesurfer
        preFS.toFS(list1, ID, change_hd,scaling[0],data_path,reference,BALSAname,diary_file,sing_fs)

    # ...................................................................................................................

    if 11 in Skip_step:
        run_cmd.msg('INFO: skip step ' + str(11), diary_file, 'OKGREEN')
    else:
        FS_surf.Wcreate(FS_dir, ID, 'lr', diary_file, sing_fs, export_fs)

        if check_visualy_each_img == True:
            FS_freeview.FSview(FS_dir, ID, 'white', sing_fs, FS_refs)

        # You can go grab a cup of coffe, it can take more than an hour...
        FS_surf.Wconfig(FS_dir, FS_refs, ID, 'lr', diary_file, sing_fs, export_fs)

    # ...................................................................................................................

    if 12 in Skip_step:
        run_cmd.msg('INFO: skip step ' + str(12), diary_file, 'OKGREEN')
    else:
        FS_surf.Pcreate(ID, 'lr', diary_file, export_fs)

        if opi(T2_img):
            f=open(new_Ref_file.replace('.nii.gz','.json'))
            jsonF = json.load(f)
            if not 'transformation of T2w image' in jsonF['Description']:
                if opi(opj(dir_prepro,ID + '_rescale.json')):
                    change_hd, _, _, _, _ = smallbrain.check(opj(dir_prepro,ID + '_rescale.json'))
                else :
                    change_hd = 0

                FS_surf.PcreateT2(T2_img,data_path,aseg_img,'lr',
                                  change_hd,diary_file,sing_fs,export_fs)

        #_12_make_pial.make_pial(FS_dir, ID, type_norm, otheranat, Hmin, Ref_file,do_surfacewith, overwrite,sing_fs, diary_file)
        if check_visualy_each_img == True:
            FS_freeview.FSview(FS_dir, ID, 'pial', sing_fs, FS_refs)

    #...................................................................................................................

    if 'flat_map' in Skip_step:
        run_cmd.msg('INFO: skip step ' + str('flat_map'), diary_file, 'OKGREEN')
    else:
        # Flat maps
        #### 4.4.1 Left hemisphere
        FS_flat.map(FS_dir, ID, 'left', cmd_tksurfer, cmd_flatten, diary_file, sing_fs, export_fs)
        #### 4.4.2 Right hemisphere
        FS_flat.map(FS_dir, ID, 'right', cmd_tksurfer, cmd_flatten, diary_file, sing_fs, export_fs)

    # ...................................................................................................................

    if 14 in Skip_step:
        run_cmd.msg('INFO: skip step ' + str(14), diary_file, 'OKGREEN')
    else:
        print(str(list_atlas))
        print(str(addatlas))
        list_atlas = templatefeat.atlas(list_atlas,addatlas)
        print(str(list_atlas))
        a = list()
        for x in list_atlas[0]:
            a.append(opj(labels_dir, ID + '_seg-' + x + '_dseg.nii.gz'))

        if opi(opj(dir_prepro, ID + '_rescale.json')):
            change_hd, _, _, _, _ = smallbrain.check(opj(dir_prepro, ID + '_rescale.json'))
        else:
            change_hd = 0

        FS_finalise.FS_finalise(FS_dir, ID, a, list_atlas[1], list_atlas[3], change_hd,
                                diary_file, sing_fs, export_fs,path_label_code, species, FS_refs,sing_afni,sing_wb)

    # ...................................................................................................................

    if 15 in Skip_step:
        run_cmd.msg('INFO: skip step ' + str(15), diary_file, 'OKGREEN')
    else:
        if opi(opj(dir_prepro, ID + '_rescale.json')):
            _, scaling, resamp, _, DATfile = smallbrain.check(opj(dir_prepro, ID + '_rescale.json'))
        else:
            scaling =''
            resamp  =''
            DATfile =''

        FS_2_WB.WB_prep(cmd_mris,FS_dir,FS_refs,data_path,ID,species,scaling,resamp,DATfile,
                       reference,list_atlas,balsa_folder,BALSAname,balsa_brainT1,
                       path_label_code,targetsuffix,listTimage,do_Nu,
                       diary_file,sing_fs,sing_wb,sing_afni,export_fs, list_transfo, type_norm,iso)

    # ...................................................................................................................

    if 16 in Skip_step:
        run_cmd.msg('INFO: skip step ' + str(16), diary_file, 'OKGREEN')
    else:

        info = backtonative.get(ID, data_path, bids_dir, Session, max_ses, targetsuffix, type_norm, BASE_SS, BASE_mask,
                                BASE_atlas_folder, study_template_atlas_folder,
                                creat_study_template, coregistration_longitudinal, reference, 'Final', species, diary_file)
        print(info)
        _16_anat_QC_SNR.anat_QC(info,type_norm, ID,volumes_dir, masks_dir, labels_dir, targetsuffix, listTimage,
                                sing_afni, diary_file)
    # ...................................................................................................................

    if 'itk_2' in Skip_step:
        run_cmd.msg('INFO: skip step ' + str('itk_2'), diary_file, 'OKGREEN')
    else:
        _200_Data_QC._itk_check_masks(volumes_dir, masks_dir, ID, type_norm,sing_itk,diary_file)

    # ...................................................................................................................

    if 'Clean' in Skip_step:
        run_cmd.msg('INFO: skip step ' + str('Clean'), diary_file, 'OKGREEN')
    else:
        _100_Data_Clean.clean(dir_prepro, volumes_dir, masks_dir, ID, diary_file)