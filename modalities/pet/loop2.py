#import
import os
import glob
import ants
import json

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

antsaff     = '0GenericAffine.mat'
antswarp    = '1Warp.nii.gz'
antsinvwarp = '1InverseWarp.nii.gz'

from Tools import run_cmd,getpath,extract_filename
from modalities.fMRI import chooseanat
from modalities.pet import studyT,MoCo,coreg2anat,preproPET

def run(ID,tracor,data_path, reference, anat_pet_same_space,IhaveanANAT,anat_raw,anat_subject,brainmask,G_mask, W_mask,V_mask,dir_transfo,templateT1w, templatemask,T_PET,do_manual, atlas_file,atlas_label,Need_check_MC,Need_check_coreg,nbvoid,studyTfile,interp_pet,Skip_step,REF, isotope, DV_time,Do_cut, duration,suv,bratio,ratio,strm,diary_file,sing_afni, sing_fs,sing_wb,sing_turku):

    Suffix = ['reorient',
              'MoCo-within',
              'MoCo-between']

    (dir_prepro_raw, dir_prepro_raw_process, dir_prepro_raw_masks, dir_prepro_raw_matrices, path_pet,
     dir_prepro_orig, dir_prepro_orig_labels, dir_prepro_orig_masks,
     dir_prepro_orig_process, dir_prepro_orig_postprocessed,
     dir_prepro_acpc, dir_prepro_acpc_labels, dir_prepro_acpc_masks,
     dir_prepro_acpc_process, dir_prepro_acpc_postprocessed,
     dir_prepro_template, dir_prepro_template_labels, dir_prepro_template_masks,
     dir_prepro_template_process, dir_prepro_template_postprocessed, dir_prepro_acpc_matrices, dir_prepro_orig_matrices)\
        = getpath.pet(data_path, reference)


    list_img = sorted(glob.glob(opj(dir_prepro_orig_process, '_'.join([ID,'trc-' + tracor, '*' + '-'.join(['desc',Suffix[0]]),'pet.nii.gz']))))
    print(list_img)

    filetype = []
    for img in list_img:
        if opi(img):
            imghd = ants.image_header_info(img)
            if imghd['nDimensions'] == 4:
                filetype.append(1)
                Tmatrices = [opj(dir_prepro_orig_matrices, '_'.join(['trc-' + tracor, 'to', 'studyTemplate', antsaff])),
                             opj(dir_prepro_orig_matrices, '_'.join(['trc-' + tracor, 'to', 'studyTemplate', antsinvwarp]))]

                refimg = studyT.automatic_crop(dir_prepro_orig_process,img,Tmatrices,studyTfile.replace('_pet.nii.gz', '_desc-crop_pet.nii.gz'),interp_pet)

                if 4 in Skip_step:
                    nl = 'skip step ' + str(1)
                    run_cmd.msg(nl, diary_file, 'OKGREEN')
                else:
                    MoCo.dynamic(img, refimg, nbvoid, Need_check_MC, diary_file, sing_afni, sing_fs)
            else:
                filetype.append(0)

    if 5 in Skip_step:
        nl = 'skip step ' + str(1)
        run_cmd.msg(nl, diary_file, 'OKGREEN')
    else:
        if list_img>1:
            for i in range(len(list_img)):
                if filetype[i]==1:
                    list_img[i] = list_img[i].replace(Suffix[0],Suffix[1])

        MoCo.static(list_img, diary_file, sing_afni, sing_fs)
        for i in range(len(list_img)):
            if filetype[i] == 0:
                list_img[i] = list_img[i].replace(Suffix[0], Suffix[2])
            elif filetype[i] == 1:
                list_img[i] = list_img[i].replace(Suffix[1], Suffix[2])


    if 6 in Skip_step:
        nl = 'skip step ' + str(1)
        run_cmd.msg(nl, diary_file, 'OKGREEN')
    else:
        for i in range(len(list_img)):
            coreg2anat.doit(list_img[i],
                            anat_raw,T_PET,brainmask,W_mask,V_mask,G_mask,
                            reference,atlas_file,atlas_label,dir_transfo,dir_prepro_orig_process, dir_prepro_orig_labels, dir_prepro_orig_masks,
            dir_prepro_orig_matrices, anat_pet_same_space, do_manual, IhaveanANAT,IhaveaCT, dir_prepro_orig_ct_matrices, dir_prepro_acpc,
            dir_prepro_acpc_process, dir_prepro_acpc_label, dir_prepro_acpc_masks, dir_prepro_template,
            dir_prepro_template_labels, dir_prepro_template_masks, BASE_template,
            dir_prepro_template_process, BASE_mask, BASE_Gmask, BASE_Wmask,
            Need_check_coreg, diary_file, sing_afni, sing_fs)


    if 7 in Skip_step:
        nl = 'skip step ' + str(1)
        run_cmd.msg(nl, diary_file, 'OKGREEN')
    else:
        dseg = sorted(glob.glob(opj(dir_prepro_orig_labels, '*_seg-' + atlas["name"] + '_dseg.nii.gz')))[0]

        vol_file     = opj(dir_prepro_orig_postprocessed, 'ROI_analysis', opb(dseg).replace('.nii.gz', '_vol_ROIs.txt'))

        for i in range(len(list_img)):
            Kinetic_file = opj(dir_prepro_orig_postprocessed, 'ROI_analysis', opb(list_img[i]).replace('.nii.gz', '_Kinetic_ROIs.txt'))

            if suv == 1 or bratio == 1:
                preproPET.suv(list_img[i], dseg, reference, brainmask, anat_subject, templateT1w, templatemask, '',
                              dir_transfo, dir_prepro_orig_labels, Do_cut, duration, ratio, bratio, diary_file, sing_afni,
                              sing_wb,
                              anat_pet_same_space)

            if strm == 1:
                preproPET.srtm(list_img[i], dseg,
                               vol_file,
                               Kinetic_file,
                               reference, brainmask, W_mask, V_mask, G_mask, anat_subject, templateT1w,
                               templatemask, '1', dir_transfo,
                               dir_prepro_orig_labels, REF, isotope, DV_time, diary_file, sing_afni, sing_turku, sing_wb,
                               anat_pet_same_space)
