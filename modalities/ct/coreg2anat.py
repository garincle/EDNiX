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
from modalities.anat.freesurfer import FS_freeview

T1interp    = 'hammingWindowedSinc'

antsaff     = '0GenericAffine.mat'
antswarp    = '1Warp.nii.gz'
antsinvwarp = '1InverseWarp.nii.gz'

resamp_interp = 'Cu'

tmp_trans = 'translation'
transCT   = 'from-ct-to-T1w'
shiftname = 'shift'
acpcname  = 'acpc'




def apply(IMG,Ref_native_raw,ctimg,Ref,Ref_native_template,ct2anat,anat2ct, ct2acpc,fname,dir_prepro_orig_matrices,to_acpc,T_CT):
    # CT to anat
    moved1 = ants.apply_transforms(fixed=Ref, moving=ctimg,
                                   transformlist=opj(dir_prepro_orig_matrices, '_'.join([fname, transCT, antsaff])),
                                   interpolator=T1interp)
    ants.image_write(moved1, ct2anat)

    dictionary = {"Sources": [IMG,
                              Ref_native_raw],
                  "Description": 'Rigid transformation of ct image to ' + T_CT + ' image (3dresample (AFNI)).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(anat2ct.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)


    # anat to ct
    moved2 = ants.apply_transforms(fixed=ctimg, moving=Ref,
                                   transformlist=opj(dir_prepro_orig_matrices, '_'.join([fname, transCT, antsaff])),
                                   whichtoinvert=[True],
                                   interpolator=T1interp)
    ants.image_write(moved2, anat2ct)
    dictionary = {"Sources": [IMG,
                              Ref_native_raw],
                  "Description": 'Rigid transformation of ' + T_CT + '  image to ct image (3dresample (AFNI)).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(anat2ct.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)

    # cT to anat in acpc space
    refacpc = ants.image_read(Ref_native_template)
    mv = ants.apply_transforms(fixed=refacpc, moving=ctimg,
                               transformlist=[to_acpc,
                                              opj(dir_prepro_orig_matrices, '_'.join([fname, transCT, antsaff]))],
                               interpolator=T1interp)
    ants.image_write(mv, ct2acpc)

    dictionary = {"Sources": [ct2anat,
                              Ref_native_template],
                  "Description": 'Rigid transformation of ct image to ' + T_CT + ' image (3dresample (AFNI)).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(anat2ct.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)


def doit(IMG,Ref_native_raw,T_CT,Ref_native_template,
         dir_transfo,dir_prepro_orig_process,dir_prepro_acpc,
         dir_prepro_acpc_process,
         dir_prepro_orig_matrices,same_space,do_manual,
         Need_check_coreg,diary_name,sing_afni,sing_fs):


    to_acpc = opj(dir_transfo,'_'.join([acpcname,antsaff]))

    fname,_ = extract_filename.get_desc(opb(IMG),'desc')

    # ct image in anat space (only for visualization)
    ct2anat = opj(dir_prepro_orig_process, '_'.join([fname,'space-anat','ct.nii.gz']))
    ct2anat_tmp = opj(dir_prepro_acpc_process, '_'.join([fname, 'space-acpc','desc-before-manual-correct', 'ct.nii.gz']))
    ct2acpc = opj(dir_prepro_acpc, '_'.join([fname, 'space-acpc', 'ct.nii.gz']))
    ct2acpc_tmp = opj(dir_prepro_acpc_process, '_'.join([fname, 'space-acpc', 'ct.nii.gz']))
    # anat image in pet space (only for visualisation)
    anat2ct  = opj(dir_prepro_acpc, '_'.join([fname, 'space-orig',T_CT + '.nii.gz']))


    list_file = [IMG,Ref_native_raw,Ref_native_template,
                    opj(dir_transfo,'_'.join([acpcname,antsaff]))]

    for fneeded in list_file:
        if not opi(fneeded):
            raise FileNotFoundError(f"Input file not found: {fneeded}")

    # ..........................................................................................................

    nl = 'Co-registration with the MRI anatomy'
    run_cmd.msg(nl, diary_name,'OkGREEN')

    # from ct to T1w
    nl = 'Transform the ct image into the MRI anatomy'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    if same_space:
        # Resample to other modality resolution
        cmd=(sing_afni + '3dresample -overwrite -master ' + IMG + ' -rmode ' + resamp_interp +
             ' -prefix ' + anat2ct + ' -input ' +  Ref_native_raw)
        run_cmd.run(cmd,diary_name)

        dictionary = {"Sources": [IMG,
                                  Ref_native_raw],
                      "Description": 'Resampling ct image to ' + T_CT + ' image (3dresample (AFNI)).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(anat2ct.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)

        cmd = (sing_afni + '3dresample -overwrite -master ' + Ref_native_raw + ' -rmode ' + resamp_interp +
               ' -prefix ' + ct2anat + ' -input ' + IMG)
        run_cmd.run(cmd, diary_name)

        dictionary = {"Sources": [IMG,
                                  Ref_native_raw],
                      "Description": 'Resampling ' + T_CT + ' image to ct image (3dresample (AFNI)).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(ct2anat.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)

    else:
        ctimg = ants.image_read(IMG)
        Ref   = ants.image_read(Ref_native_raw)

        ants.registration(fixed=Ref, moving=ctimg,
                          type_of_transform='Translation',
                          outprefix=opj(dir_prepro_orig_matrices, tmp_trans + '_'))
        ants.registration(fixed=Ref, moving=ctimg, type_of_transform='Rigid',
                          initial_transform=opj(dir_prepro_orig_matrices, '_'.join([tmp_trans, antsaff])),
                          outprefix=opj(dir_prepro_orig_matrices,'_'.join([fname,transCT,''])))

        apply(IMG,Ref_native_raw,ctimg,Ref,Ref_native_template,ct2anat,anat2ct, ct2acpc,fname,dir_prepro_orig_matrices,to_acpc,T_CT)

        if do_manual:
            if do_manual == 'orig':
                shutil.copy(ct2anat, ct2anat_tmp)

                FS_freeview.check_coreg(ct2anat, Ref_native_template, '','','ct', sing_fs)

                newref = ants.image_read(ct2anat)
                # overwrite the previous matrices
                ants.registration(fixed=newref, moving=ctimg,
                                  type_of_transform='Translation',
                                  outprefix=opj(dir_prepro_orig_matrices, tmp_trans + '_'))

                ants.registration(fixed=newref, moving=ctimg, type_of_transform='Rigid',
                                  initial_transform=opj(dir_prepro_orig_matrices, '_'.join([tmp_trans, antsaff])),
                                  outprefix=opj(dir_prepro_orig_matrices, '_'.join([fname, transCT, ''])))


                apply(IMG,Ref_native_raw,ctimg,Ref,Ref_native_template,ct2anat,anat2ct, ct2acpc,fname,dir_prepro_orig_matrices,to_acpc,T_CT)

            elif do_manual == 'acpc':

                shutil.copy(ct2acpc,ct2acpc_tmp)
                # do the manual correction, save it as it is.
                FS_freeview.check_coreg(ct2acpc,Ref_native_template,'acpc','ct',sing_fs)

                dummy  = ants.image_read(ct2acpc)
                ct_tmp = ants.image_read(ct2acpc_tmp)

                # go backward
                mv = ants.apply_transforms(fixed=ct_tmp, moving=dummy,
                                           transformlist=[to_acpc],
                                           interpolator=T1interp,whichtoinvert= [True])
                # overwrite the previous matrices
                ants.registration(fixed=mv, moving=ctimg,
                                  type_of_transform='Translation',
                                  outprefix=opj(dir_prepro_orig_matrices, tmp_trans + '_'))

                ants.registration(fixed=mv, moving=ctimg, type_of_transform='Rigid',
                                  initial_transform=opj(dir_prepro_orig_matrices, '_'.join([tmp_trans, antsaff])),
                                  outprefix=opj(dir_prepro_orig_matrices, '_'.join([fname, transCT, ''])))

                apply(IMG,Ref_native_raw,ctimg, Ref, Ref_native_template, ct2anat, anat2ct, ct2acpc, fname, dir_prepro_orig_matrices, to_acpc,T_CT)

    if Need_check_coreg:

        FS_freeview.check_coreg(ct2anat, Ref_native_template, '', '', 'ct', sing_fs)


