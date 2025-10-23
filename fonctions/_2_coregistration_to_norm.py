#import
import os
import ants
import json

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists

from Tools import run_cmd
from fonctions.extract_filename import extract_filename
from fonctions import _2b_fix_orient, _2a_correct_img


def coregist_to_norm(correction_direction, dir_prepro_fmap, dir_prepro_orig_process, RS, RS_map, nb_run, recordings, REF_int, list_map,
                     deoblique, orientation, DwellT, n_for_ANTS, overwrite,sing_afni,sing_fsl,dmap,dbold,config_fd,diary_file):


    nl = '##  Working on step ' + str(2) + '(function: _2_coregistration_to_norm).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    if recordings != 'very_old':
        
        nl = 'INFO: DwellT is equal to ' + str(DwellT) + ' please check!!'
        run_cmd.msg(nl, diary_file, 'OKGREEN')

    if recordings == 'old':

        i       = 0
        r       = REF_int
        root    = extract_filename(RS_map[i])
        root_RS = extract_filename(RS[r])

        topup_f = open(opj(dir_prepro_orig_process, '4topup.txt'), "w")
        topup_f.write(dmap + ' \n')
        topup_f.write(dbold + ' \n')
        topup_f.close()
        topup_file = [opj(dir_prepro_orig_process, '4topup.txt'), config_fd]
        fMRI_runMean_n4Bias = opj(dir_prepro_orig_process, root + '_space-func_desc-runMean_n4Bias.nii.gz')
        fMRI_runMean_fieldmap_rads = opj(dir_prepro_fmap, root + '_space-func_desc-runMean_fieldmap_rads.nii.gz')
        fMRI_runMean_unwarpped = opj(dir_prepro_orig_process, root + '_space-func_desc-runMean_unwarped.nii.gz')

        ### 1.0 Start correct img
        _2a_correct_img.correct_img(dir_prepro_orig_process, RS, list_map, RS_map, i, r, recordings, overwrite, sing_afni, sing_fsl, topup_file, diary_file)

        command = (sing_fsl + 'fugue -i ' + fMRI_runMean_n4Bias +
                   ' --dwell=' + str(DwellT) + ' --loadfmap=' + fMRI_runMean_fieldmap_rads +
                   ' --unwarpdir=' + correction_direction + ' -u ' + fMRI_runMean_unwarpped)
        run_cmd.run(command, diary_file)

        dictionary = {"Sources": [fMRI_runMean_n4Bias,
                                  fMRI_runMean_fieldmap_rads],
                      "Description": 'Distortion correction (fugue from FSL).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_prepro_orig_process, root_RS + '_xdtrf_mean_preWARP.json'), "w") as outfile:
            outfile.write(json_object)


        ### 2.0 Start fix_orient

        imgO = '_xdtr_mean_deob_ref_fudge.nii.gz'
        imgI = '_xdtrf_mean_preWARP.nii.gz'

        _2b_fix_orient.fix_orient(imgO, imgI, dir_prepro_orig_process, root_RS, deoblique, orientation,
                                            overwrite, sing_afni, diary_file)

    elif recordings == '2_mapdir':
        i = 0
        r = REF_int
        root = extract_filename(RS_map[i])
        root_RS = extract_filename(RS[r])
        topup_f = open(opj(dir_prepro_orig_process, '4topup.txt'), "w")
        topup_f.write(dmap + ' \n')
        topup_f.write(dbold + ' \n')
        topup_f.close()
        topup_file = [opj(dir_prepro_orig_process, '4topup.txt'), config_fd]

        ### 1.0 Start correct img

        _2a_correct_img.correct_img(dir_prepro_orig_process, RS, list_map, RS_map, i, r, recordings, overwrite, sing_afni, sing_fsl, topup_file, diary_file)

        command = (sing_fsl + 'fugue -i ' + fMRI_runMean_n4Bias +
                   ' --dwell=' + str(DwellT) + ' --loadfmap=' + fMRI_runMean_fieldmap_rads +
                   ' --unwarpdir=' + correction_direction + ' -u ' + fMRI_runMean_unwarpped)
        run_cmd.run(command, diary_file)

        dictionary = {"Sources": [fMRI_runMean_n4Bias,
                                  fMRI_runMean_fieldmap_rads],
                      "Description": 'Distortion correction (fugue from FSL).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_prepro_orig_process, root_RS + '_xdtrf_mean_preWARP.json'), "w") as outfile:
            outfile.write(json_object)

        ### 2.0 Start fix_orient

        imgO = '_xdtr_mean_deob_ref_fudge.nii.gz'
        imgI = '_xdtrf_mean_preWARP.nii.gz'

        _2b_fix_orient.fix_orient(imgO, imgI, dir_prepro_orig_process, root_RS, deoblique, orientation, overwrite, sing_afni, diary_file)

    elif recordings == 'new':
        topup_f = open(opj(dir_prepro_orig_process, '4topup.txt'), "w")
        topup_f.write(dmap + ' \n')
        topup_f.write(dbold + ' \n')
        topup_f.close()
        topup_file = [opj(dir_prepro_orig_process, '4topup.txt'), config_fd]

        for i, r in zip(range(0, int(len(list_map))), range(0, int(nb_run))):
            root = extract_filename(RS_map[i])
            root_RS = extract_filename(RS[r])

            for f in os.listdir(opj(dir_prepro_orig_process,'tmp')):
                os.remove(os.path.join(opj(dir_prepro_orig_process,'tmp'), f))

            ### 1.0 Start correct img

            _2a_correct_img.correct_img(dir_prepro_orig_process, RS, list_map, RS_map, i, r, recordings, overwrite, sing_afni, sing_fsl, topup_file)

            command = (sing_fsl + 'fugue -i ' + fMRI_runMean_n4Bias +
                       ' --dwell=' + str(DwellT) + ' --loadfmap=' + fMRI_runMean_fieldmap_rads +
                       ' --unwarpdir=' + correction_direction + ' -u ' + fMRI_runMean_unwarpped)
            run_cmd.run(command, diary_file)

            dictionary = {"Sources": [fMRI_runMean_n4Bias,
                                      fMRI_runMean_fieldmap_rads],
                          "Description": 'Distortion correction (fugue from FSL).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_prepro_orig_process, root_RS + '_xdtrf_mean_preWARP.json'), "w") as outfile:
                outfile.write(json_object)
            ### 2.0 Start fix_orient

            imgO = '_xdtr_mean_deob_ref_fudge.nii.gz'
            imgI = '_xdtrf_mean_preWARP.nii.gz'

            _2b_fix_orient.fix_orient(imgO, imgI, dir_prepro_orig_process, root_RS, deoblique, orientation, overwrite, sing_afni, diary_file)

    elif recordings == 'very_old':
        # ### 1.0 Normalization between runs
        r = REF_int
        root_RS = extract_filename(RS[r])
        command = (sing_afni + '3dcopy ' + fMRI_runMean_n4Bias +
                   ' ' + opj(dir_prepro_orig_process, root_RS + '_xdtrf_mean_preWARP.nii.gz') + overwrite)
        run_cmd.run(command, diary_file)

        dictionary = {"Sources": fMRI_runMean_n4Bias,
                      "Description": 'Copy.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_prepro_orig_process, root_RS + '_xdtrf_mean_preWARP.json'), "w") as outfile:
            outfile.write(json_object)

        ### 2.0 Start fix_orient

        imgO = '_xdtr_mean_deob_ref_fudge.nii.gz'
        imgI = '_xdtrf_mean_preWARP.nii.gz'

        _2b_fix_orient.fix_orient(imgO, imgI, dir_prepro_orig_process, root_RS, deoblique, orientation, overwrite, sing_afni, diary_file)

    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### 
    ###                                              fix header problems                                                                ###
    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### 

    # 1.0 first you need to fix the header problems and potentially fit anat images

    for imgO, imgI in zip(['_xdtr_mean_deob.nii.gz', '_xdtr_deob.nii.gz'], ['_xdtr_mean.nii.gz', '_xdtr.nii.gz']):
            for r in range(0, int(nb_run)):
                root_RS = extract_filename(RS[r])

                ### 2.0 Start fix_orient
                _2b_fix_orient.fix_orient(imgO, imgI, dir_prepro_orig_process, root_RS, deoblique, orientation, overwrite, sing_afni, diary_file)


    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### 
    ###                                          co-registration of each run to the norm                                                ###
    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ###

    for r in range(0, int(nb_run)):
        root_RS     = extract_filename(RS[r])
        root_RS_ref = extract_filename(RS[REF_int])
        REF = ants.image_read(opj(dir_prepro_orig_process, root_RS_ref + '_xdtr_mean_deob_ref_fudge.nii.gz'))

        if not root_RS == root_RS_ref:  # do not process ref...

            # 1.0 calculate co-registration mean image ref to mean image func
            IMG = ants.image_read(opj(dir_prepro_orig_process, root_RS + '_xdtr_mean_deob.nii.gz'))
            mTx = ants.registration(fixed=REF, moving=IMG,
                                    type_of_transform='SyNCC',
                                    initial_transform=None,
                                    outprefix=opj(dir_prepro_orig_process, root_RS + '_xdtr_mean_warp_'))

            ##  Apply the transformation  to the mean image
            moved = ants.apply_transforms(fixed=REF, moving=IMG,
                                          transformlist=mTx['fwdtransforms'],
                                          interpolator=n_for_ANTS)
            ants.image_write(moved, opj(dir_prepro_orig_process, root_RS + '_xdtr_mean_warp.nii.gz'), ri=False)

            dictionary = {"Sources": [opj(dir_prepro_orig_process, root_RS + '_xdtr_mean_deob.nii.gz'),
                                      opj(dir_prepro_orig_process, root_RS_ref + '_xdtr_mean_deob_ref_fudge.nii.gz')],
                          "Description": 'Coregistration (ANTspy).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_prepro_orig_process, root_RS + '_xdtr_mean_warp.json'), "w") as outfile:
                outfile.write(json_object)

            # 2.0 apply to all the volume in the func (_xdtr_deob)
            FUNC = ants.image_read(opj(dir_prepro_orig_process, root_RS + '_xdtr_deob.nii.gz'))
            moved = ants.apply_transforms(fixed=REF, moving=FUNC,
                                          transformlist=mTx['fwdtransforms'],
                                          interpolator=n_for_ANTS, imagetype=3)
            ants.image_write(moved, opj(dir_prepro_orig_process,root_RS + '_xdtrf_2ref.nii.gz'), ri=False)

            dictionary = {"Sources": [opj(dir_prepro_orig_process, root_RS + '_xdtr_mean_deob.nii.gz'),
                                      opj(dir_prepro_orig_process, root_RS_ref + '_xdtrf_2ref.nii.gz')],
                          "Description": 'Coregistration (ANTspy).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_prepro_orig_process, root_RS + '_xdtrf_2ref.json'), "w") as outfile:
                outfile.write(json_object)

        else:
            command = (sing_afni + '3dcopy ' + opj(dir_prepro_orig_process, root_RS_ref + '_xdtr_deob.nii.gz') +
                       ' ' + opj(dir_prepro_orig_process,root_RS + '_xdtrf_2ref.nii.gz') + overwrite)
            run_cmd.run(command, diary_file)

            dictionary = {"Sources": opj(dir_prepro_orig_process, root_RS_ref + '_xdtr_deob.nii.gz'),
                          "Description": 'Copy.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_prepro_orig_process, root_RS + '_xdtrf_2ref.json'), "w") as outfile:
                outfile.write(json_object)