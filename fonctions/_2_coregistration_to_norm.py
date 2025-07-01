#import
import os
import subprocess
from fonctions.extract_filename import extract_filename
import ants
import datetime
import json
import fonctions
import fonctions._2b_fix_orient
import fonctions._2a_correct_img
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

def coregist_to_norm(correction_direction, dir_fMRI_Refth_RS_prepro1, RS, RS_map, nb_run, recordings, REF_int, list_map,
                     deoblique, orientation, DwellT, n_for_ANTS, overwrite,s_bind,afni_sif,fsl_sif,dmap,dbold,config_fd,diary_file):


    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  Working on step ' + str(2) + '(function: _2_coregistration_to_norm).  ##'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')
    diary.write(f'\n')
    diary.close()

    if recordings != 'very_old':

        diary = open(diary_file, "a")
        nl = 'INFO: DwellT is equal to ' + str(DwellT) + ' please check!!'
        print(bcolors.OKGREEN + nl + bcolors.ENDC)
        diary.write(f'\n{nl}')
        diary.write(f'\n')
        diary.close()

    if recordings == 'old':

        i       = 0
        r       = REF_int
        root    = extract_filename(RS_map[i])
        root_RS = extract_filename(RS[r])

        topup_f = open(opj(dir_fMRI_Refth_RS_prepro1, '4topup.txt'), "w")
        topup_f.write(dmap + ' \n')
        topup_f.write(dbold + ' \n')
        topup_f.close()
        topup_file = [opj(dir_fMRI_Refth_RS_prepro1, '4topup.txt'), config_fd]

        ### 1.0 Start correct img

        fonctions._2a_correct_img.correct_img(dir_fMRI_Refth_RS_prepro1, RS, list_map, RS_map, i, r, recordings, overwrite, s_bind, afni_sif, fsl_sif, topup_file, diary_file)

        command = 'singularity run' + s_bind + fsl_sif + 'fugue -i ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean' + '.nii.gz') + \
        ' --dwell=' + str(DwellT) + ' --loadfmap=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_rads' + '.nii.gz') + \
        ' --unwarpdir=' + correction_direction + ' -u ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP' + '.nii.gz')
        nl = spgo(command)
        print(nl)
        diary = open(diary_file, "a")
        diary.write(f'\n{nl}')
        diary.write(f'\n')
        diary.close()

        dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean' + '.nii.gz'),
                                  opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_rads' + '.nii.gz')],
                      "Description": 'Distortion correction (fugue from FSL).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP.json'), "w") as outfile:
            outfile.write(json_object)


        ### 2.0 Start fix_orient

        imgO = '_xdtr_mean_deob_ref_fudge.nii.gz'
        imgI = '_xdtrf_mean_preWARP.nii.gz'

        fonctions._2b_fix_orient.fix_orient(imgO, imgI, dir_fMRI_Refth_RS_prepro1, root_RS, deoblique, orientation,
                                            overwrite, s_bind, afni_sif, diary_file)

    elif recordings == '2_mapdir':
        i = 0
        r = REF_int
        root = extract_filename(RS_map[i])
        root_RS = extract_filename(RS[r])
        topup_f = open(opj(dir_fMRI_Refth_RS_prepro1, '4topup.txt'), "w")
        topup_f.write(dmap + ' \n')
        topup_f.write(dbold + ' \n')
        topup_f.close()
        topup_file = [opj(dir_fMRI_Refth_RS_prepro1, '4topup.txt'), config_fd]

        ### 1.0 Start correct img

        fonctions._2a_correct_img.correct_img(dir_fMRI_Refth_RS_prepro1, RS, list_map, RS_map, i, r, recordings, overwrite, s_bind, afni_sif, fsl_sif, topup_file, diary_file)

        command = 'singularity run' + s_bind + fsl_sif + 'fugue -i ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean' + '.nii.gz') + \
        ' --dwell=' + str(DwellT) + ' --loadfmap=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_rads' + '.nii.gz') + \
        ' --unwarpdir=' + correction_direction + ' -u ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP' + '.nii.gz')
        nl = spgo(command)
        print(nl)
        diary = open(diary_file, "a")
        diary.write(f'\n{nl}')
        diary.write(f'\n')
        diary.close()

        dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean' + '.nii.gz'),
                                  opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_rads' + '.nii.gz')],
                      "Description": 'Distortion correction (fugue from FSL).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP.json'), "w") as outfile:
            outfile.write(json_object)

        ### 2.0 Start fix_orient

        imgO = '_xdtr_mean_deob_ref_fudge.nii.gz'
        imgI = '_xdtrf_mean_preWARP.nii.gz'

        fonctions._2b_fix_orient.fix_orient(imgO, imgI, dir_fMRI_Refth_RS_prepro1, root_RS, deoblique, orientation, overwrite, s_bind, afni_sif, diary_file)

    elif recordings == 'new':
        topup_f = open(opj(dir_fMRI_Refth_RS_prepro1, '4topup.txt'), "w")
        topup_f.write(dmap + ' \n')
        topup_f.write(dbold + ' \n')
        topup_f.close()
        topup_file = [opj(dir_fMRI_Refth_RS_prepro1, '4topup.txt'), config_fd]

        for i, r in zip(range(0, int(len(list_map))), range(0, int(nb_run))):
            root = extract_filename(RS_map[i])
            root_RS = extract_filename(RS[r])

            for f in os.listdir(opj(dir_fMRI_Refth_RS_prepro1,'tmp')):
                os.remove(os.path.join(opj(dir_fMRI_Refth_RS_prepro1,'tmp'), f))

            ### 1.0 Start correct img

            fonctions._2a_correct_img.correct_img(dir_fMRI_Refth_RS_prepro1, RS, list_map, RS_map, i, r, recordings, overwrite, s_bind, afni_sif, fsl_sif, topup_file)

            command = 'singularity run' + s_bind + fsl_sif + 'fugue -i ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean' + '.nii.gz') + \
            ' --dwell=' + str(DwellT) + ' --loadfmap=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_rads' + '.nii.gz') + \
            ' --unwarpdir=' + correction_direction + ' -u ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP' + '.nii.gz')
            nl = spgo(command)
            print(nl)
            diary = open(diary_file, "a")
            diary.write(f'\n{nl}')
            diary.write(f'\n')
            diary.close()
            dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean' + '.nii.gz'),
                                      opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_rads' + '.nii.gz')],
                          "Description": 'Distortion correction (fugue from FSL).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP.json'), "w") as outfile:
                outfile.write(json_object)
            ### 2.0 Start fix_orient

            imgO = '_xdtr_mean_deob_ref_fudge.nii.gz'
            imgI = '_xdtrf_mean_preWARP.nii.gz'

            fonctions._2b_fix_orient.fix_orient(imgO, imgI, dir_fMRI_Refth_RS_prepro1, root_RS, deoblique, orientation, overwrite, s_bind, afni_sif, diary_file)

    elif recordings == 'very_old':
        # ### 1.0 Normalization between runs
        r = REF_int
        root_RS = extract_filename(RS[r])
        command = 'singularity run' + s_bind + afni_sif + '3dcopy ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean' + '.nii.gz') + \
                  ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP.nii.gz') + overwrite
        nl = spgo(command)
        print(nl)
        diary = open(diary_file, "a")
        diary.write(f'\n{nl}')
        diary.write(f'\n')
        diary.close()

        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean' + '.nii.gz'),
                      "Description": 'Copy.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP.json'), "w") as outfile:
            outfile.write(json_object)

        ### 2.0 Start fix_orient

        imgO = '_xdtr_mean_deob_ref_fudge.nii.gz'
        imgI = '_xdtrf_mean_preWARP.nii.gz'

        fonctions._2b_fix_orient.fix_orient(imgO, imgI, dir_fMRI_Refth_RS_prepro1, root_RS, deoblique, orientation, overwrite, s_bind, afni_sif, diary_file)

    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### 
    ###                                              fix header problems                                                                ###
    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### 

    # 1.0 first you need to fix the header problems and potentially fit anat images

    for imgO, imgI in zip(['_xdtr_mean_deob.nii.gz', '_xdtr_deob.nii.gz'], ['_xdtr_mean.nii.gz', '_xdtr.nii.gz']):
            for r in range(0, int(nb_run)):
                root_RS = extract_filename(RS[r])

                ### 2.0 Start fix_orient
                fonctions._2b_fix_orient.fix_orient(imgO, imgI, dir_fMRI_Refth_RS_prepro1, root_RS, deoblique, orientation, overwrite, s_bind, afni_sif, diary_file)


    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### 
    ###                                          co-registration of each run to the norm                                                ###
    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ###

    for r in range(0, int(nb_run)):
        root_RS = extract_filename(RS[r])
        root_RS_ref = extract_filename(RS[REF_int])
        REF = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1, root_RS_ref + '_xdtr_mean_deob_ref_fudge.nii.gz'))

        if not root_RS == root_RS_ref:  # do not process ref...

            # 1.0 calculate co-registration mean image ref to mean image func
            IMG = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_deob.nii.gz'))
            mTx = ants.registration(fixed=REF, moving=IMG,
                                    type_of_transform='SyNCC',
                                    initial_transform=None,
                                    outprefix=opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp_'))

            ##  Apply the transformation  to the mean image
            moved = ants.apply_transforms(fixed=REF, moving=IMG,
                                          transformlist=mTx['fwdtransforms'],
                                          interpolator=n_for_ANTS)
            ants.image_write(moved, opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp.nii.gz'), ri=False)

            dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_deob.nii.gz'),
                                      opj(dir_fMRI_Refth_RS_prepro1, root_RS_ref + '_xdtr_mean_deob_ref_fudge.nii.gz')],
                          "Description": 'Coregistration (ANTspy).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp.json'), "w") as outfile:
                outfile.write(json_object)

            # 2.0 apply to all the volume in the func (_xdtr_deob)
            FUNC = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_deob.nii.gz'))
            moved = ants.apply_transforms(fixed=REF, moving=FUNC,
                                          transformlist=mTx['fwdtransforms'],
                                          interpolator=n_for_ANTS, imagetype=3)
            ants.image_write(moved, opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref.nii.gz'), ri=False)

            dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_deob.nii.gz'),
                                      opj(dir_fMRI_Refth_RS_prepro1, root_RS_ref + '_xdtrf_2ref.nii.gz')],
                          "Description": 'Coregistration (ANTspy).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref.json'), "w") as outfile:
                outfile.write(json_object)

        else:
            command = 'singularity run' + s_bind + afni_sif + '3dcopy ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS_ref + '_xdtr_deob.nii.gz') + \
                      ' ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref.nii.gz') + overwrite
            nl = spgo(command)
            print(nl)
            diary = open(diary_file, "a")
            diary.write(f'\n{nl}')
            diary.write(f'\n')
            diary.close()

            dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root_RS_ref + '_xdtr_deob.nii.gz'),
                          "Description": 'Copy.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref.json'), "w") as outfile:
                outfile.write(json_object)