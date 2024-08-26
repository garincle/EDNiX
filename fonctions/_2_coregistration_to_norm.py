
#import
import os
import subprocess
import glob
from fonctions.extract_filename import extract_filename
import ants

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput


def coregist_to_norm(anat_func_same_space, dir_prepro, correction_direction, dir_fMRI_Refth_RS_prepro1, RS, RS_map, nb_run, recordings, REF_int, list_map, study_fMRI_Refth, IgotbothT1T2, path_anat, otheranat,
    ID, Session, deoblique_exeption1, deoblique_exeption2, deoblique, orientation, DwellT, n_for_ANTS, overwrite,s_bind,afni_sif,fsl_sif,dmap,dbold,config_f):

    import fonctions._2a_correct_img
    import fonctions._2a_correct_img_newP
    import fonctions._2b_fix_orient

    if recordings != 'very_old':
        print('INFO: DwellT is equal to ' + str(DwellT) + ' please check!!' )

    if recordings == 'old':
        i = 0
        r = REF_int
        root = extract_filename(RS_map[i])
        root_RS = extract_filename(RS[r])

        topup_f = open(opj(dir_fMRI_Refth_RS_prepro1, '4topup.txt'), "w")
        topup_f.write(dmap + ' \n')
        topup_f.write(dbold + ' \n')
        topup_f.close()
        topup_file = [opj(dir_fMRI_Refth_RS_prepro1, '4topup.txt'), config_f]

        ###start correct img
        fonctions._2a_correct_img.correct_img(dir_fMRI_Refth_RS_prepro1, RS, list_map, RS_map, study_fMRI_Refth, i, r, overwrite,s_bind,afni_sif,fsl_sif,topup_file)

        command = 'singularity run' + s_bind + fsl_sif + 'fugue -i ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean' + '.nii.gz') + \
        ' --dwell=' + str(DwellT) + ' --loadfmap=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_rads' + '.nii.gz') + \
        ' --unwarpdir=' + correction_direction + ' -u ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP' + '.nii.gz')
        spco([command], shell=True)

        imgO = '_xdtr_mean_deob_ref_fudge.nii.gz'
        imgI = '_xdtrf_mean_preWARP.nii.gz'
        ###start fix_orient
        fonctions._2b_fix_orient.fix_orient(anat_func_same_space, imgO, imgI, dir_prepro, IgotbothT1T2, path_anat,
                                            otheranat, dir_fMRI_Refth_RS_prepro1, root_RS, ID, Session,
                                            deoblique_exeption1, deoblique_exeption2, deoblique, orientation, overwrite,s_bind,afni_sif)

    elif recordings == '2_mapdir':
        i = 0
        r = REF_int
        root = extract_filename(RS_map[i])
        root_RS = extract_filename(RS[r])
        topup_f = open(opj(dir_fMRI_Refth_RS_prepro1, '4topup.txt'), "w")
        topup_f.write(dmap + ' \n')
        topup_f.write(dbold + ' \n')
        topup_f.close()
        topup_file = [opj(dir_fMRI_Refth_RS_prepro1, '4topup.txt'), config_f]


        ###start correct img
        fonctions._2a_correct_img_newP.correct_img(dir_fMRI_Refth_RS_prepro1, RS, list_map, RS_map, study_fMRI_Refth, r, overwrite,s_bind,afni_sif,fsl_sif,topup_file)

        command = 'singularity run' + s_bind + fsl_sif + 'fugue -i ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean' + '.nii.gz') + \
        ' --dwell=' + str(DwellT) + ' --loadfmap=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_rads' + '.nii.gz') + \
        ' --unwarpdir=' + correction_direction + ' -u ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP' + '.nii.gz')
        spco([command], shell=True)

        imgO = '_xdtr_mean_deob_ref_fudge.nii.gz'
        imgI = '_xdtrf_mean_preWARP.nii.gz'
        ###start fix_orient
        fonctions._2b_fix_orient.fix_orient(anat_func_same_space, imgO, imgI, dir_prepro, IgotbothT1T2, path_anat,
                                            otheranat, dir_fMRI_Refth_RS_prepro1, root_RS, ID, Session,
                                            deoblique_exeption1, deoblique_exeption2, deoblique, orientation, overwrite,s_bind,afni_sif)

    elif recordings == 'new':

        topup_f = open(opj(dir_fMRI_Refth_RS_prepro1, '4topup.txt'), "w")
        topup_f.write(dmap + ' \n')
        topup_f.write(dbold + ' \n')
        topup_f.close()
        topup_file = [opj(dir_fMRI_Refth_RS_prepro1, '4topup.txt'), config_f]

        for i, r in zip(range(0, int(len(list_map))), range(0, int(nb_run))):
            root = extract_filename(RS_map[i])
            root_RS = extract_filename(RS[r])

            for f in os.listdir(opj(dir_fMRI_Refth_RS_prepro1,'tmp')):
                os.remove(os.path.join(opj(dir_fMRI_Refth_RS_prepro1,'tmp'), f))

            ###start correct img
            fonctions._2a_correct_img.correct_img(dir_fMRI_Refth_RS_prepro1, RS, list_map, RS_map, study_fMRI_Refth, i, r, overwrite,s_bind,afni_sif,fsl_sif,topup_file)

            command = 'singularity run' + s_bind + fsl_sif + 'fugue -i ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean' + '.nii.gz') + \
            ' --dwell=' + str(DwellT) + ' --loadfmap=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_rads' + '.nii.gz') + \
            ' --unwarpdir=' + correction_direction + ' -u ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP' + '.nii.gz')
            spco([command], shell=True)

            imgO = '_xdtr_mean_deob_ref_fudge.nii.gz'
            imgI = '_xdtrf_mean_preWARP.nii.gz'
            ###start fix_orient
            fonctions._2b_fix_orient.fix_orient(anat_func_same_space, imgO, imgI, dir_prepro, IgotbothT1T2, path_anat, otheranat, dir_fMRI_Refth_RS_prepro1, root_RS, ID, Session, deoblique_exeption1, deoblique_exeption2, deoblique, orientation, overwrite,
                                                s_bind,afni_sif)

    elif recordings == 'very_old':
        # norm between runs
        r = REF_int
        root_RS = extract_filename(RS[r])
        command = 'singularity run' + s_bind + afni_sif + '3dcopy ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean' + '.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP.nii.gz') + overwrite
        spco([command], shell=True)
        imgO = '_xdtr_mean_deob_ref_fudge.nii.gz'
        imgI = '_xdtrf_mean_preWARP.nii.gz'
        ###start fix_orient
        fonctions._2b_fix_orient.fix_orient(anat_func_same_space, imgO, imgI, dir_prepro, IgotbothT1T2, path_anat,
                                            otheranat, dir_fMRI_Refth_RS_prepro1, root_RS, ID, Session,
                                            deoblique_exeption1, deoblique_exeption2, deoblique, orientation, overwrite,s_bind,afni_sif)

    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### 
    ### ### #### ###### ### #### ###### ### #### ### fix header problems ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ##
    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### 

    #first you need to fix the header problems and potentially fit anat images
    for imgO, imgI in zip(['_xdtr_mean_deob.nii.gz', '_xdtr_deob.nii.gz'], ['_xdtr_mean.nii.gz', '_xdtr.nii.gz']):
            for r in range(0, int(nb_run)):
                root_RS = extract_filename(RS[r])
                ###start fix_orient
                fonctions._2b_fix_orient.fix_orient(anat_func_same_space, imgO, imgI, dir_prepro, IgotbothT1T2, path_anat, otheranat, dir_fMRI_Refth_RS_prepro1, root_RS, ID, Session, deoblique_exeption1, deoblique_exeption2, deoblique, orientation, overwrite,s_bind,afni_sif)

    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### 
    ### ### #### ### coregistration of each run to the norm ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ##
    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ###


    for r in range(0, int(nb_run)):
        root_RS = extract_filename(RS[r])
        root_RS_ref = extract_filename(RS[REF_int])
        REF = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1, root_RS_ref + '_xdtr_mean_deob_ref_fudge.nii.gz'))


        if not recordings == 'very_old':  # do not process ref...

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

            # 2.0 apply to all the volume in the func (_xdtr_deob)
            FUNC = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_deob.nii.gz'))
            moved = ants.apply_transforms(fixed=REF, moving=FUNC,
                                          transformlist=mTx1['fwdtransforms'],
                                          interpolator=n_for_ANTS, imagetype=3)
            ants.image_write(moved, opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref.nii.gz'), ri=False)


        else:
            command = 'singularity run' + s_bind + afni_sif + '3dcopy ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS_ref + '_xdtr_deob.nii.gz') + \
                      ' ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref.nii.gz') + overwrite
            spco([command], shell=True)
