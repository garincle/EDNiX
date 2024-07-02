    
### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### 
### ### #### ### creation and coregistration to the to the norm (opj(dir_fMRI_Refth_RS_prepro1,RS[REF_int].replace('.nii.gz','_xdtrf_mean.nii.gz'))  ???)
### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ###
    

### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### 
### ### #### ### creation the to the norm (opj(dir_fMRI_Refth_RS_prepro1,RS[REF_int].replace('.nii.gz','_xdtrf_mean.nii.gz')) XXX ???)#
### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### 
### ### #### ### ### ### #### ### ### ### #### ### ###  with ### ### #### ### ### ### #### ### ### ### #### ### ### ### ### ### #### ### 
#################################################### Method "new" and "old"" ###########################################################
#import
import os
import subprocess
import glob
import shutil
import sys


#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

from fonctions.extract_filename import extract_filename

def coregist_to_norm(TR, anat_func_same_space, dir_prepro, correction_direction, dir_fMRI_Refth_RS_prepro1, RS, RS_map, nb_run, recordings, REF_int, list_map, study_fMRI_Refth, IgotbothT1T2, path_anat, otheranat,
    ID, Session, deoblique_exeption1, deoblique_exeption2, deoblique, orientation, DwellT, n_for_ANTS, overwrite):

    import fonctions._2a_correct_img
    import fonctions._2a_correct_img_newP
    import fonctions._2b_fix_orient

    if recordings == 'old':
        i = 0
        r = REF_int
        root = extract_filename(RS_map[i])
        root_RS = extract_filename(RS[r])
        ###start correct img
        fonctions._2a_correct_img.correct_img(dir_fMRI_Refth_RS_prepro1, RS, list_map, RS_map, study_fMRI_Refth, i, r, overwrite)

        command = 'fugue -i ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean' + '.nii.gz') + \
        ' --dwell=' + str(DwellT) + ' --loadfmap=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_rads' + '.nii.gz') + \
        ' --unwarpdir=' + correction_direction + ' -u ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP' + '.nii.gz')
        spco([command], shell=True)

        imgO = '_xdtr_mean_deob_ref_fudge.nii.gz'
        imgI = '_xdtrf_mean_preWARP.nii.gz'
        ###start fix_orient
        fonctions._2b_fix_orient.fix_orient(anat_func_same_space, imgO, imgI, dir_prepro, IgotbothT1T2, path_anat,
                                            otheranat, dir_fMRI_Refth_RS_prepro1, root_RS, ID, Session,
                                            deoblique_exeption1, deoblique_exeption2, deoblique, orientation, overwrite)

    elif recordings == '2_mapdir':
        i = 0
        r = REF_int
        root = extract_filename(RS_map[i])
        root_RS = extract_filename(RS[r])
        ###start correct img
        fonctions._2a_correct_img_newP.correct_img(dir_fMRI_Refth_RS_prepro1, RS, list_map, RS_map, study_fMRI_Refth, r, overwrite)

        command = 'fugue -i ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean' + '.nii.gz') + \
        ' --dwell=' + str(DwellT) + ' --loadfmap=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_rads' + '.nii.gz') + \
        ' --unwarpdir=' + correction_direction + ' -u ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP' + '.nii.gz')
        spco([command], shell=True)

        imgO = '_xdtr_mean_deob_ref_fudge.nii.gz'
        imgI = '_xdtrf_mean_preWARP.nii.gz'
        ###start fix_orient
        fonctions._2b_fix_orient.fix_orient(anat_func_same_space, imgO, imgI, dir_prepro, IgotbothT1T2, path_anat,
                                            otheranat, dir_fMRI_Refth_RS_prepro1, root_RS, ID, Session,
                                            deoblique_exeption1, deoblique_exeption2, deoblique, orientation, overwrite)

    elif recordings == 'new':
        for i, r in zip(range(0, int(len(list_map))), range(0, int(nb_run))):
            root = extract_filename(RS_map[i])
            root_RS = extract_filename(RS[r])

            for f in os.listdir(opj(dir_fMRI_Refth_RS_prepro1,'tmp')):
                os.remove(os.path.join(opj(dir_fMRI_Refth_RS_prepro1,'tmp'), f))

            ###start correct img
            fonctions._2a_correct_img.correct_img(dir_fMRI_Refth_RS_prepro1, RS, list_map, RS_map, study_fMRI_Refth, i, r, overwrite)

            command = 'fugue -i ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean' + '.nii.gz') + \
            ' --dwell=' + str(DwellT) + ' --loadfmap=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_rads' + '.nii.gz') + \
            ' --unwarpdir=' + correction_direction + ' -u ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP' + '.nii.gz')
            spco([command], shell=True)

            imgO = '_xdtr_mean_deob_ref_fudge.nii.gz'
            imgI = '_xdtrf_mean_preWARP.nii.gz'
            ###start fix_orient
            fonctions._2b_fix_orient.fix_orient(anat_func_same_space, imgO, imgI, dir_prepro, IgotbothT1T2, path_anat, otheranat, dir_fMRI_Refth_RS_prepro1, root_RS, ID, Session, deoblique_exeption1, deoblique_exeption2, deoblique, orientation, overwrite)

    elif recordings == 'very_old':
        # norm between runs
        r = REF_int
        root_RS = extract_filename(RS[r])
        command = '3dcopy ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean' + '.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP.nii.gz') + overwrite
        spco([command], shell=True)
        imgO = '_xdtr_mean_deob_ref_fudge.nii.gz'
        imgI = '_xdtrf_mean_preWARP.nii.gz'
        ###start fix_orient
        fonctions._2b_fix_orient.fix_orient(anat_func_same_space, imgO, imgI, dir_prepro, IgotbothT1T2, path_anat,
                                            otheranat, dir_fMRI_Refth_RS_prepro1, root_RS, ID, Session,
                                            deoblique_exeption1, deoblique_exeption2, deoblique, orientation, overwrite)

    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### 
    ### ### #### ###### ### #### ###### ### #### ### fix header problems ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ##
    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### 

    #first you need to fix the header problems and potentially fit anat images
    for imgO, imgI in zip(['_xdtr_mean_deob.nii.gz', '_xdtr_deob.nii.gz'], ['_xdtr_mean.nii.gz', '_xdtr.nii.gz']):
            for r in range(0, int(nb_run)):
                root_RS = extract_filename(RS[r])
                ###start fix_orient
                fonctions._2b_fix_orient.fix_orient(anat_func_same_space, imgO, imgI, dir_prepro, IgotbothT1T2, path_anat, otheranat, dir_fMRI_Refth_RS_prepro1, root_RS, ID, Session, deoblique_exeption1, deoblique_exeption2, deoblique, orientation, overwrite)

    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### 
    ### ### #### ### coregistration of each run to the norm ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ##
    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### 

    for r in range(0, int(nb_run)):
        root_RS = extract_filename(RS[r])
        root_RS_ref = extract_filename(RS[REF_int])

        for f in os.listdir(opj(dir_fMRI_Refth_RS_prepro1,'tmp')):
            os.remove(os.path.join(opj(dir_fMRI_Refth_RS_prepro1,'tmp'), f))

        if not recordings == 'very_old':  # do not process ref...

            # 1.0 calculate co-registration mean image ref to mean image func
            command = 'antsRegistration -d 3 --float 0 --verbose 1 -w [0.01,0.99] -n ' + n_for_ANTS + \
            ' -o [' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp_') + ',' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp.nii.gz') + ']' + \
            ' -t Rigid[0.1] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
            ' -m MI[' + opj(dir_fMRI_Refth_RS_prepro1, root_RS_ref + '_xdtr_mean_deob_ref_fudge.nii.gz') + ',' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_deob.nii.gz') + ',1,32,Regular,0.2]' + \
            ' -t Affine[0.1] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
            ' -m MI[' + opj(dir_fMRI_Refth_RS_prepro1, root_RS_ref + '_xdtr_mean_deob_ref_fudge.nii.gz') + ',' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_deob.nii.gz') + ',1,32,Regular,0.2]' + \
            ' -t Syn[0.1,3,0] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
            ' -m CC[' + opj(dir_fMRI_Refth_RS_prepro1, root_RS_ref + '_xdtr_mean_deob_ref_fudge.nii.gz') + ',' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_deob.nii.gz') + ',1,4,Regular,0.2]'
            spco([command], shell=True)

            # 2.0 apply to all the volume in the func (_xdtr_deob)
            command = '3dTsplit4D' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'tmp','dummy.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_deob.nii.gz')
            spco([command], shell=True)
            WW = sorted(glob.glob(opj(dir_fMRI_Refth_RS_prepro1,'tmp','*.nii.gz')))
            for j in range(0, len(WW)):
                if j < 10:
                    command = 'antsApplyTransforms -d 3 -e 3 -i ' + WW[j] + ' -n NearestNeighbor' + \
                    ' -r ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp.nii.gz') + \
                    ' -o ' + opj(dir_fMRI_Refth_RS_prepro1,'tmp','warp1.00' + str(j) + '.nii.gz') + \
                    ' -t ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp_0GenericAffine.mat') + \
                    ' -t ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp_1Warp.nii.gz')
                elif j < 100 and j >= 10:
                    command = 'antsApplyTransforms -d 3 -e 3 -i ' + WW[j] + ' -n NearestNeighbor' + \
                    ' -r ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp.nii.gz') + \
                    ' -o ' + opj(dir_fMRI_Refth_RS_prepro1,'tmp', 'warp1.0' + str(j) + '.nii.gz') + \
                    ' -t ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp_0GenericAffine.mat') + \
                    ' -t ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp_1Warp.nii.gz')
                else:
                    command = 'antsApplyTransforms -d 3 -e 3 -i ' + WW[j] + ' -n NearestNeighbor' + \
                    ' -r ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp.nii.gz') + \
                    ' -o ' + opj(dir_fMRI_Refth_RS_prepro1,'tmp', 'warp1.' + str(j) + '.nii.gz') + \
                    ' -t ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp_0GenericAffine.mat') + \
                    ' -t ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp_1Warp.nii.gz')
                spco([command], shell=True)

            list_name = sorted(glob.glob(opj(dir_fMRI_Refth_RS_prepro1,'tmp','warp1.*.nii.gz')))
            name = ''
            for m in range(0,len(list_name)):
                name = name + ' ' + list_name[m]

            ##concat in one volume
            command = '3dTcat' + overwrite + ' ' + name + ' -tr ' + str(TR) + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref.nii.gz')
            spco([command], shell=True)

        else:
            command = '3dcopy ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS_ref + '_xdtr_deob.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref.nii.gz') + overwrite
            spco([command], shell=True)
