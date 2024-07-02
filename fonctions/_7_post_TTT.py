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



def signal_regression(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_RS_ICA_native,
    nb_run, RS, blur, TR, melodic_prior_post_TTT, extract_exterior_CSF, extract_WM, do_not_correct_signal, band, extract_Vc, overwrite):

    for i in range(0, int(nb_run)):
        root_RS = extract_filename(RS[i])

        if melodic_prior_post_TTT == True:
            command = '3dcalc' + overwrite + ' -a ' + opj(dir_RS_ICA_native, root_RS + '_norm_final_clean.nii.gz') + \
            ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz') + ' -expr "a"'
            spco([command], shell=True)

        else: 
            command = '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref_RcT_masked.nii.gz') + \
            ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz') + ' -expr "a"'
            spco([command], shell=True)

    # 4.2 Get the SVD values from the masks


    for i in range(0, int(nb_run)):
        root_RS = extract_filename(RS[i])

        if extract_exterior_CSF == True:
            command = '3dmaskSVD' + overwrite + ' -polort 2 -vnorm -mask ' + opj(dir_fMRI_Refth_RS_prepro1,'exterior_ligne.nii.gz') + \
            ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz') + \
            ' > ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS_NonB.1D')
            spco([command], shell=True)

        command = '3dBrickStat -count -non-zero ' + opj(dir_fMRI_Refth_RS_prepro1,'Wmask.nii.gz')
        count = spgo([command])

        command = '3dmaskSVD' + overwrite + ' -polort 2 -vnorm -mask ' + opj(dir_fMRI_Refth_RS_prepro1,'Wmask.nii.gz') + \
        ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz') + \
        ' > ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS_Wc.1D')
        spco([command], shell=True)

        command = '3dROIstats -nzvoxels -mask ' + opj(dir_fMRI_Refth_RS_prepro1,'Vmask.nii.gz') + ' ' +  opj(dir_fMRI_Refth_RS_prepro1,'Vmask.nii.gz')
        count = os.path.splitext(spgo([command]))[1][4:]
        if count == '0[?]':
            countVmask = 0
        else:
            countVmask = int(count.split('\t')[1])

        if countVmask>10:
            command = '3dmaskSVD' + overwrite + ' -polort 2 -vnorm -mask ' + opj(dir_fMRI_Refth_RS_prepro1,'Vmask.nii.gz') + \
            ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz') + \
            ' > ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS_Vc.1D')
            spco([command], shell=True)
        
        command = '3dTstat' + overwrite + ' -cvarinv -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS_tsnr1.nii.gz') + \
        ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz')
        spco([command], shell=True)
        command = '3dTstat' + overwrite + ' -cvar -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS_cvar.nii.gz') + \
        ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz')
        spco([command], shell=True)
        command = '3dTstat' + overwrite + ' -tsnr -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS_tsnr2.nii.gz') + \
        ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz')
        spco([command], shell=True)
        command = '3dTstat' + overwrite + ' -stdev -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS_stdev.nii.gz') + \
        ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz')
        spco([command], shell=True)


        ############################### ############################### ############################### 
        ############################### Corrections of the signal  ####################################
        ############################### ############################### ############################### 

        # 5.0 Regress out most of the noise from the data: bandpass filter, motion correction white mater noise and cbf noise , plus drift and derivatives
        # after filtering : blur within the mask and normalise the data.

        # T1 equilibrium
        command = '3dinfo -nv ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz')
        NumberofTR = str(int(spco([command], shell=True)))

        # create bandpass regressors (instead of using 3dBandpass, say)
        command = '1dBport' + overwrite + ' -nodata ' + NumberofTR + ' ' + str(TR) + ' -band ' + band + ' -invert -nozero > ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + 'bandpass_rall.1D')
        spco([command], shell=True)

        if ope(dir_RS_ICA_native) == False:
            os.makedirs(dir_RS_ICA_native)

    if do_not_correct_signal == False:
        for i in range(0, int(nb_run)):
            root_RS = extract_filename(RS[i])
            command = '3dDeconvolve -input ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz') + \
            ' -mask ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz') + \
            ' -ortvec ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + 'bandpass_rall.1D') + ' bandpass_rall'  + \
            ' -ortvec ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_demean.1D') + ' mot_demean'  + \
            ' -ortvec ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_deriv.1D') + ' mot_deriv' + \
            ' -censor ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_censor.1D') + \
            ' -polort A -float ' +                                                  \
            ' -num_stimts 0 ' + overwrite +                                                      \
            ' -fout -tout '  +                                                          \
            ' -x1D ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + 'X.xmat.1D ') +                                                       \
            ' -xjpeg ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + 'X.jpg')  +                                                    \
            ' -x1D_uncensored ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + 'X.nocensor.xmat.1D') +                                 \
            ' -fitts ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + 'Xfittssubj') +                                                 \
            ' -errts ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + 'errts') +                                                         \
            ' -x1D_stop ' +                                                          \
            ' -bucket ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + 'statssubj')

            if extract_exterior_CSF == True:
                print('t')
                command = command + ' -ortvec ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS_NonB.1D') + ' residual_norm_NonB '

            elif extract_WM == True:
                print('t')
                command = command + ' -ortvec ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS_Wc.1D') + ' residual_norm_Wc '

            elif extract_Vc == True:
                command = command + ' -ortvec ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS_Vc.1D') + ' residual_norm_Vc '

            spco([command], shell=True)

            command = '3dTproject -polort 0' + overwrite + ' -input ' + \
            opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz') + \
            ' -censor ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_censor.1D') + \
            ' -cenmode ZERO -ort ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + 'X.nocensor.xmat.1D') + \
            ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz') + ' -blur ' + str(blur)
            spco([command], shell=True)

    else:
        for i in range(0, int(nb_run)):
            root_RS = extract_filename(RS[i])
            command = '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz') + \
                      ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz') + ' -expr "a"'
            spco([command], shell=True)
