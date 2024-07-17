import os
import subprocess
import ants

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

from fonctions.extract_filename import extract_filename

##### XXX add ICA or DL to visualize pre-processing effect

def preprocess_data(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, RS, list_RS, nb_run, T1_eq, overwrite,s_bind,afni_sif):


    if ope(dir_fMRI_Refth_RS_prepro1) == False:
        os.makedirs(dir_fMRI_Refth_RS_prepro1)

    if ope(os.path.join(dir_fMRI_Refth_RS_prepro1, 'tmp')) == False:
        os.makedirs(os.path.join(dir_fMRI_Refth_RS_prepro1, 'tmp'))


    for i in range(0, int(nb_run)):

        print('work on ' + str(dir_fMRI_Refth_RS_prepro1) + ' run ' + str(i))

        root = extract_filename(RS[i])
        #copy func imag
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + list_RS[i] + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '.nii.gz') + ' -expr "a"' + overwrite
        spco([command], shell=True)

        #command = '3drefit -space ORIG ' + opj(dir_fMRI_Refth_RS_prepro1, root + '.nii.gz')
        #spco([command], shell=True)

        base_fMRI = opj(dir_fMRI_Refth_RS_prepro1, root + '.nii.gz')
        # Clean bad volumes
        if ope(opj(opd(list_RS[i]), root, '.txt'))==True:
            # Open the file in read mode
            with open(opj(opd(list_RS[i]), root, '.txt'), 'r') as file:
                # Read the first line and convert it to an integer
                cut_low = int(file.readline().strip())

                # Read the second line and convert it to an integer
                cut_high = int(file.readline().strip())

            command = 'singularity run' + s_bind + afni_sif + '3dTcat -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_x0.nii.gz') + ' ' + base_fMRI + \
            '[' + str(cut_low) + '-' + str(cut_high-1) + ']' + overwrite
            spco([command], shell=True)
            base_fMRI = opj(dir_fMRI_Refth_RS_prepro1, root + '_x0.nii.gz')

        # T1 equilibrium
        cmd = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -nv ' + base_fMRI
        nb_vol = int(spgo(cmd).split('\n')[-1])


        command = 'singularity run' + s_bind + afni_sif + '3dTcat -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_x.nii.gz') + ' ' + base_fMRI + \
        '[' + str(T1_eq) + '-' + str(nb_vol-1) + ']' + overwrite
        spco([command], shell=True)
        
        # Despiking
        command = 'singularity run' + s_bind + afni_sif + '3dDespike -NEW -nomask' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz') + \
        ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_x.nii.gz')
        spco([command], shell=True)
        
        # slice-timing correction -heptic!!!!!!
        # better be sure that the Dicom are in ascending mode ??? old option -TR ' + str(TR) + 's -tpattern XXX -slice ' + str(nslice-1) XXX
        command = 'singularity run' + s_bind + afni_sif + '3dTshift -heptic'  + overwrite + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz ') + \
        ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz')
        spco([command], shell=True)

        command = 'singularity run' + s_bind + afni_sif + '3dTstat' + overwrite + ' -mean -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt_mean.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz')
        spco([command], shell=True)

        #outlier fraction for each volume
        command = '3dToutcount' + overwrite + ' -automask -fraction -polort 4 -legendre ' + \
        opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz ') + ' > ' + \
        opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt_outcount.r$run.1D')
        spco([command], shell=True)

        # realignment intra-run (volreg)
        '''
        command = 'singularity run' + s_bind + fsl_sif + 'mcflirt -in ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz ') + \
        ' -out ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr.nii.gz ') + \
        ' -mats -plots -reffile ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt_mean.nii.gz') + ' -rmsrel -rmsabs -spline_final'
        spco([command], shell=True)
        # output mvt image in png (deplacement in absolute and relative)
        command = 'singularity run' + s_bind + fsl_sif + 'fsl_tsplot -i ' + opj(dir_fMRI_Refth_RS_prepro1,root + '_xdtr.nii.gz_abs.rms') + \
        ',' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr.nii.gz_rel.rms') + \
        ' -t "MCFLIRT estimated mean displacement (mm)" -u 1 -w 640 -h 144 -a absolute,relative' + \
        ' -o ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr.png')
        spco([command], shell=True)

        '''

        # register each volume to the base image
        command = 'singularity run' + s_bind + afni_sif + '3dvolreg' + overwrite + ' -verbose -zpad 1 -base ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt_mean.nii.gz') + \
        ' -1Dfile ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_dfile.1D') + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr.nii.gz ') + \
        ' -cubic' + \
        ' -1Dmatrix_save ' + opj(dir_fMRI_Refth_RS_prepro1, root + '.aff12.1D') + ' ' + \
        opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz')
        spco([command], shell=True)

        # censoring # see ex 10 in 1d_tool
        command = 'singularity run' + s_bind + afni_sif + '1d_tool.py' + overwrite + ' -infile ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_dfile.1D') + \
        ' -derivative -censor_prev_TR -collapse_cols euclidean_norm' + \
        ' -moderate_mask -1.2 1.2 -show_censor_count' + \
        ' -write_censor ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_censor.1D') + \
        ' -write_CENSORTR ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_CENSORTR.txt')
        spco([command], shell=True)

        # compute motion magnitude time series: the Euclidean norm
        # (sqrt(sum squares)) of the motion parameter derivatives
        command = 'singularity run' + s_bind + afni_sif + '1d_tool.py' + overwrite + ' -infile ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_dfile.1D') + ' -set_nruns 1' + \
        ' -derivative -collapse_cols euclidean_norm ' + \
        '-write ' + opj(dir_fMRI_Refth_RS_prepro1, root + 'motion_enorm.1D')
        spco([command], shell=True)

        #writing regressors # get the first derivative
        command = 'singularity run' + s_bind + afni_sif + '1d_tool.py' + overwrite + ' -infile ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_dfile.1D') + \
        ' -derivative -write ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_deriv.1D')
        spco([command], shell=True)

        #writing regressors # get the demean
        command = 'singularity run' + s_bind + afni_sif + '1d_tool.py' + overwrite + ' -infile ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_dfile.1D') + \
        ' -demean -write ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_demean.1D')
        spco([command], shell=True)

        command = 'singularity run' + s_bind + afni_sif + '3dTstat' + overwrite + ' -mean -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_mean_preWARP.nii.gz') + \
        ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr.nii.gz')
        spco([command], shell=True)

        #BiasFieldCorrection
        IMG = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_mean_preWARP.nii.gz'))
        N4 = ants.n4_bias_field_correction(IMG,
                                           shrink_factor=4,
                                           convergence={'iters': [50, 50, 50, 50], 'tol': 1e-07},
                                           spline_param=200)
        ants.image_write(N4, opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_mean.nii.gz'), ri=False)


    # for distortion correction    


