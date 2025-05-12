import os
import subprocess
import ants
import datetime
import json
import numpy as np
import nibabel as nib
import nilearn

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


# Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

spco = subprocess.check_output
spgo = subprocess.getoutput
from fonctions.extract_filename import extract_filename

##### XXX add ICA or DL to visualize pre-processing effect

def preprocess_data(dir_fMRI_Refth_RS_prepro1, RS, list_RS, nb_run, T1_eq, TR, Slice_timing_info,
                    overwrite, s_bind, afni_sif, diary_file):
    if ope(dir_fMRI_Refth_RS_prepro1) == False:
        os.makedirs(dir_fMRI_Refth_RS_prepro1)

    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  Working on step ' + str(1) + '(function: _1_fMRI_preTTT_in_fMRIspace).  ##'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    for i in range(0, int(nb_run)):
        nl = 'work on ' + str(dir_fMRI_Refth_RS_prepro1) + ' run ' + str(i +1)
        print(bcolors.OKGREEN + nl + bcolors.ENDC)
        diary.write(f'\n{nl}')
        root = extract_filename(RS[i])

        # copy func imag
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + list_RS[i] + ' -prefix ' + \
                  opj(dir_fMRI_Refth_RS_prepro1, root + '.nii.gz') + ' -expr "a"' + overwrite
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": list_RS[i],
                      "Description": 'Copy.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root + '.json'), "w") as outfile:
            outfile.write(json_object)

        base_fMRI = opj(dir_fMRI_Refth_RS_prepro1, root + '.nii.gz')

        # Clean bad volumes
        if ope(opj(opd(list_RS[i]), root, '.txt')) == True:
            # Open the file in read mode
            with open(opj(opd(list_RS[i]), root, '.txt'), 'r') as file:
                # Read the first line and convert it to an integer
                cut_low = int(file.readline().strip())

                # Read the second line and convert it to an integer
                cut_high = int(file.readline().strip())

            command = 'singularity run' + s_bind + afni_sif + '3dTcat -prefix ' + \
                      opj(dir_fMRI_Refth_RS_prepro1,root + '_x0.nii.gz') + ' ' + base_fMRI + \
                      '[' + str(cut_low) + '-' + str(cut_high - 1) + ']' + overwrite
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            base_fMRI = opj(dir_fMRI_Refth_RS_prepro1, root + '_x0.nii.gz')
            dictionary = {"Sources": [base_fMRI,
                                      opj(opd(list_RS[i]), root, '.txt')],
                          "Description": 'Remove volumes.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_x0.json'), "w") as outfile:
                outfile.write(json_object)

        img = nib.load(base_fMRI)
        nb_vol = img.shape[-1]  # For 4D data, last dimension is time

        command = 'singularity run' + s_bind + afni_sif + '3dTcat -prefix ' + \
                  opj(dir_fMRI_Refth_RS_prepro1,root + '_x.nii.gz') + ' ' + base_fMRI + \
                  '[' + str(T1_eq) + '-' + str(nb_vol - 1) + ']' + overwrite
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": base_fMRI,
                      "Description": 'Remove first volumes.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_x.json'), "w") as outfile:
            outfile.write(json_object)

        # Despiking
        command = ('singularity run' + s_bind + afni_sif + '3dDespike -NEW -nomask' + overwrite + ' -prefix ' + \
            opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz') + \
            ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_x.nii.gz'))
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_x.nii.gz'),
                      "Description": 'Despiking.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.json'), "w") as outfile:
            outfile.write(json_object)

        # slice-timing correction -heptic!!!!!!
        if Slice_timing_info == 'Auto':
            if opi(opj(dir_fMRI_Refth_RS_prepro1, 'stc.txt')):
                command = 'singularity run' + s_bind + afni_sif + '3dTshift -wsinc9' + overwrite + \
                          ' -TR ' + str(TR) + ' -tpattern @' + opj(dir_fMRI_Refth_RS_prepro1, 'stc.txt') + \
                          ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz ') + \
                          ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz')
                nl = spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
                dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz'),
                              "Description": 'Slice timing correction.', }
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.json'), "w") as outfile:
                    outfile.write(json_object)

            else:
                cmd = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -slice_timing ' + base_fMRI
                nl = spgo(cmd).split('\n')
                STC = nl[-1].split('|')
                STC = list(map(float, STC))
                if np.sum(STC) > 0:
                    # means that AFNI has access to the slice timing in the nifti header
                    command = 'singularity run' + s_bind + afni_sif + '3dTshift -wsinc9' + overwrite + \
                              ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz ') + \
                              ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz')
                    nl = spgo(command)
                    diary.write(f'\n{nl}')
                    print(nl)
                    dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz'),
                                  "Description": 'Slice timing correction.', }
                    json_object = json.dumps(dictionary, indent=2)
                    with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.json'), "w") as outfile:
                        outfile.write(json_object)
                else:
                    nl = "WARNING: Slice Timing not found, this will be particularly DANGEROUS, you SHOULD PROVIDE MANUALLY ONE!"
                    print(bcolors.WARNING + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    diary_file_WARNING = opj(opd(opd(dir_fMRI_Refth_RS_prepro1)), 'MAJOR_WARNING.txt')
                    if not opi(diary_file_WARNING):
                        diary_file_WARNING_file = open(diary_file_WARNING, "w")
                        diary_file_WARNING_file.write(f'\n{nl}')
                    else:
                        diary_file_WARNING_file = open(diary_file_WARNING, "a")
                        diary_file_WARNING_file.write(f'\n{nl}')
                    diary_file_WARNING_file.close()

                    command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz') + \
                              ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz') + ' -expr "a"' + overwrite
                    nl = spgo(command)
                    diary.write(f'\n{nl}')
                    print(nl)
                    dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz'),
                                  "Description": 'copy.', }
                    json_object = json.dumps(dictionary, indent=2)

                    with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.json'), "w") as outfile:
                        outfile.write(json_object)
        elif isinstance(Slice_timing_info, list) == False and Slice_timing_info.split(' ')[0] == '-tpattern':
            command = 'singularity run' + s_bind + afni_sif + '3dTshift -wsinc9' + overwrite + \
                      ' -TR ' + str(TR) + ' ' + Slice_timing_info + \
                      ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz ') + \
                      ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz')
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz'),
                          "Description": 'Slice timing correction.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.json'), "w") as outfile:
                outfile.write(json_object)
        elif isinstance(Slice_timing_info, list) == True:
            command = 'singularity run' + s_bind + afni_sif + '3dTshift -wsinc9' + overwrite + \
                      ' -TR ' + str(TR) + ' -tpattern @' + opj(dir_fMRI_Refth_RS_prepro1, 'stc.txt') + \
                      ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz ') + \
                      ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz')
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz'),
                          "Description": 'Slice timing correction.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.json'), "w") as outfile:
                outfile.write(json_object)
        else:
            nl = 'ERROR : please check Slice_timing_info, this variable is not define as it should'
            diary.write(f'\n{nl}')
            raise ValueError(bcolors.FAIL + nl + bcolors.ENDC)

        command = 'singularity run' + s_bind + afni_sif + '3dTstat' + overwrite + ' -mean -prefix ' + \
                  opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt_mean.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz')
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz'),
                      "Description": 'Mean image.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt_mean.json'), "w") as outfile:
            outfile.write(json_object)

        # outlier fraction for each volume
        command = 'singularity run' + s_bind + afni_sif + '3dToutcount' + overwrite + ' -automask -fraction -polort 4 -legendre ' + \
                  opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz ') + ' > ' + \
                  opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt_outcount.r$run.1D')
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        # realignment intra-run (volreg)
        # register each volume to the base image
        command = 'singularity run' + s_bind + afni_sif + '3dvolreg' + overwrite + ' -verbose -zpad 1 -base ' + \
                  opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt_mean.nii.gz') + \
                  ' -1Dfile ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_dfile.1D') + \
                  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr.nii.gz ') + \
                  ' -cubic' + \
                  ' -twodup ' + \
                  ' -1Dmatrix_save ' + opj(dir_fMRI_Refth_RS_prepro1, root + '.aff12.1D') + ' ' + \
                  opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz')
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        #with stat values
        extracted_data = nib.load(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr.nii.gz')).get_fdata()
        labeled_img2 = nilearn.image.new_img_like(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz'), extracted_data, copy_header=True )
        labeled_img2.to_filename(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr.nii.gz'))

        dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz'),
                                  opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt_mean.nii.gz')],
                      "Description": 'Rigid realignment (3dVolreg from AFNI).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr.json'), "w") as outfile:
            outfile.write(json_object)

        # censoring # see ex 10 in 1d_tool
        command = 'singularity run' + s_bind + afni_sif + '1d_tool.py' + overwrite + ' -infile ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_dfile.1D') + \
                  ' -derivative -censor_prev_TR -collapse_cols euclidean_norm' + \
                  ' -moderate_mask -1.2 1.2 -show_censor_count' + \
                  ' -write_censor ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_censor.1D') + \
                  ' -write_CENSORTR ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_CENSORTR.txt')
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        # compute motion magnitude time series: the Euclidean norm
        # (sqrt(sum squares)) of the motion parameter derivatives
        command = 'singularity run' + s_bind + afni_sif + '1d_tool.py' + overwrite + ' -infile ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_dfile.1D') + ' -set_nruns 1' + \
                  ' -derivative -collapse_cols euclidean_norm ' + \
                  '-write ' + opj(dir_fMRI_Refth_RS_prepro1, root + 'motion_enorm.1D')
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        # writing regressors # get the first derivative
        command = 'singularity run' + s_bind + afni_sif + '1d_tool.py' + overwrite + ' -infile ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_dfile.1D') + \
                  ' -derivative -write ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_deriv.1D')
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        # writing regressors get demean
        command = 'singularity run' + s_bind + afni_sif + '1d_tool.py' + overwrite + ' -infile ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_dfile.1D') + \
                  ' -demean -write ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_demean.1D')
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        command = 'singularity run' + s_bind + afni_sif + '3dTstat' + overwrite + ' -mean -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_mean_preWARP.nii.gz') + \
        ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr.nii.gz')
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr.nii.gz'),
                      "Description": 'Mean image.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_mean_preWARP.json'), "w") as outfile:
            outfile.write(json_object)

        # BiasFieldCorrection
        IMG = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_mean_preWARP.nii.gz'))
        N4 = ants.n4_bias_field_correction(IMG,
                                           shrink_factor=4,
                                           convergence={'iters': [50, 50, 50, 50], 'tol': 1e-07},
                                           spline_param=200)
        ants.image_write(N4, opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_mean.nii.gz'), ri=False)
        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_mean_preWARP.nii.gz'),
                      "Description": 'Bias field correction (N4).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_mean.json'), "w") as outfile:
            outfile.write(json_object)

    diary.write(f'\n')
    diary.close()