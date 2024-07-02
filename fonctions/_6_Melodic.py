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


############################to do 
############################the end (ICA and normalisation to the same space)

def Melodic_correct(dir_RS_ICA_native_PreTT, dir_RS_ICA_native, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, 
    nb_ICA_run, nb_run, RS, TfMRI, overwrite):

    #   6.0 first pass : 'individual runs' ICA analysis
    if not ope(dir_RS_ICA_native_PreTT): os.mkdir(dir_RS_ICA_native_PreTT)
    if not ope(dir_RS_ICA_native): os.mkdir(dir_RS_ICA_native)
    # Normalized: use fsl_regfilt to remove noises.

    for i in range(0, int(nb_run)):
        root_RS = extract_filename(RS[i])

        command = 'melodic -i ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref_RcT_masked.nii.gz') + \
        ' -o ' + opj(dir_RS_ICA_native, root_RS + '_residual.ica') + \
        ' --tr=' + str(TR) + \
        ' --bgimage=' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
        ' --nobet --bgthreshold=10 -a symm -v --report --guireport=' + opj(dir_RS_ICA_native,'report.html') + \
        ' --mask=' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz') + ' --Oall -d ' + str(nb_ICA_run) + ' --mmthresh=0.66'
        spco([command], shell=True)


        if ope(opj(dir_RS_ICA_native, root_RS + '_residual.ica', 'melodic_select.txt'))==True:
            with open(opj(dir_RS_ICA_native, root_RS + '_residual.ica'), 'melodic_select.txt') as readtxt:
                readtxt_f = readtxt.read()
            print('found txt file of the ICA to remove!')

        else:
            command = 'gnome-terminal -- itksnap -g ' +  opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz') + \
            ' -o ' + opj(dir_RS_ICA_native, root_RS + '_residual.ica', 'melodic_IC.nii.gz')
            spco([command], shell=True)

            f = []
            for nb in list(range(1, nb_ICA+1)):
                while True:
                    Your_choice = input('Is ' + str(nb) + ' an artefact ?')
                    if Your_choice == '' or not Your_choice in ['y','n']: 
                        print('Please answer with y or n !') 
                    else: 
                        break 
                        
                if Your_choice=="y":
                    f.append(nb-1)
                    print("to be removed")
                elif Your_choice=="n":
                    print("let's keep it")

            f = ','.join(str(e) for e in f)
            #open text file
            with open(opj(dir_RS_ICA_native, root_RS + '_residual.ica', 'melodic_select.txt'), "w") as file:  
                file.write(f)  

            with open(opj(dir_RS_ICA_native, root_RS + '_residual.ica', 'melodic_select.txt')) as readtxt:
                readtxt_f = readtxt.read()

        print('processing fsl regflirt !! Close itksnap please !')

        #-a switch on aggressive filtering (full instead of partial regression)
        command = 'fsl_regfilt -a -i ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref_RcT_masked.nii.gz') + \
        ' -d ' + opj(dir_RS_ICA_native, root_RS + '_residual.ica','melodic_mix') + \
        ' -o ' + opj(dir_RS_ICA_native, root_RS + '_norm_final_clean.nii.gz') + ' -f ' + str(readtxt_f)
        spco([command], shell=True)

        #DATA.append(opj(dir_RS_ICA_native, root_RS + '_norm_final_clean.nii.gz'))


    DATA = sorted(glob.glob(opj(dir_RS_ICA_native, root_RS + '*_norm_final_clean.nii.gz')))
    name  = ''
    DATA2 = ''
    for m in range(0,len(DATA)):
        name  = name  + ',' + DATA[m]
        DATA2 = DATA2 + ' ' + DATA[m]

    command = 'melodic -i ' + name[1:] + \
    ' -o ' + opj(dir_RS_ICA_native_PreTT, root_RS + '_residual.ica') + \
    ' --tr=' + str(TR) + \
    ' --bgimage=' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
    ' --nobet --bgthreshold=10 -a symm -v --report --guireport=' + opj(dir_RS_ICA_native_PreTT,'report.html') + \
    ' --mask=' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz') + ' --Oall -d ' + str(nb_ICA_run) + ' --mmthresh=0.66'
    spco([command], shell=True)
