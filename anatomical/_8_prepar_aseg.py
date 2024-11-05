import os
import subprocess
import ants
from nilearn import plotting

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


def prepar_aseg(IgotbothT1T2, Ref_file, labels_dir, volumes_dir, masks_dir, dir_transfo, BASE_SS_mask, BASE_SS_coregistr, Aseg_refLR, Aseg_ref,
                type_norm, ID, transfo_concat,w2inv_fwd, dir_prepro, list_atlases, check_visualy_each_img, n_for_ANTS, otheranat, overwrite, bids_dir, Session,
                s_bind,afni_sif,itk_sif):

    print(bcolors.OKGREEN + "INFO: template mask is " + BASE_SS_mask + bcolors.ENDC)
    print(bcolors.OKGREEN + "INFO: template is " + BASE_SS_coregistr + bcolors.ENDC)
    print(bcolors.OKGREEN + "INFO: anat Ref_file is " + Ref_file + bcolors.ENDC)

    ################################################################################################
    ########################## aseg and atlases into the native space ##############################
    ################################################################################################

    ####save the old "Ref_file"    ####test on the template (atlas or sty) to see if it works
    command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + opj(volumes_dir,ID + type_norm + '_brain_step_1.nii.gz') + ' -expr "a" -prefix ' + Ref_file
    spco([command], shell=True)

    brain_img  = ants.image_read(Ref_file)
    MSK = ants.image_read(BASE_SS_mask)
    IMG  = ants.image_read(BASE_SS_coregistr)
    TRANS = ants.apply_transforms(fixed=brain_img, moving=MSK,
                                  transformlist=transfo_concat, interpolator='nearestNeighbor',whichtoinvert=w2inv_fwd)
    ants.image_write(TRANS, opj(masks_dir,'brain_mask_in_anat_DC.nii.gz'), ri=False)
    TRANS = ants.apply_transforms(fixed=brain_img, moving=IMG,
                                  transformlist=transfo_concat, interpolator=n_for_ANTS,whichtoinvert=w2inv_fwd)
    ants.image_write(TRANS, opj(dir_prepro,'template_in_anat_DC.nii.gz'), ri=False)
    command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + Ref_file + ' -b ' + opj(dir_prepro,'template_in_anat_DC.nii.gz') + ' -expr "step(a)*b" -prefix ' + opj(dir_prepro,'template_in_anat_DC.nii.gz')
    spco([command], shell=True)

    command = 'singularity run' + s_bind + afni_sif + '3dresample' + \
          ' -master ' + opj(dir_prepro, ID + '_acpc_test_QC_' + type_norm + '.nii.gz') + \
          ' -prefix ' + opj(masks_dir, type_norm + 'brain_mask_final_QCrsp.nii.gz') + \
          ' -input ' + opj(masks_dir,'brain_mask_in_anat_DC.nii.gz') + ' -overwrite'
    spco([command], shell=True)

    if IgotbothT1T2 == True:
        command = 'singularity run' + s_bind + afni_sif + '3dresample' + \
                  ' -master ' + opj(dir_prepro, ID + '_acpc_test_QC_' + otheranat + '.nii.gz') + \
                  ' -prefix ' + opj(masks_dir, otheranat + 'brain_mask_final_QCrsp.nii.gz') + \
                  ' -input ' + opj(masks_dir,'brain_mask_in_anat_DC.nii.gz') + ' -overwrite'
        spco([command], shell=True)

    #### apply to all atlases
    if list_atlases:
        for atlas in list_atlases:
            if ope(atlas):
                IMG = ants.image_read(atlas)
                if opb(atlas) in ['Vmask.nii.gz', 'Wmask.nii.gz', 'Gmask.nii.gz']:
                    dir_out = masks_dir
                else:
                    dir_out = labels_dir

                TRANS = ants.apply_transforms(fixed=brain_img, moving=IMG,
                                              transformlist=transfo_concat, interpolator='nearestNeighbor',whichtoinvert=w2inv_fwd)
                ants.image_write(TRANS, opj(dir_out,type_norm + opb(atlas)), ri=False)
                command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + Ref_file + ' -b ' + opj(dir_out,type_norm + opb(atlas)) + ' -expr "step(a)*b" -prefix ' + opj(labels_dir,type_norm + opb(atlas))
                spco([command], shell=True)
                print(bcolors.OKGREEN + 'INFO: done with atlas in subject space! you should check that' + bcolors.ENDC)
            else:
                print(bcolors.WARNING + 'WARNING: no atlas found for path: ' + str(atlas) + bcolors.ENDC)


        if not os.path.exists(bids_dir + '/QC/FinalQC'):
            os.mkdir(bids_dir + '/QC/FinalQC')
        try:
            display = plotting.plot_anat(Ref_file, display_mode='mosaic', dim=4)
            display.add_contours(opj(dir_prepro,'template_in_anat_DC.nii.gz'),
            linewidths=.2, colors=['red'])
            display.savefig(bids_dir + '/QC/FinalQC/' + ID + '_' + str(Session) + '_' + type_norm + '_final_template_to_anat.png')
            # Don't forget to close the display
            display.close()
        except:
            display = plotting.plot_anat(Ref_file, display_mode='mosaic', dim=4)
            display.savefig(bids_dir + '/QC/FinalQC/' + ID + '_' + str(Session) + '_' + type_norm + '_final_template_to_anat.png')
            # Don't forget to close the display
            display.close()

    #### apply to asegLR
    if ope(Aseg_refLR):
        print(bcolors.OKGREEN + 'INFO: found Aseg_refLR:' + str(Aseg_refLR) + bcolors.ENDC)
        IMG = ants.image_read(Aseg_refLR)
        TRANS = ants.apply_transforms(fixed=brain_img, moving=IMG,
                                      transformlist=transfo_concat, interpolator='nearestNeighbor',whichtoinvert=w2inv_fwd)
        ants.image_write(TRANS, opj(labels_dir, type_norm + 'Aseg_refLR.nii.gz'), ri=False)
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + Ref_file + ' -b ' + opj(labels_dir, type_norm + 'Aseg_refLR.nii.gz') + ' -expr "step(a)*b" -prefix ' + opj(labels_dir, type_norm + 'Aseg_refLR.nii.gz')
        spco([command], shell=True)
    else:
        print(bcolors.WARNING + 'WARNING: no Aseg_refLR found' + bcolors.ENDC)

    if ope(Aseg_ref):
        #### apply to aseg
        IMG = ants.image_read(Aseg_ref)
        TRANS = ants.apply_transforms(fixed=brain_img, moving=IMG,
                                      transformlist=transfo_concat, interpolator='nearestNeighbor',whichtoinvert=w2inv_fwd)
        ants.image_write(TRANS, opj(labels_dir, type_norm + 'aseg.nii.gz'), ri=False)
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + Ref_file + ' -b ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "step(a)*b" -prefix ' + opj(labels_dir, type_norm + 'aseg.nii.gz')
        spco([command], shell=True)
    else:
        print(bcolors.WARNING + 'WARNING: no Aseg_ref found' + bcolors.ENDC)

    if ope(opj(labels_dir, type_norm + 'aseg.nii.gz')):
        #### check the result if you are not a fan of auto "like Simon Clavagnier
        if check_visualy_each_img == True:
            #  Check for manual correction of the seg file!!!
            def run_command_and_wait(command):
                print(bcolors.OKGREEN + 'INFO: Running command:', command + bcolors.ENDC)
                result = subprocess.run(command, shell=True)
                if result.returncode == 0:
                    print(bcolors.OKGREEN + 'INFO: Command completed successfully.' + bcolors.ENDC)
                else:
                    print(bcolors.WARNING + 'WARNING: Command failed with return code:', result.returncode, bcolors.ENDC)
            # Example usage
            print(bcolors.WARNING + 'WARNING: check aseg and if you correct it save it where it as ' + opj(labels_dir,'aseg_manual.nii.gz') +  + bcolors.ENDC)
            command = ('singularity run' + s_bind + itk_sif + 'itksnap -g ' + Ref_file + ' -s ' + opj(labels_dir, type_norm + 'aseg.nii.gz'))
            run_command_and_wait(command)

            if ope(opj(labels_dir, type_norm + 'aseg_manual.nii.gz')):
                ######New coregistration!!!!
                manual_IMG = ants.image_read(opj(labels_dir, type_norm + 'aseg_manual.nii.gz'))

                mTx = ants.registration(fixed=manual_IMG, moving=IMG,
                                        type_of_transform='SyN',
                                        outprefix=opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_'));

                TRANS = ants.apply_transforms(fixed=manual_IMG, moving=IMG,
                                              transformlist=mTx['fwdtransform'], interpolator='nearestNeighbor')
                ants.image_write(TRANS, opj(labels_dir,'aseg_manual_for_transform.nii.gz'), ri=False)

                command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + opj(labels_dir, type_norm + 'aseg_manual_for_transform.nii.gz') + ' -expr "a" -prefix ' + opj(labels_dir, type_norm + 'aseg.nii.gz')
                spco([command], shell=True)

    if ope(opj(labels_dir, type_norm + 'aseg.nii.gz')):
        ####################################################################################
        ########################## modif aseg to seg file ##################################
        ####################################################################################
        ####################################################################################
        #define the aseg number
        # for wm.mgz
        aseg_edit    = ',4,43,11,12,26,28,50,51,58,60,31,63,4,43,77,10,49'
        wm_aseg_edit = ',2,41,16,13,52'
        filled_edit   = ',4,43,11,12,26,28,50,51,58,60,31,63,4,43,77,10,49,2,41,13,52'

        CSF_edit = ',4,5,43,44,14,15,24,72,30,62,31,63'
        WM_edit  = ',2,41,85,7,46,251,252,253,254,255'
        GM_edit  = ',3,42,8,47,10,49,11,12,13,16,17,18,26,28,50,51,52,53,54,58,60'

        new_min_WM = 80
        new_max_WM = 110

        ####################################################################################

        Sig_max = brain_img.max()
        normT1 = (brain_img/Sig_max)*110

        #creat CSF file
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "amongst(a' + CSF_edit + ')" -prefix ' + opj(labels_dir, type_norm + 'CSF_select.nii.gz')
        spco([command], shell=True)
        img_CSF = ants.image_read(opj(labels_dir, type_norm + 'CSF_select.nii.gz'))
        CSF     = ants.mask_image(normT1, img_CSF, 1)

        #creat WM file
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "amongst(a' + WM_edit + ')" -prefix ' + opj(labels_dir, type_norm + 'WM_select.nii.gz')
        spco([command], shell=True)
        img_WM = ants.image_read(opj(labels_dir, type_norm + 'WM_select.nii.gz'))
        WM     = ants.mask_image(normT1, img_WM, 1)

        #creat GM file
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "amongst(a' + GM_edit + ')" -prefix ' + opj(labels_dir, type_norm + 'GM_select.nii.gz')
        spco([command], shell=True)

        img_GM  = ants.image_read(opj(labels_dir, type_norm + 'GM_select.nii.gz'))
        GM      = ants.mask_image(normT1, img_GM, 1)

        #attribute new value to your indiv image
        if 'T1' in type_norm:
            new_min_CSF = 1
            new_max_CSF = 50
            new_min_GM  = 40
            new_max_GM  = 90

        elif 'T2' in type_norm:
            new_min_CSF = 10
            new_max_CSF = 50
            new_min_GM  = 120
            new_max_GM  = 190

        nCSF  = (CSF - CSF.min()) / int(CSF.max() - CSF.min())
        pCSF  = nCSF * int(new_max_CSF - new_min_CSF ) + int(new_min_CSF)
        oCSF  = ants.mask_image(pCSF, img_CSF, 1)

        nGM = (GM - GM.min()) / int(GM.max() - GM.min())
        pGM = nGM * int(new_max_GM - new_min_GM ) + int(new_min_GM)
        oGM = ants.mask_image(pGM, img_GM, 1)

        nWM = (WM - WM.min()) / int(WM.max() - WM.min())
        pWM = nWM * int(new_max_WM - new_min_WM ) + int(new_min_WM)
        oWM = ants.mask_image(pWM, img_WM, 1)

        nii = oCSF + oGM + oWM
        ants.image_write(nii, Ref_file.replace('.nii.gz', '_norm_' + type_norm + '.nii.gz'), ri=False)

        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "amongst(a' + aseg_edit + ')*250" -prefix ' + opj(labels_dir, type_norm + 'aseg_edit.nii.gz')
        spco([command], shell=True)
        ASEG = ants.image_read(opj(labels_dir, type_norm + 'aseg_edit.nii.gz'))

        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "amongst(a' + wm_aseg_edit + ')" -prefix ' + opj(labels_dir, type_norm + 'wm_edit.nii.gz')
        spco([command], shell=True)
        img_WM2  = ants.image_read(opj(labels_dir, type_norm + 'wm_edit.nii.gz'))
        WM_seg = ants.mask_image(pWM, img_WM2, 1)
        WM_FS  = ASEG + WM_seg
        ants.image_write(WM_FS, opj(labels_dir, type_norm + 'wm.nii.gz'), ri=False)

        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "amongst(a' + filled_edit + ')" -prefix ' + opj(labels_dir, type_norm + 'filled_mask.nii.gz')
        spco([command], shell=True)

        img_filled = ants.image_read(opj(labels_dir, type_norm + 'filled_mask.nii.gz'))
        img_filled = ants.iMath(img_filled,operation='GetLargestComponent')
        nii = ants.image_clone(img_filled)
        ants.image_write(nii, opj(labels_dir, type_norm + 'filled.nii.gz'), ri=False)
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'Aseg_refLR.nii.gz') + ' -b ' + opj(labels_dir, type_norm + 'filled.nii.gz') + ' -expr "b*a" -prefix ' + opj(labels_dir, type_norm + 'filled.nii.gz') + ' -overwrite'
        spco([command], shell=True)

