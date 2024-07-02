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
import ants

def prepar_aseg(Ref_file, labels_dir, volumes_dir, masks_dir, dir_transfo, BASE_SS_mask, BASE_SS_coregistr, Aseg_refLR, Aseg_ref, type_norm, ID, transfo_concat, dir_prepro, list_atlases, creat_sutdy_template, dir_out, diratlas_orig, check_visualy_each_img, n_for_ANTS, overwrite):
    print("BASE_SS_mask is " + BASE_SS_mask)
    print("BASE_SS_coregistr is " + BASE_SS_coregistr)
    print("Ref_file is " + Ref_file)

    ################################################################################################
    ########################## aseg and atlases into the native space ##############################
    ################################################################################################
    ####test on the template (atlas or sty) to see if it works
    ####save the old "Ref_file"
    command = '3dcalc -overwrite -a ' +  Ref_file + ' -expr "a" -prefix ' + opj(volumes_dir, ID + type_norm + '_brain_coregistr_QC.nii.gz')
    spco([command], shell=True)
    print(command)

    command = 'antsApplyTransforms -d 3 -i ' + BASE_SS_mask + \
    ' -r ' + Ref_file + \
    ' -o ' + opj(masks_dir,'brain_mask_in_anat_DC.nii.gz') + \
    transfo_concat + \
    ' -n NearestNeighbor'
    spco([command], shell=True)

    command = 'antsApplyTransforms -d 3 -i ' + BASE_SS_coregistr + \
    ' -r ' + Ref_file + \
    ' -o ' + opj(dir_prepro,'template_in_anat_DC.nii.gz') + \
    transfo_concat + \
    ' -n ' + n_for_ANTS
    spco([command], shell=True)
    print(command)

    command = '3dcalc -overwrite -a ' + Ref_file + ' -b ' + opj(dir_prepro,'template_in_anat_DC.nii.gz') + ' -expr "step(a)*b" -prefix ' + opj(dir_prepro,'template_in_anat_DC.nii.gz')
    spco([command], shell=True)

    #### apply to all atlases
    for atlas in list_atlases:
        if creat_sutdy_template==True:
            diratlas = dir_out
        else:
            diratlas = diratlas_orig
        command = 'antsApplyTransforms -d 3 -i ' + atlas + \
        ' -r ' + Ref_file + \
        ' -o ' + opj(labels_dir,type_norm + opb(atlas)) + \
        transfo_concat + \
        ' -n NearestNeighbor'
        spco([command], shell=True)

        command = '3dcalc -overwrite -a ' + Ref_file + ' -b ' + opj(labels_dir,type_norm + opb(atlas)) + ' -expr "step(a)*b" -prefix ' + opj(labels_dir,type_norm + opb(atlas))
        spco([command], shell=True)

    #### apply to asegLR
    command = 'antsApplyTransforms -d 3 -i ' + Aseg_refLR + \
        ' -r ' + Ref_file + \
        ' -o ' + opj(labels_dir, type_norm + 'Aseg_refLR.nii.gz') + \
        transfo_concat + \
        ' -n NearestNeighbor'
    spco([command], shell=True)

    command = '3dcalc -overwrite -a ' + Ref_file + ' -b ' + opj(labels_dir, type_norm + 'Aseg_refLR.nii.gz') + ' -expr "step(a)*b" -prefix ' + opj(labels_dir, type_norm + 'Aseg_refLR.nii.gz')
    spco([command], shell=True)   

    #### apply to aseg
    command = 'antsApplyTransforms -d 3 -i ' + Aseg_ref + \
    ' -r ' + Ref_file + \
    ' -o ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + \
    transfo_concat + \
    ' -n NearestNeighbor'
    spco([command], shell=True)

    command = '3dcalc -overwrite -a ' + Ref_file + ' -b ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "step(a)*b" -prefix ' + opj(labels_dir, type_norm + 'aseg.nii.gz')
    spco([command], shell=True)   

    #### check the result if you are not a fan of auto "like Simon Clavagnier =)"
    if check_visualy_each_img == True:
        #  Check for manual correction of the seg file!!!"
        # if manual correction save file as opj(labels_dir,'aseg_manual.nii.gz')!!!!!!!!!!!!!
        try:
            command = 'freeview -v ' +  Ref_file + ' ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ':colormap=heat:opacity=0.2:visible=1'
            spco([command], shell=True)
        except:
            pass

        while True:
            Your_choice = input('Enter "C" if you have corrected yourself the image or "OK" to continu: ')
            if Your_choice == '' or not Your_choice in ['C','OK']: 
                print('Please answer with C or OK!') 
            else: 
                break 

        if Your_choice=="OK":
            print("lucky boy!!!! =)")

        elif Your_choice=="C":

            ######New coregistration!!!! 
            command = 'antsRegistration -d 3 --float 0 --verbose 1 -u 1 -w [0.05,0.95]  -n ' + n_for_ANTS + \
                ' -o [' + opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_') + ',' + opj(labels_dir,'aseg_manual_for_transform.nii.gz') + ']' + \
                ' -t Affine[0.1] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
                ' -m MI[' + opj(labels_dir, type_norm + 'aseg_manual.nii.gz') + ',' + Aseg_ref + ',1,32,Regular,0.2]' + \
                ' -t Syn[0.1,3,0] -f 8x4x2x1 -s 3x2x1x0vox -c [500x200x50x10,1e-6,10]' + \
                ' -m CC[' + opj(labels_dir,type_norm + 'aseg_manual.nii.gz') + ',' + Aseg_ref + ',1,4,Regular,0.2]'
            spco([command], shell=True)

            command = '3dcalc -overwrite -a ' + opj(labels_dir, type_norm + 'aseg_manual_for_transform.nii.gz') + ' -expr "a" -prefix ' + opj(labels_dir, type_norm + 'aseg.nii.gz')
            spco([command], shell=True)

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
    ####??? prkoi pas N4?
    brain_img = ants.image_read(Ref_file)

    Sig_max = brain_img.max()
    normT1 = (brain_img/Sig_max)*110

    #creat CSF file
    command = '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "amongst(a' + CSF_edit + ')" -prefix ' + opj(labels_dir, type_norm + 'CSF_select.nii.gz')
    spco([command], shell=True)
    img_CSF = ants.image_read(opj(labels_dir, type_norm + 'CSF_select.nii.gz'))
    CSF     = ants.mask_image(normT1, img_CSF, 1)

    #creat WM file
    command = '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "amongst(a' + WM_edit + ')" -prefix ' + opj(labels_dir, type_norm + 'WM_select.nii.gz')
    spco([command], shell=True)
    img_WM = ants.image_read(opj(labels_dir, type_norm + 'WM_select.nii.gz'))
    WM     = ants.mask_image(normT1, img_WM, 1)

    #creat GM file
    command = '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "amongst(a' + GM_edit + ')" -prefix ' + opj(labels_dir, type_norm + 'GM_select.nii.gz')
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

    command = '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "amongst(a' + aseg_edit + ')*250" -prefix ' + opj(labels_dir, type_norm + 'aseg_edit.nii.gz')
    spco([command], shell=True)
    ASEG = ants.image_read(opj(labels_dir, type_norm + 'aseg_edit.nii.gz'))

    command = '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "amongst(a' + wm_aseg_edit + ')" -prefix ' + opj(labels_dir, type_norm + 'wm_edit.nii.gz')
    spco([command], shell=True)
    img_WM2  = ants.image_read(opj(labels_dir, type_norm + 'wm_edit.nii.gz'))
    WM_seg = ants.mask_image(pWM, img_WM2, 1)
    WM_FS  = ASEG + WM_seg
    ants.image_write(WM_FS, opj(labels_dir, type_norm + 'wm.nii.gz'), ri=False)

    command = '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'aseg.nii.gz') + ' -expr "amongst(a' + filled_edit + ')" -prefix ' + opj(labels_dir, type_norm + 'filled_mask.nii.gz')
    spco([command], shell=True)


    img_filled = ants.image_read(opj(labels_dir, type_norm + 'filled_mask.nii.gz'))
    img_filled = ants.iMath(img_filled,operation='GetLargestComponent')
    closed = ants.morphology(img_filled, operation='close', radius=2, mtype='binary', shape='ball')
    nii = ants.image_clone(img_filled)
    #nii = nii * 0
    #nii[0:143,:,:] = closed[0:143,:,:]*127
    #nii[143:,:,:]  = closed[143:,:,:]*255
    ants.image_write(nii, opj(labels_dir, type_norm + 'filled.nii.gz'), ri=False)
    #ants.image_write(thresh, filled_3, ri=False)
    command = '3dcalc' + overwrite + ' -a ' + opj(labels_dir, type_norm + 'Aseg_refLR.nii.gz') + ' -b ' + opj(labels_dir, type_norm + 'filled.nii.gz') + ' -expr "b*a" -prefix ' + opj(labels_dir, type_norm + 'filled.nii.gz') + ' -overwrite'
    spco([command], shell=True)

