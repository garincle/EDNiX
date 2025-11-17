import os
import fnmatch

opj = os.path.join
opd = os.path.dirname
opb = os.path.basename

def get(specie,path_ATLAS,FS_tools,path_BALSA,reference,atlasname, REFname,type_norm,kind, list_atlas):

    # by default variables
    FS_refs       = opj(FS_tools, 'standard_mesh_atlases_macaque')
    BALSAname     = ''  # to get a surface to surface co-registration
    balsa_folder  = ''
    balsa_brainT1 = ''

    # search for the relevant folder
    for dirpath, dirnames, filenames in os.walk(opj(path_ATLAS, 'atlas')):
        folder_name = fnmatch.filter(dirnames, specie)
        if folder_name:
            path_ref = opj(dirpath, folder_name[0])

    path_label_code = opj(path_ref, 'label_code')
    if 'T1' in type_norm:
        suffix_template = 'T1w'
    if 'T2' in type_norm:
        suffix_template = 'T2w'

    # references ########################################################################################################

    BASE_folder        = opj(path_ref, reference, 'volumes')
    BASE_atlas_folder  = opj(BASE_folder, 'labels')
    BASE_masks_folder  = opj(BASE_folder, 'masks')
    BASE_priors_folder = opj(BASE_folder, 'priors')
    BASE_template      = opj(BASE_folder, specie + '_space-acpc_desc-template_' + suffix_template + '.nii.gz')
    BASE_SS            = opj(BASE_folder, specie + '_space-acpc_desc-SS_' + suffix_template + '.nii.gz')
    BASE_mask          = opj(BASE_masks_folder,  specie + '_mask.nii.gz')
    BASE_Gmask         = opj(BASE_masks_folder,  specie + '_desc-Gray_mask.nii.gz')
    BASE_Wmask         = opj(BASE_masks_folder,  specie + '_desc-erod-White_mask.nii.gz')
    BASE_Vmask         = opj(BASE_masks_folder,  specie + '_desc-erod-Vent_mask.nii.gz')
    CSF                = opj(BASE_priors_folder, specie + '_label-CSF_probseg.nii.gz')
    GM                 = opj(BASE_priors_folder, specie + '_label-GM_probseg.nii.gz')
    WM                 = opj(BASE_priors_folder, specie + '_label-WM_probseg.nii.gz')
    Aseg_ref           = opj(BASE_atlas_folder,  specie + '_seg-4FS_dseg.nii.gz')

    if specie == 'Human':
        FS_refs          = opj(FS_tools, 'standard_mesh_atlases')
        BALSAname        = 'S1200'
        balsa_folder     = opj(path_BALSA, BALSAname, 'modified')
        balsa_templateT1 = opj(balsa_folder, 'volumes', 'S1200_AverageT1w_restore.nii.gz')
        balsa_templateT2 = opj(balsa_folder, 'volumes', 'S1200_AverageT2w_restore.nii.gz')
        balsa_brainT1    = opj(balsa_folder, 'volumes', 'S1200_desc-SS_T1w.nii.gz')
        balsa_gray       = opj(balsa_folder, 'volumes', 'masks', 'S1200_desc-dilat-Gray_mask.nii.gz') # doesn't exist yet
        balsa_label      = opj(balsa_folder, 'volumes', 'labels')

        if reference == '':
            reference = BALSAname
            BASE_folder        = opj(balsa_folder, 'volumes')
            BASE_atlas_folder  = opj(BASE_folder, 'labels')
            BASE_masks_folder  = opj(BASE_folder, 'masks')
            BASE_priors_folder = opj(BASE_folder, 'priors')
            BASE_SS = balsa_brainT1
            BASE_mask = opj(BASE_masks_folder, BALSAname + '_mask.nii.gz')
            BASE_Gmask = opj(BASE_masks_folder, BALSAname + '_desc-Gray_mask.nii.gz')
            BASE_Wmask = opj(BASE_masks_folder, BALSAname + '_desc-erod-White_mask.nii.gz')
            BASE_Vmask = opj(BASE_masks_folder, BALSAname + '_desc-erod-Vent_mask.nii.gz')
            CSF = opj(BASE_priors_folder, BALSAname + '_label-CSF_probseg.nii.gz')
            GM = opj(BASE_priors_folder, BALSAname + '_label-GM_probseg.nii.gz')
            WM = opj(BASE_priors_folder, BALSAname + '_label-WM_probseg.nii.gz')
            Aseg_ref = opj(BASE_atlas_folder, BALSAname + '_seg-4FS_dseg.nii.gz')

    
    elif specie == 'Chimpanzee':
        FS_refs          = opj(FS_tools, 'standard_mesh_atlases_chimp')
        BALSAname        = 'CY29'
        balsa_folder     = opj(path_BALSA, BALSAname, 'modified')
        balsa_templateT1 = opj(balsa_folder, 'volumes', 'CY29_res-0.8_desc-template_T1w.nii.gz')
        balsa_templateT2 = opj(balsa_folder, 'volumes', 'CY29_res-0.8_desc-template_T2w.nii.gz')
        balsa_brainT1    = opj(balsa_folder, 'volumes', 'CY29_res-0.8_desc-SS_T1w.nii.gz')
        balsa_gray       = opj(balsa_folder, 'volumes', 'masks', 'CY29_desc-dilat-Gray_mask.nii.gz')  # doesn't exist yet
        balsa_label      = opj(balsa_folder, 'volumes', 'labels')

        if reference == '':
            reference = BALSAname
            BASE_folder = opj(balsa_folder, 'volumes')
            BASE_atlas_folder = opj(BASE_folder, 'labels')
            BASE_masks_folder = opj(BASE_folder, 'masks')
            BASE_priors_folder = opj(BASE_folder, 'priors')
            BASE_SS = balsa_brainT1
            BASE_mask = opj(BASE_masks_folder, BALSAname + '_mask.nii.gz')
            BASE_Gmask = opj(BASE_masks_folder, BALSAname + '_desc-Gray_mask.nii.gz')
            BASE_Wmask = opj(BASE_masks_folder, BALSAname + '_desc-erod-White_mask.nii.gz')
            BASE_Vmask = opj(BASE_masks_folder, BALSAname + '_desc-erod-Vent_mask.nii.gz')
            CSF = opj(BASE_priors_folder, BALSAname + '_label-CSF_probseg.nii.gz')
            GM = opj(BASE_priors_folder, BALSAname + '_label-GM_probseg.nii.gz')
            WM = opj(BASE_priors_folder, BALSAname + '_label-WM_probseg.nii.gz')
            Aseg_ref = opj(BASE_atlas_folder, BALSAname + '_seg-4FS_dseg.nii.gz')

    elif specie == 'Macaque':
        FS_refs          = opj(FS_tools, 'standard_mesh_atlases_macaque')
        BALSAname        = 'MY19'
        balsa_folder     = opj(path_BALSA, BALSAname,'modified')
        balsa_templateT1 = opj(balsa_folder, 'volumes', 'MY19_res-0.5mm_template_T1w.nii.gz')
        balsa_templateT2 = opj(balsa_folder, 'volumes', 'MY19_res-0.5mm_template_T2w.nii.gz')
        balsa_brainT1    = opj(balsa_folder, 'volumes', 'MY19_desc-SS_T1w.nii.gz')
        #ref_imc         = opj(balsa_folder, 'MY19_desc-crop_T1w.nii.gz')
        balsa_gray       = opj(balsa_folder, 'volumes', 'masks', 'MY19_desc-dilat-Gray_mask.nii.gz')
        balsa_label      = opj(balsa_folder, 'volumes', 'labels')

        if reference == '':
            reference = BALSAname
            BASE_folder = opj(balsa_folder, 'volumes')
            BASE_atlas_folder = opj(BASE_folder, 'labels')
            BASE_masks_folder = opj(BASE_folder, 'masks')
            BASE_priors_folder = opj(BASE_folder, 'priors')
            BASE_SS = balsa_brainT1
            BASE_mask = opj(BASE_masks_folder, BALSAname + '_mask.nii.gz')
            BASE_Gmask = opj(BASE_masks_folder, BALSAname + '_desc-Gray_mask.nii.gz')
            BASE_Wmask = opj(BASE_masks_folder, BALSAname + '_desc-erod-White_mask.nii.gz')
            BASE_Vmask = opj(BASE_masks_folder, BALSAname + '_desc-erod-Vent_mask.nii.gz')
            CSF = opj(BASE_priors_folder, BALSAname + '_label-CSF_probseg.nii.gz')
            GM = opj(BASE_priors_folder, BALSAname + '_label-GM_probseg.nii.gz')
            WM = opj(BASE_priors_folder, BALSAname + '_label-WM_probseg.nii.gz')
            Aseg_ref = opj(BASE_atlas_folder, BALSAname + '_seg-4FS_dseg.nii.gz')

    # for projection of volumetric results into template surfaces

    liste_fsaverage = ['MY19', 'CY29', '1200']

    if any(reference in word for word in liste_fsaverage):
        surfpath = opj(balsa_folder, 'surfaces')
        surfname = 'fsaverage_LR_32k'
        LRside = ['L', 'R']
        roi = 'atlasroi.32k_fs_LR'
    else:
        surfpath = opj(path_ref, reference, 'surfaces')
        surfname = 'Native_resol'
        LRside = ['l', 'r']
        roi = 'roi'
    surf_dir = opj(surfpath,surfname)

    # for PET studies
    # if you need a reference image
    # usually  the reference is the cerebellum (different irrigation,different types of receptors)
    REF = ''
    if not atlasname == '':
        name  = atlasname.split('-')[0]
        level = atlasname.split('-')[1]

        if name == 'BAL2013':
            if REFname == 'Cer':
                REF = [29, 30]
            Cortex_r    = [31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69]
            Cortex_l    = [32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70]
            Subcortex_r = [1, 3, 5, 7, 9,  11, 13, 15, 17, 21, 23, 29, 19]
            Subcortex_l = [2, 4, 6, 8, 10, 12, 14, 16, 18, 22, 24, 30, 20]

            listfig = [[Cortex_l, Cortex_r, Subcortex_l, Subcortex_r],
                       ['contro_Cort', 'ipsi_Cort', 'contro_subC', 'ipsi_subC']]

        elif name == 'BAL2020':
            if REFname == 'Cer':
                REF = [29, 30]
            Cortex_r = [31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69]
            Cortex_l = [32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70]
            Subcortex_r = [1, 3, 5, 7, 9, 11, 13, 15, 17, 21, 23, 29]
            Subcortex_l = [2, 4, 6, 8, 10, 12, 14, 16, 18, 22, 24, 30]
            thala_r = [95, 97, 99, 101, 103, 105, 107, 109, 19]
            thala_l = [96, 98, 100, 102, 104, 106, 108, 110, 20]

            listfig = [[Cortex_l, Cortex_r, Subcortex_l, Subcortex_r, thala_l, thala_r],
                       ['contro_Cort','ipsi_Cort','contro_subC','ipsi_subC','contro_thal','ipsi_thal']]

        elif name == 'BAL2023':
            if REFname == 'Cer':
                REF = [29, 30]

            Cortex_r    = [31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69]
            Cortex_l    = [32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70]
            Subcortex_r = [1, 3, 5, 7, 9,  11, 13, 15, 17, 21, 23, 29]
            Subcortex_l = [2, 4, 6, 8, 10, 12, 14, 16, 18, 22, 24, 30]
            thala_r     = [95, 97, 99,  101, 103, 105, 107, 109, 19]
            thala_l     = [96, 98, 100, 102, 104, 106, 108, 110, 20]

            listfig = [[Cortex_l, Cortex_r, Subcortex_l, Subcortex_r, thala_l, thala_r],
                       ['contro_Cort', 'ipsi_Cort', 'contro_subC', 'ipsi_subC', 'contro_thal', 'ipsi_thal']]

            atlasnii   = opj(BASE_atlas_folder, reference + '_seg-' + name + '_dseg.nii.gz')
            atlas_ctab = opj(BASE_folder, 'Rois', 'v_2023','fascicularis_cluster_info.txt')

        elif name == 'BAL2024':
            if REFname == 'Cer':
                REF = [29, 30]

            Cortex_r    = [31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69]
            Cortex_l    = [32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70]
            Subcortex_r = [1, 3, 5, 7, 9,  11, 13, 15, 113, 17, 21, 23, 29, 117]
            Subcortex_l = [2, 4, 6, 8, 10, 12, 14, 16, 114, 18, 22, 24, 30, 117]
            thala_r     = [19,95, 97, 99,  101, 103, 105, 107, 109]
            thala_l     = [20,96, 98, 100, 102, 104, 106, 108, 110]
            listfig = [[Cortex_l, Cortex_r, Subcortex_l, Subcortex_r, thala_l, thala_r],
                       ['contro_Cort', 'ipsi_Cort', 'contro_subC', 'ipsi_subC', 'contro_thal', 'ipsi_thal']]

            atlasnii   = opj(BASE_atlas_folder, reference + '_seg-' + name + '_desc-Gray_dseg.nii.gz')
            atlas_ctab = opj(BASE_folder, 'Rois', name + '_gray.csv')

        elif name == 'D99v2':
            if REFname == 'Cer':
                REF = [108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 125]
        elif name == 'CIVM':
            if REFname == 'Cer':
                print('This atlas is not adapted : there is no subdivision between gray and white matter')
                REF = [4]
        elif name == 'INIA19':
            if REFname == 'Cer':
                REF = [224, 225, 226, 227, 228, 229, 230, 231, 235, 236, 238, 240, 241, 242, 244, 245, 246, 247, 248, 250,
                   252, 253, 504, 505,
                   1224, 1225, 1226, 1227, 1228, 1229, 1230, 1231, 1235, 1236, 1238, 1240, 1241, 1242, 1244, 1245, 1246,
                   1247, 1248, 1250, 1252, 1253, 1504, 1505]
        elif name == 'ONPRC18':
            if REFname == 'Cer':
                REF = [57]
                print('This atlas is not adapted : there is no subdivision between gray and white matter')
        elif name == 'SARM':
            if REFname == 'Cer':
                if level == 1: REF = [203]
                elif level == 2: REF = [240]
                elif level == 3: REF = [244, 250]
                elif level == 4: REF = [244, 251, 267, 275]
                elif level == 5: REF = [244, 251, 267, 275, 246]
                elif level == 6: REF = [252, 253, 269, 270, 271, 272, 273, 274, 276, 277, ]
        elif atlasname == 'CHARM':
            if REFname == 'Cer':
                print('This atlas is not adapted : there is no cerebellum here')
                REF = ''
        elif atlasname == 'EDNIxCSCLR' or atlasname == 'EDNIxCSC':
            if REFname == 'Cer':
                REF=[7]

        else:
            REF = ''

    if kind == 'anat':
        return [FS_refs,path_ref,reference,
                balsa_folder,BALSAname,balsa_brainT1,
                BASE_folder,BASE_atlas_folder, BASE_template,BASE_SS, BASE_mask, BASE_Gmask, BASE_Wmask, BASE_Vmask,
                CSF, GM, WM, Aseg_ref,
                list_atlas, path_label_code]

    elif kind == 'func':
        return [BASE_folder, BASE_mask, BASE_Gmask, BASE_Wmask, BASE_Vmask]
    elif kind == 'pet':
        return [BASE_folder, BASE_mask,REF,atlasnii,atlas_ctab,surf_dir,LRside,roi]
    elif kind == 'RS_ana':
        return [BASE_SS, BASE_Gmask,atlasnii, BASE_mask,surf_dir,LRside,roi,atlas_ctab,listfig]


def atlas(list_atlas,add_atlas):
    if not add_atlas == '':
        list_atlas[0].append(add_atlas)
    if add_atlas == 'BAL2024':
        list_atlas[1].append('ctab')
        list_atlas[2] = list_atlas[2] + [1]
        list_atlas[3] = list_atlas[3] + [1]

    return list_atlas




