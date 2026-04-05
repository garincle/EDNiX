import os
opn = os.path.normpath
opj = os.path.join
opb = os.path.basename
opd = os.path.dirname

from Tools import Load_subject_with_BIDS
from atlases import templatefeat

def get(MAIN_PATH,bids_dir,allinfo_study_c, species,list_to_keep,list_to_remove,reference,type_norm,MNIBcorrect_indiv, masking_img, atlas_followers):

    # reference templates
    if reference == '':
        reference = 'EDNiX'

    reftemplate_path = opj(opd(MAIN_PATH), "Atlases_library")
    fs_tools = opj(reftemplate_path, 'freesurfer_refs')
    path_BALSA = opj(reftemplate_path, 'BALSA')

    (FS_refs, template_dir, reference,
     balsa_folder, BALSAname, balsa_brainT1,
     BASE_folder,BASE_atlas_folder, BASE_template, BASE_SS, BASE_mask, BASE_Gmask, BASE_WBGmask, BASE_Wmask, BASE_Vmask,
     CSF, GM, WM, Aseg_ref,
     list_atlas, path_label_code) = templatefeat.get(species, reftemplate_path, fs_tools, path_BALSA, reference,
                                                     '', '',type_norm,'anat', atlas_followers)


    all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max = Load_subject_with_BIDS.load_data_bids(
        allinfo_study_c, bids_dir, list_to_keep, list_to_remove)

    # by default values

    if MNIBcorrect_indiv == '':
        MNIBcorrect_indiv = 'N4'  # "N4" or "N3"

    if masking_img == '':
        masking_img = type_norm # "N4" or "N3"

    return (FS_refs, template_dir, reference,
     balsa_folder, BALSAname, balsa_brainT1,
     BASE_atlas_folder, BASE_template, BASE_SS, BASE_mask, BASE_Gmask, BASE_WBGmask, BASE_Wmask, BASE_Vmask,
     CSF, GM, WM, Aseg_ref,
     list_atlas, path_label_code,all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max,
            fs_tools,reftemplate_path,MNIBcorrect_indiv, masking_img)