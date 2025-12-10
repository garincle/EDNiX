import os
import glob

opj = os.path.join
opb = os.path.basename
opr = os.path.relpath
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import getpath
from Tools import run_cmd

def retrieve(ID,data_path, Session,anat_func_same_space,use_erode_V_func_masks,use_erode_WM_func_masks,TfMRI,diary_file):

    if anat_func_same_space == True:

        # The anatomy organization folders

        path_anat, dir_transfo, FS_dir, dir_prepro, dir_native, volumes_dir, labels_dir, masks_dir = getpath.anat(
            data_path, '', '', '', '', 'native')

    else:
        _, _, _, _, _, volumesdir, _, _ = getpath.anat(data_path, '', '', '', '', 'native')
        branch_path = opr(volumesdir, start=data_path)

        template_anat_for_fmri = glob.glob(opj(opd(data_path),
                                               'ses-' + str(Session), branch_path, '*_desc-SS_' + TfMRI + '.nii*'))
        if len(template_anat_for_fmri) > 1:
            nl = "WARNING: we found multiple anat template for this animal, THIS IS NOT NORMAL, please check what happened"
            run_cmd.msg(nl, diary_file, 'WARNING')

        elif len(template_anat_for_fmri) == 0:
            nl = "WARNING: we haven't found an anat template for this animal, in this SESSION, please check that this is what you want!! also check that do_anat_to_func = True  anat_func_same_space = False"
            run_cmd.msg(nl, diary_file, 'WARNING')

            nl = 'INFO: Lets scanned the whole subject folder'
            run_cmd.msg(nl, diary_file, 'OKGREEN')

            template_anat_for_fmri = glob.glob(opj(opd(data_path), '**', branch_path, '*_desc-SS_' + TfMRI + '.nii*'))

            if not template_anat_for_fmri ==[]:

                nb = len(opd(data_path).split('/'))+1
                new_sess = template_anat_for_fmri[0].split('/')[nb]

                if len(template_anat_for_fmri) == 1:
                    nl = ('WARNING: we found one anat template for this animal in this folder : ' + new_sess +
                          ' , let\'s continue with it')

                elif len(template_anat_for_fmri) > 1:
                    nl = "WARNING: we found multiple anat template for this animal, we will choose the first one by DEFAULT!"
                run_cmd.msg(nl, diary_file, 'WARNING')

                path_anat, dir_transfo, FS_dir, dir_prepro, dir_native, volumes_dir, labels_dir, masks_dir = getpath.anat(
                    opj(opd(data_path), new_sess), '', '', '', '','native')


            else:
                nl = ("ERROR: We haven't found any anat template for this animal! We can't continue ! please provide a valid link for at least one anat image! "
                         "current link is :") + str(template_anat_for_fmri)
                raise Exception(run_cmd.error(nl, diary_file))
        else:
            path_anat, dir_transfo, FS_dir, dir_prepro, dir_native, volumes_dir, labels_dir, masks_dir = getpath.anat(
                opj(opd(data_path), 'ses-' + str(Session)), '', '', '', '',
                'native')

            nl = 'INFO: We found this image as template: ' + str(template_anat_for_fmri)
            run_cmd.msg(nl, diary_file, 'OKGREEN')

    anat_subject = opj(volumes_dir, ID + '_space-acpc_desc-template_' + TfMRI + '.nii.gz')
    brainmask    = opj(masks_dir,   ID + '_space-acpc_mask.nii.gz')
    G_mask       = opj(masks_dir,   ID + '_desc-Gray_mask.nii.gz')
    if use_erode_V_func_masks == True:
        V_mask = opj(masks_dir, ID + '_desc-erod-Vent_mask.nii.gz')
    else:
        V_mask = opj(masks_dir, ID + '_desc-Vent_mask.nii.gz')
    if use_erode_WM_func_masks == True:
        W_mask = opj(masks_dir, ID + '_desc-erod-White_mask.nii.gz')
    else:
        W_mask = opj(masks_dir, ID + '_desc-White_mask.nii.gz')

    return anat_subject,brainmask,G_mask,V_mask,W_mask,dir_transfo, FS_dir, dir_prepro, volumes_dir, labels_dir, masks_dir


def create(folderforTemplate_Anat,diary_file):

    # The anatomy organization folders
    dir_transfo = ''
    dir_prepro  = ''
    volumes_dir = ''
    labels_dir  = ''
    masks_dir   = ''
    FS_dir      =''

    anat_subject = opj(folderforTemplate_Anat, 'template.nii.gz')
    brainmask = opj(folderforTemplate_Anat,    'brainmask.nii.gz')
    V_mask = opj(folderforTemplate_Anat, 'Vmask.nii.gz')
    W_mask = opj(folderforTemplate_Anat, 'Wmask.nii.gz')
    G_mask = opj(folderforTemplate_Anat, 'Gmask.nii.gz')

    nl = ' No individual anat template will be used for this animal'
    run_cmd.msg(nl, diary_file, 'OKGREEN')

    return anat_subject, brainmask, G_mask, V_mask, W_mask, dir_transfo, FS_dir, dir_prepro, volumes_dir, labels_dir, masks_dir