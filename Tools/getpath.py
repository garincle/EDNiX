#import
import os

opj = os.path.join

liste_balsa = ['S1200', 'MY19', 'CY29']

def rawdata(data_path):
    raw_dir = opj(data_path, 'anat')
    return raw_dir


def anat(data_path,reference,BALSAname,creat_study_template,coregistration_longitudinal,param):

    path_anat   = opj(data_path, 'anat')

    dir_transfo = opj(path_anat, 'matrices')
    FS_dir      = opj(path_anat, 'freesurfer')
    dir_prepro  = opj(path_anat, 'preprocessing')
    dir_native  = opj(path_anat, 'native')

    volumes_dir = opj(dir_native, 'volumes')
    labels_dir  = opj(volumes_dir, 'labels')
    masks_dir   = opj(volumes_dir, 'masks')

    wb_balsa_dir    = ''
    wb_balsa_vol    = ''
    wb_balsa_labels = ''
    wb_balsa_masks  = ''

    wb_studytemplate_dir    = ''
    wb_studytemplate_vol    = ''
    wb_studytemplate_labels = ''
    wb_studytemplate_masks  = ''

    wb_refsession_dir    = ''
    wb_refsession_vol    = ''
    wb_refsession_labels = ''
    wb_refsession_masks  = ''

    if creat_study_template == True:
        wb_studytemplate_dir    = opj(path_anat, 'templates', 'studyTemplate')
        wb_studytemplate_vol    = opj(wb_studytemplate_dir, 'volumes')
        wb_studytemplate_labels = opj(wb_studytemplate_vol, 'labels')
        wb_studytemplate_masks  = opj(wb_studytemplate_vol, 'masks')

    if coregistration_longitudinal == True:
        wb_refsession_dir    = opj(path_anat, 'templates', 'refSession')
        wb_refsession_vol    = opj(wb_refsession_dir, 'volumes')
        wb_refsession_labels = opj(wb_refsession_vol, 'labels')
        wb_refsession_masks  = opj(wb_refsession_vol, 'masks')

    if any(BALSAname in word for word in liste_balsa):
        wb_balsa_dir    = opj(path_anat, 'templates', BALSAname)
        wb_balsa_vol    = opj(wb_balsa_dir, 'volumes')
        wb_balsa_labels = opj(wb_balsa_vol, 'labels')
        wb_balsa_masks  = opj(wb_balsa_vol, 'masks')

    if not reference == BALSAname:
        wb_template_dir    = opj(path_anat, 'templates', reference)
        wb_template_vol    = opj(wb_template_dir, 'volumes')
        wb_template_labels = opj(wb_template_vol, 'labels')
        wb_template_masks  = opj(wb_template_vol, 'masks')
    else:
        wb_template_dir    = wb_balsa_dir
        wb_template_vol    = wb_balsa_vol
        wb_template_labels = wb_balsa_labels
        wb_template_masks  = wb_balsa_masks

    if param == 'all':
        return (path_anat,dir_transfo,FS_dir,dir_prepro,dir_native,volumes_dir,labels_dir,masks_dir,
            wb_template_dir,wb_template_vol,wb_template_labels,wb_template_masks,
            wb_balsa_dir,wb_balsa_vol,wb_balsa_labels,wb_balsa_masks,
            wb_studytemplate_dir,wb_studytemplate_vol,wb_studytemplate_labels,wb_studytemplate_masks,
            wb_refsession_dir,wb_refsession_vol,wb_refsession_labels,wb_refsession_masks)

    elif param == 'native':
        return (path_anat, dir_transfo, FS_dir, dir_prepro, dir_native, volumes_dir, labels_dir, masks_dir)

    elif param == 'template':
        return (path_anat, dir_transfo, FS_dir, dir_prepro, dir_native, volumes_dir, labels_dir, masks_dir,
                wb_template_dir, wb_template_vol, wb_template_labels, wb_template_masks,
                wb_balsa_dir, wb_balsa_vol, wb_balsa_labels, wb_balsa_masks)

    elif param == 'studyTemplate':
        return (path_anat, dir_transfo, FS_dir, dir_prepro, dir_native, volumes_dir, labels_dir, masks_dir,
                wb_studytemplate_dir,wb_studytemplate_vol,wb_studytemplate_labels,wb_studytemplate_masks)

    elif param == 'refSession':
        return (path_anat, dir_transfo, FS_dir, dir_prepro, dir_native, volumes_dir, labels_dir, masks_dir,
                wb_refsession_dir, wb_refsession_vol, wb_refsession_labels, wb_refsession_masks)


def surf(path_anat,reference,BALSAname):
    dir_native_resol = opj(path_anat, 'native', 'surfaces', 'Native_resol')
    if BALSAname == '':
        dir_native_32    = ''
        dir_balsa_resol  = ''
        dir_balsa_32     = ''
        dir_balsa_64     = ''

    else :
        dir_native_32   = opj(path_anat, 'native', 'surfaces', 'fsaverage_LR_32k')
        dir_balsa_resol = opj(path_anat, 'templates', BALSAname, 'surfaces', 'Native_resol')
        dir_balsa_32    = opj(path_anat, 'templates', BALSAname, 'surfaces', 'fsaverage_LR_32k')
        dir_balsa_64    = opj(path_anat, 'templates', BALSAname, 'surfaces', 'fsaverage_164k')

    return dir_native_resol,dir_native_32,dir_balsa_resol,dir_balsa_32,dir_balsa_64

def stytemplate(data_path,reference,BALSAname):

    dir_acpc    = opj(data_path, 'derivatives','acpc')
    dir_prepro  = opj(data_path, 'derivatives','tmp')
    dir_transfo = opj(dir_acpc, 'matrices')
    volumes_dir = opj(dir_acpc, 'volumes')
    labels_dir  = opj(volumes_dir, 'labels')
    masks_dir   = opj(volumes_dir, 'masks')

    wb_balsa_dir    = ''
    wb_balsa_vol    = ''
    wb_balsa_labels = ''
    wb_balsa_masks  = ''

    if not reference == BALSAname:
        wb_template_dir    = opj(data_path, 'templates', reference)
        wb_template_vol    = opj(wb_template_dir, 'volumes')
        wb_template_labels = opj(wb_template_vol, 'labels')
        wb_template_masks  = opj(wb_template_vol, 'masks')
    else:
        wb_template_dir    = wb_balsa_dir
        wb_template_vol    = wb_balsa_vol
        wb_template_labels = wb_balsa_labels
        wb_template_masks  = wb_balsa_masks

    return (dir_acpc,dir_prepro,dir_transfo, volumes_dir, labels_dir, masks_dir,
            wb_template_dir,wb_template_vol,wb_template_labels,wb_template_masks,
            wb_balsa_dir,wb_balsa_vol,wb_balsa_labels,wb_balsa_masks)

def func(data_path,reference):
    #path_rawfunc = opj(rawdata_path, 'func')
    path_func    = opj(data_path, 'func')
    #dir_rawfmap     = opj(rawdata_path, 'fmap')
    dir_fmap = opj(data_path, 'fmap')
    dir_prepro_fmap = opj(dir_fmap, 'preprocessing')

    dir_prepro_orig     = opj(path_func, 'orig-func')
    dir_prepro_orig_labels = opj(dir_prepro_orig, 'labels')
    dir_prepro_orig_masks = opj(dir_prepro_orig, 'masks')
    dir_prepro_orig_process = opj(dir_prepro_orig, 'preprocessing')
    dir_prepro_orig_rs = opj(dir_prepro_orig, 'postprocessed_rs')
    dir_prepro_orig_task = opj(dir_prepro_orig, 'postprocessed_task')

    dir_prepro_acpc     = opj(path_func, 'acpc-anat')
    dir_prepro_acpc_labels = opj(dir_prepro_acpc, 'labels')
    dir_prepro_acpc_masks = opj(dir_prepro_acpc, 'masks')
    dir_prepro_acpc_process = opj(dir_prepro_orig, 'preprocessing')
    dir_prepro_acpc_rs = opj(dir_prepro_acpc, 'postprocessed_rs')
    dir_prepro_acpc_task = opj(dir_prepro_acpc, 'postprocessed_task')

    dir_prepro_template = opj(path_func, 'templates',reference)
    dir_prepro_template_labels = opj(dir_prepro_orig, 'labels')
    dir_prepro_template_masks = opj(dir_prepro_orig, 'masks')
    dir_prepro_template_process = opj(dir_prepro_orig, 'preprocessing')
    dir_prepro_template_rs = opj(dir_prepro_orig, 'postprocessed_rs')
    dir_prepro_template_task = opj(dir_prepro_orig, 'postprocessed_task')

    return (path_func,dir_fmap, dir_prepro_fmap, dir_prepro_orig, dir_prepro_orig_labels, dir_prepro_orig_masks,
            dir_prepro_orig_process, dir_prepro_orig_rs, dir_prepro_orig_task,
            dir_prepro_acpc, dir_prepro_acpc_labels, dir_prepro_acpc_masks,
            dir_prepro_acpc_process, dir_prepro_acpc_rs, dir_prepro_acpc_task,
            dir_prepro_template, dir_prepro_template_labels, dir_prepro_template_masks,
            dir_prepro_template_process, dir_prepro_template_rs, dir_prepro_template_task)


def pet(data_path, reference, BALSAname):
    print('blabla')


