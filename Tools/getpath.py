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
    dir_rs_prepro = opj(path_func, 'preprocessing')
    dir_prepro_orig     = opj(dir_rs_prepro, 'des-space-orig')
    dir_prepro_acpc     = opj(dir_rs_prepro, 'des-space-acpc')
    dir_prepro_template = opj(dir_rs_prepro, 'des-space-template',reference)
    dir_residuals          = opj(path_func, 'residuals')
    dir_residuals_acpc     = opj(dir_residuals, 'acpc')
    dir_residuals_template = opj(dir_residuals, reference)

    return (path_func,dir_fmap,dir_rs_prepro,dir_prepro_orig,dir_prepro_acpc,dir_prepro_template,
            dir_residuals,dir_residuals_acpc,dir_residuals_template)


def pet(data_path, reference, BALSAname):
    print('blabla')


