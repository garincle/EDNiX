import os

opj = os.path.join
opd = os.path.dirname

from Tools import run_cmd


AFNI_sif     = 'afni_ub24_latest.sif'
FSL_sif      = 'fsl_6.0.5.1-cuda9.1.sif'
FS_sif       = 'freesurfer_NHP.sif'
WB_sif       = 'connectome_workbench_1.5.0-freesurfer-update.sif'
ITK_sif      = 'itksnap_5.0.9.sif'
turku_sif    = 'tpcclib.sif'
synStrip_sif = 'synthstrip.1.5.sif'

def load_requirement(MAIN_PATH,reftemplate_path,bids_dir,singularity):

    if singularity == 'no':
        sing_afni     = ''
        sing_fsl      = ''
        sing_fs       = ''
        sing_wb       = ''
        sing_itk      = ''
        sing_turku    = ''
        sing_synStrip = ''
        Unetpath = opj(MAIN_PATH,'Tool_library/Singularity/NHP-BrainExtraction/UNet_Model/models/')

    elif singularity == 'yes':
        s_path = opj(opd(MAIN_PATH), 'Tool_library', 'Singularity')
        s_bind = '--bind ' + ','.join([bids_dir,reftemplate_path,s_path])

        sing_func1  = 'singularity run ' + s_bind
        sing_func2  = 'singularity exec ' + s_bind

        sing_afni      = sing_func1 + ' ' + opj(s_path, AFNI_sif)     + ' '
        sing_fsl       = sing_func1 + ' ' + opj(s_path, FSL_sif)      + ' '
        sing_fs        = sing_func1 + ' ' + opj(s_path, FS_sif)       + ' '
        sing_itk       = sing_func1 + ' ' + opj(s_path, ITK_sif)      + ' '
        sing_wb        = sing_func1 + ' ' + opj(s_path, WB_sif)       + ' '
        sing_turku     = sing_func2 + ' ' + opj(s_path, turku_sif)    + ' '
        sing_synStrip  = sing_func1 + ' ' + opj(s_path, synStrip_sif) + ' '

        Unetpath = s_path

    return sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, sing_turku,sing_synStrip,Unetpath


def FS(fs_tools,FS_dir,diary_name,sing_fs):

    cmd = sing_fs + 'mris_convert --version'
    FS_version = run_cmd.get(cmd,diary_name)

    FS_v = FS_version[0].decode('utf-8').split(' ')
    if FS_v[2] == 'stable6':
        cmd_tksurfer = 'tksurfer '
        cmd_flatten  = 'mris_flatten '
        cmd_mris     = 'mris_convert'
    else:
        if sing_fs == '':
            cmd_tksurfer = 'tksurfer_v6 '
            cmd_flatten  = 'mris_flatten_v6 '
            cmd_mris     = 'mris_convert_v6'
        else:
            cmd_tksurfer = opj(fs_tools, 'tksurfer_v6 ')
            cmd_flatten  = opj(fs_tools, 'mris_flatten_v6 ')
            cmd_mris     = opj(fs_tools, 'mris_convert_v6')

    os.environ['SUBJECTS_DIR'] = FS_dir
    cmd = 'echo $SUBJECTS_DIR'
    run_cmd.msg(cmd,diary_name,'ENDC')

    if sing_fs != '':
        export_fs = 'export SINGULARITYENV_SUBJECTS_DIR="' + FS_dir + '";' + sing_fs
    else:
        export_fs = ''


    return cmd_tksurfer,cmd_flatten,cmd_mris,export_fs


