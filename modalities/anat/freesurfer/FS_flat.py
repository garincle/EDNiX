#import
import os
opj  = os.path.join

from Tools import run_cmd

def map(FS_dir,Subname,Side,cmd_tksurfer,cmd_flatten,diary_name,sing_fs,export_fs):

    '''
    - draw the circle around the carpus callosum and the midbrain (left click and Cut closed line)
    - draw the line to cut down the calcarine and so on (frontal pole, above and below insula, temporal pole and the POS), finish with the cut line
        ==>  left click outside the midline circle and press fill uncut area.
    - save the patch as ?h.full.patch.3d
    - Once this is done, close tksurfer and press quit
    - better Check. You may have to repeat the cut proceeding according to the distortions you gonna get.
    '''

    if Side == 'left':
        flat = 'lh'
    elif Side == 'right':
        flat = 'rh'
    else:
        print('Sorry, try again')

    nl = 'Create the ' + Side +  ' flat surface'
    run_cmd.msg(nl, diary_name,  'OKGREEN')

    cmd = export_fs + cmd_tksurfer + Subname + ' ' + flat + ' inflated -curv -gray'
    run_cmd.do(cmd,diary_name)

    cmd = sing_fs + cmd_flatten + '-w 50 ' + opj(FS_dir, Subname, 'surf', flat + '.full.patch.3d') + ' ' + \
    opj(FS_dir, Subname, 'surf',flat + '.full.flatten.patch.3d')
    run_cmd.run(cmd,diary_name)

    cmd = export_fs + cmd_tksurfer + Subname + ' ' + flat + ' inflated -curv -gray -patch ' + flat + '.full.flatten.patch.3d'
    run_cmd.do(cmd,diary_name)

    m_dir = os.getcwd()
    #os.remove(opj(m_dir,flat + '.full.flatten.patch.3d.out'))
