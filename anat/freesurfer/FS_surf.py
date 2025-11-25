#import
import os
import shutil
import ants
import glob
import json

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
opi = os.path.isfile

from Tools import run_cmd, get_orientation
from Tools import getpath

from anat.freesurfer import mgz2ants
from anat.freesurfer import preFS
from anat.freesurfer import smallbrain

surface = ['Left','Right']
Hmin    = ['l','r']
TIF     = 'h.average.curvature.filled.buckner40.tif'
GCS     = 'h.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs'


def Wcreate(FS_dir, animal,LR,diary_name,sing_fs,export_fs):

    nl = 'Create the white surface'
    run_cmd.msg(nl, diary_name,'HEADER')

    if LR =='l':
        h = [0]
    elif LR== 'r':
        h = [1]
    elif LR =='lr':
        h=range(2)

    filled = mgz2ants.read(opj(FS_dir, animal, 'mri', 'filled.mgz'),diary_name,sing_fs)
    aseg   = mgz2ants.read(opj(FS_dir, animal, 'mri', 'aseg.mgz'),diary_name,sing_fs)

    value = [int(filled[aseg == 10].min()),
             int(filled[aseg == 49].max())]

    for H in h:
        # print(value[H])
        # 1) Create a coarse white surface. It takes only into account the filled.mgz image
        cmd = sing_fs + 'mri_pretess ' + opj(FS_dir,animal,'mri','filled.mgz') + ' ' + str(value[H]) + ' ' + opj(FS_dir,animal,'mri','brain.mgz') + ' ' + \
        opj(FS_dir,animal,'mri','filled' + str(value[H]) + '.mgz')
        run_cmd.run(cmd,diary_name)
        cmd = (sing_fs + 'mri_tessellate '              + opj(FS_dir,animal,'mri', 'filled' + str(value[H]) + '.mgz') +
               ' ' + str(value[H]) + ' ' + opj(FS_dir, animal, 'surf', Hmin[H] + 'h.orig.nofix'))
        run_cmd.run(cmd,diary_name)
        cmd = sing_fs + 'mris_extract_main_component ' + opj(FS_dir,animal,'surf', Hmin[H] + 'h.orig.nofix') + ' ' + opj(FS_dir, animal, 'surf', Hmin[H] + 'h.orig.nofix')
        run_cmd.run(cmd,diary_name)
        cmd = sing_fs + 'mris_smooth -nw '             + opj(FS_dir,animal,'surf', Hmin[H] + 'h.orig.nofix') + ' ' + opj(FS_dir, animal, 'surf', Hmin[H] + 'h.smoothwm.nofix')
        run_cmd.run(cmd,diary_name)
        cmd = sing_fs + 'mris_inflate -no-save-sulc '  + opj(FS_dir, animal, 'surf', Hmin[H] + 'h.smoothwm.nofix') + ' ' + \
        opj(FS_dir,animal,'surf',Hmin[H] + 'h.inflated.nofix')
        run_cmd.run(cmd,diary_name)
        cmd = sing_fs + 'mris_sphere -q '             + opj(FS_dir,animal,'surf',Hmin[H] + 'h.inflated.nofix') + ' ' + opj(FS_dir,animal,'surf',Hmin[H] + 'h.qsphere.nofix')
        run_cmd.run(cmd,diary_name)
        
        # 2) Automatic Correction
        shutil.copyfile(opj(FS_dir,animal,'surf', Hmin[H] + 'h.orig.nofix'),    opj(FS_dir,animal,'surf', Hmin[H] + 'h.orig'))
        shutil.copyfile(opj(FS_dir,animal,'surf', Hmin[H] + 'h.inflated.nofix'),opj(FS_dir,animal,'surf', Hmin[H] + 'h.inflated'))
    
        cmd = export_fs + 'mris_fix_topology -mgz -sphere qsphere.nofix -ga ' + animal + ' ' + Hmin[H] + 'h'
        run_cmd.do(cmd,diary_name)
        cmd = sing_fs + 'mris_euler_number ' + opj(FS_dir,animal,'surf',Hmin[H] + 'h.orig')
        run_cmd.run(cmd,diary_name)
        cmd = sing_fs + 'mris_remove_intersection ' + opj(FS_dir,animal,'surf', Hmin[H] + 'h.orig') + ' ' + opj(FS_dir,animal,'surf', Hmin[H] + 'h.orig')
        run_cmd.run(cmd,diary_name)
        
        # 3) Fine white surface;  The coarse white surface is modified according to the 3 images: "wm", "aseg" and "brain".
        cmd = export_fs + 'mris_make_surfaces -aseg aseg -whiteonly -noaparc -mgz -T1 brain ' + animal + ' ' + Hmin[H] + 'h'
        run_cmd.do(cmd,diary_name)
        
        os.remove(opj(FS_dir,animal,'surf', Hmin[H] + 'h.inflated'))
        os.remove(opj(FS_dir,animal,'mri', 'filled' + str(value[H]) + '.mgz'))
        
        nl = surface[H] + ' surface done!'
        run_cmd.msg(nl, diary_name,'OKGREEN')

        
def Wconfig(FS_dir, FS_refs,animal,LR,diary_name,sing_fs,export_fs):

    nl = 'Create the sphere surface,  the spheretosphere and the curvature files '
    run_cmd.msg(nl, diary_name,'HEADER')

    if LR =='l':
        h = [0]
    elif LR== 'r':
        h = [1]
    elif LR =='lr':
        h=range(2)


    for H in h:

        # just in case.......
        #if opi(opj(FS_dir,animal,'surf',Hmin[H] + 'h.white_backup')):
        #    os.remove(opj(FS_dir,animal,'surf',Hmin[H] + 'h.white_backup'))
        shutil.copyfile(opj(FS_dir,animal,'surf',Hmin[H] + 'h.white'),opj(FS_dir,animal,'surf',Hmin[H] + 'h.white_backup'))
        
        nl = '# 1) Inflation and curvature left surface\n'
        run_cmd.msg(nl, diary_name,'OKGREEN')
        cmd = sing_fs + 'mris_smooth -n 3 -nw ' + opj(FS_dir,animal,'surf', Hmin[H] + 'h.white') + ' ' + opj(FS_dir,animal,'surf', Hmin[H] + 'h.smoothwm')
        run_cmd.run(cmd,diary_name)
        cmd = sing_fs + 'mris_inflate ' + opj(FS_dir,animal,'surf',Hmin[H] + 'h.smoothwm') + ' ' + opj(FS_dir,animal,'surf', Hmin[H] + 'h.inflated')
        run_cmd.run(cmd,diary_name)
        
        nl = '# 2) Curvature measurements (essential for the pial surface)\n'
        run_cmd.msg(nl, diary_name,'OKGREEN')

        cmd = sing_fs + 'mris_curvature -w ' + opj(FS_dir, animal, 'surf', Hmin[H] + 'h.white')
        run_cmd.run(cmd,diary_name)
        
        cmd = sing_fs + 'mris_curvature -thresh .999 -n -a 5 -w -distances 10 5 ' + opj(FS_dir,animal,'surf', Hmin[H] + 'h.inflated')
        run_cmd.run(cmd,diary_name)
        
        #      means thresholding curvature at 99.90% level
        #      normalizing curvature values.
        #      averaging curvature patterns 5 times.
        #      sampling 10 neighbors out to a distance of 10mm (10mm true for human what about animals data  ? I chose 5mm from the size factor)
        #                  lh.orig.H = mean curvature
        #	               lh.orig.K = gaussian curvature
        cmd = (export_fs + 'mris_curvature_stats -m --writeCurvatureFiles -G -o ' + opj(FS_dir,animal,'stats',Hmin[H] + 'h.curv.stats') +
               ' -F smoothwm ' + animal + ' ' + Hmin[H] + 'h curv sulc')
        run_cmd.do(cmd,diary_name)
        
        nl = '# 3) Sphere registration.\n'
        run_cmd.msg(nl, diary_name,'OKGREEN')

        #   3.1) rather time consuming
        cmd = sing_fs + 'mris_sphere ' + opj(FS_dir,animal,'surf',Hmin[H] + 'h.inflated') + ' ' + opj(FS_dir,animal,'surf', Hmin[H] + 'h.sphere')
        run_cmd.run(cmd,diary_name)

        # 3.2) VERY time consuming
        DIR = os.getcwd()
        os.chdir(opj(FS_dir,animal,'surf'))
        cmd = (export_fs + 'mris_register -curv -dist 5 -max_degrees 68 ' + Hmin[H] + 'h.sphere ' +
               opj(FS_refs, Hmin[H] + TIF) + ' ' + Hmin[H] + 'h.sphere.reg')
        run_cmd.do(cmd,diary_name)
        os.chdir(DIR)
        
        # 3.3) store the deformations
        cmd = (sing_fs + 'mris_jacobian ' + opj(FS_dir,animal,'surf',Hmin[H] + 'h.white') +
               ' ' + opj(FS_dir,animal,'surf', Hmin[H] + 'h.sphere.reg') +
               ' ' + opj(FS_dir,animal,'surf', Hmin[H] + 'h.jacobian_white'))
        run_cmd.run(cmd,diary_name)
        
        # 3.4 define what is the cortex: if thi step is not done the T2 correction 
        # (beware it has the tendency to overestimate the cortical thickness and to cross the middline)
        
        cmd = (export_fs + 'mris_ca_label -l ' + opj(FS_dir,animal,'label',Hmin[H] + 'h.cortex.label') +
               ' -aseg ' + opj(FS_dir,animal,'mri','aseg.mgz') +
               ' ' + animal + ' ' + Hmin[H] + 'h' +
               ' ' + opj(FS_dir,animal,'surf',Hmin[H] + 'h.sphere.reg') +
               ' ' + opj(FS_refs, Hmin[H] + GCS) +
               ' ' + Hmin[H] + 'h.aparc.annot')
        run_cmd.do(cmd,diary_name)

        nl = surface[H] + ' surface done!'
        run_cmd.msg(nl, diary_name,'OKGREEN')

def Pcreate(animal,LR,diary_name,export_fs):

    nl = 'Create the pial surface'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    if LR =='l':
        h = [0]
    elif LR== 'r':
        h = [1]
    elif LR =='lr':
        h=range(2)

    for H in h:

        cmd = (export_fs + 'mris_make_surfaces -white NOWRITE -aseg aseg -orig white -noaparc -mgz -T1 brain ' +
               animal + ' ' + Hmin[H] + 'h')
        run_cmd.do(cmd,diary_name)

        nl = surface[H] + ' surface done!'
        run_cmd.msg(nl, diary_name,'OKGREEN')

def PcreateT2(T2_img,dir_path,aseg,LR,change_hd,diary_name,sing_fs,export_fs):

    nl = 'Prepare for the T2w image and modify the pial surface accordingly'
    run_cmd.msg(nl, diary_name,'HEADER')

    if LR =='l':
        h = [0]
    elif LR== 'r':
        h = [1]
    elif LR =='lr':
        h=range(2)

    (_, _, FS_dir, dir_prepro, _, _, labels_dir, _,) = getpath.anat(dir_path, 'reference', '',
                                                                                         False,False,'native')

    animal  = opb(T2_img).split('_')[0]

    _, _, fwdFS_cmd, _, _ = get_orientation.use_ants(T2_img, sing_fs)

    preFS.prepa_T2img(animal, T2_img, dir_prepro, aseg, diary_name)

    new_T2_img = opj(dir_prepro, animal + '_desc-norm_T2w.nii.gz')

    if change_hd == 1:

        _, spacing, _, new_size, _ = smallbrain.check(opj(dir_prepro, animal + '_rescale.json'))

        ref  = opj(dir_prepro, animal + '_desc-norm_T1w_resamp-' + str(spacing[0]) +'.nii.gz')

        brain_img = ants.image_read(ref)
        new_orig  = brain_img.origin

        T2_file = opj(dir_prepro, animal + '_desc-norm_T2w_resamp-' + str(spacing) +'.nii.gz')
        resamp_img = ants.image_read(new_T2_img)
        ants.set_spacing(resamp_img, [spacing,spacing,spacing])
        ants.set_origin(resamp_img, [new_orig[0], new_orig[1], new_orig[2]])
        resamp_img = ants.resample_image(resamp_img, (new_size, new_size, new_size), False, 1)
        resamp_img[resamp_img < 0] = 0
        ants.image_write(resamp_img, T2_file)

        dictionary = {"Sources": [new_T2_img,
                                  ref],
                      "Description": 'reformat image with a voxel size of ' + str(
                          new_size) + 'and a resampling of ' + str(spacing) + ' (Antspy)'}
        json_object = json.dumps(dictionary, indent=2)
        with open(T2_file.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)

    else :
        T2_file = new_T2_img

    cmd = sing_fs + 'mri_convert ' + fwdFS_cmd + T2_file + ' ' + opj(FS_dir, animal, 'mri', 'T2.mgz')
    run_cmd.run(cmd,diary_name)

    for H in h:
        #if opi(opj(FS_dir,animal,'surf',Hmin[H] + 'h.pial.T1')):
        #    os.remove(opj(FS_dir,animal,'surf',Hmin[H] + 'h.pial.T1'))
        shutil.copyfile(opj(FS_dir,animal,'surf',Hmin[H] + 'h.pial'),opj(FS_dir,animal,'surf',Hmin[H] + 'h.pial.T1'))

        cmd = (export_fs + 'mris_make_surfaces -aseg aseg -mgz -orig white -nowhite -orig_white white -orig_pial pial -T2dura '
               +  opj(FS_dir, animal, 'mri', 'T2') + ' -T1 brain ' + animal + ' ' + Hmin[H] + 'h')
        run_cmd.do(cmd,diary_name)
        nl = surface[H] + ' surface done!'
        run_cmd.msg(nl, diary_name,'OKGREEN')

