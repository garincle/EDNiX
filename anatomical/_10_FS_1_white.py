#import
import os
import subprocess
import shutil

#Path to the excels files and data structure
opj = os.path.join
ope = os.path.exists
spco = subprocess.check_output

surface = ['Left','Right']
Hmin    = ['l','r']

def White_create(FS_dir, monkey):
    
    value   = [255, 127]
    
    # white matter volumes with weigth accorgind to the subcortical structures
    command = 'mri_edit_wm_with_aseg ' + opj(FS_dir,monkey,'mri','wm.seg.mgz') + ' ' + opj(FS_dir,monkey,'mri','brain.mgz') + \
    ' ' + opj(FS_dir,monkey,'mri','aseg.mgz')  + ' '  + opj(FS_dir,monkey,'mri','wm.asegedit.mgz')
    spco([command], shell=True)
    
    command = 'mri_pretess ' + opj(FS_dir,monkey,'mri','wm.asegedit.mgz') + ' wm ' + opj(FS_dir,monkey,'mri','brain.mgz') + \
    ' ' + opj(FS_dir,monkey,'mri','wm.mgz')
    spco([command], shell=True)
    
    for H in range(0,2):
        # 1) Create a coarse white surface. It takes only into account the filled.mgz image
        command = 'mri_pretess ' + opj(FS_dir,monkey,'mri','filled.mgz') + ' ' + str(value[H]) + ' ' + opj(FS_dir,monkey,'mri','brain.mgz') + ' ' + \
        opj(FS_dir,monkey,'mri','filled' + str(value[H]) + '.mgz')
        spco([command], shell=True)
        command = 'mri_tessellate '              + opj(FS_dir,monkey,'mri', 'filled' + str(value[H]) + '.mgz') + \
        ' ' + opj(FS_dir, monkey, 'surf', Hmin[H] + 'h.orig.nofix')
        spco([command], shell=True)
        command = 'mris_extract_main_component ' + opj(FS_dir,monkey,'surf', Hmin[H] + 'h.orig.nofix') + ' ' + opj(FS_dir, monkey, 'surf', Hmin[H] + 'h.orig.nofix')
        spco([command], shell=True)
        command = 'mris_smooth -nw '             + opj(FS_dir,monkey,'surf', Hmin[H] + 'h.orig.nofix') + ' ' + opj(FS_dir, monkey, 'surf', Hmin[H] + 'h.smoothwm.nofix')
        spco([command], shell=True)
        command = 'mris_inflate -no-save-sulc '  + opj(FS_dir, monkey, 'surf', Hmin[H] + 'h.smoothwm.nofix') + ' ' + \
        opj(FS_dir,monkey,'surf',Hmin[H] + 'h.inflated.nofix')
        spco([command], shell=True)
        command = 'mris_sphere -q  '             + opj(FS_dir,monkey,'surf',Hmin[H] + 'h.inflated.nofix') + ' ' + opj(FS_dir,monkey,'surf',Hmin[H] + 'h.qsphere.nofix')
        spco([command], shell=True)
        
        # 2) Automatic Correction
        shutil.copyfile(opj(FS_dir,monkey,'surf', Hmin[H] + 'h.orig.nofix'),    opj(FS_dir,monkey,'surf', Hmin[H] + 'h.orig'))
        shutil.copyfile(opj(FS_dir,monkey,'surf', Hmin[H] + 'h.inflated.nofix'),opj(FS_dir,monkey,'surf', Hmin[H] + 'h.inflated'))
    
        command = 'mris_fix_topology -mgz -sphere qsphere.nofix -ga ' + monkey + ' ' + Hmin[H] + 'h'
        spco([command], shell=True)
        command = 'mris_euler_number ' + opj(FS_dir,monkey,'surf',Hmin[H] + 'h.orig')
        spco([command], shell=True)
        command = 'mris_remove_intersection ' + opj(FS_dir,monkey,'surf', Hmin[H] + 'h.orig') + ' ' + opj(FS_dir,monkey,'surf', Hmin[H] + 'h.orig')
        spco([command], shell=True)
        
        # 3) Fine white surface;  The coarse white surface is modified according to the 3 images: "wm", "aseg" and "brain".
        command = 'mris_make_surfaces -aseg aseg -whiteonly -noaparc -mgz -T1 brain ' + monkey + ' ' + Hmin[H] + 'h'
        spco([command], shell=True)
        
        os.remove(opj(FS_dir,monkey,'surf', Hmin[H] + 'h.inflated'))
        os.remove(opj(FS_dir,monkey,'mri', 'filled' + str(value[H]) + '.mgz'))
        
        print(surface[H] + ' surface done!')

        
def White_more(FS_dir, monkey, FS_buckner40):

    # This script needs to get acces the two following FS files : 
    #                                                ?h.average.curvature.filled.buckner40.tif
    #                                                ?curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs
    
    TIF       = 'h.average.curvature.filled.buckner40.tif'
    GCS       = 'h.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs'
    
    for H in range(0,2):
        
        # just in case.......
        shutil.copyfile(opj(FS_dir,monkey,'surf',Hmin[H] + 'h.white'),opj(FS_dir,monkey,'surf',Hmin[H] + 'h.white_backup'))
        
        # 1) Inflation and curvature left surface
        command = 'mris_smooth -n 3 -nw ' + opj(FS_dir,monkey,'surf', Hmin[H] + 'h.white') + ' ' + opj(FS_dir,monkey,'surf', Hmin[H] + 'h.smoothwm')
        spco([command], shell=True)
        command = 'mris_inflate ' + opj(FS_dir,monkey,'surf',Hmin[H] + 'h.smoothwm') + ' ' + opj(FS_dir,monkey,'surf', Hmin[H] + 'h.inflated')
        spco([command], shell=True)
        
        # 2) Curvature measurements (Vital for the pial surface)
        command = 'mris_curvature -w ' + opj(FS_dir, monkey, 'surf', Hmin[H] + 'h.white')
        spco([command], shell=True)
        
        command = 'mris_curvature -thresh .999 -n -a 5 -w -distances 10 5 ' + opj(FS_dir,monkey,'surf', Hmin[H] + 'h.inflated')
        spco([command], shell=True)
        #      means thresholding curvature at 99.90% level
        #      normalizing curvature values.
        #      averaging curvature patterns 5 times.
        #      sampling 10 neighbors out to a distance of 10mm (10mm true for human what about monkeys data  ? I chose 5mm from the size factor)
        #                  lh.orig.H = mean curvature
        #	               lh.orig.K = gaussian curvature
        command = 'mris_curvature_stats -m --writeCurvatureFiles -G -o ' + opj(FS_dir,monkey,'stats',Hmin[H] + 'h.curv.stats') + \
        ' -F smoothwm ' + monkey + ' ' + Hmin[H] + 'h curv sulc'
        spco([command], shell=True)
        
        # 3) Sphere registration.
        
        #   3.1) rather time consuming
        command = 'mris_sphere ' + opj(FS_dir,monkey,'surf',Hmin[H] + 'h.inflated') + ' ' + opj(FS_dir,monkey,'surf', Hmin[H] + 'h.sphere')
        spco([command], shell=True)
        
        # 3.2) VERY time consuming (This program registers a surface to an average surface template.) XXXX to check!!!
        ##### adapt to each species if surface template exists!!!! XXXXXX

        DIR = os.getcwd()
        os.chdir(opj(FS_dir,monkey,'surf'))
        command = 'mris_register -curv -dist 5 -max_degrees 68 ' + Hmin[H] + 'h.sphere ' + opj(FS_buckner40, Hmin[H] + TIF) + \
        ' ' + Hmin[H] + 'h.sphere.reg'
        spco([command], shell=True)
        os.chdir(DIR)
        
        # 3.3) store the deformations
        command = 'mris_jacobian ' + opj(FS_dir,monkey,'surf',Hmin[H] + 'h.white') + ' ' + \
        opj(FS_dir,monkey,'surf', Hmin[H] + 'h.sphere.reg') + ' ' + opj(FS_dir,monkey,'surf', Hmin[H] + 'h.jacobian_white')
        spco([command], shell=True)
        
        # 3.4 define what is the cortex: if thi step is not done the T2 correction 
        # (beware it has the tendency to overestimate the cortical thickness and to cross the middline)
        command = 'mris_ca_label -l ' + opj(FS_dir,monkey,'label',Hmin[H] + 'h.cortex.label') + \
        ' -aseg ' + opj(FS_dir,monkey,'mri','aseg.mgz') + \
        ' ' + monkey + ' ' + Hmin[H] + 'h' + \
        ' ' + opj(FS_dir,monkey,'surf',Hmin[H] + 'h.sphere.reg') + \
        ' ' + opj('/','usr','local','freesurfer', '7.4.1', 'average', Hmin[H] + GCS) + \
        ' ' + Hmin[H] + 'h.aparc.annot'
        spco([command], shell=True)
        
        print(surface[H] + ' surface done!')