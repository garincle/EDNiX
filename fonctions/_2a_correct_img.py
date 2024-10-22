import os
import subprocess
import nibabel as nb
import numpy as np
from math import pi
from fonctions.extract_filename import extract_filename
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

def correct_img(dir_fMRI_Refth_RS_prepro1, RS, list_map, RS_map, study_fMRI_Refth, i, r, overwrite,s_bind,afni_sif,fsl_sif,topup_file):

	root = extract_filename(RS_map[i])
	root_RS = extract_filename(RS[r])

	#copy map imag in new location
	command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + list_map[i] + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, RS_map[i]) + ' -expr "a"'
	spco([command], shell=True)

	#mean of the map img
	command = 'singularity run' + s_bind + afni_sif + '3dTstat' + overwrite + ' -mean -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean_pre.nii.gz') + ' ' + \
			  opj(dir_fMRI_Refth_RS_prepro1, RS_map[i])
	spco([command], shell=True)

	# register each volume to the base image
	command = 'singularity run' + s_bind + afni_sif + '3dvolreg' + overwrite + ' -verbose -zpad 1 -base ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean_pre.nii.gz') + \
			  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_map_align.nii.gz') + \
			  ' -cubic ' + \
			  opj(dir_fMRI_Refth_RS_prepro1, RS_map[i])
	spco([command], shell=True)

	'''
	# realignment intra-run
	command = 'singularity run' + s_bind + fsl_sif + 'mcflirt -in ' + opj(dir_fMRI_Refth_RS_prepro1, RS_map[i]) + \
	' -out ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_map_align.nii.gz') + \
	' -mats -plots -reffile ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean_pre.nii.gz') + ' -rmsrel -rmsabs -spline_final'
	spco([command], shell=True)
	'''
	#mean of the map img to ref img
	command = 'singularity run' + s_bind + afni_sif + '3dTstat' + overwrite + ' -mean -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1, RS_map[i])
	spco([command], shell=True)

	#####################FRANCK????????????????????????????,
	#dell? ==> do the job ?? opj(dir_fMRI_Refth_RS_prepro1, RS[r].replace('.nii.gz','_xdtr_mean.nii.gz'))
	#command = '3dTcat -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,RS[int(ref_nb)-1].replace('.nii.gz','_fMRI_Ref.nii.gz')) + ' ' + opj(dir_fMRI_Refth_RS_prepro1, RS[int(ref_nb)-1]) + '[0-9]'
	#spco([command], shell=True)
	#command = '3dTstat -mean -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, RS_map[i].replace('.nii.gz','_fMRI_Ref_mean.nii.gz')) + \
	#' ' + opj(dir_fMRI_Refth_RS_prepro1,RS[int(ref_nb)-1].replace('.nii.gz','_fMRI_Ref.nii.gz'))
	#spco([command], shell=True)
	#os.remove(opj(dir_fMRI_Refth_RS_prepro1,RS[int(ref_nb)-1].replace('.nii.gz','_fMRI_Ref.nii.gz')))
	###resemple anat to func #XXX change opj(dir_prepro, ID + '_mprage_reorient_NU.nii.gz') for opj(dir_prepro,ID + '_mprage_reorient_NU.nii.gz')

	####### a foutu la merde sans raison!!!!!
	#command = '3dresample -master ' + opj(dir_fMRI_Refth_RS_prepro1, RS[r].replace('.nii.gz','_xdtr_mean.nii.gz')) + \
	#' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1, RS_map[i].replace('.nii.gz','_map_mean_reso.nii.gz')) + \
	#' -input ' + opj(dir_fMRI_Refth_RS_prepro1, RS_map[i].replace('.nii.gz','_map_mean.nii.gz'))
	#spco([command], shell=True)

	command = 'singularity run' + s_bind + afni_sif + '3dTcat' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_se.nii.gz') + \
	' ' + opj(dir_fMRI_Refth_RS_prepro1, opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean.nii.gz')) + \
	' ' + opj(dir_fMRI_Refth_RS_prepro1, opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean.nii.gz'))
	spco([command], shell=True)

	#####correct image for topup (i.e. remove the slices that do not fit topup requirement)
	#fslroi <input> <output> <xmin> <xsize> <ymin> <ysize> <zmin> <zsize>
	#command = 'fslroi ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_se.nii.gz') + \
	#' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_se1.nii.gz') + \
	#' 0 -1 0 -1 1 132'
	#' 1 77 0 -1 0 -1'
	#spco([command], shell=True)

	#https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;67dcb45c.1209
	#I agree with Matt that you probably have an odd number of voxels in one direction (usually in the slice direction). 
	#In topup, for various reasons, the images dimensions has to be an integer multiple of each sub-sampling level one uses. 
	#We usually just throw away the top or bottom slice (provided it is outside the brain) in these cases.
	"""
	if 'x' in correction_direction:
		intofencod = 0

	elif 'y' in correction_direction:
		intofencod = 1

	elif 'z' in correction_direction:
		intofencod = 2
	"""
	#### zeropad?? add a slice instead of removing!!!
	im = nb.load(opj(dir_fMRI_Refth_RS_prepro1, root + '_se.nii.gz'))
	imdata = im.get_fdata()
	s = imdata.shape
	dests = np.array(s)
	hdr = im.header.copy()
	hdr.set_data_shape(imdata.shape)
	for b, d in enumerate(s):
		if b < 3:
			if (d % 2) == 0:
				print(bcolors.OKGREEN + "{0} est paire, no need to remove a slice" + bcolors.ENDC)
			else:
				print(bcolors.OKGREEN + "{0} est impaire, we will have to remove a slice" + bcolors.ENDC)
				imdata = imdata.take(range(d - 1), axis=b)
	
	nb.Nifti1Image(imdata, im.affine, hdr).to_filename(opj(dir_fMRI_Refth_RS_prepro1, root + '_se1.nii.gz'))

	### se_map don't change but 1 -1
	### b02b0 don't change 

	command = 'singularity run' + s_bind + fsl_sif + 'topup --imain=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_se1.nii.gz') + \
	' --datain=' + topup_file[0] + \
	' --config=' + topup_file[1] + \
	' --fout=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap.nii.gz') + \
	' --iout=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_unwarped.nii.gz')
	spco([command], shell=True)

	##### for fugue
	command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap.nii.gz') + ' -mul ' + str(2*pi) + \
	' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_rads.nii.gz')
	spco([command], shell=True)

	command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_unwarped.nii.gz') + \
	' -Tmean ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_mag.nii.gz')
	spco([command], shell=True)

