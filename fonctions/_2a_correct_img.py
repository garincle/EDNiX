import os
import subprocess
import nibabel as nb
import numpy as np
from math import pi
from fonctions.extract_filename import extract_filename
import datetime
import json

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

def correct_img(dir_fMRI_Refth_RS_prepro1, RS, list_map, RS_map, i, r, recordings,
				overwrite,s_bind,afni_sif,fsl_sif,topup_file,diary_file):

	ct = datetime.datetime.now()
	diary = open(diary_file, "a")
	diary.write(f'\n{ct}')
	nl = '##  Working on step ' + str(2) + '(function: _2a_correct_img).  ##'
	print(bcolors.OKGREEN + nl + bcolors.ENDC)
	diary.write(f'\n{nl}')

	if recordings == '2_mapdir':
		for z in [0, 1]:
			root = extract_filename(RS_map[z])

			# copy map imag in new location
			command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + list_map[z] + \
					  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, RS_map[z]) + ' -expr "a"'
			nl = spgo(command)
			diary.write(f'\n{nl}')
			print(nl)
			dictionary = {"Sources": list_map[z],
						  "Description": 'Copy.', }
			json_object = json.dumps(dictionary, indent=2)
			with open(opj(dir_fMRI_Refth_RS_prepro1, RS_map[z][:-7] + '.json'), "w") as outfile:
				outfile.write(json_object)

			# mean of the map img
			command = 'singularity run' + s_bind + afni_sif + '3dTstat' + overwrite + ' -mean -prefix ' + \
					  opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean_pre' + str(z) + '.nii.gz') + ' ' + \
					  opj(dir_fMRI_Refth_RS_prepro1, RS_map[z])
			nl = spgo(command)
			diary.write(f'\n{nl}')
			print(nl)
			dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, RS_map[z]),
						  "Description": '4D Mean (3dTstat, AFNI).', }
			json_object = json.dumps(dictionary, indent=2)
			with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean_pre' + str(z) + '.json'), "w") as outfile:
				outfile.write(json_object)

			# register each volume to the base image
			command = 'singularity run' + s_bind + afni_sif + '3dvolreg' + overwrite + ' -verbose -zpad 1 -base ' + \
					  opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean_pre' + str(z) + '.nii.gz') + \
					  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_map_align' + str(z) + '.nii.gz') + \
					  ' -cubic ' + \
					  opj(dir_fMRI_Refth_RS_prepro1, RS_map[z])
			nl = spgo(command)
			diary.write(f'\n{nl}')
			print(nl)
			dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1, RS_map[z]),
									  opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean_pre' + str(z) + '.nii.gz')],
						  "Description": 'Rigid Realignment (3dvolreg,AFNI).', }
			json_object = json.dumps(dictionary, indent=2)
			with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_map_align' + str(z) + '.json'), "w") as outfile:
				outfile.write(json_object)

			# mean of the map img to ref img
			command = 'singularity run' + s_bind + afni_sif + '3dTstat' + overwrite + ' -mean -prefix ' + \
					  opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean' + str(z) + '.nii.gz') + ' ' + \
					  opj(dir_fMRI_Refth_RS_prepro1, root + '_map_align' + str(z) + '.nii.gz')
			nl = spgo(command)
			diary.write(f'\n{nl}')
			print(nl)
			dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_map_align' + str(z) + '.nii.gz'),
						  "Description": '4D Mean (3dTstat, AFNI).', }
			json_object = json.dumps(dictionary, indent=2)
			with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean' + str(z) + '.json'), "w") as outfile:
				outfile.write(json_object)


		root = extract_filename(RS_map[0])
		root1 = extract_filename(RS_map[1])

		command = 'singularity run' + s_bind + afni_sif + '3dTcat' + overwrite + \
				  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_se.nii.gz') + \
				  ' ' + opj(dir_fMRI_Refth_RS_prepro1, root1 + '_map_mean' + str(1) + '.nii.gz') + \
				  ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean' + str(0) + '.nii.gz')
		nl = spgo(command)
		diary.write(f'\n{nl}')
		print(nl)
		dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1, root1 + '_map_mean' + str(1) + '.nii.gz'),
								  opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean' + str(0) + '.nii.gz')],
					  "Description": 'concatenation in 4D (3dTcat,AFNI).', }
		json_object = json.dumps(dictionary, indent=2)
		with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_se.json'), "w") as outfile:
			outfile.write(json_object)

	else:
		root    = extract_filename(RS_map[i])
		root_rs = extract_filename(RS[r])

		#copy map imag in new location
		command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + list_map[i] + \
			  	' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, RS_map[i]) + ' -expr "a"'
		nl = spgo(command)
		diary.write(f'\n{nl}')
		print(nl)
		dictionary = {"Sources": list_map[i],
					  "Description": 'Copy.', }
		json_object = json.dumps(dictionary, indent=2)
		with open(opj(dir_fMRI_Refth_RS_prepro1, RS_map[i][:-7] + '.json'), "w") as outfile:
			outfile.write(json_object)

		# mean of the map img
		command = 'singularity run' + s_bind + afni_sif + '3dTstat' + overwrite + ' -mean -prefix ' \
			  	+ opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean_pre.nii.gz') + ' ' + \
			  	opj(dir_fMRI_Refth_RS_prepro1, RS_map[i])
		nl = spgo(command)
		diary.write(f'\n{nl}')
		print(nl)
		dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, RS_map[i]),
					  "Description": '4D Mean (3dTstat, AFNI).', }
		json_object = json.dumps(dictionary, indent=2)
		with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean_pre.json'), "w") as outfile:
			outfile.write(json_object)

		# register each volume to the base image
		command = 'singularity run' + s_bind + afni_sif + '3dvolreg' + overwrite + ' -verbose -zpad 1 -base ' + \
			  	opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean_pre.nii.gz') + \
			  	' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_map_align.nii.gz') + \
			  	' -cubic ' + \
			  	opj(dir_fMRI_Refth_RS_prepro1, RS_map[i])
		nl = spgo(command)
		diary.write(f'\n{nl}')
		print(nl)
		dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1, RS_map[i]),
								  opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean_pre.nii.gz')],
					  "Description": 'Rigid Realignment (3dvolreg,AFNI).', }
		json_object = json.dumps(dictionary, indent=2)
		with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_map_align.json'), "w") as outfile:
			outfile.write(json_object)

		# mean of the map img to ref img
		command = 'singularity run' + s_bind + afni_sif + '3dTstat' + overwrite + ' -mean -prefix ' + \
			  	opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean.nii.gz') + ' ' + \
			  	opj(dir_fMRI_Refth_RS_prepro1, root + '_map_align.nii.gz')
		nl = spgo(command)
		diary.write(f'\n{nl}')
		print(nl)
		dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_map_align.nii.gz'),
					  "Description": '4D Mean (3dTstat, AFNI).', }
		json_object = json.dumps(dictionary, indent=2)
		with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean.json'), "w") as outfile:
			outfile.write(json_object)


		command = 'singularity run' + s_bind + afni_sif + '3dTcat' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_se.nii.gz') + \
		' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean.nii.gz') + \
		' ' + opj(dir_fMRI_Refth_RS_prepro1, root_rs + '_xdtr_mean.nii.gz')
		nl = spgo(command)
		diary.write(f'\n{nl}')
		print(nl)
		dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1, root + '_map_mean.nii.gz'),
								  opj(dir_fMRI_Refth_RS_prepro1, root_rs + '_xdtr_mean.nii.gz')],
					  "Description": 'concatenation in 4D (3dTcat,AFNI).', }
		json_object = json.dumps(dictionary, indent=2)
		with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_se.json'), "w") as outfile:
			outfile.write(json_object)











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

	#### zeropad ?? add a slice instead of removing!!!

	im = nb.load(opj(dir_fMRI_Refth_RS_prepro1, root + '_se.nii.gz'))
	imdata = im.get_fdata()
	s = imdata.shape

	hdr = im.header.copy()
	hdr.set_data_shape(imdata.shape)
	for b, d in enumerate(s):
		if b < 3:
			if (d % 2) == 0:
				nl = "{0} is even, no need to remove a slice"

			else:
				nl = "{0} is odd, we will have to remove a slice"
				imdata = imdata.take(range(d - 1), axis=b)

			print(bcolors.OKGREEN + nl + bcolors.ENDC)
			diary.write(f'\n{nl}')

	nb.Nifti1Image(imdata, im.affine, hdr).to_filename(opj(dir_fMRI_Refth_RS_prepro1, root + '_se1.nii.gz'))
	dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_se.nii.gz'),
				  "Description": 'Make sure that the number of voxels are even in each dimension. (Nifti1Image, nilearn)', }
	json_object = json.dumps(dictionary, indent=2)
	with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_se1.json'), "w") as outfile:
		outfile.write(json_object)


	### se_map don't change but 1 -1
	### b02b0 don't change 

	command = 'singularity run' + s_bind + fsl_sif + 'topup --imain=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_se1.nii.gz') + \
	' --datain=' + topup_file[0] + \
	' --config=' + topup_file[1] + \
	' --fout=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap.nii.gz') + \
	' --iout=' + opj(dir_fMRI_Refth_RS_prepro1, root + '_unwarped.nii.gz')
	nl = spgo(command)
	diary.write(f'\n{nl}')
	print(nl)
	dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1, root + '_se1.nii.gz'),
							  topup_file[0],
							  topup_file[1]],
				  "Description": 'Create fieldmaps. (topup, FSL)', }
	json_object = json.dumps(dictionary, indent=2)
	with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap.json'), "w") as outfile:
		outfile.write(json_object)

	##### for fugue
	command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap.nii.gz') + ' -mul ' + str(2*pi) + \
	' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_rads.nii.gz')
	nl = spgo(command)
	diary.write(f'\n{nl}')
	print(nl)
	dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap.nii.gz'),
				  "Description": 'conversion into rads. (fslmaths, FSL)', }
	json_object = json.dumps(dictionary, indent=2)
	with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_rads.json'), "w") as outfile:
		outfile.write(json_object)


	command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_unwarped.nii.gz') + \
	' -Tmean ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_mag.nii.gz')
	nl = spgo(command)
	diary.write(f'\n{nl}')
	print(nl)
	dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_unwarped.nii.gz'),
				  "Description": '4D Mean (fslmaths, FSL).', }
	json_object = json.dumps(dictionary, indent=2)
	with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_fieldmap_mag.json'), "w") as outfile:
		outfile.write(json_object)

	diary.write(f'\n')
	diary.close()


