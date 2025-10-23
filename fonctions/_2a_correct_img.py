import os
import nibabel as nb
from math import pi
import json

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname

from Tools import run_cmd
from fonctions.extract_filename import extract_filename

def correct_img(dir_prepro_orig_process, dir_prepro_fmap, fMRI_runMean_n4Bias, RS, list_map, RS_map, i, r, recordings,
				overwrite,sing_afni,sing_fsl,topup_file,diary_file):

	nl = '##  Working on step ' + str(2) + '(function: _2a_correct_img).  ##'
	run_cmd.msg(nl, diary_file, 'HEADER')

	if recordings == '2_mapdir':
		for z in [0, 1]:
			root = extract_filename(RS_map[z])
			fmap_mean = opj(dir_prepro_fmap, root + '_fmap_mean' + str(z) + '.nii.gz')
			fmap_mean0 = opj(dir_prepro_fmap, root + '_fmap_mean0.nii.gz')
			fmap_mean1 = opj(dir_prepro_fmap, root + '_fmap_mean1.nii.gz')

			fmri_to_fmap_align = opj(dir_prepro_fmap, root + '_fmri_to_fmap_align' + str(z) + '.nii.gz')
			fmri_to_fmap_align_Mean = opj(dir_prepro_fmap, root + '_fmri_to_fmap_align_Mean' + str(z) + '.nii.gz')
			fmri_to_fmap_concat = opj(dir_prepro_fmap, root + '_fmri_to_fmap_concat' + str(z) + '.nii.gz')

			# mean of the map img
			command = (sing_afni + '3dTstat' + overwrite + ' -mean -prefix ' +
					   fmap_mean + ' ' +
					   list_map[z])
			run_cmd.run(command, diary_file)

			dictionary = {"Sources": list_map[z],
						  "Description": '4D Mean (3dTstat, AFNI).',
						  "Command": command,}
			json_object = json.dumps(dictionary, indent=3)

			with open(opj(dir_prepro_orig_process, root + '_map_mean_pre' + str(z) + '.json'), "w") as outfile:
				outfile.write(json_object)

			# register each volume to the base image
			command = (sing_afni + '3dvolreg' + overwrite + ' -verbose -zpad 1 -base ' +
					   fmap_mean +
					   ' -prefix ' + fmri_to_fmap_align +
					   ' -cubic ' + opj(dir_prepro_orig_process, RS_map[z]))
			run_cmd.run(command, diary_file)

			dictionary = {"Sources": [opj(dir_prepro_orig_process, RS_map[z]),
									  fmap_mean],
						  "Description": 'Rigid Realignment (3dvolreg,AFNI).', "Command": command,}
			json_object = json.dumps(dictionary, indent=3)
			with open(fmri_to_fmap_align.replace('.nii.gz','json'), "w") as outfile:
				outfile.write(json_object)

			# mean of the map img to ref img
			command = (sing_afni + '3dTstat' + overwrite + ' -mean -prefix ' +
					   fmri_to_fmap_align_Mean + ' ' +
					   fmri_to_fmap_align)
			run_cmd.run(command, diary_file)

			dictionary = {"Sources": fmri_to_fmap_align,
						  "Description": '4D Mean (3dTstat, AFNI).', "Command": command,}
			json_object = json.dumps(dictionary, indent=3)
			with open(fmri_to_fmap_align_Mean.replace('.nii.gz','json'), "w") as outfile:
				outfile.write(json_object)

		root = extract_filename(RS_map[0])
		root1 = extract_filename(RS_map[1])

		command = (sing_afni + '3dTcat' + overwrite +
				   ' -prefix ' + fmri_to_fmap_concat +
				   ' ' + fmap_mean1 +
				   ' ' + fmap_mean0)
		run_cmd.run(command, diary_file)

		dictionary = {"Sources": [fmap_mean1,
								  fmap_mean0],
					  "Description": 'concatenation in 4D (3dTcat,AFNI).', "Command": command,}
		json_object = json.dumps(dictionary, indent=3)
		with open(fmri_to_fmap_concat.replace('.nii.gz','json'), "w") as outfile:
			outfile.write(json_object)

	else:
		root    = extract_filename(RS_map[i])
		root_rs = extract_filename(RS[r])
		
		fmap_mean = opj(dir_prepro_fmap, root + '_fmap_mean.nii.gz')
		fmri_to_fmap_align = opj(dir_prepro_fmap, root + '_fmri_to_fmap_align.nii.gz')
		fmri_to_fmap_align_Mean = opj(dir_prepro_fmap, root + '_fmri_to_fmap_align_Mean.nii.gz')
		fmri_to_fmap_concat = opj(dir_prepro_fmap, root + '_fmri_to_fmap_concat.nii.gz')
		
		# mean of the map img
		command = (sing_afni + '3dTstat' + overwrite + ' -mean -prefix '
				   + fmap_mean + ' ' +
				   opj(dir_prepro_orig_process, RS_map[i]))
		run_cmd.run(command, diary_file)

		dictionary = {"Sources": opj(dir_prepro_orig_process, RS_map[i]),
					  "Description": '4D Mean (3dTstat, AFNI).', "Command": command,}
		json_object = json.dumps(dictionary, indent=3)
		with open(fmap_mean.replace('.nii.gz','json'), "w") as outfile:
			outfile.write(json_object)

		# register each volume to the base image
		command = (sing_afni + '3dvolreg' + overwrite + ' -verbose -zpad 1 -base ' +
				   fmap_mean +
				   ' -prefix ' + fmri_to_fmap_align +
				   ' -cubic ' + opj(dir_prepro_orig_process, RS_map[i]))
		run_cmd.run(command, diary_file)

		dictionary = {"Sources": [opj(dir_prepro_orig_process, RS_map[i]),
								  fmap_mean],
					  "Description": 'Rigid Realignment (3dvolreg,AFNI).', "Command": command,}
		json_object = json.dumps(dictionary, indent=3)
		
		with open(fmri_to_fmap_align.replace('.nii.gz','json'), "w") as outfile:
			outfile.write(json_object)

		# mean of the map img to ref img
		command = (sing_afni + '3dTstat' + overwrite + ' -mean -prefix ' +
				   fmri_to_fmap_align_Mean +
				   ' ' + fmri_to_fmap_align)
		run_cmd.run(command, diary_file)

		dictionary = {"Sources": fmri_to_fmap_align,
					  "Description": '4D Mean (3dTstat, AFNI).', "Command": command,}
		json_object = json.dumps(dictionary, indent=3)
		with open(fmri_to_fmap_align_Mean.replace('.nii.gz','json'), "w") as outfile:
			outfile.write(json_object)

		command = (sing_afni + '3dTcat' + overwrite + ' -prefix ' + fmri_to_fmap_concat +
				   ' ' + fmri_to_fmap_align_Mean +
				   ' ' + fMRI_runMean_n4Bias)
		run_cmd.run(command, diary_file)

		dictionary = {"Sources": [fmri_to_fmap_align_Mean,
								  fMRI_runMean_n4Bias],
					  "Description": 'concatenation in 4D (3dTcat,AFNI).', "Command": command,}
		json_object = json.dumps(dictionary, indent=3)
		with open(fmri_to_fmap_concat.replace('.nii.gz','json'), "w") as outfile:
			outfile.write(json_object)

	fmri_to_fmap_even = opj(dir_prepro_fmap, root + '_fmri_to_fmap_even.nii.gz')
	fMRI_runMean_fieldmap = opj(dir_prepro_fmap, root + '_space-func_desc-runMean_fieldmap.nii.gz')
	fMRI_runMean_fieldmap_rads = opj(dir_prepro_fmap, root + '_space-func_desc-runMean_fieldmap_rads.nii.gz')
	fMRI_runMean_fieldmap_mag = opj(dir_prepro_fmap, root + '_space-func_desc-runMean_fieldmap_mag.nii.gz')
	fMRI_b0_distortion_corrected = opj(dir_prepro_fmap, root + '_b0_distortion_corrected.nii.gz')
	# make sure that the number of voxels are even in each dimension
	im = nb.load(fmri_to_fmap_concat)
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

			run_cmd.msg(nl, diary_file, 'OKGREEN')

	nb.Nifti1Image(imdata, im.affine, hdr).to_filename(fmri_to_fmap_even)
	dictionary = {"Sources": fmri_to_fmap_concat,
				  "Description": 'Make sure that the number of voxels are even in each dimension. (Nifti1Image, nilearn)', "Command": command,}
	json_object = json.dumps(dictionary, indent=3)
	with open(fmri_to_fmap_even.replace('.nii.gz','json'), "w") as outfile:
		outfile.write(json_object)

	### se_map don't change but 1 -1
	### b02b0 don't change
	command = (sing_fsl + 'topup --imain=' + fmri_to_fmap_even +
			   ' --datain=' + topup_file[0] +
			   ' --config=' + topup_file[1] +
			   ' --fout=' + fMRI_runMean_fieldmap +
			   ' --iout=' + fMRI_b0_distortion_corrected)
	run_cmd.run(command, diary_file)

	dictionary = {"Sources": [fmri_to_fmap_even,
							  topup_file[0],
							  topup_file[1]],
				  "Description": 'Create fieldmaps. (topup, FSL)', "Command": command,}
	json_object = json.dumps(dictionary, indent=3)
	with open(fMRI_runMean_fieldmap.replace('.nii.gz','json'), "w") as outfile:
		outfile.write(json_object)

	##### for fugue
	command = (sing_fsl + 'fslmaths ' + fMRI_runMean_fieldmap + ' -mul ' + str(2*pi) +
			   ' ' + fMRI_runMean_fieldmap_rads)
	run_cmd.run(command, diary_file)

	dictionary = {"Sources": fMRI_runMean_fieldmap,
				  "Description": 'conversion into rads. (fslmaths, FSL)', "Command": command,}
	json_object = json.dumps(dictionary, indent=3)
	with open(fMRI_runMean_fieldmap_rads.replace('.nii.gz','json'), "w") as outfile:
		outfile.write(json_object)

	command = (sing_fsl + 'fslmaths ' + fMRI_b0_distortion_corrected +
			   ' -Tmean ' + fMRI_runMean_fieldmap_mag)
	run_cmd.run(command, diary_file)

	dictionary = {"Sources": fMRI_b0_distortion_corrected,
				  "Description": '4D Mean (fslmaths, FSL).', "Command": command,}
	json_object = json.dumps(dictionary, indent=3)
	with open(fMRI_runMean_fieldmap_mag.replace('.nii.gz','json'), "w") as outfile:
		outfile.write(json_object)



