import pandas as pd
import glob
import os
import subprocess

spco = subprocess.check_output
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists

output_atlas = '/srv/projects/easymribrain/data/Atlas/13_Atlas_project/Atlases_V2/'
input_atlas = '/srv/projects/easymribrain/data/Atlas/13_Atlas_project/Atlases_V2/'
orig_atlases = '/srv/projects/easymribrain/data/Atlas/13_Atlas_project/0_Atlas_modify/Atlas/'
WM_subspace = '/srv/projects/easymribrain/data/Atlas/13_Atlas_project/WM_GM/'
atlas_exel = '/srv/projects/easymribrain/data/Atlas/13_Atlas_project/Classiff/Classif_VDual_newGlass.xlsx'

species = [
    'Pig',
    'Mouse_lemur',
    'Marmoset',
    'Macaque',
    'Chimpanzee',
    'CatinDog',
    'DoginCat',
    'Mouse',
    'RatWHS',
    'Human']

species_list = ['Pig', 'Mouse_lemur', 'Marmoset', 'Macaque', 'Chimpanzee', 'DoginCat', 'CatinDog', 'Mouse', 'RatWHS', 'Human']
imposeWM_all = [
    'No',
    'No',
    'No',
    'YES',
    'YES',
    'YES',
    'YES',
    'No',
    'No',
    'No']

middle_value = ['-0.929', '0.091','0.075','0.125','0.5','-0.352','-0.289','-0.05','0','0']

MAIN_PATH = r'/srv/projects/easymribrain/code/EDNiX/'  # Update this path
atlas_folder = '/srv/projects/easymribrain/data/Atlas/'
s_bind = f' --bind {atlas_folder},{MAIN_PATH} '
s_path = opj(MAIN_PATH, 'Tool_library', 'Singularity')
afni_sif = opj(s_path, 'afni_ub24_latest.sif ')

###################### re-order altas based on excel file

for animal in species:
    template = glob.glob(output_atlas + animal + '/template*.nii.gz')
    command = 'singularity run' + s_bind + afni_sif + '3dresample -master ' + str(template[0]) + ' -prefix ' + output_atlas + animal + '/atlas_resample_template.nii.gz -input ' + output_atlas + animal + '/atlas.nii.gz -overwrite'
    spco(command, shell=True)
    if animal == 'DoginCat':
        animal_excel = pd.read_excel(atlas_exel, sheet_name='CatinDog')
    else:
        animal_excel = pd.read_excel(atlas_exel, sheet_name=animal)

    for numberlvl in [1, 2, 3, 4]:
        list_new_label = list(animal_excel['Number_level_' + str(numberlvl)])
        list_old_label = list(animal_excel['Number_level_orig'])
        string_build_atlas = str('')
        for new_label, old_label in zip(list_new_label, list_old_label):
            string_build_atlas = string_build_atlas + '(' + str(int(new_label)) + '*(equals(a,' + str(
                int(old_label)) + ')))+'
        string_build_atlas2 = "'" + string_build_atlas[:-1] + "'"
        if not os.path.exists(output_atlas + animal): os.mkdir(output_atlas + animal)

        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + output_atlas + animal + '/atlas_resample_template.nii.gz ' + ' -expr ' + string_build_atlas2 + ' -prefix ' + output_atlas + animal + '/atlaslvl' + str(
            numberlvl) + '.nii.gz -overwrite'
        spco(command, shell=True)


for animal, imposeWM in zip(species, imposeWM_all):
    if animal in species:
        print('1')

        if not os.path.exists(output_atlas + animal): os.mkdir(output_atlas + animal)
        numberlvl = 1
        # add WM
        if imposeWM == "YES":
            if os.path.exists(orig_atlases + animal + '/WM.nii.gz'):
                atlas1 = orig_atlases + animal + '/WM.nii.gz'
                atlas2 = input_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'
                atlas3 = output_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'

                command = 'singularity run' + s_bind + afni_sif + '3dresample -overwrite -master ' + atlas2 + ' -input ' + atlas1 + ' -prefix ' + atlas1
                spco(command, shell=True)

                if os.path.exists(atlas3):
                    command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + atlas3 + ' -b ' + atlas1 + ' -expr ' + '"ifelse(step(b)*6,step(b)*6,a)"' + ' -prefix ' + atlas3 + ' -overwrite'
                    spco(command, shell=True)

                if not os.path.exists(atlas3):
                    command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + atlas2 + ' -b ' + atlas1 + ' -expr ' + '"ifelse(step(b)*6,step(b)*6,a)"' + ' -prefix ' + atlas3 + ' -overwrite'
                    spco(command, shell=True)

        # Tronc
        if os.path.exists(orig_atlases + animal + '/tronc.nii.gz'):
            atlas1 = orig_atlases + animal + '/tronc.nii.gz'
            atlas2 = input_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'
            atlas3 = output_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'

            command = 'singularity run' + s_bind + afni_sif + '3dresample -overwrite -master ' + atlas2 + ' -input ' + atlas1 + ' -prefix ' + atlas1
            spco(command, shell=True)

            if os.path.exists(atlas3):
                command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + atlas3 + ' -b ' + atlas1 + ' -expr ' + '"ifelse(step(b)*13,step(b)*13,a)"' + ' -prefix ' + atlas3 + ' -overwrite'
                spco(command, shell=True)

            if not os.path.exists(atlas3):
                command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + atlas2 + ' -b ' + atlas1 + ' -expr ' + '"ifelse(step(b)*13,step(b)*13,a)"' + ' -prefix ' + atlas3 + ' -overwrite'
                spco(command, shell=True)

        # add WM_cerereb_missing
        print('4')
        if os.path.exists(orig_atlases + animal + '/WMcb.nii.gz'):
            print('WMcb')
            atlas1 = orig_atlases + animal + '/WMcb.nii.gz'
            atlas2 = input_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'
            atlas3 = output_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'

            command = 'singularity run' + s_bind + afni_sif + '3dresample -overwrite -master ' + atlas2 + ' -input ' + atlas1 + ' -prefix ' + atlas1
            spco(command, shell=True)

            if os.path.exists(atlas3):
                command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + atlas3 + ' -b ' + atlas1 + ' -expr ' + '"ifelse(step(b)*5,step(b)*5,a)"' + ' -prefix ' + atlas3 + ' -overwrite'
                spco(command, shell=True)

            if not os.path.exists(atlas3):
                command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + atlas2 + ' -b ' + atlas1 + ' -expr ' + '"ifelse(step(b)*5,step(b)*5,a)"' + ' -prefix ' + atlas3 + ' -overwrite'
                spco(command, shell=True)

        # add Ventricules
        if os.path.exists(orig_atlases + animal + '/Ventricules.nii.gz'):
            atlas1 = orig_atlases + animal + '/Ventricules.nii.gz'
            atlas2 = input_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'
            atlas3 = output_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'

            command = 'singularity run' + s_bind + afni_sif + '3dresample -overwrite -master ' + atlas2 + ' -input ' + atlas1 + ' -prefix ' + atlas1
            spco(command, shell=True)

            if os.path.exists(atlas3):
                command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + atlas3 + ' -b ' + atlas1 + ' -expr ' + '"ifelse(step(b)*4,step(b)*4,a)"' + ' -prefix ' + atlas3 + ' -overwrite'
                spco(command, shell=True)

            if not os.path.exists(atlas3):
                command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + atlas2 + ' -b ' + atlas1 + ' -expr ' + '"ifelse(step(b)*4,step(b)*4,a)"' + ' -prefix ' + atlas3 + ' -overwrite'
                spco(command, shell=True)

        # add Ventri 4
        if os.path.exists(orig_atlases + animal + '/Ventri3.nii.gz'):
            atlas1 = orig_atlases + animal + '/Ventri3.nii.gz'
            atlas2 = input_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'
            atlas3 = output_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'

            command = 'singularity run' + s_bind + afni_sif + '3dresample -overwrite -master ' + atlas2 + ' -input ' + atlas1 + ' -prefix ' + atlas1
            spco(command, shell=True)

            if os.path.exists(atlas3):
                command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + atlas3 + ' -b ' + atlas1 + ' -expr ' + '"ifelse(step(b)*2,step(b)*2,a)"' + ' -prefix ' + atlas3 + ' -overwrite'
                spco(command, shell=True)

            if not os.path.exists(atlas3):
                command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + atlas2 + ' -b ' + atlas1 + ' -expr ' + '"ifelse(step(b)*2,step(b)*2,a)"' + ' -prefix ' + atlas3 + ' -overwrite'
                spco(command, shell=True)

        # add WM
        if imposeWM == "No":
            if os.path.exists(orig_atlases + animal + '/WM.nii.gz'):
                atlas1 = orig_atlases + animal + '/WM.nii.gz'
                atlas2 = input_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'
                atlas3 = output_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'

                command = 'singularity run' + s_bind + afni_sif + '3dresample -overwrite -master ' + atlas2 + ' -input ' + atlas1 + ' -prefix ' + atlas1
                spco(command, shell=True)

                if os.path.exists(atlas3):
                    command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + atlas3 + ' -b ' + atlas1 + ' -expr ' + '"ifelse(a,a,step(b)*6)"' + ' -prefix ' + atlas3 + ' -overwrite'
                    spco(command, shell=True)

                if not os.path.exists(atlas3):
                    command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + atlas2 + ' -b ' + atlas1 + ' -expr ' + '"ifelse(a,a,step(b)*6)"' + ' -prefix ' + atlas3 + ' -overwrite'
                    spco(command, shell=True)

        # add CSF remaining

        orig_brain_mask = orig_atlases + animal + '/brain_mask.nii.gz'
        atlas1 = orig_atlases + animal + '/Mask_dilate.nii.gz'
        atlas2 = input_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'
        atlas3 = output_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'

        ####
        templates = glob.glob(output_atlas + animal + '/template*.nii.gz')

        numberlvl = 2
        atlas2 = input_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'
        atlas3 = output_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite' + ' -a ' + atlas2 + ' -expr "a"' + ' -prefix ' + atlas3
        spco(command, shell=True)
        numberlvl = 3
        atlas2 = input_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'
        atlas3 = output_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite' + ' -a ' + atlas2 + ' -expr "a"' + ' -prefix ' + atlas3
        spco(command, shell=True)
        numberlvl = 4
        atlas2 = input_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'
        atlas3 = output_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz'
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite' + ' -a ' + atlas2 + ' -expr "a"' + ' -prefix ' + atlas3
        spco(command, shell=True)

        atlas2 = input_atlas + animal + '/atlas.nii.gz'
        atlas3 = output_atlas + animal + '/atlas.nii.gz'
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite' + ' -a ' + atlas2 + ' -expr "a"' + ' -prefix ' + atlas3
        spco(command, shell=True)

        atlas2 = input_atlas + animal + '/atlaslvl1.nii.gz'
        atlas3 = output_atlas + animal + '/Wmask.nii.gz'
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite' + ' -a ' + atlas2 + ' -expr "equals(a,6)"' + ' -prefix ' + atlas3
        spco(command, shell=True)
        command = 'singularity run' + s_bind + afni_sif + '3dresample -master ' + str(
            templates[0]) + ' -prefix ' + atlas3 + ' -input ' + atlas3 + ' -overwrite'
        spco(command, shell=True)

        atlas2 = input_atlas + animal + '/atlaslvl1.nii.gz'
        atlas3 = output_atlas + animal + '/Gmask.nii.gz'
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite' + ' -a ' + atlas2 + ' -expr "equals(a,8) + equals(a,9) + equals(a,10) + equals(a,11) + equals(a,12)"' + ' -prefix ' + atlas3
        spco(command, shell=True)
        command = 'singularity run' + s_bind + afni_sif + '3dresample -master ' + str(
            templates[0]) + ' -prefix ' + atlas3 + ' -input ' + atlas3 + ' -overwrite'
        spco(command, shell=True)

        atlas2 = input_atlas + animal + '/brain_mask.nii.gz'
        atlas3 = output_atlas + animal + '/brain_mask.nii.gz'
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite' + ' -a ' + atlas2 + ' -expr "a"' + ' -prefix ' + atlas3
        spco(command, shell=True)
        command = 'singularity run' + s_bind + afni_sif + '3dresample -master ' + str(
            templates[0]) + ' -prefix ' + atlas3 + ' -input ' + atlas3 + ' -overwrite'
        spco(command, shell=True)

        for template in templates:
            command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + atlas2 + ' -b ' + str(
                template) + ' -expr "step(a)*b" -prefix ' + str(template)[:-7] + '_SS.nii.gz -overwrite'
            spco(command, shell=True)


###LR

for animal, middle_val in zip(species_list, middle_value):
	for numberlvl in [1, 2, 3, 4]:

		command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + output_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz' + ' -expr "ispositive(x-' + middle_val + ')*a"' + ' -prefix ' + output_atlas + animal + '/atlaslvl' + str(numberlvl) + '_L.nii.gz -overwrite -nscale'
		spco(command, shell=True)

		command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + output_atlas + animal + '/atlaslvl' + str(numberlvl) + '.nii.gz' + ' -expr "isnegative(x-' + middle_val + ')*a"' + ' -prefix ' + output_atlas + animal + '/atlaslvl' + str(numberlvl) + '_R.nii.gz -overwrite -nscale'
		spco(command, shell=True)

		command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + output_atlas + animal + '/atlaslvl' + str(numberlvl) + '_R.nii.gz' + ' -expr "step(a)*1000"' + ' -prefix ' + output_atlas + animal + '/atlaslvl' + str(numberlvl) + '_R1000.nii.gz -overwrite -nscale'
		spco(command, shell=True)

		command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -b ' + output_atlas + animal + '/atlaslvl' + str(numberlvl) + '_R1000.nii.gz ' + '-a ' + output_atlas + animal + '/atlaslvl' + str(numberlvl) + '_R.nii.gz ' + '-expr "b+a"'+ ' -prefix ' + output_atlas + animal + '/atlaslvl' + str(numberlvl) + '_R.nii.gz -overwrite -nscale'
		spco(command, shell=True)

		command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -b ' + output_atlas + animal + '/atlaslvl' + str(numberlvl) + '_L.nii.gz ' + '-a ' + output_atlas + animal + '/atlaslvl' + str(numberlvl) + '_R.nii.gz ' + '-expr "b+a"'+ ' -prefix ' + output_atlas + animal + '/atlaslvl' + str(numberlvl) + '_LR.nii.gz -overwrite -nscale'
		spco(command, shell=True)

		os.remove(output_atlas + animal + '/atlaslvl' + str(numberlvl) + '_L.nii.gz')
		os.remove(output_atlas + animal + '/atlaslvl' + str(numberlvl) + '_R.nii.gz')
		os.remove(output_atlas + animal + '/atlaslvl' + str(numberlvl) + '_R1000.nii.gz')


for animal in species_list:
	for numberlvl in [1, 2, 3, 4]:

		if numberlvl==1:

			list_old_label = [6, 1006, 8, 1008, 9, 1009, 10, 1010, 7, 1007, 13, 1013, 1, 1001, 4, 1004, 5, 1005, 2, 1002]
			list_new_label = [2, 41, 3, 42, 3, 42, 3, 42, 8, 47, 16, 16, 24, 24, 4, 43, 7, 46, 14, 14]

			string_build_atlas = str('')
			for new_label, old_label in zip(list_new_label, list_old_label):
				string_build_atlas = string_build_atlas + '(' + str(int(new_label)) + '*(equals(a,' + str(int(old_label)) + ')))+'
			string_build_atlas2 = "'" + string_build_atlas[:-1] + "'"
			if not os.path.exists(output_atlas + animal) : os.mkdir(output_atlas + animal)

			command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + input_atlas + animal + '/atlaslvl' + str(numberlvl) + '_LR.nii.gz' + ' -expr ' + string_build_atlas2 + ' -prefix ' + output_atlas + animal + '/atlas_forSEG_final.nii.gz -overwrite'
			spco(command, shell=True)

		if numberlvl==2:
			list_old_label = [27, 1027, 12, 1012]
			list_new_label = [17, 53, 10, 49]

			string_build_atlas = str('')
			for new_label, old_label in zip(list_new_label, list_old_label):
				string_build_atlas = string_build_atlas + '(' + str(int(new_label)) + '*(equals(a,' + str(int(old_label)) + ')))+'
			string_build_atlas2 = "'" + string_build_atlas[:-1] + "'"
			if not os.path.exists(output_atlas + animal) : os.mkdir(output_atlas + animal)

			command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + input_atlas + animal + '/atlaslvl' + str(numberlvl) + '_LR.nii.gz' + ' -expr ' + string_build_atlas2 + ' -prefix ' + output_atlas + animal + '/atlas_forSEG' + str(numberlvl) + '_ADD.nii.gz -overwrite'
			spco(command, shell=True)

			command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + output_atlas + animal + '/atlas_forSEG_final.nii.gz' + ' -b ' + output_atlas + animal + '/atlas_forSEG' + str(numberlvl) + '_ADD.nii.gz' + ' -expr ' + '"ifelse(b,b,a)"' + ' -prefix ' + output_atlas + animal + '/atlas_forSEG_final.nii.gz -overwrite'
			spco(command, shell=True)

		if numberlvl==3:
			list_old_label = [77, 1077, 79, 1079, 80, 1080, 76, 1076, 75, 1075, 75, 1075, 75, 1075]
			list_new_label = [13, 52, 18, 54, 28, 60, 3, 42, 11, 50, 12, 51, 26, 58]

			string_build_atlas = str('')
			for new_label, old_label in zip(list_new_label, list_old_label):
				string_build_atlas = string_build_atlas + '(' + str(int(new_label)) + '*(equals(a,' + str(int(old_label)) + ')))+'
			string_build_atlas2 = "'" + string_build_atlas[:-1] + "'"
			if not os.path.exists(output_atlas + animal) : os.mkdir(output_atlas + animal)

			command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + input_atlas + animal + '/atlaslvl' + str(numberlvl) + '_LR.nii.gz' + ' -expr ' + string_build_atlas2 + ' -prefix ' + output_atlas + animal + '/atlas_forSEG' + str(numberlvl) + '_ADD.nii.gz -overwrite'
			spco(command, shell=True)

			command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + output_atlas + animal + '/atlas_forSEG_final.nii.gz' + ' -b ' + output_atlas + animal + '/atlas_forSEG' + str(numberlvl) + '_ADD.nii.gz' + ' -expr ' + '"ifelse(b,b,a)"' + ' -prefix ' + output_atlas + animal + '/atlas_forSEG_final.nii.gz -overwrite'
			spco(command, shell=True)

		if numberlvl==4:
			list_old_label = [163, 1163, 164, 1164, 165, 1165]
			list_new_label = [11, 50, 12, 51, 26, 58]

			string_build_atlas = str('')
			for new_label, old_label in zip(list_new_label, list_old_label):
				string_build_atlas = string_build_atlas + '(' + str(int(new_label)) + '*(equals(a,' + str(int(old_label)) + ')))+'
			string_build_atlas2 = "'" + string_build_atlas[:-1] + "'"
			if not os.path.exists(output_atlas + animal) : os.mkdir(output_atlas + animal)

			command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + input_atlas + animal + '/atlaslvl' + str(numberlvl) + '_LR.nii.gz' + ' -expr ' + string_build_atlas2 + ' -prefix ' + output_atlas + animal + '/atlas_forSEG' + str(numberlvl) + '_ADD.nii.gz -overwrite'
			spco(command, shell=True)

			command = 'singularity run' + s_bind + afni_sif + '3dresample -overwrite -master ' + input_atlas + animal + '/atlaslvl' + str(numberlvl) + '_LR.nii.gz' + \
			' -input ' + output_atlas + animal + '/atlas_forSEG_final.nii.gz' + \
			' -prefix ' + output_atlas + animal + '/atlas_forSEG_final.nii.gz'

			command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + output_atlas + animal + '/atlas_forSEG_final.nii.gz' + ' -b ' + output_atlas + animal + '/atlas_forSEG' + str(numberlvl) + '_ADD.nii.gz' + ' -expr ' + '"ifelse(b,b,a)"' + ' -prefix ' + output_atlas + animal + '/atlas_forSEG_final.nii.gz -overwrite'
			spco(command, shell=True)

		list_old_label = [2, 41, 10, 49, 4, 43]
		list_new_label = [255, 127, 255, 127, 255, 127]

		string_build_atlas = str('')
		for new_label, old_label in zip(list_new_label, list_old_label):
			string_build_atlas = string_build_atlas + '(' + str(int(new_label)) + '*(equals(a,' + str(int(old_label)) + ')))+'
		string_build_atlas2 = "'" + string_build_atlas[:-1] + "'"
		if not os.path.exists(output_atlas + animal) : os.mkdir(output_atlas + animal)

		#command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + output_atlas + animal + '/atlas_forSEG_final.nii.gz ' + ' -expr ' + string_build_atlas2 + ' -prefix ' + output_atlas + animal + '/WM_mask_LR_freesurfer.nii.gz -overwrite'
		#spco(command, shell=True)

	os.remove(output_atlas + animal + '/atlas_forSEG2_ADD.nii.gz')
	os.remove(output_atlas + animal + '/atlas_forSEG3_ADD.nii.gz')
	os.remove(output_atlas + animal + '/atlas_forSEG4_ADD.nii.gz')

print("replace mice by the Seg copy file!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

###question of Aseg_refLR
def get_middle_value_for_species(species_list, value_list, target_species):
    for i, specie in enumerate(species_list):
        if specie == target_species:
            return value_list[i]
    return None

for species, middle_val in zip(species_list, middle_value):
	define_center = get_middle_value_for_species(species_list, middle_value, species)
	### if it doesn't exists let's make it!!! but you need an Aseg_ref!!
	####atlases files
	Aseg_ref    = opj(input_atlas,species,'atlas_forSEG_final.nii.gz')
	Aseg_refLR  = opj(input_atlas,species,'atlas_forSEG_final_LR.nii.gz')

	for numberlvl in [1, 2, 3, 4]:
		if not ope(Aseg_refLR):
		    command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + Aseg_ref + ' -expr "step(ispositive(x-' + define_center + ')*a)" -prefix ' + opd(
		        Aseg_refLR) + '/Aseg_ref_L.nii.gz'
		    spco([command], shell=True)
		    command = 'singularity run' + s_bind + afni_sif + '3dcalc -a ' + opd(
		        Aseg_refLR) + '/Aseg_ref_L.nii.gz -b ' + Aseg_ref + ' -expr "ifelse(a, a*255,step(b)*127)" -prefix ' + Aseg_refLR
		    spco([command], shell=True)
