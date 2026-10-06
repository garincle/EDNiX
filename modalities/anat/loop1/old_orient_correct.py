if deoblique == 'header':

    command = (sing_afni + '3drefit -deoblique' + overwrite + ' -orient ' + orientation +
               ' ' + anatname1 + Timage + '.nii.gz')
    run_cmd.run(command, diary_file)

    shutil.copyfile(anatname1 + Timage + '.nii.gz', anatname2 + Timage + '.nii.gz')

    dictionary = {"Sources": anatname1 + Timage + '.nii.gz',
                  "Description": 'Correction of the fields orientation.', }
    json_object = json.dumps(dictionary, indent=2)
    with open(anatname2 + Timage + '.json', "w") as outfile:
        outfile.write(json_object)

elif deoblique == 'WARPbaboon':

    command = (sing_afni + '3drefit' + overwrite + ' -duporigin ' + list_anat[0] + ' -orient ' + orientation +
               ' ' + anatname1 + Timage + '.nii.gz')
    run_cmd.run(command, diary_file)

    # reorient the fields according to the json file
    command = (sing_afni + '3dWarp ' + overwrite + ' -deoblique -wsinc5 -prefix ' +
               anatname2 + Timage + '.nii.gz' + ' ' + anatname1 + Timage + '.nii.gz')
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": anatname1 + Timage + '.nii.gz',
                  "Description": 'Correction of the fields orientation .', }
    json_object = json.dumps(dictionary, indent=2)
    with open(anatname2 + Timage + '.json', "w") as outfile:
        outfile.write(json_object)

elif deoblique == 'WARP':
    if check == True:
        command = (sing_afni + '3drefit -deoblique' + overwrite + ' -orient ' + orientation +
                   ' ' + anatname1 + Timage + '.nii.gz')
        run_cmd.run(command, diary_file)

        # reorient the fields according to the json file
        command = (sing_afni + '3dWarp' + overwrite +
                   ' -deoblique -NN -prefix ' + anatname2 + Timage + '.nii.gz' +
                   ' ' + anatname1 + Timage + '.nii.gz')
        run_cmd.run(command, diary_file)

        dictionary = {"Sources": anatname1 + Timage + '.nii.gz',
                      "Description": 'Correction of the fields orientation .', }
        json_object = json.dumps(dictionary, indent=2)
        with open(anatname2 + Timage + '.json', "w") as outfile:
            outfile.write(json_object)

        # command = (sing_afni + '3drefit' + overwrite + ' -deoblique -orient ' + orientation +
        # ' ' + opj(dir_prepro, ID + '_mprage_reorient' + Timage + '.nii.gz'))
        # run_cmd.run(command, diary_file)

    else:
        # reorient the fields according to the json file
        command = (sing_afni + '3dWarp' + overwrite +
                   ' -deoblique -NN -prefix ' + anatname2 + Timage + '.nii.gz' +
                   ' ' + anatname1 + Timage + '.nii.gz')
        run_cmd.run(command, diary_file)

        command = (sing_afni + '3drefit' + overwrite + ' -keepcen -orient ' + orientation +
                   ' -duporigin ' + anatname2 + Timage + '.nii.gz' +
                   ' ' + anatname2 + Timage + '.nii.gz')
        run_cmd.run(command, diary_file)

        dictionary = {"Sources": anatname1 + Timage + '.nii.gz',
                      "Description": 'Correction of the fields orientation .', }
        json_object = json.dumps(dictionary, indent=2)
        with open(anatname2 + Timage + '.json', "w") as outfile:
            outfile.write(json_object)

elif deoblique == 'no_deoblique':  # do nothing
    nl = 'exeption1'
    run_cmd.msg(nl, diary_file, 'OKGREEN')

    shutil.copyfile(anatname1 + Timage + '.nii.gz',
                    anatname2 + Timage + '.nii.gz')

    dictionary = {"Sources": anatname1 + Timage + '.nii.gz',
                  "Description": 'Copy .', }
    json_object = json.dumps(dictionary, indent=2)
    with open(anatname2 + Timage + '.json', "w") as outfile:
        outfile.write(json_object)

elif deoblique == 'deob_WO_orient':  # do nothing

    command = (sing_afni + '3drefit -deoblique' + overwrite +
               ' ' + anatname1 + Timage + '.nii.gz')
    run_cmd.run(command, diary_file)

    shutil.copyfile(anatname1 + Timage + '.nii.gz',
                    anatname2 + Timage + '.nii.gz')

    dictionary = {"Sources": anatname1 + Timage + '.nii.gz',
                  "Description": 'Correction of the fields orientation .', }
    json_object = json.dumps(dictionary, indent=2)
    with open(anatname2 + Timage + '.json', "w") as outfile:
        outfile.write(json_object)

elif deoblique == 'WARP_without_3drefit':  # re-alineate
    nl = 'exeption2'
    ### add something else not usefull!! XX
    run_cmd.msg(nl, diary_file, 'OKGREEN')

    command = (sing_afni + '3dWarp' + overwrite + ' -deoblique -NN -prefix ' + anatname2 + Timage + '.nii.gz' +
               ' ' + anatname1 + Timage + '.nii.gz')
    run_cmd.run(command, diary_file)

    nl = 'need to realign or use just one'
    run_cmd.msg(nl, diary_file, 'OKGREEN')

    dictionary = {"Sources": anatname1 + Timage + '.nii.gz',
                  "Description": 'deoblique .', }
    json_object = json.dumps(dictionary, indent=2)
    with open(anatname2 + Timage + '.json', "w") as outfile:
        outfile.write(json_object)