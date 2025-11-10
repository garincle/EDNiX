import os
import glob
import json
import ants
import nibabel as nib

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd
from Tools import check_nii
from Tools.extract_filename import extract_filename
from Tools import get_orientation
from anatomical.loop1 import myelin

##########################################
########  define orientation    ##########
##########################################

def correct_orient(BIDStype,
                   listTimage,
                   path_rawanat,
                   path_anat, 
                   ID, 
                   Session, 
                   otheranat,
                   listTimage_orig,
                   type_norm,
                   MNIBcorrect_indiv,
                   orientation, 
                   dir_prepro,
                   dir_transfo,
                   IgotbothT1T2,
                   force_myelin_same_space,
                   list_transfo,
                   overwrite,
                   animalP,
                   humanP,
                   sing_afni,sing_fs, sing_wb,
                   diary_file):

    nl = 'Run  anatomical._1_correct_orient.correct_orient'
    run_cmd.msg(nl, diary_file,'HEADER')
    if force_myelin_same_space== True:
        nl = 'force_myelin_same_space== True, I hope that you know what you are doing. Your two anat images, better being in the same space, otherwise it will mess everything up!!'
        run_cmd.msg(nl,diary_file,'WARNING')

    anatname1 = opj(dir_prepro, ID + '_space-raw_desc-template_')
    anatname2 = opj(dir_prepro, ID + '_space-raw_desc-reorient_')

    if MNIBcorrect_indiv == 'N3':
        anatname3 = opj(dir_prepro, ID + '_space-raw_desc-n3Bias_')
    elif MNIBcorrect_indiv == 'N4':
        anatname3 = opj(dir_prepro, ID + '_space-raw_desc-n4Bias_')

    raw_myelin         = anatname1 + 'T1wdividedbyT2w.nii.gz'
    reoriented_myelin  = anatname2 + 'T1wdividedbyT2w.nii.gz'

    if not ope(path_anat):
        os.makedirs(path_anat)

    if len(listTimage)>1:
        if len(humanP) != len(listTimage):
            if len(humanP) == 1:
                humanP = humanP*len(listTimage)
                run_cmd.msg('it is assumed that the subject position in the scanner remained the same in every image ', diary_file, 'WARNING')
            else:
                nl = 'check the parameters: ther must be a position set for each type of images. If you do not know, leave it empty'
                run_cmd.msg(nl, diary_file, 'WARNING')

        if len(animalP) != len(listTimage):
            if len(animalP) == 1:
                animalP = animalP * len(listTimage)
                run_cmd.msg('it is assumed that the subject position in the scanner remained the same in every image ',
                            diary_file, 'WARNING')
            else:
                nl = 'check the parameters: ther must be a position set for each type of images. If you do not know, leave it empty'
                run_cmd.msg(nl, diary_file, 'WARNING')
    print(listTimage_orig)
    for Im in range(len(listTimage_orig)):
        # get variables from json
        if BIDStype == 1:
            list_anat = sorted(glob.glob(opj(path_rawanat, 'sub-' + ID + '_ses-' + str(Session) + '_*' + listTimage_orig[Im] + '.nii*')))
            print(list_anat)
        elif BIDStype == 2:
            list_anat = sorted(glob.glob(opj(path_rawanat, 'sub-' + ID + '_' + str(listTimage_orig[Im]) + '.nii*')))
        else:
            nl = 'WARNING: BIDStype was not 1 or 2 so it must be define by a string that you provided'
            run_cmd.msg(nl, diary_file, 'WARNING')
            list_anat = sorted(glob.glob(opj(path_rawanat, BIDStype.format(ID=ID, Session=Session, Timage=listTimage_orig[Im]))))

        if len(list_anat)>0:
            if humanP[0] == '':
                if opi(list_anat[0].replace('.nii.gz','json')):
                    f = open(list_anat[0].replace('.nii.gz','json'))
                    info = json.load(f)
                    try:
                        humanP[Im] = info['PatientPosition']
                        run_cmd.msg('the subject was scanned with the following position parameter : ' + humanPosition,
                                    diary_file, 'OKGREEN')
                    except:
                        humanP[Im] = humanP[Im]

        if len(list_anat)==1:
            nl = 'INFO: We found only one anat images for this session'
            run_cmd.msg(nl, diary_file, 'OKGREEN')

            nl = 'list of the anat found is:' + str(list_anat)
            run_cmd.msg(nl, diary_file, 'OKGREEN')

            img = nib.load(list_anat[0])
            dims = img.shape

            if humanP[0] == '':
                if opi(list_anat[0].replace('.nii.gz','json')):
                    f = open(list_anat[0].replace('.nii.gz','json'))
                    info = json.load(f)
                    try:
                        humanPosition = info['PatientPosition']
                        run_cmd.msg('the subject was scanned with the following position parameter : ' + humanPosition,
                                    diary_file, 'OKGREEN')
                    except:
                        humanPosition = humanP[0]

            if len(dims) > 3:
                os.makedirs(opj(dir_prepro,'tmp'))
                img = ants.image_read(list_anat[0])
                unmerged = ants.ndimage_to_list(img)
                motion_corrected = list()
                for i in range(len(unmerged)):
                    mtx = ants.registration(fixed=unmerged[0], moving=unmerged[i], type_of_transform='Rigid')
                    new_img = ants.apply_transforms(fixed=unmerged[0], moving=unmerged[i],
                                                    transformlist=mtx['fwdtransforms'],
                                                    interpolator='hammingWindowedSinc')
                    #motion_corrected.append(mtx['warpedmovout']) this way we cannot set the interpolator
                    motion_corrected.append(new_img)
                motCorr = ants.list_to_ndimage(img, motion_corrected)
                new_img = ants.get_average_of_timeseries(motCorr)
                ants.image_write(new_img,  anatname1 + listTimage[Im] + '.nii.gz', ri=False)
                dictionary = {"Sources": list_anat[0],
                              "Description": 'motion correction and merge.', }
                json_object = json.dumps(dictionary, indent=2)
                with open(anatname1 + listTimage[Im] + '.json', "w") as outfile:
                        outfile.write(json_object)

            else:
                command = (sing_afni + '3dcalc -a ' + list_anat[0] +
                           ' -prefix ' + anatname1 + listTimage[Im] + '.nii.gz' + ' -expr "a"' + overwrite)
                run_cmd.do(command, diary_file)

                dictionary = {"Sources": list_anat[0],
                              "Description": 'Nothing done.', }
                json_object = json.dumps(dictionary, indent=2)
                with open(anatname1 + listTimage[Im] + '.json', "w") as outfile:
                    outfile.write(json_object)


        elif len(list_anat)>1:

            nl = 'INFO: We found ' + str(len(list_anat)) + ' anat images for this session'
            run_cmd.msg(nl, diary_file, 'OKGREEN')

            nl = 'list of the anat found is:' + str(list_anat)
            run_cmd.msg(nl, diary_file, 'OKGREEN')

            ###################################
            ###########remove norm anat #########
            ###################################

            for i, anat in enumerate(list_anat):
                if opi(opj(opd(list_anat[i]), extract_filename(list_anat[i]) + '.json')):
                    read_json = opj(opd(list_anat[i]), extract_filename(list_anat[i]) + '.json')
                else:
                    nl = "WARNING: No .json associated to the the anat image !! you might want to check that!"
                    run_cmd.msg(nl, diary_file, 'WARNING')

                    ff = open(read_json)
                    anat_T1 = json.load(ff)
                    try:
                        ImageType = anat_T1["ImageType"]

                        if 'NORM' in ImageType:
                            nl = 'Removing ' + extract_filename(list_anat[i]) + ' as it is a NORM'
                            run_cmd.msg(nl, diary_file, 'WARNING')

                            list_anat.pop(i)
                            # After removing the file, no need to increment i since the next file takes the place of the removed one
                    except:
                        nl = "WARNING No ImageType, could not check if NORM image was included in your BIDS, you might want to check that!"
                        run_cmd.msg(nl, diary_file, 'WARNING')

            if len(list_anat) > 1:
                ANAT= ' '.join(list_anat)
                command = (sing_fs + 'mri_robust_template --mov ' + ANAT +
                           ' --template ' + anatname1 + listTimage[Im] + '.nii.gz' +
                           ' --satit --inittp 1 --fixtp --noit --iscale --average 0')
                run_cmd.run(command,diary_file)

                dictionary = {"Sources": ANAT,
                              "Description": 'Co-registation and Average (mri_robust_template from Freesurfer)) .', }
                json_object = json.dumps(dictionary, indent=2)
                with open(anatname1 + listTimage[Im] + '.json', "w") as outfile:
                    outfile.write(json_object)
            else:

                command = (sing_afni + '3dcalc -a ' + list_anat[0] +
                           ' -prefix ' + anatname1 + listTimage[Im] + '.nii.gz' + ' -expr "a"' + overwrite)
                run_cmd.do(command, diary_file)

                dictionary = {"Sources": list_anat[0],
                              "Description": 'Nothing done.', }
                json_object = json.dumps(dictionary, indent=2)
                with open(anatname1 + listTimage[Im] + '.json', "w") as outfile:
                    outfile.write(json_object)

        else:
            nl = 'HELP: anat is detected with BIDStype=' + str(BIDStype) + ' but the command: '
            run_cmd.msg(nl, diary_file, 'FAIL')
            if BIDStype == 1:
                list_anat = "glob.glob(" + str(opj(path_rawanat, 'sub-' + ID + '_ses-' + str(Session) + '_*' + listTimage_orig[Im] + '.nii*'))
            elif BIDStype == 2:
                list_anat = '"glob.glob(' + str(opj(path_rawanat,'sub-' + ID + '_' + listTimage_orig[Im] + '.nii*'))
            else:
                list_anat = "glob.glob(" + str(opj(path_rawanat, BIDStype.format(ID=ID, Session=Session, Timage=listTimage_orig[Im])))
            nl = str(list_anat) + ' failed'
            raise Exception(run_cmd.error(nl, diary_file))


    if IgotbothT1T2 == True:
        # beware that it doesn't make sense if the T2 image is a FLAIR sequence
        if opi(anatname1 + type_norm + '.nii.gz') and opi(anatname1 + otheranat + '.nii.gz'):
            if force_myelin_same_space==False:
                check = myelin.check(anatname1 + type_norm + '.nii.gz', anatname1 + otheranat + '.nii.gz',
                                     sing_afni, diary_file)
            else:
                check = True
            #### if T2 or T1 not mentionned in the name then it's a mess... XXX
            for Timage in listTimage:
                if   'T1' in Timage: T1 = Timage
                elif 'T2' in Timage: T2 = Timage

            if check == True or force_myelin_same_space==True:
                myelin.divT1T2v2(T1, T2, type_norm, anatname1,
                                 raw_myelin, diary_file, sing_wb, sing_afni)

    ListimageMyeline = listTimage.copy()
    if opi(raw_myelin):
        ListimageMyeline.append('T1wdividedbyT2w')
        animalP.append(animalP[0])
        humanP.append(humanP[0])


    for Im in range(len(ListimageMyeline)):

        if animalP[Im] == '':
            animalP[Im] = 'humanlike'

        cmd = sing_afni + '3dinfo -orient ' + anatname1 + ListimageMyeline[Im] + '.nii.gz'
        msg, _ = run_cmd.get(cmd, diary_file)
        orient = msg.decode("utf-8").split('\n')[-2]

        cmd = sing_afni + '3dinfo -is_oblique ' + anatname1 + ListimageMyeline[Im] + '.nii.gz'
        msg, _ = run_cmd.get(cmd, diary_file)
        obli = msg.decode("utf-8").split('\n')[-2]

        if obli == '1':
            cmd = sing_afni + '3dWarp -overwrite -deoblique -prefix ' + anatname2 + ListimageMyeline[Im] + '.nii.gz' + ' ' + anatname1 + ListimageMyeline[Im] + '.nii.gz'
            run_cmd.run(cmd, diary_file)
            desc = 'Correction of the obliquity.'

            cmd = sing_afni + '3dinfo -orient ' + anatname2 + ListimageMyeline[Im] + '.nii.gz'
            msg, _ = run_cmd.get(cmd, diary_file)

            orient = msg.decode("utf-8").split('\n')[-2]

        else:
            cmd = sing_afni + '3dcalc -overwrite -a ' + anatname1 + ListimageMyeline[Im] + '.nii.gz' + ' -prefix ' + anatname2 + ListimageMyeline[Im] + '.nii.gz' + ' -expr "a"'
            run_cmd.do(cmd, diary_file)
            desc = 'Copy .'


        if not orientation == '':
            cmd = sing_afni + '3drefit -overwrite -orient ' + orientation + ' ' + anatname2 + ListimageMyeline[Im] + '.nii.gz'
            run_cmd.run(cmd, diary_file)
        else:
            if not animalP[Im] == 'humanlike':
                neworient = get_orientation.getreal(humanP[Im], animalP[Im], orient)
                run_cmd.msg('the new orientation will be : ' + neworient, diary_file, 'OKGREEN')

                cmd = sing_afni + '3drefit -overwrite -orient ' + neworient + ' ' + anatname2 + ListimageMyeline[Im] + '.nii.gz'
                run_cmd.run(cmd, diary_file)
                desc = 'Correction of the fields orientation.'

        dictionary = {"Sources": anatname1 + ListimageMyeline[Im] + '.nii.gz',
                      "Description": desc, }
        json_object = json.dumps(dictionary, indent=2)
        with open(anatname2 + ListimageMyeline[Im] + '.json', "w") as outfile:
            outfile.write(json_object)


    if IgotbothT1T2 == True:
        # beware that it doesn't make sense if the T2 image is a FLAIR sequence:
        if not opi(raw_myelin):
                nl =  'The validity of that result may be questionable ....'
                run_cmd.msg(nl, diary_file, 'WARNING')

                myelin.newversion(anatname2 + T1 + '.nii.gz', anatname2 + T2 + '.nii.gz',type_norm,list_transfo,
                                 reoriented_myelin, dir_transfo,diary_file,sing_afni,sing_wb)

        #### ensure that both image are now in the same space, with the same header
        check_nii.resamp(anatname2 + otheranat + '.nii.gz', anatname2 + type_norm + '.nii.gz', 'anat', '', '',
                         diary_file,
                         sing_wb)

    for Timage in listTimage:
        # Apply N4 bias field correction
        BRAIN = ants.image_read(anatname2 + Timage + '.nii.gz')
        if MNIBcorrect_indiv == 'N4':
            BC = ants.n4_bias_field_correction(BRAIN, shrink_factor=2, convergence={'iters': [50, 50, 30, 20], 'tol': 1e-7})
            nl = 'N4 Bias Field Correction done'

        elif MNIBcorrect_indiv == 'N3':
            BC = ants.n3_bias_field_correction2(BRAIN, convergence={'iters': 150, 'tol': 0.000})
            nl = 'N3 Bias Field Correction done'

        # Write the output image
        ants.image_write(BC, anatname3 + Timage + '.nii.gz', ri=False)

        dictionary = {"Sources": anatname2 + Timage + '.nii.gz',
                      "Description": 'B1 bias field non uniformity correction (' + MNIBcorrect_indiv + ' from ANTspy).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(anatname2 + Timage + '.json', "w") as outfile:
            outfile.write(json_object)

        run_cmd.msg(nl, diary_file, 'OKGREEN')





