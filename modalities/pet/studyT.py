import os
import shutil
import subprocess
import glob
import ants
import json

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

spco = subprocess.check_output
spgo = subprocess.getoutput

from Tools import run_cmd,getpath,extract_filename


def make_template(study_template_atlas_folder,transfomat,tracor,all_ID, all_data_path, reference,diary_file):

    nl = 'Run make_template'
    run_cmd.msg(nl, diary_file,'HEADER')

    template_list = []
    mat_list      = []
    fname=''

    for ID,  data_path in zip(all_ID, all_data_path):
        (_, _, _, _, path_pet,
         _, _, _, dir_prepro_orig_process, _,
         _, _, _, _, _, _, _, _, _, _, _, dir_prepro_orig_matrices) = getpath.pet(data_path, reference)

        list_img = sorted(glob.glob(opj(dir_prepro_orig_process,
                                        '_'.join([ID, 'trc-' + tracor, '*desc-reorient', 'pet.nii.gz']))))
        if list_img:
            for i in range(len(list_img)):
                imghd = ants.image_header_info(list_img[i])
                if imghd['nDimensions'] == 4:
                    list_img[i] = list_img[i].replace('desc-reorient', 'desc-reorient-mean')

            template_list.append(list_img[0])
            mat_list.append(dir_prepro_orig_matrices)
        else :
            nl = ('check the folder: ' + dir_prepro_orig_process + 'there are no available file named ' +
                  '_'.join(['xxx','desc-reorient', 'pet.nii.gz']))
            run_cmd.msg(nl, diary_file, 'WARNING')

    if template_list:
        os.mkdir(study_template_atlas_folder,exist_ok=True)
        os.mkdir(opj(study_template_atlas_folder,'pet'), exist_ok=True)
        os.mkdir(opj(study_template_atlas_folder, 'pet','matrices'), exist_ok=True)

        nii = ants.build_template(initial_template=None,
                                  image_list=template_list,
                                  iterations=5,
                                  gradient_step=0.2,
                                  blending_weight=0.75,
                                  weights=None,
                                  useNoRigid=True,
                                  type_of_transform=transfomat,
                                  output_dir=opj(study_template_atlas_folder, 'pet','matrices'))

        fname = opj(study_template_atlas_folder, 'studyTemplate' +
                                  str(len(template_list)) + '_trc-' + tracor + '_pet.nii.gz')
        ants.image_write(nii, fname)


        for img,matdir in zip(template_list,mat_list):

            reg = ants.registration(fixed=nii,
                                    moving=img,
                                    type_of_transform=transfomat,
                                    outprefix=opj(matdir,'trc-' + tracor +  "_to_studyTemplate_"))

    return fname

def manual_crop(dir,image,suffix,diary_name,sing_fs):

    new_image = opb(image.replace('_pet.nii.gz', '_'+ suffix + '_pet.nii.gz'))

    shutil.copyfile(image,opj(dir,new_image))
    cmd = (sing_fs + 'freeview -v ' + opj(dir,new_image) + ' -subtitle "CROP-THE-IMAGE ! "')
    run_cmd.do(cmd, diary_name)
    print('CROP : open Tools/Crop Volume')
    print('For brain image : keep the skin around the skull, the eyeballs (if available), remove most of the external soft tissues '
          '(ears and tip of the nose) as well as the lower jaw')
    print('Once satisfied, do : Apply and Save as :' + opj(dir,new_image))

    dictionary = {"Sources": image,
                  "Description": 'Manually Cropped image (Freeview, Freesurfer).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir,new_image.replace('.nii.gz', '.json')), "w") as outfile:
        outfile.write(json_object)


def automatic_crop(dir_prepro_orig_process,file,Tmatrices,refcrop,interp_pet):

    fname, descname = extract_filename.get_desc(opb(file), 'desc')

    newdir = opj(dir_prepro_orig_process,'forMC')

    imghd = ants.image_header_info(file)
    if imghd['nDimensions'] == 4:
        imgtype = 3
    else:
        imgtype = 0

    if descname :
        img_crop = opj(newdir, '_'.join([fname, 'desc-' + descname.split('-')[-1] + '-crop', 'pet.nii.gz']))
    else:
        img_crop = opj(newdir, '_'.join([fname, 'desc-crop', 'pet.nii.gz']))

    img = ants.image_read(file)
    ref = ants.image_read(refcrop)
    moved = ants.apply_transforms(fixed=img, moving=ref, transformlist=Tmatrices, interpolator=interp_pet,imagetype=imgtype)
    moved[moved > 0.1] = 1

    # apply crop
    cropped_img = ants.crop_image(image=img, label_image=moved, label=1)

    ants.image_write(cropped_img, img_crop)
    dictionary = {"Sources": [file,
                              refcrop,
                              Tmatrices],
                  "Description": 'Manually Cropped image (Freeview, Freesurfer).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(img_crop.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)

    return img_crop