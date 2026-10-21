#import
import os
import shutil
import ants
import json

opj = os.path.join
opd = os.path.dirname
opb = os.path.basename
ope = os.path.exists
opi = os.path.isfile


# for afni:
deob_interp   = '-wsinc5'
resamp_interp = 'Cu'
# for ants:
antsaff     = '0GenericAffine.mat'
T1interp    = 'hammingWindowedSinc'
Boldinterp  = 'lanczosWindowedSinc'
labelinterp = 'genericLabel'



from Tools import run_cmd,extract_filename, get_orientation

def reorient_apply(input, output, humanPosition,animalPosition, suffix, diary_name, sing_afni):
  nl = ' reorientation'
  run_cmd.msg(nl, diary_name,'OKGREEN')

  # Check if file exists
  if not opi(input):
    raise FileNotFoundError(f"Input file not found: {input}")

  # case very unlikely from pet image reconstruction
  cmd = (sing_afni + '3dinfo -is_oblique ' + input)
  msg, _ = run_cmd.get(cmd, diary_name)
  obli = msg.decode("utf-8").split('\n')[-2]

  if int(obli) == 1:
    cmd = (sing_afni + '3dWarp -deoblique ' + deob_interp + ' -prefix ' + output + ' ' + input)
    run_cmd.run(cmd, diary_name)
    desc = 'Deoblique and Reorientation of the header according to standard view (3dWarp and 3drefit from AFNI).'

  else:
    shutil.copy(input, output)
    desc = 'Reorientation of the header according to standard view (3drefit from AFNI).'

  orient_img, _, _, _ = get_orientation.use_afni(output, sing_afni)
  print(orient_img)

  if not animalPosition == 'na':
    new_orient = get_orientation.getreal(humanPosition, animalPosition, orient_img)
    print(new_orient)
    cmd = (sing_afni + '3drefit -orient ' + new_orient + ' ' + output)
    run_cmd.run(cmd, diary_name)

  dictionary = {"Sources": input,
                "Description": desc, }
  json_object = json.dumps(dictionary, indent=2)
  with open(output.replace('.nii.gz', '.json'), "w") as outfile:
    outfile.write(json_object)

  imghd = ants.image_header_info(output)
  if imghd['nDimensions'] == 4:
    directory = opd(output)
    newoutput = extract_filename.change(directory, output, suffix, 'mean')
    cmd = (sing_afni + '3dTstat -mean -prefix ' + newoutput + ' ' + output)
    run_cmd.run(cmd, diary_name)
    dictionary = {"Sources": output,
                  "Description": 'mean image (3dTstat from AFNI).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(newoutput.replace('.nii.gz', '.json'), "w") as outfile:
      outfile.write(json_object)


def change(input, input2, output, suffix, diary_name, sing_afni):
  nl = ' reorientation'
  run_cmd.msg(nl, diary_name,'OKGREEN')

  shutil.copy(input, output)
  desc = 'Reorientation of the header according to standard view (3drefit from AFNI).'

  orient_img, _, _, _ = get_orientation.use_afni(input2, sing_afni)

  cmd = (sing_afni + '3drefit -orient ' + orient_img + ' ' + output)
  run_cmd.run(cmd, diary_name)

  dictionary = {"Sources": input,
                "Description": desc, }
  json_object = json.dumps(dictionary, indent=2)
  with open(output.replace('.nii.gz', '.json'), "w") as outfile:
    outfile.write(json_object)

  imghd = ants.image_header_info(output)
  if imghd['nDimensions'] == 4:
    directory = opd(output)
    newoutput = extract_filename.change(directory, output, suffix, 'mean')
    cmd = (sing_afni + '3dTstat -mean -prefix ' + newoutput + ' ' + output)
    run_cmd.run(cmd, diary_name)
    dictionary = {"Sources": output,
                  "Description": 'mean image (3dTstat from AFNI).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(newoutput.replace('.nii.gz', '.json'), "w") as outfile:
      outfile.write(json_object)



def doit(new_dir,img,humanPosition,animalPosition,Ref_native_raw,Need_check,pethd,diary_name,sing_afni,sing_fs):

  nl = 'REORIENT THE FILE'
  run_cmd.msg(nl, diary_name,'OKGREEN')

  for fneeded in [img, Ref_native_raw]:
    if not opi(fneeded):
      raise FileNotFoundError(f"Input file not found: {fneeded}")

  fname,_ = extract_filename.get_desc(opb(img),'desc')

  step2 = opj(new_dir, '_'.join([fname, 'desc-reorient-mean', 'pet.nii.gz']))
  step3 = opj(new_dir, '_'.join([fname, 'desc-reorient', 'pet.nii.gz']))

  # reorient the fields according to the animalPosition
  if not animalPosition == 'na':
    if not pethd == 'wrong':
      reorient_apply(img,step3,humanPosition,animalPosition, 'desc-reorient', diary_name, sing_afni)
    else:
      change(img, Ref_native_raw,step3, 'desc-reorient', diary_name, sing_afni)

  else:
    shutil.copyfile(img,step3)

  if opi(img.replace('.nii.gz','.sif')):
    shutil.copyfile(img.replace('.nii.gz','.sif'), step3.replace('.nii.gz','.sif'))

  imghd = ants.image_header_info(img)
  if imghd['nDimensions'] == 4:
    to_see = step2
  else:
    to_see = step3

  print('better check the orientation to be correct')
  if Need_check ==1:
    nl = 'check of the orientation'
    run_cmd.msg(nl, diary_name,'OKGREEN')
    cmd = (sing_fs + 'freeview -v ' + to_see + ' ' + Ref_native_raw + ':colormap=heat:opacity=0.5:visible=1')
    run_cmd.do(cmd, diary_name)

