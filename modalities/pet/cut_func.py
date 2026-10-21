import os
import json
import numpy as np

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile


from Tools import filename
from pet import run_cut

def cut(func_dir,func_raw_dir,functype,template,petfile,petjson,runlist,atlas,atlas_label,diary_name,sing_afni, sing_wb):


    func_data_dir = opj(func_dir,'residuals',template)

    pet_dir = opd(petfile)

    f= open(opj(func_dir,'Settings.json'),'r')
    info  = json.load(f)
    T1_eq = int(info["Remove the first n volumes "])
    Subname = filename.get_subname(petfile)


    for i in range(len(runlist)):
        funcjson = opj(func_raw_dir,'_'.join([Subname,'task-rest','run-0' + str(runlist[i]),'bold.json' ]))
        runname = '_'.join([Subname,'task-rest','run-0' + str(runlist[i]), 'space-' + template,functype,'bold.nii.gz' ])
        funcfile = opj(func_data_dir,runname)
        if opi(funcfile):
            print(funcfile)
            start,stop = run_cut.getframes(funcfile,T1_eq,funcjson,petjson)
            run_cut.cut(opj(pet_dir,'templates',template),funcfile, petfile,atlas,start,stop,atlas_label, diary_name, sing_afni, sing_wb)