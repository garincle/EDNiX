import os
import glob
import subprocess
import pandas as pd
import re
import os
import shutil
import os
import glob
import subprocess
import re
# Precompiled binary linux_ubuntu_16_64: May  4 2020 (Version AFNI_20.1.06 'Otho')
# Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
spco = subprocess.check_output
opd = os.path.dirname

sub_id = 'Litchi'
for ses_id in ['1', '2']:
    ############# variables
    out_analysis = '/scratch2/EDNiX/Macaque/Test' + ses_id + '_Litchy'
    DICOMdir = '/home/common/benhalab/CASCAD/Macaque/MRI/Sociomonk/Test' + ses_id + '_Litchy/sub-Litchi_ses-0' + ses_id + '/dicom'

    in_DICOM = DICOMdir
    out_code = out_analysis + '/code'
    config_file = out_code + '/dcm2bids_config.json'
    config_file_orig = '/home/common/benhalab/CASCAD/Macaque/MRI/Sociomonk/Test1_Litchy/code/dcm2bids_config.json'

    if not os.path.exists(out_analysis): os.mkdir((out_analysis))
    if not os.path.exists(in_DICOM): os.mkdir((in_DICOM))
    # Copies only the file contents (no metadata)

    ##creat  your scaffold
    spco(['dcm2bids_scaffold', '-o', out_analysis, '--force'])

    ##creat configuration file
    spco(['dcm2bids_helper', '-d', in_DICOM, '-o', out_analysis + '/sourcedata/', '--force'])

    out_analysis = '/scratch2/EDNiX/Macaque/Test' + ses_id + '_Litchy/Bids'

    shutil.copyfile(config_file_orig, config_file)
    subprocess.check_output([
        'dcm2bids', '-d', DICOMdir, '-p', sub_id, '-s', ses_id,
        '-c', config_file, '-o', out_analysis,
        '--auto_extract_entities', '--bids_validate', '--force_dcm2bids'])
