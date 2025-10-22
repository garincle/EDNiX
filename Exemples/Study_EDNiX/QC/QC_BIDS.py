import os
import warnings
from QC_group import QC_func_raw
from QC_group import QC_func_matrix
from QC_group import QC_anat_raw
# Suppress specific warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Define path operations
opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists

# Define paths
BIDS_DIR = "/srv/projects/easymribrain/scratch/Rat/BIDS_Gd/"
OUTPUT_DIR = opj(BIDS_DIR, "group_qc_analysis_raw")
os.makedirs(OUTPUT_DIR, exist_ok=True)

qc_df = QC_anat_raw.load_qc_data(BIDS_DIR)
qc_df.to_csv(opj(OUTPUT_DIR, "qc_data_combined.csv"), index=False)
QC_anat_raw.create_comprehensive_report(qc_df, OUTPUT_DIR)

qc_df = QC_func_raw.load_qc_data(BIDS_DIR)
qc_df.to_csv(opj(OUTPUT_DIR, "qc_data_combined.csv"), index=False)
QC_func_raw.create_comprehensive_report(qc_df, OUTPUT_DIR)

qc_df = QC_func_matrix.load_qc_data(BIDS_DIR)
qc_df.to_csv(opj(OUTPUT_DIR, "qc_data_combined.csv"), index=False)
QC_func_matrix.create_comprehensive_report(qc_df, OUTPUT_DIR)