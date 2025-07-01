import random
from pathlib import Path
import os
import glob
import shutil
from os.path import join as opj, dirname as opd, basename, exists


# === CONFIGURATION ===
base_dir = "/srv/projects/easymribrain/data/MRI"
output_base = "/srv/projects/easymribrain/scratch/Training_trainSs_UNet"
train_img_dir = os.path.join(output_base, "anatT1")
train_mask_dir = os.path.join(output_base, "masks")
val_img_dir = os.path.join(output_base, "valid")
val_mask_dir = os.path.join(output_base, "valid_masks")

val_ratio = 0.2  # 20% for validation

# === CREATE OUTPUT DIRS ===
for d in [train_img_dir, train_mask_dir, val_img_dir, val_mask_dir]:
    os.makedirs(d, exist_ok=True)

# === FIND ANATOMICAL PAIRS ===
anat_pairs = []
for mask_path in glob.glob(os.path.join(base_dir, "**", "*brain_mask_final_QCrsp.nii.gz"), recursive=True):
    mask_path = Path(mask_path)
    # Anatomical image pattern in same subj/session directory
    anat_imgs = list(mask_path.parents[4].glob("**/*_acpc_test_QC_*.nii.gz"))
    if anat_imgs:
        anat_pairs.append(("anat", anat_imgs[0], mask_path))

# === FIND FUNCTIONAL PAIRS ===
func_pairs = []
for mask_path in glob.glob(os.path.join(base_dir, "**", "maskDilat_Allineate_in_func.nii.gz"), recursive=True):
    mask_path = Path(mask_path)
    func_img = mask_path.parent / "Mean_Image_test.nii.gz"
    if func_img.exists():
        func_pairs.append(("func", func_img, mask_path))

# === COMBINE AND SHUFFLE ===
all_pairs = anat_pairs + func_pairs
random.shuffle(all_pairs)

# === SPLIT ===
n_total = len(all_pairs)
n_val = int(n_total * val_ratio)
val_set = all_pairs[:n_val]
train_set = all_pairs[n_val:]

# Counters
anat_train = sum(1 for kind, *_ in train_set if kind == "anat")
func_train = sum(1 for kind, *_ in train_set if kind == "func")
anat_val   = sum(1 for kind, *_ in val_set   if kind == "anat")
func_val   = sum(1 for kind, *_ in val_set   if kind == "func")

# === COPY FUNCTION ===
def copy_and_rename(pairs, img_outdir, mask_outdir, prefix):
    for i, (kind, img_path, mask_path) in enumerate(pairs):
        tag = f"{prefix}_{kind}_{i:04d}"
        img_dst = os.path.join(img_outdir, f"{tag}.nii.gz")
        mask_dst = os.path.join(mask_outdir, f"{tag}_mask.nii.gz")
        shutil.copy(img_path, img_dst)
        shutil.copy(mask_path, mask_dst)

# === COPY FILES ===
copy_and_rename(train_set, train_img_dir, train_mask_dir, "train")
copy_and_rename(val_set, val_img_dir, val_mask_dir, "val")

# === SUMMARY ===
print(f"✅ Total samples processed: {n_total}")
print(f"   • Training: {len(train_set)} (anat: {anat_train}, func: {func_train})")
print(f"   • Validation: {len(val_set)} (anat: {anat_val}, func: {func_val})")



# --- Step 4: Build training command ---
afni_sif = '/srv/projects/easymribrain/code/EDNiX/Tool_library/Singularity/'
# --- Step 1: Define folders ---
folders = {
    '-trmsk': '/srv/projects/easymribrain/scratch/Training_trainSs_UNet/masks',
    '-out': '/srv/projects/easymribrain/scratch/Training_trainSs_UNet/result',
    '-vt1w': '/srv/projects/easymribrain/scratch/Training_trainSs_UNet/valid',
    '-vmsk': '/srv/projects/easymribrain/scratch/Training_trainSs_UNet/valid_masks',
    '-trt1w': '/srv/projects/easymribrain/scratch/Training_trainSs_UNet/anatT1'}
# Create folders
for path in folders.values():
    os.makedirs(path, exist_ok=True)

import subprocess
spgo = subprocess.getoutput
command = (
    'python3 ' + opj(opd(afni_sif), 'NHP-BrainExtraction', 'UNet_Model', 'trainSs_UNet.py') +
    ' -trt1w ' + folders['-trt1w'] +
    ' -trmsk ' + folders['-trmsk'] +
    ' -out ' + folders['-out'] +
    ' -vt1w ' + folders['-vt1w'] +
    ' -vmsk ' + folders['-vmsk'] +
    ' -init ' + opj(opd(afni_sif), 'NHP-BrainExtraction', 'UNet_Model', 'models', 'model-20_macaque-epoch'))


nl = spgo(command)
print("\n✅ Training command ready:")
print(command)