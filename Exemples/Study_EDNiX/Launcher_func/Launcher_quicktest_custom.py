import time
# quick test launcher: prints pipeline-like step markers
for step, fn in ((1, '_1_fMRI_preTTT_in_fMRIspace'), (3, '_3_mask_fMRI')):
    print(f'##  Working on step {step} (function: {fn}).  ##', flush=True); time.sleep(3)
