"""
Run selection -- launch part 1 for the EDNiX datasets (config_ednix.py).

  RUN_MAIN     main selection: leave-one-species-out consensus, AFNI 3dttest++ group
               maps of all / kept / not-kept units, surface renders
               -> <OUTPUT_ROOT>/07_state/selection_state.pkl
  RUN_CONTROL  same selection towards ONE common target (LEAVE_ONE_OUT=False), Python
               only, regions of the main run -> <OUTPUT_ROOT>_control_common_target/07_state

Requires AFNI and Connectome Workbench (singularity images, see Tools.Load_EDNiX_requirement).
The paper statistics and figures are run afterwards from Exemples/Study_EDNiX/analysis
(they only read 07_state). To check that a new run gives the same selection as an
earlier one: compare_states.py <old OUTPUT_ROOT> <new OUTPUT_ROOT>.
"""
import sys

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Statistics.Group_fMRI.run_selection import profiles as sri
from Statistics.Group_fMRI.run_selection import config_ednix as cfg
from Statistics.Group_fMRI.run_selection.selection import run_selection

RUN_MAIN = True
RUN_CONTROL = True

if __name__ == "__main__":
    df_qc_by_label = cfg.load_qc_by_label()
    if RUN_MAIN:
        run_selection(sri, cfg.BIDS_CONFIGS, config=cfg.SELECTION_CONFIG,
                      df_qc_by_label=df_qc_by_label)
    if RUN_CONTROL:
        regions = cfg.read_common_regions(cfg.SELECTION_CONFIG["OUTPUT_ROOT"])
        if regions is None:
            print("  [CONTROL] no common_regions.xlsx from the main run -- control skipped "
                  "(run RUN_MAIN first)")
        else:
            run_selection(sri, cfg.BIDS_CONFIGS,
                          config=dict(cfg.CONTROL_SELECTION_CONFIG,
                                      COMMON_REGIONS_OVERRIDE=regions),
                          df_qc_by_label=df_qc_by_label)
