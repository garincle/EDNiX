# run_selection -- cross-species consensus selection of seed-based runs

Selects, in every dataset, the units (runs, sessions or subjects) whose regional
seed-connectivity profile is closest to the profile shared by the other species,
then builds the AFNI group maps of the kept and not-kept units. It needs AFNI,
Connectome Workbench (surface renders) and EDNiX-preprocessed BIDS datasets.

Part 1 of the former `Exemples/Study_EDNiX/analysis/Fig5` pipeline, moved here unchanged
on 2026-10-27 so that it can be rerun on its own (for example with new species). Every
moved function and constant is AST-identical to its original. The only change is where
the modules are imported from.

## Run it

```
python Statistics/Group_fMRI/run_selection/launcher.py
```

- Datasets and settings are in `config_ednix.py`.
- The output goes under `SELECTION_CONFIG["OUTPUT_ROOT"]`, in the tree below.
- To add a species, add one entry to `BIDS_CONFIGS` and its tip in the time tree to
  `selection.CONFIG["SPECIES_TO_TREE_TIP"]`.

## Same result as before?

Each new run moves the previous `07_state` to `_previous_runs/<date>/`. Then run:

```
python Statistics/Group_fMRI/run_selection/compare_states.py
```

It compares the new state with the most recent archived one, entry by entry: profiles,
common regions, scores of every round, kept units and group-map profiles. Exit code 0
means the selection is the same.

## Steps (selection._run_selection)

1. Profile of every run: `profiles.compute_all_run_profiles`.
   - Region values of the seed map at `ATLAS_LEVEL`.
   - `PROFILE_METHOD` is either "percentiles" (P80 and P20 per region, `TAIL_Q=0.2`) or
     "median_r".
2. Optional run QC (`RUN_QC_ENABLE`), then aggregation to `UNIT` (`units.py`).
3. STEP 1: AFNI 3dttest++ group map of all units of each dataset (`group_maps.py`).
4. Common regions:
   - regions present at the atlas level in every dataset's atlas and sampled by every
     group map;
   - the region that contains the seed is found and flagged (`SEED_REGION_EXCLUDE_*`).
5. F22: reliability of the three profile constructions, across runs of the same animal
   (`reliability.py`).
6. Region invariance across species. The `INVARIANCE_VALUE` most invariant regions are
   used by the score (`consensus.py`).
7. Round 1:
   - each dataset is ranked against the consensus of the other species
     (`LEAVE_ONE_OUT`), by Pearson r on the scaled profiles;
   - the `N_SELECT` best units are kept.
8. Round 2: re-ranking against the dataset's own kept units.
9. Rounds 3 to `N_ITERATIVE_ROUNDS`: the cross-species target is rebuilt from the kept
   units at every round; the last round is the final selection.
10. STEP 3 / 3b:
    - AFNI group maps of the kept and not-kept units, size-matched;
    - jpg previews and surface renders on one colour scale (`surfaces.py`).
11. `07_state/selection_state.pkl` (with a readable `.txt` beside it): everything the
    statistics need, without NIfTI.

## Files

| file | content |
|---|---|
| `launcher.py` | runs the main selection and the common-target control |
| `config_ednix.py` | datasets, profile method, selection settings, QC table |
| `selection.py` | `CONFIG`, steps, `_run_selection`, `run_selection`, `save_state` / `load_state` |
| `consensus.py` | consensus, iterative selection, invariance, kept vs not-kept balance, QC validation |
| `profiles.py` | regional profiles of seed maps (former `SBA_Region_Profile_Index.py`) |
| `group_maps.py` | AFNI group seed maps (former `SBA_group_maps.py`) |
| `outlier_qc.py` | AFNI helpers for kept / outlier maps (former `_Group_anal_Outlier_QC.py`) |
| `units.py` | run to unit aggregation, run QC, group-map profiles (former `fig5_units.py`) |
| `reliability.py` | arousal state of each run, F22 profile reliability |
| `surfaces.py` | surface renders of all / kept / not-kept maps on one scale |
| `compare_states.py` | is a new selection identical to an earlier one? |

## Output tree (`OUTPUT_ROOT/`)

- `01_all_units/`: STEP 1 maps
- `02_selection/`: selection, common regions, invariance, unit counts, balance, F22
- `03_group_maps/`, `03b_group_maps_final/`: kept and not-kept maps, rounds 1 and final
- `04_jpg/`, `05_surfaces/` (and their final-round versions `04b_jpg_final/`,
  `05b_surfaces_final/`): previews and surface renders
- `07_state/`: state read by the statistics (part 2)
- `methods.txt`, `00_log.txt`
