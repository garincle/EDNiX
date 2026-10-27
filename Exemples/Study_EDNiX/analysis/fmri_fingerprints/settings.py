"""
Statistics / figure settings of the fMRI fingerprint analyses (part 2), applied on top of
the saved selection (Statistics/Group_fMRI/run_selection, 07_state). Only keys of
run_selection.selection.FIGURE_KEYS: they never change the selection.

Shared by the paper (paper_stats.py) and by the exploration pipeline
(sensitivity_analysis/fmri_selection/Launcher_exploration.py). Values unchanged from
Launcher_exploration.BASE_CONFIG (V15, 2026-10-27).
"""
FIGURE_SETTINGS = dict(
    # phylogeny / F6: profiles of the KEPT group maps of the final round
    PHYLOGENY_USE_FINAL_ROUND=True,
    RUN_MULTILEVEL_F6=False,          # F6 at the other atlas levels (exploration)
    RUN_PHYLO_COMPARATIVE=True,       # Mantel / Blomberg's K / PGLS (F11-F13, exploration)
    # fingerprints (F1/F1c, F6ab, F15-F18, F20)
    RUN_FINGERPRINTS=True,
    REFERENCE_SPECIES="Human",        # reference of F6b / F15b (Figure 4a) / F18
    # print-size panels (89 / 120 / 183 mm, 5-7 pt, TrueType PDF, PNG 600 dpi)
    RUN_PAPER_FIGURES=True,
    BRAIN_MASS_XLSX="/scratch2/EDNiX/Herculano-Houzel_mammalian_brain_dataset_complete.xlsx",
    # Figure 4b: pairs compared. "other_dataset" = same species in ANOTHER dataset vs
    # different species (pairs of one dataset share protocol, scanner and often the
    # animal: left out); "all" = every same-species pair
    WITHIN_ACROSS_PAIRS="other_dataset",
    # Figure 4b, test shown: "wilcoxon" (rank-sum on the pairs, one-sided, optimistic
    # because pairs share runs; Clément 2026-10-26: kept) or "permutation" (species
    # permuted across datasets); both are reported
    WITHIN_ACROSS_TEST="wilcoxon",
    # 6d (exploration): one panel per regional index, Human + Macaque vs the others
    FIG6D=dict(
        regions=("posterior_parietal_cortex", "prefrontal_cortex_(pfc)",
                 "somatosensory_cortex", "visual_striate_cortex", "striatum",
                 "motor_and_premotor"),
        group=("Human", "Macaque")),
    # Figure 5 (surface maps): relative brain size per species, datasets per row
    SURFACE_ZOOM={"Human": 1.0, "Macaque": 1.0, "Dog": 1.0, "Marmoset": 0.9,
                  "Mouselemur": 0.9, "Rat": 0.9, "Mouse": 0.9},
    SURFACE_DATASETS_PER_ROW=2,
    # run-level r with the kept map and with the selection criterion (Supp 9b, Supp 12)
    RUN_DISTANCE_TO_REF=True,
)

# Control run (one common target, LEAVE_ONE_OUT=False): no AFNI kept maps
CONTROL_FIGURE_SETTINGS = dict(FIGURE_SETTINGS, PHYLOGENY_USE_FINAL_ROUND=False,
                               RUN_DISTANCE_TO_REF=False)
