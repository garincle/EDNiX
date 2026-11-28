"""Figure 4 (new layout, 2026-11-27) and its supplements computed on the 3D group maps.

  4a  phylogenetic signal of the all-runs maps (P80/P20 profile, retrosplenial seed):
      K_mult and Mantel, exact P over the 7! species permutations
  4b  common profile + dataset profiles (all-runs maps by default, option: selected)
  4c  dataset x dataset similarity, raw and LOSO residual
  harmonisation: amplitude per dataset and within / between species amplitude ratio,
      seed maps and correlation matrices, with its size null
  supplements: selected vs random (3D maps and run-vector mean), noise bins (K,
      identifiability, trend), amplitude-aware phylogeny, surfaces per group

Entry point: run.py (python -m Paper.analysis.fig4_maps.run ROOT OUT_DIR); the full night
run for one invariance value is Paper/run_selection/run_final.py.
"""
