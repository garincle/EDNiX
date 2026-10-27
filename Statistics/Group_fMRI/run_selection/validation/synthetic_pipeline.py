"""Synthetic end-to-end run of the Fig5 pipeline (selection + statistics + figures), no AFNI.

    python synthetic_pipeline.py <repo_root> <out_dir> <old|new>

Runs the pipeline of the checkout <repo_root> on synthetic seed maps (10 datasets of the
paper, 14 regions x 2 hemispheres, ~20 runs each): compute_all_run_profiles is replaced
by a generator that writes small NIfTI maps and profiles them with the REAL
compute_region_profile; the AFNI steps are off. The data do not depend on the checkout,
so two checkouts can be compared:

    git worktree add /tmp/ednix_old <old commit>
    python synthetic_pipeline.py /tmp/ednix_old /tmp/out_old old    # 'old' = before the move
    python synthetic_pipeline.py <this repo>    /tmp/out_new new
    python ../compare_states.py /tmp/out_old/Fig_5_synth /tmp/out_new/Fig_5_synth
    python compare_outputs.py   /tmp/out_old/Fig_5_synth /tmp/out_new/Fig_5_synth

2026-10-27, 0707c94 (before) vs the run_selection package: same state, 58/58 tables and
texts identical.
"""
import os
import sys
import shutil


def main():
    """Run from the command line only (the test suite imports every module)."""
    repo, out_root, which = sys.argv[1], sys.argv[2], sys.argv[3]
    sys.path.insert(0, repo)
    os.environ["MPLBACKEND"] = "Agg"

    import numpy as np
    import pandas as pd
    import nibabel as nib

    if which == "old":
        from Exemples.Study_EDNiX.analysis.Fig5 import SBA_Region_Profile_Index as sri
    else:
        from Statistics.Group_fMRI.run_selection import profiles as sri
    try:
        from Exemples.Study_EDNiX.analysis.sensitivity_analysis.fmri_selection import exploration_pipeline as pipe
    except ImportError:                        # checkouts before 2026-10-27
        from Exemples.Study_EDNiX.analysis.Fig5 import fig5_pipeline as pipe

    REGIONS = ["prefrontal_cortex_(pfc)", "motor_and_premotor", "somatosensory_cortex",
               "posterior_parietal_cortex", "visual_striate_cortex", "extrastriate_cortex",
               "auditory_cortex", "insula", "cingulate_cortex", "retrosplenial_cortex",
               "striatum", "thalamus", "hippocampus", "amygdala"]
    DATASETS = [("Dog", "Dog", "BIDS_Boch_K9", 8), ("Macaque_BenHamed", "Macaque", "BIDS_BenHamed", 9),
                ("Macaque_ZhuGarin", "Macaque", "BIDS_Zhu_Garin", 8), ("Rat", "Rat", "BIDS_Grandjean", 8),
                ("Mouse", "Mouse", "BIDS_Grandjean2", 10), ("Marmoset_Tian", "Marmoset", "BIDS_Tian", 7),
                ("Marmoset_Hori", "Marmoset", "BIDS_Hori", 8), ("Mouselemur", "Mouselemur", "BIDS_Garin", 7),
                ("Human_Castrillon", "Human", "BIDS_Castrillon", 9), ("Human_Babayan", "Human", "BIDS_Babayan", 8)]
    SHAPE = (12, 12, 14)

    # fixed absolute paths (they are written in the outputs): same for every checkout
    import tempfile
    SCRATCH = os.path.join(tempfile.gettempdir(), "ednix_synthetic_pipeline")
    work = os.path.join(SCRATCH, "_inputs")
    os.makedirs(work, exist_ok=True)
    label_path = os.path.join(work, "labels.txt")
    with open(label_path, "w") as f:
        f.write("# synthetic ITK-SNAP labels\n")
        lid = 1
        for r in REGIONS:
            for h in ("L", "R"):
                f.write(f'{lid} 10 10 10 1 1 1 "{h}_{r}"\n')
                lid += 1
    # atlas: 28 labels, each a block of voxels
    atlas = np.zeros(SHAPE, dtype=np.int16)
    flat = atlas.reshape(-1)
    n_lab = 2 * len(REGIONS)
    per = flat.size // n_lab
    for i in range(n_lab):
        flat[i * per:(i + 1) * per] = i + 1
    atlas_path = os.path.join(work, "atlas.nii.gz")
    nib.save(nib.Nifti1Image(atlas, np.eye(4)), atlas_path)

    master = np.random.default_rng(12345)
    base = master.normal(0, 0.25, len(REGIONS))              # shared mammalian pattern
    species_dev = {sp: master.normal(0, 0.12, len(REGIONS)) for sp in
                   sorted({d[1] for d in DATASETS})}
    dataset_dev = {d[0]: master.normal(0, 0.06, len(REGIONS)) for d in DATASETS}
    RUNS = {}
    for lab, sp, bids, n_sub in DATASETS:
        rng = np.random.default_rng([d[0] for d in DATASETS].index(lab) + 1)
        runs = []
        for s in range(n_sub):
            animal_dev = rng.normal(0, 0.05, len(REGIONS))
            n_ses = 2 if s % 3 == 0 else 1
            for ses in range(1, n_ses + 1):
                for run in range(1, 3):
                    good = rng.random() > 0.3
                    noise = rng.normal(0, 0.08 if good else 0.3, len(REGIONS))
                    pattern = base + species_dev[sp] + dataset_dev[lab] + animal_dev + noise
                    runs.append((f"s{s:02d}", str(ses), str(run), pattern, rng.integers(0, 2 ** 31)))
        RUNS[lab] = runs


    def fake_compute_all_run_profiles(species, bids_dir, label_df, classification_path,
                                      category="All", seeds=None, atlas_level=None, **kw):
        lab = os.path.basename(bids_dir)
        profiles = {}
        mdir = os.path.join(work, lab)
        os.makedirs(mdir, exist_ok=True)
        for subj, ses, run, pattern, seed in RUNS[lab]:
            rng = np.random.default_rng(int(seed))
            vox = np.zeros(SHAPE)
            fl = vox.reshape(-1)
            for i in range(n_lab):
                fl[i * per:(i + 1) * per] = np.clip(
                    pattern[i // 2] + rng.normal(0, 0.25, per), -0.99, 0.99)
            p = os.path.join(mdir, f"sub-{subj}_ses-{ses}_run-{run}.nii.gz")
            if not os.path.isfile(p):
                nib.save(nib.Nifti1Image(vox.astype(np.float32), np.eye(4)), p)
            prof = sri.compute_region_profile(p, atlas_path, label_df, atlas_level=None)
            profiles[(subj, ses, run, seeds[0])] = prof
        return profiles, pd.DataFrame()


    sri.compute_all_run_profiles = fake_compute_all_run_profiles

    bids_configs = [dict(label=lab, species=sp, bids_dir=f"/synthetic/{lab}", label_path=label_path,
                         seeds=["BA_26_29_30__retrosplenial"],
                         state=("isoflurane" if lab in ("Rat", "Mouselemur", "Macaque_ZhuGarin")
                                else "awake"),
                         reference_subject=None, reference_session=None, reference_run=None)
                    for lab, sp, bids, n in DATASETS]

    qc_rows = []
    qrng = np.random.default_rng(7)
    for lab, sp, bids, n in DATASETS:
        for subj, ses, run, *_ in RUNS[lab]:
            qc_rows.append(dict(label=lab, bids_dir=bids, subject=subj, session=ses, run=run,
                                func_TSNR=qrng.normal(40, 8), func_mean_fd=abs(qrng.normal(.1, .04)),
                                func_gcor=abs(qrng.normal(.05, .02))))
    qc = pd.DataFrame(qc_rows)
    df_qc_by_label = {lab: g.drop(columns="label").copy() for lab, g in qc.groupby("label")}

    root = os.path.join(SCRATCH, "Fig_5_synth")
    if os.path.isdir(root):
        shutil.rmtree(root)
    config = dict(
        OUTPUT_ROOT=root, PROFILE_METHOD="percentiles", TAIL_Q=0.2, UNIT="run", ATLAS_LEVEL=3,
        RUN_QC_ENABLE=False, N_SELECT=6, INVARIANCE_MODE="top_fraction", INVARIANCE_VALUE=0.8,
        SIMILARITY="correlation", ZSCORE="scale", LEAVE_ONE_OUT=True, RUN_ROUND2_SELECTION=True,
        RUN_ITERATIVE_SELECTION=True, N_ITERATIVE_ROUNDS=5,
        SEED_REGION_EXCLUDE=["retrosplenial_cortex"], SEED_REGION_EXCLUDE_SELECTION=False,
        SEED_REGION_EXCLUDE_ANALYSES=True,
        RUN_ALL_UNITS_MAPS=False, RUN_GROUP_MAPS=False, RENDER_SURFACE=False, COPY_JPG=False,
        TREE_PATH=None, RUN_PHYLO_COMPARATIVE=False, RUN_MULTILEVEL_F6=False,
        RUN_FINGERPRINTS=True, REFERENCE_SPECIES="Human", RUN_PAPER_FIGURES=True,
        BRAIN_MASS_G={"Human": 1350, "Macaque": 90, "Marmoset": 8, "Mouselemur": 1.7,
                      "Rat": 2, "Mouse": 0.45, "Dog": 80},
        WITHIN_ACROSS_PAIRS="other_dataset", WITHIN_ACROSS_TEST="wilcoxon",
        RUN_DISTANCE_TO_REF=True, DISTANCE_REF_GROUP="kept", N_PERM_PHYLOGENY=500,
        N_DRAWS_BALANCED=50, PHYLOGENY_USE_FINAL_ROUND=False,
    )
    part = sys.argv[4] if len(sys.argv) > 4 else "both"
    pipe.run(sri, bids_configs, config=config, df_qc_by_label=df_qc_by_label, part=part)
    dst = os.path.join(out_root, "Fig_5_synth")
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    shutil.move(root, dst)
    print("SYNTH DONE", dst)


if __name__ == "__main__":
    main()
