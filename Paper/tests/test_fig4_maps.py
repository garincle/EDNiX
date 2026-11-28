"""New Figure 4 (Paper/analysis/fig4_maps) end to end on synthetic run maps: real
profiling (profiles.compute_region_profile, percentiles P80/P20), real statistics and
drawing, no AFNI / wb_command. 10 datasets of the 7 paper species, 14 runs each."""
import os

import nibabel as nib
import numpy as np
import pandas as pd
import pytest

SPECIES = [("Dog", "Dog"), ("Macaque_BenHamed", "Macaque"), ("Macaque_ZhuGarin", "Macaque"),
           ("Rat", "Rat"), ("Mouse", "Mouse"), ("Marmoset_Tian", "Marmoset"),
           ("Marmoset_Hori", "Marmoset"), ("Mouselemur", "Mouselemur"),
           ("Human_Castrillon", "Human"), ("Human_Babayan", "Human")]
SHAPE = (8, 8, 16)
N_RUNS, N_SEL = 14, 5


def _build(tmp, monkeypatch):
    from Statistics.Group_fMRI.run_selection import profiles as sri
    from Statistics.Group_fMRI.run_selection import units as fu
    from Paper.analysis.fmri_fingerprints import fingerprints as ff
    from Paper.analysis.fmri_fingerprints import profile_stats as fa
    rng = np.random.default_rng(3)
    regions = list(ff.REGION_SHORT)[:8]
    # atlas: 16 labels (L/R x 8 regions), one slab of 8 x 8 x 1 voxels each
    atlas = np.zeros(SHAPE, np.int16)
    lines = ['0 0 0 0 0 0 0 "Clear Label"']
    for k, r in enumerate(regions):
        for h, side in enumerate(("L", "R")):
            lid = 2 * k + h + 1
            atlas[:, :, 2 * k + h] = lid
            lines.append(f'{lid} 255 0 0 1 1 1 "{side}_{r}"')
    aff = np.diag([2.0, 2.0, 2.0, 1.0])
    atlas_path = os.path.join(tmp, "atlas.nii.gz")
    nib.save(nib.Nifti1Image(atlas, aff), atlas_path)
    label_path = os.path.join(tmp, "labels.txt")
    with open(label_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    monkeypatch.setattr(fu, "find_matching_atlas_volume", lambda *a, **k: atlas_path)
    sri.set_profile_method("percentiles", 0.2, keep_sketch=False)
    sri.set_tail_reference(regions)
    value_col = sri.profile_value_col("percentiles", 0.2)
    common = sri.channel_names(regions, "percentiles", 0.2)
    label_df = sri.parse_itk_snap_label_file(label_path)

    base = rng.normal(0, 0.3, len(regions))
    sp_pat = {s: base + rng.normal(0, 0.15, len(regions)) for s in sorted({s for _, s in SPECIES})}
    vox_reg = (atlas - 1) // 2
    c = dict(cfg=dict(VALUE_COL=value_col, ZSCORE="scale", ATLAS_LEVEL=3,
                      REFERENCE_SPECIES="Human", SURFACE_COMMON_PCT=90.0),
             species_of_label={l: s for l, s in SPECIES}, common=common,
             divergence=dict(fa.DIVERGENCE_TIME_MYA), region_mask=np.ones(len(common), bool),
             profiles_by_label={}, members_by_label={}, selected_units={}, final_criterion={},
             group_map_profiles={"all": {}, "kept_final": {}}, bids_configs=[],
             ana_species_profiles={}, final_round=10)
    idx_rows = []
    for d, (lab, sp) in enumerate(SPECIES):
        gain = 0.6 + 0.8 * rng.random()                    # dataset amplitude (protocol)
        mdir = os.path.join(tmp, lab, "all")
        os.makedirs(os.path.join(mdir, "input_runs"))
        map_all = os.path.join(mdir, f"{lab}_ttest-stat_fisher_cmap.nii.gz")
        maps, profs, noise = [], {}, []
        for i in range(N_RUNS):
            sd = 0.05 + 0.25 * i / N_RUNS                  # later runs noisier
            vol = gain * sp_pat[sp][vox_reg] + rng.normal(0, sd, SHAPE)
            vol[atlas == 0] = 0
            sub = f"{d:02d}{i:02d}" if sp != "Dog" else f"dog{i}"
            ses = "1"
            p = os.path.join(mdir, "input_runs", f"{i + 1:02d}_sub-{sub}_ses-{ses}_run-1.nii.gz")
            nib.save(nib.Nifti1Image(vol.astype(np.float32), aff), p)
            maps.append(vol)
            key = (sub, ses if sp != "Dog" else "1", "1", "seed")
            profs[key] = sri.compute_region_profile(p, atlas_path, label_df, atlas_level=3)
            noise.append(sd)
        mean = np.mean(maps, axis=0).astype(np.float32)
        nib.save(nib.Nifti1Image(mean, aff), map_all)
        nib.save(nib.Nifti1Image((atlas > 0).astype(np.int8), aff),
                 os.path.join(mdir, "group_mask_rs.nii.gz"))
        keys = list(profs)
        sel = [keys[j] for j in np.argsort(noise)[:N_SEL]]
        c["profiles_by_label"][lab] = profs
        c["members_by_label"][lab] = {k: [k] for k in keys}
        c["selected_units"][lab] = {k: profs[k] for k in sel}
        c["final_criterion"][lab] = sri.compute_region_profile(map_all, atlas_path, label_df,
                                                               atlas_level=3)
        c["group_map_profiles"]["all"][lab] = c["final_criterion"][lab]
        c["bids_configs"].append(dict(label=lab, species=sp, bids_dir=tmp, label_path=label_path))
        idx_rows.append(dict(label=lab, group="all", found=True, map_path=map_all))
    c["all_units_index"] = pd.DataFrame(idx_rows)
    c["map_index_final"] = pd.DataFrame(columns=["label", "group", "found", "map_path"])
    return c


def test_fig4_maps_end_to_end(tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    from Paper.analysis.fig4_maps import run as f4run
    os.makedirs(tmp_path / "data")
    c = _build(str(tmp_path / "data"), monkeypatch)
    out = str(tmp_path / "Figure4_maps")
    res = f4run.run(c, out, pool="all", n_draws=4, n_perm=2, amp_draws=20, n_bins=3,
                    do_matrices=False, do_surfaces=False, title="synthetic")
    for f in ("Figure4.png", "Figure4.pdf", "Figure4_harmonisation.png",
              "Supp_Figure4_noise_bins.png", "Supp_Figure4_bins_profiles.png",
              "Figure4_stats.xlsx", "Figure4_stats.txt"):
        assert os.path.isfile(os.path.join(out, f)), f
    ph = res["phylogeny"]
    t = ph["table"]
    a = t[(t.analysis == "3D maps") & (t.group == "all runs")].set_index("measure")
    assert set(a.index) == {"K_raw", "K_resLOSO", "M_raw", "M_res", "M_resLOSO", "M_align"}
    assert a["p_species"].between(1 / 5040 - 1e-9, 1).all()
    # the all-runs profile rebuilt from the run maps is the stored all-runs map profile
    txt = open(os.path.join(out, "Figure4_stats.txt")).read()
    assert "r = 1.0000" in txt
    # every dataset has a sibling-free / sibling profile; 4b uses the all-runs maps
    assert res["profiles"]["pool"] == "all" and len(res["profiles"]["vecs"]) == len(SPECIES)
    am = res["amplitude"]
    assert np.isfinite(am["ratio"]["all"]) and np.isfinite(am["ratio"]["selected"])
    assert am["rat_sel"].shape == (20,)
    sheets = pd.read_excel(os.path.join(out, "Figure4_stats.xlsx"), sheet_name=None)
    assert {"phylogeny_3D", "noise_bins", "amplitude_seed", "amplitude_per_dataset"} <= set(sheets)


def test_selected_option(tmp_path, monkeypatch):
    from Paper.analysis.fig4_maps import analyses as fa4
    os.makedirs(tmp_path / "data")
    c = _build(str(tmp_path / "data"), monkeypatch)
    ctx = fa4.Context(c, n_bins=3)
    ph = fa4.phylogeny(ctx, 2)
    pr_all, pr_sel = fa4.profiles(ctx, ph, "all"), fa4.profiles(ctx, ph, "selected")
    lab = SPECIES[0][0]
    assert not np.allclose(pr_all["vecs"][lab], pr_sel["vecs"][lab])


def test_surface_maps_resolve_to_their_own_file(tmp_path, monkeypatch):
    """supp_surfaces must hand render_common_scale one distinct file per group (the
    renderer resolves a path to the first '*_oversampled' map of its folder)."""
    from Paper.analysis.fig4_maps import analyses as fa4
    from Paper.analysis.fig4_maps import figures as ff4
    from Statistics.Group_fMRI.run_selection import surfaces as srf
    os.makedirs(tmp_path / "data")
    c = _build(str(tmp_path / "data"), monkeypatch)
    ctx = fa4.Context(c, n_bins=3)
    ph = fa4.phylogeny(ctx, 1)
    seen = {}

    def fake_render(sri, out_dir, maps, species_of_label, pct=90.0):
        for lab, groups in maps.items():
            seen[lab] = {g: nib.load(srf._oversampled(p)).get_fdata() for g, p in groups.items()}
        return {}, {}
    monkeypatch.setattr(srf, "render_common_scale", fake_render)
    ff4.supp_surfaces(ctx, ph, {"P_bins": []}, str(tmp_path / "out"))
    lab = SPECIES[0][0]
    assert set(seen[lab]) == {"all", "selected", "bin1", "bin2", "bin3"}
    assert not np.allclose(seen[lab]["selected"], seen[lab]["all"])
    assert not np.allclose(seen[lab]["bin1"], seen[lab]["bin3"])
