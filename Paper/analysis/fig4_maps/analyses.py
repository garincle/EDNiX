"""Analyses of the new Figure 4 and its supplements, on one finished selection.

    ctx = Context(c)                    # run table, run maps, scores, bins
    ph  = phylogeny(ctx, n_draws)       # Fig. 4a + Supp: all runs / selected / random
    pr  = profiles(ctx, ph, pool)       # Fig. 4b, 4c: dataset profiles of the 3D maps
    nb  = noise_bins(ctx, ph, n_perm)   # Supp: K and identifiability per noise bin
    am  = amplitude(ctx, n_draws)       # harmonisation: seed-map regional means
    mx  = matrices(ctx, n_draws)        # harmonisation + phylogeny: correlation matrices

Every P value is either exact over the 7! species permutations (P_species) or a
size control: (1 + #null >= observed) / (n + 1) against random draws of as many runs
per dataset as selected (selected) or against random partitions of each dataset into
bins of the same sizes (noise bins). Seeds are fixed, the analysis is deterministic.
"""
import os
import time

import numpy as np
import pandas as pd

from Paper.analysis.fig4_maps import data as fd
from Paper.analysis.fig4_maps import stats as fs

SEED = 0
PATTERN = ("K_raw", "K_resLOSO", "M_raw", "M_res", "M_resLOSO", "M_align")
AMPLITUDE = ("K_raw", "K_resSub", "K_resReg", "K_amp", "M_euc", "M_eucRes", "M_corr")


class Log:
    def __init__(self):
        self.lines = []

    def __call__(self, s=""):
        print(s, flush=True)
        self.lines.append(str(s))

    def text(self):
        return "\n".join(self.lines) + "\n"


class Context:
    """Everything the analyses share, built once from paper_stats.Inputs(...).fmri()."""

    def __init__(self, c, n_bins=5, log=None):
        from Paper.analysis.fmri_fingerprints import phylo
        from Statistics.Group_fMRI.run_selection import surfaces as srf
        self.c, self.log, self.B = c, log or Log(), int(n_bins)
        self.rt = rt = fd.RunTable(c)
        self.labels, self.sp_of_ds, self.sets = rt.labels, rt.sp_of_ds, rt.sets
        self.sps = phylo.order_species(sorted(set(rt.sp_of_ds)))
        self.Dp, self.Cphy = phylo.phylo_matrices(self.sps)
        self.T = fs.Tests(self.Dp, self.Cphy)
        self.maps = srf.map_paths_by_dataset(c["all_units_index"], c["map_index_final"])
        self.score, n_mask = fd.scores(c, rt)
        L = self.log
        L(f"  species ({len(self.sps)}): {', '.join(self.sps)}; {rt.vecs.shape[1]} profile values; "
          f"{int(rt.selected.sum())} selected / {len(rt.selected)} runs; noise score = r with the "
          f"round-{c.get('final_round', '?')} criterion on {n_mask}/{len(c['common'])} values")
        self.dm, self.mapkey, self.idx_all, self.idx_sel, self.ranked = {}, {}, {}, {}, {}
        t0 = time.time()
        for d in self.sets:
            lab = self.labels[d]
            if "all" not in self.maps.get(lab, {}):
                raise RuntimeError(f"{lab}: no all-runs group map in the selection state")
            dm = self.dm[d] = fd.DatasetMaps(c, lab, self.maps[lab]["all"])
            rows = rt.rows(d)
            m2v, miss, extra, amb = fd.pair_runs({i: rt.keys[i] for i in rows}, dm.paths)
            self.mapkey.update(m2v)
            have = np.array([i in m2v for i in rows])
            self.idx_all[d] = rows[have]
            self.idx_sel[d] = rows[have & rt.selected[rows]]
            n_sel = int(rt.selected[rows].sum())
            L(f"  {lab}: {len(dm.paths)} run maps, {len(rows)} run vectors, {len(m2v)} paired "
              f"({len(miss)} without map, {len(amb)} ambiguous, {len(extra)} maps without vector), "
              f"{len(self.idx_sel[d])}/{n_sel} selected paired; {len(dm.vox)} voxels "
              f"({dm.mask_src})")
            if len(self.idx_sel[d]) < n_sel or len(self.idx_sel[d]) < 2:
                raise RuntimeError(
                    f"{lab}: {n_sel - len(self.idx_sel[d])} of {n_sel} selected run(s) have no run "
                    f"map in {os.path.dirname(self.maps[lab]['all'])}/input_runs/. Keys -- run "
                    f"vectors: {[rt.keys[i] for i in miss[:3]]}; maps: {(extra or list(dm.paths))[:3]}")
            ok = self.idx_all[d][np.isfinite(self.score[self.idx_all[d]])]
            if len(ok) < 2 * self.B:
                raise RuntimeError(f"{lab}: {len(ok)} scored runs, fewer than 2 per bin (B = {self.B})")
            self.ranked[d] = ok[np.argsort(-self.score[ok], kind="stable")]   # closest first
        self.bins = {d: np.array_split(self.ranked[d], self.B) for d in self.sets}
        self.sizes = {d: [len(x) for x in self.bins[d]] for d in self.sets}
        L(f"  run maps indexed in {time.time() - t0:.0f} s")

    # profiles of groups of runs ------------------------------------------------------
    def keys_of(self, d, rows):
        return list(self.dm[d].paths) if rows is None else [self.mapkey[i] for i in rows]

    def prof3d(self, rows_by_ds):
        """{d: z-scored P80/P20 profile of the mean map of the rows}; rows None = every
        run map of input_runs/ (= the all-runs figure map)."""
        return {d: fs.z(self.dm[d].profile(self.keys_of(d, r))) for d, r in rows_by_ds.items()}

    def profvec(self, rows_by_ds):
        """Same runs, profile = z(mean of the per-run vectors) (former Fig. 4d)."""
        return {d: fs.z(np.nanmean(self.rt.vecs[r], axis=0)) for d, r in rows_by_ds.items()}

    def draw(self, rng):
        return {d: rng.choice(self.idx_all[d], len(self.idx_sel[d]), replace=False)
                for d in self.sets}

    def partition(self, rng):
        return {d: np.split(rng.permutation(self.ranked[d]), np.cumsum(self.sizes[d])[:-1])
                for d in self.sets}


def _stored(c, group, lab):
    from Paper.analysis.fmri_fingerprints import profile_stats as fa
    v = c["group_map_profiles"].get(group, {}).get(lab)
    if v is None:
        return None
    return fa._ordered_common_values(v[1] if isinstance(v, tuple) else v, c["common"],
                                     c["cfg"]["VALUE_COL"])


def _row(analysis, group, measure, value, p_species=np.nan, null=None, side="high"):
    r = dict(analysis=analysis, group=group, measure=measure, value=value, p_species=p_species)
    if null is not None and len(null):
        v = np.asarray(null, float)
        r.update(null_median=float(np.median(v)), null_p5=float(np.percentile(v, 5)),
                 null_p95=float(np.percentile(v, 95)),
                 p_size=fs.p_high(v, value) if side == "high" else fs.p_low(v, value),
                 p_size_side=side)
    return r


# ── Fig. 4a: phylogenetic signal of the 3D maps ─────────────────────────────────────────
def phylogeny(ctx, n_draws=100):
    L, T, sp, sps = ctx.log, ctx.T, ctx.sp_of_ds, ctx.sps
    L("\n== phylogenetic signal of the 3D group maps (P80/P20, all runs / selected / random)")
    L("  check: profile rebuilt from the run maps vs stored profile of the figure map")
    for d in ctx.sets:
        lab = ctx.labels[d]
        for g, rows in (("all", None), ("kept_final", ctx.idx_sel[d])):
            a, b = ctx.dm[d].profile(ctx.keys_of(d, rows)), _stored(ctx.c, g, lab)
            if b is None:
                L(f"    {lab:<20}{g:<11} no stored profile")
                continue
            ok = np.isfinite(a) & np.isfinite(b)
            L(f"    {lab:<20}{g:<11} r = {fs.corr(a, b):.4f}, max |diff| = "
              f"{np.max(np.abs(a[ok] - b[ok])):.2e}")
    P_all = ctx.prof3d({d: None for d in ctx.sets})
    Cg, Cs = fs.commons(P_all, sp, sps)
    P_sel = ctx.prof3d(ctx.idx_sel)
    obs = {"all runs": T.pattern(fs.species_profiles(P_all, sp, sps, Cg, Cs), keep_null=True),
           "selected": T.pattern(fs.species_profiles(P_sel, sp, sps, Cg, Cs))}
    vec = {"all runs": T.pattern(fs.species_profiles(ctx.profvec(ctx.idx_all), sp, sps, Cg, Cs)),
           "selected": T.pattern(fs.species_profiles(ctx.profvec(ctx.idx_sel), sp, sps, Cg, Cs))}
    rng = np.random.default_rng(SEED + 2)
    draws, draws_vec, draw_profiles, t0 = [], [], [], time.time()
    for i in range(n_draws):
        pick = ctx.draw(rng)
        P = ctx.prof3d(pick)
        draw_profiles.append(P)
        draws.append(T.pattern(fs.species_profiles(P, sp, sps, Cg, Cs), with_p=False))
        draws_vec.append(T.pattern(fs.species_profiles(ctx.profvec(pick), sp, sps, Cg, Cs),
                                   with_p=False))
        if (i + 1) % 20 == 0:
            print(f"    draw {i + 1}/{n_draws} ({(time.time() - t0) / (i + 1):.1f} s each)", flush=True)
    rows = []
    L(f"\n  {'measure':<11}{'all runs (P)':>18}{'selected (P)':>18}{'random med [5-95]':>24}"
      f"{'sel vs rand P':>15}")
    for k in PATTERN:
        a, s = obs["all runs"][k], obs["selected"][k]
        rr = np.array([x[k][0] for x in draws])
        L(f"  {k:<11}{a[0]:>+9.3f} ({a[1]:.4f}){s[0]:>+11.3f} ({s[1]:.4f}){np.median(rr):>+10.3f} "
          f"[{np.percentile(rr, 5):+.3f}, {np.percentile(rr, 95):+.3f}]{fs.p_high(rr, s[0]):>11.3f}")
        rows.append(_row("3D maps", "all runs", k, a[0], a[1]))
        rows.append(_row("3D maps", "selected", k, s[0], s[1], rr))
        rv = np.array([x[k][0] for x in draws_vec])
        rows.append(_row("run-vector mean", "all runs", k, vec["all runs"][k][0], vec["all runs"][k][1]))
        rows.append(_row("run-vector mean", "selected", k, vec["selected"][k][0],
                         vec["selected"][k][1], rv))
    L("  (run-vector mean = the former Fig. 4d aggregate, same runs; in the table only)")
    return dict(P_all=P_all, P_sel=P_sel, Cg=Cg, Cs=Cs, obs=obs, draws=draws,
                draw_profiles=draw_profiles, table=pd.DataFrame(rows))


# ── Fig. 4b / 4c: dataset profiles and similarity ───────────────────────────────────────
def profiles(ctx, ph, pool="all"):
    """Dataset vectors for Fig. 4b (common + per-dataset radars) and 4c (similarity):
    z-scored P80/P20 profile of the all-runs (pool='all', default, as 4a) or selected
    (pool='selected') 3D map of every dataset, keyed by dataset label."""
    from Paper.analysis.fmri_fingerprints import fingerprints as ff
    P = ph["P_all"] if pool == "all" else ph["P_sel"]
    vecs = {ctx.labels[d]: P[d] for d in ctx.sets}
    sol = ctx.c["species_of_label"]
    dist = ff.dataset_distances(vecs, sol, ctx.c["divergence"], ctx.c["cfg"]["REFERENCE_SPECIES"])
    order = ff.order_labels(vecs, sol, ctx.c["divergence"])
    D = ff.distance_matrix(vecs, order)
    res = {ctx.labels[d]: fs.resid(P[d], ph["Cs"][ctx.sp_of_ds[d]]) for d in ctx.sets}
    D_res = ff.distance_matrix(res, order)
    Dsp = fs.dist_corr(fs.species_means(P, ctx.sp_of_ds, ctx.sps))
    Dsp_res = fs.dist_corr(fs.species_profiles(P, ctx.sp_of_ds, ctx.sps, ph["Cg"], ph["Cs"])["resLOSO"])
    ctx.log(f"\n== Fig. 4b/c dataset profiles from the {pool}-runs 3D maps")
    for r in dist.itertuples():
        ctx.log(f"  {r.label:<20}r with LOSO consensus {r.r_consensus:+.3f}")
    return dict(pool=pool, vecs=vecs, dist=dist, D=D, D_res=D_res, Dsp=Dsp, Dsp_res=Dsp_res)


# ── Supp: noise bins ────────────────────────────────────────────────────────────────────
def noise_bins(ctx, ph, n_perm=100):
    """K / Mantel and species identifiability of the mean maps of B bins of runs ranked
    by their r with the selection criterion (bin 1 = closest, holds the selected runs).
    Null: random partitions of every dataset into B bins of the same sizes (bin labels
    permuted within dataset): per-bin size control and the trend over bins."""
    L, T, sp, sps, B = ctx.log, ctx.T, ctx.sp_of_ds, ctx.sps, ctx.B
    Cg, Cs = ph["Cg"], ph["Cs"]
    L(f"\n== noise bins (B = {B}, {n_perm} random partitions; residuals vs the general "
      f"LOSO common profile of the all-runs maps)")
    for b in range(B):
        L(f"  bin {b + 1}: " + ", ".join(
            f"{ctx.labels[d]} {np.median(ctx.score[ctx.bins[d][b]]):+.2f} "
            f"({int(ctx.rt.selected[ctx.bins[d][b]].sum())} sel)" for d in ctx.sets))
    P_bins = [ctx.prof3d({d: ctx.bins[d][b] for d in ctx.sets}) for b in range(B)]
    obs = [T.pattern(fs.species_profiles(P, sp, sps, Cg, Cs)) for P in P_bins]
    ident = {"all": fs.identifiability(ph["P_all"], sp, sps, Cs),
             "selected": fs.identifiability(ph["P_sel"], sp, sps, Cs)}
    for b in range(B):
        ident[f"bin{b + 1}"] = fs.identifiability(P_bins[b], sp, sps, Cs)
    rng = np.random.default_rng(SEED + 5)
    null, null_ident, t0 = [], [], time.time()
    for p in range(n_perm):
        part = ctx.partition(rng)
        row = []
        for b in range(B):
            P = ctx.prof3d({d: part[d][b] for d in ctx.sets})
            row.append(T.pattern(fs.species_profiles(P, sp, sps, Cg, Cs), with_p=False))
            null_ident.append(fs.identifiability(P, sp, sps, Cs))
        null.append(row)
        if (p + 1) % 10 == 0:
            print(f"    partition {p + 1}/{n_perm} ({(time.time() - t0) / (p + 1):.0f} s each)",
                  flush=True)
    null_sel_ident = [fs.identifiability(P, sp, sps, Cs) for P in ph["draw_profiles"]]
    rows = []
    for k in PATTERN:
        ko = np.array([o[k][0] for o in obs])
        kn = np.array([[x[k][0] for x in part] for part in null])
        pool = kn.ravel()
        L(f"\n  {k} by bin: " + "  ".join(f"{b + 1}: {ko[b]:+.3f} (P {obs[b][k][1]:.3f}, size "
                                         f"{fs.p_high(pool, ko[b]):.2f})" for b in range(B)))
        for b in range(B):
            rows.append(_row("noise bins", f"bin{b + 1}", k, ko[b], obs[b][k][1], pool))
        rho, slope = fs.trend(ko)
        tn = np.array([fs.trend(r) for r in kn])
        for name, s, col in (("trend rho", rho, 0), ("trend slope", slope, 1)):
            p2 = float((1 + np.sum(np.abs(tn[:, col]) >= abs(s) - 1e-12)) / (n_perm + 1))
            L(f"  {name:<12}{s:+.3f}  P two-sided {p2:.3f}  P decrease {fs.p_low(tn[:, col], s):.3f}")
            rows.append(dict(analysis="noise bins", group=name, measure=k, value=s,
                             p_two_sided=p2, p_decrease=fs.p_low(tn[:, col], s)))
    for g, o in ident.items():
        nl = null_sel_ident if g == "selected" else null_ident if g.startswith("bin") else None
        for k in fs.IDENT_SCALARS:
            rows.append(_row("identifiability", g, k, o[k], np.nan,
                             None if nl is None else [z_[k] for z_ in nl]))
    L("\n  species identifiability (raw / residual): gap = r same species - r between species")
    for g, o in ident.items():
        L(f"  {g:<10}gap {o['gap_raw']:+.3f} / {o['gap_res']:+.3f}   sibling hits "
          f"{o['sibling_raw']}/{o['n_sibling']} / {o['sibling_res']}/{o['n_sibling']}")
    return dict(P_bins=P_bins, obs=obs, null=null, ident=ident, null_ident=null_ident,
                null_sel_ident=null_sel_ident, table=pd.DataFrame(rows))


# ── harmonisation: amplitude of the seed maps ───────────────────────────────────────────
def _groups_linear(ctx, X, all_rows):
    """Groups of the linear analyses: X {row: vector}; all_rows {d: matrix of every run}."""
    def G(idx):
        return {d: np.nanmean([X[i] for i in idx[d]], axis=0) for d in idx}
    groups = {"all": {d: np.nanmean(all_rows[d], axis=0) for d in ctx.sets}}
    return G, groups


def _amplitude_like(ctx, name, X, all_rows, idx_all, idx_sel, ranked, n_draws, seed):
    L, T, sp, sps, B = ctx.log, ctx.T, ctx.sp_of_ds, ctx.sps, ctx.B
    G, groups = _groups_linear(ctx, X, all_rows)
    bins = {d: np.array_split(ranked[d], B) for d in ctx.sets}
    sizes = {d: [len(x) for x in bins[d]] for d in ctx.sets}
    groups["selected"] = G(idx_sel)
    for b in range(B):
        groups[f"bin{b + 1}"] = G({d: bins[d][b] for d in ctx.sets})
    _, Cs = fs.commons(groups["all"], sp, sps, zscore=False)
    obs = {g: T.amplitude(Gg, sp, sps, Cs) for g, Gg in groups.items()}
    ratio = {g: fs.amp_ratio(Gg, sp) for g, Gg in groups.items()}
    rng = np.random.default_rng(seed)
    null_sel, null_bin, rat_sel, rat_bin, amp_rand = [], [], [], [], []
    for _ in range(n_draws):
        Gr = G({d: rng.choice(idx_all[d], len(idx_sel[d]), replace=False) for d in ctx.sets})
        null_sel.append(T.amplitude(Gr, sp, sps, Cs, with_p=False))
        rat_sel.append(fs.amp_ratio(Gr, sp))
        amp_rand.append({d: fs.rms(Gr[d]) for d in ctx.sets})
        part = {d: np.split(rng.permutation(ranked[d]), np.cumsum(sizes[d])[:-1]) for d in ctx.sets}
        Gb = G({d: part[d][rng.integers(B)] for d in ctx.sets})
        null_bin.append(T.amplitude(Gb, sp, sps, Cs, with_p=False))
        rat_bin.append(fs.amp_ratio(Gb, sp))
    rows = []
    for g in groups:
        nl = null_sel if g == "selected" else null_bin if g.startswith("bin") else None
        for k in AMPLITUDE:
            rows.append(_row(name, g, k, obs[g][k][0], obs[g][k][1],
                             None if nl is None else [z_[k][0] for z_ in nl]))
        rn = rat_sel if g == "selected" else rat_bin if g.startswith("bin") else None
        rows.append(_row(name, g, "amp_ratio_within_between", ratio[g], np.nan, rn, side="low"))
    amp = pd.DataFrame({g: {ctx.labels[d]: fs.rms(Gg[d]) for d in ctx.sets}
                        for g, Gg in groups.items()})
    rnd = np.array([[a[d] for d in ctx.sets] for a in amp_rand])
    amp["random_median"] = np.median(rnd, axis=0)
    amp["random_p5"] = np.percentile(rnd, 5, axis=0)
    amp["random_p95"] = np.percentile(rnd, 95, axis=0)
    amp["selected_over_random"] = amp["selected"] / amp["random_median"]
    amp.insert(0, "species", [sp[d] for d in ctx.sets])
    tab = pd.DataFrame(rows)
    L(f"\n  {'group':<10}" + "".join(f"{k:>17}" for k in AMPLITUDE) + f"{'amp ratio':>11}")
    for g in groups:
        L(f"  {g:<10}" + "".join(f"{obs[g][k][0]:>+9.3f} ({obs[g][k][1]:.3f})" for k in AMPLITUDE)
          + f"{ratio[g]:>11.2f}")
    L("  amplitude ratio within / between species (size null, P_low: lower = more species-like "
      "than random runs)")
    for g, rn in (("selected", rat_sel),) + tuple((f"bin{b + 1}", rat_bin) for b in range(B)):
        v = np.asarray(rn)
        L(f"    {g:<10}{ratio[g]:.2f}  null {np.median(v):.2f} [{np.percentile(v, 5):.2f}, "
          f"{np.percentile(v, 95):.2f}]  P_low {fs.p_low(v, ratio[g]):.3f}   (all runs {ratio['all']:.2f})")
    return dict(groups=groups, obs=obs, ratio=ratio, rat_sel=np.asarray(rat_sel),
                rat_bin=np.asarray(rat_bin), amp=amp, table=tab)


def amplitude(ctx, n_draws=1000):
    """Harmonisation on the seed maps: mean Fisher z of every homologous region of
    c['common'] (bilateral, group-mask voxels, NOT normalised); a group = mean of the
    run regional means (linear: = regional means of the group mean map)."""
    regions = fd.regions_of_common(ctx.c["common"])
    ctx.log(f"\n== amplitude of the seed maps ({len(regions)} homologous regions, mean Fisher z, "
            f"not normalised; {n_draws} null draws)")
    X, all_rows = {}, {}
    for d in ctx.sets:
        rm = ctx.dm[d].region_means(regions)
        all_rows[d] = np.vstack(list(rm.values()))
        for i in ctx.idx_all[d]:
            X[i] = rm[ctx.mapkey[i]]
    out = _amplitude_like(ctx, "amplitude seed maps", X, all_rows, ctx.idx_all, ctx.idx_sel,
                          ctx.ranked, n_draws, SEED + 13)
    out["regions"] = regions
    ctx.log("\n  amplitude (RMS of the regional mean z) per dataset")
    ctx.log(out["amp"].to_string(float_format=lambda x: f"{x:.3f}"))
    return out


def matrices(ctx, n_draws=1000):
    """Same as amplitude() on the region x region correlation matrices of the runs
    (nodes = the homologous regions, edges present in every dataset), runs paired with
    the seed runs so that 'selected' and the bins are those of the seed analysis."""
    L, rt = ctx.log, ctx.rt
    regions = fd.regions_of_common(ctx.c["common"])
    cb_of = {b["label"]: b for b in ctx.c["bids_configs"]}
    L(f"\n== correlation matrices ({len(regions)} homologous nodes, {n_draws} null draws)")
    raw, idx_all, idx_sel, ranked = {}, {}, {}, {}
    for d in ctx.sets:
        lab = ctx.labels[d]
        runs = fd.matrix_runs(ctx.c, cb_of[lab], regions, log=L)
        rows = rt.rows(d)
        m2v, miss, extra, _ = fd.pair_runs({i: rt.keys[i] for i in rows}, runs)
        for i, k in m2v.items():
            raw[i] = runs[k]
        idx_all[d] = np.array(sorted(m2v), int)
        idx_sel[d] = idx_all[d][rt.selected[idx_all[d]]]
        L(f"  {lab}: {len(runs)} run matrices, {len(m2v)}/{len(rows)} seed runs paired, "
          f"{len(idx_sel[d])}/{int(rt.selected[rows].sum())} selected paired")
        if len(idx_sel[d]) < 2:
            raise RuntimeError(f"{lab}: fewer than 2 selected runs have a matrix; keys e.g. seed "
                               f"{[rt.keys[i] for i in rows[:3]]}, matrix {list(runs)[:3]}")
        ok = idx_all[d][np.isfinite(ctx.score[idx_all[d]])]
        ranked[d] = ok[np.argsort(-ctx.score[ok], kind="stable")]
    edges = sorted(set.intersection(*[set(e for i in idx_all[d] for e in raw[i]) for d in ctx.sets]))
    X = {i: np.array([raw[i].get(e, np.nan) for e in edges]) for i in raw}
    full = np.all([np.isfinite(np.nanmean([X[i] for i in idx_all[d]], axis=0)) for d in ctx.sets],
                  axis=0)
    edges = [e for e, k in zip(edges, full) if k]
    X = {i: v[full] for i, v in X.items()}
    n_intra, n_hom = sum("|intra" in e for e in edges), sum("|homotopic" in e for e in edges)
    L(f"  {len(edges)} edges in every dataset ({n_intra} intra, {len(edges) - n_intra - n_hom} "
      f"inter, {n_hom} homotopic)")
    all_rows = {d: np.vstack([X[i] for i in idx_all[d]]) for d in ctx.sets}
    out = _amplitude_like(ctx, "correlation matrices", X, all_rows, idx_all, idx_sel, ranked,
                          n_draws, SEED + 17)
    out["edges"] = edges
    L("\n  amplitude (RMS of the edge Fisher z) per dataset")
    L(out["amp"].to_string(float_format=lambda x: f"{x:.3f}"))
    return out
