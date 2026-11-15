"""
paper_claims_registry.py -- one entry per number of the manuscript: where it is written,
which file and cell hold the value computed by the analysis, which function writes that
cell, and which datasets it must come from. Checked by paper_claims.py.

Entry keys
    id        short unique name (figure.subject.quantity)
    figure    figure / section of the paper
    what      the quantity, in words
    prefix    text just before the number in the manuscript, typed as in Word (spaces and
              minus signs are normalised: '-' matches '−'); must be unique enough
    k         0 = first number after prefix, 1 = second... (ICC(C,1), citations and figure
              calls are not counted)
    within    optional text that the paragraph must also contain (disambiguation)
    source    dict(file=<relative to PAPER_DIR>, sheet=, where={column: value}, column=)
              or dict(file=, regex=, group=) for a .txt
              or dict(compute=function(src) -> float, file=<main file, for staleness>)
    transform optional function applied to the computed value (unit change, abs...)
    tol       optional absolute tolerance (default: half a unit of the last written digit;
              10% for numbers written with ≈ or ~)
    producer  (file relative to the repository, function, literal of the line that writes
              the value) -> file:line in the report
    data      dataset names the value must come from; checked as literals in `data_files`
              (default: the producer file)

Coverage: Fig. 2 (validation paragraph), Fig. 3 / Supp. Fig. 9 (longitudinal), Supp. Fig. 7
(FreeSurfer), Supp. Fig. 8 (longitudinal vs cross-sectional), Fig. 4 / Supp. Figs. 11-13 (run
selection), Fig. 6 (FDG and multimodal fingerprints; prefixes of the V5 text, to update when
the Fig. 6 paragraph is rewritten). Numbers without an entry appear in the 'untraced' sheet.
"""
import numpy as np

A = 'Exemples/Study_EDNiX/analysis/'
LONG = A + 'library/longitudinal.py'
HAP = A + 'library/human_assoc_primary.py'
VALID = A + 'validation_morphometry.py'
EXTRA = A + 'paper_extra_stats.py'
MORPHO = A + 'library/morphometry_figure.py'
FS = A + 'validation/supp_freesurfer_comparison.py'
ABL = A + 'library/longitudinal_ablation.py'
PDATA = A + 'paper_data.py'

F3 = 'Figure_3/Figure_3_stats.xlsx'
S7 = 'Supp_Fig_09/Supp_Fig_09_stats.xlsx'
VAL = 'Validation/validation_morphometry.xlsx'
F2 = 'Figure_2/Figure_2_stats.xlsx'
XS = 'Supplementary_stats/extra_stats.txt'
S14 = 'Supp_Fig_07/Supp_Fig_07_stats.xlsx'
S15 = 'Supp_Fig_08/Supp_Fig_08_stats.xlsx'


# ── helpers ─────────────────────────────────────────────────────────────────
def f3(species, modality, column, kind='trajectory'):
    return dict(file=F3, sheet='panels', column=column,
                where=dict(species=species, modality=modality, kind=kind))


def val(sheet, species, modality, column):
    return dict(file=VAL, sheet=sheet, column=column, where=dict(species=species, modality=modality))


def rcv(src, modality, species=None, region=None):
    """Robust CV (IQR / 1.349 / median, %) of the Figure 2 table (species x lobe)."""
    t = src.sheet(F2, 'values_by_species')
    t = t[t['modality'].isin(modality if isinstance(modality, list) else [modality])]
    if species:
        t = t[t['species'] == species]
    if region:
        t = t[t['region'] == region]
    return 100 * (t['q75'] - t['q25']) / 1.349 / t['median']


def rcv_stat(modality, how, species=None, exclude=None):
    def f(src):
        v = rcv(src, modality, species)
        if exclude:
            t = src.sheet(F2, 'values_by_species')
            t = t[t['modality'].isin(modality if isinstance(modality, list) else [modality])]
            if species:
                t = t[t['species'] == species]
            keep = ~((t['species'] == exclude[0]) & (t['region'] == exclude[1])).to_numpy()
            v = v[keep]
        return float(getattr(np, how)(v))
    return dict(compute=f, file=F2)


def s15(species, modality, column):
    return dict(file=S15, sheet='summary', column=column,
                where=dict(species=species, modality=modality))


def s15_range(column, how):
    return dict(compute=lambda src: float(getattr(np, how)(src.sheet(S15, 'summary')[column])),
                file=S15)


def s14(measure, column):
    return dict(file=S14, sheet='stats', column=column, where=dict(measure=measure))


cm2 = (lambda v: v / 100.0)          # mm² -> cm²

# ── Fig. 4 / Supp. 11-13 (run selection) and Fig. 6 (FDG) ───────────────────
SEL = A + 'fmri_fingerprints/selection_validation.py'
PSTATS = A + 'paper_stats.py'
MMF = A + 'library/multimodal_figure.py'
PETD = A + 'library/pet_datasets.py'
CFG_SEL = 'Statistics/Group_fMRI/run_selection/config_ednix.py'
F4 = 'Figure_4/Figure_4_stats.xlsx'
S11 = 'Supp_Fig_11/Supp_Fig_11_stats.xlsx'
S13 = 'Supp_Fig_13/Supp_Fig_13_stats.xlsx'
F6 = 'Figure_6/Figure_6_stats.xlsx'


def f4(sheet, column, **where):
    return dict(file=F4, sheet=sheet, column=column, where=where)


def pair_level(pool, cls):
    return f4('b_pair_levels', 'median', pool=pool, pair_class=cls)


def _random_median(file, sheet, column, **where):
    def f(src):
        t = src.sheet(file, sheet)
        for k, v in where.items():
            t = t[t[k] == v]
        return float(np.median(t.loc[t['pool'] == 'random', column]))
    return dict(compute=f, file=file)


def _share_random_ge_selected(file, sheet, column, **where):
    """One-sided P of 'selected vs random': (1 + #random >= selected) / (n random + 1), as in
    selection_validation.species_signal."""
    def f(src):
        t = src.sheet(file, sheet)
        for k, v in where.items():
            t = t[t[k] == v]
        sel = float(t.loc[t['pool'] == 'selected', column].iloc[0])
        r = t.loc[t['pool'] == 'random', column].to_numpy(float)
        return float((1 + np.sum(r >= sel - 1e-12)) / (len(r) + 1))
    return dict(compute=f, file=file)


def f6_n(dataset):
    return dict(compute=lambda src: float((src.sheet(F6, 'B_bratio')['dataset'] == dataset).sum()), file=F6)


def f6_cons(modality, column='r'):
    return dict(file=F6, sheet='D_summary_stats', column=column,
                where=dict(analysis='conservation', modality=modality))


CLAIMS = [
    # ═════════════ Fig. 2 -- morphometry validation (Results "Morphometry across seven species")
    dict(id='F2.human.area', figure='Fig. 2', what='human isocortical pial area (cm²), median of subjects',
         prefix='isocortical pial surface area was ', source=val('isocortex_by_species', 'Human', 'surface', 'median'),
         transform=cm2, tol=50, producer=(VALID, 'cross_sectional', "_summary(df, col, ['species'])"),
         data=['BIDS_Park'], data_files=[PDATA]),
    dict(id='F2.macaque.area', figure='Fig. 2', what='macaque isocortical pial area (cm²)',
         prefix='against 1,946 cm² for the midthickness neocortex [NEW-REF 3]; ',
         source=val('isocortex_by_species', 'Macaque', 'surface', 'median'), transform=cm2, tol=5,
         producer=(VALID, 'cross_sectional', "_summary(df, col, ['species'])")),
    dict(id='F2.marmoset.area', figure='Fig. 2', what='marmoset isocortical pial area (cm²)',
         prefix='reconstructed by CIVET-Macaque [NEW-REF 4]; and ',
         source=val('isocortex_by_species', 'Marmoset', 'surface', 'median'), transform=cm2,
         producer=(VALID, 'cross_sectional', "_summary(df, col, ['species'])")),
    dict(id='F2.human.thick', figure='Fig. 2', what='human mean thickness (mm)', prefix='Human cortical thickness (',
         source=val('isocortex_by_species', 'Human', 'thickness', 'median'),
         producer=(VALID, 'cross_sectional', "_summary(df, col, ['species'])")),
    dict(id='F2.mouse.thick', figure='Fig. 2', what='mouse mean thickness (mm)', prefix='and mouse thickness (',
         source=val('isocortex_by_species', 'Mouse', 'thickness', 'median'),
         producer=(VALID, 'cross_sectional', "_summary(df, col, ['species'])")),
    dict(id='F2.mouse.thick.diff', figure='Fig. 2', what='mouse thickness vs Lerch 2008 (% difference, absolute)',
         prefix='was within ', source=val('literature', 'Mouse', 'thickness', 'diff_pct'),
         transform=abs, producer=(VALID, 'literature_table', 'diff_pct=')),
    dict(id='F2.macaque.thick', figure='Fig. 2', what='macaque mean thickness (mm)', prefix='(macaque ',
         source=val('isocortex_by_species', 'Macaque', 'thickness', 'median'),
         producer=(VALID, 'cross_sectional', "_summary(df, col, ['species'])")),
    dict(id='F2.marmoset.thick', figure='Fig. 2', what='marmoset mean thickness (mm)', prefix='; marmoset ',
         within='Thickness deviated more', source=val('isocortex_by_species', 'Marmoset', 'thickness', 'median'),
         producer=(VALID, 'cross_sectional', "_summary(df, col, ['species'])")),
    dict(id='F2.rcv.sv.min', figure='Fig. 2', what='robust CV, surface + volume, min over species x lobe (%)',
         prefix='relative to the median) was ', source=rcv_stat(['surface', 'volume'], 'min'),
         producer=(MORPHO, 'figure2_table', "q25=")),
    dict(id='F2.rcv.sv.max', figure='Fig. 2', what='robust CV, surface + volume, max excluding marmoset insula (%)',
         prefix='relative to the median) was ', k=1,
         source=rcv_stat(['surface', 'volume'], 'max', exclude=('Marmoset', 'Insular_lobe')),
         producer=(MORPHO, 'figure2_table', "q25=")),
    dict(id='F2.rcv.sv.median', figure='Fig. 2', what='robust CV, surface + volume, median (%)',
         prefix='for surface area and volume (median ', source=rcv_stat(['surface', 'volume'], 'median'),
         producer=(MORPHO, 'figure2_table', "q25=")),
    dict(id='F2.rcv.th.min', figure='Fig. 2', what='robust CV, thickness, min (%)', prefix='(median 9%) and ',
         source=rcv_stat('thickness', 'min'), producer=(MORPHO, 'figure2_table', "q25=")),
    dict(id='F2.rcv.th.max', figure='Fig. 2', what='robust CV, thickness, max (%)', prefix='(median 9%) and ', k=1,
         source=rcv_stat('thickness', 'max'), producer=(MORPHO, 'figure2_table', "q25=")),
    dict(id='F2.rcv.th.median', figure='Fig. 2', what='robust CV, thickness, median (%)',
         prefix='for thickness (median ', source=rcv_stat('thickness', 'median'),
         producer=(MORPHO, 'figure2_table', "q25=")),
    dict(id='F2.rcv.insula', figure='Fig. 2', what='robust CV, marmoset insula surface (%)',
         prefix='marmoset insula (', source=dict(compute=lambda s: float(rcv(s, 'surface', 'Marmoset', 'Insular_lobe').iloc[0]), file=F2),
         producer=(MORPHO, 'figure2_table', "q25=")),
    dict(id='F2.human.rcv.min', figure='Fig. 2', what='human surface robust CV, min over lobes (%)',
         prefix='(human surface area, ', source=rcv_stat('surface', 'min', species='Human'),
         producer=(MORPHO, 'figure2_table', "q25=")),
    dict(id='F2.human.rcv.max', figure='Fig. 2', what='human surface robust CV, max over lobes (%)',
         prefix='(human surface area, ', k=1, source=rcv_stat('surface', 'max', species='Human'),
         producer=(MORPHO, 'figure2_table', "q25=")),
    dict(id='F2.mouse.rcv.min', figure='Fig. 2', what='mouse surface robust CV, min (%) [assumed surface area]',
         prefix='in plane (mouse, ', source=rcv_stat('surface', 'min', species='Mouse'),
         producer=(MORPHO, 'figure2_table', "q25=")),
    dict(id='F2.mouse.rcv.max', figure='Fig. 2', what='mouse surface robust CV, max (%) [assumed surface area]',
         prefix='in plane (mouse, ', k=1, source=rcv_stat('surface', 'max', species='Mouse'),
         producer=(MORPHO, 'figure2_table', "q25=")),
    dict(id='F2.rat.rcv.min', figure='Fig. 2', what='rat surface robust CV, min (%) [assumed surface area]',
         prefix='; rat, ', source=rcv_stat('surface', 'min', species='Rat'), producer=(MORPHO, 'figure2_table', "q25=")),
    dict(id='F2.rat.rcv.max', figure='Fig. 2', what='rat surface robust CV, max (%) [assumed surface area]',
         prefix='; rat, ', k=1, source=rcv_stat('surface', 'max', species='Rat'),
         producer=(MORPHO, 'figure2_table', "q25=")),
    dict(id='F2.fold', figure='Fig. 2', what='human / mouse isocortical area ratio', prefix='negligible against the ',
         source=dict(compute=lambda s: s.cell(VAL, 'isocortex_by_species', dict(species='Human', modality='surface'), 'median')
                     / s.cell(VAL, 'isocortex_by_species', dict(species='Mouse', modality='surface'), 'median'), file=VAL),
         tol=150, producer=(VALID, 'cross_sectional', "_summary(df, col, ['species'])")),
    dict(id='F2.scaling.R', figure='Fig. 2', what='area vs volume, Pearson R, 35 species x lobe medians',
         prefix='strongly correlated (R = ', source=dict(file=XS, regex=r'Pearson R = ([-\d.]+); species-block'),
         producer=(EXTRA, 'fig2_scaling', 'r = stats.pearsonr(d.s, d.v)[0]')),
    dict(id='F2.scaling.ci.lo', figure='Fig. 2', what='species-block bootstrap CI low', prefix='95% CI ',
         within='strongly correlated', source=dict(file=XS, regex=r'95% CI \[([-\d.]+),'),
         producer=(EXTRA, 'fig2_scaling', 'lo, hi = np.percentile(boot')),
    dict(id='F2.scaling.ci.hi', figure='Fig. 2', what='species-block bootstrap CI high', prefix='95% CI ', k=1,
         within='strongly correlated', source=dict(file=XS, regex=r'95% CI \[[-\d.]+, ([-\d.]+)\] \('),
         producer=(EXTRA, 'fig2_scaling', 'lo, hi = np.percentile(boot')),
    dict(id='F2.scaling.p', figure='Fig. 2', what='species-block permutation P', prefix='species-block permutation P ',
         source=dict(file=XS, regex=r'one-sided P = ([-\d.e]+);'), producer=(EXTRA, 'fig2_scaling', 'p_perm =')),
    dict(id='F2.scaling.R7', figure='Fig. 2', what='7 species totals, Pearson R', prefix='on the seven species totals R = ',
         source=dict(file=XS, regex=r'species totals \(sum of the five lobes\): Pearson R = ([-\d.]+)'),
         producer=(EXTRA, 'fig2_scaling', 'r7 = stats.pearsonr')),
    dict(id='F2.scaling.rho7', figure='Fig. 2', what='7 species totals, Spearman rho', prefix='Spearman ρ = ',
         within='seven species totals', source=dict(file=XS, regex=r'Spearman rho = ([-\d.]+), exact'),
         producer=(EXTRA, 'fig2_scaling', 'rho7 = stats.spearmanr')),
    dict(id='F2.scaling.p7', figure='Fig. 2', what='7 species totals, exact P', prefix='exact P = ',
         within='seven species totals', source=dict(file=XS, regex=r'exact one-sided P = ([-\d.e]+)'),
         producer=(EXTRA, 'fig2_scaling', 'p7 = float')),

    # ═════════════ Fig. 3 -- macaque development (BIDS_Zhu_Garin)
    dict(id='F3.mac.n', figure='Fig. 3', what='macaques with >= 2 sessions', prefix='In ',
         within='rhesus macaques scanned every 3 months', source=f3('Macaque', 'thickness', 'n_subjects'),
         producer=(LONG, 'figure3_stats', "n_subjects=rdf['subject'].nunique()"), data=['BIDS_Zhu_Garin']),
    dict(id='F3.mac.age.min', figure='Fig. 3', what='youngest age (years)', prefix='every 3 months from ',
         source=val('annual_change', 'Macaque', 'thickness', 'age_min'),
         producer=(VALID, 'annual_change', 'age_min='), data=['BIDS_Zhu_Garin'], data_files=[LONG]),
    dict(id='F3.mac.age.max', figure='Fig. 3', what='oldest age (years)', prefix='every 3 months from ', k=1,
         source=val('annual_change', 'Macaque', 'thickness', 'age_max'), producer=(VALID, 'annual_change', 'age_max=')),
    dict(id='F3.mac.thick.slope', figure='Fig. 3e', what='thickness LME slope (mm/yr)', prefix='(linear slope ',
         within='rhesus macaques', source=f3('Macaque', 'thickness', 'lme_slope_per_year'),
         producer=(LONG, 'figure3_stats', 'lme_slope_per_year=')),
    dict(id='F3.mac.thick.se', figure='Fig. 3e', what='thickness LME slope SE', prefix='(linear slope ', k=1,
         within='rhesus macaques', source=f3('Macaque', 'thickness', 'lme_se_per_year'),
         producer=(LONG, 'figure3_stats', 'lme_se_per_year=')),
    dict(id='F3.mac.thick.p', figure='Fig. 3e', what='thickness LME slope P', prefix='mm per year, P ',
         within='rhesus macaques', source=f3('Macaque', 'thickness', 'lme_p_slope'),
         producer=(LONG, 'figure3_stats', 'lme_p_slope=')),
    dict(id='F3.mac.thick.pnl', figure='Fig. 3e', what='thickness non-linearity LRT P',
         prefix='spline versus linear, likelihood-ratio P ', source=f3('Macaque', 'thickness', 'gamm_p_nonlinear'),
         producer=(LONG, 'figure3_stats', 'gamm_p_nonlinear=')),
    dict(id='F3.mac.area.slope', figure='Fig. 3d', what='surface LME slope (mm²/yr, left)', prefix='surface area increased modestly (',
         source=f3('Macaque', 'surface', 'lme_slope_per_year'), producer=(LONG, 'figure3_stats', 'lme_slope_per_year=')),
    dict(id='F3.mac.area.se', figure='Fig. 3d', what='surface LME slope SE', prefix='surface area increased modestly (', k=1,
         source=f3('Macaque', 'surface', 'lme_se_per_year'), producer=(LONG, 'figure3_stats', 'lme_se_per_year=')),
    dict(id='F3.mac.area.p', figure='Fig. 3d', what='surface LME slope P', prefix='mm² per year, P = ',
         within='rhesus macaques', source=f3('Macaque', 'surface', 'lme_p_slope'),
         producer=(LONG, 'figure3_stats', 'lme_p_slope=')),
    dict(id='F3.mac.area.pnl', figure='Fig. 3d', what='surface non-linearity LRT P', prefix='non-linear, P = ',
         within='rhesus macaques', source=f3('Macaque', 'surface', 'gamm_p_nonlinear'),
         producer=(LONG, 'figure3_stats', 'gamm_p_nonlinear=')),
    dict(id='F3.mac.vol.p', figure='Fig. 3f', what='volume LME slope P [text does not say which P]',
         prefix='volume did not change significantly (P = ', source=f3('Macaque', 'volume', 'lme_p_slope'),
         producer=(LONG, 'figure3_stats', 'lme_p_slope=')),
    dict(id='F3.mac.area.apc', figure='Fig. 3d', what='surface % per year (pooled LME slope / mean first value)',
         prefix='surface area grew by ', source=val('annual_change', 'Macaque', 'surface', 'apc_pooled'),
         producer=(VALID, 'annual_change', 'apc_pooled=')),
    ] + [
    dict(id=f'F3.mac.ap.{mod}.{q}', figure=f'Fig. 3{L}', what=f'macaque association vs primary, {mod}: {q}',
         prefix=pre, k=k, within='rhesus macaques',
         source=f3('Macaque', mod, col, kind='assoc_vs_primary'),
         producer=(LONG, 'figure3_stats', key))
    for mod, L, pre in (('surface', 'g', 'primary-dominated cortex in surface area ('),
                        ('thickness', 'h', '), thickness ('), ('volume', 'i', ') and volume ('))
    for q, k, col, key in (('assoc %', 0, 'assoc_pct', "assoc_pct=st['assoc_pct']"),
                           ('primary %', 1, 'primary_pct', "primary_pct=st['primary_pct']"),
                           ('t', 3, 't', "t=st['t']"))
    ] + [
    dict(id='F3.mac.ap.surface.dz', figure='Fig. 3g', what="macaque assoc vs primary surface, Cohen's d_z",
         prefix="Cohen's d_z = ", within='rhesus macaques', source=f3('Macaque', 'surface', 'dz', kind='assoc_vs_primary'),
         producer=(LONG, 'figure3_stats', "dz=st['dz']")),

    # ═════════════ Fig. 3 -- human ageing (BIDS_Park)
    dict(id='F3.hum.n', figure='Fig. 3', what='humans with >= 2 sessions', prefix='In the opposite direction, in ',
         source=f3('Human', 'surface', 'n_subjects'), producer=(LONG, 'figure3_stats', "n_subjects=rdf['subject'].nunique()"),
         data=['BIDS_Park']),
    dict(id='F3.hum.area.slope', figure='Fig. 3a', what='surface LME slope (mm²/yr, left)', prefix='isocortical surface area (',
         within='Dallas Lifespan', source=f3('Human', 'surface', 'lme_slope_per_year'), producer=(LONG, 'figure3_stats', 'lme_slope_per_year=')),
    dict(id='F3.hum.area.se', figure='Fig. 3a', what='surface LME slope SE', prefix='isocortical surface area (', k=1,
         within='Dallas Lifespan', source=f3('Human', 'surface', 'lme_se_per_year'), producer=(LONG, 'figure3_stats', 'lme_se_per_year=')),
    dict(id='F3.hum.area.p', figure='Fig. 3a', what='surface LME slope P', prefix='mm² per year, P ',
         within='Dallas Lifespan', source=f3('Human', 'surface', 'lme_p_slope'), producer=(LONG, 'figure3_stats', 'lme_p_slope=')),
    dict(id='F3.hum.vol.slope', figure='Fig. 3c', what='volume LME slope (mm³/yr, left)', prefix='and volume (',
         within='Dallas Lifespan', source=f3('Human', 'volume', 'lme_slope_per_year'), producer=(LONG, 'figure3_stats', 'lme_slope_per_year=')),
    dict(id='F3.hum.vol.se', figure='Fig. 3c', what='volume LME slope SE', prefix='and volume (', k=1,
         within='Dallas Lifespan', source=f3('Human', 'volume', 'lme_se_per_year'), producer=(LONG, 'figure3_stats', 'lme_se_per_year=')),
    dict(id='F3.hum.vol.p', figure='Fig. 3c', what='volume LME slope P', prefix='mm³ per year, P ',
         within='Dallas Lifespan', source=f3('Human', 'volume', 'lme_p_slope'), producer=(LONG, 'figure3_stats', 'lme_p_slope=')),
    dict(id='F3.hum.thick.pnl', figure='Fig. 3b', what='thickness non-linearity LRT P', prefix='thickness declined non-linearly (likelihood-ratio P = ',
         source=f3('Human', 'thickness', 'gamm_p_nonlinear'), producer=(LONG, 'figure3_stats', 'gamm_p_nonlinear=')),
    dict(id='F3.hum.thick.slope', figure='Fig. 3b', what='thickness LME slope (mm/yr)', prefix='linear slope ',
         within='Dallas Lifespan', source=f3('Human', 'thickness', 'lme_slope_per_year'), producer=(LONG, 'figure3_stats', 'lme_slope_per_year=')),
    dict(id='F3.hum.thick.p', figure='Fig. 3b', what='thickness LME slope P', prefix='mm per year, P = ',
         within='Dallas Lifespan', source=f3('Human', 'thickness', 'lme_p_slope'), producer=(LONG, 'figure3_stats', 'lme_p_slope=')),
    dict(id='F3.hum.apc.area', figure='Fig. 3a', what='surface loss % per year (pooled, absolute)', prefix='correspond to a loss of ≈',
         source=val('annual_change', 'Human', 'surface', 'apc_pooled'), transform=abs, tol=0.05,
         producer=(VALID, 'annual_change', 'apc_pooled=')),
    dict(id='F3.hum.apc.vol', figure='Fig. 3c', what='volume loss % per year (pooled, absolute)', prefix='correspond to a loss of ≈', k=1,
         source=val('annual_change', 'Human', 'volume', 'apc_pooled'), transform=abs, tol=0.05,
         producer=(VALID, 'annual_change', 'apc_pooled=')),
    dict(id='S7.vol.assoc', figure='Supp. Fig. 7c', what='human association % change, volume', prefix='lost more volume than primary-dominated cortex (',
         source=dict(file=S7, sheet='panels', where=dict(modality='volume'), column='assoc_pct'), producer=(HAP, 'supp6_stats', "assoc_pct=st['assoc_pct']")),
    dict(id='S7.vol.prim', figure='Supp. Fig. 7c', what='human primary % change, volume', prefix='lost more volume than primary-dominated cortex (', k=1,
         source=dict(file=S7, sheet='panels', where=dict(modality='volume'), column='primary_pct'), producer=(HAP, 'supp6_stats', "primary_pct=st['primary_pct']")),
    dict(id='S7.vol.t', figure='Supp. Fig. 7c', what='human assoc vs primary volume, paired t', prefix='paired t(9) = ',
         source=dict(file=S7, sheet='panels', where=dict(modality='volume'), column='t'), producer=(HAP, 'supp6_stats', "t=st['t']")),
    dict(id='S7.vol.p', figure='Supp. Fig. 7c', what='human assoc vs primary volume, P', prefix='paired t(9) = ', k=1,
         source=dict(file=S7, sheet='panels', where=dict(modality='volume'), column='p'), producer=(HAP, 'supp6_stats', "p=st['p']")),
    dict(id='S7.area.p', figure='Supp. Fig. 7a', what='human assoc vs primary surface, P', prefix='non-significant trend for surface area (P = ',
         source=dict(file=S7, sheet='panels', where=dict(modality='surface'), column='p'), producer=(HAP, 'supp6_stats', "p=st['p']")),
    dict(id='S7.thick.p', figure='Supp. Fig. 7b', what='human assoc vs primary thickness, P', prefix='none for thickness (P = ',
         source=dict(file=S7, sheet='panels', where=dict(modality='thickness'), column='p'), producer=(HAP, 'supp6_stats', "p=st['p']")),
    dict(id='F3.legend.mac.sessions', figure='Fig. 3 legend', what='macaque sessions in the trajectories', prefix='n = 8 macaques (',
         source=f3('Macaque', 'thickness', 'n_sessions'), producer=(LONG, 'figure3_stats', 'n_sessions=len(rdf)')),
    dict(id='F3.legend.hum.sessions', figure='Fig. 3 legend', what='human sessions in the trajectories', prefix='and 10 humans (',
         source=f3('Human', 'surface', 'n_sessions'), producer=(LONG, 'figure3_stats', 'n_sessions=len(rdf)')),

    # ═════════════ Supp. Fig. 7 -- EDNiX vs FreeSurfer, human (BIDS_Park); PR #77
    dict(id='S14.n', figure='Supp. Fig. 7', what='subjects with both FreeSurfer and EDNiX values', prefix='(latest session; n = ',
         source=s14('Surface area', 'n'), producer=(FS, 'agreement', 'return dict(n=n'), data=['BIDS_Park']),
    ] + [
    dict(id=f'S14.{tag}.{q}', figure='Supp. Fig. 7', what=f'{meas}: {q}', prefix=pre, k=k, within='FreeSurfer recon-all',
         source=s14(meas, col), producer=(FS, 'agreement', key))
    for meas, tag, pre in (('Surface area', 'area', 'agreement between EDNiX and FreeSurfer was r = '),
                           ('Volume', 'vol', ') for surface area, r = '),
                           ('Thickness', 'thick', ') for grey-matter volume and r = '))
    for q, k, col, key in (('Pearson r', 0, 'pearson_r', 'pearson_r=r'), ('ICC(C,1)', 1, 'icc_c1', 'icc_c1=ic'))
    ] + [
    dict(id=f'S14.{tag}.bias', figure='Supp. Fig. 7', what=f'{meas}: Bland-Altman bias (%)',
         prefix='(EDNiX - FreeSurfer, % of the mean) of ', k=k, source=s14(meas, 'bias_pct'),
         producer=(FS, 'agreement', 'bias_pct=bias'))
    for k, (meas, tag) in enumerate((('Surface area', 'area'), ('Volume', 'vol'), ('Thickness', 'thick')))
    ] + [
    dict(id='S14.thick.vertexdef.bias', figure='Supp. Fig. 7', what='thickness bias with the EDNiX definition on FreeSurfer surfaces (%)',
         prefix='on the FreeSurfer surfaces, the thickness offset was ',
         source=s14('Thickness (FreeSurfer surfaces, EDNiX definition)', 'bias_pct'), producer=(FS, 'fs_measures', 'thick_v.append')),

    # ═════════════ Supp. Fig. 8 -- longitudinal vs cross-sectional processing
    dict(id='S15.ref.diff', figure='Supp. Fig. 8', what='reference session |difference|, max over species x modality (%)',
         prefix='gave the same values (median absolute difference ', source=s15_range('ref_abs_diff_pct_median', 'max'),
         producer=(ABL, 'compare_modes', "out['ref_abs_diff_pct_median']"), data=['BIDS_Park', 'BIDS_Zhu_Garin']),
    dict(id='S15.sigma.min', figure='Supp. Fig. 15g', what='residual SD ratio cross/long, min', prefix='residual standard deviation was ',
         source=s15_range('sigma_e_ratio', 'min'), producer=(ABL, 'compare_modes', "out['sigma_e_ratio'] =")),
    dict(id='S15.sigma.max', figure='Supp. Fig. 15g', what='residual SD ratio cross/long, max', prefix='residual standard deviation was ', k=1,
         source=s15_range('sigma_e_ratio', 'max'), producer=(ABL, 'compare_modes', "out['sigma_e_ratio'] =")),
    dict(id='S15.neq.min', figure='Supp. Fig. 8', what='sample-size equivalence (SE ratio²), min', prefix='would need ',
         within='cross-sectional mode', source=s15_range('n_equivalent_ratio', 'min'), producer=(ABL, 'compare_modes', "out['n_equivalent_ratio']")),
    dict(id='S15.neq.max', figure='Supp. Fig. 8', what='sample-size equivalence (SE ratio²), max', prefix='would need ', k=1,
         within='cross-sectional mode', source=s15_range('n_equivalent_ratio', 'max'), producer=(ABL, 'compare_modes', "out['n_equivalent_ratio']")),
    dict(id='S15.delta.pmin', figure='Supp. Fig. 15h', what='slope difference long - cross, smallest bootstrap P',
         prefix='slope difference between modes, smallest P = ', source=s15_range('delta_slope_p', 'min'),
         producer=(ABL, 'compare_modes', "out[f'{name}_p'] = _boot_p")),
] + [
    # ═════════════ Fig. 4 -- run selection (Results "A cross-species criterion for selecting resting-state runs")
    dict(id='F4.sib.same', figure='Fig. 4c', what='mean distance 1 - r, same-species dataset pairs',
         prefix='mean distance 1 - r = ', source=f4('c_sibling_test', 'mean_same'),
         producer=(PSTATS, '_sibling_tables', 'sibling_test=')),
    dict(id='F4.sib.cross', figure='Fig. 4c', what='mean distance 1 - r, different-species dataset pairs',
         prefix='mean distance 1 - r = ', k=1, source=f4('c_sibling_test', 'mean_cross'),
         producer=(PSTATS, '_sibling_tables', 'sibling_test=')),
    dict(id='F4.sib.p', figure='Fig. 4c', what='sibling test, one-sided label-permutation P',
         prefix='permutation of the species labels across datasets, P = ', source=f4('c_sibling_test', 'p_one_sided'),
         producer=(PSTATS, '_sibling_tables', 'sibling_test=')),
] + [
    dict(id=f'F4.b.{tag}.{pool}', figure='Fig. 4b', what=f'held-out run-to-run r, {cls}, {pool} (median)',
         prefix=prefix, k=k, source=pair_level(pool, cls), producer=(SEL, 'heldout', 'levels = pd.DataFrame(rows)'))
    for tag, cls, prefix in (('ds', 'same dataset', 'median correlation between two runs of the same dataset '),
                             ('sp', 'same species, other dataset', 'of the same species in different datasets '),
                             ('diff', 'different species', 'and of different species '))
    for k, pool in ((0, 'selected'), (1, 'random'))
] + [
    dict(id='F4.b.pmax', figure='Fig. 4b', what='largest P of the three pair classes (selected vs random)',
         prefix='(each P ≤ ', source=dict(compute=lambda src: float(src.sheet(F4, 'b_pair_class_p')
                                                                   ['p_selected_vs_random'].max()), file=F4),
         tol=0.0005, producer=(SEL, 'heldout', 'p_cls = np.array')),
    # Supp. Fig. 11b: raw species gap
    dict(id='S11b.gap.sel', figure='Supp. Fig. 11b', what='raw species gap, selected runs',
         prefix='than for 10 random runs per dataset (', source=dict(file=S11, sheet='b_species_gap_raw',
                                                                   column='gap', where=dict(pool='selected')),
         producer=(SEL, 'raw_gap_table', 'gap=[')),
    dict(id='S11b.gap.rand', figure='Supp. Fig. 11b', what='raw species gap, random draws (median)',
         prefix='than for 10 random runs per dataset (', k=1,
         source=_random_median(S11, 'b_species_gap_raw', 'gap'), producer=(SEL, 'raw_gap_table', 'gap=[')),
    dict(id='S11b.gap.p', figure='Supp. Fig. 11b', what='raw species gap, selected vs random P',
         prefix='than for 10 random runs per dataset (', k=2,
         source=dict(file=S11, sheet='b_species_gap_raw', column='p_selected_vs_random', where=dict(pool='selected')),
         producer=(SEL, 'species_signal', 'p_vs_random=float((1 + np.sum(rr')),
    # Fig. 4d: species identity (residual profiles)
    dict(id='F4.d.sel', figure='Fig. 4d', what='species gap after removal of the common profile, selected runs',
         prefix='same-species minus different-species correlation ', source=f4('d_species_signal', 'gap', pool='selected'),
         producer=(PSTATS, 'figure4', 'd_species_signal=')),
    dict(id='F4.d.p', figure='Fig. 4d', what='species gap, selected runs, label-permutation P',
         prefix='permutation of the species labels, P = ',
         source=f4('d_species_signal', 'p_label_permutation', pool='selected'), producer=(PSTATS, 'figure4', 'd_species_signal=')),
    dict(id='F4.d.rand', figure='Fig. 4d', what='species gap, random draws (median)',
         prefix='comparable to that of random runs (median ', source=_random_median(F4, 'd_species_signal', 'gap'),
         producer=(PSTATS, 'figure4', 'd_species_signal=')),
    dict(id='F4.d.pvs', figure='Fig. 4d', what='species gap, selected vs random P',
         prefix='comparable to that of random runs (median ', k=1,
         source=_share_random_ge_selected(F4, 'd_species_signal', 'gap'), producer=(SEL, 'species_signal', 'p_vs_random=float(')),
    # Supp. Fig. 13c: units kept per state
    dict(id='S13c.marmoset.awake', figure='Supp. Fig. 13c', what='Marmoset_Hori: units kept awake',
         prefix='(marmoset: ', source=dict(file=S13, sheet='c_kept_by_state', column='kept_awake',
                                           where=dict(label='Marmoset_Hori')), producer=(PSTATS, '_kept_by_state', 'kept_awake=ka')),
    dict(id='S13c.marmoset.anaesth', figure='Supp. Fig. 13c', what='Marmoset_Hori: units kept anaesthetised',
         prefix='(marmoset: ', k=1, source=dict(file=S13, sheet='c_kept_by_state', column='kept_anaesthetised',
                                                where=dict(label='Marmoset_Hori')), producer=(PSTATS, '_kept_by_state', 'kept_anaesthetised=kn')),
    # Fig. 4e: phylogenetic gradient
    dict(id='F4.e.all.rho', figure='Fig. 4e', what='Spearman Mantel rho, all runs (raw profiles)',
         prefix='Spearman Mantel ρ = ', source=f4('e_phylo_mantel', 'rho', profiles='raw', pool='all runs'),
         producer=(PSTATS, 'figure4', 'e_phylo_mantel=')),
    dict(id='F4.e.all.p', figure='Fig. 4e', what='Mantel P, all runs', prefix='Spearman Mantel ρ = ', k=1,
         source=f4('e_phylo_mantel', 'p', profiles='raw', pool='all runs'), producer=(PSTATS, 'figure4', 'e_phylo_mantel=')),
    dict(id='F4.e.rand.rho', figure='Fig. 4e', what='Mantel rho, random draws (median)',
         prefix='nor with 10 random runs per dataset (median ρ = ',
         source=_random_median(F4, 'e_phylo_mantel', 'rho', profiles='raw'), producer=(PSTATS, 'figure4', 'e_phylo_mantel=')),
    dict(id='F4.e.sel.rho', figure='Fig. 4e', what='Mantel rho, selected runs', prefix='With the selected runs, ρ was ',
         source=f4('e_phylo_mantel', 'rho', profiles='raw', pool='selected'), producer=(PSTATS, 'figure4', 'e_phylo_mantel=')),
    dict(id='F4.e.pvs', figure='Fig. 4e', what='Mantel rho, selected vs random P',
         prefix='share of random draws at least as large, P = ',
         source=_share_random_ge_selected(F4, 'e_phylo_mantel', 'rho', profiles='raw'),
         producer=(SEL, 'species_signal', "p_vs_random=float((1 + np.sum(rr >= m_sel[k][0]")),
    dict(id='F4.e.delta', figure='Fig. 4e', what='delta rho selected - all runs', prefix='Δρ = ',
         source=f4('e_phylo_delta', 'delta', profiles='raw'), producer=(SEL, 'mantel_delta', 'delta = rho(Ds, ident)')),
    dict(id='F4.e.delta.p', figure='Fig. 4e', what='delta rho, one-sided P (shared permutation)', prefix='Δρ = ', k=1,
         source=f4('e_phylo_delta', 'p', profiles='raw'), producer=(SEL, 'mantel_delta', 'p=float(np.mean(null')),
    dict(id='F4.e.ho.sel', figure='Fig. 4e', what='held-out Mantel rho, selected (mean over splits)',
         prefix='the gain held on held-out regional values (', source=f4('e_phylo_heldout', 'rho', pool='selected'),
         producer=(SEL, 'phylo_heldout', 'selected=')),
    dict(id='F4.e.ho.rand', figure='Fig. 4e', what='held-out Mantel rho, random (median)',
         prefix='the gain held on held-out regional values (', k=1,
         source=_random_median(F4, 'e_phylo_heldout', 'rho'), producer=(SEL, 'phylo_heldout', 'random=')),
    dict(id='F4.e.ho.p', figure='Fig. 4e', what='held-out Mantel, selected vs random P',
         prefix='the gain held on held-out regional values (', k=2,
         source=f4('e_phylo_heldout', 'p', pool='selected'), producer=(SEL, 'phylo_heldout', 'p=')),
    dict(id='F4.e.sel.p', figure='Fig. 4e', what='Mantel P of the selected runs against zero', prefix='(Mantel P = ',
         source=f4('e_phylo_mantel', 'p', profiles='raw', pool='selected'), producer=(PSTATS, 'figure4', 'e_phylo_mantel=')),

    # ═════════════ Fig. 6 -- FDG (V5 text; update the prefixes with the paragraph)
    dict(id='F6.n.merida', figure='Fig. 6', what='BIDS_Merida subjects with every region', prefix='BIDS_Merida, n = ',
         source=f6_n('Human_Merida'), producer=(PETD, 'dataset_matrix', "pivot_table(index=['dataset', 'subject']")),
    dict(id='F6.n.castrillon', figure='Fig. 6', what='BIDS_Castrillon subjects with every region',
         prefix='BIDS_Castrillon, n = ', source=f6_n('Human_Castrillon'),
         producer=(PETD, 'dataset_matrix', "pivot_table(index=['dataset', 'subject']")),
    dict(id='F6.n.tremblay', figure='Fig. 6', what='BIDS_Tremblay macaques', prefix='BIDS_Tremblay, n = ',
         source=f6_n('Macaque_Tremblay'), producer=(PETD, 'dataset_matrix', "pivot_table(index=['dataset', 'subject']")),
    dict(id='F6.b.r', figure='Fig. 6b', what='human datasets agreement, Pearson r of the regional means',
         prefix='agreed closely (Pearson r = ', source=dict(file=F6, sheet='C_agreement', column='pearson_r', where={}),
         producer=(PETD, 'agreement', "'pearson_r'")),
    dict(id='F6.b.ccc', figure='Fig. 6b', what="human datasets agreement, Lin's CCC",
         prefix="Lin's concordance correlation coefficient = ",
         source=dict(file=F6, sheet='C_agreement', column='lin_ccc', where={}), producer=(PETD, 'agreement', "'lin_ccc'")),
    dict(id='F6.g.vol', figure='Fig. 6g', what='cross-species Spearman r, volume profile',
         prefix='cross-species Spearman r of the mean profiles: volume ', source=f6_cons('morphometry'),
         producer=(A + 'library/multimodal_profiles.py', 'summary_stats', 'conservation')),
    dict(id='F6.g.fc', figure='Fig. 6g', what='cross-species Spearman r, positive coupling profile',
         prefix='cross-species Spearman r of the mean profiles: volume ', k=1, source=f6_cons('pos'),
         producer=(A + 'library/multimodal_profiles.py', 'summary_stats', 'conservation')),
    dict(id='F6.g.fdg', figure='Fig. 6g', what='cross-species Spearman r, FDG profile',
         prefix='cross-species Spearman r of the mean profiles: volume ', k=2, source=f6_cons('metabolism'),
         producer=(A + 'library/multimodal_profiles.py', 'summary_stats', 'conservation')),
]


# ═════════════ Methods statements checked against the code
METHODS = [
    dict(id='M.long.reference', quote="the subject's latest session is treated as a within-subject reference",
         kind='contains', file='modalities/anat/_loop2.py', func='run', snippet='if Session != max_ses:', expected='if Session != max_ses:'),
    dict(id='M.long.spline.df', quote='a cubic B-spline of age (4 d.f.)', kind='default',
         file='Statistics/Group_fMRI/EDNiX_longitudinal.py', func='fit_gamm', arg='spline_df', expected=4),
    dict(id='M.long.ml', quote='fitted by maximum likelihood', kind='contains',
         file='Statistics/Group_fMRI/EDNiX_longitudinal.py', func='fit_lme', snippet='reml=False', expected='reml=False'),
    dict(id='M.long.min.sessions', quote='Only subjects with ≥ 2 sessions were included', kind='constant',
         file=LONG, name='MIN_SESSIONS', expected=2),
    dict(id='M.long.level', quote='Isocortical (atlas level 1', kind='constant', file=LONG, name='ATLAS_LEVEL', expected=1,
         note='Fig. 3 a-f use atlas level 1 (Isocortex)'),
    dict(id='M.long.hemisphere', quote='left hemisphere', kind='default', file=LONG, func='figure3_stats',
         arg='hemisphere', expected='left',
         note='Fig. 3 a-f use the LEFT hemisphere only: the Methods / legend must say so (quote absent = not stated)'),
    dict(id='M.val.level', quote='summed and thickness averaged over the five isocortical lobes of the level-2 atlas',
         kind='contains', file=VALID, func='cross_sectional', snippet='load_morphometry(level=2)', expected='load_morphometry(level=2)',
         note='validation_morphometry.py reads level-1 Isocortex, not the sum of the five level-2 lobes'),
    dict(id='M.val.rcv', quote='the interquartile range divided by 1.349', kind='contains', file=VALID, func='_summary',
         snippet='/ 1.349', expected='/ 1.349'),
    dict(id='M.unit.sessions', quote='Morphometric values were averaged across sessions', kind='contains', file=PDATA,
         func='per_subject', snippet='mean', expected='mean'),
    dict(id='M.fs.thickness', quote='FreeSurfer averages the shortest distance from each surface to the other',
         kind='contains', file=FS, func='fs_measures', snippet="f'{h}.thickness'", expected="f'{h}.thickness'",
         note='FreeSurfer thickness = ?h.thickness (Fischl & Dale 2000)'),
    dict(id='M.s15.boot', quote='2,000 resamples of subjects', kind='constant', file=ABL, name='N_BOOT', expected=2000),
    # run selection (Fig. 4)
    dict(id='M.sel.nselect', quote='the 10 best are kept in each dataset', kind='contains', file=CFG_SEL, func=None,
         snippet='N_SELECT=10,', expected='N_SELECT=10,'),
    dict(id='M.sel.rounds', quote='for ten rounds', kind='contains', file=CFG_SEL, func=None,
         snippet='N_ITERATIVE_ROUNDS=10,', expected='N_ITERATIVE_ROUNDS=10,'),
    dict(id='M.sel.invariance', quote='26 of 32 regional values on which species agree most', kind='contains',
         file=CFG_SEL, func=None, snippet='INVARIANCE_VALUE=0.8,', expected='INVARIANCE_VALUE=0.8,',
         note='26 of 32 = the 80 % most invariant values'),
    dict(id='M.sel.scale', quote='divided by its root-mean-square, without centring', kind='contains', file=CFG_SEL,
         func=None, snippet='ZSCORE="scale"', expected='ZSCORE="scale"'),
    dict(id='M.sel.loso', quote='leave-one-species-out', kind='contains', file=CFG_SEL, func=None,
         snippet='LEAVE_ONE_OUT=True,', expected='LEAVE_ONE_OUT=True,'),
    dict(id='M.sel.splits', quote='In each of 50 random splits', kind='constant', file=SEL, name='N_SPLITS', expected=50),
    dict(id='M.sel.draws', quote='against 200 draws of 10 runs', kind='constant', file=SEL, name='N_DRAWS', expected=200),
    # FDG (Fig. 6)
    dict(id='M.pet.minvox', quote='regions with fewer than 10 PET voxels', kind='constant', file=MMF, name='MIN_VOXELS',
         expected=10),
    dict(id='M.pet.perm', quote='10,000 permutations; Hedges', kind='contains', file=MMF, func=None,
         snippet="TEST, N_PERM, ALPHA = 'permutation', 10000, 0.05", expected="TEST, N_PERM, ALPHA = 'permutation', 10000, 0.05"),
    dict(id='M.pet.boot', quote='5,000 bootstrap resamples of the subjects', kind='default', file=PETD, func='agreement',
         arg='n_boot', expected=5000),
    dict(id='M.pet.castrillon.window', quote='mean of the frames 20–45 min after injection', kind='contains',
         file='Exemples/Study_EDNiX/Launcher_pet/Launcher_fdg_Castrillon.py', func=None,
         snippet='SUV_window = [20, 45]', expected='SUV_window = [20, 45]'),
    dict(id='M.pet.tremblay.window', quote='dynamic acquisition under propofol, mean of all frames', kind='contains',
         file='Exemples/Study_EDNiX/Launcher_pet/Launcher_fdg_Macaque_imagina.py', func=None,
         snippet='SUV_window = None', expected='SUV_window = None'),
]
