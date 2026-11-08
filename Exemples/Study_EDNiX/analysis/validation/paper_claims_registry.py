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

Coverage: Fig. 2 (validation paragraph), Fig. 3 / Supp. Fig. 7 (longitudinal), Supp. Fig. 14
(FreeSurfer), Supp. Fig. 15 (longitudinal vs cross-sectional). Figures 4-6 are not
registered yet: their numbers appear in the 'untraced' sheet of the report.
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
S7 = 'Supp_Fig_07/Supp_Fig_07_stats.xlsx'
VAL = 'Validation/validation_morphometry.xlsx'
F2 = 'Figure_2/Figure_2_stats.xlsx'
XS = 'Supplementary_stats/extra_stats.txt'
S14 = 'Supp_Fig_14/Supp_Fig_14_stats.xlsx'
S15 = 'Supp_Fig_15/Supp_Fig_15_stats.xlsx'


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

    # ═════════════ Supp. Fig. 14 -- EDNiX vs FreeSurfer, human (BIDS_Park); PR #77
    dict(id='S14.n', figure='Supp. Fig. 14', what='subjects with both FreeSurfer and EDNiX values', prefix='(latest session; n = ',
         source=s14('Surface area', 'n'), producer=(FS, 'agreement', 'return dict(n=n'), data=['BIDS_Park']),
    ] + [
    dict(id=f'S14.{tag}.{q}', figure='Supp. Fig. 14', what=f'{meas}: {q}', prefix=pre, k=k, within='FreeSurfer recon-all',
         source=s14(meas, col), producer=(FS, 'agreement', key))
    for meas, tag, pre in (('Surface area', 'area', 'agreement between EDNiX and FreeSurfer was r = '),
                           ('Volume', 'vol', ') for surface area, r = '),
                           ('Thickness', 'thick', ') for grey-matter volume and r = '))
    for q, k, col, key in (('Pearson r', 0, 'pearson_r', 'pearson_r=r'), ('ICC(C,1)', 1, 'icc_c1', 'icc_c1=ic'))
    ] + [
    dict(id=f'S14.{tag}.bias', figure='Supp. Fig. 14', what=f'{meas}: Bland-Altman bias (%)',
         prefix='(EDNiX - FreeSurfer, % of the mean) of ', k=k, source=s14(meas, 'bias_pct'),
         producer=(FS, 'agreement', 'bias_pct=bias'))
    for k, (meas, tag) in enumerate((('Surface area', 'area'), ('Volume', 'vol'), ('Thickness', 'thick')))
    ] + [
    dict(id='S14.thick.vertexdef.bias', figure='Supp. Fig. 14', what='thickness bias with the EDNiX definition on FreeSurfer surfaces (%)',
         prefix='on the FreeSurfer surfaces, the thickness offset was ',
         source=s14('Thickness (FreeSurfer surfaces, EDNiX definition)', 'bias_pct'), producer=(FS, 'fs_measures', 'thick_v.append')),

    # ═════════════ Supp. Fig. 15 -- longitudinal vs cross-sectional processing
    dict(id='S15.ref.diff', figure='Supp. Fig. 15', what='reference session |difference|, max over species x modality (%)',
         prefix='gave the same values (median absolute difference ', source=s15_range('ref_abs_diff_pct_median', 'max'),
         producer=(ABL, 'compare_modes', "out['ref_abs_diff_pct_median']"), data=['BIDS_Park', 'BIDS_Zhu_Garin']),
    dict(id='S15.sigma.min', figure='Supp. Fig. 15g', what='residual SD ratio cross/long, min', prefix='residual standard deviation was ',
         source=s15_range('sigma_e_ratio', 'min'), producer=(ABL, 'compare_modes', "out['sigma_e_ratio'] =")),
    dict(id='S15.sigma.max', figure='Supp. Fig. 15g', what='residual SD ratio cross/long, max', prefix='residual standard deviation was ', k=1,
         source=s15_range('sigma_e_ratio', 'max'), producer=(ABL, 'compare_modes', "out['sigma_e_ratio'] =")),
    dict(id='S15.neq.min', figure='Supp. Fig. 15', what='sample-size equivalence (SE ratio²), min', prefix='would need ',
         within='cross-sectional mode', source=s15_range('n_equivalent_ratio', 'min'), producer=(ABL, 'compare_modes', "out['n_equivalent_ratio']")),
    dict(id='S15.neq.max', figure='Supp. Fig. 15', what='sample-size equivalence (SE ratio²), max', prefix='would need ', k=1,
         within='cross-sectional mode', source=s15_range('n_equivalent_ratio', 'max'), producer=(ABL, 'compare_modes', "out['n_equivalent_ratio']")),
    dict(id='S15.delta.pmin', figure='Supp. Fig. 15h', what='slope difference long - cross, smallest bootstrap P',
         prefix='slope difference between modes, smallest P = ', source=s15_range('delta_slope_p', 'min'),
         producer=(ABL, 'compare_modes', "out[f'{name}_p'] = _boot_p")),
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
]
