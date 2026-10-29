"""
validation_morphometry.py -- do EDNiX morphometric values agree with the literature, and is
the between-subject dispersion acceptable?  (numbers for the "validation" paragraphs of the
paper, Results: morphometry and longitudinal).

Writes <PAPER_DIR>/Validation/:
  validation_morphometry.xlsx  sheets
      isocortex_by_species   level-1 isocortex, bilateral, one value per subject:
                             n, median, q25, q75, robust CV (IQR / 1.349 / median)
      isocortex_by_dataset   same, per dataset (between-dataset agreement within species)
      consistency            per species: median of area x thickness / volume (1 = the three
                             measures are mutually consistent)
      literature             EDNiX median vs published value, % difference (sources inline)
      annual_change          per species x modality: annualised % change per subject
                             (within-subject OLS slope / first-session value x 100),
                             median [IQR] across subjects, and pooled (LME slope / mean
                             first-session value x 100), next to published rates
      herculano              rows of HERCULANO_XLSX matching the EDNiX species, columns that
                             look like cortical surface area or thickness (for a manual check)
  validation_morphometry.txt   the same, as sentences

    python validation_morphometry.py
"""
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_data as pdat

opj = os.path.join
OUT_DIR = opj(pdat.PAPER_DIR, 'Validation')
COLS = {'surface': 'surface_area_mm2', 'volume': 'volume_mm3', 'thickness': 'thickness_mm'}
UNITS = {'surface': 'mm2', 'volume': 'mm3', 'thickness': 'mm'}

# Published values. value in the unit of UNITS; "surface_kind" says which surface the
# reference measured (EDNiX area is summed on the PIAL surface, which is larger than the
# midthickness surface by construction).
LITERATURE = [
    dict(species='Human', modality='thickness', value=2.5,
         source='Fischl & Dale 2000 PNAS 97:11050 (overall average ~2.5 mm, range 1-4.5)'),
    dict(species='Human', modality='thickness', value=2.7,
         source='Hayashi et al. 2021 NeuroImage 229:117726 (median of the group thickness map)'),
    dict(species='Macaque', modality='thickness', value=2.4,
         source='Hayashi et al. 2021 NeuroImage 229:117726 (median of the group thickness map)'),
    dict(species='Marmoset', modality='thickness', value=1.6,
         source='Hayashi et al. 2021 NeuroImage 229:117726 (median of the group thickness map)'),
    dict(species='Mouse', modality='thickness', value=0.89,
         source='Lerch et al. 2008 NeuroImage 41:243 (ex vivo MRI, whole-cortex mean)'),
    dict(species='Human', modality='surface', value=2 * 97300.0,
         source='Van Essen et al. 2012 Cereb Cortex 22:2241 (973 cm2 per hemisphere, '
                'MIDTHICKNESS surface; EDNiX = pial)'),
]

# Published annual % change (older adults, Storsve et al. 2014 J Neurosci 34:8488, mean
# annual percentage change across the cortex, 23-87 y) for the human ageing comparison.
AGEING_RATES = {'surface': -0.19, 'thickness': -0.35, 'volume': -0.51}


def _summary(df, col, keys):
    g = df.groupby(keys)[col]
    out = pd.DataFrame({'n': g.size(), 'median': g.median(), 'q25': g.quantile(0.25),
                        'q75': g.quantile(0.75)}).reset_index()
    out['robust_cv_pct'] = 100 * (out['q75'] - out['q25']) / 1.349 / out['median']
    return out


def cross_sectional():
    morph = pdat.load_morphometry(level=1)
    by_sp, by_ds, iso = [], [], {}
    for mod, col in COLS.items():
        df = morph.get(mod, pd.DataFrame())
        if df is None or df.empty:
            print(f"  [WARN] no {mod} data")
            continue
        df = pdat.select_regions(df, ['Isocortex'])
        if df.empty:
            print(f"  [WARN] no Isocortex row for {mod}")
            continue
        iso[mod] = df
        by_sp.append(_summary(df, col, ['species']).assign(modality=mod, unit=UNITS[mod]))
        by_ds.append(_summary(df, col, ['species', 'bids_dir']).assign(modality=mod,
                                                                        unit=UNITS[mod]))
    by_sp = pd.concat(by_sp, ignore_index=True) if by_sp else pd.DataFrame()
    by_ds = pd.concat(by_ds, ignore_index=True) if by_ds else pd.DataFrame()

    cons = pd.DataFrame()
    if all(m in iso for m in COLS):
        keys = ['species', 'bids_dir', 'subject']
        m = (iso['surface'][keys + [COLS['surface']]]
             .merge(iso['thickness'][keys + [COLS['thickness']]], on=keys)
             .merge(iso['volume'][keys + [COLS['volume']]], on=keys))
        m['area_x_thickness_over_volume'] = (m[COLS['surface']] * m[COLS['thickness']]
                                             / m[COLS['volume']])
        cons = (m.groupby('species')['area_x_thickness_over_volume']
                .agg(n='size', median='median', q25=lambda s: s.quantile(0.25),
                     q75=lambda s: s.quantile(0.75)).reset_index())
    return by_sp, by_ds, cons


def literature_table(by_sp):
    rows = []
    for ref in LITERATURE:
        hit = by_sp[(by_sp['species'] == ref['species']) & (by_sp['modality'] == ref['modality'])]
        if hit.empty:
            continue
        med = float(hit['median'].iloc[0])
        rows.append(dict(species=ref['species'], modality=ref['modality'], unit=UNITS[ref['modality']],
                         ednix_median=med, ednix_q25=float(hit['q25'].iloc[0]),
                         ednix_q75=float(hit['q75'].iloc[0]), published=ref['value'],
                         diff_pct=100 * (med - ref['value']) / ref['value'], source=ref['source']))
    return pd.DataFrame(rows)


def annual_change():
    from Exemples.Study_EDNiX.analysis.library import longitudinal as lg
    from Statistics.Group_fMRI.EDNiX_longitudinal import prepare_longitudinal, fit_lme
    rows = []
    for sp, morph, age_df, age_col, _ in lg.load_figure3_data():
        for mod, col in COLS.items():
            df = morph.get(mod, pd.DataFrame()) if morph else pd.DataFrame()
            if df is None or df.empty:
                continue
            m = prepare_longitudinal(df, age_df, age_col, atlas_level=lg.ATLAS_LEVEL,
                                     hemisphere='left', regions=['Isocortex'],
                                     min_sessions=lg.MIN_SESSIONS, label=f'valid/{sp}_{mod}')
            if m.empty:
                continue
            m = m[m['region'] == 'Isocortex'].dropna(subset=[col, age_col])
            per_subj, base = [], []
            for _, g in m.groupby('subject'):
                g = g.sort_values(age_col)
                if g[age_col].nunique() < 2:
                    continue
                slope = np.polyfit(g[age_col].to_numpy(float), g[col].to_numpy(float), 1)[0]
                b = float(g[col].iloc[0])
                base.append(b)
                per_subj.append(100 * slope / b)
            if not per_subj:
                continue
            per_subj = np.asarray(per_subj)
            lme = fit_lme(m, col, age_col, lg.MIN_SESSIONS)
            # slope_per_unit = per year (fit_lme's 'slope' is per SD of age)
            lme_slope = float(lme.get('slope_per_unit', np.nan)) if isinstance(lme, dict) else np.nan
            rows.append(dict(
                species=sp, modality=mod, n_subjects=len(per_subj),
                n_sessions=int(len(m)), age_min=float(m[age_col].min()),
                age_max=float(m[age_col].max()),
                baseline_mean=float(np.mean(base)),
                apc_median=float(np.median(per_subj)),
                apc_q25=float(np.percentile(per_subj, 25)),
                apc_q75=float(np.percentile(per_subj, 75)),
                n_subjects_negative=int((per_subj < 0).sum()),
                lme_slope_per_year=lme_slope,
                apc_pooled=100 * lme_slope / float(np.mean(base)) if np.isfinite(lme_slope) else np.nan,
                published_apc_ageing=AGEING_RATES.get(mod) if sp == 'Human' else np.nan))
    return pd.DataFrame(rows)


def herculano():
    path = pdat.HERCULANO_XLSX
    if not os.path.isfile(path):
        print(f"  [WARN] {path} not found")
        return pd.DataFrame()
    names = {'Human': ['homo sapiens', 'human'], 'Macaque': ['macaca', 'macaque', 'rhesus'],
             'Marmoset': ['callithrix', 'marmoset'], 'Mouselemur': ['microcebus', 'mouse lemur'],
             'Mouse': ['mus musculus', 'mouse'], 'Rat': ['rattus', 'rat'],
             'Dog': ['canis', 'dog']}
    out = []
    for sheet, df in pd.read_excel(path, sheet_name=None).items():
        if df.empty:
            continue
        txt_cols = [c for c in df.columns if df[c].dtype == object]
        val_cols = [c for c in df.columns
                    if re.search(r'thick|area|surface|\bA_?T\b|\bT\b', str(c), re.I)]
        if not val_cols or not txt_cols:
            continue
        lab = df[txt_cols].astype(str).agg(' | '.join, axis=1).str.lower()
        for sp, keys in names.items():
            pat = '|'.join(r'\b' + re.escape(k) + r'\b' for k in keys)
            hit = df[lab.str.contains(pat, regex=True)]
            if hit.empty:
                continue
            h = hit[txt_cols[:2] + val_cols].copy()
            h.insert(0, 'ednix_species', sp)
            h.insert(0, 'sheet', sheet)
            out.append(h)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def _sentences(by_sp, lit, cons, apc):
    L = ['VALIDATION OF EDNiX MORPHOMETRY', '=' * 60, '',
         'Isocortex (atlas level 1), bilateral, one value per subject (sessions averaged).', '']
    for _, r in by_sp.sort_values(['modality', 'species']).iterrows():
        L.append(f"  {r['species']:<11} {r['modality']:<9} n={int(r['n']):>3}  median "
                 f"{r['median']:.4g} {r['unit']} [IQR {r['q25']:.4g}-{r['q75']:.4g}]  "
                 f"robust CV {r['robust_cv_pct']:.1f}%")
    if not by_sp.empty:
        L += ['', 'Robust CV range by modality:']
        for mod, g in by_sp.groupby('modality'):
            L.append(f"  {mod}: {g['robust_cv_pct'].min():.1f}-{g['robust_cv_pct'].max():.1f}% "
                     f"(median {g['robust_cv_pct'].median():.1f}%)")
    if not cons.empty:
        L += ['', 'Internal consistency (area x thickness / volume, median per species):']
        for _, r in cons.iterrows():
            L.append(f"  {r['species']:<11} {r['median']:.2f} [{r['q25']:.2f}-{r['q75']:.2f}]")
    if not lit.empty:
        L += ['', 'EDNiX vs literature:']
        for _, r in lit.iterrows():
            L.append(f"  {r['species']:<9} {r['modality']:<9} EDNiX {r['ednix_median']:.4g} vs "
                     f"{r['published']:.4g} {r['unit']} ({r['diff_pct']:+.0f}%) -- {r['source']}")
    if not apc.empty:
        L += ['', 'Annualised % change of left isocortex (per subject: OLS slope / first value):']
        for _, r in apc.iterrows():
            pub = (f"; published ageing {r['published_apc_ageing']:+.2f}%/y (Storsve 2014)"
                   if np.isfinite(r['published_apc_ageing']) else '')
            L.append(f"  {r['species']:<8} {r['modality']:<9} n={r['n_subjects']} "
                     f"({r['age_min']:.1f}-{r['age_max']:.1f} y): median {r['apc_median']:+.2f}%/y "
                     f"[IQR {r['apc_q25']:+.2f} to {r['apc_q75']:+.2f}], "
                     f"{r['n_subjects_negative']}/{r['n_subjects']} subjects negative; pooled "
                     f"{r['apc_pooled']:+.2f}%/y{pub}")
    return '\n'.join(L) + '\n'


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    by_sp, by_ds, cons = cross_sectional()
    lit = literature_table(by_sp) if not by_sp.empty else pd.DataFrame()
    try:
        apc = annual_change()
    except Exception as e:  # longitudinal inputs live on other paths; report, do not stop
        print(f"  [WARN] annual change not computed: {e}")
        apc = pd.DataFrame()
    try:
        hh = herculano()
    except Exception as e:
        print(f"  [WARN] Herculano-Houzel table not read: {e}")
        hh = pd.DataFrame()
    xlsx = opj(OUT_DIR, 'validation_morphometry.xlsx')
    with pd.ExcelWriter(xlsx) as w:
        for name, t in (('isocortex_by_species', by_sp), ('isocortex_by_dataset', by_ds),
                        ('consistency', cons), ('literature', lit), ('annual_change', apc),
                        ('herculano', hh)):
            (t if not t.empty else pd.DataFrame({'note': ['no data']})).to_excel(
                w, sheet_name=name, index=False)
    txt = _sentences(by_sp, lit, cons, apc)
    with open(opj(OUT_DIR, 'validation_morphometry.txt'), 'w', encoding='utf-8') as fh:
        fh.write(txt)
    print(txt)
    print(f"  [saved] {xlsx}")


if __name__ == '__main__':
    main()
