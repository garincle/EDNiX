"""
Longitudinal analyses that are NOT panels of Figure 3 (sensitivity / exploration).

  - run_longitudinal per species x modality: <region>_trajectory.png, _pct_change.png,
    stats.xlsx under <OUT_DIR>/<Species>/<modality>/
  - stand-alone association-vs-primary PNGs (Fig3D_assoc_vs_primary_*), both species
    (the human ones are Supplementary Figure 6, drawn by
    library/human_assoc_primary.py)
  - longitudinal_summary.xlsx + longitudinal_report.txt
  - % change from baseline figure (former FigureS3_pct_change)

Figure 3 itself: library/longitudinal.py (run_paper.py).
"""
import os
import sys

import pandas as pd

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Statistics.Group_fMRI.EDNiX_longitudinal import run_longitudinal
from Exemples.Study_EDNiX.analysis.library.longitudinal import (
    OUT_DIR, ATLAS_LEVEL, REGIONS, MIN_SESSIONS, MACAQUE_BIDS_DIRS, MACAQUE_EXCEL,
    HUMAN_BIDS_DIRS, HUMAN_AGE_PREFIX, PRIMARY_REGIONS, ASSOCIATION_REGIONS,
    _build_atlas_paths, load_morphometry, format_macaque_age, format_human_age,
    run_fig3d_assoc_vs_primary, generate_report, plot_figure3_supp_pct_change)

opj = os.path.join

if __name__ == '__main__':
    os.makedirs(OUT_DIR, exist_ok=True)
    all_stats = []

    # ══════════════════════════════════════════════════════════════════════════
    # MACAQUE
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 64)
    print("  MACAQUE")
    print("=" * 64)

    atlas_paths, _ = _build_atlas_paths(['Macaque'])

    # Load ALL regions needed for both main and Fig 3D
    mac_morph = load_morphometry({'Macaque': MACAQUE_BIDS_DIRS}, atlas_paths,
                                  regions_of_interest=REGIONS)
    mac_age = format_macaque_age(MACAQUE_EXCEL)

    # Main longitudinal analysis - only surface for macaque
    for mod, metric_col, metric_label in [
        ('surface', 'surface_area_mm2', 'Surface area (mm²)'),
        ('thickness', 'thickness_mm', 'Cortical thickness (mm)'),
        ('volume', 'volume_mm3', 'Volume (mm³)'),
    ]:
        df_mod = mac_morph.get(mod, pd.DataFrame())
        if df_mod is None or df_mod.empty:
            print(f"  [Macaque/{mod}] no data")
            continue

        stats = run_longitudinal(
            label=f'Macaque_{mod}',
            df_morph=df_mod,
            df_age=mac_age,
            metric_col=metric_col,
            metric_label=metric_label,
            age_col='age_years',
            age_label='Age (years)',
            out_dir=opj(OUT_DIR, 'Macaque', mod),
            atlas_level=ATLAS_LEVEL,
            hemisphere='left',
            regions=['Isocortex'],
            min_sessions=MIN_SESSIONS,
        )
        if stats is not None and not stats.empty:
            stats['species'] = 'Macaque'
            stats['modality'] = mod
            all_stats.append(stats)

    # ── FIG 3D: Association vs Primary Cortex (Macaque) ──────────────────────
    print("\n" + "=" * 64)
    print("  FIG 3D — Association vs Primary Cortex (Macaque)")
    print("=" * 64)
    for data_key in ['surface', 'volume', 'thickness']:
        fig3d_stats, fig3d_path = run_fig3d_assoc_vs_primary(
            mac_morph, mac_age,
            out_dir=opj(OUT_DIR, 'Macaque', 'fig3d'),
            primary_regions=PRIMARY_REGIONS,
            association_regions=ASSOCIATION_REGIONS,
            species_label="Macaque",
            data_key=data_key,
        )
        if fig3d_stats:
            all_stats.append(pd.DataFrame([fig3d_stats]))

    # ══════════════════════════════════════════════════════════════════════════
    # HUMAN
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 64)
    print("  HUMAN")
    print("=" * 64)

    atlas_paths, _ = _build_atlas_paths(['Human'])

    hum_morph = load_morphometry({'Human': HUMAN_BIDS_DIRS}, atlas_paths,
                                  regions_of_interest=REGIONS)
    hum_age = format_human_age(HUMAN_BIDS_DIRS, age_prefix=HUMAN_AGE_PREFIX)

    # Main longitudinal analysis - all modalities for human
    for mod, metric_col, metric_label in [
        ('surface', 'surface_area_mm2', 'Surface area (mm²)'),
        ('thickness', 'thickness_mm', 'Cortical thickness (mm)'),
        ('volume', 'volume_mm3', 'Volume (mm³)'),
    ]:
        df_mod = hum_morph.get(mod, pd.DataFrame())
        if df_mod is None or df_mod.empty:
            print(f"  [Human/{mod}] no data")
            continue

        stats = run_longitudinal(
            label=f'Human_{mod}',
            df_morph=df_mod,
            df_age=hum_age,
            metric_col=metric_col,
            metric_label=metric_label,
            age_col='age',
            age_label='Age (years)',
            out_dir=opj(OUT_DIR, 'Human', mod),
            atlas_level=ATLAS_LEVEL,
            hemisphere='left',
            regions=['Isocortex'],
            min_sessions=MIN_SESSIONS,
        )
        if stats is not None and not stats.empty:
            stats['species'] = 'Human'
            stats['modality'] = mod
            all_stats.append(stats)

    # ── FIG 3D: Association vs Primary Cortex (Human) ──────────────────────────
    print("\n" + "=" * 64)
    print("  FIG 3D — Association vs Primary Cortex (Human)")
    print("=" * 64)
    for data_key in ['surface', 'thickness', 'volume']:
        fig3d_stats_human, fig3d_path_human = run_fig3d_assoc_vs_primary(
            hum_morph, hum_age,
            out_dir=opj(OUT_DIR, 'Human', 'fig3d'),
            primary_regions=PRIMARY_REGIONS,
            association_regions=ASSOCIATION_REGIONS,
            species_label="Human",
            data_key=data_key,
        )
        if fig3d_stats_human:
            all_stats.append(pd.DataFrame([fig3d_stats_human]))

    # ── Combined summary Excel ────────────────────────────────────────────────
    if all_stats:
        summary = pd.concat(all_stats, ignore_index=True, sort=False)
        out_xlsx = opj(OUT_DIR, 'longitudinal_summary.xlsx')
        summary.to_excel(out_xlsx, index=False)
        print(f"\n  Summary Excel → {out_xlsx}")

        # ── Generate clean report ──────────────────────────────────────────────
        generate_report(
            all_stats, OUT_DIR, mac_age, hum_age,
            mac_surface=mac_morph.get('surface') if mac_morph else None,
            hum_surface=hum_morph.get('surface') if hum_morph else None,
            mac_thickness=mac_morph.get('thickness') if mac_morph else None,
            hum_thickness=hum_morph.get('thickness') if hum_morph else None,
            mac_volume=mac_morph.get('volume') if mac_morph else None,
            hum_volume=hum_morph.get('volume') if hum_morph else None,
        )

    # ── % change from baseline (former Figure S3) ────────────────────────────
    fig3_rows = [
        ('Human',   hum_morph, hum_age, 'age',       'Human'),
        ('Macaque', mac_morph, mac_age, 'age_years', 'Macaque'),
    ]
    plot_figure3_supp_pct_change(fig3_rows, out_dir=opj(OUT_DIR, 'sensitivity'))

    print("\n  Done.")