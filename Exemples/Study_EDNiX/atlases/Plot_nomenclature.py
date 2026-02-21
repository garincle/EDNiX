"""
EDNiX Comparative Surface Atlas — Hierarchical Sunburst (matplotlib)
Run:  python ednix_sunburst_matplotlib.py
Out:  ednix_sunburst.png  /  ednix_sunburst.svg
"""
import numpy as np
import pandas as pd
import matplotlib;

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import matplotlib.patheffects as pe

color_data = [
    ["Frontal L.", 2, "Frontal L.", 221, 0, 0, ""],
    ["Parietal L.", 2, "Parietal L.", 230, 230, 0, ""],
    ["Occipital L.", 2, "Occipital L.", 0, 230, 0, ""],
    ["Temporal L.", 2, "Temporal L.", 0, 179, 246, ""],
    ["Insular L.", 2, "Insular L.", 100, 179, 46, ""],
    ["Paleocortex", 2, "Paleocortex", 10, 50, 200, ""],
    ["Archicortex", 2, "Archicortex", 20, 20, 20, ""],
    ["Periarchicortex", 2, "Periarchicortex", 0, 189, 246, ""],
    ["Lateral_ventricules", 2, "Lateral ventricules", 98, 99, 103, ""],
    ["Cerebral_nuclei", 2, "Cerebral nuclei", 111, 0, 222, ""],
    ["Cortical_subplate", 2, "Cortical subplate", 210, 133, 0, ""],
    ["Diencephalon", 2, "Diencephalon", 10, 20, 0, ""],
    ["White_matter", 2, "White matter", 230, 230, 0, ""],
    ["Cerebellum", 2, "Cerebellum", 162, 167, 169, ""],
    ["Cerebellum_White", 2, "Cerebellum White", 98, 120, 103, ""],
    ["Brainstem", 2, "Brainstem", 169, 67, 127, ""],
    ["3rde_ventricules", 2, "3rde ventricules", 38, 72, 67, ""],
    ["CSF", 2, "CSF", 98, 0, 103, ""],
    ["oFC", 3, "oFC", 90, 0, 0, "Frontal L."],
    ["oPFC", 3, "oPFC", 190, 0, 0, "Frontal L."],
    ["dlPFC", 3, "dlPFC", 150, 0, 0, "Frontal L."],
    ["mPFC", 3, "mPFC", 130, 0, 0, "Frontal L."],
    ["Motor_premotor", 3, "Motor/premotor", 110, 0, 0, "Frontal L."],
    ["vlPFC", 3, "vlPFC", 170, 0, 0, "Frontal L."],
    ["Somatosensory_Cx", 3, "Somatosensory Cx", 230, 160, 0, "Parietal L."],
    ["Post_parietal_Cx", 3, "Post. parietal Cx", 230, 200, 0, "ParietalL. "],
    ["PMC", 3, "PMC", 229, 230, 0, "Parietal L."],
    ["Striate_Cx", 3, "Striate Cx", 70, 0, 0, "Occipital L."],
    ["Extra_striate_Cx", 3, "Extra striate Cx", 0, 70, 0, "Occipital L."],
    ["Auditory_Cx", 3, "Auditory Cx", 0, 30, 246, "Temporal L."],
    ["MIPT", 3, "MIPT", 0, 90, 246, "Temporal L."],
    ["Vent_temporal", 3, "Vent. temporal", 0, 150, 246, "Temporal L."],
    ["Insula", 3, "Insula", 0, 210, 246, "Insular L."],
    ["Septum", 3, "Septum", 113, 65, 150, "Cerebral_nuclei"],
    ["Olfactory_Cx", 3, "Olfactory Cx", 116, 116, 13, "Paleocortex"],
    ["Olfactory_bulb", 3, "Olfactory bulb", 170, 195, 154, "Paleocortex"],
    ["Hippocampal_form", 3, "Hippocampal form.", 116, 60, 13, "Archicortex"],
    ["Subiculum", 3, "Subiculum", 116, 60, 140, "Archicortex"],
    ["Periarchicortex_L3", 3, "Periarchicortex", 0, 189, 246, "Periarchicortex"],
    ["Striatum", 3, "Striatum", 111, 0, 196, "Cerebral_nuclei"],
    ["Basal_forebrain", 3, "Basal forebrain", 113, 65, 65, "Cerebral_nuclei"],
    ["Pallidum", 3, "Pallidum", 111, 0, 176, "Cerebral_nuclei"],
    ["Claustrum", 3, "Claustrum", 205, 133, 0, "Cortical_subplate"],
    ["Amygdala", 3, "Amygdala", 235, 133, 0, "Cortical_subplate"],
    ["Hypothalamus", 3, "Hypothalamus", 26, 128, 127, "Diencephalon"],
    ["Thalamus", 3, "Thalamus", 26, 176, 127, "Diencephalon"],
    ["Epithalamus", 3, "Epithalamus", 16, 30, 127, "Diencephalon"],
    ["Subthalamus", 3, "Subthalamus", 26, 30, 127, "Diencephalon"],
    ["Tectum", 3, "Tectum", 169, 67, 80, "Brainstem"],
    ["Tegmentum", 3, "Tegmentum", 169, 67, 100, "Brainstem"],
    ["Cerebral_ped", 3, "Cerebral ped.", 168, 67, 127, "Brainstem"],
    ["Pons", 3, "Pons", 205, 81, 149, "Brainstem"],
    ["Medulla", 3, "Medulla", 205, 81, 193, "Brainstem"],
    ["Colliculus", 3, "Colliculus", 169, 67, 127, "Brainstem"],
    ["Gustatory_Cx", 4, "Gustatory Cx", 171, 220, 154, "oFC"],
    ["oProiso", 4, "oProiso", 50, 0, 0, "oFC"],
    ["BA_10", 4, "BA 10", 120, 0, 0, "oPFC"],
    ["BA_11", 4, "BA 11", 220, 0, 0, "oPFC"],
    ["BA_13", 4, "BA 13", 230, 0, 0, "oPFC"],
    ["BA_12_47", 4, "BA 12/47", 240, 0, 0, "oPFC"],
    ["IFJ", 4, "IFJ", 40, 0, 0, "dlPFC"],
    ["BA_9", 4, "BA 9", 179, 0, 0, "dlPFC"],
    ["BA_46", 4, "BA 46", 140, 0, 0, "dlPFC"],
    ["BA_8", 4, "BA 8", 129, 0, 0, "dlPFC"],
    ["BA_55", 4, "BA 55", 160, 0, 0, "dlPFC"],
    ["BA_32", 4, "BA 32", 80, 0, 0, "mPFC"],
    ["BA_14", 4, "BA 14", 92, 0, 0, "mPFC"],
    ["BA_24", 4, "BA 24", 100, 0, 0, "mPFC"],
    ["BA_25", 4, "BA 25", 110, 10, 0, "mPFC"],
    ["BA_33", 4, "BA 33", 169, 0, 0, "mPFC"],
    ["BA_6", 4, "BA 6", 60, 0, 0, "Motor_premotor"],
    ["BA_4_M1", 4, "BA 4 (M1)", 70, 30, 0, "Motor_premotor"],
    ["BA_47", 4, "BA 47", 189, 0, 0, "vlPFC"],
    ["BA_45", 4, "BA 45", 200, 0, 0, "vlPFC"],
    ["BA_44", 4, "BA 44", 209, 0, 0, "vlPFC"],
    ["BA_1_3_S1", 4, "BA 1-3 (S1)", 229, 160, 0, "Somatosensory_Cx"],
    ["BA_43", 4, "BA 43", 230, 170, 0, "Somatosensory_Cx"],
    ["S2", 4, "S2", 230, 180, 0, "Somatosensory_Cx"],
    ["BA_5", 4, "BA 5", 230, 190, 0, "Post_parietal_Cx"],
    ["BA_40", 4, "BA 40", 229, 200, 0, "Post_parietal_Cx"],
    ["BA_7", 4, "BA 7", 230, 210, 0, "Post_parietal_Cx"],
    ["BA_39", 4, "BA 39", 230, 220, 0, "Post_parietal_Cx"],
    ["PGM", 4, "PGM", 180, 0, 0, "PMC"],
    ["BA_23", 4, "BA 23", 220, 230, 0, "PMC"],
    ["BA_7m", 4, "BA 7m", 230, 150, 0, "PMC"],
    ["BA_31", 4, "BA 31", 230, 250, 0, "PMC"],
    ["BA_17_V1", 4, "BA 17 (V1)", 0, 110, 0, "Striate_Cx"],
    ["BA_18_V2", 4, "BA 18 (V2)", 0, 160, 0, "Extra_striate_Cx"],
    ["BA_19", 4, "BA 19", 10, 230, 0, "Extra_striate_Cx"],
    ["BA_41_42_A1", 4, "BA 41/42 (A1)", 0, 20, 246, "Auditory_Cx"],
    ["Belt_A2", 4, "Belt (A2)", 0, 45, 246, "Auditory_Cx"],
    ["Parabelt_A3", 4, "Parabelt (A3)", 0, 60, 246, "Auditory_Cx"],
    ["RTp", 4, "RTp", 0, 75, 246, "Auditory_Cx"],
    ["BA_22", 4, "BA 22", 0, 80, 246, "Auditory_Cx"],
    ["BA_21", 4, "BA 21", 0, 105, 246, "MIPT"],
    ["BA_20", 4, "BA 20", 0, 120, 246, "MIPT"],
    ["BA_38", 4, "BA 38", 0, 135, 246, "MIPT"],
    ["Parahippocampal", 4, "Parahippocampal", 0, 180, 246, "Vent_temporal"],
    ["BA_37", 4, "BA 37", 10, 210, 246, "Vent_temporal"],
    ["Ri", 4, "Ri", 60, 210, 246, "Insula"],
    ["Insular_Cx", 4, "Insular Cx", 0, 225, 246, "Insula"],
    ["BA_52", 4, "BA 52", 0, 240, 246, "Insula"],
    ["CA", 4, "CA", 116, 59, 15, "Hippocampal_form"],
    ["DG", 4, "DG", 116, 60, 70, "Hippocampal_form"],
    ["BA_28", 4, "BA 28", 116, 60, 210, "Hippocampal_form"],
    ["BA_35_36", 4, "BA 35/36", 0, 165, 246, "Periarchicortex_L3"],
    ["Presubiculum", 4, "Presubiculum", 115, 60, 140, "Periarchicortex_L3"],
    ["BA_26_29_30", 4, "BA 26/29/30", 230, 240, 0, "Periarchicortex_L3"],
    ["Caudate", 4, "Caudate", 111, 0, 186, "Striatum"],
    ["Putamen", 4, "Putamen", 111, 0, 200, "Striatum"],
    ["N_accumbens", 4, "N. accumbens", 111, 0, 156, "Striatum"],
    ["Globus_pallidus", 4, "Globus pallidus", 113, 65, 113, "Pallidum"],
    ["Bed_n_ST", 4, "Bed n. stria term.", 115, 60, 13, "Amygdala"],
    ["Pulvinar", 4, "Pulvinar", 26, 50, 127, "Thalamus"],
    ["Geniculate_n", 4, "Geniculate n.", 26, 128, 0, "Thalamus"],
    ["Thalamic_nuclei", 4, "Thalamic nuclei", 0, 0, 127, "Thalamus"],
    ["Sup_colliculus", 4, "Sup. colliculus", 169, 57, 127, "Tectum"],
    ["Inf_colliculus", 4, "Inf. colliculus", 168, 67, 100, "Tectum"],
    ["PAG", 4, "PAG", 20, 67, 20, "Tegmentum"],
    ["Red_n", 4, "Red n.", 169, 100, 2, "Tegmentum"],
    ["Substantia_nigra", 4, "Substantia nigra", 50, 67, 50, "Tegmentum"],
    ["Crus_cerebri", 4, "Crus cerebri", 169, 200, 200, "Cerebral_ped"], ]

ABBREV = {
    "Paleocortex": "Paleo", "White matter": "WM",
    "Archicortex": "Archi", "Periarchicortex": "Peri", "Cerebral nuclei": "CN",
    "Cortical subplate": "CS", "Diencephalon": "Dien.", "Cerebellum": "Cbl.",
    "Brainstem": "BS", "Lateral ventricules": "LV", "Cerebellum White": "CbW",
    "3rde ventricules": "3V", "Cortical White matter": "CWM",
    "Motor/premotor": "M/PM", "Somatosensory Cx": "SmCx",
    "Post. parietal Cx": "PPC", "Striate Cx": "V1 Cx", "Extra striate Cx": "V2+",
    "Auditory Cx": "AC", "Vent. temporal": "VT", "Olfactory Cx": "OC",
    "Olfactory bulb": "OB", "Hippocampal form.": "HF", "Basal forebrain": "BF",
    "Cerebral ped.": "CP", "Hypothalamus": "Hyp.", "Epithalamus": "Epi.",
    "Subthalamus": "Sub.", "Gustatory Cx": "GC", "BA 1-3 (S1)": "S1",
    "BA 4 (M1)": "M1", "BA 41/42 (A1)": "A1", "Belt (A2)": "A2",
    "Parabelt (A3)": "A3", "Parahippocampal": "PHC", "Insular Cx": "InsCx",
    "N. accumbens": "NAc", "Globus pallidus": "GP", "Bed n. stria term.": "BST",
    "Geniculate n.": "GN", "Thalamic nuclei": "TN", "Sup. colliculus": "SC",
    "Inf. colliculus": "IC", "Substantia nigra": "SN", "Crus cerebri": "CC",
    "BA 17 (V1)": "V1", "BA 18 (V2)": "V2", "BA 12/47": "12/47",
    "BA 26/29/30": "26-30", "BA 35/36": "35/36", "Insular": "Ins"}

# Create a list of labels that should only appear at level 4
LEVEL4_ONLY_LABELS = []

df = pd.DataFrame(color_data, columns=["ID", "Level", "Label", "R", "G", "B", "Parent"])
df['Label'] = df['Label'].str.replace('_', ' ')


def propagate_missing(df):
    extra = []
    for level in [2, 3]:
        at = df[df['Level'] == level]
        below = df[df['Level'] == level + 1]
        for _, row in at.iterrows():
            # Skip propagation for labels that should only appear at level 4
            if row['Label'] in LEVEL4_ONLY_LABELS:
                continue
            if not (below['Parent'] == row['ID']).any():
                extra.append([row['ID'] + f'_L{level + 1}', level + 1,
                              row['Label'], row['R'], row['G'], row['B'], row['ID']])
    if extra:
        return pd.concat([df, pd.DataFrame(extra, columns=df.columns)], ignore_index=True)
    return df


df = propagate_missing(df)
df = propagate_missing(df)

nodes = {row['ID']: {'label': row['Label'], 'level': row['Level'],
                     'rgb': (row['R'], row['G'], row['B']),
                     'parent': row['Parent'], 'children': []}
         for _, row in df.iterrows()}
for nid, n in nodes.items():
    if n['parent'] and n['parent'] in nodes:
        nodes[n['parent']]['children'].append(nid)

roots = [nid for nid, n in nodes.items() if n['level'] == 2]


def count_leaves(nid):
    n = nodes[nid]
    if not n['children']: return 1
    return sum(count_leaves(c) for c in n['children'])


for nid in nodes:
    nodes[nid]['leaves'] = count_leaves(nid)

R_BANDS = {2: (0.22, 0.52), 3: (0.54, 0.80), 4: (0.82, 1.13)}


def lum(r, g, b): return (0.299 * r + 0.587 * g + 0.114 * b) / 255


# Create figure with sunburst taking most space
fig = plt.figure(figsize=(28, 22), facecolor='white')

# Sunburst axes - properly enlarged to fill space
ax = fig.add_axes([0.02, 0.08, 0.96, 0.86])  # [left, bottom, width, height]
ax.set_aspect('equal')
ax.axis('off')
# Expand the data limits to make the sunburst larger relative to the axes
scale_factor = 0.80  # Reduce this to make the sunburst appear larger
ax.set_xlim(-1.6 * scale_factor, 1.6 * scale_factor)
ax.set_ylim(-1.6 * scale_factor, 1.6 * scale_factor)

wedge_data = []


def collect(nid_list, level, a0, a_total):
    tot = sum(nodes[n]['leaves'] for n in nid_list)
    a = a0
    for nid in nid_list:
        n = nodes[nid]
        da = (n['leaves'] / tot) * a_total
        wedge_data.append({'nid': nid, 'level': level,
                           'a0': a, 'a1': a + da, 'amid': a + da / 2, 'da': da})
        if n['children'] and level < 4:
            collect(n['children'], level + 1, a, da)
        a += da


collect(roots, 2, np.pi / 2, 2 * np.pi)

for w in wedge_data:
    r, g, b = nodes[w['nid']]['rgb']
    r0, r1 = R_BANDS[w['level']]
    ax.add_patch(Wedge((0, 0), r1, np.degrees(w['a0']), np.degrees(w['a1']),
                       width=r1 - r0, facecolor=(r / 255, g / 255, b / 255),
                       edgecolor='white', linewidth=0.7, zorder=w['level']))

FS = {2: 12, 3: 10, 4: 9}
MIN_DA = {2: 0.00, 3: 0.04, 4: 0.055}
used_abbrevs = {}

for w in wedge_data:
    n = nodes[w['nid']]
    lv = w['level']
    da = w['da']
    mid = w['amid']
    r, g, b = n['rgb']
    fg = 'white' if lum(r, g, b) < 0.2 else '#111111'
    label = n['label']
    if da < MIN_DA[lv]:
        continue
    r0, r1 = R_BANDS[lv]
    rm = (r0 + r1) / 2
    wx = rm * np.cos(mid)
    wy = rm * np.sin(mid)
    char_w = 0.052
    arc_len = da * rm
    abbr = ABBREV.get(label, label)
    if len(label) * char_w <= arc_len * 0.92:
        disp = label
    else:
        disp = abbr
        if abbr != label:
            used_abbrevs[abbr] = label
    txt = ax.text(wx, wy, disp, ha='center', va='center',
                  fontsize=FS[lv], color=fg,
                  fontweight='bold',
                  fontfamily='DejaVu Sans', zorder=lv + 10)

ax.add_patch(plt.Circle((0, 0), 0.21, color='white', zorder=30))
ax.text(0, 0.055, 'EDNiX', ha='center', va='center', fontsize=18,
        fontweight='bold', color='#222', zorder=35)
ax.text(0, -0.075, 'Atlases', ha='center', va='center', fontsize=12,
        fontweight='bold', color='#666', zorder=35)

# Adjust title positions
ax.text(0, 1.3,
        'EDNiX Comparative Surface Atlas — Hierarchical Parcellation',
        ha='center', fontsize=16, fontweight='bold', color='#111')

ax.text(0, 1.2,
        'LVL2: Lobar territories  ·  LVL3: functional areas  ·  LVL4: cytoarchitectonic areas',
        ha='center', fontsize=15, color='#555', fontstyle='italic',
        fontweight='bold')

# Legend - Compact but readable at the bottom
ax2 = fig.add_axes([0.05, 0.02, 0.9, 0.12])  # Slightly taller for better spacing
ax2.axis('off')
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1)

# Legend header - increased space below
ax2.text(0.02, 0.80, 'Abbreviations', fontsize=15, fontweight='bold',
         color='#111', va='top', transform=ax2.transAxes)
ax2.plot([0.15, 0.98], [0.72, 0.72], color='#bbb', lw=2,
         transform=ax2.transAxes)

# Legend entries - tight grid with minimal spacing between abbreviation and full name
entries = sorted(used_abbrevs.items())
n_ent = len(entries)

if n_ent > 0:
    # Calculate grid layout with tighter spacing
    n_cols = 6  # More columns for tighter layout
    n_rows = (n_ent + n_cols - 1) // n_cols

    col_width = 0.15  # Narrower columns
    row_height = 0.10  # Tighter row spacing

    for i, (abbr, full) in enumerate(entries):
        col = i // n_rows
        row = i % n_rows

        x_pos = 0.05 + col * col_width
        y_pos = 0.62 - row * row_height

        # Abbreviation (bold) - closer to full name
        ax2.text(x_pos, y_pos + 0.02, abbr, fontsize=12, fontweight='bold',
                 color='#222', ha='left', va='center', transform=ax2.transAxes)
        # Full name - closer to abbreviation (reduced from 0.09 to 0.05)
        ax2.text(x_pos + 0.05, y_pos - 0.02, full, fontsize=11, color='#444',
                 ha='left', va='center', transform=ax2.transAxes)

plt.savefig('ednix_sunburst.png', dpi=200, bbox_inches='tight',
            facecolor='white', pad_inches=0.15)
plt.savefig('ednix_sunburst.svg', format='svg', bbox_inches='tight',
            facecolor='white', pad_inches=0.15)

print(f"Saved! {len(used_abbrevs)} abbreviations in legend.")