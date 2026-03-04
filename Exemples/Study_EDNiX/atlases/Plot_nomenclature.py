"""
EDNiX Comparative Surface Atlas — Dual Hierarchical Sunburst (matplotlib)
Run:  python ednix_sunburst_matplotlib.py
Out:  ednix_sunburst.png  /  ednix_sunburst.svg
"""
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import matplotlib.patheffects as pe
# Color data from EDNIxCSCLR_l.ctab
color_data = [
    # Level 2
    ["Frontal_lobe",          2, "Frontal lobe",         220,  50,  50, "", "cortical"],
    ["Insular_lobe", 2, "Insular lobe", 180, 100, 200, "", "cortical"],
    ["Temporal_lobe", 2, "Temporal lobe", 50, 100, 220, "", "cortical"],
    ["Paleocortex", 2, "Paleocortex", 150, 150, 200, "", "cortical"],
    ["Archicortex", 2, "Archicortex", 130, 130, 130, "", "cortical"],
    ["Cortical_subplate", 2, "Cortical subplate", 200, 150, 100, "", "cortical"],
    ["Periarchicortex", 2, "Periarchicortex", 140, 170, 140, "", "cortical"],
    ["Parietal_lobe",         2, "Parietal lobe",        255, 215,   0, "", "cortical"],
    ["Occipital_lobe",        2, "Occipital lobe",        50, 180,  50, "", "cortical"],
    ["Cerebral_nuclei",       2, "Cerebral nuclei",      140,  70, 140, "", "subcortical"],
    ["Diencephalon",          2, "Diencephalon",          38,  72,  67, "", "subcortical"],  # Changed to 3V green
    ["Brainstem",             2, "Brainstem",            255, 215,   0, "", "subcortical"],  # Changed to yellow
    ["Cerebellum", 2, "Cerebellum", 162, 167, 169, "", "subcortical"],
    ["Cerebellum_White", 2, "Cerebellum White", 160, 160, 160, "", "subcortical"],  # Changed to gray
    ["White_matter",          2, "White matter",         192, 192, 192, "", "subcortical"],  # Changed to light gray
    ["Lateral_ventricules",   2, "Lateral ventricules",  128, 128, 128, "", "subcortical"],  # Changed to medium gray
    ["3rde_ventricules",      2, "3rde ventricules",     150, 150, 150, "", "subcortical"], # Changed to medium gray
    ["CSF",                   2, "CSF",                   98,   0, 103, "", "subcortical"],

    # Level 3
    ["oFC",              3, "oFC (orb. frontal Cx)",     200,  40,  40, "Frontal_lobe", "cortical"],
    ["oPFC",             3, "oPFC (orb. PFC)",           190,  45,  45, "Frontal_lobe", "cortical"],
    ["dlPFC",            3, "dlPFC",                     210,  50,  50, "Frontal_lobe", "cortical"],
    ["mPFC",             3, "mPFC",                      180,  55,  55, "Frontal_lobe", "cortical"],
    ["Motor_premotor",   3, "Motor/premotor",            170,  60,  60, "Frontal_lobe", "cortical"],
    ["vlPFC",            3, "vlPFC",                     160,  65,  65, "Frontal_lobe", "cortical"],
    ["Somatosensory_Cx", 3, "Somatosensory Cx",         240, 200,  30, "Parietal_lobe", "cortical"],
    ["Post_parietal_Cx", 3, "Post. parietal Cx",        245, 190,  20, "Parietal_lobe", "cortical"],
    ["PMC",              3, "PMC",                       250, 180,  10, "Parietal_lobe", "cortical"],
    ["Striate_Cx",       3, "Visual striate Cx",          30, 150,  30, "Occipital_lobe", "cortical"],  # Updated to match reference
    ["Extra_striate_Cx", 3, "Visual extra-striate Cx",   60, 170,  60, "Occipital_lobe", "cortical"],
    ["Auditory_Cx",      3, "Auditory Cx",                30,  85, 235, "Temporal_lobe", "cortical"],  # Updated to match reference
    ["MIPT",             3, "MIPT",                       50, 110, 230, "Temporal_lobe", "cortical"],
    ["Vent_temporal",    3, "Vent. temporal",             60, 130, 220, "Temporal_lobe", "cortical"],
    ["Insula",           3, "Insula & lat. sulcus",      170,  90, 210, "Insular_lobe", "cortical"],
    ["Olfactory_Cx",     3, "Olfactory Cx",              140, 140, 190, "Paleocortex", "cortical"],
    ["Olfactory_bulb",   3, "Olfactory bulb",            150, 150, 180, "Paleocortex", "cortical"],
    ["Hippocampal_form", 3, "Hippocampal form.",         120, 120, 120, "Archicortex", "cortical"],
    ["Subiculum",        3, "Subiculum",                 125, 125, 125, "Archicortex", "cortical"],
    ["Periarchicortex_L3",3,"Periarchicortex",           135, 165, 135, "Periarchicortex", "cortical"],
    ["Septum",           3, "Septum",                    130,  60, 130, "Cerebral_nuclei", "subcortical"],
    ["Striatum",         3, "Striatum",                  145,  60, 145, "Cerebral_nuclei", "subcortical"],
    ["Basal_forebrain",  3, "Basal forebrain",           135,  70, 135, "Cerebral_nuclei", "subcortical"],
    ["Pallidum",         3, "Pallidum",                  125,  80, 125, "Cerebral_nuclei", "subcortical"],
    ["Claustrum",        3, "Claustrum",                 195, 145,  95, "Cortical_subplate", "cortical"],
    ["Amygdala",         3, "Amygdala",                  190, 140, 110, "Cortical_subplate", "cortical"],
    ["Hypothalamus",     3, "Hypothalamus",               38,  72,  67, "Diencephalon", "subcortical"],  # Changed to 3V green
    ["Thalamus",         3, "Thalamus",                   70, 100,  80, "Diencephalon", "subcortical"],  # Lighter green gradient
    ["Epithalamus",      3, "Epithalamus",                90, 120,  90, "Diencephalon", "subcortical"],  # Lighter green gradient
    ["Subthalamus",      3, "Subthalamus",               110, 140, 100, "Diencephalon", "subcortical"],  # Lighter green gradient
    ["Tectum",           3, "Tectum",                    255, 215,   0, "Brainstem", "subcortical"],  # Yellow
    ["Tegmentum",        3, "Tegmentum",                 255, 200,  50, "Brainstem", "subcortical"],  # Yellow-orange gradient
    ["Cerebral_ped",     3, "Cerebral ped.",             255, 185,  75, "Brainstem", "subcortical"],  # Yellow-orange gradient
    ["Pons",             3, "Pons",                      255, 170, 100, "Brainstem", "subcortical"],  # Yellow-orange gradient
    ["Medulla",          3, "Medulla",                   245, 155, 125, "Brainstem", "subcortical"],  # Yellow-orange gradient
    ["Colliculus",       3, "Colliculus",                255, 200,  25, "Brainstem", "subcortical"],  # Yellow

    # Level 4 - cortical
    ["Gustatory_Cx",    4, "Gustatory Cx",              195,  50, 40, "oFC", "cortical"],
    ["oProiso",         4, "oProiso",                   180,  40,  40, "oFC", "cortical"],
    ["BA_10",           4, "BA 10",                     185,  45,  45, "oPFC", "cortical"],
    ["BA_11",           4, "BA 11",                     190,  50,  50, "oPFC", "cortical"],
    ["BA_13",           4, "BA 13",                     195,  55,  55, "oPFC", "cortical"],
    ["BA_12_47",        4, "BA 12/47",                  200,  60,  60, "oPFC", "cortical"],
    ["IFJ",             4, "IFJ",                       205,  50,  50, "dlPFC", "cortical"],
    ["BA_9",            4, "BA 9",                      210,  55,  55, "dlPFC", "cortical"],
    ["BA_46",           4, "BA 46",                     215,  60,  60, "dlPFC", "cortical"],
    ["BA_8",            4, "BA 8",                      200,  65,  65, "dlPFC", "cortical"],
    ["BA_55",           4, "BA 55",                     195,  70,  70, "dlPFC", "cortical"],
    ["BA_32",           4, "BA 32",                     185,  60,  60, "mPFC", "cortical"],
    ["BA_14",           4, "BA 14",                     180,  65,  65, "mPFC", "cortical"],
    ["BA_24",           4, "BA 24",                     175,  70,  70, "mPFC", "cortical"],
    ["BA_25",           4, "BA 25",                     160, 40, 35, "mPFC", "cortical"],
    ["BA_33",           4, "BA 33",                     165, 50, 40, "mPFC", "cortical"],
    ["BA_6",            4, "BA 6",                      170,  55,  55, "Motor_premotor", "cortical"],
    ["BA_4_M1",         4, "BA 4 (M1)",                 165,  65,  50, "Motor_premotor", "cortical"],
    ["BA_47",           4, "BA 47",                     190,  70,  70, "vlPFC", "cortical"],
    ["BA_45",           4, "BA 45",                     185,  75,  75, "vlPFC", "cortical"],
    ["BA_44",           4, "BA 44",                     180,  80,  80, "vlPFC", "cortical"],
    ["BA_1_3_S1",       4, "BA 1-3 (S1)",               235, 195,  40, "Somatosensory_Cx", "cortical"],
    ["BA_43",           4, "BA 43",                     240, 190,  35, "Somatosensory_Cx", "cortical"],
    ["S2",              4, "S2",                        245, 185,  30, "Somatosensory_Cx", "cortical"],
    ["BA_5",            4, "BA 5",                      250, 180,  25, "Post_parietal_Cx", "cortical"],
    ["BA_40",           4, "BA 40",                     245, 175,  20, "Post_parietal_Cx", "cortical"],
    ["BA_7",            4, "BA 7",                      240, 170,  15, "Post_parietal_Cx", "cortical"],
    ["BA_39",           4, "BA 39",                     235, 165,  10, "Post_parietal_Cx", "cortical"],
    ["PGM",             4, "PGM",                       230, 175,   5, "PMC", "cortical"],
    ["BA_23",           4, "BA 23",                     225, 180,  15, "PMC", "cortical"],
    ["BA_7m",           4, "BA 7m",                     235, 185,  20, "PMC", "cortical"],
    ["BA_31",           4, "BA 31",                     240, 190,  25, "PMC", "cortical"],
    ["BA_17_V1",        4, "BA 17 (V1)",                 30, 150,  30, "Striate_Cx", "cortical"],
    ["BA_18_V2",        4, "BA 18 (V2)",                 45, 160,  45, "Extra_striate_Cx", "cortical"],
    ["BA_19",           4, "BA 19",                      60, 170,  60, "Extra_striate_Cx", "cortical"],
    ["BA_41_42_A1",     4, "BA 41/42 (A1)",              30,  85, 235, "Auditory_Cx", "cortical"],
    ["Belt_A2",         4, "Belt (A2)",                  40,  95, 230, "Auditory_Cx", "cortical"],
    ["Parabelt_A3",     4, "Parabelt (A3)",              50, 105, 225, "Auditory_Cx", "cortical"],
    ["RTp",             4, "RTp",                        60, 115, 220, "Auditory_Cx", "cortical"],
    ["BA_22",           4, "BA 22",                      70, 125, 215, "Auditory_Cx", "cortical"],
    ["BA_21",           4, "BA 21",                      55, 115, 225, "MIPT", "cortical"],
    ["BA_20",           4, "BA 20",                      65, 125, 220, "MIPT", "cortical"],
    ["BA_38",           4, "BA 38",                      75, 135, 215, "MIPT", "cortical"],
    ["Parahippocampal", 4, "Parahippocampal",            85, 145, 210, "Vent_temporal", "cortical"],
    ["BA_37",           4, "BA 37",                      95, 155, 205, "Vent_temporal", "cortical"],
    ["Ri",              4, "Ri",                        160,  95, 200, "Insula", "cortical"],
    ["Insular_Cx",      4, "Insular Cx",                175, 100, 205, "Insula", "cortical"],
    ["BA_52",           4, "BA 52",                     165, 105, 195, "Insula", "cortical"],
    ["BA_27",           4, "BA 27",                     145, 145, 190, "Olfactory_Cx", "cortical"],
    ["Olfactory_tubercle",4,"Olfactory tubercle",       150, 150, 185, "Olfactory_Cx", "cortical"],
    ["OB",              4, "Olfactory bulb (OB)",       155, 155, 180, "Olfactory_bulb", "cortical"],
    ["CA",              4, "CA",                        115, 115, 115, "Hippocampal_form", "cortical"],
    ["DG",              4, "DG",                        120, 120, 120, "Hippocampal_form", "cortical"],
    ["BA_28",           4, "BA 28",                     125, 125, 125, "Hippocampal_form", "cortical"],
    ["BA_35_36",        4, "BA 35/36",                  135, 165, 135, "Periarchicortex_L3", "cortical"],
    ["Presubiculum",    4, "Presubiculum",              130, 160, 130, "Periarchicortex_L3", "cortical"],
    ["BA_26_29_30",     4, "BA 26/29/30",               135, 155, 145, "Periarchicortex_L3", "cortical"],

    # Level 4 - subcortical
    ["Caudate",         4, "Caudate",                   150,  65, 150, "Striatum", "subcortical"],
    ["Putamen",         4, "Putamen",                   145,  70, 145, "Striatum", "subcortical"],
    ["N_accumbens",     4, "N. accumbens",              140,  75, 140, "Striatum", "subcortical"],
    ["Globus_pallidus", 4, "Globus pallidus",           115,  80, 115, "Pallidum", "subcortical"],
    ["Bed_n_ST",        4, "Bed n. ST",                 200, 160,  100, "Amygdala", "cortical"],  # Updated to match reference
    ["Pulvinar",        4, "Pulvinar",                   90, 120,  90, "Thalamus", "subcortical"],  # Green gradient
    ["Geniculate_n",    4, "Geniculate n.",             110, 140, 100, "Thalamus", "subcortical"],  # Green gradient
    ["Thalamic_nuclei", 4, "Thalamic nuclei",           130, 160, 110, "Thalamus", "subcortical"],  # Green gradient
    ["Sup_colliculus",  4, "Sup. colliculus",           255, 215,   0, "Tectum", "subcortical"],  # Yellow
    ["Inf_colliculus",  4, "Inf. colliculus",           255, 200,  50, "Tectum", "subcortical"],  # Yellow-orange
    ["PAG",             4, "PAG",                       255, 185,  75, "Tegmentum", "subcortical"],  # Yellow-orange
    ["Red_n",           4, "Red n.",                    255, 170, 100, "Tegmentum", "subcortical"],  # Yellow-orange
    ["Subst_nigra",     4, "Subst. nigra",              245, 155, 125, "Tegmentum", "subcortical"],  # Yellow-orange
    ["Crus_cerebri",    4, "Crus cerebri",              255, 185,  75, "Cerebral_ped", "subcortical"],  # Yellow-orange
]

ABBREV = {
    "Paleocortex": "Paleo", "White matter": "WM",
    "Archicortex": "Archi", "Periarchicortex": "Peri",
    "Cerebral nuclei": "CN", "Cortical subplate": "CS",
    "Diencephalon": "Dien.", "Cerebellum": "Cbl.",
    "Brainstem": "BS", "Lateral ventricules": "LV",
    "Cerebellum White": "CbW", "3rde ventricules": "3V",
    "Motor/premotor": "M/PM", "Somatosensory Cx": "SmCx",
    "Post. parietal Cx": "PPC", "Visual striate Cx": "V1",
    "Visual extra-striate Cx": "V2+", "Auditory Cx": "AC",
    "Vent. temporal": "VT", "Olfactory Cx": "OC",
    "Olfactory bulb": "OB", "Hippocampal form.": "HF",
    "Basal forebrain": "BF", "Cerebral ped.": "CP",
    "Hypothalamus": "Hyp.", "Epithalamus": "Epi.",
    "Subthalamus": "Sub.", "Gustatory Cx": "GC",
    "BA 1-3 (S1)": "S1", "BA 4 (M1)": "M1",
    "BA 41/42 (A1)": "A1", "Belt (A2)": "A2",
    "Parabelt (A3)": "A3", "Parahippocampal": "PHC",
    "Insular Cx": "InsCx", "N. accumbens": "NAc",
    "Globus pallidus": "GP", "Bed n. ST": "BST",
    "Geniculate n.": "GN", "Thalamic nuclei": "TN",
    "Sup. colliculus": "SC", "Inf. colliculus": "IC",
    "Subst. nigra": "SN", "Crus cerebri": "CC",
    "BA 17 (V1)": "V1", "BA 18 (V2)": "V2",
    "BA 12/47": "12/47", "BA 26/29/30": "26-30",
    "BA 35/36": "35/36", "Insula & lat. sulcus": "Ins&LS",
    "oFC (orb. frontal Cx)": "oFC", "oPFC (orb. PFC)": "oPFC",
    "Olfactory bulb (OB)": "OB", "Olfactory tubercle": "OT",
}

df = pd.DataFrame(color_data, columns=["ID","Level","Label","R","G","B","Parent","Type"])
df['Label'] = df['Label'].str.replace('_', ' ')

def propagate_missing(df):
    extra = []
    for level in [2, 3]:
        at = df[df['Level'] == level]
        below = df[df['Level'] == level + 1]
        for _, row in at.iterrows():
            if not (below['Parent'] == row['ID']).any():
                extra.append([row['ID']+f'_L{level+1}', level+1,
                               row['Label'], row['R'], row['G'], row['B'],
                               row['ID'], row['Type']])
    if extra:
        return pd.concat([df, pd.DataFrame(extra, columns=df.columns)], ignore_index=True)
    return df

df = propagate_missing(df)
df = propagate_missing(df)

# Split into cortical and subcortical
df_cortical = df[df['Type'] == 'cortical'].copy()
df_subcortical = df[df['Type'] == 'subcortical'].copy()

def build_tree(df_subset):
    nodes = {row['ID']: {'label': row['Label'], 'level': row['Level'],
                         'rgb': (row['R'], row['G'], row['B']),
                         'parent': row['Parent'], 'children': []}
             for _, row in df_subset.iterrows()}
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
    return nodes, roots

nodes_cort, roots_cort = build_tree(df_cortical)
nodes_sub, roots_sub = build_tree(df_subcortical)

R_BANDS = {2: (0.22, 0.52), 3: (0.54, 0.80), 4: (0.82, 1.13)}
def lum(r, g, b): return (0.299*r + 0.587*g + 0.114*b) / 255

# Create dual sunburst figure
fig = plt.figure(figsize=(32, 18), facecolor='white')

# Left sunburst (cortical)
ax1 = fig.add_axes([0.02, 0.12, 0.5, 0.8])
ax1.set_aspect('equal'); ax1.axis('off')
ax1.set_xlim(-1.2, 1.2); ax1.set_ylim(-1.2, 1.2)

# Right sunburst (subcortical)
ax2 = fig.add_axes([0.5, 0.12, 0.5, 0.8])
ax2.set_aspect('equal'); ax2.axis('off')
ax2.set_xlim(-1.2, 1.2); ax2.set_ylim(-1.2, 1.2)

FS = {2: 13, 3: 13, 4: 12}
MIN_DA = {2: 0.00, 3: 0.04, 4: 0.055}
used_abbrevs = {}

def draw_sunburst(ax, nodes, roots, title):
    wedge_data = []
    def collect(nid_list, level, a0, a_total):
        tot = sum(nodes[n]['leaves'] for n in nid_list)
        a = a0
        for nid in nid_list:
            n = nodes[nid]
            da = (n['leaves'] / tot) * a_total
            wedge_data.append({'nid': nid, 'level': level,
                               'a0': a, 'a1': a+da, 'amid': a+da/2, 'da': da})
            if n['children'] and level < 4:
                collect(n['children'], level+1, a, da)
            a += da
    collect(roots, 2, np.pi/2, 2*np.pi)

    # Draw wedges
    for w in wedge_data:
        r, g, b = nodes[w['nid']]['rgb']
        r0, r1 = R_BANDS[w['level']]
        ax.add_patch(Wedge((0,0), r1, np.degrees(w['a0']), np.degrees(w['a1']),
                           width=r1-r0, facecolor=(r/255, g/255, b/255),
                           edgecolor='white', linewidth=0.7, zorder=w['level']))

    # Draw labels
    for w in wedge_data:
        n = nodes[w['nid']]
        lv, da, mid = w['level'], w['da'], w['amid']
        r, g, b = n['rgb']
        fg = 'white' if lum(r, g, b) < 0.35 else '#111111'
        label = n['label']
        if da < MIN_DA[lv]: continue
        r0, r1 = R_BANDS[lv]
        rm = (r0 + r1) / 2
        wx, wy = rm * np.cos(mid), rm * np.sin(mid)
        char_w, arc_len = 0.052, da * rm
        abbr = ABBREV.get(label, label)
        if len(label) * char_w <= arc_len * 0.92:
            disp = label
        else:
            disp = abbr
            if abbr != label: used_abbrevs[abbr] = label
        txt = ax.text(wx, wy, disp, ha='center', va='center',
                      fontsize=FS[lv], color=fg, fontweight='bold',
                      fontfamily='DejaVu Sans', zorder=lv+10)

    # Center
    ax.add_patch(plt.Circle((0,0), 0.21, color='white', zorder=30))
    ax.text(0, 0, title, ha='center', va='center', fontsize=25,
            fontweight='bold', color='#222', zorder=35)

draw_sunburst(ax1, nodes_cort, roots_cort, 'Cortical\nareas')
draw_sunburst(ax2, nodes_sub, roots_sub, 'Non-cortical\nareas')

# Main title
fig.text(0.5, 0.9, 'EDNiX Comparative Surface Atlas — Hierarchical Parcellation',
         ha='center', fontsize=18, fontweight='bold', color='#111')
fig.text(0.5, 0.88, 'LVL2: Lobar territories  ·  LVL3: Functional areas  ·  LVL4: Cytoarchitectonic areas',
         ha='center', fontsize=14, color='#555', fontstyle='italic', fontweight='bold')

# Legend
ax_leg = fig.add_axes([0.02, 0.02, 0.98, 0.10])
ax_leg.axis('off'); ax_leg.set_xlim(0, 1); ax_leg.set_ylim(0, 1)
ax_leg.text(0.02, 0.9, 'Abbreviations', fontsize=16, fontweight='bold',
            color='#111', va='top')
ax_leg.plot([0.15, 0.98], [0.85, 0.85], color='#bbb', lw=2)
entries = sorted(used_abbrevs.items())
if entries:
    n_cols, n_rows = 6, (len(entries) + 5) // 6
    col_width, row_height = 0.15, 0.10
    for i, (abbr, full) in enumerate(entries):
        col, row = i // n_rows, i % n_rows
        x_pos, y_pos = 0.05 + col * col_width, 0.62 - row * row_height
        ax_leg.text(x_pos, y_pos+0.02, abbr, fontsize=13, fontweight='bold',
                    color='#222', ha='left', va='center')
        ax_leg.text(x_pos+0.05, y_pos-0.02, full, fontsize=13, color='#444',
                    ha='left', va='center')

plt.savefig('ednix_sunburst.png', dpi=200, bbox_inches='tight',
            facecolor='white', pad_inches=0.15)
plt.savefig('ednix_sunburst.svg', format='svg', bbox_inches='tight',
            facecolor='white', pad_inches=0.15)
print(f"Saved! {len(used_abbrevs)} abbreviations in legend.")