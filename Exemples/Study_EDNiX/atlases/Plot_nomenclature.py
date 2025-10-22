import pandas as pd
import plotly.express as px
import plotly.colors as pc
import colorsys

# Load data
file_path = '/srv/projects/easymribrain/data/Atlas/13_Atlas_project/Classiff/Legende_VDualvf2_formatrix.xlsx'
df = pd.read_excel(file_path, sheet_name="Legend_2023")
df = df.fillna('')

# Create display labels
def create_display_labels(row):
    labels = []
    prev_label = ""
    for col in ['NEWLVL1', 'NEWLVL2', 'NEWLVL3', 'NEWLVL4']:
        if row[col] and str(row[col]) != prev_label:
            labels.append(str(row[col]))
            prev_label = str(row[col])
    return '|'.join(labels)

df['display_label'] = df.apply(create_display_labels, axis=1)

# Define base colors for NEWLVL1
lvl1_names = df['NEWLVL1'].unique()
#base_colors = px.colors.qualitative.Plotly
'''
base_colors = [
    '#ff9896', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#98df8a', '#ffbb78', '#9edae5', '#aec7e8', '#c5b0d5',
    '#c49c94', '#f7b6d2', '#1f77b4', '#dbdb8d', '#c7c7c7']
base_color_map = {name: base_colors[i % len(base_colors)] for i, name in enumerate(lvl1_names)}
'''

base_color_map = {'CSF': '#ff9896', '3rde ventricules': '#ff7f0e', '4th ventricules': '#ffbb78', 'Lateral ventricules': '#2ca02c',
                  'Cerebellum White': '#9467bd', 'Cortical White ': '#8c564b', 'Cerebellum': '#e377c2', 'Brain stem': '#7f7f7f', 'Allocortex': '#bcbd22',
                  'Subcortical areas': '#17becf', 'Diencephalon': '#98df8a', 'Isocortex': '#d62728', 'Periallocortex ': '#9edae5'}

# Utility to lighten a color
def lighten_color(color, factor):
    rgb = pc.hex_to_rgb(color)
    hls = colorsys.rgb_to_hls(*[v/255 for v in rgb])
    lighter_rgb = colorsys.hls_to_rgb(hls[0], min(1, hls[1] + factor*(1 - hls[1])), hls[2])
    return pc.label_rgb([int(c * 255) for c in lighter_rgb])

# Assign hierarchical gradient color
def compute_color(row):
    lvl1 = row['NEWLVL1']
    depth = 0
    for col in ['NEWLVL2', 'NEWLVL3', 'NEWLVL4']:
        if row[col]: depth += 1
    base = base_color_map.get(lvl1, '#cccccc')
    return lighten_color(base, 0.15 * depth)

df['color'] = df.apply(compute_color, axis=1)

# Build the plot
fig = px.sunburst(
    df,
    path=['NEWLVL1', 'NEWLVL2', 'NEWLVL3', 'NEWLVL4'],
    color='color',
    color_discrete_map={k: k for k in df['color'].unique()},
    width=1400,
    height=1400,
    branchvalues='total',
    hover_data={'display_label': True}
)

# Customize layout
fig.update_layout(
    title={
        'y':0.95,
        'x':0.5,
        'xanchor': 'center',
        'yanchor': 'top',
        'font': {'size': 24, 'family': "Arial"}
    },
    margin=dict(t=150, l=20, r=20, b=20),
    uniformtext=dict(minsize=14, mode='hide'),
    plot_bgcolor='white'
)

fig.update_traces(
    textinfo="label",
    insidetextorientation='radial',
    marker=dict(
        line=dict(color='white', width=0.5),
    ),
    hovertemplate='<b>%{label}</b><br>Path: %{customdata[0]}<extra></extra>',
    textfont=dict(size=12)
)

# Save output
output_path = "/srv/projects/easymribrain/data/Atlas/13_Atlas_project/Classiff/NEWLVL_hierarchy_gradient.png"
fig.write_image(output_path, scale=2, engine="kaleido")
print(f"Visualisation avec gradient sauvegardée à : {output_path}")
