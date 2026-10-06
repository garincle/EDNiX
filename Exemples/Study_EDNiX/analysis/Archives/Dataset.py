import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import matplotlib.transforms as transforms
import re

# ── Load ───────────────────────────────────────────────────────────────────────
df = pd.read_excel('/home/common/benhalab/CASCAD/EDNiX/databse_ednix_study.xlsx')
df.columns = df.columns.str.strip()
df = df.rename(columns={'Modality ': 'Modality'})

# ── Colours ────────────────────────────────────────────────────────────────────
COL_AWAKE  = '#27ae60'
COL_ANESTH = '#2980b9'
modality_colors = {'rs-fMRI': '#c0392b', 'PET FDG': '#b7770d'}

taxon_map = {
    'Homo sapiens':          ('Primate',    0, 0),
    'Macaca mulatta':        ('Primate',    0, 1),
    'Callithrix jacchus':    ('Primate',    0, 2),
    'Microcebus murinus':    ('Primate',    0, 3),
    'Canis lupus':           ('Carnivore',  1, 0),
    'Rattus norvegicus':     ('Rodent',     2, 0),
    'Mus musculus':          ('Rodent',     2, 1),
    'Phyllostomus discolor': ('Chiroptera', 3, 0),
}
taxon_band  = {'Primate':'#D6EAF8','Carnivore':'#FADBD8','Rodent':'#D5F5E3','Chiroptera':'#E8DAEF'}
taxon_label = {'Primate':'#1a5276','Carnivore':'#922b21','Rodent':'#1e8449','Chiroptera':'#6c3483'}

df['Taxon']      = df['Species'].map(lambda s: taxon_map.get(s,('Other',99,0))[0])
df['taxon_rank'] = df['Species'].map(lambda s: taxon_map.get(s,('Other',99,0))[1])
df['sp_rank']    = df['Species'].map(lambda s: taxon_map.get(s,('Other',99,0))[2])
df = df.sort_values(['taxon_rank','sp_rank','datasetname']).reset_index(drop=True)

def parse_n(val):
    if pd.isna(val): return None
    nums = [int(x) for x in re.findall(r'\d+', str(val))]
    return nums[-1] if nums else None
df['N_val'] = df['N'].apply(parse_n)

def parse_anesth(val):
    s = str(val).strip().rstrip('\xa0').lower()
    if s == 'awake':
        return [('awake', COL_AWAKE)]
    elif 'awake' in s and ('iso' in s or '/' in s):
        return [('awake', COL_AWAKE), ('isoflurane', COL_ANESTH)]
    else:
        label = str(val).strip().rstrip('\xa0')
        return [(label, COL_ANESTH)]
df['anesth_pills'] = df['Anesthesia'].apply(parse_anesth)

def parse_modality(val):
    s = str(val).strip()
    mods = []
    if 'rs-fMRI' in s: mods.append('rs-fMRI')
    if 'PET FDG' in s: mods.append('PET FDG')
    return mods if mods else [s]
df['mod_list'] = df['Modality'].apply(parse_modality)

def short_cite(paper):
    if pd.isna(paper): return '—'
    s = str(paper).strip()
    if s in ['not published','Unpublished dataset']: return s
    m_author = re.match(r'([A-Z][a-z]+)', s)
    m_year   = re.search(r'(20\d{2}|19\d{2})', s)
    parts = s.split('.')
    title_part = parts[1].strip() if len(parts) > 1 else ''
    words = title_part.split()[:9]
    title_short = ' '.join(words) + ('…' if len(title_part.split()) > 9 else '')
    author = m_author.group(1) if m_author else ''
    year   = m_year.group(1)   if m_year   else ''
    if author and year:
        return f"{author} et al., {year} — {title_short}"
    return s[:75] + '…'
df['cite'] = df['Paper'].apply(short_cite)

# ══════════════════════════════════════════════════════════════════════════════
#  LAYOUT
# ══════════════════════════════════════════════════════════════════════════════
# Key idea: row height = content height + padding
# We compute the natural content height first, then set ROW_H to fit

FS_HEAD  = 8.5    # column headers
FS_BODY  = 9.0    # dataset, species, N, reference
FS_PILL  = 8.0    # pill text
FS_TAXON = 8.5    # taxon group label
FS_LEG   = 7.5    # legend

# Pill dimensions — compact
PILL_H  = 0.24    # inches
PILL_A  = 0.50    # opacity
PILL_GAP= 0.12    # gap between side-by-side pills

# Row height: enough for a pill + small margin
ROW_H   = PILL_H + 0.22
HEAD_H  = 0.42
FOOT_H  = 0.42
NROWS   = len(df)

# Left margin = taxon bar width
MARGIN_L = 0.50
MARGIN_R = 0.10
FIG_W    = 15
FIG_H    = HEAD_H + NROWS * ROW_H + FOOT_H + 0.10

fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')
ax  = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, FIG_W)
ax.set_ylim(0, FIG_H)
ax.axis('off')

def row_yc(i):
    top = FIG_H - HEAD_H
    return top - (i + 0.5) * ROW_H
def row_ytop(i): return row_yc(i) + ROW_H/2
def row_ybot(i): return row_yc(i) - ROW_H/2

head_yc = FIG_H - HEAD_H/2

# Columns (x_left, width) in inches — auto-sized to content
col_defs = {
    'dataset':  (MARGIN_L,        2),
    'species':  (MARGIN_L+2,   2.5),
    'N':        (MARGIN_L+4.5,   0.55),
    'anesth':   (MARGIN_L+5,   2.10),
    'modality': (MARGIN_L+7,   2.10),
    'ref':      (MARGIN_L+5,   FIG_W - MARGIN_L - 2 - MARGIN_R),
}
def col_xc(k): return col_defs[k][0] + col_defs[k][1]/2
def col_xl(k): return col_defs[k][0]
def col_w(k):  return col_defs[k][1]

# ── Primitives ─────────────────────────────────────────────────────────────────
def draw_pill(x, y, w, h, color, alpha=PILL_A, zorder=4):
    pad = min(h*0.30, 0.05)
    p = FancyBboxPatch((x-w/2, y-h/2), w, h,
                       boxstyle=f'round,pad={pad}',
                       facecolor=color, edgecolor='none',
                       alpha=alpha, transform=ax.transData, zorder=zorder)
    ax.add_patch(p)

def pill_w(label, fs):
    # Empirical: ~0.055 inch/char at fs=8
    return len(label) * (fs/8) * 0.055 + 0.16

def text_pill(x, y, label, bg, fs=FS_PILL, max_w=None, ph=PILL_H):
    w = pill_w(label, fs)
    if max_w: w = min(w, max_w)
    draw_pill(x, y, w, ph, bg)
    ax.text(x, y, label, ha='center', va='center',
            fontsize=fs, color='#111', fontweight='bold',
            fontfamily='Liberation Serif',
            transform=ax.transData, zorder=5)
    return w

def draw_bg(x0, y0, x1, y1, color, alpha=1.0, zorder=1, r=0.04):
    p = FancyBboxPatch((x0,y0), x1-x0, y1-y0,
                       boxstyle=f'round,pad={r}',
                       facecolor=color, edgecolor='none',
                       alpha=alpha, transform=ax.transData, zorder=zorder)
    ax.add_patch(p)

def side_by_side(xc, yc, items, col_width, fs=FS_PILL, ph=PILL_H):
    """Pills side-by-side, centred in column."""
    widths = [pill_w(lbl, fs) for lbl, _ in items]
    max_w_each = (col_width - PILL_GAP*(len(items)-1) - 0.08) / len(items)
    widths = [min(w, max_w_each) for w in widths]
    total = sum(widths) + PILL_GAP*(len(items)-1)
    x0 = xc - total/2
    for (label, color), w in zip(items, widths):
        cx = x0 + w/2
        draw_pill(cx, yc, w, ph, color)
        ax.text(cx, yc, label, ha='center', va='center',
                fontsize=fs, color='#111', fontweight='bold',
                fontfamily='Liberation Serif',
                transform=ax.transData, zorder=5)
        x0 += w + PILL_GAP

# ══════════════════════════════════════════════════════════════════════════════
#  HEADER ROW
# ══════════════════════════════════════════════════════════════════════════════
draw_bg(0.08, head_yc-HEAD_H/2+0.025, FIG_W-0.08, head_yc+HEAD_H/2-0.025,
        '#1a1a2e', r=0.06, zorder=2)
for key, label in [('dataset','Dataset'),('species','Species'),('N','N'),
                   ('anesth','Anesthesia'),('modality','Modality'),('ref','Reference')]:
    ax.text(col_xc(key), head_yc, label,
            ha='center', va='center', fontsize=FS_HEAD, fontweight='bold',
            color='white', fontfamily='Liberation Serif',
            transform=ax.transData, zorder=5)

# ══════════════════════════════════════════════════════════════════════════════
#  TAXON BANDS + LEFT-MARGIN LABELS
# ══════════════════════════════════════════════════════════════════════════════
taxon_groups = {}
for i, row in df.iterrows():
    t = row['Taxon']
    taxon_groups.setdefault(t, {'start': i, 'end': i})
    taxon_groups[t]['end'] = i

BAR_X0, BAR_X1 = 0.08, 0.44   # taxon label bar in left margin

for taxon, info in taxon_groups.items():
    y_top = row_ytop(info['start']) + 0.01
    y_bot = row_ybot(info['end'])   - 0.01
    tc    = taxon_label[taxon]

    # Full-width band (light)
    draw_bg(BAR_X1, y_bot, FIG_W-0.08, y_top,
            taxon_band[taxon], alpha=0.40, r=0.04)

    # Left margin coloured bar
    draw_bg(BAR_X0, y_bot+0.03, BAR_X1-0.02, y_top-0.03,
            tc, alpha=0.15, r=0.04, zorder=2)

    # Taxon label — rotated, centred in bar
    ax.text((BAR_X0+BAR_X1)/2 - 0.01, (y_top+y_bot)/2, taxon,
            ha='center', va='center',
            fontsize=FS_TAXON, fontweight='bold', color=tc, rotation=90,
            fontfamily='Liberation Serif',
            transform=ax.transData, zorder=4)

# Separator lines
for i in range(1, NROWS):
    taxon_break = df.iloc[i]['Taxon'] != df.iloc[i-1]['Taxon']
    lw  = 1.0 if taxon_break else 0.35
    col = '#aaa' if taxon_break else '#e0e0e0'
    ax.plot([BAR_X1, FIG_W-0.10], [row_ytop(i)]*2,
            color=col, lw=lw, transform=ax.transData, zorder=3)

# ══════════════════════════════════════════════════════════════════════════════
#  DATA ROWS
# ══════════════════════════════════════════════════════════════════════════════
for i, row in df.iterrows():
    yc = row_yc(i)

    ax.text(col_xc('dataset'), yc, row['datasetname'],
            ha='center', va='center', fontsize=FS_BODY,
            color='#1a1a2e', fontfamily='Liberation Mono',
            transform=ax.transData, zorder=4)

    ax.text(col_xc('species'), yc, row['Species'],
            ha='center', va='center', fontsize=FS_BODY,
            color='#222', style='italic',
            fontfamily='Liberation Serif',
            transform=ax.transData, zorder=4)

    n = row['N_val']
    ax.text(col_xc('N'), yc, str(n) if n is not None else '—',
            ha='center', va='center', fontsize=FS_BODY+1,
            color='#1a1a2e', fontweight='bold',
            transform=ax.transData, zorder=4)

    anesth_items = row['anesth_pills']
    if len(anesth_items) == 1:
        text_pill(col_xc('anesth'), yc,
                  anesth_items[0][0], anesth_items[0][1],
                  max_w=col_w('anesth')-0.08)
    else:
        side_by_side(col_xc('anesth'), yc, anesth_items, col_w('anesth'))

    mods = row['mod_list']
    mod_items = [(m, modality_colors.get(m,'#546E7A')) for m in mods]
    if len(mod_items) == 1:
        text_pill(col_xc('modality'), yc,
                  mod_items[0][0], mod_items[0][1],
                  max_w=col_w('modality')-0.08)
    else:
        side_by_side(col_xc('modality'), yc, mod_items, col_w('modality'))

    ax.text(col_xl('ref') + col_w('ref')/2, yc, row['cite'],
            ha='center', va='center', fontsize=FS_BODY-1, color='#444',
            style='italic', fontfamily='Liberation Serif',
            transform=ax.transData, zorder=4)

# ══════════════════════════════════════════════════════════════════════════════
#  LEGEND — compact, bottom strip
# ══════════════════════════════════════════════════════════════════════════════
LEGY  = FOOT_H / 2
LPH   = 0.20    # legend pill height — smaller than body pills
LFS   = FS_LEG

draw_bg(0.08, 0.06, FIG_W-0.08, FOOT_H-0.06, '#F4F6F7', r=0.05, zorder=1)

lx = 0.38
ax.text(lx, LEGY, 'Anesthesia:', ha='left', va='center',
        fontsize=LFS+0.5, fontweight='bold', color='#333',
        fontfamily='Liberation Serif', transform=ax.transData, zorder=5)
lx += 1.05

for label, col in [('awake', COL_AWAKE), ('anesthetised', COL_ANESTH)]:
    pw = pill_w(label, LFS)
    draw_pill(lx+pw/2, LEGY, pw, LPH, col, zorder=4)
    ax.text(lx+pw/2, LEGY, label, ha='center', va='center',
            fontsize=LFS, color='#111', fontweight='bold',
            fontfamily='Liberation Serif', transform=ax.transData, zorder=5)
    lx += pw + 0.12

lx += 0.18
ax.plot([lx, lx], [LEGY-0.12, LEGY+0.12], color='#ccc', lw=0.7,
        transform=ax.transData, zorder=3)
lx += 0.20

ax.text(lx, LEGY, 'Modality:', ha='left', va='center',
        fontsize=LFS+0.5, fontweight='bold', color='#333',
        fontfamily='Liberation Serif', transform=ax.transData, zorder=5)
lx += 0.88

for label, col in [('rs-fMRI', modality_colors['rs-fMRI']),
                   ('PET FDG', modality_colors['PET FDG'])]:
    pw = pill_w(label, LFS)
    draw_pill(lx+pw/2, LEGY, pw, LPH, col, zorder=4)
    ax.text(lx+pw/2, LEGY, label, ha='center', va='center',
            fontsize=LFS, color='#111', fontweight='bold',
            fontfamily='Liberation Serif', transform=ax.transData, zorder=5)
    lx += pw + 0.12

# ── Outer frame ────────────────────────────────────────────────────────────────
outer = FancyBboxPatch((0.05, 0.04), FIG_W-0.10, FIG_H-0.08,
                        boxstyle='round,pad=0.04',
                        facecolor='none', edgecolor='#bbb', linewidth=0.8,
                        transform=ax.transData, zorder=6)
ax.add_patch(outer)

# ── Save ───────────────────────────────────────────────────────────────────────
for ext, dpi in [('pdf',150),('png',260)]:
    path = f'/scratch2/EDNiX/ednix_table.{ext}'
    fig.savefig(path, dpi=dpi, bbox_inches='tight', facecolor='white')
    print(f'Saved {path}')