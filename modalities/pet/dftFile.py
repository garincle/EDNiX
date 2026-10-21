"""
DFT files (Turku PET Centre TAC format) of the regional time-activity curves.

Format (https://www.turkupetcentre.net/petanalysis/format_tpc_dft.html), TAB separated here:
    line 1: 'DFT' + region names
    line 2: study id (<= 8 characters) + hemisphere / second name ('.' = none)
    line 3: concentration unit + planes ('All')
    line 4: 'Times (sec)' (frame start and end) + region volumes in mm3 (= uL)
    then one line per frame: start, end, one value per region ('.' = missing)
"""
import numpy as np


def unit_label(unit):
    """Unit string accepted by tpcclib (kBq/mL, Bq/mL, nCi/mL...)."""
    if not unit:
        return 'kBq/mL'
    return str(unit).replace('cc', 'mL').replace('ml', 'mL')


def write(path, t_start, t_end, values, names, volumes_mm3, unit='kBq/mL', study='EDNiX'):
    values = np.atleast_2d(np.asarray(values, dtype=float))
    if values.shape[0] != len(t_start):
        values = values.T
    if values.shape != (len(t_start), len(names)):
        raise ValueError(f'DFT {path}: {values.shape} values for {len(t_start)} frames x {len(names)} regions')
    if len(volumes_mm3) != len(names):
        raise ValueError(f'DFT {path}: {len(volumes_mm3)} volumes for {len(names)} regions')
    rows = ['\t'.join(['DFT'] + [str(n) for n in names]),
            '\t'.join([str(study)[:8]] + ['.'] * len(names)),
            '\t'.join([unit_label(unit)] + ['All'] * len(names)),
            '\t'.join(['Times (sec)'] + [f'{float(v):.6g}' for v in volumes_mm3])]
    for j in range(len(t_start)):
        vals = ['.' if not np.isfinite(v) else f'{v:.6g}' for v in values[j]]
        rows.append('\t'.join([f'{float(t_start[j]):g}', f'{float(t_end[j]):g}'] + vals))
    with open(path, 'w') as f:
        f.write('\n'.join(rows) + '\n')
    return path


def read(path):
    """(frame starts, frame ends, values frames x regions, names, volumes) of a DFT written in 'Times' mode
    (weights columns, added by dftweigh, are dropped)."""
    with open(path) as f:
        lines = [l.rstrip('\n') for l in f if l.strip() and not l.startswith('#')]
    tab = '\t' in lines[0]
    split = (lambda l: l.split('\t')) if tab else (lambda l: l.split())
    names = split(lines[0])[1:]
    head4 = split(lines[3])
    if not head4[0].startswith('Times'):
        raise ValueError(f'{path}: only DFT files with frame start and end times are read')
    vols = head4[1:] if tab or not head4[1].startswith('(') else head4[2:]
    data = [split(l) for l in lines[4:]]
    keep = [i for i, n in enumerate(names) if n.lower() != 'weight']
    arr = np.array([[np.nan if x == '.' else float(x) for x in row[2:2 + len(names)]] for row in data])
    t_start = np.array([float(r[0]) for r in data])
    t_end = np.array([float(r[1]) for r in data])
    vol = np.array([float(v) if v not in ('.', '') else np.nan for v in vols[:len(names)]])
    return t_start, t_end, arr[:, keep], [names[i] for i in keep], vol[keep]


def roi(file, tac_df, timing, unit, study='EDNiX'):
    """DFT of every region of a TAC table (roi_signal.tac: columns = labels, attrs 'names', 'volumes_mm3')."""
    names = tac_df.attrs['names']
    return write(file, timing['start'], timing['end'], tac_df.values, names, tac_df.attrs['volumes_mm3'], unit, study)


def ref(file, data4d, labels_vol, labels, voxel_mm3, timing, unit, ref_name, study='EDNiX'):
    """
    Reference TAC = mean concentration over the union of the voxels of `labels` (the engineer averaged the
    regional means, which weighs a small region as much as a large one).
    """
    sel = np.isin(labels_vol, labels)
    if not sel.any():
        raise ValueError(f'reference region {ref_name}: none of the labels {labels} is in the atlas in PET space')
    tac = data4d[sel].mean(axis=0)
    return write(file, timing['start'], timing['end'], tac[:, None], [ref_name], [sel.sum() * voxel_mm3], unit, study), int(sel.sum())


def mask(file, data4d, msk, voxel_mm3, timing, unit, ref_name, study='EDNiX'):
    """Reference TAC = mean concentration in a mask (the engineer wrote the SUM of the activity)."""
    if not msk.any():
        raise ValueError(f'reference mask {ref_name} is empty in PET space')
    tac = data4d[msk].mean(axis=0)
    return write(file, timing['start'], timing['end'], tac[:, None], [ref_name], [msk.sum() * voxel_mm3], unit, study), int(msk.sum())
