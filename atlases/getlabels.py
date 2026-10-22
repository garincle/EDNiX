"""
Reference region of the kinetic models (SRTM, Logan) of the PET pipeline.

The region is given in the launcher (`reference_signal`) or chosen interactively (PET engineer's
windows, tkinter) when it is 'unknown'; the interactive choice is saved and reused for the other
subjects of the study.

reference_signal
----------------
{'type': 'atlas', 'atlas': 'EDNIxCSCLR', 'level': 0, 'labels': [7, 1007], 'name': 'Cer'}
{'type': 'atlas', 'atlas': 'EDNIxCSCLR', 'level': 0, 'labels': 'Cerebell', 'name': 'Cer'}   regular expression
{'type': 'atlas', 'atlas': 'EDNIxCSCLR', 'level': 0, 'labels': 'default',  'name': 'Cer'}   table cerebellum()
{'type': 'mask',  'name': 'white'}             white | csf | white+csf | gray | brain   (masks in PET space)
'unknown'  (or {'type': 'unknown'})            interactive choice
"""
import os
import re
import json

import numpy as np

opj = os.path.join
opi = os.path.isfile

MASK_NAMES = ['white', 'csf', 'white+csf', 'gray', 'brain']


def label_table(path_label_code, atlas_name):
    """label_id / region_name of an atlas (label files of <species>/label_code)."""
    from Tools.atlas_labels import parse_label_file
    for fn in (atlas_name + '_label.txt', atlas_name + '_StatsLUT.txt'):
        if opi(opj(path_label_code, fn)):
            return parse_label_file(opj(path_label_code, fn))
    raise FileNotFoundError(f'no label file for {atlas_name} in {path_label_code}')


def cerebellum(atlas):
    """Cerebellum labels of the atlases used with EDNiX (table of the PET engineer; the commas missing in
    ['BAL2013','BAL2020''BAL2023''BAL2024'] and the INIA19 entry stored in REF instead of ref are fixed)."""
    name, level = atlas['name'], atlas.get('level', 0)
    ref = None
    if name in ['BAL2013', 'BAL2020', 'BAL2023', 'BAL2024']:
        ref = [29, 30]
    elif name == 'D99v2':
        ref = [108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 125]
    elif name == 'CIVM':
        ref = [4]
        print('This atlas is not adapted: there is no subdivision between gray and white matter')
    elif name == 'INIA19':
        ref = [224, 225, 226, 227, 228, 229, 230, 231, 235, 236, 238, 240, 241, 242, 244, 245, 246, 247, 248, 250,
               252, 253, 504, 505,
               1224, 1225, 1226, 1227, 1228, 1229, 1230, 1231, 1235, 1236, 1238, 1240, 1241, 1242, 1244, 1245, 1246,
               1247, 1248, 1250, 1252, 1253, 1504, 1505]
    elif name == 'ONPRC18':
        ref = [57]
        print('This atlas is not adapted: there is no subdivision between gray and white matter')
    elif name in ['CHARM-SARM', 'SARM']:
        ref = {0: [203], 1: [240], 2: [244, 250], 3: [244, 251, 267, 275], 4: [244, 251, 267, 275, 246],
               5: [252, 253, 269, 270, 271, 272, 273, 274, 276, 277]}.get(level)
    elif name == 'EDNIxCSCLR':
        ref = [7, 1007]   # level 0 only
    elif name == 'EDNIxCSC':
        ref = [7]
    if ref is None:
        raise ValueError(f'no cerebellum definition for atlas {name} level {level}: give the labels')
    return ref


def find_labels(label_df, pattern, present=None):
    """Labels whose name matches the regular expression `pattern` (case-insensitive), restricted to `present`."""
    sel = label_df[label_df['region_name'].str.contains(pattern, flags=re.IGNORECASE, regex=True)]
    labels = [int(l) for l in sel['label_id']]
    if present is not None:
        labels = [l for l in labels if l in set(present)]
    return labels


# ------------------------------------------------------------------------------------------------------------------
# Interactive choice (PET engineer). tkinter is imported only here, so that the module imports on a server.
# ------------------------------------------------------------------------------------------------------------------
def getval(present_labels, result, label_df):
    """Window: multiple choice of the regions of reference (present in the atlas) and of their name."""
    import tkinter as tk
    from tkinter import ttk
    df = label_df[label_df['label_id'].isin(present_labels)]
    names_list = df['region_name'].tolist()

    root = tk.Tk()
    root.title('Select the region(s) of reference (ex: cerebellum)')
    tk.Label(root, text='Choose one or multiple regions:').grid(row=0, column=0)
    listbox = tk.Listbox(root, selectmode=tk.MULTIPLE, width=50, height=25)
    for name in names_list:
        listbox.insert(tk.END, name)
    listbox.grid(row=1, column=0)
    tk.Label(root, text='Name of the reference:').grid(row=0, column=1)
    ref_var = tk.StringVar(value='Cer')
    ttk.Entry(root, textvariable=ref_var).grid(row=1, column=1, sticky='n')

    def validation():
        chosen = [names_list[i] for i in listbox.curselection()]
        result['labels'] = sorted(int(l) for l in df[df['region_name'].isin(chosen)]['label_id'])
        result['name'] = ref_var.get()
        root.destroy()

    tk.Button(root, text='Validate', command=validation).grid(row=2, column=1)
    root.mainloop()
    return result


def _choose(title, names_list):
    import tkinter as tk
    out = {}
    root = tk.Tk()
    root.title(title)
    tk.Label(root, text='Choose one:').grid(row=0, column=0)
    listbox = tk.Listbox(root)
    for name in names_list:
        listbox.insert(tk.END, name)
    listbox.grid(row=1, column=0)
    listbox.select_set(0)

    def validation():
        out['v'] = listbox.get(listbox.curselection())
        root.destroy()

    tk.Button(root, text='Validate', command=validation).grid(row=2, column=1)
    root.mainloop()
    return out.get('v')


def getmask(result):
    result['name'] = _choose('Select the mask of reference (ex: white)', MASK_NAMES)
    return result


def gettype(atlas, present_labels, label_df):
    """Interactive: reference from the atlas (choice of regions) or from a tissue mask."""
    result = {'type': _choose('Select the reference type', ['atlas', 'mask']), 'name': None, 'labels': None}
    if result['type'] == 'atlas':
        result.update(atlas=atlas['name'], level=atlas.get('level', 0))
        result = getval(present_labels, result, label_df)
    elif result['type'] == 'mask':
        result = getmask(result)
    return result


def resolve(reference_signal, atlas, present_labels, label_df, saved_json):
    """
    Reference region as {'type', 'name', 'atlas', 'level', 'labels'} (atlas) or {'type': 'mask', 'name'}.
    atlas: {'name', 'level'} used for an interactive or 'default' choice; present_labels: labels of that
    level in PET space. The interactive choice is written to saved_json and reused when it exists.
    """
    ref = reference_signal
    if ref in (None, '', 'unknown') or (isinstance(ref, dict) and ref.get('type') in (None, 'unknown')):
        if saved_json and opi(saved_json):
            with open(saved_json) as f:
                return json.load(f)
        ref = gettype(atlas, present_labels, label_df)
        if not ref.get('type') or (ref['type'] == 'atlas' and not ref.get('labels')):
            raise ValueError('no reference region chosen')
        if saved_json:
            os.makedirs(os.path.dirname(saved_json), exist_ok=True)
            with open(saved_json, 'w') as f:
                json.dump(ref, f, indent=2)
        return ref
    ref = dict(ref)
    if ref['type'] == 'mask':
        if ref.get('name') not in MASK_NAMES:
            raise ValueError(f"reference mask {ref.get('name')!r}: expected one of {MASK_NAMES}")
        return ref
    if ref['type'] != 'atlas':
        raise ValueError(f"reference_signal type {ref['type']!r}: expected 'atlas', 'mask' or 'unknown'")
    ref.setdefault('atlas', atlas['name'])
    ref.setdefault('level', atlas.get('level', 0))
    ref.setdefault('name', 'ref')
    labels = ref.get('labels')
    if labels == 'default':
        ref['labels'] = cerebellum({'name': ref['atlas'], 'level': ref['level']})
    elif isinstance(labels, str):
        ref['labels'] = find_labels(label_df, labels, present_labels)
    else:
        ref['labels'] = [int(l) for l in np.atleast_1d(labels)]
    if not ref['labels']:
        raise ValueError(f'reference region {ref}: no label selected')
    return ref
