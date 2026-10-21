import ants
import pandas as pd
import tkinter as tk
from tkinter import ttk

def getval(filenii,result,level,filecsv):

    df = pd.read_csv(filecsv,sep='\t')

    img = ants.image_read(filenii)
    if img.dimension == 4:
        tmp = ants.ndimage_to_list(img)
        values= tmp[level].unique()
    else:
        values= img.unique()

    df_filtered = df[df['ROIs'].isin(values)]
    names_list  = df_filtered['Names'].tolist()

    # --- GUI ---
    root = tk.Tk()
    root.title("Select the relevant regions of reference (ex: Cerebellum")

    label = tk.Label(root, text="Choose one or multiple regions :")
    label.grid(row=0, column=0)
    #label.pack()

    # multi-selection Listbox
    listbox = tk.Listbox(root, selectmode=tk.MULTIPLE, width=50)
    for name in names_list:
        listbox.insert(tk.END, name)
    #listbox.pack()
    listbox.grid(row=1, column=0)

    ref_var = tk.StringVar(value='Cer')
    ref_menu = ttk.Entry(root, textvariable=ref_var)
    ref_menu.grid(row=1, column=1)

    # --- get ---
    def validation():
        selected_indices = listbox.curselection()
        selected_names = [names_list[i] for i in selected_indices]

        # Mapping name -> ROIs
        selected_rois = df_filtered[df_filtered['Names'].isin(selected_names)]['ROIs'].tolist()
        result["rois"] = sorted(selected_rois)
        result["Refname"] = ref_var.get()
        root.destroy()

    btn = tk.Button(root, text="Validate", command=validation)
    btn.grid(row=2, column=1)
    #btn.pack()

    root.mainloop()

    return result


def cerebellum(atlas):
    if atlas["name"] in ['BAL2013','BAL2020''BAL2023''BAL2024']:
        ref = [29,30]

    elif atlas["name"] == 'D99v2':
        ref = [108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 125]
    elif atlas["name"] == 'CIVM':
        ref = [4]
        print('This atlas is not adapted : there is no subdivision between gray and white matter')

    elif atlas["name"] == 'INIA19':
        REF = [224, 225, 226, 227, 228, 229, 230, 231, 235, 236, 238, 240, 241, 242, 244, 245, 246, 247, 248, 250,252,
               253, 504, 505,
               1224, 1225, 1226, 1227, 1228, 1229, 1230, 1231, 1235, 1236, 1238, 1240, 1241, 1242,1244, 1245, 1246,1247,
               1248, 1250, 1252, 1253, 1504, 1505]

    elif atlas["name"] == 'ONPRC18':
        ref = [57]
        print('This atlas is not adapted : there is no subdivision between gray and white matter')

    elif atlas["name"] in ['CHARM-SARM','SARM']:
        if atlas["level"] == 0: ref = [203]
        elif atlas["level"] == 1: ref = [240]
        elif atlas["level"] == 2: ref = [244, 250]
        elif atlas["level"] == 3: ref = [244, 251, 267, 275]
        elif atlas["level"] == 4: ref = [244, 251, 267, 275, 246]
        elif atlas["level"] == 5: ref = [252, 253, 269, 270, 271, 272, 273, 274, 276, 277]

    elif atlas["name"] == 'EDNIxCSCLR':
        ref=[7,1007] # so far only present in atlas["level"] == 0
    elif atlas["name"] == 'EDNIxCSC':
        ref = [7]

    return ref

def gettype(filenii, level, filecsv):

    result = {"type": None, "Refname": None, "rois": None}

    def gui(result):

        names_list = ['atlas','mask']

        # --- GUI ---
        root = tk.Tk()
        root.title("Select the reference type(ex: white")

        label = tk.Label(root, text="Choose one:")
        label.grid(row=0, column=0)
        # label.pack()

        # selection
        listbox = tk.Listbox(root)
        for name in names_list:
            listbox.insert(tk.END, name)
        # listbox.pack()
        listbox.grid(row=1, column=0)

        listbox.select_set(0)

        # --- get ---
        def validation():
            # Mapping name -> ROIs
            result["type"] = listbox.get(listbox.curselection())
            root.destroy()

        btn = tk.Button(root, text="Validate", command=validation)
        btn.grid(row=2, column=1)
        # btn.pack()

        root.mainloop()

        return result

    result= gui(result)

    if result["type"] == 'atlas':
        result = getval(filenii, result, level, filecsv)
    elif result["type"] == 'mask':
        result = getmask(result)


    return result


def getmask(result):

    names_list = ['white','csf','white+csf','gray','brain']

    # --- GUI ---
    root = tk.Tk()
    root.title("Select the mask for reference (ex: white")

    label = tk.Label(root, text="Choose one:")
    label.grid(row=0, column=0)
    # label.pack()

    # multi-selection Listbox
    listbox = tk.Listbox(root)
    for name in names_list:
        listbox.insert(tk.END, name)
    # listbox.pack()
    listbox.grid(row=1, column=0)
    listbox.select_set(0)

    # --- get ---
    def validation():
        # Mapping name -> ROIs
        result["Refname"] = listbox.get(listbox.curselection())
        result["rois"] = [1]

        root.destroy()

    btn = tk.Button(root, text="Validate", command=validation)
    btn.grid(row=2, column=1)
    # btn.pack()

    root.mainloop()

    return result