"""
Atlas label utilities used by the pipeline (moved from
Exemples/Study_EDNiX/analysis/ednix_bids_tools.py so that the core package
does not depend on study scripts).
"""
import os
import re
import pandas as pd

opj = os.path.join


def _linux_path(p):
    return str(p).replace("\\", "/")

def parse_label_file(label_path):
    """
    Parse an EDNiX atlas label / StatsLUT file.

    Supported formats:
      StatsLUT:  ID  NAME  R  G  B  A      (one entry per line)
      label:     NAME\\nID R G B A\\n...   (alternating name / colour lines)

    Returns DataFrame with columns:
        label_id, region_name, base_name, hemisphere, R, G, B, A
    where hemisphere is 'left' | 'right' | 'bilateral' and
    base_name has the L_/R_ prefix stripped.
    """
    with open(label_path, encoding="utf-8-sig") as f:
        lines = [l.strip() for l in f if l.strip()]

    def _hemi(name):
        if   re.match(r"^[Ll]_", name): return "left",     re.sub(r"^[Ll]_", "", name)
        elif re.match(r"^[Rr]_", name): return "right",    re.sub(r"^[Rr]_", "", name)
        else:                            return "bilateral", name

    # Detect format: StatsLUT if first non-empty token of first line is a digit
    rows = []
    first_parts = lines[0].split()
    is_statslut = len(first_parts) >= 2 and first_parts[0].isdigit()

    if is_statslut:
        for line in lines:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                lid = int(parts[0])
            except ValueError:
                continue
            name = parts[1]
            if name.lower() in ("unknown", "???"):
                continue
            try:
                r, g, b, a = int(parts[2]), int(parts[3]), int(parts[4]), int(parts[5])
            except (IndexError, ValueError):
                r = g = b = a = 0
            hemi, base = _hemi(name)
            rows.append(dict(label_id=lid, region_name=name, base_name=base,
                             hemisphere=hemi, R=r, G=g, B=b, A=a))
    else:
        # Alternating name / "ID R G B A" lines
        i = 0
        while i < len(lines) - 1:
            name_line  = lines[i]
            color_line = lines[i + 1]
            parts = color_line.split()
            if len(parts) >= 5 and parts[0].isdigit():
                lid = int(parts[0])
                if lid != 0 and name_line not in ("???", "Unknown"):
                    r, g, b, a = int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])
                    hemi, base = _hemi(name_line)
                    rows.append(dict(label_id=lid, region_name=name_line,
                                     base_name=base, hemisphere=hemi,
                                     R=r, G=g, B=b, A=a))
                i += 2
            else:
                i += 1

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df[df["label_id"] != 0].reset_index(drop=True)
    return df

def get_atlas_label_path(atlas_library_root, species_path_fragment, atlas_name,
                          prefer_statslut=False):
    label_dir = opj(atlas_library_root, species_path_fragment, "label_code")
    priority = (
        [f"{atlas_name}_StatsLUT.txt", f"{atlas_name}_label.txt"]
        if prefer_statslut
        else [f"{atlas_name}_label.txt", f"{atlas_name}_StatsLUT.txt"]
    )
    for fn in priority:
        c = opj(label_dir, fn)
        if os.path.exists(c):
            return _linux_path(c)
    raise FileNotFoundError(f"No label file for atlas '{atlas_name}' in {label_dir}")

def find_species_path(atlas_library_root, species_name):
    atlas_root = opj(atlas_library_root, "atlas")
    matches = []
    for dirpath, dirnames, _ in os.walk(atlas_root):
        if os.path.basename(dirpath).lower() == species_name.lower():
            matches.append(_linux_path(os.path.relpath(dirpath, atlas_library_root)))
            dirnames.clear()
    if not matches:
        raise FileNotFoundError(f"No species folder '{species_name}' under {atlas_root}")
    if len(matches) > 1:
        raise ValueError(f"Ambiguous species '{species_name}':\n" +
                         "\n".join(f"  {m}" for m in matches))
    return matches[0]
