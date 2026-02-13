"""
fill_all_ednix_labels.py

Scans /home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/ for:
  - *.label.gii  (single hemisphere, e.g. Human.r.EDNiXCSCLR.native.label.gii)
  - *.dlabel.nii (both hemispheres CIFTI, e.g. Dog.EDNIxCSC.dlabel.nii)

Strategy: for each cluster of unlabeled (==0) vertices, check its boundary:
  - Single-label boundary → reconstruction artifact → fill with that label (NN)
  - Mixed-label boundary  → intentional empty zone  → keep as 0
  - No boundary           → isolated / unreachable  → keep as 0

Output: *filled.label.gii / *filled.dlabel.nii in the same directory.
"""

import os
import re
import glob
import nibabel as nib
import numpy as np
from scipy.spatial import cKDTree


# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────
ATLAS_ROOT    = '/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/'
PATTERNS_GII  = ['*EDNiXCSC.native.label.gii',  '*EDNiXCSCLR.native.label.gii']
PATTERNS_NII  = ['*EDNIxCSC.dlabel.nii',         '*EDNIxCSCLR.dlabel.nii',
                 '*EDNiXCSC.dlabel.nii',          '*EDNiXCSCLR.dlabel.nii']
DRY_RUN       = False   # Set True to preview without saving


# ─────────────────────────────────────────────────────────────────────────────
# FILE DISCOVERY
# ─────────────────────────────────────────────────────────────────────────────
def find_files(root, patterns):
    found = []
    for pattern in patterns:
        found.extend(glob.glob(os.path.join(root, '**', pattern), recursive=True))
    return sorted(set(f for f in found if 'filled' not in os.path.basename(f)))


# Cache: surf_path -> n_vertices
_surf_vertex_cache = {}

def get_surf_n_vertices(surf_path):
    if surf_path not in _surf_vertex_cache:
        try:
            s = nib.load(surf_path)
            _surf_vertex_cache[surf_path] = len(s.darrays[0].data)
        except Exception:
            _surf_vertex_cache[surf_path] = -1
    return _surf_vertex_cache[surf_path]


def find_surface(label_path, hemi, n_vertices=None):
    """
    Find midthickness surface for a given hemisphere ('l' or 'r').
    If n_vertices is given, searches the whole species subtree for a surface
    with exactly that vertex count (for dlabel.nii resolution matching).
    Otherwise searches same dir / sibling dirs by name.
    """
    label_dir = os.path.dirname(label_path)

    # Same directory — always check first
    candidates = glob.glob(os.path.join(label_dir, f'*.{hemi}.*midthickness*.surf.gii'))
    if candidates:
        if n_vertices is None:
            return candidates[0]
        for c in candidates:
            if get_surf_n_vertices(c) == n_vertices:
                return c

    # Search species subtree (go up to find the species root)
    # Walk up until we find a directory whose parent is 'atlas'
    search_root = label_dir
    for _ in range(6):  # max 6 levels up
        parent = os.path.dirname(search_root)
        if os.path.basename(parent) == 'atlas' or parent == search_root:
            break
        search_root = parent

    candidates = glob.glob(
        os.path.join(search_root, '**', f'*.{hemi}.*midthickness*.surf.gii'),
        recursive=True
    )

    if n_vertices is not None:
        # Match by vertex count exactly
        for c in sorted(candidates):
            if get_surf_n_vertices(c) == n_vertices:
                return c
        return None  # No surface with matching vertex count found
    else:
        if candidates:
            res_folder = os.path.basename(label_dir)
            same_res = [c for c in candidates if res_folder in c]
            return same_res[0] if same_res else candidates[0]

    return None


# ─────────────────────────────────────────────────────────────────────────────
# MESH HELPERS
# ─────────────────────────────────────────────────────────────────────────────
def build_adjacency(triangles, n_vertices):
    adjacency = [set() for _ in range(n_vertices)]
    for tri in triangles:
        for i in range(3):
            for j in range(3):
                if i != j:
                    adjacency[tri[i]].add(tri[j])
    return adjacency


def find_clusters(labels, adjacency):
    visited = set()
    clusters = []
    for v in np.where(labels == 0)[0]:
        if v in visited:
            continue
        cluster = []
        stack = [v]
        while stack:
            node = stack.pop()
            if node in visited or labels[node] != 0:
                continue
            visited.add(node)
            cluster.append(node)
            stack.extend(adjacency[node])
        clusters.append(cluster)
    return clusters


def get_boundary_labels(cluster, adjacency, labels):
    boundary = set()
    for v in cluster:
        for nb in adjacency[v]:
            if labels[nb] != 0:
                boundary.add(int(labels[nb]))
    return boundary


# ─────────────────────────────────────────────────────────────────────────────
# CORE FILL ALGORITHM
# ─────────────────────────────────────────────────────────────────────────────
def fill_hemi(original_labels, adjacency, coords, level_idx, hemi_name):
    """
    Apply boundary-homogeneity fill to one hemisphere / one level.
    Returns filled labels array and stats.
    """
    labels      = original_labels.copy()
    n_unlabeled = (labels == 0).sum()

    print(f"    [{hemi_name}] Unlabeled: {n_unlabeled}")

    if n_unlabeled == 0:
        return labels, 0, 0

    clusters = find_clusters(labels, adjacency)
    sizes    = sorted([len(c) for c in clusters], reverse=True)

    print(f"    [{hemi_name}] {len(clusters)} clusters | Top 10: {sizes[:10]}")
    print(f"    [{hemi_name}] "
          f"< 10: {sum(1 for s in sizes if s < 10):4d} | "
          f"< 50: {sum(1 for s in sizes if s < 50):4d} | "
          f"<100: {sum(1 for s in sizes if s < 100):4d} | "
          f">=100: {sum(1 for s in sizes if s >= 100):4d} | "
          f">=500: {sum(1 for s in sizes if s >= 500):4d}")

    to_fill        = []
    to_keep        = []
    single_count   = 0
    mixed_count    = 0
    no_bound_count = 0

    for c in clusters:
        bl = get_boundary_labels(c, adjacency, original_labels)
        if len(bl) == 0:
            no_bound_count += 1
            to_keep.extend(c)
        elif len(bl) == 1:
            single_count += 1
            to_fill.extend(c)
        else:
            mixed_count += 1
            to_keep.extend(c)

    print(f"    [{hemi_name}] Single-label (fill): {single_count:4d} clusters ({len(to_fill)} vertices)")
    print(f"    [{hemi_name}] Mixed-label  (keep): {mixed_count:4d} clusters")
    print(f"    [{hemi_name}] No boundary  (keep): {no_bound_count:4d} clusters")

    if len(to_fill) > 0:
        labeled_mask    = original_labels != 0
        tree            = cKDTree(coords[labeled_mask])
        labeled_indices = np.where(labeled_mask)[0]
        to_fill_arr     = np.array(to_fill)
        _, indices      = tree.query(coords[to_fill_arr], k=1)
        labels[to_fill_arr] = original_labels[labeled_indices[indices]]

    assert np.all(labels[original_labels != 0] == original_labels[original_labels != 0]), \
        f"ERROR [{hemi_name}]: original labels modified!"

    print(f"    [{hemi_name}] Remaining unlabeled: {(labels == 0).sum()}")
    return labels, len(to_fill), len(to_keep)


# ─────────────────────────────────────────────────────────────────────────────
# PROCESS .label.gii (single hemisphere)
# ─────────────────────────────────────────────────────────────────────────────
def process_label_gii(label_path):
    print(f"\n{'='*70}")
    print(f"[GII] {os.path.basename(label_path)}")
    print(f"      {label_path}")

    # Detect hemisphere
    m = re.match(r'^.+\.(l|r)\.', os.path.basename(label_path))
    if not m:
        print("  ⚠ Cannot detect hemisphere from filename, skipping.")
        return False
    hemi = m.group(1)

    surf_path = find_surface(label_path, hemi)
    if surf_path is None:
        print(f"  ⚠ No midthickness surface found for hemi={hemi}, skipping.")
        return False
    print(f"  Surface: {os.path.basename(surf_path)}")

    label_img = nib.load(label_path)
    surf      = nib.load(surf_path)
    coords    = surf.darrays[0].data
    triangles = surf.darrays[1].data
    n_verts   = len(label_img.darrays[0].data)

    if len(coords) != n_verts:
        print(f"  ⚠ Vertex mismatch (label={n_verts}, surf={len(coords)}), skipping.")
        return False

    print(f"  Building adjacency ({n_verts} vertices)...")
    adjacency    = build_adjacency(triangles, n_verts)
    total_filled = 0

    for i, darray in enumerate(label_img.darrays):
        print(f"\n  --- Level {i} ---")
        filled, n_filled, _ = fill_hemi(
            darray.data.copy(), adjacency, coords, i, hemi.upper()
        )
        label_img.darrays[i].data = filled
        total_filled += n_filled

    out_name = os.path.basename(label_path).replace('.label.gii', '.filled.label.gii')
    out_path = os.path.join(os.path.dirname(label_path), out_name)

    if not DRY_RUN:
        nib.save(label_img, out_path)
        print(f"\n  ✓ Saved: {out_name}  ({total_filled} vertices filled)")
    else:
        print(f"\n  [DRY RUN] Would save: {out_path}")

    return True


# ─────────────────────────────────────────────────────────────────────────────
# PROCESS .dlabel.nii (both hemispheres CIFTI)
# ─────────────────────────────────────────────────────────────────────────────
def process_dlabel_nii(label_path):
    print(f"\n{'='*70}")
    print(f"[NII] {os.path.basename(label_path)}")
    print(f"      {label_path}")

    cifti     = nib.load(label_path)
    data      = cifti.get_fdata(dtype=np.float32)   # shape (n_levels, n_vertices)
    brain_ax  = cifti.header.get_axis(1)             # BrainModelAxis

    # Parse hemispheres from CIFTI brain model axis
    hemis = {}
    for name, slc, bm in brain_ax.iter_structures():
        if 'LEFT' in name.upper():
            hemi_key = 'l'
        elif 'RIGHT' in name.upper():
            hemi_key = 'r'
        else:
            continue
        # bm.nvertices = size of the full native surface
        # bm.vertex    = indices of CIFTI vertices in that surface (subset)
        n_native = list(bm.nvertices.values())[0]
        hemis[hemi_key] = {
            'slice':    slc,
            'n_cifti':  len(bm.vertex),   # vertices in CIFTI
            'n_native': n_native,          # vertices in native surface
            'vertex':   bm.vertex.copy(),  # mapping CIFTI -> native indices
        }

    if not hemis:
        print("  ⚠ No cortical structures found in CIFTI, skipping.")
        return False

    print(f"  Hemispheres: " +
          ", ".join(
              f"{h.upper()} ({v['n_cifti']} cifti / {v['n_native']} native)"
              for h, v in hemis.items()
          ))

    # Load surfaces for each hemisphere — match by native vertex count
    for hemi, info in hemis.items():
        surf_path = find_surface(label_path, hemi, n_vertices=info['n_native'])
        if surf_path is None:
            print(f"  ⚠ No midthickness surface with {info['n_native']} vertices "
                  f"for hemi={hemi}, skipping whole file.")
            return False
        info['surf_path'] = surf_path
        print(f"  Surface {hemi.upper()}: {os.path.basename(surf_path)} "
              f"({info['n_native']} vertices)")

        surf      = nib.load(surf_path)
        coords    = surf.darrays[0].data
        triangles = surf.darrays[1].data

        print(f"  Building adjacency {hemi.upper()} ({info['n_native']} vertices)...")
        info['coords']    = coords
        info['adjacency'] = build_adjacency(triangles, info['n_native'])

    # Process each level x each hemisphere
    n_levels     = data.shape[0]
    filled_data  = data.copy()
    total_filled = 0

    for i in range(n_levels):
        print(f"\n  --- Level {i} ---")
        for hemi, info in hemis.items():
            slc            = info['slice']
            vertex_indices = info['vertex']    # CIFTI -> native surface indices
            n_native       = info['n_native']
            adjacency      = info['adjacency']
            coords         = info['coords']

            # Get CIFTI labels for this level/hemi
            cifti_labels = data[i, slc].astype(np.int32)

            # Map CIFTI labels onto full native surface (unlabeled = 0)
            native_labels = np.zeros(n_native, dtype=np.int32)
            native_labels[vertex_indices] = cifti_labels

            # Run fill on the full native surface
            filled_native, n_filled, _ = fill_hemi(
                native_labels, adjacency, coords, i, hemi.upper()
            )

            # Map filled labels back to CIFTI
            filled_cifti = filled_native[vertex_indices]
            filled_data[i, slc] = filled_cifti
            total_filled += n_filled

    # Rebuild CIFTI with filled data
    new_img = nib.Cifti2Image(filled_data, header=cifti.header, nifti_header=cifti.nifti_header)

    out_name = os.path.basename(label_path).replace('.dlabel.nii', '.filled.dlabel.nii')
    out_path = os.path.join(os.path.dirname(label_path), out_name)

    if not DRY_RUN:
        nib.save(new_img, out_path)
        print(f"\n  ✓ Saved: {out_name}  ({total_filled} vertices filled)")
    else:
        print(f"\n  [DRY RUN] Would save: {out_path}")

    return True


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
def main():
    print(f"Atlas root : {ATLAS_ROOT}")
    if DRY_RUN:
        print("*** DRY RUN MODE — no files will be saved ***")

    gii_files = find_files(ATLAS_ROOT, PATTERNS_GII)
    nii_files = find_files(ATLAS_ROOT, PATTERNS_NII)

    print(f"\nFound {len(gii_files)} .label.gii file(s)")
    print(f"Found {len(nii_files)} .dlabel.nii file(s)")

    success, failed = 0, 0

    for path in gii_files:
        try:
            ok = process_label_gii(path)
            success += ok
            failed  += not ok
        except Exception as e:
            print(f"  ✗ ERROR: {e}")
            import traceback; traceback.print_exc()
            failed += 1

    for path in nii_files:
        try:
            ok = process_dlabel_nii(path)
            success += ok
            failed  += not ok
        except Exception as e:
            print(f"  ✗ ERROR: {e}")
            import traceback; traceback.print_exc()
            failed += 1

    print(f"\n{'='*70}")
    print(f"Done! {success} processed, {failed} failed/skipped.")


if __name__ == '__main__':
    main()