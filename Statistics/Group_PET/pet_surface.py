"""
Surface figures of the PET group maps, with the function of the seed-based fMRI figures
(Tools/niftotoWBsurface.nifti_to_surface_wb: wb_command volume-to-surface mapping, dense scalar, scene).

The surfaces are those of the reference template of the atlas library (<atlas_lib>/<species_fragment>/
EDNiX/surfaces/Native_resol): use them with maps computed in that template space (pet_group with
template = the `reference` of the PET launcher), not in a study template.

Projection: 'ribbon' by default (-ribbon-constrained between the white and pial surfaces). A PET voxel is
large compared with the cortical thickness; trilinear sampling on the midthickness surface ('vox') mixes the
CSF outside the pial surface and the white matter into the cortical value, the ribbon projection averages
the voxels of the cortical ribbon only (Glasser et al. 2013, HCP minimal preprocessing). The analysis mask of
pet_group (relative threshold, smoothing inside the brain) already limits the CSF signal: vertices whose
ribbon is outside the mask receive 0 (not displayed with display_zero=False).
"""
import os

import numpy as np
import nibabel as nib

opj = os.path.join
opb = os.path.basename
ope = os.path.exists

# map key of pet_group.run_group_analysis -> (palette, symmetric scale)
DEFAULT_MAPS = {
    'zstat_thresholded': ('JET256', True),
    'zstat': ('JET256', True),
    'percent_difference': ('JET256', True),
    'mean_all': ('JET256', False),
}
# other maps (rendered when asked with maps=[...]): group mean images of pet_group.run_group_mean
PALETTES = dict(DEFAULT_MAPS, sd_all=('videen_style', False), cv_all=('videen_style', False))


def render_pet_surfaces(result, species, atlas_lib, species_fragment, sing_wb, scene_template=None,
                        scene_label='New Scene 1', maps=None, projection='ribbon', surface_type='inflated',
                        output_dir=None, png_width=1200, png_height=400, delete_surf_file=True, ranges=None):
    """
    Projects maps of a group result onto the surfaces of the species and renders them.

    result          : dict returned by pet_group.run_group_analysis (or {'maps': {key: path}, 'output_dir': ...})
    species, atlas_lib, species_fragment, sing_wb, scene_template, scene_label : as the SBA figures
                      (Exemples/Study_EDNiX/analysis/Figure_5_SBA.py); the scene template must use
                      '<species>.thickness.dscalar.nii' as placeholder (<species>_thickness.scene of the atlas
                      library). Without scene, only the .dscalar.nii files are written.
    maps            : keys of result['maps'] to render (default: those of DEFAULT_MAPS present; a result of
                      pet_group.run_group_mean: 'mean_all', plus 'sd_all' / 'cv_all' if asked)
    projection      : 'ribbon' (default), 'ribbon_dil', 'enclosing', 'vox'... (niftotoWBsurface)
    ranges          : {key: (neg_max, pos_max)} fixed colour scale (e.g. the same bratio range for every species of
                      a cross-species figure); default: 1st / 99th percentile of each map (symmetric: 99th of |x|)
    Returns {key: dscalar path} (the PNG has the same name, .png, when a scene is given).
    """
    from Tools import niftotoWBsurface
    out_dir = output_dir or opj(result['output_dir'], 'surfaces')
    os.makedirs(out_dir, exist_ok=True)
    keys = maps or [k for k in DEFAULT_MAPS if k in result['maps']]
    done = {}
    for k in keys:
        path = result['maps'].get(k)
        if not path or not ope(path):
            print(f'  [surface] {k}: no map, skipped', flush=True)
            continue
        data = np.asarray(nib.load(path).dataobj, float)
        nz = data[np.isfinite(data) & (data != 0)]
        if nz.size == 0:
            print(f'  [surface] {opb(path)} is empty (no significant cluster): skipped', flush=True)
            continue
        palette, symmetric = PALETTES.get(k, ('ROY-BIG-BL', True))
        if ranges and k in ranges:
            neg_max, pos_max = (float(v) for v in ranges[k])
        elif symmetric:
            vmax = float(np.percentile(np.abs(nz), 99))
            pos_max, neg_max = vmax, -vmax
        else:
            pos_max, neg_max = float(np.percentile(nz, 99)), float(np.percentile(nz, 1))
        niftotoWBsurface.nifti_to_surface_wb(
            nifti_path=path, output_dir=out_dir, species_name=species, atlas_lib_path=atlas_lib,
            species_fragment=species_fragment, sing_wb=sing_wb, surface_type=surface_type, projection=projection,
            palette_name=palette, pos_max=pos_max, neg_max=neg_max, colorbar_mode='fixed',
            scene_template=scene_template, scene_label=scene_label, png_width=png_width, png_height=png_height,
            delete_surf_file=delete_surf_file, display_zero=False, display_negative=symmetric)
        stem = opb(path).replace('.nii.gz', '')
        done[k] = opj(out_dir, stem + '.dscalar.nii')
    return done
