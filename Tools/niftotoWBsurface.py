import os
import re
import glob
import shutil
import struct
import numpy as np
import nibabel as nib
import subprocess

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output


# ─────────────────────────────────────────────────────────────────────────────
def _patch_cifti_palette_xml(nii_path, replacements: dict):
    """
    Patch PaletteColorMapping fields directly in the NIfTI2/CIFTI binary XML
    extension (ecode=32).  This is the ONLY reliable way to set fields that
    wb_command -cifti-palette does not expose as CLI flags, such as:
        ShowTickMarksSelected, NumericSubivisions, PrecisionDigits,
        DisplayZeroData, ColorBarValuesMode, NumericFormatMode.

    The XML stores the PaletteColorMapping as HTML-escaped text inside a
    <Value> node, so tags appear as &lt;TAG&gt;value&lt;/TAG&gt;.
    Both escaped and plain forms are handled.

    Parameters
    ----------
    nii_path     : str  - Path to .dscalar.nii (modified in-place)
    replacements : dict - {tag_name: new_value}, e.g.
                          {'ShowTickMarksSelected': 'true', 'PaletteName': 'JET256'}
    """
    with open(nii_path, 'rb') as f:
        raw = f.read()

    hdr_end    = 540          # NIfTI2 header size
    ext_start  = hdr_end + 4  # skip 4-byte extension flag
    esize      = struct.unpack_from('<i', raw, ext_start)[0]
    ecode      = struct.unpack_from('<i', raw, ext_start + 4)[0]
    assert ecode == 32, f"Expected CIFTI ecode=32, got {ecode}"

    content_start = ext_start + 8
    content_end   = ext_start + esize
    xml_str       = raw[content_start:content_end].decode('utf-8', errors='replace')

    for tag, value in replacements.items():
        for open_tag, close_tag in [
            (f'&lt;{tag}&gt;', f'&lt;/{tag}&gt;'),
            (f'<{tag}>',       f'</{tag}>'),
        ]:
            pattern     = re.escape(open_tag) + r'[^<&]*' + re.escape(close_tag)
            replacement = f'{open_tag}{value}{close_tag}'
            xml_str, n  = re.subn(pattern, replacement, xml_str)
            if n:
                break   # found in one form, stop

    new_bytes = xml_str.encode('utf-8')
    max_len   = esize - 8
    assert len(new_bytes) <= max_len, \
        f"Patched XML too long: {len(new_bytes)} > {max_len}"
    padded = new_bytes + b'\x00' * (max_len - len(new_bytes))

    with open(nii_path, 'wb') as f:
        f.write(raw[:content_start] + padded + raw[content_end:])


# ─────────────────────────────────────────────────────────────────────────────
def render_stat_maps(results_dir, results_subdir, stat_suffixes,
                     regions_of_interest, surface_cfg):
    """
    For each region of interest and each stat suffix, find the output nifti
    and render it to surface.

    Parameters
    ----------
    results_dir         : base Results/ directory
    results_subdir      : subdirectory name, e.g. 'Grp_SBA_3dTTEST'
    stat_suffixes       : list of filename substrings to glob for
    regions_of_interest : list of region folder names to process
    surface_cfg         : dict of kwargs forwarded to nifti_to_surface_wb
    """
    for region in regions_of_interest:
        region_dir = opj(results_dir, results_subdir, region)
        if not os.path.isdir(region_dir):
            print(f"  [WARN] Result dir not found, skipping: {region_dir}")
            continue
        for suffix in stat_suffixes:
            matches = glob.glob(opj(region_dir, f'*{suffix}*.nii.gz'))
            if not matches:
                print(f"  [WARN] No nifti found for '{suffix}' in {region_dir}")
                continue
            for nifti_path in matches:
                print(f"\n  [surface] Rendering: {opb(nifti_path)}")
                nifti_to_surface_wb(nifti_path=nifti_path,
                                    output_dir=region_dir,
                                    **surface_cfg)


# ─────────────────────────────────────────────────────────────────────────────
def nifti_to_surface_wb(
    nifti_path: str,
    output_dir: str,
    species_name: str,
    atlas_lib_path: str,
    species_fragment: str,
    sing_wb: str,
    surface_type: str    = 'midthickness',
    projection: str      = 'enclosing',
    palette_name: str    = 'JET256',
    pos_max: float       = None,
    neg_max: float       = None,
    colorbar_mode: str   = 'fixed',
    scene_template: str  = None,
    scene_label: str     = 'Map0',
    png_width: int       = 500,
    png_height: int      = 100,
    map_prefix: str      = None,
    overwrite: bool      = True,
    delete_surf_file: bool = True,
    # Colorbar display options
    precision_digits: int      = 2,
    numeric_subdivisions: int  = 0,
    show_tick_marks: bool      = False,
    colorbar_values_mode: str  = 'DATA',
    numeric_format_mode: str   = 'AUTO',
    display_zero: bool         = True,
    display_positive: bool     = True,
    display_negative: bool     = True,
):
    """
    Project any 3D or 4D NIfTI volume onto a species surface using wb_command,
    produce .dscalar.nii files (one per volume/component), and render PNGs
    via wb_command -show-scene if a scene template is provided.

    Parameters
    ----------
    nifti_path           : Path to the input .nii.gz (3D or 4D)
    output_dir           : Where outputs go (.func.gii, .dscalar.nii, .png)
    species_name         : Species label, e.g. 'Mouselemur'
    atlas_lib_path       : Root of the atlas library
    species_fragment     : Relative path to species in atlas lib,
                           e.g. 'atlas/mammals/primates/Mouselemur'
    sing_wb              : Singularity + wb_command prefix string
    surface_type         : Surface mesh (only for projection='enclosing'):
                           'midthickness'|'inflated'|'very_inflated'|
                           'pial'|'white'|'flat'
    projection           : 'enclosing'|'ribbon'|'ribbon_dil'|
                           'vox'|'vox_dil'|'vox_pial'|'vox_inflated'
    palette_name         : wb_command palette name (default: 'JET256')
    pos_max              : Upper bound (meaning depends on colorbar_mode)
    neg_max              : Lower bound (meaning depends on colorbar_mode)
    colorbar_mode        : 'fixed'|'percent'|'abs_pct'|'full'
    scene_template       : Path to .scene template.
                           Must use '{species}.thickness.dscalar.nii' as placeholder.
    scene_label          : Scene name in .scene file (default: 'Map0')
    png_width/height     : PNG dimensions in pixels
    map_prefix           : Output filename prefix (default: NIfTI stem)
    overwrite            : Overwrite existing outputs (default: True)
    delete_surf_file     : Delete Native_resol copy after rendering
    precision_digits     : Decimal digits on colorbar labels (default: 2)
    numeric_subdivisions : Number of colorbar tick subdivisions (default: 0)
    show_tick_marks      : Show tick marks on colorbar (default: False)
    colorbar_values_mode : 'DATA' or 'PERCENTAGE' (default: 'DATA')
    numeric_format_mode  : 'AUTO'|'Decimal'|'Scientific' (default: 'AUTO')
    display_zero         : Show zero values on brain (default: True)
    display_positive     : Show positive values on brain (default: True)
    display_negative     : Show negative values on brain (default: True)

    Returns
    -------
    dscalar_paths : list of str
    png_paths     : list of str
    """

    # ── Internal helper ───────────────────────────────────────────────────────
    def _run(label, cmd):
        print(f"  [{label}] {cmd}", flush=True)
        out = spco(cmd, shell=True, stderr=subprocess.STDOUT)
        if out:
            print(f"  [{label}] {out.decode(errors='replace')[:300]}", flush=True)

    # ── Validate inputs ───────────────────────────────────────────────────────
    valid_modes = ('fixed', 'percent', 'abs_pct', 'full')
    if colorbar_mode not in valid_modes:
        raise ValueError(f"colorbar_mode must be one of {valid_modes}, got '{colorbar_mode}'")

    valid_proj = ('enclosing', 'ribbon', 'ribbon_dil',
                  'vox', 'vox_dil', 'vox_pial', 'vox_inflated')
    if projection not in valid_proj:
        raise ValueError(f"projection must be one of {valid_proj}, got '{projection}'")

    # ── Resolve surface paths ─────────────────────────────────────────────────
    surf_dir = opj(atlas_lib_path, species_fragment, 'EDNiX', 'surfaces', 'Native_resol')
    if not ope(surf_dir):
        raise FileNotFoundError(f"Surface directory not found: {surf_dir}")

    def _surf(h, stype):
        p = opj(surf_dir, f'{species_name}.{h}.{stype}.surf.gii')
        if not ope(p):
            raise FileNotFoundError(f"Surface file not found: {p}")
        return p

    surf_mid   = {h: _surf(h, 'midthickness') for h in ('l', 'r')}
    surf_pial  = {h: _surf(h, 'pial')         for h in ('l', 'r')}
    surf_white = {h: _surf(h, 'white')        for h in ('l', 'r')}
    surf_inf   = {h: _surf(h, 'inflated')     for h in ('l', 'r')}

    surf_main = {
        'enclosing'   : {h: _surf(h, surface_type) for h in ('l', 'r')},
        'ribbon'      : surf_mid,
        'ribbon_dil'  : surf_mid,
        'vox'         : surf_mid,
        'vox_dil'     : surf_mid,
        'vox_pial'    : surf_pial,
        'vox_inflated': surf_inf,
    }[projection]

    # ── Load NIfTI ────────────────────────────────────────────────────────────
    img    = nib.load(nifti_path)
    shape  = img.shape
    n_vols = shape[3] if len(shape) == 4 else 1
    print(f"  NIfTI shape: {shape}  ->  {n_vols} vol(s) | projection: {projection}",
          flush=True)

    os.makedirs(output_dir, exist_ok=True)

    stem   = opb(nifti_path).replace('.nii.gz', '').replace('.nii', '')
    prefix = map_prefix if map_prefix else stem
    hemis  = ['l', 'r']

    # ── Compute palette scale values ──────────────────────────────────────────
    data     = img.get_fdata()
    data_max = float(np.nanmax(data))
    data_min = float(np.nanmin(data))

    if colorbar_mode == 'fixed':
        _pos = pos_max if pos_max is not None else data_max
        _neg = neg_max if neg_max is not None else data_min
        print(f"  Palette [fixed]:   neg={_neg:.4f}  ->  pos={_pos:.4f}", flush=True)
        scale_mode    = 'MODE_USER_SCALE'
        scale_flags   = f'-pos-user 0 {_pos:.6f} -neg-user {_neg:.6f} 0'
        user_values   = f'{_neg:.6f} 0.000000 0.000000 {_pos:.6f}'
        pct_values    = '98.000000 2.000000 4.000000 96.000000'
        abs_values    = '20.000000 80.000000'

    elif colorbar_mode == 'percent':
        pct_pos = pos_max if pos_max is not None else 98.0
        pct_neg = neg_max if neg_max is not None else 2.0
        print(f"  Palette [percent]: lower={pct_neg:.2f}%  upper={pct_pos:.2f}%", flush=True)
        scale_mode    = 'MODE_AUTO_SCALE_PERCENTAGE'
        scale_flags   = f'-pos-percent 0 {pct_pos:.2f} -neg-percent {pct_neg:.2f} 0'
        user_values   = f'{data_min:.6f} 0.000000 0.000000 {data_max:.6f}'
        pct_values    = f'{pct_pos:.6f} {pct_neg:.6f} {pct_neg:.6f} {pct_pos:.6f}'
        abs_values    = '20.000000 80.000000'

    elif colorbar_mode == 'abs_pct':
        pct = pos_max if pos_max is not None else 98.0
        print(f"  Palette [abs_pct]: clip at {pct:.2f}%", flush=True)
        scale_mode    = 'MODE_AUTO_SCALE_ABSOLUTE_PERCENTAGE'
        scale_flags   = f'-pos-percent 0 {pct:.2f} -neg-percent 0 {pct:.2f}'
        user_values   = f'{data_min:.6f} 0.000000 0.000000 {data_max:.6f}'
        pct_values    = '98.000000 2.000000 4.000000 96.000000'
        pct_lo        = 100.0 - pct
        abs_values    = f'{pct_lo:.6f} {pct:.6f}'

    elif colorbar_mode == 'full':
        print(f"  Palette [full]:    {data_min:.4f}  ->  {data_max:.4f}", flush=True)
        scale_mode    = 'MODE_AUTO_SCALE'
        scale_flags   = ''
        user_values   = f'{data_min:.6f} 0.000000 0.000000 {data_max:.6f}'
        pct_values    = '98.000000 2.000000 4.000000 96.000000'
        abs_values    = '20.000000 80.000000'

    disp_flags   = (
        f'-disp-pos {"true" if display_positive else "false"} '
        f'-disp-neg {"true" if display_negative else "false"} '
        f'-disp-zero {"true" if display_zero else "false"}'
    )
    palette_flags = f'-palette-name {palette_name} {scale_flags} {disp_flags}'.strip()

    # Fields that -cifti-palette does NOT expose as CLI flags —
    # must be patched directly in the binary CIFTI XML header.
    cifti_xml_patches = {
        'ScaleMode'                      : scale_mode,
        'AutoScalePercentageValues'      : pct_values,
        'AutoScaleAbsolutePercentageValues': abs_values,
        'UserScaleValues'                : user_values,
        'PaletteName'                    : palette_name,
        'DisplayPositiveData'            : 'true' if display_positive else 'false',
        'DisplayZeroData'                : 'true' if display_zero     else 'false',
        'DisplayNegativeData'            : 'true' if display_negative else 'false',
        'NumericFormatMode'              : numeric_format_mode,
        'PrecisionDigits'                : str(precision_digits),
        'NumericSubivisions'             : str(numeric_subdivisions),   # note WB typo
        'ColorBarValuesMode'             : colorbar_values_mode,
        'ShowTickMarksSelected'          : 'true' if show_tick_marks else 'false',
    }

    # tick_str used by Patch D (m_showTickMarksSelected direct nodes in scene)
    tick_str = 'true' if show_tick_marks else 'false'

    # ── Prepare Native_resol copy for scene rendering ─────────────────────────
    surf_copy_dir = opj(output_dir, 'Native_resol')
    if scene_template is not None:
        if overwrite and ope(surf_copy_dir):
            shutil.rmtree(surf_copy_dir)
            print(f"  [copy] Removed existing Native_resol copy.", flush=True)
        if not ope(surf_copy_dir):
            print(f"  [copy] Copying Native_resol -> {surf_copy_dir} ...", flush=True)
            shutil.copytree(surf_dir, surf_copy_dir)
            print(f"  [copy] Done.", flush=True)

    dscalar_paths = []
    png_paths     = []

    # ── Loop over volumes ─────────────────────────────────────────────────────
    for vol_idx in range(n_vols):

        map_name     = f'{prefix}_vol{vol_idx}' if n_vols > 1 else prefix
        dscalar_path = opj(output_dir, f'{map_name}.dscalar.nii')

        if not overwrite and ope(dscalar_path):
            print(f"  [skip] {opb(dscalar_path)} already exists.", flush=True)
            dscalar_paths.append(dscalar_path)
            continue

        print(f"\n  === Volume {vol_idx+1}/{n_vols} : {map_name} ===", flush=True)

        if n_vols > 1:
            vol_nii  = nib.Nifti1Image(img.dataobj[..., vol_idx], img.affine, img.header)
            vol_path = opj(output_dir, f'_tmp_{map_name}.nii.gz')
            nib.save(vol_nii, vol_path)
        else:
            vol_path = nifti_path

        func_giis = {}

        # ── Step 1 : volume-to-surface mapping ────────────────────────────────
        for h in hemis:
            out_gii      = opj(output_dir, f'{map_name}.{h}.func.gii')
            func_giis[h] = out_gii
            if overwrite and ope(out_gii):
                os.remove(out_gii)

            if projection == 'enclosing':
                cmd = (f'{sing_wb} wb_command -volume-to-surface-mapping '
                       f'{vol_path} {surf_main[h]} {out_gii} -enclosing')
            elif projection == 'ribbon':
                cmd = (f'{sing_wb} wb_command -volume-to-surface-mapping '
                       f'{vol_path} {surf_pial[h]} {out_gii} '
                       f'-ribbon-constrained {surf_white[h]} {surf_pial[h]}')
            elif projection == 'ribbon_dil':
                cmd = (f'{sing_wb} wb_command -volume-to-surface-mapping '
                       f'{vol_path} {surf_pial[h]} {out_gii} '
                       f'-ribbon-constrained {surf_white[h]} {surf_pial[h]} '
                       f'-dilate-missing 10')
            elif projection == 'vox':
                cmd = (f'{sing_wb} wb_command -volume-to-surface-mapping '
                       f'{vol_path} {surf_mid[h]} {out_gii} -trilinear')
            elif projection == 'vox_dil':
                cmd = (f'{sing_wb} wb_command -volume-to-surface-mapping '
                       f'{vol_path} {surf_mid[h]} {out_gii} '
                       f'-trilinear -dilate-missing 10')
            elif projection == 'vox_pial':
                cmd = (f'{sing_wb} wb_command -volume-to-surface-mapping '
                       f'{vol_path} {surf_pial[h]} {out_gii} -trilinear')
            elif projection == 'vox_inflated':
                cmd = (f'{sing_wb} wb_command -volume-to-surface-mapping '
                       f'{vol_path} {surf_inf[h]} {out_gii} -trilinear')
            _run(f'vol-to-surf {h}', cmd)

        # ── Step 2 : set map names ─────────────────────────────────────────────
        for h in hemis:
            _run(f'set-map-name {h}',
                 f'{sing_wb} wb_command -set-map-names {func_giis[h]} -map 1 {map_name}')

        # ── Step 3 : create dense scalar ──────────────────────────────────────
        if overwrite and ope(dscalar_path):
            os.remove(dscalar_path)
        _run('cifti-create',
             f'{sing_wb} wb_command -cifti-create-dense-scalar '
             f'{dscalar_path} '
             f'-left-metric  {func_giis["l"]} '
             f'-right-metric {func_giis["r"]}')

        # ── Step 4 : set map name on dscalar ──────────────────────────────────
        _run('set-map-name dscalar',
             f'{sing_wb} wb_command -set-map-names {dscalar_path} -map 1 {map_name}')

        # ── Step 5 : apply palette via wb_command (scale/name only) ───────────
        # -cifti-palette handles ScaleMode + PaletteName but silently resets
        # ShowTickMarksSelected, NumericSubivisions, DisplayZeroData etc. to defaults.
        _run('cifti-palette',
             f'{sing_wb} wb_command -cifti-palette '
             f'{dscalar_path} {scale_mode} {dscalar_path} {palette_flags}')

        # ── Step 5b : patch CIFTI binary XML directly ─────────────────────────
        # Fix fields that -cifti-palette cannot set via CLI.
        # This writes directly into the NIfTI2 extension XML blob.
        print(f"  [cifti-xml-patch] Applying {len(cifti_xml_patches)} field(s)...",
              flush=True)
        _patch_cifti_palette_xml(dscalar_path, cifti_xml_patches)

        dscalar_paths.append(dscalar_path)

        # ── Step 6 : render PNG via scene ─────────────────────────────────────
        if scene_template is not None:
            if not ope(scene_template):
                print(f"  [WARN] Scene template not found: {scene_template}", flush=True)
            else:
                with open(scene_template) as f:
                    scene_txt = f.read()

                # Patch A: dscalar filename
                old_dscalar_base = f'{species_name}.thickness.dscalar.nii'
                new_dscalar_base = opb(dscalar_path)
                scene_txt = scene_txt.replace(old_dscalar_base, new_dscalar_base)
                scene_txt = scene_txt.replace('thickness.dscalar.nii', new_dscalar_base)

                # Patch B: fix hardcoded /volumes path
                atlas_ednix = opj(atlas_lib_path, species_fragment, 'EDNiX')
                m = re.search(r'(/[^\s<>"]+/volumes)(?=[</"])', scene_txt)
                if m:
                    old_vol = m.group(1)
                    new_vol = opj(atlas_ednix, 'volumes')
                    scene_txt = scene_txt.replace(old_vol, new_vol)
                    print(f"  [scene] vol path: {old_vol} -> {new_vol}", flush=True)

                # Patch C: remove ALL savedPaletteColorMappingArray blocks.
                # When the scene overrides palette via savedPaletteColorMapping,
                # WB ignores ShowTickMarksSelected / DisplayZeroData / PrecisionDigits
                # from the file and resets them to defaults.
                # By removing these blocks (setting Length=0), WB falls back to
                # reading the palette exclusively from the dscalar file itself —
                # which we already patched correctly in Step 5b.
                n_removed = [0]
                def _remove_saved_palette(match):
                    n_removed[0] += 1
                    indent = re.search(r'(\s*)<ObjectArray', match.group(0))
                    ind = indent.group(1) if indent else '                                        '
                    return f'{ind}<ObjectArray Type="class" Name="savedPaletteColorMappingArray" Length="0">\n{ind}</ObjectArray>'

                scene_txt = re.sub(
                    r'\s*<ObjectArray Type="class" Name="savedPaletteColorMappingArray".*?</ObjectArray>',
                    _remove_saved_palette, scene_txt, flags=re.DOTALL)
                print(f"  [scene] Removed {n_removed[0]} savedPaletteColorMappingArray block(s) "
                      f"(palette read from file instead)", flush=True)

                # Patch D: m_showTickMarksSelected direct XML nodes in
                # all m_colorBar Annotation objects
                scene_txt = re.sub(
                    r'(<Object Type="boolean" Name="m_showTickMarksSelected">)'
                    r'(?:true|false)(</Object>)',
                    rf'\g<1>{tick_str}\g<2>', scene_txt)

                # Write patched scene
                scene_out = opj(surf_copy_dir, f'{map_name}.scene')
                if overwrite and ope(scene_out):
                    os.remove(scene_out)
                with open(scene_out, 'w') as f:
                    f.write(scene_txt)

                # Copy dscalar (with both name variants) into surf_copy_dir
                for link_name in (new_dscalar_base,
                                  f'{species_name}.{new_dscalar_base}'):
                    dst = opj(surf_copy_dir, link_name)
                    if overwrite and ope(dst):
                        os.remove(dst)
                    if not ope(dst):
                        shutil.copy2(dscalar_path, dst)

                # Render PNG
                png_out          = opj(output_dir, f'{map_name}.png')
                sing_wb_headless = sing_wb.replace('vglrun ', '')
                _run('show-scene',
                     f'cd {surf_copy_dir} && '
                     f'{sing_wb_headless} wb_command -show-scene '
                     f'{opb(scene_out)} "{scene_label}" '
                     f'../{opb(png_out)} {png_width} {png_height} -use-window-size')
                png_paths.append(png_out)

        # ── Cleanup temp volume ───────────────────────────────────────────────
        if n_vols > 1 and ope(vol_path):
            os.remove(vol_path)

    # ── Cleanup Native_resol copy ─────────────────────────────────────────────
    if scene_template is not None and ope(surf_copy_dir) and delete_surf_file:
        shutil.rmtree(surf_copy_dir)
        print(f"\n  [cleanup] Removed Native_resol copy: {surf_copy_dir}", flush=True)

    print(f"\n  Done. {len(dscalar_paths)} .dscalar.nii | {len(png_paths)} .png produced.",
          flush=True)
    return dscalar_paths, png_paths