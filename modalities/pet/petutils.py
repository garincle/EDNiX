"""
Shared helpers of the PET and CT pipelines: file names, provenance sidecars, frame timing,
isotopes and decay, activity units, subject data for SUV, duration-weighted averages and
small image operations (nibabel / numpy, no external program).

Conventions
-----------
- Times are in seconds, relative to the BIDS `TimeZero` of the PET sidecar.
- `timing` is a dict with numpy arrays `start`, `duration`, `end`, `mid` (one value per frame).
- Activity concentrations are converted to Bq/mL when a quantitative value is computed (SUV).
"""
import os
import re
import json
import glob
import math
import datetime

import numpy as np
import nibabel as nib

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
opi = os.path.isfile
ope = os.path.exists

# ---------------------------------------------------------------------------------------------------------------------
# Isotopes
# ---------------------------------------------------------------------------------------------------------------------
# Physical half-lives (s), NNDC / Laboratoire National Henri Becquerel values.
HALF_LIFE_S = {'F-18': 6586.2,      # 109.77 min
               'C-11': 1221.8,      # 20.364 min
               'O-15': 122.24,
               'N-13': 597.9,       # 9.965 min
               'Ga-68': 4062.6,     # 67.71 min
               'Cu-64': 45720.0,    # 12.70 h
               'Zr-89': 282276.0,   # 78.41 h
               'Rb-82': 75.45,
               'I-124': 360806.0,   # 4.176 d
               'Br-76': 58320.0}    # 16.2 h

# Radionuclide of the tracers used with EDNiX (list of the PET engineer, completed).
KNOWN_TRACER = {'FDG': 'F-18', 'RACLO': 'C-11', 'RAC': 'C-11', 'BUPRE': 'C-11', '2FNQ1P': 'F-18',
                'H2O': 'O-15', 'PE2I': 'C-11', 'FDOPA': 'F-18', 'PIB': 'C-11', 'UCBJ': 'C-11',
                'FMZ': 'C-11', 'FLB457': 'C-11', 'DASB': 'C-11', 'MDL100907': 'C-11'}


def normalise_isotope(code):
    """'F18', '18F', 'F-18', '^18^Fluorine', 'Fluorine-18' -> 'F-18' (None if not recognised)."""
    if code is None:
        return None
    s = str(code).strip()
    if s in HALF_LIFE_S:
        return s
    names = {'fluorine': 'F', 'carbon': 'C', 'oxygen': 'O', 'nitrogen': 'N', 'gallium': 'Ga',
             'copper': 'Cu', 'zirconium': 'Zr', 'rubidium': 'Rb', 'iodine': 'I', 'bromine': 'Br'}
    low = s.lower()
    for name, sym in names.items():
        if name in low:
            m = re.search(r'(\d+)', s)
            if m:
                cand = f'{sym}-{m.group(1)}'
                return cand if cand in HALF_LIFE_S else None
    m = re.fullmatch(r'\^?(\d+)\^?([A-Za-z]{1,2})', s) or re.fullmatch(r'([A-Za-z]{1,2})-?(\d+)', s)
    if m:
        a, b = m.groups()
        sym, num = (b, a) if a.isdigit() else (a, b)
        cand = f'{sym.capitalize()}-{num}'
        return cand if cand in HALF_LIFE_S else None
    return None


def isotope_from_half_life(half_life_s, tol=0.05):
    """Isotope whose half-life is within `tol` (relative) of `half_life_s`, else None."""
    try:
        h = float(half_life_s)
    except (TypeError, ValueError):
        return None
    best = min(HALF_LIFE_S, key=lambda k: abs(HALF_LIFE_S[k] - h) / HALF_LIFE_S[k])
    return best if abs(HALF_LIFE_S[best] - h) / HALF_LIFE_S[best] <= tol else None


def decay_constant(isotope):
    if isotope not in HALF_LIFE_S:
        raise ValueError(f'unknown isotope {isotope!r}: expected one of {sorted(HALF_LIFE_S)}')
    return math.log(2) / HALF_LIFE_S[isotope]


# ---------------------------------------------------------------------------------------------------------------------
# Units
# ---------------------------------------------------------------------------------------------------------------------
_PREFIX = {'': 1.0, 'k': 1e3, 'M': 1e6, 'G': 1e9, 'm': 1e-3, 'u': 1e-6, 'µ': 1e-6, 'n': 1e-9, 'p': 1e-12}
CI_TO_BQ = 3.7e10


def _activity_factor(unit):
    """Factor converting an activity in `unit` (Bq, kBq, MBq, GBq, mCi, uCi, nCi...) to Bq."""
    u = unit.strip()
    m = re.fullmatch(r'([kMGmuµnp]?)(Bq|Ci)', u)
    if not m:
        raise ValueError(f'activity unit not recognised: {unit!r}')
    pre, base = m.groups()
    return _PREFIX[pre] * (CI_TO_BQ if base == 'Ci' else 1.0)


def concentration_to_bq_per_ml(unit):
    """Factor converting an image in `unit` (e.g. 'Bq/mL', 'kBq/cc', 'nCi/mL') to Bq/mL; None if unknown."""
    if not unit:
        return None
    u = str(unit).replace(' ', '').replace('cc', 'mL').replace('ml', 'mL')
    m = re.fullmatch(r'([kMGmuµnp]?(?:Bq|Ci))/([mµu]?L)', u)
    if not m:
        return None
    act, vol = m.groups()
    vol_ml = {'mL': 1.0, 'L': 1000.0, 'uL': 1e-3, 'µL': 1e-3}[vol]
    return _activity_factor(act) / vol_ml


def activity_to_mbq(value, unit):
    return float(value) * _activity_factor(unit or 'MBq') / 1e6


# ---------------------------------------------------------------------------------------------------------------------
# File names and provenance
# ---------------------------------------------------------------------------------------------------------------------
def strip_ext(path):
    b = opb(path)
    for ext in ('.nii.gz', '.nii', '.json', '.sif', '.dft', '.tsv', '.txt'):
        if b.endswith(ext):
            return b[:-len(ext)]
    return b


def base_name(path, suffix='pet'):
    """BIDS raw name without extension, suffix and derivative entities:
    'sub-01_ses-1_trc-FDG_run-1_desc-reorient_pet.nii.gz' -> 'sub-01_ses-1_trc-FDG_run-1'."""
    parts = strip_ext(path).split('_')
    if parts and parts[-1] == suffix:
        parts = parts[:-1]
    return '_'.join(p for p in parts if p.split('-')[0] not in ('desc', 'space', 'seg'))


def session_base(run_base):
    """Name shared by the runs of one session (run entity removed)."""
    return '_'.join(p for p in run_base.split('_') if not p.startswith('run-'))


def out_name(folder, base, desc=None, space=None, suffix='pet', ext='.nii.gz'):
    """<folder>/<base>[_space-<space>][_desc-<desc>]_<suffix><ext> (BIDS entity order)."""
    name = base
    if space:
        name += '_space-' + space
    if desc:
        name += '_desc-' + desc
    return opj(folder, name + '_' + suffix + ext)


def sidecar(path):
    for ext in ('.nii.gz', '.nii'):
        if path.endswith(ext):
            return path[:-len(ext)] + '.json'
    return path + '.json'


def write_json(path, sources, description, **extra):
    """Provenance sidecar written next to every output image (as in the rest of EDNiX)."""
    d = {"Sources": sources, "Description": description}
    d.update(extra)
    with open(sidecar(path), 'w') as f:
        json.dump(d, f, indent=2, default=str)


def read_json(path):
    if path and opi(path):
        with open(path) as f:
            return json.load(f)
    return {}


def raw_sidecar(raw_img):
    """BIDS sidecar of a raw image (same name, .json)."""
    return read_json(sidecar(raw_img))


def list_raw_images(path_mod, pattern, suffix='pet'):
    """Raw images of one session: files of <path_mod> matching `pattern` and ending with _<suffix>.nii[.gz]."""
    out = []
    for f in sorted(glob.glob(opj(path_mod, pattern))):
        if opi(f) and re.search(r'_' + suffix + r'\.nii(\.gz)?$', opb(f)):
            out.append(f)
    return out


# ---------------------------------------------------------------------------------------------------------------------
# Images
# ---------------------------------------------------------------------------------------------------------------------
def load(path):
    img = nib.load(path)
    return img, np.asarray(img.dataobj, dtype=np.float32)


def save_like(data, ref_img, path):
    """Save `data` with the geometry of `ref_img` (nibabel image), float32, no scaling."""
    hdr = ref_img.header.copy()
    hdr.set_data_dtype(np.float32)
    data = np.asarray(data, dtype=np.float32)
    hdr.set_data_shape(data.shape)
    if data.ndim == 3 and len(hdr.get_zooms()) > 3:
        hdr.set_zooms(hdr.get_zooms()[:3])
    out = nib.Nifti1Image(data, ref_img.affine, hdr)
    out.header['scl_slope'] = 1
    out.header['scl_inter'] = 0
    nib.save(out, path)
    return path


def n_frames(path):
    shape = nib.load(path).shape
    return int(shape[3]) if len(shape) > 3 and shape[3] > 1 else 1


def is_dynamic(path):
    return n_frames(path) > 1


def voxel_volume_ml(path):
    return float(np.prod(nib.load(path).header.get_zooms()[:3])) / 1000.0


def mask_bool(path, ref=None):
    """Boolean mask of a mask image (values > 0.5); checks the grid against `ref` (path) if given."""
    m_img, m = load(path)
    if m.ndim > 3:
        m = m[..., 0]
    if ref is not None:
        r = nib.load(ref)
        if r.shape[:3] != m.shape or not np.allclose(r.affine, m_img.affine, atol=1e-3):
            raise ValueError(f'mask {path} is not on the grid of {ref}')
    return m > 0.5


def head_radius_mm(data, zooms):
    """Radius (mm) of the sphere with the volume of the head (voxels above 10 % of the 99th percentile)."""
    thr = 0.1 * np.percentile(data[np.isfinite(data)], 99)
    vol = float((data > thr).sum()) * float(np.prod(zooms[:3]))
    return (3.0 * vol / (4.0 * math.pi)) ** (1.0 / 3.0) if vol > 0 else 50.0


# ---------------------------------------------------------------------------------------------------------------------
# Frame timing
# ---------------------------------------------------------------------------------------------------------------------
def frame_timing(info, nframes):
    """
    Frame timing from a BIDS-PET sidecar (FrameTimesStart, FrameDuration, s). The engineer's
    FrameReferenceTime (mid-frame times) is used when present, otherwise the frame middle.
    Returns None for a static image without timing.
    """
    start = info.get('FrameTimesStart')
    dur = info.get('FrameDuration')
    if start is None or dur is None:
        if nframes == 1:
            return None
        raise ValueError('dynamic PET without FrameTimesStart / FrameDuration in its JSON sidecar')
    start = np.atleast_1d(np.asarray(start, dtype=float))
    dur = np.atleast_1d(np.asarray(dur, dtype=float))
    if len(start) != nframes or len(dur) != nframes:
        raise ValueError(f'the JSON sidecar describes {len(start)} frames ({len(dur)} durations) '
                         f'but the image has {nframes}')
    if np.any(dur <= 0) or np.any(np.diff(start) < 0):
        raise ValueError('frame timing not increasing or non-positive frame duration')
    mid = info.get('FrameReferenceTime')
    mid = np.asarray(mid, dtype=float) if mid is not None and len(mid) == nframes else start + dur / 2.0
    return dict(start=start, duration=dur, end=start + dur, mid=mid)


def timing_to_json(timing):
    return {"FrameTimesStart": timing['start'].tolist(), "FrameDuration": timing['duration'].tolist(),
            "FrameReferenceTime": timing['mid'].tolist()}


def load_timing(json_file):
    info = read_json(json_file)
    if not info:
        return None
    return frame_timing(info, len(info.get('FrameTimesStart', [])) or 1)


def window_frames(timing, window_min):
    """Indices of the frames whose middle falls in window_min = [t0, t1] (minutes); None or 'all' = every frame."""
    n = len(timing['mid'])
    if window_min in (None, 'all', ''):
        return np.arange(n)
    t0, t1 = [60.0 * float(x) for x in window_min]
    idx = np.where((timing['mid'] >= t0) & (timing['mid'] <= t1))[0]
    if idx.size == 0:
        raise ValueError(f'no frame has its middle in the window {window_min} min '
                         f'(frames span {timing["start"][0] / 60:.1f}-{timing["end"][-1] / 60:.1f} min)')
    return idx


def weighted_mean(data4d, timing, frames=None):
    """Time average of the activity concentration over `frames` (weights = frame durations)."""
    if data4d.ndim == 3:
        return data4d
    if frames is None:
        frames = np.arange(data4d.shape[3])
    if timing is None:
        return data4d[..., frames].mean(axis=3)
    w = timing['duration'][frames]
    return np.tensordot(data4d[..., frames], w / w.sum(), axes=([3], [0]))


def write_mean(img_path, out_path, timing=None, window_min=None, description=None):
    """Duration-weighted mean (3D) of a dynamic image (copy for a static one)."""
    img, d = load(img_path)
    if d.ndim == 3:
        mean = d
        frames = None
    else:
        frames = window_frames(timing, window_min) if timing is not None else np.arange(d.shape[3])
        mean = weighted_mean(d, timing, frames)
    save_like(mean, img, out_path)
    write_json(out_path, img_path, description or 'Time average of the frames (weights = frame durations).',
               Frames=None if frames is None else [int(i) for i in frames])
    return out_path


# ---------------------------------------------------------------------------------------------------------------------
# Subject data for SUV (BIDS: pet sidecar + participants.tsv / sessions.tsv)
# ---------------------------------------------------------------------------------------------------------------------
_WEIGHT_COLS = ('body_weight', 'BodyWeight', 'weight', 'weight_kg')
_DOSE_COLS = ('injected_radioactivity', 'InjectedRadioactivity', 'tracor_qte_MBq', 'dose_MBq')
_GLU_MMOL_COLS = ('blood_glucose', 'glucose_mmol_L', 'glycemia_mmol_L')
_GLU_GL_COLS = ('blood_concentration', 'glucose_g_L')
GLUCOSE_MMOL_PER_G = 1.0 / 0.18016   # glucose, 180.16 g/mol


def _tsv_row(tsv, ID, Session):
    """Row of a participants.tsv / sessions.tsv for this subject (and session when the file has one)."""
    import pandas as pd
    if not opi(tsv):
        return None
    df = pd.read_csv(tsv, sep='\t', dtype=str)
    sub_cols = [c for c in ('participant_id', 'id', 'subject') if c in df.columns]
    if sub_cols:
        c = sub_cols[0]
        df = df[df[c].str.replace('sub-', '', regex=False) == str(ID)]
    ses_cols = [c for c in ('session_id', 'session_name', 'session') if c in df.columns]
    if ses_cols:
        c = ses_cols[0]
        df = df[df[c].str.replace('ses-', '', regex=False) == str(Session)]
    return None if df.empty else df.iloc[0]


def _first_number(row, cols):
    if row is None:
        return None, None
    for c in cols:
        if c in row.index:
            try:
                v = float(row[c])
            except (TypeError, ValueError):
                continue
            if np.isfinite(v):
                return v, c
    return None, None


def subject_values(bids_dir, ID, Session, pet_info):
    """
    Body weight (kg), injected radioactivity (MBq) and blood glucose (mmol/L) of one session.
    Sources, by priority: PET sidecar (InjectedRadioactivity[Units]), sub-<ID>_sessions.tsv,
    participants.tsv (BIDS columns, and the engineer's columns weight_kg / tracor_qte_MBq /
    blood_concentration in g/L). Missing values are None.
    """
    rows = [_tsv_row(opj(bids_dir, 'sub-' + str(ID), 'sub-' + str(ID) + '_sessions.tsv'), ID, Session),
            _tsv_row(opj(bids_dir, 'participants.tsv'), ID, Session)]
    src = {}
    weight = dose = glu = None
    for r, name in zip(rows, ('sessions.tsv', 'participants.tsv')):
        if weight is None:
            weight, c = _first_number(r, _WEIGHT_COLS)
            if weight is not None:
                src['weight'] = f'{name}:{c}'
        if dose is None:
            dose, c = _first_number(r, _DOSE_COLS)
            if dose is not None:
                src['dose'] = f'{name}:{c}'
        if glu is None:
            glu, c = _first_number(r, _GLU_MMOL_COLS)
            if glu is not None:
                src['glucose'] = f'{name}:{c} (mmol/L)'
            else:
                g, c = _first_number(r, _GLU_GL_COLS)
                if g is not None:
                    glu = g * GLUCOSE_MMOL_PER_G
                    src['glucose'] = f'{name}:{c} (g/L -> mmol/L)'
    if pet_info.get('InjectedRadioactivity') is not None:
        dose = activity_to_mbq(pet_info['InjectedRadioactivity'],
                               pet_info.get('InjectedRadioactivityUnits') or pet_info.get('InjectedRadioActivityUnits'))
        src['dose'] = 'PET sidecar: InjectedRadioactivity'
    if pet_info.get('BodyWeight') is not None and weight is None:
        weight = float(pet_info['BodyWeight'])
        src['weight'] = 'PET sidecar: BodyWeight'
    return dict(weight_kg=weight, dose_MBq=dose, glucose_mmol_L=glu, sources=src)


def suv_scale(pet_info, values, isotope):
    """
    Factor f such that SUV = f x image (body-weight SUV, g/mL, density 1 g/mL; Boellaard 2009, Kinahan &
    Fletcher 2010). The injected radioactivity is decayed to the decay-correction reference time of the
    image (BIDS ImageDecayCorrectionTime and InjectionStart, s relative to TimeZero); without that
    information the image is assumed decay-corrected to the injection time. Returns (f, report) or
    (None, reason).
    """
    w, d = values.get('weight_kg'), values.get('dose_MBq')
    if w is None or d is None:
        return None, 'body weight or injected radioactivity unknown: SUV not computed'
    unit = pet_info.get('Units')
    k = concentration_to_bq_per_ml(unit)
    report = {'weight_kg': w, 'dose_MBq': d, 'image_units': unit}
    if k is None:
        k = 1.0
        report['warning'] = f'image units {unit!r} not recognised: assumed Bq/mL'
    if pet_info.get('ImageDecayCorrected') is False:
        report['warning_decay'] = ('image NOT decay-corrected: it is decay-corrected to the injection time '
                                   'frame by frame before the SUV (decay_correct_frames)')
        dose_ref = d
    else:
        t_ref, t_inj = pet_info.get('ImageDecayCorrectionTime'), pet_info.get('InjectionStart')
        if t_ref is not None and t_inj is not None and isotope in HALF_LIFE_S:
            dose_ref = d * math.exp(-decay_constant(isotope) * (float(t_ref) - float(t_inj)))
            report['dose_decayed_to_reference_MBq'] = dose_ref
        else:
            dose_ref = d
            report['assumption'] = 'image assumed decay-corrected to the injection time'
    f = k * (w * 1000.0) / (dose_ref * 1e6)
    report['factor'] = f
    return f, report


def decay_correct_frames(data4d, timing, isotope, t_injection):
    """Decay correction of a non-decay-corrected dynamic image to the injection time
    (frame-averaged factor: lambda*dt / (1 - exp(-lambda*dt)) * exp(lambda*(t_start - t_inj)))."""
    lam = decay_constant(isotope)
    dt = timing['duration']
    f = np.exp(lam * (timing['start'] - t_injection)) * (lam * dt) / (1.0 - np.exp(-lam * dt))
    return data4d * f.astype(np.float32)


# ---------------------------------------------------------------------------------------------------------------------
# Motion parameters of an ANTs rigid / affine matrix
# ---------------------------------------------------------------------------------------------------------------------
def rigid_parameters(mat_file, radius_mm):
    """(tx, ty, tz) mm, (rx, ry, rz) deg (ZYX Euler) and mean displacement (mm) of the points of a sphere of
    radius `radius_mm` around the centre of rotation, for an ANTs transform x' = R (x - c) + c + t."""
    import ants
    tr = ants.read_transform(mat_file)
    p = np.asarray(tr.parameters, dtype=float)
    R, t = p[:9].reshape(3, 3), p[9:12]
    rx = math.degrees(math.atan2(R[2, 1], R[2, 2]))
    ry = math.degrees(math.asin(max(-1.0, min(1.0, -R[2, 0]))))
    rz = math.degrees(math.atan2(R[1, 0], R[0, 0]))
    golden = math.pi * (3.0 - math.sqrt(5.0))
    n = 200
    k = np.arange(n)
    z = 1 - 2 * (k + 0.5) / n
    r = np.sqrt(1 - z * z)
    pts = radius_mm * np.stack([r * np.cos(golden * k), r * np.sin(golden * k), z], axis=1)
    disp = np.linalg.norm(pts @ R.T + t - pts, axis=1).mean()
    return t, np.array([rx, ry, rz]), float(disp)


def now():
    return datetime.datetime.now().isoformat(timespec='seconds')
