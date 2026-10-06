"""Run signal_regression end to end with AFNI / FSL commands emulated in numpy."""
import os, re, shlex, sys, shutil, tempfile, numpy as np, nibabel as nib
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import modalities.fMRI._7_post_TTT as N

def _load(p): return np.asarray(nib.load(p).dataobj, dtype=np.float32)
def _save(a, ref, p): nib.save(nib.Nifti1Image(a.astype(np.float32), nib.load(ref).affine), p)
def fourier_columns(n, dt, stop):
    """cos / sin regressors (no constant, no sine at Nyquist) of the Fourier frequencies f for which stop(f) is True."""
    t = np.arange(n); cols = []
    for k in range(1, n // 2 + 1):
        if stop(k / (n * dt)):
            cols.append(np.cos(2 * np.pi * k * t / n))
            if not (n % 2 == 0 and k == n // 2): cols.append(np.sin(2 * np.pi * k * t / n))
    return np.column_stack(cols) if cols else np.zeros((n, 0))

def _ts_bytes(ts): return ('\n'.join(' '.join(f'{v:.6f}' for v in np.atleast_1d(r)) for r in ts) + '\n').encode()

CALLS = []
def emulate(cmd, *a, **k):
    CALLS.append(cmd.split()[0])
    tok = shlex.split(cmd.replace('>', ' > '))
    prog = tok[0]
    arg = lambda f: tok[tok.index(f) + 1]
    if prog == '3dTstat': return b'', b''
    if prog == '3dmaskSVD':
        if not os.path.exists(arg('-mask')):
            return b'', b'** FATAL ERROR: cannot open mask'      # what AFNI does, logged by run_cmd
        d = _load(tok[-1]); m = _load(arg('-mask')) > 0
        return _ts_bytes(d[m].mean(axis=0)), b''
    if prog == '1dBport':
        # -band a b -invert: the regressors of the frequencies OUTSIDE the band (to be projected out)
        n = int(arg('-nodata')); dt = float(tok[tok.index('-nodata') + 2]); a, b = float(arg('-band')), float(tok[tok.index('-band') + 2])
        return _ts_bytes(fourier_columns(n, dt, lambda f: not (a <= f <= b))), b''
    if prog == '3dDeconvolve':
        orts = [tok[i + 1] for i, x in enumerate(tok) if x == '-ortvec']
        if any(os.path.getsize(o) <= 5 for o in orts):
            raise RuntimeError('3dDeconvolve: empty -ortvec file')
        cols = [np.atleast_2d(np.loadtxt(o)).reshape(-1, 1) if np.loadtxt(o).ndim == 1 else np.loadtxt(o) for o in orts]
        n = cols[0].shape[0]
        X = np.column_stack([np.ones(n), np.linspace(-1, 1, n)] + cols)
        np.savetxt(arg('-x1D_uncensored'), X); np.savetxt(arg('-x1D'), X); return b'', b''
    if prog == '3dTproject':
        d = _load(arg('-input')); n = d.shape[3]
        orts = [tok[i + 1] for i, x in enumerate(tok) if x == '-ort']
        dt = N._header_tr(arg('-input'))                       # 3dTproject reads the TR in the header
        filt = []
        if '-passband' in tok:
            a, b = float(arg('-passband')), float(tok[tok.index('-passband') + 2]); filt.append(fourier_columns(n, dt, lambda f: not (a <= f <= b)))
        if '-stopband' in tok:
            a, b = float(arg('-stopband')), float(tok[tok.index('-stopband') + 2]); filt.append(fourier_columns(n, dt, lambda f: a <= f <= b))
        X = np.column_stack([np.ones(n)] + filt + [np.loadtxt(o).reshape(n, -1) for o in orts])
        keep = np.loadtxt(arg('-censor')).astype(bool)
        Y = d.reshape(-1, n).T
        B = np.linalg.lstsq(X[keep], Y[keep], rcond=None)[0]
        R = Y - X @ B; R[~keep] = 0
        _save(R.T.reshape(d.shape), arg('-input'), arg('-prefix')); return b'', b''
    if prog == '1dcat':
        out = tok[tok.index('>') + 1]; files = tok[1:tok.index('>')]
        np.savetxt(out, np.column_stack([np.loadtxt(f).reshape(len(np.loadtxt(files[0])), -1) for f in files])); return b'', b''
    if prog == 'fslmeants':
        d = _load(arg('-i')); m = _load(arg('-m')) > 0; np.savetxt(arg('-o'), d[m].mean(axis=0)); return b'', b''
    if prog == 'fslmaths':
        if '-Tmean' in tok: _save(_load(tok[1]).mean(axis=3), tok[1], tok[-1])
        else: _save(_load(tok[1]) + _load(arg('-add'))[..., None], tok[1], tok[-1])
        return b'', b''
    if prog == '3dcalc':
        shutil.copy(arg('-a'), arg('-prefix')); return b'', b''
    raise RuntimeError('not emulated: ' + prog)

MSG = []


def install(monkeypatch):
    """Replace the command runner of step 7 by the emulator, for one test only."""
    monkeypatch.setattr(N.run_cmd, 'run', lambda cmd, diary, *a, **k: emulate(cmd))
    monkeypatch.setattr(N.run_cmd, 'do', lambda cmd, diary, *a, **k: emulate(cmd))
    monkeypatch.setattr(N.run_cmd, 'get', lambda cmd, diary, *a, **k: emulate(cmd))
    monkeypatch.setattr(N.run_cmd, 'msg', lambda m, d, c='ENDC': MSG.append((c, m)))

def make_subject(T=200, tr=2.0, seed=0, censor=True, header_tr=None):
    rng = np.random.default_rng(seed); root = tempfile.mkdtemp()
    dirs = {k: os.path.join(root, k) for k in ('proc', 'ica', 'masks', 'raw', 'post', 'mat')}
    for p in dirs.values(): os.makedirs(p)
    A = np.eye(4) * 2; A[3, 3] = 1
    shape = (16, 16, 10)
    brain = np.zeros(shape, bool); brain[2:14, 2:14, 1:9] = True
    wm = np.zeros(shape, bool); wm[5:11, 5:11, 3:7] = True
    csf = np.zeros(shape, bool); csf[7:9, 7:9, 4:6] = True
    motion = np.cumsum(rng.normal(0, 0.02, (T, 6)), axis=0)
    noise_ts = rng.normal(size=(T, 2))
    data = 1000 + rng.normal(0, 5, shape + (T,)) + (motion[:, :1].T * 50)[None, None] + \
        brain[..., None] * (noise_ts[:, 0] * 8)[None, None, None]
    data = data * brain[..., None]
    rs = 'sub-01_run-1_bold'
    for name in ('_space-acpc-func_desc-fMRI_run_inRef.nii.gz', '_space-acpc-func_desc-fMRI_run_inRef_SS.nii.gz'):
        img = nib.Nifti1Image(data.astype('f4'), A)
        img.header.set_zooms((2, 2, 2, tr if header_tr is None else header_tr)); img.header.set_xyzt_units('mm', 'sec')
        nib.save(img, os.path.join(dirs['proc'], rs + name))
    for n, m in (('mask_ref', brain), ('Wmask', wm), ('Vmask', csf)):
        nib.save(nib.Nifti1Image(m.astype('u1'), A), os.path.join(dirs['masks'], n + '.nii.gz'))
    keep = np.ones(T, int)
    if censor: keep[[10, 11, 50]] = 0
    np.savetxt(os.path.join(dirs['raw'], rs + '_space-func_desc-censor.1D'), keep, fmt='%d')
    np.savetxt(os.path.join(dirs['raw'], rs + '_space-func_desc-demean.1D'), motion - motion.mean(0))
    np.savetxt(os.path.join(dirs['raw'], rs + '_space-func_desc-deriv.1D'), np.vstack([np.zeros((1, 6)), np.diff(motion, axis=0)]))
    np.savetxt(os.path.join(dirs['mat'], rs + '_space-func_desc-motion_correction.1D'), motion)
    return dirs, rs, keep

def run(method, normalize='zscore', T=200, tr=2.0, band='0.01 0.1', GS=True, Vc=True, NonB=False, censor=True, header_tr=None, **kw):
    MSG.clear()
    dirs, rs, keep = make_subject(T, tr, censor=censor, header_tr=header_tr)
    N.signal_regression(dirs['proc'], dirs['ica'], dirs['masks'], dirs['raw'], 1, [rs + '.nii.gz'], 0, tr, 'Skip',
                        NonB, True, normalize, 'sub-01', method, False, band, Vc, GS,
                        dirs['post'], dirs['mat'], ' -overwrite', '', '', 'diary.txt', **kw)
    res = os.path.join(dirs['post'], rs + '_space-acpc-func_desc-fMRI_residual.nii.gz')
    out = _load(res) if os.path.exists(res) else None
    return out, keep, [m for c, m in MSG if c in ('FAIL', 'WARNING')]
