"""Reset cleaners (anat, PET, CT), QC reports (subject / all subjects / dataset + Excel) and surface QC, on a
synthetic BIDS tree."""
import json
import os

import numpy as np
import pytest

from Tools import QC_report as Q
from Tools import data_reset as dr
from modalities.anat import _000_Data_Reset as anat_reset
from modalities.pet import _000_Data_Reset as pet_reset

PNG = b'\x89PNG\r\n\x1a\n'


def _touch(p, data=b'x'):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(data)
    return p


def _nii(p, shape, zooms, tr=None):
    import nibabel as nib
    img = nib.Nifti1Image(np.zeros(shape, np.float32), np.diag(list(zooms[:3]) + [1]))
    img.header.set_zooms(zooms)
    if tr is not None:
        img.header.set_xyzt_units('mm', 'sec')
    p.parent.mkdir(parents=True, exist_ok=True)
    nib.save(img, str(p))


@pytest.fixture
def bids(tmp_path):
    b = tmp_path / 'BIDS'
    ses = b / 'sub-01' / 'ses-1'
    # raw BIDS
    _nii(ses / 'anat' / 'sub-01_ses-1_T1w.nii.gz', (8, 8, 8), (1, 1, 1))
    (ses / 'anat' / 'sub-01_ses-1_T1w.json').write_text('{"MagneticFieldStrength": 3}')
    _nii(ses / 'func' / 'sub-01_ses-1_task-rest_bold.nii.gz', (4, 4, 4, 10), (3, 3, 3, 2.0), tr=2.0)
    (ses / 'func' / 'sub-01_ses-1_task-rest_bold.json').write_text('{"RepetitionTime": 0.002}')   # ms in s field
    # EDNiX anat outputs + manual files
    a = ses / 'anat'
    _touch(a / '01 session 12026-11-10.txt', b'ERROR: something failed\nWARNING: hmm\n')
    _touch(a / 'preprocessing' / '01_space-acpc_desc-64_T1w.nii.gz')
    _touch(a / 'native' / 'volumes' / 'masks' / '01_final_mask.nii.gz')
    _touch(a / 'native' / 'volumes' / 'masks' / '01_final_mask_2.nii.gz')
    _touch(a / 'native' / 'volumes' / 'masks' / '01_backtonative_mask.nii.gz')
    _touch(a / 'native' / 'volumes' / 'masks' / '01_space-acpc_mask.nii.gz')
    _touch(a / 'native' / 'volumes' / '01_space-acpc_desc-SS_T1w.nii.gz')
    _touch(a / 'native' / 'surfaces' / 'Native_resol' / 'Exemple1.scene')
    _touch(a / 'native' / 'surfaces' / 'Native_resol' / '01.l.white.surf.gii')
    _touch(a / 'freesurfer' / '01' / 'mri' / 'segmentation.mgz')
    _touch(a / 'freesurfer' / '01' / 'mri' / 'brain.mgz')
    _touch(a / 'QC_anat' / 'T1w_QC_summary.png', PNG)
    (a / 'QC_anat' / '01_T1w_QC_values.json').write_text(json.dumps({'cnr': 2.5, 'cjv': 0.4}))
    # func outputs
    f = ses / 'func'
    _touch(f / 'raw' / 'masks' / 'manual_final_mask.nii.gz')
    _touch(f / 'QC' / 'sub-01_ses-1_task-rest_QC_summary.png', PNG)
    (f / 'QC' / 'sub-01_ses-1_task-rest_QC_values.json').write_text(json.dumps({'mean_fd': 0.1, 'TSNR': 50.0}))
    # PET
    p = ses / 'pet'
    _touch(p / 'sub-01_ses-1_trc-FDG_pet.nii.gz')
    _touch(p / 'sub-01_ses-1_trc-FDG_pet.json', b'{}')
    _touch(p / 'Prepro_Settings.json', b'{}')
    _touch(p / 'orig' / 'preprocessing' / 'forMC' / 'sub-01_ses-1_desc-crop_pet.nii.gz')
    _touch(p / 'orig' / 'preprocessing' / 'forMC' / 'sub-01_ses-1_desc-cropbox_mask.nii.gz')
    _touch(p / 'orig' / 'matrices' / 'x_0GenericAffine.mat')
    # dataset QC folders, names as written by the different QC writers
    _touch(b / 'QC' / 'skullstrip_step1' / '01_1_T1w_skullstriped.png', PNG)
    _touch(b / 'QC' / 'meanIMG_in_template' / '01meanIMG_in_template.png', PNG)
    _touch(b / 'QC' / 'meanIMG_in_template' / '011meanIMG_in_template.png', PNG)
    _touch(b / 'sub-011' / 'ses-1' / 'anat' / 'sub-011_ses-1_T1w.nii.gz')
    return b


def test_assign_images_longest_id_and_session(bids):
    ss = Q.sessions(str(bids))
    imgs = [p for v in Q._qc_check_images(str(bids)).values() for p in v]
    got = Q.assign_images(imgs, ss)
    names = {k: sorted(os.path.basename(p) for p in v) for k, v in got.items()}
    assert names[('01', '1')] == ['01_1_T1w_skullstriped.png', '01meanIMG_in_template.png']
    assert names[('011', '1')] == ['011meanIMG_in_template.png']


def test_anat_reset_keeps_raw_and_manual_files(bids):
    dp = str(bids / 'sub-01' / 'ses-1')
    root, to_delete, kept = anat_reset.plan_anat(dp)
    deleted = {os.path.relpath(p, root) for p in to_delete}
    kept_names = {os.path.basename(p) for p in kept}
    assert {'sub-01_ses-1_T1w.nii.gz', 'sub-01_ses-1_T1w.json', '01_final_mask.nii.gz', '01_final_mask_2.nii.gz',
            '01_backtonative_mask.nii.gz', 'segmentation.mgz', 'Exemple1.scene'} <= kept_names
    assert os.path.join('native', 'volumes', 'masks', '01_space-acpc_mask.nii.gz') in deleted
    assert '01 session 12026-11-10.txt' in deleted
    n = dr.execute('test', root, to_delete, kept, str(bids / 'diary.txt'), dry_run=False, skip_confirmation=True,
                   verbose=False)
    assert n == len(to_delete)
    assert (bids / 'sub-01' / 'ses-1' / 'anat' / 'native' / 'volumes' / 'masks' / '01_final_mask.nii.gz').exists()
    assert not (bids / 'sub-01' / 'ses-1' / 'anat' / 'preprocessing').exists()      # emptied folders pruned


def test_anat_reset_surfaces_scope(bids):
    root, to_delete, kept = anat_reset.plan_anat(str(bids / 'sub-01' / 'ses-1'), scope='surfaces')
    rel = {os.path.relpath(p, root) for p in to_delete}
    assert rel == {os.path.join('native', 'surfaces', 'Native_resol', '01.l.white.surf.gii'),
                   os.path.join('freesurfer', '01', 'mri', 'brain.mgz')}


def test_pet_reset_keeps_crop_and_optionally_matrices(bids):
    dp = str(bids / 'sub-01' / 'ses-1')
    root, to_delete, kept = pet_reset.plan_emission(dp)
    kn = {os.path.basename(p) for p in kept}
    assert kn == {'sub-01_ses-1_trc-FDG_pet.nii.gz', 'sub-01_ses-1_trc-FDG_pet.json',
                  'sub-01_ses-1_desc-crop_pet.nii.gz'}
    assert any(p.endswith('Prepro_Settings.json') for p in to_delete)
    _, to_delete2, kept2 = pet_reset.plan_emission(dp, keep_registration=True)
    assert any(p.endswith('x_0GenericAffine.mat') for p in kept2)


def test_reset_is_skipped_for_unattended_workers(monkeypatch):
    monkeypatch.setenv('EDNIX_SKIP_RESET', '1')
    assert dr.run_batch([dict(label='x', root='/nonexistent', to_delete=['/nonexistent/a'], kept=[],
                              diary_file='')], dry_run=False) == 0


def test_qc_build_three_levels(bids):
    pd = pytest.importorskip('pandas')
    out = Q.build(str(bids), surface=False, verbose=False)
    page = (bids / 'QC' / 'subjects' / 'sub-01_ses-1.html').read_text()
    assert '01_1_T1w_skullstriped.png' in page and 'mean_fd' in page and '1 error line' in page
    assert (bids / 'QC' / 'QC_all_subjects.html').exists() and (bids / 'QC' / 'QC_dataset.html').exists()
    acq = pd.read_excel(out['xlsx'], sheet_name='acquisition')
    func = acq[acq['modality'] == 'func'].iloc[0]
    assert 'out of range' in func['TR_check'] and func['TR_used_s'] == pytest.approx(2.0)
    assert func['duration_min'] == pytest.approx(10 * 2.0 / 60, abs=0.01)
    assert set(pd.ExcelFile(out['xlsx']).sheet_names) >= {'overview', 'subjects', 'acquisition', 'status', 'func_QC'}


def _sphere(r, center, n=24):
    th, ph = np.meshgrid(np.linspace(0.1, np.pi - 0.1, n), np.linspace(0, 2 * np.pi, n, endpoint=False), indexing='ij')
    v = np.c_[(r * np.sin(th) * np.cos(ph)).ravel(), (r * np.sin(th) * np.sin(ph)).ravel(), (r * np.cos(th)).ravel()]
    f = []
    for i in range(n - 1):
        for j in range(n):
            a, b, c, d = i * n + j, i * n + (j + 1) % n, (i + 1) * n + j, (i + 1) * n + (j + 1) % n
            f += [[a, b, c], [b, d, c]]
    return v + center, np.asarray(f, np.int32)


def test_surface_qc_on_synthetic_spheres(tmp_path):
    nib = pytest.importorskip('nibabel')
    pytest.importorskip('matplotlib')
    sp, n = 'S', 48
    mri, surf = tmp_path / sp / 'mri', tmp_path / sp / 'surf'
    mri.mkdir(parents=True)
    surf.mkdir(parents=True)
    ijk = np.indices((n, n, n)).reshape(3, -1).T - (n - 1) / 2
    rad = np.linalg.norm(ijk, axis=1).reshape(n, n, n)
    aff = np.diag([-1., 1., 1., 1.])
    filled = (rad <= 10).astype(np.uint8) * 127
    aseg = np.where(rad <= 10, 2, np.where(rad <= 14, 3, 0)).astype(np.int32)
    for name, arr in (('orig.mgz', (rad <= 14) * 100.), ('filled.mgz', filled), ('aseg.mgz', aseg)):
        nib.save(nib.MGHImage(np.asarray(arr, np.float32 if name == 'orig.mgz' else np.int32), aff), str(mri / name))
    tkr = nib.load(str(mri / 'orig.mgz')).header.get_vox2ras_tkr()
    center = (tkr @ np.r_[(n - 1) / 2, (n - 1) / 2, (n - 1) / 2, 1])[:3]
    for h in ('lh', 'rh'):
        nib.freesurfer.write_geometry(str(surf / f'{h}.white'), *_sphere(10, center))
        nib.freesurfer.write_geometry(str(surf / f'{h}.pial'), *_sphere(14, center))
        v = _sphere(10, center)[0]
        nib.freesurfer.write_morph_data(str(surf / f'{h}.thickness'), np.full(len(v), 4.0, np.float32))
    m = Q.surface_qc(str(tmp_path), sp, str(tmp_path / 'qc.png'), real_mm_per_fs_mm=0.25)
    assert (tmp_path / 'qc.png').exists() and (tmp_path / 'qc.json').exists()
    assert m['thick_median_mm'] == pytest.approx(1.0)              # 4 FS mm x 0.25
    assert m['white_within_1vox'] > 0.9 and m['pial_within_1vox'] > 0.9
    assert m['has_sphere_reg'] is False


def test_run_all_ednix_splits_launchers_by_subject(tmp_path):
    import importlib.util
    import subprocess
    import sys
    from conftest import REPO
    spec = importlib.util.spec_from_file_location('run_all_ednix', REPO / 'Exemples' / 'Study_EDNiX' / 'run_all_ednix.py')
    R = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(R)
    bids = tmp_path / 'BIDS'
    for s in ('a', 'b', 'c'):
        for ses in ('1', '2'):
            (bids / f'sub-{s}' / f'ses-{ses}' / 'anat').mkdir(parents=True)
    out = tmp_path / 'out.txt'
    launcher = tmp_path / 'Launcher_anat_X.py'
    launcher.write_text(f"bids_dir = {str(bids)!r}\nlist_to_keep = []\nlist_to_remove = [('c', '2')]\n"
                        f"creat_study_template = False\n"
                        f"open({str(out)!r}, 'a').write(repr(sorted(list_to_keep)) + '\\n')\n")
    L = R.parse_launcher(launcher, 'anat')
    assert L['bids'] == str(bids) and L['remove'] == [('c', '2')] and not L['template']
    jobs = R.make_jobs([L], tmp_path)
    assert len(jobs) == 3
    for j in jobs:
        subprocess.run([sys.executable, j['script']], check=True)
    got = sorted(out.read_text().split('\n')[:-1])
    assert got == ["[('a', '1'), ('a', '2')]", "[('b', '1'), ('b', '2')]", "[('c', '1')]"]
    # dependencies: func and pet of the same BIDS wait for the anat
    F = dict(L, kind='func', name='Launcher_func_X')
    P = dict(L, kind='pet', name='Launcher_pet_X')
    jobs = [dict(launcher=L), dict(launcher=F), dict(launcher=P)]
    assert R.dependencies(jobs) == {0: set(), 1: {0}, 2: {0}}
    # a study template keeps the launcher in one process
    T = dict(L, template=True)
    assert len(R.make_jobs([T], tmp_path)) == 1


def test_build_library_V2_discovery_and_manifest(tmp_path, monkeypatch):
    """Library V2: template discovery (freesurfer folders ignored) and checksum verification."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        'build_library_V2', os.path.join(os.path.dirname(__file__), '..', '..', 'Exemples', 'Study_EDNiX', 'atlases',
                                         'build_library_V2.py'))
    lib = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(lib)
    root = tmp_path / 'Atlases_library_V2'
    for sp in ('Rat', 'Macaque'):
        vol = root / 'atlas' / 'group' / sp / 'EDNiX' / 'volumes'
        vol.mkdir(parents=True)
        (vol / f'{sp}_space-acpc_desc-SS_T1w.nii.gz').write_bytes(b'x')
    fs = root / 'atlas' / 'group' / 'Rat' / 'EDNiX' / 'freesurfer' / 'Rat'
    fs.mkdir(parents=True)  # a folder named like a species under freesurfer/ must not be taken as a template
    found = lib.find_templates(str(root))
    assert [(s, r) for s, r, _ in found] == [('Macaque', 'EDNiX'), ('Rat', 'EDNiX')]
    assert [s for s, _, _ in lib.find_templates(str(root), ['Rat'])] == ['Rat']
    import json as _json
    (root / 'EDNiX_library_manifest.json').write_text(_json.dumps({'files': lib.checksums(str(root))}))
    assert lib.verify(str(root)) == []
    (root / 'atlas' / 'group' / 'Rat' / 'EDNiX' / 'volumes' / 'Rat_space-acpc_desc-SS_T1w.nii.gz').write_bytes(b'y')
    assert any(d.startswith('changed:') for d in lib.verify(str(root)))
