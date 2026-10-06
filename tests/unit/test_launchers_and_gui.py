import sys
import pytest

from conftest import REPO

sys.path.insert(0, str(REPO / 'GUI'))
import ednix_gui_core as G  # noqa: E402

LAUNCHERS = sorted((REPO / 'Exemples' / 'Study_EDNiX').glob('Launcher_*/*.py'))


@pytest.mark.parametrize('launcher', [p for p in LAUNCHERS if 'Launcher_pet' not in str(p)], ids=lambda p: p.name)
def test_example_launchers_only_use_known_parameters(launcher):
    # preprocess_* accept **kwargs: a misspelled parameter would be silently ignored
    assert G.LauncherModel.check_calls(launcher.read_text(encoding='utf-8'), REPO) == []


@pytest.mark.parametrize('launcher', LAUNCHERS, ids=lambda p: p.name)
def test_gui_rewrites_launchers_without_changing_them(launcher):
    m = G.LauncherModel(launcher)
    src, n = m.render({f.key: f.value for f in m.fields})
    assert n == 0 and src == m.src


def test_label_file_formats(tmp_path):
    (tmp_path / 'wb.txt').write_text('L_A\n162 200 30 40 255\nR_A\n1162 200 30 40 255\n')
    (tmp_path / 'lut.txt').write_text('0 Unknown 0 0 0 0\n162 L_A 200 30 40 0\n1162 R_A 200 30 40 0\n')
    (tmp_path / 'itk.txt').write_text('    0 0 0 0 0 0 0 "Clear"\n  162 200 30 40 1 1 1 "L_A"\n 1162 200 30 40 1 1 1 "R_A"\n')
    for f in ('wb.txt', 'lut.txt', 'itk.txt'):
        labels = G.parse_label_text(tmp_path / f)
        assert labels[162]['name'] == 'L_A' and labels[1162]['hemi'] == 'right' and labels[162]['rgb'] == [200, 30, 40]


def test_pipeline_docstrings_document_most_parameters():
    docs = G.pipeline_docs(REPO)
    assert len(docs.get('preprocess_fmri', {})) > 50 and len(docs.get('preprocess_anat', {})) > 25


def test_qc_scan_table_and_path_safety(tmp_path):
    import pandas as pd
    qc = tmp_path / 'QC' / 'skullstrip_step1'
    qc.mkdir(parents=True)
    (qc / 'sub-01.png').write_bytes(b'\x89PNG')
    pd.DataFrame({'subject': [f'{i:02d}' for i in range(10)],
                  'mean_FD': [0.05, 0.06, 0.05, 0.04, 0.05, 0.06, 0.05, 0.05, 0.04, 0.9]}).to_csv(tmp_path / 'QC' / 'm.csv', index=False)
    scan = G.qc_scan(tmp_path)
    assert scan['steps'][0]['images'] == ['skullstrip_step1/sub-01.png'] and scan['tables'] == ['m.csv']
    t = G.qc_table(tmp_path, 'm.csv')
    assert t['n_flagged_rows'] == 1 and t['flags'][9][1]
    with pytest.raises(ValueError):
        G.qc_path(tmp_path, '../../etc/passwd')


def test_run_progress_is_read_from_the_log(tmp_path):
    log = tmp_path / 'run.log'
    log.write_text('##  Working on step 1 (function: _1_fMRI_preTTT_in_fMRIspace).  ##\nx\n'
                   '##  Working on step 7 (function: _7_post_TTT).  ##\n')
    rm = G.RunManager(tmp_path / 'state')
    assert rm.progress({'log': str(log)}) == dict(step='7', function='7_post_TTT', steps_seen=2)
