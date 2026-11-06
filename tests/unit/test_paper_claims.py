"""validation/paper_claims.py: text -> number -> stats cell -> code line, on a tiny manuscript."""
import pandas as pd
import pytest

docx = pytest.importorskip('docx')
from Exemples.Study_EDNiX.analysis.validation import paper_claims as PC

F3 = PC.opj('Figure_3', 'Figure_3_stats.xlsx')
LONG = 'Exemples/Study_EDNiX/analysis/library/longitudinal.py'


def _manuscript(path, text):
    d = docx.Document()
    for t in ('Results', text, 'Discussion', 'Methods', 'Isocortical trajectories, left hemisphere.'):
        d.add_paragraph(t)
    d.save(path)


def _stats(paper):
    (paper / 'Figure_3').mkdir(parents=True)
    pd.DataFrame([dict(species='Macaque', modality='thickness', kind='trajectory',
                       lme_slope_per_year=0.0212, lme_se_per_year=0.0049, lme_p_slope=0.0004,
                       n_subjects=8)]).to_excel(paper / F3, sheet_name='panels', index=False)


def _claim(id_, prefix, column, k=0):
    return dict(id=id_, prefix=prefix, k=k, producer=(LONG, 'figure3_stats', f'{column}='),
                source=dict(file=F3, sheet='panels', column=column,
                            where=dict(species='Macaque', modality='thickness', kind='trajectory')))


def test_claims_ok_mismatch_pending(tmp_path):
    ms = tmp_path / 'paper.docx'
    _manuscript(ms, 'In 8 macaques, thickness increased (linear slope 0.021 ± 0.006 mm per year, '
                    'P < 0.001). The cross-check gave r = [XX: r].')
    _stats(tmp_path)
    reg = [_claim('slope', '(linear slope ', 'lme_slope_per_year'),
           _claim('se', '(linear slope ', 'lme_se_per_year', k=1),
           _claim('p', 'mm per year, P ', 'lme_p_slope'),
           dict(_claim('pending', 'gave r = ', 'n_subjects')),
           _claim('gone', 'a sentence that was rewritten ', 'n_subjects')]
    rep = PC.check(str(ms), str(tmp_path), registry=reg, methods=[])
    st = rep['claims'].set_index('id')['status'].to_dict()
    assert st == dict(slope='OK', se='MISMATCH', p='OK', pending='PENDING', gone='NOT_IN_TEXT')
    code = rep['claims'].set_index('id').loc['slope', 'code']
    assert code.startswith(LONG + ':') and int(code.rsplit(':', 1)[1]) > 1000
    # the numbers registered are covered; '8' is not
    assert '8' in set(rep['untraced']['number'])
    assert '0.021' not in set(rep['untraced']['number'])


def test_methods_check(tmp_path):
    ms = tmp_path / 'paper.docx'
    _manuscript(ms, 'Nothing.')
    meth = [dict(id='hemi', quote='left hemisphere', kind='default', file=LONG, func='figure3_stats',
                 arg='hemisphere', expected='left'),
            dict(id='lvl', quote='Isocortical trajectories', kind='constant', file=LONG, name='ATLAS_LEVEL',
                 expected=2),
            dict(id='absent', quote='never written', kind='constant', file=LONG, name='MIN_SESSIONS', expected=2)]
    m = PC.check(str(ms), str(tmp_path), registry=[], methods=meth)['methods'].set_index('id')['status']
    assert m.to_dict() == dict(hemi='OK', lvl='FAIL', absent='NOT_IN_METHODS')


def test_registry_is_well_formed():
    from Exemples.Study_EDNiX.analysis.validation import paper_claims_registry as R
    ids = [c['id'] for c in R.CLAIMS]
    assert len(ids) == len(set(ids))
    for c in R.CLAIMS:
        assert c['prefix'] and 'source' in c and len(c.get('producer', ('', '', ''))) == 3
