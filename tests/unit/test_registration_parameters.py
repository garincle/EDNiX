from modalities.anat import transfoparams as T
from modalities.anat.load_transfo_parameters import build_transfos


def _lookup(blocks, name):
    return [b for b in blocks if b['name'] == name]


def test_every_step_block_exists_and_is_independent():
    _, blocks = T.set('do_it_for_me', True, True, 'bet2_ANTS', 'bet2_ANTS', 'bet2_ANTS', '', [])
    assert len({id(b) for b in blocks}) == len(blocks)
    for name in T.PIPELINE_NAMES:
        assert len(_lookup(blocks, name)) == 1
    assert _lookup(blocks, 'coreg')[0]['type_of_transform'] == 'SyNCC'


def test_build_transfos_names_are_the_ones_the_pipeline_looks_up():
    lt = build_transfos(align={'type_of_transform': 'Rigid', 'affmetric': 'MI', 'affmetricT': 'MI'},
                        ss_coarse={'type_of_transform': 'SyN'})
    _, blocks = T.set('do_as_I_said', False, False, '', '', '', 'Ants', lt)
    assert _lookup(blocks, 'SS1')[0]['type_of_transform'] == 'SyN'
