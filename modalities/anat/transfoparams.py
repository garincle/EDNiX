
default_intramod_type = 'SyNCC' # should be "SyN","SyNCC","SynBold","Rigid" or "Affine"
default_intermod_type = 'SyN'   # should be "SyN","SyNCC","SynBold","Rigid" or "Affine"
default_affmetric     = 'MI'    # should be "mattes" (means MI for mutual information),"GC" (global correlation) or "meansquare"
default_affmetricT    = 'mattes' # should be "mattes" (means MI for mutual information),"GC" (global correlation) or "meansquare"
default_interp        = 'hammingWindowedSinc' # should be 'nearestNeighbour','linear','hammingWindowedSinc'


# Names the pipeline looks up in list_transfo, and legacy aliases.
PIPELINE_NAMES = ('align', 'between', 'coreg', 'SS1', 'SS2', 'SS3', 'stdyT')
LEGACY_NAMES   = {'intra': 'between', 'ss_coarse': 'SS1', 'ss_fine': 'SS2', 'ss_template': 'SS3'}
KEYS           = ("interpol", "type_of_transform", "affmetric", "affmetricT")


def _defaults_for(name, coreg_block):
    """Default parameters of a block; SS1/SS2/SS3/stdyT inherit from 'coreg'."""
    if name in ('align', 'between'):
        return dict(zip(KEYS, (default_interp, 'Rigid', default_affmetric, default_affmetricT)))
    if name == 'coreg':
        return dict(zip(KEYS, (default_interp, default_intramod_type, default_affmetric, default_affmetricT)))
    return {k: coreg_block[k] for k in KEYS}


def set(transfo_message, IgotbothT1T2, creat_study_template, brain_skullstrip_1, brain_skullstrip_2,
        template_skullstrip, Align_img_to_template, list_transfo):
    """
    Return (Align_img_to_template, list_transfo) with one independent dict per
    pipeline step (align, between, coreg, SS1, SS2, SS3, stdyT), every empty
    field filled with its default. Legacy names (intra, ss_coarse, ss_fine,
    ss_template) are renamed; unknown names raise.
    """
    if Align_img_to_template == '':
        Align_img_to_template = 'Ants'

    if transfo_message == 'do_it_for_me':
        user = {}
    elif transfo_message == 'do_as_I_said':
        user = {}
        for d in list_transfo:
            name = LEGACY_NAMES.get(d["name"], d["name"])
            if name not in PIPELINE_NAMES:
                raise ValueError(f'Unknown transformation block "{d["name"]}". '
                                 f'Valid names: {PIPELINE_NAMES} (legacy: {tuple(LEGACY_NAMES)})')
            if name in user:
                raise ValueError(f'Transformation block "{name}" defined twice.')
            user[name] = {k: d.get(k, '') for k in KEYS}
    else:
        raise ValueError('transfo_message must be "do_it_for_me" or "do_as_I_said"')

    coreg = _defaults_for('coreg', None)
    coreg.update({k: v for k, v in user.get('coreg', {}).items() if v != ''})

    out = []
    for name in PIPELINE_NAMES:
        block = coreg.copy() if name == 'coreg' else _defaults_for(name, coreg)
        if name != 'coreg':
            block.update({k: v for k, v in user.get(name, {}).items() if v != ''})
        block["name"] = name
        out.append(block)

    return Align_img_to_template, out
