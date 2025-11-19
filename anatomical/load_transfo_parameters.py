def build_transfos(
        align=None,
        intra=None,
        coreg=None,
        ss_coarse=None,
        ss_fine=None,
        ss_template=None,
        stdyT=None
):
    """
    Build all registration parameter blocks for the EDNiX pipeline.

    Each argument overrides the default parameters for that step.
    Only provide the keys you want to modify.

    Parameters
    ----------
    align : dict
        Native → template rigid alignment (ACPC-like). Typically: Rigid + MI.
    intra : dict
        Intra-subject alignment: e.g., T1 → T2, or M0 → T1 for some vendors.
    coreg : dict
        Coregistration of normalized anatomy to the high-resolution template.
    ss_coarse : dict
        Skull-stripping coarse pass (large mask).
    ss_fine : dict
        Skull-stripping fine pass (refined mask).
    ss_template : dict
        Skull stripping applied to the *study-level template*.
    stdyT : dict
        Registration of the study template → reference population template.

    Valid keys inside each block:
        'interpol'          : interpolation ("nearestNeighbor", "linear", "hammingWindowedSinc", ...)
        'type_of_transform' : ANTs transform ("Rigid", "Affine", "SyN", "SyNCC", ...)
        'affmetric'         : metric for affine stage ("MI", "mattes", "GC", "meansquare")
        'affmetricT'        : metric for nonlinear stage (optional, depends on type)

    Returns
    -------
    list of dict
        Ordered pipeline steps:
            [
                align,
                intra,
                coreg,
                ss_coarse,
                ss_fine,
                ss_template,
                stdyT
            ]
    """

    def block(name, defaults, override):
        cfg = defaults.copy()
        if override:
            cfg.update(override)
        cfg["name"] = name
        return cfg

    # --- Default sets ---
    base = {"interpol": "", "type_of_transform": "", "affmetric": "", "affmetricT": ""}

    align_def = {"interpol": "", "type_of_transform": "Rigid", "affmetric": "mattes", "affmetricT": "mattes"}

    coreg_def = {
        "interpol": "hammingWindowedSinc",
        "type_of_transform": "SyNCC",
        "affmetric": "MI",
        "affmetricT": ""
    }

    # --- Assemble blocks ---
    b_align      = block("align",       align_def, align)
    b_intra      = block("intra",       base,      intra)
    b_coreg      = block("coreg",       coreg_def, coreg)
    b_ss_coarse  = block("ss_coarse",   base,      ss_coarse)
    b_ss_fine    = block("ss_fine",     base,      ss_fine)
    b_ss_tpl     = block("ss_template", base,      ss_template)
    b_stdyT      = block("stdyT",       base,      stdyT)

    return [
        b_align,
        b_intra,
        b_coreg,
        b_ss_coarse,
        b_ss_fine,
        b_ss_tpl,
        b_stdyT
    ]