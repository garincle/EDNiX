"""
ednix_atlas_regions.py
======================
EDNiX atlas region hierarchy — replaces the Legend xlsx as authoritative source.

STRUCTURE
---------
  §1  HIERARCHY LISTS         — by level and modality (surface / volume)
  §2  PARENT MAPS             — child → parent lookups at each transition
  §3  FUNCTIONAL GROUPINGS    — primary sensory/motor, associative, limbic, etc.
  §4  HELPER FUNCTIONS        — get_level_regions(), get_children(), get_parent(),
                                get_region_to_level(), get_functional_group()

USAGE
-----
    from ednix_atlas_regions import (
        get_level_regions,
        get_region_to_level,
        get_children,
        get_parent,
        FUNCTIONAL_GROUPS,
        get_functional_group,
        ABBREV,
    )

    # Region names at each level (cortical surface only)
    lvl1 = get_level_regions(1, cortical_only=True)
    lvl3 = get_level_regions(3, cortical_only=True)

    # Dict used by Stage 0 — base_name → atlas_level
    rtl = get_region_to_level()          # all regions (cortical + subcortical)
    rtl = get_region_to_level(cortical_only=True)

    # Hierarchy navigation
    children = get_children("Frontal_lobe")   # → ['oFC', 'PFC', 'Motor_premotor']
    parent   = get_parent("oFC")              # → 'Frontal_lobe'

    # Functional grouping
    group = get_functional_group("BA_1_3_S1") # → 'primary_sensorimotor'
"""

# ═══════════════════════════════════════════════════════════════════════════════
# §1  HIERARCHY LISTS
# ═══════════════════════════════════════════════════════════════════════════════

# ── Level 1 ───────────────────────────────────────────────────────────────────

level1_cortical = [
    "Isocortex",
    "Allocortex",
    "Periallocortex",
]

level1_subcortical = [
    "Subcortical_areas",
    "Diencephalon",
    "Cortical_White_matter",
]

# ── Level 2 ───────────────────────────────────────────────────────────────────

level2_cortical = [
    "Frontal_lobe",
    "Insular_lobe",
    "Temporal_lobe",
    "Parietal_lobe",
    "Occipital_lobe",
    "Paleocortex",
    "Archicortex",
    "Cortical_subplate",
    "Periarchicortex",
]

level2_subcortical = [
    "Cerebral_nuclei",
    "Diencephalon",
    "Brainstem",
    "Cerebellum",
    "Cerebellum_White",
    "White_matter",
    "Lateral_ventricules",
    "3rde_ventricules",
    "CSF",
]

# ── Level 3 ───────────────────────────────────────────────────────────────────

level3_cortical = [
    # Frontal
    "oFC",
    "PFC",
    "Motor_premotor",
    # Parietal
    "Somatosensory_Cx",
    "Post_parietal_Cx",
    "PMC",
    # Occipital
    "Striate_Cx",
    "Extra_striate_Cx",
    # Temporal
    "Auditory_Cx",
    "MIPT",
    "Vent_temporal",
    # Insular
    "Insula",
    # Allocortex / Paleocortex
    "Olfactory_Cx",
    "Olfactory_bulb",
    "Hippocampal_form",
    "Subiculum",
    "Periarchicortex_L3",
    # Cortical subplate
    "Claustrum",
    "Amygdala",
]

level3_subcortical = [
    "Septum",
    "Striatum",
    "Basal_forebrain",
    "Pallidum",
    "Hypothalamus",
    "Thalamus",
    "Epithalamus",
    "Subthalamus",
    "Tectum",
    "Tegmentum",
    "Cerebral_ped",
    "Pons",
    "Medulla",
    "Colliculus",
]

# ── Level 4 ───────────────────────────────────────────────────────────────────

level4_cortical = [
    # oFC
    "Gustatory_Cx",
    "oProiso",
    # PFC
    "BA_10", "BA_11", "BA_13", "BA_12_47", "IFJ",
    "BA_9", "BA_46", "BA_8", "BA_55",
    "BA_32", "BA_14", "BA_24", "BA_25", "BA_33",
    "BA_47", "BA_45", "BA_44",
    # Motor/premotor
    "BA_6", "BA_4_M1",
    # Somatosensory
    "BA_1_3_S1", "BA_43", "S2",
    # Post-parietal
    "BA_5", "BA_40", "BA_7", "BA_39",
    # PMC (posteromedial cortex)
    "PGM", "BA_23", "BA_7m", "BA_31",
    # Striate / extra-striate
    "BA_17_V1",
    "BA_18_V2", "BA_19",
    # Auditory
    "BA_41_42_A1", "Belt_A2", "Parabelt_A3", "RTp", "BA_22",
    # MIPT / inferotemporal
    "BA_21", "BA_20", "BA_38",
    # Ventral temporal
    "Parahippocampal", "BA_37",
    # Insula
    "Ri", "Insular_Cx", "BA_52",
    # Olfactory / paleocortex
    "BA_27", "Olfactory_tubercle",
    # Olfactory bulb
    "OB",
    # Hippocampal formation
    "CA", "DG", "BA_28",
    # Periarchicortex
    "BA_35_36", "Presubiculum", "BA_26_29_30",
    # Cortical subplate / amygdala complex
    "Bed_n_ST",
]

level4_subcortical = [
    # Striatum
    "Caudate", "Putamen", "N_accumbens",
    # Pallidum
    "Globus_pallidus",
    # Thalamus
    "Pulvinar", "Geniculate_n", "Thalamic_nuclei",
    # Tectum
    "Sup_colliculus", "Inf_colliculus",
    # Tegmentum
    "PAG", "Red_n", "Subst_nigra",
    # Cerebral peduncle
    "Crus_cerebri",
]

# ═══════════════════════════════════════════════════════════════════════════════
# §2  PARENT MAPS  (child → parent at next coarser level)
# ═══════════════════════════════════════════════════════════════════════════════

# lvl2 → lvl1
_PARENT_L2 = {
    # cortical → Isocortex
    "Frontal_lobe":       "Isocortex",
    "Insular_lobe":       "Isocortex",
    "Temporal_lobe":      "Isocortex",
    "Parietal_lobe":      "Isocortex",
    "Occipital_lobe":     "Isocortex",
    # cortical → Allocortex / Periallocortex
    "Paleocortex":        "Allocortex",
    "Archicortex":        "Allocortex",
    "Cortical_subplate":  "Periallocortex",
    "Periarchicortex":    "Periallocortex",
    # subcortical → Subcortical_areas / Diencephalon
    "Cerebral_nuclei":    "Subcortical_areas",
    "Brainstem":          "Subcortical_areas",
    "Cerebellum":         "Subcortical_areas",
    "Cerebellum_White":   "Subcortical_areas",
    "White_matter":       "Cortical_White_matter",
    "Lateral_ventricules":"Subcortical_areas",
    "3rde_ventricules":   "Subcortical_areas",
    "CSF":                "Subcortical_areas",
    "Diencephalon":       "Diencephalon",   # maps to itself at lvl1
}

# lvl3 → lvl2
_PARENT_L3 = {
    # Frontal
    "oFC":              "Frontal_lobe",
    "PFC":              "Frontal_lobe",
    "Motor_premotor":   "Frontal_lobe",
    # Parietal
    "Somatosensory_Cx": "Parietal_lobe",
    "Post_parietal_Cx": "Parietal_lobe",
    "PMC":              "Parietal_lobe",
    # Occipital
    "Striate_Cx":       "Occipital_lobe",
    "Extra_striate_Cx": "Occipital_lobe",
    # Temporal
    "Auditory_Cx":      "Temporal_lobe",
    "MIPT":             "Temporal_lobe",
    "Vent_temporal":    "Temporal_lobe",
    # Insular
    "Insula":           "Insular_lobe",
    # Paleocortex
    "Olfactory_Cx":     "Paleocortex",
    "Olfactory_bulb":   "Paleocortex",
    # Archicortex
    "Hippocampal_form": "Archicortex",
    "Subiculum":        "Archicortex",
    # Periarchicortex
    "Periarchicortex_L3": "Periarchicortex",
    # Cortical subplate
    "Claustrum":        "Cortical_subplate",
    "Amygdala":         "Cortical_subplate",
    # Subcortical
    "Septum":           "Cerebral_nuclei",
    "Striatum":         "Cerebral_nuclei",
    "Basal_forebrain":  "Cerebral_nuclei",
    "Pallidum":         "Cerebral_nuclei",
    "Hypothalamus":     "Diencephalon",
    "Thalamus":         "Diencephalon",
    "Epithalamus":      "Diencephalon",
    "Subthalamus":      "Diencephalon",
    "Tectum":           "Brainstem",
    "Tegmentum":        "Brainstem",
    "Cerebral_ped":     "Brainstem",
    "Pons":             "Brainstem",
    "Medulla":          "Brainstem",
    "Colliculus":       "Brainstem",
}

# lvl4 → lvl3
_PARENT_L4 = {
    # oFC
    "Gustatory_Cx":       "oFC",
    "oProiso":            "oFC",
    # PFC
    "BA_10":  "PFC", "BA_11": "PFC", "BA_13": "PFC", "BA_12_47": "PFC",
    "IFJ":    "PFC", "BA_9":  "PFC", "BA_46": "PFC", "BA_8":     "PFC",
    "BA_55":  "PFC", "BA_32": "PFC", "BA_14": "PFC", "BA_24":    "PFC",
    "BA_25":  "PFC", "BA_33": "PFC", "BA_47": "PFC", "BA_45":    "PFC",
    "BA_44":  "PFC",
    # Motor/premotor
    "BA_6":         "Motor_premotor",
    "BA_4_M1":      "Motor_premotor",
    # Somatosensory
    "BA_1_3_S1":    "Somatosensory_Cx",
    "BA_43":        "Somatosensory_Cx",
    "S2":           "Somatosensory_Cx",
    # Post-parietal
    "BA_5":  "Post_parietal_Cx", "BA_40": "Post_parietal_Cx",
    "BA_7":  "Post_parietal_Cx", "BA_39": "Post_parietal_Cx",
    # PMC
    "PGM":   "PMC", "BA_23": "PMC", "BA_7m": "PMC", "BA_31": "PMC",
    # Striate / extra-striate
    "BA_17_V1":     "Striate_Cx",
    "BA_18_V2":     "Extra_striate_Cx",
    "BA_19":        "Extra_striate_Cx",
    # Auditory
    "BA_41_42_A1":  "Auditory_Cx",
    "Belt_A2":      "Auditory_Cx",
    "Parabelt_A3":  "Auditory_Cx",
    "RTp":          "Auditory_Cx",
    "BA_22":        "Auditory_Cx",
    # MIPT / inferotemporal
    "BA_21": "MIPT", "BA_20": "MIPT", "BA_38": "MIPT",
    # Ventral temporal
    "Parahippocampal": "Vent_temporal",
    "BA_37":           "Vent_temporal",
    # Insula
    "Ri":         "Insula",
    "Insular_Cx": "Insula",
    "BA_52":      "Insula",
    # Olfactory / paleocortex
    "BA_27":              "Olfactory_Cx",
    "Olfactory_tubercle": "Olfactory_Cx",
    # Olfactory bulb
    "OB": "Olfactory_bulb",
    # Hippocampal formation
    "CA": "Hippocampal_form", "DG": "Hippocampal_form", "BA_28": "Hippocampal_form",
    # Periarchicortex
    "BA_35_36":    "Periarchicortex_L3",
    "Presubiculum":"Periarchicortex_L3",
    "BA_26_29_30": "Periarchicortex_L3",
    # Cortical subplate / amygdala
    "Bed_n_ST": "Amygdala",
    # Subcortical
    "Caudate":        "Striatum",
    "Putamen":        "Striatum",
    "N_accumbens":    "Striatum",
    "Globus_pallidus":"Pallidum",
    "Pulvinar":       "Thalamus",
    "Geniculate_n":   "Thalamus",
    "Thalamic_nuclei":"Thalamus",
    "Sup_colliculus": "Tectum",
    "Inf_colliculus": "Tectum",
    "PAG":            "Tegmentum",
    "Red_n":          "Tegmentum",
    "Subst_nigra":    "Tegmentum",
    "Crus_cerebri":   "Cerebral_ped",
}

# ═══════════════════════════════════════════════════════════════════════════════
# §3  FUNCTIONAL GROUPINGS
# ═══════════════════════════════════════════════════════════════════════════════
# Four orthogonal axes, each a dict {group_name: [region_list]}
# Regions not listed in a given axis get group "other" / "unclassified".

# ── 3a  Primary vs Integrative (unimodal → heteromodal gradient) ──────────────
FUNCTIONAL_GROUPS = {

    # PRIMARY unimodal sensory/motor cortices
    "primary_sensorimotor": [
        "BA_17_V1",         # primary visual
        "BA_41_42_A1",      # primary auditory
        "BA_1_3_S1",        # primary somatosensory
        "BA_4_M1",          # primary motor
        "OB",               # primary olfactory relay
        "Striate_Cx",       # lvl3 alias for V1 cluster
        "Somatosensory_Cx", # lvl3 alias
        "Motor_premotor",   # lvl3 alias (includes M1)
        "Auditory_Cx",      # lvl3 alias
        "Olfactory_bulb",   # lvl3 alias
    ],

    # UNIMODAL ASSOCIATION — one modality, higher processing
    "unimodal_association": [
        "BA_18_V2", "BA_19",            # visual association
        "Belt_A2", "Parabelt_A3", "RTp","BA_22",  # auditory belt/parabelt
        "S2", "BA_43",                  # secondary somatosensory
        "BA_6",                         # premotor
        "BA_5",                         # somatosensory association
        "Extra_striate_Cx",             # lvl3 alias
    ],

    # HETEROMODAL / HIGH-ORDER ASSOCIATION
    "heteromodal_association": [
        # Parietal
        "BA_40", "BA_7", "BA_39",       # IPL, SPL, angular gyrus
        "Post_parietal_Cx",             # lvl3 alias
        # Temporal
        "BA_21", "BA_20", "BA_38",      # MTG, ITG, temporal pole
        "BA_37",                        # fusiform / posterior IT
        "Parahippocampal",
        "MIPT",                         # lvl3 alias
        "Vent_temporal",                # lvl3 alias
        # Prefrontal lateral
        "BA_9", "BA_46", "BA_8",        # DLPFC
        "BA_45", "BA_44", "BA_47",      # IFG
        "IFJ",
        # Insula
        "Ri", "Insular_Cx", "BA_52",
        "Insula",                       # lvl3 alias
    ],

    # TRANSMODAL / DEFAULT MODE
    "transmodal_default_mode": [
        # Medial prefrontal
        "BA_10", "BA_11", "BA_14",      # frontal pole / vmPFC
        "BA_25",                        # subgenual ACC
        "BA_32", "BA_24",               # ACC
        "BA_33",
        # Posteromedial
        "PGM", "BA_23", "BA_7m", "BA_31",
        "PMC",                          # lvl3 alias (precuneus / RSC)
        # Angular gyrus (overlap with heteromodal — listed here for DMN emphasis)
        # Medial temporal
        "CA", "DG", "BA_28",            # hippocampus
        "Hippocampal_form",             # lvl3 alias
        "BA_35_36", "Presubiculum", "BA_26_29_30",
        "Periarchicortex_L3",           # lvl3 alias
        # Temporal pole
        "BA_38",
    ],

    # LIMBIC / PARALIMBIC
    "limbic": [
        "Amygdala", "Bed_n_ST",
        "Hippocampal_form",             # lvl3 alias
        "Subiculum",
        "Olfactory_Cx", "BA_27", "Olfactory_tubercle",
        "Paleocortex",                  # lvl2 alias
        "Archicortex",                  # lvl2 alias
        "Periarchicortex",              # lvl2 alias
        "Cortical_subplate",            # lvl2 alias
        "Claustrum",
        "Basal_forebrain",
        "Septum",
        "Hypothalamus",
    ],

    # FRONTAL EXECUTIVE / CONTROL
    "frontal_executive": [
        "BA_10", "BA_11", "BA_13",      # frontal pole / OFC
        "BA_12_47", "BA_55",
        "oFC", "PFC",                   # lvl3 aliases
        "oProiso", "Gustatory_Cx",
        "BA_32", "BA_24",               # ACC (also DMN)
        "BA_9", "BA_46", "BA_8",        # DLPFC (overlap heteromodal)
        "BA_44", "BA_45", "BA_47",      # IFG / Broca
        "IFJ",
        "Frontal_lobe",                 # lvl2 alias
    ],

    # SUBCORTICAL CIRCUITS
    "subcortical_basal_ganglia": [
        "Caudate", "Putamen", "N_accumbens",
        "Globus_pallidus",
        "Striatum", "Pallidum",         # lvl3 aliases
        "Subst_nigra", "Red_n",         # midbrain
    ],

    "subcortical_thalamus": [
        "Pulvinar", "Geniculate_n", "Thalamic_nuclei",
        "Thalamus",                     # lvl3 alias
        "Epithalamus", "Subthalamus",
    ],

    "subcortical_brainstem_cerebellum": [
        "Tectum", "Tegmentum", "Cerebral_ped", "Pons", "Medulla",
        "Sup_colliculus", "Inf_colliculus",
        "PAG", "Red_n", "Subst_nigra", "Crus_cerebri",
        "Cerebellum", "Cerebellum_White",
        "Brainstem",                    # lvl2 alias
    ],
}

# ── 3b  Cortical type axis (agranular → granular / Mesulam) ───────────────────
CORTICAL_TYPE = {
    "agranular": [
        # Motor, premotor, anterior cingulate — no layer IV
        "BA_4_M1", "BA_6", "Motor_premotor",
        "BA_24", "BA_25", "BA_33", "BA_32",
        "Insular_Cx", "Ri",             # anterior insula
    ],
    "dysgranular": [
        # Transitional — thin layer IV
        "BA_43", "S2", "Somatosensory_Cx",
        "BA_52", "Insula",
        "Olfactory_Cx", "Paleocortex",
        "Amygdala", "Cortical_subplate",
    ],
    "granular_eulaminate": [
        # Full 6-layer — primary + unimodal + heteromodal
        "BA_1_3_S1", "BA_17_V1", "BA_41_42_A1",
        "BA_18_V2", "BA_19", "BA_22", "Belt_A2", "Parabelt_A3",
        "BA_5", "BA_7", "BA_40", "BA_39",
        "BA_9", "BA_46", "BA_8", "BA_45", "BA_44",
        "BA_21", "BA_20", "BA_37",
        "PGM", "BA_23", "BA_7m", "BA_31",
        "BA_10", "BA_11",
    ],
    "paralimbic": [
        # Mesocortex / periallocortex
        "CA", "DG", "BA_28", "Subiculum",
        "BA_35_36", "Presubiculum", "BA_26_29_30",
        "Periarchicortex_L3", "Periarchicortex",
        "Archicortex", "Hippocampal_form",
    ],
}

# ── 3c  Sensory modality axis ─────────────────────────────────────────────────
SENSORY_MODALITY = {
    "visual": [
        "BA_17_V1", "BA_18_V2", "BA_19",
        "Striate_Cx", "Extra_striate_Cx",
        "Pulvinar", "Geniculate_n",
        "Sup_colliculus",
    ],
    "auditory": [
        "BA_41_42_A1", "Belt_A2", "Parabelt_A3", "RTp", "BA_22",
        "Auditory_Cx",
        "Inf_colliculus",
    ],
    "somatomotor": [
        "BA_1_3_S1", "BA_43", "S2",
        "BA_4_M1", "BA_6",
        "BA_5",
        "Somatosensory_Cx", "Motor_premotor",
    ],
    "olfactory": [
        "OB", "Olfactory_bulb",
        "Olfactory_Cx", "BA_27", "Olfactory_tubercle",
        "Paleocortex",
    ],
    "gustatory": [
        "Gustatory_Cx",
        "BA_43",                         # secondary taste (overlap somatomotor)
    ],
    "multimodal": [
        "Insula", "Ri", "Insular_Cx", "BA_52",
        "BA_38", "BA_39", "BA_40",
        "Post_parietal_Cx",
        "Thalamic_nuclei",
    ],
}

# ── 3d  Association / integration axis (simple primary→transmodal gradient) ───
# Useful for plotting a single scalar "integration level" per region.
INTEGRATION_LEVEL = {
    # 1 = primary unimodal
    1: ["BA_17_V1", "BA_41_42_A1", "BA_1_3_S1", "BA_4_M1", "OB",
        "Striate_Cx", "Somatosensory_Cx", "Motor_premotor", "Auditory_Cx",
        "Olfactory_bulb"],
    # 2 = unimodal association
    2: ["BA_18_V2", "BA_19", "Belt_A2", "Parabelt_A3", "RTp", "BA_22",
        "S2", "BA_43", "BA_6", "BA_5",
        "Extra_striate_Cx"],
    # 3 = heteromodal association
    3: ["BA_40", "BA_7", "BA_39", "BA_9", "BA_46", "BA_8",
        "BA_21", "BA_20", "BA_37", "BA_45", "BA_44", "BA_47", "IFJ",
        "Ri", "Insular_Cx", "BA_52", "Parahippocampal",
        "Post_parietal_Cx", "MIPT", "Vent_temporal", "Insula"],
    # 4 = transmodal / default mode
    4: ["BA_10", "BA_11", "BA_14", "BA_25", "BA_32", "BA_24",
        "PGM", "BA_23", "BA_7m", "BA_31",
        "CA", "DG", "BA_28", "BA_35_36", "Presubiculum", "BA_26_29_30",
        "BA_38",
        "PMC", "Hippocampal_form", "Periarchicortex_L3"],
}

# ═══════════════════════════════════════════════════════════════════════════════
# §4  HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def get_level_regions(atlas_level: int, cortical_only: bool = False,
                      subcortical_only: bool = False) -> list:
    """
    Return region names at a given atlas level (1–4).

    Parameters
    ----------
    atlas_level     : 1–4
    cortical_only   : if True, return only cortical regions
    subcortical_only: if True, return only subcortical regions

    Returns
    -------
    list of str  (region base names, no hemisphere prefix)
    """
    _map = {
        1: (level1_cortical,  level1_subcortical),
        2: (level2_cortical,  level2_subcortical),
        3: (level3_cortical,  level3_subcortical),
        4: (level4_cortical,  level4_subcortical),
    }
    if atlas_level not in _map:
        raise ValueError(f"atlas_level must be 1–4, got {atlas_level}")
    cortical, subcortical = _map[atlas_level]
    if cortical_only:   return list(cortical)
    if subcortical_only: return list(subcortical)
    return list(cortical) + list(subcortical)


def get_region_to_level(cortical_only: bool = False,
                        subcortical_only: bool = False) -> dict:
    """
    Return {base_name: atlas_level} for all regions.
    Regions that appear at multiple levels are assigned their COARSEST level.
    This is what Stage 0 needs to know which dlabel map to query.

    Compatible drop-in for atlas_config['region_to_level'].
    """
    rtl = {}
    for lvl in range(1, 5):
        for r in get_level_regions(lvl, cortical_only=cortical_only,
                                   subcortical_only=subcortical_only):
            if r not in rtl:   # first occurrence = coarsest
                rtl[r] = lvl
    return rtl


def get_parent(region: str) -> str | None:
    """Return the parent region at the next coarser atlas level, or None."""
    for lut in (_PARENT_L4, _PARENT_L3, _PARENT_L2):
        if region in lut:
            return lut[region]
    return None


def get_children(region: str, atlas_level: int | None = None) -> list:
    """
    Return all immediate children of a region.

    Parameters
    ----------
    region      : parent region name
    atlas_level : if given, children are restricted to that level;
                  otherwise all finer-level children are returned.
    """
    # Build reverse lookup lazily
    _luts = [(_PARENT_L2, 2), (_PARENT_L3, 3), (_PARENT_L4, 4)]
    children = []
    for lut, child_lvl in _luts:
        if atlas_level is not None and child_lvl != atlas_level:
            continue
        for child, parent in lut.items():
            if parent == region:
                children.append(child)
    return children


def get_functional_group(region: str) -> str | None:
    """
    Return the FUNCTIONAL_GROUPS key for a region, or None if unclassified.
    If a region belongs to multiple groups (e.g. BA_38 in DMN and heteromodal),
    the first match in dict order is returned.
    """
    for group, members in FUNCTIONAL_GROUPS.items():
        if region in members:
            return group
    return None


def get_integration_level(region: str) -> int | None:
    """Return the integration level (1–4) for a region, or None."""
    for lvl, members in INTEGRATION_LEVEL.items():
        if region in members:
            return lvl
    return None


def get_cortical_type(region: str) -> str | None:
    """Return the cortical laminar type for a region, or None."""
    for ctype, members in CORTICAL_TYPE.items():
        if region in members:
            return ctype
    return None


def get_sensory_modality(region: str) -> str | None:
    """Return the primary sensory modality label, or None."""
    for mod, members in SENSORY_MODALITY.items():
        if region in members:
            return mod
    return None


def build_atlas_config_from_py(cortical_only: bool = False) -> dict:
    """
    Build an atlas_config-compatible dict from this module —
    drop-in replacement for build_atlas_config(label_path, legend_path).

    Returns the same structure expected by _process_one_subject:
        {
          'region_to_level': {base_name: int},
          'name_to_ids':     {} (empty — IDs come from dlabel at runtime),
          'legend':          {'level1': [...], ..., 'region_to_level': {...}},
        }

    Usage in run_wb_extraction (no legend xlsx needed):
        from ednix_atlas_regions import build_atlas_config_from_py
        atlas_config = build_atlas_config_from_py()
    """
    rtl = get_region_to_level(cortical_only=cortical_only)
    legend = {f"level{i}": get_level_regions(i, cortical_only=cortical_only)
              for i in range(1, 5)}
    legend["region_to_level"] = rtl
    return {
        "label_df":        None,
        "id_to_info":      {},
        "name_to_ids":     {},
        "legend":          legend,
        "region_to_level": rtl,
    }


# ── Abbreviations (for plot tick labels) ──────────────────────────────────────
ABBREV = {
    "Paleocortex": "Paleo", "White matter": "WM",
    "Archicortex": "Archi", "Periarchicortex": "Peri",
    "Cerebral nuclei": "CN", "Cortical subplate": "CS",
    "Diencephalon": "Dien.", "Cerebellum": "Cbl.",
    "Brainstem": "BS", "Lateral ventricules": "LV",
    "Cerebellum White": "CbW", "3rde ventricules": "3V",
    "Motor/premotor": "M/PM", "Somatosensory Cx": "SmCx",
    "Post. parietal Cx": "PPC", "Visual striate Cx": "V1",
    "Visual extra-striate Cx": "V2+", "Auditory Cx": "AC",
    "Vent. temporal": "VT", "Olfactory Cx": "OC",
    "Olfactory bulb": "OB", "Hippocampal form.": "HF",
    "Basal forebrain": "BF", "Cerebral ped.": "CP",
    "Hypothalamus": "Hyp.", "Epithalamus": "Epi.",
    "Subthalamus": "Sub.", "Gustatory Cx": "GC",
    "BA 1-3 (S1)": "S1", "BA 4 (M1)": "M1",
    "BA 41/42 (A1)": "A1", "Belt (A2)": "A2",
    "Parabelt (A3)": "A3", "Parahippocampal": "PHC",
    "Insular Cx": "InsCx", "N. accumbens": "NAc",
    "Globus pallidus": "GP", "Bed n. ST": "BST",
    "Geniculate n.": "GN", "Thalamic nuclei": "TN",
    "Sup. colliculus": "SC", "Inf. colliculus": "IC",
    "Subst. nigra": "SN", "Crus cerebri": "CC",
    "BA 17 (V1)": "V1", "BA 18 (V2)": "V2",
    "BA 12/47": "12/47", "BA 26/29/30": "26-30",
    "BA 35/36": "35/36", "Insula & lat. sulcus": "Ins&LS",
    "oFC (orb. frontal Cx)": "oFC",
    "Prefrontal Cx": "PFC",
    "Olfactory bulb (OB)": "OB", "Olfactory tubercle": "OT",
    # Key base-name abbreviations (for axis labels)
    "Post_parietal_Cx": "PPC", "Somatosensory_Cx": "SmCx",
    "Motor_premotor": "M/PM", "Extra_striate_Cx": "V2+",
    "Striate_Cx": "V1cx", "Vent_temporal": "VT",
    "Olfactory_Cx": "OC", "Olfactory_bulb": "OB",
    "Hippocampal_form": "HF", "Periarchicortex_L3": "Peri",
    "Cortical_subplate": "CS", "Cerebral_nuclei": "CN",
    "Frontal_lobe": "FL", "Parietal_lobe": "PL",
    "Temporal_lobe": "TL", "Occipital_lobe": "OL",
    "Insular_lobe": "Ins", "Isocortex": "Iso",
    "Allocortex": "Allo", "Periallocortex": "Peri",
    "Subcortical_areas": "Subc.", "Cortical_White_matter": "WM",
}


# ═══════════════════════════════════════════════════════════════════════════════
# §5  SURFACE-SAFE REGION LISTS
# ═══════════════════════════════════════════════════════════════════════════════
# Regions that have a genuine pial cortical surface representation.
# Deep structures (hippocampus CA/DG, amygdala, OB, claustrum, white matter,
# ventricles) are excluded — wb_command surface-vertex-areas would return
# zero or garbage for them.
#
# Exclusion rationale per region:
#   Archicortex        → CA/DG are rolled deep into hippocampal sulcus
#   Cortical_subplate  → amygdala + claustrum, entirely subcortical
#   Hippocampal_form   → CA fields, no pial sheet
#   Subiculum          → deep to hippocampal fissure
#   Amygdala           → deep nuclei
#   Claustrum          → deep to insular cortex
#   Olfactory_bulb     → a bulb structure, not a cortical sheet
#   OB (lvl4)          → same
#   CA, DG, BA_28      → hippocampus proper / entorhinal (deep collateral sulcus)
#   Bed_n_ST           → deep amygdala / BNST
#
# Note on borderline cases:
#   Periarchicortex / Periarchicortex_L3 → parahippocampal gyrus DOES have
#     a gyral surface; included but may be missing in some species.
#   Paleocortex → piriform cortex has a thin surface in many species; included.
#   Olfactory_Cx (lvl3) → piriform + anterior olfactory; has surface.
#   BA_27, Olfactory_tubercle → medial temporal surface; included for primates.

# ── Level 1 ───────────────────────────────────────────────────────────────────
level1_surface = [
    "Isocortex",
    # Allocortex and Periallocortex excluded at lvl1:
    # they aggregate deep structures (hippocampus, amygdala) alongside
    # surface-mapped ones; use lvl2/3 for finer control.
]

# ── Level 2 ───────────────────────────────────────────────────────────────────
level2_surface = [
    "Frontal_lobe",
    "Parietal_lobe",
    "Temporal_lobe",
    "Occipital_lobe",
    "Insular_lobe",
    "Paleocortex",        # piriform / olfactory cortex — has surface
    "Periarchicortex",    # parahippocampal gyrus — has gyral surface
    # EXCLUDED: "Archicortex"       (CA/DG — deep)
    # EXCLUDED: "Cortical_subplate" (amygdala + claustrum — deep)
]

# ── Level 3 ───────────────────────────────────────────────────────────────────
level3_surface = [
    # Frontal
    "oFC",
    "PFC",
    "Motor_premotor",
    # Parietal
    "Somatosensory_Cx",
    "Post_parietal_Cx",
    "PMC",
    # Occipital
    "Striate_Cx",
    "Extra_striate_Cx",
    # Temporal
    "Auditory_Cx",
    "MIPT",
    "Vent_temporal",
    # Insular
    "Insula",
    # Paleocortex
    "Olfactory_Cx",       # piriform / anterior olfactory — surface-mapped
    # Periarchicortex
    "Periarchicortex_L3", # parahippocampal / entorhinal gyrus — has surface
    # EXCLUDED: "Olfactory_bulb"    (not a cortical sheet)
    # EXCLUDED: "Hippocampal_form"  (CA/DG — deep)
    # EXCLUDED: "Subiculum"         (deep to hippocampal fissure)
    # EXCLUDED: "Amygdala"          (subcortical nuclei)
    # EXCLUDED: "Claustrum"         (deep to insula)
]

# ── Level 4 ───────────────────────────────────────────────────────────────────
level4_surface = [
    # oFC
    "Gustatory_Cx",
    "oProiso",
    # PFC
    "BA_10", "BA_11", "BA_13", "BA_12_47", "IFJ",
    "BA_9",  "BA_46", "BA_8",  "BA_55",
    "BA_32", "BA_14", "BA_24", "BA_25", "BA_33",
    "BA_47", "BA_45", "BA_44",
    # Motor / premotor
    "BA_6", "BA_4_M1",
    # Somatosensory
    "BA_1_3_S1", "BA_43", "S2",
    # Post-parietal
    "BA_5", "BA_40", "BA_7", "BA_39",
    # PMC (posteromedial cortex / precuneus / RSC)
    "PGM", "BA_23", "BA_7m", "BA_31",
    # Striate / extra-striate
    "BA_17_V1",
    "BA_18_V2", "BA_19",
    # Auditory
    "BA_41_42_A1", "Belt_A2", "Parabelt_A3", "RTp", "BA_22",
    # MIPT / inferotemporal
    "BA_21", "BA_20", "BA_38",
    # Ventral temporal
    "Parahippocampal", "BA_37",
    # Insula
    "Ri", "Insular_Cx", "BA_52",
    # Olfactory / medial temporal surface
    "BA_27", "Olfactory_tubercle",
    # Periarchicortex — parahippocampal / entorhinal gyrus surface
    "BA_35_36", "Presubiculum", "BA_26_29_30",
    # EXCLUDED: "OB"        (olfactory bulb — not a sheet)
    # EXCLUDED: "CA"        (hippocampus CA — deep)
    # EXCLUDED: "DG"        (dentate gyrus — deep)
    # EXCLUDED: "BA_28"     (entorhinal — deep collateral sulcus, unreliable)
    # EXCLUDED: "Bed_n_ST"  (deep amygdala / BNST)
]

# Flat set for fast membership test
_SURFACE_SAFE = (
    set(level1_surface) | set(level2_surface) |
    set(level3_surface) | set(level4_surface)
)

# Regions excluded from surface extraction (documented for reference)
SURFACE_EXCLUDED = {
    # lvl1
    "Allocortex":        "aggregates deep structures at lvl1",
    "Periallocortex":    "aggregates deep structures at lvl1",
    # lvl2
    "Archicortex":       "CA/DG — deep, no pial sheet",
    "Cortical_subplate": "amygdala + claustrum — subcortical nuclei",
    # lvl3
    "Olfactory_bulb":    "bulb structure, not a cortical sheet",
    "Hippocampal_form":  "CA fields — deep to hippocampal sulcus",
    "Subiculum":         "deep to hippocampal fissure",
    "Amygdala":          "subcortical nuclei",
    "Claustrum":         "deep to insular cortex",
    # lvl4
    "OB":       "olfactory bulb — not a cortical sheet",
    "CA":       "hippocampus CA — deep",
    "DG":       "dentate gyrus — deep",
    "BA_28":    "entorhinal cortex — deep collateral sulcus, unreliable",
    "Bed_n_ST": "deep amygdala / bed nucleus stria terminalis",
}


def get_surface_regions(atlas_level: int) -> list:
    """
    Return surface-safe region names at a given atlas level.
    These are the only regions that should be passed to Stage 0
    surface/thickness extraction.

    Parameters
    ----------
    atlas_level : 1–4

    Returns
    -------
    list of str
    """
    _map = {
        1: level1_surface,
        2: level2_surface,
        3: level3_surface,
        4: level4_surface,
    }
    if atlas_level not in _map:
        raise ValueError(f"atlas_level must be 1–4, got {atlas_level}")
    return list(_map[atlas_level])


def is_surface_safe(region: str) -> bool:
    """Return True if the region has a pial cortical surface representation."""
    return region in _SURFACE_SAFE


def get_surface_region_to_level() -> dict:
    """
    Like get_region_to_level() but restricted to surface-safe regions only.
    Drop-in for atlas_config['region_to_level'] when running Stage 0
    surface/thickness extraction.
    """
    rtl = {}
    for lvl in range(1, 5):
        for r in get_surface_regions(lvl):
            if r not in rtl:
                rtl[r] = lvl
    return rtl


def build_surface_atlas_config() -> dict:
    """
    Build atlas_config restricted to surface-safe regions.
    Use this instead of build_atlas_config_from_py() for Stage 0
    surface/thickness extraction to avoid wasted wb_command calls
    on deep structures that have no vertex coverage.

    Usage:
        from ednix_atlas_regions import build_surface_atlas_config
        atlas_config = build_surface_atlas_config()
        # pass to _process_one_subject via run_wb_extraction
    """
    rtl = get_surface_region_to_level()
    legend = {f"level{i}": get_surface_regions(i) for i in range(1, 5)}
    legend["region_to_level"] = rtl
    return {
        "label_df":        None,
        "id_to_info":      {},
        "name_to_ids":     {},
        "legend":          legend,
        "region_to_level": rtl,
    }


# ── Quick self-test ────────────────────────────────────────────────────────────
if __name__ == "__main__":
    rtl = get_region_to_level()
    print(f"Total regions: {len(rtl)}")
    from collections import Counter
    dist = Counter(rtl.values())
    for lvl in sorted(dist):
        print(f"  level {lvl}: {dist[lvl]} regions")

    print("\nChildren of 'Frontal_lobe' at lvl3:")
    print(" ", get_children("Frontal_lobe", atlas_level=3))

    print("\nChildren of 'PFC' at lvl4:")
    print(" ", get_children("PFC", atlas_level=4))

    print("\nParent chain for 'BA_4_M1':")
    r = "BA_4_M1"
    chain = [r]
    while True:
        p = get_parent(r)
        if p is None: break
        chain.append(p); r = p
    print(" ", " → ".join(chain))

    print("\nFunctional group for 'BA_1_3_S1':", get_functional_group("BA_1_3_S1"))
    print("Functional group for 'BA_10':", get_functional_group("BA_10"))
    print("Integration level for 'BA_17_V1':", get_integration_level("BA_17_V1"))
    print("Integration level for 'BA_10':", get_integration_level("BA_10"))

    print("\nbuild_atlas_config_from_py()['region_to_level'] sample:")
    cfg = build_atlas_config_from_py()
    sample = list(cfg["region_to_level"].items())[:8]
    for k, v in sample:
        print(f"  {k}: lvl{v}")

    print("\n=== Surface-safe regions ===")
    srtl = get_surface_region_to_level()
    from collections import Counter as _Counter
    sdist = _Counter(srtl.values())
    for lvl in sorted(sdist):
        regions = get_surface_regions(lvl)
        print(f"  level {lvl}: {sdist[lvl]} regions")

    print(f"\nTotal surface-safe: {len(srtl)}")
    print(f"Total excluded:     {len(SURFACE_EXCLUDED)}")
    print("\nExcluded regions:")
    for r, reason in SURFACE_EXCLUDED.items():
        print(f"  {r:25s} — {reason}")

    print("\nis_surface_safe checks:")
    for r in ["Isocortex", "CA", "Hippocampal_form", "BA_17_V1", "Amygdala", "Insula"]:
        print(f"  {r:25s} → {is_surface_safe(r)}")

# ═══════════════════════════════════════════════════════════════════════════════
# §6  FLAT-ATLAS SPECIES — DLABEL NAME ALIASES
# ═══════════════════════════════════════════════════════════════════════════════
# For species whose EDNIxCSC dlabel uses long descriptive region names
# (Rat, Dog, Mouse, Mouselemur, Marmoset, Bat) rather than short codes,
# all 4 dlabel maps are IDENTICAL — there is no hierarchy encoded in map number.
# This alias table maps each dlabel base name (after stripping l_/r_ prefix)
# to the corresponding short region name used in the rest of the EDNiX pipeline.
#
# Usage: pass to _build_dlabel_name_map via atlas_config["dlabel_aliases"]
# Effect: dlabel regions are renamed on-the-fly before level filtering,
#         so level_regions matching works correctly for these species.

DLABEL_NAME_ALIASES = {
    # ── Level 1 ─────────────────────────────────────────────────────────────
    "Isocortex":                    "Isocortex",
    "Allocortex":                   "Allocortex",
    "Periallocortex":               "Periallocortex",
    "Subcortical_areas":            "Subcortical_areas",
    "Cortical_White_matter":        "Cortical_White_matter",
    "Diencephalon":                 "Diencephalon",
    # ── Level 2 ─────────────────────────────────────────────────────────────
    "Frontal_lobe":                 "Frontal_lobe",
    "Insular_lobe":                 "Insular_lobe",
    "Temporal_lobe":                "Temporal_lobe",
    "Parietal_lobe":                "Parietal_lobe",
    "Occipital_lobe":               "Occipital_lobe",
    "Paleocortex":                  "Paleocortex",
    "Periarchicortex":              "Periarchicortex",
    "Archicortex":                  "Archicortex",
    "Cortical_subplate":            "Cortical_subplate",
    "Cerebral_nuclei":              "Cerebral_nuclei",
    "Brainstem":                    "Brainstem",
    "Cerebellum":                   "Cerebellum",
    "Cerebellum_White":             "Cerebellum_White",
    # ── Level 3 ─────────────────────────────────────────────────────────────
    "Orbital_frontal_cortex_(oFC)":                              "oFC",
    "Orbital_PFC_(oPFC)":                                        "oFC",
    "Prefrontal_cortex_(PFC)":                                   "PFC",
    "Medial_prefrontal_cortex_(mPFC)":                           "PFC",
    "Dorsolateral_prefrontal_cortex_(dlPFC)":                    "PFC",
    "Ventrolateral_prefrontal_cortex":                           "PFC",
    "Motor_and_premotor":                                        "Motor_premotor",
    "Somatosensory_cortex":                                      "Somatosensory_Cx",
    "Posterior_parietal_cortex":                                 "Post_parietal_Cx",
    "Posterior_medial_cortex_(PMC)":                             "PMC",
    "Visual_striate_cortex":                                     "Striate_Cx",
    "Visual_pre_and_extra_striate_cortex":                       "Extra_striate_Cx",
    "Auditory_cortex_(Superior_temporal)":                       "Auditory_Cx",
    "Middle_Temporal_Inferior_temporal_Temporal_pole_(MIPT)":    "MIPT",
    "Temporal_association_area":                                 "MIPT",
    "Ventral_areas_of_the_temporal_lobe_(vent_Temp)":            "Vent_temporal",
    "Insula_and_others_in_lateral_sulcus":                       "Insula",
    "Olfactory_cortex":                                          "Olfactory_Cx",
    "Olfactory_bulb":                                            "Olfactory_bulb",
    "Olfactory_bulb_(OB)":                                       "Olfactory_bulb",
    "Hippocampal_formation":                                     "Hippocampal_form",
    "Subiculum":                                                 "Subiculum",
    "Parahippocampal_complex":                                   "Periarchicortex_L3",
    "Amygdala":                                                  "Amygdala",
    "Claustrum":                                                 "Claustrum",
    "Septum":                                                    "Septum",
    "Striatum":                                                  "Striatum",
    "Basal_forebrain":                                           "Basal_forebrain",
    "Pallidum":                                                  "Pallidum",
    "Hypothalamus":                                              "Hypothalamus",
    "Thalamus":                                                  "Thalamus",
    "Epithalamus":                                               "Epithalamus",
    "Subthalamus":                                               "Subthalamus",
    "Tectum":                                                    "Tectum",
    "Tegmentum":                                                 "Tegmentum",
    "Midbrain":                                                  "Tegmentum",
    "Cerebral_Peduncles":                                        "Cerebral_ped",
    "Pons":                                                      "Pons",
    "Medulla":                                                   "Medulla",
    # ── Level 4 cortical ────────────────────────────────────────────────────
    "Gustatory_cortex":                                          "Gustatory_Cx",
    "Orbital_proisocortex_and_preiallocortex_(oProiso/preialloC)":"oProiso",
    "BA_10_(Orbital_cortex)":                                    "BA_10",
    "BA_11":                                                     "BA_11",
    "BA_13":                                                     "BA_13",
    "BA_12(/47)":                                                "BA_12_47",
    "IFJ":                                                       "IFJ",
    "BA_9":                                                      "BA_9",
    "BA_46":                                                     "BA_46",
    "BA_8":                                                      "BA_8",
    "BA_55":                                                     "BA_55",
    "BA_32":                                                     "BA_32",
    "BA_14":                                                     "BA_14",
    "BA_24":                                                     "BA_24",
    "BA_25":                                                     "BA_25",
    "BA_33":                                                     "BA_33",
    "BA_6_(supplementary_premotor_area)":                        "BA_6",
    "BA_4_(primary_motor_cortex_M1)":                            "BA_4_M1",
    "BA_47":                                                     "BA_47",
    "BA_45":                                                     "BA_45",
    "BA_44":                                                     "BA_44",
    "BA_1_3_(primary_somatosensory_cortex_S1)":                 "BA_1_3_S1",
    "BA_43":                                                     "BA_43",
    "Secondary_somatosensory_cortex_(S2)":                       "S2",
    "BA_5":                                                      "BA_5",
    "BA_40":                                                     "BA_40",
    "BA_7(a,b,ab)":                                             "BA_7",
    "BA_7(m)":                                                   "BA_7m",
    "BA_39_(angular_gyrus)":                                     "BA_39",
    "Prostriate_area_Parietal_area_medial_(PGM)":                "PGM",
    "BA_23":                                                     "BA_23",
    "BA_31":                                                     "BA_31",
    "BA_17_(V1_cuneus_striate_cortex)":                          "BA_17_V1",
    "BA_18_(V2)":                                                "BA_18_V2",
    "BA_19_(V3_V4_V6)":                                          "BA_19",
    "BA_41/42_(primary_auditory_A1_PEG)":                        "BA_41_42_A1",
    "Ectosylvian_area_(AEG_PEG)":                                "BA_41_42_A1",
    "Belt_(secondary_auditory_A2)":                              "Belt_A2",
    "Parabelt_(tertiary_auditory_cortex_A3)":                    "Parabelt_A3",
    "Polar_rostrotemporal_cortex_(RTp)":                         "RTp",
    "BA_22":                                                     "BA_22",
    "BA_21_(STS_TPJ_MT/V5)":                                    "BA_21",
    "BA_20_(ITS_TE)":                                            "BA_20",
    "BA_38_(TG)":                                                "BA_38",
    "BA_37_(Fusiform)":                                          "BA_37",
    "Temporoparietal_transitional_area_Retroinsula_(Ri)":        "Ri",
    "Insular_cortex":                                            "Insular_Cx",
    "BA_52_(Parainsular_cortex)":                                "BA_52",
    "BA_27_(Pre_)Piriform/amygdaloid_cortex":                    "BA_27",
    "Olfactory_tubercle":                                        "Olfactory_tubercle",
    "Cornu_ammonis_(CA)":                                        "CA",
    "Dentate_gyrus_(DG)":                                        "DG",
    "BA_28_(entorhinal)":                                        "BA_28",
    "BA_35/36_(Perihinal_ectorhinal_cortex)":                    "BA_35_36",
    "Presubiculum":                                              "Presubiculum",
    "BA_26/29/30_(retrosplenial)":                               "BA_26_29_30",
    "Bed_Nucleus_of_the_Stria_Terminalis_(BedNofST)":            "Bed_n_ST",
    # ── Level 4 subcortical ──────────────────────────────────────────────────
    "Caudate":                                                   "Caudate",
    "Putamen":                                                   "Putamen",
    "Nucleus_accumbens":                                         "N_accumbens",
    "Globus_pallidus":                                           "Globus_pallidus",
    "Pulvinar":                                                  "Pulvinar",
    "Geniculate_nuclei":                                         "Geniculate_n",
    "Thalamic_nuclei":                                           "Thalamic_nuclei",
    "Superior_colliculus":                                       "Sup_colliculus",
    "Inferior_colliculus":                                       "Inf_colliculus",
    "Periaqueductal_gray":                                       "PAG",
    "Red_nucleus":                                               "Red_n",
    "Substantia_nigra":                                          "Subst_nigra",
    "Crus_cerebri":                                              "Crus_cerebri",
    # ── Intentionally excluded (not in hierarchy / non-surface) ─────────────
    # 3rd_ventricules, 4th_ventricules, CSF, Lateral_ventricules,
    # Pineal_gland, Habenula, Diagonal_Band_of_Broca, Medial_Septal_Nucleus,
    # Nucleus_Basalis_of_Meynert, Substantia_Innominata, unknown_*
}

# Pre-compute normalized alias lookup for fast runtime use
import re as _re
def _norm(s):
    s = _re.sub(r"^[lLrR]_", "", str(s))
    s = s.lower()
    s = _re.sub(r"[\s_]+", "_", s)
    s = _re.sub(r"[(),.\-/\\]", "", s)
    return _re.sub(r"_+", "_", s).strip("_")

DLABEL_NAME_ALIASES_NORM = {_norm(k): v for k, v in DLABEL_NAME_ALIASES.items()}


def build_flat_atlas_config(surface_only: bool = True) -> dict:
    """
    Build atlas_config for flat-atlas species (Rat, Dog, Mouse, Mouselemur,
    Marmoset, Bat) whose EDNIxCSC dlabel has identical content in all 4 maps
    and uses long descriptive region names.

    This config:
    1. Provides DLABEL_NAME_ALIASES_NORM so _build_dlabel_name_map can rename
       dlabel regions to their short codes on-the-fly.
    2. Provides region_to_level with the correct hierarchy assignments so
       level-filtered extraction works correctly.
    3. Since all 4 maps are identical, extraction at each level pulls only the
       regions assigned to that level — no duplication across levels.

    Parameters
    ----------
    surface_only : if True (default), uses surface-safe regions only.

    Usage in run_wb_extraction:
        from ednix_atlas_regions import build_flat_atlas_config
        run_wb_extraction(..., atlas_config_override=build_flat_atlas_config())
    """
    if surface_only:
        rtl = get_surface_region_to_level()
    else:
        rtl = get_region_to_level()

    legend = {f"level{i}": (get_surface_regions(i) if surface_only
                             else get_level_regions(i))
              for i in range(1, 5)}
    legend["region_to_level"] = rtl

    return {
        "label_df":           None,
        "id_to_info":         {},
        "name_to_ids":        {},
        "legend":             legend,
        "region_to_level":    rtl,
        "dlabel_aliases":     DLABEL_NAME_ALIASES_NORM,  # ← key addition
        "is_flat_atlas":      True,
    }
