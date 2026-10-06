"""
ednix_checkpoint.py
===================
Checkpoint detection for both the EDNiX fMRI and anatomical pipelines.

Takes the user's existing Skip_step and list_to_remove as inputs — no new
configuration needed. Detects what is already done per subject, combines
with the user's intended skips, groups subjects by their resume point, and
returns (list_to_keep, Skip_step) pairs ready to pass directly to the
existing launcher functions unchanged.

────────────────────────────────────────────────────────────────────────────
How it works
────────────────────────────────────────────────────────────────────────────

User already expresses intent in the launcher:

    Skip_step    = [10, 12, 13, 14, 'itk_1', 'itk_2', 'itk_3', 'Clean']
    list_to_remove = [('Trinity', '6'), ('Trinity', '3')]

The checkpoint function:
  1. Respects list_to_remove  — those subjects are never included
  2. Respects Skip_step       — those steps are never run for anyone
  3. Detects completed steps  — already-done steps are added to each
                                subject's skip list
  4. Groups subjects          — subjects at the same resume point share
                                a Skip_step and are launched together

────────────────────────────────────────────────────────────────────────────
fMRI usage
────────────────────────────────────────────────────────────────────────────
    from ednix_checkpoint import fmri_resume_groups

    Skip_step      = [10, 12, 'itk_1', 'itk_2', 'Clean']
    list_to_remove = [('Trinity', '6'), ('Trinity', '3')]

    groups = fmri_resume_groups(
        allinfo_study_c, bids_dir,
        reference      = 'EDNiX',
        resting_or_task= 'resting',
        type_norm      = 'T1w',
        ICA_cleaning   = 'Skip',
        endfmri        = '*_task-rest_*.nii.gz',
        Skip_step      = Skip_step,
        list_to_remove = list_to_remove)

    for list_to_keep, Skip_step_group in groups:
        modalities.fMRI._0_Pipeline_launcher.preprocess_data(
            Skip_step_group, MAIN_PATH, bids_dir, ...,
            list_to_keep=list_to_keep,
            list_to_remove=list_to_remove, ...)

────────────────────────────────────────────────────────────────────────────
Anat usage
────────────────────────────────────────────────────────────────────────────
    from ednix_checkpoint import anat_resume_groups

    Skip_step      = [10, 11, 14, 15, 'flat_map', 'itk_3']
    list_to_remove = [('Trinity', '6')]

    groups = anat_resume_groups(
        allinfo_study_c, bids_dir,
        reference      = 'EDNiX',
        type_norm      = 'T1w',
        species        = 'Dog',
        Skip_step      = Skip_step,
        list_to_remove = list_to_remove)

    for list_to_keep, Skip_step_group in groups:
        anat._0_Pipeline_launcher.preprocess_anat(
            Skip_step_group, MAIN_PATH, bids_dir, ...,
            list_to_keep=list_to_keep,
            list_to_remove=list_to_remove, ...)
"""

import os
import glob

opj = os.path.join
ope = os.path.exists
opi = os.path.isfile
opb = os.path.basename


# =============================================================================
# fMRI step definitions
# =============================================================================

def _all_runs(rs_roots, directory, pattern):
    """True if output pattern exists for every run."""
    return all(opi(opj(directory, r + pattern)) for r in rs_roots)


FMRI_STEPS = [
    {
        'step_id': 1,
        'name':    'Step 1 — fMRI preprocessing (motion correction, STC, N4)',
        'check':   lambda p, rs, ica, tn: _all_runs(
            rs, p['dir_prepro_raw_process'],
            '_space-func_desc-runMean_n4Bias.nii.gz'),
    },
    {
        'step_id': 2,
        'name':    'Step 2 — Coregistration to norm (distortion correction)',
        'check':   lambda p, rs, ica, tn: _all_runs(
            rs, p['dir_prepro_raw_process'],
            '_space-func_desc-fMRI_run_inRef.nii.gz'),
    },
    {
        'step_id': 3,
        'name':    'Step 3 — Functional brain masking',
        # Output: {ID}_fMRI_mask.nii.gz — named after subject ID, not run root.
        # ID is not directly available here so we check for any *_fMRI_mask.nii.gz
        # in the masks directory (there is exactly one per subject/session).
        'check':   lambda p, rs, ica, tn: bool(
            glob.glob(opj(p['dir_prepro_raw_masks'], '*_fMRI_mask.nii.gz'))),
    },
    {
        'step_id': 'itk_1',
        'name':    'itk_1 — Manual mask QC step 1 (ITK-SNAP)',
        'check':   lambda p, rs, ica, tn: bool(
            glob.glob(opj(p['dir_prepro_raw_masks'], '*_fMRI_mask.nii.gz'))),
    },
    {
        'step_id': 5,
        'name':    'Step 5 — Anat-to-func registration',
        'check':   lambda p, rs, ica, tn: _all_runs(
            rs, p['dir_prepro_orig_process'],
            '_space-acpc-func_desc-fMRI_run_inRef.nii.gz'),
    },
    {
        'step_id': 'itk_2',
        'name':    'itk_2 — Anat-func overlap QC (ITK-SNAP)',
        'check':   lambda p, rs, ica, tn: _all_runs(
            rs, p['dir_prepro_orig_process'],
            '_space-acpc-func_desc-fMRI_run_inRef.nii.gz'),
    },
    {
        'step_id': 6,
        'name':    'Step 6 — ICA / Melodic',
        'check':   lambda p, rs, ica, tn: (
            ica == 'Skip' or
            ope(opj(p['dir_RS_ICA_native'], 'melodic_IC.nii.gz'))),
    },
    {
        'step_id': 7,
        'name':    'Step 7 — Signal regression',
        # acpc-func/postprocessed_rs/{root}_space-acpc-func_desc-fMRI_residual.nii.gz
        'check':   lambda p, rs, ica, tn: _all_runs(
            rs, p['dir_prepro_orig_postprocessed'],
            '_space-acpc-func_desc-fMRI_residual.nii.gz'),
    },
    {
        'step_id': 8,
        'name':    'Step 8 — fMRI to anat space',
        # acpc-anat/postprocessed_rs/{root}_space-acpc-anat_desc-fMRI_residual.nii.gz
        'check':   lambda p, rs, ica, tn: _all_runs(
            rs, p['dir_prepro_acpc_postprocessed'],
            '_space-acpc-anat_desc-fMRI_residual.nii.gz'),
    },
    {
        'step_id': 9,
        'name':    'Step 9 — Coregistration to template space',
        # templates/EDNiX/postprocessed_rs/{root}_space-template_desc-fMRI_residual.nii.gz
        'check':   lambda p, rs, ica, tn: _all_runs(
            rs, p['dir_prepro_template_postprocessed'],
            '_space-template_desc-fMRI_residual.nii.gz'),
    },
    {
        'step_id': 10,
        'name':    'Step 10 — Functional connectivity matrix',
        # acpc-func/Stats/Correl_matrix/{atlas}/{atlas}_*_correlation_matrix.csv
        'check':   lambda p, rs, ica, tn: bool(
            glob.glob(opj(p['path_func'], 'acpc-func', 'Stats',
                          'Correl_matrix', '*', '*', '*_correlation_matrix.csv'))),
    },
    {
        'step_id': 11,
        'name':    'Step 11 — Seed-based analysis',
        # templates/EDNiX/Stats/SBA/ contains region subdirectories when done
        'check':   lambda p, rs, ica, tn: bool(
            glob.glob(opj(p['dir_prepro_template'], 'Stats', 'SBA', '*', '*.nii.gz'))),
    },
    {
        'step_id': 12,
        'name':    'Step 12 — fMRI QC report',
        # func/QC/{root}_QC_summary.png
        'check':   lambda p, rs, ica, tn: (
            opi(opj(p['path_func'], 'QC',
                    rs[0] + '_QC_summary.png')) if rs else False),
    },
    {
        'step_id': 'itk_3',
        'name':    'itk_3 — Not applicable for fMRI (anat-only step)',
        # itk_3 is an anat-only step. Recognized here so it is preserved in
        # Skip_step when the user passes it, instead of being silently dropped.
        'check':   lambda p, rs, ica, tn: True,
    },
    {
        'step_id': 'Clean',
        'name':    'Clean — Intermediate file cleanup',
        # Clean is a post-processing step — always mark as complete since
        # it is safe to skip (cleaning up intermediate files is optional).
        # If the user has it in Skip_step it will be caught there first.
        'check':   lambda p, rs, ica, tn: True,
    },
]


# =============================================================================
# Anat step definitions
# =============================================================================

ANAT_STEPS = [
    # ---- Loop 1 -------------------------------------------------------------
    {
        'step_id': 1,
        'loop':    1,
        'name':    'Step 1 — Orientation correction + N4/N3 bias field',
        # Try all known desc- variants and both T1w/T2w
        # 'Bias' covers _space-raw_desc-Bias_ seen in some species (e.g. Bat)
        'check':   lambda p, ID, tn, sp: any(
            opi(opj(p['dir_prepro'],
                    ID + '_space-raw_desc-' + d + '_' + t + '.nii.gz'))
            for d in ['n4Bias', 'n3Bias', 'reorient', 'Bias']
            for t in ['T1w', 'T2w']),
    },
    {
        'step_id': 2,
        'loop':    1,
        'name':    'Step 2 — Skull-strip step 1 + ACPC alignment',
        # SS-step1 is deleted by Clean. Use acpc_0GenericAffine.mat as the
        # primary proxy (always survives) combined with any SS variant present.
        # acpc mat alone is sufficient since it's only produced after ACPC alignment.
        'check':   lambda p, ID, tn, sp: (
            opi(opj(p['dir_transfo'], 'acpc_0GenericAffine.mat')) and
            any(opi(opj(p['volumes_dir'],
                     ID + '_space-acpc_desc-' + desc + '_' + t + '.nii.gz'))
                for desc in ['SS-step1', 'SS-step2', 'SS', 'template']
                for t in [tn, 'T1w', 'T2w'])),
    },
    {
        'step_id': 'itk_1',
        'loop':    1,
        'name':    'itk_1 — Manual mask QC step 1 (ITK-SNAP)',
        # desc-step1_mask is deleted by Clean — use acpc mat as proxy
        # (only exists after step 2 which follows itk_1)
        'check':   lambda p, ID, tn, sp: (
            opi(opj(p['masks_dir'], ID + '_desc-step1_mask.nii.gz')) or
            opi(opj(p['dir_transfo'], 'acpc_0GenericAffine.mat'))),
    },
    # ---- Study template steps (loop 1 break — only when creat_study_template=True)
    # sp = species string (unused by these checks — bids_dir not available here
    # so we check via path_anat going up to find sty_template)
    {
        'step_id': 3,
        'loop':    1,
        'name':    'Step 3 — Study template creation',
        # Output: bids_dir/sty_template/templates/{reference}/volumes/studyTemplate_space-{reference}.nii.gz
        # sp = (species, bids_dir, reference) when creat_study_template=True
        # sp = species string otherwise — check is safely False (no sty_template)
        'check':   lambda p, ID, tn, sp: (
            opi(opj(sp[1], 'sty_template', 'templates', sp[2], 'volumes',
                    'studyTemplate_space-' + sp[2] + '.nii.gz'))
            if isinstance(sp, tuple) else False),
    },
    {
        'step_id': 4,
        'loop':    1,
        'name':    'Step 4 — Study template skull-stripping',
        'check':   lambda p, ID, tn, sp: (
            opi(opj(sp[1], 'sty_template', 'derivatives', 'acpc', 'volumes',
                    'masks', 'studyTemplate_mask.nii.gz'))
            if isinstance(sp, tuple) else False),
    },
    {
        'step_id': 7,
        'loop':    1,
        'name':    'Step 7 — Study template to atlas registration',
        # Output: sty_template/derivatives/acpc/volumes/labels/studyTemplate_seg-*_dseg.nii.gz
        # Also: sty_template/derivatives/acpc/matrices/studyTemplate_to_{ref}_Final_0GenericAffine.mat
        'check':   lambda p, ID, tn, sp: (
            opi(opj(sp[1], 'sty_template', 'derivatives', 'acpc', 'matrices',
                    'studyTemplate_to_' + sp[2] + '_Final_0GenericAffine.mat')) and
            bool(glob.glob(opj(sp[1], 'sty_template', 'derivatives', 'acpc',
                               'volumes', 'labels', 'studyTemplate_seg-*_dseg.nii.gz')))
            if isinstance(sp, tuple) else False),
    },
    # ---- Loop 2 -------------------------------------------------------------
    {
        'step_id': 5,
        'loop':    2,
        'name':    'Step 5 — Skull-strip step 2 / individual template brain',
        # SS-step2 and desc-step2_mask are deleted by Clean.
        # Use space-acpc_desc-template_ and desc-Bias_ which survive cleaning.
        'check':   lambda p, ID, tn, sp: (
            any(opi(opj(p['volumes_dir'],
                     ID + '_space-acpc_desc-template_' + t + '.nii.gz'))
                for t in [tn, 'T1w', 'T2w']) or
            any(opi(opj(p['volumes_dir'],
                     ID + '_space-acpc_desc-SS-step2_' + t + '.nii.gz'))
                for t in [tn, 'T1w', 'T2w'])),
    },
    {
        'step_id': 'itk_2',
        'loop':    2,
        'name':    'itk_2 — Manual mask QC step 2 (ITK-SNAP)',
        # Same fallback — use template file if SS-step2 was cleaned
        'check':   lambda p, ID, tn, sp: (
            any(opi(opj(p['volumes_dir'],
                    ID + '_space-acpc_desc-template_' + t + '.nii.gz'))
                for t in [tn, 'T1w', 'T2w']) or
            any(opi(opj(p['volumes_dir'],
                    ID + '_space-acpc_desc-SS-step2_' + t + '.nii.gz'))
                for t in [tn, 'T1w', 'T2w'])),
    },
    {
        'step_id': 6,
        'loop':    2,
        'name':    'Step 6 — Registration to reference template',
        'check':   lambda p, ID, tn, sp: any(
            opi(opj(p['dir_transfo'], f))
            for f in [
                'native_to_EDNiX_Final_0GenericAffine.mat',
                'native_to_EDNiX_Final_1Warp.nii.gz',
                'native_to_studyTemplate_Final_0GenericAffine.mat',
                'native_to_refSession_Final_0GenericAffine.mat',
            ]),
    },
    # ---- Loop 3 -------------------------------------------------------------
    {
        'step_id': 8,
        'loop':    3,
        'name':    'Step 8 — Backtonative (labels/masks to ACPC space)',
        'check':   lambda p, ID, tn, sp: (
            any(opi(opj(p['volumes_dir'],
                     ID + '_space-acpc_desc-SS_' + t + '.nii.gz'))
                for t in [tn, 'T1w', 'T2w']) and
            opi(opj(p['masks_dir'], ID + '_space-acpc_mask.nii.gz'))),
    },
    {
        'step_id': 'itk_3',
        'loop':    3,
        'name':    'itk_3 — Final mask QC / backtonative mask (ITK-SNAP)',
        'check':   lambda p, ID, tn, sp: any(
            opi(opj(p['volumes_dir'],
                    ID + '_space-acpc_desc-SS_' + t + '.nii.gz'))
            for t in [tn, 'T1w', 'T2w']),
    },
    {
        'step_id': 10,
        'loop':    3,
        'name':    'Step 10 — FreeSurfer recon-all preparation',
        # Step 9 = prepa_img (always runs) + toFS (only when FS enabled).
        # prepa_img outputs: desc-norm T1w/T2w, seg-wm, seg-filled in labels_dir.
        # toFS outputs: orig.mgz + brain.mgz in FS_dir.
        # If FS_dir does not exist the user intentionally skipped FS (steps 10-16
        # should be in their Skip_step) — prepa_img outputs alone mark step 9 done.
        # prepa_img auto-detects T1w/T2w from the input filename so we check both.
        'check':   lambda p, ID, tn, sp: (
            any(opi(opj(p['dir_prepro'],
                        ID + '_desc-norm_' + t + '.nii.gz'))
                for t in [tn, 'T1w', 'T2w']) and
            opi(opj(p['labels_dir'], ID + '_seg-wm_dseg.nii.gz')) and
            opi(opj(p['labels_dir'], ID + '_seg-filled_dseg.nii.gz')) and
            (not ope(p['FS_dir']) or (
                opi(opj(p['FS_dir'], ID, 'mri', 'orig.mgz')) and
                opi(opj(p['FS_dir'], ID, 'mri', 'brain.mgz'))))),
    },
    {
        'step_id': 11,
        'loop':    3,
        'name':    'Step 11 — White surface construction',
        'check':   lambda p, ID, tn, sp: (
            opi(opj(p['FS_dir'], ID, 'surf', 'lh.white')) and
            opi(opj(p['FS_dir'], ID, 'surf', 'rh.white'))),
    },
    {
        'step_id': 12,
        'loop':    3,
        'name':    'Step 12 — Pial surface construction',
        'check':   lambda p, ID, tn, sp: (
            opi(opj(p['FS_dir'], ID, 'surf', 'lh.pial')) and
            opi(opj(p['FS_dir'], ID, 'surf', 'rh.pial'))),
    },
    {
        'step_id': 'flat_map',
        'loop':    3,
        'name':    'flat_map — Flat map construction',
        'check':   lambda p, ID, tn, sp: (
            opi(opj(p['FS_dir'], ID, 'surf', 'lh.flat')) and
            opi(opj(p['FS_dir'], ID, 'surf', 'rh.flat'))),
    },
    {
        'step_id': 14,
        'loop':    3,
        'name':    'Step 14 — FreeSurfer finalise + atlas label mapping',
        'check':   lambda p, ID, tn, sp: (
            opi(opj(p['FS_dir'], ID, 'surf', 'lh.sphere.reg')) and
            opi(opj(p['FS_dir'], ID, 'surf', 'rh.sphere.reg'))),
    },
    {
        'step_id': 15,
        'loop':    3,
        'name':    'Step 15 — Connectome Workbench conversion',
        # Native_resol lives inside surfaces/ subdirectory of dir_native
        'check':   lambda p, ID, tn, sp: opi(opj(
            p['dir_native'], 'surfaces', 'Native_resol',
            ID + '_native_LR.wb.spec')),
    },
    {
        'step_id': 16,
        'loop':    3,
        'name':    'Step 16 — Anatomical QC SNR report',
        # Output: {path_anat}/QC_anat/{type_norm}_QC_summary.png
        # Written by _16_anat_QC_SNR.create_qc_figure() for all species.
        # Check both tn and T1w/T2w variants since type may differ per subject.
        'check':   lambda p, ID, tn, sp: any(
            opi(opj(p['path_anat'], 'QC_anat', t + '_QC_summary.png'))
            for t in [tn, 'T1w', 'T2w']),
    },
    {
        'step_id': 'Clean',
        'loop':    3,
        'name':    'Clean — Intermediate file cleanup',
        # Always mark as complete — safe to skip, cleanup is optional.
        'check':   lambda p, ID, tn, sp: True,
    },
]


# =============================================================================
# Path builders
# =============================================================================

def _fmri_paths(data_path, reference, resting_or_task):
    from Tools import getpath
    # getpath.func returns 24 values — unpack all explicitly
    (dir_prepro_raw, dir_prepro_raw_process, dir_prepro_raw_masks,
     dir_prepro_raw_matrices, path_func, dir_fmap, dir_prepro_fmap,
     dir_prepro_orig, dir_prepro_orig_labels, dir_prepro_orig_masks,
     dir_prepro_orig_process, dir_prepro_orig_postprocessed,
     dir_prepro_acpc, dir_prepro_acpc_labels, dir_prepro_acpc_masks,
     dir_prepro_acpc_process, dir_prepro_acpc_postprocessed,
     dir_prepro_template, dir_prepro_template_labels, dir_prepro_template_masks,
     dir_prepro_template_process, dir_prepro_template_postprocessed,
     dir_prepro_acpc_matrices, dir_prepro_orig_matrices
     ) = getpath.func(data_path, reference, resting_or_task)
    # All paths come directly from getpath.func() — no derivation needed
    # dir_prepro_orig_postprocessed = acpc-func/postprocessed_rs  (step 7)
    # dir_prepro_acpc_postprocessed = acpc-anat/postprocessed_rs  (step 8)
    # dir_prepro_template_postprocessed = templates/EDNiX/postprocessed_rs (step 9)
    # path_func = func/  (for QC and Stats)
    return {
        'dir_prepro_raw_process':            dir_prepro_raw_process,
        'dir_prepro_raw_masks':              dir_prepro_raw_masks,
        'dir_prepro_raw_matrices':           dir_prepro_raw_matrices,
        'dir_prepro_orig_process':           dir_prepro_orig_process,
        'dir_prepro_orig_postprocessed':     dir_prepro_orig_postprocessed,
        'dir_prepro_acpc_process':           dir_prepro_acpc_process,
        'dir_prepro_acpc_masks':             dir_prepro_acpc_masks,
        'dir_prepro_acpc_postprocessed':     dir_prepro_acpc_postprocessed,
        'dir_prepro_template_process':       dir_prepro_template_process,
        'dir_prepro_template_postprocessed': dir_prepro_template_postprocessed,
        'dir_prepro_template':               dir_prepro_template,
        'path_func':                         path_func,
        'dir_RS_ICA_native':                 opj(dir_prepro_orig_process, 'Melodic'),
    }


def _anat_paths(data_path, reference):
    from Tools import getpath
    (path_anat, dir_transfo, FS_dir, dir_prepro,
     dir_native, volumes_dir, labels_dir, masks_dir,
     *_) = getpath.anat(data_path, reference, '', False, False, 'native')
    return {
        'path_anat':   path_anat,
        'dir_prepro':  dir_prepro,
        'dir_transfo': dir_transfo,
        'FS_dir':      FS_dir,
        'dir_native':  dir_native,
        'volumes_dir': volumes_dir,
        'labels_dir':  labels_dir,
        'masks_dir':   masks_dir,
    }


def _rs_roots(data_path, endfmri):
    from modalities.fMRI.extract_filename import extract_filename
    files = sorted(glob.glob(opj(data_path, 'func', endfmri)))
    return [extract_filename(f) for f in files]


def _is_removed(ID, Session, list_to_remove):
    """Return True if (ID, Session) is in list_to_remove."""
    return any(
        str(ID) == str(r[0]) and str(Session) == str(r[1])
        for r in list_to_remove)


# =============================================================================
# Core detection
# =============================================================================

def _detect(all_ID, all_Session, all_data_path,
            steps, path_fn, check_fn,
            user_Skip_step, list_to_remove):
    """
    Core detection loop shared by fMRI and anat.

    For each subject:
      effective_skip = user_Skip_step (never run for anyone)
                     + completed_steps (already done for this subject)

    A step needs to run only if:
      - not in user_Skip_step, AND
      - not already complete on disk
    """
    status = []

    for ID, Session, data_path in zip(all_ID, all_Session, all_data_path):
        key = f'sub-{ID}_ses-{Session}'

        if _is_removed(ID, Session, list_to_remove):
            status.append({
                'ID': ID, 'Session': Session, 'key': key,
                'removed': True, 'complete': False,
                'steps': {}, 'effective_skip': [],
                'needs_run': [], 'error': None,
            })
            continue

        entry = {
            'ID':             ID,
            'Session':        Session,
            'key':            key,
            'removed':        False,
            'complete':       False,
            'steps':          {},
            'effective_skip': [],
            'needs_run':      [],
            'error':          None,
        }

        try:
            paths = path_fn(data_path)

            for step in steps:
                sid = step['step_id']

                if sid in user_Skip_step:
                    entry['steps'][sid] = 'skipped_by_user'
                    entry['effective_skip'].append(sid)
                    continue

                try:
                    done = check_fn(step, paths, data_path)
                except Exception as e:
                    done = False
                    entry['error'] = f'Check failed step {sid}: {e}'

                entry['steps'][sid] = done
                if done:
                    entry['effective_skip'].append(sid)
                else:
                    entry['needs_run'].append(sid)

            entry['complete'] = len(entry['needs_run']) == 0

        except Exception as e:
            entry['error'] = str(e)

        status.append(entry)

    return status


# =============================================================================
# Grouping
# =============================================================================

def _group_by_resume_point(status, steps):
    """
    Group incomplete subjects by their first step that still needs to run.

    Each group shares a Skip_step = intersection of all subjects' effective
    skips so no subject reruns a step it already completed.
    Ordered least-to-most advanced so dependencies are always met.

    Returns list of (list_to_keep, Skip_step) tuples.
    list_to_keep contains (ID, Session) tuples matching the launcher format.
    """
    all_step_ids = [s['step_id'] for s in steps]
    step_index   = {sid: i for i, sid in enumerate(all_step_ids)}

    groups = {}
    for s in status:
        if s['removed'] or s['complete'] or s['error']:
            continue
        first_needed = s['needs_run'][0] if s['needs_run'] else None
        if first_needed is None:
            continue
        groups.setdefault(first_needed, []).append(s)

    result = []
    for first_step, subjects in sorted(
            groups.items(),
            key=lambda x: step_index.get(x[0], 999)):

        # (ID, Session) tuples — Session kept as original type to match
        # allinfo_study_c and load_data_bids expectations
        list_to_keep = [(s['ID'], s['Session']) for s in subjects]

        # Skip_step = steps every subject in this group has completed or
        # had user-skipped — intersection ensures no subject reruns anything
        all_skips   = [set(s['effective_skip']) for s in subjects]
        shared_skip = list(all_skips[0].intersection(*all_skips[1:])
                           if len(all_skips) > 1 else all_skips[0])

        # Preserve original step order
        Skip_step = [sid for sid in all_step_ids if sid in shared_skip]

        result.append((list_to_keep, Skip_step))

    return result


# =============================================================================
# Report printer
# =============================================================================

def _print_report(status, steps, title, user_Skip_step, list_to_remove):
    step_labels = {s['step_id']: s['name'] for s in steps}
    step_ids    = [s['step_id'] for s in steps]
    loop_map    = {s['step_id']: s.get('loop') for s in steps}

    print(f"\n{'='*80}")
    print(f"EDNiX {title} Checkpoint Report")
    print(f"{'='*80}")
    print(f"  User Skip_step    : {user_Skip_step}")
    print(f"  list_to_remove    : {list_to_remove}")

    n_complete   = sum(1 for s in status if not s['removed'] and s['complete'])
    n_incomplete = sum(1 for s in status
                       if not s['removed'] and not s['complete'] and not s['error'])
    n_removed    = sum(1 for s in status if s['removed'])
    n_errors     = sum(1 for s in status if s['error'])
    print(f"  Complete: {n_complete}   Incomplete: {n_incomplete}   "
          f"Removed: {n_removed}   Errors: {n_errors}\n")

    for s in status:
        label = f"sub-{s['ID']} ses-{s['Session']}"

        if s['removed']:
            print(f"  {label}  ✕  REMOVED (list_to_remove)")
            continue
        if s['error']:
            print(f"  {label}  ⚠  ERROR: {s['error']}")
            continue
        if s['complete']:
            n_user_skip = sum(1 for v in s['steps'].values()
                              if v == 'skipped_by_user')
            print(f"  {label}  ✓  COMPLETE"
                  + (f"  ({n_user_skip} user-skipped)" if n_user_skip else ""))
            continue

        n_done      = sum(1 for v in s['steps'].values() if v is True)
        n_user_skip = sum(1 for v in s['steps'].values()
                          if v == 'skipped_by_user')
        n_todo      = len(s['needs_run'])
        n_tot       = len(s['steps'])

        print(f"  {label}  ✗  INCOMPLETE  "
              f"[{n_done} done / {n_user_skip} user-skipped / "
              f"{n_todo} to run / {n_tot} total]")
        print(f"      Resume from : "
              f"{step_labels.get(s['needs_run'][0], s['needs_run'][0])}")

        current_loop = None
        for sid in step_ids:
            loop = loop_map.get(sid)
            if loop and loop != current_loop:
                current_loop = loop
                print(f"      ── Loop {loop} ──")
            val = s['steps'].get(sid)
            if val is True:
                symbol, note = '✓', ''
            elif val == 'skipped_by_user':
                symbol, note = '—', '  (user Skip_step)'
            elif val is False:
                symbol, note = '✗', ''
            else:
                symbol, note = '?', '  (not checked)'
            marker = '  ← TO RUN' if sid == s['needs_run'][0] else ''
            print(f"        {symbol} {step_labels.get(sid, str(sid))}"
                  f"{marker}{note}")
        print()

    print(f"{'='*80}\n")


def _print_groups(groups, steps, title):
    step_labels  = {s['step_id']: s['name'] for s in steps}
    all_step_ids = [s['step_id'] for s in steps]

    print(f"  Launcher groups ({title}):")
    if not groups:
        print("    All subjects complete or removed — nothing to run.\n")
        return
    for i, (ids, skip) in enumerate(groups):
        needs_run = [s for s in all_step_ids if s not in skip]
        resume    = needs_run[0] if needs_run else 'done'
        print(f"\n  Group {i+1}  ({len(ids)} subject(s))")
        print(f"    list_to_keep = {ids}")
        print(f"    Skip_step    = {skip}")
        print(f"    Will run     : {needs_run}")
        print(f"    Resume from  : {step_labels.get(resume, resume)}")
    print()


# =============================================================================
# Public API — fMRI
# =============================================================================

def fmri_resume_groups(allinfo_study_c, bids_dir,
                       reference, resting_or_task, type_norm,
                       Skip_step=None,
                       list_to_remove=None,
                       ICA_cleaning='Skip',
                       endfmri='*_task-rest_*.nii.gz',
                       verbose=True):
    """
    Detect fMRI pipeline status and return launcher-ready (list_to_keep,
    Skip_step) groups.

    Parameters
    ----------
    allinfo_study_c : pandas.DataFrame
    bids_dir : str
    reference : str
    resting_or_task : str
    type_norm : str
    Skip_step : list — steps never to run (your launcher's Skip_step)
    list_to_remove : list of (ID, Session) tuples
    ICA_cleaning : str
    endfmri : str
    verbose : bool

    Returns
    -------
    groups : list of (list_to_keep, Skip_step) tuples

    Example
    -------
    groups = fmri_resume_groups(
        allinfo_study_c, bids_dir, 'EDNiX', 'resting', 'T1w',
        Skip_step=['itk_1', 'itk_2', 'Clean'],
        list_to_remove=[('Trinity', '6')])

    for list_to_keep, Skip_step_group in groups:
        modalities.fMRI._0_Pipeline_launcher.preprocess_data(
            Skip_step_group, MAIN_PATH, bids_dir, ...,
            list_to_keep=list_to_keep,
            list_to_remove=list_to_remove, ...)
    """
    from Tools import Load_subject_with_BIDS

    if Skip_step      is None: Skip_step      = []
    if list_to_remove is None: list_to_remove = []

    all_ID, all_Session, all_data_path, *_ = \
        Load_subject_with_BIDS.load_data_bids(allinfo_study_c, bids_dir, [], [])

    def path_fn(dp):
        return _fmri_paths(dp, reference, resting_or_task)

    def check_fn(step, paths, dp):
        rs = _rs_roots(dp, endfmri)
        return step['check'](paths, rs, ICA_cleaning, type_norm)

    status = _detect(all_ID, all_Session, all_data_path,
                     FMRI_STEPS, path_fn, check_fn,
                     Skip_step, list_to_remove)

    if verbose:
        _print_report(status, FMRI_STEPS, 'fMRI', Skip_step, list_to_remove)

    groups = _group_by_resume_point(status, FMRI_STEPS)

    if verbose:
        _print_groups(groups, FMRI_STEPS, 'fMRI')

    return groups


# =============================================================================
# Public API — Anat
# =============================================================================

def anat_resume_groups(allinfo_study_c, bids_dir,
                       reference, type_norm, species,
                       Skip_step=None,
                       list_to_remove=None,
                       creat_study_template=False,
                       verbose=True):
    """
    Detect anat pipeline status and return launcher-ready (list_to_keep,
    Skip_step) groups.

    Parameters
    ----------
    allinfo_study_c : pandas.DataFrame
    bids_dir : str
    reference : str
    type_norm : str
    species : str
    Skip_step : list — steps never to run (your launcher's Skip_step)
    list_to_remove : list of (ID, Session) tuples
    verbose : bool

    Returns
    -------
    groups : list of (list_to_keep, Skip_step) tuples

    Example
    -------
    groups = anat_resume_groups(
        allinfo_study_c, bids_dir, 'EDNiX', 'T1w', 'Dog',
        Skip_step=[10, 11, 'flat_map', 'itk_3'],
        list_to_remove=[('Trinity', '6')],
        creat_study_template=False)  # auto-skips steps 3 and 4 when False

    for list_to_keep, Skip_step_group in groups:
        anat._0_Pipeline_launcher.preprocess_anat(
            Skip_step_group, MAIN_PATH, bids_dir, ...,
            list_to_keep=list_to_keep,
            list_to_remove=list_to_remove, ...)
    """
    from Tools import Load_subject_with_BIDS

    if Skip_step      is None: Skip_step      = []
    if list_to_remove is None: list_to_remove = []

    # Steps 3 and 4 are dataset-level (study template), not per-subject.
    # They run once across all subjects before per-subject processing continues.
    #
    # creat_study_template=False: auto-skip 3 and 4 for everyone.
    # creat_study_template=True:  check if the template already exists.
    #   - If it exists: skip 3 and 4 for everyone (already done).
    #   - If it doesn't: keep 3 and 4 as needed — the launcher will run them
    #     for the whole group before continuing per-subject steps.
    if not creat_study_template:
        for sid in [3, 4, 7]:
            if sid not in Skip_step:
                Skip_step = list(Skip_step) + [sid]
    else:
        # Check if study template already exists (dataset-level check)
        sty_template_done = (
            opi(opj(bids_dir, 'sty_template', 'templates', reference,
                    'volumes', 'studyTemplate_space-' + reference + '.nii.gz')) and
            opi(opj(bids_dir, 'sty_template', 'derivatives', 'acpc',
                    'volumes', 'masks', 'studyTemplate_mask.nii.gz')))
        sty_template_step7_done = (
            sty_template_done and
            opi(opj(bids_dir, 'sty_template', 'derivatives', 'acpc', 'matrices',
                    'studyTemplate_to_' + reference + '_Final_0GenericAffine.mat')) and
            bool(glob.glob(opj(bids_dir, 'sty_template', 'derivatives', 'acpc',
                               'volumes', 'labels', 'studyTemplate_seg-*_dseg.nii.gz'))))
        if sty_template_done:
            for sid in [3, 4]:
                if sid not in Skip_step:
                    Skip_step = list(Skip_step) + [sid]
        if sty_template_step7_done:
            if 7 not in Skip_step:
                Skip_step = list(Skip_step) + [7]
        # If not done, steps stay in needs_run — launcher will build them

    all_ID, all_Session, all_data_path, *_ = \
        Load_subject_with_BIDS.load_data_bids(allinfo_study_c, bids_dir, [], [])

    status = []
    for ID, Session, data_path in zip(all_ID, all_Session, all_data_path):
        key = f'sub-{ID}_ses-{Session}'

        if _is_removed(ID, Session, list_to_remove):
            status.append({
                'ID': ID, 'Session': Session, 'key': key,
                'removed': True, 'complete': False,
                'steps': {}, 'effective_skip': [],
                'needs_run': [], 'error': None,
            })
            continue

        entry = {
            'ID': ID, 'Session': Session, 'key': key,
            'removed': False, 'complete': False,
            'steps': {}, 'effective_skip': [], 'needs_run': [], 'error': None,
        }

        try:
            paths = _anat_paths(data_path, reference)

            for step in ANAT_STEPS:
                sid = step['step_id']

                if sid in Skip_step:
                    entry['steps'][sid] = 'skipped_by_user'
                    entry['effective_skip'].append(sid)
                    continue

                try:
                    done = step['check'](paths, ID, type_norm, (species, bids_dir, reference))
                except Exception as e:
                    done = False
                    entry['error'] = f'Check failed step {sid}: {e}'

                entry['steps'][sid] = done
                if done:
                    entry['effective_skip'].append(sid)
                else:
                    entry['needs_run'].append(sid)

            entry['complete'] = len(entry['needs_run']) == 0

        except Exception as e:
            entry['error'] = str(e)

        status.append(entry)

    if verbose:
        _print_report(status, ANAT_STEPS, 'Anatomical',
                      Skip_step, list_to_remove)

    groups = _group_by_resume_point(status, ANAT_STEPS)

    if verbose:
        _print_groups(groups, ANAT_STEPS, 'Anatomical')

    return groups
