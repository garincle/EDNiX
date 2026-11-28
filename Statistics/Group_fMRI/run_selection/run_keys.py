"""Run identifier of a BOLD / SBA-map file name: the one rule shared by
profiles._extract_run_id, reliability.run_token and
group_maps.collect_bold_and_means.

Why a shared rule. The run id used to be the run-<n> token alone. Two files
of the same subject/session that differ only by acq- (Babayan: acq-AP /
acq-PA, same run number) or by a non-'rest' task- label (Hori: restawake /
restanest) then got the SAME key, and in the profile dictionaries keyed by
(subject, session, run, seed) one silently overwrote the other, whereas the
group maps (built from the file list) kept both. The id now carries those
labels as a prefix:

    sub-1_ses-1_task-rest_acq-AP_run-01_bold      -> 'acqAP-01'
    sub-1_ses-1_task-restawake_run-1_bold         -> 'restawake-1'
    sub-1_ses-1_task-restawake_acq-AP_run-2_bold  -> 'restawake-acqAP-2'
    sub-1_task-rest_run-RL-4_bold                 -> 'RL-4'   (unchanged)
    sub-1_ses-1_task-rest_run-3_bold              -> '3'      (unchanged)

Files with no acq- label and a plain 'rest' task (every other paper dataset)
keep exactly their previous id, so their selections, maps and caches are
unaffected.
"""
import os
import re

__all__ = ["core_run", "variant_tags", "with_variants", "run_id"]

_RLLR = re.compile(r"run-(RL|LR)-(\d+)")
_RUN = re.compile(r"run-(\d+)")
_TASK = re.compile(r"task-([A-Za-z0-9]+)")
_ACQ = re.compile(r"acq-([A-Za-z0-9]+)")


def core_run(filename):
    """'RL-4' / 'LR-2' for direction-tagged runs, else the digits of run-<n>
    (kept as written, '01' stays '01'), else None."""
    name = os.path.basename(str(filename))
    m = _RLLR.search(name)
    if m:
        return f"{m.group(1)}-{m.group(2)}"
    m = _RUN.search(name)
    return m.group(1) if m else None


def variant_tags(filename):
    """Labels that distinguish two files with the same run number: the task-
    label when it is not plain 'rest', then 'acq' + the acq- label."""
    name = os.path.basename(str(filename))
    tags = []
    t = _TASK.search(name)
    if t and t.group(1).lower() != "rest":
        tags.append(t.group(1))
    a = _ACQ.search(name)
    if a:
        tags.append("acq" + a.group(1))
    return tags


def with_variants(filename, core):
    """`core` (a caller's legacy run id) prefixed by the variant labels of
    `filename`; `core` unchanged when the file has none. A missing core
    (None / 'NA') leaves the labels alone ('acqAP')."""
    tags = variant_tags(filename)
    if not tags:
        return core
    if core is None or str(core) == "NA":
        return "-".join(tags)
    return "-".join(tags + [str(core)])


def run_id(filename, default="NA"):
    """Run id of a file: core_run + variant labels; `default` when the name
    has neither a run token nor a variant label."""
    core = core_run(filename)
    out = with_variants(filename, core)
    return default if out is None else out
