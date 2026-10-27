"""
Profile reliability (F22) and arousal-state resolution of each run.

StateResolver   state of a run (awake / anaesthetised) from the dataset config
profile_reliability / compare_methods   how reproducible each profile construction
                (median_r, percentiles, fraction_posneg) is across runs of one animal

Moved VERBATIM from Exemples/Study_EDNiX/analysis/Fig5/fig5_fingerprints.py (2026-10-27);
fig5_fingerprints imports these names back.
"""
import os
import re
import glob
import itertools
import contextlib

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats as st

from Statistics.Group_fMRI.run_selection import consensus as fa


opj = os.path.join


opb = os.path.basename


def run_token(filename):
    """Run id from a file name -- same rule as
    SBA_Region_Profile_Index._extract_run_id (kept in sync by a unit test)."""
    m = re.search(r"run-(RL|LR)-(\d+)", filename)
    if m:
        return f"{m.group(1)}-{m.group(2)}"
    m = re.search(r"run-(\d+)", filename)
    return m.group(1) if m else "NA"


def task_of(filename):
    m = re.search(r"task-([A-Za-z0-9]+)", filename)
    return m.group(1) if m else None


class StateResolver:
    """State of one unit key (subject, session, run, seed) from the rule in a
    BIDS_CONFIGS entry. With rule {"from": "task"} the raw BIDS func/ folder of
    that subject/session is listed once and cached."""

    def __init__(self, bids_config):
        self.rule = bids_config.get("state")
        self.bids_dir = bids_config.get("bids_dir")
        self._tasks = {}
        self.animal_re = (self.rule.get("animal") if isinstance(self.rule, dict) else None)

    def animal(self, key):
        """Animal identity of a unit. Default: the BIDS subject. With
        rule["animal"] = regex, the regex is removed from the subject label
        (case-insensitive result), e.g. Mouse "jgrAesISOc11L" -> "c11l": the
        state is encoded in the subject label but the animal is the same."""
        subj = str(key[0])
        if not self.animal_re:
            return subj
        return (re.sub(self.animal_re, "", subj) or subj).lower()

    def _task_pairs(self, subj, ses):
        key = (subj, ses)
        if key not in self._tasks:
            files = []
            for d in (opj(self.bids_dir, f"sub-{subj}", f"ses-{ses}", "func"),
                      opj(self.bids_dir, f"sub-{subj}", "func")):
                if os.path.isdir(d):
                    files = sorted(glob.glob(opj(d, "*_bold.nii*")))
                    if files:
                        break
            self._tasks[key] = [(run_token(opb(f)), task_of(opb(f))) for f in files]
        return self._tasks[key]

    def state(self, key):
        subj, ses, run = str(key[0]), str(key[1]), str(key[2])
        rule = self.rule
        if rule is None:
            return "unknown"
        if isinstance(rule, str):
            return rule
        src, mapping = rule.get("from"), rule.get("map")
        if src == "subject":
            items = mapping.items() if isinstance(mapping, dict) else mapping
            for pattern, state in items:          # first match wins: MEDISO before ISO
                if pattern in subj:
                    return state
            return "unknown"
        if src == "task":
            tasks = {t for r, t in self._task_pairs(subj, ses)
                     if t and (run == "ALL" or r == run)}
            states = {mapping.get(t, "unknown") for t in tasks}
            if not states:
                return "unknown"
            if len(states) == 1:
                return states.pop()
            return "mixed" if run == "ALL" else "ambiguous"
        raise ValueError(f"state rule {rule!r}: 'from' must be 'task' or 'subject'")


def discriminability(D, animal):
    """Bridgeford et al. 2021 (PLoS Comput Biol 17:e1009279): probability that
    a run is closer to another run of the SAME animal than to a run of a
    different animal. D = distance matrix (runs x runs), animal = label per
    run. Mean over every ordered same-animal pair (i, j) of the fraction of
    runs k of other animals with D[i, k] > D[i, j] (ties count 1/2). 0.5 =
    chance, 1 = every run identifies its animal."""
    animal = np.asarray(animal)
    num, den = 0.0, 0
    for i in range(len(animal)):
        same = np.flatnonzero((animal == animal[i]) & (np.arange(len(animal)) != i))
        other = np.flatnonzero(animal != animal[i])
        if not same.size or not other.size:
            continue
        dij = D[i, same][:, None]
        dik = D[i, other][None, :]
        num += float(np.sum((dik > dij).mean(axis=1) + 0.5 * (dik == dij).mean(axis=1)))
        den += same.size
    return num / den if den else np.nan


def _method_matrix(sri, profs, keys, common, method, norm, q):
    """Runs x values for one profile method, computed from the raw profiles
    (which carry the voxel sketch)."""
    if method == "median_r":
        names, vcol = list(common), "median_r"
        rows = [fa._ordered_common_values(profs[k], names, vcol) for k in keys]
    else:
        names, vcol = sri.channel_names(common, method, q), sri.profile_value_col(method, q)
        rows = [fa._ordered_common_values(sri.to_channel_profile(profs[k], method,
                                                                 reference=common, q=q),
                                          names, vcol) for k in keys]
    return np.vstack([fa._normalize_values(np.asarray(r, float), norm) for r in rows])


def profile_reliability(sri, runs_by_label, common, resolvers, norm, q=0.2,
                        methods=("median_r", "fraction_posneg", "percentiles"),
                        n_perm=500, seed=0):
    """Per dataset with several runs per animal and per profile method, on
    EVERY run (before any selection): median r inter-run (same session),
    same animal, different animals, the animal gap, and the discriminability
    (with its permutation p against chance: animal labels shuffled across
    runs). Returns a long DataFrame (one row per dataset x method)."""
    rng = np.random.default_rng(seed)
    rows = []
    for lab, profs in runs_by_label.items():
        keys = sorted(profs, key=lambda k: tuple(str(x) for x in k))
        if len(keys) < 6:
            continue
        skeys = [tuple(str(x) for x in k) for k in keys]
        res = resolvers.get(lab)
        animal = np.array([res.animal(k) if res else k[0] for k in skeys])
        session = np.array([f"{k[0]}|{k[1]}" for k in skeys])
        if len(set(animal)) < 2 or len(set(animal)) == len(animal):
            continue                                  # no repeated runs per animal
        C = np.where(animal[:, None] != animal[None, :], 2,
                     np.where(session[:, None] != session[None, :], 1, 0))
        iu = np.triu_indices(len(keys), 1)
        perms = [rng.permutation(animal) for _ in range(n_perm)]
        for method in methods:
            try:
                V = _method_matrix(sri, profs, keys, common, method, norm, q)
            except ValueError as e:
                print(f"  [F22] {lab} / {method}: {e}")
                continue
            ok = np.isfinite(V).all(axis=0)
            if ok.sum() < 3:
                continue
            R = np.corrcoef(V[:, ok])
            r, c = R[iu], C[iu]
            med = {k: float(np.median(r[c == v])) if np.any(c == v) else np.nan
                   for k, v in (("inter_run", 0), ("inter_session", 1), ("inter_animal", 2))}
            same = r[c <= 1]
            dist = 1.0 - R
            d_obs = discriminability(dist, animal)
            null = np.array([discriminability(dist, a) for a in perms])
            rows.append(dict(
                label=lab, method=method, n_runs=len(keys), n_animals=len(set(animal)),
                n_values=int(ok.sum()), r_inter_run=med["inter_run"],
                r_inter_session=med["inter_session"], r_inter_animal=med["inter_animal"],
                r_same_animal=float(np.median(same)) if same.size else np.nan,
                gap_animal=(float(np.median(same)) - med["inter_animal"]) if same.size
                else np.nan,
                discriminability=d_obs,
                p_discriminability=(1 + int(np.sum(null >= d_obs - 1e-12))) / (n_perm + 1)))
    return pd.DataFrame(rows)


def compare_methods(table, measure="discriminability", reference="median_r"):
    """Paired comparison of each method with `reference` over datasets: mean
    difference and exact two-sided sign-flip p (2^n assignments)."""
    out = {}
    if table.empty:
        return out
    w = table.pivot(index="label", columns="method", values=measure)
    if reference not in w:
        return out
    for m in w.columns:
        if m == reference:
            continue
        d = (w[m] - w[reference]).dropna().values
        if len(d) < 2:
            continue
        signs = np.array(list(itertools.product((1, -1), repeat=len(d))))
        null = np.abs((signs * d).mean(axis=1))
        out[m] = dict(mean_diff=float(d.mean()), n=len(d),
                      p=float(np.mean(null >= abs(d.mean()) - 1e-12)),
                      n_better=int(np.sum(d > 0)))
    return out
