"""Supp. Fig. 16 (awake vs anaesthetised): statistics and drawing on synthetic profiles."""
import numpy as np
import pandas as pd

from Paper.analysis import paper_stats as pst
from Paper.analysis.fmri_fingerprints import fingerprints as ff
from Paper.analysis.fmri_fingerprints import state_maps_figure as smf


def _synthetic():
    rng = np.random.default_rng(0)
    common = [f"{r}|q80" for r in ff.REGION_SHORT] + [f"{r}|q20" for r in ff.REGION_SHORT]
    n = len(common)
    sol = {"Mouse": "Mouse", "Rat": "Rat", "Human_X": "Human"}
    group = {lab: rng.normal(0, 1, n) for lab in sol}
    rows, vecs = [], {}
    for j, (subj, state) in enumerate([("AWa1", "awake"), ("AWa2", "awake"), ("AWa3", "awake"),
                                       ("ISOa1", "isoflurane"), ("ISOa2", "isoflurane")]):
        key = (subj, "1", "1", "seed")
        vecs.setdefault("Mouse", {})[key] = group["Mouse"] + rng.normal(0, 0.3, n)
        rows.append(dict(label="Mouse", species="Mouse", subject=subj, session="1", run="1",
                         key=key, animal=subj[-2:], state=state,
                         state_class=ff.state_class(state), kept=True))
    return pd.DataFrame(rows), vecs, group, sol, common


def test_compute_counts_and_pairs():
    runs, vecs, group, sol, common = _synthetic()
    res = smf.compute(runs, vecs, group, sol)
    s = res["summary"].iloc[0]
    assert (s.n_runs_awake, s.n_runs_anesthetized) == (3, 2)
    assert s.n_animals_both_states == 2
    kinds = res["pairs"]["kind"].value_counts()
    # 5 animal x state profiles: 10 pairs = 3+1 same state, 2 same animal, 4 other
    assert kinds[smf.PAIR_KINDS[0]] == 4 and kinds[smf.PAIR_KINDS[1]] == 2
    assert kinds[smf.PAIR_KINDS[2]] == 4
    assert 0.5 < s.r_states <= 1


def test_draw_without_maps(tmp_path):
    runs, vecs, group, sol, common = _synthetic()
    res = smf.compute(runs, vecs, group, sol)
    base = smf.draw(res, common, None, None, str(tmp_path / "S16"))
    assert (tmp_path / "S16.png").is_file() and (tmp_path / "S16_legend.txt").is_file()
    assert base


def test_section_registered():
    assert "supp16" in pst.SECTIONS and pst.FOLDER["supp16"] == "Supp_Fig_16"
