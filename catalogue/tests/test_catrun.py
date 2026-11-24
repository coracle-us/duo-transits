"""catrun on synthetic light curves: it reproduces the frozen round-3 harness for a single transit,
and the duotransit table keeps the true period."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "catalogue"))
import catrun  # noqa: E402

h3 = catrun.h3
P, T0, DUR, DEPTH = 9.3, 1403.1, 0.2, 0.006


def _sectors(seed=3):
    rng = np.random.default_rng(seed)
    out = []
    for s, start in ((5, 1400.0), (6, 1428.0), (32, 2130.0), (33, 2158.0)):
        t = np.concatenate([np.arange(start, start + 12.5, 1 / 144), np.arange(start + 14, start + 26.5, 1 / 144)])
        f = 1 + 2e-3 * np.sin(2 * np.pi * t / 6.1) + rng.normal(0, 6e-4, t.size)
        k = np.round((t - T0) / P)
        f[np.abs(t - T0 - k * P) < DUR / 2] -= DEPTH
        out.append((s, "SPOC", t, f))
    return out


def _def(**kw):
    d = dict(row_id="x", tic=1, arm="A2", t_cat=T0 + 0.004, dur_cat_d=DUR, depth_cat=DEPTH, depth_cat_err=np.nan,
             truth_P=P, truth_t0=T0, mstar=1.0, rstar=1.0, multi=False, t2_cat=np.nan)
    d.update(kw)
    return d


@pytest.fixture(scope="module")
def searched():
    sec = _sectors()
    d = _def()
    return sec, d, catrun.search_row(d, sec)


def test_single_matches_frozen_harness(searched):
    sec, d, (row, ev, ctx) = searched
    ref = h3.score_trial(d, sec)
    assert row["frac_lnP_excluded"] == ref["frac_lnP_excluded"]
    assert row["ref_measured"] == ref["ref_measured"] and row["ref_snr"] == ref["ref_snr"]
    cols, tabs = catrun.finish_row(d, ctx, ev, {e["event_id"] for e in ev if e["passes"]})
    assert sum(e["passes"] for e in ev) == ref["n_vetted"] > 0
    for k in ("call", "n_alias_ok", "n_alias", "t2", "frac_alias_excluded"):
        assert cols[k] == ref[k], k
    n = int(round(abs(cols["t2"] - d["t_cat"]) / P))
    assert tabs[0][3][n - 1][2] != "excluded"


def test_no_accepted_event_is_no_repeat(searched):
    _, d, (row, ev, ctx) = searched
    cols, tabs = catrun.finish_row(d, ctx, ev, set())
    assert cols["call"] == "no repeat" and not tabs
    assert 0 < row["frac_lnP_excluded"] <= row["frac_lnP_excluded_testable"] or row["p_longest_testable_d"] > 1000


def test_excluded_runs_cover_the_excluded_fraction(searched):
    _, d, (row, ev, ctx) = searched
    assert ctx["excluded_runs"] and all(a <= b for a, b in ctx["excluded_runs"])
    assert not any(a <= P <= b for a, b in ctx["excluded_runs"])       # the true period is not excluded


def test_duotransit_keeps_true_period_and_third_event():
    sec = _sectors()
    k2 = 80                                                             # a transit in the second pair of sectors
    d = _def(t2_cat=T0 + k2 * P + 0.003)
    row, ev, ctx = catrun.search_row(d, sec)
    assert any(e["at_t2_cat"] for e in ev)
    acc = {e["event_id"] for e in ev if e["passes"]}
    cols, tabs = catrun.finish_row(d, ctx, ev, acc)
    kind, t1, t2, rows = tabs[0]
    assert kind == "duo" and cols["second_recovered"]
    assert rows[k2 - 1][2] != "excluded"
    assert cols["n_alias_ok"] < cols["n_alias"]
    if "t3" in cols:                                                    # a third transit was found and accepted
        assert rows[k2 - 1][3] in ("yes", "unconstrained")
        assert tabs[1][0] == "ref_third"


def test_stated_period_check(searched):
    _, d, (row, ev, ctx) = searched
    exc, nep = catrun.period_excluded(ctx, d["t_cat"], P)
    assert not exc and nep > 5
    assert catrun.period_excluded(ctx, d["t_cat"], 13.37)[0]             # a wrong period with clean epochs
    cols, tabs = catrun.finish_row(d, ctx, ev, {e["event_id"] for e in ev if e["passes"]})
    assert catrun.alias_status(ctx, tabs, P)["single"] != "excluded"
    assert catrun.alias_status(ctx, tabs, 13.37)["single"] == "not an alias"
