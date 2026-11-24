"""Synthetic checks of the calibration harness. No network, no real light curves."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import monocal as mc  # noqa: E402

CAD = 10 / 1440
DUR, DEPTH, P_TRUE, T0 = 4 / 24, 5000e-6, 20.0, 5.0


def raw_sectors(seed=1, starts=(0.0, 400.0, 800.0), noise=300e-6):
    """Three 27 d sectors with a slow trend and white noise, no transit."""
    rng = np.random.default_rng(seed)
    out = []
    for i, s0 in enumerate(starts):
        t = np.arange(s0, s0 + 27, CAD)
        f = 1 + 0.003 * np.sin(2 * np.pi * (t - s0) / 9.0) + rng.normal(0, noise, t.size)
        out.append((i + 1, "QLP", t, f))
    return out


def with_box_transits(sectors, t0=T0, P=P_TRUE, depth=DEPTH, dur=DUR, only=None):
    out = []
    for s, a, t, f in sectors:
        f = f.copy()
        _, ep = mc.epochs(t0, P, t.min(), t.max())
        for e in ep if only is None else only:
            f[np.abs(t - e) < dur / 2] -= depth
        out.append((s, a, t, f))
    return out


def test_duration_matches_the_pipeline_at_b_zero_and_shrinks_with_b():
    d0 = mc.transit_duration_d(50.0, 1.0, 1.0, 0.1, 0.0)
    assert d0 == pytest.approx(mc.core.max_duration_d(50.0, 1.0, 1.0, 0.1), rel=1e-3)
    assert mc.transit_duration_d(50.0, 1.0, 1.0, 0.1, 0.8) < d0


def test_batman_model_has_the_asked_depth_and_duration():
    t = np.arange(0, 10, CAD)
    ar = mc.core.a_over_rstar(50.0, 1.0, 1.0)
    m = mc.transit_flux(t, 5.0, 50.0, 0.1, ar, 0.0, (0.0, 0.0), CAD)
    assert 1 - m.min() == pytest.approx(0.01, rel=0.02)
    width = (m < 1 - 1e-6).sum() * CAD
    assert width == pytest.approx(mc.transit_duration_d(50.0, 1.0, 1.0, 0.1, 0.0), abs=3 * CAD)


def test_mask_windows_removes_only_the_windows():
    t = np.arange(0, 10, 0.1)
    t2, _ = mc.mask_windows(t, t, [5.0], 0.25)
    assert t.size - t2.size == 5 and not np.any(np.abs(t2 - 5.0) <= 0.25)


def test_excluded_periods_uses_every_epoch_but_the_reference():
    step, c0 = 0.02, 0.0
    excl = np.zeros(5001, bool)            # timeline 0..100 d
    excl[2495:2506] = True                 # flat data proven around t = 50
    ex = mc.excluded_periods(0.0, [50.0, 25.0, 30.0, 100.0, 10.0], c0, step, excl, 0.0, 100.0)
    assert ex.tolist() == [True, True, False, False, True]
    # chunking gives the same answer
    per = np.exp(np.linspace(0, np.log(100), 3000))
    a = mc.excluded_periods(0.0, per, c0, step, excl, 0.0, 100.0)
    b = mc.excluded_periods(0.0, per, c0, step, excl, 0.0, 100.0, max_elems=500)
    assert (a == b).all() and a.any() and not a.all()


def test_visible_repeats_are_found_and_the_true_period_survives():
    r = mc.run_trial(with_box_transits(raw_sectors()), T0, DUR, DEPTH, truth_P=P_TRUE)
    assert r["status"] == "ok" and r["ref_snr"] > 10
    assert r["n_clean_repeats"] >= 3 and r["repeat_detected"]
    assert r["truth_in_alias_set"] and not r["truth_excluded_single"] and not r["truth_excluded_grid"]


def test_a_lone_transit_excludes_periods_that_land_in_flat_data_only():
    sec = with_box_transits(raw_sectors(), only=[T0])
    # 20 d would put a transit at t = 25 in flat data; 500 d puts every other epoch in gaps
    wrong = mc.run_trial(sec, T0, DUR, DEPTH, truth_P=20.0)
    assert wrong["call"] == "no repeat" and wrong["truth_excluded_single"] and wrong["truth_excluded_grid"]
    right = mc.run_trial(sec, T0, DUR, DEPTH, truth_P=500.0)
    assert not right["truth_excluded_single"] and not right["truth_excluded_grid"]
    assert right["n_clean_repeats"] == 0 and 0.2 < right["frac_lnP_excluded"] < 1.0


def test_a_shallow_reference_excludes_little():
    """At low SNR, half the depth minus 3 sigma is below zero, so only a strong brightening at a
    predicted epoch excludes. Far less of ln P is excluded than for the deep transit above."""
    sec = with_box_transits(raw_sectors(noise=2000e-6), only=[T0], depth=600e-6)
    r = mc.run_trial(sec, T0, DUR, 600e-6, truth_P=500.0)
    assert not r["truth_excluded_single"] and r["frac_lnP_excluded"] < 0.1


def test_limb_darkening_comes_from_the_claret_grid():
    u1, u2, how = mc.limb_dark(5750.0, 1.0, 1.0)      # log g 4.44 -> grid point 4.5
    assert (u1, u2) == (0.3075, 0.2743) and how == "claret"
    cool = mc.limb_dark(3000.0, 0.3, 0.3)
    assert "below the grid" in cool[2]
    assert "assumed" in mc.limb_dark(float("nan"), float("nan"), 1.0)[2]
    assert mc.limb_dark(4000.0, 0.6, 0.6)[0] > u1      # cooler stars are more limb darkened


def test_scoring_counts_rates_and_judges_the_thresholds():
    import pandas as pd
    import score
    n = 100
    t = pd.DataFrame(dict(status="ok", ref_snr=[10.0] * n, truth_excluded_single=[True] * 4 + [False] * 96,
                          truth_n_clean_epochs=[1] * 50 + [0] * 50, n_clean_repeats=[1] * 40 + [0] * 60,
                          repeat_detected=[True] * 30 + [False] * 70, call=["narrowed"] * 20 + ["no repeat"] * 80,
                          truth_in_alias_set=[False] * 1 + [True] * 19 + [None] * 80, frac_lnP_excluded=0.5))
    r = score.score(t)
    g = lambda m, b="all": r[(r.metric == m) & (r.snr_bin == b)].iloc[0]
    assert (g("false_exclusion").k, g("false_exclusion").n) == (4, 100)
    assert g("false_exclusion", "8-12").n == 100 and g("false_exclusion", ">20").n == 0
    assert (g("false_exclusion_given_a_clean_epoch").k, g("false_exclusion_given_a_clean_epoch").n) == (4, 50)
    assert (g("repeat_recovery").k, g("repeat_recovery").n) == (30, 40)
    assert (g("wrong_period").k, g("wrong_period").n) == (1, 20)
    lines, ok = score.verdict(r, "B", True)
    assert not ok and any("wrong period" in x and "FAIL" in x for x in lines)   # 1 of 20 = 5% > 3%
    assert score.upper95(0, 59) == pytest.approx(0.0495, abs=1e-3)


def test_an_event_before_the_reference_still_gives_an_alias_table():
    """Reference = the last visible transit; the repeats are earlier. (Found by the independent check.)"""
    r = mc.run_trial(with_box_transits(raw_sectors()), 825.0, DUR, DEPTH, truth_P=P_TRUE, truth_t0=T0)
    assert r["t2"] < 825.0 and r["n_aliases"] > 0 and r["truth_in_alias_set"]
