"""Round 3 rule: with every round-3 switch off it is round 2's alias table; the window grows as stated."""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "round3"))
sys.path.insert(0, str(HERE / "round2"))
import harness as h2  # noqa: E402
import rule3  # noqa: E402


def _lc(P=7.0, t0=3.0, dur=0.15, depth=0.004, seed=1, drift=0.0):
    rng = np.random.default_rng(seed)
    t = np.concatenate([np.arange(0, 27, 1 / 144), np.arange(700, 754, 1 / 144)])
    f = 1 + rng.normal(0, 3e-4, t.size)
    k = np.round((t - t0) / P)
    f[np.abs(t - t0 - k * P) < dur / 2] -= depth
    return t, f


def _args(t, f, dur=0.15, depth=0.004):
    return dict(box=h2.Boxes(t, f), tmin=t.min(), tmax=t.max(), dur=dur, depth0=depth, derr=2e-4,
                sig=h2.core.duration_scale_noise(t, f, dur), ref_measured=True, multi=False, pmin=np.nan)


def test_switches_off_is_round2():
    t, f = _lc()
    a = _args(t, f)
    for t2 in (3.0 + 7.0 * 2, 3.0 + 7.0 * 3 + 0.01, 710.0):
        new = rule3.alias_table3(3.0, t2, 0.0, r=rule3.ROUND2, **a)
        old = h2.alias_table(3.0, t2, rule=h2.FROZEN, **a)
        assert [(n, s) for n, _, s in new] == [(n, s) for n, _, s in old]


def test_drift_window_keeps_the_true_alias():
    """Second transit timed 25 min late, two orbits on: the true alias drifts off the far transits."""
    t, f = _lc()
    a = _args(t, f)
    t2 = 3.0 + 14.0 + 25 / 1440
    old = rule3.alias_table3(3.0, t2, 0.0, r=rule3.ROUND2, **a)
    new = rule3.alias_table3(3.0, t2, np.hypot(15, 15) / 1440, r=rule3.Rule3(c_drift=2.0, cap_dur=2.0), **a)
    assert old[1][2] == "excluded" and new[1][2] != "excluded"
    assert sum(s == "excluded" for _, _, s in new) >= 0.8 * sum(s == "excluded" for _, _, s in old)


def test_vetting_adds_the_sector_box_cut():
    ev = dict(snr=12.0, depth_ratio=1.0, dur_ratio=1.0, anti_snr=2.0, edge_dist_d=2.0, frac=1.0, sec_box_snr=8.5)
    assert rule3.vet3(ev, rule3.ROUND2, 20.0)
    assert not rule3.vet3(ev, rule3.Rule3(sec_box_snr_min=9.0), 20.0)
    assert rule3.vet3(dict(ev, sec_box_snr=9.5), rule3.Rule3(sec_box_snr_min=9.0), 20.0)
