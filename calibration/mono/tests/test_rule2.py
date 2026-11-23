"""Round 2 rule: the round-1 setting must reproduce st2core's verdicts, and each new part must act."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "round2"))
import monocal as mc  # noqa: E402
import rule2  # noqa: E402
from rule2 import ROUND1, Rule  # noqa: E402

core = mc.core


def _lc(seed=1, noise=3e-4):
    rng = np.random.default_rng(seed)
    t = np.concatenate([np.arange(0, 25, 10 / 1440), np.arange(700, 725, 10 / 1440)])
    return t, 1 + rng.normal(0, noise, t.size)


def _boxes(t, f, ep, dur, m=8):
    off = np.arange(-m, m + 1) * dur / 8
    c = np.atleast_1d(ep)[:, None] + off[None, :]
    d, fr, _ = core.box_stats(t, f, c.ravel(), dur)
    e = core.edge_distance(core.gap_edges(t), c.ravel())
    return d.reshape(c.shape), fr.reshape(c.shape), e.reshape(c.shape)


def test_round1_rule_matches_st2core():
    t, f = _lc()
    f[(t > 10) & (t < 10.2)] -= 2e-3
    f[(t > 705) & (t < 705.2)] += 3e-3          # a brightening
    dur, d0, e0 = 0.2, 2e-3, 2e-4
    ep = np.arange(1.0, 724, 0.37)
    v = core.epoch_verdicts(ep, t, f, dur, d0, e0)
    d, fr, ed = _boxes(t, f, ep, dur)
    sig = core.duration_scale_noise(t, f, dur)
    mine = rule2.epoch_excludes(d, fr, ed, sig, d0, e0, ROUND1)
    assert (mine == (v["verdict"] == "excludes")).all() and mine.any()


def test_floor_stops_a_brightening_from_excluding():
    t, f = _lc()
    f[(t > 705) & (t < 705.2)] += 3e-3
    dur, d0, e0 = 0.2, 2e-3, 2e-4
    d, fr, ed = _boxes(t, f, [705.1], dur)
    sig = core.duration_scale_noise(t, f, dur)
    assert rule2.epoch_excludes(d, fr, ed, sig, d0, e0, ROUND1)[0]
    assert not rule2.epoch_excludes(d, fr, ed, sig, d0, e0, Rule(floor_sig=3))[0]


def test_weak_reference_and_sensitivity():
    t, f = _lc()
    dur = 0.2
    d, fr, ed = _boxes(t, f, [12.0], dur)
    sig = core.duration_scale_noise(t, f, dur)
    assert rule2.epoch_excludes(d, fr, ed, sig, 2e-3, 2e-4, ROUND1)[0]
    assert not rule2.epoch_excludes(d, fr, ed, sig, 2e-3, 2e-4, Rule(weak_ref_excludes=False), ref_measured=False)[0]
    # a reference too shallow for the noise: the threshold is negative, so nothing may exclude
    assert not rule2.epoch_excludes(d, fr, ed, sig, 1e-4, 1e-5, Rule(min_thr_sig=0))[0]


def test_timing_tolerance_sees_a_shifted_transit():
    t, f = _lc()
    dur, d0, e0 = 0.2, 2e-3, 2e-4
    f[np.abs(t - 12.19) < dur / 2] -= d0        # the transit arrives 0.95 durations late
    d, fr, ed = _boxes(t, f, [12.0], dur)
    sig = core.duration_scale_noise(t, f, dur)
    assert rule2.epoch_excludes(d, fr, ed, sig, d0, e0, ROUND1)[0]
    assert not rule2.epoch_excludes(d, fr, ed, sig, d0, e0, Rule(w_dur=1.0))[0]
    assert rule2.epoch_excludes(*_boxes(t, f, [5.0], dur), sig, d0, e0, Rule(w_dur=1.0))[0]


def test_timeline_agrees_with_epochs():
    t, f = _lc()
    dur, d0, e0 = 0.2, 2e-3, 2e-4
    f[np.abs(t - 12.15) < dur / 2] -= d0
    step = dur / 8
    c = np.arange(t.min(), t.max() + step, step)
    dep, fr, _ = core.box_stats(t, f, c, dur)
    keep = fr >= 0.5
    cl = (fr[keep] >= core.CLEAN_FRAC) & (core.edge_distance(core.gap_edges(t), c[keep]) >= core.CLEAN_EDGE_D)
    sig = core.duration_scale_noise(t, f, dur)
    for rule in (ROUND1, Rule(floor_sig=3, w_dur=0.5), Rule(w_dur=1.0)):
        tl = rule2.timeline_excludes(dep[keep], fr[keep], cl, np.flatnonzero(keep), c.size, sig, d0, e0, rule)
        i = np.arange(100, c.size - 100, 7)
        ep = rule2.epoch_excludes(*_boxes(t, f, c[i], dur), sig, d0, e0, rule)
        assert (tl[i] == ep).all()


def test_vetting():
    ev = dict(snr=9.0, depth_ratio=1.6, dur_ratio=1.0, anti_snr=3.0, local_anti_snr=4.5,
              edge_dist_d=3.0, frac=1.0)
    assert rule2.vet(ev, ROUND1)
    assert not rule2.vet(ev, Rule(snr_min=10))
    assert not rule2.vet(ev, Rule(depth_ratio_max=1.5))
    assert not rule2.vet(ev, Rule(local_anti_max=4))
    assert rule2.best_event([ev, dict(ev, snr=20.0, frac=0.6)], ROUND1)["snr"] == 9.0


def test_harness_boxes_and_alias_table_match_st2core_under_round1_rule():
    import harness
    t, f = _lc(seed=3)
    dur, d0, e0 = 0.2, 2e-3, 2e-4
    for tc in (10.0, 710.0):
        f[np.abs(t - tc) < dur / 2] -= d0
    box = harness.Boxes(t, f)
    c = np.linspace(1, 720, 500)
    d, fr, ed = box(c, dur)
    d2, fr2, _ = core.box_stats(t, f, c, dur)
    assert np.allclose(d, d2, equal_nan=True) and np.allclose(fr, fr2)
    assert np.allclose(ed, core.edge_distance(core.gap_edges(t), c))
    sig = core.duration_scale_noise(t, f, dur)
    mine = harness.alias_table(10.0, 710.0, box, t.min(), t.max(), dur, d0, e0, sig, ROUND1, True, False, 5.0)
    ref = core.alias_table(10.0, 710.0, dur, d0, e0, t, f, p_floor=1.0, pmin_density=5.0)
    assert [m[2] for m in mine] == [r["status"] for r in ref]
    assert {"excluded", "allowed"} <= {m[2] for m in mine}
