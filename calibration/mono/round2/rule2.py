"""
Round 2 of the single-transit calibration: the rule with its tunable parts (prereg_round2.yaml,
tuning.what_may_change R1-R4). ROUND1 reproduces the rule that failed in round 1.

Everything here works on arrays of box depth, coverage and edge distance, so the same functions
score a stored summary (tune.py) and a light curve (the frozen harness).

No file input or output. Usage: imported by dump.py, tune.py and the tests.
"""
from dataclasses import dataclass

import numpy as np

CLEAN_FRAC, CLEAN_EDGE_D, MIN_FRAC = 0.8, 0.5, 0.5   # st2core's clean-epoch rule, not tunable


@dataclass(frozen=True)
class Rule:
    # R1 floor: an excluding box must not be a brightening. depth >= -floor_sig * (epoch noise).
    floor_sig: float = np.inf
    # R1 sensitivity: the threshold 0.5 depth0 - 3 sigma must itself be >= min_thr_sig * (epoch noise).
    min_thr_sig: float = -np.inf
    # R2: may a reference the pipeline could not measure (catalogue depth used) exclude?
    weak_ref_excludes: bool = True
    # R3: every covered box within +-w_dur durations of the epoch must also be below the threshold.
    w_dur: float = 0.0
    w_dur_multi: float = -1.0      # for stars with more than one planet candidate; -1 = same as w_dur
    # R4 vetting of a second event (round 1: 7, 2, 2, 1.5, off)
    snr_min: float = 7.0
    depth_ratio_max: float = 2.0   # accepted range is 1/x to x
    dur_ratio_max: float = 2.0
    anti_ratio: float = 1.5
    local_anti_max: float = np.inf  # veto if a brightening within 1 d has SNR above this
    # NOT in the approved R4 list (prereg_round2.yaml 1.0); explored on the tuning set only, and
    # usable in a frozen rule only if Jon approves a new version of the pre-registration:
    require_clean: bool = False     # the event itself must be clean (coverage and edge distance)


ROUND1 = Rule()


def threshold(err, depth0, depth0_err):
    return 0.5 * depth0 - 3 * np.hypot(err, depth0_err)


def box_excludes(depth, frac, sig, depth0, depth0_err, rule, ref_measured=True):
    """True where a covered box is below the threshold and passes the R1 and R2 conditions.
    Coverage and edge distance are judged separately (clean)."""
    depth, frac = np.asarray(depth, float), np.asarray(frac, float)
    err = sig / np.sqrt(np.clip(frac, 1e-3, None))
    thr = threshold(err, depth0, depth0_err)
    ok = (depth < thr) & (depth >= -rule.floor_sig * err) & (thr >= rule.min_thr_sig * err)
    if not rule.weak_ref_excludes and not ref_measured:
        ok &= False
    return ok & (frac >= MIN_FRAC)


def clean(frac, edge_d):
    return (np.asarray(frac) >= CLEAN_FRAC) & (np.asarray(edge_d) >= CLEAN_EDGE_D)


def w_steps(rule, step_dur_frac, multi=False):
    w = rule.w_dur_multi if (multi and rule.w_dur_multi >= 0) else rule.w_dur
    return int(round(w / step_dur_frac))


def epoch_excludes(depth, frac, edge_d, sig, depth0, depth0_err, rule, ref_measured=True, multi=False,
                   step_dur_frac=0.125):
    """Verdict for epochs given boxes at offsets around each one.

    depth, frac, edge_d: arrays (n_epochs, 2m+1), column m is the predicted time and neighbouring
    columns are step_dur_frac durations apart. An epoch excludes if its central box is clean and
    excludes, and every covered box within +-w is also below its own threshold (R3).
    """
    depth, frac, edge_d = (np.atleast_2d(x) for x in (depth, frac, edge_d))
    m = depth.shape[1] // 2
    ex = box_excludes(depth, frac, sig, depth0, depth0_err, rule, ref_measured)
    out = ex[:, m] & clean(frac[:, m], edge_d[:, m])
    w = w_steps(rule, step_dur_frac, multi)
    if w > m:
        raise ValueError("offsets do not reach the timing tolerance")
    if w:
        err = sig / np.sqrt(np.clip(frac, 1e-3, None))
        deep = (frac >= MIN_FRAC) & (depth >= threshold(err, depth0, depth0_err))
        out &= ~deep[:, m - w:m + w + 1].any(axis=1)
    return out


def timeline_excludes(depth, frac, is_clean, covered_idx, n, sig, depth0, depth0_err, rule,
                      ref_measured=True, multi=False, step_dur_frac=0.125):
    """The same verdict on a regular time grid (step = step_dur_frac durations): bool array of
    length n. depth, frac, is_clean are given for the covered grid points covered_idx only."""
    ex = np.zeros(n, bool)
    ex[covered_idx] = box_excludes(depth, frac, sig, depth0, depth0_err, rule, ref_measured) & is_clean
    w = w_steps(rule, step_dur_frac, multi)
    if w:
        err = sig / np.sqrt(np.clip(frac, 1e-3, None))
        deep = np.zeros(n, bool)
        deep[covered_idx] = depth >= threshold(err, depth0, depth0_err)
        c = np.concatenate([[0], np.cumsum(deep)])
        i = np.arange(n)
        ex &= (c[np.minimum(i + w + 1, n)] - c[np.maximum(i - w, 0)]) == 0
    return ex


def vet(ev, rule):
    """R4: does a second event pass the search cuts and the vetting?"""
    if "error" in ev or ev.get("snr") is None:
        return False
    return bool(ev["snr"] >= rule.snr_min
                and 1 / rule.depth_ratio_max <= ev["depth_ratio"] <= rule.depth_ratio_max
                and 1 / rule.dur_ratio_max <= ev["dur_ratio"] <= rule.dur_ratio_max
                and ev["snr"] >= rule.anti_ratio * max(ev["anti_snr"], 1e-9)
                and ev.get("local_anti_snr", 0.0) <= rule.local_anti_max
                and (not rule.require_clean
                     or (ev["edge_dist_d"] >= CLEAN_EDGE_D and ev["frac"] >= CLEAN_FRAC)))


def best_event(events, rule):
    """As st2_pipeline.alias_events: the highest-SNR clean vetted event, else the highest-SNR one."""
    ok = [e for e in events if vet(e, rule)]
    if not ok:
        return None
    return max(ok, key=lambda e: (e["edge_dist_d"] >= CLEAN_EDGE_D and e["frac"] >= CLEAN_FRAC, e["snr"]))
