"""
Round 3 harness: one trial scored exactly under the round-3 rule (rule3.Rule3), from light curves.

A trial is prepared exactly as in round 2 (round2/harness.py, frozen): same light curves, masking,
injection, detrending, reference, event search and null run. Its single-transit exclusion columns
(truth_excluded, frac_lnP_excluded, ...) come from round2/harness.score unchanged, and that
function's round-2 call is kept as r2_* columns. The round-3 columns come from rule3.call3.

Columns added by round 3:
  n_vetted, false_event, repeat_detected, t2, call, n_alias_ok, n_alias, truth_in_alias_set,
  event_true, frac_alias_excluded, frac_alias_excluded_r2 (round 2's table through the catalogued
  time and the same event), and for arm B the same with suffix _own (the host's own null-run
  events removed first).

Inputs  : trial definitions (CSV), MAST light curves (cache), gaia_dr4/candidates_frozen.csv
Outputs : rows CSV
Usage   : python harness3.py --defs defs.csv --out rows.csv [--workers 12]
          python harness3.py --make-defs a|b --sample s.csv [--seed N --noise 0.1] --out defs.csv
"""
import argparse
import multiprocessing as mp
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "round2"))
import harness as h2  # noqa: E402
import dump3  # noqa: E402
import rule3  # noqa: E402
from rule3 import Rule3  # noqa: E402

mc, core, dump = h2.mc, h2.core, h2.dump
CACHE = dump3.CACHE

# The rule proposed for the freeze (FROZEN_RULE3.md). Not frozen until APPROVALS.md says so.
FROZEN3 = Rule3(c_drift=2.0, sig_floor_min=15.0, cap_dur=2.0, sec_box_snr_min=9.0)


def own_removed(events, null_events, dur):
    nt = np.array([e["tc"] for e in null_events])
    return [e for e in events
            if e.get("at_truth", False) or not (nt.size and np.min(np.abs(nt - e["tc"])) < dur)]


def score_trial(d, sectors, rule=FROZEN3):
    tic = int(d["tic"])
    prep, null_events = dump3.prep_trial(d, sectors, tic)
    if prep is None:
        return dict(d, status="reference not in data")
    r2 = h2.score(prep, d["t_cat"], d["truth_P"], d["truth_t0"], d["mstar"], d["rstar"], bool(d["multi"]),
                  h2.FROZEN, null_events)
    call_keys = ("n_vetted", "false_event", "repeat_detected", "t2", "call", "n_alias_ok", "truth_in_alias_set")
    out = dict(d)
    for k, v in r2.items():
        base = k[:-4] if k.endswith("_own") else k
        out[("r2_" + k) if base in call_keys else k] = v
    s = dump.summarize(prep, d["t_cat"], d["truth_P"], d["truth_t0"], grid=False)
    T, F, dur = prep["T"], prep["F"], prep["dur_ref"]
    ev = [e for e in s["events"] if "error" not in e]
    for e in ev:
        e["sec_box_snr"] = dump3.sector_box_snr(prep["det"], e, dur)
    tr_, fr_ = prep["det"][s["ref_sector"]]
    ref = core.refine_event(tr_, fr_, d["t_cat"], d["dur_cat_d"])
    tr = dict(d=d, dur=dur, depth0=prep["depth_ref"], derr=prep["derr_ref"], sig=s["sig"],
              ref_measured=s["ref_measured"], ref_snr=s["ref_snr"],
              ref_tc=ref["tc"] if ref else np.nan, ref_dur=ref["dur_d"] if ref else np.nan,
              tmin=float(T.min()), tmax=float(T.max()),
              pmin=core.pmin_from_duration(dur, d["mstar"], d["rstar"], prep["depth_ref"]),
              pmin_cat=core.pmin_from_duration(d["dur_cat_d"], d["mstar"], d["rstar"], prep["depth_ref"]))
    box = h2.Boxes(T, F)
    by_tc = {e["tc"]: e for e in ev}

    def one(events):
        o = rule3.call3(events, tr, rule, box, d["truth_P"])
        o["frac_alias_excluded_r2"] = np.nan
        if o["call"] in ("resolved", "narrowed", "all excluded"):
            t1, t2 = sorted([d["t_cat"], o["t2"]])
            rr = rule3.alias_table3(t1, t2, 0.0, box, tr["tmin"], tr["tmax"], dur, tr["depth0"],
                                    by_tc[o["t2"]]["depth_err"], tr["sig"], tr["ref_measured"],
                                    bool(d["multi"]), tr["pmin"], rule3.ROUND2)
            o["frac_alias_excluded_r2"] = float(np.mean([x[2] == "excluded" for x in rr])) if rr else np.nan
        return o

    out.update(one(ev))
    if null_events is not None:
        out.update({k + "_own": v for k, v in one(own_removed(ev, null_events, dur)).items()})
    return out


def run_star(defs, cache=CACHE, rule=FROZEN3):
    sectors, _ = mc.load_star(int(defs[0]["tic"]), Path(cache))
    if not sectors:
        return [dict(d, status="no light curve") for d in defs]
    return [score_trial(d, sectors, rule) for d in defs]


def add_noise(defs, noise, seed):
    """Arm B: Gaussian error of `noise` (fraction) on the catalogued depth and duration, one generator
    per star (numpy seed [seed, tic]), two draws per definition in file order (depth, then duration).
    The injected values stay in inj_depth and inj_dur_d."""
    out = []
    for tic, g in pd.DataFrame(defs).groupby("tic", sort=False):
        rng = np.random.default_rng([seed, int(tic)])
        for d in g.to_dict("records"):
            d["inj_depth"], d["inj_dur_d"] = d["depth_cat"], d["dur_cat_d"]
            d["depth_cat"] *= max(0.2, 1 + noise * rng.normal())
            d["dur_cat_d"] *= max(0.2, 1 + noise * rng.normal())
            out.append(d)
    return out


def _star(args):
    fn, a = args
    try:
        return fn(*a)
    except Exception as e:
        return dict(error=f"{type(e).__name__}: {str(e)[:300]}", args=str(a[0])[:200])


def _pool(fn, arglist, workers):
    with mp.get_context("spawn").Pool(workers) as pool:
        for i, r in enumerate(pool.imap(_star, [(fn, a) for a in arglist]), 1):
            print(f"[{i}/{len(arglist)}]", flush=True)   # progress only, no outcome
            yield r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--defs")
    ap.add_argument("--make-defs", choices=["a", "b"])
    ap.add_argument("--sample")
    ap.add_argument("--seed", type=int, default=20261203)
    ap.add_argument("--noise", type=float, default=0.0)
    ap.add_argument("--noise-seed", type=int, default=20261204)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cache", default=CACHE)
    ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args()
    if a.make_defs:
        items = pd.read_csv(a.sample).to_dict("records")
        fn = h2.make_defs_a if a.make_defs == "a" else h2.make_defs_b
        args = [(it, a.cache) if a.make_defs == "a" else (it, a.cache, a.seed) for it in items]
        defs, notes, errs = [], [], []
        for it, r in zip(items, _pool(fn, args, a.workers)):
            if isinstance(r, dict):
                errs.append(r)
                notes.append(dict(tic=it["tic"], note="error"))
                continue
            defs += r[0]
            notes.append(dict(tic=it["tic"], note=r[1]))
        if a.noise and defs:
            defs = add_noise(defs, a.noise, a.noise_seed)
        pd.DataFrame(defs).to_csv(a.out, index=False)
        pd.DataFrame(notes).to_csv(Path(a.out).with_suffix(".notes.csv"), index=False)
        if errs:
            pd.DataFrame(errs).to_csv(Path(a.out).with_suffix(".errors.csv"), index=False)
        return
    df = pd.read_csv(a.defs)
    groups = [(g.to_dict("records"), a.cache) for _, g in df.groupby("tic", sort=False)]
    rows, errs = [], []
    for r in _pool(run_star, groups, a.workers):
        (errs.append(r) if isinstance(r, dict) else rows.extend(r))
    pd.DataFrame(rows).to_csv(a.out, index=False)
    if errs:
        pd.DataFrame(errs).to_csv(Path(a.out).with_suffix(".errors.csv"), index=False)


if __name__ == "__main__":
    main()
