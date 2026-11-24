"""
Round 3 tuning: store, for every already-seen trial, what the alias table and the vetting need, so
that rule variants can be scored without detrending and searching again. No rule is applied here.

A trial is prepared exactly as in round 2 (round2/harness.run_star, files frozen). Stored per trial:
the definition, the scalars of the reference, the events (with at_truth and the round-3 features),
the null-run events (arm B), and - only if some event has SNR >= 7 - the detrended T and F.

Round-3 event features (T5): sec_box_snr = event depth over the robust scatter (1.4826 MAD) of
same-duration box depths in the event's own sector (boxes every half duration, coverage >= 0.8).

Inputs  : trial definitions (CSV, as round2/harness.py writes them), MAST light curves (cache)
Outputs : <out>/<tic>.pkl, one per star (calibration/mono/work/, gitignored)
Usage   : python calibration/mono/round3/dump3.py --defs a.csv [b.csv ...] --out calibration/mono/work/round3_tuning
          [--noise 0.1 --noise-seed 20261205]   (10% Gaussian error on catalogued depth and duration, arm B)
"""
import argparse
import multiprocessing as mp
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "round2"))
import harness as h2  # noqa: E402  (frozen round-2 harness)

mc, rb, core, dump = h2.mc, h2.rb, h2.core, h2.dump
CACHE = "/Users/mcbee/Documents/GitHub/coracle/duo-transits/salinas-second-transits/work/lc_cache"


def sector_box_snr(det, ev, dur):
    t, f = det.get(ev.get("sector"), (None, None))
    if t is None or t.size < 50:
        return np.nan
    cen = np.arange(t.min() + dur / 2, t.max() - dur / 2, dur / 2)
    d, fr, _ = core.box_stats(t, f, cen, dur)
    d = d[(fr >= 0.8) & (np.abs(cen - ev["tc"]) > dur)]
    if d.size < 10:
        return np.nan
    mad = 1.4826 * np.median(np.abs(d - np.median(d)))
    return float(ev["depth"] / mad) if mad > 0 else np.nan


def prep_trial(d, sectors, tic):
    """As round2/harness.run_star, up to the prepared trial."""
    derr = d.get("depth_cat_err", np.nan)
    null_events = None
    if d["arm"] in ("A1", "A2"):
        sec = sectors
        if d["arm"] == "A1":
            Traw = np.concatenate([t for _, _, t, _ in sectors])
            _, ep = mc.epochs(d["truth_t0"], d["truth_P"], Traw.min(), Traw.max())
            ref = ep[np.argmin(np.abs(ep - d["t_cat"]))]
            others = ep[np.abs(ep - ref) > d["dur_cat_d"] / 2]
            sec = [(s, a, *mc.mask_windows(t, f, others, 1.5 * d["dur_cat_d"])) for s, a, t, f in sectors]
            sec = [x for x in sec if x[2].size >= 50]
    else:
        base = [(s, a, *mc.mask_windows(t, f, h2._events_of(tic), rb.EVENT_MASK_D)) for s, a, t, f in sectors]
        base = [x for x in base if x[2].size >= 50]
        p = dict(t0=d["inj_t0"], P=d["inj_P"], k=d["inj_k"], ar=d["inj_ar"], b=d["inj_b"])
        sec = rb.inject(base, p, (d["u1"], d["u2"]))
        null = mc.prepare(base, d["t_cat"], d["dur_cat_d"], d["depth_cat"])
        null_events = [] if null is None else [
            e for e in null["events"] if "error" not in e and abs(e["tc"] - d["t_cat"]) >= d["dur_cat_d"]]
    return mc.prepare(sec, d["t_cat"], d["dur_cat_d"], d["depth_cat"], derr), null_events


def star(args):
    defs, out, noise, seed = args
    tic = int(defs[0]["tic"])
    path = Path(out) / f"{tic}.pkl"
    if path.exists():
        return tic, "cached"
    try:
        sectors, _ = mc.load_star(tic, Path(CACHE))
        trials = []
        if sectors:
            rng = np.random.default_rng([seed, tic])
            for d in defs:
                d = dict(d)
                if noise and d["arm"] == "B":
                    d["depth_cat"] *= max(0.2, 1 + noise * rng.normal())
                    d["dur_cat_d"] *= max(0.2, 1 + noise * rng.normal())
                prep, null_events = prep_trial(d, sectors, tic)
                if prep is None:
                    continue
                s = dump.summarize(prep, d["t_cat"], d["truth_P"], d["truth_t0"])
                T, F, dur = prep["T"], prep["F"], prep["dur_ref"]
                ev = [e for e in s["events"] if "error" not in e]
                for e in ev:
                    e["sec_box_snr"] = sector_box_snr(prep["det"], e, dur)
                rs = s.get("ref_sector")
                tr, fr = prep["det"][rs]
                ref = core.refine_event(tr, fr, d["t_cat"], d["dur_cat_d"])
                keepTF = any(e.get("snr") is not None and e["snr"] >= 7 for e in ev)
                trials.append(dict(
                    d=d, events=ev, null_events=null_events, dur=dur, depth0=prep["depth_ref"],
                    derr=prep["derr_ref"], sig=s["sig"], ref_measured=s["ref_measured"], ref_snr=s["ref_snr"],
                    ref_tc=ref["tc"] if ref else np.nan, ref_dur=ref["dur_d"] if ref else np.nan,
                    tmin=float(T.min()), tmax=float(T.max()), n_clean_repeats=s["n_clean_repeats"],
                    pmin=core.pmin_from_duration(dur, d["mstar"], d["rstar"], prep["depth_ref"]),
                    pmin_cat=core.pmin_from_duration(d["dur_cat_d"], d["mstar"], d["rstar"], prep["depth_ref"]),
                    T=T if keepTF else None, F=F if keepTF else None))
        pickle.dump(trials, open(path, "wb"))
        return tic, f"{len(trials)} trials"
    except Exception as e:
        return tic, f"ERROR {type(e).__name__}: {str(e)[:200]}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--defs", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--noise", type=float, default=0.0)
    ap.add_argument("--noise-seed", type=int, default=20261205)
    a = ap.parse_args()
    Path(a.out).mkdir(parents=True, exist_ok=True)
    df = pd.concat([pd.read_csv(f).assign(src=Path(f).name) for f in a.defs])
    groups = [(g.to_dict("records"), a.out, a.noise, a.noise_seed) for _, g in df.groupby("tic")]
    with mp.get_context("spawn").Pool(a.workers) as pool:
        for i, (tic, msg) in enumerate(pool.imap_unordered(star, groups), 1):
            print(f"[{i}/{len(groups)}] {tic} {msg}", flush=True)


if __name__ == "__main__":
    main()
