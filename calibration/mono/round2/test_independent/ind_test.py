#!/usr/bin/env python3
"""Independent re-derivation of the round-2 TEST-SET per-trial outcomes.

Adapted from ../xcheck/independent/ind2.py (an earlier blind re-implementation), re-read against
FROZEN_RULE.md, prereg_round2.yaml (1.1) and the pipeline (st2core.py, st2_pipeline.py) only,
without sight of round2/harness.py, rule2.py, score2.py, tune.py, dump.py, run_test.py,
TUNING_LOG.md, tuning/, monocal.py, the round-1 drivers or the primary test outputs.

Changes from ind2.py (see AMBIGUITIES.md):
  1. `call`, `n_alias_ok`, `truth_in_alias_set` (and the _own versions) use the alias event's own
     depth error in the alias table, as FROZEN_RULE.md "Calls" says. ind2.py reported that as the
     `_everr` variant and the reference-error table as the primary; here the reference-error table
     is the `_referr` variant.
  2. --defs / --out arguments, up to 12 workers, rows keyed by the definition's row number,
     a per-star checkpoint (jsonl), the sector/author list of each star in `sectors`.

Input   a defs CSV (test/arm_a_defs.csv or test/arm_b_defs.csv), gaia_dr4/candidates_frozen.csv
Output  --out rows CSV (one row per definition)
Usage   python ind_test.py --defs ../test/arm_a_defs.csv --out ind_a_rows.csv [--workers 12]
"""
import json
import argparse
import multiprocessing as mp
import sys
import time
import traceback
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "salinas-second-transits"))
import st2core as core  # noqa: E402
import st2_pipeline as pl  # noqa: E402

warnings.filterwarnings("ignore")
LC_CACHE = Path("/Users/mcbee/Documents/GitHub/coracle/duo-transits/salinas-second-transits/work/lc_cache")
MIN_PTS = 50
MASK_D = 0.75
# R4
SNR_MIN, DR_LO, DR_HI, DUR_LO, DUR_HI, ANTI = 8.0, 1 / 1.7, 1.7, 0.5, 2.0, 1.5
CLEAN_FRAC, CLEAN_EDGE, COVERED = 0.8, 0.5, 0.5


def retry(fn, *a, tries=4, wait=8, **k):
    for i in range(tries):
        try:
            return fn(*a, **k)
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(wait * (i + 1))


def load_star(tic):
    sr, prods = retry(pl.lc_products, tic)
    best = pl.best_per_sector(prods)
    secs, errs = [], []
    for s in sorted(best):
        if s in pl.SPECIAL_SECTORS:
            continue
        row = best[s]
        try:
            t, f = pl.load_sector_lc(sr, row, LC_CACHE)
        except Exception as e:
            errs.append(f"S{s} load: {str(e)[:100]}")
            continue
        if t.size < MIN_PTS:
            continue
        secs.append(dict(sector=int(s), author=row["author"], t=np.asarray(t, float), f=np.asarray(f, float)))
    return secs, errs


def cut(secs, times, half):
    """Remove points within `half` of each time; drop sectors left with fewer than MIN_PTS."""
    times = np.sort(np.asarray(times, float))
    out = []
    for s in secs:
        t = s["t"]
        keep = np.ones(t.size, bool)
        for e in times[(times >= t.min() - half) & (times <= t.max() + half)]:
            keep &= np.abs(t - e) > half
        if keep.sum() >= MIN_PTS:
            out.append(dict(s, t=t[keep], f=s["f"][keep]))
    return out


def inject(secs, d):
    import batman
    pm = batman.TransitParams()
    pm.t0, pm.per, pm.rp, pm.a = d["inj_t0"], d["inj_P"], d["inj_k"], d["inj_ar"]
    pm.inc = float(np.degrees(np.arccos(d["inj_b"] / d["inj_ar"])))
    pm.ecc, pm.w, pm.u, pm.limb_dark = 0.0, 90.0, [d["u1"], d["u2"]], "quadratic"
    out = []
    for s in secs:
        cad = core.cadence_days(s["t"])
        ss = max(1, int(round(cad / (2.0 / 1440.0))))
        m = batman.TransitModel(pm, np.ascontiguousarray(s["t"]), supersample_factor=ss, exp_time=cad).light_curve(pm)
        out.append(dict(s, f=s["f"] * m))
    return out


def search(secs, t_cat, dur_cat, depth_cat, depth_cat_err):
    """Reference rule, detrending and event search; returns the trial's working set."""
    log = []
    det, ref, ref_i = {}, None, None
    for i, s in enumerate(secs):
        if s["t"].min() - 0.5 <= t_cat <= s["t"].max() + 0.5:
            try:
                det[i] = pl.detrend(s["t"], s["f"], dur_cat)
            except Exception as e:
                log.append(f"S{s['sector']} ref detrend: {str(e)[:80]}")
                det[i] = None
                break
            ref = core.refine_event(det[i][0], det[i][1], t_cat, dur_cat)
            ref_i = i
            break
    measured = bool(ref and ref["snr"] >= 5 and 0.33 < ref["depth"] / depth_cat < 3)
    if measured:
        depth0, dur, derr = ref["depth"], ref["dur_d"], ref["depth_err"]
    else:
        depth0, dur = depth_cat, dur_cat
        derr = depth_cat_err if (depth_cat_err is not None and np.isfinite(depth_cat_err)) else 0.2 * depth_cat
    parts, events = [], []
    for i, s in enumerate(secs):
        if i not in det:
            try:
                det[i] = pl.detrend(s["t"], s["f"], dur)
            except Exception as e:
                log.append(f"S{s['sector']} detrend: {str(e)[:80]}")
                continue
        if det[i] is None or det[i][0].size == 0:
            continue
        t, f = det[i]
        parts.append((t, f))
        try:
            evs = pl.sector_events(t, f, s["sector"], s["author"], dur, depth0)
        except Exception as e:
            log.append(f"S{s['sector']} events: {str(e)[:80]}")
            continue
        for ev in evs:
            if i == ref_i and abs(ev["tc"] - t_cat) < dur:
                continue
            events.append(ev)
    if not parts:
        raise RuntimeError("no detrended data")
    T = np.concatenate([p[0] for p in parts])
    F = np.concatenate([p[1] for p in parts])
    o = np.argsort(T, kind="stable")
    T, F = T[o], F[o]
    return dict(ref=ref, measured=measured, depth0=float(depth0), dur=float(dur), derr=float(derr),
                T=T, F=F, sig=float(core.duration_scale_noise(T, F, dur)), edges=core.gap_edges(T),
                events=events, log=log)


def vet(ev):
    return bool(ev["snr"] >= SNR_MIN and DR_LO <= ev["depth_ratio"] <= DR_HI
                and DUR_LO <= ev["dur_ratio"] <= DUR_HI
                and ev["snr"] >= ANTI * max(ev["anti_snr"], 1e-9)
                and ev["frac"] >= CLEAN_FRAC and ev["edge_dist_d"] >= CLEAN_EDGE)


def boxes(R, c, derr):
    depth, frac, _ = core.box_stats(R["T"], R["F"], c, R["dur"])
    err = R["sig"] / np.sqrt(np.maximum(frac, 0.001))
    thr = 0.5 * R["depth0"] - 3 * np.hypot(err, derr)
    return depth, frac, err, thr


def clean_at(R, c):
    _, frac, _ = core.box_stats(R["T"], R["F"], c, R["dur"])
    return (frac >= CLEAN_FRAC) & (core.edge_distance(R["edges"], c) >= CLEAN_EDGE)


def excludes(R, c, w, derr=None):
    """R1-R3 at each predicted epoch c."""
    c = np.atleast_1d(np.asarray(c, float))
    out = np.zeros(c.size, bool)
    if not R["measured"] or c.size == 0:
        return out
    derr = R["derr"] if derr is None else derr
    depth, frac, err, thr = boxes(R, c, derr)
    edge = core.edge_distance(R["edges"], c)
    with np.errstate(invalid="ignore"):
        cand = (frac >= CLEAN_FRAC) & (edge >= CLEAN_EDGE) & (depth >= -3 * err) & (depth < thr)
    idx = np.flatnonzero(cand)
    for j in range(-w, w + 1):
        if idx.size == 0:
            break
        if j == 0:
            continue
        d, fr, _, th = boxes(R, c[idx] + j * R["dur"] / 8.0, derr)
        with np.errstate(invalid="ignore"):
            bad = (fr >= COVERED) & ~(d < th)
        idx = idx[~bad]
    out[idx] = True
    return out


def truth_excluded(R, t_cat, P, w):
    tmin, tmax = R["T"].min(), R["T"].max()
    k = np.arange(int(np.floor((tmin - t_cat) / P)) - 1, int(np.ceil((tmax - t_cat) / P)) + 2)
    ep = t_cat + k * P
    ep = ep[(k != 0) & (ep >= tmin) & (ep <= tmax)]
    return bool(excludes(R, ep, w).any())


def lnp_fraction(R, t_cat, w):
    T, dur = R["T"], R["dur"]
    tmin, tmax = float(T.min()), float(T.max())
    dln = dur / (4.0 * max(tmax - tmin, dur))
    n_p = int(np.floor(np.log(1000.0) / dln)) + 1
    if not R["measured"]:
        return 0.0
    dtg = dur / 8.0
    ng = int(np.floor((tmax - tmin) / dtg)) + 1
    tg = tmin + np.arange(ng + 1) * dtg
    ex = excludes(R, tg, w)
    if not ex.any():
        return 0.0
    P = np.exp(np.arange(n_p) * dln)
    out = np.zeros(n_p, bool)
    for sign in (1.0, -1.0):
        room = (tmax - t_cat) if sign > 0 else (t_cat - tmin)
        k = 1
        while k <= room + dtg:
            m = int(np.searchsorted(P, (room + dtg) / k, side="right"))
            if m == 0:
                break
            idx = np.rint((t_cat + sign * k * P[:m] - tmin) / dtg).astype(np.int64)
            ok = (idx >= 0) & (idx <= ng)
            hit = np.zeros(m, bool)
            hit[ok] = ex[idx[ok]]
            out[:m] |= hit
            k += 1
    return float(out.mean())


def alias_call(R, t_cat, t_ev, w, P_true, mstar, rstar, derr=None):
    derr = R["derr"] if derr is None else derr
    T = R["T"]
    tmin, tmax = float(T.min()), float(T.max())
    t1, t2 = sorted([float(t_cat), float(t_ev)])
    dT = t2 - t1
    eps, owner, Ps = [], [], []
    for n in range(1, 3001):
        P = dT / n
        if P < 1.0:
            break
        k = np.arange(int(np.floor((tmin - t1) / P)) - 1, int(np.ceil((tmax - t1) / P)) + 2)
        k = k[(k != 0) & (k != n)]
        eps.append(t1 + k * P)
        owner.append(np.full(k.size, n - 1))
        Ps.append(P)
    na = len(Ps)
    if na == 0:
        return dict(call="all excluded", n_alias_ok=0, truth_in=False, status=[])
    eps, owner = np.concatenate(eps), np.concatenate(owner)
    ex = excludes(R, eps, w, derr)
    depth, frac, err, _ = boxes(R, eps, derr)
    edge = core.edge_distance(R["edges"], eps)
    with np.errstate(invalid="ignore"):
        sup = (frac >= CLEAN_FRAC) & (edge >= CLEAN_EDGE) & (np.abs(depth - R["depth0"]) < 3 * np.hypot(err, derr))
    n_ex = np.bincount(owner[ex], minlength=na)
    n_su = np.bincount(owner[sup], minlength=na)
    pmin = core.pmin_from_duration(R["dur"], mstar, rstar, R["depth0"])
    status = []
    for i in range(na):
        if n_ex[i]:
            status.append("excluded")
        elif n_su[i]:
            status.append("supported")
        elif np.isfinite(pmin) and Ps[i] < pmin:
            status.append("allowed_below_pmin")
        else:
            status.append("allowed")
    n_ok = sum(s in ("allowed", "supported") for s in status)
    n = int(round(dT / P_true))
    truth_in = bool(1 <= n <= na and abs(dT - n * P_true) < R["dur"] / 2 and status[n - 1] != "excluded")
    return dict(call="resolved" if n_ok == 1 else "narrowed" if n_ok > 1 else "all excluded",
                n_alias_ok=int(n_ok), truth_in=truth_in, status=status)


def true_epochs(R, d):
    """True transit epochs other than the reference, over the data span (one period of margin)."""
    P, t0 = d["truth_P"], d["truth_t0"]
    tmin, tmax = R["T"].min(), R["T"].max()
    k = np.arange(int(np.floor((tmin - t0) / P)) - 1, int(np.ceil((tmax - t0) / P)) + 2)
    ep = t0 + k * P
    return ep[np.abs(ep - d["t_cat"]) > R["dur"] / 2]


def score_events(R, d, w, events, suffix, cache):
    other = true_epochs(R, d)

    def at_true(tc):
        return bool(other.size and np.min(np.abs(other - tc)) < R["dur"] / 2)

    passing = [e for e in events if vet(e)]
    out = {"n_vetted" + suffix: len(passing),
           "false_event" + suffix: any(not at_true(e["tc"]) for e in passing),
           "repeat_detected" + suffix: any(at_true(e["tc"]) for e in passing)}
    if not passing:
        out.update({"t2" + suffix: np.nan, "call" + suffix: "no repeat", "n_alias_ok" + suffix: 0,
                    "truth_in_alias_set" + suffix: False})
        out.update({"call_referr" + suffix: "no repeat", "n_alias_ok_referr" + suffix: 0,
                    "truth_in_alias_set_referr" + suffix: False})
        return out
    best = max(passing, key=lambda e: e["snr"])
    key = float(best["tc"])
    if key not in cache:
        cache[key] = (alias_call(R, d["t_cat"], key, w, d["truth_P"], d["mstar"], d["rstar"]),
                      alias_call(R, d["t_cat"], key, w, d["truth_P"], d["mstar"], d["rstar"], derr=best["depth_err"]))
    a, b = cache[key]   # a: reference depth error (variant); b: the alias event's own (FROZEN_RULE)
    out.update({"t2" + suffix: key, "t2_snr" + suffix: float(best["snr"]),
                "call" + suffix: b["call"], "n_alias_ok" + suffix: b["n_alias_ok"],
                "truth_in_alias_set" + suffix: b["truth_in"],
                "call_referr" + suffix: a["call"], "n_alias_ok_referr" + suffix: a["n_alias_ok"],
                "truth_in_alias_set_referr" + suffix: a["truth_in"]})
    return out


def do_trial(d, raw, cat_times):
    arm = d["arm"]
    w = 8 if d["multi"] else 4
    derr_cat = d["depth_cat_err"] if np.isfinite(d["depth_cat_err"]) else None
    null = None
    if arm == "A2":
        secs = raw
    elif arm == "A1":
        P, t0, dc = d["truth_P"], d["truth_t0"], d["dur_cat_d"]
        lo = min(s["t"].min() for s in raw)
        hi = max(s["t"].max() for s in raw)
        k = np.arange(int(np.floor((lo - t0) / P)) - 1, int(np.ceil((hi - t0) / P)) + 2)
        ep = t0 + k * P
        i_ref = int(np.argmin(np.abs(ep - d["t_cat"])))
        other = np.delete(ep, i_ref)
        other = other[np.abs(other - ep[i_ref]) > dc / 2]
        secs = cut(raw, other, 1.5 * dc)
    else:
        null = cut(raw, cat_times, MASK_D)
        secs = inject(null, d)
    R = search(secs, d["t_cat"], d["dur_cat_d"], d["depth_cat"], derr_cat)
    out = dict(_order=d["_order"], arm=arm, tic=d["tic"], t_cat=d["t_cat"],
               sectors=" ".join(f"{s['sector']}:{s['author']}" for s in secs),
               ref_snr=R["ref"]["snr"] if R["ref"] else np.nan, ref_measured=R["measured"],
               depth_ref=R["depth0"], dur_ref_d=R["dur"],
               truth_excluded=truth_excluded(R, d["t_cat"], d["truth_P"], w),
               frac_lnP_excluded=lnp_fraction(R, d["t_cat"], w))
    cache = {}
    out.update(score_events(R, d, w, R["events"], "", cache))
    other = true_epochs(R, d)
    out["n_clean_repeats"] = int(clean_at(R, other).sum()) if other.size else 0
    notes = list(R["log"])
    if arm == "B":
        N = search(null, d["t_cat"], d["dur_cat_d"], d["depth_cat"], derr_cat)
        own = [e for e in N["events"] if abs(e["tc"] - d["t_cat"]) >= d["dur_cat_d"]]
        own_t = np.array([e["tc"] for e in own], float)

        def at_true(tc):
            return bool(other.size and np.min(np.abs(other - tc)) < R["dur"] / 2)

        kept = [e for e in R["events"]
                if not (own_t.size and np.min(np.abs(own_t - e["tc"])) < R["dur"] and not at_true(e["tc"]))]
        out.update(score_events(R, d, w, kept, "_own", cache))
        nv = sum(vet(e) for e in own)
        out.update(null_event=bool(nv > 0), n_null_events=len(own), n_null_vetted=int(nv),
                   null_ref_measured=N["measured"])
        notes += ["null " + x for x in N["log"]]
    out["notes"] = "; ".join(notes)
    return out


def do_star(args):
    tic, defs, cat_times = args
    t0 = time.time()
    rows = []
    try:
        raw, errs = load_star(tic)
        if not raw:
            raise RuntimeError("no light curve: " + "; ".join(errs))
    except Exception as e:
        return tic, [dict(_order=d["_order"], arm=d["arm"], tic=tic, t_cat=d["t_cat"], notes=f"FAILED load: {type(e).__name__}: {str(e)[:200]}")
                     for d in defs], time.time() - t0
    for d in defs:
        try:
            r = do_trial(d, raw, cat_times)
            if errs:
                r["notes"] = "; ".join(filter(None, [r["notes"]] + errs))
        except Exception as e:
            r = dict(_order=d["_order"], arm=d["arm"], tic=tic, t_cat=d["t_cat"],
                     notes=f"FAILED: {type(e).__name__}: {str(e)[:200]} | {traceback.format_exc()[-600:]}")
        rows.append(r)
    return tic, rows, time.time() - t0


COLS = ["arm", "tic", "toi", "draw", "multi", "t_cat", "truth_P", "sectors", "ref_snr", "ref_measured",
        "depth_ref", "dur_ref_d", "truth_excluded", "frac_lnP_excluded", "n_vetted", "false_event",
        "repeat_detected", "t2", "t2_snr", "call", "n_alias_ok", "truth_in_alias_set", "n_clean_repeats",
        "n_vetted_own", "false_event_own", "repeat_detected_own", "t2_own", "t2_snr_own", "call_own",
        "n_alias_ok_own", "truth_in_alias_set_own", "null_event", "n_null_events", "n_null_vetted",
        "null_ref_measured",
        # variant: alias table with the reference depth error instead of the event's
        "call_referr", "n_alias_ok_referr", "truth_in_alias_set_referr",
        "call_referr_own", "n_alias_ok_referr_own", "truth_in_alias_set_referr_own", "notes"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--defs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--tics", default="")
    a = ap.parse_args()
    defs = pd.read_csv(a.defs)
    defs["multi"] = defs["multi"].astype(str).str.strip().str.lower() == "true"
    defs["_order"] = np.arange(len(defs))
    for c in ("toi", "draw"):
        if c not in defs:
            defs[c] = np.nan
    sel = defs
    if a.tics:
        sel = defs[defs["tic"].isin([int(x) for x in a.tics.split(",")])]
    ck = Path(a.out).with_suffix(".ckpt.jsonl")
    rows = []
    if ck.exists():
        rows = [json.loads(x) for x in ck.read_text().splitlines() if x.strip()]
    done = {r["_order"] for r in rows}
    cat = pd.read_csv(ROOT / "gaia_dr4" / "candidates_frozen.csv")
    jobs = []
    for tic, g in sel[~sel["_order"].isin(done)].groupby("tic", sort=False):
        ct = cat.loc[cat["tic"] == tic, ["t1_btjd", "t2_btjd"]].to_numpy(float).ravel()
        jobs.append((int(tic), g.to_dict("records"), ct[np.isfinite(ct)]))
    print(f"{len(done)} trials from checkpoint, {len(jobs)} stars to run", flush=True)
    with mp.get_context("spawn").Pool(min(a.workers, 12)) as pool, open(ck, "a") as fh:
        for i, (tic, rr, dt) in enumerate(pool.imap_unordered(do_star, jobs), 1):
            print(f"[{i}/{len(jobs)}] TIC {tic}: {len(rr)} trials ({dt:.0f} s)", flush=True)
            for r in rr:
                fh.write(json.dumps(r, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n")
            fh.flush()
            rows += rr
    df = pd.DataFrame(rows).drop(columns=["arm", "tic", "t_cat"], errors="ignore")
    df = sel[["_order", "arm", "tic", "toi", "draw", "multi", "t_cat", "truth_P"]].merge(df, on="_order", how="left")
    df = df.sort_values("_order")
    for c in COLS:
        if c not in df:
            df[c] = np.nan
    for c in [c for c in COLS if c.startswith("n_")]:
        df[c] = pd.to_numeric(df[c]).astype("Int64")
    df[COLS].to_csv(a.out, index=False)
    print(f"wrote {len(df)} rows to {a.out}")


if __name__ == "__main__":
    main()
