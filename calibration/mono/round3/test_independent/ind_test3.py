#!/usr/bin/env python3
"""Independent re-derivation of the ROUND-3 TEST-SET per-trial outcomes.

A copy of ../xcheck/independent/ind3.py (a blind re-implementation of the round-3 per-trial
scoring), re-read against round3/FROZEN_RULE3.md, round2/FROZEN_RULE.md, prereg_round3.yaml,
prereg_round2.yaml and the pipeline (st2core.py, st2_pipeline.py). Not read: the round-3 harness,
rule, scorer, drivers, tuning log or primary test outputs, nor those of rounds 1 and 2.

Changes from ind3.py (see AMBIGUITIES.md):
  1. truth_excluded (round 2's R1-R3 at the true period's epochs, reference depth error) and
     n_clean_repeats are restored from round2/test_independent/ind_test.py; the score needs them.
  2. The leftover "highest SNR among clean events, else highest SNR" and its t2_boxclean
     diagnostic are removed: the frozen text says only "the passing event with the highest SNR".
  3. Diagnostics for the recovery analysis (not part of the rule): what R4 alone would have
     detected (repeat_detected_r4, false_event_r4), why the true repeat was lost (true_ev_*),
     the product and cadence of the reference sector and of the clean repeat epochs.
  4. ROOT for this directory; the arm/tic/t_cat of a failed row are kept.

Input   a defs CSV (../test/arm_a_defs.csv, ../test/arm_b_defs.csv), gaia_dr4/candidates_frozen.csv,
        the shared light-curve cache
Output  --out rows CSV (one row per definition) and a .ckpt.jsonl beside it
Usage   python ind_test3.py --defs ../test/arm_a_defs.csv --out ind_a_rows.csv [--workers 12]
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
# R4 (round 2)
SNR_MIN, DR_LO, DR_HI, DUR_LO, DUR_HI, ANTI = 8.0, 1 / 1.7, 1.7, 0.5, 2.0, 1.5
CLEAN_FRAC, CLEAN_EDGE, COVERED = 0.8, 0.5, 0.5
# T5 (round 3)
SEC_BOX_SNR_MIN, SEC_BOX_MIN_N = 9.0, 10
# T1 (round 3)
SIG_DT_D = float(np.hypot(15.0, 15.0)) / 1440.0   # 21.2 min
H_CAP = 16


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
    parts, events, secdata, secmeta = [], [], {}, []
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
        secdata[int(s["sector"])] = (t, f)
        secmeta.append((float(t.min()), float(t.max()), int(s["sector"]), s["author"],
                        float(core.cadence_days(t)) * 1440.0))
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
                events=events, secdata=secdata, secmeta=secmeta, ref_i=ref_i,
                ref_meta=(secs[ref_i]["author"], float(core.cadence_days(secs[ref_i]["t"])) * 1440.0)
                if ref_i is not None else ("", np.nan), log=log)


def sec_box_snr(R, ev):
    """T5: event depth over the scatter of same-duration boxes in the event's own sector."""
    sd = R["secdata"].get(int(ev["sector"]))
    if sd is None:
        return np.nan
    t, f = sd
    dur = R["dur"]
    cen = np.arange(t.min() + dur / 2, t.max() - dur / 2, dur / 2)
    if cen.size == 0:
        return np.nan
    d, frac, _ = core.box_stats(t, f, cen, dur)
    keep = (frac >= 0.8) & (np.abs(cen - ev["tc"]) > dur)
    if keep.sum() < SEC_BOX_MIN_N:
        return np.nan
    d = d[keep]
    sc = 1.4826 * np.median(np.abs(d - np.median(d)))
    if not np.isfinite(sc) or sc <= 0:
        return np.nan
    return float(ev["depth"] / sc)


def r4_fail(ev):
    """Which R4 conditions an event fails (diagnostic)."""
    why = []
    if not ev["snr"] >= SNR_MIN:
        why.append("snr")
    if not DR_LO <= ev["depth_ratio"] <= DR_HI:
        why.append("depth")
    if not DUR_LO <= ev["dur_ratio"] <= DUR_HI:
        why.append("dur")
    if not ev["snr"] >= ANTI * max(ev["anti_snr"], 1e-9):
        why.append("anti")
    if not (ev["frac"] >= CLEAN_FRAC and ev["edge_dist_d"] >= CLEAN_EDGE):
        why.append("clean")
    return why


def vet_r4(ev):
    return bool(ev["snr"] >= SNR_MIN and DR_LO <= ev["depth_ratio"] <= DR_HI
                and DUR_LO <= ev["dur_ratio"] <= DUR_HI
                and ev["snr"] >= ANTI * max(ev["anti_snr"], 1e-9)
                and ev["frac"] >= CLEAN_FRAC and ev["edge_dist_d"] >= CLEAN_EDGE)


def vet3(R, ev):
    """Round 3: R4 and sec_box_snr >= 9 (undefined fails)."""
    if not vet_r4(ev):
        return False
    if "sec_box_snr" not in ev:
        ev["sec_box_snr"] = sec_box_snr(R, ev)
    s = ev["sec_box_snr"]
    return bool(np.isfinite(s) and s >= SEC_BOX_SNR_MIN)


def boxes(R, c, derr):
    depth, frac, _ = core.box_stats(R["T"], R["F"], c, R["dur"])
    err = R["sig"] / np.sqrt(np.maximum(frac, 0.001))
    thr = 0.5 * R["depth0"] - 3 * np.hypot(err, derr)
    return depth, frac, err, thr


def excludes(R, c, h, derr):
    """Exclusion at each predicted epoch c with its own window half-width h (steps of dur/8)."""
    c = np.atleast_1d(np.asarray(c, float))
    h = np.broadcast_to(np.asarray(h, np.int64), c.shape)
    out = np.zeros(c.size, bool)
    if not R["measured"] or c.size == 0:
        return out
    depth, frac, err, thr = boxes(R, c, derr)
    edge = core.edge_distance(R["edges"], c)
    with np.errstate(invalid="ignore"):
        cand = (frac >= CLEAN_FRAC) & (edge >= CLEAN_EDGE) & (depth >= -3 * err) & (depth < thr)
    idx = np.flatnonzero(cand)
    j = 1
    while idx.size and j <= int(h[idx].max()):
        for sign in (-1.0, 1.0):
            sel = idx[h[idx] >= j]
            if sel.size == 0:
                continue
            d, fr, _, th = boxes(R, c[sel] + sign * j * R["dur"] / 8.0, derr)
            with np.errstate(invalid="ignore"):
                bad = (fr >= COVERED) & ~(d < th)
            if bad.any():
                idx = np.setdiff1d(idx, sel[bad], assume_unique=True)
        j += 1
    out[idx] = True
    return out


def clean_at(R, c):
    _, frac, _ = core.box_stats(R["T"], R["F"], c, R["dur"])
    return (frac >= CLEAN_FRAC) & (core.edge_distance(R["edges"], c) >= CLEAN_EDGE)


def truth_excluded(R, t_cat, P, w):
    """Round 2 use (i): any epoch t_cat + k P_true (k != 0) inside the data span excludes."""
    tmin, tmax = R["T"].min(), R["T"].max()
    k = np.arange(int(np.floor((tmin - t_cat) / P)) - 1, int(np.ceil((tmax - t_cat) / P)) + 2)
    ep = t_cat + k * P
    ep = ep[(k != 0) & (ep >= tmin) & (ep <= tmax)]
    return bool(excludes(R, ep, w, R["derr"]).any())


def alias_table(R, t_cat, t_ev, w0, derr, mstar, rstar, round3=True):
    """Statuses of the aliases P_n = dT/n (n = 1 ..). round3=False: h_k = w0, no cap (round 2)."""
    T, dur = R["T"], R["dur"]
    tmin, tmax = float(T.min()), float(T.max())
    t1, t2 = sorted([float(t_cat), float(t_ev)])
    dT = t2 - t1
    eps, owner, kk, Ps = [], [], [], []
    for n in range(1, 3001):
        P = dT / n
        if P < 1.0:
            break
        k = np.arange(int(np.floor((tmin - t1) / P)) - 1, int(np.ceil((tmax - t1) / P)) + 2)
        k = k[(k != 0) & (k != n)]
        eps.append(t1 + k * P)
        kk.append(k)
        owner.append(np.full(k.size, n - 1))
        Ps.append(P)
    na = len(Ps)
    if na == 0:
        return dT, []
    eps, owner, kk = np.concatenate(eps), np.concatenate(owner), np.concatenate(kk)
    nn = owner + 1
    if round3:
        far = np.minimum(np.abs(kk), np.abs(kk - nn))
        h = w0 + np.ceil(2.0 * far * (SIG_DT_D / nn) / (dur / 8.0) - 1e-9).astype(np.int64)
        use = h <= H_CAP
    else:
        h = np.full(eps.size, w0, np.int64)
        use = np.ones(eps.size, bool)
    ex = np.zeros(eps.size, bool)
    ex[use] = excludes(R, eps[use], h[use], derr)
    depth, frac, err, _ = boxes(R, eps, derr)
    edge = core.edge_distance(R["edges"], eps)
    with np.errstate(invalid="ignore"):
        sup = use & (frac >= CLEAN_FRAC) & (edge >= CLEAN_EDGE) & (np.abs(depth - R["depth0"]) < 3 * np.hypot(err, derr))
    n_ex = np.bincount(owner[ex], minlength=na)
    n_su = np.bincount(owner[sup], minlength=na)
    pmin = core.pmin_from_duration(dur, mstar, rstar, R["depth0"])
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
    return dT, status


def alias_call(R, d, t_ev, w0, derr):
    dT, st = alias_table(R, d["t_cat"], t_ev, w0, derr, d["mstar"], d["rstar"], round3=True)
    _, st2 = alias_table(R, d["t_cat"], t_ev, w0, derr, d["mstar"], d["rstar"], round3=False)
    na = len(st)
    n_ok = sum(s in ("allowed", "supported") for s in st)
    P_true = d["truth_P"]
    n = int(round(dT / P_true))
    truth_in = bool(1 <= n <= na and abs(dT - n * P_true) < R["dur"] / 2 and st[n - 1] != "excluded")
    return dict(call="resolved" if n_ok == 1 else "narrowed" if n_ok > 1 else "all excluded",
                n_alias_ok=int(n_ok), n_alias=na, truth_in=truth_in,
                frac=(sum(s == "excluded" for s in st) / na) if na else np.nan,
                frac_r2=(sum(s == "excluded" for s in st2) / na) if na else np.nan,
                n_excl=sum(s == "excluded" for s in st), n_supp=sum(s == "supported" for s in st),
                n_below=sum(s == "allowed_below_pmin" for s in st),
                ok_n=" ".join(str(i + 1) for i, s in enumerate(st) if s in ("allowed", "supported"))[:400])


def true_epochs(R, d):
    """True transit epochs other than the reference, over the data span (one period of margin)."""
    P, t0 = d["truth_P"], d["truth_t0"]
    tmin, tmax = R["T"].min(), R["T"].max()
    k = np.arange(int(np.floor((tmin - t0) / P)) - 1, int(np.ceil((tmax - t0) / P)) + 2)
    ep = t0 + k * P
    return ep[np.abs(ep - d["t_cat"]) > R["dur"] / 2]


def score_events(R, d, w0, events, suffix, cache):
    other = true_epochs(R, d)

    def at_true(tc):
        return bool(other.size and np.min(np.abs(other - tc)) < R["dur"] / 2)

    passing = [e for e in events if vet3(R, e)]
    r4 = [e for e in events if vet_r4(e)]
    out = {"n_r4" + suffix: len(r4),
           "false_event_r4" + suffix: any(not at_true(e["tc"]) for e in r4),
           "repeat_detected_r4" + suffix: any(at_true(e["tc"]) for e in r4),
           "n_vetted" + suffix: len(passing),
           "false_event" + suffix: any(not at_true(e["tc"]) for e in passing),
           "repeat_detected" + suffix: any(at_true(e["tc"]) for e in passing)}
    if not passing:
        out.update({"t2" + suffix: np.nan, "call" + suffix: "no repeat", "n_alias_ok" + suffix: 0,
                    "n_alias" + suffix: 0, "truth_in_alias_set" + suffix: False,
                    "event_true" + suffix: False, "frac_alias_excluded" + suffix: np.nan,
                    "frac_alias_excluded_r2" + suffix: np.nan})
        return out
    best = max(passing, key=lambda e: e["snr"])
    key = float(best["tc"])
    if key not in cache:
        cache[key] = alias_call(R, d, key, w0, float(best["depth_err"]))
    a = cache[key]
    out.update({"t2" + suffix: key, "t2_snr" + suffix: float(best["snr"]),
                "t2_sec_box_snr" + suffix: float(best["sec_box_snr"]),
                "call" + suffix: a["call"], "n_alias_ok" + suffix: a["n_alias_ok"],
                "n_alias" + suffix: a["n_alias"], "truth_in_alias_set" + suffix: a["truth_in"],
                "event_true" + suffix: at_true(key),
                "frac_alias_excluded" + suffix: a["frac"], "frac_alias_excluded_r2" + suffix: a["frac_r2"],
                "n_alias_excluded" + suffix: a["n_excl"], "n_alias_supported" + suffix: a["n_supp"],
                "n_alias_below_pmin" + suffix: a["n_below"], "alias_ok_n" + suffix: a["ok_n"]})
    return out


def do_trial(d, raw, cat_times):
    arm = d["arm"]
    w0 = 8 if d["multi"] else 4
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
    other = true_epochs(R, d)
    cl = clean_at(R, other) if other.size else np.zeros(0, bool)
    prods = set()
    for ep in other[cl]:
        for lo, hi, sec, au, cad in R["secmeta"]:
            if lo <= ep <= hi:
                prods.add(f"{au}:{cad:.0f}")
    tev = [e for e in R["events"] if other.size and np.min(np.abs(other - e["tc"])) < R["dur"] / 2]
    tinfo = dict(n_true_ev=len(tev))
    if tev:
        for e in tev:
            if vet_r4(e) and "sec_box_snr" not in e:
                e["sec_box_snr"] = sec_box_snr(R, e)
        tr4 = [e for e in tev if vet_r4(e)]
        tb = max(tr4 or tev, key=lambda e: (e.get("sec_box_snr", -np.inf) if tr4 else e["snr"]))
        tinfo.update(true_ev_r4=bool(tr4), true_ev_snr=float(tb["snr"]),
                     true_ev_sbs=float(tb.get("sec_box_snr", np.nan)),
                     true_ev_fail="+".join(r4_fail(tb)), true_ev_author=tb["author"],
                     true_ev_depth_ratio=float(tb["depth_ratio"]), true_ev_dur_ratio=float(tb["dur_ratio"]),
                     true_ev_anti=float(tb["anti_snr"]))
    out = dict(_order=d["_order"], arm=arm, tic=d["tic"], t_cat=d["t_cat"],
               truth_excluded=truth_excluded(R, d["t_cat"], d["truth_P"], w0),
               n_clean_repeats=int(cl.sum()), repeat_products=" ".join(sorted(prods)),
               ref_author=R["ref_meta"][0], ref_cad_min=R["ref_meta"][1],
               ref_depth_err=R["derr"], sig=R["sig"], **tinfo,
               sectors=" ".join(f"{s['sector']}:{s['author']}" for s in secs),
               ref_snr=R["ref"]["snr"] if R["ref"] else np.nan, ref_measured=R["measured"],
               depth_ref=R["depth0"], dur_ref_d=R["dur"], n_events=len(R["events"]))
    cache = {}
    out.update(score_events(R, d, w0, R["events"], "", cache))
    notes = list(R["log"])
    if arm == "B":
        N = search(null, d["t_cat"], d["dur_cat_d"], d["depth_cat"], derr_cat)
        own = [e for e in N["events"] if abs(e["tc"] - d["t_cat"]) >= d["dur_cat_d"]]
        own_t = np.array([e["tc"] for e in own], float)

        def at_true(tc):
            return bool(other.size and np.min(np.abs(other - tc)) < R["dur"] / 2)

        kept = [e for e in R["events"]
                if not (own_t.size and np.min(np.abs(own_t - e["tc"])) < R["dur"] and not at_true(e["tc"]))]
        out.update(score_events(R, d, w0, kept, "_own", cache))
        out.update(n_null_events=len(own), null_ref_measured=N["measured"],
                   n_null_r4=sum(vet_r4(e) for e in own), n_null_vetted=sum(vet3(N, e) for e in own))
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


CALL = ["n_vetted", "false_event", "repeat_detected", "t2", "call", "n_alias_ok", "n_alias",
        "truth_in_alias_set", "event_true", "frac_alias_excluded", "frac_alias_excluded_r2"]
EXTRA = ["t2_snr", "t2_sec_box_snr", "n_r4", "false_event_r4", "repeat_detected_r4", "n_alias_excluded", "n_alias_supported", "n_alias_below_pmin",
         "alias_ok_n"]
DIAG = ["n_true_ev", "true_ev_r4", "true_ev_snr", "true_ev_sbs", "true_ev_fail", "true_ev_author",
        "true_ev_depth_ratio", "true_ev_dur_ratio", "true_ev_anti", "repeat_products", "ref_author",
        "ref_cad_min", "ref_depth_err", "sig", "n_null_r4", "n_null_vetted"]
COLS = (["tic", "arm", "toi", "draw", "ref_snr", "ref_measured", "truth_excluded", "n_clean_repeats"] + CALL + [c + "_own" for c in CALL]
        # diagnostics, not asked for
        + ["multi", "t_cat", "truth_P", "depth_ref", "dur_ref_d", "sectors", "n_events"] + EXTRA
        + [c + "_own" for c in EXTRA] + DIAG + ["n_null_events", "null_ref_measured", "notes"])


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
    df = sel[["_order", "arm", "tic", "toi", "draw", "multi", "t_cat", "truth_P"]].merge(
        pd.DataFrame(rows).drop(columns=["arm", "tic", "t_cat"], errors="ignore"),
        on="_order", how="left").sort_values("_order")
    for c in COLS:
        if c not in df:
            df[c] = np.nan
    for c in [c for c in COLS if c.startswith("n_")]:
        df[c] = pd.to_numeric(df[c]).astype("Int64")
    df[COLS].to_csv(a.out, index=False)
    print(f"wrote {len(df)} rows to {a.out}")


if __name__ == "__main__":
    main()
