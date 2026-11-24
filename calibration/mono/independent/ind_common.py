#!/usr/bin/env python3
"""Shared pieces of the BLIND independent re-derivation of the single-transit calibration.

Purpose
    Re-implements, from calibration/mono/prereg_mono_calibration.yaml (v1.1) only, the harness
    around the pipeline being calibrated (salinas-second-transits/st2core.py, st2_pipeline.py):
    loading a star, running the single-transit search on (possibly modified) raw sectors, the
    single-transit excluded set, and the ln P excluded fraction. Written without sight of
    calibration/mono/{monocal,run_arm_a,run_arm_b,runner,score}.py or any result of the first run.

Inputs
    MAST light curves through pl.lc_products / pl.load_sector_lc (cached in LC_CACHE),
    TIC parameters through pl.tic_info.

Outputs
    work/raw/<tic>.npz   raw normalised sectors of a star (so reruns do not touch MAST)
    work/tic/<tic>.json  TIC parameters

Usage
    imported by ind_arm_a.py and ind_arm_b.py; not run directly.
"""
import json
import sys
import time
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "salinas-second-transits"))
sys.path.insert(0, str(ROOT / "salinas-second-transits" / "round2"))
import st2core as core  # noqa: E402
import st2_pipeline as pl  # noqa: E402

warnings.filterwarnings("ignore")

WORK = HERE / "work"
LC_CACHE = Path("/Users/mcbee/Documents/GitHub/coracle/duo-transits/salinas-second-transits/work/lc_cache")
# the argparse defaults of pl.main
CUTS = SimpleNamespace(snr_min=7.0, depth_lo=0.5, depth_hi=2.0, dur_lo=0.5, dur_hi=2.0, anti_ratio=1.5)
P_FLOOR = 1.0
T_SCATTER_D = 15.0 / 1440.0
LNP_LO, LNP_HI = 1.0, 1000.0


def retry(fn, *a, tries=4, wait=8, **k):
    for i in range(tries):
        try:
            return fn(*a, **k)
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(wait * (i + 1))


def tic_params(tic):
    p = WORK / "tic" / f"{tic}.json"
    if p.exists():
        return json.loads(p.read_text())
    info = retry(pl.tic_info, tic)
    info = {k: (None if not np.isfinite(v) else float(v)) for k, v in info.items()}
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(info))
    return info


def load_star_raw(tic):
    """Raw normalised sectors (best product per sector, special sectors skipped).

    Returns (list of dict(sector, author, exptime, t, f), list of load errors).
    """
    p = WORK / "raw" / f"{tic}.npz"
    if p.exists():
        d = np.load(p, allow_pickle=False)
        meta = json.loads(str(d["meta"]))
        secs = [dict(sector=m["sector"], author=m["author"], exptime=m["exptime"],
                     t=d[f"t{m['sector']}"], f=d[f"f{m['sector']}"]) for m in meta["sectors"]]
        return secs, meta["errors"]
    sr, prods = retry(pl.lc_products, tic)
    best = pl.best_per_sector(prods)
    secs, errors = [], []
    for s, row in sorted(best.items()):
        if s in pl.SPECIAL_SECTORS:
            continue
        try:
            t, f = retry(pl.load_sector_lc, sr, row, LC_CACHE, tries=2)
        except Exception as e:
            errors.append(f"S{s}: {str(e)[:120]}")
            continue
        o = np.argsort(t)
        t, f = np.asarray(t[o], float), np.asarray(f[o], float)
        if t.size < 50:
            errors.append(f"S{s}: only {t.size} points")
            continue
        secs.append(dict(sector=int(s), author=row["author"], exptime=float(row["exptime"]), t=t, f=f))
    p.parent.mkdir(parents=True, exist_ok=True)
    arrs = {}
    for s in secs:
        arrs[f"t{s['sector']}"] = s["t"]
        arrs[f"f{s['sector']}"] = s["f"]
    meta = dict(sectors=[dict(sector=s["sector"], author=s["author"], exptime=s["exptime"]) for s in secs],
                errors=errors)
    np.savez_compressed(p, meta=json.dumps(meta), **arrs)
    return secs, errors


def concat(parts):
    T = np.concatenate([p[0] for p in parts])
    F = np.concatenate([p[1] for p in parts])
    o = np.argsort(T)
    return T[o], F[o]


def clean_mask(T, ep, dur_d, F=None):
    """Clean epochs (st2core's definition): coverage >= CLEAN_FRAC, >= CLEAN_EDGE_D from a gap/edge."""
    depth, frac, _ = core.box_stats(T, np.ones_like(T) if F is None else F, ep, dur_d)
    edge = core.edge_distance(core.gap_edges(T), ep)
    return (frac >= core.CLEAN_FRAC) & (edge >= core.CLEAN_EDGE_D) & np.isfinite(depth)


def epochs_in_span(t_ref, P, tmin, tmax):
    k = np.arange(int(np.floor((tmin - t_ref) / P)) - 1, int(np.ceil((tmax - t_ref) / P)) + 2)
    return k, t_ref + k * P


def lnp_excluded_fraction(t_cat, T, F, dur, depth, derr, sig, edges):
    """Fraction of ln P in [1, 1000] d that is in the single-transit excluded set.

    Spec: periods log-spaced with a step that moves the farthest epoch by a quarter duration;
    epochs looked up (nearest point) on a time grid of one eighth of a duration. Instead of looping
    over periods, each run of excluding grid times is mapped to the period-grid index ranges whose
    epoch k rounds into it (identical result, far fewer operations).
    """
    tmin, tmax = float(T.min()), float(T.max())
    dtg = dur / 8.0
    tg = np.arange(tmin, tmax + dtg, dtg)
    v = core.epoch_verdicts(tg, T, F, dur, depth, derr, sig, edges)
    ex = v["verdict"] == "excludes"
    t_far = max(tmax - t_cat, t_cat - tmin)
    dln = dur / (4.0 * t_far)
    n_p = int(np.floor(np.log(LNP_HI / LNP_LO) / dln)) + 1
    j = np.flatnonzero(ex)
    if j.size == 0:
        return 0.0, n_p, 0.0
    brk = np.flatnonzero(np.diff(j) > 1)
    starts = j[np.concatenate([[0], brk + 1])]
    ends = j[np.concatenate([brk, [j.size - 1]])]
    lo_all, hi_all = [], []
    for a, b in zip(tg[starts] - dtg / 2 - t_cat, tg[ends] + dtg / 2 - t_cat):
        pieces = [(a, b)] if a >= 0 else [(-b, -a)] if b <= 0 else [(0.0, b), (0.0, -a)]
        for lo, hi in pieces:
            if hi < LNP_LO:
                continue
            k = np.arange(max(1, int(np.ceil(lo / LNP_HI))), int(np.floor(hi / LNP_LO)) + 1)
            if not k.size:
                continue
            plo = np.maximum(lo / k, LNP_LO)
            phi = np.minimum(hi / k, LNP_HI)
            ilo = np.ceil(np.log(plo / LNP_LO) / dln).astype(np.int64)
            ihi = np.minimum(np.floor(np.log(phi / LNP_LO) / dln).astype(np.int64), n_p - 1)
            ok = ihi >= ilo
            lo_all.append(ilo[ok])
            hi_all.append(ihi[ok] + 1)
    if not lo_all:
        return 0.0, n_p, float(ex.mean())
    diff = (np.bincount(np.concatenate(lo_all), minlength=n_p + 1)
            - np.bincount(np.concatenate(hi_all), minlength=n_p + 1))
    excl = np.cumsum(diff[:n_p]) > 0
    return float(excl.mean()), n_p, float(ex.mean())


def lnp_excluded_fraction_brute(t_cat, T, F, dur, depth, derr, sig, edges, max_periods=40000):
    """Literal loop over (a thinned set of) grid periods; used only to test the fast version."""
    tmin, tmax = float(T.min()), float(T.max())
    dtg = dur / 8.0
    tg = np.arange(tmin, tmax + dtg, dtg)
    ex = core.epoch_verdicts(tg, T, F, dur, depth, derr, sig, edges)["verdict"] == "excludes"
    t_far = max(tmax - t_cat, t_cat - tmin)
    dln = dur / (4.0 * t_far)
    n_p = int(np.floor(np.log(LNP_HI / LNP_LO) / dln)) + 1
    idx = np.unique(np.linspace(0, n_p - 1, min(n_p, max_periods)).astype(int))
    out = np.zeros(idx.size, bool)
    for m, i in enumerate(idx):
        P = LNP_LO * np.exp(i * dln)
        k, ep = epochs_in_span(t_cat, P, tmin, tmax)
        ep = ep[(k != 0) & (ep >= tmin - dtg / 2) & (ep < tmax + dtg / 2)]
        jj = np.rint((ep - tg[0]) / dtg).astype(int)
        jj = jj[(jj >= 0) & (jj < tg.size)]
        out[m] = ex[jj].any()
    return float(out.mean())


def run_trial(secs, t_cat, depth_cat, dur_cat, depth_err_cat, P_true, t_ref_true, dur_true,
              mstar, rstar, want_lnp=True):
    """The pipeline's single-transit treatment of one (possibly modified) star.

    secs: raw sectors. t_cat, depth_cat (fraction), dur_cat (d), depth_err_cat: what a catalogue
    would hand the pipeline. P_true, t_ref_true, dur_true: the truth, used only for scoring.
    """
    out = dict(status="ok", t_cat=float(t_cat), truth_P=float(P_true), n_sectors=len(secs))
    # ---- reference event, as pl.search_one: the sector holding t_cat, detrended at the catalogue duration
    ref, ref_i, det = None, None, {}
    for i, s in enumerate(secs):
        if s["t"].min() - 0.5 <= t_cat <= s["t"].max() + 0.5:
            try:
                det[i] = pl.detrend(s["t"], s["f"], dur_cat)
            except Exception as e:
                out["ref_error"] = str(e)[:120]
                continue
            ref = core.refine_event(det[i][0], det[i][1], t_cat, dur_cat)
            ref_i = i
            break
    if ref and ref["snr"] >= 5 and 0.33 < ref["depth"] / depth_cat < 3:
        depth_ref, dur_ref, derr_ref, ref_used = ref["depth"], ref["dur_d"], ref["depth_err"], "box"
    else:
        depth_ref, dur_ref = depth_cat, dur_cat
        derr_ref = depth_err_cat if (depth_err_cat is not None and np.isfinite(depth_err_cat)) else 0.2 * depth_cat
        ref_used = "catalogue"
    out.update(ref_snr=float(ref["snr"]) if ref else np.nan, ref_box_depth_ppm=ref["depth"] * 1e6 if ref else np.nan,
               ref_box_dur_h=ref["dur_d"] * 24 if ref else np.nan, ref_used=ref_used,
               depth_ref_ppm=depth_ref * 1e6, dur_ref_h=dur_ref * 24, depth_ref_err_ppm=derr_ref * 1e6,
               ref_sector=secs[ref_i]["sector"] if ref_i is not None else None,
               ref_author=secs[ref_i]["author"] if ref_i is not None else None,
               ref_cadence_min=float(core.cadence_days(secs[ref_i]["t"]) * 1440) if ref_i is not None else np.nan,
               ref_exptime_s=secs[ref_i]["exptime"] if ref_i is not None else np.nan)

    # ---- detrend every other sector at the reference duration (pl.search_sector / all_lc_for_tic)
    parts, events, sector_errors = [], [], []
    for i, s in enumerate(secs):
        if i not in det:
            try:
                det[i] = pl.detrend(s["t"], s["f"], dur_ref)
            except Exception as e:
                sector_errors.append(f"S{s['sector']}: {str(e)[:80]}")
                continue
        t, f = det[i]
        if t.size < 20:
            continue
        parts.append((t, f))
        # ---- search; the reference sector too, minus one duration either side of the reference
        if i == ref_i:
            m = np.abs(t - t_cat) >= dur_ref
            ts, fs = t[m], f[m]
        else:
            ts, fs = t, f
        if ts.size < 20:
            continue
        try:
            evs = pl.sector_events(ts, fs, s["sector"], s["author"], dur_ref, depth_ref)
        except Exception as e:
            sector_errors.append(f"S{s['sector']} search: {str(e)[:80]}")
            continue
        for ev in evs:
            ev["passes"] = bool(pl.passes(ev, CUTS))
            ev["clean"] = bool(ev["edge_dist_d"] >= core.CLEAN_EDGE_D and ev["frac"] >= core.CLEAN_FRAC)
            events.append(ev)
    out["sector_errors"] = "; ".join(sector_errors)
    if not parts:
        out["status"] = "no usable data"
        return out
    T, F = concat(parts)
    sig = core.duration_scale_noise(T, F, dur_ref)
    edges = core.gap_edges(T)
    out.update(noise_dur_ppm=sig * 1e6, baseline_d=float(np.ptp(T)), n_points=int(T.size))

    # ---- truth bookkeeping: other true transits at clean epochs of the analysed data
    k, ep = epochs_in_span(t_ref_true, P_true, T.min(), T.max())
    other = ep[k != 0]
    cl = clean_mask(T, other, dur_true, F) if other.size else np.zeros(0, bool)
    clean_repeats = other[cl]
    out["n_clean_repeats"] = int(cl.sum())
    out["n_repeats_any_data"] = int((core.box_stats(T, F, other, dur_true)[1] > 0).sum()) if other.size else 0

    # ---- single-transit excluded set: the true period at its exact value
    kk, epc = epochs_in_span(t_cat, P_true, T.min(), T.max())
    epc = epc[kk != 0]
    v = core.epoch_verdicts(epc, T, F, dur_ref, depth_ref, derr_ref, sig, edges)
    exm = v["verdict"] == "excludes"
    out["truth_excluded"] = bool(exm.any())
    out["n_truth_epochs_clean"] = int(np.isin(v["verdict"], ["excludes", "supports", "shallow", "deeper"]).sum())
    out["n_truth_epochs_excluding"] = int(exm.sum())
    out["n_truth_epochs_supporting"] = int((v["verdict"] == "supports").sum())
    out["truth_excluding_epochs"] = ";".join(
        f"{e:.4f}:{d*1e6:.0f}+/-{er*1e6:.0f}ppm" for e, d, er in zip(epc[exm], v["depth"][exm], v["depth_err"][exm]))
    if want_lnp:
        fr, n_p, fr_t = lnp_excluded_fraction(t_cat, T, F, dur_ref, depth_ref, derr_ref, sig, edges)
        out.update(frac_lnP_excluded=fr, n_period_grid=n_p)
    else:
        out.update(frac_lnP_excluded=np.nan, n_period_grid=0)

    # ---- repeats: passing events, best clean passing event (pl.alias_events), alias table
    passing = [e for e in events if e["passes"]]
    out["n_events"] = len(events)
    out["n_passing_events"] = len(passing)
    out["n_clean_passing_events"] = sum(e["clean"] for e in passing)
    tol = 0.5 * max(dur_ref, dur_true)

    def is_true_repeat(tc):
        return bool(other.size and np.min(np.abs(other - tc)) < tol)

    out["repeat_detected"] = bool(any(is_true_repeat(e["tc"]) for e in passing))
    out["n_passing_true"] = sum(is_true_repeat(e["tc"]) for e in passing)
    out["n_passing_false"] = len(passing) - out["n_passing_true"]
    out["passing_events"] = ";".join(
        f"S{e['sector']}:{e['tc']:.4f}:snr{e['snr']:.1f}:{'clean' if e['clean'] else 'notclean'}:"
        f"{'true' if is_true_repeat(e['tc']) else 'false'}" for e in passing)
    pmin = core.pmin_from_duration(dur_ref, mstar, rstar, depth_ref)
    out["pmin_circ_d"] = float(pmin) if np.isfinite(pmin) else np.nan
    if not passing:
        out.update(call="no repeat", truth_in_alias_set=False, n_alias_ok=0)
        return out
    best = sorted(passing, key=lambda e: (e["clean"], e["snr"]), reverse=True)[0]
    t1, t2 = sorted([float(t_cat), float(best["tc"])])
    at = core.alias_table(t1, t2, dur_ref, depth_ref, best["depth_err"], T, F,
                          p_floor=P_FLOOR, pmin_density=pmin)
    status = {r["n"]: r["status"] for r in at}
    n_ok = sum(s in ("allowed", "supported") for s in status.values())
    dt = t2 - t1
    n = int(round(dt / P_true))
    on_alias = n >= 1 and abs(dt - n * P_true) < 0.5 * dur_ref
    truth_status = status.get(n, "not in table") if on_alias else "event off the true ephemeris"
    out.update(call="resolved" if n_ok == 1 else "narrowed" if n_ok > 1 else "all excluded",
               n_alias_ok=int(n_ok), n_aliases=len(at),
               n_alias_below_pmin=sum(s == "allowed_below_pmin" for s in status.values()),
               n_alias_supported=sum(s == "supported" for s in status.values()),
               truth_in_alias_set=bool(on_alias and n in status and status[n] != "excluded"),
               truth_alias_status=truth_status, truth_alias_n=n,
               alias_event_tc=float(best["tc"]), alias_event_snr=float(best["snr"]),
               alias_event_sector=int(best["sector"]), alias_event_clean=bool(best["clean"]),
               alias_event_is_true_repeat=is_true_repeat(best["tc"]),
               alias_event_depth_ratio=float(best["depth_ratio"]), alias_event_dur_ratio=float(best["dur_ratio"]))
    return out


def jsonable(o):
    if isinstance(o, dict):
        return {k: jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, (float, np.floating)):
        return None if not np.isfinite(o) else float(o)
    return o
