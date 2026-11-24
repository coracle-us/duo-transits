"""
Round 2 tuning, step 1: run the tuning-set stars once and store, for every trial, what a rule
variant needs (box depths on a time grid, boxes around the true period's epochs, every event).

Tuning set (prereg_round2.yaml, tuning.data): round 1's arm A stars (A1 and A2) and arm B hosts
(the 4 round-1 draws, 8 new low-depth draws per host, and a search with nothing injected for each
draw). No star outside round 1's samples is loaded.

Inputs  : calibration/mono/arm_a_sample.csv, arm_b_hosts.csv, gaia_dr4/candidates_frozen.csv,
          MAST light curves (cached)
Outputs : calibration/mono/work/round2_tuning/arm_a/<tic>.pkl, arm_b/<tic>.pkl (not in git)
Usage   : python calibration/mono/round2/dump.py --arm a|b [--workers 12] [--limit N]
"""
import argparse
import pickle
import sys
import multiprocessing as mp
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import monocal as mc  # noqa: E402
import run_arm_a as ra  # noqa: E402
import run_arm_b as rb  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
core, pl = mc.core, mc.pl
STEP_FRAC = 0.125            # grid step in durations
M_OFF = 8                    # boxes at -1 to +1 duration around each true epoch
EXTRA_SEED = 20261124        # the 8 new draws per arm-B host
EXTRA_DRAWS = 8
EXTRA_DEPTH = (500e-6, 5000e-6)


def summarize(prep, t_cat, truth_P=None, truth_t0=None, grid=True):
    """Arrays and events of one trial; no rule is applied here."""
    T, F, dur = prep["T"], prep["F"], prep["dur_ref"]
    sig = core.duration_scale_noise(T, F, dur)
    edges = core.gap_edges(T)
    step = dur * STEP_FRAC
    s = dict(prep["out"], t_cat=t_cat, derr_ref=prep["derr_ref"], sig=sig, step=step,
             ref_cadence_min=prep["ref_cadence_d"] * 1440, tmin=float(T.min()), tmax=float(T.max()))
    if grid:
        centers = np.arange(T.min(), T.max() + step, step)
        depth, frac, _ = core.box_stats(T, F, centers, dur)
        keep = frac >= 0.5
        s.update(c0=float(centers[0]), n_grid=int(centers.size), g_idx=np.flatnonzero(keep).astype(np.int32),
                 g_depth=depth[keep], g_frac=frac[keep].astype(np.float32),
                 g_clean=(frac[keep] >= core.CLEAN_FRAC) & (core.edge_distance(edges, centers[keep]) >= core.CLEAN_EDGE_D))
    events = []
    for e in prep["events"]:
        e = dict(e)
        if "error" not in e:
            td, fd = prep["det"][e["sector"]]
            cen, _, snr, _, _ = core.box_search(td, fd, dur)
            near = (np.abs(cen - e["tc"]) < 1.0) & (np.abs(cen - e["tc"]) > dur)
            v = snr[near]
            e["local_anti_snr"] = float(-np.nanmin(v)) if np.isfinite(v).any() else 0.0
        events.append(e)
    if truth_P is not None:
        t0 = t_cat if truth_t0 is None else truth_t0
        _, ep = mc.epochs(t0, truth_P, T.min(), T.max())
        other = ep[np.abs(ep - t_cat) > dur / 2]
        cl = mc.clean_mask(T, ep, dur) & (np.abs(ep - t_cat) > dur / 2)
        s["n_clean_repeats"] = int(cl.sum())
        for e in events:
            if "error" not in e:
                e["at_truth"] = bool(other.size and np.min(np.abs(other - e["tc"])) < dur / 2)
        k1, ep1 = mc.epochs(t_cat, truth_P, T.min(), T.max())
        ep1 = ep1[k1 != 0]
        off = np.arange(-M_OFF, M_OFF + 1) * step
        c = ep1[:, None] + off[None, :]
        d, fr, _ = core.box_stats(T, F, c.ravel(), dur)
        s.update(truth_P=truth_P, tr_depth=d.reshape(c.shape), tr_frac=fr.reshape(c.shape),
                 tr_edge=core.edge_distance(edges, c.ravel()).reshape(c.shape))
    s["events"] = events
    return s


def star_a(item, cache):
    """Arm A exactly as run_arm_a.run_star chooses the reference; summaries instead of scores."""
    tic = int(item["tic"])
    rng = np.random.default_rng([ra.SEED, tic])
    sectors, _ = mc.load_star(tic, Path(cache))
    info = pl.tic_info(tic)
    if not sectors:
        return dict(tic=tic, skipped="no light curve", trials=[])
    dur, depth, P0 = item["duration_h"] / 24, item["depth_ppm"] * 1e-6, item["period_d"]
    T, F = [], []
    for _, _, t, f in sectors:
        try:
            td, fd = pl.detrend(t, f, dur)
            T.append(td)
            F.append(fd)
        except Exception:
            pass
    T, F = np.concatenate(T), np.concatenate(F)
    o = np.argsort(T)
    perr = item.get("period_err")
    half = max(5 * (perr if perr is not None and np.isfinite(perr) else 0.0), 2e-4 * P0)
    eph = ra._s10().refine_ephemeris(T[o], F[o], P0, item["t0"], dur, half)
    P, t0 = eph["P"], eph["t0"]
    Traw = np.sort(np.concatenate([t for _, _, t, _ in sectors]))
    _, ep = mc.epochs(t0, P, Traw.min(), Traw.max())
    cl = mc.clean_mask(Traw, ep, dur)
    if cl.sum() < ra.MIN_CLEAN:
        return dict(tic=tic, skipped="under 3 clean transits", trials=[])
    ref = float(ep[cl][0])
    t_cat = ref + rng.normal(0, ra.T_CAT_ERR_D)
    others = ep[np.abs(ep - ref) > dur / 2]
    hidden = [(s, a, *mc.mask_windows(t, f, others, 1.5 * dur)) for s, a, t, f in sectors]
    hidden = [x for x in hidden if x[2].size >= 50]
    derr = item.get("depth_err")
    derr = derr * 1e-6 if derr is not None and np.isfinite(derr) else np.nan
    mstar, rstar = info.get("mstar", np.nan), info.get("rstar", np.nan)
    trials = []
    for arm, sec in (("A1", hidden), ("A2", sectors)):
        prep = mc.prepare(sec, t_cat, dur, depth, derr)
        if prep is None:
            continue
        s = summarize(prep, t_cat, P, t0)
        s.update(arm=arm, tic=tic, toi=item["toi"], n_toi_on_star=int(item["n_toi_on_star"]), tmag=item["tmag"],
                 depth_cat=depth, pmin=core.pmin_from_duration(prep["dur_ref"], mstar, rstar, prep["depth_ref"]))
        trials.append(s)
    return dict(tic=tic, skipped="", trials=trials)


def star_b(item, cache):
    """Arm B: the 4 round-1 draws (same seed, so the same planets), 8 new shallower draws, and for
    every draw a search of the star with nothing injected, using that draw's depth and duration."""
    tic = int(item["tic"])
    events = np.asarray(item.get("events", []), float)
    sectors, _ = mc.load_star(tic, Path(cache))
    info = pl.tic_info(tic)
    if not sectors:
        return dict(tic=tic, skipped="no light curve", trials=[])
    sectors = [(s, a, *mc.mask_windows(t, f, events, rb.EVENT_MASK_D)) for s, a, t, f in sectors]
    sectors = [x for x in sectors if x[2].size >= 50]
    mstar, rstar = info.get("mstar", np.nan), info.get("rstar", np.nan)
    u1, u2, _ = mc.limb_dark(info.get("teff", np.nan), mstar, rstar)
    if not (np.isfinite(mstar) and np.isfinite(rstar) and mstar > 0 and rstar > 0):
        mstar, rstar = 1.0, 1.0
    trials = []
    plan = [("B", np.random.default_rng([rb.SEED, tic]), rb.DRAWS, rb.DEPTH_RANGE),
            ("Bx", np.random.default_rng([EXTRA_SEED, tic]), EXTRA_DRAWS, EXTRA_DEPTH)]
    for arm, rng, n, drange in plan:
        old = rb.DEPTH_RANGE
        rb.DEPTH_RANGE = drange
        try:
            for d in range(n):
                p = rb.draw(rng, sectors, mstar, rstar)
                if p is None:
                    continue
                t_cat = p["t_ref"] + rng.normal(0, rb.T_CAT_ERR_D)
                prep = mc.prepare(rb.inject(sectors, p, (u1, u2)), t_cat, p["dur_d"], p["depth"])
                if prep is None:
                    continue
                s = summarize(prep, t_cat, p["P"], p["t0"])
                s.update(arm=arm, tic=tic, draw=d, n_toi_on_star=1, tmag=info.get("tmag", np.nan),
                         depth_cat=p["depth"], inj_b=p["b"],
                         pmin=core.pmin_from_duration(prep["dur_ref"], mstar, rstar, prep["depth_ref"]))
                null = mc.prepare(sectors, t_cat, p["dur_d"], p["depth"])
                if null is not None:
                    ns = summarize(null, t_cat, grid=False)
                    s["null_events"] = [e for e in ns["events"]
                                        if "error" not in e and abs(e["tc"] - t_cat) >= p["dur_d"]]
                trials.append(s)
        finally:
            rb.DEPTH_RANGE = old
    return dict(tic=tic, skipped="", trials=trials)


def _call(args):
    fn, item, cache, path = args
    try:
        r = fn(item, cache)
    except Exception as e:
        r = dict(tic=int(item["tic"]), error=f"{type(e).__name__}: {str(e)[:300]}", trials=[])
    Path(path).write_bytes(pickle.dumps(r))
    return r["tic"], len(r["trials"]), r.get("error", r.get("skipped", ""))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm", choices=["a", "b"], required=True)
    ap.add_argument("--cache", default=str(ROOT / "calibration/mono/work/lc_cache"))
    ap.add_argument("--out", default=str(ROOT / "calibration/mono/work/round2_tuning"))
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--limit", type=int)
    a = ap.parse_args()
    if a.arm == "a":
        items = pd.read_csv(ROOT / "calibration/mono/arm_a_sample.csv").to_dict("records")
        fn = star_a
    else:
        cand = pd.read_csv(ROOT / "gaia_dr4/candidates_frozen.csv")
        items = []
        for tic in pd.read_csv(ROOT / "calibration/mono/arm_b_hosts.csv").tic.astype(int):
            c = cand[cand.tic == tic]
            items.append(dict(tic=tic, events=pd.concat([c.t1_btjd, c.t2_btjd]).dropna().tolist()))
        fn = star_b
    items = items[:a.limit] if a.limit else items
    out = Path(a.out) / f"arm_{a.arm}"
    out.mkdir(parents=True, exist_ok=True)
    todo = [(fn, it, a.cache, str(out / f"{int(it['tic'])}.pkl")) for it in items
            if not (out / f"{int(it['tic'])}.pkl").exists()]
    print(f"{len(items) - len(todo)} done already, {len(todo)} to run", flush=True)
    with mp.get_context("spawn").Pool(a.workers) as pool:
        for i, (tic, n, note) in enumerate(pool.imap_unordered(_call, todo), 1):
            print(f"[{i}/{len(todo)}] TIC {tic} trials={n} {note}", flush=True)


if __name__ == "__main__":
    main()
