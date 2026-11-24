#!/usr/bin/env python3
"""Arm B (injections) of the blind independent re-derivation.

Purpose
    For each of the 150 hosts in arm_b_hosts.csv: mask the catalogued events (+-0.75 d around every
    t1_btjd / t2_btjd in gaia_dr4/candidates_frozen.csv), then 4 draws of a batman transit train
    (P log-uniform 20-1000 d, (Rp/Rs)^2 log-uniform 500-20000 ppm, b uniform 0-0.9, circular,
    duration from P, TIC density and b, quadratic limb darkening from the Claret 2017 table at the
    nearest grid point) multiplied into the RAW normalised sectors, i.e. before detrending. The
    phase is redrawn until a transit falls in clean data; the first such transit is the reference
    and its time, scattered by 15 min, is the catalogued time handed to the pipeline.

Inputs
    ../arm_b_hosts.csv, ../data/claret2017_table25_Z0_xi2.csv, gaia_dr4/candidates_frozen.csv,
    MAST light curves and TIC (see ind_common.py).

Outputs
    work/arm_b/<tic>.json  one file per host: 4 trials, or one 'skipped' record.

Usage
    python ind_arm_b.py [--workers 4] [--limit N] [--tics 1,2,3] [--force]
"""
import argparse
import json
import multiprocessing as mp
import time
import traceback

import numpy as np
import pandas as pd

import ind_common as ic
from ind_common import core

SEED = 20261123
OUT = ic.WORK / "arm_b"
N_DRAWS = 4
MASK_D = 0.75
MAX_PHASE_DRAWS = 5000
TEFF_SUN = 5772.0


def limb_darkening(teff, logg, table):
    """Nearest grid point: nearest Teff first, then the nearest log g the grid has at that Teff."""
    cool = teff < 3500
    te = table["Teff"].unique()
    t_near = te[np.argmin(np.abs(te - max(teff, 3500.0)))]
    g = table[table["Teff"] == t_near]
    r = g.iloc[int(np.argmin(np.abs(g["logg"].to_numpy() - logg)))]
    return float(r["aLSM"]), float(r["bLSM"]), float(r["Teff"]), float(r["logg"]), bool(cool)


def t14_days(P, aRs, k, b):
    inc = np.arccos(b / aRs)
    x = np.sqrt((1 + k) ** 2 - b ** 2) / (aRs * np.sin(inc))
    return P / np.pi * np.arcsin(min(x, 1.0)), inc


def inject(secs, t_ref, P, k, aRs, inc_rad, u, T14):
    import batman
    pm = batman.TransitParams()
    pm.t0, pm.per, pm.rp, pm.a, pm.inc = t_ref, P, k, aRs, np.degrees(inc_rad)
    pm.ecc, pm.w, pm.u, pm.limb_dark = 0.0, 90.0, list(u), "quadratic"
    out, max_depth, n_inj = [], 0.0, 0
    for s in secs:
        t = s["t"]
        ph = (t - t_ref + P / 2) % P - P / 2
        if not np.any(np.abs(ph) < T14):
            out.append(s)
            continue
        cad = core.cadence_days(t)
        m = batman.TransitModel(pm, t, supersample_factor=7, exp_time=cad).light_curve(pm)
        max_depth = max(max_depth, float(1 - m.min()))
        n_inj += 1
        out.append(dict(s, f=s["f"] * m))
    return out, max_depth, n_inj


def do_host(r):
    tic = int(r["tic"])
    p = OUT / f"{tic}.json"
    t_start = time.time()
    base = dict(tic=tic, tmag=float(r["Tmag"]), tmag_bin=r["tmag_bin"], arm="B")
    try:
        secs, errs = ic.load_star_raw(tic)
        info = ic.tic_params(tic)
        if not secs:
            raise RuntimeError("no light curve: " + "; ".join(errs))
        mstar, rstar, teff = info.get("mstar"), info.get("rstar"), info.get("teff")
        solar = mstar is None or rstar is None or mstar <= 0 or rstar <= 0
        # solar values for the injection; the pipeline's minimum-period limit gets what the TIC has
        m_pipe, r_pipe = (np.nan, np.nan) if solar else (mstar, rstar)
        if solar:
            mstar, rstar = 1.0, 1.0
        teff_missing = teff is None
        if teff_missing:
            teff = TEFF_SUN
        logg = 4.438 + np.log10(mstar) - 2 * np.log10(rstar)
        table = pd.read_csv(ic.HERE.parent / "data" / "claret2017_table25_Z0_xi2.csv")
        u1, u2, ld_teff, ld_logg, cool = limb_darkening(teff, logg, table)
        base.update(load_errors="; ".join(errs), solar_values=bool(solar), teff_missing=bool(teff_missing),
                    cool_star_3500=cool, teff=teff, logg=logg, mstar=mstar, rstar=rstar,
                    u1=u1, u2=u2, ld_teff=ld_teff, ld_logg=ld_logg)
        # mask the catalogued events
        cat = pd.read_csv(ic.ROOT / "gaia_dr4" / "candidates_frozen.csv")
        times = cat.loc[cat["tic"] == tic, ["t1_btjd", "t2_btjd"]].to_numpy(float).ravel()
        times = times[np.isfinite(times)]
        masked = []
        for s in secs:
            keep = np.ones(s["t"].size, bool)
            for tc in times:
                keep &= np.abs(s["t"] - tc) > MASK_D
            if keep.sum() >= 50:
                masked.append(dict(s, t=s["t"][keep], f=s["f"][keep]))
        base["n_catalogued_times_masked"] = int(times.size)
        if not masked:
            raise RuntimeError("no data left after masking")
        Traw = np.sort(np.concatenate([s["t"] for s in masked]))
        rng = np.random.default_rng([SEED, tic])
        res = []
        for d in range(N_DRAWS):
            n_param_redraws = 0
            while True:
                P = float(np.exp(rng.uniform(np.log(20), np.log(1000))))
                depth = float(np.exp(rng.uniform(np.log(500e-6), np.log(20000e-6))))
                b = float(rng.uniform(0, 0.9))
                k = np.sqrt(depth)
                aRs = float(core.a_over_rstar(P, mstar, rstar))
                if aRs > 1 + k + 0.05:
                    break
                n_param_redraws += 1  # orbit inside the star
            T14, inc = t14_days(P, aRs, k, b)
            t_ref, n_ph = None, 0
            while n_ph < MAX_PHASE_DRAWS:
                n_ph += 1
                ph = Traw[0] + rng.uniform(0, P)
                kk, ep = ic.epochs_in_span(ph, P, Traw[0], Traw[-1])
                cl = ic.clean_mask(Traw, ep, T14)
                if cl.any():
                    t_ref = float(ep[np.where(cl)[0][0]])
                    break
            trial = dict(base, draw=d, inj_P=P, inj_depth_ppm=depth * 1e6, inj_b=b, inj_duration_h=T14 * 24,
                         inj_k=k, inj_aRs=aRs, n_phase_draws=n_ph, n_param_redraws=n_param_redraws)
            if t_ref is None:
                res.append(dict(trial, status="skipped: no clean transit in %d phase draws" % MAX_PHASE_DRAWS))
                continue
            t_cat = t_ref + rng.normal(0, ic.T_SCATTER_D)
            inj, max_depth, n_inj = inject(masked, t_ref, P, k, aRs, inc, (u1, u2), T14)
            tr = ic.run_trial(inj, t_cat, depth, T14, np.nan, P, t_ref, T14, m_pipe, r_pipe)
            res.append(dict(trial, t_ref_true=t_ref, t_cat_offset_min=(t_cat - t_ref) * 1440,
                            inj_model_max_depth_ppm=max_depth * 1e6, n_sectors_injected=n_inj, **tr))
        p.write_text(json.dumps(ic.jsonable(res), indent=1))
        return tic, "ok", time.time() - t_start
    except Exception as e:
        res = [dict(base, status=f"error: {type(e).__name__}: {str(e)[:200]}",
                    traceback=traceback.format_exc()[-1500:])]
        p.write_text(json.dumps(ic.jsonable(res), indent=1))
        return tic, res[0]["status"], time.time() - t_start


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--tics", default="")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    s = pd.read_csv(ic.HERE.parent / "arm_b_hosts.csv")
    if a.tics:
        s = s[s["tic"].isin([int(x) for x in a.tics.split(",")])]
    if a.limit:
        s = s.head(a.limit)
    todo = [r for r in s.to_dict("records") if a.force or not (OUT / f"{int(r['tic'])}.json").exists()]
    print(f"arm B: {len(todo)} of {len(s)} hosts to do", flush=True)
    with mp.get_context("spawn").Pool(a.workers, maxtasksperchild=10) as pool:
        for i, (tic, st, dt) in enumerate(pool.imap_unordered(do_host, todo), 1):
            print(f"[{i}/{len(todo)}] TIC {tic}: {st} ({dt:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
