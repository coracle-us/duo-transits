#!/usr/bin/env python3
"""Arm A (confirmed planets) of the blind independent re-derivation: A1 hidden, A2 visible.

Purpose
    For each of the 120 TOIs in arm_a_sample.csv: refine the TOI ephemeris on the star's own
    detrended data (s10_epoch_rule_test.refine_ephemeris, half-width max(5 period_err, 2e-4 P)),
    take the first clean transit as the reference, scatter its time by 15 min (one draw per star,
    used by both arms), then run the single-transit pipeline twice:
      A2  nothing removed
      A1  every other transit window of the planet (+-1.5 TOI durations) removed from the raw
          light curve before detrending.

Inputs
    ../arm_a_sample.csv; MAST light curves and TIC (see ind_common.py).

Outputs
    work/arm_a/<tic>.json  one file per star: [A1 trial, A2 trial] or one 'skipped' record.

Usage
    python ind_arm_a.py [--workers 4] [--limit N] [--tics 1,2,3] [--force]
"""
import argparse
import json
import multiprocessing as mp
import time
import traceback

import numpy as np
import pandas as pd

import ind_common as ic
from ind_common import core, pl

SEED = 20261122
OUT = ic.WORK / "arm_a"


def do_star(r):
    tic = int(r["tic"])
    p = OUT / f"{tic}.json"
    t_start = time.time()
    base = dict(tic=tic, toi=float(r["toi"]), n_toi_on_star=int(r["n_toi_on_star"]), depth_bin=r["depth_bin"],
                tmag=float(r["tmag"]), toi_period_d=float(r["period_d"]), toi_depth_ppm=float(r["depth_ppm"]),
                toi_duration_h=float(r["duration_h"]))
    try:
        secs, errs = ic.load_star_raw(tic)
        info = ic.tic_params(tic)
        if not secs:
            raise RuntimeError("no light curve: " + "; ".join(errs))
        mstar, rstar = info.get("mstar"), info.get("rstar")
        # arm A injects nothing, so mass and radius only feed the pipeline's circular-orbit minimum
        # period; as in the pipeline, a star without them simply gets no such limit
        no_mr = mstar is None or rstar is None
        if no_mr:
            mstar, rstar = np.nan, np.nan
        base.update(load_errors="; ".join(errs), no_tic_mass_or_radius=no_mr)
        dur_d, d0 = float(r["duration_h"]) / 24.0, float(r["depth_ppm"]) * 1e-6
        d0_err = float(r["depth_err"]) * 1e-6 if np.isfinite(r["depth_err"]) else np.nan
        T, F = ic.concat([pl.detrend(s["t"], s["f"], dur_d) for s in secs])
        import s10_epoch_rule_test as s10
        half = max(5 * np.nan_to_num(float(r["period_err"])), 2e-4 * float(r["period_d"]))
        eph = s10.refine_ephemeris(T, F, float(r["period_d"]), float(r["t0"]), dur_d, half)
        P, t0 = eph["P"], eph["t0"]
        base.update(period_offset_d=eph["P_offset_d"], at_grid_edge=bool(eph["at_grid_edge"]),
                    stack_depth_ppm=eph["depth"] * 1e6)
        k, ep = ic.epochs_in_span(t0, P, T.min(), T.max())
        clean = ic.clean_mask(T, ep, dur_d, F)
        base["n_clean_transits"] = int(clean.sum())
        if clean.sum() < 3:
            res = [dict(base, arm="A", status="skipped: under 3 transits in clean data", truth_P=P)]
            p.write_text(json.dumps(ic.jsonable(res), indent=1))
            return tic, res[0]["status"], time.time() - t_start
        i_ref = int(np.where(clean)[0][0])
        rng = np.random.default_rng([SEED, tic])
        t_ref = float(ep[i_ref])
        t_cat = t_ref + rng.normal(0, ic.T_SCATTER_D)
        base.update(t_ref_true=t_ref, t_cat_offset_min=(t_cat - t_ref) * 1440)
        others = np.delete(ep, i_ref)
        res = []
        for arm in ("A2", "A1"):
            if arm == "A1":
                use = []
                for s in secs:
                    keep = np.ones(s["t"].size, bool)
                    lo = np.searchsorted(others, s["t"].min() - 2 * dur_d)
                    hi = np.searchsorted(others, s["t"].max() + 2 * dur_d)
                    for e in others[lo:hi]:
                        keep &= np.abs(s["t"] - e) > 1.5 * dur_d
                    if keep.sum() >= 50:
                        use.append(dict(s, t=s["t"][keep], f=s["f"][keep]))
            else:
                use = secs
            tr = ic.run_trial(use, t_cat, d0, dur_d, d0_err, P, t_ref, dur_d, mstar, rstar)
            res.append(dict(base, arm=arm, **tr))
        res.sort(key=lambda x: x["arm"])
        p.write_text(json.dumps(ic.jsonable(res), indent=1))
        return tic, "ok", time.time() - t_start
    except Exception as e:
        res = [dict(base, arm="A", status=f"error: {type(e).__name__}: {str(e)[:200]}",
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
    s = pd.read_csv(ic.HERE.parent / "arm_a_sample.csv")
    if a.tics:
        s = s[s["tic"].isin([int(x) for x in a.tics.split(",")])]
    if a.limit:
        s = s.head(a.limit)
    todo = [r for r in s.to_dict("records") if a.force or not (OUT / f"{int(r['tic'])}.json").exists()]
    print(f"arm A: {len(todo)} of {len(s)} stars to do", flush=True)
    with mp.get_context("spawn").Pool(a.workers, maxtasksperchild=10) as pool:
        for i, (tic, st, dt) in enumerate(pool.imap_unordered(do_star, todo), 1):
            print(f"[{i}/{len(todo)}] TIC {tic}: {st} ({dt:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
