"""
Arm B of the single-transit calibration: inject planets into raw light curves of frozen hosts,
before detrending, and run each through monocal.run_trial.

Inputs  : calibration/mono/arm_b_hosts.csv (samples.py), gaia_dr4/candidates_frozen.csv (catalogued
          event times to mask), MAST light curves
Outputs : <out>/arm_b/<tic>.json per star; <out>/arm_b_trials.csv and arm_b_stars.csv
Usage   : python calibration/mono/run_arm_b.py [--workers 4]
          python calibration/mono/run_arm_b.py --tics 1,2,3 --timing-only   (the pilot)

--timing-only discards outcomes in memory; it is the pilot the pre-registration allows.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import monocal as mc  # noqa: E402
import runner  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
SEED = 20261123
DRAWS = 4
P_RANGE, DEPTH_RANGE, B_MAX = (20.0, 1000.0), (500e-6, 20000e-6), 0.9
EVENT_MASK_D = 0.75           # +-0.75 d around each catalogued time (version 1.1, item 4)
T_CAT_ERR_D = 15 / 1440       # scatter of the catalogued time (version 1.1, item 6)


def draw(rng, sectors, mstar, rstar, max_tries=100):
    """One planet whose ephemeris puts at least one transit in clean raw data."""
    T = np.sort(np.concatenate([t for _, _, t, _ in sectors]))
    P = float(np.exp(rng.uniform(*np.log(P_RANGE))))
    depth = float(np.exp(rng.uniform(*np.log(DEPTH_RANGE))))
    b = float(rng.uniform(0, B_MAX))
    k = np.sqrt(depth)
    dur = mc.transit_duration_d(P, mstar, rstar, k, b)
    for _ in range(max_tries):
        t0 = float(T.min() + rng.uniform(0, P))
        _, ep = mc.epochs(t0, P, T.min(), T.max())
        cl = mc.clean_mask(T, ep, dur)
        if cl.any():
            return dict(P=P, depth=depth, b=b, k=k, dur_d=dur, t0=t0, t_ref=float(ep[cl][0]),
                        ar=float(mc.core.a_over_rstar(P, mstar, rstar)))
    return None


def inject(sectors, p, u):
    out = []
    for s, a, t, f in sectors:
        m = mc.transit_flux(t, p["t0"], p["P"], p["k"], p["ar"], p["b"], u, mc.core.cadence_days(t))
        out.append((s, a, t, f * m))
    return out


def run_host(item, cache, timing_only=False):
    tic = int(item["tic"])
    events = np.asarray(item.get("events", []), float)
    rng = np.random.default_rng([SEED, tic])
    t_a = time.time()
    sectors, errors = mc.load_star(tic, Path(cache))
    info = mc.pl.tic_info(tic)
    meta = dict(tic=tic, n_sectors=len(sectors), n_load_errors=len(errors),
                n_points=int(sum(t.size for _, _, t, _ in sectors)), load_s=time.time() - t_a, n_trials=0)
    if not sectors:
        return [], meta
    sectors = [(s, a, *mc.mask_windows(t, f, events, EVENT_MASK_D)) for s, a, t, f in sectors]
    sectors = [x for x in sectors if x[2].size >= 50]
    mstar, rstar = info.get("mstar", np.nan), info.get("rstar", np.nan)
    u1, u2, ld_how = mc.limb_dark(info.get("teff", np.nan), mstar, rstar)
    fallback = not (np.isfinite(mstar) and np.isfinite(rstar) and mstar > 0 and rstar > 0)
    if fallback:
        mstar, rstar = 1.0, 1.0
    rows, secs = [], []
    for d in range(DRAWS):
        p = draw(rng, sectors, mstar, rstar)
        if p is None:
            continue
        t_b = time.time()
        t_cat = p["t_ref"] + rng.normal(0, T_CAT_ERR_D)
        r = mc.run_trial(inject(sectors, p, (u1, u2)), t_cat, p["dur_d"], p["depth"], mstar, rstar,
                         truth_P=p["P"], truth_t0=p["t0"])
        secs.append(time.time() - t_b)
        if not timing_only:
            rows.append(dict(arm="B", tic=tic, draw=d, tmag=info.get("tmag", np.nan), stellar_fallback=fallback,
                             limb_dark=ld_how, u1=u1, u2=u2, **{f"inj_{k}": v for k, v in p.items()},
                             t_cat=t_cat, **r))
    meta.update(n_trials=len(secs), trial_s=float(np.mean(secs)) if secs else np.nan)
    return rows, meta


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tics", help="comma list of TICs (default: arm_b_hosts.csv)")
    ap.add_argument("--cache", default=str(ROOT / "calibration/mono/work/lc_cache"))
    ap.add_argument("--out", default=str(ROOT / "calibration/mono/results"))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--timing-only", action="store_true")
    a = ap.parse_args()
    tics = ([int(x) for x in a.tics.split(",")] if a.tics
            else pd.read_csv(ROOT / "calibration/mono/arm_b_hosts.csv").tic.astype(int).tolist())
    cand = pd.read_csv(ROOT / "gaia_dr4/candidates_frozen.csv")
    items = []
    for tic in tics:
        c = cand[cand.tic == tic]
        items.append(dict(tic=tic, events=pd.concat([c.t1_btjd, c.t2_btjd]).dropna().tolist()))
    out = Path(a.out)
    Path(a.cache).mkdir(parents=True, exist_ok=True)
    sub = out / ("arm_b_timing" if a.timing_only else "arm_b")
    runner.run(run_host, items, sub, a.workers, cache=a.cache, timing_only=a.timing_only)
    rows, metas = runner.collect(sub)
    metas.to_csv(out / ("arm_b_timing.csv" if a.timing_only else "arm_b_stars.csv"), index=False)
    if not a.timing_only:
        rows.to_csv(out / "arm_b_trials.csv", index=False)


if __name__ == "__main__":
    main()
