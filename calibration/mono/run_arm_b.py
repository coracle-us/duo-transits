"""
Arm B of the single-transit calibration: inject planets into raw light curves of frozen hosts,
before detrending, and run each through monocal.run_trial.

Inputs  : gaia_dr4/candidates_frozen.csv (catalogued event times to mask), MAST light curves
Outputs : <out>/arm_b_trials.csv  (one row per trial; NOT written with --timing-only)
          <out>/arm_b_timing.csv  (seconds per star and per trial; no outcomes)
Usage   : python calibration/mono/run_arm_b.py --tics 1,2,3 --timing-only
          python calibration/mono/run_arm_b.py --hosts hosts.csv        (after the sample is fixed)

--timing-only is the pilot the pre-registration allows before the sample size is fixed: outcomes
are discarded in memory and never written or printed.

Not final until PROPOSED_v1.1.md is approved: LIMB_DARK, EVENT_MASK_D and T_CAT_ERR_D below.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import monocal as mc  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
SEED = 20261123
DRAWS = 4
P_RANGE, DEPTH_RANGE, B_MAX = (20.0, 1000.0), (500e-6, 20000e-6), 0.9
LIMB_DARK = (0.35, 0.25)      # placeholder; v1.0 says "from the TIC Teff and logg"
EVENT_MASK_D = 0.75           # placeholder; v1.0 says +-1.5 durations, the frozen list has no durations
T_CAT_ERR_D = 15 / 1440       # scatter of the catalogued time, as in s10_epoch_rule_test.py


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


def inject(sectors, pl_):
    out = []
    for s, a, t, f in sectors:
        m = mc.transit_flux(t, pl_["t0"], pl_["P"], pl_["k"], pl_["ar"], pl_["b"], LIMB_DARK,
                            mc.core.cadence_days(t))
        out.append((s, a, t, f * m))
    return out


def run_host(tic, cache, events, timing_only):
    rng = np.random.default_rng([SEED, int(tic)])
    t_a = time.time()
    sectors, errors = mc.load_star(tic, cache)
    info = mc.pl.tic_info(tic)
    t_load = time.time() - t_a
    trials, timing = [], dict(tic=int(tic), n_sectors=len(sectors), n_load_errors=len(errors),
                              n_points=int(sum(t.size for _, _, t, _ in sectors)), load_s=t_load)
    if len(sectors) == 0:
        return trials, {**timing, "trial_s": np.nan, "n_trials": 0}
    sectors = [(s, a, *mc.mask_windows(t, f, events, EVENT_MASK_D)) for s, a, t, f in sectors]
    mstar, rstar = info.get("mstar", np.nan), info.get("rstar", np.nan)
    fallback = not (np.isfinite(mstar) and np.isfinite(rstar) and mstar > 0 and rstar > 0)
    if fallback:
        mstar, rstar = 1.0, 1.0
    secs = []
    for d in range(DRAWS):
        p = draw(rng, sectors, mstar, rstar)
        if p is None:
            continue
        t_b = time.time()
        t_cat = p["t_ref"] + rng.normal(0, T_CAT_ERR_D)
        r = mc.run_trial(inject(sectors, p), t_cat, p["dur_d"], p["depth"], mstar, rstar,
                         truth_P=p["P"], truth_t0=p["t0"])
        secs.append(time.time() - t_b)
        if not timing_only:
            trials.append(dict(tic=int(tic), draw=d, tmag=info.get("tmag", np.nan),
                               stellar_fallback=fallback, **{f"inj_{k}": v for k, v in p.items()},
                               t_cat=t_cat, **r))
        del r
    return trials, {**timing, "trial_s": float(np.mean(secs)) if secs else np.nan, "n_trials": len(secs)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tics", help="comma list of TICs")
    ap.add_argument("--hosts", help="CSV with a tic column")
    ap.add_argument("--cache", default=str(ROOT / "calibration/mono/work/lc_cache"))
    ap.add_argument("--out", default=str(ROOT / "calibration/mono/results"))
    ap.add_argument("--timing-only", action="store_true")
    a = ap.parse_args()
    tics = [int(x) for x in a.tics.split(",")] if a.tics else pd.read_csv(a.hosts).tic.astype(int).tolist()
    cand = pd.read_csv(ROOT / "gaia_dr4/candidates_frozen.csv")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    Path(a.cache).mkdir(parents=True, exist_ok=True)
    trials, timing = [], []
    for i, tic in enumerate(tics, 1):
        c = cand[cand.tic == tic]
        ev = pd.concat([c.t1_btjd, c.t2_btjd]).dropna().to_numpy()
        try:
            tr, tm = run_host(tic, Path(a.cache), ev, a.timing_only)
        except Exception as e:
            tr, tm = [], dict(tic=tic, error=str(e)[:200])
        trials += tr
        timing.append(tm)
        print(f"[{i}/{len(tics)}] TIC {tic}: " + ", ".join(
            f"{k}={v:.1f}" if isinstance(v, float) else f"{k}={v}" for k, v in tm.items() if k != "tic"),
            flush=True)
        pd.DataFrame(timing).to_csv(out / "arm_b_timing.csv", index=False)
        if not a.timing_only:
            pd.DataFrame(trials).to_csv(out / "arm_b_trials.csv", index=False)


if __name__ == "__main__":
    main()
