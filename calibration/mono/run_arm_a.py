"""
Arm A of the single-transit calibration: confirmed planets treated as single transits.

For each star of arm_a_sample.csv: the TOI ephemeris is refined on the star's own detrended data
(s10_epoch_rule_test.refine_ephemeris, as in the 2026-11-01 test), the first clean transit becomes
the catalogued single transit, and the star is run twice through monocal.run_trial:
  A1  every other transit window (+-1.5 durations) removed from the raw flux;
  A2  nothing removed.

Inputs  : calibration/mono/arm_a_sample.csv (samples.py), MAST light curves
Outputs : <out>/arm_a/<tic>.json per star; <out>/arm_a_trials.csv and arm_a_stars.csv
Usage   : python calibration/mono/run_arm_a.py [--workers 4]
"""
import argparse
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import monocal as mc  # noqa: E402
import runner  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
SEED = 20261122
T_CAT_ERR_D = 15 / 1440
MIN_CLEAN = 3


def _s10():
    spec = importlib.util.spec_from_file_location(
        "s10", ROOT / "salinas-second-transits/round2/s10_epoch_rule_test.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def run_star(item, cache):
    tic = int(item["tic"])
    rng = np.random.default_rng([SEED, tic])
    t_a = time.time()
    sectors, errors = mc.load_star(tic, Path(cache))
    info = mc.pl.tic_info(tic)
    meta = dict(tic=tic, toi=item["toi"], n_sectors=len(sectors), n_load_errors=len(errors),
                load_s=time.time() - t_a, n_trials=0, skipped="")
    if not sectors:
        return [], {**meta, "skipped": "no light curve"}
    dur, depth, P0 = item["duration_h"] / 24, item["depth_ppm"] * 1e-6, item["period_d"]
    T, F = [], []
    for _, _, t, f in sectors:
        try:
            td, fd = mc.pl.detrend(t, f, dur)
            T.append(td)
            F.append(fd)
        except Exception:
            pass
    T, F = np.concatenate(T), np.concatenate(F)
    o = np.argsort(T)
    perr = item.get("period_err")
    half = max(5 * (perr if perr is not None and np.isfinite(perr) else 0.0), 2e-4 * P0)
    eph = _s10().refine_ephemeris(T[o], F[o], P0, item["t0"], dur, half)
    P, t0 = eph["P"], eph["t0"]
    meta.update(period_d=P, period_offset_d=eph["P_offset_d"], at_grid_edge=bool(eph["at_grid_edge"]),
                stack_depth_ppm=eph["depth"] * 1e6)
    Traw = np.sort(np.concatenate([t for _, _, t, _ in sectors]))
    _, ep = mc.epochs(t0, P, Traw.min(), Traw.max())
    cl = mc.clean_mask(Traw, ep, dur)
    meta["n_clean_transits"] = int(cl.sum())
    if cl.sum() < MIN_CLEAN:
        return [], {**meta, "skipped": f"under {MIN_CLEAN} clean transits"}
    ref = float(ep[cl][0])
    t_cat = ref + rng.normal(0, T_CAT_ERR_D)
    others = ep[np.abs(ep - ref) > dur / 2]
    hidden = [(s, a, *mc.mask_windows(t, f, others, 1.5 * dur)) for s, a, t, f in sectors]
    hidden = [x for x in hidden if x[2].size >= 50]
    mstar, rstar = info.get("mstar", np.nan), info.get("rstar", np.nan)
    rows = []
    for arm, sec in (("A1", hidden), ("A2", sectors)):
        r = mc.run_trial(sec, t_cat, dur, depth, mstar, rstar, truth_P=P, truth_t0=t0)
        rows.append(dict(arm=arm, tic=tic, toi=item["toi"], depth_bin=item["depth_bin"],
                         n_toi_on_star=item["n_toi_on_star"], tmag=item["tmag"], truth_P=P, truth_t0=t0,
                         t_cat=t_cat, depth_cat=depth, dur_cat_d=dur, at_grid_edge=bool(eph["at_grid_edge"]),
                         stack_depth_ppm=eph["depth"] * 1e6, **r))
    meta["n_trials"] = len(rows)
    return rows, meta


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache", default=str(ROOT / "calibration/mono/work/lc_cache"))
    ap.add_argument("--out", default=str(ROOT / "calibration/mono/results"))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, help="first N stars only (smoke test)")
    a = ap.parse_args()
    s = pd.read_csv(ROOT / "calibration/mono/arm_a_sample.csv")
    items = s.head(a.limit).to_dict("records") if a.limit else s.to_dict("records")
    out = Path(a.out)
    Path(a.cache).mkdir(parents=True, exist_ok=True)
    runner.run(run_star, items, out / "arm_a", a.workers, cache=a.cache)
    rows, metas = runner.collect(out / "arm_a")
    rows.to_csv(out / "arm_a_trials.csv", index=False)
    metas.to_csv(out / "arm_a_stars.csv", index=False)


if __name__ == "__main__":
    main()
