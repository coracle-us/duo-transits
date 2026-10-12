"""
Shared K2-singles analysis code, rewritten after the Fable 5.1 adversarial review (2026-10-12,
k2singles/review_fable/; findings summarised in k2singles/README.md). It replaces run_pilot.py and
run_full.py, which stay in the repo for the record.

Fixes relative to run_full.py:
  C1  no monkeypatching: TESS and K2 light curves are loaded explicitly per target, and the runner
      uses max_tasks_per_child=1 (the old closure chain mixed other stars' K2 data into later targets)
  C2  robust_periods(): every repeat can anchor trial periods. Each period reports how many repeats
      it explains; a single spurious repeat no longer deletes the true period
  M4  no K2 'P_min' cut. K2 data are now in the search, so K2 flat windows exclude periods directly
      (the old cut ignored gaps)
  M5  'excluded' windows need cover >= 0.8 (V-shaped transits sampled only on their wings were
      being called flat)
  M6  K2 repeats need n_in >= 4 and beta <= 2. K2 cadences flagged 'possible thruster firing'
      (bit 524288) are masked in addition to lightkurve's default mask. The HARD mask is deliberately
      NOT used: EPIC 211924561's C5 transit sits in 'rolling band' cadences.
  min_allowed_P / allowed_ranges are not reported (M3: grid artefacts). frac_logP_excluded is kept.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "singles" / "scripts"))
import alias_pipeline as ap  # noqa: E402
import pilot  # noqa: E402

CACHE = ROOT / "k2singles" / "cache"
K2_FOLDERS = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "102", "111", "112", "12", "13", "14", "15", "16",
              "17", "18", "19"]
THRUSTER_POSSIBLE = 524288
EXCL_MIN_COVER = 0.8
K2_MIN_NIN, K2_MAX_BETA = 4, 2.0
PRODUCTS = ap.PRODUCT_RANK + ["OWN"]


def _get(url, dest, tries=3):
    import time
    import urllib.error
    import urllib.request
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                dest.write_bytes(r.read())
            return True
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return False
            time.sleep(5 * (k + 1))
        except Exception:  # noqa: BLE001
            time.sleep(5 * (k + 1))
    raise RuntimeError(f"download failed after {tries} tries: {url}")


def k2_lightcurves(epic):
    """All K2 campaigns of this EPIC (K2 mission PDCSAP), cached as npz. The mask is
    DEFAULT + 'possible thruster firing'. Network failures raise; nothing is cached on failure."""
    from astropy.io import fits
    from lightkurve.utils import KeplerQualityFlags
    cdir = CACHE / str(epic)
    (cdir / "dl").mkdir(parents=True, exist_ok=True)
    s = str(epic)
    if not (cdir / "checked_v2.txt").exists():
        for c in K2_FOLDERS:
            fn = cdir / "dl" / f"ktwo{s}-c{int(c):02d}_llc.fits"
            if fn.exists():
                continue
            url = (f"https://archive.stsci.edu/missions/k2/lightcurves/c{c}/{s[:4]}00000/{s[4:6]}000/"
                   f"ktwo{s}-c{int(c):02d}_llc.fits")
            _get(url, fn)
        for old in cdir.glob("k2_c*.npz"):
            old.unlink()
        mask = KeplerQualityFlags.DEFAULT_BITMASK | THRUSTER_POSSIBLE
        for fn in sorted((cdir / "dl").glob("ktwo*_llc.fits")):
            camp = int(fn.stem.split("-c")[1].split("_")[0])
            with fits.open(fn) as h:
                d = h[1].data
                t = np.asarray(d["TIME"], float) + 2454833 - 2457000
                f = np.asarray(d["PDCSAP_FLUX"], float)
                q = np.asarray(d["SAP_QUALITY"], int)
            ok = np.isfinite(t) & np.isfinite(f) & ((q & mask) == 0)
            if ok.sum() > 50:
                np.savez(cdir / f"k2_c{camp:03d}.npz", t=t[ok], f=f[ok] / np.nanmedian(f[ok]))
        (cdir / "checked_v2.txt").write_text("all K2 campaign folders checked (v2 mask)")
    parts = []
    for fn in sorted(cdir.glob("k2_c*.npz")):
        camp = int(fn.stem.split("_c")[1])
        z = np.load(fn)
        o = np.argsort(z["t"])
        t, f = z["t"][o], z["f"][o]
        sec = 1000 + camp
        seg = sec * 1000 + np.concatenate([[0], np.cumsum(np.diff(t) > 0.5)])
        parts.append(pd.DataFrame({"time": t, "flux": f, "err": np.nan, "sector": sec, "product": "K2PDC",
                                   "segment": seg}))
    return pd.concat(parts, ignore_index=True) if parts else None


def load_combined(tic, epic):
    tess = ap.load_lightcurves(tic, products=PRODUCTS)
    k2 = k2_lightcurves(epic)
    lc = tess if k2 is None else pd.concat([tess, k2], ignore_index=True)
    lc = lc.sort_values("time").reset_index(drop=True)
    # guard against the C1 failure mode: duplicate timestamps within a sector
    dup = lc.duplicated(subset=["sector", "time"]).sum()
    if dup:
        raise RuntimeError(f"{dup} duplicate timestamps within sectors")
    return lc, k2


def grid(lc, D, delta, epoch):
    g = pilot.grid_windows(lc, D, delta, epoch)
    if "cover" in g:
        weak = (g.verdict == "excluded") & (g.cover.fillna(0) < EXCL_MIN_COVER)
        g.loc[weak, "verdict"] = "unconstrained"
    return g


def find_repeats(g, delta, D):
    reps = []
    if "snr" not in g:
        return reps
    for _, r in g[~g.near_event].sort_values("snr", ascending=False).iterrows():
        if not pilot.robust(r, delta):
            continue
        if r.sector >= 1000 and (r.n_in < K2_MIN_NIN or r.beta > K2_MAX_BETA):
            continue
        if all(abs(r.T_pred - x["T"]) > D for x in reps):
            reps.append(dict(T=round(r.T_pred, 4), sector=int(r.sector), product=r["product"],
                             depth_ppm=round(r.depth * 1e6), snr=round(r.snr, 1), n_in=int(r.n_in),
                             beta=round(float(r.beta), 2)))
    return sorted(reps, key=lambda x: x["T"])


def robust_periods(g, epoch, reps, D, p_min=5.0):
    """Trial periods anchored on EVERY repeat: P = (T_anchor - epoch)/n. For each, the repeats
    consistent with a linear ephemeris (|O-C| < D/3) are collected iteratively; the period is
    refitted from the epoch plus those repeats, and kept if no predicted epoch falls in a flat window.
    Returns a list of dicts sorted by support (number of repeats explained), then P. Periods
    within 0.1% are merged."""
    out = []
    if not reps:
        return out
    Ts = np.array([x["T"] for x in reps])
    for a, Ta in enumerate(Ts):
        dT = abs(Ta - epoch)
        for n in range(1, int(dT / p_min) + 1):
            P = dT / n
            used = [a]
            for _ in range(2):
                cyc = np.round((Ts[used] - epoch) / P)
                A = np.vstack([np.ones(len(used) + 1), np.concatenate([[0], cyc])]).T
                y = np.concatenate([[epoch], Ts[used]])
                T0, Pf = np.linalg.lstsq(A, y, rcond=None)[0]
                allc = np.round((Ts - T0) / Pf)
                res = np.abs(Ts - (T0 + allc * Pf))
                new = [i for i in range(len(Ts)) if res[i] < D / 3]
                if a not in new:
                    break
                used = sorted(set(new))
                P = Pf
            if a not in used:
                continue
            cyc = np.round((Ts[used] - epoch) / P)
            if len(set(cyc)) < len(used) or 0 in cyc:
                continue
            if pilot.excluded_by_flat(g, epoch, P, D):
                continue
            out.append(dict(P=float(P), n_support=len(used), support_T=[float(Ts[i]) for i in used]))
    out.sort(key=lambda d: (-d["n_support"], d["P"]))
    merged = []
    for d in out:
        if all(abs(d["P"] - m["P"]) / m["P"] > 1e-3 for m in merged):
            merged.append(d)
    return merged


def frac_excluded(g, epoch, D):
    return pilot.period_constraints(g, epoch, D)["frac_logP_excluded"]
