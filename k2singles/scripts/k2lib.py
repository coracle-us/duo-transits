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

v3 (2026-10-12, after the second Fable review, k2singles/review_fable2/REPORT.md):
  P1  the red-noise factor was silently set to 1 whenever the 0.75-d baseline held fewer than three
      in-transit-length bins, i.e. always for D >~ 0.5 d at K2 cadence. Granulation on giants then
      made most windows 'ambiguous' and alias lists far too long. test_window_k2 uses a baseline of
      max(0.75 d, 2 D) on each side; when beta still cannot be measured, the window cannot exclude.
      A dip below a quarter of the expected depth no longer blocks exclusion, and very strong
      exclusions (d + 3 sigma < delta/4) are allowed up to beta 5.
  P3  K2 repeats: n_in >= 4 and (beta <= 2 or SNR >= 15); the SNR already includes beta. Baseline
      sides: >= 0.25 x 0.75 d of cadences each and >= 0.75 d in total (was 0.5 x 0.75 d each side,
      which threw away a 63-sigma S92 eclipse with a short one-sided baseline).
  P2  partial events at sector/campaign edges are still not tested (known limitation).
  P4  repeat times are grid-quantised (step D/4); periods for the note come from dedicated fits.
  The singles project's alias_pipeline.py / pilot.py are NOT changed.
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


BETA_MAX_EXCL_STRONG = 5.0
SIDE_FRAC_EACH, SIDE_TOTAL = 0.25, 1.0   # x 0.75 d of cadences
K2_MIN_SNR_HIGH_BETA = 15.0


def red_noise_beta_k2(r, n_bin):
    """As ap.red_noise_beta, but NaN (not 1.0) when there are fewer than three bins."""
    if n_bin <= 1:
        return 1.0
    nb = len(r) // n_bin
    if nb < 3:
        return np.nan
    rb = r[: nb * n_bin].reshape(nb, n_bin).mean(axis=1)
    sw = np.std(r, ddof=1)
    if sw == 0:
        return np.nan
    return float(max(1.0, np.std(rb, ddof=1) / (sw / np.sqrt(n_bin)) * np.sqrt(nb / (nb - 1))))


def test_window_k2(lc_t, lc_f, T, D, delta, lc_seg):
    """ap.test_window (sigma_T = 0) with the v3 changes described in the module docstring."""
    base_w = max(ap.BASE_W, 2 * D)
    half = D / 2
    mask_half = half + 0.02
    i0, i1 = np.searchsorted(lc_t, [T - mask_half - base_w, T + mask_half + base_w])
    t, f, sg = lc_t[i0:i1], lc_f[i0:i1], lc_seg[i0:i1]
    oot = np.abs(t - T) > mask_half
    anchored = np.isin(sg, np.unique(sg[oot])) if oot.any() else np.zeros(len(t), bool)
    t, f, sg, oot = t[anchored], f[anchored], sg[anchored], oot[anchored]
    res = dict(T_pred=T, sigma_T=0.0, n_window=len(t), base_w=base_w)
    if len(t) == 0:
        return {**res, "verdict": "no-data"}
    n_left, n_right = np.sum(oot & (t < T)), np.sum(oot & (t > T))
    res.update(n_oot_left=int(n_left), n_oot_right=int(n_right))
    if n_left + n_right < 10:
        flat_any = np.sum(np.abs(t - T) < ap.FLAT * half)
        return {**res, "verdict": "unconstrained" if flat_any else "no-data"}
    deg = 2 if (n_left >= 10 and n_right >= 10) else (1 if min(n_left, n_right) >= 5 else 0)
    segs = np.unique(sg)
    X = np.column_stack([(t - T) ** p for p in range(1, deg + 1)] + [(sg == u).astype(float) for u in segs])
    coef, *_ = np.linalg.lstsq(X[oot], f[oot], rcond=None)
    r = f / (X @ coef) - 1.0
    r_oot = r[oot]
    sig_w = 1.4826 * np.median(np.abs(r_oot - np.median(r_oot)))
    cad = np.median(np.diff(t)) if len(t) > 1 else 0.0
    n_flat_full = max(1, int(round(ap.FLAT * D / cad))) if cad > 0 else 1
    beta = red_noise_beta_k2(r_oot, n_flat_full)
    b_det = 1.0 if np.isnan(beta) else beta
    inn = np.abs(t - T) < ap.FLAT * half
    n_in = int(inn.sum())
    cover = n_in * cad / (ap.FLAT * D) if cad > 0 else 0.0
    res.update(beta=b_det, beta_measured=not np.isnan(beta), sigma_white=sig_w, baseline_deg=deg, cadence=cad)
    if n_in < ap.MIN_IN or cover < ap.MIN_COVER:
        return {**res, "verdict": "unconstrained"}
    d = -np.mean(r[inn])
    s = b_det * sig_w / np.sqrt(n_in)
    res.update(offset=0.0, depth=d, depth_err=s, snr=d / s, n_in=n_in, cover=min(cover, 1.0))
    if d > ap.DEPTH_FRAC * delta and d / s > ap.SNR_MIN:
        return {**res, "verdict": "supported"}
    strong = d + ap.K_SIG * s < 0.25 * delta
    excl = (not np.isnan(beta) and d + ap.K_SIG * s < ap.DEPTH_FRAC * delta
            and abs(d) < max(ap.K_SIG * s, 0.25 * delta)
            and (beta <= ap.BETA_MAX_EXCL or (strong and beta <= BETA_MAX_EXCL_STRONG)))
    return {**res, "verdict": "excluded" if excl else "ambiguous"}


def grid_windows_k2(lc, D, delta, epoch):
    """pilot.grid_windows with test_window_k2."""
    t, f, sg = lc.time.to_numpy(), lc.flux.to_numpy(), lc.segment.to_numpy()
    T_all = np.arange(t.min(), t.max(), D / 4)
    j = np.searchsorted(t, T_all)
    j0, j1 = np.clip(j - 1, 0, len(t) - 1), np.clip(j, 0, len(t) - 1)
    near = np.minimum(np.abs(t[j0] - T_all), np.abs(t[j1] - T_all)) < D / 2
    rows = []
    for T in T_all[near]:
        r = test_window_k2(t, f, T, D, delta, sg)
        k = min(np.searchsorted(t, T), len(t) - 1)
        r["sector"] = int(lc.sector.iloc[k])
        r["product"] = lc["product"].iloc[k]
        rows.append(r)
    g = pd.DataFrame(rows)
    g["near_event"] = (g.T_pred - epoch).abs() < D
    return g


def robust_k2(r, delta):
    """pilot.robust with the v3 baseline-side rule."""
    if r.get("verdict") != "supported":
        return False
    per_bw = ap.BASE_W / r["cadence"]
    return (r["snr"] >= pilot.DETECT_SNR and pilot.DMIN * delta <= r["depth"] <= pilot.DMAX * delta
            and r["cover"] >= ap.DETECT_COVER
            and min(r["n_oot_left"], r["n_oot_right"]) >= SIDE_FRAC_EACH * per_bw
            and r["n_oot_left"] + r["n_oot_right"] >= SIDE_TOTAL * per_bw)


def grid(lc, D, delta, epoch, v3=True):
    if v3:
        return grid_windows_k2(lc, D, delta, epoch)
    g = pilot.grid_windows(lc, D, delta, epoch)
    if "cover" in g:
        weak = (g.verdict == "excluded") & (g.cover.fillna(0) < EXCL_MIN_COVER)
        g.loc[weak, "verdict"] = "unconstrained"
    return g


def find_repeats(g, delta, D, v3=True):
    reps = []
    if "snr" not in g:
        return reps
    for _, r in g[~g.near_event].sort_values("snr", ascending=False).iterrows():
        if not (robust_k2(r, delta) if v3 else pilot.robust(r, delta)):
            continue
        if r.sector >= 1000 and (r.n_in < K2_MIN_NIN or (r.beta > K2_MAX_BETA
                                                         and not (v3 and r.snr >= K2_MIN_SNR_HIGH_BETA))):
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
