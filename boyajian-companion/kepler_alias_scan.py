#!/usr/bin/env python3
"""Kepler alias scan for the Boyajian's star (KIC 8462852) TESS single transit.

Downloads Kepler long-cadence light curves (Q0-Q17) from MAST's public AWS bucket
(stpubdata), then, for every trial period P, predicts earlier transits
t_n = T0 + n*P that land in the Kepler baseline and asks whether the data allow
a 1.1%-deep, ~21 h transit there.

Criterion ("additive deficit", default): flux deficits add, so if a transit occurred,
every cadence in the transit core (|t - t_n| < T23/2) must sit at least 0.8% below
the local quiet level (95th percentile within +-3 d). P is excluded if any predicted
core contains >= 3 cadences with a smaller deficit. Transits falling in data gaps or
inside dips at least that deep survive (the same idea as Madurga-Favieres et al. 2026,
who masked Kepler dips deeper than 0.1%).

Usage:  python kepler_alias_scan.py 770 845        # scan 770-845 d in 0.01 d steps
        python kepler_alias_scan.py 770 845 --mask 0.002   # looser 'masked dip' rule
Outputs: prints surviving ranges; writes kepler_alias_survivors.csv
"""
import argparse, io, re, sys, urllib.request
import numpy as np
from astropy.io import fits

BUCKET = "https://stpubdata.s3.amazonaws.com/"
PREFIX = "kepler/public/lightcurves/0084/008462852/"
T0_BJD = 2458730.164          # Madurga-Favieres et al. 2026, Table 4 (TESS S15)
T0 = T0_BJD - 2454833.0       # in BKJD
T23 = 0.50                    # d, full-depth part of the transit (our trapezoid fit gives ~0.40)

def kepler_lc():
    keys = re.findall(r"<Key>([^<]*llc\.fits)</Key>",
                      urllib.request.urlopen(f"{BUCKET}?list-type=2&prefix={PREFIX}").read().decode())
    T, F = [], []
    for k in sorted(keys):
        d = fits.open(io.BytesIO(urllib.request.urlopen(BUCKET + k).read()))[1].data
        m = (d["SAP_QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"]) & np.isfinite(d["TIME"])
        fl = d["PDCSAP_FLUX"][m]
        T.append(d["TIME"][m]); F.append(fl / np.nanmedian(fl))
    t, f = np.concatenate(T), np.concatenate(F)
    o = np.argsort(t)
    return t[o], f[o]

def quiet_level(t, f, half=3.0, pct=95):
    grid = np.arange(t.min(), t.max() + 0.1, 0.1); base = np.full(grid.size, np.nan)
    for i, g in enumerate(grid):
        a, b = np.searchsorted(t, [g - half, g + half])
        if b - a > 50: base[i] = np.percentile(f[a:b], pct)
    ok = np.isfinite(base)
    return np.interp(t, grid[ok], base[ok])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pmin", type=float); ap.add_argument("pmax", type=float)
    ap.add_argument("--step", type=float, default=0.01)
    ap.add_argument("--mask", type=float, default=None,
                    help="use the masking rule instead: cadences more than MASK below quiet count as hidden")
    a = ap.parse_args()
    t, f = kepler_lc()
    deficit = quiet_level(t, f) - f
    visible = deficit < (a.mask if a.mask is not None else 0.008)
    def status(tc):
        i, j = np.searchsorted(t, [tc - T23 / 2, tc + T23 / 2])
        if visible[i:j].sum() >= 3: return "excluded"
        return "hidden" if j - i >= 3 else "gap"
    Ps = np.round(np.arange(a.pmin, a.pmax + 1e-9, a.step), 4)
    surv = []
    for P in Ps:
        st = [(n, status(T0 + n * P)) for n in range(-12, 0) if 120.0 < T0 + n * P < 1591.5]
        if all(s != "excluded" for _, s in st): surv.append((P, st))
    rows, i = [], 0
    while i < len(surv):
        j = i
        while j + 1 < len(surv) and abs(surv[j + 1][0] - surv[j][0] - a.step) < 1e-6: j += 1
        rows.append((surv[i][0], surv[j][0], surv[i][1])); i = j + 1
    with open("kepler_alias_survivors.csv", "w") as fh:
        fh.write("p_min,p_max,epochs\n")
        for lo, hi, st in rows:
            print(f"{lo:.2f}-{hi:.2f}  {st}")
            fh.write(f"{lo:.2f},{hi:.2f},\"{st}\"\n")

if __name__ == "__main__":
    main()
