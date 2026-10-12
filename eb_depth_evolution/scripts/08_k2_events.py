"""K2 C0: fit the single primary (EVEREST, K2SFF, own TPF after local detrend) and set limits on the secondary by scanning
phase 0.20-0.45 relative to the K2 primary with a fixed-shape template (generalized Gaussian, w = 0.114 d (2021 TESS secondary width)
scaled, p=2.8)."""
import numpy as np, glob, pandas as pd, sys
from astropy.io import fits
sys.path.insert(0, "scripts"); from common import fit_event, gg
off = 2454833 - 2457000
ev = fits.open(glob.glob("data/k2/**/hlsp_everest*.fits", recursive=True)[0])[1].data
sff = fits.open(glob.glob("data/k2/**/hlsp_k2sff*.fits", recursive=True)[0])[1].data
own = pd.read_csv("data/phot/k2_own_tpf.csv")
sets = {"EVEREST": (ev["TIME"]+off, ev["FLUX"], ev["QUALITY"]), "K2SFF": (sff["T"]+off, sff["FCOR"], np.zeros(len(sff))),
        "own TPF": (own.btjd.values, own.flux.values, own.q.values)}
rows = []
for k, (t, f, q) in sets.items():
    m = np.isfinite(t) & np.isfinite(f) & (t > -232)
    t, f = t[m].astype(float), f[m].astype(float); f = f/np.median(f)
    e = np.full(len(f), np.std(np.diff(f))/np.sqrt(2))
    if k == "own TPF":   # remove thruster-fire outliers by sigma clipping a running median
        rm = pd.Series(f).rolling(9, center=True, min_periods=1).median().values
        good = np.abs(f - rm) < 4*np.std(f - rm); t, f, e = t[good], f[good], e[good]
    r = fit_event(t, f, e, -211.557, win=0.8)
    rows.append(dict(product=k, kind="P", tc=r["tc"], e_tc=r["e_tc"], depth=r["D"], e_depth=r["e_D"], w=r["w"], p=r["p"], beta=r["beta"], rms=r["rms"], n=r["n"]))
    # integration-time effect: K2 long cadence 29.4 min. Generalised Gaussian fit to binned data underestimates true depth slightly;
    # estimate by convolving best model and refitting
    # secondary scan
    if k in ("EVEREST", "K2SFF"):
        best = []
        for ph in np.arange(0.20, 0.46, 0.0025):
            for tsec in [-211.557 + ph*21.77465, -211.557 - (1-ph)*21.77465]:
                win = np.abs(t - tsec) < 0.8
                if win.sum() < 30: continue
                tw, fw = t[win], f[win]
                for wsec in (0.114, 0.145):
                    s = 1 - gg(tw, tsec, 1.0, wsec, 2.7)        # unit-depth template
                    A = np.vstack([np.ones_like(tw), tw - tsec, -s]).T
                    c, *_ = np.linalg.lstsq(A, fw, rcond=None)
                    res = fw - A @ c; sig = np.std(res)
                    cov = np.linalg.inv(A.T @ A)*sig**2
                    best.append((ph, tsec, wsec, c[2], np.sqrt(cov[2, 2]), win.sum()))
        b = pd.DataFrame(best, columns=["phase", "t", "w", "depth", "err", "n"])
        b.to_csv(f"data/k2_secondary_scan_{k}.csv", index=False)
        for ts in sorted(set(np.round(b.t.values // 5))):
            pass
        b["epoch"] = np.where(b.t < -215, "E-1 (-226)", "E0 (-204)")
        for epn, g in b.groupby("epoch"):
            i = g.depth.idxmax()
            print(f"{k} secondary {epn}: phase 0.20-0.45 max depth {g.depth.max()*1e6:.0f} ppm at phase {g.loc[i,'phase']:.3f}; "
                  f"typical err {g.err.median()*1e6:.0f} ppm; 3-sigma upper limit over scan {(g.depth.max()+3*g.err.median())*1e6:.0f} ppm")
print(pd.DataFrame(rows).to_string(index=False))
pd.DataFrame(rows).to_csv("data/k2_primary_fits.csv", index=False)
