"""Scan WASP (2004-2008) and KELT (2006-2008) for primary / secondary eclipses at any time offset delta in [-1.5, +1.5] d from the
linear TESS ephemeris (P=21.77465, T0=2488.59020; secondary at phase 0.33).  At each delta, fit flux = c_night + D*template, where
template is a unit-depth generalized-Gaussian (w=0.075 d, p=2.7, the K2 2014 primary shape; narrow = conservative), with a free
zero point per night.  Report D, sigma_D, and number of in-eclipse points (|t-tc|<w) for each delta."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, glob
P, T0 = 21.77465, 2488.59020
def kelt():
    out = []
    for f in glob.glob("data/kelt/*_raw_lc.tbl"):
        a = np.array([l.split() for l in open(f) if l.strip() and l.strip()[0].isdigit()], float)
        out.append(pd.DataFrame(dict(hjd=a[:, 0], mag=a[:, 1], src="east" if "east" in f else "west")))
    d = pd.concat(out); d["dm"] = d.mag - d.groupby("src").mag.transform("median")
    mad = 1.48*np.median(np.abs(d.dm)); d = d[np.abs(d.dm) < 0.15]    # remove gross outliers only (>15%)
    d = d[d.dm > -5*mad]
    d["flux"] = 10**(-0.4*d.dm); d["btjd"] = d.hjd - 2457000; return d[["btjd", "flux", "src"]]
def wasp():
    d = pd.read_csv("data/wasp_clean2.csv"); return d[["btjd", "flux"]].assign(src="cam" + d.cam.astype(str))
w, p = 0.075, 2.7
res = []
for name, d in [("WASP", wasp()), ("KELT", kelt())]:
    d = d.sort_values("btjd").reset_index(drop=True)
    d["night"] = np.floor(d.btjd + 0.5).astype(int)   # local night (site longitudes ~ -18 and -112 deg; noon break)
    sig = 1.48*np.median(np.abs(d.flux - d.groupby("night").flux.transform("median")))
    for kind, ph in [("P", 0.0), ("S", 0.33)]:
        for delta in np.arange(-1.5, 1.5001, 0.005):
            E = np.round((d.btjd - T0 - ph*P - delta)/P)
            tc = T0 + (E + ph)*P + delta
            tmpl = np.exp(-np.abs((d.btjd - tc)/w)**p)
            nin = int((tmpl > 0.5).sum())
            if nin < 3: res.append((name, kind, delta, np.nan, np.nan, nin, 0)); continue
            # only nights with in-eclipse points matter for D; others only constrain their own zero point
            nights_in = d.night[tmpl > 0.05].unique()
            m = d.night.isin(nights_in)
            dd, tt = d[m], tmpl[m]
            # design matrix: one column per night + template
            nid = pd.factorize(dd.night)[0]
            A = np.zeros((len(dd), nid.max() + 2)); A[np.arange(len(dd)), nid] = 1; A[:, -1] = -tt
            c, *_ = np.linalg.lstsq(A, dd.flux.values, rcond=None)
            r = dd.flux.values - A @ c; s = max(sig, 1.48*np.median(np.abs(r)))
            cov = np.linalg.pinv(A.T @ A)*s**2
            res.append((name, kind, delta, c[-1], np.sqrt(cov[-1, -1]), nin, len(nights_in)))
R = pd.DataFrame(res, columns=["survey", "kind", "delta_d", "depth", "err", "n_in", "n_nights"])
R.to_csv("data/ground_eclipse_scan.csv", index=False)
fig, axs = plt.subplots(2, 2, figsize=(14, 7), sharex=True)
for i, s in enumerate(["WASP", "KELT"]):
    for j, k in enumerate(["P", "S"]):
        g = R[(R.survey == s) & (R.kind == k)]; ax = axs[i, j]
        ax.errorbar(g.delta_d, g.depth, g.err, fmt=".", ms=2, lw=0.5)
        ax.axhline(0, color="k"); ax.set_ylim(-0.05, 0.25); ax.set_title(f"{s} {k}: fitted depth vs time offset from TESS linear ephemeris")
        ax.set_xlabel("delta (d)")
        for lev, col in [(0.061, "orange"), (0.20, "r")]: ax.axhline(lev, color=col, ls=":")
plt.tight_layout(); plt.savefig("figures/ground_eclipse_scan.png", dpi=80)
for (s, k), g in R.groupby(["survey", "kind"]):
    cov = g.err.notna()
    print(f"{s} {k}: covered offsets {cov.mean()*100:.0f}% of [-1.5,1.5] d; max depth/err {np.nanmax(g.depth/g.err):.1f} at delta={g.delta_d.values[np.nanargmax(g.depth/g.err)]:.3f} "
          f"(D={g.depth.values[np.nanargmax(g.depth/g.err)]:.4f}); median 3-sigma limit {np.nanmedian(g.depth+3*g.err):.4f}")
    for lo, hi in [(-0.3, 0.3), (-0.6, 0.6), (-1.0, 1.0)]:
        h = g[(g.delta_d >= lo) & (g.delta_d <= hi)]
        print(f"     window [{lo},{hi}] d: covered {h.err.notna().mean()*100:.0f}%; worst 3-sigma limit among covered {np.nanmax(h.depth+3*h.err):.4f}")
