"""Step 6: collect every eclipse window from every data source into one table.

Sources:
  own  : own TESScut aperture photometry (r = 2.0 px; other radii kept as columns)
  spoc : TESS-SPOC PDCSAP (S43-45)
  qlp  : QLP SAP_FLUX (S43-45, S71-72) (raw SAP; we normalise locally ourselves)
  k2   : K2 PDCSAP, K2SFF FCOR, EVEREST FLUX (C0)
All times in BTJD = BJD_TDB - 2457000 (K2 BKJD - 2167).
Output: data/windows.csv and fig_windows.png
"""
import os, glob, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(H, "data")
BADBITS = 1 | 2 | 4 | 8 | 16 | 32 | 128 | 1024

P0, T00, DSEC = 21.7765, 2488.60, 7.20   # rough, refined later
HW = 1.0                                  # half window (d)

def season(t):
    if t < 0: return 2014
    if t < 3000: return 2021
    if t < 4000: return 2023
    return 2026

rows = []
def add(src, t, f, sector, texp, extra=None):
    t = np.asarray(t, float); f = np.asarray(f, float)
    ok = np.isfinite(t) & np.isfinite(f)
    t, f = t[ok], f[ok]
    ex = {k: np.asarray(v)[ok] for k, v in (extra or {}).items()}
    nmin, nmax = int(np.floor((t.min() - T00) / P0)) - 1, int(np.ceil((t.max() - T00) / P0)) + 1
    for n in range(nmin, nmax + 1):
        for typ, tc in (("P", T00 + n * P0), ("S", T00 + n * P0 + DSEC)):
            m = np.abs(t - tc) < HW
            if m.sum() < 10:
                continue
            for j in np.where(m)[0]:
                r = dict(src=src, sector=sector, texp=texp, season=season(t[j]), epoch=n, type=typ,
                         tpred=tc, time=t[j], flux=f[j])
                for k, v in ex.items():
                    r[k] = v[j]
                rows.append(r)

# own TESScut photometry
for fpath in sorted(glob.glob(os.path.join(D, "phot_s*.csv"))):
    d = pd.read_csv(fpath)
    sec = int(os.path.basename(fpath)[6:10])
    d = d[(d.quality & BADBITS) == 0]
    texp = 600.0 if sec < 56 else 200.0
    add("own", d.time, d["flux_r2.0"], sec, texp,
        {"f15": d["flux_r1.5"], "f25": d["flux_r2.5"], "f30": d["flux_r3.0"], "bkg": d.bkg})

for fpath in sorted(glob.glob(os.path.join(D, "mastDownload", "HLSP", "*", "*.fits"))):
    name = os.path.basename(fpath)
    h = fits.open(fpath); dd = h[1].data
    if "tess-spoc" in name:
        sec = h[0].header["SECTOR"]
        q = dd["QUALITY"]; m = (q & BADBITS) == 0
        add("spoc", dd["TIME"][m], dd["PDCSAP_FLUX"][m], sec, 600.0)
    elif "qlp" in name:
        sec = h[0].header["SECTOR"]
        q = dd["QUALITY"]; m = q == 0
        add("qlp", dd["TIME"][m], dd["SAP_FLUX"][m], sec, 600.0 if sec < 56 else 200.0)
    elif "k2sff" in name:
        add("k2sff", dd["T"] - 2167.0, dd["FCOR"], 0, 1765.5)
    elif "everest" in name:
        m = dd["QUALITY"] == 0
        add("everest", dd["TIME"][m] - 2167.0, dd["FLUX"][m], 0, 1765.5)
for fpath in glob.glob(os.path.join(D, "mastDownload", "K2", "*", "ktwo*llc.fits")):
    dd = fits.open(fpath)[1].data
    m = dd["SAP_QUALITY"] == 0
    add("k2pdc", dd["TIME"][m] - 2167.0, dd["PDCSAP_FLUX"][m], 0, 1765.5)

W = pd.DataFrame(rows)
W.to_csv(os.path.join(D, "windows.csv"), index=False)
g = W.groupby(["src", "season", "epoch", "type"])
print(g.size().to_string())

# plot all windows
keys = sorted(g.groups.keys(), key=lambda k: (k[2], k[3], k[0]))
ev = sorted(set((k[2], k[3]) for k in keys))
fig, axes = plt.subplots(int(np.ceil(len(ev) / 4)), 4, figsize=(18, 3.2 * np.ceil(len(ev) / 4)))
axes = axes.ravel()
for ax, (n, typ) in zip(axes, ev):
    for k in keys:
        if (k[2], k[3]) != (n, typ):
            continue
        d = g.get_group(k)
        oot = np.abs(d.time - d.tpred) > 0.5
        norm = np.nanmedian(d.flux[oot]) if oot.sum() > 5 else np.nanmax(d.flux)
        ax.plot(d.time - d.tpred, d.flux / norm, ".", ms=2, label=k[0])
    ax.set_title("epoch %d %s (tpred %.2f)" % (n, typ, d.tpred.iloc[0]), fontsize=9)
    ax.legend(fontsize=6)
for ax in axes[len(ev):]:
    ax.axis("off")
plt.tight_layout()
plt.savefig(os.path.join(H, "fig_windows.png"), dpi=65)
