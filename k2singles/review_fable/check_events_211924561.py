"""Independent look at the four claimed events of EPIC 211924561 / TIC 238606535 using only the
pipeline products already on disk (no download): K2 mission PDCSAP + SAP from
k2singles/cache/211924561/dl/, TESS SPOC 2-min PDCSAP/SAP/background/centroids from
data/raw/lightcurves/238606535/.  Local quadratic detrend, box depth over the central 70 % of a 3-h
window, and the same statistic at the two 70.247-d-only epochs that fall just outside data (to show
they are untested) and at the m=3/m=4 excluding epochs (to show they are flat).

Output: k2singles/review_fable/events_211924561.png and printed table.
Run:    .venv/Scripts/python.exe k2singles/review_fable/check_events_211924561.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
D = 3.0 / 24


def k2(camp):
    fn = ROOT / "k2singles" / "cache" / "211924561" / "dl" / f"ktwo211924561-c{camp:02d}_llc.fits"
    with fits.open(fn) as h:
        d = h[1].data
        t = np.asarray(d["TIME"], float) + 2454833 - 2457000
        return dict(t=t, pdc=np.asarray(d["PDCSAP_FLUX"], float), sap=np.asarray(d["SAP_FLUX"], float),
                    q=(np.asarray(d["SAP_QUALITY"], int) & 1130799), bkg=np.full(len(t), np.nan),
                    c1=np.asarray(d["MOM_CENTR1"], float), c2=np.asarray(d["MOM_CENTR2"], float), cad=0.0204)


def tess(sec):
    fn = ROOT / "data" / "raw" / "lightcurves" / "238606535" / f"s{sec:03d}_SPOC_120s.csv.gz"
    d = pd.read_csv(fn)
    return dict(t=d.time.values, pdc=d.pdcsap_flux.values, sap=d.sap_flux.values, q=d.quality.values.astype(int),
                bkg=d.sap_bkg.values, c1=d.mom_centr1.values, c2=d.mom_centr2.values, cad=2 / 1440)


def local(dat, T, col="pdc", half=0.6):
    t, f, q = dat["t"], dat[col], dat["q"]
    m = (np.abs(t - T) < half) & np.isfinite(t) & np.isfinite(f) & (q == 0)
    if m.sum() < 10:
        return None
    tt, ff = t[m], f[m]
    oot = np.abs(tt - T) > D / 2 + 0.03
    if oot.sum() < 8:
        return None
    base = np.polyval(np.polyfit(tt[oot] - T, ff[oot], 2), tt - T)
    r = ff / base - 1
    inn = np.abs(tt - T) < 0.7 * D / 2
    sig = 1.4826 * np.median(np.abs(r[oot] - np.median(r[oot])))
    n_in = inn.sum()
    depth = -np.mean(r[inn]) if n_in else np.nan
    return dict(tt=tt - T, r=r, depth_ppt=1e3 * depth, err_ppt=1e3 * sig / np.sqrt(max(n_in, 1)), n_in=int(n_in),
                cover=n_in * dat["cad"] / (0.7 * D))


events = [("K2 C5", "K5", 179.8978), ("K2 C16", "K16", 1163.3554), ("TESS S46", "S46", 2568.2983), ("TESS S72", "S72", 3270.7646)]
P1, P2, T0 = 140.49401, 70.24701, 179.89774
extra = [("m=3 C18", "K18", 1257.0185), ("m=3 S44", "S44", 2521.4646), ("m=4 C5", "K5", 144.7742), ("m=4 C16", "K16", 1128.2323),
         ("m=4 S45", "S45", 2533.1724), ("P70 C16 start?", "K16", T0 + 13 * P2), ("P70 C18 end?", "K18", T0 + 16 * P2),
         ("P70 S44 start?", "S44", T0 + 33 * P2)]
rows = []
fig, axs = plt.subplots(3, 4, figsize=(16, 9))
for ax, (lab, s, T) in zip(axs.flat, events + extra):
    dat = k2(int(s[1:])) if s[0] == "K" else tess(int(s[1:]))
    for col, c, ms in (("sap", "0.7", 3), ("pdc", "C0", 4)):
        L = local(dat, T, col)
        if L is None:
            ax.text(0.5, 0.5, "no data at epoch", transform=ax.transAxes, ha="center")
            if col == "pdc":
                rows.append(dict(event=lab, T=T, depth_ppt=np.nan, note="no data"))
            continue
        ax.plot(L["tt"] * 24, L["r"] * 1e3, ".", color=c, ms=ms, label=col)
        if col == "pdc":
            rows.append(dict(event=lab, T=round(T, 4), depth_ppt=round(L["depth_ppt"], 2), err_ppt=round(L["err_ppt"], 2),
                             n_in=L["n_in"], cover=round(L["cover"], 2)))
    ax.axvspan(-1.5, 1.5, color="r", alpha=0.08)
    ax.set_ylim(-3.5, 2); ax.set_xlim(-14, 14)
    ax.set_title(f"{lab}  T={T:.3f}", fontsize=9)
    ax.set_xlabel("h from T"); ax.set_ylabel("ppt")
axs.flat[0].legend(fontsize=7)
fig.suptitle("EPIC 211924561: SAP (grey) and PDCSAP (blue) at the 4 claimed events (row 1), the m=3/m=4 excluding epochs (row 2) "
             "and the 70.247-d-only epochs nearest the data (row 3)", fontsize=9)
fig.tight_layout()
fig.savefig(OUT / "events_211924561.png", dpi=90)
tab = pd.DataFrame(rows)
print(tab.to_string(index=False))

# TESS background / centroid behaviour during the two TESS events
print("\nTESS SAP background and centroid, in-transit vs out-of-transit (same +-0.6 d window):")
for lab, s, T in events[2:]:
    dat = tess(int(s[1:]))
    m = (np.abs(dat["t"] - T) < 0.6) & (dat["q"] == 0)
    inn = np.abs(dat["t"] - T) < 0.7 * D / 2
    for key in ("bkg", "c1", "c2"):
        v = dat[key]
        a, b = np.nanmean(v[m & inn]), np.nanmean(v[m & ~inn])
        sd = np.nanstd(v[m & ~inn]) / np.sqrt((m & inn).sum())
        print(f"  {lab} {key}: in-out = {a - b:+.4f} (+-{sd:.4f})  [oot mean {b:.3f}]")
