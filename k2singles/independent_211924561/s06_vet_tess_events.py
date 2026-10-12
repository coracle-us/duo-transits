"""Zoom on each TESS candidate in every reduction + diagnostics (bkg, centroids, quality/momentum dumps)."""
import os
import numpy as np
from astropy.io import fits
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lcio

CANDS = [(44, 2509.81), (45, 2539.52), (45, 2541.05), (46, 2568.30), (72, 3270.75), (72, 3275.46)]
HW = 1.0  # half-window in days


def binned(t, f, dt):
    e = np.arange(t.min(), t.max() + dt, dt)
    i = np.digitize(t, e)
    tb = np.array([t[i == k].mean() for k in np.unique(i)])
    fb = np.array([np.nanmean(f[i == k]) for k in np.unique(i)])
    return tb, fb


for sec, tc in CANDS:
    rows = []
    f2 = lcio._find(f"tess*-s{sec:04d}-0000000238606535-*-s_lc.fits")
    d = fits.getdata(f2, 1)
    t = d["TIME"]; w = np.isfinite(t) & (np.abs(t - tc) < HW)
    q = d["QUALITY"]
    dumps = t[w & ((q & 32) > 0 | ((q & 1024) > 0))]  # bit6 reaction wheel desat, bit 11 impulsive? keep 32/1024 info
    dumps = t[np.isfinite(t) & ((q & 32) > 0) & (np.abs(t - tc) < HW)]
    qbad = t[np.isfinite(t) & (q != 0) & (np.abs(t - tc) < HW)]
    good = w & (q == 0)
    for col in ["PDCSAP_FLUX", "SAP_FLUX", "SAP_BKG", "MOM_CENTR1", "MOM_CENTR2"]:
        y = d[col][good].astype(float)
        rows.append((f"SPOC2min {col}", t[good], y / np.nanmedian(y) if "FLUX" in col and col != "SAP_BKG" else y - np.nanmedian(y)))
    for col in ["PDCSAP_FLUX", "SAP_FLUX"]:
        try:
            dd = fits.getdata(lcio._find(f"hlsp_tess-spoc_tess_phot_0000000238606535-s{sec:04d}_tess_v1_lc.fits"), 1)
            g = (dd["QUALITY"] == 0) & (np.abs(dd["TIME"] - tc) < HW) & np.isfinite(dd[col])
            rows.append((f"TESS-SPOC {col}", dd["TIME"][g], dd[col][g] / np.nanmedian(dd[col][g])))
        except FileNotFoundError:
            pass
    dq = fits.getdata(lcio._find(f"hlsp_qlp_tess_ffi_s{sec:04d}-0000000238606535_tess_v01_llc.fits"), 1)
    qc = "KSPSAP_FLUX" if "KSPSAP_FLUX" in dq.columns.names else "DET_FLUX"
    for col in [qc, "SAP_FLUX", "SAP_BKG"]:
        g = (dq["QUALITY"] == 0) & (np.abs(dq["TIME"] - tc) < HW) & np.isfinite(dq[col])
        y = dq[col][g].astype(float)
        rows.append((f"QLP {col}", dq["TIME"][g], y / np.nanmedian(y) if col != "SAP_BKG" else y - np.nanmedian(y)))
    dt_ = fits.getdata(lcio._find(f"hlsp_tars_tess_ffi_s{sec:04d}-0000000238606535_tess_v01_lc.fits"), 1)
    g = (np.abs(dt_["TIME"] - tc) < HW) & np.isfinite(dt_["FLUX"])
    rows.append(("TARS FLUX", dt_["TIME"][g], dt_["FLUX"][g] / np.nanmedian(dt_["FLUX"][g])))

    fig, ax = plt.subplots(len(rows), 1, figsize=(9, 1.7 * len(rows)), sharex=True)
    for a, (lab, x, y) in zip(ax, rows):
        a.plot(x, y, ".", ms=2, color="0.6")
        if len(x) > 10:
            xb, yb = binned(x, y, 30 / 1440.)
            a.plot(xb, yb, "o", ms=3, color="C0")
        a.set_ylabel(lab, fontsize=7, rotation=0, ha="right")
        a.axvspan(tc - 0.06, tc + 0.06, color="r", alpha=0.08)
        for dmp in dumps:
            a.axvline(dmp, color="m", lw=0.8)
        if "FLUX" in lab and "BKG" not in lab:
            lo, hi = np.nanpercentile(yb if len(x) > 10 else y, [1, 99]); a.set_ylim(lo - 0.0015, hi + 0.0015)
    ax[0].set_title(f"S{sec} candidate near {tc}  (magenta = momentum dumps; flagged cadences n={len(qbad)})")
    plt.tight_layout()
    plt.savefig(os.path.join(lcio.HERE, f"fig_vet_S{sec}_{tc:.2f}.png"), dpi=90)
    plt.close()
    print(sec, tc, "dumps within 1 d:", np.round(np.unique(np.round(dumps, 3)), 3), " nonzero-q cadences:", len(qbad),
          "  range", (np.round(qbad.min(), 3), np.round(qbad.max(), 3)) if len(qbad) else "")
