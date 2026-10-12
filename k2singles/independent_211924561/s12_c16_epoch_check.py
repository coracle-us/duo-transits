"""Close look at K2 C16 around BTJD 1163.36 (the epoch predicted by every P = 140.494/m): all reductions and flags."""
import os
import numpy as np
from astropy.io import fits
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lcio
from lightkurve.utils import KeplerQualityFlags as K

TP = 1163.356
ev = fits.getdata(lcio._find("hlsp_everest_k2_llc_211924561-c16_kepler_v2.0_lc.fits"), 1)
sap = fits.getdata(lcio._find("ktwo211924561-c16_llc.fits"), 1)
sff = fits.getdata(lcio._find("hlsp_k2sff_k2_lightcurve_211924561-c16_kepler_v1_llc.fits"), 1)
tev = ev["TIME"] - 2167; tsap = sap["TIME"] - 2167; tsff = sff["T"] - 2167
w = np.abs(tsap - TP) < 0.5
print("K2SFF cadences within +-1 d:", np.sum(np.abs(tsff - TP) < 1.0), " nearest K2SFF times:", np.round(np.sort(tsff[np.abs(tsff - TP) < 2.5]), 3)[[0, -1]] if np.any(np.abs(tsff - TP) < 2.5) else None)
gap = tsff[(tsff > TP - 3) & (tsff < TP + 3)]
dg = np.diff(gap); ig = np.argmax(dg); print("K2SFF largest gap near epoch: %.3f -> %.3f" % (gap[ig], gap[ig + 1]))
medev = np.nanmedian(ev["FLUX"][np.abs(tev - TP) < 1.5])
meds = np.nanmedian(sap["SAP_FLUX"][w])
print("  BTJD     SAPq  flags                                    SAP[ppt]  EVq      EVEREST[ppt]")
for i in np.where(w)[0]:
    c = sap["CADENCENO"][i]
    j = np.where(ev["CADN"] == c)[0]
    print("%9.4f %8d %-45s %7.2f | %8d %7.2f" % (tsap[i], sap["SAP_QUALITY"][i], ",".join(x[:12] for x in K.decode(sap["SAP_QUALITY"][i]))[:45],
          1e3 * (sap["SAP_FLUX"][i] / meds - 1), ev["QUALITY"][j[0]] if len(j) else -1, 1e3 * (ev["FLUX"][j[0]] / medev - 1) if len(j) else np.nan))

fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
ww = np.abs(tev - TP) < 3
good = ww & (ev["QUALITY"] == 0)
ax[0].plot(tev[ww], 1e3 * (ev["FLUX"][ww] / medev - 1), ".", color="0.7", label="EVEREST all")
ax[0].plot(tev[good], 1e3 * (ev["FLUX"][good] / medev - 1), "o", ms=3, label="EVEREST q=0")
ax[0].set_ylim(-4, 2); ax[0].legend(fontsize=7)
ws = np.abs(tsff - TP) < 3
ax[1].plot(tsff[ws], 1e3 * (sff["FCOR"][ws] / np.nanmedian(sff["FCOR"][ws]) - 1), "o", ms=3); ax[1].set_ylim(-4, 2); ax[1].set_ylabel("K2SFF")
w3 = np.abs(tsap - TP) < 3
ax[2].plot(tsap[w3], sap["SAP_FLUX"][w3] / np.nanmedian(sap["SAP_FLUX"][w3]), ".", ms=3); ax[2].set_ylabel("SAP")
for a in ax:
    a.axvspan(TP - 0.062, TP + 0.062, color="r", alpha=0.15)
plt.tight_layout(); plt.savefig(os.path.join(lcio.HERE, "fig_c16_predicted_epoch.png"), dpi=100)
