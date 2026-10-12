"""Cadence-level table around the C5 event: flags and flux in every reduction (no quality masking)."""
import numpy as np
from astropy.io import fits
import lcio

T0 = 179.89
sap = fits.getdata(lcio._find("ktwo211924561-c05_llc.fits"), 1)
ev = fits.getdata(lcio._find("hlsp_everest_k2_llc_211924561-c05_kepler_v2.0_lc.fits"), 1)
sff = fits.getdata(lcio._find("hlsp_k2sff_k2_lightcurve_211924561-c05_kepler_v1_llc.fits"), 1)
k2sc = fits.getdata(lcio._find("hlsp_k2sc_k2_llc_211924561-c05_kepler_v2_lc.fits"), 1)

t = sap["TIME"] - 2167
w = np.where(np.abs(t - T0) < 0.45)[0]
sff_t = sff["T"] - 2167
ev_t = ev["TIME"] - 2167
k_t = k2sc["time"] - 2167
med_sap = np.nanmedian(sap["SAP_FLUX"][w])
med_ev = np.nanmedian(ev["FLUX"][(np.abs(ev_t - T0) < 1.5)])
med_sff = np.nanmedian(sff["FCOR"][(np.abs(sff_t - T0) < 1.5)])
med_k = np.nanmedian(k2sc["flux"][(np.abs(k_t - T0) < 1.5)])
print(" cad      BTJD     SAPq      SAP[ppt]  bkg    cen1     cen2     pos1    pos2  | EVq   EVEREST  | K2SFF  mov | K2SC q mflag")
for i in range(w[0], w[-1] + 1):
    c = sap["CADENCENO"][i]
    j = np.where(ev["CADN"] == c)[0]
    k = np.where(sff["CADENCENO"] == c)[0]
    m = np.where(k2sc["cadence"] == c)[0]
    s = "%6d %9.4f %8d %8.2f %7.1f %8.3f %8.3f %7.3f %7.3f" % (
        c, t[i], sap["SAP_QUALITY"][i], 1e3 * (sap["SAP_FLUX"][i] / med_sap - 1), sap["SAP_BKG"][i],
        sap["MOM_CENTR1"][i], sap["MOM_CENTR2"][i], sap["POS_CORR1"][i], sap["POS_CORR2"][i])
    s += " | %7d %7.2f" % ((ev["QUALITY"][j[0]], 1e3 * (ev["FLUX"][j[0]] / med_ev - 1)) if len(j) else (-1, np.nan))
    s += " | %7.2f %3d" % ((1e3 * (sff["FCOR"][k[0]] / med_sff - 1), sff["MOVING"][k[0]]) if len(k) else (np.nan, -1))
    s += " | %7.2f %d %d" % ((1e3 * (k2sc["flux"][m[0]] / med_k - 1), k2sc["quality"][m[0]], k2sc["mflags"][m[0]]) if len(m) else (np.nan, -1, -1))
    print(s)
