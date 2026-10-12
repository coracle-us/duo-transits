"""K2 C0 light curves of EPIC 202126877: SAP, EVEREST, K2SFF, and own large-aperture TPF photometry. Plots full campaign with
predicted primary/secondary times (from the K2 eclipse at BTJD -211.557 and phase 0.33)."""
import numpy as np, glob, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, pandas as pd
from astropy.io import fits
off = 2454833 - 2457000   # BKJD -> BTJD
ev = fits.open(glob.glob("data/k2/**/hlsp_everest*.fits", recursive=True)[0])[1].data
sff = fits.open(glob.glob("data/k2/**/hlsp_k2sff*.fits", recursive=True)[0])[1].data
sap = fits.open(glob.glob("data/k2/**/ktwo*llc.fits", recursive=True)[0])[1].data
tpf = fits.open(glob.glob("data/k2/**/ktwo*lpd-targ.fits.gz", recursive=True)[0])
td = tpf[1].data; img = td["FLUX"]; tt = td["TIME"]; q = td["QUALITY"]
med = np.nanmedian(img[np.isfinite(tt)], axis=0)
print("TPF shape", img.shape)
# own photometry: big aperture (all pixels above 3x background) and bkg from the faintest 30% pixels
bk = np.nanpercentile(med, 30)
ap = med > bk + 0.02*np.nanmax(med)
print("aperture pixels", ap.sum())
bkg = np.array([np.nanmedian(im[med < np.nanpercentile(med, 30)]) for im in img])
fl = np.array([np.nansum(im[ap]) for im in img]) - bkg*ap.sum()
own = pd.DataFrame(dict(btjd=tt+off, flux=fl/np.nanmedian(fl), q=q))
own.to_csv("data/phot/k2_own_tpf.csv", index=False)
series = {"SAP": (sap["TIME"]+off, sap["SAP_FLUX"]), "EVEREST": (ev["TIME"]+off, ev["FLUX"]), "K2SFF": (sff["T"]+off, sff["FCOR"]),
          "own TPF": (own.btjd.values, own.flux.values)}
fig, axs = plt.subplots(len(series), 1, figsize=(15, 11), sharex=True)
tp = -211.557; P = 21.77465
for ax, (k, (t, f)) in zip(axs, series.items()):
    m = np.isfinite(t) & np.isfinite(f); f = f/np.nanmedian(f[m])
    ax.plot(t[m], f[m], "k.", ms=1); ax.set_ylabel(k); ax.set_ylim(0.85, 1.05)
    for n in range(-3, 2):
        ax.axvline(tp + n*P, color="b", alpha=.4)
        for ph in [0.30, 0.33, 0.36]: ax.axvline(tp + (n+ph)*P, color="g", alpha=.3)
plt.tight_layout(); plt.savefig("figures/k2_full_campaign.png", dpi=80)
