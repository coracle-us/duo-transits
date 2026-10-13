"""Step 5: own aperture photometry on all TESScut cutouts (same aperture definition every sector).

Aperture: circle of radius R (px) centred on the target's WCS position; pixel included if its
centre lies within R. Background: per-frame median of 'dark' pixels (> 6 px from the target,
time-median flux below the 50th percentile of those pixels). Outputs data/phot_<sector>.csv
with fluxes for R = 1.5, 2.0, 2.5, 3.0 px.
"""
import os, glob, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(H, "data")
RA, DEC = 97.98012, 25.48207
RADII = [1.5, 2.0, 2.5, 3.0]

files = sorted(glob.glob(os.path.join(D, "tesscut", "*.fits")))
fig, axes = plt.subplots(len(files), 2, figsize=(14, 2.6 * len(files)), gridspec_kw=dict(width_ratios=[1, 5]))
for k, f in enumerate(files):
    h = fits.open(f)
    sector = h[0].header["SECTOR"]
    d = h[1].data
    t = np.array(d["TIME"], float)
    cube = np.array(d["FLUX"], float)
    bkg0 = np.array(d["FLUX_BKG"], float)
    q = np.array(d["QUALITY"], int)
    w = WCS(h[2].header)
    x0, y0 = w.all_world2pix(RA, DEC, 0)
    ny, nx = cube.shape[1:]
    yy, xx = np.mgrid[0:ny, 0:nx]
    rr = np.hypot(xx - x0, yy - y0)
    med = np.nanmedian(cube, axis=0)
    far = rr > 6
    thr = np.nanpercentile(med[far], 50)
    bmask = far & (med < thr)
    bkg = np.nanmedian(cube[:, bmask], axis=1)
    out = dict(time=t, quality=q, bkg=bkg, pipe_bkg_in_ap=np.nansum(bkg0[:, rr <= 2.0], axis=1))
    for R in RADII:
        ap = rr <= R
        fl = np.nansum(cube[:, ap] - bkg[:, None], axis=1)
        out["flux_r%.1f" % R] = fl
        out["npix_r%.1f" % R] = ap.sum()
    # flux-weighted centroid in r<=2.5
    ap = rr <= 2.5
    im = cube[:, ap] - bkg[:, None]
    out["cx"] = np.nansum(im * xx[ap], axis=1) / np.nansum(im, axis=1)
    out["cy"] = np.nansum(im * yy[ap], axis=1) / np.nansum(im, axis=1)
    df = pd.DataFrame(out)
    df = df[np.isfinite(df.time)]
    df.to_csv(os.path.join(D, "phot_s%04d.csv" % sector), index=False)
    print("S%d target px (%.2f, %.2f), npix r2=%d, nbkg=%d, median bkg/px=%.1f, pipeline bkg mean=%.1f, flux r2=%.0f"
          % (sector, x0, y0, (rr <= 2).sum(), bmask.sum(), np.nanmedian(bkg), np.nanmean(bkg0), np.nanmedian(df["flux_r2.0"])))
    ax = axes[k, 0]
    ax.imshow(np.log10(np.clip(med, 1, None)), origin="lower")
    ax.add_patch(plt.Circle((x0, y0), 2.0, fill=False, color="r"))
    ax.contour(bmask, levels=[0.5], colors="w", linewidths=0.5)
    ax.set_title("S%d" % sector)
    ax = axes[k, 1]
    fl = df["flux_r2.0"].values
    ax.plot(df.time, fl / np.nanmedian(fl), "k.", ms=1)
    bad = df.quality.values != 0
    ax.plot(df.time[bad], (fl / np.nanmedian(fl))[bad], "r.", ms=1)
    ax.set_ylim(0.6, 1.1)
plt.tight_layout()
plt.savefig(os.path.join(H, "fig_tesscut_phot.png"), dpi=70)
