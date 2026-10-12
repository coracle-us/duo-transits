"""Own fixed-aperture photometry on TESScut 21x21 cutouts, all sectors, several aperture radii.
Background: per-cadence median of pixels > 6 px from the star, after masking the brightest 10% of them.
Aperture centre: flux-weighted centroid of the median image in a 5x5 box around the peak (fixed per sector).
Outputs data/phot/tess_<sector>.csv with columns time, q, bkg, f_r1.0, f_r1.5, f_r2.0, f_r2.5, cx, cy (per-cadence centroid in r=2 ap)."""
import glob, os, numpy as np, pandas as pd
from astropy.io import fits
os.makedirs("data/phot", exist_ok=True)
radii = [1.0, 1.5, 2.0, 2.5]
for f in sorted(glob.glob("data/tesscut/tesscut/*.fits")):
    sec = os.path.basename(f).split("-")[1]
    h = fits.open(f); d = h[1].data
    t, fl, q = d["TIME"], d["FLUX"].astype(float), d["QUALITY"]
    good = np.isfinite(t) & (np.nansum(fl, axis=(1, 2)) != 0)
    t, fl, q = t[good], fl[good], q[good]
    ny, nx = fl.shape[1:]
    med = np.nanmedian(fl[q == 0], axis=0)
    y0, x0 = 10, 10
    sub = med[y0-2:y0+3, x0-2:x0+3]; sub = sub - np.nanmin(sub)
    yy, xx = np.mgrid[-2:3, -2:3]
    cy = y0 + np.nansum(yy*sub)/np.nansum(sub); cx = x0 + np.nansum(xx*sub)/np.nansum(sub)
    Y, X = np.mgrid[0:ny, 0:nx]
    R = np.hypot(Y-cy, X-cx)
    bmask = R > 6
    bvals = fl[:, bmask]
    thr = np.nanpercentile(med[bmask], 90)
    keep = med[bmask] < thr
    bkg = np.nanmedian(bvals[:, keep], axis=1)
    out = {"time": t, "q": q, "bkg": bkg}
    for r in radii:
        ap = R <= r
        out[f"f_r{r}"] = np.nansum(fl[:, ap], axis=1) - bkg*ap.sum()
    ap2 = R <= 2.0
    fs = fl - bkg[:, None, None]
    w = np.where(ap2, fs, 0)
    out["cx"] = (w*X).sum(axis=(1,2))/w.sum(axis=(1,2)); out["cy"] = (w*Y).sum(axis=(1,2))/w.sum(axis=(1,2))
    pd.DataFrame(out).to_csv(f"data/phot/tess_s{int(sec[1:])}.csv", index=False)
    print(sec, f"centre=({cx:.2f},{cy:.2f})", len(t), "npix:", [int((R<=r).sum()) for r in radii])
