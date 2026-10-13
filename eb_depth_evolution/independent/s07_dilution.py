"""Step 7: estimate the third-light fraction in our r=2.0 px TESScut aperture.

Fit each sector's median image with a Gaussian-PSF scene model: stars at TIC positions
(Tmag < 17.5 within 14 px), relative fluxes fixed by Tmag, free overall scale, PSF sigma,
constant sky.  Then:
  contam = sum(non-target flux in aperture) / (target flux in aperture)
  bias   = (our dark-pixel background estimate - fitted sky) * Npix  (over-subtraction)
  l3_eff = (contam_flux - bias) / measured aperture flux      (signed third light)
"""
import os, glob, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
from scipy.special import erf

H = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(H, "data")
RA, DEC = 97.98012, 25.48207
tic = pd.read_csv(os.path.join(D, "tic_neighbours.csv"))
# need a wider list for stars up to 14 px (~300")
from astroquery.mast import Catalogs
from astropy.coordinates import SkyCoord
import astropy.units as u
cache = os.path.join(D, "tic_wide.csv")
if not os.path.exists(cache):
    t = Catalogs.query_region(SkyCoord(RA, DEC, unit="deg"), radius=330 * u.arcsec, catalog="TIC")
    t = t[["ID", "ra", "dec", "Tmag"]].to_pandas()
    t.to_csv(cache, index=False)
tic = pd.read_csv(cache)
tic = tic[np.isfinite(tic.Tmag) & (tic.Tmag < 17.5)].reset_index(drop=True)
itgt = int(np.where(tic.ID == 86491361)[0][0])

def pix_gauss(xx, yy, x0, y0, s):
    fx = 0.5 * (erf((xx + 0.5 - x0) / (np.sqrt(2) * s)) - erf((xx - 0.5 - x0) / (np.sqrt(2) * s)))
    fy = 0.5 * (erf((yy + 0.5 - y0) / (np.sqrt(2) * s)) - erf((yy - 0.5 - y0) / (np.sqrt(2) * s)))
    return fx * fy

res = []
for f in sorted(glob.glob(os.path.join(D, "tesscut", "*.fits"))):
    h = fits.open(f)
    sec = h[0].header["SECTOR"]
    cube = h[1].data["FLUX"]; q = h[1].data["QUALITY"]
    med = np.nanmedian(cube[q == 0], axis=0)
    w = WCS(h[2].header)
    xs, ys = w.all_world2pix(tic.ra.values, tic.dec.values, 0)
    ny, nx = med.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    inn = (xs > -3) & (xs < nx + 2) & (ys > -3) & (ys < ny + 2)
    rel = 10 ** (-0.4 * (tic.Tmag.values - tic.Tmag.values[itgt]))
    # restrict fit to central 17x17 region to limit edge effects
    x0, y0 = xs[itgt], ys[itgt]
    reg = (np.abs(xx - x0) <= 8) & (np.abs(yy - y0) <= 8)

    def model(p, return_parts=False):
        scale, s, sky, dx, dy = p
        img = np.full(med.shape, sky)
        parts = {}
        for j in np.where(inn)[0]:
            g = scale * rel[j] * pix_gauss(xx, yy, xs[j] + dx, ys[j] + dy, s)
            img += g
            if return_parts:
                parts[j] = g
        return (img, parts) if return_parts else img

    p0 = [np.nansum(med[np.hypot(xx - x0, yy - y0) <= 2.5]) * 1.2, 0.8, np.nanpercentile(med, 10), 0, 0]
    r = least_squares(lambda p: ((model(p) - med)[reg]) / np.sqrt(np.abs(med[reg]) + 50), p0,
                      bounds=([0, 0.3, -1e4, -1, -1], [1e7, 3, 1e4, 1, 1]))
    img, parts = model(r.x, True)
    rr = np.hypot(xx - x0, yy - y0)
    ap = rr <= 2.0
    tgt = parts[itgt][ap].sum()
    oth = sum(parts[j][ap].sum() for j in parts if j != itgt)
    # our background estimator
    far = rr > 6
    thr = np.nanpercentile(med[far], 50)
    bmask = far & (med < thr)
    bk_ours = np.nanmedian(med[bmask])
    bias = (bk_ours - r.x[2]) * ap.sum()
    meas = np.nansum(med[ap]) - bk_ours * ap.sum()
    l3_eff = (oth - bias) / meas
    resid_rms = np.std((med - img)[reg]) / np.median(med[reg])
    print("S%-5d sigma=%.2f px sky=%.1f ours=%.1f  target-in-ap=%.0f other-in-ap=%.0f (contam %.3f)"
          "  bkg-bias flux=%.0f  measured=%.0f  l3_eff=%.3f  model/meas tgt=%.3f"
          % (sec, r.x[1], r.x[2], bk_ours, tgt, oth, oth / (tgt + oth), bias, meas, l3_eff, (tgt + oth - bias) / meas))
    res.append(dict(sector=sec, psf_sigma=r.x[1], sky_fit=r.x[2], sky_ours=bk_ours, tgt_in_ap=tgt, other_in_ap=oth,
                    contam_frac=oth / (tgt + oth), bkg_bias_flux=bias, measured=meas, l3_eff=l3_eff))
pd.DataFrame(res).to_csv(os.path.join(H, "dilution.csv"), index=False)
