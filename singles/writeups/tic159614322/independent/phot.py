"""Own aperture photometry on TESScut FFI cutouts with a local background from pixels away from the star."""
import glob, os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

RA, DEC = 220.00977, 65.484041
DIRS = [r"C:\git\duo-transits\data\raw\tesscut\159614322",
        r"C:\git\duo-transits\data\raw\tesscut\159614322\tesscut",
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "tesscut_dl")]


def find_cutout(sector):
    for d in DIRS:
        g = glob.glob(os.path.join(d, f"tess-s{sector:04d}-*astrocut.fits"))
        if g:
            return g[0]
    return None


def photometry(sector, radii=(1.0, 1.5, 2.0, 2.5, 3.0), bkg_rmin=4.5):
    fn = find_cutout(sector)
    h = fits.open(fn)
    d = h[1].data
    w = WCS(h[2].header)
    x0, y0 = w.all_world2pix([[RA, DEC]], 0)[0]  # column (x), row (y) zero-based
    t = np.asarray(d["TIME"], dtype=np.float64)  # native byte order (batman needs it)
    cube = d["FLUX"].astype(float)  # (n, 15, 15) [row, col]
    q = d["QUALITY"]
    ny, nx = cube.shape[1:]
    yy, xx = np.mgrid[0:ny, 0:nx]
    r = np.hypot(xx - x0, yy - y0)
    med_img = np.nanmedian(cube, axis=0)
    # background pixels: far from the star and not bright (exclude pixels above 60th percentile of far pixels)
    far = r > bkg_rmin
    thr = np.nanpercentile(med_img[far], 60)
    bmask = far & (med_img < thr)
    bkg = np.nanmedian(cube[:, bmask], axis=1)
    sub = cube - bkg[:, None, None]
    out = dict(time=t, quality=q, bkg_per_pix=bkg, x0=x0, y0=y0, nbkg=int(bmask.sum()), file=os.path.basename(fn),
               med_img=med_img, bmask=bmask)
    for rad in radii:
        ap = r <= rad
        out[f"ap{rad}"] = np.nansum(sub[:, ap], axis=1)
        out[f"npix{rad}"] = int(ap.sum())
    # flux-weighted centroid within 2.5 px
    cm = r <= 2.5
    wsum = np.nansum(sub[:, cm], axis=1)
    out["cx"] = np.nansum(sub[:, cm] * xx[cm], axis=1) / wsum
    out["cy"] = np.nansum(sub[:, cm] * yy[cm], axis=1) / wsum
    # alternative background: annulus 4.5-7 px, sigma-clipped mean
    ann = (r > 4.5) & (r < 7)
    vals = cube[:, ann]
    medv = np.nanmedian(vals, axis=1)
    mad = 1.4826 * np.nanmedian(np.abs(vals - medv[:, None]), axis=1)
    clip = np.where(np.abs(vals - medv[:, None]) < 3 * mad[:, None], vals, np.nan)
    out["bkg_ann"] = np.nanmean(clip, axis=1)
    ap = r <= 2.0
    out["ap2.0_annbkg"] = np.nansum(cube[:, ap] - out["bkg_ann"][:, None], axis=1)
    return out
