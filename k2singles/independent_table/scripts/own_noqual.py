"""Own TESScut aperture photometry WITHOUT the quality mask (to look at flagged cadences
near gaps). Usage: own_noqual.py EPIC SECTOR TSTART TEND [radius]"""
import sys, os, glob, re, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import lib
from targets import T
e, sec, a, b = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
rad = float(sys.argv[5]) if len(sys.argv) > 5 else 1.5
r = T[T.epic == e].iloc[0]
fn = [f for f in lib.tesscut_files(r.ra, r.dec) if f"tess-s{sec:04d}" in f][0]
h = fits.open(fn); d = h[1].data
t = d["TIME"]; F = d["FLUX"].astype(float); q = d["QUALITY"]
x0, y0 = WCS(h[2].header).all_world2pix([[r.ra, r.dec]], 0)[0]
yy, xx = np.mgrid[0:F.shape[1], 0:F.shape[2]]; rr = np.hypot(xx - x0, yy - y0)
med = np.nanmedian(F, axis=0); outer = rr > 3.5; bm = outer & (med <= np.nanpercentile(med[outer], 40))
bkg = np.nanmedian(F[:, bm], axis=1)
fl = np.nansum((F - bkg[:, None, None])[:, rr <= rad], axis=1)
m = (t > a) & (t < b) & np.isfinite(fl) & (fl != 0)
ref = np.nanmedian(fl[m])
for ti, fi, qi, bi in zip(t[m], fl[m] / ref, q[m], bkg[m]):
    print(f"{ti:.4f} {fi:.4f} q={qi} bkg={bi:.1f}")
