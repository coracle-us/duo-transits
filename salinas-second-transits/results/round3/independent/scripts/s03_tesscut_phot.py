"""Simple aperture photometry on TESScut cutouts, with barycentric re-correction to the
target direction.

TESScut TIME is barycentred for the FFI reference direction (TIMECORR column holds the
correction that was applied). We compute the correction for the target direction with
astropy (geocentre; TESS-geocentre offset contributes <~1 s to the *difference*) and apply
  t_new = t + (bc_target - TIMECORR).
The aperture is a circle centred on the target WCS position; radius is chosen to minimise
the 1-h binned scatter of the detrended light curve. Background = per-cadence median of
the faintest 40% of pixels outside r>3.5 px.
"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.coordinates import SkyCoord, get_body_barycentric
import astropy.units as u
from common import DATA, TARGETS, detrend

C = 299792.458  # km/s


def bary_corr(t_btjd, ra, dec):
    """Barycentric light-travel-time correction (days) for geocentre toward (ra, dec)."""
    tt = Time(t_btjd + 2457000, format="jd", scale="tdb")
    pos = get_body_barycentric("earth", tt)
    xyz = np.vstack([pos.x.to(u.km).value, pos.y.to(u.km).value, pos.z.to(u.km).value])
    n = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs").cartesian.xyz.value
    return (n @ xyz) / C / 86400.0


def phot_one(fn, ra, dec):
    with fits.open(fn) as h:
        hd0 = h[0].header
        sec = int(hd0["SECTOR"])
        d = h[1].data
        t = np.array(d["TIME"], float)
        tc = np.array(d["TIMECORR"], float)
        cube = np.array(d["FLUX"], float)
        q = np.array(d["QUALITY"])
        w = WCS(h[2].header)
    x0, y0 = w.all_world2pix(ra, dec, 0)
    ny, nx = cube.shape[1:]
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(xx - x0, yy - y0)
    good = np.isfinite(t) & (q == 0) & np.all(np.isfinite(cube.reshape(len(t), -1)), axis=1)
    if good.sum() < 100:
        return sec, None, None, None, None, None, 0.0
    t, tc, cube = t[good], tc[good], cube[good]
    med = np.nanmedian(cube, axis=0)
    outer = r > 3.5
    thr = np.nanpercentile(med[outer], 40)
    bmask = outer & (med <= thr)
    bkg = np.nanmedian(cube[:, bmask], axis=1)
    # reject cutouts where the target is not on science pixels (e.g. S11 for TIC 241550976:
    # the cutout is all ~0 e/s)
    ix, iy = int(round(float(x0))), int(round(float(y0)))
    peak = med[max(0, iy - 1):iy + 2, max(0, ix - 1):ix + 2].max() - np.nanmedian(med[bmask])
    if not (0 <= ix < nx and 0 <= iy < ny) or peak < 50:
        return sec, None, None, None, None, None, peak
    sub = cube - bkg[:, None, None]
    best = None
    for rad in [1.0, 1.5, 2.0, 2.5, 3.0]:
        ap = r <= rad
        f = sub[:, ap].sum(axis=1)
        ok = np.isfinite(f) & (f > 0)
        if ok.sum() < 100:
            continue
        fn_ = f[ok] / np.nanmedian(f[ok])
        flat, _ = detrend(t[ok], fn_, 1.0)
        # 1-h binned scatter
        nb = max(1, int(round((1 / 24) / np.nanmedian(np.diff(t)))))
        k = (len(flat) // nb) * nb
        s = np.nanstd(flat[:k].reshape(-1, nb).mean(axis=1))
        if best is None or s < best[0]:
            best = (s, rad, ok, fn_, flat)
    s, rad, ok, f, flat = best
    t, tc = t[ok], tc[ok]
    bc = bary_corr(t, ra, dec)
    dt = bc - tc
    tnew = t + dt
    pt = np.abs(np.diff(flat))
    e = np.full(len(f), 1.4826 * np.nanmedian(pt) / np.sqrt(2))
    return sec, tnew, f, e, rad, dt, s


if __name__ == "__main__":
    st = TARGETS
    os.makedirs(os.path.join(DATA, "tesscut_phot"), exist_ok=True)
    for tic, p in st.items():
        files = glob.glob(os.path.join(DATA, "tesscut", "**", f"*.fits"), recursive=True)
        for fn in sorted(files):
            hd = fits.getheader(fn, 0)
            if abs(hd.get("RA_OBJ", -99) - p["ra"]) > 0.01 or abs(hd.get("DEC_OBJ", -99) - p["dec"]) > 0.01:
                continue
            sec, t, f, e, rad, dt, s = phot_one(fn, p["ra"], p["dec"])
            if t is None:
                print(f"TIC {tic} S{sec}: target not on science pixels (peak {s:.1f} e/s) - skipped")
                stale = os.path.join(DATA, "tesscut_phot", f"tic{tic}_s{sec:04d}.npz")
                if os.path.exists(stale):
                    os.remove(stale)
                continue
            np.savez(os.path.join(DATA, "tesscut_phot", f"tic{tic}_s{sec:04d}.npz"),
                     sector=sec, t=t, f=f, e=e, rad=rad)
            print(f"TIC {tic} S{sec}: n={len(t)} r_ap={rad} 1h-scatter={s*1e6:.0f} ppm "
                  f"bary re-correction {dt.min()*1440:+.2f}..{dt.max()*1440:+.2f} min")
