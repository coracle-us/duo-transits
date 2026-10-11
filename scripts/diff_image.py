"""
Vetting: locate the source of a transit signal with TESS FFI difference images.

For each requested sector: download a TESScut FFI cutout (cached in data/raw/tesscut/),
average the out-of-transit and in-transit frames around every transit of the given
ephemeris in that sector, and form diff = out - in (positive where flux was lost).
Each frame stack is locally detrended per pixel with a linear fit to out-of-transit frames
of that transit. We then:
  - fit the diff image with a PRF-like 2-D Gaussian (sigma ~ 0.75 px) on a free position and
    report its offset from the target (Gaia DR3 position via the cutout WCS),
  - compare with Gaia neighbours,
  - make a figure: mean image, diff image, diff/err, with Gaia stars overplotted.

Usage:  python scripts/diff_image.py TIC P T0 dur_h SECTOR [SECTOR ...]
Output: results/vetting/<TIC>_diffimg_s<SECTOR>.png and a summary line per sector.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "tesscut"
OUT = ROOT / "results" / "vetting"
CAT = ROOT / "data" / "catalog"
SIZE = 15


def get_cutout(tic, ra, dec, sector):
    d = RAW / str(tic)
    d.mkdir(parents=True, exist_ok=True)
    hits = list(d.glob(f"*s{sector:04d}*.fits"))
    if hits:
        return hits[0]
    import lightkurve as lk
    sr = lk.search_tesscut(f"{ra} {dec}", sector=sector)
    if len(sr) == 0:
        return None
    tpf = sr[0].download(cutout_size=SIZE, download_dir=str(d))
    return Path(tpf.path)


def gaussian(p, yy, xx):
    a, x0, y0, s, c = p
    return a * np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * s ** 2)) + c


def analyse(tic, P, T0, D, sector, ra, dec, gaia, out_dir=None):
    path = get_cutout(tic, ra, dec, sector)
    if path is None:
        return None
    with fits.open(path) as h:
        tab = h[1].data
        wcs = WCS(h[2].header)
        t = tab["TIME"]
        flux = tab["FLUX"]
        q = tab["QUALITY"]
    good = np.isfinite(t) & (q == 0) & np.all(np.isfinite(flux.reshape(len(t), -1)), axis=1)
    t, flux = t[good], flux[good]
    ks = np.arange(np.ceil((t.min() - T0) / P), np.floor((t.max() - T0) / P) + 1)
    diffs, errs, used = [], [], []
    for k in ks:
        Tk = T0 + k * P
        dt = t - Tk
        inn = np.abs(dt) < 0.35 * D
        out = (np.abs(dt) > D / 2 + 0.05) & (np.abs(dt) < D / 2 + 0.6)
        if inn.sum() < 2 or out.sum() < 6 or not ((dt[out] < 0).any() and (dt[out] > 0).any()):
            continue
        # per-pixel linear trend from out-of-transit frames
        A = np.column_stack([np.ones(out.sum()), dt[out]])
        Y = flux[out].reshape(out.sum(), -1)
        coef, *_ = np.linalg.lstsq(A, Y, rcond=None)
        pred_in = np.column_stack([np.ones(inn.sum()), dt[inn]]) @ coef
        resid_out = Y - A @ coef
        d = (pred_in - flux[inn].reshape(inn.sum(), -1)).mean(axis=0)
        e = resid_out.std(axis=0) * np.sqrt(1 / inn.sum() + 1 / out.sum())
        diffs.append(d)
        errs.append(e)
        used.append(round(Tk, 3))
    if not diffs:
        return dict(sector=sector, n_transits=0)
    w = 1 / np.array(errs) ** 2
    diff = (np.sum(np.array(diffs) * w, axis=0) / np.sum(w, axis=0)).reshape(SIZE, SIZE)
    derr = (1 / np.sqrt(np.sum(w, axis=0))).reshape(SIZE, SIZE)
    mean_img = np.nanmedian(flux, axis=0)
    # target and neighbour pixel positions
    tx, ty = wcs.all_world2pix(ra, dec, 0)
    gx, gy = wcs.all_world2pix(gaia.ra.values, gaia.dec.values, 0)
    yy, xx = np.mgrid[0:SIZE, 0:SIZE]
    iy, ix = np.unravel_index(np.argmax(diff), diff.shape)
    r = least_squares(lambda p: ((gaussian(p, yy, xx) - diff) / derr).ravel(),
                      [diff.max(), ix, iy, 0.8, 0.0],
                      bounds=([0, 0, 0, 0.4, -np.inf], [np.inf, SIZE - 1, SIZE - 1, 2.0, np.inf]))
    try:
        cov = np.linalg.inv(r.jac.T @ r.jac) * max(1, np.sum(r.fun ** 2) / (r.fun.size - 5))
        ex, ey = np.sqrt(cov[1, 1]), np.sqrt(cov[2, 2])
    except np.linalg.LinAlgError:
        ex = ey = np.nan
    fx, fy = r.x[1], r.x[2]
    # Out-of-transit centroid of the target itself (same Gaussian model, 5x5 box around the
    # WCS position): offsets are measured from THIS, so WCS/pixel-convention errors cancel
    # (Bryson et al. 2013 style).
    bx0, by0 = int(round(float(tx))), int(round(float(ty)))
    sl = (slice(max(by0 - 2, 0), by0 + 3), slice(max(bx0 - 2, 0), bx0 + 3))
    m_img = mean_img[sl]
    myy, mxx = yy[sl], xx[sl]
    rm = least_squares(lambda p: (gaussian(p, myy, mxx) - m_img).ravel(),
                       [m_img.max() - np.median(mean_img), bx0, by0, 0.8, np.median(mean_img)],
                       bounds=([0, bx0 - 2, by0 - 2, 0.4, -np.inf], [np.inf, bx0 + 2, by0 + 2, 2.0, np.inf]))
    cx, cy = rm.x[1], rm.x[2]
    wcs_off = float(np.hypot(cx - tx, cy - ty))
    # Offsets are measured from the target's WCS (Gaia) position. Revision 2026-10-11: the
    # earlier version measured from the out-of-transit image centroid, which LOCKS ONTO A
    # BRIGHTER NEIGHBOUR when one sits within ~2-3 px (TIC 265466589: false 'off-target').
    # The WCS agrees with the image centroid to <= 0.17 px for isolated targets (validated on
    # TOI-5542, TOI-2589, TIC 1504460, 14445414), so it is the safer reference. The centroid
    # offset is kept as a diagnostic; > 0.5 px flags a crowded field.
    off = np.hypot(fx - tx, fy - ty)
    off_err = np.hypot(ex, ey)
    off_vs_centroid = float(np.hypot(fx - cx, fy - cy))
    # nearest Gaia star to the fitted source position
    dist = np.hypot(gx - fx, gy - fy)
    j = np.argmin(dist)
    res = dict(sector=sector, n_transits=len(used), transits=used, wcs_vs_image_centroid_px=wcs_off,
               crowded_reference=wcs_off > 0.5, offset_from_image_centroid_px=off_vs_centroid,
               snr_peak=float(diff.max() / derr.ravel()[np.argmax(diff)]),
               offset_from_target_px=float(off), offset_err_px=float(off_err),
               offset_from_target_arcsec=float(off * 21), offset_sigma=float(off / off_err) if off_err else None,
               nearest_gaia_to_source=dict(source_id=int(gaia.source_id.iloc[j]),
                                           G=float(gaia.phot_g_mean_mag.iloc[j]),
                                           sep_from_target_arcsec=float(gaia.sep.iloc[j]),
                                           dist_from_fit_px=float(dist[j])))
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.3))
    for ax, img, lab in zip(axs, [mean_img, diff, diff / derr],
                            ["median image", "difference (out - in)", "difference / error"]):
        im = ax.imshow(img, origin="lower", cmap="viridis")
        plt.colorbar(im, ax=ax, fraction=0.046)
        bright = gaia.phot_g_mean_mag < gaia.phot_g_mean_mag.iloc[0] + 5
        ax.scatter(gx[bright], gy[bright], s=40 * 10 ** (-0.2 * (gaia.phot_g_mean_mag[bright]
                                                                  - gaia.phot_g_mean_mag.iloc[0])),
                   facecolors="none", edgecolors="w", lw=0.8)
        ax.plot(tx, ty, "r+", ms=14, mew=2)
        ax.plot(fx, fy, "mx", ms=10, mew=2)
        ax.set_title(lab, fontsize=9)
        ax.set_xlim(-0.5, SIZE - 0.5)
        ax.set_ylim(-0.5, SIZE - 0.5)
    fig.suptitle(f"TIC {tic} S{sector}: {len(used)} transits; red + target, magenta x fitted "
                 f"source (offset {off*21:.1f}\" +/- {off_err*21:.1f}\"); circles Gaia (G<target+5)",
                 fontsize=9)
    fig.tight_layout()
    out_dir = Path(out_dir) if out_dir else OUT
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{tic}_diffimg_s{sector}.png", dpi=85)
    plt.close(fig)
    return res


def gaia_field(ra, dec):
    """Gaia DR3 sources (G < 19) covering the cutout, nearest first (row 0 = target)."""
    from astroquery.gaia import Gaia
    q = f"""SELECT source_id, ra, dec, phot_g_mean_mag,
            DISTANCE(POINT({ra},{dec}), POINT(ra,dec))*3600 AS sep
            FROM gaiadr3.gaia_source
            WHERE 1=CONTAINS(POINT(ra,dec), CIRCLE({ra},{dec},{SIZE*21/2*1.5/3600}))
            AND phot_g_mean_mag < 19 ORDER BY sep"""
    return Gaia.launch_job(q).get_results().to_pandas()


def main():
    tic, P, T0, dur_h = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    sectors = [int(s) for s in sys.argv[5:]]
    c = pd.read_csv(CAT / "candidates.csv").set_index("TIC").loc[tic]
    gaia = gaia_field(c.RA, c.Dec)
    tgt = gaia.iloc[0]
    for s in sectors:
        r = analyse(tic, P, T0, dur_h / 24, s, tgt.ra, tgt.dec, gaia)
        print(r)


if __name__ == "__main__":
    main()
