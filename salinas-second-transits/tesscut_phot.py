"""TESScut FFI cutouts: simple aperture photometry, difference images, centroids.

For sectors with no light curve product, and for pixel-level vetting. Times are
BTJD as in the cutout files (SPOC FFI barycentric correction per target).

  t, f, e, info = sector_lightcurve(tic, ra, dec, sector, cache_dir)
"""
from pathlib import Path

import numpy as np


def cutout(ra, dec, sector, cache_dir, size=15):
    """Download (or reuse) a size x size TESScut cutout; returns the TessTargetPixelFile."""
    import lightkurve as lk
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    hits = sorted(cache_dir.glob(f"*s{sector:04d}*{ra:.6f}*{dec:.6f}*{size}x{size}*.fits"))
    if hits:
        return lk.read(hits[0])
    sr = lk.search_tesscut(f"{ra} {dec}", sector=int(sector))
    if len(sr) == 0:
        raise RuntimeError(f"no TESScut data for sector {sector}")
    tpf = sr[0].download(cutout_size=size, download_dir=str(cache_dir), quality_bitmask="none")
    return tpf


def target_pixel(tpf, ra, dec):
    """(col, row) of the target in cutout pixel coordinates (0-based)."""
    x, y = tpf.wcs.all_world2pix([[ra, dec]], 0)[0]
    return float(x), float(y)


def apertures(tpf, ra, dec, n_sigma=3.0, max_radius=3.0):
    """Target aperture, background mask and distance (px) from the target to the
    nearest pixel off the science array.

    Aperture: pixels n_sigma above background, within max_radius px of the
    target. Background: faintest 30% of on-chip pixels more than 4 px away.
    Off-chip pixels (cutouts at a detector edge are padded with virtual pixels
    near zero with almost no noise) are never used.
    """
    img = np.nanmedian(tpf.flux.value, axis=0)
    valid = np.isfinite(img) & (np.nanstd(tpf.flux.value, axis=0) > 0.2)
    x0, y0 = target_pixel(tpf, ra, dec)
    yy, xx = np.indices(img.shape)
    r = np.hypot(xx - x0, yy - y0)
    iy, ix = int(round(y0)), int(round(x0))
    on = 0 <= iy < img.shape[0] and 0 <= ix < img.shape[1] and valid[iy, ix]
    # distance to the nearest off-chip pixel; negative (to the nearest on-chip one) if off-chip
    edge_px = (float(r[~valid].min()) if (~valid).any() else np.inf) if on else -float(r[valid].min())
    far = (r > 4) & valid
    bkg_mask = far & (img <= np.nanpercentile(img[far], 30))
    bkg = np.nanmedian(img[bkg_mask])
    sig = 1.4826 * np.nanmedian(np.abs(img[bkg_mask] - bkg))
    ap = (img - bkg > n_sigma * sig) & (r <= max_radius) & valid
    if not ap.any() or not on or not ap[iy, ix]:
        ap = (r <= 1.5) & valid  # faint, blended or at the edge: fixed 1.5 px radius
    return ap, bkg_mask, edge_px


def sector_lightcurve(tic, ra, dec, sector, cache_dir, size=15, bkg_nsig=3.0):
    """Background-subtracted aperture light curve, normalised.

    Background per cadence is the median of the background mask (SPOC
    background-subtracts some FFI releases but not others; this works for
    both). Keeps quality == 0 cadences whose background is within bkg_nsig
    robust sigma of the sector median: FFI quality flags miss most scattered
    light, which leaves ramps and fake dips at the start and end of orbits.
    info['edge_px'] < ~2 means part of the star falls off the detector.
    """
    tpf = cutout(ra, dec, sector, cache_dir, size)
    ap, bkg_mask, edge_px = apertures(tpf, ra, dec)
    if edge_px < 0:
        raise RuntimeError(f"target falls {-edge_px:.1f} px off the detector's science pixels "
                           f"(camera {tpf.camera} CCD {tpf.ccd}); no usable photometry")
    fl = tpf.flux.value
    ok = (tpf.quality == 0) & np.isfinite(tpf.time.value)
    bkg = np.nanmedian(fl[:, bkg_mask], axis=1)
    if np.nanmax(np.abs(tpf.flux_bkg.value)) > 0:  # SPOC subtracted it: cut on theirs
        cut = np.nanmedian(tpf.flux_bkg.value[:, ap], axis=1) + bkg
    else:
        cut = bkg
    f = np.nansum(fl[:, ap], axis=1) - bkg * ap.sum()
    e = np.sqrt(np.nansum(tpf.flux_err.value[:, ap] ** 2, axis=1))
    t = tpf.time.value
    m = ok & np.isfinite(f) & (f > 0) & np.isfinite(cut)
    b0 = np.nanmedian(cut[m])
    bs = 1.4826 * np.nanmedian(np.abs(cut[m] - b0))
    low_bkg = cut < b0 + bkg_nsig * bs
    n_scat = int((m & ~low_bkg).sum())
    m &= low_bkg
    med = np.nanmedian(f[m])
    info = dict(tic=tic, sector=int(sector), n_ap=int(ap.sum()), edge_px=edge_px,
                cadence_s=float(np.nanmedian(np.diff(t[m])) * 86400),
                n_high_bkg_dropped=n_scat, aperture=ap, tpf=tpf, ok=m)
    return t[m], f[m] / med, e[m] / med, info


def difference_image(tpf, tc, dur_d, ok=None):
    """Mean out-of-transit minus mean in-transit image (positive where flux drops),
    out-of-transit taken from 0.75-2.5 durations either side of tc."""
    t = tpf.time.value
    if ok is None:
        ok = (tpf.quality == 0) & np.isfinite(t)
    dt = np.abs(t - tc)
    inn = ok & (dt < 0.4 * dur_d)
    out = ok & (dt > 0.75 * dur_d) & (dt < 2.5 * dur_d)
    if inn.sum() < 3 or out.sum() < 6:
        return None, None, int(inn.sum()), int(out.sum())
    fl = tpf.flux.value
    img_out = np.nanmean(fl[out], axis=0)
    diff = img_out - np.nanmean(fl[inn], axis=0)
    return diff, img_out, int(inn.sum()), int(out.sum())


def centroid(img, mask=None):
    """Flux-weighted centroid (col, row) of the positive part of img, within mask."""
    w = np.where(np.isfinite(img), img, 0.0).clip(min=0)
    if mask is not None:
        w = w * mask
    yy, xx = np.indices(img.shape)
    s = w.sum()
    return (float((w * xx).sum() / s), float((w * yy).sum() / s)) if s > 0 else (np.nan, np.nan)
