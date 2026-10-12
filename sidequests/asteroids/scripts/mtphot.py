"""
Moving-target photometry of an asteroid from a TESScut moving-target cutout.

Steps (see README for the reasoning):
 1. Cutout: astroquery.mast.Tesscut.download_cutouts(objectname=<Horizons id>, moving_target=True,
    sector=S, size=SIZE). Cached in cache/tesscut/ under a safe name. The table gives, per FFI,
    a SIZExSIZE pixel stamp centred on the Horizons position (TGT_X, TGT_Y within the stamp,
    TGT_RA/TGT_DEC), and HDU 2 carries the full-CCD SIP WCS, from which the CCD position of
    every stamp corner is recovered (corner = round(ccd_xy(TGT_RA, TGT_DEC) - TGT_XY); checked to
    be integer to < 0.25 px).
 2. Per-frame background = median of the faintest 40% of stamp pixels (the convention of
    singles/scripts/tesscut_photometry.photometry; that function writes into data/raw/ so it is
    not called here, only its rule is reused).
 3. Static scene: every stamp pixel is mapped to its CCD pixel; for each CCD pixel the median of
    its background-subtracted values over all good frames in which the asteroid was > R_EX px away
    is the static-sky model. Subtracting it removes the stars the asteroid passes.
 4. Aperture photometry on the difference stamp at (TGT_X, TGT_Y): pixel centres within R_AP;
    a residual sky from the median of the annulus R_IN..R_OUT is removed. (Options for a soft
    aperture and re-centring exist but are off; see run() and README.)
 5. Flags per frame: QUALITY != 0; stamp not entirely on science pixels (flag 16); fewer than NMIN static samples for an aperture pixel;
    star contamination (static-scene flux in the aperture, 'star_ap', see period.py for the cut);
    noisy frames (annulus robust rms > 2.5 x the sector median: scattered light, glints).
 6. Ephemeris (JPL Horizons, observer = TESS, code 500@-95): r, Delta, alpha, V, interpolated to
    the frame times; light-time corrected times t_lt = t - Delta/c.
Output: data/lightcurves/<safe id>_s<SECTOR>.csv.gz

Reused from the repo (imported): scripts/diff_image.py `gaussian` (2-D Gaussian PRF model) for
the stacked co-moving difference-image centroid check.

Usage: python mtphot.py OBJ_ID SECTOR [SIZE]     (OBJ_ID in Horizons form, e.g. '624;')
"""
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import diff_image as di  # noqa: E402

CACHE = HERE / "cache"
LCDIR = HERE / "data" / "lightcurves"
SIZE = 15
R_EX = 3.5     # px: exclusion radius around the asteroid for the static-scene model
R_AP = 1.6     # px: aperture radius (pixel centres)
R_IN, R_OUT = 3.5, 7.0
NMIN = 8       # min. static-scene samples per aperture pixel
C_AU_PER_DAY = 173.1446327


def safe(obj):
    return obj.replace(";", "").replace(" ", "_").replace("/", "_")


def fetch_cutout(obj, sector, size=SIZE):
    d = CACHE / "tesscut"
    d.mkdir(parents=True, exist_ok=True)
    out = d / f"{safe(obj)}_s{sector:03d}_{size}.fits"
    if out.exists():
        return out
    from astroquery.mast import Tesscut
    tmp = d / f"tmp_{safe(obj)}_{sector}"
    tmp.mkdir(exist_ok=True)
    import time
    for attempt in range(7):  # TESScut answers 429/502 under load: back off and retry
        try:
            m = Tesscut.download_cutouts(objectname=obj, moving_target=True, sector=sector, size=size,
                                         path=str(tmp))
            break
        except Exception as e:  # noqa: BLE001
            if attempt == 6 or not any(c in str(e) for c in ("429", "502", "503", "504", "imed out")):
                raise
            wait = min(60 * 2 ** attempt, 600)
            print(f"  {obj} S{sector}: {str(e)[:40]} -> retry in {wait}s", flush=True)
            time.sleep(wait)
    paths = [Path(p) for p in m["Local Path"]]
    if not paths:
        raise RuntimeError(f"no cutout for {obj} S{sector}")
    if len(paths) > 1:  # object crossed more than one camera/CCD: keep the longest
        paths.sort(key=lambda p: p.stat().st_size, reverse=True)
        print(f"  {len(paths)} cutouts (several CCDs?), keeping the largest", flush=True)
    shutil.move(str(paths[0]), out)
    shutil.rmtree(tmp, ignore_errors=True)
    return out


def fetch_ephem(obj, t_btjd, step_h=3):
    """Horizons observer table from TESS (500@-95): r, delta, alpha, V at the frame times."""
    from astroquery.jplhorizons import Horizons
    f = CACHE / "ephem" / f"{safe(obj)}_{t_btjd.min():.2f}_{t_btjd.max():.2f}.csv"
    f.parent.mkdir(parents=True, exist_ok=True)
    if not f.exists():
        jd0, jd1 = t_btjd.min() + 2457000 - 0.2, t_btjd.max() + 2457000 + 0.2
        idt = obj.rstrip(";")
        h = Horizons(id=idt, id_type="smallbody", location="500@-95",
                     epochs={"start": f"JD{jd0:.4f}", "stop": f"JD{jd1:.4f}", "step": f"{step_h}h"})
        e = h.ephemerides(quantities="1,9,19,20,24", extra_precision=True)
        cols = {"datetime_jd": "jd", "RA": "ra", "DEC": "dec", "V": "V", "r": "r", "delta": "delta",
                "alpha": "alpha"}
        df = e.to_pandas()
        if "V" not in df.columns and "Tmag" in df.columns:  # comets
            df["V"] = df["Tmag"]
        df = df[[c for c in cols if c in df.columns]].rename(columns=cols)
        df.to_csv(f, index=False)
    df = pd.read_csv(f)
    # Horizons epochs are UT; BTJD is TDB (offset ~69 s, irrelevant at a 3-h interpolation grid)
    x = df.jd.values - 2457000
    return {c: np.interp(t_btjd, x, df[c].values) for c in ["r", "delta", "alpha", "V", "ra", "dec"]}


def load_cutout(path):
    with fits.open(path) as h:
        d = h[1].data
        wcs = WCS(h[2].header)
        hdr = h[0].header
        tab = dict(t=np.array(d["TIME"], float), flux=np.array(d["FLUX"], float),
                   q=np.array(d["QUALITY"]), tx=np.array(d["TGT_X"], float), ty=np.array(d["TGT_Y"], float),
                   ra=np.array(d["TGT_RA"], float), dec=np.array(d["TGT_DEC"], float))
        tab.update(camera=hdr.get("CAMERA"), ccd=hdr.get("CCD"), sector=hdr.get("SECTOR"))
    X, Y = wcs.all_world2pix(tab["ra"], tab["dec"], 0)
    cx, cy = X - tab["tx"], Y - tab["ty"]
    tab["corner_frac_max"] = float(max(np.nanmax(np.abs(cx - np.round(cx))), np.nanmax(np.abs(cy - np.round(cy)))))
    tab["cx"], tab["cy"] = np.round(cx).astype(int), np.round(cy).astype(int)
    tab["ccd_x"], tab["ccd_y"] = X, Y
    return tab


def photometry(tab, dx=0.0, dy=0.0, soft=False):
    """dx, dy: aperture offset from the Horizons position (px; from the stacked centroid).
    soft=True: aperture weight clip(R_AP + 0.5 - d, 0, 1), a linear 1-px taper approximating
    exact pixel overlap; soft=False: hard pixel-centre aperture (the first version)."""
    t, flux, q = tab["t"], tab["flux"], tab["q"]
    n, S, _ = flux.shape
    yy, xx = np.mgrid[0:S, 0:S]
    good = np.isfinite(t) & (q == 0) & np.all(np.isfinite(flux.reshape(n, -1)), axis=1)
    # the whole stamp must lie on science pixels (cols 44-2091, rows 0-2047, 0-based); TESScut
    # keeps following the target into the collateral columns (found on (12198) in validation)
    onchip = (tab["cx"] >= 44) & (tab["cx"] + S - 1 <= 2091) & (tab["cy"] >= 0) & (tab["cy"] + S - 1 <= 2047)
    good &= onchip
    # per-frame background (faintest-40% median), stamp-level
    fl = flux.reshape(n, -1)
    p40 = np.nanpercentile(fl, 40, axis=1)
    bkg = np.array([np.nanmedian(r[r <= p]) if np.isfinite(p) else np.nan for r, p in zip(fl, p40)])
    sub = flux - bkg[:, None, None]
    # CCD pixel ids
    W = 4096
    X = tab["cx"][:, None, None] + xx[None]
    Y = tab["cy"][:, None, None] + yy[None]
    pid = (Y * W + X)
    dist = np.hypot(xx[None] - tab["tx"][:, None, None], yy[None] - tab["ty"][:, None, None])
    use = good[:, None, None] & (dist > R_EX) & np.isfinite(sub)
    ids, vals = pid[use], sub[use]
    o = np.argsort(ids, kind="stable")
    ids, vals = ids[o], vals[o]
    uid, start, counts = np.unique(ids, return_index=True, return_counts=True)
    med = np.array([np.median(vals[s:s + c]) for s, c in zip(start, counts)])
    # look-up of the static model for every stamp pixel
    k = np.searchsorted(uid, pid.ravel())
    k = np.clip(k, 0, len(uid) - 1)
    found = uid[k] == pid.ravel()
    static = np.where(found, med[k], np.nan).reshape(pid.shape)
    nsamp = np.where(found, counts[k], 0).reshape(pid.shape)
    diff = sub - static
    dap = np.hypot(xx[None] - (tab["tx"] + dx)[:, None, None], yy[None] - (tab["ty"] + dy)[:, None, None])
    wap = np.clip(R_AP + 0.5 - dap, 0, 1) if soft else (dap <= R_AP).astype(float)
    ap = wap > 0
    ann = (dist > R_IN) & (dist <= R_OUT)
    F = np.full(n, np.nan)
    sky = np.full(n, np.nan)
    rms = np.full(n, np.nan)
    star_ap = np.full(n, np.nan)
    nmin_ap = np.zeros(n, int)
    npix = np.zeros(n, int)
    werr = np.ones(n)
    for i in range(n):
        a, b = ap[i], ann[i] & np.isfinite(diff[i])
        if b.sum() < 10:
            continue
        av = diff[i][b]
        sky[i] = np.median(av)
        rms[i] = 1.4826 * np.median(np.abs(av - sky[i]))
        w = wap[i][a]
        F[i] = np.nansum(w * (diff[i][a] - sky[i]))
        npix[i] = a.sum()
        werr[i] = np.sqrt(np.sum(w ** 2))
        star_ap[i] = np.nansum(w * static[i][a])
        nmin_ap[i] = nsamp[i][a].min()
    err = rms * werr
    noisy = rms > 2.5 * np.nanmedian(rms[good])
    flag = (~good).astype(int) * 1 + (nmin_ap < NMIN) * 2 + noisy * 4 + ~np.isfinite(F) * 8 + (~onchip) * 16
    return dict(F=F, err=err, sky=sky, rms=rms, bkg=bkg, star_ap=star_ap, nmin_ap=nmin_ap, flag=flag,
                diff=diff, good=good)


def centroid_check(tab, ph, nbin=200):
    """Stack co-moving difference stamps (shifted so the target is at the stamp centre to the
    nearest pixel) and fit the 2-D Gaussian of diff_image.py: offset of the stacked source from
    the Horizons position, in px."""
    ok = ph["flag"] == 0
    if ok.sum() < 20:
        return dict(cen_dx=np.nan, cen_dy=np.nan, cen_amp=np.nan)
    S = ph["diff"].shape[1]
    c = S // 2
    stack = []
    for i in np.where(ok)[0]:
        sx, sy = int(round(c - tab["tx"][i])), int(round(c - tab["ty"][i]))
        stack.append(np.roll(np.roll(ph["diff"][i] - ph["sky"][i], sy, 0), sx, 1))
        stack[-1] = (stack[-1], tab["tx"][i] + sx - c, tab["ty"][i] + sy - c)
    img = np.nanmedian(np.array([s[0] for s in stack]), axis=0)
    fx = np.mean([s[1] for s in stack])
    fy = np.mean([s[2] for s in stack])
    yy, xx = np.mgrid[0:S, 0:S]
    sl = (slice(c - 3, c + 4), slice(c - 3, c + 4))
    r = least_squares(lambda p: (di.gaussian(p, yy[sl], xx[sl]) - img[sl]).ravel(),
                      [max(float(img[c, c]), 1e-3), c, c, 0.8, 0.0],
                      bounds=([0, c - 3, c - 3, 0.3, -np.inf], [np.inf, c + 3, c + 3, 3.0, np.inf]))
    return dict(cen_dx=float(r.x[1] - (c + fx)), cen_dy=float(r.x[2] - (c + fy)), cen_amp=float(r.x[0]),
                cen_sigma=float(r.x[3]), stack_img=img)


def run(obj, sector, size=SIZE):
    path = fetch_cutout(obj, sector, size)
    tab = load_cutout(path)
    # Adopted: hard pixel-centre aperture at the Horizons position (version 1). A soft-edged
    # aperture and re-centring on the stacked source were tried during validation
    # (compare_aperture.py -> data/aperture_comparison.csv): they helped some objects and broke the
    # period of (698), so they are NOT used; the stacked offset is kept as a diagnostic only.
    ph = photometry(tab, soft=False)
    cen0 = centroid_check(tab, ph)
    dx = dy = 0.0
    cen = cen0
    cen["ap_dx"], cen["ap_dy"] = dx, dy
    eph = fetch_ephem(obj, tab["t"][np.isfinite(tab["t"])])
    t = tab["t"]
    tt = np.where(np.isfinite(t), t, np.nanmedian(t))
    eph = {k: np.interp(tt, tab["t"][np.isfinite(t)], v) for k, v in eph.items()}
    df = pd.DataFrame(dict(time=t, t_lt=t - eph["delta"] / C_AU_PER_DAY, flux=ph["F"], flux_err=ph["err"],
                           flag=ph["flag"], star_ap=ph["star_ap"], nmin_ap=ph["nmin_ap"], sky_res=ph["sky"],
                           rms=ph["rms"], bkg=ph["bkg"], quality=tab["q"], tgt_x=tab["tx"], tgt_y=tab["ty"],
                           ccd_x=tab["ccd_x"], ccd_y=tab["ccd_y"],
                           r=eph["r"], delta=eph["delta"], alpha=eph["alpha"], V=eph["V"]))
    LCDIR.mkdir(parents=True, exist_ok=True)
    out = LCDIR / f"{safe(obj)}_s{sector:03d}.csv.gz"
    df.to_csv(out, index=False, compression="gzip", float_format="%.6f")
    meta = dict(obj=obj, sector=sector, camera=tab["camera"], ccd=tab["ccd"], n_frames=len(t),
                n_good=int((ph["flag"] == 0).sum()), corner_frac_max=round(tab["corner_frac_max"], 3),
                cen_dx=round(cen0["cen_dx"], 3), cen_dy=round(cen0["cen_dy"], 3),
                ap_dx=round(dx, 3), ap_dy=round(dy, 3),
                cen_sigma=round(cen.get("cen_sigma", np.nan), 3),
                cen_amp=round(cen["cen_amp"], 2), median_flux=float(np.nanmedian(ph["F"][ph["flag"] == 0])),
                V_med=round(float(np.nanmedian(eph["V"])), 2), track_px=round(float(np.hypot(
                    np.ptp(tab["ccd_x"]), np.ptp(tab["ccd_y"]))), 1))
    if "stack_img" in cen:
        np.save(LCDIR / f"{safe(obj)}_s{sector:03d}_stack.npy", cen["stack_img"])
    return out, meta


if __name__ == "__main__":
    obj, sector = sys.argv[1], int(sys.argv[2])
    size = int(sys.argv[3]) if len(sys.argv) > 3 else SIZE
    out, meta = run(obj, sector, size)
    print(out)
    print(meta)
