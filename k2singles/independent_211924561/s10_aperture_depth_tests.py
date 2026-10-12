"""Depth vs aperture (core vs halo) for the K2 C5 event and both TESS events, plus TESS difference images.

K2: sum pixels in an aperture, then decorrelate the summed flux against a quadratic in (POS_CORR1,2)
plus a linear time term fitted out of event (|t-tc|>0.12 d) within +-1.5 d.
TESS 2-min TPFs: sum pixels, quadratic time baseline out of event (|t-tc|>0.1 d) within +-0.8 d.
Depth = mean deficit in the event core (|t-tc|<0.035 d) relative to the aperture flux.
If the signal is on the target, depth(core) ~ depth(halo); a contaminant offset by a few pixels
makes the halo or one side much deeper.
"""
import os, json
import numpy as np
from astropy.io import fits
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lcio

out = {}


def depth_k2(fl, t, p1, p2, mask, tc, oot, w, inn):
    y = np.nansum(fl[:, mask], axis=1)
    X = lambda s: np.vstack([np.ones(s.sum()), p1[s], p2[s], p1[s] ** 2, p2[s] ** 2, p1[s] * p2[s], t[s] - tc]).T
    c, *_ = np.linalg.lstsq(X(oot), y[oot], rcond=None)
    r = y[w] / (X(w) @ c) - 1
    rin = y[inn] / (X(inn) @ c) - 1
    roo = y[oot] / (X(oot) @ c) - 1
    return -np.mean(rin) * 1e3, np.std(roo) / np.sqrt(inn.sum()) * 1e3, r


# ---------- K2
TC = 179.90
d = fits.getdata(lcio._find("ktwo211924561-c05_lpd-targ.fits.gz"), 1)
t = d["TIME"] - 2167.0; fl = d["FLUX"].astype(float); q = d["QUALITY"]; p1, p2 = d["POS_CORR1"], d["POS_CORR2"]
ok = np.isfinite(t) & np.isfinite(p1) & np.isfinite(p2) & ((q & (1 << 20)) == 0)
w = ok & (np.abs(t - TC) < 1.5); inn = w & (np.abs(t - TC) < 0.035); oot = w & (np.abs(t - TC) > 0.12)
img = np.nanmedian(fl[w], axis=0)
iy, ix = np.unravel_index(np.nanargmax(img), img.shape)
yy, xx = np.indices(img.shape)
rr = np.hypot(yy - iy, xx - ix)
fin = np.isfinite(img)
masks = {"core r<=1": fin & (rr <= 1.0), "ring 1<r<=2.3": fin & (rr > 1.0) & (rr <= 2.3), "outer r>2.3": fin & (rr > 2.3),
         "left half": fin & (xx < ix) & (rr <= 2.3), "right half": fin & (xx > ix) & (rr <= 2.3),
         "lower half": fin & (yy < iy) & (rr <= 2.3), "upper half": fin & (yy > iy) & (rr <= 2.3), "all r<=2.3": fin & (rr <= 2.3)}
out["K2_C5"] = {}
for k, m in masks.items():
    dep, err, _ = depth_k2(fl, t, p1, p2, m, TC, oot, w, inn)
    out["K2_C5"][k] = dict(depth_ppt=round(dep, 3), err_ppt=round(err, 3), npix=int(m.sum()), flux_frac=round(float(np.nansum(img[m]) / np.nansum(img[fin])), 3))
    print("K2", k, out["K2_C5"][k])


# ---------- TESS
def tess(sector, tc):
    f = lcio._find(f"tess*-s{sector:04d}-0000000238606535-*-s_tp.fits")
    h = fits.open(f)
    d = h[1].data
    t = d["TIME"]; fl = d["FLUX"].astype(float); q = d["QUALITY"]
    ok = np.isfinite(t) & (q == 0) & np.all(np.isfinite(fl.reshape(len(t), -1)), axis=1)
    w = ok & (np.abs(t - tc) < 0.8); inn = w & (np.abs(t - tc) < 0.035); oot = w & (np.abs(t - tc) > 0.1)
    img = np.nanmedian(fl[w], axis=0)
    # per-pixel quadratic baseline, deficit image
    ny, nx = img.shape
    dimg = np.zeros_like(img); eimg = np.zeros_like(img)
    for a in range(ny):
        for b in range(nx):
            y = fl[:, a, b]
            c = np.polyfit(t[oot] - tc, y[oot], 2)
            r = y - np.polyval(c, t - tc)
            dimg[a, b] = -np.mean(r[inn]); eimg[a, b] = np.std(r[oot]) / np.sqrt(inn.sum())
    apmask = h[2].data
    ap = (apmask & 2) > 0
    iy, ix = np.unravel_index(np.nanargmax(img), img.shape)
    yy, xx = np.indices(img.shape)
    rr = np.hypot(yy - iy, xx - ix)
    res = {}
    for k, m in {"SPOC aperture": ap, "core r<=1": rr <= 1, "ring 1<r<=2.3": (rr > 1) & (rr <= 2.3), "outer 2.3<r<=4": (rr > 2.3) & (rr <= 4)}.items():
        y = np.nansum(fl[:, m], axis=1)
        c = np.polyfit(t[oot] - tc, y[oot], 2)
        r = y / np.polyval(c, t - tc) - 1
        res[k] = dict(depth_ppt=round(-np.mean(r[inn]) * 1e3, 3), err_ppt=round(np.std(r[oot]) / np.sqrt(inn.sum()) * 1e3, 3), npix=int(m.sum()))
        print(f"S{sector}", k, res[k])
    # centroid of the deficit image vs centroid of the direct image (flux-weighted, within r<=2.3)
    m = rr <= 2.3
    cx_img = np.sum(xx[m] * img[m]) / np.sum(img[m]); cy_img = np.sum(yy[m] * img[m]) / np.sum(img[m])
    dpos = np.clip(dimg, 0, None)
    cx_d = np.sum(xx[m] * dpos[m]) / np.sum(dpos[m]); cy_d = np.sum(yy[m] * dpos[m]) / np.sum(dpos[m])
    # bootstrap uncertainty on deficit centroid using per-pixel errors
    rng = np.random.default_rng(1)
    cs = []
    for _ in range(500):
        dd = np.clip(dimg + rng.standard_normal(dimg.shape) * eimg, 0, None)
        cs.append((np.sum(xx[m] * dd[m]) / np.sum(dd[m]), np.sum(yy[m] * dd[m]) / np.sum(dd[m])))
    cs = np.array(cs)
    res["centroid_offset_pix"] = dict(dx=round(cx_d - cx_img, 3), dy=round(cy_d - cy_img, 3), sx=round(float(np.std(cs[:, 0])), 3), sy=round(float(np.std(cs[:, 1])), 3))
    print(f"S{sector} deficit-centroid minus image-centroid (pix):", res["centroid_offset_pix"])
    snr = dimg / eimg
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    im = ax[0].imshow(img, origin="lower"); ax[0].contour(ap, levels=[0.5], colors="w"); ax[0].set_title(f"S{sector} median image"); plt.colorbar(im, ax=ax[0])
    vm = np.nanmax(np.abs(dimg)); im = ax[1].imshow(dimg, origin="lower", cmap="RdBu_r", vmin=-vm, vmax=vm); ax[1].contour(ap, levels=[0.5], colors="k")
    ax[1].set_title("deficit image (oot - in)"); plt.colorbar(im, ax=ax[1])
    im = ax[2].imshow(snr, origin="lower", cmap="RdBu_r", vmin=-5, vmax=5); ax[2].contour(ap, levels=[0.5], colors="k"); ax[2].set_title("deficit S/N"); plt.colorbar(im, ax=ax[2])
    plt.tight_layout(); plt.savefig(os.path.join(lcio.HERE, f"fig_tess_diffimg_S{sector}.png"), dpi=100); plt.close()
    return res


out["TESS_S46"] = tess(46, 2568.30)
out["TESS_S72"] = tess(72, 3270.75)
json.dump(out, open(os.path.join(lcio.HERE, "aperture_tests.json"), "w"), indent=1)
