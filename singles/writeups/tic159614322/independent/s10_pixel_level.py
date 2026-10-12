"""Pixel-level localisation of the S21 dip vs the known S14/S41 events.
1) Per-pixel fit: pixel(t) = a + b*(t-tc) [+ c*(t-tc)^2] - D_pix * S(t), using a local background-subtracted cube.
   The depth map D_pix is the 'difference image'.  Compare it with the median image (on-target expectation:
   D map proportional to the target's PSF, i.e. fractional depth roughly uniform over the core pixels).
2) Empirical calibration of the in-transit centroid-shift statistic: evaluate the same statistic (one-sided
   baseline, as forced by the orbit gap) at many random trial times in S14-S23 30-min data.
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from model import shape
from phot import find_cutout, RA, DEC
from lcio import OUT

E = json.load(open(os.path.join(OUT, "epochs.json")))
T14 = E["T14_d"]
TC21 = json.load(open(os.path.join(OUT, "highn_density_ephem.json")))["tc_S21"]["value"]


def cube_for(sector):
    h = fits.open(find_cutout(sector))
    d = h[1].data
    w = WCS(h[2].header)
    x0, y0 = w.all_world2pix([[RA, DEC]], 0)[0]
    t = np.asarray(d["TIME"], dtype=np.float64)
    cube = np.asarray(d["FLUX"], dtype=np.float64)
    q = np.asarray(d["QUALITY"])
    yy, xx = np.mgrid[0:cube.shape[1], 0:cube.shape[2]]
    r = np.hypot(xx - x0, yy - y0)
    med = np.nanmedian(cube, axis=0)
    far = r > 4.5
    bm = far & (med < np.nanpercentile(med[far], 60))
    bkg = np.nanmedian(cube[:, bm], axis=1)
    return t, cube - bkg[:, None, None], q, x0, y0, xx, yy, r


def depth_map(sector, tc, tmin, tmax, deg=1):
    t, cube, q, x0, y0, xx, yy, r = cube_for(sector)
    cad = np.median(np.diff(t)) * 86400
    m = (t > tmin) & (t < tmax) & (q == 0) & np.all(np.isfinite(cube), axis=(1, 2))
    tt, cc = t[m], cube[m]
    S = shape(tt, tc, cad)
    X = np.vstack([(tt - tc) ** i for i in range(deg + 1)] + [-S]).T
    Y = cc.reshape(len(tt), -1)
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ beta
    rms = resid.std(axis=0, ddof=X.shape[1])
    sig = np.sqrt(np.linalg.inv(X.T @ X)[-1, -1]) * rms
    D = beta[-1].reshape(cube.shape[1:])
    sD = sig.reshape(cube.shape[1:])
    base = beta[0].reshape(cube.shape[1:])
    return dict(D=D, sD=sD, base=base, x0=x0, y0=y0, xx=xx, yy=yy, r=r)


def summarize(dm, label):
    D, sD, base, r, xx, yy = dm["D"], dm["sD"], dm["base"], dm["r"], dm["xx"], dm["yy"]
    core = r <= 2.5
    # least-squares scale of D onto the baseline image within core (on-target model: D = delta * base_target)
    w = 1 / sD[core] ** 2
    delta = np.sum(w * D[core] * base[core]) / np.sum(w * base[core] ** 2)
    chi_on = np.sum(w * (D[core] - delta * base[core]) ** 2)
    # centroid of depth map vs centroid of image, both within r<=3
    c3 = r <= 3.0
    cxD = np.sum(D[c3] * xx[c3]) / np.sum(D[c3]); cyD = np.sum(D[c3] * yy[c3]) / np.sum(D[c3])
    cxI = np.sum(base[c3] * xx[c3]) / np.sum(base[c3]); cyI = np.sum(base[c3] * yy[c3]) / np.sum(base[c3])
    # bootstrap-free error via linear propagation for the depth-map centroid
    Ds = np.sum(D[c3])
    ex = np.sqrt(np.sum(((xx[c3] - cxD) * sD[c3]) ** 2)) / Ds
    ey = np.sqrt(np.sum(((yy[c3] - cyD) * sD[c3]) ** 2)) / Ds
    snr_tot = np.sum(D[core]) / np.sqrt(np.sum(sD[core] ** 2))
    out = dict(label=label, frac_depth_core=float(delta), chi2_on_target_model=float(chi_on), n_core=int(core.sum()),
               depthmap_centroid=[float(cxD), float(cyD)], image_centroid=[float(cxI), float(cyI)],
               target_pix=[float(dm["x0"]), float(dm["y0"])],
               offset_D_minus_image_pix=[float(cxD - cxI), float(cyD - cyI)], offset_err_pix=[float(ex), float(ey)],
               total_snr_core=float(snr_tot))
    print(f"{label}: frac depth (core fit) {delta*1e3:.2f} ppt, chi2(on-target PSF model, {core.sum()} px) = {chi_on:.1f}; "
          f"depth-map centroid - image centroid = ({cxD-cxI:+.3f}±{ex:.3f}, {cyD-cyI:+.3f}±{ey:.3f}) px; SNR {snr_tot:.1f}")
    return out


res = {}
maps = {}
maps["S14 known"] = depth_map(14, E["TC1"], E["TC1"] - 0.8, E["TC1"] + 0.8)
maps["S41 known"] = depth_map(41, E["TC2"], E["TC2"] - 0.8, E["TC2"] + 0.8)
maps["S21 dip"] = depth_map(21, TC21, TC21 - 1.0, 1883.90)
maps["S21 dip (quad)"] = depth_map(21, TC21, TC21 - 1.0, 1883.90, deg=2)
# same one-sided set-up applied to the known S14 event (truncate post-ingress data at the same relative time)
maps["S14 known, one-sided"] = depth_map(14, E["TC1"], E["TC1"] - 1.0, E["TC1"] + (1883.90 - TC21))
for k, v in maps.items():
    res[k] = summarize(v, k)

fig, axs = plt.subplots(2, len(maps), figsize=(4 * len(maps), 8))
for i, (k, v) in enumerate(maps.items()):
    axs[0, i].imshow(v["base"], origin="lower"); axs[0, i].plot(v["x0"], v["y0"], "r+")
    axs[0, i].set_title(f"{k}: baseline image", fontsize=8)
    snr = v["D"] / v["sD"]
    im = axs[1, i].imshow(snr, origin="lower", cmap="RdBu_r", vmin=-8, vmax=8)
    axs[1, i].plot(v["x0"], v["y0"], "k+")
    axs[1, i].set_title(f"{k}: depth map / sigma", fontsize=8)
plt.colorbar(im, ax=axs[1, -1])
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_pixel_depth_maps.png"), dpi=90)

# --- empirical calibration of the one-sided depth-map centroid offset at random times (no transit) ---
rng = np.random.default_rng(5)
trials = []
for s in [14, 21, 22, 23]:  # 30-min cutouts available locally
    t, cube, q, *_ = cube_for(s)
    ok = (q == 0)
    for _ in range(60):
        tc = rng.uniform(t.min() + 1.2, t.max() - 0.3)
        if min(abs(tc - E["TC1"]), abs(tc - TC21)) < 1.5:
            continue
        sel = ok & (t > tc - 1.0) & (t < tc + 0.14)
        if sel.sum() < 40 or np.sum(sel & (np.abs(t - tc) < 0.17)) < 12:
            continue
        try:
            dm = depth_map(s, tc, tc - 1.0, tc + 0.14)
        except Exception:
            continue
        D, base, r, xx, yy, sD = dm["D"], dm["base"], dm["r"], dm["xx"], dm["yy"], dm["sD"]
        core = r <= 2.5
        w = 1 / sD[core] ** 2
        delta = np.sum(w * D[core] * base[core]) / np.sum(w * base[core] ** 2)
        chi = np.sum(w * (D[core] - delta * base[core]) ** 2)
        trials.append((s, tc, delta, chi))
trials = np.array(trials)
res["random_trials"] = dict(n=int(len(trials)), chi2_on_model_median=float(np.median(trials[:, 3])),
                            chi2_on_model_p90=float(np.percentile(trials[:, 3], 90)),
                            chi2_on_model_p99=float(np.percentile(trials[:, 3], 99)),
                            frac_depth_abs_p95=float(np.percentile(np.abs(trials[:, 2]), 95)))
print("random one-sided trials (no transit expected):", res["random_trials"])
json.dump(res, open(os.path.join(OUT, "pixel_level.json"), "w"), indent=1)
