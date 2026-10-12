"""Pixel-level test of the K2 C5 event.

For every pixel: model out-of-event flux within +-1.5 d as a quadratic in (POS_CORR1, POS_CORR2) plus a
linear time term (fit on cadences with |t - tc| > 0.12 d and no thruster firing), predict the in-event
cadences, and form a deficit image (mean predicted - observed over the 4 deepest cadences).
A real transit on the target: deficit image ~ proportional to the mean image (deficit/flux constant
~2e-3 over the PSF). A rolling-band / additive artefact: deficit roughly uniform per row, present in faint
pixels too. Also builds my own aperture light curve from the decorrelated pixels.
"""
import os
import numpy as np
from astropy.io import fits
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lcio

TC = 179.90
f = lcio._find("ktwo211924561-c05_lpd-targ.fits.gz")
h = fits.open(f)
d = h[1].data
t = d["TIME"] - 2167.0
fl = d["FLUX"].astype(float)
q = d["QUALITY"]
p1, p2 = d["POS_CORR1"], d["POS_CORR2"]
ok = np.isfinite(t) & np.all(np.isfinite(fl.reshape(len(t), -1)) | ~np.isfinite(fl[0].reshape(-1)), axis=1) \
    & np.isfinite(p1) & np.isfinite(p2) & ((q & (1 << 20)) == 0)
w = ok & (np.abs(t - TC) < 1.5)
inn = w & (np.abs(t - TC) < 0.035)    # 179.865..179.935 -> deepest 3-4 cadences
oot = w & (np.abs(t - TC) > 0.12)
print("n in", inn.sum(), "n oot", oot.sum(), "in-event times", np.round(t[inn], 4))

ny, nx = fl.shape[1:]
mean_img = np.nanmedian(fl[w], axis=0)
deficit = np.full((ny, nx), np.nan)
resid_lc = np.zeros((w.sum(), ny, nx))
X = lambda s: np.vstack([np.ones(s.sum()), p1[s], p2[s], p1[s] ** 2, p2[s] ** 2, p1[s] * p2[s], t[s] - TC]).T
for iy in range(ny):
    for ix in range(nx):
        y = fl[:, iy, ix]
        if not np.isfinite(mean_img[iy, ix]):
            continue
        c, *_ = np.linalg.lstsq(X(oot), y[oot], rcond=None)
        pred_w = X(w) @ c
        resid_lc[:, iy, ix] = y[w] - pred_w
        deficit[iy, ix] = np.mean((X(inn) @ c) - y[inn])

frac = deficit / mean_img
ap = mean_img > np.nanpercentile(mean_img, 80)
print("total deficit in brightest-20% aperture (e-/s):", np.nansum(deficit[ap]), " flux:", np.nansum(mean_img[ap]),
      " frac depth (ppt):", 1e3 * np.nansum(deficit[ap]) / np.nansum(mean_img[ap]))
faint = mean_img < np.nanpercentile(mean_img, 40)
print("faint pixels: total deficit", np.nansum(deficit[faint]), "e-/s over", faint.sum(), "pixels; per pixel",
      np.nanmean(deficit[faint]), "+-", np.nanstd(deficit[faint]) / np.sqrt(faint.sum()))
# correlation of the deficit with the PSF (mean image): linear fit deficit = a*img + b
m = np.isfinite(deficit) & np.isfinite(mean_img)
A = np.vstack([mean_img[m], np.ones(m.sum())]).T
(a, b0), *_ = np.linalg.lstsq(A, deficit[m], rcond=None)
res = deficit[m] - A @ np.array([a, b0])
print(f"deficit = {a*1e3:.3f}e-3 * image + {b0:.2f} e-/s/pix ; resid rms {np.std(res):.2f}; corr coef {np.corrcoef(mean_img[m], deficit[m])[0,1]:.3f}")
# row-wise mean deficit in faint pixels (rolling-band signature)
for iy in range(ny):
    row = faint[iy] & np.isfinite(deficit[iy])
    if row.sum():
        pass
# my own aperture LC
lc = np.nansum(resid_lc[:, ap], axis=1) / np.nansum(mean_img[ap])
tw = t[w]
np.save(os.path.join(lcio.HERE, "cache", "k2_pixel_aperture_lc.npy"), np.vstack([tw, lc]))

fig, ax = plt.subplots(1, 4, figsize=(18, 4.2))
im = ax[0].imshow(mean_img, origin="lower"); ax[0].set_title("C5 median image"); plt.colorbar(im, ax=ax[0])
vm = np.nanmax(np.abs(deficit))
im = ax[1].imshow(deficit, origin="lower", cmap="RdBu_r", vmin=-vm, vmax=vm); ax[1].set_title("in-event deficit (e-/s)"); plt.colorbar(im, ax=ax[1])
ax[1].contour(ap, levels=[0.5], colors="k", linewidths=0.8)
ax[2].plot(mean_img[m], deficit[m], "o", ms=3); xx = np.linspace(0, np.nanmax(mean_img), 10); ax[2].plot(xx, a * xx + b0, "r-")
ax[2].set_xlabel("pixel flux (e-/s)"); ax[2].set_ylabel("deficit (e-/s)"); ax[2].set_title("deficit vs PSF: slope %.2f ppt" % (a * 1e3))
ax[3].plot(tw, 1e3 * lc, "o", ms=3); ax[3].axvline(TC, color="r", alpha=0.3); ax[3].set_title("own pixel-decorrelated aperture LC (ppt)")
plt.tight_layout(); plt.savefig(os.path.join(lcio.HERE, "fig_k2_pixel_test.png"), dpi=100)
