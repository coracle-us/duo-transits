"""
Leads triage, step 3b: pixel difference images for EPIC 245929407 (off-target call from the light-curve centroid,
Gaia neighbour 4.8" away) and EPIC 211894612 (on target by centroid; EB by size). Rules fixed before the first run.

K2 (target pixel files, as k2singles/independent_211924561/s09_k2_pixels.py):
  per pixel, within +-max(1.5 d, 2D) of the event, fit out-of-event cadences (|t-T| > D/2 + 0.1 d, no thruster
  flag) with a quadratic in (POS_CORR1, POS_CORR2) plus a linear time term; deficit image = mean(predicted -
  observed) over in-event cadences (|t-T| < 0.35 D); error from the out-of-event residual scatter.
  (a) free Gaussian fit to the deficit -> fitted source position; offset from the target's Gaia position (TPF WCS)
      and from the out-of-event image centroid (same Gaussian, WCS-free).
  (b) template fit: the median image is fitted with Gaussians at the Gaia positions of the target and every
      neighbour with G < G_target + 5 inside the stamp (common shift, common width, free amplitudes, background);
      the deficit is then fitted with the same templates (amplitudes and background free, linear). The depth each
      star would need is deficit amplitude / image amplitude.
  On target if the fitted source is within 3 sigma (and 1 px) of the target and the template fit gives the target a
  deficit > 3 sigma with every neighbour's deficit consistent with zero (< 3 sigma) or implying a depth > 100%.
  Off target if the fitted source is > 3 sigma and > 1 px from the target, or the target's template deficit is
  < 3 sigma while a neighbour's is > 3 sigma at a depth <= 100%.
TESS (TESScut FFI cutouts): the difference image of scripts/diff_image.py analyse() (per-pixel linear trend from
  out-of-transit frames, diff = predicted - observed in transit), fitted with the same Gaussian, but started at the
  peak of diff/error instead of the peak of diff. Changed after the first run: in S44 of 211894612 the largest
  absolute difference is a saturated G = 9.8 star on the cutout edge (diff/err ~ 3), and the fit started there,
  while the target pixel has diff/err ~ 70. Offset from the target's Gaia position (cutout WCS), 21"/px.

Output: k2singles/leads_triage/pixels.csv, figs/pix_*.png
Usage:  python k2singles/leads_triage/pixels_diff.py
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from astropy.io import fits  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import diff_image as di  # noqa: E402

PIX = HERE / "data" / "pixels"
FIG = HERE / "figs"
K2_SCALE = 3.98
# events: (EPIC, TIC, mission, campaign/sector, T_BTJD, D_days); D from the search's refined fit
EVENTS = [
    ("245929407", 49733422, "K2", 12, 760.123, 3.87 / 24),
    ("245929407", 49733422, "TESS", 69, 3196.343, 3.87 / 24),
    ("211894612", 14678155, "K2", 5, 182.855, 11.94 / 24),
    ("211894612", 14678155, "K2", 18, 1297.244, 10.78 / 24),
    ("211894612", 14678155, "TESS", 44, 2523.029, 8.78 / 24),
]


def gauss(x0, y0, s, yy, xx):
    return np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * s ** 2))


def gaia_cone(ra, dec, r_arcsec):
    from astroquery.gaia import Gaia
    q = f"""SELECT source_id, ra, dec, pmra, pmdec, phot_g_mean_mag,
            DISTANCE(POINT({ra},{dec}), POINT(ra,dec))*3600 AS sep FROM gaiadr3.gaia_source
            WHERE 1=CONTAINS(POINT(ra,dec), CIRCLE({ra},{dec},{r_arcsec / 3600})) ORDER BY sep"""
    return Gaia.launch_job(q).get_results().to_pandas()


def k2_event(epic, camp, T, D, gaia):
    f = next((PIX / "k2" / epic).rglob(f"ktwo{epic}-c{camp:02d}_lpd-targ.fits*"))
    import lightkurve as lk
    wcs = lk.read(str(f), quality_bitmask="none").wcs
    with fits.open(f) as h:
        d = h[1].data
        t = d["TIME"] - 2167.0
        fl = d["FLUX"].astype(float)
        q = d["QUALITY"]
        p1, p2 = d["POS_CORR1"], d["POS_CORR2"]
        date = h[0].header.get("DATE-OBS") or h[1].header.get("DATE-OBS")
    ny, nx = fl.shape[1:]
    pixok = np.isfinite(np.nanmedian(fl, axis=0))
    ok = np.isfinite(t) & np.isfinite(p1) & np.isfinite(p2) & ((q & (1 << 20)) == 0) & \
        np.all(np.isfinite(fl[:, pixok]), axis=1)
    B = max(1.5, 2 * D)
    w = ok & (np.abs(t - T) < D / 2 + B)
    inn = w & (np.abs(t - T) < 0.35 * D)
    oot = w & (np.abs(t - T) > D / 2 + 0.1)
    X = lambda s: np.vstack([np.ones(s.sum()), p1[s], p2[s], p1[s] ** 2, p2[s] ** 2, p1[s] * p2[s], t[s] - T]).T
    mean_img = np.nanmedian(fl[oot], axis=0)
    deficit = np.full((ny, nx), np.nan)
    derr = np.full((ny, nx), np.nan)
    for iy in range(ny):
        for ix in range(nx):
            if not pixok[iy, ix]:
                continue
            y = fl[:, iy, ix]
            c, *_ = np.linalg.lstsq(X(oot), y[oot], rcond=None)
            r = y[oot] - X(oot) @ c
            deficit[iy, ix] = np.mean(X(inn) @ c - y[inn])
            derr[iy, ix] = np.std(r) * np.sqrt(1 / inn.sum() + 1 / oot.sum())
    # Gaia positions at the K2 epoch (proper motion from J2016.0)
    ep = float(date[:4]) + (float(date[5:7]) - 0.5) / 12 if date else 2016.0
    dt_yr = ep - 2016.0
    gra = gaia.ra + gaia.pmra.fillna(0) / 3.6e6 / np.cos(np.radians(gaia.dec)) * dt_yr
    gde = gaia.dec + gaia.pmdec.fillna(0) / 3.6e6 * dt_yr
    gx, gy = wcs.all_world2pix(gra.values, gde.values, 0)
    yy, xx = np.mgrid[0:ny, 0:nx]
    m = pixok & np.isfinite(deficit)
    tx, ty = float(gx[0]), float(gy[0])
    # (a) free Gaussian on the deficit, and on the median image for the WCS-free reference
    def free_fit(img, err, x0, y0):
        res = least_squares(lambda p: ((p[0] * gauss(p[1], p[2], p[3], yy, xx) + p[4] - img) / err)[m],
                            [max(np.nanmax(img[m]), 1e-3), x0, y0, 0.8, 0.0],
                            bounds=([0, -1, -1, 0.3, -np.inf], [np.inf, nx, ny, 3.0, np.inf]))
        J = res.jac
        cov = np.linalg.pinv(J.T @ J) * max(1.0, np.sum(res.fun ** 2) / (res.fun.size - 5))
        return res.x, np.sqrt(np.diag(cov))
    iy0, ix0 = np.unravel_index(np.nanargmax(np.where(m, deficit, -np.inf)), deficit.shape)
    pdef, edef = free_fit(deficit, derr, float(ix0), float(iy0))
    pimg, _ = free_fit(mean_img, np.sqrt(np.abs(mean_img)) + 1.0, tx, ty)
    off_t = np.hypot(pdef[1] - tx, pdef[2] - ty)
    off_c = np.hypot(pdef[1] - pimg[1], pdef[2] - pimg[2])
    e_off = float(np.hypot(edef[1], edef[2]))
    # (b) template fit: target + neighbours inside the stamp
    inside = (gx > -1) & (gx < nx) & (gy > -1) & (gy < ny) & (gaia.phot_g_mean_mag < gaia.phot_g_mean_mag.iloc[0] + 5)
    idx = np.where(inside)[0]
    if 0 not in idx:
        idx = np.r_[0, idx]
    k = len(idx)

    def model_img(p):
        dx, dy, s, bg = p[:4]
        amps = p[4:]
        return sum(a * gauss(gx[i] + dx, gy[i] + dy, s, yy, xx) for a, i in zip(amps, idx)) + bg
    eimg = np.sqrt(np.abs(mean_img)) + 1.0
    flux0 = np.nansum(mean_img[m])
    a0 = [flux0 * 10 ** (-0.4 * (gaia.phot_g_mean_mag.iloc[i] - gaia.phot_g_mean_mag.iloc[0])) / 4 for i in idx]
    rimg = least_squares(lambda p: ((model_img(p) - mean_img) / eimg)[m], [0, 0, 0.9, 0] + a0,
                         bounds=([-2, -2, 0.3, -np.inf] + [0] * k, [2, 2, 3.0, np.inf] + [np.inf] * k))
    dx, dy, s = rimg.x[:3]
    A_img = rimg.x[4:]
    T_ = np.column_stack([gauss(gx[i] + dx, gy[i] + dy, s, yy, xx)[m] for i in idx] + [np.ones(m.sum())])
    Wt = 1 / derr[m]
    coef, *_ = np.linalg.lstsq(T_ * Wt[:, None], deficit[m] * Wt, rcond=None)
    cov = np.linalg.pinv((T_ * Wt[:, None]).T @ (T_ * Wt[:, None]))
    resid = (deficit[m] - T_ @ coef) * Wt
    cov *= max(1.0, np.sum(resid ** 2) / (m.sum() - k - 1))
    ecoef = np.sqrt(np.diag(cov))
    stars = []
    for j, i in enumerate(idx):
        stars.append(dict(source_id=str(int(gaia.source_id.iloc[i])), sep=round(float(gaia.sep.iloc[i]), 1),
                          dG=round(float(gaia.phot_g_mean_mag.iloc[i] - gaia.phot_g_mean_mag.iloc[0]), 2),
                          deficit_sigma=round(float(coef[j] / ecoef[j]), 1),
                          implied_depth=round(float(coef[j] / A_img[j]), 4) if A_img[j] > 0 else None))
    ap = mean_img > np.nanpercentile(mean_img[m], 80)
    depth_ap = float(np.nansum(deficit[ap]) / np.nansum(mean_img[ap]))
    out = dict(EPIC=epic, mission="K2", sector=camp, T=T, n_in=int(inn.sum()), n_oot=int(oot.sum()),
               aperture_depth=round(depth_ap, 5),
               fit_offset_from_target_arcsec=round(off_t * K2_SCALE, 2),
               fit_offset_from_image_centroid_arcsec=round(off_c * K2_SCALE, 2),
               fit_offset_err_arcsec=round(e_off * K2_SCALE, 2), fit_offset_sigma=round(off_t / e_off, 1) if e_off else None,
               wcs_vs_image_centroid_arcsec=round(float(np.hypot(pimg[1] - tx, pimg[2] - ty)) * K2_SCALE, 2),
               template_psf_sigma_px=round(float(s), 2), template_shift_arcsec=round(float(np.hypot(dx, dy)) * K2_SCALE, 2),
               stars=stars)
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
    for a, img, lab in zip(ax, [mean_img, deficit, deficit / derr], ["median out-of-event image", "deficit (e-/s)", "deficit / error"]):
        im = a.imshow(img, origin="lower", cmap="viridis")
        plt.colorbar(im, ax=a, fraction=0.046)
        a.scatter(gx[idx[1:]], gy[idx[1:]], s=60, facecolors="none", edgecolors="w")
        a.plot(tx, ty, "r+", ms=14, mew=2)
        a.plot(pdef[1], pdef[2], "mx", ms=10, mew=2)
        a.set_title(lab, fontsize=9)
    fig.suptitle(f"EPIC {epic} K2 C{camp} event {T:.3f}: red + target (Gaia), magenta x fitted deficit, circles Gaia "
                 f"neighbours G<target+5; offset {off_t * K2_SCALE:.1f}\" +- {e_off * K2_SCALE:.1f}\"", fontsize=9)
    fig.tight_layout()
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / f"pix_{epic}_k2c{camp}.png", dpi=85)
    plt.close(fig)
    return out


def verdict_k2(r):
    tgt = r["stars"][0]
    others = r["stars"][1:]
    nb_real = [s for s in others if s["deficit_sigma"] >= 3 and s["implied_depth"] is not None and s["implied_depth"] <= 1]
    off_fit = r["fit_offset_sigma"] is not None and r["fit_offset_sigma"] > 3 and r["fit_offset_from_target_arcsec"] > K2_SCALE
    on_fit = (r["fit_offset_sigma"] is not None and r["fit_offset_sigma"] <= 3) or r["fit_offset_from_target_arcsec"] <= K2_SCALE
    if off_fit or (tgt["deficit_sigma"] < 3 and nb_real):
        return "off-target"
    if on_fit and tgt["deficit_sigma"] >= 3 and not nb_real:
        return "on-target"
    return "ambiguous"


def tess_event(tic, sec, T, D, gaia):
    path = next((PIX / "tesscut" / str(tic)).glob(f"*s{sec:04d}*.fits"))
    from astropy.wcs import WCS
    with fits.open(path) as h:
        tab = h[1].data
        wcs = WCS(h[2].header)
        t, flux, q = tab["TIME"], tab["FLUX"], tab["QUALITY"]
    good = np.isfinite(t) & (q == 0) & np.all(np.isfinite(flux.reshape(len(t), -1)), axis=1)
    t, flux = t[good], flux[good]
    n = flux.shape[1]
    dt = t - T
    inn = np.abs(dt) < 0.35 * D
    out = (np.abs(dt) > D / 2 + 0.05) & (np.abs(dt) < D / 2 + 0.6)
    A = np.column_stack([np.ones(out.sum()), dt[out]])
    Y = flux[out].reshape(out.sum(), -1)
    coef, *_ = np.linalg.lstsq(A, Y, rcond=None)
    diff = ((np.column_stack([np.ones(inn.sum()), dt[inn]]) @ coef) - flux[inn].reshape(inn.sum(), -1)).mean(axis=0)
    derr = (Y - A @ coef).std(axis=0) * np.sqrt(1 / inn.sum() + 1 / out.sum())
    diff, derr = diff.reshape(n, n), derr.reshape(n, n)
    tx, ty = (float(v) for v in wcs.all_world2pix(gaia.ra.iloc[0], gaia.dec.iloc[0], 0))
    gx, gy = wcs.all_world2pix(gaia.ra.values, gaia.dec.values, 0)
    yy, xx = np.mgrid[0:n, 0:n]
    iy, ix = np.unravel_index(np.argmax(diff / derr), diff.shape)
    r = least_squares(lambda p: ((di.gaussian(p, yy, xx) - diff) / derr).ravel(),
                      [max(float(diff[iy, ix]), 1e-6), float(np.clip(ix, 0.01, n - 1.01)), float(np.clip(iy, 0.01, n - 1.01)), 0.8, 0.0],
                      bounds=([0, 0, 0, 0.4, -np.inf], [np.inf, n - 1, n - 1, 2.0, np.inf]))
    cov = np.linalg.pinv(r.jac.T @ r.jac) * max(1, np.sum(r.fun ** 2) / (r.fun.size - 5))
    fx, fy = r.x[1], r.x[2]
    off, eoff = float(np.hypot(fx - tx, fy - ty)), float(np.hypot(np.sqrt(cov[1, 1]), np.sqrt(cov[2, 2])))
    j = int(np.argmin(np.hypot(gx - fx, gy - fy)))
    bx, by = int(round(tx)), int(round(ty))
    o, e = off * 21, eoff * 21
    res = dict(EPIC=None, mission="TESS", sector=sec, T=T, n_in=int(inn.sum()), n_oot=int(out.sum()),
               target_pixel_diff_sigma=round(float(diff[by, bx] / derr[by, bx]), 1),
               peak_diff_sigma=round(float(np.max(diff / derr)), 1),
               fit_offset_from_target_arcsec=round(o, 2), fit_offset_err_arcsec=round(e, 2),
               fit_offset_sigma=round(o / e, 1) if e > 0 else None,
               nearest_gaia_to_source=dict(source_id=str(int(gaia.source_id.iloc[j])), G=round(float(gaia.phot_g_mean_mag.iloc[j]), 2),
                                           sep_from_target_arcsec=round(float(gaia.sep.iloc[j]), 1),
                                           dist_from_fit_arcsec=round(float(np.hypot(gx[j] - fx, gy[j] - fy)) * 21, 1)))
    # a significant offset under half a TESS pixel cannot tell the target from a star a few arcsec away: "unresolved"
    # (label added after the first run; 245929407's S69 fit lies 1.7" from its 4.8" neighbour)
    res["verdict"] = ("off-target" if (o / e > 3 and o > 21) else "on-target" if o / e <= 3
                      else "unresolved (< 0.5 px)" if o <= 21 / 2 else "ambiguous")
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.3))
    for ax, img, lab in zip(axs, [np.nanmedian(flux, axis=0).reshape(n, n), diff, diff / derr],
                            ["median image", "difference (out - in)", "difference / error"]):
        im = ax.imshow(img, origin="lower", cmap="viridis")
        plt.colorbar(im, ax=ax, fraction=0.046)
        b = gaia.phot_g_mean_mag < gaia.phot_g_mean_mag.iloc[0] + 5
        ax.scatter(gx[b], gy[b], s=30, facecolors="none", edgecolors="w", lw=0.8)
        ax.plot(tx, ty, "r+", ms=14, mew=2)
        ax.plot(fx, fy, "mx", ms=10, mew=2)
        ax.set_xlim(-0.5, n - 0.5)
        ax.set_ylim(-0.5, n - 0.5)
        ax.set_title(lab, fontsize=9)
    fig.suptitle(f"TIC {tic} S{sec} event {T:.3f}: red + target, magenta x fitted source (offset {o:.1f}\" +- {e:.1f}\"); "
                 "circles Gaia G<target+5", fontsize=9)
    fig.tight_layout()
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / f"pix_{tic}_s{sec}.png", dpi=85)
    plt.close(fig)
    return res


def main():
    cat = pd.read_csv(HERE.parent / "data" / "k2_singles.csv", dtype={"EPIC": str}).drop_duplicates("EPIC").set_index("EPIC")
    look = pd.read_csv(HERE / "catalog_lookup.csv", dtype={"EPIC": str, "gaia_source": str}).set_index("EPIC")
    rows = []
    for epic, tic, mis, sec, T, D in EVENTS:
        r = cat.loc[epic]
        g = gaia_cone(r.ra, r.dec, 120 if mis == "K2" else 15 * 21 / 2 * 1.5)
        sid = look.loc[epic, "gaia_source"]
        g = pd.concat([g[g.source_id.astype(str) == sid], g[g.source_id.astype(str) != sid]], ignore_index=True)
        g["sep"] = np.hypot((g.ra - g.ra.iloc[0]) * np.cos(np.radians(g.dec.iloc[0])), g.dec - g.dec.iloc[0]) * 3600
        if mis == "K2":
            out = k2_event(epic, sec, T, D, g)
            out["verdict"] = verdict_k2(out)
        else:
            gt = g[g.phot_g_mean_mag < 19].reset_index(drop=True)
            out = tess_event(tic, sec, T, D, gt)
            out["EPIC"] = epic
        print(json.dumps(out, default=str))
        rows.append(out)
    df = pd.DataFrame(rows)
    df["stars"] = df.get("stars", pd.Series(dtype=object)).apply(lambda v: json.dumps(v) if isinstance(v, list) else "")
    df.to_csv(HERE / "pixels.csv", index=False)


if __name__ == "__main__":
    main()
