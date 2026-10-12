"""
Is the Sector 21 dip of TIC 159614322 a third transit?

Context: the 179.93-d and 89.97-d aliases (n = 4, 8 of dT = 719.72 d) predict a transit at
BTJD 1883.76 +/- ~0.02. TESS S21 data end at the orbit-1 downlink gap (TESS-SPOC 1883.97),
just after the predicted egress (1883.94), so the alias pipeline called the window
'ambiguous' (no baseline after it). Checks here:

  1. Timing/depth fit: transit shape fixed from the S14+S41 joint fit (b=0 solution at P=89.97),
     free mid-time, depth scale and a linear pre-ingress baseline; per product (TESS-SPOC PDC
     and SAP, QLP). Delta-chi2 vs a no-transit model with the same baseline freedom.
  2. Pixel level (TESScut S21): target-aperture light curve, background, and 5 comparison stars;
     difference image (pre-ingress minus in-transit) with a Gaussian-fit source position vs the
     target's WCS position.
  3. False-alarm check: the same in/out test placed 0.4 d before the end of every orbit in all
     sectors of this star (14-16, 21-23, 41, 48, 49, 75, 76): how often does a 'drop' of this size
     appear just before a downlink gap?

Output: s21_check.json, fig_s21.png (diagnostic), next to this script.
Run:    python singles/writeups/tic159614322/s21_check.py
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import diff_image as di  # noqa: E402
import vet_candidate as vc  # noqa: E402

TIC, RA, DEC = 159614322, 220.00977, 65.484041
T_PRED, P8 = 1883.758, 89.9651
D = 8.773 / 24
OUT_LO, OUT_HI = 1883.00, 1883.55     # after the 1882.97 momentum dump, before ingress
IN_LO, IN_HI = 1883.62, 1883.96


def shape():
    a = json.load(open(HERE / "analysis.json"))
    b = a["best_b_P90"]
    return b["rp"], b["b"], b["aR"]


def fit_product(p, rp, b, aR, col="flux"):
    t = p["time"]
    m = (t > OUT_LO) & (t < 1884.2)
    t, f = t[m], p[col][m]
    ok = np.isfinite(f)
    t, f = t[ok], f[ok] / np.nanmedian(f[ok][t[ok] < OUT_HI])
    pre = t < OUT_HI
    sig = np.std(f[pre] - np.polyval(np.polyfit(t[pre], f[pre], 1), t[pre]))

    def model(x):
        T0, s, c0, c1 = x
        tr = vc.model(t, T0, P8, rp, b, aR, p["exptime"])
        return (c0 + c1 * (t - 1883.5)) * (1 - s * (1 - tr))

    def r_tr(x):
        return (f - model(x)) / sig

    best = None
    for T0g in np.arange(1883.55, 1883.95, 0.02):
        s_ = least_squares(r_tr, [T0g, 1.0, 1.0, 0.0],
                           bounds=([1883.3, 0.0, 0.9, -1], [1884.2, 3.0, 1.1, 1]))
        if best is None or s_.cost < best.cost:
            best = s_
    # no transit, same baseline freedom (linear) -- plus a version allowing a quadratic
    nt = least_squares(lambda x: (f - (x[0] + x[1] * (t - 1883.5))) / sig, [1, 0])
    nt2 = least_squares(lambda x: (f - (x[0] + x[1] * (t - 1883.5) + x[2] * (t - 1883.5) ** 2)) / sig, [1, 0, 0])
    J = best.jac
    cov = np.linalg.inv(J.T @ J) * max(2 * best.cost / (len(t) - 4), 1)
    return dict(product=p["product"], column=col, n=int(len(t)), sigma_ppm=float(sig * 1e6),
                T0=float(best.x[0]), T0_err=float(np.sqrt(cov[0, 0])),
                depth_scale=float(best.x[1]), depth_scale_err=float(np.sqrt(cov[1, 1])),
                chi2_transit=float(2 * best.cost), chi2_flat_linear=float(2 * nt.cost),
                chi2_flat_quadratic=float(2 * nt2.cost),
                dchi2_vs_linear=float(2 * nt.cost - 2 * best.cost),
                dchi2_vs_quadratic=float(2 * nt2.cost - 2 * best.cost),
                data=(t.tolist(), f.tolist()), model=model(best.x).tolist())


def pixel_level():
    path = di.get_cutout(TIC, RA, DEC, 21)
    with fits.open(path) as h:
        tab, wcs = h[1].data, WCS(h[2].header)
    t, fl, bkg, q = tab["TIME"], tab["FLUX"].astype(float), tab["FLUX_BKG"].astype(float), tab["QUALITY"]
    good = np.isfinite(t) & (q == 0) & np.all(np.isfinite(fl.reshape(len(t), -1)), axis=1)
    t, fl, bkg = t[good], fl[good], bkg[good].reshape(good.sum(), -1).mean(1)
    tx, ty = wcs.all_world2pix(RA, DEC, 0)
    g = di.gaia_field(RA, DEC)
    gx, gy = wcs.all_world2pix(g.ra.values, g.dec.values, 0)
    ny, nx = fl.shape[1:]

    def ap_lc(x, y, r=1):
        xi, yi = int(round(float(x))), int(round(float(y)))
        if not (r <= xi < nx - r and r <= yi < ny - r):
            return None
        return fl[:, yi - r:yi + r + 1, xi - r:xi + r + 1].sum(axis=(1, 2))

    out = (t > OUT_LO) & (t < OUT_HI)
    inn = (t > IN_LO) & (t < IN_HI)

    def drop(lc):
        o = lc[out]
        return float(1 - lc[inn].mean() / o.mean()), float(o.std() / o.mean() / np.sqrt(inn.sum()))

    res = dict(n_in=int(inn.sum()), n_out=int(out.sum()))
    tgt = ap_lc(tx, ty)
    res["target_drop"] = drop(tgt)
    res["background_change_e_s_per_pix"] = float(bkg[inn].mean() - bkg[out].mean())
    res["background_level_e_s_per_pix"] = float(bkg[out].mean())
    # comparison stars: brightest Gaia stars in the cutout, > 3 px from the target
    comps = []
    order = np.argsort(g.phot_g_mean_mag.values)
    for i in order:
        if np.hypot(gx[i] - tx, gy[i] - ty) < 3 or g.phot_g_mean_mag.values[i] > 14:
            continue
        lc = ap_lc(gx[i], gy[i])
        if lc is None:
            continue
        d, e = drop(lc)
        comps.append(dict(G=float(g.phot_g_mean_mag.values[i]), dx=float(gx[i] - tx), dy=float(gy[i] - ty),
                          drop=d, err=e))
        if len(comps) == 6:
            break
    res["comparison_stars"] = comps
    # difference image: pre-ingress (linear-free; constant) minus in-transit
    diff = fl[out].mean(0) - fl[inn].mean(0)
    err = fl[out].std(0) * np.sqrt(1 / inn.sum() + 1 / out.sum())
    yy, xx = np.mgrid[:ny, :nx]
    sel = (np.abs(xx - tx) <= 3) & (np.abs(yy - ty) <= 3)
    p0 = [diff[sel].max(), tx, ty, 0.8, 0]
    fit = least_squares(lambda p: ((di.gaussian(p, yy, xx) - diff) / np.maximum(err, 1e-3))[sel], p0,
                        bounds=([0, tx - 3, ty - 3, 0.3, -np.inf], [np.inf, tx + 3, ty + 3, 3, np.inf]))
    J = fit.jac
    cov = np.linalg.inv(J.T @ J) * max(2 * fit.cost / (sel.sum() - 5), 1)
    off_px = np.hypot(fit.x[1] - tx, fit.x[2] - ty)
    off_err = np.sqrt(cov[1, 1] + cov[2, 2])
    pk = int(np.argmax(np.where(sel, diff, -np.inf)))
    res["diff_image"] = dict(offset_px=float(off_px), offset_arcsec=float(off_px * 21), offset_err_arcsec=float(off_err * 21),
                             peak_snr=float((diff / err).flat[pk]), fit=[float(v) for v in fit.x])
    res["_img"] = dict(diff=diff.tolist(), tx=float(tx), ty=float(ty),
                       gx=gx[g.phot_g_mean_mag.values < 15].tolist(), gy=gy[g.phot_g_mean_mag.values < 15].tolist())
    res["_lc"] = dict(t=t.tolist(), target=(tgt / tgt[out].mean()).tolist(), bkg=bkg.tolist())
    return res


def false_alarm(prods):
    """Same in/out test 0.4 d before the end of every contiguous orbit (gap > 0.5 d follows)."""
    rows = []
    for p in prods:
        if p["product"] not in ("TESS-SPOC", "SPOC_120", "QLP"):
            continue
        t, f = p["time"], p["flux"]
        ok = np.isfinite(f)
        t, f = t[ok], f[ok]
        gaps = np.where(np.diff(t) > 0.5)[0]
        ends = list(t[gaps]) + [t[-1]]
        for te in ends:
            inn = (t > te - 0.34) & (t <= te)
            out = (t > te - 0.96) & (t < te - 0.41)
            if inn.sum() < 5 or out.sum() < 5:
                continue
            o = f[out]
            d = 1 - f[inn].mean() / o.mean()
            e = o.std() / o.mean() / np.sqrt(inn.sum())
            rows.append(dict(sector=int(p["sector"]), product=p["product"], t_end=float(te),
                             drop_ppm=float(d * 1e6), snr=float(d / e)))
    return rows


def main():
    rp, b, aR = shape()
    prods = vc.all_products(TIC)
    s21 = [p for p in prods if p["sector"] == 21]
    fits_ = []
    for p in s21:
        for col in ("flux", "sap_flux"):
            r = fit_product(p, rp, b, aR, col)
            fits_.append(r)
    pix = pixel_level()
    fa = false_alarm(prods)
    res = dict(T_pred_n4_n8=T_PRED, shape=dict(rp=rp, b=b, aR=aR),
               fits=[{k: v for k, v in r.items() if k not in ("data", "model")} for r in fits_],
               pixel={k: v for k, v in pix.items() if not k.startswith("_")},
               false_alarm_orbit_ends=fa)
    (HERE / "s21_check.json").write_text(json.dumps(res, indent=1))

    fig = plt.figure(figsize=(11, 7))
    for i, r in enumerate(fits_):
        ax = fig.add_subplot(3, 2, i + 1)
        t, f = np.array(r["data"][0]), np.array(r["data"][1])
        ax.plot(t, (f - 1) * 1e3, "k.", ms=3)
        ax.plot(t, (np.array(r["model"]) - 1) * 1e3, "-", color="#d9682b")
        ax.axvline(T_PRED, color="g", ls=":")
        ax.set_title(f"S21 {r['product']} {r['column']}: T0={r['T0']:.3f}+/-{r['T0_err']:.3f}, "
                     f"scale {r['depth_scale']:.2f}+/-{r['depth_scale_err']:.2f}, dchi2 {r['dchi2_vs_linear']:.0f}", fontsize=7.5)
    ax = fig.add_subplot(3, 2, 5)
    lc = pix["_lc"]
    t = np.array(lc["t"])
    m = (t > 1882.5) & (t < 1884.2)
    ax.plot(t[m], (np.array(lc["target"])[m] - 1) * 1e3, "k.", ms=3, label="target 3x3 px (TESScut)")
    bk = np.array(lc["bkg"])
    ax.plot(t[m], (bk[m] - bk[m].mean()) / 10, "b.", ms=2, label="background /10 (e/s/px)")
    ax.axvspan(OUT_LO, OUT_HI, color="grey", alpha=0.15)
    ax.axvspan(IN_LO, IN_HI, color="orange", alpha=0.2)
    ax.legend(fontsize=6)
    ax = fig.add_subplot(3, 2, 6)
    im = pix["_img"]
    ax.imshow(np.array(im["diff"]), origin="lower", cmap="magma")
    ax.plot(im["tx"], im["ty"], "o", mfc="none", mec="cyan", ms=12)
    ax.plot(im["gx"], im["gy"], "s", mfc="none", mec="w", ms=5)
    fx, fy = pix["diff_image"]["fit"][1:3]
    ax.plot(fx, fy, "x", color="lime", ms=8)
    ax.set_xlim(im["tx"] - 5, im["tx"] + 5)
    ax.set_ylim(im["ty"] - 5, im["ty"] + 5)
    ax.set_title(f"S21 diff image: source {pix['diff_image']['offset_arcsec']:.1f}\" +/- "
                 f"{pix['diff_image']['offset_err_arcsec']:.1f}\" from target, SNR {pix['diff_image']['peak_snr']:.1f}", fontsize=7.5)
    fig.tight_layout()
    fig.savefig(HERE / "fig_s21_diagnostic.png", dpi=130)
    print(json.dumps({k: v for k, v in res.items() if k != "false_alarm_orbit_ends"}, indent=1))
    fa_s = sorted(fa, key=lambda r: -r["drop_ppm"])
    print("orbit-end drops (top 8 of", len(fa), "):", [(r["sector"], r["product"], round(r["drop_ppm"]), round(r["snr"], 1)) for r in fa_s[:8]])


if __name__ == "__main__":
    main()
