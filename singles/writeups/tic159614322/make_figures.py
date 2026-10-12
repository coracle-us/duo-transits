"""
Figures for the TIC 159614322 (CTOI 159614322.01) write-up.
  fig_transits.png  the three transits (S14, S21, S41) with one transit shape; S21 shows
                    TESS-SPOC SAP, QLP and our own TESScut photometry, which end at the
                    orbit-1 downlink gap just before the predicted egress
  fig_aliases.png   predicted epochs of P = dT/n (n = 1..12) against TESS coverage
  fig_maps.png      per-pixel transit-amplitude maps (S14, S21, S41): each pixel's light curve
                    fitted with a quadratic baseline + the transit shape, so slow pointing drift
                    goes into the baseline; source position from a Gaussian fit
Also writes maps.json (source offsets).
Run:  python singles/writeups/tic159614322/analysis.py, s21_check.py, ephemeris.py, then this.
"""
import glob
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
import alias_pipeline as ap  # noqa: E402
import diff_image as di  # noqa: E402
import vet_candidate as vc  # noqa: E402

TIC, RA, DEC = 159614322, 220.00977, 65.484041
INK, MUTED, ACCENT, GOOD, BAD = "#1d2b3a", "#5d6b7c", "#d9682b", "#2b8a5e", "#b8323a"
EPH = json.load(open(HERE / "ephemeris.json"))
SH = json.load(open(HERE / "analysis.json"))["best_b_P90"]
TR = EPH["transits"]
P8 = EPH["ephemerides"]["n8"]["P"]


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=8, colors=MUTED)


def own_s21():
    """Background-subtracted 2-px-radius aperture photometry on the S21 TESScut cutout."""
    f = glob.glob(str(ROOT / "data/raw/tesscut/159614322/**/*s0021*.fits"), recursive=True)[0]
    with fits.open(f) as h:
        tab, wcs = h[1].data, WCS(h[2].header)
    t, fl, q = tab["TIME"], tab["FLUX"].astype(float), tab["QUALITY"]
    ok = np.isfinite(t) & (q == 0) & np.all(np.isfinite(fl.reshape(len(t), -1)), axis=1)
    t, fl = t[ok], fl[ok]
    tx, ty = [float(v) for v in wcs.all_world2pix(RA, DEC, 0)]
    yy, xx = np.mgrid[:fl.shape[1], :fl.shape[2]]
    r = np.hypot(xx - tx, yy - ty)
    bkg = np.median(fl[:, r > 5.5], axis=1)
    lc = (fl[:, r <= 2.0] - bkg[:, None]).sum(1)
    base = (t > 1883.0) & (t < 1883.55)
    return t, lc / lc[base].mean(), bkg


def fig_transits():
    prods = vc.all_products(TIC)
    get = lambda s, pr: next(q for q in prods if q["sector"] == s and q["product"] == pr)  # noqa: E731
    fig, axs = plt.subplots(1, 3, figsize=(10, 3.3), sharey=True)
    tt = np.linspace(-14, 14, 600)
    D = 8.773 / 24
    for ax, r, lab in zip(axs, TR, ("Sector 14 (Jul 2019)", "Sector 21 (Jan 2020)", "Sector 41 (Jul 2021)")):
        Tm = r["T"]
        if r["sector"] == 21:
            for pr, col, mk, off, name in (("TESS-SPOC", INK, "o", 0, "TESS-SPOC SAP"), ("QLP", "#3d7ab8", "s", -6, "QLP")):
                p = get(21, pr)
                c = "sap_flux" if pr == "TESS-SPOC" else "flux"
                m = (p["time"] > Tm - 0.8) & (p["time"] < 1884.0) & np.isfinite(p[c])
                base = m & (p["time"] < 1883.55)
                f = p[c][m] / np.median(p[c][base])
                ax.plot((p["time"][m] - Tm) * 24, (f - 1) * 1e3 + off, mk, ms=3, color=col, label=name)
            t, f, bkg = own_s21()
            m = (t > Tm - 0.8) & (t < 1884.0)
            ax.plot((t[m] - Tm) * 24, (f[m] - 1) * 1e3 - 12, "^", ms=3, color=GOOD, label="own TESScut")
            for off in (0, -6, -12):
                ax.plot(tt, (vc.model(Tm + tt / 24, Tm, P8, SH["rp"], SH["b"], SH["aR"], 120) - 1) * 1e3 * r["depth_scale"] + off,
                        "-", color=ACCENT, lw=1.3)
            ax.axvspan((1883.975 - Tm) * 24, 14, color="#e6e9ee", zorder=0)
            ax.text((1883.975 - Tm) * 24 + 0.3, 1.2, "downlink\ngap", fontsize=7, color=MUTED, va="top")
            ax.legend(fontsize=6.5, loc="lower right", frameon=True, framealpha=0.9)
        else:
            pr = "TESS-SPOC" if r["sector"] == 14 else "SPOC_120"
            p = get(r["sector"], pr)
            w = vc.window(p, Tm, D)
            t, f = w[0], w[1]
            x = (t - Tm) * 24
            if p["exptime"] >= 600:
                ax.plot(x, (f - 1) * 1e3, "o", ms=3, color=INK)
            else:
                ax.plot(x, (f - 1) * 1e3, ".", ms=1.5, color="#9fb3c8")
                bins = np.arange(-14, 14.01, 0.5)
                idx = np.digitize(x, bins)
                ax.plot([x[idx == i].mean() for i in np.unique(idx)], [((f - 1) * 1e3)[idx == i].mean() for i in np.unique(idx)],
                        "o", ms=3, color=INK)
            ax.plot(tt, (vc.model(Tm + tt / 24, Tm, P8, SH["rp"], SH["b"], SH["aR"], 120) - 1) * 1e3 * r["depth_scale"],
                    "-", color=ACCENT, lw=1.3)
        ax.set_title(f"{lab}\nBTJD {Tm:.3f} +/- {r['T_err']:.3f}", fontsize=8.5, color=INK)
        ax.set_xlim(-14, 14)
        ax.set_xlabel("hours from mid-transit", fontsize=8.5)
        style(ax)
    axs[0].set_ylim(-17, 2.5)
    axs[0].set_ylabel("relative flux (ppt; S21 curves offset)", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(HERE / "fig_transits.png", dpi=170)
    plt.close(fig)


def fig_aliases():
    lc = ap.load_lightcurves(TIC)
    spans = lc.groupby("sector").time.agg(["min", "max"])
    T1, T3 = TR[0]["T"], TR[2]["T"]
    dT = T3 - T1
    am = json.load(open(HERE / "analysis.json"))["alias_map"]
    fig, ax = plt.subplots(figsize=(10, 3.6))
    for s, (a, b) in spans.iterrows():
        if a > 3500:
            continue
        ax.axvspan(a, b, color="#e6e9ee", zorder=0)
        ax.text((a + b) / 2, 12.85 + 0.35 * (list(spans.index).index(s) % 2), f"S{s}", fontsize=6, ha="center", color=MUTED)
    for n in range(1, 13):
        P = dT / n
        ks = np.arange(0, n + 1)
        ks = np.arange(np.floor((1650 - T1) / P), np.ceil((3450 - T1) / P) + 1)
        for k in ks:
            T = T1 + k * P
            if not 1650 < T < 3450:
                continue
            v = [r for r in am if r["n"] == n and abs(r["T"] - T) < 0.3]
            if abs(T - T1) < 0.3 or abs(T - T3) < 0.3:
                ax.plot(T, n, "o", color=INK, ms=6)
            elif abs(T - TR[1]["T"]) < 0.3:
                ax.plot(T, n, "o", color=ACCENT, ms=7)
            elif v and v[0]["verdict"] == "excluded":
                ax.plot(T, n, "x", color=BAD, ms=7, mew=2)
            elif v:
                ax.plot(T, n, "o", mfc="none", mec=MUTED, ms=5)
            else:
                ax.plot(T, n, ".", color="#b9c2cc", ms=3)
        alive = n in (4, 8)
        ax.text(3465, n, f"{P:7.2f} d", fontsize=7.5, va="center", color=GOOD if alive else MUTED,
                fontweight="bold" if alive else "normal")
    ax.axvline(TR[1]["T"], color=ACCENT, lw=0.8, ls=":")
    ax.set_ylim(13.5, 0.4)
    ax.set_xlim(1650, 3600)
    ax.set_yticks(range(1, 13))
    ax.set_ylabel("n  (P = dT / n)", fontsize=8.5)
    ax.set_xlabel("BTJD", fontsize=8.5)
    ax.plot([], [], "o", color=INK, label="known transits (S14, S41)")
    ax.plot([], [], "o", color=ACCENT, label="S21 transit (ingress + 9 h)")
    ax.plot([], [], "x", color=BAD, mew=2, label="flat TESS data: excluded")
    ax.plot([], [], ".", color="#b9c2cc", label="no data")
    ax.legend(fontsize=7, loc="lower center", bbox_to_anchor=(0.45, 1.0), frameon=False, ncol=4)
    style(ax)
    fig.tight_layout()
    fig.savefig(HERE / "fig_aliases.png", dpi=170)
    plt.close(fig)


def amp_map(sec, T0, lo, hi, exptime):
    f = glob.glob(str(ROOT / f"data/raw/tesscut/159614322/**/*s{sec:04d}*.fits"), recursive=True)[0]
    with fits.open(f) as h:
        tab, wcs = h[1].data, WCS(h[2].header)
    t, fl, q = tab["TIME"], tab["FLUX"].astype(float), tab["QUALITY"]
    ok = np.isfinite(t) & (q == 0) & np.all(np.isfinite(fl.reshape(len(t), -1)), axis=1) & (t > lo) & (t < hi)
    t, fl = t[ok], fl[ok]
    tx, ty = [float(v) for v in wcs.all_world2pix(RA, DEC, 0)]
    ny, nx = fl.shape[1:]
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(xx - tx, yy - ty)
    bkg = np.median(fl[:, r > 5.5], axis=1)
    F = (fl - bkg[:, None, None]).reshape(len(t), -1)
    tr = 1 - vc.model(t, T0, P8, SH["rp"], SH["b"], SH["aR"], exptime)
    A = np.column_stack([np.ones(len(t)), t - T0, (t - T0) ** 2, tr / tr.max()])
    coef, *_ = np.linalg.lstsq(A, F, rcond=None)
    resid = F - A @ coef
    covA = np.linalg.inv(A.T @ A)
    c = -coef[3].reshape(ny, nx)
    ce = (resid.std(0) * np.sqrt(covA[3, 3])).reshape(ny, nx)
    sel = r <= 3.5
    fit = least_squares(lambda p: ((di.gaussian(p, yy, xx) - c) / np.maximum(ce, 1e-3))[sel],
                        [c[sel].max(), tx, ty, 0.8, 0],
                        bounds=([0, tx - 3, ty - 3, 0.3, -np.inf], [np.inf, tx + 3, ty + 3, 3, np.inf]))
    J = fit.jac
    cov = np.linalg.inv(J.T @ J) * max(2 * fit.cost / (sel.sum() - 5), 1)
    pk = np.max(np.where(sel, c / ce, -np.inf))
    return dict(c=c, tx=tx, ty=ty, fx=fit.x[1], fy=fit.x[2], n=int(len(t)),
                off=float(np.hypot(fit.x[1] - tx, fit.x[2] - ty) * 21), off_err=float(np.sqrt(cov[1, 1] + cov[2, 2]) * 21),
                snr=float(pk), wcs=wcs)


def fig_maps():
    g = di.gaia_field(RA, DEC)
    cfg = ((14, TR[0]["T"], TR[0]["T"] - 0.95, TR[0]["T"] + 0.95, 1800),
           (21, TR[1]["T"], 1883.0, 1883.89, 1800),
           (41, TR[2]["T"], TR[2]["T"] - 0.95, TR[2]["T"] + 0.95, 600))
    fig, axs = plt.subplots(1, 3, figsize=(8.4, 3.0))
    out = {}
    for ax, (sec, T0, lo, hi, ex) in zip(axs, cfg):
        m = amp_map(sec, T0, lo, hi, ex)
        gx, gy = m["wcs"].all_world2pix(g.ra.values, g.dec.values, 0)
        ax.imshow(np.clip(m["c"], 0, None), origin="lower", cmap="magma", interpolation="nearest")
        br = (g.phot_g_mean_mag < 16).values
        ax.plot(gx[br], gy[br], "s", ms=6, mfc="none", mec="white", mew=1)
        ax.plot(m["tx"], m["ty"], "o", ms=13, mfc="none", mec="#7fd4ff", mew=2)
        ax.plot(m["fx"], m["fy"], "+", ms=10, color="lime", mew=1.5)
        ax.set_xlim(m["tx"] - 4, m["tx"] + 4)
        ax.set_ylim(m["ty"] - 4, m["ty"] + 4)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"S{sec}: source {m['off']:.1f}\" +/- {m['off_err']:.1f}\"\nfrom target (SNR {m['snr']:.1f})",
                     fontsize=8, color=INK)
        out[f"S{sec}"] = dict(offset_arcsec=m["off"], offset_err_arcsec=m["off_err"], peak_snr=m["snr"], n_frames=m["n"])
    fig.tight_layout()
    fig.subplots_adjust(top=0.8)
    fig.savefig(HERE / "fig_maps.png", dpi=170)
    plt.close(fig)
    (HERE / "maps.json").write_text(json.dumps(out, indent=1))
    print(out)


if __name__ == "__main__":
    fig_transits()
    fig_aliases()
    fig_maps()
    print("done")
