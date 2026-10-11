"""
Figures for the TIC 147035057 (CTOI 147035057.01) write-up.
  fig_eclipses.png  - the three TESS eclipses (S7, S33, S108) with the common model
  fig_s108.png      - Sector 108 (own TESScut photometry) with every alias that predicted an
                      eclipse in it; only P = 37.315 d shows one
  fig_diffimg.png   - difference images at the S33 and S108 eclipses
Run:  python singles/writeups/tic147035057/make_figures.py
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

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
import alias_pipeline as ap  # noqa: E402
import diff_image as di  # noqa: E402
import vet_candidate as vc  # noqa: E402

OUT = Path(__file__).resolve().parent
VET = ROOT / "singles" / "results" / "vetting"
TIC, RA, DEC = 147035057, 104.274396, -16.478476
INK, MUTED, ACCENT = "#1d2b3a", "#5d6b7c", "#d9682b"
D = 4.53 / 24


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=8, colors=MUTED)


def eph():
    return json.load(open(VET / "147035057_ephemeris.json"))


def fig_eclipses():
    e = eph()
    r = json.load(open(VET / "147035057_n1.json"))
    rp, b, aR, Psh = r["fit"]["rp"], r["fit"]["b"], r["fit"]["aR"], r["fit"]["P"]
    prods = vc.all_products(TIC)
    fig, axs = plt.subplots(1, 3, figsize=(9.5, 3.0), sharey=True)
    for ax, Tm, sec, pref in zip(axs, e["epochs"], (7, 33, 108), ("TESS-SPOC", "SPOC_120", "OWN")):
        p = next(q for q in prods if q["sector"] == sec and q["product"] == pref)
        w = vc.window(p, Tm, D)
        t, f = w[0], w[1]
        x = (t - Tm) * 24
        sparse = p["exptime"] >= 600  # 10/30-min FFI cadence: show every point prominently
        ax.plot(x, (f - 1) * 1e3, "o" if sparse else ".", ms=3.5 if sparse else 2,
                color=INK if sparse else "#9fb3c8")
        bins = np.arange(-14, 14.01, 0.5)
        idx = np.digitize(x, bins)
        bx = [x[idx == i].mean() for i in np.unique(idx) if (idx == i).sum() > 1]
        by = [((f - 1) * 1e3)[idx == i].mean() for i in np.unique(idx) if (idx == i).sum() > 1]
        ax.plot(bx, by, "o", ms=3.5, color=INK)
        tt = np.linspace(-14, 14, 500)
        ax.plot(tt, (vc.model(Tm + tt / 24, Tm, Psh, rp, b, aR, 120) - 1) * 1e3, "-", color=ACCENT, lw=1.6)
        lab = {"TESS-SPOC": "TESS-SPOC 30-min", "SPOC_120": "SPOC 2-min", "OWN": "own TESScut 200-s"}[pref]
        yr = {7: "Jan 2019", 33: "Jan 2021", 108: "Sep 2026"}[sec]
        ax.set_title(f"Sector {sec} ({yr})\nBTJD {Tm:.4f}  [{lab}]", fontsize=8.5, color=INK)
        ax.set_xlim(-14, 14)
        ax.set_ylim(-5.5, 2)
        ax.set_xlabel("hours from mid-eclipse", fontsize=8.5)
        style(ax)
    axs[0].set_ylabel("relative flux (ppt)", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(OUT / "fig_eclipses.png", dpi=170)
    plt.close(fig)


def fig_s108():
    e = eph()
    own = ap.load_lightcurves(TIC, products=["OWN"])
    o = own[own.sector == 108]
    Tc1, dT = 1515.0652, 2224.0514 - 1515.0652
    fig, ax = plt.subplots(figsize=(9.5, 2.6))
    ax.plot(o.time, (o.flux - 1) * 1e3, ".", ms=1, color="#9fb3c8")
    tb = np.arange(o.time.min(), o.time.max(), 1 / 24)
    idx = np.digitize(o.time, tb)
    ax.plot([o.time[idx == i].mean() for i in np.unique(idx)],
            [((o.flux - 1) * 1e3)[idx == i].mean() for i in np.unique(idx)], "-", color=INK, lw=0.8)
    for n, verdict in ((12, "inconclusive"), (13, "inconclusive"), (14, "flat: excluded"), (19, "eclipse")):
        P = dT / n
        k = round((4305 - Tc1) / P)
        for kk in (k - 1, k, k + 1):
            T = Tc1 + kk * P
            if o.time.min() < T < o.time.max():
                col = ACCENT if n == 19 else MUTED
                ax.axvline(T, color=col, lw=1.2 if n == 19 else 0.8, ls="-" if n == 19 else "--")
                ax.text(T + 0.2, -4.6, f"P = {P:.2f} d\n{verdict}", fontsize=7.5, color=col, va="top")
    ax.set_ylim(-7, 4)
    ax.set_xlabel("BTJD", fontsize=8.5)
    ax.set_ylabel("relative flux (ppt)", fontsize=8.5)
    style(ax)
    fig.tight_layout()
    fig.savefig(OUT / "fig_s108.png", dpi=170)
    plt.close(fig)


def diff_panel(ax, sector, Tm, title):
    path = di.get_cutout(TIC, RA, DEC, sector)
    with fits.open(path) as h:
        tab, wcs = h[1].data, WCS(h[2].header)
    t, fl, q = tab["TIME"], tab["FLUX"].astype(float), tab["QUALITY"]
    good = np.isfinite(t) & (q == 0) & np.all(np.isfinite(fl.reshape(len(t), -1)), axis=1)
    t, fl = t[good], fl[good]
    dt = t - Tm
    inn = np.abs(dt) < 0.35 * D
    out = (np.abs(dt) > D / 2 + 0.05) & (np.abs(dt) < D / 2 + 0.6)
    A = np.column_stack([np.ones(out.sum()), dt[out]])
    coef, *_ = np.linalg.lstsq(A, fl[out].reshape(out.sum(), -1), rcond=None)
    pred = np.column_stack([np.ones(inn.sum()), dt[inn]]) @ coef
    diff = (pred - fl[inn].reshape(inn.sum(), -1)).mean(0).reshape(fl.shape[1:])
    tx, ty = wcs.all_world2pix(RA, DEC, 0)
    g = di.gaia_field(RA, DEC)
    gx, gy = wcs.all_world2pix(g.ra.values, g.dec.values, 0)
    ax.imshow(np.clip(diff, 0, None), origin="lower", cmap="magma", interpolation="nearest")
    br = (g.phot_g_mean_mag < 15).values
    ax.plot(gx[br], gy[br], "s", ms=6, mfc="none", mec="white", mew=1)
    ax.plot(tx, ty, "o", ms=12, mfc="none", mec="#7fd4ff", mew=2)
    ax.set_xlim(tx - 4, tx + 4)
    ax.set_ylim(ty - 4, ty + 4)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=8.5, color=INK)


def fig_diffimg():
    e = eph()
    fig, axs = plt.subplots(1, 2, figsize=(6.0, 3.1))
    diff_panel(axs[0], 33, e["epochs"][1], "S33 difference image\n(source 2.9\" +/- 1.0\" from target)")
    diff_panel(axs[1], 108, e["epochs"][2], "S108 difference image\n(source 2.2\" +/- 0.8\" from target)")
    fig.tight_layout()
    fig.savefig(OUT / "fig_diffimg.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    fig_eclipses()
    fig_s108()
    fig_diffimg()
    print("done")
