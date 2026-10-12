"""EPIC 211390677 / TIC 366568535: own TESScut aperture photometry for S44/S45/S46 (QLP-only sectors)
to test whether the QLP S45 dip at BTJD ~2547.8-2548.0 is real and on-target. Also pixel-level: light
curves of the target aperture vs an annulus, and a difference image (in-dip minus out-of-dip).
Everything is written into review_fable2/ (nothing in data/raw is touched)."""
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astroquery.mast import Tesscut
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
TIC, RA, DEC = 366568535, 127.703812, 11.49457
T_DIP, D = 2547.85, 39.5 / 24
cdir = HERE / "tesscut_366568535"
cdir.mkdir(exist_ok=True)
coord = SkyCoord(RA, DEC, unit="deg")
out = open(HERE / "tesscut_s45_366568535.txt", "w")


def say(*a):
    print(*a)
    print(*a, file=out)


for sec in (44, 45, 46, 72):
    hits = list(cdir.glob(f"*s00{sec}*.fits"))
    if not hits:
        tabs = Tesscut.download_cutouts(coordinates=coord, size=15, sector=sec, path=str(cdir))
        hits = list(cdir.glob(f"*s00{sec}*.fits"))
    if not hits:
        say(f"S{sec}: no cutout")
        continue
    with fits.open(hits[0]) as h:
        tab, wcs = h[1].data, WCS(h[2].header)
    t, flux, q = tab["TIME"], tab["FLUX"].astype(float), tab["QUALITY"]
    good = np.isfinite(t) & (q == 0) & np.isfinite(flux).all(axis=(1, 2))
    t, flux = t[good], flux[good]
    tx, ty = (float(x) for x in wcs.all_world2pix(RA, DEC, 0))
    yy, xx = np.mgrid[0:flux.shape[1], 0:flux.shape[2]]
    rr = np.hypot(xx - tx, yy - ty)
    med = np.nanmedian(flux, axis=0)
    bkg = np.nanmedian(flux[:, med <= np.nanpercentile(med, 40)], axis=1)
    say(f"S{sec}: {len(t)} good cadences, BTJD {t.min():.2f}-{t.max():.2f}; target at pixel ({tx:.2f},{ty:.2f})")
    for rad in (1.0, 1.6, 2.5):
        ap = rr <= rad
        f = np.nansum(flux[:, ap] - bkg[:, None], axis=1)
        f /= np.nanmedian(f)
        # box test at T_DIP with linear baseline from +-(D/2+0.1 .. D/2+1.5)
        inn = np.abs(t - T_DIP) < 0.35 * D
        oot = (np.abs(t - T_DIP) > 0.5 * D + 0.1) & (np.abs(t - T_DIP) < 0.5 * D + 1.5)
        if inn.sum() > 5 and oot.sum() > 20:
            c = np.polyfit(t[oot] - T_DIP, f[oot], 1)
            r = f / np.polyval(c, t - T_DIP) - 1
            sig = 1.4826 * np.median(np.abs(r[oot] - np.median(r[oot])))
            d = -np.mean(r[inn])
            say(f"   aperture r<={rad} px ({ap.sum()} px): depth at {T_DIP} = {d*1e3:.2f} ppt, oot rms {sig*1e3:.2f} ppt/cad, "
                f"SNR(white) = {d/(sig/np.sqrt(inn.sum())):.1f}")
        else:
            say(f"   aperture r<={rad} px: dip window not covered in S{sec} (in={inn.sum()}, oot={oot.sum()})")
        if sec == 45 and rad == 1.6:
            fig, ax = plt.subplots(2, 1, figsize=(10, 6))
            ax[0].plot(t, f, ".", ms=2); ax[0].set_title(f"TIC {TIC} S45 own TESScut r<=1.6 px"); ax[0].axvline(T_DIP, color="r", lw=0.5)
            m = np.abs(t - T_DIP) < 3
            ax[1].plot(t[m], f[m], ".", ms=3); ax[1].axvspan(T_DIP - D / 2, T_DIP + D / 2, alpha=0.2, color="r")
            fig.tight_layout(); fig.savefig(HERE / "tesscut_s45_366568535_lc.png", dpi=100)
    if sec == 72:
        T72 = 151.7809 + 17 * 184.3087  # = 3285.03; identical for P = 92.154 (k=34) and 61.436 (k=51)
        for rad in (1.0, 1.6, 2.5):
            ap = rr <= rad
            f = np.nansum(flux[:, ap] - bkg[:, None], axis=1); f /= np.nanmedian(f)
            base = (t > T72 - 0.82 - 2.5) & (t < T72 - 0.82 - 0.1)
            c = np.polyfit(t[base] - T72, f[base], 1)
            r = f / np.polyval(c, t - T72) - 1
            sig = 1.4826 * np.median(np.abs(r[base] - np.median(r[base])))
            say(f"   S72 (ends {t.max():.2f}): predicted eclipse centre {T72:.2f}, ingress from {T72-0.82:.2f}. r<={rad}px: baseline rms {sig*1e3:.2f} ppt;"
                f" mean residual in [T1, end] = {np.mean(r[t > T72-0.82])*1e3:+.2f} ppt (n={(t > T72-0.82).sum()}),"
                f" in [T72-0.4, end] = {np.mean(r[t > T72-0.4])*1e3:+.2f} ppt (n={(t > T72-0.4).sum()}); last 0.5 d before T1: {np.mean(r[(t > T72-1.32)&(t < T72-0.82)])*1e3:+.2f} ppt")
        fig, ax = plt.subplots(1, 1, figsize=(10, 4))
        ap = rr <= 1.6; f = np.nansum(flux[:, ap] - bkg[:, None], axis=1); f /= np.nanmedian(f)
        m = t > T72 - 6
        ax.plot(t[m], f[m], ".", ms=3); ax.axvline(T72, color="r"); ax.axvspan(T72 - 0.82, T72 + 0.82, alpha=0.2, color="r")
        ax.set_title("TIC 366568535 S72 own TESScut r<=1.6 px; red = predicted eclipse for P=184.31/92.15/61.44")
        fig.tight_layout(); fig.savefig(HERE / "tesscut_s72_366568535_lc.png", dpi=100)
    if sec == 45:
        # difference image: mean out-of-dip (both sides) minus mean in-dip
        inn = np.abs(t - T_DIP) < 0.35 * D
        oot = (np.abs(t - T_DIP) > 0.5 * D + 0.1) & (np.abs(t - T_DIP) < 0.5 * D + 1.2)
        diff = np.nanmean(flux[oot], axis=0) - np.nanmean(flux[inn], axis=0)
        noise = np.nanstd(flux[oot], axis=0) * np.sqrt(1 / oot.sum() + 1 / inn.sum())
        snr = diff / noise
        iy, ix = np.unravel_index(np.nanargmax(snr), snr.shape)
        # centroid of positive diff within 2 px of peak
        w = np.clip(diff, 0, None) * (np.hypot(xx - ix, yy - iy) <= 2)
        cx, cy = np.sum(w * xx) / w.sum(), np.sum(w * yy) / w.sum()
        say(f"   S45 difference image: peak SNR {np.nanmax(snr):.1f} at pixel ({ix},{iy}); flux-weighted centroid ({cx:.2f},{cy:.2f});"
            f" offset from target = {np.hypot(cx-tx, cy-ty):.2f} px = {np.hypot(cx-tx, cy-ty)*21:.0f} arcsec")
        # also: per-pixel depth (ppt) for the 3x3 around target vs around peak
        fig, axs = plt.subplots(1, 3, figsize=(12, 4))
        axs[0].imshow(np.log10(np.clip(med, 1, None)), origin="lower"); axs[0].plot(tx, ty, "r+"); axs[0].set_title("median frame (log)")
        axs[1].imshow(diff, origin="lower"); axs[1].plot(tx, ty, "r+"); axs[1].set_title("oot - in-dip")
        axs[2].imshow(snr, origin="lower"); axs[2].plot(tx, ty, "r+"); axs[2].plot(cx, cy, "wx"); axs[2].set_title("SNR")
        fig.tight_layout(); fig.savefig(HERE / "tesscut_s45_366568535_diff.png", dpi=100)
        # neighbours: pixel light curves of the 8 neighbours around target, depth each
        say("   per-pixel depth (ppt of that pixel's own flux) in 5x5 around the target:")
        for dy in range(2, -3, -1):
            row = []
            for dx in range(-2, 3):
                px, py = int(round(tx)) + dx, int(round(ty)) + dy
                fp = flux[:, py, px] - bkg
                c = np.polyfit(t[oot] - T_DIP, fp[oot], 1)
                r = fp / np.polyval(c, t - T_DIP) - 1
                row.append(f"{-np.mean(r[inn])*1e3:7.1f}")
            say("     ", " ".join(row))
# QLP comparison from the repo's cached file
qlp = pd.read_csv(HERE.parents[1] / "data" / "raw" / "lightcurves" / str(TIC) / "s045_QLP_600s.csv.gz")
say(f"\nQLP S45 columns: {list(qlp.columns)}")
q0 = qlp[(qlp.quality == 0) & np.isfinite(qlp.flux)]
m = np.abs(q0.time - T_DIP) < 3
say(f"QLP S45 points within +-3 d of dip: {m.sum()}; min flux in window {q0.flux[m].min()/np.nanmedian(q0.flux):.4f}")
out.close()
