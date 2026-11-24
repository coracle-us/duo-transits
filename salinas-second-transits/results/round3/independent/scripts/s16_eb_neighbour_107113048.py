"""TIC 107113048: the G=12.74 Gaia star 97" away (DR3 3463196954429755904) is an eclipsing binary whose
eclipse in S90 overlaps transit E29. Checks that E29 is a real transit of the target:
  * 3x3-pixel light curves of the target and the binary in every 21x21 / S101 cutout (plot);
  * the binary's eclipse centre, depth and FWHM near each eclipse seen, and the change in the target's
    3x3 flux at the binary's mid-eclipse (leakage);
  * the E29 time predicted from E0, E15 and E36 alone, vs the observed E29 time and the binary's eclipse.
Outputs: out/eb_neighbour_107113048.txt, out/eb_neighbour_107113048.png"""
import sys, os, glob, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import TARGETS, OUT, DATA

TIC, EB = 107113048, 3463196954429755904
tg = TARGETS[TIC]
nb = [n for n in json.load(open(os.path.join(OUT, f"gaia_{TIC}.json")))["neighbours"] if n["source_id"] == EB][0]
eph = json.load(open(os.path.join(OUT, f"ephem_{TIC}.json"))); P, T0 = eph["P"], eph["T0_E0"]
files = sorted(glob.glob(os.path.join(DATA, "tesscut_diff", "**", f"*{tg['ra']:.6f}*21x21*.fits"), recursive=True) +
               glob.glob(os.path.join(DATA, "tesscut", "**", f"*s0101*{tg['ra']:.6f}*.fits"), recursive=True))
L = [f"Eclipsing-binary neighbour Gaia DR3 {EB}: G={nb['G']:.2f}, {nb['sep']:.1f}\" from TIC {TIC}"]


def lcs(fn):
    h = fits.open(fn); d = h[1].data; w = WCS(h[2].header)
    ok = (d["QUALITY"] == 0) & np.isfinite(d["TIME"]); t = d["TIME"][ok]; c = d["FLUX"][ok]
    bkg = np.nanmedian(c.reshape(len(t), -1), axis=1)
    r = {}
    for lab, ra, de in [("target", tg["ra"], tg["dec"]), ("binary", nb["ra"], nb["dec"])]:
        x, y = w.all_world2pix(ra, de, 0); ix, iy = int(round(float(x))), int(round(float(y)))
        f = np.nansum(c[:, iy - 1:iy + 2, ix - 1:ix + 2] - bkg[:, None, None], axis=(1, 2))
        r[lab] = f / np.nanmedian(f)
    return int(h[0].header["SECTOR"]), t, r


fig, axs = plt.subplots(len(files), 1, figsize=(10, 2.8 * len(files)))
for ax, fn in zip(axs, files):
    sec, t, r = lcs(fn)
    k = np.round((np.median(t) - T0) / P); tc = T0 + k * P
    ax.plot(t, r["target"], ".", ms=1, label="target 3x3")
    ax.plot(t, r["binary"] - 0.05, ".", ms=1, label=f"G{nb['G']:.2f} @{nb['sep']:.0f}\" 3x3 (-0.05)")
    ax.axvline(tc, color="r", lw=0.5); ax.set_ylim(0.8, 1.05); ax.legend(fontsize=7)
    ax.set_title(f"S{sec}: red line = target transit E{int(k)} ({tc:.2f})", fontsize=9)
    # binary eclipses in this sector: runs of points > 5% below median
    low = r["binary"] < 0.95
    if not low.any():
        L.append(f"S{sec}: no binary eclipse in the cutout data")
        continue
    groups = np.split(np.where(low)[0], np.where(np.diff(t[low]) > 0.5)[0] + 1)
    for gidx in groups:
        e0 = np.median(t[gidx])
        oot = (np.abs(t - e0) > 0.6) & (np.abs(t - e0) < 2.0)
        win = np.abs(t - e0) < 0.6
        if oot.sum() < 20 or win.sum() < 20:
            L.append(f"S{sec}: binary eclipse near {e0:.2f} at a data edge - not measured")
            continue
        n = {lab: r[lab] / np.polyval(np.polyfit(t[oot] - e0, r[lab][oot], 2), t - e0) for lab in r}
        defi = np.clip(1 - n["binary"][win], 0, None)
        ce = np.sum(t[win] * defi) / np.sum(defi)
        mid = np.abs(t - ce) < 0.03
        dep = 1 - np.median(n["binary"][mid]); fw = t[win][defi > dep / 2]
        if (fw.max() - fw.min()) * 24 > 15:
            L.append(f"S{sec}: broad dip near {ce:.2f} (FWHM {(fw.max()-fw.min())*24:.0f} h) seen in both stars' pixels "
                     f"(target {(np.median(n['target'][mid])-1)*100:+.1f}%) - scattered-light systematic, not an eclipse")
            continue
        bins = n["target"][oot][:oot.sum() // 30 * 30].reshape(-1, 30).mean(1)
        L.append(f"S{sec}: binary eclipse centre {ce:.4f} BTJD, depth {dep*100:.1f}% in its 3x3, FWHM "
                 f"{(fw.max()-fw.min())*24:.1f} h; target 3x3 flux at binary mid-eclipse {(np.median(n['target'][mid])-1)*1e6:+.0f} ppm "
                 f"(1-h scatter {np.std(bins)*1e6:.0f} ppm); nearest target transit E{int(k)} is {(ce-tc)*24:+.2f} h away")
fig.tight_layout(); fig.savefig(os.path.join(OUT, f"eb_neighbour_{TIC}.png"), dpi=70)
# E29 predicted from the other three transits
import csv
rows = [r for r in csv.DictReader(open(os.path.join(OUT, f"times_{TIC}.csv"))) if r["tc"] and r["timing_source"] == "1"]
E = np.array([int(r["epoch"]) for r in rows]); T = np.array([float(r["tc"]) for r in rows]); S = np.array([float(r["tc_err"]) for r in rows])
o = E != 29
A = np.vstack([np.ones(o.sum()), E[o]]).T; W = np.diag(1 / S[o] ** 2); C = np.linalg.inv(A.T @ W @ A); p = C @ A.T @ W @ T[o]
v = np.array([1, 29]); pred, pe = v @ p, np.sqrt(v @ C @ v)
obs, oe = T[E == 29][0], S[E == 29][0]
L.append(f"E29 predicted from E{', E'.join(map(str, E[o]))}: {pred:.5f} ± {pe*1440:.1f} min (P = {p[1]:.6f} ± {np.sqrt(C[1,1]):.6f} d); "
         f"observed {obs:.5f} ± {oe*1440:.1f} min -> O-C {(obs-pred)*1440:+.1f} min ({(obs-pred)/np.hypot(pe, oe):+.1f} sigma)")
txt = "\n".join(L)
open(os.path.join(OUT, f"eb_neighbour_{TIC}.txt"), "w").write(txt + "\n"); print(txt)
