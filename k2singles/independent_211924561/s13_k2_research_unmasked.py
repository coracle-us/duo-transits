"""Re-run the single-event search on K2 C5/C16/C18 WITHOUT EVEREST's outlier bits.

Lesson from s12: EVEREST sets outlier bits (2^22, 2^24) on in-transit cadences, so QUALITY==0 hides
transits. Here: EVEREST keeping everything except thruster firings (2^20); and K2-pipeline SAP with
only hard spacecraft flags removed (attitude tweak, safe mode, coarse point, earth point, desat,
manual exclude, thruster firing), detrended against roll with a running median of 0.75 d.
"""
import os, json
import numpy as np
from astropy.io import fits
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lcio

HARD = 1 | 2 | 4 | 8 | 32 | 128 | (1 << 20)


def everest(c):
    d = fits.getdata(lcio._find(f"hlsp_everest_k2_llc_211924561-c{c:02d}_kepler_v2.0_lc.fits"), 1)
    q = (d["QUALITY"] & (1 << 20)) == 0
    return lcio._norm(d["TIME"] - 2167.0, np.asarray(d["FLUX"], float), None, q)


def sap(c):
    d = fits.getdata(lcio._find(f"ktwo211924561-c{c:02d}_llc.fits"), 1)
    q = (d["SAP_QUALITY"] & HARD) == 0
    return lcio._norm(d["TIME"] - 2167.0, np.asarray(d["SAP_FLUX"], float), None, q)


def box_search(t, r, durations=(0.08, 0.10, 0.12, 0.14), step=0.005):
    sig_pt = 1.4826 * np.median(np.abs(r - np.median(r)))
    cad = np.median(np.diff(t))
    out = []
    for D in durations:
        nb = max(int(round(D / cad)), 1)
        nbins = len(r) // nb
        means = r[: nbins * nb].reshape(nbins, nb).mean(1)
        beta = max(1.0, np.std(means) / (sig_pt / np.sqrt(nb)))
        for tc in np.arange(t.min(), t.max(), step):
            inb = np.abs(t - tc) < D / 2
            n = inb.sum()
            if n < max(3, 0.6 * D / cad):
                continue
            dep = -np.mean(r[inb])
            out.append((tc, D, dep, dep / (beta * sig_pt / np.sqrt(n)), n, beta))
    return np.array(out), sig_pt


res = {}
fig, axes = plt.subplots(6, 1, figsize=(16, 13))
i = 0
for c in (5, 16, 18):
    for nm, fn in (("EVEREST_noOutlierMask", everest), ("SAP_hardflags", sap)):
        d = fn(c)
        t, f = d["t"], d["f"]
        o = np.argsort(t); t, f = t[o], f[o]
        flat, trend = lcio.flatten(t, f, window=0.75)
        r = flat - 1
        g = np.isfinite(r); t, r = t[g], r[g]
        # SAP retains the K2 roll sawtooth: remove it with a 6-cadence running median of residuals? no -
        # instead clip nothing and rely on depth+duration; report top events
        s, sig = box_search(t, r)
        cand = []
        for k in np.argsort(-s[:, 3]):
            if any(abs(s[k, 0] - x[0]) < 0.3 for x in cand):
                continue
            cand.append(tuple(s[k]))
            if len(cand) >= 6:
                break
        key = f"C{c}_{nm}"
        res[key] = [dict(tc=round(x[0], 3), dur=x[1], depth_ppt=round(1e3 * x[2], 2), snr=round(x[3], 1)) for x in cand]
        print(key, f"rms={sig*1e6:.0f}ppm span {t.min():.2f}-{t.max():.2f}")
        for x in res[key]:
            print("    ", x)
        ax = axes[i]; i += 1
        ax.plot(t, 1e3 * r, ".", ms=1.5, color="k"); ax.set_ylim(-3.5, 1.5)
        ax.text(0.005, 0.05, key, transform=ax.transAxes)
        for x in cand:
            ax.axvline(x[0], color="r" if x[3] > 7 else "orange", alpha=0.4)
plt.tight_layout(); plt.savefig(os.path.join(lcio.HERE, "fig_k2_research_unmasked.png"), dpi=80)
json.dump(res, open(os.path.join(lcio.HERE, "search_k2_unmasked.json"), "w"), indent=1)
