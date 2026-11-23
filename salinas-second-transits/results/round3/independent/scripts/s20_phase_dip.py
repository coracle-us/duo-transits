"""Diagnostics for a dip found by s17 at a fixed orbital phase (not the transit).

Uses the QLP light curves (all sectors), which carry three apertures (DET_FLUX_SML, DET_FLUX, DET_FLUX_LAG;
KSPSAP_* in older files),
the raw SAP_FLUX, the background SAP_BKG and the centroids SAP_X/SAP_Y. For the dip and, as calibration,
for the transit itself:
 - depth in each aperture, fitted with the same trapezoid + quadratic-per-orbit model as s17, errors scaled
   by the empirical noise of that quantity over all phases (a dip from a neighbour gets deeper in the
   larger aperture; one on the target keeps the same depth, as the transit does);
 - the same fit on the background and on the centroids (in units of their scatter);
 - for each orbit: depth alone, and the distance from the dip to the nearest data gap or sector edge;
 - depth at the half-period phases (dip phase +- 0.5), to test for a signal at P/2.
Usage: s20_phase_dip.py TIC PHASE DURATION_HOURS. Outputs out/phasedip_<tic>.txt and .png."""
import sys, os, glob, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from astropy.io import fits
from common import OUT, DATA
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

tic, PH, DH = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
D = DH / 24
eph = json.load(open(os.path.join(OUT, f"ephem_{tic}.json")))
sh = json.load(open(os.path.join(OUT, f"shape_{tic}.json")))
P, T0, k = eph["P"], eph["T0_E0"], sh["rp"]
Dtr = P / np.pi * np.arcsin(np.sqrt((1 + k) ** 2 - sh["b"] ** 2) / sh["a"])
L = []
pr = L.append
BIN = 10 / 1440
COLS = ["SAP_FLUX", "DET_FLUX_SML", "DET_FLUX", "DET_FLUX_LAG", "SAP_BKG", "SAP_X", "SAP_Y"]

files = sorted(glob.glob(os.path.join(DATA, "lk", "**", f"*qlp*{tic:016d}*_llc.fits"), recursive=True))
segs = []
for fn in files:
    with fits.open(fn) as h:
        d = h[1].data
        sec = int(h[0].header["SECTOR"])
        t = np.array(d["TIME"], float)
        # older QLP files (S10, S37) name the detrended fluxes KSPSAP_*
        src = {c: (c if c in d.names else c.replace("DET_", "KSPSAP_")) for c in COLS}
        ok = np.isfinite(t) & (d["QUALITY"] == 0)
        for c in COLS:
            ok &= np.isfinite(np.array(d[src[c]], float))
        t = t[ok]
        cols = {c: np.array(d[src[c]], float)[ok] for c in COLS}
    i = np.floor((t - t[0]) / BIN).astype(int)
    n = np.bincount(i)
    g = n > 0
    tb = np.bincount(i, t)[g] / n[g]
    cb = {}
    for c in COLS:
        v = np.bincount(i, cols[c])[g] / n[g]
        cb[c] = v / np.median(v) if "FLUX" in c else v - np.median(v)
    segs.append(dict(sector=sec, t=tb, **cb))
pr(f"TIC {tic}: QLP sectors {[s['sector'] for s in segs]}; dip at phase {PH:.4f}, {DH:.2f} h; "
   f"transit {Dtr*24:.2f} h")
T = np.concatenate([s["t"] for s in segs])
S = np.concatenate([np.full(len(s["t"]), s["sector"]) for s in segs])
V = {c: np.concatenate([s[c] for s in segs]) for c in COLS}


def trap(dt, Dd):
    return np.clip((Dd / 2 - np.abs(dt)) / (Dd * k / (1 + k)), 0, 1)


def fit(y, phi, Dd, sel=None, nmin=3):
    """Depth (positive = dip) with a quadratic baseline per orbit over +-2Dd; white error from residuals."""
    x = (T - T0) / P - phi
    cyc = np.round(x)
    dt = (x - cyc) * P
    w = np.abs(dt) < 2 * Dd
    if sel is not None:
        w &= sel
    FR, SR = [], []
    for c in np.unique(cyc[w] * 1000 + S[w]):
        m = w & (cyc * 1000 + S == c)
        s = trap(dt[m], Dd)
        if (s > 0.5).sum() < 5 or (dt[m] < -Dd / 2).sum() < 5 or (dt[m] > Dd / 2).sum() < 5:
            continue
        u = dt[m] / Dd
        X = np.vstack([np.ones_like(u), u, u**2]).T
        Q, _ = np.linalg.qr(X)
        FR.append(y[m] - Q @ (Q.T @ y[m])); SR.append(-(s - Q @ (Q.T @ s)))
    if len(FR) < nmin:
        return None
    fr, sr = np.concatenate(FR), np.concatenate(SR)
    d = (sr @ fr) / (sr @ sr)
    sig = np.std(fr - d * sr)
    return d, sig / np.sqrt(sr @ sr), len(FR)


def noise_scale(y, Dd, avoid):
    z = []
    for phi in np.arange(0.02, 0.98, Dd / 4 / P):
        if min(abs(phi - a) for a in avoid) * P < Dd:
            continue
        r = fit(y, phi, Dd)
        if r:
            z.append(r[0] / r[1])
    z = np.array(z)
    return max(1.0, 1.4826 * np.median(np.abs(z - np.median(z)))), z


pr("depth by quantity (fluxes in ppm; background in e/s; centroids in milli-pixels), error x noise scale:")
pr(f"{'quantity':14s} {'transit (phase 0)':>26s}   {'dip':>26s}  noise scale")
res = {}
for c in COLS:
    ns, z = noise_scale(V[c], D, [0.0, 1.0, PH])
    rt, rd = fit(V[c], 0.0, Dtr), fit(V[c], PH, D)
    sc = 1e6 if "FLUX" in c else (1e3 if c in ("SAP_X", "SAP_Y") else 1.0)
    res[c] = (rt, rd, ns)
    pr(f"{c:14s} {rt[0]*sc:+10.1f} +- {rt[1]*sc*ns:8.1f} ({rt[0]/rt[1]/ns:+6.1f})   "
       f"{rd[0]*sc:+10.1f} +- {rd[1]*sc*ns:8.1f} ({rd[0]/rd[1]/ns:+6.1f})  {ns:.2f}")
for c in ("DET_FLUX_SML", "DET_FLUX_LAG"):
    rt, rd, _ = res[c]
    rt0, rd0, _ = res["DET_FLUX"]
    pr(f"ratio {c}/DET_FLUX: transit {rt[0]/rt0[0]:.3f}, dip {rd[0]/rd0[0]:.2f} +- "
       f"{abs(rd[0]/rd0[0])*np.hypot(rd[1]/rd[0], rd0[1]/rd0[0]):.2f} (errors correlated; indicative)")

# half-period phases
for ph2 in [(PH + 0.5) % 1, (PH + 1 / 3) % 1, (PH + 2 / 3) % 1]:
    r = fit(V["SAP_FLUX"], ph2, D)
    ns = res["SAP_FLUX"][2]
    pr(f"SAP_FLUX at phase {ph2:.4f} (P/2, P/3 tests): {r[0]*1e6:+.0f} +- {r[1]*1e6*ns:.0f} ppm")

# per orbit: depth and gap distance
pr("each orbit (SAP_FLUX), with the distance from the dip centre to the nearest gap > 0.3 d or sector edge:")
cyc = np.round((T - T0) / P - PH)
ns = res["SAP_FLUX"][2]
for c in np.unique(cyc):
    tc = T0 + (c + PH) * P
    r = fit(V["SAP_FLUX"], PH, D, sel=(cyc == c), nmin=1)
    s_ = int(np.median(S[cyc == c]))
    ts = T[S == s_]
    gaps = np.concatenate([[ts[0]], ts[:-1][np.diff(ts) > 0.3], ts[1:][np.diff(ts) > 0.3], [ts[-1]]])
    dg = np.min(np.abs(gaps - tc))
    txt = f"{r[0]*1e6:+6.0f} +- {r[1]*1e6*ns:4.0f} ppm" if r else "not enough data"
    pr(f"  E{int(c):4d} S{s_:3d} BTJD {tc:.2f}: {txt}; nearest gap/edge {dg:.2f} d")

# ---------- plot ----------
fig, ax = plt.subplots(2, 2, figsize=(12, 8))


def fold(y, phi, Dd):
    x = (T - T0) / P - phi
    cyc = np.round(x)
    dt = (x - cyc) * P
    w = np.abs(dt) < 2 * Dd
    xs, ys = [], []
    for c in np.unique(cyc[w] * 1000 + S[w]):
        m = w & (cyc * 1000 + S == c)
        oot = m & (np.abs(dt) > Dd / 2)
        if oot.sum() < 10:
            continue
        p = np.polyfit(dt[oot], y[oot], 2)
        xs.append(dt[m] * 24); ys.append(y[m] - np.polyval(p, dt[m]))
    xs, ys = np.concatenate(xs), np.concatenate(ys)
    bins = np.arange(-2 * Dd * 24, 2 * Dd * 24 + 1e-6, Dd * 24 / 6)
    i = np.digitize(xs, bins)
    B = np.array([(xs[i == j].mean(), ys[i == j].mean(), ys[i == j].std() / np.sqrt((i == j).sum()))
                  for j in range(1, len(bins)) if (i == j).sum() > 5])
    return B


for c, col in zip(["DET_FLUX_SML", "DET_FLUX", "DET_FLUX_LAG"], ["C0", "C1", "C2"]):
    for j, (phi, Dd) in enumerate([(0.0, Dtr), (PH, D)]):
        B = fold(V[c], phi, Dd)
        ax[0, j].errorbar(B[:, 0], B[:, 1] * 1e6, B[:, 2] * 1e6, fmt="o-", ms=3, lw=0.7, color=col, label=c)
ax[0, 0].set_title("transit (calibration)"); ax[0, 1].set_title(f"dip at phase {PH:.3f}")
for a in ax[0]:
    a.set_xlabel("hours"); a.set_ylabel("flux - 1 (ppm)"); a.legend(fontsize=7); a.axhline(0, color="gray", lw=0.5)
B = fold(V["SAP_BKG"], PH, D)
ax[1, 0].errorbar(B[:, 0], B[:, 1], B[:, 2], fmt="o-", ms=3)
ax[1, 0].set_title("background at the dip"); ax[1, 0].set_xlabel("hours"); ax[1, 0].set_ylabel("SAP_BKG - median")
for c, col in zip(["SAP_X", "SAP_Y"], ["C3", "C4"]):
    B = fold(V[c], PH, D)
    ax[1, 1].errorbar(B[:, 0], B[:, 1] * 1e3, B[:, 2] * 1e3, fmt="o-", ms=3, color=col, label=c)
ax[1, 1].set_title("centroids at the dip"); ax[1, 1].set_xlabel("hours"); ax[1, 1].set_ylabel("milli-pixel")
ax[1, 1].legend(fontsize=7)
fig.suptitle(f"TIC {tic}: QLP diagnostics for the phase-{PH:.3f} dip")
fig.tight_layout()
fig.savefig(os.path.join(OUT, f"phasedip_{tic}.png"), dpi=110)
txt = "\n".join(L)
open(os.path.join(OUT, f"phasedip_{tic}.txt"), "w").write(txt + "\n")
print(txt)
