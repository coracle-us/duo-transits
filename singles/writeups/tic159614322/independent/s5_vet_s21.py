"""Vet the dip at the predicted epoch BTJD ~1883.758 (S21, end of orbit 1).
- own TESScut photometry with several apertures and two background estimators
- fixed-shape depth fit at predicted time, and free mid-time scan (chi2 vs tc)
- comparison with the known S14 event measured the same way
- background and centroid behaviour
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from phot import photometry
from lcio import load, OUT
from model import shape

E = json.load(open(os.path.join(OUT, "epochs.json")))
T14 = E["T14_d"]
TPRED = [r for r in E["epochs"] if r["n"] == 4 and r["k"] == 1][0]
TC_PRED, SIG_PRED = TPRED["tc"], TPRED["sig_tc"]
TC1, TC2 = E["TC1"], E["TC2"]


def depth_fit(t, f, tc, cad, deg=1, win=(None, None)):
    S = shape(t, tc, cad)
    x = t - tc
    X = np.vstack([x ** i for i in range(deg + 1)] + [-S]).T
    beta, *_ = np.linalg.lstsq(X, f, rcond=None)
    r = f - X @ beta
    rms = np.std(r, ddof=X.shape[1])
    cov = np.linalg.inv(X.T @ X) * rms ** 2
    return beta[-1], np.sqrt(cov[-1, -1]), np.sum(r ** 2), rms


def scan(t, f, cad, centre, half=0.6, step=0.002, deg=1):
    grid = np.arange(centre - half, centre + half, step)
    chi = []
    A = []
    for tc in grid:
        S = shape(t, tc, cad)
        if np.sum(S > 0.5) < 2:
            chi.append(np.nan); A.append(np.nan); continue
        a, sa, c2, rms = depth_fit(t, f, tc, cad, deg)
        chi.append(c2); A.append(a)
    chi = np.array(chi); A = np.array(A)
    # null model
    x = t - centre
    X0 = np.vstack([x ** i for i in range(deg + 1)]).T
    b0, *_ = np.linalg.lstsq(X0, f, rcond=None)
    r0 = f - X0 @ b0
    rms0 = np.std(r0)
    return grid, chi, A, np.sum(r0 ** 2), rms0


def best_tc(grid, chi, rms):
    c = chi / rms ** 2
    i = np.nanargmin(c)
    ok = np.isfinite(c) & (c < c[i] + 1)
    return grid[i], grid[ok].min(), grid[ok].max(), c


results = {"tc_pred": TC_PRED, "sig_pred": SIG_PRED, "ingress_pred": TC_PRED - T14 / 2, "egress_pred": TC_PRED + T14 / 2}

# ----- products -----
prod = {}
for p in ["TESS-SPOC", "QLP"]:
    for which in ["flux", "sap"]:
        d = load(21, p, which, strict_quality=True)
        m = (d.time > TC_PRED - 1.3) & (d.time < 1884.0)
        prod[f"{p}_{which}"] = (d.time[m].values, d.f[m].values / np.median(d.f[m].values), 1800)
# QLP including flagged cadences (the last ~3 h before the gap are flagged 4096)
d = load(21, "QLP", "flux", strict_quality=False)
m = (d.time > TC_PRED - 1.3) & (d.time < 1884.0) & (d.f > 0.9)
prod["QLP_flux_allq"] = (d.time[m].values, d.f[m].values / np.median(d.f[m].values), 1800)

# ----- own photometry -----
ph21 = photometry(21)
ph14 = photometry(14)
ph41 = photometry(41)
results["phot_files"] = [ph21["file"], ph14["file"], ph41["file"]]
results["star_pix_s21"] = [float(ph21["x0"]), float(ph21["y0"])]
for key in ["ap1.0", "ap1.5", "ap2.0", "ap2.5", "ap3.0", "ap2.0_annbkg"]:
    t, f, q = ph21["time"], ph21[key], ph21["quality"]
    m = (t > TC_PRED - 1.3) & (t < 1884.0) & (q == 0) & np.isfinite(f)
    prod[f"own_{key}"] = (t[m], f[m] / np.median(f[m][t[m] < TC_PRED - T14 / 2]), 1800)

# known-event depths in own photometry (S14 30-min FFI, S41 10-min FFI) for calibration
results["known_depths_own"] = {}
for sec, ph, tck, cad in [(14, ph14, TC1, 1800), (41, ph41, TC2, 600)]:
    for key in ["ap1.0", "ap1.5", "ap2.0", "ap2.5", "ap3.0", "ap2.0_annbkg"]:
        t, f, q = ph["time"], ph[key], ph["quality"]
        m = (np.abs(t - tck) < 1.0) & (q == 0) & np.isfinite(f)
        a, sa, c2, rms = depth_fit(t[m], f[m] / np.median(f[m]), tck, cad, deg=2)
        results["known_depths_own"][f"S{sec}_{key}"] = [float(a), float(sa)]

# ----- fits -----
results["fits"] = {}
for k, (t, f, cad) in prod.items():
    a, sa, c2, rms = depth_fit(t, f, TC_PRED, cad, deg=1)
    grid, chi, A, c0, rms0 = scan(t, f, cad, TC_PRED, half=0.5, step=0.002)
    tb, tlo, thi, c = best_tc(grid, chi, rms)
    dchi_pred = (c0 - np.nanmin(chi)) / rms ** 2
    entry = dict(depth_at_pred=float(a), sig=float(sa), snr=float(a / sa), rms=float(rms), npts=int(t.size),
                 tc_best=float(tb), tc_lo=float(tlo), tc_hi=float(thi), dchi2_vs_flat=float(dchi_pred),
                 depth_at_best=float(A[np.nanargmin(chi)]))
    results["fits"][k] = entry
    print(f"{k:22s} depth@pred={a*1e3:5.2f}±{sa*1e3:4.2f} ppt (SNR {a/sa:4.1f})  tc_best={tb:.4f} [{tlo:.4f},{thi:.4f}]  "
          f"Δtc={(tb-TC_PRED)*24:+.2f} h  Δχ²(flat-best)={dchi_pred:.1f}")
print("known event depths in own photometry:")
for k, v in results["known_depths_own"].items():
    print(f"   {k:20s} {v[0]*1e3:5.2f}±{v[1]*1e3:4.2f} ppt")

# ----- plots -----
fig, axs = plt.subplots(5, 1, figsize=(10, 15), sharex=True)
ax = axs[0]
off = 0
for k in ["TESS-SPOC_flux", "TESS-SPOC_sap", "QLP_flux_allq", "QLP_flux"]:
    t, f, cad = prod[k]
    ax.plot(t - TC_PRED, f + off, ".-", lw=0.4, ms=4, label=k)
    off -= 0.005
ax.set_title(f"S21 pipeline products around predicted Tc={TC_PRED:.3f}±{SIG_PRED:.3f} (n=4/8); shaded = predicted T14")
ax.legend(fontsize=7)
ax = axs[1]
off = 0
for k in ["own_ap1.0", "own_ap1.5", "own_ap2.0", "own_ap2.5", "own_ap3.0", "own_ap2.0_annbkg"]:
    t, f, cad = prod[k]
    ax.plot(t - TC_PRED, f + off, ".-", lw=0.4, ms=4, label=k)
    off -= 0.004
ax.set_title("Own TESScut photometry, apertures r=1..3 px (offsets)")
ax.legend(fontsize=7)
t = ph21["time"]
m = (t > TC_PRED - 1.3) & (t < 1884.0)
axs[2].plot(t[m] - TC_PRED, ph21["bkg_per_pix"][m], ".", label="bkg per pix (median, far pixels)")
axs[2].plot(t[m] - TC_PRED, ph21["bkg_ann"][m], ".", label="bkg per pix (annulus, clipped mean)")
axs[2].legend(fontsize=7); axs[2].set_ylabel("e-/s/pix")
axs[3].plot(t[m] - TC_PRED, ph21["cx"][m] - np.nanmedian(ph21["cx"][m]), ".", label="centroid col")
axs[3].plot(t[m] - TC_PRED, ph21["cy"][m] - np.nanmedian(ph21["cy"][m]), ".", label="centroid row")
axs[3].legend(fontsize=7); axs[3].set_ylabel("pix")
for k in ["TESS-SPOC_flux", "QLP_flux_allq", "own_ap2.0"]:
    tt, ff, cad = prod[k]
    grid, chi, A, c0, rms0 = scan(tt, ff, cad, TC_PRED, half=0.5)
    axs[4].plot(grid - TC_PRED, (chi - c0) / rms0 ** 2, label=k)
axs[4].axvspan(-SIG_PRED, SIG_PRED, color="r", alpha=0.3)
axs[4].set_ylabel("Δχ² vs flat (mid-time scan)"); axs[4].legend(fontsize=7)
for a in axs[:4]:
    a.axvspan(-T14 / 2, T14 / 2, color="r", alpha=0.07)
axs[4].set_xlabel(f"t - {TC_PRED:.3f} [d]")
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig_vet_s21.png"), dpi=100)

# image with apertures/background mask
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
ax[0].imshow(np.log10(np.clip(ph21["med_img"], 1, None)), origin="lower")
ax[0].plot(ph21["x0"], ph21["y0"], "r+")
ax[0].set_title("S21 median image (log)")
ax[1].imshow(ph21["bmask"], origin="lower"); ax[1].set_title("background pixels")
plt.savefig(os.path.join(OUT, "fig_s21_image.png"), dpi=90)
json.dump(results, open(os.path.join(OUT, "vet_s21.json"), "w"), indent=1)
