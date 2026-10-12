"""Own TESScut aperture photometry (local background from far pixels) at every predicted epoch that is
partially covered, flagged (scattered light) or covered by only one product.  Uses ALL cadences with finite flux
(flags recorded), fits a fixed-shape transit at the predicted time with a quadratic local baseline, and also
reports results for each aperture and both background estimators.
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import shape
from phot import photometry
from lcio import load, OUT

E = json.load(open(os.path.join(OUT, "epochs.json")))
T14 = E["T14_d"]
targets = [  # (tc, sector, n list, note)
    (1900.116, 22, [11], "S22 start; pipeline flux flagged (4096/32768)"),
    (1909.463, 22, [7], "S22 ~90% covered"),
    (1919.744, 22, [10], "S22 full"),
    (1943.735, 23, [3, 6, 9, 12], "S23; TESS-SPOC PDCSAP NaN (4096), QLP only"),
    (2639.464, 49, [10], "S49 start; SPOC/TESS-SPOC flagged 4096"),
    (2663.455, 49, [3, 6, 9, 12], "S49 full"),
    (3339.555, 75, [11], "S75 start"),
    (3359.184, 75, [10], "S75; TESS-SPOC PDCSAP NaN (4096)"),
]
# refine predicted times from epochs.json
for i, (tc, s, ns, note) in enumerate(targets):
    best = min(E["epochs"], key=lambda r: abs(r["tc"] - tc))
    targets[i] = (best["tc"], s, ns, note)

cache = {}
out = []
fig, axs = plt.subplots(len(targets), 1, figsize=(10, 3.0 * len(targets)))
for ax, (tc, s, ns, note) in zip(axs, targets):
    if s not in cache:
        cache[s] = photometry(s)
    ph = cache[s]
    t = ph["time"]
    cad = np.median(np.diff(t)) * 86400
    win = np.abs(t - tc) < 1.0
    entry = dict(tc=tc, sector=s, ns=ns, note=note, cutout=ph["file"], fits={})
    S_all = shape(t[win], tc, cad) if win.sum() else np.array([])
    nin = int(np.sum(S_all > 0.05)) if win.sum() else 0
    entry["n_in_transit_cadences_all"] = nin
    entry["n_in_transit_q0"] = int(np.sum((S_all > 0.05) & (ph["quality"][win] == 0))) if win.sum() else 0
    entry["flags_in_transit"] = {str(k): int(v) for k, v in zip(*np.unique(ph["quality"][win][S_all > 0.05], return_counts=True))} if nin else {}
    off = 0
    for key in ["ap1.5", "ap2.0", "ap2.5", "ap2.0_annbkg"]:
        f = ph[key]
        m = win & np.isfinite(f) & (f > 0)
        if m.sum() < 10:
            continue
        tt, ff = t[m], f[m] / np.median(f[m])
        # drop >5-sigma outliers
        dev = np.abs(ff - np.median(ff))
        k = dev < 5 * 1.4826 * np.median(dev)
        tt, ff = tt[k], ff[k]
        S = shape(tt, tc, cad)
        if np.sum(S > 0.05) < 3:
            entry["fits"][key] = None
        else:
            x = tt - tc
            X = np.vstack([np.ones_like(x), x, x ** 2, -S]).T
            b, *_ = np.linalg.lstsq(X, ff, rcond=None)
            r = ff - X @ b
            rms = np.std(r, ddof=4)
            cov = np.linalg.inv(X.T @ X) * rms ** 2
            # chi2 penalty of forcing a 3.0 ppt transit
            X0 = X[:, :3]
            b0, *_ = np.linalg.lstsq(X0, ff + 0.003 * S, rcond=None)
            chi_ref = np.sum((ff + 0.003 * S - X0 @ b0) ** 2) / rms ** 2 - np.sum(r ** 2) / rms ** 2
            entry["fits"][key] = dict(depth=float(b[-1]), sig=float(np.sqrt(cov[-1, -1])), rms=float(rms),
                                      n_in=int(np.sum(S > 0.05)), dchi2_3ppt=float(chi_ref))
        tb = np.arange(tc - 1, tc + 1, 1 / 48)
        idx = np.digitize(tt, tb)
        xb = np.array([tt[idx == j].mean() for j in np.unique(idx)])
        yb = np.array([ff[idx == j].mean() for j in np.unique(idx)])
        ax.plot(tt - tc, ff + off, ",", color="0.6")
        ax.plot(xb - tc, yb + off, ".", ms=4, label=key)
        mod = 1 - 0.003 * shape(np.linspace(tc - 1, tc + 1, 600), tc, 120)
        ax.plot(np.linspace(-1, 1, 600), mod + off, "k-", lw=0.6)
        off -= 0.006
    qf = ph["quality"][win] != 0
    if qf.any():
        ax.plot(t[win][qf] - tc, np.full(qf.sum(), 1.004), "r|", ms=4, label="flagged cadence")
    ax.axvspan(-T14 / 2, T14 / 2, color="r", alpha=0.07)
    ax.set_xlim(-1, 1); ax.set_ylim(0.97, 1.006)
    ax.set_title(f"own phot S{s}: predicted {tc:.3f} (n={ns}) {note}", fontsize=9)
    ax.legend(fontsize=6, ncol=5, loc="lower left")
    out.append(entry)
    print(f"tc={tc:.3f} S{s} n={ns}  in-transit cadences all={nin} q0={entry['n_in_transit_q0']} flags={entry['flags_in_transit']}")
    for k2, v in entry["fits"].items():
        if v:
            print(f"     {k2:14s} depth={v['depth']*1e3:6.2f}±{v['sig']*1e3:4.2f} ppt  n_in={v['n_in']}  Δχ²(3ppt forced)={v['dchi2_3ppt']:.1f}")
        else:
            print(f"     {k2:14s} no in-transit data")
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig_gap_epochs_ownphot.png"), dpi=90)
json.dump(out, open(os.path.join(OUT, "gap_epochs_ownphot.json"), "w"), indent=1)
