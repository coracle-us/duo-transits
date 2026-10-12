"""(a) Mid-time of the S21 event from ingress-side data only (exclude the orbit-end background ramp).
(b) Empirical false-alarm rate: slide the fixed-shape transit (full or partial) over all FFI light curves of this
    star (TESS-SPOC and QLP, q=0) and count dips as significant as the S21 one, away from the known events.
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import shape
from lcio import load, OUT
from phot import photometry

E = json.load(open(os.path.join(OUT, "epochs.json")))
T14 = E["T14_d"]
TC_PRED = [r for r in E["epochs"] if r["n"] == 4 and r["k"] == 1][0]["tc"]
res = {}


def fit_scan(t, f, cad, lo, hi, step=0.001, deg=1, depth=None):
    grid = np.arange(lo, hi, step)
    chi = np.full(grid.size, np.nan)
    A = np.full(grid.size, np.nan)
    for i, tc in enumerate(grid):
        S = shape(t, tc, cad)
        if np.sum(S > 0.5) < 2:
            continue
        x = t - TC_PRED
        cols = [x ** j for j in range(deg + 1)]
        if depth is None:
            X = np.vstack(cols + [-S]).T
            b, *_ = np.linalg.lstsq(X, f, rcond=None)
            r = f - X @ b
            A[i] = b[-1]
        else:
            X = np.vstack(cols).T
            b, *_ = np.linalg.lstsq(X, f + depth * S, rcond=None)
            r = f + depth * S - X @ b
        chi[i] = np.sum(r ** 2)
    return grid, chi, A


# (a) ingress-only timing, data t < 1883.90 (before background ramp), several products
series = {}
for p in ["TESS-SPOC", "QLP"]:
    d = load(21, p, "flux", strict_quality=False)
    m = (d.time > TC_PRED - 1.2) & (d.time < 1883.90) & (d.f > 0.9)
    series[p + "_flux"] = (d.time[m].values, d.f[m].values / np.median(d.f[m].values), 1800)
    d = load(21, p, "sap", strict_quality=False)
    m = (d.time > TC_PRED - 1.2) & (d.time < 1883.90) & (d.f > 0.9)
    series[p + "_sap"] = (d.time[m].values, d.f[m].values / np.median(d.f[m].values), 1800)
ph = photometry(21)
for key in ["ap1.5", "ap2.0", "ap2.5"]:
    t, f, q = ph["time"], ph[key], ph["quality"]
    m = (t > TC_PRED - 1.2) & (t < 1883.90) & (q == 0)
    series["own_" + key] = (t[m], f[m] / np.median(f[m]), 1800)
res["ingress_timing"] = {}
print("Ingress-side timing (data before 1883.90):")
for k, (t, f, cad) in series.items():
    for depth in [None, 0.0030]:
        g, c, A = fit_scan(t, f, cad, TC_PRED - 0.15, TC_PRED + 0.15, depth=depth)
        rms = np.sqrt(np.nanmin(c) / (t.size - 3))
        cc = c / rms ** 2
        i = np.nanargmin(cc)
        ok = cc < cc[i] + 1
        lab = "free depth" if depth is None else "depth fixed 3.0 ppt"
        res["ingress_timing"][f"{k} ({lab})"] = dict(tc=float(g[i]), lo=float(g[ok].min()), hi=float(g[ok].max()),
                                                    depth=float(A[i]) if depth is None else 0.003)
        print(f"  {k:16s} {lab:20s} tc={g[i]:.4f} [{g[ok].min():.4f},{g[ok].max():.4f}]  Δ={(g[i]-TC_PRED)*24:+.2f} h"
              + (f" depth={A[i]*1e3:.2f} ppt" if depth is None else ""))

# (b) empirical false-alarm search
print("\nEmpirical search for transit-like dips over all FFI data (q=0):")
hits = []
allstat = []
total_days = 0
for p in ["TESS-SPOC", "QLP"]:
    for s in [14, 15, 16, 21, 22, 23, 41, 48, 49, 75, 76]:
        d = load(s, p, "flux", strict_quality=True)
        if d is None:
            continue
        t, f = d.time.values, d.f.values
        cad = d.cad.iloc[0]
        total_days += t.size * cad / 86400
        # trial mid-times every 0.02 d; also allow events partially outside the data
        for tc in np.arange(t.min() - 0.15, t.max() + 0.15, 0.02):
            m = np.abs(t - tc) < 1.0
            if m.sum() < 20:
                continue
            tt, ff = t[m], f[m]
            S = shape(tt, tc, cad)
            nin = np.sum(S > 0.5)
            if nin < max(5, int(0.15 * T14 * 86400 / cad)):
                continue
            # require some out-of-transit baseline on at least one side
            if np.sum(S == 0) < 10:
                continue
            x = tt - tc
            X = np.vstack([np.ones_like(x), x, -S]).T
            b, *_ = np.linalg.lstsq(X, ff, rcond=None)
            r = ff - X @ b
            rms = np.std(r, ddof=3)
            cov = np.linalg.inv(X.T @ X) * rms ** 2
            snr = b[-1] / np.sqrt(cov[-1, -1])
            known = min(abs(tc - E["TC1"]), abs(tc - E["TC2"]), abs(tc - TC_PRED)) < 0.5
            allstat.append((p, s, tc, b[-1], snr, known))
            if snr > 6 and b[-1] > 0.0015 and not known:
                hits.append((p, s, round(tc, 3), round(b[-1] * 1e3, 2), round(snr, 1)))
allstat = np.array([(a[2], a[3], a[4], a[5]) for a in allstat], dtype=float)
# merge hits into distinct events (same product, within 0.3 d)
events = []
for h in sorted(hits, key=lambda x: (x[0], x[2])):
    if events and events[-1][0] == h[0] and abs(events[-1][2] - h[2]) < 0.3:
        if h[4] > events[-1][4]:
            events[-1] = h
    else:
        events.append(h)
res["fap_search"] = dict(total_days_both_products=total_days, n_events=len(events), events=events)
print(f"  data searched: {total_days:.0f} product-days; distinct dips with depth>1.5 ppt & SNR>6 away from known/S21: {len(events)}")
for e in events:
    print("   ", e)
# compare with SNR of S21 event in same search
s21 = [a for a in allstat if abs(a[0] - TC_PRED) < 0.05]
res["s21_in_search"] = [list(map(float, a)) for a in s21]
print("  S21 event in the same search (tc, depth, snr):", [(round(a[0], 3), round(a[1] * 1e3, 2), round(a[2], 1)) for a in s21])

fig, ax = plt.subplots(figsize=(10, 4))
ax.hist(allstat[allstat[:, 3] == 0, 2], bins=100, histtype="step", label="all trial epochs (not known)")
for a in s21:
    ax.axvline(a[2], color="r", lw=0.5)
ax.set_yscale("log"); ax.set_xlabel("depth SNR"); ax.legend()
plt.savefig(os.path.join(OUT, "fig_fap_snr_hist.png"), dpi=100)
json.dump(res, open(os.path.join(OUT, "ingress_fap.json"), "w"), indent=1, default=str)
