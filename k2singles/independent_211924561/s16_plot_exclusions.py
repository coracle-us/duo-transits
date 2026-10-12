"""Plot the epochs that exclude P=46.83 (m=3) and P=35.12 (m=4), with the expected transit overlaid."""
import os, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lcio, tfit
import s14_alias_analysis as A  # re-runs the analysis (fast enough) and exposes the cached datasets

a = json.load(open(os.path.join(lcio.HERE, "alias_table.json")))
rows = []
for p in a["periods"][2:4]:
    for e in p["excluding_epochs"]:
        for ds in e["ds"]:
            rows.append((p["P"], e["tc"], ds))
fig, ax = plt.subplots(len(rows), 1, figsize=(7, 1.9 * len(rows)))
for axx, (P, tc, ds) in zip(ax, rows):
    c = A.cache[ds]
    w = np.abs(c["t"] - tc) < 0.6
    t, f = c["t"][w], c["f"][w]
    x = t - tc
    base = np.polyval(np.polyfit(x, f, 2), x)
    axx.plot(x * 24, 1e3 * (f / base - 1), ".", ms=2, color="0.5")
    xx = np.linspace(-0.6, 0.6, 1000)
    axx.plot(xx * 24, 1e3 * (tfit.model(tc + xx, tc, A.k_, A.b_, A.aR_, c["u"], c["exp"]) - 1), "r-")
    axx.set_ylim(-3.5, 1.5); axx.set_xlim(-14, 14)
    axx.text(0.01, 0.06, f"P={P:.3f}  {ds}  predicted {tc:.3f}", transform=axx.transAxes, fontsize=7)
ax[-1].set_xlabel("hours from predicted mid-time")
plt.tight_layout(); plt.savefig(os.path.join(lcio.HERE, "fig_exclusions_m3_m4.png"), dpi=90)
