"""Quick look at the two known events in every product (raw normalised flux, PDCSAP/QLP, and SAP)."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lcio import load, OUT
import os

events = [(14, 1703.83), (41, 2423.55)]
fig, axes = plt.subplots(2, 2, figsize=(13, 8))
for i, (s, t0) in enumerate(events):
    for j, which in enumerate(["flux", "sap"]):
        ax = axes[i, j]
        off = 0
        for p in ["SPOC", "TESS-SPOC", "QLP"]:
            d = load(s, p, which)
            if d is None:
                continue
            m = np.abs(d.time - t0) < 2.0
            x, y = d.time[m].values, d.f[m].values
            if p == "SPOC":  # bin to 30 min for display
                bins = np.arange(t0 - 2, t0 + 2, 1 / 48)
                idx = np.digitize(x, bins)
                xb = np.array([x[idx == k].mean() for k in np.unique(idx)])
                yb = np.array([y[idx == k].mean() for k in np.unique(idx)])
                ax.plot(x - t0, y + off, ",", color="0.7")
                x, y = xb, yb
            ax.plot(x - t0, y + off, ".", ms=4, label=p)
            off -= 0.006
        ax.axvline(0, color="r", lw=0.5)
        ax.set_title(f"S{s} {which} around {t0}")
        ax.legend(fontsize=7)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig_known_events_quicklook.png"), dpi=110)
