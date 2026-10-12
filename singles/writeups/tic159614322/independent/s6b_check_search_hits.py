"""Plot the single-product dips flagged by the empirical search, in both products, to see whether they are product systematics."""
import numpy as np, os, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lcio import load, OUT
regions = [(16, 1744.5, 1764.0), (75, 3343.5, 3361.5)]
fig, axs = plt.subplots(len(regions), 1, figsize=(13, 7))
for ax, (s, lo, hi) in zip(axs, regions):
    for off, p in [(0, "TESS-SPOC"), (-0.008, "QLP")]:
        for strict, mk in [(False, "x"), (True, ".")]:
            d = load(s, p, "flux", strict_quality=strict)
            m = (d.time > lo) & (d.time < hi)
            ax.plot(d.time[m], d.f[m] + off, mk, ms=2 if strict else 3, color=("C0" if p == "TESS-SPOC" else "C2") if strict else "0.7",
                    label=f"{p} {'q=0' if strict else 'flagged'}")
    ax.set_ylim(0.982, 1.006); ax.legend(fontsize=7, ncol=4); ax.set_title(f"S{s}")
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_search_hits.png"), dpi=90)
