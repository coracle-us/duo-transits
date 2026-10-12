"""Plot every C5 reduction (full campaign + zoom on the event)."""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lcio

T0 = 179.89
prods = {
    "EVEREST": lcio.k2_everest(5),
    "K2SFF": lcio.k2_k2sff(5),
    "K2SC": lcio.k2_k2sc(5),
    "SAP": lcio.k2_sap(5),
    "PDC": lcio.k2_pdc(5),
}
fig, ax = plt.subplots(len(prods), 2, figsize=(14, 2.6 * len(prods)))
for i, (k, d) in enumerate(prods.items()):
    ax[i, 0].plot(d["t"], d["f"], ".", ms=2)
    ax[i, 0].axvline(T0, color="r", alpha=0.3)
    ax[i, 0].set_ylabel(k)
    lo, hi = np.nanpercentile(d["f"], [0.5, 99.5]); ax[i, 0].set_ylim(lo - 0.002, hi + 0.002)
    w = np.abs(d["t"] - T0) < 1.5
    ax[i, 1].plot(d["t"][w], d["f"][w], "o", ms=3)
    ax[i, 1].axvline(T0, color="r", alpha=0.3)
    print(k, "n=%d  rms(ppm, global pt-to-pt/sqrt2)=%.0f" % (len(d["t"]), 1e6 * np.nanmedian(np.abs(np.diff(d["f"]))) * 1.4826 / np.sqrt(2)))
plt.tight_layout()
plt.savefig(os.path.join(lcio.HERE, "fig_k2_c5_products.png"), dpi=110)
