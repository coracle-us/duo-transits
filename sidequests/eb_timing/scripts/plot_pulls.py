"""
Pull histograms before (formal Jacobian errors from run_etv) and after (noise model +
wild cluster bootstrap, calibrated) for the B26 NCVZ systems: L1 (P3 and A, calibration set)
and D1 (P3 only, held out). Input results/calibration.csv; output figs/pulls.png.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
c = pd.read_csv(HERE / "results" / "calibration.csv")
L1, D1 = c[c.group == "L1"], c[c.group == "D1"]
old = np.concatenate([L1.pullP_old.dropna(), L1.pullA_old.dropna()])
new = np.concatenate([L1.pullP_new_cal.dropna(), L1.pullA_new_cal.dropna()])
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
bins = np.linspace(-6, 6, 25)
x = np.linspace(-6, 6, 200)
for a, z, title in ((ax[0], old, "before: formal errors"), (ax[1], new, "after: noise model + bootstrap")):
    a.hist(np.clip(z, -5.9, 5.9), bins=bins, color="0.6", label=f"L1 P3+A (n={len(z)})")
    a.plot(x, len(z) * (bins[1] - bins[0]) * np.exp(-x ** 2 / 2) / np.sqrt(2 * np.pi), "k-", lw=1, label="N(0,1)")
    a.set_title(f"{title}\nRMS {np.sqrt(np.mean(z**2)):.2f}, robust sd {1.4826*np.median(np.abs(z)):.2f}", fontsize=9)
    a.set_xlabel("(ours - B26) / sigma  (clipped at +-6)")
zd = D1.pullP_new_cal.dropna()
ax[1].hist(np.clip(zd, -5.9, 5.9), bins=bins, histtype="step", color="C3", lw=1.5, label=f"D1 P3, held out (n={len(zd)})")
ax[1].legend(fontsize=7)
ax[0].legend(fontsize=7)
fig.savefig(HERE / "figs" / "pulls.png", dpi=100, bbox_inches="tight")
print("ok")
