import glob, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
P, T0 = 21.77465, 2488.59020
fs = sorted(glob.glob("data/phot/tess_s*.csv"), key=lambda s:int(s.split("_s")[1][:-4]))
fig, axs = plt.subplots(len(fs), 1, figsize=(14, 2.6*len(fs)))
for ax, f in zip(axs, fs):
    d = pd.read_csv(f); m = d.q == 0
    for r, c in [("1.5", "k"), ("2.5", "r")]:
        y = d[f"f_r{r}"]/np.nanmedian(d[f"f_r{r}"][m])
        ax.plot(d.time[m], y[m] + (0 if r=="1.5" else -0.0), ".", ms=1, color=c, label=f"r={r}")
    ax2 = ax.twinx(); ax2.plot(d.time[m], d.bkg[m], "-", color="c", lw=0.5, alpha=0.6)
    for k in range(-5, 100):
        for ph, col in [(0, "b"), (0.33, "g")]:
            tt = T0 + (k+ph)*P
            if d.time.min()-1 < tt < d.time.max()+1: ax.axvline(tt, color=col, alpha=0.3)
    ax.set_title(f.split("/")[-1]); ax.set_ylim(0.6, 1.05)
plt.tight_layout(); plt.savefig("figures/tess_raw_lightcurves.png", dpi=90)
