"""Search TESS own photometry (r=1.0 and r=1.5) and SPOC/QLP for extra dips (tertiary eclipses) outside the binary eclipses:
mask +-0.6 d around each primary/secondary, detrend with a 1.5-d running median, list runs of >=6 consecutive cadences below 3 sigma."""
import numpy as np, pandas as pd, glob, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ev = pd.read_csv("data/eclipse_times.csv")
P = 21.7748; Tp = 2488.58926
ecl = []
for k in range(-10, 120):
    ecl += [Tp + k*P, Tp + (k + 0.330)*P]
ecl = np.array(ecl)
fig, axs = plt.subplots(6, 1, figsize=(15, 14))
for ax, fn in zip(axs, sorted(glob.glob("data/phot/tess_s*.csv"), key=lambda s: int(s.split("_s")[1][:-4]))):
    d = pd.read_csv(fn); d = d[d.q == 0]
    t = d.time.values; f = (d["f_r1.0"]/np.nanmedian(d["f_r1.0"])).values; b = d.bkg.values
    m = np.min(np.abs(t[:, None] - ecl[None, :]), axis=1) > 0.6
    dt = np.median(np.diff(t)); win = int(1.5/dt) | 1
    trend = pd.Series(np.where(m, f, np.nan)).rolling(win, center=True, min_periods=win//4).median().values
    r = f/trend - 1; s = 1.48*np.nanmedian(np.abs(r[m]))
    low = (r < -3*s) & m
    # runs
    runs = []; i = 0
    while i < len(low):
        if low[i]:
            j = i
            while j + 1 < len(low) and low[j + 1]: j += 1
            if j - i + 1 >= 6: runs.append((t[i], t[j], r[i:j+1].min()))
            i = j + 1
        else: i += 1
    print(fn, f"sigma={s*1e3:.2f} ppt; candidate dips:", [(round(a, 2), round(b_, 2), round(c*100, 2)) for a, b_, c in runs])
    ax.plot(t[m], r[m]*100, "k.", ms=1); ax.set_ylim(-3, 2); ax.set_title(fn); ax.set_ylabel("resid %")
    ax2 = ax.twinx(); ax2.plot(t, b, "c-", lw=0.4)
plt.tight_layout(); plt.savefig("figures/tess_residuals_extra_events.png", dpi=80)
