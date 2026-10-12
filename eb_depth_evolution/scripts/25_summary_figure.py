"""Summary figure: (top) primary and secondary eclipse depth vs year from K2, TESS and the 2004-2008 ground-survey limits, with the
depth implied by linear / quadratic i(t) through the base physical model; (bottom) inclination vs year."""
import numpy as np, pandas as pd, json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
T = pd.read_csv("results_depth_events.csv")
C = pd.read_csv("data/depth_vs_inclination.csv")
R = pd.read_csv("data/inclination_trend.csv"); b = R[(R.variant == "base") & (R.ph14 == 0.32)].iloc[0]
ti = np.array([2014.356, 2021.808, 2023.866, 2026.041]); ii = np.array([b.i2014, b.i2021, b.i2023, b.i2026])
lin = np.polyfit(ti - 2020, ii, 1); quad = np.polyfit(ti - 2020, ii, 2)
yrs = np.linspace(2003, 2031, 600)
fig, ax = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
for kind, col, lab in [("P", "C0", "primary"), ("S", "C2", "secondary")]:
    g = T[(T.kind == kind) & (T.coverage > 0.8)]
    ax[0].plot(g.year, g.depth_modelfree * 100, "o", color=col, label=f"{lab} (K2 / TESS, measured)")
    for nm, c, ls in [("linear i(t)", lin, "-"), ("quadratic i(t)", quad, ":")]:
        iy = np.polyval(c, yrs - 2020); iy = np.where(iy > 90, 180 - iy, iy)
        ax[0].plot(yrs, np.interp(iy, C.inc, C["d%s_tess" % kind]) * 100, ls, color=col, lw=1, alpha=0.7,
                   label=f"{lab}: model, {nm}" if kind == "P" or nm == "linear i(t)" else None)
ax[0].errorbar([2004.7, 2007.0], [4.5, 4.5], yerr=[[1.5, 1.5], [0, 0]], uplims=True, fmt="v", color="C0", label="primary 3-sigma limit, WASP+KELT")
ax[0].errorbar([2007.3], [2.3], yerr=[[0.8], [0]], uplims=True, fmt="v", color="C2", label="secondary limit, WASP+KELT (median)")
ax[0].errorbar([2014.4], [0.03], yerr=[[0.0], [0]], fmt="v", color="C2", ms=7)
ax[0].text(2014.6, 1.0, "K2 secondary < 0.03%", color="C2", fontsize=8)
ax[0].set_ylabel("eclipse depth (%)"); ax[0].legend(fontsize=7, loc="upper left"); ax[0].set_ylim(-1, 42)
ax[0].set_title("EPIC 202126877 / TIC 86491361: eclipse depth history")
ax[1].plot(ti, ii, "ko", label="physical-model inclination (K2 + TESS seasons)")
ax[1].plot(yrs, np.polyval(lin, yrs - 2020), "k-", lw=1, label=f"linear: {lin[0]:.3f} deg/yr")
ax[1].plot(yrs, np.polyval(quad, yrs - 2020), "k:", lw=1, label="quadratic")
ax[1].axhline(b.iP, color="C0", ls="--", lw=0.8, label=f"primary eclipses vanish below i={b.iP:.2f}")
ax[1].axhline(b.iS, color="C2", ls="--", lw=0.8, label=f"secondary eclipses vanish below i={b.iS:.2f}")
ax[1].axhline(90, color="gray", lw=0.5)
ax[1].set_ylim(82.5, 91); ax[1].set_ylabel("inclination (deg)"); ax[1].set_xlabel("year"); ax[1].legend(fontsize=7, loc="lower right")
plt.tight_layout(); plt.savefig("figures/depth_history_summary.png", dpi=110)
print("lin", lin, "quad", quad)
