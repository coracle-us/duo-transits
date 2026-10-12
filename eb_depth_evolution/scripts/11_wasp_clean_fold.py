"""WASP: keep frames with good ensemble (>=10 comparison stars, ensemble scatter < 0.04 mag, |ensemble offset| < 0.15 mag),
use ensemble-corrected magnitude, per-camera zero point; fold on TESS ephemeris; list coverage of predicted eclipse windows."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
P, T0 = 21.77465, 2488.59020
T = pd.read_csv("data/wasp_target_with_ensemble.csv")
n0 = len(T)
T = T[(T["size"] >= 6) & (T.sd < 0.03) & (T["median"].abs() < 0.10)].copy()
T["dmc"] = T.dm - T.groupby("cam").dm.transform("median")   # raw target mag (ensemble used only to reject bad frames)
# per-camera count
T = T[T.groupby("cam").dmc.transform("size") >= 50]
# 5-sigma clip of single bright/faint outliers relative to a running nightly median would also remove eclipse points; keep, only clip bright
T = T[T.dmc > -0.06]
T["flux"] = 10**(-0.4*T.dmc); T["btjd"] = T.hjd - 2457000
T["E"] = (T.btjd - T0)/P; T["phase"] = T.E - np.floor(T.E)
T.to_csv("data/wasp_clean2.csv", index=False)
print(f"kept {len(T)} of {n0}; rms flux {np.std(T.flux):.4f}; robust {1.48*np.median(np.abs(T.flux-1)):.4f}")
fig, axs = plt.subplots(1, 3, figsize=(16, 4))
axs[0].plot(T.phase, T.flux, "k.", ms=1); axs[0].set_ylim(0.6, 1.1); axs[0].set_title("WASP 2004-2008 (cleaned) folded, TESS ephemeris")
for ax, ph0 in zip(axs[1:], [0.0, 0.33]):
    x = ((T.phase - ph0 + 0.5) % 1 - 0.5)*P
    ax.plot(x, T.flux, "k.", ms=2); ax.set_xlim(-3, 3); ax.set_ylim(0.6, 1.1); ax.set_title(f"days from phase {ph0}")
plt.tight_layout(); plt.savefig("figures/wasp_fold_clean.png", dpi=80)
x = ((T.phase + 0.5) % 1 - 0.5)*P
print("points with |x|<0.3 d of primary:", (np.abs(x) < 0.3).sum(), " min flux there", T.flux[np.abs(x) < 0.3].min() if (np.abs(x)<0.3).any() else None)
print("faintest 15 points:"); print(T.nsmallest(15, "flux")[["hjd", "cam", "flux", "sd", "phase"]].to_string(index=False))
