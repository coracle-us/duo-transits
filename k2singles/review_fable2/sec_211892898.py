"""EPIC 211892898: look for a secondary eclipse at the Gaia-predicted secondary conjunction (nu+w=270)
in the K2 data, and phase-fold at P = 72.5086 d. Also report the eclipse depth/duration expectations."""
import glob
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
P, T0 = 72.508568, 189.331573
T_SEC = [206.315, 1149.251, 1294.320]  # from ecl_211892898.txt (Gaia SB1, secondary conjunction), +-0.5 d
out = open(HERE / "sec_211892898.txt", "w")


def say(*a):
    print(*a); print(*a, file=out)


parts = []
for fn in sorted(glob.glob(str(ROOT / "k2singles" / "cache" / "211892898" / "k2_c*.npz"))):
    z = np.load(fn); parts.append(pd.DataFrame(dict(t=z["t"], f=z["f"], camp=int(Path(fn).stem.split("_c")[1]))))
k2 = pd.concat(parts).sort_values("t").reset_index(drop=True)
t, f = k2.t.to_numpy(), k2.f.to_numpy()

# detrend each campaign with a running median (window 5 d) excluding +-1.5 d around primary eclipses
ph = ((t - T0) / P + 0.5) % 1 - 0.5
mask_pri = np.abs(ph * P) < 1.6
fd = np.empty_like(f)
for c in k2.camp.unique():
    m = (k2.camp == c).to_numpy()
    tt, ff = t[m], f[m]
    base = np.empty_like(ff)
    for i in range(len(tt)):
        w = (np.abs(tt - tt[i]) < 2.5) & ~mask_pri[m]
        base[i] = np.median(ff[w]) if w.sum() > 20 else np.nan
    fd[m] = ff / base
ok = np.isfinite(fd)
t, fd, ph = t[ok], fd[ok], ph[ok]

say(f"Phase-folded at P={P}, T0={T0}; primary eclipse depth (phase 0, |dt|<0.5 d): {(1-np.mean(fd[np.abs(ph*P)<0.5]))*1e3:.2f} ppt")
say("Secondary-conjunction windows (Gaia prediction +-0.5 d; test box +-0.6 d, baseline +-1..3 d):")
for Ts in T_SEC:
    inn = np.abs(t - Ts) < 0.6
    oot = (np.abs(t - Ts) > 1.2) & (np.abs(t - Ts) < 3.0)
    if inn.sum() < 10 or oot.sum() < 20:
        say(f"  {Ts}: not covered (in {inn.sum()}, oot {oot.sum()})"); continue
    d = np.mean(fd[oot]) - np.mean(fd[inn])
    sig = np.std(fd[oot]) / np.sqrt(inn.sum())
    say(f"  {Ts}: depth = {d*1e3:+.2f} +- {sig*1e3:.2f} ppt  (n_in={inn.sum()})")
# scan: binned folded light curve away from primary
say("\nFolded light curve, 1-d bins in time-from-primary (ppt, mean +- sem):")
dtp = ph * P
for lo in np.arange(-36, 36, 3.0):
    m = (dtp >= lo) & (dtp < lo + 3)
    if m.sum() > 20:
        say(f"  dt {lo:+6.1f}..{lo+3:+6.1f} d: {(1-np.mean(fd[m]))*1e3:+.2f} +- {np.std(fd[m])/np.sqrt(m.sum())*1e3:.2f}  (n={m.sum()})")
fig, ax = plt.subplots(2, 1, figsize=(10, 6))
ax[0].plot(dtp, (fd - 1) * 1e3, ".", ms=2); ax[0].set_xlabel("days from primary eclipse"); ax[0].set_ylabel("ppt")
for Ts in T_SEC:
    ax[0].axvline(((Ts - T0) / P + 0.5) % 1 * P - 0.5 * P, color="r", lw=0.5)
ax[0].set_ylim(-9, 4)
ax[1].plot(dtp, (fd - 1) * 1e3, ".", ms=2); ax[1].set_xlim(-3, 3); ax[1].set_ylim(-9, 2); ax[1].set_xlabel("days from primary")
fig.tight_layout(); fig.savefig(HERE / "sec_211892898.png", dpi=100)

# physical sanity: mass function and expectations
K, e, Porb = 22.747, 0.4801, 72.5336
fM = 1.0361e-7 * K**3 * Porb * (1 - e**2) ** 1.5
say(f"\nGaia SB1 mass function f(M) = {fM:.4f} Msun")
for M1 in (1.0, 1.2, 1.5):
    from scipy.optimize import brentq
    M2 = brentq(lambda m: m**3 / (M1 + m) ** 2 - fM, 0.01, 5)
    a = ((M1 + M2) * (Porb / 365.25) ** 2) ** (1 / 3) * 215.03
    say(f"  M1={M1}: M2(sin i=1) = {M2:.3f} Msun; a = {a:.1f} Rsun; R1=8.74 Rsun -> R1/a = {8.74/a:.3f};"
        f" central-transit duration ~ {Porb/np.pi*8.74/a*24*np.sqrt(1-e**2)/(1+e*np.sin(np.deg2rad(201.53))):.0f} h (at nu+w=90)")
say("Observed: depth ~6.2-6.9 ppt (=> R2/R1 ~ 0.08 ignoring limb darkening -> R2 ~ 0.7 Rsun), T14 ~ 51-55 h.")
out.close()
