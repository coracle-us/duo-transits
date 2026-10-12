"""DASCH DR7 (Harvard plates, ~1890-1990) light curve of APASS_J063155.2+252856 = EPIC 202126877, from the public Starglass API.
Quality cuts: calibrated magcal_magdep present, magcal_local_rms < 0.3, limiting_mag_local > mag + 0.75, aflags == 0, reject_flag == 0,
plate_center_dist_deg < 15, no blend (blended_mag empty). Look for a faint tail (eclipses) by decade; scan for periodicity near 21.77 d."""
import json, numpy as np, pandas as pd, io, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from astropy.time import Time
raw = json.load(open("data/dasch/lc_raw.json"))
d = pd.read_csv(io.StringIO("\n".join(raw)))
print("rows", len(d))
d = d[d.magcal_magdep.notna()]
print("with magcal_magdep", len(d))
q = d[(d.magcal_local_rms < 0.3) & (d.limiting_mag_local > d.magcal_magdep + 0.75) & (d.aflags == 0) & (d.reject_flag == 0)
      & (d.plate_center_dist_deg < 15) & (d.blended_mag.isna())].copy()
print("after cuts", len(q))
q["year"] = Time(q.date_jd, format="jd").decimalyear
q["dm"] = q.magcal_magdep - np.median(q.magcal_magdep)
sig = 1.48*np.median(np.abs(q.dm))
print(f"median mag {np.median(q.magcal_magdep):.3f}; robust sigma {sig:.3f}")
q["decade"] = (q.year//10*10).astype(int)
summ = q.groupby("decade").dm.agg(["size", lambda x: 1.48*np.median(np.abs(x-np.median(x))), lambda x: (x > 0.3).sum(), lambda x: (x > 0.4).sum()])
summ.columns = ["n", "rsig", "n_faint>0.3", "n_faint>0.4"]; print(summ)
# expected number of points in a ~30%-deep part of the eclipse if eclipsing like 2023-26: fraction of time with depth > 0.25 ~ 0.17 d / 21.77 d (P+S)
print("expected fraction of time with depth>0.25 mag-equivalent (2026-like geometry, primary only): ~", round(0.17/21.77, 4))
q.to_csv("data/dasch/dasch_clean.csv", index=False)
faint = q[q.dm > 0.3].sort_values("date_jd")
print("faint points (>0.3 mag below median):"); print(faint[["date_jd", "year", "series", "plate_number", "magcal_magdep", "magcal_magdep_rms", "limiting_mag_local", "dm"]].to_string(index=False))
fig, ax = plt.subplots(2, 1, figsize=(12, 7))
ax[0].plot(q.year, q.magcal_magdep, "k.", ms=2); ax[0].invert_yaxis(); ax[0].set_title("DASCH DR7 EPIC 202126877 (quality-cut)")
ax[1].hist(q.dm, bins=np.arange(-1, 1.01, 0.025), log=True); ax[1].set_xlabel("mag - median")
plt.tight_layout(); plt.savefig("figures/dasch_lightcurve.png", dpi=80)
