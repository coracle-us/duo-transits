"""Step 16: overlay of all primaries and secondaries, aligned on measured mid-times (own photometry
for TESS, K2SFF for 2014), coloured by season.  -> fig_overlay.png"""
import os, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ebdata

H = os.path.dirname(os.path.abspath(__file__))
E = pd.read_csv(os.path.join(H, "per_eclipse_nominal.csv")).set_index("ecl")
k2 = json.load(open(os.path.join(H, "k2_results_nominal.json")))
df = ebdata.load("own", seasons=[2021, 2023, 2026])
col = {2014: "k", 2021: "C0", 2023: "C1", 2026: "C3"}
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for ecl, d in df.groupby("ecl"):
    if ecl not in E.index:
        continue
    a = ax[0] if ecl[0] == "P" else ax[1]
    a.plot((d.time - E.loc[ecl, "tmid"]) * 24, d.fn, ".", ms=2, color=col[int(d.season.iloc[0])])
W = pd.read_csv(os.path.join(H, "data", "windows.csv"))
w = W[(W.src == "k2sff") & (W.epoch == -124) & (W.type == "P")]
tm = k2["primary_by_pipeline"]["k2sff"]["tmin"]
oot = np.abs(w.time - tm) > 0.3
ax[0].plot((w.time - tm) * 24, w.flux / np.median(w.flux[oot]), "o", ms=4, color="k")
for s, c in col.items():
    ax[0].plot([], [], "o", color=c, label=str(s) + (" (K2, Kepler band)" if s == 2014 else ""))
ax[0].legend(); ax[0].set_xlim(-9, 9); ax[0].set_title("Primary eclipses"); ax[0].set_xlabel("hours from mid-eclipse")
ax[1].set_xlim(-9, 9); ax[1].set_title("Secondary eclipses (none in 2014; 2026 not covered)"); ax[1].set_xlabel("hours from mid-eclipse")
ax[0].set_ylabel("normalised flux (own r=2 px aperture)")
plt.tight_layout(); plt.savefig(os.path.join(H, "fig_overlay.png"), dpi=80)
