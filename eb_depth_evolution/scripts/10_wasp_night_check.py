"""Check suspicious WASP nights against ~29 comparison stars within 3 arcmin observed on the same frames.
For each comparison star, compute the per-frame deviation from its own season median; flag frames where the ensemble median
deviation is large (clouds, twilight). Then plot the target on suspicious nights together with the ensemble."""
import numpy as np, pandas as pd, glob, sys, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def rd(f):
    d = pd.read_csv(f); d.columns = ["hjd", "cam", "mag", "err"]; return d
tgt = rd("data/wasp/J063155.22+252855.4.csv")
comps = []
for f in glob.glob("data/wasp/comps/*.csv"):
    if "J063155.22" in f or "J063154.08" in f: continue
    d = rd(f)
    if len(d) < 1000 or 1.48*np.median(np.abs(d.mag - d.mag.median())) > 0.04: continue
    d["dm"] = d.mag - d.groupby("cam").mag.transform("median"); d["star"] = f.split("comps")[1][1:-4]
    d["magref"] = d.mag.median()
    comps.append(d)
C = pd.concat(comps)
C["key"] = C.hjd.round(5)
ens = C.groupby("key").dm.agg(["median", "size", lambda x: 1.48*np.median(np.abs(x - np.median(x)))]).rename(columns={"<lambda_0>": "sd"})
tgt["key"] = tgt.hjd.round(5)
tgt["dm"] = tgt.mag - tgt.groupby("cam").mag.transform("median")
T = tgt.merge(ens, left_on="key", right_index=True, how="left")
T["dm_corr"] = T.dm - T["median"]
T.to_csv("data/wasp_target_with_ensemble.csv", index=False)
print("frames matched", T["median"].notna().sum(), "of", len(T))
nights = [int(a) for a in sys.argv[1:]] or [2454110, 2454109, 2454087, 2454066]
fig, axs = plt.subplots(len(nights), 1, figsize=(10, 3*len(nights)))
for ax, n in zip(np.atleast_1d(axs), nights):
    m = (T.hjd > n) & (T.hjd < n + 1.2)
    ax.plot(T.hjd[m] - n, T.dm[m], "ko", ms=3, label="target dm")
    ax.plot(T.hjd[m] - n, T["median"][m], "r^", ms=3, label="ensemble median dm (comps <13.5 mag)")
    ax.invert_yaxis(); ax.legend(fontsize=7); ax.set_title(f"WASP night {n}")
    print(f"night {n}:"); print(T[m][["hjd", "cam", "mag", "dm", "median", "size", "sd", "dm_corr"]].to_string(index=False))
plt.tight_layout(); plt.savefig("figures/wasp_suspicious_nights.png", dpi=80)
