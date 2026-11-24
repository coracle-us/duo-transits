"""Plot raw normalised flux at the predicted epochs of P/2 and P/3 (where in data) and at
every predicted epoch of the adopted P, for a visual check of the alias rejections.
Usage: s08_alias_plots.py TIC   (reads out/ephem_<tic>.json)"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import TARGETS, OUT, load_lc_files, load_tesscut_phot
from s07_aliases import depth_at, RANK

tic = int(sys.argv[1])
dur, dep = TARGETS[tic]["dur"], TARGETS[tic]["depth"]
eph = json.load(open(os.path.join(OUT, f"ephem_{tic}.json")))
sh = json.load(open(os.path.join(OUT, f"shape_{tic}.json")))
P, T0 = eph["P"], eph["T0_E0"]
lcs = [lc for lc in load_lc_files(tic) + load_tesscut_phot(tic) if len(lc["t"]) > 100]
tmin = min(lc["t"].min() for lc in lcs); tmax = max(lc["t"].max() for lc in lcs)
panels = []
for div, lab in [(1, "P"), (2, "P/2"), (3, "P/3")]:
    p = P / div
    for k in np.arange(np.ceil((tmin - T0) / p), np.floor((tmax - T0) / p) + 1):
        if div > 1 and int(k) % div == 0:
            continue  # coincides with a P epoch
        tc = T0 + k * p
        m = depth_at(lcs, sh, tc, dur)
        if m is None or m[2] < 0.5:
            continue
        panels.append((lab, tc, m))
nc = 5
nr = int(np.ceil(len(panels) / nc))
fig, axs = plt.subplots(nr, nc, figsize=(18, 2.6 * nr), squeeze=False)
for ax, (lab, tc, m) in zip(axs.flat, panels):
    s, e, cov, sec, src = m
    lc = [l for l in lcs if l["sector"] == sec and l["source"] == src][0]
    w = np.abs(lc["t"] - tc) < 1.5 * dur
    t, f = lc["t"][w], lc["f"][w]
    oot = np.abs(t - tc) > 0.6 * dur
    c = np.polyfit(t[oot] - tc, f[oot], 2) if oot.sum() > 5 else [0, 0, np.median(f)]
    fn = f / np.polyval(c, t - tc)
    ax.plot((t - tc) * 24, fn, ".", ms=1.5, color="0.6")
    nb = max(1, int(round((20 / 1440) / lc["exptime"])))
    k = (len(t) // nb) * nb
    if k:
        ax.plot((t[:k].reshape(-1, nb).mean(1) - tc) * 24, fn[:k].reshape(-1, nb).mean(1), "k.", ms=3)
    ax.axvspan(-dur * 12, dur * 12, color="C1" if lab == "P" else "C3", alpha=0.12)
    ax.axhline(1 - dep, color="r", lw=0.5)
    ax.set_ylim(1 - 2 * dep, 1 + dep)
    ax.set_title(f"{lab}: Tc={tc:.3f} S{sec} {src}  s={s:.2f}±{e:.2f}", fontsize=8)
for ax in list(axs.flat)[len(panels):]:
    ax.axis("off")
fig.suptitle(f"TIC {tic}: adopted P={P:.5f} d; shaded = expected transit window")
fig.tight_layout(); fig.savefig(os.path.join(OUT, f"aliases_{tic}.png"), dpi=65)
print(f"{len(panels)} panels -> out/aliases_{tic}.png")
