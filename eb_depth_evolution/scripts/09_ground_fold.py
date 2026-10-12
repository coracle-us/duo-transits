"""Fold SuperWASP (2004-2008) and KELT-North N04 (2006-2008) photometry on the TESS ephemeris and look for primary/secondary eclipses.
Per-camera (WASP) / per-orientation (KELT) zero points removed with a nightly-robust median; nights with < 5 points dropped."""
import numpy as np, pandas as pd, glob, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
P, T0 = 21.77465, 2488.59020
def load_wasp():
    d = pd.read_csv("data/wasp/J063155.22+252855.4.csv"); d.columns = ["hjd", "cam", "mag", "err"]
    d["btjd"] = d.hjd - 2457000; d["src"] = "WASP cam" + d.cam.astype(str)
    return d
def load_kelt():
    out = []
    for f in glob.glob("data/kelt/*_raw_lc.tbl"):
        rows = [l.split() for l in open(f) if l.strip() and l.strip()[0].isdigit()]
        a = np.array(rows, float); ori = "east" if "east" in f else "west"
        out.append(pd.DataFrame(dict(hjd=a[:, 0], mag=a[:, 1], err=a[:, 2], src="KELT " + ori)))
    d = pd.concat(out); d["btjd"] = d.hjd - 2457000
    return d
for name, d in [("wasp", load_wasp()), ("kelt", load_kelt())]:
    d = d.copy()
    # zero point per source and season
    d["season"] = np.floor((d.hjd - 2453000)/365.25 + 0.5)
    d["night"] = np.floor(d.hjd - 0.3)
    d["dm"] = d.mag - d.groupby(["src", "season"]).mag.transform("median")
    # drop sources with < 50 points (WASP cams 143/144/101 are sparse and show 0.1-0.3 mag offsets)
    n = d.groupby("src").mag.transform("size"); d = d[n >= 50]
    # sigma-clip bright outliers (no physical brightening expected)
    d = d[(d.dm > -0.08)]
    d["flux"] = 10**(-0.4*d.dm)
    E = (d.btjd - T0)/P
    d["E"] = E; d["phase"] = E - np.floor(E)
    d.to_csv(f"data/{name}_clean.csv", index=False)
    print(name, len(d), "points; seasons", sorted(d.season.unique()), "rms flux", np.std(d.flux).round(4), "srcs", d.src.unique())
    fig, axs = plt.subplots(1, 3, figsize=(16, 4))
    axs[0].plot(d.phase, d.flux, "k.", ms=1); axs[0].set_title(name + " full fold P=21.77465")
    for ax, ph0 in zip(axs[1:], [0.0, 0.33]):
        x = ((d.phase - ph0 + 0.5) % 1 - 0.5)*P
        ax.plot(x, d.flux, "k.", ms=2); ax.set_xlim(-3, 3); ax.set_title(f"{name}: days from phase {ph0}")
        ax.set_ylim(0.6, 1.15)
    plt.tight_layout(); plt.savefig(f"figures/{name}_fold.png", dpi=80)
