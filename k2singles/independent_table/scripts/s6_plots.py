"""Overview + zoom plots per star, and a merged list of candidate dips across products."""
import sys, os, pickle, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import lib
from targets import T

SCAN = os.path.join(lib.ROOT, "data", "scan")
PL = os.path.join(lib.ROOT, "plots")
PREF = ["K2pdc", "K2ev", "K2sff", "spoc", "tspoc", "qlp", "tglc", "own", "eleanor", "qlpk", "tars", "t16", "cdips", "K2sap"]


def clusters(scan, dur):
    rows = []
    for x in scan:
        m = x["det"]
        if not m.any(): continue
        t0 = x["t0"][m]; snr = (x["depth"] / x["err"])[m]; d = x["depth"][m]
        groups = np.split(np.arange(len(t0)), np.where(np.diff(t0) > dur)[0] + 1)
        for g in groups:
            i = g[np.argmax(snr[g])]
            rows.append(dict(t0=t0[i], depth=d[i], snr=snr[i], src=x["src"], kind=x["kind"]))
    if not rows: return pd.DataFrame(columns=["t0", "depth", "snr", "src", "kind"]), []
    df = pd.DataFrame(rows).sort_values("t0").reset_index(drop=True)
    cl = np.cumsum(np.r_[0, np.diff(df.t0) > dur])
    df["cl"] = cl
    merged = []
    for c, g in df.groupby("cl"):
        merged.append(dict(t0=np.median(g.t0), nprod=g.src.nunique(), kinds=",".join(sorted(set(g.kind))),
                           maxsnr=g.snr.max(), med_depth_ppm=np.median(g.depth) * 1e6,
                           srcs="; ".join(f"{s}:{d*1e6:.0f}/{q:.0f}" for s, d, q in zip(g.src, g.depth, g.snr))))
    return df, merged


def group_key(s):
    return ("K2" if s["kind"].startswith("K2") else "TESS", str(s["camp"]))


def main(only=None):
    for _, r in T.iterrows():
        e = int(r.epic)
        if only and e not in only: continue
        fn = os.path.join(SCAN, f"{e}.pkl")
        if not os.path.exists(fn): continue
        P = pickle.load(open(fn, "rb"))
        segs, scan = P["segs"], P["scan"]
        dur = r.dur / 24; dep = r.depth * 1e-6
        df, merged = clusters(scan, dur)
        pd.DataFrame(merged).to_csv(os.path.join(SCAN, f"{e}_dips.csv"), index=False)
        # overview: one panel per campaign/sector with the preferred product
        groups = {}
        for s in segs:
            groups.setdefault(group_key(s), []).append(s)
        keys = sorted(groups, key=lambda k: np.median([np.median(s["t"]) for s in groups[k]]))
        n = len(keys)
        fig, axes = plt.subplots(n, 1, figsize=(14, 2.2 * n + 0.5), squeeze=False)
        for ax, k in zip(axes[:, 0], keys):
            ss = sorted(groups[k], key=lambda s: PREF.index(s["kind"]) if s["kind"] in PREF else 99)
            s = ss[0]
            tb, fb = lib.bin_series(s["t"], s["f"], 1 / 24)
            ax.plot(s["t"], s["f"], ",", color="0.7")
            ax.plot(tb, fb, ".", ms=2, color="k")
            ax.axhline(1 - dep, color="r", lw=0.5, ls="--")
            lo, hi = np.nanpercentile(s["f"], [0.5, 99.5])
            ax.set_ylim(min(lo, 1 - 1.3 * dep), hi + 0.3 * (hi - lo))
            for m in merged:
                if s["t"][0] - 1 < m["t0"] < s["t"][-1] + 1:
                    ax.axvline(m["t0"], color="b" if m["nprod"] > 1 else "orange", alpha=0.5)
            ax.text(0.01, 0.85, s["src"] + "  (others: " + ", ".join(x["kind"] for x in ss[1:]) + ")", transform=ax.transAxes, fontsize=7)
        axes[-1, 0].set_xlabel("BTJD")
        fig.suptitle(f"EPIC {e} / TIC {int(r.tic)}  depth {r.depth:.0f} ppm  dur {r.dur} h  K2 epoch {r.t0}")
        fig.tight_layout()
        fig.savefig(os.path.join(PL, f"{e}_overview.png"), dpi=90)
        plt.close(fig)
        # zooms: catalogue epoch + every merged cluster
        def interesting(m):
            ks = set(m["kinds"].split(","))
            if ks <= {"K2sap", "tars", "t16", "cdips"}:
                return m["maxsnr"] > 40 and "K2sap" in ks and m["med_depth_ppm"] > 0.3 * r.depth
            return True
        times = [r.t0] + [m["t0"] for m in merged if abs(m["t0"] - r.t0) > dur and interesting(m)]
        for f_ in os.listdir(PL):
            if f_.startswith(f"{e}_zoom_"): os.remove(os.path.join(PL, f_))
        for i, tc in enumerate(times):
            zoom(e, segs, tc, dur, dep, f"{e}_zoom_{i:02d}_{tc:.2f}.png")


def zoom(e, segs, tc, dur, dep, name, halfwin=None):
    hw = halfwin or max(3 * dur, 1.0)
    inwin = [s for s in segs if np.any(np.abs(s["t"] - tc) < hw)]
    if not inwin: return
    fig, axes = plt.subplots(len(inwin), 1, figsize=(9, 1.8 * len(inwin) + 0.6), squeeze=False, sharex=True)
    for ax, s in zip(axes[:, 0], inwin):
        m = np.abs(s["t"] - tc) < hw
        ax.plot(s["t"][m] - tc, s["f"][m], ".", ms=2, color="0.5")
        cad = lib.cadence(s["t"])
        if cad < 0.015:
            tb, fb = lib.bin_series(s["t"][m], s["f"][m], 0.5 / 24)
            ax.plot(tb - tc, fb, "k.", ms=4)
        ax.axvspan(-dur / 2, dur / 2, color="y", alpha=0.2)
        med = np.nanmedian(s["f"][m])
        ax.axhline(med * (1 - dep), color="r", lw=0.5, ls="--")
        ax.text(0.01, 0.05, s["src"], transform=ax.transAxes, fontsize=7)
    axes[-1, 0].set_xlabel(f"BTJD - {tc:.3f}")
    fig.suptitle(f"EPIC {e} around {tc:.3f}")
    fig.tight_layout()
    fig.savefig(os.path.join(PL, name), dpi=80)
    plt.close(fig)


if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or None)
