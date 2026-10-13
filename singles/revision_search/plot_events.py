"""
Vetting plots for extra events from the unconstrained search.

For each TIC: one row per event (the catalogued event first, then the extra events in
extra_events_auto.csv that are not category d, or all with --all), one column per product,
raw normalised flux with a local quadratic baseline divided out (fitted outside the event).

Usage: python singles/revision_search/plot_events.py TIC [TIC ...] [--all] [--times T1,T2,...]
Output: singles/revision_search/data/plots/<TIC>_events.png
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import search_any as sa  # noqa: E402

OUT = HERE / "data" / "plots"
PRODS = ["SPOC_120", "TESS-SPOC", "QLP", "OWN"]


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    show_all = "--all" in sys.argv
    extra_times = []
    for a in sys.argv[1:]:
        if a.startswith("--times="):
            extra_times = [float(x) for x in a.split("=")[1].split(",")]
    OUT.mkdir(parents=True, exist_ok=True)
    ex = pd.read_csv(HERE / "extra_events_auto.csv")
    for tic in map(int, args):
        j = json.loads((HERE / "data" / "per_star" / f"{tic}.json").read_text())
        e = ex[ex.TIC == tic]
        if not show_all:
            e = e[e.category != "d"]
        rows = [("catalogued", j["event_T"], j["dur_h"] / 24, j["depth_ppm"] * 1e-6)]
        rows += [(f"{r.category} {r.category_reason[:22]}", r.T_btjd, r.dur_fit_h / 24, r.depth_fit_ppm * 1e-6)
                 for r in e.itertuples()]
        rows += [("manual", T, j["dur_h"] / 24, j["depth_ppm"] * 1e-6) for T in extra_times]
        rows = rows[:16]
        prods, _ = sa.load_products(tic)
        fig, axs = plt.subplots(len(rows), len(PRODS), figsize=(13, 2.1 * len(rows)), squeeze=False)
        for i, (lab, T, D, dep) in enumerate(rows):
            W = max(1.0, 2.5 * D)
            for k, key in enumerate(PRODS):
                ax = axs[i, k]
                for sec, av in prods.items():
                    if key not in av:
                        continue
                    p = av[key]
                    m = np.abs(p["t"] - T) < W
                    if m.sum() < 5:
                        continue
                    tt, ff = p["t"][m], p["f"][m]
                    oot = np.abs(tt - T) > D / 2 + 0.05
                    if oot.sum() > 10:
                        ff = ff / np.polyval(np.polyfit(tt[oot] - T, ff[oot], 2), tt - T)
                    ax.plot((tt - T) * 24, (ff - 1) * 1e3, ".", ms=1.5, color="0.4")
                    nb = max(1, int(round(0.5 / 24 / np.median(np.diff(tt)))))
                    if nb > 1 and len(tt) > 2 * nb:
                        k2 = len(tt) // nb * nb
                        ax.plot((tt[:k2].reshape(-1, nb).mean(1) - T) * 24, (ff[:k2].reshape(-1, nb).mean(1) - 1) * 1e3,
                                "o", ms=2.5, color="C3")
                    ax.set_title(f"{key} s{sec}", fontsize=7)
                ax.axvspan(-D * 12, D * 12, color="C0", alpha=0.1)
                ax.axhline(-dep * 1e3, color="C0", lw=0.6, ls="--")
                ax.set_ylim(-3 * max(dep, 3e-4) * 1e3, 1.5 * max(dep, 3e-4) * 1e3)
                ax.tick_params(labelsize=6)
                if k == 0:
                    ax.set_ylabel(f"{lab}\nT={T:.2f}", fontsize=6)
        fig.suptitle(f"TIC {tic} ({j['cls']}); x: hours from T, y: ppt", fontsize=9)
        fig.tight_layout()
        fig.savefig(OUT / f"{tic}_events.png", dpi=75)
        plt.close(fig)
        print(OUT / f"{tic}_events.png")


if __name__ == "__main__":
    main()
