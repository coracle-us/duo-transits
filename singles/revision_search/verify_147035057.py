"""
Independent check (main session) of the agent's claim that TIC 147035057 has a deep eclipse at phase ~0.239
of P = 37.31497 d (S6 and S87), i.e. that the catalogued 3.2 ppt event is the secondary.
Uses the project's own light-curve loader, all products, and prints the depth at the predicted times.
Output: singles/revision_search/verify_147035057.png and stdout
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "singles" / "scripts"))
import alias_pipeline as ap  # noqa: E402

TIC = 147035057
P, T0 = 37.31497, 1515.0661          # BTJD; catalogued (shallow) event ephemeris, singles/writeups/tic147035057
PH = 0.239


def main():
    lc = ap.load_lightcurves(TIC, products=ap.PRODUCT_RANK + ["OWN"])
    t, f = lc.time.to_numpy(), lc.flux.to_numpy()
    n_all = np.arange(int((t.min() - T0) / P) - 2, int((t.max() - T0) / P) + 2)
    fig, axs = plt.subplots(1, 4, figsize=(14, 3.2))
    k = 0
    for n in n_all:
        for ph, lab in ((PH, "deep?"), (0.0, "catalogued")):
            tc = T0 + (n + ph) * P
            m = np.abs(t - tc) < 0.6
            if m.sum() < 20:
                continue
            for (sec, prod), g in lc[m].groupby(["sector", "product"]):
                tt, ff = g.time.to_numpy(), g.flux.to_numpy()
                oot = np.abs(tt - tc) > 0.25
                base = np.median(ff[oot]) if oot.sum() > 5 else np.median(ff)
                inn = np.abs(tt - tc) < 0.06
                dep = 1 - np.median(ff[inn]) / base if inn.sum() > 2 else np.nan
                print(f"n={n:4d} {lab:10s} T={tc:9.3f} S{int(sec):3d} {prod:10s} depth={dep * 1e3:7.2f} ppt  (n_in={inn.sum()})")
                if lab == "deep?" and k < 4 and np.isfinite(dep):
                    axs[k].plot((tt - tc) * 24, (ff / base - 1) * 1e3, ".", ms=2)
                    axs[k].set_title(f"S{int(sec)} {prod} phase {ph}", fontsize=9)
                    axs[k].set_xlabel("h")
                    k += 1
    axs[0].set_ylabel("ppt")
    fig.tight_layout()
    fig.savefig(Path(__file__).with_name("verify_147035057.png"), dpi=110)


if __name__ == "__main__":
    main()
