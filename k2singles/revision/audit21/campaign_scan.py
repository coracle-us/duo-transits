"""
Campaign-wide companion to audit21.py: for each of the 21 non-supported catalogued K2 events, run
k2lib.test_window_k2 (catalogued depth and duration) on a grid of step D/8 across the WHOLE campaign in three
products (PDC_all, EVEREST without K2-hard-flagged cadences, K2SFF) and list the deepest 'supported' windows (distinct by > D).
Purpose: catch catalogue epoch errors larger than the +-1.5 d plot window, and catalogue depths that belong to a
different event (e.g. primary vs secondary eclipse).
Output: k2singles/revision/audit21/campaign_scan.csv (top 3 dips per event x product) and
plots/<EPIC>_campaign.png (whole-campaign overview of the three products).
Usage:  python k2singles/revision/audit21/campaign_scan.py
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import audit21 as A  # noqa: E402

warnings.filterwarnings("ignore")


def main():
    rows = []
    for _, r in A.targets().iterrows():
        epic, camp, T = int(r.EPIC), int(r.campaign), float(r.epoch_btjd)
        D, delta = r.dur_h / 24, r.depth_ppm * 1e-6
        # load_products windows +-4 d; reload the full campaign by temporarily widening the window
        prods = full_products(epic, camp, T, D)
        overview_plot(epic, camp, T, D, delta, prods)
        for name, (t, f, flag) in prods.items():
            sg = A._seg(t)
            grid = np.arange(t.min() + D, t.max() - D, D / 8)
            res = [A.k2lib.test_window_k2(t, f, x, D, delta, sg) for x in grid]
            dep = np.array([x.get("depth", np.nan) for x in res])
            snr = np.array([x.get("snr", np.nan) for x in res])
            sup = np.array([x["verdict"] == "supported" for x in res])
            order = np.argsort(-np.where(sup, dep, -np.inf))
            picked = []
            for k in order:
                if not sup[k] or len(picked) == 3:
                    break
                if all(abs(grid[k] - grid[j]) > D for j in picked):
                    picked.append(k)
            for rank, k in enumerate(picked):
                rows.append(dict(EPIC=epic, campaign=camp, product=name, rank=rank + 1, T_dip=grid[k],
                                 offset_d=grid[k] - T, ratio=dep[k] / delta, snr=snr[k], cat_epoch=T,
                                 cat_depth_ppm=r.depth_ppm, dur_h=r.dur_h))
            if not picked:
                rows.append(dict(EPIC=epic, campaign=camp, product=name, rank=0, T_dip=np.nan, offset_d=np.nan,
                                 ratio=np.nanmax(dep) / delta, snr=np.nan, cat_epoch=T, cat_depth_ppm=r.depth_ppm,
                                 dur_h=r.dur_h))
        print(epic, "done", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(HERE / "campaign_scan.csv", index=False)
    pd.set_option("display.width", 200)
    print(d.round(3).to_string(index=False))


def overview_plot(epic, camp, T, D, delta, prods):
    """Whole-campaign view of each product, catalogued epoch (red) and depth (orange) marked."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(len(prods), 1, figsize=(11, 2.3 * len(prods)), sharex=True, squeeze=False)
    for ax, (name, (t, f, _)) in zip(axs[:, 0], prods.items()):
        ax.plot(t, f, ".", ms=1.5, color="k")
        ax.axvline(T, color="tab:red", lw=0.8, alpha=0.6)
        ax.axhline(1 - delta, color="tab:orange", ls="--", lw=0.8)
        lo, hi = np.nanpercentile(f, [0.2, 99.8])
        ax.set_ylim(min(lo, 1 - 1.3 * delta), hi)
        ax.set_ylabel(name, fontsize=8)
    axs[-1, 0].set_xlabel("BTJD")
    axs[0, 0].set_title(f"EPIC {epic} C{camp} whole campaign; red = catalogued epoch {T:.3f}; orange = catalogued "
                        f"depth {delta * 1e6:.0f} ppm", fontsize=9)
    fig.tight_layout()
    fig.savefig(A.PLOTS / f"{epic}_campaign.png", dpi=90)
    plt.close(fig)


def full_products(epic, camp, T, D):
    """PDC_all, EVEREST (all finite cadences) and K2SFF over the whole campaign, each normalised by its median."""
    from astropy.io import fits
    out = {}
    cdir = A.ROOT / "k2singles" / "cache" / str(epic) / "dl"
    for fn in sorted(cdir.glob(f"ktwo{epic}-c{camp:02d}*_llc.fits")):
        with fits.open(fn) as h:
            d = h[1].data
            t = np.asarray(d["TIME"], float) + A.BK
            f = np.asarray(d["PDCSAP_FLUX"], float)
        ok = np.isfinite(t) & np.isfinite(f)
        out["PDC_all"] = (t[ok], f[ok] / np.median(f[ok]), np.zeros(ok.sum(), bool))
    hl = A.DATA / "mastDownload" / "HLSP"
    for p in hl.glob(f"hlsp_everest_k2_llc_{epic}-c{camp:02d}*/*.fits"):
        with fits.open(p) as h:
            a = h["ARRAYS"].data
            t = np.asarray(a["TIME"], float) + A.BK
            f = np.asarray(a["FLUX"], float)
            q = np.asarray(a["QUALITY"], np.int64)
        from lightkurve.utils import KeplerQualityFlags
        hard = (q & (KeplerQualityFlags.DEFAULT_BITMASK | A.k2lib.THRUSTER_POSSIBLE)) != 0
        ok = np.isfinite(t) & np.isfinite(f) & ~hard
        out["EVEREST"] = (t[ok], f[ok] / np.median(f[ok]), np.zeros(ok.sum(), bool))
    for p in hl.glob(f"hlsp_k2sff_k2_lightcurve_{epic}-c{camp:02d}*/*.fits"):
        with fits.open(p) as h:
            a = h["BESTAPER"].data
            t = np.asarray(a["T"], float) + A.BK
            f = np.asarray(a["FCOR"], float)
        ok = np.isfinite(t) & np.isfinite(f)
        out["K2SFF"] = (t[ok], f[ok] / np.median(f[ok]), np.zeros(ok.sum(), bool))
    for k, (t, f, fl) in list(out.items()):
        o = np.argsort(t)
        out[k] = (t[o], f[o], fl[o])
    return out


if __name__ == "__main__":
    main()
