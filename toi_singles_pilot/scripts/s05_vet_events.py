"""
Vet a TOI's catalogued event and its repeats: dedicated trapezoid fits (k2singles/scripts/eclipse_times.fit_event)
in EVERY available product per sector (best product plus QLP / TESS-SPOC / SPOC alternatives, fetched into
data/lightcurves), depth/T14 consistency, an optional linear ephemeris for a trial period, and a plot.

Usage:  python toi_singles_pilot/scripts/s05_vet_events.py SET TOI T1 T2 ... [--P P_trial]
        (T1 = catalogued event; SET = pilot | validation)
Output: toi_singles_pilot/results/<SET>/vet_<TOI>.json, toi_singles_pilot/plots/vet_<TOI>.png
"""
import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common  # noqa: E402
import s01_run  # noqa: E402
import eclipse_times as et  # noqa: E402  (k2singles/scripts, on sys.path via common)


def load_all_products(tic):
    """Every product file per sector, each normalised and segmented like load_lightcurves."""
    out = {}
    for f in (common.LCD / str(tic)).glob("s*.csv.gz"):
        m = re.fullmatch(r"s(\d+)_(.+)_(\d+)s\.csv\.gz", f.name)
        if not m:
            continue
        sec, author, exp = int(m.group(1)), m.group(2), int(m.group(3))
        df = pd.read_csv(f, usecols=["time", "flux", "quality"])
        df = df[np.isfinite(df.time) & np.isfinite(df.flux) & (df.quality == 0)].sort_values("time")
        if len(df) < 50:
            continue
        out[(sec, f"{author}_{exp}")] = (df.time.to_numpy(), (df.flux / np.nanmedian(df.flux)).to_numpy())
    return out


def main():
    a = argparse.ArgumentParser()
    a.add_argument("set")
    a.add_argument("toi")
    a.add_argument("times", nargs="+", type=float)
    a.add_argument("--P", type=float)
    a.add_argument("--no-fetch", action="store_true")
    o = a.parse_args()
    row = next(r for r in s01_run.targets(o.set) if r["TOI"] == o.toi)
    tic, D = int(row["TIC"]), row["dur_h"] / 24
    if not o.no_fetch:
        lc0 = common.ap.load_lightcurves(tic, products=common.PRODUCTS)
        secs = {int(lc0.sector.iloc[int(np.argmin(np.abs(lc0.time.to_numpy() - T)))]) for T in o.times}
        common.fetch_lightcurves(tic, extra_authors=("QLP", "TESS-SPOC", "SPOC"), only_sectors=secs)
    prods = load_all_products(tic)
    sibs = common.sibling_ephemerides(tic, o.toi)
    fits = []
    for T in o.times:
        for (sec, prod), (t, f) in sorted(prods.items()):
            if not (t.min() - 0.3 < T < t.max() + 0.3):
                continue
            lc = pd.DataFrame(dict(time=t, flux=f))
            lc, _ = common.mask_siblings(lc, sibs)
            r = et.fit_event(lc.time.to_numpy(), lc.flux.to_numpy(), T, D)
            if r is None:
                continue
            r.update(guess=T, sector=sec, product=prod)
            fits.append(r)
            print(f"{T:10.3f} s{sec} {prod:12s} T={r['T']:.4f}+/-{r['sigma_T']:.4f} depth {r['depth'] * 1e6:6.0f}"
                  f"+/-{r['sigma_depth'] * 1e6:4.0f} ppm T14 {r['T14_h']:.2f} h fin {r['fin']:.2f} beta {r['beta']:.2f}")
    res = dict(TOI=o.toi, TIC=tic, toi_depth_ppm=row["depth_ppm"], toi_dur_h=row["dur_h"],
               fits=[{k: v for k, v in x.items() if not k.startswith("_")} for x in fits])
    if o.P and fits:
        best = {}
        for x in fits:
            if x["guess"] not in best or x["sigma_T"] < best[x["guess"]]["sigma_T"]:
                best[x["guess"]] = x
        Ts = np.array([b["T"] for b in best.values()])
        s = np.array([b["sigma_T"] for b in best.values()])
        n = np.round((Ts - Ts[0]) / o.P)
        A = np.vstack([np.ones_like(n), n]).T
        W = np.diag(1 / s ** 2)
        cov = np.linalg.inv(A.T @ W @ A)
        T0, P = cov @ A.T @ W @ Ts
        rs = Ts - (T0 + n * P)
        chi2 = float(np.sum((rs / s) ** 2))
        dof = len(Ts) - 2
        sc = np.sqrt(chi2 / dof) if dof > 0 and chi2 > dof else 1.0
        res["ephemeris"] = dict(P=float(P), sigma_P=float(np.sqrt(cov[1, 1]) * sc), T0=float(T0),
                                sigma_T0=float(np.sqrt(cov[0, 0]) * sc), cycles=n.astype(int).tolist(),
                                oc_min=(rs * 1440).round(1).tolist(), chi2=chi2, dof=dof)
        print("ephemeris", res["ephemeris"])
    out = common.PILOT / "results" / o.set
    (out / f"vet_{o.toi}.json").write_text(json.dumps(res, indent=1, default=float))
    k = len(fits)
    if k:
        ncol = min(k, 4)
        nrow = int(np.ceil(k / ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 2.8 * nrow), squeeze=False)
        for ax, x in zip(axs.ravel(), fits):
            ax.plot((x["_tt"] - x["T"]) * 24, (x["_ff"] - 1) * 1e6, ".", ms=2, color="0.5")
            xh = (x["_tt"] - x["T"]) * 24
            edges = np.arange(xh.min(), xh.max() + 0.5, 0.5)
            idx = np.digitize(xh, edges)
            bx = [xh[idx == i].mean() for i in np.unique(idx) if (idx == i).sum() >= 2]
            by = [(x["_ff"][idx == i].mean() - 1) * 1e6 for i in np.unique(idx) if (idx == i).sum() >= 2]
            ax.plot(bx, by, "o", ms=3, color="b")
            ax.plot(xh, (x["_model"] - 1) * 1e6, "r-", lw=1)
            lim = 4 * max(row["depth_ppm"], 1.5 * x["depth"] * 1e6)
            ax.set_ylim(-lim, lim / 2)
            ax.set_title(f"s{x['sector']} {x['product']} T={x['T']:.3f}\n{x['depth'] * 1e6:.0f} ppm, "
                         f"T14 {x['T14_h']:.1f} h", fontsize=8)
            ax.set_xlabel("h from mid", fontsize=8)
        for ax in axs.ravel()[k:]:
            ax.axis("off")
        fig.suptitle(f"TOI-{o.toi} (TIC {tic}); TOI depth {row['depth_ppm']:.0f} ppm, T14 {row['dur_h']:.1f} h",
                     fontsize=9)
        fig.tight_layout()
        fig.savefig(common.PILOT / "plots" / f"vet_{o.toi}.png", dpi=100)
        plt.close(fig)


if __name__ == "__main__":
    main()
