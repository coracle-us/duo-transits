"""
Shape screen for every pipeline repeat ("repeats need a consistent depth and duration").

For the catalogued event and each repeat: a trapezoid fit (k2singles eclipse_times.fit_event) in the best
product (the one the search used). A repeat PASSES if
  - the fit converged with sigma_T < D/4 and sigma_depth < depth/3 (a real, measurable dip)
  - depth within 0.5-2x the event's fitted depth, T14 within 0.67-1.5x the event's fitted T14
  - data on both sides: >= 0.3 d (or >= T14/2) of cadences before ingress and after egress within the fit window
  - no gap > 0.1 d inside the fitted T14
Periods are then recomputed with k2lib.robust_periods from the PASSING repeats only.
Usage:  python toi_singles_pilot/scripts/s08_screen_repeats.py SET [TOI ...]
Output: toi_singles_pilot/results/<SET>/screen_<TOI>.json and toi_singles_pilot/screen_<SET>.csv,
        toi_singles_pilot/plots/screen_<TOI>.png (binned event + repeats with fits)
"""
import json
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
import eclipse_times as et  # noqa: E402


def fit(lc, T, D):
    m = np.abs(lc.time - T) < 1.5 * D + 0.3
    d = lc[m]
    if len(d) < 10:
        return None
    sec = d.sector.mode().iloc[0]
    d = d[d.sector == sec]
    r = et.fit_event(d.time.to_numpy(), d.flux.to_numpy(), T, D)
    if r is None:
        return None
    t = r["_tt"]
    half = r["T14_h"] / 48
    side = max(0.3, half)
    r["n_before"] = int(np.sum((t < r["T"] - half) & (t > r["T"] - half - side)))
    r["n_after"] = int(np.sum((t > r["T"] + half) & (t < r["T"] + half + side)))
    inn = np.sort(t[np.abs(t - r["T"]) < half])
    r["max_gap_in"] = float(np.max(np.diff(inn))) if len(inn) > 2 else 99.0
    r["cad"] = float(np.median(np.diff(t)))
    r["sector"] = int(sec)
    return r


def screen(which, toi):
    out = common.PILOT / "results" / which
    j = json.loads((out / f"{toi}.json").read_text())
    if j.get("status") != "ok":
        return None
    tic, D = j["TIC"], j["dur_h"] / 24
    lc = common.ap.load_lightcurves(tic, products=common.PRODUCTS)
    lc, _ = common.mask_siblings(lc, common.sibling_ephemerides(tic, toi))
    ev = fit(lc, j["event_T"], D)
    res = dict(TOI=toi, TIC=tic, event=None, repeats=[])
    if ev is None:
        res["note"] = "event fit failed"
        return res, j, lc, ev, []
    res["event"] = {k: v for k, v in ev.items() if not k.startswith("_")}
    fits = []
    for rp in j.get("repeats") or []:
        r = fit(lc, rp["T"], D)
        row = dict(T_grid=rp["T"], sector=rp["sector"], grid_depth_ppm=rp["depth_ppm"], grid_snr=rp["snr"])
        if r is None:
            row.update(passed=False, why="fit failed")
        else:
            min_side = min(max(0.3, r["T14_h"] / 48) / r["cad"] * 0.5, 30)
            why = []
            if not (r["sigma_T"] < D / 4 and r["sigma_depth"] < r["depth"] / 3):
                why.append("unconstrained fit")
            if not (0.5 <= r["depth"] / ev["depth"] <= 2):
                why.append(f"depth ratio {r['depth'] / ev['depth']:.2f}")
            if not (0.67 <= r["T14_h"] / ev["T14_h"] <= 1.5):
                why.append(f"T14 ratio {r['T14_h'] / ev['T14_h']:.2f}")
            if r["n_before"] < min_side or r["n_after"] < min_side:
                why.append(f"one-sided ({r['n_before']}/{r['n_after']} pts)")
            if r["max_gap_in"] > 0.1:
                why.append(f"gap {r['max_gap_in']:.2f} d in transit")
            row.update(passed=not why, why="; ".join(why), T_fit=r["T"], sigma_T=r["sigma_T"],
                       depth_ppm=round(r["depth"] * 1e6), sigma_depth_ppm=round(r["sigma_depth"] * 1e6),
                       T14_h=round(r["T14_h"], 2))
            fits.append((row, r))
        res["repeats"].append(row)
    passing = [dict(T=x["T_fit"], sector=x["sector"]) for x in res["repeats"] if x.get("passed")]
    g = pd.read_csv(out / f"{tic}_{toi}_grid.csv.gz")
    g["near_event"] = (g.T_pred - j["event_T"]).abs() < D
    pers = common.k2lib.robust_periods(g, ev["T"], passing, D)
    res["periods_screened"] = pers
    res["n_passing"] = len(passing)
    return res, j, lc, ev, fits


def plot(res, j, ev, fits, path):
    panels = [("event", ev, None)] + [(f"s{row['sector']} {'PASS' if row['passed'] else 'fail'}", r, row)
                                     for row, r in fits]
    k = len(panels)
    ncol = min(k, 4)
    nrow = int(np.ceil(k / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 2.7 * nrow), squeeze=False)
    for ax, (lab, r, row) in zip(axs.ravel(), panels):
        xh = (r["_tt"] - r["T"]) * 24
        ax.plot(xh, (r["_ff"] - 1) * 1e6, ".", ms=1, color="0.75")
        edges = np.arange(xh.min(), xh.max() + 0.5, 0.5)
        idx = np.digitize(xh, edges)
        u = [i for i in np.unique(idx) if (idx == i).sum() >= 2]
        ax.plot([xh[idx == i].mean() for i in u], [(r["_ff"][idx == i].mean() - 1) * 1e6 for i in u], "o", ms=2.5,
                color="b")
        ax.plot(xh, (r["_model"] - 1) * 1e6, "r-", lw=1)
        lim = 4 * ev["depth"] * 1e6
        ax.set_ylim(-lim, lim / 2)
        ax.set_title(f"{lab} T={r['T']:.2f}\n{r['depth'] * 1e6:.0f} ppm T14 {r['T14_h']:.1f} h"
                     + (f"\n{row['why'][:40]}" if row and row["why"] else ""), fontsize=7)
    for ax in axs.ravel()[k:]:
        ax.axis("off")
    fig.suptitle(f"TOI-{res['TOI']} (TIC {res['TIC']}); TOI {j['depth_ppm']:.0f} ppm, {j['dur_h']:.1f} h", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=90)
    plt.close(fig)


def main():
    which = sys.argv[1]
    out = common.PILOT / "results" / which
    tois = sys.argv[2:] or [f.stem for f in out.glob("*.json") if not f.stem.startswith(("vet_", "screen_", "diffs_"))]
    rows = []
    for toi in sorted(tois):
        s = screen(which, toi)
        if s is None:
            continue
        res, j, lc, ev, fits = s
        (out / f"screen_{toi}.json").write_text(json.dumps(res, indent=1, default=float))
        if ev is not None:
            plot(res, j, ev, fits, common.PILOT / "plots" / f"screen_{toi}.png")
        pers = res.get("periods_screened") or []
        top = max([p["n_support"] for p in pers], default=0)
        best = sorted(p["P"] for p in pers if p["n_support"] == top)
        rows.append(dict(TOI=toi, event_depth_ppm=round(ev["depth"] * 1e6) if ev else None,
                         event_T14_h=round(ev["T14_h"], 2) if ev else None,
                         n_repeats=len(res["repeats"]), n_passing=res.get("n_passing", 0),
                         passing=" ".join(f"{x['T_fit']:.3f}(s{x['sector']},{x['depth_ppm']}ppm,{x['T14_h']}h)"
                                          for x in res["repeats"] if x.get("passed")),
                         failing=" | ".join(f"{x['T_grid']:.2f}(s{x['sector']}): {x['why']}"
                                            for x in res["repeats"] if not x.get("passed")),
                         max_support=top, n_top=len(best),
                         top_periods=" ".join(f"{p:.3f}" for p in best[:15]) + (" ..." if len(best) > 15 else "")))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(common.PILOT / f"screen_{which}.csv", index=False)


if __name__ == "__main__":
    main()
