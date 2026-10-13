"""
Ephemeris check for stars whose confirmed events (catalogued event, known repeats and extra
events of categories a-c, fine-fit centres) may share one period.

Events are folded at a trial period P_guess (e.g. from the spacing of two similar events) and
split into phase clusters (gaps > 0.3 d); each cluster gets cycle numbers
round((T - T_first)/P_guess) and a linear ephemeris is fitted. Reported per family: P, T0,
rms and max O-C; a family is 'ok' if max |O-C| < 0.05 d and no two events share a cycle. The
phase of the second family is given relative to the first.

Usage: python singles/revision_search/period_check.py TIC:P_guess[+] [...]
       ('+' also uses extra events rejected only for being in a single reduction)
Output: singles/revision_search/data/period_check.csv (+ printed)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def events_for(tic, ex, include_single_reduction=False):
    j = json.loads((HERE / "data" / "per_star" / f"{tic}.json").read_text())
    rows = []
    Dc = j["dur_h"] / 24
    for ev in j["events"]:
        known = any(abs(ev["T"] - Tk) < (ev["D"] + Dc) / 2 + 0.1 for Tk in [j["event_T"]] + j["known_repeats"])
        e = ex[(ex.TIC == tic) & ((ex.T_btjd - ev["T"]).abs() < 1e-3)]
        ok = known or (len(e) and (e.category.iloc[0] != "d"
                                   or (include_single_reduction and e.category_reason.iloc[0] == "single-reduction")))
        if ok:
            rows.append((ev.get("T_fit", ev["T"]), ev.get("depth_fit", ev["depth"]), ev.get("D_fit", ev["D"]),
                         ev["product"], "known" if known else e.category.iloc[0]))
    return pd.DataFrame(rows, columns=["T", "depth", "D", "product", "kind"]).sort_values("T")


def families(ev, Pg, gap=0.3):
    """Phase clusters at the trial period: events folded (in days, mod P_guess) and split at
    gaps > 0.3 d on the circle."""
    T = ev["T"].to_numpy()
    x = (T - T[0]) % Pg
    o = np.argsort(x)
    xs = x[o]
    d = np.diff(np.concatenate([xs, [xs[0] + Pg]]))
    if (d > gap).sum() == 0:
        return [ev]
    start = (d.argmax() + 1) % len(xs)  # rotate so the largest gap is at the end
    xs_r = np.roll(o, -start)
    xr = (x[xs_r] - x[xs_r[0]]) % Pg
    cut = np.where(np.diff(xr) > gap)[0]
    groups = np.split(xs_r, cut + 1)
    return sorted([ev.iloc[g] for g in groups], key=lambda f: -len(f))


def fit(fam, Pg):
    T = fam["T"].to_numpy()
    n = np.round((T - T[0]) / Pg)
    if len(np.unique(n)) < len(n):
        return dict(ok=False, why="cycle collision", n=len(T))
    if len(T) == 1:
        return dict(ok=True, P=Pg, T0=float(T[0]), n=1)
    A = np.column_stack([np.ones_like(n), n])
    (T0, P), *_ = np.linalg.lstsq(A, T, rcond=None)
    oc = T - (T0 + n * P)
    return dict(ok=bool(np.max(np.abs(oc)) < 0.05), P=float(P), T0=float(T0), n=len(T),
                rms_oc_min=float(np.sqrt(np.mean(oc ** 2)) * 1440), max_oc_min=float(np.max(np.abs(oc)) * 1440))


def main():
    ex = pd.read_csv(HERE / "extra_events_auto.csv")
    out = []
    for arg in sys.argv[1:]:
        flag = arg.endswith("+")
        tic, Pg = arg.rstrip("+").split(":")
        tic, Pg = int(tic), float(Pg)
        ev = events_for(tic, ex, include_single_reduction=flag)
        fams = families(ev, Pg)
        res = [fit(f, Pg) for f in fams]
        P1, T1 = res[0].get("P", Pg), res[0].get("T0", ev["T"].iloc[0])
        desc = []
        for f, r in zip(fams, res):
            ph = ((r.get("T0", f["T"].iloc[0]) - T1) / P1) % 1.0
            desc.append(dict(n=len(f), depth_ppm=round(float(f.depth.median()) * 1e6), dur_h=round(float(f.D.median()) * 24, 1),
                             phase=round(ph, 4), ok=r["ok"], P=round(r.get("P", np.nan), 5),
                             max_oc_min=round(r.get("max_oc_min", np.nan), 1), T=" ".join(f"{t:.3f}" for t in f["T"])))
        row = dict(TIC=tic, P_guess=Pg, includes_single_reduction=flag, n_events=len(ev), n_clusters=len(fams),
                   n_in_clusters_ge2=int(sum(len(f) for f in fams if len(f) >= 2)), P_fit=round(P1, 5),
                   T0_fit=round(T1, 4), clusters=json.dumps(desc))
        print(row)
        out.append(row)
    pd.DataFrame(out).to_csv(HERE / "data" / "period_check.csv", index=False)


if __name__ == "__main__":
    main()
