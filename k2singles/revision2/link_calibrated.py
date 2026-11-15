"""
AAS82225 revision: POST-HOC linked-period search, calibrated on the inverted control (PLAN.md, "Things decided after
first results", 2026-11-15). The pre-registered result (search_periodic.py) stays the primary one.

Why: the pre-registered single-event threshold came out at SNR 53 (K2) and 60 (TESS). The control events that set it
are large BRIGHTENINGS (bad cadences, glints; up to 150x the stellar flux) that become dips when the flux is
inverted and have no dimming counterpart, so the inverted control overstates single false dips. A period needs
several events of the same depth and duration on one ephemeris, which artefacts rarely give; so calibrate that.

Events: first-stage events (SNR >= 8, search_any rules) that are eligible (n_reductions >= 2 or n_alt_tested = 0)
and not common-mode; the same for control events (common-mode computed among control events).
Linking: for every pair of events with depth ratio 0.5-2 and duration ratio 0.5-2, every P = dT/n >= 1 d. Support =
the number of events with depth and duration within 0.5-2x of the pair's mean and within max(0.1 d, D/4) of the
ephemeris. A period is kept unless a clean predicted epoch is flat (box depth < 0.5 x pair depth - 3 sigma), as in
search_periodic.link. Per star: max_support over kept periods.
Calibration: N* = the smallest support with at most 0.05 control stars per star at or above it (pooled over stars).
Stars with max_support >= N* are detections; their periods with support >= N* are reported.
Output: k2singles/revision2/data/link_calibrated.json, link_calibrated.csv
Usage:  python k2singles/revision2/link_calibrated.py [--workers N]
"""
import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import search_periodic as sp  # noqa: E402

DATA = HERE / "data"
FA_PER_STAR = 0.05


def cm_flags(reps, key):
    rows = [(r["EPIC"], e["sector"], e["T"]) for r in reps for e in r[key]]
    df = pd.DataFrame(rows, columns=["EPIC", "sector", "T"])
    flag = set()
    for _, g in df.groupby("sector"):
        for epic, T in zip(g.EPIC, g["T"]):
            near = g[(g.EPIC != epic) & (np.abs(g["T"] - T) < sp.CM_WINDOW)]
            if near.EPIC.nunique() >= 2:
                flag.add((epic, round(T, 5)))
    return flag


def usable(r, key, cm):
    return [e for e in r[key] if sp.eligible(e) and (r["EPIC"], round(e["T"], 5)) not in cm]


def link_star(evs, lc, grids):
    """Pairwise linking. Returns list of kept periods with support >= 2 (P, support, T list)."""
    if len(evs) < 2:
        return []
    T = np.array([e["T_fit"] for e in evs])
    d = np.array([e["depth_fit"] for e in evs])
    D = np.array([e["D_fit"] for e in evs])
    snr = np.array([e["snr"] for e in evs])
    kept = {}
    for i in range(len(evs)):
        for j in range(i + 1, len(evs)):
            if not (0.5 <= d[i] / d[j] <= 2 and 0.5 <= D[i] / D[j] <= 2):
                continue
            dT = abs(T[j] - T[i])
            if dT < sp.P_MIN:
                continue
            dm, Dm = (d[i] + d[j]) / 2, (D[i] + D[j]) / 2
            like = (d / dm >= 0.5) & (d / dm <= 2) & (D / Dm >= 0.5) & (D / Dm <= 2)
            tol = max(0.1, Dm / 4)
            t0 = min(T[i], T[j])
            ns = np.arange(1, int(dT / sp.P_MIN) + 1)
            Ps = dT / ns
            ph = (T[like][None, :] - t0) / Ps[:, None]
            on = np.abs(ph - np.round(ph)) * Ps[:, None] <= tol
            sup = on.sum(axis=1)
            key_D = round(Dm * 24 * 4) / 4 / 24
            if key_D not in grids:
                grids[key_D] = sp.box_grid(lc, key_D) if lc is not None else (np.array([]),) * 4
            Tg, dg, sg, cg = grids[key_D]
            d_err = dm / float(np.median(snr[[i, j]]))
            for P, s, onrow in zip(Ps, sup, on):
                if s < 2:
                    continue
                k = round(P, 3)
                if k in kept and kept[k]["support"] >= s:
                    continue
                excl = False
                if len(Tg):
                    k0, k1 = np.ceil((Tg[0] - t0) / P), np.floor((Tg[-1] - t0) / P)
                    Tk = t0 + np.arange(k0, k1 + 1) * P
                    if len(Tk):
                        Tk = Tk[np.min(np.abs(Tk[:, None] - T[like][onrow][None, :]), axis=1) > tol]
                    if len(Tk):
                        idx = np.clip(np.searchsorted(Tg, Tk), 1, len(Tg) - 1)
                        idx = np.where(np.abs(Tg[idx - 1] - Tk) < np.abs(Tg[idx] - Tk), idx - 1, idx)
                        near = np.abs(Tg[idx] - Tk) <= key_D / 8
                        clean = near & (cg[idx] >= 0.8)
                        excl = bool(np.any(clean & (dg[idx] < 0.5 * dm - 3 * np.hypot(sg[idx], d_err))))
                if not excl:
                    kept[k] = dict(P=float(P), support=int(s), T=[round(x, 4) for x in sorted(T[like][onrow])],
                                   depth_ppm=round(dm * 1e6, 1), D_h=round(Dm * 24, 2))
    return sorted(kept.values(), key=lambda x: (-x["support"], -x["P"]))


def run(r, cm_real, cm_ctl):
    lc = sp.primary_lcs(r)
    grids = {}
    real = link_star(usable(r, "events", cm_real), lc, grids)
    # control: same linking on the inverted light curve's events; flat-epoch test on the inverted flux
    lci = None if lc is None else dict(lc, f=2.0 - lc["f"])
    ctl = link_star(usable(r, "control_events", cm_ctl), lci, {})
    return dict(EPIC=r["EPIC"], TIC=r["TIC"], real=real[:200], control=ctl[:50],
                max_support=max((p["support"] for p in real), default=0),
                max_support_control=max((p["support"] for p in ctl), default=0))


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--workers", type=int, default=12)
    args = a.parse_args()
    reps = [json.loads(f.read_text()) for f in sorted((DATA / "per_star").glob("*.json"))]
    reps = [r for r in reps if r.get("status") == "ok"]
    cm_real, cm_ctl = cm_flags(reps, "events"), cm_flags(reps, "control_events")
    res = []
    with ProcessPoolExecutor(args.workers, max_tasks_per_child=4) as ex:
        futs = [ex.submit(run, r, cm_real, cm_ctl) for r in reps]
        for k, fu in enumerate(as_completed(futs), 1):
            x = fu.result()
            res.append(x)
            print(f"[{k}/{len(reps)}] EPIC {x['EPIC']} max support {x['max_support']} (control {x['max_support_control']})", flush=True)
    n = len(res)
    ctl = np.array([x["max_support_control"] for x in res])
    allowed = int(np.floor(FA_PER_STAR * n))
    N_star = next(N for N in range(2, 50) if (ctl >= N).sum() <= allowed)
    rows = []
    for x in res:
        det = x["max_support"] >= N_star
        best = [p for p in x["real"] if p["support"] >= N_star]
        rows.append(dict(EPIC=x["EPIC"], TIC=x["TIC"], max_support=x["max_support"],
                         max_support_control=x["max_support_control"], detected=det, n_periods=len(best),
                         periods=" ".join(f"{p['P']:.4f}({p['support']})" for p in best[:12]),
                         depth_ppm=best[0]["depth_ppm"] if best else None, D_h=best[0]["D_h"] if best else None))
    out = dict(N_star=N_star, n_stars=n, allowed_control_stars=allowed,
               control_max_support_hist={int(k): int(v) for k, v in zip(*np.unique(ctl, return_counts=True))},
               real_max_support_hist={int(k): int(v) for k, v in zip(*np.unique([x["max_support"] for x in res], return_counts=True))},
               stars=res)
    (DATA / "link_calibrated.json").write_text(json.dumps(out, indent=1, default=float))
    pd.DataFrame(rows).sort_values(["detected", "max_support"], ascending=False).to_csv(DATA / "link_calibrated.csv", index=False)
    print(f"N* = {N_star} (control stars with max support >= N*: {(ctl >= N_star).sum()} of {n}, allowed {allowed})")
    print("control max-support histogram:", out["control_max_support_hist"])
    print("real max-support histogram:   ", out["real_max_support_hist"])
    print(f"detected stars: {sum(r['detected'] for r in rows)}")


if __name__ == "__main__":
    main()
