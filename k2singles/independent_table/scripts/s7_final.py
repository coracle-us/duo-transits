"""Final per-star synthesis: refit every vetted event, derive periods with the exclusion
machinery, check P/2..P/4 exclusions, odd/even depths, phase-0.5 secondary, and write
results.csv + data/final_details.json.  Event lists below were vetted by eye (plots/)."""
import sys, os, json, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
import lib, s5_periods as s5
from fit_events import fit_all
from epoch_summary import summary
from targets import T

K2K = {"K2pdc", "K2ev", "K2sff"}
TK = {"spoc", "tspoc", "qlp", "own", "tglc", "eleanor"}

# vetted events: (approx time, label). label 'P' = primary (used for period), 'S' = secondary,
# 'p' = partial / supporting primary (not used in the period fit)
EVENTS = {
    248811085: [(948.38, "P"), (2573.633, "P"), (3263.135, "P")],
    246696483: [(855.59, "P"), (2495.796, "P"), (3224.776, "P")],
    247967714: [(845.078, "P"), (2515.718, "P"), (3258.224, "P")],
    246089124: [(774.470, "P"), (3215.229, "P"), (3806.929, "P")],
    211892898: [(189.33, "P"), (1131.94, "P"), (1276.96, "P")],
    211503363: [(197.518, "P"), (2532.72, "P"), (3270.14, "P")],
    211390677: [(151.786, "P"), (1257.637, "P"), (2547.823, "P"), (3284.95, "p")],
    211411112: [(178.953, "P"), (1117.876, "P"), (2557.556, "P")],
    211489484: [(192.621, "P"), (1261.983, "P"), (1508.757, "P"), (2249.087, "P"), (2578.118, "P"),
                (3236.180, "P"), (208.605, "S"), (1277.968, "S"), (3252.17, "S")],
    202071902: [(-208.427, "P"), (2517.165, "P"), (3250.097, "P"), (3273.002, "p")],
    202072917: [(-213.370, "P"), (1489.165, "P"), (2496.787, "P"), (2531.533, "P"), (4060.341, "P")],
    230887315: [(660.31, "P")],
    212152316: [(207.194, "P")],
    202126877: [(-211.560, "P"), (2488.584, "P"), (2510.368, "P"), (2532.147, "P"), (3250.695, "P"),
                (3272.484, "P"), (4056.378, "P"), (2495.826, "S"), (2517.592, "S"), (3236.092, "S"),
                (3257.846, "S"), (3279.610, "S")],
}


def combine(e, tc, durg=None, depg=None):
    rows = fit_all(e, tc, durg=durg, depg=depg, kinds=K2K | TK, quiet=True)
    rows = [(s, f) for s, f in rows if np.isfinite(f["e_t0"]) and f["e_t0"] < 0.1 and f["depth"] > 0]
    if not rows:
        return None
    t0s = np.array([f["t0"] for _, f in rows]); ets = np.array([f["e_t0"] for _, f in rows])
    t0 = np.median(t0s)
    et = max(np.median(ets), np.std(t0s) if len(t0s) > 1 else 0)
    best = min(rows, key=lambda x: x[1]["e_depth"])
    return dict(t0=t0, e_t0=et, n_prod=len(rows), best_src=best[0], depth=best[1]["depth"],
                e_depth=best[1]["e_depth"], T14_h=best[1]["T14"] * 24,
                per_prod={s: (round(f["t0"], 5), round(f["depth"] * 1e6)) for s, f in rows})


def run(only=None, extra=None):
    ev = dict(EVENTS)
    if extra: ev.update(extra)
    out = {}
    for _, r in T.iterrows():
        e = int(r.epic)
        if only and e not in only: continue
        if e not in ev: continue
        res = dict(epic=e, tic=int(r.tic), events=[])
        for tc, lab in ev[e]:
            c = combine(e, tc)
            if c is None:
                res["events"].append(dict(t_guess=tc, label=lab, fit="failed")); continue
            c["label"] = lab
            res["events"].append(c)
        prim = [(x["t0"], x["e_t0"]) for x in res["events"] if x.get("label") == "P" and "t0" in x]
        if len(prim) >= 2:
            a = s5.analyse(e, prim, Pmin=1.0, verbose=False)
            ok = [c for c in a["cands"] if not c["excluded"]]
            res["allowed"] = [(c["P"], c["eP"], c["n"]) for c in ok]
            res["n_commensurate"] = len(a["cands"])
            if ok:
                Pb = max(c["P"] for c in ok)
                T0 = prim[0][0]
                res["submultiples"] = {}
                for m in (2, 3, 4):
                    h = s5.why_excluded(e, T0, Pb / m, show=False)
                    res["submultiples"][m] = [(k, round(tc, 3), p[:2]) for k, tc, p in h[:3]]
                # epochs with data for the longest allowed P and at phase 0.5
                res["epochs"] = [(k, round(tc, 3), src, None if not np.isfinite(d) else round(d * 1e6), None if not np.isfinite(s) else round(s * 1e6), ex, det)
                                 for k, tc, src, d, s, cov, ex, det in summary(e, T0, Pb, 0.0, show=False)]
                res["phase05"] = [(k, round(tc, 3), src, None if not np.isfinite(d) else round(d * 1e6), None if not np.isfinite(s) else round(s * 1e6), ex, det)
                                  for k, tc, src, d, s, cov, ex, det in summary(e, T0, Pb, 0.5, show=False)]
        elif len(prim) == 1:
            a = s5.analyse(e, prim, Pmin=1.0, verbose=False)
            res["allowed_windows"] = [(round(x, 3), round(y, 3)) for x, y in a["allowed"]]
        out[e] = res
        print(json.dumps({k: v for k, v in res.items() if k in ("epic", "allowed", "submultiples")}, default=str))
    return out


if __name__ == "__main__":
    o = run([int(a) for a in sys.argv[1:]] or None)
    fn = os.path.join(lib.ROOT, "data", "final_details.json")
    old = json.load(open(fn)) if os.path.exists(fn) else {}
    old.update({str(k): v for k, v in o.items()})
    json.dump(old, open(fn, "w"), indent=1, default=str)
