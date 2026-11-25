"""
tess_wide step 3 (PLAN.md sections 2, 4, 5.4, 7): single-event search of every downloaded star.

For each star (tess_wide/data/lc_s<NNN>/<TIC>.npz):
  light curve  SAP_FLUX where QUALITY == 0 and finite; the first 0.5 d of each orbit (ORBITID) dropped;
               normalised by its median; segments split at gaps > 0.5 d (search_any.SEG_GAP)
  dumps        times with QUALITY bits 32, 1 or 4 (search_any.DUMP_BITS)
  search       search_any.search (imported, unchanged rules) at durations 2, 3, 5, 8, 13 h, on the flux
               (real) and on 2 - flux (control)
  flags        per event: in-event SAP_BKG change and SAP_X / SAP_Y shift against the flanks, in sigma
               (PLAN 5.4: flags only)
Injection mode (PLAN 7): one batman transit per star on 3,000 stars drawn with seed 20261125, injected into
SAP_FLUX before anything else; real sign only.

Output (outside git): tess_wide/data/s<NNN>_events_{real,ctrl,inj}.csv, s<NNN>_inj_truth.csv, s<NNN>_search.json
Usage:  python tess_wide/s3_search.py --sector 104 --mode real|inj [--workers 10]
"""
import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "singles" / "revision_search"))
import search_any as sa  # noqa: E402
sys.path.insert(0, str(HERE))
import inject  # noqa: E402

DATA = HERE / "data"
DURS = np.array([2, 3, 5, 8, 13]) / 24.0      # PLAN 4
ORBIT_TRIM = 0.5                              # d, PLAN 2
N_INJ, SEED = 3000, 20261125                  # PLAN 7
DEPTH_RANGE, DUR_RANGE_H, B_MAX = (300e-6, 15000e-6), (2.0, 13.0), 0.9


def load(npz):
    z = np.load(npz)
    t, f, q, orb = z["time"], z["sap"].astype(float), z["quality"], z["orbit"]
    aux = dict(bkg=z["bkg"].astype(float), x=z["x"].astype(float), y=z["y"].astype(float))
    meta = dict(camera=int(z["camera"]), ccd=int(z["ccd"]))
    return t, f, q, orb, aux, meta


def prepare(t, f, q, orb, aux):
    ok = np.isfinite(t)
    dumps = np.unique(np.round(t[ok & ((q & sa.DUMP_BITS) > 0)], 4))
    m = ok & np.isfinite(f) & (q == 0)
    for o in np.unique(orb[m]):
        mo = m & (orb == o)
        if mo.any():
            m &= ~(mo & (t - t[mo].min() < ORBIT_TRIM))
    t, f = t[m], f[m]
    a = {k: v[m] for k, v in aux.items()}
    if len(t) < 50:
        return None, dumps, a
    o = np.argsort(t)
    t, f, a = t[o], f[o] / np.nanmedian(f[o]), {k: v[o] for k, v in a.items()}
    seg = np.concatenate([[0], np.cumsum(np.diff(t) > sa.SEG_GAP)])
    return dict(t=t, f=f, seg=seg), dumps, a


def shift_sigma(t, v, T, D):
    """In-event mean of v minus the flank mean, in units of its error (flanks: D/2 + D/4 to D/2 + max(0.75, 2D))."""
    v = np.asarray(v, float)
    ok = np.isfinite(v)
    dt = np.abs(t - T)
    i = ok & (dt < 0.35 * D)
    w = ok & (dt > D / 2 + D / 4) & (dt < D / 2 + max(0.75, 2 * D))
    if i.sum() < 3 or w.sum() < 10:
        return np.nan, np.nan
    sd = np.std(v[w])
    d = v[i].mean() - v[w].mean()
    err = sd * np.sqrt(1 / i.sum() + 1 / w.sum())
    return (d / err if err > 0 else np.nan), d


def run_one(spec):
    tic, npz, sector, inj = spec
    t, f, q, orb, aux, meta = load(npz)
    truth = None
    if inj is not None:
        m = inject.model(t, inj)
        f = f * m
        truth = dict(depth_true=float(1 - m[np.isfinite(t)].min()))
    lc, dumps, a = prepare(t, f, q, orb, aux)
    out = []
    if lc is None:
        return tic, out, truth
    used = {sector: "QLP"}
    lc["seg"] = sector * 1000 + lc["seg"]
    for sign, kind in ((1.0, "inj" if inj is not None else "real"), (-1.0, "ctrl")):
        if inj is not None and sign < 0:
            continue
        evs, _ = sa.search(lc, used, dumps, DURS, sign)
        for e in evs:
            row = dict(tic=tic, kind=kind, T=e["T"], D=e["D"], depth=e["depth"], snr=e["snr"],
                       depth_med=e.get("depth_med"), beta=e.get("beta"), durs_detected=e.get("durs_detected"),
                       dump_in_window=e.get("dump_in_window"), **meta)
            for c in ("bkg", "x", "y"):
                s, d = shift_sigma(lc["t"], a[c], e["T"], e["D"])
                row[f"{c}_sig"], row[f"{c}_shift"] = s, d
            out.append(row)
    return tic, out, truth


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", type=int, default=104)
    ap.add_argument("--mode", choices=["real", "inj"], default="real")
    ap.add_argument("--workers", type=int, default=10)
    a = ap.parse_args()
    s = a.sector
    lcdir = DATA / f"lc_s{s:04d}"
    files = sorted(p for p in lcdir.glob("*.npz") if ".tmp" not in p.name)
    specs = []
    truth_rows = []
    if a.mode == "real":
        specs = [(int(p.stem), p, s, None) for p in files]
    else:
        rng = np.random.default_rng(SEED)
        pick = rng.choice(len(files), size=min(N_INJ, len(files)), replace=False)
        for i in pick:
            p = files[i]
            z = np.load(p)
            tt, qq = z["time"], z["quality"]
            g = tt[np.isfinite(tt) & (qq == 0) & np.isfinite(z["sap"])]
            inj = dict(depth_draw=float(np.exp(rng.uniform(*np.log(DEPTH_RANGE)))),
                       dur_h=float(np.exp(rng.uniform(*np.log(DUR_RANGE_H)))),
                       b=float(rng.uniform(0, B_MAX)), T=float(rng.uniform(g.min(), g.max())))
            specs.append((int(p.stem), p, s, inj))
            truth_rows.append(dict(tic=int(p.stem), **inj))
    t0 = time.time()
    rows, truths, n = [], {}, 0
    with ProcessPoolExecutor(a.workers) as ex:
        for tic, out, truth in ex.map(run_one, specs, chunksize=20):
            rows += out
            if truth:
                truths[tic] = truth
            n += 1
            if n % 5000 == 0:
                print(f"{n}/{len(specs)} stars, {time.time() - t0:.0f} s", flush=True)
    el = time.time() - t0
    ev = pd.DataFrame(rows)
    if a.mode == "real":
        for kind in ("real", "ctrl"):
            ev[ev.kind == kind].to_csv(DATA / f"s{s:04d}_events_{kind}.csv", index=False)
    else:
        ev.to_csv(DATA / f"s{s:04d}_events_inj.csv", index=False)
        tr = pd.DataFrame(truth_rows)
        tr["depth_true"] = tr.tic.map(lambda x: truths.get(x, {}).get("depth_true", np.nan))
        tr.to_csv(DATA / f"s{s:04d}_inj_truth.csv", index=False)
    js = DATA / f"s{s:04d}_search.json"
    summ = json.loads(js.read_text()) if js.exists() else {}
    summ[a.mode] = dict(stars=len(specs), seconds=round(el), workers=a.workers,
                        events={k: int((ev.kind == k).sum()) for k in ev.kind.unique()} if len(ev) else {})
    js.write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ[a.mode]))


if __name__ == "__main__":
    main()
