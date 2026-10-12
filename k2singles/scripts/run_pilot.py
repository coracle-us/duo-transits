"""
K2 singles x TESS pilot: run the single-transit pipeline on the 50 planet-like K2 single transits
with TESS coverage (selection pre-registered in k2singles/README.md).

Per target:
  - singles/scripts/pilot.run_target with the K2 epoch (BTJD) as the reference event. It runs the
    repeat search over all TESS data (pipeline products, plus our own TESScut photometry for sectors
    >= 80 that lack them), flat-window period constraints, periods implied by any repeat, and
    difference images at repeats.
  - K2 campaign coverage (MAST observation t_min/t_max for the EPIC's K2 light curves). The K2
    light curve showed only one transit, so P must exceed the longer of (epoch - campaign start,
    campaign end - epoch), ignoring gaps. Candidate periods below that are dropped and logged.
Resumable: targets with an existing JSON are skipped.
Output: k2singles/results/pilot/<TIC>.json (+ window plots, diff images), summary.csv
Usage:  python k2singles/scripts/run_pilot.py [n_workers]
"""
import json
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "singles" / "scripts"))
sys.path.insert(0, str(ROOT / "scripts"))
OUT = ROOT / "k2singles" / "results" / "pilot"


def select():
    d = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv").drop_duplicates("EPIC")
    pl = d[(d.depth_ppm <= 20000) & (~d.comment.fillna("").str.contains("secondary|Deep", case=False))]
    return pl[(pl.n_tess_observed > 0) & (pl.tess_snr_per_transit >= 7)].copy()


def k2_coverage(epic, epoch_btjd):
    """Start/end (BTJD) of the K2 campaign that contains the epoch. Fields re-observed in later
    campaigns (C5/C16/C18, C6/C17, ...) have several intervals years apart; only the containing
    one bounds P without gaps. (Bug fixed 2026-10-12: the first version took min/max over all
    campaigns, giving spurious P_min ~1000 d.)"""
    import lightkurve as lk
    try:
        t = lk.search_lightcurve(f"EPIC {epic}", mission="K2").table
    except Exception:  # noqa: BLE001
        return None, None, []
    iv = sorted({(round(float(a) + 2400000.5 - 2457000, 2), round(float(b) + 2400000.5 - 2457000, 2))
                 for a, b in zip(t["t_min"], t["t_max"])})
    cont = [(a, b) for a, b in iv if a - 0.5 <= epoch_btjd <= b + 0.5]
    if not cont:
        return None, None, iv
    a, b = min(x[0] for x in cont), max(x[1] for x in cont)
    return a, b, iv


def work(row_dict):
    import pilot
    pilot.OWN_MIN_SECTOR = 80
    row = pd.Series(row_dict)
    tic = int(row.TIC)
    f = OUT / f"{tic}.json"
    if f.exists():
        return tic, "exists"
    try:
        rep = pilot.run_target(row, OUT, True)
        tmin, tmax, iv = k2_coverage(int(row.EPIC), float(row.epoch_btjd))
        rep["EPIC"] = int(row.EPIC)
        rep["k2_campaign"] = str(row.campaign)
        rep["k2_tmin_btjd"], rep["k2_tmax_btjd"] = tmin, tmax
        rep["k2_intervals_btjd"] = iv
        if tmin is not None:
            pmin = max(float(row.epoch_btjd) - tmin, tmax - float(row.epoch_btjd))
            rep["P_min_from_K2"] = pmin
            cp = rep.get("candidate_periods") or []
            rep["candidate_periods_all"] = cp
            rep["candidate_periods"] = [x for x in cp if x[1] >= pmin]
            rep["n_candidate_periods"] = len(rep["candidate_periods"])
    except Exception:  # noqa: BLE001
        rep = dict(TIC=tic, EPIC=int(row.EPIC), status="error", error=traceback.format_exc()[-600:])
    f.write_text(json.dumps(rep, indent=1, default=float))
    return tic, (f"repeats {rep.get('n_repeats')} periods {[round(p, 2) for _, p in (rep.get('candidate_periods') or [])][:5]} "
                 f"Pmin_K2 {rep.get('P_min_from_K2')} own {rep.get('own_photometry_sectors')}")


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    OUT.mkdir(parents=True, exist_ok=True)
    t = select()
    rows = []
    for r in t.itertuples():
        rows.append(dict(TIC=int(r.TIC), EPIC=int(r.EPIC), campaign=r.campaign, RA=float(r.ra), Dec=float(r.dec),
                         Tmag=float(r.Tmag), depth_ppm=float(r.depth_ppm), dur_h=float(r.dur_h),
                         epoch_btjd=float(r.epoch_btjd), CTOI=None))
    todo = [r for r in rows if not (OUT / f"{r['TIC']}.json").exists()]
    print(f"pilot targets: {len(rows)}; to do: {len(todo)}", flush=True)
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = [ex.submit(work, r) for r in todo]
        for i, fu in enumerate(as_completed(futs), 1):
            tic, msg = fu.result()
            print(f"[{i}/{len(todo)}] TIC {tic}: {msg}", flush=True)
    summ = []
    for r in rows:
        f = OUT / f"{r['TIC']}.json"
        if f.exists():
            j = json.loads(f.read_text())
            summ.append(dict(TIC=r["TIC"], EPIC=r["EPIC"], Tmag=r["Tmag"], depth_ppm=r["depth_ppm"], dur_h=r["dur_h"],
                             status=j.get("status", "ok"), n_sectors=j.get("n_sectors"), n_repeats=j.get("n_repeats"),
                             n_periods=j.get("n_candidate_periods"), P_min_K2=j.get("P_min_from_K2"),
                             frac_logP_excluded=j.get("frac_logP_excluded"),
                             periods=" ".join(f"{p:.2f}" for _, p in (j.get("candidate_periods") or [])[:6])))
    pd.DataFrame(summ).to_csv(OUT / "summary.csv", index=False)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
