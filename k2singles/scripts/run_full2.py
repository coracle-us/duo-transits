"""
K2 singles x TESS, full run v2 (after the Fable review). It uses k2lib (explicit TESS + K2 loading,
robust period inference, stricter exclusion and K2-repeat rules) and one task per worker process.
The selection is unchanged from run_full.py: n_tess_observed >= 1, tess_snr_per_transit >= 7,
depth < 1e6 ppm. Each EVENT is one job, keyed EPIC_<i>, so a star with two catalogued events keeps
both. Catalogue values are used as-is for every target, including EPIC 211924561.
TIC-match quality is flagged: |Tmag - Kp| > 1 mag, or more than one TIC within 3".
Output: k2singles/results/full2/<EPIC>_<i>.json (+ plots, diff images), summary.csv. Resumable.
Usage:  python k2singles/scripts/run_full2.py [n_workers]
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
sys.path.insert(0, str(HERE))
OUT = ROOT / "k2singles" / "results" / "full2"


def select():
    d = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv")
    d["event_idx"] = d.groupby("EPIC").cumcount()
    s = d[(d.n_tess_observed > 0) & (d.tess_snr_per_transit >= 7) & (d.depth_ppm < 1e6)].copy()
    s["match_flag"] = ((s.Tmag - s.Kp).abs() > 1) | (s.n_tic_3as > 1)
    return s


def work(row):
    import k2lib
    import pilot
    pilot.OWN_MIN_SECTOR = 80
    key = f"{row['EPIC']}_{row['event_idx']}"
    f = OUT / f"{key}.json"
    if f.exists():
        return key, "exists"
    rep = dict(EPIC=row["EPIC"], event_idx=row["event_idx"], TIC=row["TIC"], Tmag=row["Tmag"], Kp=row["Kp"],
               match_flag=bool(row["match_flag"]), depth_ppm=row["depth_ppm"], dur_h=row["dur_h"],
               epoch_btjd=row["epoch_btjd"], k2_campaign=row["campaign"], catalog_comment=row["comment"],
               TOI=row["TOI"])
    try:
        tic, epic = int(row["TIC"]), int(row["EPIC"])
        D, delta, epoch = row["dur_h"] / 24, row["depth_ppm"] * 1e-6, float(row["epoch_btjd"])
        import fetch_data as fd
        if not (fd.LC_DIR / str(tic) / "search_results.csv").exists():
            fd._fetch_one(0, tic)
        rep["own_photometry_sectors"] = pilot.add_own_photometry(tic, row["RA"], row["Dec"])
        lc, k2 = k2lib.load_combined(tic, epic)
        rep["k2_campaigns_loaded"] = [] if k2 is None else sorted({int(s) - 1000 for s in k2.sector})
        rep["tess_sectors_loaded"] = sorted({int(s) for s in lc.sector if s < 1000})
        g = k2lib.grid(lc, D, delta, epoch)
        reps = k2lib.find_repeats(g, delta, D)
        rep["repeats"] = reps
        rep["n_repeats"] = len(reps)
        rep["periods"] = k2lib.robust_periods(g, epoch, reps, D)
        rep["n_periods"] = len(rep["periods"])
        rep["max_support"] = max([p["n_support"] for p in rep["periods"]], default=0)
        rep["frac_logP_excluded"] = k2lib.frac_excluded(g, epoch, D)
        rep["n_windows"] = int(len(g))
        rep["n_excluded_windows"] = int((g.verdict == "excluded").sum())
        tess_reps = [x for x in reps if x["sector"] < 1000]
        rep["diff_images"] = pilot.compute_diffs(tic, pd.Series(row), lc[lc.sector < 1000], epoch, tess_reps, D, OUT)
        pilot.plot_windows(tic, lc, epoch, reps, D, delta, OUT)
        rep["status"] = "ok"
    except Exception:  # noqa: BLE001
        rep.update(status="error", error=traceback.format_exc()[-600:])
    f.write_text(json.dumps(rep, indent=1, default=float))
    return key, (f"{rep.get('status')} repeats {rep.get('n_repeats')} periods "
                 f"{[(round(p['P'], 2), p['n_support']) for p in (rep.get('periods') or [])[:4]]}")


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    OUT.mkdir(parents=True, exist_ok=True)
    s = select()
    rows = [dict(EPIC=int(r.EPIC), event_idx=int(r.event_idx), TIC=int(r.TIC), Tmag=float(r.Tmag), Kp=float(r.Kp),
                 match_flag=bool(r.match_flag), RA=float(r.ra), Dec=float(r.dec), depth_ppm=float(r.depth_ppm),
                 dur_h=float(r.dur_h), epoch_btjd=float(r.epoch_btjd), campaign=str(r.campaign),
                 comment=str(r.comment), TOI=str(r.TOI) if pd.notna(r.TOI) else "", CTOI=None)
            for r in s.itertuples()]
    todo = [r for r in rows if not (OUT / f"{r['EPIC']}_{r['event_idx']}.json").exists()]
    print(f"events: {len(rows)}; to do: {len(todo)}; TIC-match flags: {int(s.match_flag.sum())}", flush=True)
    with ProcessPoolExecutor(max_workers=n_workers, max_tasks_per_child=1) as ex:
        futs = [ex.submit(work, r) for r in todo]
        for i, fu in enumerate(as_completed(futs), 1):
            key, msg = fu.result()
            print(f"[{i}/{len(todo)}] {key}: {msg}", flush=True)
    summ = []
    for r in rows:
        fn = OUT / f"{r['EPIC']}_{r['event_idx']}.json"
        if fn.exists():
            j = json.loads(fn.read_text())
            reps = j.get("repeats") or []
            summ.append(dict(EPIC=r["EPIC"], event_idx=r["event_idx"], TIC=r["TIC"], TOI=r["TOI"], Tmag=r["Tmag"],
                             match_flag=r["match_flag"], depth_ppm=r["depth_ppm"], dur_h=r["dur_h"], comment=r["comment"],
                             status=j.get("status"), k2_campaigns=" ".join(map(str, j.get("k2_campaigns_loaded") or [])),
                             n_repeats=len(reps), n_repeats_k2=sum(1 for x in reps if x["sector"] >= 1000),
                             n_periods=j.get("n_periods"), max_support=j.get("max_support"),
                             frac_logP_excluded=j.get("frac_logP_excluded"),
                             periods=" ".join(f"{p['P']:.3f}({p['n_support']})" for p in (j.get("periods") or [])[:6])))
    pd.DataFrame(summ).to_csv(OUT / "summary.csv", index=False)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
