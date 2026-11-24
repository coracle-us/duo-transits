"""
The catalogue run (prereg_catalogue.yaml). Refuses to read any candidate's light curve until
catalogue/APPROVALS.md records Jon's step-2 approval (the freeze).

  search : steps 1 and 2 for every row of the definitions - reference, excluded periods, events.
  tables : step 4, after vetting - alias tables through the events both vetters accepted
           (the verdicts file), and the final table. Recomputes steps 1 and 2 from the same light
           curves, so nothing but the verdicts passes between the two stages.

--check DEFS runs both stages on other definitions (calibration stars only: it refuses any star of
catalogue/defs.csv) with every passing event taken as accepted; this is the driver check.

Inputs  : catalogue/defs.csv (build_defs.py), MAST light curves (cache),
          <out>/vetting/verdicts.csv for the tables stage (event_id, vetter_1, vetter_2)
Outputs : <out>/rows_search.csv, events.csv, excluded_periods.csv            (search)
          <out>/catalogue.csv, alias_tables.csv                               (tables)
Usage   : python catalogue/run_catalogue.py search|tables [--out catalogue/results/v0] [--workers 12]
          python catalogue/run_catalogue.py tables --check DEFS --out DIR
"""
import argparse
import multiprocessing as mp
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import catrun  # noqa: E402

CACHE = catrun.dump3.CACHE
EV_COLS = ["event_id", "row_id", "sector", "author", "tc", "dur_d", "depth", "depth_err", "snr", "frac",
           "depth_ratio", "dur_ratio", "anti_snr", "local_anti_snr", "sec_box_snr", "edge_dist_d",
           "noise_ppm", "passes", "at_t2_cat", "at_other_row"]


def frozen():
    p = HERE / "APPROVALS.md"
    return p.exists() and "Step 2 of 2" in p.read_text()


def star(args):
    defs, cache, stage, accepted = args
    tic = int(defs[0]["tic"])
    try:
        sectors, _ = catrun.mc.load_star(tic, Path(cache))
        if not sectors:
            return dict(rows=[dict(row_id=d["row_id"], status="no light curve") for d in defs])
        rows, events, excl, tables = [], [], [], []
        times = {d["row_id"]: [t for t in (d["t_cat"], d.get("t2_cat", np.nan)) if np.isfinite(t)] for d in defs}
        for d in defs:
            other = [t for r, ts in times.items() if r != d["row_id"] for t in ts]
            row, ev, ctx = catrun.search_row(d, sectors, other)
            if ctx is not None:
                excl += [dict(row_id=d["row_id"], p_lo_d=a, p_hi_d=b) for a, b in ctx["excluded_runs"]]
                if stage == "tables":
                    ok = {e["event_id"] for e in ev if e["passes"]} if accepted is None else accepted
                    cols, tabs = catrun.finish_row(d, ctx, ev, ok)
                    row.update(cols)
                    for kind, t1, t2, rr in tabs:
                        tables += [dict(row_id=d["row_id"], table=kind, t1=t1, t2=t2, n=n, period_d=P,
                                        status=st, fits_third=f) for n, P, st, f in rr]
            rows.append(row)
            events += [{k: e.get(k) for k in EV_COLS} for e in ev]
        return dict(rows=rows, events=events, excl=excl, tables=tables)
    except Exception as e:
        return dict(rows=[dict(row_id=d["row_id"], status=f"error: {type(e).__name__}: {str(e)[:200]}")
                          for d in defs])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["search", "tables"])
    ap.add_argument("--defs", default=str(HERE / "defs.csv"))
    ap.add_argument("--check")
    ap.add_argument("--out", default=str(HERE / "results/v0"))
    ap.add_argument("--cache", default=CACHE)
    ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args()
    cand = set(pd.read_csv(HERE / "defs.csv").tic)
    if a.check:
        defs = pd.read_csv(a.check)
        if set(defs.tic) & cand:
            sys.exit("--check refuses definitions on a star of catalogue/defs.csv. Nothing was loaded.")
        accepted = None
    else:
        if not frozen():
            sys.exit("The catalogue run is not frozen: catalogue/APPROVALS.md has no step-2 approval. "
                     "Nothing was loaded.")
        defs = pd.read_csv(a.defs)
        accepted = None
        if a.stage == "tables":
            v = pd.read_csv(Path(a.out) / "vetting/verdicts.csv")
            accepted = set(v.event_id[(v.vetter_1 == "accept") & (v.vetter_2 == "accept")])
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    groups = [(g.to_dict("records"), a.cache, a.stage, accepted) for _, g in defs.groupby("tic", sort=False)]
    res = dict(rows=[], events=[], excl=[], tables=[])
    with mp.get_context("spawn").Pool(a.workers) as pool:
        for i, r in enumerate(pool.imap(star, groups), 1):
            print(f"[{i}/{len(groups)}]", flush=True)   # progress only, no outcome
            for k in res:
                res[k] += r.get(k, [])
    rows = defs.merge(pd.DataFrame(res["rows"]), on="row_id", how="left", validate="1:1")
    if a.stage == "search":
        rows.to_csv(out / "rows_search.csv", index=False)
        pd.DataFrame(res["events"], columns=EV_COLS).to_csv(out / "events.csv", index=False)
        pd.DataFrame(res["excl"], columns=["row_id", "p_lo_d", "p_hi_d"]).to_csv(out / "excluded_periods.csv", index=False)
    else:
        rows.to_csv(out / "catalogue.csv", index=False)
        pd.DataFrame(res["tables"], columns=["row_id", "table", "t1", "t2", "n", "period_d", "status",
                                             "fits_third"]).to_csv(out / "alias_tables.csv", index=False)


if __name__ == "__main__":
    main()
