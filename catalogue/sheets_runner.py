"""
Run the frozen sheet generators (make_sheets.star_v1_sheets, make_sheets.star_candidates) star by
star with a socket timeout, retries and per-star checkpoints. Written after the freeze, because the
frozen make_sheets.main hung on MAST with no timeout (2026-11-23). It changes no sheet: the stars,
the picks (make_sheets.pick_v1, same seed) and the sheet code are the frozen ones.

Inputs  : as make_sheets.py; for vetter-check, <out>/vetter_check/events_listed.csv (written by the
          frozen make_sheets.py before it hung)
Outputs : <out>/vetter_check/sheets/, <out>/vetter_check/key.csv; <out>/vetting/sheets/, <out>/vetting/key.csv
Usage   : python catalogue/sheets_runner.py vetter-check|candidates [--out catalogue/results/v0] [--workers 4]
"""
import argparse
import json
import multiprocessing as mp
import socket
import sys
import time
from pathlib import Path

import pandas as pd

socket.setdefaulttimeout(120)          # also in every spawned worker, which imports this module
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_sheets as ms  # noqa: E402


def one(job):
    mode, name, args, ckdir = job
    ck = Path(ckdir) / f"{name}.json"
    if ck.exists():
        return json.loads(ck.read_text())
    fn = ms.star_v1_sheets if mode == "vetter-check" else ms.star_candidates
    for attempt in range(4):
        key = fn(args)
        if not any("error" in k and k["error"] == k["error"] and k["error"] for k in key):
            ck.write_text(json.dumps(key, default=str))
            return key
        time.sleep(30 * (attempt + 1))
    return key                           # still failing: not checkpointed, reported in key.csv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["candidates", "vetter-check"])
    ap.add_argument("--out", default=str(HERE / "results/v0"))
    ap.add_argument("--cache", default=ms.CACHE)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    if not ms.frozen():
        sys.exit("not frozen")
    out = Path(a.out)
    if a.mode == "vetter-check":
        base = out / "vetter_check"
        defs = ms.round3_defs()
        groups = [g.to_dict("records") for _, g in defs.groupby("tic", sort=False)]
        listed = pd.read_csv(base / "events_listed.csv")
        ok = listed[listed["error"].isna()] if "error" in listed else listed
        picks = ms.pick_v1(ok.astype({"at_truth": bool}))
        by = {t: g.to_dict("records") for t, g in picks.groupby("tic")}
        jobs = [(a.mode, str(int(g[0]["tic"])), (g, by[int(g[0]["tic"])], base / "sheets", a.cache), None)
                for g in groups if int(g[0]["tic"]) in by]
    else:
        base = out / "vetting"
        defs = pd.read_csv(HERE / "defs.csv").set_index("row_id", drop=False)
        rows = pd.read_csv(out / "rows_search.csv").set_index("row_id", drop=False)
        ev = pd.read_csv(out / "events.csv")
        ev = ev[ev.passes & ~ev.at_t2_cat]
        ev["tic"] = ev.row_id.map(defs.tic)
        jobs = [(a.mode, str(int(t)), (int(t), {r: defs.loc[r].to_dict() for r in g.row_id.unique()},
                                       {r: rows.loc[r].to_dict() for r in g.row_id.unique()},
                                       g.to_dict("records"), base / "sheets", a.cache), None)
                for t, g in ev.groupby("tic")]
    ckdir = base / "_checkpoints"
    ckdir.mkdir(parents=True, exist_ok=True)
    jobs = [(m, n, x, str(ckdir)) for m, n, x, _ in jobs]
    key = []
    with mp.get_context("spawn").Pool(a.workers) as pool:
        for i, r in enumerate(pool.imap_unordered(one, jobs), 1):
            print(f"[{i}/{len(jobs)}]", flush=True)
            key += r
    pd.DataFrame(key).to_csv(base / "key.csv", index=False)


if __name__ == "__main__":
    main()
