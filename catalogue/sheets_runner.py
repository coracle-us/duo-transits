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
import pickle
import subprocess
from multiprocessing.pool import ThreadPool
import socket
import sys
import time
from pathlib import Path

import pandas as pd

STAR_TIMEOUT = 1500
socket.setdefaulttimeout(120)          # also in every spawned worker, which imports this module
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_sheets as ms  # noqa: E402
import numpy as np  # noqa: E402
from astroquery.mast import Catalogs  # noqa: E402


def _tic_by_cone(catalog=None, coordinates=None, radius=None, Tmag=None, **kw):
    """Stand-in for Catalogs.query_criteria(catalog="Tic", coordinates, radius, Tmag=[lo, hi]) as
    sheets.neighbours calls it. From 2026-11-24 about 11:35 CST MAST's filtered TIC service stopped
    answering (the request hangs); the cone service answers. Same stars: every TIC entry within
    `radius` degrees, kept if lo <= Tmag <= hi."""
    if coordinates is None:               # e.g. the look-up by TIC ID in pl.tic_info: unchanged
        return _ORIG_QUERY_CRITERIA(catalog=catalog, radius=radius, **kw) if radius is not None else _ORIG_QUERY_CRITERIA(catalog=catalog, **kw)
    from astropy.coordinates import SkyCoord
    ra, dec = (float(x) for x in str(coordinates).split())
    # sheets.neighbours uses only the 10 nearest of the list, so a smaller cone that already holds
    # 20 such stars gives the same 10. Crowded fields do not answer at 0.5 deg.
    for rad in (0.12, 0.5):
        rad = min(rad, radius)
        r = Catalogs.query_region(SkyCoord(ra, dec, unit="deg"), radius=rad, catalog="TIC")
        tm = np.ma.filled(np.ma.asarray(r["Tmag"], dtype=float), np.nan)
        r = r[(tm >= Tmag[0]) & (tm <= Tmag[1])]
        if len(r) >= 20 or rad >= radius:
            return r


_ORIG_QUERY_CRITERIA = Catalogs.query_criteria
Catalogs.query_criteria = _tic_by_cone   # in every worker process, which imports this module


def one(job):
    mode, name, args, ckdir = job
    ck = Path(ckdir) / f"{name}.json"
    if ck.exists():
        return json.loads(ck.read_text())
    fn = ms.star_v1_sheets if mode == "vetter-check" else ms.star_candidates
    for attempt in range(2):
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
    jobdir = base / "_jobs"
    jobdir.mkdir(exist_ok=True)
    jobs = [(m, n, x, str(ckdir)) for m, n, x, _ in jobs]

    def guarded(job):
        """One star in its own process, killed after STAR_TIMEOUT s (a MAST read can hang for hours
        whatever the socket timeout) and tried again, up to 3 times."""
        ck = ckdir / f"{job[1]}.json"
        if ck.exists():
            return json.loads(ck.read_text())
        jf = jobdir / f"{job[1]}.pkl"
        jf.write_bytes(pickle.dumps(job))
        for _ in range(3):
            try:
                subprocess.run([sys.executable, "-W", "ignore", __file__, "--one", str(jf)], timeout=STAR_TIMEOUT + 300 * len(job[2][1] if job[0] == "vetter-check" else job[2][3]),
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except subprocess.TimeoutExpired:
                pass
            if ck.exists():
                return json.loads(ck.read_text())
        return [dict(sheet_id="", tic=job[1], error="no sheets after 3 tries (timeout or error)")]

    key = []
    with ThreadPool(a.workers) as pool:
        for i, r in enumerate(pool.imap_unordered(guarded, jobs), 1):
            print(f"[{i}/{len(jobs)}]", flush=True)
            key += r
    pd.DataFrame(key).to_csv(base / "key.csv", index=False)


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--one":
        one(pickle.loads(Path(sys.argv[2]).read_bytes()))
    else:
        main()
