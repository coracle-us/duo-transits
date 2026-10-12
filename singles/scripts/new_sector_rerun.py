"""
Monthly re-run: re-analyse the analysed single-transit CTOIs that have NEW TESS data, and flag
anything that changed.

What counts as new:
  - a sector that observed the star (tess-point), has ended at least 7 days ago (TESS orbit-times
    table, refreshed from tess.mit.edu), and has no SPOC/TESS-SPOC/QLP light curve in our cache.
    That includes sectors where we so far have only our own TESScut photometry ('OWN'), so pipeline
    products replace it when they appear. This generalises the weekly S107 recheck.
Per star with pending sectors:
  1. scripts/fetch_data._fetch_one: re-search MAST, download only new pipeline products;
  2. pilot.run_target: repeat search, periods, own TESScut photometry for sectors >= 100 that still
     lack pipeline products, difference images. Output goes to singles/results/rerun_<YYYY-MM-DD>/;
  3. compare with the star's previous JSON (latest of pilot/batch2/batch3/full/rerun_*):
     NEW REPEAT     a repeat more than D from every earlier repeat and from the event
     LOST REPEAT    an earlier repeat no longer found (e.g. an OWN-photometry repeat not confirmed
                    by the pipeline product)
     PERIODS        a change in the candidate-period list
     EVENT          a change in the event verdict
Stars whose file set did not change (data not on MAST yet) are listed as 'no new data yet', and
their run JSON is removed.
Nothing is emailed or posted. The report says which flagged stars need vetting by hand.

Output: singles/results/rerun_<date>/<TIC>.json, changes.csv, report.md
Usage:  python singles/scripts/new_sector_rerun.py [n_workers] [--max N]
Then:   python singles/scripts/classify_final.py   (it includes rerun_* directories, latest wins)
"""
import json
import re
import sys
import traceback
import urllib.request
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "scripts"))
RES = ROOT / "singles" / "results"
DATA = ROOT / "singles" / "data"
LC = ROOT / "data" / "raw" / "lightcurves"
ORBIT_URL = "https://tess.mit.edu/public/files/TESS_orbit_times.csv"
ORBIT_LATEST = DATA / "tess_orbit_times_latest.csv"
ORBIT_RAW = ROOT / "data" / "raw" / "tess_orbit_times.csv"
PIPE = ("SPOC", "TESS-SPOC", "QLP")
BASE_RUNS = ["pilot", "batch2", "batch3", "full"]


def orbit_table():
    try:
        with urllib.request.urlopen(ORBIT_URL, timeout=60) as r:
            txt = r.read().decode()
        if "Sector" in txt.splitlines()[0]:
            ORBIT_LATEST.write_text(txt)
    except Exception as e:  # noqa: BLE001
        print(f"orbit table refresh failed ({e}); using cached copy", flush=True)
    f = ORBIT_LATEST if ORBIT_LATEST.exists() else ORBIT_RAW
    o = pd.read_csv(f, comment="#")
    o["end"] = pd.to_datetime(o["End of Orbit"], utc=True)
    return o.groupby("Sector").end.max(), f


def cached(tic):
    pipe, own = set(), set()
    for f in (LC / str(tic)).glob("s*.csv.gz"):
        m = re.fullmatch(r"s(\d+)_(.+)_(\d+)s\.csv\.gz", f.name)
        if not m:
            continue
        (own if m.group(2) == "OWN" else pipe if m.group(2) in PIPE else set()).add(int(m.group(1)))
    return pipe, own


def file_set(tic):
    return {f.name for f in (LC / str(tic)).glob("s*.csv.gz")}


def previous_json(tic, exclude):
    runs = BASE_RUNS + sorted(p.name for p in RES.glob("rerun_*") if p.is_dir() and p.name != exclude)
    prev = None
    for r in runs:
        f = RES / r / f"{tic}.json"
        if f.exists():
            prev = (r, json.loads(f.read_text()))
    return prev


def work(args):
    row, outdir = args
    import fetch_data as fd
    import pilot
    tic = int(row["TIC"])
    before = file_set(tic)
    try:
        fd._fetch_one(0, tic)
        rep = pilot.run_target(pd.Series(row), Path(outdir), True)
    except Exception:  # noqa: BLE001
        return tic, dict(TIC=tic, status="error", error=traceback.format_exc()[-600:]), sorted(file_set(tic) - before)
    return tic, rep, sorted(file_set(tic) - before)


def compare(prev, new, D):
    out = []
    if prev is None:
        return ["no previous run"]
    pr = prev[1]
    old_T = [pr.get("event_T", 0)] + [x["T"] for x in pr.get("repeats") or []]
    new_T = [x["T"] for x in new.get("repeats") or []]
    for x in new.get("repeats") or []:
        if all(abs(x["T"] - t) > D for t in old_T):
            out.append(f"NEW REPEAT T={x['T']:.3f} S{x['sector']} {x['product']} depth {x['depth_ppm']} ppm SNR {x['snr']}")
    for x in pr.get("repeats") or []:
        if all(abs(x["T"] - t) > D for t in new_T):
            out.append(f"LOST REPEAT T={x['T']:.3f} S{x['sector']} {x['product']}")
    op = [round(p, 3) for _, p in (pr.get("candidate_periods") or [])]
    np_ = [round(p, 3) for _, p in (new.get("candidate_periods") or [])]
    if op != np_:
        out.append(f"PERIODS {len(op)} -> {len(np_)}: {np_[:6]}")
    if pr.get("event_verdict") != new.get("event_verdict"):
        out.append(f"EVENT {pr.get('event_verdict')} -> {new.get('event_verdict')}")
    return out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    n_workers = int(args[0]) if args else 3
    nmax = int(sys.argv[sys.argv.index("--max") + 1]) if "--max" in sys.argv else None
    today = date.today().isoformat()
    outdir = RES / f"rerun_{today}"
    outdir.mkdir(parents=True, exist_ok=True)
    ends, orbit_file = orbit_table()
    released = {int(s) for s, e in ends.items() if e < datetime.now(timezone.utc) - timedelta(days=7)}
    last_listed = int(max(ends.index))
    import tess_stars2px as tp
    cls = pd.read_csv(RES / "final_classification.csv")
    cat = pd.read_csv(DATA / "ctoi_singles.csv").drop_duplicates("TIC").set_index("TIC")
    todo, pending_info = [], {}
    for tic in cls.TIC.astype(int):
        if tic not in cat.index:
            continue
        r = cat.loc[tic]
        o = tp.tess_stars2px_function_entry(tic, float(r.RA), float(r.Dec))
        observed = {int(s) for s in o[3]}
        pipe, own = cached(tic)
        pending = sorted((observed & released) - pipe)
        pending = [s for s in pending if s <= 200]
        if pending:
            row = r.to_dict()
            row["TIC"] = tic
            todo.append((row, str(outdir)))
            pending_info[tic] = dict(pending=pending, own_only=sorted(set(pending) & own))
    n_pending_total = len(todo)
    if nmax:
        todo = todo[:nmax]
    print(f"orbit table: {orbit_file.name}, last sector listed {last_listed}; released sectors up to "
          f"{max(released)}; stars with pending sectors: {n_pending_total} (running {len(todo)})", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = [ex.submit(work, a) for a in todo]
        for i, fu in enumerate(as_completed(futs), 1):
            tic, rep, new_files = fu.result()
            f = outdir / f"{tic}.json"
            info = pending_info[tic]
            if rep.get("status") == "error":
                rows.append(dict(TIC=tic, status="error", pending=info["pending"], new_files="", changes=rep["error"][-200:]))
                print(f"[{i}/{len(todo)}] TIC {tic}: ERROR", flush=True)
                continue
            if not new_files:
                rows.append(dict(TIC=tic, status="no new data yet", pending=info["pending"], new_files="", changes=""))
                for g in outdir.glob(f"{tic}_*.png"):
                    g.unlink()
                continue
            f.write_text(json.dumps(rep, indent=1, default=float))
            ch = compare(previous_json(tic, outdir.name), rep, rep.get("dur_h", 5) / 24)
            rows.append(dict(TIC=tic, status="rerun", pending=info["pending"], new_files=" ".join(new_files),
                             changes=" | ".join(ch)))
            print(f"[{i}/{len(todo)}] TIC {tic}: new {new_files} -> {ch if ch else 'no change'}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(outdir / "changes.csv", index=False)
    flagged = df[df.changes.fillna("").str.contains("NEW REPEAT|LOST REPEAT|PERIODS|EVENT")] if len(df) else df
    fc = cls.set_index("TIC")
    lines = [f"# New-sector re-run {today}", "",
             f"- Orbit table: `{orbit_file.name}`; last sector listed {last_listed}; "
             f"sectors released (ended more than 7 d ago) up to S{max(released)}.",
             f"- Stars with pending sectors: {len(todo)}; re-run with new files: {int((df.status == 'rerun').sum()) if len(df) else 0}; "
             f"no new data yet: {int((df.status == 'no new data yet').sum()) if len(df) else 0}; "
             f"errors: {int((df.status == 'error').sum()) if len(df) else 0}.",
             "", "## Flagged changes (vet by hand before believing anything)", ""]
    if len(flagged):
        for _, r in flagged.iterrows():
            c0 = fc.loc[r.TIC] if r.TIC in fc.index else None
            prev = f"{c0.cls}, source {c0.source}, Gaia binary {c0.gaia_binary_flag}" if c0 is not None else ""
            lines.append(f"- **TIC {r.TIC}** ({prev}): {r.changes}")
    else:
        lines.append("None.")
    if (outdir / "..").exists() and len(todo) and last_listed <= max(released):
        lines += ["", f"Note: the orbit table lists sectors only to S{last_listed}. If it was not refreshed, "
                  "later sectors are invisible to this check."]
    (outdir / "report.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
