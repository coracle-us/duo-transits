"""
Independent re-derivation of the catalogue run's search stage (steps 1 and 2 of
catalogue/prereg_catalogue.yaml) for every row of catalogue/defs.csv.

Written without reading catalogue/catrun.py, run_catalogue.py, report.py, the tests, the driver
check outputs or anything under catalogue/results/ (other than this script's own folder).
Shared, frozen pipeline functions are imported and called (monocal, round2/harness, rule2, dump,
round3/dump3, rule3, st2core, st2_pipeline); the glue, the R4/T5 test, sec_box_snr, the excluding
time grid and the R5 lookup are written here from the frozen rule texts, with the shared versions
called beside them as a cross-check where one exists.

Inputs  : catalogue/defs.csv, MAST light curves through mc.load_star (shared cache)
Outputs : catalogue/results/v0/independent/stars/<tic>.json (checkpoint, one per star),
          ind_rows.csv, ind_events.csv, ind_log.txt
Usage   : caffeinate -i .venv/bin/python -W ignore catalogue/independent/ind_search.py --run [--workers 6]
          .venv/bin/python -W ignore catalogue/independent/ind_search.py --star TIC   (one star, used by --run)
          .venv/bin/python -W ignore catalogue/independent/ind_search.py --collect
"""
import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "catalogue/results/v0/independent"
STARS = OUT / "stars"
DEFS = ROOT / "catalogue/defs.csv"
CACHE = "/Users/mcbee/Documents/GitHub/coracle/duo-transits/salinas-second-transits/work/lc_cache"
TIMEOUT_S = 900
TRIES = 3

# R5 table of round2/FROZEN_RULE.md: (label, k, n); p = (k + 0.5) / (n + 1)
R5 = [("unmeasured", 0, 507), ("5-8", 0, 211), ("8-12", 0, 202), ("12-20", 0, 192), ("over 20", 1, 511)]


def r5_bin(ref_measured, ref_snr):
    if not ref_measured or not np.isfinite(ref_snr) or ref_snr < 5:
        i = 0
    elif ref_snr < 8:
        i = 1
    elif ref_snr < 12:
        i = 2
    elif ref_snr < 20:
        i = 3
    else:
        i = 4
    label, k, n = R5[i]
    return label, (k + 0.5) / (n + 1)


def own_passes_r4(e, dur_unused=None):
    """R4 as written in round2/FROZEN_RULE.md."""
    return bool(e["snr"] >= 8.0
                and 1 / 1.7 <= e["depth_ratio"] <= 1.7
                and 0.5 <= e["dur_ratio"] <= 2.0
                and e["snr"] >= 1.5 * e["anti_snr"]
                and e["frac"] >= 0.8 and e["edge_dist_d"] >= 0.5)


def own_sec_box_snr(core, det, e, dur):
    """T5's sec_box_snr as written in round3/FROZEN_RULE3.md."""
    if e["sector"] not in det:
        return np.nan
    t, f = det[e["sector"]]
    cen = np.arange(t.min() + dur / 2, t.max() - dur / 2, dur / 2)
    if cen.size == 0:
        return np.nan
    d, fr, _ = core.box_stats(t, f, cen, dur)
    d = d[(fr >= 0.8) & (np.abs(cen - e["tc"]) > dur)]
    if d.size < 10:
        return np.nan
    scat = 1.4826 * np.median(np.abs(d - np.median(d)))
    return float(e["depth"] / scat) if scat > 0 else np.nan


def own_excluding_grid(core, T, F, dur, depth0, derr, sig, measured, w):
    """R1-R3 on the time grid of step dur/8 starting at min(T), written from the rule text:
    one box per grid point and per offset, no cumulative-sum shortcut. Returns (c0, step, mask)."""
    step = dur / 8
    cen = np.arange(T.min(), T.max() + step, step)
    n = cen.size
    out = np.zeros(n, bool)
    if not measured:
        return float(cen[0]), step, out
    edges = core.gap_edges(T)
    d, fr, _ = core.box_stats(T, F, cen, dur)
    err = sig / np.sqrt(np.maximum(fr, 0.001))
    thr = 0.5 * depth0 - 3 * np.hypot(err, derr)
    clean = (fr >= 0.8) & (core.edge_distance(edges, cen) >= 0.5)
    cand = clean & (d >= -3 * err) & (d < thr)
    idx = np.flatnonzero(cand)
    ok = np.ones(idx.size, bool)
    for j in range(-w, w + 1):
        cj = cen[idx] + j * step
        dj, fj, _ = core.box_stats(T, F, cj, dur)
        ej = sig / np.sqrt(np.maximum(fj, 0.001))
        tj = 0.5 * depth0 - 3 * np.hypot(ej, derr)
        with np.errstate(invalid="ignore"):
            ok &= ~((fj >= 0.5) & (dj >= tj))
    out[idx[ok]] = True
    return float(cen[0]), step, out


def run_star(tic, accept_partial):
    sys.path.insert(0, str(ROOT / "calibration/mono/round3"))
    sys.path.insert(0, str(ROOT / "calibration/mono/round2"))
    import harness as h2
    import dump3
    import rule3
    import harness3
    mc, core, dump, rule2 = h2.mc, h2.core, h2.dump, h2.rule2

    defs = pd.read_csv(DEFS)
    defs = defs[defs.tic == tic]
    assert len(defs), tic
    sectors, errors = mc.load_star(int(tic), Path(CACHE))
    if errors and not accept_partial:
        print("LOAD ERRORS", errors, flush=True)
        sys.exit(3)
    rows, events = [], []
    used = ";".join(f"{s}:{a}" for s, a, _, _ in sectors)
    for d in defs.to_dict("records"):
        row = dict(row_id=d["row_id"], tic=int(tic), catalogue=d["catalogue"], sectors=used,
                   n_load_errors=len(errors), load_errors=" | ".join(errors))
        if not sectors:
            row["status"] = "no light curve" if not errors else "load error"
            rows.append(row)
            continue
        multi = bool(d["multi"])
        prep = mc.prepare(sectors, float(d["t_cat"]), float(d["dur_cat_d"]), float(d["depth_cat"]),
                          float(d["depth_cat_err"]))
        if prep is None:
            row["status"] = "reference not in data"
            rows.append(row)
            continue
        o = prep["out"]
        T, F, dur = prep["T"], prep["F"], prep["dur_ref"]
        depth0, derr, meas = prep["depth_ref"], prep["derr_ref"], bool(o["ref_measured"])
        s = dump.summarize(prep, float(d["t_cat"]))
        sig = s["sig"]
        # shared path, as round2/harness.score
        ex = rule2.timeline_excludes(s["g_depth"], s["g_frac"], s["g_clean"], s["g_idx"], s["n_grid"], sig,
                                     depth0, derr, h2.FROZEN, meas, multi)
        per = mc.period_grid(np.ptp(T), dur)
        frac = float(mc.excluded_periods(float(d["t_cat"]), per, s["c0"], s["step"], ex, T.min(), T.max()).mean())
        # own excluding grid, same period folding
        w = 8 if multi else 4
        c0, step, ex_own = own_excluding_grid(core, T, F, dur, depth0, derr, sig, meas, w)
        frac_own = float(mc.excluded_periods(float(d["t_cat"]), per, c0, step, ex_own, T.min(), T.max()).mean())
        label, p = r5_bin(meas, o["ref_snr"])
        ev = [e for e in s["events"] if "error" not in e]
        n_pass = 0
        for e in ev:
            sb = own_sec_box_snr(core, prep["det"], e, dur)
            sb_shared = dump3.sector_box_snr(prep["det"], e, dur)
            r4 = own_passes_r4(e)
            ok = bool(r4 and np.isfinite(sb) and sb >= 9.0)
            e2 = dict(e, sec_box_snr=sb_shared)
            shared_ok = bool(rule3.vet3(e2, harness3.FROZEN3, o["ref_snr"]))
            n_pass += ok
            events.append(dict(row_id=d["row_id"], tic=int(tic), sector=int(e["sector"]), author=e["author"],
                               tc=e["tc"], snr=e["snr"], depth=e["depth"], depth_err=e["depth_err"],
                               dur_d=e["dur_d"], depth_ratio=e["depth_ratio"], dur_ratio=e["dur_ratio"],
                               anti_snr=e["anti_snr"], frac=e["frac"], edge_dist_d=e["edge_dist_d"],
                               sec_box_snr=sb, sec_box_snr_shared=sb_shared, passes_r4=r4, passes=ok,
                               passes_shared=shared_ok))
        row.update(status="ok", ref_measured=meas, ref_snr=o["ref_snr"], depth_ref=depth0, dur_ref_d=dur,
                   derr_ref=derr, ref_sector=int(o["ref_sector"]), ref_author=o["ref_author"],
                   n_sectors=o["n_sectors"], n_points=o["n_points"], baseline_d=o["baseline_d"], sig=sig,
                   multi=multi, w_steps=w, frac_lnP_excluded=frac, frac_lnP_own=frac_own,
                   n_grid_periods=int(per.size), snr_bin=label, p_false_exclusion=p,
                   n_events=len(ev), n_passing=int(n_pass),
                   n_sector_errors=sum("error" in e for e in s["events"]))
        rows.append(row)

    def clean(v):
        if isinstance(v, (np.bool_, bool)):
            return bool(v)
        if isinstance(v, (np.integer,)):
            return int(v)
        if isinstance(v, (np.floating, float)):
            return None if not np.isfinite(v) else float(v)
        return v
    STARS.mkdir(parents=True, exist_ok=True)
    tmp = STARS / f"{tic}.json.tmp"
    tmp.write_text(json.dumps(dict(rows=[{k: clean(v) for k, v in r.items()} for r in rows],
                                   events=[{k: clean(v) for k, v in e.items()} for e in events])))
    tmp.rename(STARS / f"{tic}.json")


def one(tic):
    path = STARS / f"{tic}.json"
    if path.exists():
        return tic, "cached"
    msg = ""
    for attempt in range(1, TRIES + 1):
        cmd = [sys.executable, "-W", "ignore", __file__, "--star", str(tic)]
        if attempt == TRIES:
            cmd.append("--accept-partial")
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_S)
            if path.exists():
                return tic, f"ok (attempt {attempt})"
            msg = f"exit {r.returncode}: {(r.stdout + r.stderr).strip()[-300:]}"
        except subprocess.TimeoutExpired:
            msg = "timeout"
        time.sleep(5 * attempt)
    return tic, "FAILED " + msg.replace("\n", " ")


def run(workers):
    STARS.mkdir(parents=True, exist_ok=True)
    tics = list(dict.fromkeys(pd.read_csv(DEFS).tic.astype(int)))
    log = open(OUT / "ind_log.txt", "a")
    with ThreadPoolExecutor(workers) as pool:
        for i, (tic, msg) in enumerate(pool.map(one, tics), 1):
            if msg != "cached":
                log.write(f"[{i}/{len(tics)}] {tic} {msg}\n")   # no outcome, only whether it ran
                log.flush()
    log.close()
    print("done", flush=True)


def collect():
    defs = pd.read_csv(DEFS)
    rows, events, missing = [], [], []
    for tic in dict.fromkeys(defs.tic.astype(int)):
        p = STARS / f"{tic}.json"
        if not p.exists():
            missing.append(tic)
            rows += [dict(row_id=r, tic=tic, status="error") for r in defs[defs.tic == tic].row_id]
            continue
        j = json.loads(p.read_text())
        rows += j["rows"]
        events += j["events"]
    R = pd.DataFrame(rows).set_index("row_id").loc[defs.row_id].reset_index()
    R.to_csv(OUT / "ind_rows.csv", index=False)
    pd.DataFrame(events).to_csv(OUT / "ind_events.csv", index=False)
    print(len(R), "rows", len(events), "events", len(missing), "stars missing", missing[:20])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--star", type=int)
    ap.add_argument("--accept-partial", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--collect", action="store_true")
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    if a.star is not None:
        run_star(a.star, a.accept_partial)
    elif a.run:
        run(a.workers)
    elif a.collect:
        collect()
