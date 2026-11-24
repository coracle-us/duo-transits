"""
Vetting sheets (sheets.py) for the catalogue run, and for the vetter check V1.

  candidates   : one sheet per passing event of <out>/events.csv, except a Hawthorn row's events at its
                 catalogued second transit (that transit is catalogued, not vetted).
  vetter-check : V1 - round 3's test trials are prepared again exactly as harness3 prepared them; every
                 passing event not at a true transit (arm A1 and arm B) gets a sheet, and an equal number
                 of passing events at a true transit (arm A2 and arm B), drawn with numpy seed 20261205.
                 Arm-B sheets show the injected planet in every product of the target.
Both refuse to start before the freeze (catalogue/APPROVALS.md). --check DEFS makes sheets for
calibration definitions on stars outside catalogue/defs.csv (a test of the generator).

Sheets carry a hashed name. The key (sheet -> event, and for V1 the truth) is written beside the
sheets folder, never inside it; vetters are given the sheets folder only.

Inputs  : catalogue/defs.csv, <out>/rows_search.csv, <out>/events.csv; calibration/mono/round3/test/*_defs.csv
Outputs : <out>/vetting/sheets/*.png|json + <out>/vetting/key.csv
          <out>/vetter_check/sheets/*.png|json + <out>/vetter_check/key.csv
Usage   : python catalogue/make_sheets.py candidates|vetter-check [--out catalogue/results/v0] [--workers 8]
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
import sheets  # noqa: E402
from run_catalogue import CACHE, frozen  # noqa: E402

mc, pl, dump, dump3, h2, h3 = catrun.mc, catrun.mc.pl, catrun.dump, catrun.dump3, catrun.h2, catrun.h3
R3 = HERE.parent / "calibration/mono/round3/test"
V1_SEED = 20261205


def _ref(row):
    return dict(depth=row["depth_ref"], depth_err=row["depth_ref_err"], dur_d=row["dur_ref_d"],
                snr=row["ref_snr"], measured=bool(row["ref_measured"]), sector=int(row["ref_sector"]))


def star_candidates(args):
    tic, defs, rows, events, out, cache = args
    try:
        sectors, _ = mc.load_star(tic, Path(cache))
        sr, prods = pl.lc_products(tic)
        info = pl.tic_info(tic)
        key = []
        for ev in events:
            d, row = defs[ev["row_id"]], rows[ev["row_id"]]
            sid = sheets.make_sheet(d, ev, _ref(row), sectors, sr, prods, info, out, cache)
            key.append(dict(sheet_id=sid, event_id=ev["event_id"], row_id=ev["row_id"]))
        return key
    except Exception as e:
        return [dict(sheet_id="", event_id=ev["event_id"], row_id=ev["row_id"],
                     error=f"{type(e).__name__}: {str(e)[:200]}") for ev in events]


def trial_events(d, sectors, tic):
    """Passing events of one round-3 trial, with at_truth, prepared as harness3.score_trial does."""
    prep, _ = dump3.prep_trial(d, sectors, tic)
    if prep is None:
        return [], None, None
    s = dump.summarize(prep, d["t_cat"], d["truth_P"], d["truth_t0"], grid=False)
    out = []
    for e in s["events"]:
        if "error" in e:
            continue
        e["sec_box_snr"] = dump3.sector_box_snr(prep["det"], e, prep["dur_ref"])
        if catrun.rule3.vet3(e, h3.FROZEN3, s["ref_snr"]):
            out.append(e)
    ref = dict(depth=prep["depth_ref"], depth_err=prep["derr_ref"], dur_d=prep["dur_ref"], snr=s["ref_snr"],
               measured=bool(s["ref_measured"]), sector=int(s["ref_sector"]))
    return out, ref, prep


def star_v1_list(args):
    """Stage 1 of V1: list the passing events of a star's trials (no sheets yet)."""
    defs, cache = args
    tic = int(defs[0]["tic"])
    try:
        sectors, _ = mc.load_star(tic, Path(cache))
        out = []
        for i, d in enumerate(defs):
            evs, ref, _ = trial_events(d, sectors, tic) if sectors else ([], None, None)
            for e in evs:
                # false events from A1 and B, true events from A2 and B (an A-star event is in A1 and A2)
                if (d["arm"] == "A1" and e["at_truth"]) or (d["arm"] == "A2" and not e["at_truth"]):
                    continue
                out.append(dict(tic=tic, arm=d["arm"], trial=d["trial"], sector=int(e["sector"]), tc=e["tc"],
                                at_truth=bool(e["at_truth"]), snr=e["snr"]))
        return out
    except Exception as e:
        return [dict(tic=tic, error=f"{type(e).__name__}: {str(e)[:200]}")]


def star_v1_sheets(args):
    defs, picks, out, cache = args
    tic = int(defs[0]["tic"])
    try:
        sectors, _ = mc.load_star(tic, Path(cache))
        sr, prods = pl.lc_products(tic)
        info = pl.tic_info(tic)
        key = []
        for d in defs:
            want = [p for p in picks if p["trial"] == d["trial"]]
            if not want:
                continue
            evs, ref, prep = trial_events(d, sectors, tic)
            inject, shown = None, sectors
            if d["arm"] == "B":
                p = dict(t0=d["inj_t0"], P=d["inj_P"], k=d["inj_k"], ar=d["inj_ar"], b=d["inj_b"])
                u = (d["u1"], d["u2"])

                def inject(s, a, t, f, p=p, u=u):
                    return f * mc.transit_flux(t, p["t0"], p["P"], p["k"], p["ar"], p["b"], u, mc.core.cadence_days(t))
                shown = [(s, a, t, inject(s, a, t, f)) for s, a, t, f in sectors]
            for w in want:
                e = next(x for x in evs if int(x["sector"]) == w["sector"] and abs(x["tc"] - w["tc"]) < 1e-6)
                e = dict(e, event_id=f"v1|{d['trial']}|s{int(e['sector'])}|{e['tc']:.4f}", at_other_row=False)
                dd = dict(d, row_id=d["trial"])
                sid = sheets.make_sheet(dd, e, ref, shown, sr, prods, info, out, cache, inject=inject, salt="v1")
                key.append(dict(sheet_id=sid, event_id=e["event_id"], arm=d["arm"], tic=tic, at_truth=w["at_truth"]))
        return key
    except Exception as e:
        return [dict(sheet_id="", tic=tic, error=f"{type(e).__name__}: {str(e)[:200]}")]


def round3_defs():
    a = pd.read_csv(R3 / "arm_a_defs.csv")
    b = pd.read_csv(R3 / "arm_b_defs.csv")
    a["trial"] = a.arm + ":" + a.toi.astype(str)
    b["trial"] = "B:" + b.tic.astype(str) + ":" + b.draw.astype(str)
    return pd.concat([a, b], ignore_index=True)


def _pool(fn, args, workers):
    with mp.get_context("spawn").Pool(workers) as pool:
        for i, r in enumerate(pool.imap(fn, args), 1):
            print(f"[{i}/{len(args)}]", flush=True)
            yield r


def pick_v1(listed, seed=V1_SEED):
    """Every false event, and as many true events drawn without replacement (all, if fewer)."""
    listed = listed.sort_values(["tic", "trial", "sector", "tc"]).reset_index(drop=True)
    false, true = listed[~listed.at_truth], listed[listed.at_truth]
    rng = np.random.default_rng(seed)
    take = true.iloc[np.sort(rng.choice(len(true), size=min(len(false), len(true)), replace=False))]
    return pd.concat([false, take], ignore_index=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["candidates", "vetter-check"])
    ap.add_argument("--out", default=str(HERE / "results/v0"))
    ap.add_argument("--check", help="vetter-check on these calibration definitions only (generator test)")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--cache", default=CACHE)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    out = Path(a.out)
    cand = set(pd.read_csv(HERE / "defs.csv").tic)
    if a.check:
        defs = pd.read_csv(a.check)
        if set(defs.tic) & cand:
            sys.exit("--check refuses definitions on a star of catalogue/defs.csv. Nothing was loaded.")
    elif not frozen():
        sys.exit("The catalogue run is not frozen: catalogue/APPROVALS.md has no step-2 approval. Nothing was loaded.")
    if a.mode == "candidates":
        defs = pd.read_csv(HERE / "defs.csv").set_index("row_id", drop=False)
        rows = pd.read_csv(out / "rows_search.csv").set_index("row_id", drop=False)
        ev = pd.read_csv(out / "events.csv")
        ev = ev[ev.passes & ~ev.at_t2_cat]
        ev["tic"] = ev.row_id.map(defs.tic)
        args = [(int(t), {r: defs.loc[r].to_dict() for r in g.row_id.unique()},
                 {r: rows.loc[r].to_dict() for r in g.row_id.unique()}, g.to_dict("records"),
                 out / "vetting/sheets", a.cache) for t, g in ev.groupby("tic")]
        key = [k for r in _pool(star_candidates, args, a.workers) for k in r]
        pd.DataFrame(key).to_csv(out / "vetting/key.csv", index=False)
        return
    if not a.check:
        defs = round3_defs()
    elif "trial" not in defs:
        defs["trial"] = defs.arm + ":" + defs.get("toi", defs.get("draw")).astype(str)
    groups = [g.to_dict("records") for _, g in defs.groupby("tic", sort=False)]
    if a.limit:
        groups = groups[:a.limit]
    vc = out / "vetter_check"
    vc.mkdir(parents=True, exist_ok=True)
    listed = pd.DataFrame([x for r in _pool(star_v1_list, [(g, a.cache) for g in groups], a.workers) for x in r])
    listed.to_csv(vc / "events_listed.csv", index=False)
    ok = listed[listed.get("error").isna()] if "error" in listed else listed
    picks = pick_v1(ok.astype({"at_truth": bool}))
    by = {t: g.to_dict("records") for t, g in picks.groupby("tic")}
    args = [(g, by[int(g[0]["tic"])], vc / "sheets", a.cache) for g in groups if int(g[0]["tic"]) in by]
    key = [k for r in _pool(star_v1_sheets, args, a.workers) for k in r]
    pd.DataFrame(key).sample(frac=1.0, random_state=V1_SEED).to_csv(vc / "key.csv", index=False)


if __name__ == "__main__":
    main()
