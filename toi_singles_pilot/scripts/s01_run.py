"""
Run the repeat search (common.run_toi) on a target set, in parallel (one task per worker process).

  python toi_singles_pilot/scripts/s01_run.py validation [n_workers]   # the 16 community-solved TOIs, no diffs
  python toi_singles_pilot/scripts/s01_run.py pilot [n_workers]        # the 20-TOI pilot, with diff images

Inputs come ONLY from toi_singles_k2/targets_merged.csv columns TOI, TIC, Tmag, ra_deg, dec_deg,
epoch_btjd, dur_h, depth_ppm (no period, comment or upload columns are read, so the validation run is
blind to the uploaded periods). Outputs: toi_singles_pilot/results/<set>/<TOI>.json, grids, plots.
"""
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
PILOT = HERE.parent
ROOT = PILOT.parent
COLS = ["TOI", "TIC", "Tmag", "ra_deg", "dec_deg", "epoch_btjd", "dur_h", "depth_ppm", "group",
        "n_new_since_update"]

PILOT_14 = ["2270.01", "2296.01", "2098.01", "1301.02", "4585.01", "5975.01", "2102.02", "2298.01", "6075.01",
            "1893.01", "289.01", "4348.01", "4355.01", "6691.01"]


def targets(which):
    d = pd.read_csv(ROOT / "toi_singles_k2" / "targets_merged.csv", usecols=COLS, dtype={"TOI": str})
    d["TOI"] = d.TOI.map(lambda x: f"{float(x):.2f}")
    if which == "validation":
        s = d[d.group == "community/comment period after last TOI update"]
    else:
        o = d[d.group == "new TESS sectors since update, no newer period"].sort_values(
            ["n_new_since_update", "Tmag"], ascending=[False, True])
        extra = [x for x in o.TOI if x not in PILOT_14][:6]
        s = o[o.TOI.isin(PILOT_14 + extra)]
    s = s.rename(columns={"ra_deg": "ra", "dec_deg": "dec"})
    return [dict(r) for _, r in s.iterrows()]


def work(row, which):
    sys.path.insert(0, str(HERE))
    import common
    out = PILOT / "results" / which
    rep = common.run_safe(row, out, do_diff=(which == "pilot"))
    return row["TOI"], rep


def main():
    which = sys.argv[1]
    nw = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    out = PILOT / "results" / which
    out.mkdir(parents=True, exist_ok=True)
    rows = targets(which)
    print(f"{which}: {len(rows)} TOIs: {[r['TOI'] for r in rows]}", flush=True)
    with ProcessPoolExecutor(max_workers=nw, max_tasks_per_child=1) as ex:
        futs = [ex.submit(work, r, which) for r in rows]
        for i, fu in enumerate(as_completed(futs), 1):
            toi, rep = fu.result()
            per = [(round(p["P"], 3), p["n_support"]) for p in (rep.get("periods") or [])[:5]]
            print(f"[{i}/{len(rows)}] TOI {toi}: {rep.get('status')} sectors {rep.get('n_sectors')} event "
                  f"{rep.get('event_verdict')} snr {rep.get('event_snr')} repeats {rep.get('n_repeats')} "
                  f"periods(top5) {per} minP {rep.get('min_allowed_P_fine')} "
                  f"{rep.get('error', '')[-300:] if rep.get('status') == 'error' else ''}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
