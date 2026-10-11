"""
Unattended run over ALL remaining single-transit CTOIs (those not in an earlier batch).

Order: expected event SNR (depth_ppm * 10^(-0.2 (Tmag - 10))), depth 200-50000 ppm, Tmag >= 7.
Per target: ExoFOP community-period check (skip + log if someone already posted a period),
then pilot.run_target (event check, repeat search, period inference, own photometry for
unprocessed sectors >= 100, difference images). One JSON per target, written as soon as it
is done, so the run is resumable: targets with an existing JSON are skipped.

Usage:  python singles/scripts/run_remaining.py [n_workers]
Output: singles/results/full/<TIC>.json, singles/results/full/skipped.csv, run log on stdout
Then:   python singles/scripts/summarize_pilot.py full
"""
import json
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "scripts"))
OUT = ROOT / "singles" / "results" / "full"
DATA = ROOT / "singles" / "data"


def work(row_dict):
    import pilot
    import select_pilot as sp
    import fetch_exofop as fx
    row = pd.Series(row_dict)
    tic = int(row.TIC)
    f = OUT / f"{tic}.json"
    if f.exists():
        return tic, "exists"
    try:
        fx.fetch(tic)
        hits = sp.community_period(tic)
        if hits:
            return tic, "skip: community period on ExoFOP: " + " | ".join(hits)[:300]
        rep = pilot.run_target(row, OUT, True)
    except Exception:  # noqa: BLE001
        rep = dict(TIC=tic, status="error", error=traceback.format_exc()[-600:])
    f.write_text(json.dumps(rep, indent=1, default=float))
    return tic, (f"repeats {rep.get('n_repeats')} periods {rep.get('candidate_periods', [])[:4]} "
                 f"diff {[(d.get('kind'), d.get('offset_arcsec'), d.get('snr_peak')) for d in rep.get('diff_images', []) or []][:3]}")


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    OUT.mkdir(parents=True, exist_ok=True)
    c = pd.read_csv(DATA / "ctoi_singles.csv").drop_duplicates("TIC")
    done = set(pd.read_csv(DATA / "done_batches.csv").TIC)
    for b in ("batch3_targets.csv",):
        if (DATA / b).exists():
            done |= set(pd.read_csv(DATA / b).TIC)
    c = c[~c.TIC.isin(done) & (c.Tmag >= 7) & (c.depth_ppm >= 200) & (c.depth_ppm <= 50000)].copy()
    c["snr_score"] = c.depth_ppm * 10 ** (-0.2 * (c.Tmag - 10))
    c = c.sort_values("snr_score", ascending=False)
    todo = [r for r in c.to_dict("records") if not (OUT / f"{int(r['TIC'])}.json").exists()]
    print(f"remaining: {len(c)} targets; to do now: {len(todo)}; workers: {n_workers}", flush=True)
    skipped = []
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = [ex.submit(work, r) for r in todo]
        for i, fu in enumerate(as_completed(futs), 1):
            tic, msg = fu.result()
            print(f"[{i}/{len(todo)}] TIC {tic}: {msg}", flush=True)
            if msg.startswith("skip"):
                skipped.append(dict(TIC=tic, reason=msg))
                pd.DataFrame(skipped).to_csv(OUT / "skipped.csv", index=False)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
