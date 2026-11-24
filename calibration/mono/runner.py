"""Run one function per star in a process pool; one JSON per star, so a run can be resumed."""
import json
import multiprocessing as mp
import time
from pathlib import Path

import pandas as pd


def _call(args):
    fn, item, kw, path = args
    t0 = time.time()
    try:
        rows, meta = fn(item, **kw)
    except Exception as e:  # recorded, never silently dropped
        rows, meta = [], dict(tic=int(item["tic"]), error=f"{type(e).__name__}: {str(e)[:300]}")
    meta["wall_s"] = time.time() - t0
    Path(path).write_text(json.dumps(dict(rows=rows, meta=meta), default=float))
    return meta


def run(fn, items, out_dir, workers, **kw):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    todo = [(fn, it, kw, str(out_dir / f"{int(it['tic'])}.json")) for it in items
            if not (out_dir / f"{int(it['tic'])}.json").exists()]
    print(f"{len(items) - len(todo)} done already, {len(todo)} to run", flush=True)
    with mp.get_context("spawn").Pool(workers) as pool:
        for i, meta in enumerate(pool.imap_unordered(_call, todo), 1):
            # progress only: no outcome is printed
            print(f"[{i}/{len(todo)}] TIC {meta.get('tic')} sectors={meta.get('n_sectors')} "
                  f"trials={meta.get('n_trials')} wall={meta['wall_s']:.0f}s {meta.get('error', '')}", flush=True)


def collect(out_dir):
    rows, metas = [], []
    for p in sorted(Path(out_dir).glob("*.json")):
        d = json.loads(p.read_text())
        rows += d["rows"]
        metas.append(d["meta"])
    return pd.DataFrame(rows), pd.DataFrame(metas)
