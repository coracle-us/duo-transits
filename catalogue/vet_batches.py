"""
Split a sheets folder into batches for one vetter, and join the batch verdicts back.
(Each vetter is several fresh agents, one per batch: catalogue/APPROVALS.md, deviation of 2026-11-23.)

Inputs  : <dir>/sheets/*.png|json
Outputs : split -> <dir>/batches/v<vetter>/b<k>/ (copies of the sheets, shuffled with numpy seed
          [20261205, vetter]); join -> <dir>/verdicts_<vetter>.csv from every b<k>/verdicts.csv
Usage   : python catalogue/vet_batches.py split|join <dir> <vetter 1|2> [--size 40]
"""
import argparse
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("mode", choices=["split", "join"])
ap.add_argument("dir")
ap.add_argument("vetter", type=int)
ap.add_argument("--size", type=int, default=40)
a = ap.parse_args()
d = Path(a.dir)
root = d / "batches" / f"v{a.vetter}"
if a.mode == "split":
    ids = sorted(p.stem for p in (d / "sheets").glob("*.png"))
    ids = list(np.random.default_rng([20261205, a.vetter]).permutation(ids))
    for k in range(0, len(ids), a.size):
        b = root / f"b{k // a.size:03d}"
        b.mkdir(parents=True, exist_ok=True)
        for s in ids[k:k + a.size]:
            for ext in ("png", "json"):
                shutil.copy(d / "sheets" / f"{s}.{ext}", b / f"{s}.{ext}")
    print(len(ids), "sheets in", len(list(root.glob("b*"))), "batches")
else:
    parts, missing = [], []
    for b in sorted(root.glob("b*")):
        f = b / "verdicts.csv"
        if not f.exists():
            missing.append(b.name)
            continue
        v = pd.read_csv(f, dtype=str).fillna("")
        want = {p.stem for p in b.glob("*.png")}
        assert set(v.sheet_id) == want and v.sheet_id.is_unique, f"{b}: verdicts do not match the sheets"
        parts.append(v[["sheet_id", "verdict", "reason", "note"]])
    assert not missing, f"no verdicts in {missing}"
    pd.concat(parts).to_csv(d / f"verdicts_{a.vetter}.csv", index=False)
    print(sum(map(len, parts)), "verdicts")
