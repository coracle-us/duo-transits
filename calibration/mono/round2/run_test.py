"""
Run the round-2 test set. Refuses to start until APPROVALS.md records Jon's step-2 approval (the
freeze), because nothing of the test set may be loaded before it.

Arm A: every star of samples/arm_a_test.csv, trials chosen as in round 1 (harness.make_defs_a).
Arm B: per Tmag bin the first 50 hosts of samples/arm_b_test_order.csv that have a light curve
       (a host without one is replaced by the next in the order), 4 injections each with numpy seed
       [20261202, tic] (harness.make_defs_b), and a null run per draw.

Inputs  : samples/arm_a_test.csv, samples/arm_b_test_order.csv, MAST light curves
Outputs : test/arm_a_defs.csv, arm_b_defs.csv, arm_a_rows.csv, arm_b_rows.csv (+ notes and errors)
Usage   : python calibration/mono/round2/run_test.py [--workers 12]   then   python score2.py
"""
import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
SEED_B, PER_BIN = 20261202, 50


def sh(*args):
    subprocess.run([sys.executable, "-W", "ignore", str(HERE / "harness.py"), *map(str, args)], check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--cache", default=str(HERE.parent / "work/lc_cache"))
    a = ap.parse_args()
    if "Step 2 of 2" not in (HERE.parent / "APPROVALS.md").read_text():
        sys.exit("The rule is not frozen: APPROVALS.md has no step-2 approval. Nothing was loaded.")
    out = HERE / "test"
    out.mkdir(exist_ok=True)
    common = ["--cache", a.cache, "--workers", a.workers]
    if not (out / "arm_a_defs.csv").exists():
        sh("--make-defs", "a", "--sample", HERE / "samples/arm_a_test.csv", "--out", out / "arm_a_defs.csv", *common)
    order = pd.read_csv(HERE / "samples/arm_b_test_order.csv")
    if not (out / "arm_b_defs.csv").exists():
        defs, notes, done = [], [], {b: 0 for b in order.tmag_bin.unique()}
        left = {b: g.sort_values("rank_in_bin") for b, g in order.groupby("tmag_bin")}
        rnd = 0
        while any(done[b] < PER_BIN and len(left[b]) for b in done):
            batch = []
            for b in done:
                need = PER_BIN - done[b]
                batch.append(left[b].head(need))
                left[b] = left[b].iloc[need:]
            batch = pd.concat(batch)
            batch.to_csv(out / f"arm_b_batch{rnd}.csv", index=False)
            sh("--make-defs", "b", "--sample", out / f"arm_b_batch{rnd}.csv", "--seed", SEED_B,
               "--out", out / f"arm_b_batch{rnd}_defs.csv", *common)
            d = pd.read_csv(out / f"arm_b_batch{rnd}_defs.csv") if (out / f"arm_b_batch{rnd}_defs.csv").stat().st_size > 5 else pd.DataFrame()
            notes.append(pd.read_csv(out / f"arm_b_batch{rnd}_defs.notes.csv").assign(batch=rnd))
            if len(d):
                defs.append(d)
                got = batch[batch.tic.isin(d.tic)].tmag_bin.value_counts()
                for b, n in got.items():
                    done[b] += int(n)
            rnd += 1
        pd.concat(defs).to_csv(out / "arm_b_defs.csv", index=False)
        pd.concat(notes).to_csv(out / "arm_b_hosts_tried.csv", index=False)
    for arm in ("a", "b"):
        if not (out / f"arm_{arm}_rows.csv").exists():
            sh("--defs", out / f"arm_{arm}_defs.csv", "--out", out / f"arm_{arm}_rows.csv", *common)


if __name__ == "__main__":
    main()
