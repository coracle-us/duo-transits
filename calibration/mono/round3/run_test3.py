"""
Run the round-3 test set. Refuses to start until APPROVALS.md records Jon's step-2 approval of
round 3 (the freeze), because nothing of the test set may be loaded before it.

Arm A: every star of samples/arm_a_test.csv, trials chosen as in rounds 1 and 2 (make_defs_a).
Arm B: per Tmag bin the next 100 hosts of samples/arm_b_test_order.csv that have a light curve (in
       the bin below Tmag 10, all that have one); 4 injections each with numpy seed [20261203, tic];
       catalogued depth and duration with a 10% Gaussian error, numpy seed [20261204, tic]; a null
       run per draw.

Inputs  : samples/arm_a_test.csv, samples/arm_b_test_order.csv, MAST light curves
Outputs : test/arm_a_defs.csv, arm_b_defs.csv, arm_a_rows.csv, arm_b_rows.csv (+ notes and errors)
Usage   : python calibration/mono/round3/run_test3.py [--workers 12]   then   python score3.py
"""
import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
SEED_B, SEED_NOISE, NOISE, PER_BIN = 20261203, 20261204, 0.1, 100


def sh(*args):
    subprocess.run([sys.executable, "-W", "ignore", str(HERE / "harness3.py"), *map(str, args)], check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args()
    text = (HERE.parent / "APPROVALS.md").read_text()
    if "## Round 3" not in text or "Step 2 of 2" not in text[text.index("## Round 3"):]:
        sys.exit("The round-3 rule is not frozen: APPROVALS.md has no step-2 approval. Nothing was loaded.")
    out = HERE / "test"
    out.mkdir(exist_ok=True)
    common = ["--workers", a.workers]
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
            f = out / f"arm_b_batch{rnd}"
            batch.to_csv(f"{f}.csv", index=False)
            sh("--make-defs", "b", "--sample", f"{f}.csv", "--seed", SEED_B, "--noise", NOISE,
               "--noise-seed", SEED_NOISE, "--out", f"{f}_defs.csv", *common)
            d = pd.read_csv(f"{f}_defs.csv") if Path(f"{f}_defs.csv").stat().st_size > 5 else pd.DataFrame()
            notes.append(pd.read_csv(f"{f}_defs.notes.csv").assign(batch=rnd))
            if len(d):
                defs.append(d)
                for b, n in batch[batch.tic.isin(d.tic)].tmag_bin.value_counts().items():
                    done[b] += int(n)
            rnd += 1
        pd.concat(defs).to_csv(out / "arm_b_defs.csv", index=False)
        pd.concat(notes).to_csv(out / "arm_b_hosts_tried.csv", index=False)
    for arm in ("a", "b"):
        if not (out / f"arm_{arm}_rows.csv").exists():
            sh("--defs", out / f"arm_{arm}_defs.csv", "--out", out / f"arm_{arm}_rows.csv", *common)


if __name__ == "__main__":
    main()
