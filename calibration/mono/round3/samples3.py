"""
Write the round-3 test samples exactly as prereg_round3.yaml (1.0) states. No light curve is read
and MAST is not queried.

Arm A: every star of round 2's arm-A pool (round2/samples2.arm_a) that was not drawn, plus every
       ExoFOP TOI dispositioned CP or KP with period above 100 d, depth 500-25000 ppm, with epoch,
       duration and depth, on a star not used in any round and not a frozen host. No draw.
Arm B: the order of round2/samples/arm_b_test_order.csv without the hosts tried in round 2
       (round2/test/arm_b_hosts_tried.csv). The hosts are taken in that order at run time.

Inputs  : data/raw/exofop_toi.csv (ExoFOP, fetched 2026-11-22), the files named above
Outputs : calibration/mono/round3/samples/arm_a_test.csv, arm_b_test_order.csv
Usage   : python calibration/mono/round3/samples3.py
"""
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / "round2"))
import samples2  # noqa: E402


def arm_a():
    drawn, pool = samples2.arm_a()
    left = pool[~pool.tic.isin(drawn.tic)]
    t = pd.read_csv(ROOT / "data/raw/exofop_toi.csv", low_memory=False)
    used = set(pool.tic)
    used |= set(pd.read_csv(ROOT / "salinas-second-transits/results/round2/epoch_rule_toi_sample.csv").tic)
    used |= set(pd.read_csv(ROOT / "calibration/mono/arm_a_sample.csv").tic)
    used |= set(pd.read_csv(ROOT / "gaia_dr4/hosts_frozen.csv").tic)
    n_toi = t.groupby("TIC ID").TOI.size()
    s = pd.DataFrame(dict(tic=t["TIC ID"], toi=t["TOI"], disposition=t["TFOPWG Disposition"],
                          period_d=t["Period (days)"], period_err=t["Period (days) err"],
                          t0=t["Epoch (BJD)"] - 2457000, duration_h=t["Duration (hours)"],
                          depth_ppm=t["Depth (ppm)"], depth_err=t["Depth (ppm) err"], tmag=t["TESS Mag"],
                          n_toi_on_star=t["TIC ID"].map(n_toi)))
    s = s[s.disposition.isin(["CP", "KP"]) & (s.period_d > 100) & s.depth_ppm.between(500, 25000)]
    s = s.dropna(subset=["period_d", "t0", "duration_h", "depth_ppm"])
    s = s[~s.tic.isin(used)].sort_values("toi").drop_duplicates("tic")
    return pd.concat([left.assign(part="round-2 pool, not drawn"), s.assign(part="period above 100 d")]) \
        .sort_values("toi").reset_index(drop=True)


def arm_b():
    order = pd.read_csv(HERE.parent / "round2/samples/arm_b_test_order.csv")
    tried = set(pd.read_csv(HERE.parent / "round2/test/arm_b_hosts_tried.csv").tic)
    return order[~order.tic.isin(tried)].reset_index(drop=True)


if __name__ == "__main__":
    a, b = arm_a(), arm_b()
    (HERE / "samples").mkdir(exist_ok=True)
    a.to_csv(HERE / "samples/arm_a_test.csv", index=False)
    b.to_csv(HERE / "samples/arm_b_test_order.csv", index=False)
    print("arm A", len(a), a.part.value_counts().to_dict(), "one-TOI", int((a.n_toi_on_star == 1).sum()))
    print("arm B left per bin", b.tmag_bin.value_counts().to_dict())
