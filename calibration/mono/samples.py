"""
Draw the two calibration samples exactly as the pre-registration states. No light curve is read.

Arm A: ExoFOP TOIs dispositioned CP or KP, period 10-100 d, not among the 83 stars of the
       2026-11-01 epoch-rule test, not a frozen host. Sorted by TOI, then up to 40 per depth bin
       (500-2000, 2000-5000, 5000-25000 ppm) with numpy seed 20261122.
Arm B: 150 frozen hosts, 50 per TESS-magnitude bin (<10, 10-12, >12), numpy seed 20261123,
       after sorting by TIC.

Inputs  : data/raw/exofop_toi.csv (ExoFOP, fetched 2026-11-22 by gaia_dr4/t16_overlap.py)
          salinas-second-transits/results/round2/epoch_rule_toi_sample.csv
          gaia_dr4/hosts_frozen.csv, gaia_dr4/results_dr3/ids.csv (Tmag from the TIC)
Outputs : calibration/mono/arm_a_sample.csv, arm_b_hosts.csv
Usage   : python calibration/mono/samples.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SEED_A, SEED_B = 20261122, 20261123


def arm_a():
    t = pd.read_csv(ROOT / "data/raw/exofop_toi.csv", low_memory=False)
    tuned = set(pd.read_csv(ROOT / "salinas-second-transits/results/round2/epoch_rule_toi_sample.csv").tic)
    hosts = set(pd.read_csv(ROOT / "gaia_dr4/hosts_frozen.csv").tic)
    n_toi = t.groupby("TIC ID").TOI.size()
    s = pd.DataFrame(dict(tic=t["TIC ID"], toi=t["TOI"], disposition=t["TFOPWG Disposition"],
                          period_d=t["Period (days)"], period_err=t["Period (days) err"],
                          t0=t["Epoch (BJD)"] - 2457000, duration_h=t["Duration (hours)"],
                          depth_ppm=t["Depth (ppm)"], depth_err=t["Depth (ppm) err"], tmag=t["TESS Mag"],
                          n_toi_on_star=t["TIC ID"].map(n_toi)))
    s = s[s.disposition.isin(["CP", "KP"]) & s.period_d.between(10, 100)]
    s = s.dropna(subset=["period_d", "t0", "duration_h", "depth_ppm"])
    s = s[~s.tic.isin(tuned) & ~s.tic.isin(hosts)].sort_values("toi").drop_duplicates("tic")
    s["depth_bin"] = pd.cut(s.depth_ppm, [500, 2000, 5000, 25000]).astype(str)
    s = s[s.depth_bin != "nan"]
    rng = np.random.default_rng(SEED_A)
    parts = []
    for _, g in s.groupby("depth_bin", sort=True):
        parts.append(g if len(g) <= 40 else g.iloc[np.sort(rng.choice(len(g), 40, replace=False))])
    return pd.concat(parts).sort_values("toi").reset_index(drop=True), s.depth_bin.value_counts().to_dict()


def arm_b():
    h = pd.read_csv(ROOT / "gaia_dr4/hosts_frozen.csv")
    ids = pd.read_csv(ROOT / "gaia_dr4/results_dr3/ids.csv", usecols=["tic", "Tmag"])
    h = h.merge(ids, on="tic", how="left").sort_values("tic")
    h["tmag_bin"] = pd.cut(h.Tmag, [-np.inf, 10, 12, np.inf], labels=["<10", "10-12", ">12"]).astype(str)
    h = h[h.tmag_bin != "nan"]
    rng = np.random.default_rng(SEED_B)
    parts = []
    for b in ["<10", "10-12", ">12"]:
        g = h[h.tmag_bin == b]
        parts.append(g if len(g) <= 50 else g.iloc[np.sort(rng.choice(len(g), 50, replace=False))])
    return pd.concat(parts).reset_index(drop=True), h.tmag_bin.value_counts().to_dict()


if __name__ == "__main__":
    a, na = arm_a()
    b, nb = arm_b()
    a.to_csv(HERE / "arm_a_sample.csv", index=False)
    b.to_csv(HERE / "arm_b_hosts.csv", index=False)
    print("arm A pool per depth bin", na, "-> drawn", a.depth_bin.value_counts().to_dict(), "total", len(a),
          "| stars with >1 TOI:", int((a.n_toi_on_star > 1).sum()))
    print("arm B pool per Tmag bin", nb, "-> drawn", b.tmag_bin.value_counts().to_dict(), "total", len(b))
