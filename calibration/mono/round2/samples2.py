"""
Draw the round-2 test samples exactly as prereg_round2.yaml (1.1) states. No light curve is read
and MAST is not queried.

Arm A: ExoFOP TOIs dispositioned CP or KP, period 5-100 d, depth 500-25000 ppm, with epoch,
       duration and depth; not among the 83 stars of the 2026-11-01 test, not in round 1's
       arm_a_sample.csv, not a frozen host. Sorted by TOI, up to 50 per depth bin, numpy seed 20261201.
Arm B: frozen hosts not drawn in round 1, without a stated period in candidates_frozen.csv, without
       an NSS orbit in gaia_dr4/results_dr3/report.csv, and without an eclipsing-binary class in
       singles/results/final_classification.csv. Sorted by TIC, then each Tmag bin is put in a random
       order with numpy seed 20261202. The first 50 of a bin are the draw; a host with no light curve
       is replaced by the next in that order (decided at run time, after the freeze).

Inputs  : data/raw/exofop_toi.csv (ExoFOP, fetched 2026-11-22), the files named above
Outputs : calibration/mono/round2/samples/arm_a_test.csv, arm_b_test_order.csv
Usage   : python calibration/mono/round2/samples2.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
SEED_A, SEED_B, PER_BIN = 20261201, 20261202, 50


def arm_a():
    t = pd.read_csv(ROOT / "data/raw/exofop_toi.csv", low_memory=False)
    seen = set(pd.read_csv(ROOT / "salinas-second-transits/results/round2/epoch_rule_toi_sample.csv").tic)
    seen |= set(pd.read_csv(ROOT / "calibration/mono/arm_a_sample.csv").tic)
    hosts = set(pd.read_csv(ROOT / "gaia_dr4/hosts_frozen.csv").tic)
    n_toi = t.groupby("TIC ID").TOI.size()
    s = pd.DataFrame(dict(tic=t["TIC ID"], toi=t["TOI"], disposition=t["TFOPWG Disposition"],
                          period_d=t["Period (days)"], period_err=t["Period (days) err"],
                          t0=t["Epoch (BJD)"] - 2457000, duration_h=t["Duration (hours)"],
                          depth_ppm=t["Depth (ppm)"], depth_err=t["Depth (ppm) err"], tmag=t["TESS Mag"],
                          n_toi_on_star=t["TIC ID"].map(n_toi)))
    s = s[s.disposition.isin(["CP", "KP"]) & s.period_d.between(5, 100)]
    s = s.dropna(subset=["period_d", "t0", "duration_h", "depth_ppm"])
    s = s[~s.tic.isin(seen) & ~s.tic.isin(hosts)].sort_values("toi").drop_duplicates("tic")
    s["depth_bin"] = pd.cut(s.depth_ppm, [500, 2000, 5000, 25000]).astype(str)
    s = s[s.depth_bin != "nan"]
    rng = np.random.default_rng(SEED_A)
    parts = []
    for _, g in s.groupby("depth_bin", sort=True):
        parts.append(g if len(g) <= PER_BIN else g.iloc[np.sort(rng.choice(len(g), PER_BIN, replace=False))])
    return pd.concat(parts).sort_values("toi").reset_index(drop=True), s


def arm_b():
    h = pd.read_csv(ROOT / "gaia_dr4/hosts_frozen.csv")
    ids = pd.read_csv(ROOT / "gaia_dr4/results_dr3/ids.csv", usecols=["tic", "Tmag"])
    h = h.merge(ids, on="tic", how="left").sort_values("tic")
    used = set(pd.read_csv(ROOT / "calibration/mono/arm_b_hosts.csv").tic)
    cand = pd.read_csv(ROOT / "gaia_dr4/candidates_frozen.csv")
    rep = pd.read_csv(ROOT / "gaia_dr4/results_dr3/report.csv")
    cls = pd.read_csv(ROOT / "singles/results/final_classification.csv")
    has_period = set(cand[cand.period_d.notna()].tic)
    has_nss = set(rep[rep.n_nss.fillna(0) > 0].tic)
    is_eb = set(cls[cls.cls.str.contains("EB", na=False)].TIC)
    h["excluded"] = np.select([h.tic.isin(used), h.tic.isin(has_period), h.tic.isin(has_nss), h.tic.isin(is_eb)],
                              ["round 1", "stated period", "NSS orbit", "EB class"], "")
    counts = h.excluded.value_counts().to_dict()
    h = h[h.excluded == ""].drop(columns="excluded")
    h["tmag_bin"] = pd.cut(h.Tmag, [-np.inf, 10, 12, np.inf], labels=["<10", "10-12", ">12"]).astype(str)
    h = h[h.tmag_bin != "nan"]
    rng = np.random.default_rng(SEED_B)
    parts = []
    for b in ["<10", "10-12", ">12"]:
        g = h[h.tmag_bin == b]
        g = g.iloc[rng.permutation(len(g))].copy()
        g["rank_in_bin"] = np.arange(1, len(g) + 1)
        parts.append(g)
    return pd.concat(parts).reset_index(drop=True), counts


if __name__ == "__main__":
    a, pool = arm_a()
    b, counts = arm_b()
    (HERE / "samples").mkdir(exist_ok=True)
    a.to_csv(HERE / "samples/arm_a_test.csv", index=False)
    b.to_csv(HERE / "samples/arm_b_test_order.csv", index=False)
    print("arm A pool", len(pool), "one-TOI", int((pool.n_toi_on_star == 1).sum()), "per bin",
          pool.depth_bin.value_counts().to_dict(), "-> drawn", len(a), a.depth_bin.value_counts().to_dict(),
          "| one-TOI drawn:", int((a.n_toi_on_star == 1).sum()))
    print("arm B hosts removed", counts, "-> eligible per bin", b.tmag_bin.value_counts().to_dict())
