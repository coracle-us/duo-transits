"""False-alarm check: run the identical detector on the sign-flipped residuals
(negative excursions) of every host-sector. Real flares only go up, so the
negative-side event rate estimates how many 'flares' noise/systematics produce.
Transits and other dips also land here, so this is an upper bound on the noise rate.
Writes false_alarm.csv.
"""
import numpy as np
import pandas as pd

import flarelib as F


def main():
    cat = F.load_catalog()
    sectors = pd.read_csv(F.SQ / "sectors.csv")
    rows = []
    for _, s in sectors.iterrows():
        lc = F.detrend(F.load_sector(int(s.tic), int(s.sector), cat))
        neg = dict(lc)
        neg["resid"] = -lc["resid"]
        fl = F.find_flares(neg)
        if len(fl):
            fl = F.auto_flags(fl, neg)
            rows.append(fl)
    fa = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    fa.to_csv(F.SQ / "false_alarm.csv", index=False)
    tot = sectors.obs_days.sum()
    passing = fa[(~fa.near_dump) & (~fa.edge) & fa.shape_ok] if len(fa) else fa
    print(f"total obs {tot:.1f} d; negative events {len(fa)}; passing auto flags {len(passing)}")
    if len(passing):
        print(passing[["tic", "sector", "t_peak", "amp_sig", "n_above", "dur_min", "ed_s"]].to_string())


if __name__ == "__main__":
    main()
