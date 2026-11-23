"""Next transits of the four periodic Salinas stars, from the round 3 ephemerides.

Usage: python round3/s15_predictions.py
Reads results/round3/oddeven_summary.csv; writes results/round3/predictions_2026-11-22.csv
(the next 3 transits after 2026-11-22 for each star, with 1-sigma timing errors from P and T0).
Times: BTJD = BJD_TDB - 2457000; UTC from astropy (TDB -> UTC; the light-travel term is already
in BJD, so the UTC column is the barycentric time on the UTC scale, good to ~8 min for planning).
"""
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time

HERE = Path(__file__).resolve().parents[1]
SUM = HERE / "results" / "round3" / "oddeven_summary.csv"
OUT = HERE / "results" / "round3" / "predictions_2026-11-22.csv"
START = Time("2026-11-22T00:00:00", scale="utc").tdb.jd - 2457000

rows = []
for r in pd.read_csv(SUM).itertuples():
    e0 = int(np.ceil((START - r.T0_btjd) / r.P_d))
    for e in range(e0, e0 + 3):
        tc = r.T0_btjd + e * r.P_d
        err_d = np.hypot(r.T0_err_d, e * r.P_err_d)
        utc = Time(tc + 2457000, format="jd", scale="tdb").utc.iso[:16]
        rows.append(dict(tic=r.tic, epoch=e, tc_btjd=round(tc, 5), tc_err_min=round(err_d * 1440, 1),
                         tc_utc=utc, t14_h=round(r.t14_h, 2), depth_ppm=round(r.depth_ppm)))
df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)
print(df.to_string(index=False))
