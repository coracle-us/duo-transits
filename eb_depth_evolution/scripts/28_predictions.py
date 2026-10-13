"""
Predicted primary eclipses of EPIC 202126877 for the 2026-27 season, as BJD_TDB and as geocentric UTC (with the
barycentric light-travel time removed). Ephemeris: the linear fit to the TESS primaries (REPORT.txt section 3,
13_oc_basic.py): P = 21.774799 d, T0 = BTJD 2488.58926. Allow +-30 min for the observed eclipse-timing variations.
Output: eb_depth_evolution/predictions_2026_27.csv
"""
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

P, T0 = 21.774799, 2488.58926
COORD = SkyCoord(97.98012 * u.deg, 25.48207 * u.deg)
GEO = EarthLocation.from_geocentric(0 * u.m, 0 * u.m, 0 * u.m)


def main():
    start = Time("2026-10-13").tdb.jd - 2457000
    n0 = int(np.ceil((start - T0) / P))
    rows = []
    for n in range(n0, n0 + 7):
        t = T0 + n * P
        tb = Time(t + 2457000, format="jd", scale="tdb", location=GEO)
        lt = tb.light_travel_time(COORD, kind="barycentric")
        rows.append(dict(cycle=n, BJD_TDB=round(t + 2457000, 4), utc_geocentric=(tb - lt).utc.iso[:16],
                         ltt_min=round(float(lt.to(u.min).value), 1)))
    df = pd.DataFrame(rows)
    df.to_csv(Path(__file__).resolve().parents[1] / "predictions_2026_27.csv", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
