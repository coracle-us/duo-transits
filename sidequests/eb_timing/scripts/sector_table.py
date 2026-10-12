"""
Per-sector observing time for TESS Sectors 1-107, measured from SPOC 2-min light curves
already on disk (data/raw/lightcurves/*/sNNN_SPOC_120s.csv.gz; read-only).

For each sector we take up to 5 files and keep the one with the most good cadences
(quality == 0, finite flux): t_start, t_end (BTJD) and good_days = N_good * 120 s.
Output: sidequests/eb_timing/data/sector_table.csv
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
LC = ROOT / "data" / "raw" / "lightcurves"
OUT = Path(__file__).resolve().parents[1] / "data" / "sector_table.csv"


def main():
    files = {}
    for f in LC.glob("*/s*_SPOC_120s.csv.gz"):
        s = int(re.match(r"s(\d+)_", f.name).group(1))
        if s <= 200:
            files.setdefault(s, []).append(f)
    rows = []
    for s in sorted(files):
        best = None
        for f in sorted(files[s], key=lambda p: -p.stat().st_size)[:5]:
            d = pd.read_csv(f, usecols=["time", "flux", "quality"])
            g = d[(d.quality == 0) & np.isfinite(d.flux) & np.isfinite(d.time)]
            if best is None or len(g) > best[0]:
                best = (len(g), d.time.min(), d.time.max())
        rows.append(dict(sector=s, t_start=best[1], t_end=best[2], good_days=best[0] * 120 / 86400))
    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False, float_format="%.4f")
    print(df.describe().T[["min", "50%", "max"]])


if __name__ == "__main__":
    main()
