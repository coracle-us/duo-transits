"""
Pilot sample: numbered Jupiter Trojans crossing TESS Sector 72 (200-s FFIs; after the Konkoly
EDR2 coverage, which ends at S48; the sector with the most Trojans on silicon in S56-S107, see
data/sector_counts_trojan_H13.5.csv).

Rule (fixed before any pilot photometry):
  - select_targets.sector_table('trojan', 72): >= 15 d on silicon (Keplerian pre-filter).
  - The 20 brightest by predicted median V, plus 10 more spread evenly in predicted V over
    18.0 <= V < 19.0 (every k-th object when sorted by V), to map the faint limit.
  - LCDB status is NOT used for selection (it is joined afterwards for the comparison).
Output: data/pilot_sample.csv
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import select_targets as st  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
SECTOR = 72


def main():
    t = st.sector_table("trojan", SECTOR)
    t = t[t.days_on_silicon >= 15].sort_values("V_med").reset_index(drop=True)
    bright = t.iloc[:20].copy()
    bright["tier"] = "brightest20"
    faint = t[(t.V_med >= 18.0) & (t.V_med < 19.0)]
    idx = np.linspace(0, len(faint) - 1, 10).round().astype(int)
    faint = faint.iloc[idx].copy()
    faint["tier"] = "V18-19"
    p = pd.concat([bright, faint])
    m = pd.read_csv(HERE / "cache" / "population_lcdb_merge.csv.gz", low_memory=False)
    m["pdes"] = m.number.astype("Int64").astype(str)
    p = p.merge(m[["pdes", "period_h", "amax", "U"]].rename(columns={"period_h": "P_lcdb_h", "amax": "amax_lcdb",
                                                                     "U": "U_lcdb"}), on="pdes", how="left")
    p["sector"] = SECTOR
    p.to_csv(HERE / "data" / "pilot_sample.csv", index=False)
    print(len(t), "Trojans with >=15 d in S72")
    print(p[["pdes", "name", "H", "V_med", "days_on_silicon", "camera", "alpha", "tier", "P_lcdb_h", "U_lcdb"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
