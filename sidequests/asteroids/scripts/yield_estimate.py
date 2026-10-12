"""
Upper bound on new Trojan/Hilda periods obtainable from S56-S107 with this method.

Validation + pilot say the method works (auto_pass) only for V <~ 16.5 (Trojans/MBAs in clean
fields) and fails at V >~ 17.5. For every sector 56-107: members on silicon >= 10 d with predicted
median V < VLIM; union over sectors; joined with LCDB (U >= 2 known or not).
Output: data/yield_estimate.csv (one row per object) and a printed summary.
Usage: python yield_estimate.py [VLIM]
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import select_targets as st  # noqa: E402

HERE = Path(__file__).resolve().parents[1]


def main():
    vlim = float(sys.argv[1]) if len(sys.argv) > 1 else 16.5
    m = pd.read_csv(HERE / "cache" / "population_lcdb_merge.csv.gz", low_memory=False)
    m["pdes"] = m.number.astype("Int64").astype(str)
    umap = {"1": 1.0, "1+": 1.3, "2-": 1.7, "2": 2.0, "2+": 2.3, "3-": 2.7, "3": 3.0}
    m["Unum"] = m.U.map(umap)
    rows = []
    for pop in ["trojan", "hilda"]:
        counts = pd.read_csv(HERE / "data" / f"sector_counts_{pop}_H13.5.csv")
        for s in counts[counts.n_any > 0].sector:
            t = st.sector_table(pop, int(s), hmax=13.5)
            t = t[(t.days_on_silicon >= 10) & (t.V_med < vlim)]
            for r in t.itertuples():
                rows.append(dict(pop=pop, pdes=r.pdes, name=r.name, H=r.H, sector=int(s), V_med=r.V_med))
            print(pop, s, len(t), flush=True)
    d = pd.DataFrame(rows).merge(m[["pdes", "period_h", "U", "Unum"]], on="pdes", how="left")
    d.to_csv(HERE / "data" / "yield_estimate.csv", index=False)
    u = d.sort_values("V_med").drop_duplicates(["pop", "pdes"])
    for pop, g in u.groupby("pop"):
        print(f"{pop}: {len(g)} distinct objects with V<{vlim} and >=10 d in some sector 56-107; "
              f"LCDB U>=2: {(g.Unum >= 2).sum()}; without: {(~(g.Unum >= 2)).sum()}")


if __name__ == "__main__":
    main()
