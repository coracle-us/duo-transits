"""
Periods stated elsewhere, to be compared with the run's result afterwards (never an input to it).

Sources per scored row:
  stated        the source catalogue's own period (defs.csv, stated_period_d)
  v2_solved     data/catalog/status.csv solved Hawthorn period (validation/v2_truth.csv)
  v3_published  arXiv:2611.03019, a period published elsewhere (validation/v3_truth.csv)
  v3_best       arXiv:2611.03019, their highest-weight period, state RESOLVED or EXCLUSION
  own           a period this project derived earlier (own_ctoi, own_salinas rows on the same star)

Inputs  : catalogue/defs.csv, catalogue/validation/v2_truth.csv, v3_truth.csv, gaia_dr4/candidates_frozen.csv
Outputs : catalogue/check_periods.csv
Usage   : python catalogue/build_check_periods.py
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent


def main():
    d = pd.read_csv(HERE / "defs.csv")
    parts = [d[d.stated_period_d > 0][["row_id", "stated_period_d"]].rename(
        columns=dict(stated_period_d="period_d")).assign(source="stated")]
    v2 = pd.read_csv(HERE / "validation/v2_truth.csv")
    parts.append(v2[["row_id", "P_true_d"]].rename(columns=dict(P_true_d="period_d")).assign(source="v2_solved"))
    v3 = pd.read_csv(HERE / "validation/v3_truth.csv")
    pub = v3[v3.P_published_d.notna()]
    parts.append(pub[["row_id", "P_published_d"]].rename(columns=dict(P_published_d="period_d")).assign(source="v3_published"))
    best = v3[v3.State.isin(["RESOLVED", "EXCLUSION"]) & v3.P_best_A_d.notna()]
    parts.append(best[["row_id", "P_best_A_d"]].rename(columns=dict(P_best_A_d="period_d")).assign(source="v3_best"))
    c = pd.read_csv(HERE.parent / "gaia_dr4/candidates_frozen.csv")
    own = c[c.catalogue.str.startswith("own") & c.period_d.notna()][["tic", "period_d"]]
    parts.append(d[["row_id", "tic"]].merge(own, on="tic")[["row_id", "period_d"]].assign(source="own"))
    out = pd.concat(parts, ignore_index=True)[["row_id", "source", "period_d"]].drop_duplicates()
    out.to_csv(HERE / "check_periods.csv", index=False)
    print(out.source.value_counts().to_dict())


if __name__ == "__main__":
    main()
