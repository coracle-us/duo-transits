"""
Measured MAST light-curve availability versus TESS magnitude, to set the magnitude cut of the
southern run.

For every Prsa-catalogue star we fetched (cache/lc/<TIC>/search_results.csv lists the MAST
products found: SPOC 120 s, TESS-SPOC, QLP), compare the sectors with at least one product
against the sectors on silicon from tess-point (data/coverage_all.csv), up to the maximum
sector of that fetch (107 south, 86 for the NCVZ test sets).
Output: data/lc_availability.csv (per star) and a binned summary printed / written to
data/lc_availability_bins.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
D = HERE / "data"


def main():
    a = pd.read_csv(D / "sample_all.csv").set_index("TIC")
    cov = pd.read_csv(D / "coverage_all.csv").fillna("").set_index("TIC").sectors
    rows = []
    for d in sorted((HERE / "cache" / "lc").iterdir()):
        tic = int(d.name)
        f = d / "search_results.csv"
        if tic not in a.index or not f.exists():
            continue
        s = pd.read_csv(f)
        smax = 86 if (a.loc[tic, "ecl_lat"] > 0) else 107
        pred = {int(x) for x in cov.loc[tic].split(";") if x and int(x) <= smax}
        have = set(s[s.author.isin(["SPOC", "TESS-SPOC", "QLP"])].sector.astype(int)) & pred
        have120 = set(s[(s.author == "SPOC") & (s.exp == 120)].sector.astype(int)) & pred
        rows.append(dict(TIC=tic, Tmag=a.loc[tic, "Tmag"], n_pred=len(pred), n_have=len(have),
                         n_spoc120=len(have120), frac=len(have) / max(len(pred), 1)))
    df = pd.DataFrame(rows)
    df.to_csv(D / "lc_availability.csv", index=False)
    bins = [0, 9, 10, 11, 12, 12.5, 13, 13.5, 14, 14.5, 16]
    df["bin"] = pd.cut(df.Tmag, bins)
    g = df.groupby("bin", observed=True).agg(n=("TIC", "size"), median_frac=("frac", "median"),
                                             min_frac=("frac", "min"), frac_ge_0p8=("frac", lambda x: np.mean(x >= 0.8)))
    g.to_csv(D / "lc_availability_bins.csv")
    print(g.round(2).to_string())


if __name__ == "__main__":
    main()
