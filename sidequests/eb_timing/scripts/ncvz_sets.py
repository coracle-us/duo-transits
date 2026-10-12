"""
Northern CVZ test sets, fixed BEFORE any timing is run on them.

Recall set: all 35 Borkovits+2026 group-L1 (certain pure-LTTE) systems, from
  data/b26_triples.csv (solution 1). Seed ephemeris: B26 P1 (only 5 of the 35 are in the
  Prsa+2022 catalogue); the epoch is arbitrary (run_etv moves phase 0 to the template minimum).
False-alarm set: Prsa+2022 EBs in the NCVZ (ecliptic latitude > 78 deg, inside the B26 'in
  and near the NCVZ' footprint) passing the same cuts as the southern sample
  (0.4 <= P <= 30 d, 0 <= morph < 0.6, Dp-2g >= 0.01) plus Tmag <= 15, and NOT in
  data/known_multiples.csv (which contains every TIC named in Mitnyan+2024 and B26, i.e.
  their 135/168 candidates and long-term-variation systems, plus 25 other multiples papers).
  50 drawn with numpy default_rng(20261012).
Baseline: S14-86 (B26), enforced by fetching/loading sectors <= 86.

Output: data/ncvz_recall_list.csv, data/ncvz_fa_list.csv, data/ncvz_ephem.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd

D = Path(__file__).resolve().parents[1] / "data"
SEED = 20261012


def main():
    b = pd.read_csv(D / "b26_triples.csv")
    L1 = b[(b.group == "L1") & (b.solution == 1)].copy()
    rec = L1[["tic", "P1", "P2", "P2_err", "A_ltte_s", "A_ltte_err_s", "e2", "e2_err"]].rename(columns={"tic": "TIC"})
    rec.to_csv(D / "ncvz_recall_list.csv", index=False)

    a = pd.read_csv(D / "sample_all.csv")
    known = set(pd.read_csv(D / "known_multiples.csv").TIC)
    pool = a[(a.ecl_lat > 78) & a.Per.between(0.4, 30) & (a.Morph >= 0) & (a.Morph < 0.6)
             & (a["Dp-2g"] >= 0.01) & (a.Tmag <= 15) & ~a.TIC.isin(known) & ~a.TIC.isin(L1.tic)]
    pool = pool.sort_values("TIC").reset_index(drop=True)
    pick = np.sort(np.random.default_rng(SEED).choice(len(pool), size=min(50, len(pool)), replace=False))
    fa = pool.loc[pick, ["TIC", "Per", "BJD0", "Tmag", "Morph", "Dp-2g", "ecl_lat"]]
    with open(D / "ncvz_fa_list.csv", "w") as fh:
        fh.write(f"# NCVZ false-alarm set: {len(fa)} of {len(pool)} eligible, rng seed {SEED}; fixed before timing\n")
        fa.to_csv(fh, index=False)
    # ephemeris seeds for run_etv --ephem (Prsa where available, else B26 P1 with epoch 2000)
    eph = pd.concat([fa[["TIC", "Per", "BJD0"]],
                     pd.DataFrame({"TIC": rec.TIC, "Per": rec.P1, "BJD0": 2000.0})])
    cat = a.set_index("TIC")
    in_cat = eph.TIC.isin(cat.index)
    eph.loc[in_cat, "Per"] = cat.loc[eph.TIC[in_cat], "Per"].values
    eph.loc[in_cat, "BJD0"] = cat.loc[eph.TIC[in_cat], "BJD0"].values
    eph.drop_duplicates("TIC").to_csv(D / "ncvz_ephem.csv", index=False)
    print(f"recall {len(rec)}; FA pool {len(pool)} -> {len(fa)}; ephem rows {eph.TIC.nunique()}")


if __name__ == "__main__":
    main()
