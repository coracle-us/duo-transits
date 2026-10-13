"""
Second, out-of-sample NCVZ false-alarm set (fixed BEFORE timing; the first set and its
generator ncvz_sets.py are left untouched so data/ncvz_fa_list.csv stays reproducible).

Footprint: the first set used ecliptic latitude > 78 deg; only 16 such stars remained eligible.
This set uses Mitnyan+2024's own input criterion for the B26 sample ("in and near the NCVZ"):
>= 8 sectors in Year 2 (S14-26), from data/coverage_all.csv. Cuts identical to the first
set: Prsa+2022 EB, 0.4 <= P <= 30 d, 0 <= morph < 0.6, Dp-2g >= 0.01, Tmag <= 15.
Exclusions:
  - the first FA set (data/ncvz_fa_list.csv) and the recall sets (data/ncvz_recall_list.csv);
  - data/known_multiples.csv;
  - every TIC in data/b26_triples.csv (all 168, parsed table by table, superscripts handled);
  - every 6-10 digit number in the LaTeX of Borkovits+2026 and Mitnyan+2024 that is a TIC in
    the Prsa catalogue (catches the 15 superscripted TICs known_triples.py missed, and any
    TIC named in the text, e.g. the long-term-variation systems).
Selection: all eligible stars if <= 50, else 50 drawn with numpy default_rng(SEED).
Output: data/ncvz_fa2_list.csv, data/ncvz_fa2_ephem.csv
"""
import io
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
D = HERE / "data"
SEED = 20261013


def latex_numbers(aid):
    with tarfile.open(fileobj=io.BytesIO((HERE / "cache" / "arxiv" / f"{aid}.bin").read_bytes())) as tf:
        txt = "\n".join(tf.extractfile(m).read().decode("utf-8", "replace") for m in tf.getmembers()
                        if m.isfile() and m.name.endswith(".tex"))
    return {int(x) for x in re.findall(r"(?<!\d)(\d{6,10})(?!\d)", txt)}


def main():
    a = pd.read_csv(D / "sample_all.csv")
    cov = pd.read_csv(D / "coverage_all.csv").fillna("").set_index("TIC").sectors
    a["n_y2"] = a.TIC.map(cov.apply(lambda s: sum(1 for x in s.split(";") if x and 14 <= int(x) <= 26)))
    cut = (a.Per.between(0.4, 30) & (a.Morph >= 0) & (a.Morph < 0.6) & (a["Dp-2g"] >= 0.01) & (a.Tmag <= 15)
           & (a.n_y2 >= 8))
    excl = set(pd.read_csv(D / "ncvz_fa_list.csv", comment="#").TIC)
    excl |= set(pd.read_csv(D / "ncvz_recall_list.csv").TIC)
    excl |= set(pd.read_csv(D / "known_multiples.csv").TIC)
    excl |= set(pd.read_csv(D / "b26_triples.csv").tic)
    tex = latex_numbers("2604.04004") | latex_numbers("2402.01486")
    excl |= tex & set(a.TIC)
    pool = a[cut & ~a.TIC.isin(excl)].sort_values("TIC").reset_index(drop=True)
    if len(pool) > 50:
        pool = pool.loc[np.sort(np.random.default_rng(SEED).choice(len(pool), 50, replace=False))]
    fa = pool[["TIC", "Per", "BJD0", "Tmag", "Morph", "Dp-2g", "ecl_lat", "n_y2"]]
    n_cut = int(cut.sum())
    with open(D / "ncvz_fa2_list.csv", "w") as fh:
        fh.write(f"# NCVZ false-alarm set 2: {len(fa)} stars; footprint >= 8 Year-2 sectors; {n_cut} pass the cuts, "
                 f"{len(excl & set(a[cut].TIC))} excluded; seed {SEED} (used only if > 50 eligible); fixed before timing\n")
        fa.to_csv(fh, index=False)
    fa[["TIC", "Per", "BJD0"]].to_csv(D / "ncvz_fa2_ephem.csv", index=False)
    print(f"pass cuts {n_cut}; excluded {len(excl & set(a[cut].TIC))}; selected {len(fa)}; "
          f"beta>78: {(fa.ecl_lat > 78).sum()}")
    print(f"B26/Mitnyan LaTeX TICs in Prsa: {len(tex & set(a.TIC))}; of those not in known_multiples: "
          f"{len((tex & set(a.TIC)) - set(pd.read_csv(D / 'known_multiples.csv').TIC))}")


if __name__ == "__main__":
    main()
