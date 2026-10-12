"""
Predicted TESS FFI coverage (Sectors 1-107) for every EB in the Prsa et al. 2022 catalog,
using tess-point (on-silicon test with SPOC collateral margins, no aberration).

Output: sidequests/eb_timing/data/coverage_all.csv  (TIC, sectors as ';'-joined string)
"""
from pathlib import Path

import numpy as np
import pandas as pd
from tess_stars2px import tess_stars2px_function_entry as t2p

D = Path(__file__).resolve().parents[1] / "data"
MAX_SECTOR = 107


def main():
    cat = pd.read_csv(D / "prsa2022_tess_ebs.csv").drop_duplicates("TIC")
    out = t2p(cat.TIC.values, cat.RAJ2000.values, cat.DEJ2000.values)
    tic, sec = out[0], out[3]
    ok = (tic > 0) & (sec <= MAX_SECTOR)
    df = pd.DataFrame({"TIC": tic[ok], "sector": sec[ok]}).drop_duplicates()
    g = df.groupby("TIC").sector.apply(lambda s: ";".join(str(x) for x in sorted(s)))
    res = pd.DataFrame({"TIC": cat.TIC.values})
    res["sectors"] = res.TIC.map(g).fillna("")
    res.to_csv(D / "coverage_all.csv", index=False)
    print(f"{len(res)} EBs; with >=1 sector: {(res.sectors != '').sum()}")


if __name__ == "__main__":
    main()
