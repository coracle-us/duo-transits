"""
AAS82225 revision: does MAST hold a SPOC data-validation (DV) report, i.e. a SPOC TCE, for each Table 1 star? SPOC writes a
DV report only for a target with a TCE. Both collections are listed: TESS (SPOC 2-min) and HLSP (TESS-SPOC FFI;
file names hlsp_tess-spoc_...). The first version (2026-11-24) listed only the TESS collection and missed two
TESS-SPOC reports (review 2, PLAN.md). Listing only; nothing is downloaded.
Output: results/literature/dv_listing.txt
Usage:  python k2singles/revision2/dv_listing.py
"""
import os

import pandas as pd
from astroquery.mast import Observations

HERE = os.path.dirname(os.path.abspath(__file__))
t = pd.read_csv(os.path.join(HERE, "..", "rnaas", "table.csv"))
NEW = list(t.EPIC)
out = []
for e in NEW:
    tic = int(t.loc[t.EPIC == e, "TIC"].iloc[0])
    obs = Observations.query_criteria(target_name=[str(tic), f"{tic:016d}"], dataproduct_type="timeseries")
    if len(obs) == 0:
        out.append(f"EPIC {e} TIC {tic}: no time series on MAST")
        continue
    p = Observations.get_product_list(obs)
    dv = [f"{f} ({int(z)} B)" for f, z in zip(p["productFilename"], p["size"]) if f.endswith("dvr.xml")]
    out.append(f"EPIC {e} TIC {tic} ({t.loc[t.EPIC == e, 'group'].iloc[0]}): collections "
               f"{sorted(set(map(str, obs['obs_collection'])))}; DV reports: {'; '.join(dv) if dv else 'none'}")
print("\n".join(out))
open(os.path.join(HERE, "results", "literature", "dv_listing.txt"), "w").write("\n".join(out) + "\n")
