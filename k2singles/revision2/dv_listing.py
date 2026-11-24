"""
AAS82225 revision: does MAST hold a SPOC data-validation (DV) report, i.e. a SPOC TCE, for each star whose period the
note calls new? SPOC writes a DV report only for a target with a TCE. Listing only; nothing is downloaded.
Output: results/literature/dv_listing.txt
Usage:  python k2singles/revision2/dv_listing.py
"""
import os

import pandas as pd
from astroquery.mast import Observations

HERE = os.path.dirname(os.path.abspath(__file__))
NEW = [246194159, 248811085, 246089124, 246696483, 211390677, 247967714, 212012030, 211503363]
t = pd.read_csv(os.path.join(HERE, "..", "rnaas", "table.csv"))
out = []
for e in NEW:
    tic = int(t.loc[t.EPIC == e, "TIC"].iloc[0])
    obs = Observations.query_criteria(target_name=str(tic), obs_collection="TESS", dataproduct_type="timeseries")
    if len(obs) == 0:
        out.append(f"EPIC {e} TIC {tic}: no SPOC time series on MAST, so no SPOC TCE or DV report")
        continue
    p = Observations.get_product_list(obs)
    dv = [f for f in p["productFilename"] if f.endswith("_dvr.xml")]
    secs = sorted({int(o["sequence_number"]) for o in obs if o["sequence_number"]})
    out.append(f"EPIC {e} TIC {tic}: SPOC sectors {secs}; DV reports: {dv if dv else 'none'}")
print("\n".join(out))
open(os.path.join(HERE, "results", "literature", "dv_listing.txt"), "w").write("\n".join(out) + "\n")
