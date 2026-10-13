"""
Count the CTOIs removed by the 'actively worked' owner exclusion in singles/scripts/build_ctoi_singles.py, for
editor comment 3 on AAS82207. It replicates the script's masks and reports, among single/duo, open CTOIs, how many
belong to each excluded ExoFOP user.
Usage:  python singles/revision/count_owner_exclusion.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ACTIVE_OWNERS = {"essack", "dragomir", "skinner", "hawthorn", "salinas"}

c = pd.read_csv(ROOT / "data" / "raw" / "exofop_ctoi.csv")
notes = (c["Notes"].astype(str) + " " + c["Tag"].astype(str)).str.lower()
per = pd.to_numeric(c["Period (days)"], errors="coerce")
is_single = per.isna() | (per <= 0) | notes.str.contains("single|duo|mono")
is_open = c["Promoted to TOI"].isna() & ~c["TFOPWG Disposition"].astype(str).isin(["FP", "FA", "CP", "KP"])
user = c["User"].astype(str).str.lower()
pool = c[is_single & is_open]
excl = pool[user.loc[pool.index].isin(ACTIVE_OWNERS)]
print(f"single/duo open CTOIs: {len(pool)}; excluded as actively worked: {len(excl)}; kept: {len(pool) - len(excl)}")
print(excl["User"].str.lower().value_counts().to_string())
