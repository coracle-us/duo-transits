"""Stage 0: load the Jia et al. (2026) BeT catalogue (China-VO r101801) to parquet.

Source files (downloaded by hand with curl, no login needed, 2026-11-22):
  https://nadc.china-vo.org/res/r101801/  ->  bet_table.csv, wise_lc.zip, bn_tess_data.zip
into ~/astro-cache/be-star-pulsations/catalogue/.

bet_table.csv columns: DR3Name, name (Gaia DR3 source_id), TIC, ra, dec, color_ptp
((W1-W2) peak-to-peak), w1_ptp, r_value (Pearson r of W1 vs W1-W2), Gmag, BPmag, RPmag,
P, Type (0 = TESS-checked, not EB; 1 EA, 2 EB, 3 EW; blank = not in their TESS table),
best_author (TESS pipeline they used). No outburst/decay dates or durations.
"""
import pandas as pd
from common import CAT_DIR, RESULTS

def load():
    d = pd.read_csv(CAT_DIR / "bet_table.csv")
    out = pd.DataFrame({
        "gaia_dr3_id": d["name"].astype("int64"),
        "tic": d["TIC"].astype("Int64"),
        "ra": d["ra"], "dec": d["dec"],
        "w1_ptp": d["w1_ptp"], "w1w2_ptp": d["color_ptp"],
        "pearson_r": d["r_value"],
        "gmag": d["Gmag"], "bpmag": d["BPmag"], "rpmag": d["RPmag"],
        "period": d["P"], "eb_type": d["Type"].astype("Int64"),
        "tess_author": d["best_author"],
    })
    return out

if __name__ == "__main__":
    c = load()
    assert len(c) == 3300, len(c)
    print("rows", len(c), "with TIC", c.tic.notna().sum(), "unique TIC", c.tic.nunique())
    print("eb_type", c.eb_type.value_counts(dropna=False).to_dict())
    print("EBs (type 1-3):", c.eb_type.isin([1, 2, 3]).sum())
    print("pearson_r range", c.pearson_r.min(), c.pearson_r.max())
    c.to_parquet(RESULTS / "catalogue.parquet", index=False)
    print("wrote", RESULTS / "catalogue.parquet")
