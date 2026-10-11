"""
Single-transit project, step 1: build the target list of stale community single-transit
candidates (CTOIs) from the ExoFOP CTOI table.

Keep a CTOI if:
  - not promoted to a TOI, and TFOPWG disposition is not FP/FA/CP/KP
  - single or duo: period missing/0, or notes/tag mention single/duo/mono
  - not owned by a group that is actively working its singles (TSTPC: essack, dragomir,
    skinner; plus hawthorn, salinas) -- see docs/problem_set_survey.md
  - has an epoch, depth and duration
Drop (with a recorded reason) if the star:
  - hosts a confirmed planet in NASA pscomppars (matched by TIC)
  - has any TOI (the TOI process is handling it)
  - is a known TESS eclipsing binary (Prsa et al. 2022, ApJS 258, 16; VizieR J/ApJS/258/16,
    Sectors 1-26; saved as data/raw/tess_eb_prsa2022.tsv)

Inputs : data/raw/exofop_ctoi.csv, data/raw/exofop_toi.csv, data/raw/nasa_pscomppars.csv
Outputs: singles/data/ctoi_singles.csv  (kept), singles/data/ctoi_singles_dropped.csv

Usage:  python singles/scripts/build_ctoi_singles.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "singles" / "data"
ACTIVE_OWNERS = {"essack", "dragomir", "skinner", "hawthorn", "salinas"}


def main():
    c = pd.read_csv(RAW / "exofop_ctoi.csv")
    toi = pd.read_csv(RAW / "exofop_toi.csv")
    ps = pd.read_csv(RAW / "nasa_pscomppars.csv")
    ps_tic = set(pd.to_numeric(ps.tic_id.str.replace("TIC ", ""), errors="coerce").dropna().astype(int))
    toi_tic = set(toi["TIC ID"].astype(int))
    eb = pd.read_csv(RAW / "tess_eb_prsa2022.tsv", sep="\t", comment="#", skiprows=[1, 2],
                     dtype=str).dropna(subset=["TIC"])
    eb = eb[eb.TIC.str.strip().str.isdigit()]
    eb_per = dict(zip(eb.TIC.str.strip().astype(int), pd.to_numeric(eb.Per, errors="coerce")))

    notes = (c["Notes"].astype(str) + " " + c["Tag"].astype(str)).str.lower()
    per = pd.to_numeric(c["Period (days)"], errors="coerce")
    is_single = per.isna() | (per <= 0) | notes.str.contains("single|duo|mono")
    is_open = c["Promoted to TOI"].isna() & ~c["TFOPWG Disposition"].astype(str).isin(["FP", "FA", "CP", "KP"])
    owner_active = c["User"].astype(str).str.lower().isin(ACTIVE_OWNERS)
    has_params = (c["Transit Epoch (BJD)"].notna() & (c["Depth ppm"] > 0) & (c["Duration (hrs)"] > 0))
    base = c[is_single & is_open & ~owner_active].copy()
    print(f"single/duo, open, not actively owned: {len(base)} CTOIs on {base['TIC ID'].nunique()} TICs")

    base["drop_reason"] = ""
    base.loc[~has_params.loc[base.index], "drop_reason"] = "missing epoch/depth/duration"
    base.loc[base["TIC ID"].isin(ps_tic) & (base.drop_reason == ""), "drop_reason"] = "host has confirmed planet (pscomppars)"
    base.loc[base["TIC ID"].isin(toi_tic) & (base.drop_reason == ""), "drop_reason"] = "star has a TOI"
    is_eb = base["TIC ID"].isin(eb_per.keys()) & (base.drop_reason == "")
    base.loc[is_eb, "drop_reason"] = [f"known TESS EB (Prsa+2022), P={eb_per[t]:.4f} d"
                                      for t in base.loc[is_eb, "TIC ID"]]
    kept = base[base.drop_reason == ""].copy()
    dropped = base[base.drop_reason != ""]

    kept["epoch_btjd"] = kept["Transit Epoch (BJD)"] - 2457000
    kept = kept.rename(columns={"TIC ID": "TIC", "TESS Mag": "Tmag", "Depth ppm": "depth_ppm",
                                "Duration (hrs)": "dur_h", "Period (days)": "period_listed"})
    cols = ["TIC", "CTOI", "User", "Tmag", "RA", "Dec", "epoch_btjd", "depth_ppm", "dur_h",
            "period_listed", "Planet Radius (R_Earth)", "Stellar Radius (R_Sun)",
            "Stellar Eff Temp (K)", "CTOI lastmod", "Notes"]
    kept = kept[cols].sort_values("Tmag")
    OUT.mkdir(parents=True, exist_ok=True)
    kept.to_csv(OUT / "ctoi_singles.csv", index=False)
    dropped.to_csv(OUT / "ctoi_singles_dropped.csv", index=False)
    print(f"kept {len(kept)} CTOIs on {kept.TIC.nunique()} TICs; dropped {len(dropped)}:")
    print(dropped.drop_reason.str.split(",").str[0].value_counts().to_string())
    print("kept by uploader:", kept.User.value_counts().head(8).to_dict())
    print("Tmag quartiles:", kept.Tmag.quantile([.25, .5, .75]).round(2).tolist(),
          "| depth_ppm median:", round(kept.depth_ppm.median()), "| dur_h median:", round(kept.dur_h.median(), 2))
    print("lastmod year:", pd.to_datetime(kept["CTOI lastmod"]).dt.year.value_counts().sort_index().to_dict())


if __name__ == "__main__":
    main()
