"""
Step 1: open single / duo / no-period / maximum-period TOIs from the local ExoFOP TOI table.

Selection (reproduces the 156-row list):
  - TFOPWG disposition not in {CP, KP, FP, FA} (PC, APC or blank)
  - AND any of
      (a) no period (Period blank or 0)
      (b) comment mentions single / duo / mono
      (c) comment mentions a maximum period ("max period", "maximum ...")
Refinements (flags, nothing dropped):
  - max_period_d: parsed from the comment ("max period of ~702.47 d"), else the catalogue period
    when the comment calls it a maximum, else NaN
  - comment_says_resolved: comment states the period is already pinned down (e.g. "confirms period is 1/46")
  - still_single_sector: N from "single ... as of sector N" / "in sN" style comments (last sector named)
  - n_named_transits: "two transits" / "three transits" / "3 transits" parsed where present
  - ecliptic latitude (for the K2 footprint filter in step 2)

Input:  data/raw/exofop_toi.csv (2026-10-11; never edited)
Output: toi_singles_k2/targets.csv
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import BarycentricTrueEcliptic, SkyCoord

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "toi_singles_k2"

WORDNUM = {"one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "single": 1}


def parse_max_period(c):
    m = re.search(r"max(?:imum)?\.?\s*(?:period)?\s*(?:between events[^0-9~]*)?(?:of|is|=)?\s*~?\s*([\d.]+)", c, re.I)
    if m:
        try:
            return float(m.group(1).rstrip("."))
        except ValueError:
            pass
    m = re.search(r"~?\s*([\d.]+)\s*d(?:ays?)?\s+is\s+(?:the\s+)?maximum period", c, re.I)
    if m:
        return float(m.group(1))
    return np.nan


def parse_n_transits(c):
    m = re.search(r"\b(one|two|three|four|five|single|\d)\s+transits?\b", c, re.I)
    if not m:
        return np.nan
    w = m.group(1).lower()
    return WORDNUM.get(w, int(w) if w.isdigit() else np.nan)


def parse_last_sector(c):
    secs = [int(x) for x in re.findall(r"(?:sector|\bs)\s*-?\s*(\d{1,3})\b", c, re.I)]
    return max(secs) if secs else np.nan


def main():
    t = pd.read_csv(ROOT / "data" / "raw" / "exofop_toi.csv")
    o = t[~t["TFOPWG Disposition"].isin(["CP", "KP", "FP", "FA"])].copy()
    c = o.Comments.fillna("")
    P = o["Period (days)"]
    nop = P.isna() | (P == 0)
    sing = c.str.contains(r"single|duo|mono", case=False)
    mx = c.str.contains(r"max(?:imum)? period|period (?:of|is) ?(?:at most|<)|maximum", case=False)
    sel = nop | sing | mx
    d = o[sel].copy()
    print(f"open TOIs: {len(o)}; selected: {len(d)} (no-period {nop.sum()}, single/duo/mono {sing.sum()}, "
          f"max-period comment {mx.sum()})")

    d["sel_noP"] = nop[sel].values
    d["sel_single_duo"] = sing[sel].values
    d["sel_maxP"] = mx[sel].values
    cc = d.Comments.fillna("")
    d["epoch_btjd"] = d["Epoch (BJD)"] - 2457000.0
    d["max_period_d"] = [parse_max_period(x) for x in cc]
    # catalogue period equals the max period when the comment says so and no number was parsed
    fill = d.max_period_d.isna() & d.sel_maxP & (d["Period (days)"] > 0)
    d.loc[fill, "max_period_d"] = d.loc[fill, "Period (days)"]
    d["comment_says_resolved"] = cc.str.contains(r"confirms period|period confirmed|true period is \d", case=False)
    d["n_named_transits"] = [parse_n_transits(x) for x in cc]
    d["last_sector_in_comment"] = [parse_last_sector(x) for x in cc]

    sc = SkyCoord(d.RA.values, d.Dec.values, unit=(u.hourangle, u.deg))
    d["ra_deg"] = sc.ra.deg
    d["dec_deg"] = sc.dec.deg
    d["ecl_lat"] = sc.transform_to(BarycentricTrueEcliptic()).lat.deg
    d["ecl_lon"] = sc.transform_to(BarycentricTrueEcliptic()).lon.deg

    cols = ["TOI", "TIC ID", "TFOPWG Disposition", "TESS Disposition", "TESS Mag", "ra_deg", "dec_deg", "ecl_lat",
            "ecl_lon", "epoch_btjd", "Epoch (BJD) err", "Period (days)", "max_period_d", "Duration (hours)",
            "Depth (ppm)", "Planet Radius (R_Earth)", "Stellar Radius (R_Sun)", "Stellar Eff Temp (K)",
            "Sectors", "Source", "Date TOI Alerted (UTC)", "Date TOI Updated (UTC)", "sel_noP", "sel_single_duo",
            "sel_maxP", "n_named_transits", "last_sector_in_comment", "comment_says_resolved", "Comments"]
    d = d[cols].rename(columns={"TIC ID": "TIC", "TFOPWG Disposition": "TFOPWG", "TESS Mag": "Tmag",
                                "Period (days)": "period_d", "Duration (hours)": "dur_h", "Depth (ppm)": "depth_ppm",
                                "Planet Radius (R_Earth)": "Rp_Re", "Stellar Radius (R_Sun)": "Rstar",
                                "Stellar Eff Temp (K)": "Teff", "Date TOI Alerted (UTC)": "alerted",
                                "Date TOI Updated (UTC)": "updated", "Epoch (BJD) err": "epoch_err"})
    d = d.sort_values("TOI").reset_index(drop=True)
    d.to_csv(OUT / "targets.csv", index=False)
    print(f"stars: {d.TIC.nunique()}; with max period parsed/assigned: {d.max_period_d.notna().sum()}; "
          f"comment says resolved: {d.comment_says_resolved.sum()}")
    print(f"|ecliptic lat| < 20 deg: {(d.ecl_lat.abs() < 20).sum()}")
    print("written", OUT / "targets.csv")


if __name__ == "__main__":
    main()
