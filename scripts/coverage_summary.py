"""
Plan step 5 (coverage): which TESS sectors each candidate has beyond the data used by
Hawthorn et al. (2024).

"Paper data" = Sectors 1-39 (Cycles 1 and 3) plus the sectors of each candidate's own Tc1/Tc2
(one candidate, TIC 188620407, has its second transit in Sector 42). Everything else is "new".

Inputs : data/raw/lightcurves/<TIC>/search_results.csv  (all MAST products, incl. not downloaded)
         data/raw/lightcurves/<TIC>/s<NNN>_<author>_<exptime>s.csv.gz  (downloaded light curves)
Output : data/catalog/coverage.csv  one row per candidate:
           paper_sectors, new_sectors, n_new_sectors, new_days (days of good-quality data in new
           sectors, best product per sector), products available per new sector
         data/catalog/coverage_by_sector.csv  one row per (TIC, sector, product) downloaded

Usage:  python scripts/coverage_summary.py
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
LC = ROOT / "data" / "raw" / "lightcurves"
CAT = ROOT / "data" / "catalog"

# Best product first; used to pick one light curve per sector for the good-data estimate
PRODUCT_RANK = ["SPOC_20", "SPOC_120", "TESS-SPOC", "QLP", "TGLC", "GSFC-ELEANOR-LITE"]
VALID_SECTOR_MAX = 200  # MAST has at least one mislabelled product ("Sector 1751")


def product_key(author, exptime):
    return f"SPOC_{int(exptime)}" if author == "SPOC" else author


def good_days(path):
    """Days of usable data: finite flux, quality == 0, summed over cadences."""
    df = pd.read_csv(path, usecols=["time", "flux", "quality"])
    ok = np.isfinite(df.time) & np.isfinite(df.flux) & (df.quality == 0)
    t = df.time[ok].to_numpy()
    if len(t) < 2:
        return 0.0, len(t)
    cad = np.median(np.diff(t))
    return float(ok.sum() * cad), int(ok.sum())


def main():
    cand = pd.read_csv(CAT / "candidates.csv")
    rows, per_sector = [], []
    for _, c in cand.iterrows():
        tdir = LC / str(c.TIC)
        paper = set(range(1, 40)) | {int(c.S1), int(c.S2)}
        sr_path = tdir / "search_results.csv"
        if not sr_path.exists():
            rows.append(dict(TIC=c.TIC, note="no search results (fetch not run/failed)"))
            continue
        sr = pd.read_csv(sr_path)
        sr["sector"] = sr.mission.str.extract(r"(\d+)").astype(float)
        bad = sr[sr.sector > VALID_SECTOR_MAX]
        sr = sr[sr.sector <= VALID_SECTOR_MAX]
        all_sec = sorted(int(s) for s in sr.sector.dropna().unique())
        new_sec = [s for s in all_sec if s not in paper]

        # downloaded files -> good data per sector, best product
        files = {}
        for f in tdir.glob("s*.csv.gz"):
            m = re.fullmatch(r"s(\d+)_(.+)_(\d+)s\.csv\.gz", f.name)
            if not m:
                continue
            sec, author, exp = int(m.group(1)), m.group(2), int(m.group(3))
            if sec > VALID_SECTOR_MAX:
                continue
            files.setdefault(sec, {})[product_key(author, exp)] = f
        new_days, new_products, missing_dl = 0.0, [], []
        for s in new_sec:
            avail = files.get(s, {})
            if not avail:
                missing_dl.append(s)
                continue
            best = next(k for k in PRODUCT_RANK if k in avail)
            d, n = good_days(avail[best])
            new_days += d
            new_products.append(f"{s}:{'+'.join(k for k in PRODUCT_RANK if k in avail)}")
            for k in avail:
                dd, nn = good_days(avail[k]) if k == best else (np.nan, np.nan)
                per_sector.append(dict(TIC=c.TIC, sector=s, product=k, best=(k == best),
                                       good_days=dd, n_good=nn))
        spoc2min_new = sorted(s for s in new_sec if "SPOC_120" in files.get(s, {}) or "SPOC_20" in files.get(s, {}))
        rows.append(dict(
            TIC=c.TIC, Tc1_sector=int(c.S1), Tc2_sector=int(c.S2),
            paper_sectors=" ".join(map(str, [s for s in all_sec if s in paper])),
            new_sectors=" ".join(map(str, new_sec)), n_new_sectors=len(new_sec),
            new_good_days=round(new_days, 1),
            new_spoc_2min_sectors=" ".join(map(str, spoc2min_new)),
            new_products=" ".join(new_products),
            new_sectors_not_downloaded=" ".join(map(str, missing_dl)),
            note=f"ignored mislabelled MAST product(s): {list(bad.mission)}" if len(bad) else "",
        ))
    out = pd.DataFrame(rows)
    out.to_csv(CAT / "coverage.csv", index=False)
    pd.DataFrame(per_sector).to_csv(CAT / "coverage_by_sector.csv", index=False)
    print(f"wrote {CAT / 'coverage.csv'}")
    print(out[["TIC", "n_new_sectors", "new_good_days", "new_sectors"]]
          .sort_values("new_good_days", ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
