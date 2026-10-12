"""
Light-curve fetcher for the EB timing side quest.

Why not scripts/fetch_data._fetch_one directly: it downloads every product of every sector
(including 20-s data, ~10x larger) into data/raw/. For EBs in the CVZ that is 40+ sectors x
4-5 products per star. Here we keep ONE product per sector, in the same file format and
naming as fetch_data (sNNN_<AUTHOR>_<exp>s.csv.gz, columns time, flux, flux_err, quality +
sap_flux when present), reusing fetch_data._col for the column conversion, and write to
sidequests/eb_timing/cache/lc/<TIC>/ (gitignored). FITS files are deleted after conversion.

Product preference per sector: SPOC 120 s > TESS-SPOC (FFI) > QLP (FFI).
Usage:  python sidequests/eb_timing/scripts/fetch_lc.py TIC [TIC ...]
"""
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import fetch_data as fd  # noqa: E402

LCDIR = HERE / "cache" / "lc"
FITS = HERE / "cache" / "fits"
RANK = [("SPOC", 120), ("TESS-SPOC", None), ("QLP", None)]
MAX_SECTOR = 107


def fetch(tic):
    import lightkurve as lk
    tdir = LCDIR / str(tic)
    tdir.mkdir(parents=True, exist_ok=True)
    sr = None
    for attempt in range(3):
        try:
            sr = lk.search_lightcurve(f"TIC {tic}")
            break
        except Exception as e:  # noqa: BLE001
            print(f"TIC {tic} search retry {attempt + 1}: {e}", flush=True)
            time.sleep(10)
    if sr is None or len(sr) == 0:
        return tic, 0, 0
    tab = sr.table.to_pandas()
    tab["sector"] = tab.mission.str.extract(r"(\d+)").astype(float)
    tab["exp"] = tab.exptime.astype(float).round().astype(int)
    tab = tab[tab.sector.notna() & (tab.sector <= MAX_SECTOR)]
    # TESS-SPOC/QLP can only appear in a sector without SPOC 120 s
    tab[["sector", "author", "exp"]].to_csv(tdir / "search_results.csv", index=False)
    n_ok = n_fail = 0
    for s, g in tab.groupby("sector"):
        pick = None
        for author, exp in RANK:
            m = g[(g.author == author) & ((g.exp == exp) if exp else True)]
            if len(m):
                pick = m.index[0]
                break
        if pick is None:
            continue
        row = tab.loc[pick]
        fname = tdir / f"s{int(s):03d}_{row.author}_{row.exp}s.csv.gz"
        if fname.exists():
            n_ok += 1
            continue
        ddir = FITS / str(tic)
        lc = None
        for attempt in range(3):
            try:
                lc = sr[int(pick)].download(quality_bitmask="none", download_dir=str(ddir))
                break
            except Exception as e:  # noqa: BLE001
                print(f"TIC {tic} s{int(s)} retry {attempt + 1}: {str(e)[:100]}", flush=True)
                shutil.rmtree(ddir, ignore_errors=True)
                time.sleep(5)
        if lc is None:
            n_fail += 1
            continue
        df = pd.DataFrame({"time": fd._col(lc.time), "flux": fd._col(lc.flux),
                           "flux_err": fd._col(lc.flux_err),
                           "quality": np.nan_to_num(fd._col(lc.quality)).astype(int)})
        if "sap_flux" in lc.columns:
            df["sap_flux"] = fd._col(lc["sap_flux"])
        tmp = fname.with_suffix(".tmp")
        df.to_csv(tmp, index=False, compression="gzip")
        tmp.replace(fname)
        shutil.rmtree(ddir, ignore_errors=True)
        n_ok += 1
    print(f"TIC {tic}: {n_ok} sectors saved, {n_fail} failed", flush=True)
    return tic, n_ok, n_fail


def main():
    tics = [int(x) for x in sys.argv[1:]]
    with ProcessPoolExecutor(max_workers=4) as ex:
        for tic, n_ok, n_fail in ex.map(fetch, tics):
            pass


if __name__ == "__main__":
    main()
