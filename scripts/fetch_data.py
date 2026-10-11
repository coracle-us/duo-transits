"""
Duotransit project - data fetch script (run on Jon's machine).

Downloads everything the cloud workspace can't reach:
  1. LaTeX source of Hawthorn et al. (arXiv:2310.17268v2) for the full candidate + alias tables
  2. NASA Exoplanet Archive confirmed planets (TIC-matched) and the ExoFOP TOI table
  3. All available TESS light curves for the 85 duotransit candidate stars

Usage (from the repo root):
    python -m pip install -r requirements.txt
    python scripts/fetch_data.py

Re-running is safe: finished targets and already-saved files are skipped; failed or partial
targets are retried. Output goes to data/raw/
"""
import io
import json
import re
import sys
import tarfile
import time
import threading
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "data" / "raw"
LC_DIR = OUT / "lightcurves"
OUT.mkdir(parents=True, exist_ok=True)
LC_DIR.mkdir(exist_ok=True)

TICS = [
    1504460, 2669681, 9844069, 13072758, 13713700, 14445414, 20904104, 25194908,
    32179255, 38138512, 39167176, 42428568, 52195587, 66439839, 67599025, 70561926,
    71028120, 77437543, 81089255, 100776118, 101824521, 107113345, 118339710,
    123763494, 124029677, 140215502, 140750416, 142278054, 145006304, 145913596,
    153838604, 156716001, 157119927, 157698565, 159490807, 161169240, 176518126,
    188620407, 193318850, 207783865, 211409161, 215402824, 218977148, 220622886,
    221915858, 224279805, 235058563, 237605045, 242241304, 251057075, 256912435,
    265465927, 265466589, 268534931, 269333648, 275878706, 279727635, 287204963,
    289840544, 290165539, 292719109, 294097549, 296737508, 300394149, 303317324,
    304339227, 306249066, 313671132, 317923092, 323295479, 332697924, 333736132,
    339399841, 349091983, 381553868, 393229954, 394340183, 396720998, 404518509,
    412386707, 437293313, 439491923, 442893494, 457649900, 466206508,
]
assert len(TICS) == 85

LOG = OUT / "fetch_log.txt"


_log_lock = threading.Lock()


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    with _log_lock:
        print(line, flush=True)
        with open(LOG, "a", encoding="utf-8") as f:
            f.write(line + "\n")


def get(url, **kw):
    r = requests.get(url, timeout=300, headers={"User-Agent": "duotransit-research/0.1"}, **kw)
    r.raise_for_status()
    return r


# ---------------------------------------------------------------- 1. paper source
def fetch_paper():
    dest = OUT / "paper_src"
    if dest.exists() and any(dest.iterdir()):
        log("paper source already present")
        return
    dest.mkdir(exist_ok=True)
    log("downloading arXiv source 2310.17268v2")
    data = get("https://arxiv.org/e-print/2310.17268v2").content
    try:
        with tarfile.open(fileobj=io.BytesIO(data)) as tf:
            for m in tf.getmembers():
                if m.isfile() and m.name.lower().endswith((".tex", ".txt", ".csv", ".dat", ".bbl")):
                    target = dest / Path(m.name).name
                    target.write_bytes(tf.extractfile(m).read())
    except tarfile.ReadError:
        (dest / "source_single_file.tex.gz").write_bytes(data)
    log(f"paper source saved: {[p.name for p in dest.iterdir()]}")


# ---------------------------------------------------------------- 2. catalogs
def fetch_catalogs():
    ps = OUT / "nasa_pscomppars.csv"
    if not ps.exists():
        log("downloading NASA Exoplanet Archive pscomppars")
        # pscomppars has no pl_refname; the period's provenance is pl_orbper_reflink
        q = ("select pl_name,hostname,tic_id,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_tranmid,"
             "pl_rade,pl_bmassj,disc_year,disc_facility,disc_refname,pl_orbper_reflink "
             "from pscomppars")
        url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
        ps.write_bytes(get(url, params={"query": q, "format": "csv"}).content)
    toi = OUT / "exofop_toi.csv"
    if not toi.exists():
        log("downloading ExoFOP TOI table")
        toi.write_bytes(get("https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv").content)
    ctoi = OUT / "exofop_ctoi.csv"
    if not ctoi.exists():
        log("downloading ExoFOP CTOI table")
        try:
            ctoi.write_bytes(get("https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv").content)
        except Exception as e:
            log(f"CTOI download failed (non-fatal): {e}")
    log("catalogs done")


# ---------------------------------------------------------------- 3. light curves
# Products kept: SPOC 2-min/20-s, TESS-SPOC FFI, QLP FFI, TGLC FFI, eleanor-lite. We keep ALL
# available products per sector (small files), so the analysis can choose.
AUTHORS = ("SPOC", "TESS-SPOC", "QLP", "TGLC", "GSFC-ELEANOR-LITE")
# Extra columns kept for vetting (raw flux, background, centroids, pointing), whichever exist
EXTRA_COLS = ("sap_flux", "pdcsap_flux", "kspsap_flux", "sap_bkg", "mom_centr1", "mom_centr2",
              "pos_corr1", "pos_corr2", "sap_x", "sap_y", "orbitid", "psf_flux", "aperture_flux",
              "cal_psf_flux", "cal_aper_flux", "background", "tglc_flags", "corr_flux", "raw_flux",
              "pca_flux", "flux_bkg")
N_WORKERS = 6


def _col(x):
    """Column -> float ndarray, masked entries -> NaN (SPOC fluxes are astropy Masked)."""
    v = x.value if hasattr(x, "value") else x
    mask = getattr(v, "mask", None)
    arr = np.asarray(getattr(v, "unmasked", v), dtype=float)
    if mask is not None:
        arr = np.where(np.asarray(mask), np.nan, arr)
    return arr


def _fetch_one(i, tic):
    import lightkurve as lk

    log(f"[{i}/{len(TICS)}] TIC {tic}")
    tdir = LC_DIR / str(tic)
    tdir.mkdir(exist_ok=True)
    sr = None
    for attempt in range(3):
        try:
            sr = lk.search_lightcurve(f"TIC {tic}")
            break
        except Exception as e:
            log(f"  TIC {tic} search retry {attempt+1}: {e}")
            time.sleep(10)
    if sr is None:
        raise RuntimeError("search failed 3 times")
    if len(sr) == 0:
        log(f"  TIC {tic}: no light curves found")
        return {"tic": tic, "status": "done", "n_products": 0, "n_failed": 0}

    tab = sr.table.to_pandas()
    keep_cols = [c for c in ["mission", "author", "exptime", "target_name", "distance"] if c in tab.columns]
    tab[keep_cols].to_csv(tdir / "search_results.csv", index=False)

    n_saved = n_failed = 0
    for j in range(len(sr)):
        author = str(tab.iloc[j]["author"])
        if author not in AUTHORS:
            continue
        mission = str(tab.iloc[j]["mission"])
        exptime = float(tab.iloc[j]["exptime"])
        sector = "".join(ch for ch in mission if ch.isdigit()) or "NA"
        if sector != "NA":
            fname = tdir / f"s{int(sector):03d}_{author}_{int(round(exptime))}s.csv.gz"
        else:
            fname = tdir / f"sNA_{author}_{j}.csv.gz"
        if fname.exists():
            n_saved += 1
            continue
        lc = None
        for attempt in range(3):
            try:
                lc = sr[j].download(quality_bitmask="none")
                break
            except Exception as e:
                log(f"  TIC {tic} download retry {attempt+1} ({author} s{sector}): {e}")
                # A corrupt cached file (e.g. from an interrupted run) is re-read on every
                # retry unless removed, so delete it to force a fresh download.
                m = re.search(r"Data product (.+?\.fits)", str(e))
                if m and Path(m.group(1)).exists():
                    Path(m.group(1)).unlink()
                    log(f"  TIC {tic} removed corrupt cache file {Path(m.group(1)).name}")
                time.sleep(10)
        if lc is None:
            n_failed += 1
            continue
        df = pd.DataFrame({
            "time": _col(lc.time),
            "flux": _col(lc.flux),
            "flux_err": _col(lc.flux_err),
            "quality": np.nan_to_num(_col(lc.quality)).astype(int),
        })
        for extra in EXTRA_COLS:
            if extra in lc.columns:
                df[extra] = _col(lc[extra])
        tmp = fname.with_suffix(".tmp")
        df.to_csv(tmp, index=False, compression="gzip")
        tmp.replace(fname)  # atomic: a killed run never leaves a half-written file
        n_saved += 1
    log(f"  TIC {tic}: saved {n_saved} light curves ({n_failed} failed)")
    status = "done" if n_failed == 0 else "partial"
    return {"tic": tic, "status": status, "n_products": n_saved, "n_failed": n_failed}


def _fetch_one_safe(i, tic):
    try:
        return _fetch_one(i, tic)
    except Exception:
        log(f"  TIC {tic} FAILED:\n" + traceback.format_exc())
        return {"tic": tic, "status": "failed", "n_products": 0, "n_failed": 0}


def fetch_lightcurves():
    inv_path = OUT / "lc_inventory.csv"
    inventory = {}
    if inv_path.exists():
        inventory = {int(r["tic"]): r for r in pd.read_csv(inv_path).to_dict("records")}
    todo = [(i, t) for i, t in enumerate(TICS, 1) if inventory.get(t, {}).get("status") != "done"]
    log(f"light curves: {len(todo)} targets to fetch, {N_WORKERS} workers")

    # Processes, not threads: lightkurve silences download progress bars by swapping the
    # process-wide sys.stdout for a devnull handle and closing it, which breaks other threads.
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = {ex.submit(_fetch_one_safe, i, tic): tic for i, tic in todo}
        for fut in as_completed(futs):
            rec = fut.result()
            inventory[rec["tic"]] = rec
            pd.DataFrame(list(inventory.values())).sort_values("tic").to_csv(inv_path, index=False)
    bad = [t for t, r in inventory.items() if r["status"] != "done"]
    log(f"light curves finished; not fully done: {bad if bad else 'none'} (re-run to retry)")


def main():
    log(f"python {sys.version.split()[0]}")
    for step in (fetch_paper, fetch_catalogs, fetch_lightcurves):
        try:
            step()
        except Exception:
            log(f"step {step.__name__} failed:\n" + traceback.format_exc())
    log("ALL DONE - tell Claude it finished")


if __name__ == "__main__":
    main()
