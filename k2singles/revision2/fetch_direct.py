"""
K2 revision (editor, 2026-11-15): fetch every light curve the blind search needs, for ALL 163 catalogue stars
(no selection on catalogued depth, duration or expected SNR), using direct archive file URLs.

Why direct URLs: on 2026-11-15 the MAST portal product-list service (used by lightkurve.search_lightcurve) failed
with "Request failed" after ~3 min, while archive file URLs and the TESScut API worked.

Per star (written in the same format as scripts/fetch_data.py, so alias_pipeline.load_lightcurves reads them):
  sectors  TESScut sector list (API) union tess-point 0.9.5; sectors > 200 are skipped, as in the note's pipeline
  SPOC     2-min light curves, filename from the per-sector MAST bulk script (cached in data/bulk/)
  TESS-SPOC, QLP   HLSP file URLs (a 404 means no product)
  OWN      our TESScut photometry for sectors >= 80 without a pipeline product (as pilot.add_own_photometry in
           the submitted note; cutouts from the TESScut API, see own_photometry)
  K2       mission PDCSAP for every campaign (k2lib.k2_lightcurves, direct URLs), and EVEREST v2 and K2SFF HLSPs,
           saved as npz in k2singles/revision2/data/hlsp/<EPIC>/
Resumable; one process per star. Status: k2singles/revision2/data/fetch_status.csv.
Usage:  python k2singles/revision2/fetch_direct.py [n_workers] [EPIC ...]
"""
import re
import sys
import time
import traceback
import urllib.error
import urllib.request
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = HERE / "data"
HLSP = DATA / "hlsp"
BULK = DATA / "bulk"
TMP = DATA / "tmp"
for p in (ROOT / "scripts", ROOT / "singles" / "scripts", ROOT / "k2singles" / "scripts"):
    sys.path.insert(0, str(p))
ARCH = "https://archive.stsci.edu"


def stars():
    d = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv")
    return d.drop_duplicates("EPIC")[["EPIC", "TIC", "ra", "dec"]].reset_index(drop=True)


def get(url, dest, tries=4):
    """Download url to dest. True on success, False on 404."""
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                data = r.read()
            tmp = dest.with_suffix(dest.suffix + ".part")
            tmp.write_bytes(data)
            tmp.replace(dest)
            return True
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return False
            time.sleep(10 * (k + 1))
        except Exception:  # noqa: BLE001
            time.sleep(10 * (k + 1))
    raise RuntimeError(f"download failed: {url}")


def tid_path(tic):
    s = f"{int(tic):016d}"
    return s, f"{s[0:4]}/{s[4:8]}/{s[8:12]}/{s[12:16]}"


def spoc_filename(sector, tic):
    """2-min LC filename for this TIC from the sector's bulk curl script, or None."""
    BULK.mkdir(parents=True, exist_ok=True)
    f = BULK / f"tesscurl_sector_{sector}_lc.sh"
    if not f.exists():
        if not get(f"{ARCH}/missions/tess/download_scripts/sector/tesscurl_sector_{sector}_lc.sh", f):
            f.write_text("")
    pat = f"-s{sector:04d}-{int(tic):016d}-"
    for line in f.read_text().splitlines():
        if pat in line and line.startswith("curl"):
            return line.split(" -o ")[1].split()[0]
    return None


def sectors(ra, dec, tesscut_only=False):
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astroquery.mast import Tesscut
    from tess_stars2px import tess_stars2px_function_entry as t2p
    secs = set()
    for k in range(3):
        try:
            tab = Tesscut.get_sectors(coordinates=SkyCoord(ra * u.deg, dec * u.deg))
            secs |= {int(s) for s in tab["sector"]}
            break
        except Exception:  # noqa: BLE001
            time.sleep(15)
    if not tesscut_only:
        r = t2p(0, ra, dec)
        secs |= {int(s) for s in r[3] if int(s) > 0}
    return sorted(s for s in secs if s <= 200)


def save_csv(fits_path, out):
    import fetch_data as fd
    import lightkurve as lk
    lc = lk.read(str(fits_path), quality_bitmask="none")
    df = pd.DataFrame({"time": fd._col(lc.time), "flux": fd._col(lc.flux), "flux_err": fd._col(lc.flux_err),
                       "quality": np.nan_to_num(fd._col(lc.quality)).astype(int)})
    for extra in fd.EXTRA_COLS:
        if extra in lc.columns:
            df[extra] = fd._col(lc[extra])
    tmp = out.with_suffix(".tmp")
    df.to_csv(tmp, index=False, compression="gzip")
    tmp.replace(out)


def tess(tic, ra, dec):
    import fetch_data as fd
    tdir = fd.LC_DIR / str(tic)
    tdir.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)
    s16, path = tid_path(tic)
    rows = []
    for sec in sectors(ra, dec):
        cands = []
        fn = spoc_filename(sec, tic)
        if fn:
            cands.append(("SPOC", 120, f"{ARCH}/missions/tess/tid/s{sec:04d}/{path}/{fn}"))
        cands.append(("TESS-SPOC", None, f"{ARCH}/hlsps/tess-spoc/s{sec:04d}/target/{path}/"
                                         f"hlsp_tess-spoc_tess_phot_{s16}-s{sec:04d}_tess_v1_lc.fits"))
        cands.append(("QLP", None, f"{ARCH}/hlsps/qlp/s{sec:04d}/{path}/hlsp_qlp_tess_ffi_s{sec:04d}-{s16}_tess_v01_llc.fits"))
        for author, exp, url in cands:
            if exp is None:  # FFI cadence by sector: 1800 s (S1-26), 600 s (S27-55), 200 s (S56+)
                exp = 1800 if sec <= 26 else 600 if sec <= 55 else 200
            out = tdir / f"s{sec:03d}_{author}_{exp}s.csv.gz"
            if out.exists():
                rows.append(dict(mission=f"TESS Sector {sec:02d}", author=author, exptime=exp))
                continue
            tmpf = TMP / f"{tic}_{sec}_{author}.fits"
            if get(url, tmpf):
                save_csv(tmpf, out)
                tmpf.unlink()
                rows.append(dict(mission=f"TESS Sector {sec:02d}", author=author, exptime=exp))
    pd.DataFrame(rows, columns=["mission", "author", "exptime"]).to_csv(tdir / "search_results.csv", index=False)
    return " ".join(r["author"] + ":" + r["mission"].split()[-1] for r in rows)


def own_photometry(tic, ra, dec):
    """Our TESScut photometry (singles/scripts/tesscut_photometry.py) for sectors >= 80 without a pipeline product,
    as pilot.add_own_photometry, but the cutout comes from the TESScut API directly: pilot uses
    lightkurve.search_tesscut, which goes through the failing portal service and returns nothing."""
    import astropy.units as u
    import diff_image as di
    import fetch_data as fd
    import tesscut_photometry as tp
    from astropy.coordinates import SkyCoord
    from astroquery.mast import Tesscut
    have = set()
    for f in (fd.LC_DIR / str(tic)).glob("s*.csv.gz"):
        m = re.fullmatch(r"s(\d+)_(.+)_(\d+)s\.csv\.gz", f.name)
        if m and m.group(2) in ("SPOC", "TESS-SPOC", "QLP", "OWN"):
            have.add(int(m.group(1)))
    added = []
    # only sectors TESScut lists: tess-point also predicts sectors TESScut does not (yet) hold, and each of those
    # cost two minutes of retries before failing (seen on 2026-11-15)
    for sec in sectors(ra, dec, tesscut_only=True):
        if not (80 <= sec <= 200) or sec in have:
            continue
        d = di.RAW / str(tic)
        d.mkdir(parents=True, exist_ok=True)
        if not list(d.glob(f"*s{sec:04d}*.fits")):
            for k in range(3):
                try:
                    Tesscut.download_cutouts(coordinates=SkyCoord(ra * u.deg, dec * u.deg), size=di.SIZE, sector=sec,
                                             path=str(d))
                    break
                except Exception:  # noqa: BLE001
                    time.sleep(20 * (k + 1))
        if list(d.glob(f"*s{sec:04d}*.fits")):
            try:
                tp.photometry(tic, sec, ra, dec)
                added.append(sec)
            except BaseException as e:  # noqa: BLE001  (tesscut_photometry raises SystemExit)
                print(f"  own photometry TIC {tic} S{sec} failed: {e}")
    return added


def hlsp(epic, campaigns):
    from astropy.io import fits
    out = HLSP / str(epic)
    out.mkdir(parents=True, exist_ok=True)
    s = str(epic)
    got = []
    for c in campaigns:
        for author, tmpl in (("EVEREST", "everest/v2/c{cc}/{a}00000/{b}/hlsp_everest_k2_llc_{e}-c{cc}_kepler_v2.0_lc.fits"),
                             ("K2SFF", "k2sff/c{cc}/{a}00000/{b}/hlsp_k2sff_k2_lightcurve_{e}-c{cc}_kepler_v1_llc.fits")):
            fn = out / f"{author}_c{c:03d}.npz"
            if fn.exists():
                got.append(fn.stem)
                continue
            tmpf = TMP / f"{epic}_{c}_{author}.fits"
            ok = get(f"{ARCH}/hlsps/" + tmpl.format(cc=f"{c:02d}", a=s[:4], b=s[4:9], e=s), tmpf)
            if not ok:
                continue
            with fits.open(tmpf) as h:
                d = h[1].data
                t = np.asarray(d["TIME"] if "TIME" in d.names else d["T"], float) + 2454833 - 2457000
                if author == "EVEREST":
                    f = np.asarray(d["FLUX"], float)
                    q = np.asarray(d["QUALITY"], int)
                    ok_ = np.isfinite(t) & np.isfinite(f) & (q == 0)
                else:
                    f = np.asarray(d["FCOR"], float)
                    ok_ = np.isfinite(t) & np.isfinite(f)
            tmpf.unlink()
            if ok_.sum() > 50:
                np.savez(fn, t=t[ok_], f=f[ok_] / np.nanmedian(f[ok_]))
                got.append(fn.stem)
    return " ".join(got)


def work(row):
    warnings.filterwarnings("ignore")
    import k2lib
    epic, tic = int(row["EPIC"]), int(row["TIC"])
    rec = dict(EPIC=epic, TIC=tic)
    try:
        rec["tess"] = tess(tic, row["ra"], row["dec"])
        rec["own_sectors"] = " ".join(map(str, own_photometry(tic, row["ra"], row["dec"])))
        k2lib.k2_lightcurves(epic)
        camps = sorted(int(f.stem.split("-c")[1].split("_")[0]) for f in (k2lib.CACHE / str(epic) / "dl").glob("ktwo*_llc.fits"))
        rec["k2pdc_campaigns"] = " ".join(map(str, camps))
        rec["hlsp"] = hlsp(epic, camps)
        rec["status"] = "done"
    except Exception:  # noqa: BLE001
        rec["status"] = "failed"
        rec["error"] = traceback.format_exc()[-1500:]
    return rec


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--own-pass":
        # second pass for stars fetched before own_photometry existed
        warnings.filterwarnings("ignore")
        st = stars()
        for r in st.to_dict("records"):
            added = own_photometry(int(r["TIC"]), r["ra"], r["dec"])
            if added:
                print(f"EPIC {r['EPIC']} TIC {r['TIC']}: OWN sectors {added}", flush=True)
        return
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    only = {int(x) for x in sys.argv[2:]}
    DATA.mkdir(parents=True, exist_ok=True)
    st = stars()
    status_f = DATA / "fetch_status.csv"
    done = set()
    if status_f.exists():
        s = pd.read_csv(status_f)
        done = set(s[s.status == "done"].EPIC)
    todo = [r for r in st.to_dict("records") if r["EPIC"] not in done and (not only or r["EPIC"] in only)]
    print(f"{len(st)} stars, {len(todo)} to fetch, {n} workers", flush=True)
    with ProcessPoolExecutor(n, max_tasks_per_child=1) as ex:
        futs = [ex.submit(work, r) for r in todo]
        for i, fu in enumerate(as_completed(futs), 1):
            r = fu.result()
            print(f"[{i}/{len(todo)}] EPIC {r['EPIC']} {r['status']} | {r.get('tess','')} | own={r.get('own_sectors','')} "
                  f"| k2={r.get('k2pdc_campaigns','')} | {r.get('hlsp','')}", flush=True)
            if r["status"] != "done":
                print(r.get("error", ""), flush=True)
            old = pd.read_csv(status_f) if status_f.exists() else pd.DataFrame(columns=["EPIC"])
            pd.concat([old[old.EPIC != r["EPIC"]], pd.DataFrame([r])]).to_csv(status_f, index=False)
    print("done", flush=True)


if __name__ == "__main__":
    main()
