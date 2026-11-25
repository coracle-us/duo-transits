"""
tess_wide step 2 (PLAN.md section 2): download the QLP light curve of every star that passed the cut, keep the
columns the search needs as a compressed .npz, and delete the FITS file.

Kept per star: TIME (float64), SAP_FLUX, SAP_BKG, SAP_X, SAP_Y (float32), QUALITY, ORBITID (int32), and the
header's CAMERA, CCD, TESSMAG. Downloads run 16 at a time from MAST; a rerun skips stars already done.

Output (outside git): tess_wide/data/lc_s<NNN>/<TIC>.npz, tess_wide/data/s<NNN>_download.json (bytes, seconds)
Usage:  python tess_wide/s2_download.py --sector 104 [--limit N]
"""
import argparse
import json
import re
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"


def urls(sector):
    sh = (DATA / f"hlsp_qlp_tess_ffi_s{sector:04d}_tess_v01_llc-fits.sh").read_text()
    return {int(m.group(1)): m.group(0) for m in
            re.finditer(rf"https://mast\.stsci\.edu/api/v0\.1/Download/file/\?uri=mast:HLSP/qlp/s{sector:04d}/"
                        rf"[0-9/]+/hlsp_qlp_tess_ffi_s{sector:04d}-0*(\d+)_tess_v01_llc\.fits", sh)}


def one(tic, url, out, tmpdir):
    if out.exists():
        return 0, 0
    tmp = tmpdir / f"{tic}.fits"
    for attempt in range(5):
        try:
            urllib.request.urlretrieve(url, tmp)
            break
        except Exception:  # noqa: BLE001
            time.sleep(5 * (attempt + 1))
    else:
        return -1, 0
    nbytes = tmp.stat().st_size
    with fits.open(tmp) as h:
        d, H = h[1].data, h[0].header
        np.savez_compressed(out.with_suffix(".tmp.npz"), time=np.asarray(d["TIME"], "f8"),
                            sap=np.asarray(d["SAP_FLUX"], "f4"), bkg=np.asarray(d["SAP_BKG"], "f4"),
                            x=np.asarray(d["SAP_X"], "f4"), y=np.asarray(d["SAP_Y"], "f4"),
                            quality=np.asarray(d["QUALITY"], "i4"), orbit=np.asarray(d["ORBITID"], "i4"),
                            camera=int(H["CAMERA"]), ccd=int(H["CCD"]), tessmag=float(H.get("TESSMAG") or np.nan))
    out.with_suffix(".tmp.npz").rename(out)
    tmp.unlink()
    return nbytes, out.stat().st_size


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", type=int, default=104)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--stars", default="", help="star list to use instead of s<NNN>_stars.csv (e.g. the partial one)")
    a = ap.parse_args()
    s = a.sector
    stars = pd.read_csv(Path(a.stars) if a.stars else DATA / f"s{s:04d}_stars.csv")
    u = urls(s)
    tics = [int(t) for t in stars.id if int(t) in u]
    if a.limit:
        tics = tics[:a.limit]
    od, tmpdir = DATA / f"lc_s{s:04d}", DATA / f"tmp_s{s:04d}"
    od.mkdir(parents=True, exist_ok=True)
    tmpdir.mkdir(parents=True, exist_ok=True)
    t0, got, kept, fail, n, newly = time.time(), 0, 0, 0, 0, 0
    with ThreadPoolExecutor(a.workers) as ex:
        for b_in, b_out in ex.map(lambda t: one(t, u[t], od / f"{t}.npz", tmpdir), tics):
            n += 1
            if b_in < 0:
                fail += 1
            got, kept = got + max(b_in, 0), kept + b_out
            newly += b_in > 0
            if n % 2000 == 0:
                el = time.time() - t0
                print(f"{n}/{len(tics)}  {got / 1e9:.2f} GB in, {kept / 1e9:.2f} GB kept, {el:.0f} s, "
                      f"{n / el:.1f} files/s, failed {fail}", flush=True)
    el = time.time() - t0
    summ = dict(sector=s, stars=len(tics), downloaded_now=newly, downloaded_bytes=got, kept_bytes=kept, failed=fail, seconds=round(el),
                kept_total_bytes=sum(p.stat().st_size for p in od.glob("*.npz")))
    prev = DATA / f"s{s:04d}_download.json"
    if prev.exists() and not a.limit:  # a resumed run: keep the first run's timing too
        summ["previous"] = json.loads(prev.read_text())
    prev.write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ))


if __name__ == "__main__":
    main()
