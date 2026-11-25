"""
tess_wide step 1 (PLAN.md section 3): TIC 8.2 rows for every star in a sector's QLP file list, then the star cut.

The QLP bulk script (MAST) gives the TIC ids. TIC rows come from the MAST TIC TAP service in chunks of 300 ids
(longer IN lists hit a 414 URI limit; cone queries are capped at 0.5 deg there), 16 requests at a time. Chunks
are cached, so a rerun resumes.

Output (under tess_wide/data/, outside git):
  s<NNN>_tic.csv     every star in the QLP list with its TIC columns
  s<NNN>_stars.csv   the stars that pass the cut, with the rule each failed star missed in s<NNN>_tic.csv
  s<NNN>_starcut.json  counts per rule
Usage:  python tess_wide/s1_starcut.py --sector 104
"""
import argparse
import io
import json
import re
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io.votable import parse_single_table

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TAP = "https://mast.stsci.edu/vo-tap/api/v0.1/tic/sync"
COLS = "id,ra,dec,tmag,rad,teff,logg,contratio,objtype,disposition"
CHUNK = 300

# PLAN.md section 3 (v1.0)
TMAG_MAX, RAD_MAX, TEFF_MIN, TEFF_MAX, LOGG_MIN, CONTRATIO_MAX = 12.0, 2.0, 3000.0, 7500.0, 3.5, 1.0


def tic_ids(sector):
    sh = DATA / f"hlsp_qlp_tess_ffi_s{sector:04d}_tess_v01_llc-fits.sh"
    if not sh.exists():
        url = f"https://archive.stsci.edu/hlsps/qlp/download_scripts/{sh.name}"
        urllib.request.urlretrieve(url, sh)
    ids = sorted({int(m) for m in re.findall(rf"s{sector:04d}-(\d{{16}})_tess", sh.read_text())})
    return ids, sh


def fetch(chunk_ids, out):
    if out.exists() and out.stat().st_size > 0:
        return out
    q = f"SELECT {COLS} FROM dbo.catalogrecord WHERE id IN ({','.join(map(str, chunk_ids))})"
    body = urllib.parse.urlencode({"LANG": "ADQL", "QUERY": q}).encode()
    for attempt in range(5):
        try:
            with urllib.request.urlopen(urllib.request.Request(TAP, data=body), timeout=300) as r:
                x = r.read()
            if b'QUERY_STATUS" value="ERROR' in x:
                raise RuntimeError(x[:400])
            tmp = out.with_suffix(".tmp")
            tmp.write_bytes(x)
            tmp.rename(out)
            return out
        except Exception as e:  # noqa: BLE001
            print("retry", out.name, attempt, str(e)[:120], flush=True)
            time.sleep(10 * (attempt + 1))
    raise RuntimeError(f"chunk failed: {out.name}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", type=int, default=104)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--partial", action="store_true",
                    help="no fetching: apply the cut to the chunks already cached (to start downloads early)")
    a = ap.parse_args()
    s = a.sector
    ids, sh = tic_ids(s)
    cdir = DATA / f"tic_chunks_s{s:04d}"
    cdir.mkdir(parents=True, exist_ok=True)
    jobs = [(ids[i:i + CHUNK], cdir / f"{i // CHUNK:05d}.vot") for i in range(0, len(ids), CHUNK)]
    t0 = time.time()
    if a.partial:
        jobs = [j for j in jobs if j[1].exists() and j[1].stat().st_size > 0]
    else:
      with ThreadPoolExecutor(a.workers) as ex:
        for k, _ in enumerate(ex.map(lambda j: fetch(*j), jobs), 1):
            if k % 100 == 0:
                print(f"{k}/{len(jobs)} chunks, {time.time() - t0:.0f} s", flush=True)
    tag = "_partial" if a.partial else ""
    parts = [parse_single_table(io.BytesIO(p.read_bytes())).to_table().to_pandas() for _, p in jobs]
    tic = pd.concat(parts, ignore_index=True).drop_duplicates("id")
    tic["objtype"] = tic.objtype.astype(str)
    tic["disposition"] = tic.disposition.fillna("").astype(str).replace({"--": "", "None": "", "nan": ""})
    # first failed rule, in PLAN.md order
    rules = [
        ("objtype", tic.objtype != "STAR"),
        ("disposition", tic.disposition.str.len() > 0),
        ("tmag", ~(tic.tmag <= TMAG_MAX)),
        ("rad", ~(tic.rad <= RAD_MAX)),
        ("teff", ~((tic.teff >= TEFF_MIN) & (tic.teff <= TEFF_MAX))),
        ("logg", tic.logg.notna() & (tic.logg < LOGG_MIN)),
        ("contratio", tic.contratio.notna() & (tic.contratio > CONTRATIO_MAX)),
    ]
    tic["cut"] = ""
    for name, bad in rules:
        tic.loc[(tic.cut == "") & bad.fillna(True).to_numpy(), "cut"] = name
    tic.to_csv(DATA / f"s{s:04d}_tic{tag}.csv", index=False)
    stars = tic[tic.cut == ""].drop(columns="cut")
    stars.to_csv(DATA / f"s{s:04d}_stars{tag}.csv", index=False)
    summ = dict(sector=s, qlp_files=len(ids), tic_rows=len(tic), missing_from_tic=len(ids) - len(tic),
                failed_first=tic.cut.value_counts().to_dict(), passed=len(stars),
                tmag_le_12_any=int((tic.tmag <= TMAG_MAX).sum()), tap_seconds=round(time.time() - t0),
                qlp_file_bytes_each=682560, est_bytes_cut_sector=int(len(stars) * 682560))
    (DATA / f"s{s:04d}_starcut{tag}.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
