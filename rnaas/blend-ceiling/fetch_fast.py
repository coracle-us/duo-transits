"""
Bulk download of the light curves make_windows.py needs (same files, same cache, same lcwin.prep), using a few large
MAST observation and product queries and direct parallel downloads instead of one lightkurve search per file.

Input:  data/windows.csv, data/sample.csv, data/eb_windows.csv, data/eb_sample.csv
Cache:  ~/Documents/GitHub/coracle/duo-transits-cache/blend-ceiling/{mast,binned}
Usage:  .venv/bin/python -I fetch_fast.py [--threads 12]
"""
import argparse
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # -I drops the script dir
import numpy as np
import pandas as pd
import requests

import lcwin
from make_windows import CACHE, D, npz_path

URL = "https://mast.stsci.edu/api/v0.1/Download/file?uri="
LOG = open(os.path.join(CACHE, "fetch_fast.log"), "a")


def say(m):
    LOG.write(m + "\n")
    LOG.flush()


def needed():
    s = pd.read_csv(os.path.join(D, "sample.csv"))
    w = pd.read_csv(os.path.join(D, "windows.csv"))
    w = w[w.toi.isin(s[s.drop_reason.isna()].toi)]
    eb = pd.read_csv(os.path.join(D, "eb_sample.csv"))
    ew = pd.read_csv(os.path.join(D, "eb_windows.csv"))
    ew = ew[ew.tic.isin(eb[eb.drop_reason.isna()].tic)]
    j = pd.concat([w, ew])[["tic", "sector", "author"]].drop_duplicates()
    j = j[[not os.path.exists(npz_path(int(r.tic), int(r.sector), r.author)) for r in j.itertuples()]]
    return j


def uris(j):
    from astroquery.mast import Observations as O
    out = []
    for author in ("SPOC", "QLP"):
        jj = j[j.author == author]
        tics = sorted(jj.tic.unique())
        want = set(zip(jj.tic, jj.sector))
        for i in range(0, len(tics), 150):
            names = [str(t) for t in tics[i:i + 150]]
            for attempt in range(4):
                try:
                    if author == "SPOC":
                        ob = O.query_criteria(obs_collection="TESS", dataproduct_type="timeseries",
                                              target_name=names, t_exptime=120)
                    else:
                        ob = O.query_criteria(obs_collection="HLSP", provenance_name="QLP", target_name=names)
                    break
                except Exception as e:
                    say(f"obs query failed: {str(e)[:100]}")
                    time.sleep(15)
            ob = ob.to_pandas()
            if author == "SPOC":
                ob = ob[ob.obs_id.str.endswith("-s")]
            ob["tic"] = ob.target_name.astype(int)
            ob["sector"] = ob.sequence_number.astype(int)
            ob = ob[[(t, s) in want for t, s in zip(ob.tic, ob.sector)]]
            for k in range(0, len(ob), 300):
                part = ob.iloc[k:k + 300]
                for attempt in range(4):
                    try:
                        p = O.get_product_list([str(x) for x in part.obsid]).to_pandas()
                        break
                    except Exception as e:
                        say(f"product query failed: {str(e)[:100]}")
                        time.sleep(15)
                if author == "SPOC":
                    p = p[p.productSubGroupDescription == "LC"]
                else:
                    p = p[p.productFilename.str.endswith("_llc.fits")]
                p = p.merge(part[["obsid", "tic", "sector"]].astype({"obsid": str}),
                            left_on=p.obsID.astype(str), right_on="obsid", how="inner")
                for r in p.itertuples():
                    out.append((int(r.tic), int(r.sector), author, r.dataURI, r.productFilename))
            say(f"{author}: {i + len(names)}/{len(tics)} TICs -> {len(out)} files so far")
    return pd.DataFrame(out, columns=["tic", "sector", "author", "uri", "fname"]).drop_duplicates(["tic", "sector", "author"])


def get(row):
    import lightkurve as lk
    tic, sector, author, uri, fname = row
    p = npz_path(tic, sector, author)
    if os.path.exists(p):
        return "cached"
    dest = os.path.join(CACHE, "mast", "fast", fname)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    for attempt in range(4):
        try:
            if not os.path.exists(dest):
                r = requests.get(URL + uri, timeout=120)
                r.raise_for_status()
                open(dest + ".part", "wb").write(r.content)
                os.replace(dest + ".part", dest)
            t, f = lcwin.prep(lk.read(dest))
            os.makedirs(os.path.dirname(p), exist_ok=True)
            np.savez_compressed(p, t=t, f=f)
            return "ok"
        except Exception as e:
            err = str(e)[:120]
            if os.path.exists(dest) and attempt >= 1:
                os.remove(dest)
            time.sleep(5 * (attempt + 1))
    return "error: " + err


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=12)
    a = ap.parse_args()
    j = needed()
    say(f"files still needed: {len(j)}")
    u = uris(j)
    u.to_csv(os.path.join(CACHE, "uris.csv"), index=False)
    miss = len(j) - len(u)
    say(f"URIs found: {len(u)}; not found on MAST: {miss}")
    st = []
    with ThreadPoolExecutor(a.threads) as ex:
        for i, s in enumerate(ex.map(get, [tuple(r) for r in u.itertuples(index=False)]), 1):
            st.append(s)
            if i % 250 == 0:
                say(f"  downloaded {i}/{len(u)}")
    s = pd.Series(st)
    say(f"done: {(s == 'ok').sum()} ok, {(s == 'cached').sum()} cached, {s.str.startswith('error').sum()} errors")


if __name__ == "__main__":
    main()
