"""
tess_wide step 6 (PLAN.md section 8): which CNN-vetted single events are already known.

Events cross-matched: real events with CNN verdict "pass" or "out_of_range" (counted separately).
Known if any of:
  TOI (any disposition), CTOI           same TIC; ExoFOP tables downloaded on the run day
  confirmed planet                      same TIC; NASA Exoplanet Archive pscomppars (TAP), run day
  TESS EB catalogue (Prsa+2022)         within 42" (VizieR J/ApJS/258/16)
  Gaia DR3 eclipsing binaries           within 42" (VizieR I/358/veb)
  VSX eclipsing types E, EA, EB, EW     within 42" (VizieR B/vsx/vsx; ELL is not eclipsing)
  Hawthorn+2024 duotransits             same TIC (data/catalog/candidates.csv)
  Salinas+2026 singles                  same TIC (salinas-second-transits/candidates.csv)
VizieR goes through the CfA mirror helper of the K2 revision (k2singles/revision2/vizier_targeted.py).

Output (outside git): tess_wide/data/s<NNN>_crossmatch.csv (names candidate stars: never commit),
                      s<NNN>_crossmatch.json (counts), catalogue snapshots in tess_wide/data/catalogs/
"""
import argparse
import datetime as dt
import json
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = HERE / "data"
CAT = DATA / "catalogs"
sys.path.insert(0, str(ROOT / "k2singles" / "revision2"))
import vizier_targeted as vt  # noqa: E402

R_ARCSEC = 42.0
VIZ = {"prsa2022": "J/ApJS/258/16", "gaia_veb": "I/358/veb", "vsx": "B/vsx/vsx"}
VSX_ECL = re.compile(r"^E($|A|B|W|[:/|+ ])")


def fetch(url, out):
    """curl, not urllib: on 2026-11-25 ExoFOP answered curl in ~20 s but never finished a urllib read."""
    if not out.exists():
        tmp = out.with_suffix(".tmp")
        subprocess.run(["curl", "-s", "-f", "--max-time", "600", "-o", str(tmp), url], check=True)
        tmp.rename(out)
    return out


def tic_sets(day):
    CAT.mkdir(parents=True, exist_ok=True)
    toi = pd.read_csv(fetch("https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv",
                            CAT / f"toi_{day}.csv"), comment="#")
    ctoi = pd.read_csv(fetch("https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv",
                             CAT / f"ctoi_{day}.csv"), comment="#")
    q = urllib.parse.urlencode({"query": "select distinct tic_id from pscomppars where tic_id is not null",
                                "format": "csv"})
    ps = pd.read_csv(fetch(f"https://exoplanetarchive.ipac.caltech.edu/TAP/sync?{q}", CAT / f"pscomppars_tic_{day}.csv"))
    haw = pd.read_csv(ROOT / "data" / "catalog" / "candidates.csv")
    sal = pd.read_csv(ROOT / "salinas-second-transits" / "candidates.csv")
    tcol = lambda d: next(c for c in d.columns if c.strip().lower() in ("tic id", "tic", "tic_id"))
    num = lambda s: set(pd.to_numeric(s.astype(str).str.replace("TIC", "").str.strip(), errors="coerce").dropna().astype(int))
    return {"toi": num(toi[tcol(toi)]), "ctoi": num(ctoi[tcol(ctoi)]),
            "confirmed": num(ps.tic_id), "hawthorn": num(haw.TIC), "salinas": num(sal.tic)}


def cone(cat, ra, dec):
    for attempt in range(4):
        v = vt.get({"-source": cat, "-c": f"{ra:.6f} {dec:+.6f}", "-c.rs": R_ARCSEC, "-out.max": 200,
                    "-out.all": 1, "-out.add": "_r"})
        if v is not None:
            tabs = [t.to_table(use_names_over_ids=True) for t, _ in vt.tables(v)]
            tabs = [t.to_pandas() for t in tabs if len(t)]
            return pd.concat(tabs, ignore_index=True) if tabs else pd.DataFrame()
        time.sleep(30 * (attempt + 1))
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", type=int, default=104)
    s = ap.parse_args().sector
    day = dt.date.today().isoformat()
    cnn = pd.read_csv(DATA / f"s{s:04d}_cnn_real.csv")
    ev = cnn[cnn.cnn.isin(["pass", "out_of_range"])].copy()
    stars = pd.read_csv(DATA / f"s{s:04d}_stars.csv").set_index("id")
    sets = tic_sets(day)
    rows = []
    for e in ev.itertuples():
        r = dict(tic=e.tic, T=e.T, cnn=e.cnn)
        for k, ids in sets.items():
            r[k] = e.tic in ids
        ra, dec = stars.loc[e.tic, "ra"], stars.loc[e.tic, "dec"]
        for k, cat in VIZ.items():
            t = cone(cat, ra, dec)
            if t is None:
                r[k] = "query failed"
                continue
            if k == "vsx" and len(t):
                tcol = next((c for c in t.columns if c.lower() == "type"), None)
                t = t[t[tcol].astype(str).str.match(VSX_ECL)] if tcol else t.iloc[0:0]
            r[k] = bool(len(t))
            if len(t):
                r[f"{k}_sep"] = float(np.nanmin(t["_r"])) if "_r" in t else np.nan
        r["known"] = any(r[k] is True for k in list(sets) + list(VIZ))
        r["query_failed"] = any(r[k] == "query failed" for k in VIZ)
        rows.append(r)
        print(e.tic, {k: r[k] for k in list(sets) + list(VIZ)}, flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(DATA / f"s{s:04d}_crossmatch.csv", index=False)
    summ = dict(run_day=day, catalog_sizes={k: len(v) for k, v in sets.items()})
    for verdict in ("pass", "out_of_range"):
        o = out[out.cnn == verdict] if len(out) else out
        summ[verdict] = dict(events=len(o), known=int(o.known.sum()) if len(o) else 0,
                             new=int((~o.known).sum()) if len(o) else 0,
                             query_failed=int(o.query_failed.sum()) if len(o) else 0,
                             by_catalogue={k: int((o[k] == True).sum()) for k in list(sets) + list(VIZ)} if len(o) else {})  # noqa: E712
    (DATA / f"s{s:04d}_crossmatch.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
