"""
Fresh novelty check for the pilot TOIs (run just before any claim).

  1. ExoFOP target JSON, re-fetched now (https://exofop.ipac.caltech.edu/tess/target.php?id=<TIC>&json),
     saved to toi_singles_pilot/data/exofop_<YYYYMMDD>/<TIC>.json. Every non-pipeline planet-parameter row
     (community uploads) is listed with period, user, date, notes; also the TOI rows' current period.
  2. NASA Exoplanet Archive pscomppars, queried live by TIC (TAP), plus pl_name/period.
Output: toi_singles_pilot/novelty_check.csv
Usage:  python toi_singles_pilot/scripts/s04_novelty_check.py [TOI ...]   (default: all pilot TOIs)
"""
import json
import sys
import time
from pathlib import Path

import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
PILOT = HERE.parent
sys.path.insert(0, str(HERE))
import s01_run  # noqa: E402

BASE = "https://exofop.ipac.caltech.edu/tess"
UA = {"User-Agent": "duotransit-research/0.1"}


def exofop(tic, jdir):
    for attempt in range(5):
        try:
            text = requests.get(f"{BASE}/target.php?id={tic}&json", timeout=120, headers=UA).text
            d = json.loads(text)
            (jdir / f"{tic}.json").write_text(text, encoding="utf-8")
            return d
        except Exception:  # noqa: BLE001
            time.sleep(5 * (attempt + 1))
    return None


def pscomppars(tic):
    q = (f"select pl_name,pl_orbper,disc_year,pl_orbper_reflink from pscomppars where tic_id='TIC {tic}'")
    r = requests.get("https://exoplanetarchive.ipac.caltech.edu/TAP/sync", params={"query": q, "format": "json"},
                     timeout=120, headers=UA)
    r.raise_for_status()
    return r.json()


def main():
    jdir = PILOT / "data" / f"exofop_{time.strftime('%Y%m%d')}"
    jdir.mkdir(parents=True, exist_ok=True)
    rows = s01_run.targets("pilot")
    if len(sys.argv) > 1:
        rows = [r for r in rows if r["TOI"] in sys.argv[1:]]
    out = []
    for r in rows:
        tic = int(r["TIC"])
        j = exofop(tic, jdir)
        ups, toi_rows, prov = [], [], None
        for p in (j or {}).get("planet_parameters", []):
            if p.get("prov"):
                prov = p.get("prov_title")
                continue
            line = f"{p.get('name') or p.get('toi')}: P={p.get('per')} by {p.get('puser')} {p.get('pdate')}"
            if str(prov).startswith("TOIs") or str(p.get("puser", "")).lower() in ("", "none"):
                toi_rows.append(line)
            else:
                ups.append(line + f" [{str(p.get('pnotes'))[:120]}]")
        try:
            ps = pscomppars(tic)
        except Exception as e:  # noqa: BLE001
            ps = [dict(error=str(e)[:100])]
        out.append(dict(TOI=r["TOI"], TIC=tic, checked_at=time.strftime("%Y-%m-%d %H:%M"),
                        exofop_ok=j is not None, n_community_uploads=len(ups),
                        community_uploads=" | ".join(ups), toi_rows=" | ".join(toi_rows),
                        n_files=len((j or {}).get("files", [])),
                        n_time_series=len((j or {}).get("time_series", [])),
                        pscomppars=" | ".join(f"{x.get('pl_name')} P={x.get('pl_orbper')}" for x in ps) if ps else ""))
        print(out[-1], flush=True)
        time.sleep(1)
    pd.DataFrame(out).to_csv(PILOT / "novelty_check.csv", index=False)


if __name__ == "__main__":
    main()
