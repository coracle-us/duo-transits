"""
Step 4 (part): ExoFOP target-page JSON for every target star, plus pscomppars by TIC.

Downloads https://exofop.ipac.caltech.edu/tess/target.php?id=<TIC>&json into toi_singles_k2/data/exofop/
(resumable; 1 s pause). Re-uses nothing from data/raw/exofop_targets (never written).
For each TOI: every planet-parameter row on the star that is NOT the SPOC/QLP TOI pipeline row
(user/community uploads), with its period, user, date and notes; and whether the uploaded period is
an integer alias of the TOI's maximum period (|n - round(n)| < 0.002).
Also: TIC in NASA pscomppars (data/raw/nasa_pscomppars.csv), TOI disposition counts.

Output: toi_singles_k2/exofop_uploads.csv (one row per upload) and columns merged into
        toi_singles_k2/targets_status.csv (one row per TOI: n_uploads, uploaded periods, pscomppars).
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "toi_singles_k2"
JDIR = OUT / "data" / "exofop"
BASE = "https://exofop.ipac.caltech.edu/tess"


def fetch(tic):
    js = JDIR / f"{tic}.json"
    if js.exists():
        return json.loads(js.read_text(encoding="utf-8"))
    for attempt in range(5):
        try:
            text = requests.get(f"{BASE}/target.php?id={tic}&json", timeout=120,
                                headers={"User-Agent": "duotransit-research/0.1"}).text
            d = json.loads(text)
            js.write_text(text, encoding="utf-8")
            time.sleep(1)
            return d
        except Exception:  # noqa: BLE001
            time.sleep(5 * (attempt + 1))
    return None


def main():
    JDIR.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(OUT / "targets.csv")
    ps = pd.read_csv(ROOT / "data" / "raw" / "nasa_pscomppars.csv")
    ps_tic = set(pd.to_numeric(ps.tic_id.astype(str).str.replace("TIC ", ""), errors="coerce").dropna().astype(int))
    ups, stat = [], []
    for r in d.itertuples():
        j = fetch(int(r.TIC))
        if j is None:
            stat.append(dict(TOI=r.TOI, TIC=r.TIC, exofop_ok=False))
            continue
        prov, my = None, []
        for p in j.get("planet_parameters", []):
            if p.get("prov"):
                prov = p.get("prov_title")
                continue
            if str(prov).startswith("TOIs") or str(p.get("puser", "")).lower() in ("", "none"):
                continue
            per = pd.to_numeric(p.get("per"), errors="coerce")
            toi_name = str(p.get("toi") or p.get("name") or "")
            Pm = r.max_period_d
            n = Pm / per if (np.isfinite(Pm) and per and per > 0) else np.nan
            row = dict(TOI=r.TOI, TIC=r.TIC, provenance=prov, name=p.get("name"), toi_field=toi_name,
                       period=per, period_err=p.get("per_e"), epoch=p.get("epoch"), depth_ppm=p.get("dep_p"),
                       dur_h=p.get("dur"), user=p.get("puser"), date=p.get("pdate"), notes=p.get("pnotes"),
                       n_alias_of_maxP=round(n, 4) if n == n else np.nan,
                       integer_alias=bool(n == n and abs(n - round(n)) < 0.002))
            ups.append(row)
            my.append(row)
        same = [u for u in my if str(r.TOI).split(".")[0] in str(u["toi_field"]) + str(u["name"])] or my
        pers = sorted({round(u["period"], 4) for u in same if u["period"] == u["period"] and u["period"] > 0})
        stat.append(dict(TOI=r.TOI, TIC=r.TIC, exofop_ok=True, n_uploads=len(my),
                         uploaded_periods=" ".join(map(str, pers)),
                         latest_upload=max([str(u["date"]) for u in my], default=""),
                         upload_users=" ".join(sorted({str(u["user"]) for u in my})),
                         n_files=len(j.get("files", [])), n_time_series=len(j.get("time_series", [])),
                         n_spectra=len(j.get("spectroscopy", [])), n_imaging=len(j.get("imaging", [])),
                         in_pscomppars=int(r.TIC) in ps_tic))
        print(r.TOI, stat[-1].get("n_uploads"), stat[-1].get("uploaded_periods"), flush=True)
    pd.DataFrame(ups).to_csv(OUT / "exofop_uploads.csv", index=False)
    s = pd.DataFrame(stat)
    s.to_csv(OUT / "targets_status.csv", index=False)
    print(f"ExoFOP JSON ok: {s.exofop_ok.sum()}/{len(s)}; with >=1 community planet-parameter upload: "
          f"{(s.n_uploads > 0).sum()}; with an uploaded period: {(s.uploaded_periods.fillna('') != '').sum()}; "
          f"in pscomppars: {s.in_pscomppars.sum()}")


if __name__ == "__main__":
    main()
