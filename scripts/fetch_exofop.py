"""
Triage step 3: per-target ExoFOP records for all 85 candidates.

ExoFOP target pages carry information the bulk TOI/CTOI tables do not: community-uploaded
planet parameters (often with a refined period and a note such as "Period constrained by
transits in sectors ..."), observing notes, and follow-up file listings.

Downloads (resumable) into data/raw/exofop_targets/:
    <TIC>.json           https://exofop.ipac.caltech.edu/tess/target.php?id=<TIC>&json
    <TIC>_obsnotes.html  https://exofop.ipac.caltech.edu/tess/edit_obsnotes.php?id=<TIC>
Writes data/catalog/exofop_planet_params.csv: every planet-parameter row per TIC, with the
alias number n = dT / P for each period (integer n => consistent with the duotransit).

Usage:  python scripts/fetch_exofop.py
"""
import json
import re
import time
from html import unescape
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "raw" / "exofop_targets"
CAT = ROOT / "data" / "catalog"
BASE = "https://exofop.ipac.caltech.edu/tess"


def get(url):
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=120, headers={"User-Agent": "duotransit-research/0.1"})
            r.raise_for_status()
            return r
        except Exception as e:
            print(f"  retry {attempt + 1}: {e}")
            time.sleep(5)
    raise RuntimeError(f"failed: {url}")


def fetch(tic):
    js = OUT / f"{tic}.json"
    if not js.exists():
        # ExoFOP intermittently returns an empty body; retry until it parses
        for attempt in range(5):
            text = get(f"{BASE}/target.php?id={tic}&json").text
            try:
                json.loads(text)
                break
            except json.JSONDecodeError:
                print(f"  TIC {tic}: empty/invalid JSON, retry {attempt + 1}")
                time.sleep(10)
        else:
            raise RuntimeError(f"TIC {tic}: no valid JSON from ExoFOP")
        js.write_text(text, encoding="utf-8")
        time.sleep(1)
    on = OUT / f"{tic}_obsnotes.html"
    if not on.exists():
        on.write_text(get(f"{BASE}/edit_obsnotes.php?id={tic}").text, encoding="utf-8")


def obsnotes_text(tic):
    html = (OUT / f"{tic}_obsnotes.html").read_text(encoding="utf-8")
    txt = unescape(re.sub(r"<[^>]+>", " ", html))
    return re.sub(r"\s+", " ", txt).strip()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cand = pd.read_csv(CAT / "candidates.csv").set_index("TIC")
    for tic in cand.index:
        fetch(tic)
    rows = []
    for tic, c in cand.iterrows():
        d = json.loads((OUT / f"{tic}.json").read_text(encoding="utf-8"))
        prov = None
        for p in d.get("planet_parameters", []):
            if "prov" in p:
                prov = p["prov_title"]
                continue
            per = pd.to_numeric(p.get("per"), errors="coerce")
            ep = pd.to_numeric(p.get("epoch"), errors="coerce")
            n = c.dT / per if per and per > 0 else float("nan")
            rows.append(dict(
                TIC=tic, provenance=prov, name=p.get("name"), toi=p.get("toi"), disp=p.get("disp"),
                period=per, period_err=p.get("per_e"), epoch_btjd=ep - 2457000 if ep > 2e6 else ep,
                depth_ppm=p.get("dep_p"), dur_h=p.get("dur"), n_alias=round(n, 4),
                is_alias=abs(n - round(n)) < 0.002 if n == n else False,
                user=p.get("puser"), date=p.get("pdate"), notes=p.get("pnotes"),
            ))
        files = d.get("files", [])
        notes = obsnotes_text(tic)
        rows.append(dict(TIC=tic, provenance="(summary)", name="", n_files=len(files),
                         n_time_series=len(d.get("time_series", [])),
                         obsnotes_len=len(notes)))
    out = pd.DataFrame(rows)
    out.to_csv(CAT / "exofop_planet_params.csv", index=False)
    up = out[(out.provenance != "(summary)") & ~out.provenance.astype(str).str.startswith("TOIs")]
    print(f"{len(up)} user-uploaded parameter rows for {up.TIC.nunique()} TICs")


if __name__ == "__main__":
    main()
