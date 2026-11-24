"""Has anyone else already reported a period for the four periodic Salinas stars?

Read-only checks, run on the date in OUT:
- ExoFOP TOI and CTOI tables (fresh download): any TOI or CTOI on these TICs.
- ExoFOP per-target JSON: TOIs, CTOIs, uploaded files, observing notes, tags.
- NASA Exoplanet Archive TAP: pscomppars (confirmed) and toi tables by TIC.

Usage: python round3/s13_novelty.py
Writes results/round3/novelty_2026-11-22.csv and prints a summary. Raw downloads go
to work/round3_novelty/ (not committed).
"""
import io
import json
import time
from pathlib import Path

import pandas as pd
import requests

HERE = Path(__file__).resolve().parents[1]
RAW = HERE / "work" / "round3_novelty"
OUT = HERE / "results" / "round3" / "novelty_2026-11-22.csv"
TICS = [275647095, 241550976, 107113048, 230133616]
UA = {"User-Agent": "duo-transits research script (read-only)"}

EXOFOP = "https://exofop.ipac.caltech.edu/tess"
TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"


def get(url, name, params=None):
    RAW.mkdir(parents=True, exist_ok=True)
    path = RAW / name
    r = requests.get(url, params=params, headers=UA, timeout=120)
    r.raise_for_status()
    path.write_bytes(r.content)
    time.sleep(1)
    return r.content


def main():
    rows = []
    toi = pd.read_csv(io.BytesIO(get(f"{EXOFOP}/download_toi.php",
                                     "exofop_toi.csv", {"sort": "toi", "output": "csv"})))
    ctoi = pd.read_csv(io.BytesIO(get(f"{EXOFOP}/download_ctoi.php",
                                      "exofop_ctoi.csv", {"sort": "ctoi", "output": "csv"})))
    print(f"ExoFOP TOI rows {len(toi)}, CTOI rows {len(ctoi)}")
    ticcol_toi = [c for c in toi.columns if c.lower().startswith("tic id")][0]
    ticcol_ctoi = [c for c in ctoi.columns if "tic id" in c.lower()][0]

    for tic in TICS:
        row = {"tic": tic}
        t = toi[toi[ticcol_toi] == tic]
        c = ctoi[ctoi[ticcol_ctoi] == tic]
        row["exofop_toi"] = ";".join(str(x) for x in t["TOI"]) if len(t) else ""
        row["exofop_toi_period"] = ";".join(str(x) for x in t["Period (days)"]) if len(t) else ""
        row["exofop_ctoi"] = ";".join(str(x) for x in c.iloc[:, 0]) if len(c) else ""

        js = json.loads(get(f"{EXOFOP}/target.php", f"exofop_target_{tic}.json",
                            {"id": tic, "json": ""}))
        row["target_tois"] = len(js.get("tois", []) or [])
        row["target_ctois"] = len(js.get("ctois", []) or [])
        files = js.get("files", []) or []
        row["target_files"] = len(files)
        row["target_file_names"] = ";".join(f.get("fname", "") for f in files)[:300]
        notes = js.get("obsnotes", []) or js.get("notes", []) or []
        row["target_notes"] = len(notes)
        row["target_tags"] = ";".join(str(x.get("tag", x)) for x in (js.get("tags", []) or []))
        for key in ["imaging", "spectroscopy", "time_series", "stellar_companions"]:
            row[f"target_{key}"] = len(js.get(key, []) or [])

        for table, cols in [("pscomppars", "pl_name,pl_orbper,disc_refname"),
                            ("toi", "toi,tfopwg_disp,pl_orbper")]:
            q = f"select {cols} from {table} where tid={tic}" if table == "toi" else \
                f"select {cols} from {table} where tic_id='TIC {tic}'"
            txt = get(TAP, f"nea_{table}_{tic}.csv", {"query": q, "format": "csv"}).decode()
            n = max(0, len(txt.strip().splitlines()) - 1)
            row[f"nea_{table}_rows"] = n
            row[f"nea_{table}"] = txt.strip().splitlines()[1] if n else ""
        rows.append(row)
        print(tic, {k: v for k, v in row.items() if v not in ("", 0)})

    pd.DataFrame(rows).to_csv(OUT, index=False)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
