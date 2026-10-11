"""
Triage step 2 (automated literature pass): for each candidate, list the papers SIMBAD
associates with the star (via its TIC identifier), published 2022 or later.
This catches published follow-up (e.g. NGTS/WINE discoveries) where the planet carries a
survey name rather than a TOI number. SIMBAD indexing lags publication, so recent arXiv
papers must still be checked by hand.

Output: data/catalog/simbad_refs.csv  (TIC, main_id, bibcode, year, journal, title)
Raw query results are cached in data/raw/simbad_refs_raw.csv (re-fetched only if missing).

Usage:  python scripts/triage_literature.py
"""
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
CAT = ROOT / "data" / "catalog"
TAP = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"


def tap(query):
    r = requests.post(TAP, data={"request": "doQuery", "lang": "adql", "format": "csv",
                                 "query": query}, timeout=300)
    r.raise_for_status()
    from io import StringIO
    return pd.read_csv(StringIO(r.text))


def main():
    cand = pd.read_csv(CAT / "candidates.csv")
    raw = RAW / "simbad_refs_raw.csv"
    if not raw.exists():
        ids = ",".join(f"'TIC {t}'" for t in cand.TIC)
        q = f"""
        SELECT i.id AS tic_ident, b.main_id, r.bibcode, r."year", r.journal, r.title
        FROM ident AS i
        JOIN basic AS b ON b.oid = i.oidref
        JOIN has_ref AS h ON h.oidref = i.oidref
        JOIN ref AS r ON r.oidbib = h.oidbibref
        WHERE i.id IN ({ids})
        """
        df = tap(q)
        df.to_csv(raw, index=False)
    df = pd.read_csv(raw)
    df["TIC"] = df.tic_ident.str.replace("TIC ", "").astype(int)
    found = set(df.TIC)
    print(f"SIMBAD: {len(found)} of {len(cand)} TICs have references; {len(df)} refs total")
    recent = df[df["year"] >= 2022].sort_values(["TIC", "year"])
    recent = recent[["TIC", "main_id", "bibcode", "year", "journal", "title"]]
    recent.to_csv(CAT / "simbad_refs.csv", index=False)
    print(f"{len(recent)} refs from 2022+ -> {CAT / 'simbad_refs.csv'}")


if __name__ == "__main__":
    main()
