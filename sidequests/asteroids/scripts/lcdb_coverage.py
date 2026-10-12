"""
How underserved are Jupiter Trojans, Hildas and the outer main belt in the LCDB?

Inputs
  - LCDB public summary (Warner, Harris & Pravec), cache/lcdb/lc_summary_pub.txt, from
    https://minplanobs.org/MPInfo/datazips/LCLIST_PUB_CURRENT.zip (fetched 2026-10-11; the file
    inside is still the 2023 Oct 1 release).
  - JPL SBDB query API: all NUMBERED asteroids of classes TJN (Jupiter Trojans), OMB, MBA
    (pdes, a, e, i, H), cached as cache/sbdb_<class>.json.

Population definitions (numbered objects only)
  - Jupiter Trojans: SBDB class TJN.
  - Hildas: a in [3.70, 4.20] au, e < 0.30, i < 20 deg (the usual 3:2-resonance box).
  - Outer main belt (non-Hilda): a in [2.82, 3.70) au (classes MBA/OMB).
  - Middle MB a in [2.50, 2.82), inner MB a < 2.50 (for contrast).

Output: data/lcdb_coverage.csv (population x H bin: N numbered, N with LCDB U>=2, U>=3-,
fraction) and data/lcdb_slow_rotators.csv (period distribution of U>=2 LCDB entries by population).
Usage: python sidequests/asteroids/scripts/lcdb_coverage.py
"""
import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
CACHE = HERE / "cache"
OUT = HERE / "data"
LCDB = CACHE / "lcdb" / "lc_summary_pub.txt"


def read_lcdb(path=LCDB):
    """Parse the fixed-width LCDB summary using the column starts in its header line."""
    lines = path.read_text(encoding="latin-1").splitlines()
    hdr = lines[3]
    names = ["NUMBER", "NAME", "DESIG", "FAM", "PERIOD", "AMIN", "AMAX", "U", "NOTES", "BIN"]
    starts = {}
    for n in names:
        starts[n] = hdr.index(n + " ") if n != "PERIOD" else hdr.index("PERIOD")
    # the 'F' flag column precedes PERIOD; the PERIOD field runs up to the 'P DESC' column
    ends = {"NUMBER": starts["NAME"], "NAME": starts["DESIG"], "DESIG": starts["FAM"],
            "FAM": starts["FAM"] + 8, "PERIOD": hdr.index("P DESC"), "AMIN": starts["AMAX"],
            "AMAX": starts["U"], "U": starts["U"] + 3, "NOTES": starts["NOTES"] + 6,
            "BIN": starts["BIN"] + 4}
    rows = []
    for ln in lines[5:]:
        if not ln.strip():
            continue
        rows.append({n: ln[starts[n]:ends[n]].strip() for n in names})
    df = pd.DataFrame(rows)
    # a '*' flag can sit in the last column of the NUMBER field (601 rows, e.g. 1143 Odysseus);
    # the first version dropped those rows (bug found while picking the pilot sample)
    df["number"] = pd.to_numeric(df.NUMBER.str.replace("*", "", regex=False).str.strip(), errors="coerce")
    df["period_h"] = pd.to_numeric(df.PERIOD, errors="coerce")
    df["amax"] = pd.to_numeric(df.AMAX, errors="coerce")
    # U codes: 1, 1+, 2-, 2, 2+, 3-, 3 ; map to numbers for thresholds
    umap = {"1": 1.0, "1+": 1.3, "2-": 1.7, "2": 2.0, "2+": 2.3, "3-": 2.7, "3": 3.0}
    df["Unum"] = df.U.map(umap)
    return df


def sbdb(cls):
    f = CACHE / f"sbdb_{cls}.json"
    if not f.exists():
        url = ("https://ssd-api.jpl.nasa.gov/sbdb_query.api?fields=pdes,name,a,e,i,H"
               f"&sb-class={cls}&sb-kind=a&sb-ns=n")
        with urllib.request.urlopen(url, timeout=600) as r:
            f.write_bytes(r.read())
    j = json.loads(f.read_text())
    d = pd.DataFrame(j["data"], columns=j["fields"])
    for c in ["a", "e", "i", "H"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["number"] = pd.to_numeric(d.pdes, errors="coerce")
    d["sbclass"] = cls
    return d


def populations():
    d = pd.concat([sbdb(c) for c in ["TJN", "OMB", "MBA"]], ignore_index=True)
    pop = np.full(len(d), "", dtype=object)
    hilda = (d.a >= 3.70) & (d.a <= 4.20) & (d.e < 0.30) & (d.i < 20)
    pop[(d.sbclass != "TJN") & (d.a < 2.50)] = "inner MB"
    pop[(d.sbclass != "TJN") & (d.a >= 2.50) & (d.a < 2.82)] = "middle MB"
    pop[(d.sbclass != "TJN") & (d.a >= 2.82) & (d.a < 3.70)] = "outer MB (a 2.82-3.70)"
    pop[(d.sbclass != "TJN") & hilda] = "Hilda"
    pop[d.sbclass == "TJN"] = "Jupiter Trojan"
    d["population"] = pop
    return d[d.population != ""]


def main():
    lc = read_lcdb()
    lc = lc[lc.number.notna()].sort_values("Unum", ascending=False).drop_duplicates("number")
    pops = populations()
    m = pops.merge(lc[["number", "period_h", "amax", "U", "Unum"]], on="number", how="left")
    m["has_U2"] = m.Unum >= 2.0
    m["has_U3m"] = m.Unum >= 2.7
    m["any_entry"] = m.U.notna() & (m.U != "")
    hb = [-10, 10, 11, 12, 13, 14, 15, 99]
    m["Hbin"] = pd.cut(m.H, hb, right=False, labels=["<10", "10-11", "11-12", "12-13", "13-14",
                                                     "14-15", ">=15"])
    rows = []
    for (p, h), g in m.groupby(["population", "Hbin"], observed=True):
        rows.append(dict(population=p, H_bin=str(h), N_numbered=len(g), N_lcdb_any=int(g.any_entry.sum()),
                         N_U_ge_2=int(g.has_U2.sum()), N_U_ge_3minus=int(g.has_U3m.sum()),
                         frac_U_ge_2=round(g.has_U2.mean(), 4)))
    for p, g in m.groupby("population"):
        rows.append(dict(population=p, H_bin="all", N_numbered=len(g), N_lcdb_any=int(g.any_entry.sum()),
                         N_U_ge_2=int(g.has_U2.sum()), N_U_ge_3minus=int(g.has_U3m.sum()),
                         frac_U_ge_2=round(g.has_U2.mean(), 4)))
    cov = pd.DataFrame(rows)
    OUT.mkdir(exist_ok=True)
    cov.to_csv(OUT / "lcdb_coverage.csv", index=False)
    print(cov[cov.H_bin == "all"].to_string(index=False))
    print()
    print(cov[cov.population.isin(["Jupiter Trojan", "Hilda"])].to_string(index=False))
    # slow rotators among reliable LCDB periods
    s = []
    for p, g in m[m.has_U2].groupby("population"):
        P = g.period_h
        s.append(dict(population=p, N_U_ge_2=len(g), median_P_h=round(P.median(), 2),
                      frac_P_gt_24h=round((P > 24).mean(), 3), frac_P_gt_100h=round((P > 100).mean(), 3)))
    s = pd.DataFrame(s)
    s.to_csv(OUT / "lcdb_slow_rotators.csv", index=False)
    print()
    print(s.to_string(index=False))
    m[["number", "name", "population", "a", "e", "i", "H", "period_h", "amax", "U"]].to_csv(
        CACHE / "population_lcdb_merge.csv.gz", index=False, compression="gzip")


if __name__ == "__main__":
    main()
