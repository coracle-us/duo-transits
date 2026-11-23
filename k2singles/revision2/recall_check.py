"""
AAS82225 revision, literature search: the pre-registered recall check (PLAN.md, "Literature classification: rules").
For each known (star, paper) pair from the submitted note, does each source find it?
  ADS     hits.csv rows with found_by in (ads, both)
  SIMBAD  hits.csv rows with found_by in (simbad, both)
  VizieR  vizier.csv rows whose catalogue is the paper's VizieR catalogue (J/<journal>/<vol>/<page>)
  VizieR fixed list  the same in vizier_targeted.csv (post-hoc, vizier_targeted.py). That list holds these papers'
          catalogues by construction, so this column only shows whether the positional query works; it is not recall.
Nothing is tuned on this; the numbers go into the revised note as the search's known gaps.
Usage:  python k2singles/revision2/recall_check.py
Output: stdout, and k2singles/revision2/data/literature/recall.csv
"""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LIT = HERE / "data" / "literature"

# bibcodes from the submitted note's reference list (k2singles/rnaas/rnaas.tex); VizieR ids follow from them
PAPERS = {
    "LaCourse & Jacobs 2018": ("2018RNAAS...2...28L", None),
    "Sandquist+2021": ("2021AJ....161...59S", "J/AJ/161/59"),
    "IJspeert+2024": ("2024A&A...691A.242I", "J/A+A/691/A242"),
    "Valizadegan+2025": ("2025AJ....170..287V", "J/AJ/170/287"),
    "Kovalev+2024": ("2024MNRAS.527..521K", "J/MNRAS/527/521"),
    "Guo+2025": ("2025ApJS..278...46G", "J/ApJS/278/46"),
    "Kovacs 2020": ("2020A&A...643A.169K", "J/A+A/643/A169"),
    "Libralato+2016": ("2016MNRAS.463.1780L", "J/MNRAS/463/1780"),
}
PAIRS = [("Sandquist+2021", 211411112), ("IJspeert+2024", 202126877), ("IJspeert+2024", 202072917),
         ("Valizadegan+2025", 202072917), ("Valizadegan+2025", 211489484), ("Kovalev+2024", 202071902),
         ("Guo+2025", 202071902), ("Kovacs 2020", 211503363), ("Libralato+2016", 211892898)]


def vizier_rows():
    # always from the cache: vizier.csv is written at the end of a run and lists only the stars that succeeded, so it
    # cannot tell an unqueried star from one with no match
    files = sorted((LIT / "vizier_cache").glob("*.json"))
    rows = [r for x in files for r in json.loads(x.read_text())]
    return pd.DataFrame(rows, columns=["EPIC", "TIC", "catalog", "title", "periods"]), {int(x.stem) for x in files}


def main():
    stars = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv").drop_duplicates("EPIC").EPIC.astype(int)
    h = pd.read_csv(LIT / "hits.csv")
    v, done = vizier_rows()
    vt, vt_done, vt_exists = None, set(), set()
    if (LIT / "vizier_targeted.csv").exists():
        vt = pd.read_csv(LIT / "vizier_targeted.csv")
        vt_done = {int(x.stem) for x in (LIT / "vizier_targeted_cache").glob("*.json")}
        c = pd.read_csv(LIT / "vizier_catalogues.csv")
        vt_exists = set(c.catalog[c.exists.astype(bool)])
    print(f"VizieR all-catalogue census: {len(done)} of {len(stars)} stars queried so far")
    out = []

    def check(paper, epic):
        bib, viz = PAPERS[paper]
        s = h[(h.EPIC == epic) & (h.bibcode == bib)]
        ads = bool(s.found_by.isin(["ads", "both"]).any())
        sim = bool(s.found_by.isin(["simbad", "both"]).any())
        if viz is None:
            vz = "no catalogue"
        elif epic not in done:
            vz = "not queried yet"
        else:
            vz = bool(v[(v.EPIC == epic) & v.catalog.astype(str).str.startswith(viz + "/")].shape[0])
        if viz is None:
            fx = "no catalogue"
        elif vt is None or epic not in vt_done:
            fx = "not queried yet"
        elif viz not in vt_exists:
            fx = "catalogue not on VizieR"
        else:
            fx = bool(vt[(vt.EPIC == epic) & vt.catalog.astype(str).str.startswith(viz + "/")].shape[0])
        out.append(dict(paper=paper, EPIC=epic, ads=ads, simbad=sim, vizier=vz, vizier_fixed_list=fx))

    for e in stars:
        check("LaCourse & Jacobs 2018", int(e))
    for p, e in PAIRS:
        check(p, e)
    r = pd.DataFrame(out)
    r.to_csv(LIT / "recall.csv", index=False)
    lj = r[r.paper == "LaCourse & Jacobs 2018"]
    print(f"LaCourse & Jacobs 2018 (all {len(lj)} stars): ADS {lj.ads.sum()}, SIMBAD {lj.simbad.sum()}, "
          f"either {(lj.ads | lj.simbad).sum()}")
    k = r[r.paper != "LaCourse & Jacobs 2018"]
    print(k.to_string(index=False))
    vz = k[k.vizier.apply(lambda x: isinstance(x, bool))]
    print(f"other known pairs ({len(k)}): ADS {k.ads.sum()}, SIMBAD {k.simbad.sum()}, "
          f"VizieR {sum(x is True for x in vz.vizier)} of {len(vz)} queried, any source "
          f"{(k.ads | k.simbad | k.vizier.apply(lambda x: x is True)).sum()}")
    fx = k[k.vizier_fixed_list.apply(lambda x: isinstance(x, bool))]
    print(f"VizieR fixed list (not recall: the list holds these catalogues): {sum(x is True for x in fx.vizier_fixed_list)} "
          f"of {len(fx)} pairs whose catalogue is on VizieR; "
          f"{(k.vizier_fixed_list == 'catalogue not on VizieR').sum()} pairs have no catalogue there")


if __name__ == "__main__":
    main()
