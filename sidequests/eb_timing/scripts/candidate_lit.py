"""
Per-candidate literature check (required for every ETV candidate before it is reported).

For each TIC (coordinates from Prsa+2022, or --ra/--dec style rows in data/ncvz_*_ephem are not
needed: the TIC catalogue at MAST is queried when the star is not in Prsa):
  1. SIMBAD (object cross-identified with 'TIC n', else the nearest within 10"): main id, object type, all identifiers, full reference list (bibcodes).
  2. Tokovinin Multiple Star Catalog (VizieR J/ApJS/235/6) within 20 arcsec.
  3. Washington Double Star catalogue (VizieR B/wds/wds) within 10 arcsec.
  4. AAVSO VSX (B/vsx/vsx) within 10 arcsec.
  5. data/known_multiples.csv (TICs harvested from 27 triple/ETV papers) and data/b26_triples.csv.
  5b. SIMBAD reference titles are scanned for the timing words below as well.
  6. ADS full text (scripts/ads_check.py: token from ADS_DEV_KEY or the Windows user
     environment, never printed): full:"<id>" OR ... over the TIC number and the SIMBAD
     identifiers of types TIC/HD/HIP/TYC/BD/CD/CPD/GSC/UCAC4/2MASS/Gaia DR3/V*/ASAS/NSV.
     'already timed' flag (the CATS-style check): any hit whose title/abstract/body highlight
     mentions O-C, minima, eclipse timing, period change/variation, light-time or third body.
Positive-control caveat (notes/ads_fulltext_check.md): ADS full text finds papers ABOUT a star
reliably, but mostly misses stars that are only rows in big tables.

Usage:  python sidequests/eb_timing/scripts/candidate_lit.py TIC [TIC ...]
Output: notes/candidate_lit/<TIC>.json and data/candidate_lit.csv (rows updated per TIC)
"""
import json
import re
import sys
import time
from pathlib import Path

import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ads_check  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
OUTD = HERE / "notes" / "candidate_lit"
OUTC = HERE / "data" / "candidate_lit.csv"
ID_TYPES = ("TIC", "HD", "HIP", "TYC", "BD", "CD", "CPD", "GSC", "UCAC4", "2MASS", "Gaia DR3", "V*", "ASAS", "NSV")
TIMING_WORDS = re.compile(r"O\s*[-−–]\s*C|minim|eclipse tim|timing|period (change|variation)|light[- ]time|"
                          r"LTTE|LITE|third bod|tertiary|period (analysis|study)|ephemeris", re.I)


def coords(tic, cat):
    if tic in cat.index:
        return SkyCoord(cat.loc[tic, "RAJ2000"] * u.deg, cat.loc[tic, "DEJ2000"] * u.deg)
    from astroquery.mast import Catalogs
    r = Catalogs.query_criteria(catalog="TIC", ID=int(tic)).to_pandas()
    return SkyCoord(r.ra.iloc[0] * u.deg, r.dec.iloc[0] * u.deg)


def viz(cat_id, c, radius):
    try:
        t = Vizier(row_limit=10, columns=["**"]).query_region(c, radius=radius * u.arcsec, catalog=cat_id)
        return t[0].to_pandas().astype(str).to_dict("records") if len(t) else []
    except Exception as e:  # noqa: BLE001
        return [{"error": str(e)[:100]}]


def check(tic, cat, tok):
    c = coords(tic, cat)
    out = {"tic": tic, "ra": c.ra.deg, "dec": c.dec.deg}
    sim = Simbad()
    sim.add_votable_fields("otype", "nbref")
    ids = [str(tic)]
    try:
        s = sim.query_object(f"TIC {tic}")           # exact cross-id first
        if s is None or not len(s):
            s = sim.query_region(c, radius=10 * u.arcsec)
            if s is not None and len(s):               # nearest object, not SIMBAD's first row
                sp = s.to_pandas()
                sp.columns = [x.lower() for x in sp.columns]
                sc = SkyCoord(sp.ra.values * u.deg, sp.dec.values * u.deg)
                s = s[[int(c.separation(sc).argmin())]]
        if s is not None and len(s):
            s = s.to_pandas()
            s.columns = [x.lower() for x in s.columns]
            mid = str(s.iloc[0]["main_id"])
            out.update(simbad_id=mid, simbad_otype=str(s.iloc[0].get("otype")), simbad_nbref=int(s.iloc[0].get("nbref") or 0))
            q = mid.replace("'", "''")
            refs = Simbad.query_tap("SELECT b.bibcode, b.title FROM basic AS o JOIN has_ref AS h ON h.oidref=o.oid "
                                    f"JOIN ref AS b ON b.oidbib=h.oidbibref WHERE o.main_id='{q}'").to_pandas()
            out["simbad_refs"] = [f"{r.bibcode} {str(r.title)[:100]}" for r in refs.itertuples()]
            out["simbad_timing_refs"] = [f"{r.bibcode} {str(r.title)[:100]}" for r in refs.itertuples()
                                         if TIMING_WORDS.search(str(r.title))]
            idt = Simbad.query_objectids(mid)
            allids = [str(x) for x in idt[idt.colnames[0]]] if idt is not None else []
            out["simbad_ids"] = allids
            for i in allids:
                i = re.sub(r"\s+", " ", i).strip()
                if i.startswith(ID_TYPES):
                    if i.startswith(("TIC", "Gaia DR3")):
                        i = i.replace("TIC ", "").replace("Gaia DR3 ", "")
                    elif i.startswith("V* "):
                        i = i[3:]                      # 'V* SX Dra' -> 'SX Dra'
                    ids.append(i)
    except Exception as e:  # noqa: BLE001
        out["simbad_err"] = str(e)[:100]
    out["msc"] = viz("J/ApJS/235/6", c, 20)
    out["wds"] = viz("B/wds/wds", c, 10)
    out["vsx"] = viz("B/vsx/vsx", c, 10)
    known = pd.read_csv(HERE / "data" / "known_multiples.csv").set_index("TIC").papers
    out["known_multiples"] = str(known.get(tic, ""))
    b26 = pd.read_csv(HERE / "data" / "b26_triples.csv")
    out["b26"] = b26[b26.tic == tic].group.tolist()
    ids = list(dict.fromkeys(ids))
    out["ads_ids"] = ids
    try:
        q = " OR ".join(f'full:"{i}"' for i in ids)
        n, docs = ads_check.query(q, tok)
        out["ads_n"] = n
        out["ads_docs"] = docs
        out["ads_timing_hits"] = [d["bibcode"] for d in docs
                                  if TIMING_WORDS.search(d["title"] + " " + " ".join(d["highlights"]))]
    except Exception as e:  # noqa: BLE001
        out["ads_err"] = str(e)[:100]
    return out


def main():
    cat = pd.read_csv(HERE / "data" / "prsa2022_tess_ebs.csv").drop_duplicates("TIC").set_index("TIC")
    tok = ads_check.token()
    OUTD.mkdir(parents=True, exist_ok=True)
    old = pd.read_csv(OUTC) if OUTC.exists() else pd.DataFrame(columns=["tic"])
    rows = {int(r["tic"]): r for r in old.to_dict("records")}
    for tic in [int(x) for x in sys.argv[1:]]:
        r = check(tic, cat, tok)
        (OUTD / f"{tic}.json").write_text(json.dumps(r, indent=1, default=str), encoding="utf-8")
        rows[tic] = dict(tic=tic, simbad_id=r.get("simbad_id"), otype=r.get("simbad_otype"),
                         n_refs=r.get("simbad_nbref"), n_msc=len([x for x in r["msc"] if "error" not in x]),
                         n_wds=len([x for x in r["wds"] if "error" not in x]),
                         vsx=";".join(x.get("Name", "") for x in r["vsx"]), known=r["known_multiples"],
                         b26=";".join(r["b26"]), ads_n=r.get("ads_n"),
                         ads_timing_hits=";".join(r.get("ads_timing_hits", [])),
                         simbad_timing_refs=" | ".join(r.get("simbad_timing_refs", [])))
        print(rows[tic], flush=True)
        time.sleep(1)
    pd.DataFrame(list(rows.values())).to_csv(OUTC, index=False)


if __name__ == "__main__":
    main()
