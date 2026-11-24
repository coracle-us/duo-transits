"""
AAS82225 revision, literature search part 2b: VizieR queries of a FIXED list of catalogues. POST-HOC (PLAN.md, "Things
decided after first results", 2026-11-22): the pre-registered all-catalogue census (vizier_all.py) could not run, because
the CfA mirror's catalogue database ("metaviz@CFA") was down most of 2026-11-22 and an all-catalogue reply takes several
minutes when it is up. vizier_all.py keeps running serially; if it finishes, it is the primary result.

Catalogues (CATALOGUES below), chosen before any of them was queried:
  the VizieR catalogues of the papers in the recall check (recall_check.PAPERS) and Wang et al. 2021 (RAA 21, 292);
  variable-star and period catalogues (VSX, Gaia DR3 variability and non-single stars, ASAS-SN, ATLAS, ZTF, WISE, TESS EBs);
  K2 variability, eclipsing-binary, rotation and planet-candidate catalogues.
Plus, for the Table 1 stars only, the VizieR catalogue J/<journal>/<vol>/<page> of every paper in table1_hits.csv
(online-only tables: Wang et al. 2021's Table 3 was found only there).
Because the list holds the recall papers' own catalogues, VizieR recall on the known pairs measures only whether the query
works, not how well a blind search would do.

Each catalogue is first checked to exist (a query without a position must return a table). A reply counts only if it has
no database error, ends in </VOTABLE> and (per star) contains the star's Gaia DR3 or 2MASS row, queried alongside as controls;
otherwise it is retried (up to 30 times, 1-5 min apart). A catalogue counts as absent only if the mirror says "not found"; one whose reply lists no tables is recorded as "no
tables in a cone search" and not queried.
Server: vizier.cfa.harvard.edu (the CDS host and the other mirrors were unreachable on 2026-11-22).
Output: data/literature/vizier_targeted.csv (one row per star and catalogue table), vizier_catalogues.csv. Resumable.
Usage:  python k2singles/revision2/vizier_targeted.py
"""
import io
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
from astropy.io.votable import parse

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "data" / "literature"
CACHE = OUT / "vizier_targeted_cache"
URL = "https://vizier.cfa.harvard.edu/viz-bin/votable"
GAIA = "I/355/gaiadr3"
TMASS = "II/246/out"  # 2MASS (epoch ~2000), the control for high proper-motion stars: EPIC 211498244 moves 346 mas/yr and
                      # its Gaia DR3 (epoch 2016) position is 5.2" from the catalogue's J2000 position
PCOLS = ("per", "period", "p", "porb", "plc", "pdays", "p_orb", "p1", "pbls", "ptess", "prot", "freq")

CATALOGUES = {
    # recall-check papers (recall_check.PAPERS) and Wang et al. 2021
    "J/AJ/161/59": "Sandquist+2021, M67 eclipsing binaries",
    "J/A+A/691/A242": "IJspeert+2024, Gaia/TESS eclipsing binaries",
    "J/AJ/170/287": "Valizadegan+2025",
    "J/MNRAS/527/521": "Kovalev+2024",
    "J/ApJS/278/46": "Guo+2025",
    "J/A+A/643/A169": "Kovacs 2020, K2 M67",
    "J/MNRAS/463/1780": "Libralato+2016, K2 M67 variables",
    "J/other/RAA/21.292": "Wang+2021, M67 variables",
    # variable stars and periods
    "B/vsx": "AAVSO VSX",
    "I/358": "Gaia DR3 variability",
    "I/357": "Gaia DR3 non-single stars",
    "II/366": "ASAS-SN variable stars",
    "J/AJ/156/241": "ATLAS variable stars",
    "J/ApJS/249/18": "ZTF periodic variables (Chen+2020)",
    "J/ApJS/237/28": "WISE periodic variables (Chen+2018)",
    "J/ApJS/258/16": "TESS eclipsing binaries (Prsa+2022)",
    # K2
    "J/MNRAS/456/2260": "K2 variable stars (Armstrong+2016)",
    "J/MNRAS/452/3561": "K2 C0 eclipsing binaries (LaCourse+2015)",
    "J/A+A/635/A43": "K2 rotation periods (Reinhold & Hekker 2020)",
    "J/ApJS/244/11": "K2 planet candidates C5-C8 (Kruse+2019)",
    "J/AJ/155/136": "K2 planet candidates (Mayo+2018)",
    "J/ApJS/222/14": "K2 planet candidates C0-C5 (Vanderburg+2016)",
    "J/ApJS/226/7": "K2 planet candidates C0-C4 (Crossfield+2016)",
    "J/A+A/594/A100": "K2 planet candidates C1-C6 (Barros+2016)",
    "J/MNRAS/461/3399": "K2 planet candidates C5-C6 (Pope+2016)",
    "J/AJ/155/21": "K2 planet candidates C5-C8 (Petigura+2018)",
    "J/AJ/161/259": "K2 planet candidates C5-C18 (Zink+2021)",
}
JOURNALS = {"AJ": "AJ", "ApJ": "ApJ", "ApJS": "ApJS", "MNRAS": "MNRAS", "PASP": "PASP", "A&A": "A+A", "A&AS": "A+AS",
            "AcA": "AcA", "PASJ": "PASJ"}


def bib_to_vizier(b):
    """2021AJ....161...59S -> J/AJ/161/59; A&A article ids (A242) kept; None for journals VizieR does not file this way."""
    m = re.fullmatch(r"\d{4}([A-Za-z&]+)\.*(\d+)([A-Z]?)\.*(\d+)[A-Z.]", b.strip())
    if m and m.group(1) == "RAA":
        return f"J/other/RAA/{m.group(2)}.{m.group(4)}"
    if not m or m.group(1) not in JOURNALS:
        return None
    return f"J/{JOURNALS[m.group(1)]}/{m.group(2)}/{m.group(3)}{m.group(4)}"


def get(params):
    """One VizieR request; returns the parsed VOTable, or None when the mirror's database is down or the reply is cut."""
    url = URL + "?" + urllib.parse.urlencode(params, doseq=True)
    try:
        with urllib.request.urlopen(url, timeout=900) as r:
            raw = r.read()
    except Exception:  # noqa: BLE001
        return None
    txt = raw.decode("utf-8", "replace")
    # only this message means the database is down; "Table or Catalog not found" replies also carry STDB_ERR / Postgres
    # syntax errors, and mean the catalogue is absent
    if "not currently reachable" in txt or not txt.rstrip().endswith("</VOTABLE>"):
        return None
    v = parse(io.BytesIO(raw), verify="ignore")
    v.raw_text = txt
    return v


def patient(params, what, ok=lambda v: True):
    """Retry until the reply is complete and passes `ok` (the mirror also returns complete-looking replies with the data
    missing: on 2026-11-22 one query for EPIC 202071902 came back empty, the same query minutes later had 4 matches)."""
    for k in range(30):
        v = get(params)
        if v is not None and ok(v):
            return v
        time.sleep(min(300, 60 * (k + 1)))
    raise RuntimeError(f"VizieR failed: {what}")


def has_rows(v, prefix):
    return any(t.name and t.name.startswith(prefix) and len(t.to_table()) for t, _ in tables(v))


def tables(v):
    for res in v.iter_resources() if hasattr(v, "iter_resources") else v.resources:
        for t in res.tables:
            yield t, str(res.description or t.description or "")


def check_catalogues(ids):
    f = OUT / "vizier_catalogues.csv"
    done = pd.read_csv(f) if f.exists() else pd.DataFrame(columns=["catalog", "exists", "status", "tables"])
    rows = done.to_dict("records")
    for c in ids:
        if c in set(done.catalog):
            continue
        # three outcomes: tables listed (queried); "Catalog not found" (absent); a clean reply with empty RESOURCE elements
        # and no tables (e.g. J/A+AS/100/243: a catalogue without positions, which a cone search cannot query)
        v = patient({"-source": c, "-c": "0 +0", "-c.rs": 1, "-out.max": 1}, f"catalogue {c}")
        names = [t.name for t, _ in tables(v) if t.name]
        status = ("exists" if names else "not found" if "Catalog not found" in v.raw_text
                  else "no tables in a cone search")
        rows.append(dict(catalog=c, exists=bool(names), status=status, tables=" ".join(names)))
        print(f"catalogue {c}: {status}{', ' + str(len(names)) + ' tables' if names else ''}", flush=True)
        pd.DataFrame(rows).to_csv(f, index=False)
    d = pd.DataFrame(rows)
    return set(d.catalog[d.exists.astype(bool)])


def one(epic, tic, ra, dec, cats):
    f = CACHE / f"{epic}.json"
    if f.exists():
        return json.loads(f.read_text())
    # Gaia DR3 and 2MASS ride along as controls: every star is in them, so a reply with neither row is a failed one
    v = patient({"-source": " ".join(sorted(cats | {GAIA, TMASS})), "-c": f"{ra:.6f} {dec:+.6f}", "-c.rs": 3, "-out.max": 5,
                 "-out.all": 1, "-out.add": "_r"}, f"EPIC {epic}", ok=lambda v: has_rows(v, GAIA) or has_rows(v, TMASS))
    rows = []
    for t, desc in tables(v):
        if t.name in (GAIA, TMASS):
            continue
        tab = t.to_table(use_names_over_ids=True)
        if not len(tab):
            continue
        pc = [c for c in tab.colnames if c.lower() in PCOLS or c.lower().startswith("per")]
        rows.append(dict(EPIC=epic, TIC=tic, catalog=t.name, title=desc[:120], n_rows=len(tab),
                         sep_arcsec=float(tab["_r"][0]) if "_r" in tab.colnames else None,
                         periods="; ".join(f"{c}={tab[c][0]}" for c in pc)))
    f.write_text(json.dumps(rows, default=str))
    return rows


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv").drop_duplicates("EPIC")
    t1 = pd.read_csv(OUT / "table1_hits.csv")
    t1cats = {e: {c for c in (bib_to_vizier(b) for b in g.bibcode) if c} for e, g in t1.groupby("EPIC")}
    ids = set(CATALOGUES) | set().union(*t1cats.values())
    ok = check_catalogues(sorted(ids))
    print(f"{len(ok)} of {len(ids)} catalogues exist on the mirror", flush=True)
    allrows = []
    for k, r in enumerate(d.itertuples(), 1):
        cats = (set(CATALOGUES) | t1cats.get(int(r.EPIC), set())) & ok
        rows = one(int(r.EPIC), int(r.TIC), float(r.ra), float(r.dec), cats)
        allrows += rows
        print(f"[{k}/{len(d)}] EPIC {r.EPIC}: {len(rows)} tables, {sum(bool(x['periods']) for x in rows)} with a period column",
              flush=True)
    pd.DataFrame(allrows).to_csv(OUT / "vizier_targeted.csv", index=False)


if __name__ == "__main__":
    main()
