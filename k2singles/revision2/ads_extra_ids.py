"""
AAS82225 revision: ADS full-text search under the Gaia and 2MASS identifiers of the stars SIMBAD does not know.
literature.py took identifiers from SIMBAD; for stars with no SIMBAD entry (n_identifiers = 2 in per_star.csv) it
searched only "EPIC n" and "TIC n". This adds, for those stars: the Gaia DR3 source_id (gaia_ids.csv, position match),
and the Gaia DR2 id and 2MASS designation from the TIC (MAST Catalogs, TIC v8.2), as phrases and bare numbers.
Output: results/literature/ads_extra_ids.csv (one row per star; new bibcodes not already in hits.csv)
Usage:  python k2singles/revision2/ads_extra_ids.py
"""
import os
import sys

import pandas as pd
from astroquery.mast import Catalogs

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import literature as L  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
LIT = os.path.join(HERE, "data", "literature")
tok = open(os.path.expanduser("~/.ads/dev_key")).read().strip()
ps = pd.read_csv(os.path.join(LIT, "per_star.csv"))
hits = pd.read_csv(os.path.join(LIT, "hits.csv"))
gid = pd.read_csv(os.path.join(LIT, "gaia_ids.csv")).set_index("EPIC")
todo = ps[ps.n_identifiers <= 2]
tic = Catalogs.query_criteria(catalog="Tic", ID=[int(x) for x in todo.TIC]).to_pandas()
tic["ID"] = tic["ID"].astype(int)
tic = tic.set_index("ID")
rows = []
for r in todo.itertuples():
    phrases = []
    g3 = gid.source_id.get(r.EPIC)
    if pd.notna(g3):
        phrases += [f"Gaia DR3 {int(g3)}", str(int(g3))]
    if r.TIC in tic.index:
        t = tic.loc[r.TIC]
        if pd.notna(t.get("GAIA")) and str(t["GAIA"]).strip():
            phrases += [f"Gaia DR2 {t['GAIA']}", str(t["GAIA"])]
        if pd.notna(t.get("TWOMASS")) and str(t["TWOMASS"]).strip():
            phrases += [f"2MASS J{t['TWOMASS']}", f"J{t['TWOMASS']}"]
    phrases = sorted({L.ads_phrase(p) for p in phrases if len(p) >= 6})
    if not phrases:
        rows.append(dict(EPIC=r.EPIC, TIC=r.TIC, phrases="", numfound=0, new_bibcodes=""))
        continue
    q = "full:(" + " OR ".join(f'"{p}"' for p in phrases) + ")"
    docs, _, n = L.ads(q, tok)
    old = set(hits.bibcode[hits.EPIC == r.EPIC])
    new = [f"{d['bibcode']} ({(d.get('title') or [''])[0][:80]})" for d in docs if d["bibcode"] not in old]
    rows.append(dict(EPIC=r.EPIC, TIC=r.TIC, phrases=" | ".join(phrases), numfound=n, new_bibcodes=" ;; ".join(new)))
    print(r.EPIC, n, len(new), flush=True)
o = pd.DataFrame(rows)
o.to_csv(os.path.join(HERE, "results", "literature", "ads_extra_ids.csv"), index=False)
print(f"{len(o)} stars searched; {int((o.numfound > 0).sum())} with any ADS hit; "
      f"{int((o.new_bibcodes.fillna('') != '').sum())} with a paper not already found")
