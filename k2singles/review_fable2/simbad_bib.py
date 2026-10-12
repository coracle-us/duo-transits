"""SIMBAD bibliography for the review targets."""
from astroquery.simbad import Simbad
names = ["EPIC 211892898", "EPIC 211390677", "EPIC 211411112", "EPIC 202072917", "EPIC 202071902", "EPIC 202126877"]
for n in names:
    try:
        t = Simbad.query_bibobj(n).to_pandas()
        ids = Simbad.query_objectids(n).to_pandas()
    except Exception as e:
        print(n, "ERR", e); continue
    print(f"\n=== {n}: {len(t)} refs")
    for _, r in t.iterrows():
        print(f"  {r.get('bibcode')}  {str(r.get('title'))[:120]}")
    print("  ids:", "; ".join(ids.iloc[:, 0].astype(str)))
