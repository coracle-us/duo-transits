"""Novelty check: SIMBAD identifiers/bibliography and VizieR cone search around EPIC 202126877."""
import warnings; warnings.filterwarnings("ignore")
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
c = SkyCoord(97.98012, 25.48210, unit="deg")
print("=== SIMBAD ===")
try:
    r = Simbad.query_region(c, radius=10*u.arcsec)
    print(r)
    if r is not None and len(r):
        name = r[0]["main_id"] if "main_id" in r.colnames else r[0]["MAIN_ID"]
        print("IDs:", Simbad.query_objectids(name))
        try:
            b = Simbad.query_bibcode(name) if hasattr(Simbad,'query_bibcode') else None
        except Exception as e: b=None
        res = Simbad.query_tap(f"SELECT bibcode FROM has_ref JOIN ident ON has_ref.oidref=ident.oidref JOIN ref ON ref.oidbib=has_ref.oidbibref WHERE ident.id='{name}'")
        print("Bibcodes:", list(res["bibcode"]))
except Exception as e:
    print("SIMBAD error", e)
print("=== VizieR cone 10 arcsec ===")
v = Vizier(columns=["*"], row_limit=50)
try:
    res = v.query_region(c, radius=10*u.arcsec, catalog=None)
    for k in res.keys():
        t = res[k]
        print(f"--- {k} ({len(t)} rows)  {t.meta.get('description','')}")
        print(t[:3])
except Exception as e:
    print("VizieR error", e)
