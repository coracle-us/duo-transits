"""VizieR all-catalogue cone search (5 arcsec) for the review targets: which published catalogues contain them?"""
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u, warnings
warnings.filterwarnings("ignore")
T = {"EPIC 211892898": (130.894296, 18.622201), "EPIC 211390677": (127.703812, 11.49457), "EPIC 211411112": (132.612448, 11.808697),
     "EPIC 202072917": (99.78832, 16.82623), "EPIC 202071902": (101.84677, 21.1342), "EPIC 202126877": (97.98012, 25.4821)}
v = Vizier(columns=["**"], row_limit=5, timeout=600)
for n, (ra, dec) in T.items():
    print(f"\n==================== {n}")
    try:
        res = v.query_region(SkyCoord(ra, dec, unit="deg"), radius=5 * u.arcsec)
    except Exception as e:
        print("ERR", e); continue
    for tab in res:
        cat = tab.meta.get("ID", tab.meta.get("name", "?"))
        desc = tab.meta.get("description", "")[:90]
        print(f"  {cat:28s} n={len(tab):2d}  {desc}")
        # print period-like columns if any
        pcols = [c for c in tab.colnames if c.lower() in ("per", "period", "p", "porb", "pvar", "prot", "f", "freq")]
        for c in pcols:
            print(f"      {c} = {list(tab[c])}")
