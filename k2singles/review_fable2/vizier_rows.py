"""Pull the specific catalogue rows that matter (full columns) for the review targets."""
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u, warnings, sys
warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
T = {"EPIC 211892898": (130.894296, 18.622201), "EPIC 211390677": (127.703812, 11.49457), "EPIC 211411112": (132.612448, 11.808697),
     "EPIC 202072917": (99.78832, 16.82623), "EPIC 202071902": (101.84677, 21.1342), "EPIC 202126877": (97.98012, 25.4821)}
want = {"EPIC 211892898": ["J/ApJS/249/22", "J/ApJ/926/191", "J/A+A/677/A21", "J/AJ/167/50", "J/AJ/154/224/table2", "J/MNRAS/442/2081"],
        "EPIC 211390677": ["J/MNRAS/522/29", "J/ApJ/858/L7", "J/ApJS/281/52", "J/A+A/705/A247", "J/ApJS/264/17", "J/ApJS/272/40", "J/ApJ/926/191", "J/A+A/677/A21", "J/ApJS/246/9", "J/A+A/690/A323"],
        "EPIC 211411112": ["J/AJ/161/190", "J/AJ/150/97", "J/MNRAS/463/1831", "J/ApJS/281/52", "J/ApJS/280/15"],
        "EPIC 202072917": ["J/MNRAS/452/3561", "J/A+A/579/A19", "J/MNRAS/529/1802", "J/AJ/170/287"],
        "EPIC 202071902": ["J/MNRAS/452/3561", "J/A+A/579/A19", "J/MNRAS/495/1252/table6", "J/MNRAS/527/521"],
        "EPIC 202126877": ["J/MNRAS/452/3561", "J/A+A/579/A19", "J/MNRAS/486/251"]}
for n, cats in want.items():
    ra, dec = T[n]
    print(f"\n==================== {n}")
    for cat in cats:
        v = Vizier(columns=["**"], row_limit=10, timeout=300, catalog=cat)
        try:
            res = v.query_region(SkyCoord(ra, dec, unit="deg"), radius=5 * u.arcsec)
        except Exception as e:
            print(f"  {cat}: ERR {e}"); continue
        for tab in res:
            print(f"  --- {tab.meta.get('ID')}  ({tab.meta.get('description','')[:80]})")
            for row in tab:
                items = [f"{c}={row[c]}" for c in tab.colnames if str(row[c]).strip() not in ("", "--", "nan", "None")]
                print("     ", "; ".join(items)[:1500])
