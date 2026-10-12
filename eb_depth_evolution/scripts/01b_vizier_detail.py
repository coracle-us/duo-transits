import warnings; warnings.filterwarnings("ignore")
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
c = SkyCoord(97.98012, 25.48210, unit="deg")
v = Vizier(columns=["**"], row_limit=50)
for cat in ["J/MNRAS/452/3561/table1","J/A+A/691/A242/obaf-eb1","J/A+A/691/A242/tessobaf","J/AJ/169/266/obs","J/A+A/684/A29/main","B/vsx/vsx","I/355/gaiadr3","J/ApJS/256/14/table4","V/162/dr11mem","J/MNRAS/506/2269/shift","II/271A/patch2","V/162/dr11m"]:
    try:
        r = v.query_region(c, radius=10*u.arcsec, catalog=cat)
        for k in r.keys():
            print("=====", k, r[k].meta.get('description',''))
            t=r[k]
            for row in t:
                for cn in t.colnames:
                    print(f"  {cn} = {row[cn]}")
                print("  --")
    except Exception as e: print(cat, "ERR", e)
