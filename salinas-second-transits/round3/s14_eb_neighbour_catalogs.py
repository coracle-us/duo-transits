"""Is the G = 12.74 eclipsing neighbour of TIC 107113048 (Gaia DR3 3463196954429755904, 97" away)
already known? Checks Gaia DR3 vari_eclipsing_binary, the TESS EB catalogue (Prsa+ 2022,
J/ApJS/258/16), VSX (B/vsx/vsx) and the TIC entry, all through VizieR (CfA mirror).

Usage: python round3/s14_eb_neighbour_catalogs.py   (prints; no files written)
"""
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier

SRC = 3463196954429755904
r = Gaia.launch_job(f"select source_id, ra, dec, phot_g_mean_mag, phot_variable_flag, ruwe, "
                    f"non_single_star from gaiadr3.gaia_source where source_id={SRC}").get_results()
print(r)
c = SkyCoord(r["ra"][0], r["dec"][0], unit="deg")
eb = Gaia.launch_job(f"select source_id, frequency, global_ranking from gaiadr3.vari_eclipsing_binary "
                     f"where source_id={SRC}").get_results()
print("Gaia DR3 vari_eclipsing_binary rows:", len(eb))
if len(eb):
    print(eb, " -> period", 1 / eb["frequency"][0], "d")
v = Vizier(columns=["**"], row_limit=20, vizier_server="vizier.cfa.harvard.edu")
for cat in ["J/ApJS/258/16", "B/vsx/vsx", "IV/39/tic82"]:
    res = v.query_region(c, radius=10 * u.arcsec, catalog=cat)
    print(cat, "->", [len(t) for t in res])
    for t in res:
        cols = [x for x in t.colnames if x.lower() in ("tic", "name", "type", "period", "per", "tmag", "_r")]
        print(t[cols] if cols else t)
