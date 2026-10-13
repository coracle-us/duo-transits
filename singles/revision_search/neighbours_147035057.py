"""
Can a neighbour produce the 5.3-5.6% primary eclipse of TIC 147035057? Lists Gaia DR3 sources within 60"
and their flux relative to the target (G band). A neighbour supplying f of the total flux can produce at most
an f-deep eclipse, so the deep eclipse must be on the target if the neighbours' total is well below 5%.
Output: stdout (quoted in singles/rnaas/numbers.md)
"""
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.mast import Catalogs

TIC = 147035057
RADIUS = 60 * u.arcsec


def main():
    t = Catalogs.query_criteria(catalog="TIC", ID=TIC)[0]
    c = SkyCoord(t["ra"], t["dec"], unit="deg")
    Gaia.ROW_LIMIT = -1
    r = Gaia.cone_search_async(c, radius=RADIUS).get_results()
    r.sort("phot_g_mean_mag")
    g0 = r["phot_g_mean_mag"][0]
    print(f"TIC {TIC}: target G = {g0:.2f}; Gaia DR3 sources within {RADIUS:.0f}: {len(r)}")
    total = 0.0
    for row in r[1:]:
        if row["phot_g_mean_mag"] is None or not float(row["phot_g_mean_mag"]) > 0:
            continue
        f = 10 ** (-0.4 * (row["phot_g_mean_mag"] - g0))
        total += f
        sep = row["dist"] * 3600
        if f > 1e-3:
            print(f"  G = {row['phot_g_mean_mag']:.2f} at {sep:5.1f}\"  flux ratio {f * 100:.2f}%")
    print(f"Total neighbour flux / target flux = {total * 100:.2f}%")


if __name__ == "__main__":
    main()
