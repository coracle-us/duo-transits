"""Gaia DR3 sources within 90 arcsec, and the target's astrophysical parameters (FLAME)."""
import os, json
from astroquery.gaia import Gaia
HERE = os.path.dirname(os.path.abspath(__file__))
q = """SELECT s.source_id, s.ra, s.dec, s.phot_g_mean_mag, s.bp_rp, s.parallax, s.ruwe,
 DISTANCE(POINT('ICRS', s.ra, s.dec), POINT('ICRS', 132.652546, 19.091365))*3600 AS sep,
 s.non_single_star, a.radius_flame, a.mass_flame, a.teff_gspphot, a.logg_gspphot, a.age_flame, a.lum_flame
 FROM gaiadr3.gaia_source AS s LEFT JOIN gaiadr3.astrophysical_parameters AS a ON s.source_id=a.source_id
 WHERE 1=CONTAINS(POINT('ICRS', s.ra, s.dec), CIRCLE('ICRS', 132.652546, 19.091365, 0.025)) ORDER BY sep"""
r = Gaia.launch_job(q).get_results()
r.pprint(max_width=250, max_lines=60)
r.write(os.path.join(HERE, "gaia_neighbours.csv"), overwrite=True)
